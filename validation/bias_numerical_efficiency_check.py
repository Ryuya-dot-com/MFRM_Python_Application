"""Paired complete G3 numerical pipelines on a separate development cohort.

One replicate of each declared design, same input in each implementation.
This checks implementation transport, not error rates or statistical power.
Source/runtime checks in the formal controller are excluded from timing;
generation, both Q fits, four screens and durable ledger writes are included.
"""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import zipfile

import numpy as np

from validation import bias_population_study as study, mml_population_study as shared
from validation import mml_grouped_fit as current_fit
from validation.mml_grouped_likelihood import GroupedMmlLikelihood
from validation.bias_population_protocol import make_registered_bias_plan, planned_bias_runs
from validation.mml_grouped_fit_efficiency_check import LIMITS, digest
from validation.mml_numerical_efficiency_check import module, paired_comparison
from validation.mml_population_study import AttemptLedger, json_bytes


def compare_screen(old, new, path='', errors=None):
    """Exact decisions/missingness/structure; bounded numerical discrepancies."""
    if errors is None: errors=[]
    if isinstance(old, dict) and isinstance(new, dict):
        if set(old)!=set(new): errors.append(dict(path=path, reason='keys_differ'))
        for key in old.keys() & new.keys(): compare_screen(old[key],new[key],path+'/'+key,errors)
    elif isinstance(old, list) and isinstance(new, list):
        if len(old)!=len(new): errors.append(dict(path=path, reason='length_differs'))
        for i,(a,b) in enumerate(zip(old,new)): compare_screen(a,b,path+'/'+str(i),errors)
    elif type(old) is float and type(new) is float:
        if not np.isclose(old,new,atol=2e-5,rtol=1e-7,equal_nan=True):
            errors.append(dict(path=path,old=old,new=new,reason='numeric_tolerance'))
    elif old!=new: errors.append(dict(path=path,old=old,new=new,reason='value_differs'))
    return errors


def numeric_fit(value):
    value=deepcopy(value)
    for key in ['point','covariance']:
        if value.get(key) is not None: value[key]=np.asarray(value[key],dtype=float)
    if 'evaluated' in value:
        for key in ['eap','sd']:
            value['evaluated'][key]=np.asarray(value['evaluated'][key],dtype=float)
    return dict(result=value)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    root=Path(__file__).resolve().parents[1]
    old_fit=module(root/'validation/grouped_fit_before.py','g3_prior_fit')
    old_engine=module(root/'validation/mml_grouped_likelihood_baseline.py','g3_prior_engine').GroupedMmlLikelihood
    versions=dict(frozen=(old_engine,old_fit),candidate=(GroupedMmlLikelihood,current_fit))
    original=make_registered_bias_plan()
    plan=dict(classification='DEVELOPMENT_ONLY',qualification_execution_allowed=False,
        study_id='DEVELOPMENT_ONLY_g3_efficiency_20261006',master_seed=2**117+2026100603,
        estimator=original['estimator'],screen=original['screen'],
        scenarios=[dict(s,replicates=1,independent_reference_replicates=[]) for s in original['scenarios']],
        scientific_inference_ready=False,scope=__doc__,fit_limits=LIMITS,
        screen_comparison=dict(atol=2e-5,rtol=1e-7,decisions_and_missingness='exact'),
        utc=datetime.now(timezone.utc).isoformat())
    sources={str(p.relative_to(root)):digest(p) for p in root.rglob('*') if p.is_file()
             and '__pycache__' not in p.parts and p.suffix in ['.py','.R','.json','.cff']}
    with zipfile.ZipFile(out/'sources.zip','x',zipfile.ZIP_DEFLATED) as z:
        for name in sources:z.write(root/name,name)
    (out/'source_manifest.json').write_bytes(json_bytes(sources))
    plan['source_manifest_sha256']=digest(out/'source_manifest.json')
    (out/'plan.json').write_bytes(json_bytes(plan));identity=digest(out/'plan.json')
    ledgers={label:AttemptLedger(out/(label+'.sqlite'),plan_sha256=identity) for label in versions}
    records=[]
    try:
        for index,run in enumerate(planned_bias_runs(plan)):
            scenario=next(s for s in plan['scenarios'] if s['id']==run['scenario_id'])
            order=list(versions) if index%2==0 else list(reversed(versions));timings={}
            for label in order:
                engine,solver=versions[label];solver.GroupedMmlLikelihood=engine
                study.fit_grouped_joint=solver.fit_grouped_joint
                study.GroupedMmlLikelihood=shared.GroupedMmlLikelihood=engine
                study.normal_rule=shared.normal_rule=solver.normal_rule
                start=time.perf_counter()
                study.execute_bias_run(ledgers[label],run,scenario,plan,verify_environment=lambda:None)
                timings[label]=time.perf_counter()-start
            key=run['run_id'];a,b=[ledgers[label] for label in versions]
            failures=[];fits=[]
            old_input,new_input=[l.get(key,'input_retained') for l in [a,b]]
            input_match=old_input==new_input
            for q in plan['estimator']['quadrature']:
                old,new=[l.get(key,f'Q{q}_finished') for l in [a,b]]
                if old is not None and new is not None:
                    fits.append(dict(q=q,**paired_comparison(numeric_fit(old),numeric_fit(new),scenario['design']['persons'])))
                else:
                    old,new=[l.get(key,f'Q{q}_failed') for l in [a,b]]
                    failures.extend(compare_screen(old,new,f'Q{q}_failed'))
            for phase in ['known_posterior','calibration_assessment','terminal',
                          'known_calibration_known_person','known_calibration_EAP',
                          'fitted_calibration_known_person','fitted_calibration_EAP']:
                failures.extend(compare_screen(a.get(key,phase),b.get(key,phase),phase))
            record=dict(run_id=key,input_match=input_match,fits=fits,screen_differences=failures,
                seconds=timings,time_ratio=timings['frozen']/timings['candidate'],
                all_checks_pass=input_match and not failures and all(f['all_checks_pass'] for f in fits))
            records.append(record);(out/(scenario['id']+'.json')).write_bytes(json_bytes(record))
            print(key,record['all_checks_pass'],round(record['time_ratio'],3),flush=True)
        summary=dict(classification=plan['classification'],scientific_inference_ready=False,
            all_checks_pass=all(r['all_checks_pass'] for r in records),records=records,
            ledger_chains={label:l.verify() for label,l in ledgers.items()})
        (out/'summary.json').write_bytes(json_bytes(summary))
        return 0 if summary['all_checks_pass'] else 1
    finally:
        for ledger in ledgers.values():ledger.close()


if __name__=='__main__':raise SystemExit(main())
