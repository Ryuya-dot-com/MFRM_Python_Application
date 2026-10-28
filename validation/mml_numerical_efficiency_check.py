"""Compare frozen, first-improved and current MML on retained development inputs.

This is an efficiency/equivalence experiment, never a population qualification
or permission to replace frozen study code. Every fixed case and failure stays.
"""
import argparse
from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import shutil
import statistics
import zipfile

import numpy as np
import pandas as pd

from validation import mml_grouped_fit as current_fit
from validation.mml_grouped_likelihood import GroupedMmlLikelihood
from validation.mml_grouped_fit_efficiency_check import LIMITS, compare, digest, timed_fit
from validation.mml_gpcm_joint_check import prepare_problem
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data
from validation.mml_population_study import json_bytes


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec); spec.loader.exec_module(value)
    return value


def paired_comparison(a, b, persons):
    result = compare(a, b, persons)
    old, new = a['result'], b['result']; checks = result['checks']
    # An identical prefit rejection is valid agreement, not a successful fit.
    checks['same_numerical_pass'] = old['numerical_pass'] == new['numerical_pass']
    del checks['both_numerical_pass']
    checks['same_recession_audit'] = old.get('recession') == new.get('recession')
    checks['same_error_type'] = old.get('error_type') == new.get('error_type')
    if old.get('point') is None and new.get('point') is None:
        for key in LIMITS: checks.pop(key, None)
    if old.get('covariance') is None and new.get('covariance') is None:
        checks.pop('covariance_scaled', None)
    result['both_numerical_pass'] = bool(old['numerical_pass'] and new['numerical_pass'])
    result['all_checks_pass'] = all(checks.values())
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('parent', type=Path); parser.add_argument('output', type=Path)
    args = parser.parse_args(); root = Path(__file__).resolve().parents[1]
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    old_fit = module(root/'validation/grouped_fit_before.py', 'prior_grouped_fit')
    versions = {
        'frozen': (module(root/'validation/mml_grouped_likelihood_baseline.py', 'frozen_grouped').GroupedMmlLikelihood, old_fit),
        'previous': (module(root/'validation/grouped_before.py', 'previous_grouped').GroupedMmlLikelihood, old_fit),
        'candidate': (GroupedMmlLikelihood, current_fit),
    }
    names = [r['id'] for r in json.loads((args.parent/'summary.json').read_bytes())['records']]
    # Retain the actual import tree AND runtime resources, not just Python.
    sources = {str(p.relative_to(root)):digest(p) for p in root.rglob('*') if p.is_file()
               and '__pycache__' not in p.parts and (p.suffix in ['.py','.R','.json','.cff'])}
    with zipfile.ZipFile(out/'sources.zip', 'x', zipfile.ZIP_DEFLATED) as archive:
        for name in sources: archive.write(root/name, name)
    (out/'source_manifest.json').write_bytes(json_bytes(sources))
    inputs = {}
    for name in names:
        for suffix in ['.csv','.json']:
            source = args.parent/(name+suffix)
            inputs[name+suffix] = digest(source); shutil.copyfile(source,out/source.name)
        metadata = json.loads((out/(name+'.json')).read_bytes())
        if inputs[name+'.csv'] != metadata['input_sha256']: raise ValueError('Input changed')
    assert old_fit.OPTIONS == current_fit.OPTIONS and old_fit.LIMITS == current_fit.LIMITS
    plan = dict(classification='DEVELOPMENT_NUMERICAL_EFFICIENCY_ONLY', scientific_inference_ready=False,
        utc=datetime.now(timezone.utc).isoformat(), cases=names, input_sha256=inputs,
        source_manifest_sha256=digest(out/'source_manifest.json'), repeats=3,
        point_comparison=dict(atol=1e-9,rtol=1e-12,points=['zero','nonzero','stress']),
        fit_comparison_limits=LIMITS, solver_options=current_fit.OPTIONS, solver_limits=current_fit.LIMITS,
        rejection_rule='Same numerical status, solver checks, recession audit and error type; compare all available points/moments; covariance only when returned by both. Held fits never count as successful fits.',
        timing_scope='One serial worker; rotating version order; cold solver starts and mandatory polish; quadrature cache persists; all elapsed times included; two frozen workers may run concurrently',
        omitted_result_fields=['evaluated.posterior'])
    (out/'plan.json').write_bytes(json_bytes(plan)); records=[]
    for name in names:
        metadata = json.loads((out/(name+'.json')).read_bytes())
        case=metadata['case']; data=pd.read_csv(out/(name+'.csv')); q=case['q']
        engines={label:engine(data,metadata['config']) for label,(engine,_) in versions.items()}
        d=engines['frozen'].dimension
        points=[np.zeros(d), np.asarray(metadata['point']), np.r_[np.linspace(-3.,3.,d-1),np.log(4.)]]
        points[-1][engines['frozen'].linear_d:-1] = np.linspace(-1.,1.,engines['frozen'].slope_d)
        nodes,weights=current_fit.normal_rule(q); point_errors={}; point_checks=[]
        for index,point in enumerate(points):
            for info in [False,True]:
                evaluated={label:e.evaluate(point,nodes,weights,information=info) for label,e in engines.items()}
                for baseline in ['frozen','previous']:
                    for key in evaluated[baseline]:
                        old,new=np.asarray(evaluated[baseline][key]),np.asarray(evaluated['candidate'][key])
                        point_checks.append(bool(np.allclose(old,new,atol=plan['point_comparison']['atol'],
                                                            rtol=plan['point_comparison']['rtol'])))
                        error=float(np.max(abs(old-new)))
                        point_errors[key]=max(point_errors.get(key,0.),error)
        if case['model']=='GPCM': _,problem=prepare_problem(data,q,rating_max=case['categories']-1)
        else:
            _,problem=prepare_app_free_sd_problem_from_data(data,person_col='Person',facet_cols=['Rater','Task','Criterion'],
                score_col='Score',rating_min=0,rating_max=case['categories']-1,model=case['model'],noncenter_facet='Criterion',
                step_facet=None if case['model']=='RSM' else 'Criterion',quadrature_points=q)
        timings={label:[] for label in versions}; comparisons=[]; statuses={label:[] for label in versions}
        for repeat in range(plan['repeats']):
            order=list(versions);order=order[repeat:]+order[:repeat];fits={}
            for label in order:
                engine,fit_module=versions[label]
                fits[label]=timed_fit(engine,data,problem,q,fit_module=fit_module)
                timings[label].append(fits[label]['seconds']);statuses[label].append(fits[label]['result']['status'])
                (out/f'{name}_{repeat}_{label}.json').write_bytes(json_bytes(fits[label]))
            comparisons.extend(dict(baseline=label,repeat=repeat,**paired_comparison(fits[label],fits['candidate'],case['persons']))
                               for label in ['frozen','previous'])
        medians={label:statistics.median(values) for label,values in timings.items()}
        record=dict(id=name,case=case,point_checks_pass=all(point_checks),point_errors=point_errors,
            median_seconds=medians,statuses=statuses,comparisons=comparisons,
            speed_ratios={label:medians[label]/medians['candidate'] for label in ['frozen','previous']},
            all_checks_pass=all(point_checks) and all(v['all_checks_pass'] for v in comparisons))
        records.append(record);(out/(name+'_comparison.json')).write_bytes(json_bytes(record))
        print(name,record['all_checks_pass'],record['speed_ratios'],flush=True)
    summary=dict(classification=plan['classification'],scientific_inference_ready=False,records=records,
        all_checks_pass=all(r['all_checks_pass'] for r in records),completed_utc=datetime.now(timezone.utc).isoformat())
    (out/'summary.json').write_bytes(json_bytes(summary))
    return 0 if summary['all_checks_pass'] else 1


if __name__=='__main__': raise SystemExit(main())
