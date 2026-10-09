"""Development transport controls: DGP -> joint fit -> actual bias screen -> R.

This small observed fixture is excluded from any prospective G3 cohort. Its
purpose is to verify targets and software transport before fixing a statistical
study. It does not estimate reliable error rates or qualify a reference law.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

import numpy as np
import pandas as pd
import streamlit_app as app
from validation.bias_population_design import contrast_matrix,generate_shifted_sample,screen_from_point,prepare_bias_problem
from validation.mfrm_recovery_design import replicate_seed
from validation.mml_population_protocol import design,make_plan
from validation.mml_population_study import AttemptLedger,json_bytes,truth_joint,compact_fit,compare_orders
from validation.mml_grouped_fit import fit_grouped_joint,normal_rule
from validation.mml_grouped_likelihood import GroupedMmlLikelihood
from validation.mml_population_controller import source_paths


def write(path,value):
    with path.open('xb') as handle:handle.write(json_bytes(value)+b'\n')


def compare_reference(python,families,reference,limits):
    comparisons=[];adjustments=[]
    if len({r['id'] for r in python})!=len(python):raise ValueError('Duplicate reference case identity')
    for py,r in zip(python,reference['cases'],strict=True):
        if py['id']!=r['id']:raise ValueError('Reference case identity changed')
        row=py['row'];status=row['Cell Status']==r['status'];errors={};missing=[]
        if r['status']=='interior':
            required=['estimate','se','lr']
            if row.get('Reference Status')=='unit_weight_conditional_lr_candidate':required+=['p','ci','ci_status']
            missing=[key for key in required if key not in r]
            for key,column in [('estimate','Bias Size'),('se','S.E.'),('lr','LR ChiSq'),('p','LR Prob.')]:
                if key in r:errors[key]=abs(row[column]-r[key])
            if 'ci' in r:
                if np.asarray(r['ci']).shape!=(2,):raise ValueError('Reference interval must have two endpoints')
                errors['ci']=float(np.max(abs(np.array([row['Profile CI Lower'],row['Profile CI Upper']])-r['ci'])))
                status &= row['Profile CI Status']==r.get('ci_status')
        comparisons.append(dict(id=py['id'],errors=errors,status_match=bool(status),missing_reference_fields=missing,
            pass_=bool(status and not missing and all(v<=limits[k] for k,v in errors.items()))))
    for family,r in zip(families,reference['adjustments'],strict=True):
        for method in ['holm','BH']:
            a=app.adjust_p_values(np.asarray(family,dtype=float),method);b=np.asarray(r[method],dtype=float)
            adjustments.append(dict(method=method,pass_=bool(a.shape==b.shape and np.allclose(a,b,rtol=0,atol=limits['adjustment'],equal_nan=True))))
    return comparisons,adjustments


def evaluate_case(out,ledger,key,spec,seed,magnitude,plan,reference_cases,python_cells,families):
    ledger.append(key,'generation_started',dict(seed=seed,design=asdict(spec),magnitude=magnitude))
    sample=generate_shifted_sample(spec,seed,contrast_matrix(len(spec.rater),len(spec.task),magnitude))
    raw=sample['data'].to_csv(index=False).encode();(out/(key+'.csv')).write_bytes(raw)
    retained={k:sample[k] for k in ['truth','declared_levels','metadata']}
    retained.update(csv_utf8=raw.decode(),csv_sha256=hashlib.sha256(raw).hexdigest())
    ledger.append(key,'input_retained',retained)
    prep,problem=prepare_bias_problem(sample['data'],spec,641)
    fits=[]
    for q in plan['estimator']['quadrature']:
        ledger.append(key,f'Q{q}_started',dict(input_sha256=retained['csv_sha256'],q=q))
        fit=compact_fit(fit_grouped_joint(sample['data'],problem,q=q));fits.append(fit)
        ledger.append(key,f'Q{q}_finished',fit)
    sensitivity=compare_orders(sample['data'],problem,fits,plan)
    point=truth_joint(problem.config,sample['truth']);grouped=GroupedMmlLikelihood(sample['data'],problem.config)
    known=[grouped.evaluate(point,*normal_rule(q)) for q in plan['estimator']['quadrature']]
    limits=plan['estimator']['quadrature_limits'];known_error=dict(nll_per_person=abs(known[0]['nll']-known[1]['nll'])/spec.persons,
        cross_score_supnorm=float(np.max(abs(known[0]['gradient']-known[1]['gradient']))),
        posterior_moment=max(float(np.max(abs(known[0][k]-known[1][k]))) for k in ['eap','sd']))
    known_ok=all(v<=limits[k] for k,v in known_error.items());eligible=sensitivity['pass_']
    lanes=[('known_calibration_known_person',point,sample['truth']['Person'],True),
        ('known_calibration_EAP',point,known[0]['eap'],known_ok),
        ('fitted_calibration_known_person',fits[0]['point'],sample['truth']['Person'],eligible),
        ('fitted_calibration_EAP',fits[0]['point'],fits[0].get('evaluated',{}).get('eap'),eligible)]
    statuses={}
    for basis,calibration,theta,ok in lanes:
        if not ok:
            statuses[basis]='numerical_calibration_unavailable';ledger.append(key,basis,dict(status=statuses[basis]));continue
        result=screen_from_point(prep,problem,calibration,theta,sample,basis=basis)
        ledger.append(key,basis,result);statuses[basis]='computed'
        reference=pd.DataFrame(result['reference_rows'])
        p=[]
        for row in result['table']:
            identifier=key+'_'+basis+'_'+row['FacetA_Level']+'_'+row['FacetB_Level'];filename=identifier+'.csv'
            data=reference.loc[reference.Rater.eq(row['FacetA_Level'])&reference.Task.eq(row['FacetB_Level'])]
            data.to_csv(out/filename,index=False,float_format='%.17g')
            reference_cases.append(dict(id=identifier,file=filename,bound=10.,excluded_extreme=row['Cell Status']=='excluded_extreme_level'))
            python_cells.append(dict(id=identifier,row=row))
            p.append(float(row['LR Prob.']) if np.isfinite(row['LR Prob.']) else None)
        families.append(p)
    result=dict(id=key,numerical_calibration_eligible=eligible,q_comparison=sensitivity,known_calibration_q_errors=known_error,
        lanes=statuses,scientific_inference_ready=False)
    ledger.append(key,'terminal',result);return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path);args=parser.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);plan=make_plan()
    study_id='DEVELOPMENT_ONLY_G3_TRANSPORT_20261004';master=923760110294
    cases=[]
    for model in ['RSM','PCM','GPCM']:
        for assignment in ['regular','sparse','unbalanced']:
            for magnitude in [0.,.8]:
                key=f'{model}_{assignment}_B{round(100*magnitude)}'
                cases.append(dict(id=key,model=model,assignment=assignment,magnitude=magnitude,
                    seed=replicate_seed(master,study_id+':'+key,1)))
    limits=dict(estimate=1e-8,se=1e-8,lr=1e-8,p=1e-8,ci=1e-7,adjustment=1e-14)
    protocol=dict(classification='DEVELOPMENT_ONLY',study_id=study_id,cases=cases,persons=48,categories=4,
        quadrature=[641,1281],numerical_contract=plan['estimator'],screen=dict(alpha=.05,min_n=5,practical_logit=.5,max_abs=10.,omit_extreme=True),
        reference_limits=limits,reference='Independent R formula/roots/CI and Holm/BH, every computed native screen cell',
        scientific_inference_ready=False,scope='Software/estimand transport control; one observed sample per condition, not an operating-characteristic qualification study')
    write(out/'protocol.json',protocol)
    with zipfile.ZipFile(out/'sources.zip','x',zipfile.ZIP_DEFLATED) as archive:
        for p in source_paths():archive.write(p,str(p.relative_to(ROOT)))
    identity=hashlib.sha256((out/'protocol.json').read_bytes()).hexdigest();ledger=AttemptLedger(out/'events.sqlite',plan_sha256=identity)
    records=[];references=[];python=[];families=[]
    for case in cases:
        try:
            spec=design(case['model'],case['assignment'],persons=48)
            record=evaluate_case(out,ledger,case['id'],spec,case['seed'],case['magnitude'],plan,references,python,families)
        except Exception as exc:
            record=dict(id=case['id'],status='execution_failed',error=dict(type=type(exc).__name__,message=str(exc)))
            ledger.append(case['id'],'terminal_failure',record)
        records.append(record);print(case['id'],record.get('lanes',record),flush=True)
    ledger_state=ledger.verify();ledger.close()
    write(out/'input.json',dict(cases=references,families=families));write(out/'python.json',python)
    with (out/'R.log').open('xb') as log:
        run=subprocess.run(['/usr/local/bin/Rscript','--vanilla',str(ROOT/'validation/bias_conditional_reference.R'),
            str(out/'input.json'),str(out/'R.json')],stdout=log,stderr=log)
    comparisons=[];adjustments=[]
    if run.returncode==0:
        reference=json.loads((out/'R.json').read_bytes())
        comparisons,adjustments=compare_reference(python,families,reference,limits)
    summary=dict(classification='DEVELOPMENT_ONLY',cases=records,planned_cases=len(cases),ledger=ledger_state,R_returncode=run.returncode,
        planned_screen_lanes=4*len(cases),declared_cell_slots=6*4*len(cases),
        required_reference_cells=len(python),comparisons=comparisons,adjustments=adjustments,
        all_transport_checks_pass=len(records)==len(cases) and all('lanes' in r and all(v=='computed' for v in r['lanes'].values()) for r in records)
            and len(comparisons)==len(python) and all(r['pass_'] for r in comparisons+adjustments),scientific_inference_ready=False)
    write(out/'summary.json',summary)
    write(out/'manifest.json',{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()})
    print(json.dumps(dict(cases=len(records),reference_cells=len(comparisons),all_transport_checks_pass=summary['all_transport_checks_pass'])),flush=True)


if __name__=='__main__':main()
