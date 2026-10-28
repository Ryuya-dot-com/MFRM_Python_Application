"""Deterministic conditional-slice controls, not simulated error-rate samples."""
from dataclasses import asdict

import numpy as np
import pandas as pd
from scipy.special import logsumexp,softmax

import streamlit_app as app
from validation.bias_population_design import prepare_bias_problem,screen_from_point
from validation.mml_population_protocol import design
from validation.mml_population_study import truth_joint


def scalar_controls():
    """Known closed-form roots, nonexistence, search bounds and weight scope."""
    base=pd.DataFrame(dict(eta=[.2]*10,score=[0]*3+[1]*7,slope=[1.]*10,weight=[1.]*10,step1=[0.]*10))
    cases=[]
    def add(name,data,status='interior',bound=10.,**expected):
        cases.append(dict(id=name,data=data,bound=bound,expected=dict(status=status,**expected)))
    for a in [.5,1.,2.]:
        add('binary_slope_'+str(a),base.assign(slope=a),estimate=np.log(7/3)/a-.2,
            se=1/np.sqrt(10*a*a*.7*.3),p_available=True,ci_status='ok')
    for y in [0,1]:add('all_'+str(y),base.assign(score=y),status='no_finite_mle',p_available=False)
    add('root_outside_range',base,status='search_range_limited',bound=.1,p_available=False)
    add('fractional_weights',base.assign(weight=.5),estimate=np.log(7/3)-.2,
        se=1/np.sqrt(5*.7*.3),p_available=False,ci_status='not available')
    add('zero_weight_row',pd.concat([base,base.iloc[:1].assign(score=1,weight=0.)],ignore_index=True),
        estimate=np.log(7/3)-.2,se=1/np.sqrt(10*.7*.3),p_available=True,ci_status='ok')
    add('zero_support',base.assign(weight=0.),status='invalid_input',p_available=False)
    bad=base.copy();bad.loc[0,'weight']=-1.
    add('negative_weight',bad,status='invalid_input',p_available=False)
    add('zero_slope',base.assign(slope=0.),status='invalid_input',p_available=False)
    add('fractional_score',base.assign(score=.5),status='invalid_input',p_available=False)
    add('zero_search_bound',base,status='invalid_input',bound=0.,p_available=False)
    add('interval_upper_outside_range',base,bound=.9,p_available=True,ci_status='ok; limited by search range')
    add('interval_lower_outside_range',base.assign(eta=-.2,score=1-base.score),bound=.9,p_available=True,ci_status='limited by search range; ok')
    add('interval_both_outside_range',base.iloc[:2].assign(eta=0.,score=[0,1]),bound=.9,p_available=True,ci_status='limited by search range')
    # With fixed thresholds, missing internal observed categories do not imply
    # nonexistence of the one-dimensional shift MLE. At b=-eta, the symmetric
    # fixed distribution has mean 1.5 and the curvature below is exact.
    variance=(2.25+.25*np.exp(.7))/(1+np.exp(.7))
    four=base.assign(score=[0]*5+[3]*5,step1=-.7,step2=0.,step3=.7)
    add('fixed_steps_missing_middle_categories',four,estimate=-.2,se=1/np.sqrt(10*variance),p_available=True,ci_status='ok')
    return cases


def native_scalar(case):
    data=case['data'];y=data.score.to_numpy();a=data.slope.to_numpy();w=data.weight.to_numpy()
    step_cols=sorted([c for c in data if c.startswith('step')],key=lambda c:int(c[4:]))
    cumulative=np.c_[np.zeros(len(data)),np.cumsum(data[step_cols].to_numpy(),axis=1)]
    k=np.arange(len(step_cols)+1)
    def logits(b):return a[:,None]*((data.eta.to_numpy()+b)[:,None]*k-cumulative)
    def probability(b):return softmax(logits(b),axis=1)
    def nll(b):
        values=logits(b)
        return float(np.sum(w*(logsumexp(values,axis=1)-values[np.arange(len(y)),y.astype(int)])))
    fit=app._fit_conditional_bias_cell(nll,probability,y,a,w,case['bound'])
    expected=case['expected'];checks=dict(status=fit['status']==expected['status'],
        p_available=bool(np.isfinite(fit['p']))==expected['p_available'])
    for name in ['estimate','se']:
        if name in expected:checks[name]=bool(np.isfinite(fit[name]) and abs(fit[name]-expected[name])<=1e-10)
    if 'ci_status' in expected:checks['ci_status']=fit['ci']['status']==expected['ci_status']
    row={'Cell Status':fit['status'],'Reference Status':fit['reference_status'],'Bias Size':fit['estimate'],
        'S.E.':fit['se'],'LR ChiSq':fit['lr'],'LR Prob.':fit['p'],'Profile CI Lower':fit['ci']['lower'],
        'Profile CI Upper':fit['ci']['upper'],'Profile CI Status':fit['ci']['status']}
    return dict(id=case['id'],row=row,expected=expected,checks=checks,all_checks_pass=all(checks.values()))


def family_control(model):
    """Fixed response patterns distinguish empty cells, extreme levels and MLEs.

    The chosen responses are deliberate controls; no recovery/coverage or
    population operating-characteristic claim is made from these rows.
    """
    spec=design(model,'regular',persons=4);levels=dict(Person=[f'P{i}' for i in range(4)],
        Rater=['R0','R1','R2'],Task=['T0','T1'],Criterion=['C0','C1','C2'])
    rows=[]
    for p in range(4):
        for r in range(3):
            for t in range(2):
                if (r,t)==(2,1):continue
                for c in range(3):
                    y=0 if r==0 else 3 if (r,t)==(1,0) else (p+c+t)%4
                    rows.append(dict(Person=f'P{p}',Rater=f'R{r}',Task=f'T{t}',Criterion=f'C{c}',Score=y))
    sample=dict(data=pd.DataFrame(rows),declared_levels=levels,metadata=dict(design=asdict(spec),classification='DETERMINISTIC_GUARD_CONTROL'),
        truth=dict(Person=np.array([-.5,-.1,.2,.4]),Rater=spec.rater,Task=spec.task,Criterion=spec.criterion,
            Step=spec.steps,LogSlope=spec.log_slopes,PopulationSD=spec.population_sd,RaterTaskEtaShift=np.zeros((3,2))))
    prep,problem=prepare_bias_problem(sample['data'],spec,31);point=truth_joint(problem.config,sample['truth'])
    result=screen_from_point(prep,problem,point,sample['truth']['Person'],sample,basis='known_calibration_known_person')
    expected=['excluded_extreme_level','excluded_extreme_level','no_finite_mle','interior','interior','unobserved_declared_cell']
    return dict(sample=sample,result=result,expected_statuses=expected,
        native_status_pass=[c['status'] for c in result['cells']]==expected,scientific_inference_ready=False)


def execute_controls(output,plan,identity,root,rscript,verifier):
    """One retained control attempt per group; no operating-rate evidence."""
    import json
    from pathlib import Path
    from validation.bias_population_check import compare_reference
    from validation.bias_population_references import _write,_digest,_family_input,_r_family,review_family
    from validation.mml_population_reference_execution import _once
    from validation.mml_population_study import AttemptLedger,json_bytes

    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True);root=Path(root)
    verifier();ledger=AttemptLedger(output/'events.sqlite',plan_sha256=identity)
    definition=plan['deterministic_controls'];cases=scalar_controls();models=definition['family_cases']
    if len(cases)!=definition['scalar_cases'] or models!=['RSM','PCM','GPCM']:
        ledger.close();raise ValueError('Required deterministic controls changed')

    def run_scalar():
        target=output/'scalar';target.mkdir(exist_ok=False)
        native=[native_scalar(c) for c in cases];request=[]
        for c in cases:
            filename=c['id']+'.csv';c['data'].to_csv(target/filename,index=False,float_format='%.17g')
            request.append(dict(id=c['id'],file=filename,bound=c['bound']))
        _write(target/'input.json',dict(cases=request,families=[],include_points=False))
        _write(target/'native.json',native)
        r=_r_family(root,rscript,target,target/'R')
        compared,_=compare_reference(native,[],r['result'],plan['independent_reference']['limits'])
        return dict(native=[dict(id=v['id'],pass_=v['all_checks_pass']) for v in native],references=compared,
            files={str(p.relative_to(output)):_digest(p) for p in target.rglob('*') if p.is_file()})

    def run_family(model):
        value=family_control(model);sample=value['sample'];retained={k:v for k,v in sample.items() if k!='data'}
        retained['csv_utf8']=sample['data'].to_csv(index=False)
        scenario=dict(design=retained['metadata']['design'],eta_shift=np.zeros((3,2)).tolist())
        target=output/model
        _family_input(target,value['result'],retained,scenario,'known_calibration_known_person',plan)
        _write(target/'native.json',dict(expected=value['expected_statuses'],
            actual=[c['status'] for c in value['result']['cells']],native_status_pass=value['native_status_pass']))
        r=_r_family(root,rscript,target,target/'R')
        # In-memory native rows contain NaNs; the file is retained independently.
        from validation.bias_population_references import _decode
        py=_decode(json.loads((target/'python.json').read_bytes()));comparison=review_family(py,r['result'],plan)
        return dict(native=[dict(id=model,pass_=value['native_status_pass'])],
            references=comparison['cells']+comparison['adjustments']+comparison['derived_adjustments'],
            files={str(p.relative_to(output)):_digest(p) for p in target.rglob('*') if p.is_file()})

    results=[]
    try:
        ledger.verify()
        groups=[('scalar',len(cases),len(cases),run_scalar)]
        groups += [(model,1,10,lambda m=model:run_family(m)) for model in models]
        for key,native_n,reference_n,work in groups:
            record=_once(ledger,key,'control',work,verifier);value=record.get('value',{})
            for name,sha in value.get('files',{}).items():
                if _digest(output/name)!=sha:raise ValueError('Retained control artifact changed')
            native=value.get('native',[]);refs=value.get('references',[])
            results.append(dict(id=key,status=record['status'],required_native=native_n,required_reference=reference_n,
                passed_native=sum(bool(r['pass_']) for r in native),passed_reference=sum(bool(r['pass_']) for r in refs),
                all_checks_pass=len(native)==native_n and len(refs)==reference_n and all(r['pass_'] for r in native+refs)))
        summary=dict(classification='DETERMINISTIC_FORMULA_AND_GUARD_CONTROLS',groups=results,
            required_native=sum(r['required_native'] for r in results),required_reference=sum(r['required_reference'] for r in results),
            passed_native=sum(r['passed_native'] for r in results),passed_reference=sum(r['passed_reference'] for r in results),
            all_checks_pass=all(r['all_checks_pass'] for r in results),ledger=ledger.verify(),scientific_inference_ready=False,
            interpretation='Fixed controls, not independent simulation replicates or LR-distribution qualification')
        target=output/'summary.json'
        if target.exists():
            if target.read_bytes()!=json_bytes(summary)+b'\n':raise ValueError('Completed control summary differs; no replacement')
        else:_write(target,summary)
        return summary
    finally:ledger.close()
