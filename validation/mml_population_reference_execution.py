"""Single-attempt independent-reference jobs, bound to a frozen study.

The sealed controller verifies sources/runtime and serializes this command
against population fitting. Numerical failures and interrupted reference jobs
remain unavailable; no replacement, reference-point reselection or refitting.
"""
import hashlib
import json
from pathlib import Path
import subprocess

from validation.mml_population_study import AttemptLedger,json_bytes
from validation.mml_population_reference_review import (
    continuous_values,continuous_review,finite_review,tam_review)


def _write(path,value):
    with path.open('xb') as handle:handle.write(json_bytes(value)+b'\n')


def _once(ledger,key,phase,work,verifier):
    finished=ledger.get(key,phase+'_finished')
    if finished is not None:return finished
    if ledger.get(key,phase+'_started') is not None:
        result=dict(status='interrupted_reference_retained',scientific_inference_ready=False)
    else:
        verifier();ledger.append(key,phase+'_started',dict(single_attempt=True))
        try:
            value=work();verifier()
            result=dict(status='completed',value=value,scientific_inference_ready=False)
        except Exception as exc:
            result=dict(status='reference_failed',error=dict(type=type(exc).__name__,message=str(exc)),scientific_inference_ready=False)
    ledger.append(key,phase+'_finished',result)
    print(key,phase,result['status'],flush=True)
    return result


def _r_one(root,rscript,output,case,*,settings=None):
    output.mkdir(parents=True,exist_ok=False)
    request=dict(classification='DEVELOPMENT_ONLY',cases=[case])
    if settings is not None:request['continuous_settings']=settings
    _write(output/'input.json',request)
    with (output/'stdout.txt').open('xb') as stdout,(output/'stderr.txt').open('xb') as stderr:
        result=subprocess.run([str(rscript),'--vanilla',str(root/'validation/mml_joint_information_reference.R'),
            str(output/'input.json'),str(output/'result.json')],cwd=root,stdout=stdout,stderr=stderr)
    _write(output/'execution.json',dict(returncode=result.returncode,
        files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir() if p.is_file()}))
    if result.returncode:raise RuntimeError('Independent R exited '+str(result.returncode))
    values=json.loads((output/'result.json').read_bytes())['cases']
    if len(values)!=1 or values[0]['id']!=case['id']:raise ValueError('R returned a different required case')
    return values[0]


def execute_bundle(bundle,plan,identity,root,rscript,verifier):
    bundle=Path(bundle);selection=json.loads((bundle/'reference_plan.json').read_bytes())
    if selection['study_sha256']!=identity:raise ValueError('Reference input belongs to another study')
    for name,sha in json.loads((bundle/'manifest.json').read_bytes()).items():
        if hashlib.sha256((bundle/name).read_bytes()).hexdigest()!=sha:raise ValueError('Reference input changed: '+name)
    output=bundle/'execution';output.mkdir(exist_ok=True)
    ledger=AttemptLedger(output/'events.sqlite',plan_sha256=identity)
    python={c['id']:c for c in json.loads((bundle/'python_finite.json').read_bytes())}
    native={c['id']:c for c in json.loads((bundle/'native_rowwise.json').read_bytes())}
    finite=json.loads((bundle/'finite_input.json').read_bytes())['cases']
    continuous=json.loads((bundle/'continuous_input.json').read_bytes())['cases']
    tam=json.loads((bundle/'tam_input/reference_input.json').read_bytes())['cases']
    by_lane=dict(finite_R={c['id']:c for c in finite},continuous={c['id']:c for c in continuous},TAM={c['id']:c for c in tam})
    reviews=[]
    try:
        ledger.verify()
        for required in selection['required_cases']:
            key,lane=required['id'],required['lane']
            if not required['available']:
                reviews.append(dict(id=key,lane=lane,status='required_input_or_point_unavailable',all_checks_pass=False));continue
            case=by_lane[lane][key]
            if lane in ['finite_R','continuous']:
                r=_once(ledger,key,lane+'_R',lambda:_r_one(root,rscript,output/(key+'_'+lane+'_R'),case,
                    settings=plan['numerical_reference']['continuous_settings'] if lane=='continuous' else None),verifier)
                if lane=='continuous':
                    value=_once(ledger,key,'continuous_python',lambda:continuous_values(case,plan),verifier)
                    work=lambda:continuous_review(case,value['value'],r['value']['continuous'],python.get(key),plan)
                    available=value['status']=='completed' and r['status']=='completed'
                else:
                    work=lambda:finite_review(case,python[key],native[key],r['value'],plan)
                    available=r['status']=='completed'
                review=_once(ledger,key,lane+'_review',work,verifier) if available else dict(status='reference_unavailable')
                reviews.append(dict(id=key,lane=lane,status=review['status'],all_checks_pass=review.get('value',{}).get('all_checks_pass',False)))
            else:
                def tam_work():
                    target=output/(key+'_TAM');log=output/(key+'_TAM_stdout.txt');error=output/(key+'_TAM_stderr.txt')
                    with log.open('xb') as stdout,error.open('xb') as stderr:
                        completed=subprocess.run([str(rscript),'--vanilla',str(root/'validation/mml_rasch_tam_refit.R'),
                            str(bundle/'tam_input'),str(target),'^'+key+'$'],cwd=root,stdout=stdout,stderr=stderr)
                    files={str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in [log,error,*target.rglob('*')] if p.is_file()}
                    return dict(returncode=completed.returncode,files=files)
                fitted=_once(ledger,key,'TAM_fit',tam_work,verifier)
                target=output/(key+'_TAM')
                for grid in ['b12_h10','b16_h05']:
                    grid_key=key+'_'+grid
                    def review_tam():
                        if fitted['status']!='completed':raise ValueError('TAM attempt unavailable')
                        for name,sha in fitted['value']['files'].items():
                            if hashlib.sha256((output/name).read_bytes()).hexdigest()!=sha:raise ValueError('TAM artifact changed')
                        value=json.loads((target/(grid_key+'.json')).read_bytes())
                        if 'error' in value:raise ValueError(value['error'])
                        design=json.loads((target/(key+'_design.json')).read_bytes())
                        return tam_review(case,design,value,python.get(key),case['point'],plan)
                    review=_once(ledger,grid_key,'TAM_review',review_tam,verifier)
                    reviews.append(dict(id=grid_key,lane=lane,status=review['status'],all_checks_pass=review.get('value',{}).get('all_checks_pass',False)))
        result=dict(run_id=selection['run_id'],required_cases=len(selection['required_cases']),comparisons=reviews,
            all_checks_pass=bool(reviews) and all(r['all_checks_pass'] for r in reviews),ledger=ledger.verify(),scientific_inference_ready=False)
        path=output/'summary.json'
        if path.exists():
            if path.read_bytes()!=json_bytes(result)+b'\n':raise ValueError('Completed reference summary differs; no replacement')
        else:_write(path,result)
        return result
    finally:ledger.close()
