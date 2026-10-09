"""Prespecified G3 family references, retaining every required failure slot.

R receives full row designs and derives empty/extreme cells itself. Formula
agreement is conditional on those transported row designs, not proof of the
LR reference law or a replacement for calibration/recovery comparisons.
"""
from collections import Counter
import hashlib
from io import StringIO
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
import streamlit_app as app

from validation.bias_population_check import compare_reference
from validation.mfrm_recovery_design import replicate_seed
from validation.mml_population_study import AttemptLedger,json_bytes
from validation.mml_population_reference_execution import _once


def _write(path,value):
    with Path(path).open('xb') as handle:handle.write(json_bytes(value)+b'\n')


def _digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _decode(value):
    if isinstance(value,dict):
        if set(value)=={'nonfinite'}:
            if value['nonfinite'] not in ['nan','inf','-inf']:raise ValueError('Invalid nonfinite evidence tag')
            return float(value['nonfinite'])
        return {k:_decode(v) for k,v in value.items()}
    if isinstance(value,list):return [_decode(v) for v in value]
    return value


def _family_input(target,result,retained,scenario,basis,plan):
    levels=retained['declared_levels'];keys=[(r,t) for r in levels['Rater'] for t in levels['Task']]
    for facet,parameter,prefix in [('Rater','rater','R'),('Task','task','T'),('Criterion','criterion','C')]:
        if levels[facet]!=[prefix+str(i) for i in range(len(scenario['design'][parameter]))]:
            raise ValueError('Declared reference levels differ from scenario')
    if result.get('basis')!=basis or [(c['rater'],c['task']) for c in result['cells']]!=keys:
        raise ValueError('Reference lane or declared targets changed')
    if any(c['DGP_shift']!=b or c['oracle_target']!=(basis=='known_calibration_known_person')
        for c,b in zip(result['cells'],np.asarray(scenario['eta_shift']).ravel(),strict=True)):
        raise ValueError('Reference cell targets changed')
    original=pd.read_csv(StringIO(retained['csv_utf8']),dtype={k:str for k in levels})
    rows=pd.DataFrame(result['reference_rows'])
    # Verify the response and assignment transport before comparing formulae.
    for column in ['Person','Rater','Task','Criterion']:
        if not np.array_equal(rows[column].to_numpy(),original[column].to_numpy()):
            raise ValueError('Reference row assignment changed: '+column)
    if not np.array_equal(rows['score'].to_numpy(),original.Score.to_numpy()-scenario['design']['rating_min']):
        raise ValueError('Reference response rows changed')
    expected={(r['FacetA_Level'],r['FacetB_Level']):r for r in result['table']}
    if len(expected)!=len(result['table']) or set(expected)-set(keys):raise ValueError('Unexpected or duplicate native cells')
    requests=[];python=[];p=[];observed_ids=[];actual_holm=[]
    for (r,t),cell in zip(keys,result['cells'],strict=True):
        present=bool(((rows.Rater==r)&(rows.Task==t)).any())
        if present!=((r,t) in expected):raise ValueError('Native observed-cell family differs from reference rows')
        identifier=r+'_'+t;row=expected.get((r,t),{'Cell Status':'unobserved_declared_cell'})
        if row['Cell Status']!=cell['status']:raise ValueError('Cell/table status mismatch')
        requests.append(dict(id=identifier,rater=r,task=t,bound=plan['screen']['max_abs']))
        python.append(dict(id=identifier,row=row))
        if present:
            observed_ids.append(identifier)
            p.append(float(row['LR Prob.']) if np.isfinite(row['LR Prob.']) else None)
            actual_holm.append(cell.get('p_holm',np.nan))
    expected_holm=app.adjust_p_values(np.asarray(p,dtype=float),'holm')
    if not np.allclose(expected_holm,np.asarray(actual_holm,dtype=float),rtol=0,atol=1e-14,equal_nan=True):
        raise ValueError('Retained Holm values differ from the complete observed-cell family')
    target.mkdir(exist_ok=False);rows.to_csv(target/'family.csv',index=False,float_format='%.17g')
    request=dict(family_file='family.csv',cases=requests,families=[p],include_points=False,
        omit_extreme=plan['screen']['omit_extreme'])
    _write(target/'input.json',request)
    _write(target/'python.json',dict(cases=python,families=[p],observed_ids=observed_ids))


def export_reference(out,ledger,scenario,replicate,plan,identity):
    """No point selection from observed fits; every declared lane remains required."""
    registered=next((s for s in plan['scenarios'] if s['id']==scenario['id']),None)
    if registered is None or json_bytes(registered)!=json_bytes(scenario):raise ValueError('Unregistered reference scenario')
    if type(replicate) is not int or replicate not in scenario['independent_reference_replicates']:
        raise ValueError('Unregistered reference replicate')
    if ledger.connection.execute('SELECT plan_sha256 FROM binding').fetchall()!=[(identity,)]:
        raise ValueError('Reference ledger belongs to another plan')
    ledger.verify();run_id=f"{scenario['id']}__{replicate:04d}"
    if ledger.get(run_id,'terminal') is None:raise ValueError('Reference input attempt is incomplete')
    started=ledger.get(run_id,'generation_started')
    seed=replicate_seed(plan['master_seed'],plan['study_id']+':'+scenario['id'],replicate)
    if started is None or started!=dict(run_id=run_id,scenario_id=scenario['id'],replicate=replicate,seed=seed):
        raise ValueError('Reference attempt seed/identity differs from the plan')
    retained=ledger.get(run_id,'input_retained');out=Path(out);out.mkdir(parents=True,exist_ok=False)
    if retained is not None and (hashlib.sha256(retained['csv_utf8'].encode()).hexdigest()!=retained['csv_sha256']
        or retained['seed']!=seed or retained['scenario_id']!=scenario['id']
        or json_bytes(retained['metadata']['design'])!=json_bytes(scenario['design'])
        or not np.array_equal(retained['truth']['RaterTaskEtaShift'],scenario['eta_shift'])):raise ValueError('Retained input changed')
    families=[];ncells=len(scenario['design']['rater'])*len(scenario['design']['task'])
    for basis in plan['scope']['lanes']:
        result=_decode(ledger.get(run_id,basis));available=False;reason='native_screen_unavailable'
        if result is not None and 'cells' in result and retained is not None:
            try:_family_input(out/basis,result,retained,scenario,basis,plan);available=True;reason='available'
            except Exception as exc:reason=f'{type(exc).__name__}: {exc}'
        families.append(dict(basis=basis,available=available,input_status=reason,required_scalar_cells=ncells,
            required_adjustments=2,required_derived_adjustments=2))
    _write(out/'selection.json',dict(study_sha256=identity,scenario_id=scenario['id'],replicate=replicate,run_id=run_id,
        input_sha256=retained['csv_sha256'] if retained is not None else None,families=families,scientific_inference_ready=False))
    _write(out/'manifest.json',{str(p.relative_to(out)):_digest(p) for p in out.rglob('*') if p.is_file()})


def _r_family(root,rscript,source,target):
    target.mkdir(parents=True,exist_ok=False)
    with (target/'stdout.txt').open('xb') as stdout,(target/'stderr.txt').open('xb') as stderr:
        run=subprocess.run([str(rscript),'--vanilla',str(root/'validation/bias_conditional_reference.R'),
            str(source/'input.json'),str(target/'result.json')],cwd=root,stdout=stdout,stderr=stderr)
    files={p.name:_digest(p) for p in target.iterdir() if p.is_file()}
    _write(target/'execution.json',dict(returncode=run.returncode,files=files))
    if run.returncode:raise RuntimeError('Independent R exited '+str(run.returncode))
    files['execution.json']=_digest(target/'execution.json')
    return dict(result=json.loads((target/'result.json').read_bytes()),files=files)


def review_family(python,reference,plan):
    limits=plan['independent_reference']['limits']
    comparisons,adjustments=compare_reference(python['cases'],python['families'],reference,limits)
    if reference.get('family_observed_ids')!=python['observed_ids']:
        raise ValueError('R-derived observed-cell family differs from retained data')
    derived=reference.get('derived_adjustments')
    if not isinstance(derived,list) or len(derived)!=1:raise ValueError('Missing independently derived R adjustments')
    values=np.asarray(python['families'][0],dtype=float);propagated=len(values)*limits['p']+limits['adjustment']
    independent=[]
    for method in ['holm','BH']:
        a=app.adjust_p_values(values,method);b=np.asarray(derived[0][method],dtype=float)
        match=bool(a.shape==b.shape and np.allclose(a,b,rtol=0,atol=propagated,equal_nan=True))
        decisions=bool(a.shape==b.shape and np.array_equal(a<plan['screen']['alpha'],b<plan['screen']['alpha']))
        independent.append(dict(method=method,pass_=match and decisions,propagated_tolerance=propagated,
            threshold_decisions_match=decisions))
    return dict(cells=comparisons,adjustments=adjustments,derived_adjustments=independent,
        all_checks_pass=all(r['pass_'] for r in comparisons+adjustments+independent),scientific_inference_ready=False)


def execute_reference_bundle(bundle,plan,identity,root,rscript,verifier):
    """Called under the controller lock. Interrupted R jobs are never rerun."""
    bundle=Path(bundle).resolve();root=Path(root);verifier()
    selection=json.loads((bundle/'selection.json').read_bytes())
    if selection['study_sha256']!=identity:raise ValueError('Reference bundle belongs to another plan')
    for name,sha in json.loads((bundle/'manifest.json').read_bytes()).items():
        path=(bundle/name).resolve()
        if not path.is_relative_to(bundle) or _digest(path)!=sha:raise ValueError('Reference input changed: '+name)
    expected=plan['scope']['lanes']
    if [r['basis'] for r in selection['families']]!=expected:raise ValueError('Required reference lanes changed')
    scenario=next((s for s in plan['scenarios'] if s['id']==selection['scenario_id']),None)
    if (scenario is None or selection['replicate'] not in scenario['independent_reference_replicates']
        or selection['run_id']!=f"{scenario['id']}__{selection['replicate']:04d}"):
        raise ValueError('Required reference run changed')
    ncells=len(scenario['design']['rater'])*len(scenario['design']['task'])
    if any([r[k] for k in ['required_scalar_cells','required_adjustments','required_derived_adjustments']]!=[ncells,2,2]
        for r in selection['families']):raise ValueError('Required reference denominator changed')
    output=bundle/'execution';output.mkdir(exist_ok=True);reviews=[]
    ledger=AttemptLedger(output/'events.sqlite',plan_sha256=identity)
    try:
        ledger.verify()
        for family in selection['families']:
            basis=family['basis'];key=selection['run_id']+'_'+basis
            required=sum(family[k] for k in ['required_scalar_cells','required_adjustments','required_derived_adjustments'])
            if not family['available']:
                reviews.append(dict(basis=basis,status='required_input_unavailable',required_slots=required,passed_slots=0,
                    all_checks_pass=False,input_status=family['input_status']));continue
            ref=_once(ledger,key,'R',lambda:_r_family(root,rscript,bundle/basis,output/basis),verifier)
            if ref['status']=='completed':
                # Check even when a previous completed comparison is reused.
                for name,sha in ref['value']['files'].items():
                    if _digest(output/basis/name)!=sha:raise ValueError('Retained R artifact changed')
                def review():
                    python=_decode(json.loads((bundle/basis/'python.json').read_bytes()))
                    return review_family(python,ref['value']['result'],plan)
                checked=_once(ledger,key,'review',review,verifier)
            else:checked=dict(status=ref['status'])
            value=checked.get('value',{});rows=value.get('cells',[])+value.get('adjustments',[])+value.get('derived_adjustments',[])
            complete=len(rows)==required
            reviews.append(dict(basis=basis,status=checked['status'],required_slots=required,
                passed_slots=sum(bool(r['pass_']) for r in rows),all_checks_pass=bool(complete and value.get('all_checks_pass')),
                reference_status=ref['status']))
        result=dict(run_id=selection['run_id'],families=reviews,status_counts=dict(Counter(r['status'] for r in reviews)),
            required_slots=sum(r['required_slots'] for r in reviews),passed_slots=sum(r['passed_slots'] for r in reviews),
            all_checks_pass=bool(reviews) and all(r['all_checks_pass'] for r in reviews),ledger=ledger.verify(),scientific_inference_ready=False)
        path=output/'summary.json'
        if path.exists():
            if path.read_bytes()!=json_bytes(result)+b'\n':raise ValueError('Completed reference summary differs; no replacement')
        else:_write(path,result)
        return result
    finally:ledger.close()
