"""Export the preselected independent-reference cases from retained bytes.

Selection depends only on the frozen scenario/replicate list. Missing inputs
or fits remain explicit required cases. No references are selected for looking
well behaved, and no reference outcome changes recovery denominators.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import streamlit_app as app

from validation.mfrm_recovery_metrics import calibration_coordinates
from validation.mml_grouped_fit import normal_rule
from validation.mml_population_study import json_bytes, prepare_retained_problem, truth_joint


CONFIG_KEYS = ['model','facet_names','facet_levels','facet_signs','noncenter_facet',
    'step_facet','slope_facet','n_cat','rating_min','rating_max']


def _write(path,value):
    with path.open('xb') as handle:handle.write(json_bytes(value)+b'\n')


def export_reference_bundle(out,ledger,scenario,replicate,plan,study_sha256):
    if replicate not in set(scenario['independent_R_replicates']+scenario['TAM_replicates']):
        raise ValueError('Reference replicate was not selected before generation')
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    run_id=f"{scenario['id']}__{replicate:04d}"
    retained=ledger.get(run_id,'input_retained')
    expected=[];finite=[];continuous=[];tam=[];native=[];python=[]
    for model in scenario['fit_models']:
        problem=None;data=None
        try:
            if retained is None:raise ValueError('Retained input unavailable')
            data,problem,invalid=prepare_retained_problem(retained,model,plan['estimator']['primary_q'])
            if invalid:raise ValueError(invalid)
            config={key:problem.config[key] for key in CONFIG_KEYS}
        except Exception as exc:
            input_error=dict(type=type(exc).__name__,message=str(exc))
        else: input_error=None
        for q in plan['estimator']['quadrature']:
            identifier=f'{run_id}_{model}_Q{q}'
            fit=ledger.get(run_id,f'{model}_Q{q}_finished')
            available=input_error is None and fit is not None and fit.get('point') is not None
            if replicate in scenario['independent_R_replicates']:
                expected.append(dict(id=identifier,lane='finite_R',available=available,input_error=input_error,
                    native_fit_status=fit.get('status') if fit else 'fit_unavailable'))
                if available:
                    point=np.array(fit['point']);nodes,weights=normal_rule(q);keep=weights>0
                    case=dict(id=identifier,config=config,data=data.to_dict('records'),point=point,
                        requested_order=q,nodes=nodes[keep],weights=weights[keep],group_rows=True,
                        point_kind='fitted',parent_input_sha256=retained['csv_sha256'])
                    finite.append(case)
                    quad=dict(nodes=nodes[keep]*np.exp(point[-1]),weights=weights[keep],sd=float(np.exp(point[-1])))
                    nll,gradient=app.mfrm_loglik_mml_value_grad(point[:-1],problem.idx,problem.config,problem.sizes,quad,include_log_sigma=True)
                    native.append(dict(id=identifier,nll=nll,gradient=gradient))
                    expanded=[row for row in calibration_coordinates(config,point) if row['kind']!='LogPopulationSD']
                    matrix=np.array([row['jacobian'] for row in expanded]);cov=fit.get('covariance')
                    value=dict(id=identifier,**fit['evaluated'],covariance=cov,expanded=[row['value'] for row in expanded],
                        expanded_jacobian=matrix,expanded_covariance=matrix@cov@matrix.T if cov is not None else None,
                        operationally_eligible=bool(fit.get('numerical_pass')))
                    python.append(value)
                if q==plan['estimator']['primary_q']:
                    expected.append(dict(id=identifier,lane='continuous',available=available,input_error=input_error))
                    if available:continuous.append(dict(case,continuous_only=True))
            if q==plan['estimator']['primary_q'] and replicate in scenario['TAM_replicates'] and model in ['RSM','PCM']:
                expected.append(dict(id=identifier,lane='TAM',available=input_error is None,input_error=input_error))
                if input_error is None:
                    tam.append(dict(id=identifier,config=config,data=data.to_dict('records'),
                        point=fit.get('point') if fit else None,native_fit_status=fit.get('status') if fit else 'unavailable'))
        if model==scenario['primary_model'] and replicate in scenario['independent_R_replicates']:
            identifier=f'{run_id}_{model}_truth'
            expected.append(dict(id=identifier,lane='continuous',available=input_error is None,input_error=input_error))
            if input_error is None:
                point=truth_joint(config,retained['truth']);nodes,weights=normal_rule(plan['estimator']['primary_q']);keep=weights>0
                continuous.append(dict(id=identifier,config=config,data=data.to_dict('records'),point=point,
                    nodes=nodes[keep],weights=weights[keep],group_rows=True,continuous_only=True,
                    point_kind='true_calibration',parent_input_sha256=retained['csv_sha256']))
    _write(out/'reference_plan.json',dict(study_sha256=study_sha256,run_id=run_id,scenario_id=scenario['id'],
        seed=retained['seed'] if retained else None,required_cases=expected,
        missing_required_cases=sum(not case['available'] for case in expected),scientific_inference_ready=False))
    # Low-level reference tools retain their own non-qualifying classification;
    # the separate sealed study identity binds the prospective selection.
    _write(out/'finite_input.json',dict(classification='DEVELOPMENT_ONLY',cases=finite))
    _write(out/'continuous_input.json',dict(classification='DEVELOPMENT_ONLY',cases=continuous,
        continuous_settings=plan['numerical_reference']['continuous_settings']))
    _write(out/'native_rowwise.json',native);_write(out/'python_finite.json',python)
    tam_parent=out/'tam_input';tam_parent.mkdir()
    _write(tam_parent/'reference_input.json',dict(classification='DEVELOPMENT_ONLY',cases=tam))
    for case in tam:
        with (tam_parent/(case['id']+'.csv')).open('xb') as handle:handle.write(retained['csv_utf8'].encode('utf-8'))
    _write(tam_parent/'manifest.json',{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in tam_parent.iterdir() if p.is_file()})
    _write(out/'manifest.json',{str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in out.rglob('*') if p.is_file()})
    return dict(run_id=run_id,required_cases=len(expected),available_cases=sum(c['available'] for c in expected))
