"""Single-attempt G3 execution over retained generated response bytes.

This is an executor component, not a registration/receipt controller. The
current draft refuses qualification execution. Development controls must use
a separately named cohort. Calibration failures do not erase known-baseline
screen lanes, and a screen failure cannot trigger another optimizer attempt.
"""
import hashlib

import numpy as np

from validation.bias_population_design import generate_shifted_sample, prepare_bias_problem, screen_from_point
from validation.mfrm_recovery_design import RecoveryDesign, replicate_seed
from validation.mml_grouped_fit import fit_grouped_joint, normal_rule
from validation.mml_grouped_likelihood import GroupedMmlLikelihood
from validation.mml_population_study import (compact_fit, compare_orders, json_bytes,
    prepare_retained_problem, truth_joint)


def _failure(status,exc):
    return dict(status=status,error=dict(type=type(exc).__name__,message=str(exc)),scientific_inference_ready=False)


def execute_bias_run(ledger,run,scenario,plan,*,verify_environment):
    development=plan['classification']=='DEVELOPMENT_ONLY' and plan['study_id'].startswith('DEVELOPMENT_ONLY_')
    if not development and not plan['qualification_execution_allowed']:
        raise ValueError('Draft G3 plan does not authorize qualification data generation')
    registered=next((s for s in plan['scenarios'] if s['id']==scenario['id']),None)
    replicate=run['replicate']
    if (registered is None or json_bytes(registered)!=json_bytes(scenario) or type(replicate) is not int
        or not 1<=replicate<=scenario['replicates']):raise ValueError('Unregistered scenario or replicate')
    expected=dict(run_id=f"{scenario['id']}__{replicate:04d}",scenario_id=scenario['id'],replicate=replicate,
        seed=replicate_seed(plan['master_seed'],plan['study_id']+':'+scenario['id'],replicate))
    if run!=expected or type(run['seed']) is not int:raise ValueError('Unregistered RunId or seed')
    run_id=run['run_id']
    if ledger.started(run_id):raise ValueError('Attempt already started; numerical replay is prohibited')
    implementation=verify_environment();ledger.append(run_id,'generation_started',run)
    if implementation is not None:ledger.append(run_id,'execution_implementation',implementation)
    try:
        spec=RecoveryDesign(**scenario['design'])
        sample=generate_shifted_sample(spec,run['seed'],scenario['eta_shift'])
        csv=sample['data'].to_csv(index=False)
        retained={k:sample[k] for k in ['truth','declared_levels','metadata']}
        retained.update(csv_utf8=csv,csv_sha256=hashlib.sha256(csv.encode()).hexdigest(),
            seed=run['seed'],scenario_id=scenario['id'])
        ledger.append(run_id,'input_retained',retained)
    except Exception as exc:
        ledger.append(run_id,'terminal',_failure('generation_failed',exc));return
    retained=ledger.get(run_id,'input_retained')
    if (retained['seed']!=run['seed'] or retained['scenario_id']!=scenario['id']
        or json_bytes(retained['metadata']['design'])!=json_bytes(scenario['design'])
        or not np.array_equal(retained['truth']['RaterTaskEtaShift'],scenario['eta_shift'])):
        raise ValueError('Retained generation differs from the declared run')
    qs=plan['estimator']['quadrature'];lanes={};fits=[]
    try:
        data,_,invalid=prepare_retained_problem(retained,spec.model,qs[0])
        if invalid:raise ValueError(invalid)
        prep,problem=prepare_bias_problem(data,spec,qs[0])
        sample=dict(retained,data=data);point=truth_joint(problem.config,retained['truth'])
        ledger.append(run_id,'preparation',dict(status='prepared',input_sha256=retained['csv_sha256']))
    except Exception as exc:
        ledger.append(run_id,'terminal',_failure('preparation_failed',exc));return

    def screen(basis,calibration,theta,available,reason):
        verify_environment()
        if not available:result=dict(status=reason,scientific_inference_ready=False)
        else:
            ledger.append(run_id,basis+'_started',dict(input_sha256=retained['csv_sha256']))
            try:
                result=screen_from_point(prep,problem,calibration,theta,sample,basis=basis,
                    omit_extreme=plan['screen']['omit_extreme'],max_abs=plan['screen']['max_abs'],
                    alpha=plan['screen']['alpha'],min_n=plan['screen']['min_positive_weight_observations'],
                    practical_logit=plan['screen']['practical_logit'])
                # Large rowwise reference transports are needed only for the
                # prespecified independent-reference subset. Other runs retain
                # inputs, calibration points and all cell/pairwise evidence.
                if replicate not in scenario['independent_reference_replicates']:
                    del result['reference_rows']
                    result['reference_storage']='Reconstruct from retained input and calibration/Person basis'
            except Exception as exc:result=_failure('screen_execution_failed',exc)
        ledger.append(run_id,basis,result);lanes[basis]=result.get('status','computed')

    screen('known_calibration_known_person',point,retained['truth']['Person'],True,None)
    verify_environment();ledger.append(run_id,'known_posterior_started',dict(point=point,quadrature=qs))
    try:
        grouped=GroupedMmlLikelihood(data,problem.config)
        known=[grouped.evaluate(point,*normal_rule(q)) for q in qs]
        errors=dict(nll_per_person=abs(known[0]['nll']-known[1]['nll'])/spec.persons,
            cross_score_supnorm=float(np.max(abs(known[0]['gradient']-known[1]['gradient']))),
            posterior_moment=max(float(np.max(abs(known[0][k]-known[1][k]))) for k in ['eap','sd']))
        known_ok=all(v<=plan['estimator']['quadrature_limits'][k] for k,v in errors.items())
        known_result=dict(status='evaluated',numerical_pass=known_ok,errors=errors,eap=known[0]['eap'],sd=known[0]['sd'])
    except Exception as exc:
        known_result=_failure('known_posterior_failed',exc);known_ok=False
    ledger.append(run_id,'known_posterior',known_result)
    screen('known_calibration_EAP',point,known_result.get('eap'),known_ok,'known_posterior_unavailable')
    for q in qs:
        verify_environment();phase=f'Q{q}'
        ledger.append(run_id,phase+'_started',dict(input_sha256=retained['csv_sha256'],q=q,
            structural_start=problem.structural_start,sigma_start=problem.sigma_start))
        try:
            fit=compact_fit(fit_grouped_joint(data,problem,q=q))
            ledger.append(run_id,phase+'_finished',fit)
        except Exception as exc:
            fit=dict(_failure('calibration_execution_failed',exc),numerical_pass=False)
            ledger.append(run_id,phase+'_failed',fit)
        fits.append(fit)
    try:comparison=compare_orders(data,problem,fits,plan)
    except Exception as exc:comparison=dict(_failure('quadrature_comparison_failed',exc),pass_=False)
    ledger.append(run_id,'calibration_assessment',comparison);eligible=comparison['pass_']
    screen('fitted_calibration_known_person',fits[0].get('point'),retained['truth']['Person'],eligible,'calibration_unavailable')
    screen('fitted_calibration_EAP',fits[0].get('point'),fits[0].get('evaluated',{}).get('eap'),eligible,'calibration_unavailable')
    ledger.append(run_id,'terminal',dict(status='completed',lanes=lanes,scientific_inference_ready=False))
