"""Prespecified numerical comparisons for the prospective population study.

All failed and unavailable comparisons remain in the required-case denominator.
Reference agreement never changes operational eligibility or sampling outcomes.
"""
import numpy as np
import pandas as pd

from validation.mfrm_recovery_metrics import calibration_coordinates
from validation.mml_gpcm_integral import GpcmIntegral
from validation.mml_grouped_likelihood import GroupedMmlLikelihood
from validation.mml_rasch_integral import RaschIntegral


def distance(a,b):
    a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
    if a.shape!=b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Reference shapes or finite values differ')
    return float(np.max(abs(a-b)))


def expanded_se(config,point,covariance):
    covariance=np.asarray(covariance,dtype=float)
    if (covariance.shape!=(len(point),len(point)) or not np.isfinite(covariance).all()
        or np.max(abs(covariance-covariance.T))>1e-10 or np.linalg.eigvalsh(covariance).min()<=0):
        raise ValueError('Full positive-definite joint covariance required')
    jacobian=np.array([r['jacobian'] for r in calibration_coordinates(config,point)])
    variable=np.any(jacobian!=0,axis=1)
    variance=np.einsum('ij,jk,ik->i',jacobian[variable],covariance,jacobian[variable])
    if np.any(variance<=0):raise ValueError('Nonpositive expanded variance')
    return np.sqrt(variance)


def finite_review(case,python,native,r,plan):
    n=len({row['Person'] for row in case['data']});limits=plan['numerical_reference']['finite_R_limits']
    errors=dict(nll_per_person=distance(python['nll'],r['nll'])/n,
        gradient_per_person=distance(python['gradient'],r['gradient'])/n,
        information_per_person=distance(python['information'],r['information'])/n,
        hessian_step_per_person=float(r['hessian_step_disagreement'])/n,
        posterior_moment=max(distance(python[k],r[k]) for k in ['eap','sd']),
        expansion=max(distance(python['expanded'],r['expanded']),distance(python['expanded_jacobian'],r['expanded_jacobian'])))
    checks={k:v<=limits[k] for k,v in errors.items()}
    checks['rowwise_R']=r['rowwise_nll_disagreement']/n<=limits['nll_per_person']
    checks['native_rowwise']=all(distance(native[k],python[k])<=plan['numerical_reference']['native_limits'][k] for k in ['nll','gradient'])
    checks['positive_definite']=bool(r['positive_definite'])
    checks['joint_covariance']=checks['expanded_se']=False
    if python['covariance'] is not None and r['covariance'] is not None:
        a,b=np.array(python['covariance']),np.array(r['covariance'])
        scale=np.sqrt(np.outer(np.diag(b),np.diag(b)))
        errors['covariance_scaled']=float(np.max(abs(a-b)/scale))
        errors['expanded_se_relative']=float(np.max(abs(expanded_se(case['config'],case['point'],a)/
            expanded_se(case['config'],case['point'],b)-1)))
        checks['joint_covariance']=errors['covariance_scaled']<=limits['covariance_scaled']
        checks['expanded_se']=errors['expanded_se_relative']<=limits['expanded_se_relative']
    return dict(errors=errors,checks=checks,all_checks_pass=all(checks.values()),scientific_inference_ready=False)


def continuous_values(case,plan):
    data=pd.DataFrame(case['data']);cfg=case['config']
    reference=(GpcmIntegral if cfg['model']=='GPCM' else RaschIntegral)(data,cfg)
    information=case['point_kind']=='fitted'
    return [reference.continuous(case['point'],**setting,information=information)
        for setting in plan['numerical_reference']['continuous_settings']]


def reference_checks(values,plan,n,*,information):
    a,b=values;limits=plan['numerical_reference']['continuous_limits']
    errors={k:distance(a[k],b[k]) for k in ['nll','gradient','eap','sd']+(['information'] if information else [])}
    checks=dict(reference_nll=errors['nll']<=limits['reference_nll'],
        reference_gradient=errors['gradient']<=limits['reference_gradient'],
        reference_moments=max(errors['eap'],errors['sd'])<=limits['reference_moments'],
        relative_mass_error=b['numeric_relative_mass_error_sum']<=limits['relative_mass_error'],
        tail_mass=b['tail_relative_mass_bound_sum']<=limits['tail_mass'],
        tail_gradient=b['tail_gradient_abs_bound']<=limits['tail_gradient'])
    if information:
        checks.update(reference_information=errors['information']/n<=limits['reference_information_per_person'],
            tail_information=b['tail_information_abs_bound']<=limits['tail_information'],
            numeric_information_allowance=b['numeric_information_error_allowance']/n<=limits['numeric_information_allowance_per_person'])
    return errors,checks


def continuous_review(case,values,r,python,plan):
    n=len({row['Person'] for row in case['data']});limits=plan['numerical_reference']['continuous_limits']
    fitted=case['point_kind']=='fitted';tight=values[-1]
    refinement,checks=reference_checks(values,plan,n,information=fitted)
    r_errors={k:distance(tight[k],r[-1][k]) for k in ['nll','eap','sd','density_log_sigma_score']}
    r_refinement={k:distance(r[0][k],r[1][k]) for k in ['nll','eap','sd','density_log_sigma_score']}
    checks.update(r_nll=max(r_errors['nll'],r_refinement['nll'])<=limits['r_nll'],
        r_moments=max(r_errors['eap'],r_errors['sd'],r_refinement['eap'],r_refinement['sd'])<=limits['r_moments'],
        r_sigma_score=max(r_errors['density_log_sigma_score'],r_refinement['density_log_sigma_score'])<=limits['r_sigma_score'],
        r_mass=r[-1]['numeric_relative_mass_error_sum']<=limits['relative_mass_error'],
        r_tail=r[-1]['tail_relative_mass_bound_sum']<=limits['tail_mass'])
    if python is None:
        python=GroupedMmlLikelihood(pd.DataFrame(case['data']),case['config']).evaluate(
            np.array(case['point']),np.array(case['nodes']),np.array(case['weights']))
    errors={k:distance(python[k],tight[k]) for k in ['nll','gradient','eap','sd']}
    checks.update(finite_nll=errors['nll']/n<=limits['nll_per_person'],
        finite_gradient=errors['gradient']<=limits['score_supnorm'],
        finite_moments=max(errors['eap'],errors['sd'])<=limits['posterior_moment'])
    if fitted:
        checks['stationary']=float(np.max(abs(tight['gradient'])))<=plan['estimator']['numerical_limits']['score']
        checks['expanded_se']=False
        if np.linalg.eigvalsh(tight['information']).min()>0 and python['covariance'] is not None:
            errors['expanded_se_relative']=float(np.max(abs(expanded_se(case['config'],case['point'],python['covariance'])/
                expanded_se(case['config'],case['point'],np.linalg.inv(tight['information']))-1)))
            checks['expanded_se']=errors['expanded_se_relative']<=limits['expanded_se_relative']
    return dict(refinement=refinement,R_error=r_errors,R_refinement=r_refinement,integration_error=errors,
        checks=checks,all_checks_pass=all(checks.values()),scientific_inference_ready=False)


def tam_review(case,design,tam,python,python_point,plan):
    cfg=case['config'];data=pd.DataFrame(case['data']);point=np.array(tam['point'])
    reference=RaschIntegral(data,cfg);limits=plan['numerical_reference']['TAM']['limits']
    item_map={tuple(row[f] for f in cfg['facet_names']):i for i,row in enumerate(design['items'])}
    indices=[item_map[tuple(row[f] for f in cfg['facet_names'])] for row in case['data']]
    response=np.array(design['response'],dtype=float)
    checks=dict(A=distance(np.array(design['A'])[indices],reference.design)<=limits['design'],
        B=distance(np.array(design['B'])[:,:,0],np.broadcast_to(reference.k,(len(design['items']),len(reference.k))))<=limits['design'],
        response=bool(np.isfinite(response).sum()==len(data) and np.array_equal(response[reference.person,indices],reference.y)),
        person_order=design['persons']==reference.person_labels,terminated=bool(tam['terminated_below_cap']))
    sigma=np.exp(point[-1]);theta=np.array(tam['theta']).ravel()
    weights=np.exp(-.5*(theta/sigma)**2);weights/=weights.sum()
    # Only value/moments use this fixed-theta grid, not its log-SD derivatives.
    fixed=GroupedMmlLikelihood(data,cfg).evaluate(point,theta/sigma,weights)
    values=[reference.continuous(point,**s,information=True) for s in plan['numerical_reference']['continuous_settings']]
    tight=values[-1];refinement,check=reference_checks(values,plan,reference.n,information=True);checks.update(check)
    errors=dict(grid_nll_per_person=abs(fixed['nll']-tam['deviance']/2)/reference.n,
        joint_score=float(np.max(abs(tight['gradient']))))
    checks.update(grid_nll=errors['grid_nll_per_person']<=limits['nll_per_person'],joint_score=errors['joint_score']<=limits['joint_score'])
    for field in ['person','fresh_person']:
        persons=pd.DataFrame(tam[field]).set_index('pid').loc[reference.person_labels]
        errors[field+'_moments']=max(distance(persons['EAP'],fixed['eap']),distance(persons['SD.EAP'],fixed['sd']))
    checks['fresh_moments']=errors['fresh_person_moments']<=limits['fresh_moments']
    checks['refit_point']=checks['refit_nll']=checks['expanded_se']=False
    if python is not None and python_point is not None:
        errors['coordinate']=distance(point,python_point)
        errors['nll_per_person']=distance(tight['nll'],python['nll'])/reference.n
        checks['refit_point']=errors['coordinate']<=limits['coordinate']
        checks['refit_nll']=errors['nll_per_person']<=limits['nll_per_person']
        if python['covariance'] is not None and np.linalg.eigvalsh(tight['information']).min()>0:
            errors['expanded_se_relative']=float(np.max(abs(expanded_se(cfg,point,np.linalg.inv(tight['information']))/
                expanded_se(cfg,python_point,python['covariance'])-1)))
            checks['expanded_se']=errors['expanded_se_relative']<=limits['expanded_se_relative']
    return dict(continuous=values,refinement=refinement,errors=errors,checks=checks,all_checks_pass=all(checks.values()),
        native_tam_se_definition='Retained native xsi SE is not treated as full joint inverse information',scientific_inference_ready=False)
