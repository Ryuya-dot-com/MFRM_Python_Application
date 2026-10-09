"""Declared cell-shift DGP and native bias-screen transport for G3 controls.

The generator changes response probabilities, not the additive calibration
model. Fitted-calibration screens describe residual departures; their scalar
intervals are not silently relabelled as intervals for the original DGP effect.
"""
import copy
import numpy as np
import pandas as pd
from scipy.special import softmax

import streamlit_app as app
from validation.mfrm_recovery_design import (STREAMS,category_probabilities,
    generate_population_sample)
from validation.mfrm_recovery_metrics import calibration_coordinates
from validation.mml_population_study import truth_joint


GENERATOR_VERSION='mfrm_bias_rater_task_shift_v1'


def prepare_bias_problem(data,design,q):
    from validation.mml_gpcm_joint_check import prepare_problem
    from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data
    rating_max=design.rating_min+len(design.steps[0])
    if design.model=='GPCM':return prepare_problem(data,q,rating_min=design.rating_min,rating_max=rating_max)
    return prepare_app_free_sd_problem_from_data(data,person_col='Person',facet_cols=['Rater','Task','Criterion'],
        score_col='Score',rating_min=design.rating_min,rating_max=rating_max,model=design.model,
        noncenter_facet='Criterion',step_facet=None if design.model=='RSM' else 'Criterion',quadrature_points=q)


def contrast_matrix(raters,tasks,magnitude):
    if (type(raters) is not int or type(tasks) is not int or min(raters,tasks)<2
        or isinstance(magnitude,(bool,np.bool_)) or not np.isfinite(magnitude)):
        raise ValueError('At least two raters/tasks and a finite effect required')
    matrix=np.zeros((raters,tasks));matrix[:2,:2]=float(magnitude)*np.array([[1,-1],[-1,1]])
    return matrix


def shifted_probabilities(design,theta,shift):
    design.validate();theta=np.asarray(theta,dtype=float);shift=np.asarray(shift,dtype=float)
    if (theta.shape!=(design.persons,) or not np.isfinite(theta).all()
        or shift.shape!=(len(design.rater),len(design.task)) or not np.isfinite(shift).all()
        or np.max(abs(shift.sum(axis=0)))>1e-12 or np.max(abs(shift.sum(axis=1)))>1e-12):
        raise ValueError('Finite native Person vector and row/column sum-zero Rater x Task shift required')
    if not np.any(shift):return category_probabilities(design,theta)
    steps=np.asarray(design.steps)
    if design.model=='RSM':steps=np.repeat(steps,len(design.criterion),axis=0)
    eta=(theta[:,None,None,None]-np.asarray(design.rater)[None,:,None,None]
        -np.asarray(design.task)[None,None,:,None]-np.asarray(design.criterion)[None,None,None,:]
        +shift[None,:,:,None])
    cumulative=np.c_[np.zeros(len(steps)),np.cumsum(steps,axis=1)]
    slope=np.exp(design.log_slopes) if design.model=='GPCM' else np.ones(len(design.criterion))
    logits=slope[None,None,None,:,None]*(eta[...,None]*np.arange(steps.shape[1]+1)-cumulative[None,None,None,:,:])
    if not np.isfinite(logits).all():raise FloatingPointError('Nonfinite declared shift logits')
    return softmax(logits,axis=-1)


def generate_shifted_sample(design,seed,shift):
    sample=generate_population_sample(design,seed)
    probability=shifted_probabilities(design,sample['theta'],shift)
    # Reuse the exact declared response stream and full-cell draw order.
    # Assignment and Person draws remain identical for null/effect controls.
    rng=np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence(seed,spawn_key=(STREAMS['response'],))))
    uniform=rng.random(probability.shape[:-1]);cdf=probability.cumsum(axis=-1);cdf[...,-1]=1.
    scores=np.sum(uniform[...,None]>=cdf,axis=-1)+design.rating_min
    data=sample['data'].copy()
    indices=[data[name].map({label:i for i,label in enumerate(levels)}).to_numpy(dtype=int)
        for name,levels in sample['declared_levels'].items()]
    data['Score']=scores[tuple(indices)];sample['data']=data
    sample['truth']['RaterTaskEtaShift']=np.array(shift,dtype=float)
    sample['metadata'].update(base_generator_version=sample['metadata']['generator_version'],generator_version=GENERATOR_VERSION,
        shift_parameterization='Additive eta shift before multiplying by Criterion slope; row/column sums zero in declared level space',
        effect_target='DGP Rater x Task effect; fitted additive calibration may absorb part of this effect')
    return sample


def screen_from_point(prep,problem,point,theta,sample,*,basis,omit_extreme=True,max_abs=10.,weights=None,
        alpha=.05,min_n=5,practical_logit=.5):
    """Use the actual application cell, pairwise and multiplicity code paths.

    A caller must choose a named calibration/Person basis. Only the fully
    known DGP basis identifies its scalar target with the generated cell shift.
    Numerical eligibility of estimated points remains a caller-owned gate.
    """
    allowed={'known_calibration_known_person','known_calibration_EAP',
        'fitted_calibration_known_person','fitted_calibration_EAP'}
    if basis not in allowed:raise ValueError('Unknown calibration/Person basis')
    point=np.asarray(point,dtype=float);theta=np.asarray(theta,dtype=float)
    if theta.shape!=(problem.config['n_person'],) or not np.isfinite(theta).all():
        raise ValueError('A finite Person vector in native sorted order is required')
    coordinates=calibration_coordinates(problem.config,point)
    if basis.startswith('known_calibration') and not np.array_equal(point,truth_joint(problem.config,sample['truth'])):
        raise ValueError('Known calibration must equal the retained native DGP truth')
    if basis.endswith('known_person') and not np.array_equal(theta,np.asarray(sample['truth']['Person'])):
        raise ValueError('Known Person values must equal the retained DGP truth')
    levels=sample['declared_levels']
    for facet,declared in levels.items():
        if list(prep['levels'][facet])!=list(declared):raise ValueError('Prepared and declared level order differs')
    for column in [*levels,'Score']:
        if not np.array_equal(prep['data'][column].to_numpy(),sample['data'][column].to_numpy()):
            raise ValueError('Prepared rows differ from retained sample')
    work=copy.deepcopy(prep)
    if weights is not None:
        weights=np.asarray(weights,dtype=float)
        if weights.shape!=(len(work['data']),) or not np.isfinite(weights).all() or np.any(weights<0):
            raise ValueError('One finite nonnegative screening weight per observed row required')
        work['data']['Weight']=weights
    params=app.expand_params(point[:-1],problem.sizes,problem.config)
    res=dict(prep=work,config=problem.config,params=params,facets={'person':pd.DataFrame({'Estimate':theta})})
    obs=app.compute_obs_table(res)
    measures=pd.DataFrame([dict(Facet=row['kind'],Level=row['level'],Estimate=row['value'],SE=np.nan)
        for row in coordinates if row['kind'] in problem.config['facet_names']])
    output=app.estimate_bias_interaction(res,dict(obs=obs,measures=measures),'Rater','Task',
        omit_extreme=omit_extreme,max_abs=max_abs)
    if '_skip_reason' in output:raise ValueError(output['_skip_reason'])
    table=output['table'];screen=app.build_dff_bias_screening_table(output,alpha=alpha,min_n=min_n,practical_logit=practical_logit)
    pairwise=app.calc_bias_pairwise(table,'Rater','Task')
    shifts=np.asarray(sample['truth']['RaterTaskEtaShift'],dtype=float)
    lookup={(r,t):float(shifts[i,j]) for i,r in enumerate(levels['Rater']) for j,t in enumerate(levels['Task'])}
    screen_lookup={(r.FacetA_Level,r.FacetB_Level):r for r in screen.itertuples(index=False)}
    oracle=basis=='known_calibration_known_person';cells=[]
    for (rater,task),truth in lookup.items():
        found=table.loc[table.FacetA_Level.eq(rater)&table.FacetB_Level.eq(task)]
        if found.empty:
            cells.append(dict(rater=rater,task=task,DGP_shift=truth,status='unobserved_declared_cell',
                p_available=False,interval_available=False,oracle_target=oracle));continue
        row=found.iloc[0];candidate=screen_lookup[rater,task]
        valid=bool(candidate.InferenceAvailable);ci=bool(row['Profile CI Status']=='ok')
        cells.append(dict(rater=rater,task=task,DGP_shift=truth,status=row['Cell Status'],
            estimate=row['Bias Size'],conditional_se=row['S.E.'],p=row['LR Prob.'],p_holm=candidate.p_holm,
            p_available=valid,interval_available=ci,oracle_target=oracle,
            oracle_error=float(row['Bias Size']-truth) if oracle and np.isfinite(row['Bias Size']) else None,
            oracle_covered=bool(row['Profile CI Lower']<=truth<=row['Profile CI Upper']) if oracle and ci else None,
            raw_rejection=bool(candidate.p<alpha) if valid else None,
            holm_rejection=bool(candidate.p_holm<alpha) if valid else None,
            operational_statistical_flag=bool(candidate.EvidenceLevel in ['Strong review','Statistical review']) if valid else None,
            review_label=candidate.EvidenceLevel,review_flag=bool(candidate.Flag),observations=int(row['ObsN']),persons=int(row['PersonN'])))
    contrasts=pairwise.to_dict('records')
    for row in contrasts:
        row.update(DGP_local_contrast=(lookup[row['Target'],row['Context1']]-lookup[row['Target'],row['Context2']])/row['Target Sign'],oracle_target=oracle)
    idx=app.build_indices(work,step_facet=problem.config['step_facet'],slope_facet=problem.config.get('slope_facet'))
    reference=work['data'][['Person','Rater','Task','Criterion']].astype(str).copy()
    reference['eta']=app.compute_eta(idx,params,problem.config,theta_override=theta)
    reference['score']=idx['score_k'];reference['weight']=idx['weight']
    reference['slope']=params['slopes'][idx['slope_idx']] if problem.config['model']=='GPCM' else 1.
    steps=np.tile(params['steps'],(len(reference),1)) if problem.config['model']=='RSM' else params['steps_mat'][idx['step_idx']]
    for k in range(steps.shape[1]):reference['step'+str(k+1)]=steps[:,k]
    return dict(basis=basis,observed_family_size=len(table),declared_cells=len(lookup),
        family_complete=len(table)==len(lookup),cells=cells,pairwise=contrasts,
        table=table.to_dict('records'),screen=screen.to_dict('records'),
        reference_rows=reference.to_dict('records'),
        scientific_inference_ready=False,
        interpretation='Original DGP shift is a scalar recovery target only in the known-calibration/known-Person lane; fitted/EAP lanes assess the operational residual screen')
