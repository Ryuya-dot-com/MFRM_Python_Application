#!/usr/bin/env python3
"""Paired missingness/weight probes on new PCM responses; no inference qualification."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import copy
from dataclasses import replace
from decimal import Decimal, localcontext
import json
from pathlib import Path
import platform
import runpy
import subprocess
import traceback

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import least_squares
from scipy.special import logsumexp, roots_hermitenorm

from mml_observed_pcm_refit import ROOT, ObservedPCMReference, discrepancies, app_evaluate
from mml_pcm_continuous_refit import dump, sha
from mfrm_app.mml_engine_v2 import StationarityContract, run_free_sd_two_stage, assess_free_sd_stationarity
from mfrm_app.mml_stationarity import JointPolishOptions, polish_joint_free_sd, information_diagnostics, audit_joint_gradient
from mfrm_app.mml_quadrature_sensitivity import QuadratureSensitivityContract, assess_quadrature_sensitivity
import streamlit_app as app

PREVIOUS = ROOT/'validation/generated/mml_observed_pcm_continuous_20260914'
GH_PREVIOUS = ROOT/'validation/generated/mml_observed_pcm_refit_20260914'
TEST_SOURCE = ROOT/'tests/test_mml_free_sd_stationarity_adapter.py'
PROTOCOL = dict(
    id='independent_observation_designs_20260914_v1', classification='OBSERVED_DEVELOPMENT_ONLY',
    scientific_inference_ready=False, qualification_eligible=False,
    question='For each new response dataset, how do independent random deletion and response weighting affect finite-GH stationarity, Q121/Q181 sensitivity and total scoring differences from a refitted continuous reference?',
    datasets=[dict(id='ordinary', seed=2026091441, sigma=1.), dict(id='wide', seed=2026091442, sigma=2.5)],
    persons=32, designs=['complete_unit', 'missing_unit', 'missing_weighted'],
    randomization='PCG64 with SeedSequence(seed).spawn(3) for responses, deletion, weights. Retain each row when independent deletion uniform >= 0.30. Uniform weight draw in {0,0.5,1.5,2}. No response edits or rejected/regenerated datasets.',
    pairing='Within each dataset all designs share responses; both missing designs share exactly the same retained rows. Across ordinary/wide both seed and generating SD differ, so this is not an isolated causal SD contrast or a replicated performance-rate study.',
    true_structural=[.4, -.2, -.1, -.25, .35, -1., -.1, -.7, .2],
    covariate='Fixed linspace(-1,2,N); theta = 0.4*x + sigma*z. Rater 1 fixed +0.35, positive rater sign; task mean 0.2 and negative sign; negative criterion sign; sum-zero PCM steps.',
    orders=[121,181], starts=[.5,3.], structural_start=[0.,0.,0.,0.,0.,-1.,0.,-1.,0.],
    continuous_starts='Own Q121 and Q181 endpoints from start SD 0.5; two raw likelihood polishes then two bounded score-norm polishes at each start. Preselected Q181-derived continuous reference, admitted only when both continuous starts pass and agree.',
    r_selection='Q181/start SD 0.5 and its refitted continuous reference, in every design cell; two integration settings.',
    decimal_selection='Q181/start SD 0.5 anchor/final, missing_weighted cells only, at 60 and 90 digits.',
    finite_parity=dict(nll=1e-8, gradient=1e-8, eap=1e-9, sd=1e-9),
    arithmetic_gradient_limit=1e-6, cross_start_limit=1e-7, replay_limit=1e-9,
    decimal_limits=dict(absolute=1e-18, relative=1e-11, precision=1e-45),
    normalization='Stationarity/Q-sensitivity per-observation normalization uses actual response-weight sum; continuous NLL loss per person uses positive-weight persons. Neither is an effective-sample-size claim.',
    workers=2,
    failure_policy='Freeze all settings/sources before generation; save all full responses, masks and weights before fitting. Retain every solver failure and unmet target. No threshold changes or best-case replacement; exceptions leave partial artifacts and a failure record.',
    scope='Two independent numerical probes, each with paired observation designs, not coverage/recovery rates or general missingness/survey-weight inference. No native TAM/ConQuest fit, public estimator or default change.',
)


def generate(condition):
    response_rng, missing_rng, weight_rng = [np.random.Generator(np.random.PCG64(s))
        for s in np.random.SeedSequence(condition['seed']).spawn(3)]
    n = PROTOCOL['persons']
    x = np.linspace(-1.,2.,n)
    z = response_rng.normal(size=n)
    theta = .4*x + condition['sigma']*z
    rater, task, criterion = [.35,-.2], [-.1,.5], [-.25,.35]
    cumulative = np.array([[0.,-1.,-1.1,0.],[0.,-.7,-.5,0.]])
    rows, probabilities = [], []
    for p in range(n):
        for r in range(2):
            for t in range(2):
                for c in range(2):
                    logits = np.arange(4)*(theta[p]+rater[r]-task[t]-criterion[c])-cumulative[c]
                    probability = np.exp(logits-logsumexp(logits))
                    y = int(response_rng.choice(4,p=probability))
                    rows.append(dict(Person=f'P{p:03d}',Rater=f'R{r+1}',Task=f'T{t+1}',Criterion=f'C{c+1}',Score=y))
                    probabilities.append(probability)
    deletion_uniform = missing_rng.random(len(rows))
    weights = weight_rng.choice([0.,.5,1.5,2.],size=len(rows))
    return dict(condition=condition,n_person=n,x=x,z=z,theta=theta,full_responses=rows,
                generating_probabilities=probabilities,deletion_uniform=deletion_uniform,
                retained=deletion_uniform >= .30,weights=weights)


def prepare(template, data, design):
    config = copy.deepcopy(template.config)
    n = data['n_person']
    config['n_person'] = n
    config['population_model']['X'] = np.array(data['x'])[:,None]
    people = [f'P{p:03d}' for p in range(n)]
    config['theta_spec'] = app.build_facet_constraint(people)
    frame = pd.DataFrame(data['full_responses'])
    keep = np.ones(len(frame),bool) if design == 'complete_unit' else np.array(data['retained'],bool)
    frame = frame.loc[keep].copy()
    frame['Weight'] = np.array(data['weights'])[keep] if design == 'missing_weighted' else 1.
    levels = dict(Person=people,**config['facet_levels'])
    for name in ['Person',*config['facet_names']]:
        frame[name] = pd.Categorical(frame[name],categories=levels[name])
    frame['score_k'] = frame['Score'].astype(int)
    idx = app.build_indices(dict(data=frame,levels=levels,facet_names=config['facet_names']),step_facet='Criterion')
    return replace(template,config=config,idx=idx,observations=float(frame['Weight'].sum()),
                   structural_start=tuple(PROTOCOL['structural_start']),structural_bounds=((None,None),)*9), np.flatnonzero(keep)


def continuous_fit(reference, start, spec):
    settings, limits = spec['continuous'], spec['continuous']['limits']
    evaluate = lambda p:reference.continuous(p,**settings['integration'])
    tight = lambda p:reference.continuous(p,**settings['tight_integration'])
    def joint(p):
        v=evaluate(p);return v['nll'],v['gradient']
    raw = lambda p,sd:evaluate(np.r_[p,np.log(sd)])['nll']
    def structural(p,sd):
        v=evaluate(np.r_[p,np.log(sd)]);return v['nll'],v['gradient'][:-1]
    point=np.array(start); initial=evaluate(point); preliminary=[]
    for _ in range(2):
        fitted=polish_joint_free_sd(point[:-1],np.exp(point[-1]),raw,structural,
            sigma_bounds=(.05,10.),joint_value_gradient=joint,options=JointPolishOptions(**settings['preliminary']))
        preliminary.append(fitted.to_dict());point=np.array(fitted.joint_coordinates)
    score_runs=[]
    for _ in range(2):
        first=point.copy(); first_nll=evaluate(first)['nll']; last=None
        def score(p):
            nonlocal last
            value=evaluate(p);last=dict(coordinates=p.copy(),nll=value['nll'],gradient=value['gradient'].copy())
            return value['gradient']
        fit=least_squares(score,first,bounds=(np.r_[np.full(9,-np.inf),np.log(.05)],np.r_[np.full(9,np.inf),np.log(10.)]),**settings['score_polish'])
        point=fit.x.copy();value=evaluate(point)
        score_runs.append(dict(initial=first,returned=point,initial_nll=first_nll,final_nll=value['nll'],
            gradient=value['gradient'],last_residual_evaluation=last,success=bool(fit.success),status=int(fit.status),
            message=str(fit.message),cost=float(fit.cost),optimality=float(fit.optimality),nfev=int(fit.nfev)))
    value,tighter=evaluate(point),tight(point)
    info=[information_diagnostics(joint,point,relative_step=h).to_dict() for h in settings['information_steps']]
    fd=[audit_joint_gradient(lambda p:evaluate(p)['nll'],joint,point,relative_step=h).to_dict() for h in settings['gradient_steps']]
    refinement=discrepancies(value,tighter)
    checks=dict(terminated=all(r['success'] for r in score_runs),gradient=max(abs(tighter['gradient']))<limits['gradient'],
        curvature=all(r['positive_definite'] and r['condition_number']<limits['information_condition']
                      and r['relative_symmetry_residual']<limits['information_symmetry']
                      and r['newton_correction_supnorm']<limits['newton_correction']
                      and r['standardized_score_supnorm']<limits['standardized_score'] for r in info),
        gradient_fd=max(r['maximum_absolute_difference'] for r in fd)<limits['gradient_fd'],
        integration=max(refinement.values())<limits['integration_refinement']
                    and tighter['numeric_relative_mass_error_sum']<limits['numeric_relative_mass_error']
                    and tighter['tail_relative_mass_bound_sum']<limits['tail_relative_mass'],
        nonworsening=tighter['nll']<=initial['nll']+limits['raw_worsening']
                     and all(r['final_nll']<=r['initial_nll']+limits['raw_worsening'] for r in score_runs),
        restart=max(abs(score_runs[-1]['returned']-score_runs[-1]['initial']))<limits['restart_coordinates'],
        sigma_interior=.05<tighter['sigma']<10.)
    return dict(coordinates=point,initial=initial,value=value,tight=tighter,preliminary=preliminary,
                score_runs=score_runs,information=info,gradient_fd=fd,refinement=refinement,checks=checks,
                local_reference_pass=bool(all(checks.values())))


def run_case(folder, template, data, design, spec):
    folder=Path(folder);problem,original_rows=prepare(template,data,design)
    reference=ObservedPCMReference(problem);count=int(np.count_nonzero(reference.person_weight))
    assert count>0 and problem.observations>0
    dump(folder/'input.json',dict(indices=problem.idx,original_rows=original_rows,x=reference.x,n_person=reference.n,
        person_weight=reference.person_weight,weight_sum=problem.observations,informative_persons=count,
        classification=PROTOCOL['classification'],scientific_inference_ready=False,r_integration=spec['continuous']['r_integration']))
    truth=np.r_[PROTOCOL['true_structural'],np.log(data['condition']['sigma'])]
    logits=reference.offset+reference.design@truth[:-1]+np.exp(truth[-1])*np.array(data['z'])[reference.person,None]*reference.k
    generator_error=float(np.max(abs(np.exp(logits-logsumexp(logits,axis=1,keepdims=True))-np.array(data['generating_probabilities'])[original_rows])))
    stationary=StationarityContract(**dict(spec['gh']['stationarity'],objective_worsening_tolerance_per_observation=1e-18/problem.observations))
    sensitivity=QuadratureSensitivityContract(**spec['gh']['sensitivity'])
    rules={}
    for q in PROTOCOL['orders']:
        z,w=roots_hermitenorm(q);rules[q]=(z,w/np.sqrt(2*np.pi))
    observations={}; runs={}; problems={}
    for q in PROTOCOL['orders']:
        model=replace(problem,quadrature_points=q);problems[q]=model
        for sd in PROTOCOL['starts']:
            key=f'q{q}_sd{sd:g}'
            result=run_free_sd_two_stage(problem.structural_start,sd,model.value,model.value_gradient,
                observations=problem.observations,structural_bounds=problem.structural_bounds,sigma_bounds=(.05,10.),
                joint_value_gradient=model.joint_value_gradient,difference_factory=model.likelihood_difference,
                constraint_residual_function=model.constraint_residual,
                anchor_gradient_limit=spec['gh']['anchor_gradient_limit'],
                preliminary_options=JointPolishOptions(**spec['gh']['preliminary']),refinement_options=JointPolishOptions(**spec['gh']['refinement']))
            dump(folder/(key+'_run.json'),result.to_dict())
            final=result.refinement if result.refinement is not None else result.preliminary
            runs[key]=final;p=np.array(final.restart_polish.joint_coordinates)
            actual=app_evaluate(model,p);finite=reference.finite(p,*rules[q]);cont=reference.continuous(p,**spec['continuous']['tight_integration'])
            parity=discrepancies(actual,finite);assess=assess_free_sd_stationarity(final,stationary)
            replay=[]
            for stage in (result.preliminary,result.refinement):
                if stage is None:continue
                for polish in (stage.primary_polish,stage.restart_polish):
                    shift=polish.objective_shift
                    delta=model.likelihood_difference(np.array(shift.anchor_coordinates)) if shift else None
                    for label,coords,nll,expected in (
                        ('initial',polish.initial_joint_coordinates,polish.initial_objective,shift.initial_difference if shift else 0.),
                        ('returned',polish.joint_coordinates,polish.final_objective,shift.final_difference if shift else 0.),
                        ('last_query',polish.last_query_coordinates,polish.last_query_raw_objective,polish.last_query_optimizer_objective if shift else 0.)):
                        point=np.array(coords);replay.append(dict(point=label,raw_error=abs(model.joint_value_gradient(point)[0]-nll),
                            difference_error=abs(delta(point)-expected) if delta else 0.))
            empty=np.flatnonzero(reference.person_weight==0)
            prior_error=max((max(float(np.max(abs(v['eap'][empty]-reference.x[empty]*p[0]))),float(np.max(abs(v['sd'][empty]-np.exp(p[-1]))))) for v in (actual,finite,cont)),default=0.) if len(empty) else 0.
            observations[key]=dict(coordinates=p,app=actual,finite=finite,continuous_at_same_point=cont,
                same_point_error=discrepancies(finite,cont),parity=parity,assessment=assess.to_dict(),replay=replay,prior_error=prior_error)
            dump(folder/(key+'.json'),observations[key])
            print(folder.name,key,'stationary',assess.stationarity_pass,'SD',np.exp(p[-1]),flush=True)
    continuous={}
    for q in PROTOCOL['orders']:
        key=f'from_q{q}';continuous[key]=continuous_fit(reference,observations[f'q{q}_sd0.5']['coordinates'],spec)
        dump(folder/(key+'.json'),continuous[key])
    cross_cont=float(np.max(abs(continuous['from_q121']['coordinates']-continuous['from_q181']['coordinates'])))
    reference_valid=all(v['local_reference_pass'] for v in continuous.values()) and cross_cont<PROTOCOL['cross_start_limit']
    target=continuous['from_q181']['tight'];target_point=continuous['from_q181']['coordinates'];limits=spec['continuous']['limits']
    comparisons={}
    for key,value in observations.items():
        p=value['coordinates'];total=discrepancies(value['app'],target,keys=('eap','sd'))
        structural=float(np.max(abs(p[:-1]-target_point[:-1])));log_sigma=abs(float(p[-1]-target_point[-1]));loss=(value['continuous_at_same_point']['nll']-target['nll'])/count
        comparisons[key]=dict(structural=structural,log_sigma=log_sigma,continuous_nll_loss_per_person=loss,total_score_error=total,
            pass_targets=bool(reference_valid and value['assessment']['StationarityPass'] and structural<limits['gh_structural']
                and log_sigma<limits['gh_log_sigma'] and total['eap']<limits['gh_eap'] and total['sd']<limits['gh_sd']
                and -limits['raw_worsening']/count<=loss<limits['gh_continuous_nll_loss_per_informative_person']))
    sensitivities={}
    for sd in PROTOCOL['starts']:
        sensitivities[str(sd)]=assess_quadrature_sensitivity(problem_digest=sha(folder/'input.json'),primary_quadrature_points=121,sensitivity_quadrature_points=181,
            primary_run=runs[f'q121_sd{sd:g}'],sensitivity_run=runs[f'q181_sd{sd:g}'],primary_value_function=problems[121].value,
            sensitivity_value_function=problems[181].value,stationarity_contract=stationary,sensitivity_contract=sensitivity).to_dict()
    cross_start={str(q):float(np.max(abs(observations[f'q{q}_sd0.5']['coordinates']-observations[f'q{q}_sd3']['coordinates']))) for q in PROTOCOL['orders']}
    r_cases=dict(gh181=dict(coordinates=observations['q181_sd0.5']['coordinates']),continuous=dict(coordinates=target_point))
    dump(folder/'r_cases.json',r_cases)
    with (folder/'r.log').open('x') as log:
        r_status=subprocess.run(['Rscript',str(Path(__file__).with_suffix('.R')),str(folder/'input.json'),str(folder/'r_cases.json'),str(folder)],stdout=log,stderr=subprocess.STDOUT).returncode
    r_comparisons={}
    if r_status==0:
        for name,value in [('gh181',observations['q181_sd0.5']['continuous_at_same_point']),('continuous',target)]:
            outputs=json.loads((folder/(name+'_r.json')).read_text())['continuous']
            difference=discrepancies(value,outputs[-1],keys=('nll','eap','sd'))
            difference['log_sigma_score']=abs(value['gradient'][-1]-outputs[-1]['log_sigma_nll_score'])
            refinement=discrepancies(*outputs,keys=('nll','eap','sd'))
            r_comparisons[name]=dict(difference=difference,refinement=refinement,passed=all(v<limits['r_'+k] for k,v in difference.items()) and max(refinement.values())<limits['integration_refinement'])
    decimal=None
    if design=='missing_weighted':
        helper=runpy.run_path(str(TEST_SOURCE))['_decimal_observed_nll'];model=problems[181]
        anchor=np.array(runs['q181_sd0.5'].primary_polish.objective_shift.anchor_coordinates);point=observations['q181_sd0.5']['coordinates'];values={}
        for digits in (60,90):
            with localcontext() as context:
                context.prec=digits;left=helper(model,anchor,digits);right=helper(model,point,digits)
                values[str(digits)]=dict(anchor=str(left),final=str(right),delta=str(right-left))
        with localcontext() as context:
            context.prec=90;exact=Decimal(values['90']['delta']);error=abs(Decimal.from_float(model.likelihood_difference(anchor)(point))-exact)
            precision=max(abs(Decimal(values['60'][k])-Decimal(values['90'][k])) for k in ('anchor','final','delta'))
            passed=error<=Decimal('1e-18')+Decimal('1e-11')*abs(exact) and precision<Decimal('1e-45')
        decimal=dict(values=values,error=str(error),precision=str(precision),passed=bool(passed));dump(folder/'decimal.json',decimal)
    checks=dict(generator=generator_error<1e-14,arithmetic=all(all(v['parity'][k]<limit for k,limit in PROTOCOL['finite_parity'].items()) for v in observations.values()),
        replay=all(max(r['raw_error'],r['difference_error'])<PROTOCOL['replay_limit'] for v in observations.values() for r in v['replay']),
        prior=all(v['prior_error']<1e-10 for v in observations.values()),continuous_reference=reference_valid,
        independent_r=r_status==0 and len(r_comparisons)==2 and all(v['passed'] for v in r_comparisons.values()),
        decimal=decimal is None or decimal['passed'],raw_reconstruction=all(v['objective_reconstruction_pass'] for v in sensitivities.values()))
    summary=dict(classification=PROTOCOL['classification'],scientific_inference_ready=False,qualification_eligible=False,
        rows=len(original_rows),weight_sum=problem.observations,informative_persons=count,generator_error=generator_error,
        observations=observations,continuous=continuous,continuous_cross_start=cross_cont,comparisons=comparisons,
        sensitivities=sensitivities,cross_start=cross_start,r_exit_code=r_status,independent_r=r_comparisons,
        checks=checks,implementation_checks_pass=bool(all(checks.values())),
        artifact_sha256={p.name:sha(p) for p in sorted(folder.iterdir()) if p.is_file()})
    dump(folder/'summary.json',summary)
    print(folder.name,'complete','checks',checks,'total accuracy',{k:v['pass_targets'] for k,v in comparisons.items()},flush=True)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True)
    out=parser.parse_args().output_dir.resolve();load=lambda p:json.loads(p.read_text())
    assert sha(PREVIOUS/'summary.json')=='106b21b62b07042a27497210253e23317ab558fe072a4ceb2fb7406fefa7ee5a'
    assert sha(GH_PREVIOUS/'summary.json')=='9204597c2128374f36793387d3cd9feaab02988a7b89c66fcdda98b39735be99'
    assert sha(GH_PREVIOUS/'protocol.json')==load(GH_PREVIOUS/'summary.json')['protocol_sha256']
    previous=load(PREVIOUS/'summary.json');previous_spec=load(PREVIOUS/'protocol.json')
    for path,digest in previous['artifact_sha256'].items():assert sha(PREVIOUS/path)==digest,path
    for path,digest in previous_spec['source_sha256'].items():assert sha(ROOT/path)==digest,path
    sources=[Path(__file__).resolve(),Path(__file__).with_suffix('.R'),TEST_SOURCE,ROOT/'validation/mml_free_sd_stationarity_adapter.py',ROOT/'mfrm_app/mml_engine_v2.py',ROOT/'mfrm_app/mml_quadrature_sensitivity.py',GH_PREVIOUS/'protocol.json',PREVIOUS/'summary.json']
    hashes=dict(previous_spec['source_sha256'],**{str(p.relative_to(ROOT)):sha(p) for p in sources})
    spec=dict(PROTOCOL,source_sha256=hashes,continuous=previous_spec,gh=load(GH_PREVIOUS/'protocol.json'),
              environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    out.mkdir(parents=True,exist_ok=False);dump(out/'protocol.json',spec)
    datasets={}
    for condition in PROTOCOL['datasets']:
        data=generate(condition);datasets[condition['id']]=data;dump(out/(condition['id']+'_full_input.json'),data)
    helpers=runpy.run_path(str(TEST_SOURCE));template=helpers['_observed_problem'](helpers['native_free_sd_result'].__wrapped__(),'PCM')
    jobs=[]
    for name,data in datasets.items():
        for design in PROTOCOL['designs']:
            folder=out/(name+'_'+design);folder.mkdir();jobs.append((folder,template,data,design,spec))
    results={};failures={}
    with ProcessPoolExecutor(max_workers=PROTOCOL['workers']) as executor:
        pending={executor.submit(run_case,*job):job[0] for job in jobs}
        for future in as_completed(pending):
            folder=pending[future]
            try:results[folder.name]=future.result()
            except Exception as exc:
                failure=dict(error=repr(exc),traceback=traceback.format_exc());failures[folder.name]=failure;dump(folder/'failure.json',failure);print(folder.name,'FAILED',repr(exc),flush=True)
    checks=dict(all_six_retained=len(results)==6 and not failures,
        implementation=all(v['implementation_checks_pass'] for v in results.values()) and len(results)==6,
        sources_unchanged=all(sha(ROOT/p)==h for p,h in hashes.items()),inference_withheld=all(not v['scientific_inference_ready'] for v in results.values()))
    dump(out/'summary.json',dict(classification=PROTOCOL['classification'],scientific_inference_ready=False,qualification_eligible=False,
        checks=checks,implementation_checks_pass=bool(all(checks.values())),results=results,failures=failures,
        protocol_sha256=sha(out/'protocol.json'),artifact_sha256={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}))
    print('Study checks:',checks,flush=True)
    return 0 if all(checks.values()) else 1


if __name__=='__main__':
    raise SystemExit(main())
