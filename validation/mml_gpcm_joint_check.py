"""Independent finite-GH GPCM joint-information development controls.

Preserves every point and failure. Uses SLSQP to keep ALL expanded log-slopes
bounded; the RSM/PCM v2 solver's support contract is not broadened. No sampling
uncertainty, continuous-integral accuracy or G4 qualification is asserted.
"""
import argparse
import copy
from dataclasses import replace
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import zipfile

import numpy as np
import pandas as pd
from scipy.optimize import LinearConstraint, minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import streamlit_app as app
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data
from validation.mml_joint_information_check import save
from validation.mml_rasch_quadrature_review import verify


def prepare_problem(data, q, *, rating_min=0, rating_max=3):
    """Input-only native coordinates; separate development scope for GPCM."""
    prep, pcm = prepare_app_free_sd_problem_from_data(data,person_col='Person',
        facet_cols=['Rater','Task','Criterion'],score_col='Score',rating_min=rating_min,rating_max=rating_max,
        model='PCM',noncenter_facet='Criterion',step_facet='Criterion',quadrature_points=q)
    config=copy.deepcopy(pcm.config)
    config.update(model='GPCM',slope_facet='Criterion',gpcm_log_slope_bounds=(-3.,3.))
    sizes=app.build_param_sizes(config)
    config['facet_regularization']=app.prepare_facet_regularization(None,config,sizes)['runtime']
    problem=replace(pcm, config=config,sizes=sizes,
        idx=copy.deepcopy(app.build_indices(prep,step_facet='Criterion',slope_facet='Criterion')),
        structural_start=tuple(pcm.structural_start)+(0.,)*sizes['log_slopes'],
        structural_bounds=tuple(app.build_optimizer_bounds(sizes,config)))
    return prep,problem


def joint_constraints(problem):
    """Append a zero SD coefficient to the native derived-log-slope constraint."""
    return [LinearConstraint(np.column_stack([np.atleast_2d(c.A),np.zeros((np.atleast_2d(c.A).shape[0],1))]),c.lb,c.ub)
            for c in app.build_optimizer_constraints(problem.sizes,problem.config)]


def expanded(problem, point):
    """Native output vector; covariance Jacobian is computed explicitly below."""
    p=app.expand_params(point[:-1],problem.sizes,problem.config)
    return np.r_[np.concatenate([p['facets'][f] for f in problem.config['facet_names']]),
                  p['steps_mat'].ravel(),p['log_slopes'],p['slopes'],np.exp(point[-1])]


def expansion_jacobian(config,point):
    """Independent exact expansion, including dependent log-slopes and exp map."""
    labels=[]; rows=[]; dimension=len(point); cursor=0
    def block(prefix,count,centered):
        nonlocal cursor
        free=count-int(centered); basis=np.zeros((count,dimension))
        basis[:free,cursor:cursor+free]=np.eye(free)
        if centered: basis[-1,cursor:cursor+free]=-1.
        cursor+=free;labels.extend(f'{prefix}:{i}' for i in range(count));rows.extend(basis)
        return basis
    for facet in config['facet_names']:
        block(facet,len(config['facet_levels'][facet]),facet!=config['noncenter_facet'])
    for level in config['facet_levels'][config['step_facet']]:
        block('step:'+str(level),config['n_cat']-1,True)
    slope_count=len(config['facet_levels'][config['slope_facet']])
    slope_basis=block('log_slope',slope_count,True)
    slopes=np.exp(slope_basis@point)
    rows.extend(slopes[:,None]*slope_basis);labels.extend(f'slope:{i}' for i in range(slope_count))
    assert cursor==dimension-1
    sd=np.zeros(dimension);sd[-1]=np.exp(point[-1]);rows.append(sd);labels.append('population_sd')
    return labels,np.array(rows)


def generate(seed,design,persons):
    """Independent four-category DGP with genuine step-shape/slope variation."""
    rng=np.random.default_rng(seed); theta=rng.normal(0.,1.2,persons)
    rater,task,criterion=[-.4,0.,.4],[-.2,.2],[-.4,.1,.7]
    steps=np.array([[-.8,-.1,.9],[-1.,.2,.8],[-.5,-.2,.7]])
    slopes=[.7,1.,1/.7]; rows=[]
    for p in range(persons):
        for r in range(3):
            for t in range(2):
                for c in range(3):
                    if design=='sparse' and (r==p%3 or t!=p%2 or c==(p+1)%3): continue
                    if design=='unbalanced' and (r,t,c)!=(p%3,p%2,p%3):
                        if rng.random()>[.25,.6,.95][r]*[.45,.8,1.][c]: continue
                    logits=slopes[c]*((theta[p]-rater[r]-task[t]-criterion[c])*np.arange(4)-np.r_[0.,np.cumsum(steps[c])])
                    prob=np.exp(logits-logits.max());prob/=prob.sum()
                    rows.append((f'P{p:03d}',f'R{r}',f'T{t}',f'C{c}',rng.choice(4,p=prob)))
    return pd.DataFrame(rows,columns=['Person','Rater','Task','Criterion','Score'])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    parser.add_argument('--refinement-from',type=Path)
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    contract=dict(classification='DEVELOPMENT_ONLY',scientific_inference_ready=False,
        persons=48,designs=['regular','sparse','unbalanced'],seeds=[2026100440,2026100441,2026100442],
        quadrature_points=[31,61],initialization='zero facets/log-slopes; zero-sum [-1,0,1] steps; sigma=1',
        optimizer=[dict(method='SLSQP',maxiter=500,ftol=1e-12),dict(method='SLSQP',maxiter=500,ftol=1e-14)],
        bounds=dict(expanded_log_slopes=[-3.,3.],sigma=[.05,10.]),
        derivative_steps=[1e-5,5e-6],inversion='direct positive-definite inverse only; no flooring',
        tolerances=dict(nll=1e-8,gradient=2e-5,information=1e-4,information_step=1e-4,covariance=1e-5,
            expanded=1e-10,expanded_jacobian=1e-8,expanded_covariance=1e-5,eap=1e-10,sd=1e-10,posterior=1e-10),
        fit_checks=dict(score_supnorm=1e-4,constraint_residual=1e-10,boundary_distance=1e-6),
        scope='same-point finite-GH joint probability/score/information/coordinate transformation; not continuous or sampling accuracy')
    parent=args.refinement_from
    if parent:
        verify(parent)
        previous=json.loads((parent/'summary.json').read_text())
        if previous['completed']!=6 or previous['failures']:
            raise ValueError('Refinement requires all six retained parent points')
        contract.update(quadrature_points=[121,241,321],
            initialization='Q61 parent point for all orders; additional Q321 from Q31 for every design',
            parent_manifest_sha256=hashlib.sha256((parent/'manifest.json').read_bytes()).hexdigest(),
            cross_start_tolerance=1e-6)
    contract['planned_cases']=12 if parent else 6
    save(out/'contract.json',contract)
    source_paths=[ROOT/'streamlit_app.py',*sorted((ROOT/'mfrm_app').glob('*.py')),
        Path(__file__),ROOT/'validation/mml_free_sd_stationarity_adapter.py',
        ROOT/'validation/mml_joint_information_check.py',ROOT/'validation/mml_joint_information_reference.R']
    with zipfile.ZipFile(out/'sources.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for path in source_paths:archive.write(path,str(path.relative_to(ROOT)))
    save(out/'source_manifest.json',{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths})
    (out/'reference.R').write_bytes((ROOT/'validation/mml_joint_information_reference.R').read_bytes())
    save(out/'runtime.json',dict(python=sys.version,packages={d.metadata['Name']:d.version for d in importlib.metadata.distributions()}))
    problems=[];prefit=[]
    for design,seed in zip(contract['designs'],contract['seeds'],strict=True):
        case_id='GPCM_'+design
        data=pd.read_csv(parent/(case_id+'.csv')) if parent else generate(seed,design,contract['persons'])
        data.to_csv(out/(case_id+'.csv'),index=False)
        prep,problem=prepare_problem(data,31)
        prefit.append(dict(id=case_id,seed=seed,rows=len(data),config=problem.config,
            start=problem.structural_start,sigma_start=problem.sigma_start,
            input_sha256=hashlib.sha256((out/(case_id+'.csv')).read_bytes()).hexdigest(),
            bounds=problem.structural_bounds,linear_constraint=joint_constraints(problem)[0].A))
        problems.append((case_id,data,problem))
    # Config contains NaN free-anchor sentinels: store a compact declared subset.
    for p in prefit:
        p['config']={k:p['config'][k] for k in ['model','facet_names','facet_levels','facet_signs','noncenter_facet','step_facet','slope_facet','n_cat','rating_min','rating_max']}
    save(out/'prefit_manifest.json',prefit)
    python_cases=[];references=[];failures=[]
    for case_id,data,base in problems:
        jobs=[(q,'from61' if parent else '') for q in contract['quadrature_points']]
        if parent: jobs.append((321,'from31'))
        for q,start_label in jobs:
            run_id=f'{case_id}_Q{q}'+('_'+start_label if start_label else '')
            try:
                problem=replace(base,quadrature_points=q)
                start=np.r_[problem.structural_start,np.log(problem.sigma_start)]
                if parent:
                    start=np.array(json.loads((parent/(case_id+'_Q'+start_label.removeprefix('from')+'.json')).read_text())['point'])
                bounds=list(problem.structural_bounds)+[tuple(np.log(problem.sigma_bounds))]
                runs=[]
                for options in contract['optimizer']:
                    fit=minimize(problem.joint_value_gradient,start,jac=True,method=options['method'],
                        bounds=bounds,constraints=joint_constraints(problem),
                        options={k:v for k,v in options.items() if k!='method'})
                    runs.append(dict(initial=start,point=fit.x,nll=fit.fun,gradient=fit.jac,
                        success=bool(fit.success),status=int(fit.status),message=str(fit.message),nit=fit.nit))
                    start=fit.x
                point=fit.x; sigma=float(np.exp(point[-1]));quad=app.gauss_hermite_normal(q,sd=sigma)
                information=[app._compute_mml_observed_information(point,problem.idx,problem.config,problem.sizes,quad,
                    rel_step=h,include_log_sigma=True) for h in contract['derivative_steps']]
                eigen=np.linalg.eigvalsh(information[-1]);positive=bool(eigen.min()>0)
                covariance=np.linalg.inv(information[-1]) if positive else None
                p=app.expand_params(point[:-1],problem.sizes,problem.config)
                posterior,ll=app._e_step_posteriors(problem.idx,problem.config,p,quad)
                person=app.compute_person_eap(problem.idx,problem.config,p,quad)
                labels,jacobian=expansion_jacobian(problem.config,point)
                bound_distance=min(np.min(3-abs(p['log_slopes'])),point[-1]-np.log(.05),np.log(10.)-point[-1])
                record=dict(id=run_id,runs=runs,point=point,nll=problem.joint_value_gradient(point)[0],
                    gradient=problem.joint_value_gradient(point)[1],information=information[-1],
                    information_step_disagreement=float(np.max(abs(information[0]-information[1]))),
                    positive_definite=positive,eigenvalues=eigen,covariance=covariance,
                    expanded=expanded(problem,point),expanded_labels=labels,expanded_jacobian=jacobian,
                    expanded_covariance=jacobian@covariance@jacobian.T if positive else None,
                    eap=person.Estimate.to_numpy(),sd=person.SD.to_numpy(),posterior=posterior,
                    constraint_residual=problem.constraint_residual(point[:-1]),boundary_distance=float(bound_distance))
                save(out/(run_id+'.json'),record);python_cases.append(record)
                config={k:problem.config[k] for k in ['model','facet_names','facet_levels','noncenter_facet','facet_signs',
                    'step_facet','slope_facet','n_cat','rating_min','rating_max']}
                references.append(dict(id=run_id,config=config,data=data.to_dict('records'),point=point,
                    nodes=app.gauss_hermite_normal(q)['nodes'],weights=quad['weights']))
                print(run_id,'Python complete; score=',np.max(abs(record['gradient'])),flush=True)
            except Exception as exc:
                failures.append(dict(id=run_id,type=type(exc).__name__,message=str(exc)))
                save(out/'failures.json',failures)
                print(run_id,'failed',repr(exc),flush=True)
    save(out/'reference_input.json',dict(classification='DEVELOPMENT_ONLY',cases=references))
    with (out/'R_execution.log').open('w') as log:
        result=subprocess.run(['/usr/local/bin/Rscript',str(out/'reference.R'),str(out/'reference_input.json'),str(out/'R_results.json')],stdout=log,stderr=subprocess.STDOUT)
    comparisons=[]
    if result.returncode==0:
        reference=json.loads((out/'R_results.json').read_text());by_id={c['id']:c for c in reference['cases']}
        for case in python_cases:
            r=by_id[case['id']]
            errors={k:float(np.max(abs(np.array(case[k])-np.array(r[k])))) for k in
                    ['nll','gradient','information','expanded','expanded_jacobian','posterior','eap','sd']}
            errors['information_step']=max(case['information_step_disagreement'],r['hessian_step_disagreement'])
            for k in ['covariance','expanded_covariance']:
                errors[k]=float(np.max(abs(np.array(case[k])-np.array(r[k])))) if case[k] is not None and r[k] is not None else None
            checks={k:v is not None and v<=contract['tolerances'][k] for k,v in errors.items()}
            fit_checks=dict(solver_terminated=case['runs'][-1]['success'],score=np.max(abs(np.array(case['gradient'])))<contract['fit_checks']['score_supnorm'],
                constraints=case['constraint_residual']<=contract['fit_checks']['constraint_residual'],
                interior=case['boundary_distance']>contract['fit_checks']['boundary_distance'])
            comparisons.append(dict(id=case['id'],errors=errors,checks=checks,fit_checks=fit_checks,
                                    all_checks_pass=all(checks.values()) and all(fit_checks.values())))
    cross_starts={}
    if parent:
        by_id={c['id']:c for c in python_cases}
        for design in contract['designs']:
            a,b=(by_id.get(f'GPCM_{design}_Q321_{label}') for label in ['from31','from61'])
            cross_starts[design]=float(np.max(abs(a['point']-b['point']))) if a is not None and b is not None else None
    summary=dict(classification='DEVELOPMENT_ONLY',scientific_inference_ready=False,planned=contract['planned_cases'],completed=len(python_cases),
        cross_start_coordinate_difference=cross_starts,
        R_returncode=result.returncode,failures=failures,comparisons=comparisons,
        all_checks_pass=len(comparisons)==contract['planned_cases'] and not failures and all(c['all_checks_pass'] for c in comparisons)
            and all(v is not None and v<contract['cross_start_tolerance'] for v in cross_starts.values()))
    save(out/'summary.json',summary)
    save(out/'manifest.json',{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file()})
    print((out/'summary.json').read_text(),flush=True)
    return 0 if summary['all_checks_pass'] else 1


if __name__=='__main__':raise SystemExit(main())
