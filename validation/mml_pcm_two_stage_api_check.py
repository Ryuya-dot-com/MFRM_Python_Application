#!/usr/bin/env python3
"""Check the optional two-stage API against frozen PCM results and raw-Q gates."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import platform
import subprocess

import mpmath as mp
import numpy as np
import scipy

from mml_pcm_continuous_refit import ROOT, dump, sha
from mml_pcm_independent_conditions import PCMIntegral, difference
from mml_pcm_conquest_history_audit import finite_gh
from mml_pcm_gh_likelihood_difference import gh_likelihood_difference, mp_gh
from mfrm_app.mml_stationarity import JointPolishOptions
from mfrm_app.mml_engine_v2 import StationarityContract, run_free_sd_two_stage, assess_free_sd_stationarity
from mfrm_app.mml_quadrature_sensitivity import QuadratureSensitivityContract, assess_quadrature_sensitivity

REFERENCE=ROOT/'validation/generated/mml_pcm_gh_two_stage_20260914'
SNAPSHOT=ROOT/'validation/generated/mml_api_pre_shift_20260914'
STATIONARITY=StationarityContract(max_projected_gradient_supnorm=1e-8,
    max_standardized_score_supnorm=1e-7,max_newton_correction_supnorm=1e-7,
    max_restart_improvement_per_observation=1e-10,max_restart_displacement=1e-7,
    max_gradient_fd_disagreement=1e-7,max_objective_value_disagreement=1e-8,
    max_information_relative_symmetry_residual=1e-6,max_information_condition_number=1e6,
    max_constraint_residual=1e-12,objective_worsening_tolerance_total=1e-18,
    objective_worsening_tolerance_per_observation=1e-18/80,sigma_boundary_tolerance=1e-6)
SENSITIVITY=QuadratureSensitivityContract(primary_quadrature_points=121,sensitivity_quadrature_points=181,
    max_structural_parameter_difference=1e-4,max_log_sigma_difference=1e-4,
    max_sensitivity_optimization_gain_per_observation=1e-6,max_negative_sensitivity_gain_per_observation=1e-10,
    max_abs_primary_back_evaluation_change_per_observation=1e-6,max_objective_reconstruction_disagreement=1e-8)
PROTOCOL=dict(id='pcm_two_stage_api_20260914_v1',classification='DEVELOPMENT_ONLY',
    qualification_eligible=False,scientific_inference_ready=False,
    question='Does the optional API preserve raw NLL, stage failures and the existing Q-sensitivity reconstruction while reproducing the validated GH solutions?',
    datasets=['ordinary_new_1','wide_new_1'],orders=[121,181],starts=[.5,3.],
    preliminary=dict(maxiter=250,gtol=1e-8,ftol=1e-15,maxls=50,log_sigma_relative_step=1e-6),
    refinement=dict(maxiter=250,gtol=1e-8,ftol=0.,maxls=50,log_sigma_relative_step=1e-6),
    anchor_gradient_limit=1e-3,stationarity=asdict(STATIONARITY),sensitivity=asdict(SENSITIVITY),
    parity_limits=dict(coordinates=1e-7,nll=1e-8,eap=1e-7,sd=1e-7),
    independent_selection='Start SD 0.5, both orders, both datasets: R at both finite rules and continuous integration; 60/90 digits at API anchor and final point',
    mp_limits=dict(absolute=1e-18,relative=1e-11,refinement=1e-45),
    r_limits=dict(nll=1e-8,eap=1e-9,sd=1e-9),
    failure_policy='Freeze before fitting; retain all failures without retuning. Old source manifests stay unchanged and are checked against the snapshot where needed.',
    scope='Optional numerical API integration on observed complete PCM; no default Streamlit wiring, native TAM/ConQuest fits or scientific inference qualification')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    out=parser.parse_args().output_dir.resolve();load=lambda p:json.loads(p.read_text())
    assert sha(REFERENCE/'summary.json')=='a58dee9cbf1c3b12f51a1088ae3614a31a37ab09f3771c5eb54764790939a929'
    old=load(REFERENCE/'summary.json');old_spec=load(REFERENCE/'protocol.json')
    for name,h in old['artifact_sha256'].items():assert sha(REFERENCE/name)==h,name
    snapshot=load(SNAPSHOT/'manifest.json')
    for name,v in snapshot['sources'].items():assert sha(SNAPSHOT/v['snapshot'])==v['original_sha256'],name
    for name,h in old_spec['source_sha256'].items():
        archived=SNAPSHOT/'source'/Path(name).relative_to(ROOT)
        assert sha(archived)==h,name
    paths=[Path(__file__).resolve(),ROOT/'mfrm_app/mml_stationarity.py',ROOT/'mfrm_app/mml_engine_v2.py',
        ROOT/'mfrm_app/mml_quadrature_sensitivity.py',SNAPSHOT/'manifest.json']+[
        Path(__file__).with_name('mml_pcm_'+f) for f in ('gh_likelihood_difference.py','likelihood_difference.py',
        'conquest_history_audit.py','independent_conditions.py','continuous_refit.py','independent_conditions.R')]
    hashes={str(p):sha(p) for p in paths};out.mkdir(parents=True,exist_ok=False)
    dump(out/'protocol.json',dict(PROTOCOL,source_sha256=hashes,reference_summary_sha256=sha(REFERENCE/'summary.json'),
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__)))
    rules={q:load(REFERENCE/'rules.json')[str(q)] for q in PROTOCOL['orders']}
    results={};mp_checks=[];r_jobs=[]
    for name in PROTOCOL['datasets']:
        folder=out/name;folder.mkdir();inp=load(REFERENCE/name/'input.json');inp['orders']=PROTOCOL['orders']
        dump(folder/'input.json',inp);problem=PCMIntegral(inp['data']);records={};runs={};evaluators={};r_cases={}
        for q in PROTOCOL['orders']:
            z,w=np.array(rules[q]['nodes']),np.array(rules[q]['weights'])
            evaluate=lambda p,z=z,w=w:finite_gh(problem,p,z,w)
            raw=lambda p,sd,evaluate=evaluate:evaluate(np.r_[p,np.log(sd)])['nll']
            def structural(p,sd,evaluate=evaluate):
                v=evaluate(np.r_[p,np.log(sd)]);return v['nll'],v['gradient'][:-1]
            def joint(p,evaluate=evaluate):
                v=evaluate(p);return v['nll'],v['gradient']
            evaluators[q]=raw
            for sd in PROTOCOL['starts']:
                key=f'q{q}_sd{sd:g}'
                run=run_free_sd_two_stage(np.zeros(23),sd,raw,structural,observations=len(problem.y),
                    joint_value_gradient=joint,difference_factory=lambda a:gh_likelihood_difference(problem,a,z,w),
                    anchor_gradient_limit=PROTOCOL['anchor_gradient_limit'],
                    preliminary_options=JointPolishOptions(**PROTOCOL['preliminary']),
                    refinement_options=JointPolishOptions(**PROTOCOL['refinement']),constraint_residual_function=lambda p:0.)
                dump(folder/(key+'_run.json'),run.to_dict())
                final=run.refinement if run.refinement is not None else run.preliminary
                p=np.array(final.restart_polish.joint_coordinates);value=evaluate(p)
                continuous=problem.evaluate(p,**inp['integration'][-1])
                expected=old['results'][name]['records'][key]
                error=dict(difference(value,expected['finite']),coordinates=float(np.max(abs(p-expected['coordinates']))))
                assessment=assess_free_sd_stationarity(final,STATIONARITY)
                record=dict(anchor_admitted=run.anchor_admitted,assessment=assessment.to_dict(),
                    parity_error=error,parity_pass=all(error[k]<v for k,v in PROTOCOL['parity_limits'].items()),
                    finite=value,continuous=continuous,integration_error=difference(value,continuous),
                    finite_by_rule={str(other):finite_gh(problem,p,np.array(rule['nodes']),np.array(rule['weights'])) for other,rule in rules.items()})
                records[key]=record;runs[key]=final;dump(folder/(key+'.json'),record)
                print(name,key,'stationarity',assessment.stationarity_pass,'parity',record['parity_pass'],flush=True)
                if sd==.5:
                    r_cases[f'q{q}']=dict(coordinates=p)
                    a=np.array(final.primary_polish.objective_shift.anchor_coordinates);precise={}
                    for digits in (60,90):
                        with mp.workdps(digits):
                            left,right=mp_gh(problem,a,z,w)['nll'],mp_gh(problem,p,z,w)['nll']
                            precise[str(digits)]=dict(anchor=mp.nstr(left,digits),final=mp.nstr(right,digits),delta=mp.nstr(right-left,digits))
                        dump(folder/f'q{q}_mp_{digits}.json',precise[str(digits)])
                    with mp.workdps(90):
                        exact=mp.mpf(precise['90']['delta']);delta=final.restart_polish.objective_shift.final_difference
                        err=abs(mp.mpf(delta)-exact);limit=PROTOCOL['mp_limits']['absolute']+PROTOCOL['mp_limits']['relative']*abs(exact)
                        precision=max(abs(mp.mpf(precise['60'][k])-mp.mpf(precise['90'][k])) for k in ('anchor','final','delta'))
                    mp_checks.append(dict(dataset=name,q=q,error=float(err),precision=float(precision),passed=bool(err<=limit and precision<PROTOCOL['mp_limits']['refinement'])))
        sensitivities={}
        for sd in PROTOCOL['starts']:
            common=dict(problem_digest=sha(folder/'input.json'),primary_quadrature_points=121,sensitivity_quadrature_points=181,
                primary_run=runs[f'q121_sd{sd:g}'],sensitivity_run=runs[f'q181_sd{sd:g}'],
                stationarity_contract=STATIONARITY,sensitivity_contract=SENSITIVITY)
            assessment=assess_quadrature_sensitivity(primary_value_function=evaluators[121],sensitivity_value_function=evaluators[181],**common)
            wrong={q:lambda p,sigma,q=q,origin=runs[f'q{q}_sd{sd:g}'].restart_polish.objective_shift.anchor_objective:
                   evaluators[q](p,sigma)-origin for q in rules}
            rejected=assess_quadrature_sensitivity(primary_value_function=wrong[121],sensitivity_value_function=wrong[181],**common)
            sensitivities[str(sd)]=dict(assessment=assessment.to_dict(),shifted_instead_of_raw_rejected=not rejected.objective_reconstruction_pass)
        dump(folder/'r_cases.json',r_cases);log=(folder/'r.log').open('x')
        process=subprocess.Popen(['Rscript',str(paths[-1]),str(folder/'input.json'),str(folder/'r_cases.json'),str(folder)],stdout=log,stderr=subprocess.STDOUT)
        r_jobs.append((name,process,log));results[name]=dict(records=records,sensitivities=sensitivities)
    independent={}
    for name,process,log in r_jobs:
        status=process.wait();log.close();comparisons={}
        dump(out/name/'r_execution.json',dict(exit_code=status))
        if status==0:
            for q in rules:
                r=load(out/name/f'q{q}_r.json');py=results[name]['records'][f'q{q}_sd0.5']
                comparisons[str(q)]=dict(continuous=difference(py['continuous'],r['continuous'][-1]),
                    **{f'q{other}':difference(py['finite_by_rule'][str(other)],r['finite'][str(other)]) for other in rules})
        independent[name]=dict(exit_code=status,comparisons=comparisons)
    rows=[r for d in results.values() for r in d['records'].values()]
    checks=dict(all_eight_retained=len(rows)==8,parity=all(r['parity_pass'] for r in rows),
        stationarity=all(r['anchor_admitted'] and r['assessment']['StationarityPass'] for r in rows),
        mp=all(c['passed'] for c in mp_checks),
        r=all(v['exit_code']==0 and len(v['comparisons'])==2 and all(e[k]<limit for c in v['comparisons'].values() for e in c.values() for k,limit in PROTOCOL['r_limits'].items()) for v in independent.values()),
        raw_reconstruction=all(v['assessment']['objective_reconstruction_pass'] and v['shifted_instead_of_raw_rejected'] for d in results.values() for v in d['sensitivities'].values()),
        inference_withheld=all(not r['assessment']['InferenceReady'] for r in rows) and all(not v['assessment']['scientific_inference_ready'] for d in results.values() for v in d['sensitivities'].values()),
        sources_unchanged=all(sha(p)==h for p,h in hashes.items()))
    dump(out/'summary.json',dict(classification=PROTOCOL['classification'],scientific_inference_ready=False,qualification_eligible=False,
        checks=checks,api_checks_pass=all(checks.values()),results=results,mp_checks=mp_checks,independent_r=independent,
        protocol_sha256=sha(out/'protocol.json'),artifact_sha256={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}))
    print('API checks:',checks,flush=True)
    return 0 if all(checks.values()) else 1


if __name__=='__main__':
    raise SystemExit(main())
