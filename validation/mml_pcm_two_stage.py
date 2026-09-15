#!/usr/bin/env python3
"""Predeclared two-stage fixed-grid PCM check on fresh, unedited data."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
import json
from pathlib import Path
import platform
import subprocess

import mpmath as mp
import numpy as np
import scipy

from mml_pcm_continuous_refit import ROOT, dump, sha
from mml_pcm_independent_conditions import PCMIntegral, simulate, difference
from mml_pcm_external_review import fixed_grid
from mml_pcm_likelihood_difference import likelihood_difference
from mml_pcm_optimizer_resolution import mp_nll
from mml_pcm_relative_optimizer import run_pair
from mfrm_app.mml_stationarity import central_difference_gradient, information_diagnostics

PREVIOUS=ROOT/'validation/generated/mml_pcm_relative_optimizer_20260914'
SIMULATION=ROOT/'validation/mml_pcm_independent_protocol.json'
PROTOCOL=dict(
    id='pcm_two_stage_20260914_v1',classification='DEVELOPMENT_ONLY',
    qualification_eligible=False,scientific_inference_ready=False,
    question='Can a run construct its own fixed anchor and attain the requested gradient on new PCM data, without hiding preliminary failures or integration errors?',
    datasets=[dict(id='ordinary_new_1',seed=2026091411,sigma=1.),
              dict(id='ordinary_new_2',seed=2026091412,sigma=1.),
              dict(id='wide_new_1',seed=2026091413,sigma=2.5),
              dict(id='wide_new_2',seed=2026091414,sigma=2.5)],
    simulation='Reuse the frozen 80-person PCM design and PCG64 inverse-CDF generator, with new seeds; no response editing, selection or redraws. Truth is archived but never used in fitting.',
    starts=[.5,1.,3.],structural_start='All 23 free structural coordinates zero',
    grids=[dict(id='b12_q241',bound=12,q=241,spacing=.1),
           dict(id='b20_q801',bound=20,q=801,spacing=.05)],
    preliminary=dict(gtol=1e-8,ftol=1e-15,maxiter=250,maxls=50,maxcor=10,maxfun=15000),
    refinement=dict(gtol=1e-8,ftol=0.,maxiter=250,maxls=50,maxcor=10,maxfun=15000),
    passes_per_stage=2,planned_trajectories=24,maximum_optimizer_calls=96,sigma_bounds=[.05,10.],
    anchor='The second preliminary return, regardless of raw success, provided it is finite, interior and has raw gradient <=1e-3. No best-start or best-iteration selection. Freeze it through refinement and restart.',
    admission_gradient=1e-3,
    refusal='For a usable finite endpoint that misses admission, omit refinement and mark its checks false. Unusable arithmetic aborts with saved inputs/prior outputs retained; do not redraw or change settings.',
    reporting='Keep preliminary success, refinement success, requested-gradient checks, local curvature, restart and integration separately. Neither a successful refinement nor an execution flag overwrites preliminary failure or grants inference.',
    diagnostic_limits=dict(gradient=1e-8,gradient_fd=1e-7,newton=1e-7,standardized_score=1e-7,
        condition=1e6,relative_symmetry=1e-6,restart_coordinates=1e-7,
        refinement_worsening=1e-18,scalar_reconstruction=1e-8,continuous_refinement=1e-8,
        continuous_error_bound=1e-8,r_nll=1e-8,r_moments=1e-9,
        mp_absolute=1e-18,mp_relative=1e-11,mp_refinement=1e-45),
    integration=[dict(bound=12,rel_tol=1e-10,abs_tol=1e-12),
                 dict(bound=14,rel_tol=1e-12,abs_tol=1e-13)],
    integration_targets=dict(nll_per_person=1e-6,eap=1e-4,sd=1e-4),
    independent_selection='Start SD 1, both grids, all four new datasets; independent R finite-grid and continuous moments, and 60/90-digit NLL at preliminary anchor and final endpoint. Retain failed endpoints too.',
    scope='Prospectively specified numerical development check, not recovery/coverage qualification. Complete unit-weight PCM, mean zero, unnormalized fixed-grid normal prior; no app/native external fits or moving GH nodes.',
)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    out=parser.parse_args().output_dir.resolve()
    load=lambda path:json.loads(path.read_text())
    assert sha(PREVIOUS/'summary.json')=='b5c3a17b91d193090e1ac810691a2d02d22e1c07dc310de1333c96cb14365eff'
    assert sha(SIMULATION)=='b732c535860e554e9b6a3028ac28c3fd41c2992d48baa87c06517e5592a76b05'
    previous=load(PREVIOUS/'summary.json')
    for name,h in previous['artifact_sha256'].items():
        assert sha(PREVIOUS/name)==h,name
    for name,h in load(PREVIOUS/'input.json')['source_sha256'].items():
        assert sha(name)==h,name
    simulation=load(SIMULATION)
    assert not ({x['seed'] for x in simulation['datasets']} & {x['seed'] for x in PROTOCOL['datasets']})
    paths=[Path(__file__).resolve(),SIMULATION,ROOT/'mfrm_app/mml_stationarity.py',
        ROOT/'mfrm_app/mml_engine_v2.py',ROOT/'validation/mml_pcm_grid_check.R']+[
        Path(__file__).with_name('mml_pcm_'+n+'.py') for n in
        ['relative_optimizer','continuous_refit','independent_conditions','external_review','likelihood_difference','optimizer_resolution']]
    hashes={str(p):sha(p) for p in paths}
    out.mkdir(parents=True,exist_ok=False)
    # Freeze settings and source hashes before generating any new responses.
    dump(out/'protocol.json',dict(PROTOCOL,simulation_design=simulation,source_sha256=hashes,
        previous_summary_sha256=sha(PREVIOUS/'summary.json'),
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__)))
    bounds=[(None,None)]*23+[tuple(np.log(PROTOCOL['sigma_bounds']))]
    limits=PROTOCOL['diagnostic_limits']; results={}; planned_r={}
    for condition in PROTOCOL['datasets']:
        name=condition['id']; folder=out/name;folder.mkdir()
        data,latent,truth=simulate(simulation,condition)
        dump(folder/'input.json',dict(classification=PROTOCOL['classification'],qualification_eligible=False,
            scientific_inference_ready=False,condition=condition,data=data,latent_theta=latent,truth=truth,
            forced_extreme_persons=[],integration=PROTOCOL['integration']))
        problem=PCMIntegral(data); records={}; r_cases={}; mp_review={}
        for config in PROTOCOL['grids']:
            theta=np.linspace(-config['bound'],config['bound'],config['q'])
            evaluate=lambda p:fixed_grid(problem,p,theta)
            for sd in PROTOCOL['starts']:
                key=f'{config["id"]}_sd{sd:g}'
                initial=np.r_[np.zeros(23),np.log(sd)]
                preliminary=run_pair(evaluate,lambda p,v:v['nll'],initial,bounds,PROTOCOL['preliminary'])
                dump(folder/(key+'_preliminary.json'),preliminary)
                anchor=np.array(preliminary[-1]['returned']);anchor_value=evaluate(anchor)
                admitted=bool(np.all(np.isfinite(anchor)) and bounds[-1][0]<anchor[-1]<bounds[-1][1]
                    and np.max(abs(anchor_value['gradient']))<=PROTOCOL['admission_gradient'])
                delta=likelihood_difference(problem,anchor,theta)
                refinement=run_pair(evaluate,lambda p,v:delta(p),anchor,bounds,PROTOCOL['refinement']) if admitted else []
                dump(folder/(key+'_runs.json'),dict(initial=initial,anchor=anchor,admitted=admitted,
                    preliminary=preliminary,refinement=refinement))
                p=np.array(refinement[-1]['returned'] if admitted else anchor)
                finite=evaluate(p)
                continuous=[problem.evaluate(p,**setting) for setting in PROTOCOL['integration']]
                error=difference(finite,continuous[-1]);refinement_error=difference(continuous[0],continuous[1])
                local=likelihood_difference(problem,p,theta)
                numerical=central_difference_gradient(local,p,relative_step=1e-6)
                def vg(a):
                    v=evaluate(a)
                    return v['nll'],v['gradient']
                info=information_diagnostics(vg,p)
                restart_displacement=float(np.max(abs(refinement[0]['returned']-p))) if admitted else None
                checks=dict(anchor_admitted=admitted,algorithm_terminated=bool(admitted and all(r['success'] for r in refinement)),
                    requested_gradient=bool(admitted and all(r['requested_gradient_pass'] for r in refinement)),
                    sigma_interior=bool(bounds[-1][0]<p[-1]<bounds[-1][1]),
                    curvature=bool(info.finite and info.positive_definite and info.condition_number<limits['condition']
                        and info.relative_symmetry_residual<limits['relative_symmetry']),
                    standardized_score=bool(info.standardized_score_supnorm<limits['standardized_score']),
                    newton=bool(info.newton_correction_supnorm<limits['newton']),
                    gradient_fd=bool(np.max(abs(numerical-finite['gradient']))<limits['gradient_fd']),
                    restart=bool(admitted and restart_displacement<limits['restart_coordinates']),
                    nonworsening=bool(admitted and all(r['returned_objective']<=delta(r['initial'])+limits['refinement_worsening'] for r in refinement)))
                compact=lambda r:{k:r[k] for k in ['success','status','message','nit','nfev','gradient_supnorm','requested_gradient_pass']}
                record=dict(config=config,start_sd=sd,anchor=anchor,coordinates=p,sigma=float(np.exp(p[-1])),
                    preliminary_results=list(map(compact,preliminary)),refinement_results=list(map(compact,refinement)),
                    preliminary_success=all(r['success'] for r in preliminary),refinement_checks=checks,
                    refinement_checks_pass=all(checks.values()),gradient_supnorm=float(np.max(abs(finite['gradient']))),
                    gradient_fd_error=float(np.max(abs(numerical-finite['gradient']))),curvature=asdict(info),
                    restart_displacement=restart_displacement,anchor_displacement=float(np.max(abs(p-anchor))),
                    anchor_to_final_nll_delta=delta(p),finite=finite,continuous=continuous[-1],
                    continuous_refinement_error=refinement_error,integration_error=error,
                    integration_target_met=bool(error['nll']/len(problem.y)<1e-6 and error['eap']<1e-4 and error['sd']<1e-4),
                    scientific_inference_ready=False)
                records[key]=record;dump(folder/(key+'.json'),record)
                print(name,key,'preliminary',record['preliminary_success'],'refinement',record['refinement_checks_pass'],
                      'gradient',record['gradient_supnorm'],'integration',record['integration_target_met'],flush=True)
                if sd==1:
                    r_cases[config['id']]=dict(coordinates=p,grids=[config],continuous=True)
                    precise={}
                    for digits in (60,90):
                        with mp.workdps(digits):
                            a,b=mp_nll(problem,anchor,theta),mp_nll(problem,p,theta)
                            precise[str(digits)]=dict(anchor=mp.nstr(a,digits),final=mp.nstr(b,digits),delta=mp.nstr(b-a,digits))
                        dump(folder/f'{config["id"]}_mp_{digits}.json',precise[str(digits)])
                    with mp.workdps(90):
                        exact=mp.mpf(precise['90']['delta'])
                        diff=abs(mp.mpf(delta(p))-exact)
                        tolerance=limits['mp_absolute']+limits['mp_relative']*abs(exact)
                        digit_error=max(abs(mp.mpf(precise['60'][k])-mp.mpf(precise['90'][k])) for k in ('anchor','final'))
                    mp_review[config['id']]=dict(error=float(diff),tolerance=float(tolerance),passed=bool(diff<=tolerance),
                        refinement=float(digit_error),refinement_passed=bool(digit_error<limits['mp_refinement']))
            print(name,config['id'],'three starts complete',flush=True)
        cross_start=[]
        for config in PROTOCOL['grids']:
            matched=[r for r in records.values() if r['config']['id']==config['id']]
            cross_start.append(dict(grid=config['id'],coordinate_range=float(np.max(np.ptp([r['coordinates'] for r in matched],axis=0))),
                eap_range=float(np.max(np.ptp([r['finite']['eap'] for r in matched],axis=0))),
                sd_range=float(np.max(np.ptp([r['finite']['sd'] for r in matched],axis=0)))))
        results[name]=dict(records=records,cross_start=cross_start,mp_review=mp_review)
        dump(folder/'r_cases.json',r_cases);planned_r[name]=folder
        dump(folder/'python.json',results[name])
    def run_r(name):
        folder=planned_r[name]
        with (folder/'r_review.log').open('x') as log:
            result=subprocess.run(['Rscript',str(ROOT/'validation/mml_pcm_grid_check.R'),
                str(folder/'input.json'),str(folder),'evaluate'],stdout=log,stderr=subprocess.STDOUT)
        return name,result.returncode
    with ThreadPoolExecutor(max_workers=4) as pool:
        statuses=dict(pool.map(run_r,planned_r))
    independent={}
    for name,folder in planned_r.items():
        cases={}
        if statuses[name]==0:
            for config in PROTOCOL['grids']:
                r=load(folder/(config['id']+'_r.json'))
                py=results[name]['records'][config['id']+'_sd1']
                cases[config['id']]=dict(finite=difference(py['finite'],r['finite'][config['id']]),
                                        continuous=difference(py['continuous'],r['continuous']))
        independent[name]=dict(exit_code=statuses[name],cases=cases)
    rows=[r for d in results.values() for r in d['records'].values()]
    executions=[v for r in rows for v in r['preliminary_results']+r['refinement_results']]
    arithmetic=dict(trajectories_retained=len(rows)==PROTOCOL['planned_trajectories'],
        continuous_refinement=all(max(r['continuous_refinement_error'].values())<limits['continuous_refinement'] for r in rows),
        continuous_error_bounds=all(r['continuous']['numeric_relative_mass_error_sum']<limits['continuous_error_bound'] and
            r['continuous']['tail_relative_mass_bound_sum']<limits['continuous_error_bound'] for r in rows),
        gradient_fd=all(r['gradient_fd_error']<limits['gradient_fd'] for r in rows),
        mp_delta=all(v['passed'] for d in results.values() for v in d['mp_review'].values()),
        mp_refinement=all(v['refinement_passed'] for d in results.values() for v in d['mp_review'].values()),
        independent_r=all(v['exit_code']==0 and len(v['cases'])==2 and all(
            c['nll']<limits['r_nll'] and max(c['eap'],c['sd'])<limits['r_moments']
            for pair in v['cases'].values() for c in pair.values()) for v in independent.values()),
        sources_unchanged=all(sha(p)==h for p,h in hashes.items()),
        historical_summary_unchanged=sha(PREVIOUS/'summary.json')=='b5c3a17b91d193090e1ac810691a2d02d22e1c07dc310de1333c96cb14365eff')
    dump(out/'summary.json',dict(classification=PROTOCOL['classification'],qualification_eligible=False,scientific_inference_ready=False,
        arithmetic_checks=arithmetic,arithmetic_checks_pass=all(arithmetic.values()),optimizer_calls=len(executions),
        all_optimizer_calls_successful=all(r['success'] for r in executions),
        all_refinement_checks_pass=all(r['refinement_checks_pass'] for r in rows),
        all_integration_targets_met=all(r['integration_target_met'] for r in rows),
        results=results,independent_r=independent,protocol_sha256=sha(out/'protocol.json'),
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}))
    print('Arithmetic:',arithmetic,'Refinement:',sum(r['refinement_checks_pass'] for r in rows),'/',len(rows),
          'Raw optimizer success:',sum(r['success'] for r in executions),'/',len(executions),flush=True)
    return 0 if all(arithmetic.values()) else 1


if __name__=='__main__':
    raise SystemExit(main())
