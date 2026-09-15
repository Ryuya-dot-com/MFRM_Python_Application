#!/usr/bin/env python3
"""Frozen moving-GH two-stage PCM experiment; preserve failures and local solutions."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import platform
import subprocess

import mpmath as mp
import numpy as np
import scipy
from scipy.special import roots_hermitenorm

from mml_pcm_continuous_refit import ROOT, dump, sha
from mml_pcm_independent_conditions import PCMIntegral, difference
from mml_pcm_conquest_history_audit import finite_gh
from mml_pcm_gh_likelihood_difference import gh_likelihood_difference, mp_gh
from mml_pcm_relative_optimizer import run_pair
from mml_pcm_two_stage import PROTOCOL as FIXED
from mfrm_app.mml_stationarity import central_difference_gradient, information_diagnostics

DATA = ROOT/'validation/generated/mml_pcm_two_stage_20260914'
ARITHMETIC = ROOT/'validation/generated/mml_pcm_gh_likelihood_difference_20260914'
INPUT_HASHES = {str(DATA/'summary.json'):'c3558407a71b08d2ea555d860904f3750f0b23b41e1e537e7465bb5cd8018c15',
               str(ARITHMETIC/'summary.json'):'359596321455a45c8b514be12cb5b3634607014283a722d7abf5a5b2fa001749'}
PROTOCOL = dict({k:FIXED[k] for k in ('datasets','starts','structural_start','preliminary','refinement',
    'passes_per_stage','sigma_bounds','anchor','admission_gradient','refusal','reporting',
    'diagnostic_limits','integration','integration_targets')},
    id='pcm_gh_two_stage_20260914_v1',classification='DEVELOPMENT_ONLY',
    qualification_eligible=False,scientific_inference_ready=False,
    question='Can self-anchored moving-GH refinement attain gtol, and do distinct starting SDs or quadrature orders still give different solutions and integration errors?',
    data_status='All four retained, previously observed 80-person datasets; no new generation, response editing, warm starts or endpoint selection',
    orders=[31,61,121,181],planned_trajectories=48,maximum_optimizer_calls=192,
    measure='Standard-normal GH weights fixed and positive; theta=exp(log-SD)*z. Exact moving-node ratio and unchanged analytic finite-GH gradient.',
    independent_selection='Starting SD 1, every dataset and fitted order, including failures: 60/90-digit NLL at preliminary anchor and final return. Independent R finite GH at all four rules for each endpoint, plus its continuous moments.',
    cross_order='At each SD-1 endpoint, evaluate all four rules without refitting; compare to continuous integration at those same coordinates',
    external_reference='Retained b20_q801_sd1 fixed-grid endpoint for each dataset, used only after fitting to compare coordinates and continuous objective; never a starting value or claimed exact MLE',
    failure_policy='Freeze inputs, source and settings before fitting. No additional retries, changed thresholds, best-start selection, or overwritten raw failure flags.',
    scope='Complete unit-weight mean-zero PCM numerical development only; no app/API changes, native TAM/ConQuest fitting, global-optimum or SE/CI/coverage qualification',
)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    out=parser.parse_args().output_dir.resolve()
    load=lambda p:json.loads(p.read_text())
    for file,h in INPUT_HASHES.items():
        file=Path(file);assert sha(file)==h,file
        for name,value in load(file)['artifact_sha256'].items():
            assert sha(file.parent/name)==value,name
        for name,value in load(file.parent/'protocol.json')['source_sha256'].items():
            assert sha(name)==value,name
    paths=[Path(__file__).resolve(),ROOT/'mfrm_app/mml_stationarity.py',ROOT/'mfrm_app/mml_engine_v2.py']+[
        Path(__file__).with_name('mml_pcm_'+name) for name in ('two_stage.py','relative_optimizer.py',
        'continuous_refit.py','independent_conditions.py','conquest_history_audit.py',
        'gh_likelihood_difference.py','likelihood_difference.py','independent_conditions.R')]
    hashes={str(p):sha(p) for p in paths}
    out.mkdir(parents=True,exist_ok=False)
    dump(out/'protocol.json',dict(PROTOCOL,input_sha256=INPUT_HASHES,source_sha256=hashes,
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__)))
    rules={}
    for q in PROTOCOL['orders']:
        z,w=roots_hermitenorm(q);w/=np.sqrt(2*np.pi)
        assert np.all(w>0) and abs(w.sum()-1)<1e-12
        rules[q]=(z,w)
    dump(out/'rules.json',{str(q):dict(nodes=z,weights=w) for q,(z,w) in rules.items()})
    bounds=[(None,None)]*23+[tuple(np.log(PROTOCOL['sigma_bounds']))]
    limits=PROTOCOL['diagnostic_limits'];targets=PROTOCOL['integration_targets']
    results={};r_jobs=[]
    for condition in PROTOCOL['datasets']:
        name=condition['id'];folder=out/name;folder.mkdir()
        inp=load(DATA/name/'input.json');inp['orders']=PROTOCOL['orders']
        dump(folder/'input.json',inp)
        problem=PCMIntegral(inp['data']);records={};mp_review={}
        for q,(z,w) in rules.items():
            evaluate=lambda p:finite_gh(problem,p,z,w)
            for sd in PROTOCOL['starts']:
                key=f'q{q}_sd{sd:g}'
                initial=np.r_[np.zeros(23),np.log(sd)]
                preliminary=run_pair(evaluate,lambda p,v:v['nll'],initial,bounds,PROTOCOL['preliminary'])
                dump(folder/(key+'_preliminary.json'),preliminary)
                anchor=np.array(preliminary[-1]['returned']);av=evaluate(anchor)
                admitted=bool(np.all(np.isfinite(anchor)) and bounds[-1][0]<anchor[-1]<bounds[-1][1]
                    and np.max(abs(av['gradient']))<=PROTOCOL['admission_gradient'])
                delta=gh_likelihood_difference(problem,anchor,z,w)
                refinement=run_pair(evaluate,lambda p,v:delta(p),anchor,bounds,PROTOCOL['refinement']) if admitted else []
                dump(folder/(key+'_runs.json'),dict(initial=initial,anchor=anchor,admitted=admitted,
                    preliminary=preliminary,refinement=refinement))
                p=np.array(refinement[-1]['returned'] if admitted else anchor);finite=evaluate(p)
                continuous=[problem.evaluate(p,**setting) for setting in PROTOCOL['integration']]
                error=difference(finite,continuous[-1])
                local=gh_likelihood_difference(problem,p,z,w)
                step=1e-6*max(1,abs(p[-1]))
                numerical=(central_difference_gradient(local,p,relative_step=1e-6)
                    if bounds[-1][0]<=p[-1]-step and p[-1]+step<=bounds[-1][1] else None)
                fd_error=float(np.max(abs(numerical-finite['gradient']))) if numerical is not None else None
                def vg(point):
                    v=evaluate(point)
                    return v['nll'],v['gradient']
                info=information_diagnostics(vg,p)
                restart=float(np.max(abs(refinement[0]['returned']-p))) if admitted else None
                reconstruction=max((abs(v['objective']+av['nll']-v['raw_nll']) for r in refinement for v in r['queries']),default=0.)
                checks=dict(anchor_admitted=admitted,algorithm_terminated=bool(admitted and all(r['success'] for r in refinement)),
                    requested_gradient=bool(admitted and all(r['requested_gradient_pass'] for r in refinement)),
                    sigma_interior=bool(bounds[-1][0]<p[-1]<bounds[-1][1]),
                    curvature=bool(info.finite and info.positive_definite and info.condition_number<limits['condition']
                        and info.relative_symmetry_residual<limits['relative_symmetry']),
                    standardized_score=bool(info.standardized_score_supnorm<limits['standardized_score']),
                    newton=bool(info.newton_correction_supnorm<limits['newton']),
                    gradient_fd=bool(fd_error is not None and fd_error<limits['gradient_fd']),
                    restart=bool(admitted and restart<limits['restart_coordinates']),
                    nonworsening=bool(admitted and all(r['returned_objective']<=delta(r['initial'])+limits['refinement_worsening'] for r in refinement)))
                compact=lambda r:{k:r[k] for k in ('success','status','message','nit','nfev','gradient_supnorm','requested_gradient_pass')}
                record=dict(q=q,start_sd=sd,anchor=anchor,coordinates=p,sigma=float(np.exp(p[-1])),
                    preliminary_results=list(map(compact,preliminary)),refinement_results=list(map(compact,refinement)),
                    preliminary_success=all(r['success'] for r in preliminary),refinement_checks=checks,
                    refinement_checks_pass=all(checks.values()),finite=finite,continuous=continuous[-1],
                    gradient_fd_error=fd_error,curvature=asdict(info),
                    restart_displacement=restart,anchor_displacement=float(np.max(abs(p-anchor))),
                    anchor_to_final_nll_delta=delta(p),scalar_reconstruction_error=reconstruction,
                    integration_error=error,continuous_refinement_error=difference(continuous[0],continuous[1]),
                    integration_target_met=bool(error['nll']/len(problem.y)<targets['nll_per_person'] and
                        error['eap']<targets['eap'] and error['sd']<targets['sd']),scientific_inference_ready=False)
                records[key]=record;dump(folder/(key+'.json'),record)
                print(name,key,'refinement',record['refinement_checks_pass'],'gradient',float(np.max(abs(finite['gradient']))),
                      'sigma',record['sigma'],'integration',record['integration_target_met'],flush=True)
                if sd==1:
                    precise={}
                    for digits in (60,90):
                        with mp.workdps(digits):
                            a,b=mp_gh(problem,anchor,z,w)['nll'],mp_gh(problem,p,z,w)['nll']
                            precise[str(digits)]=dict(anchor=mp.nstr(a,digits),final=mp.nstr(b,digits),delta=mp.nstr(b-a,digits))
                        dump(folder/f'q{q}_mp_{digits}.json',precise[str(digits)])
                    with mp.workdps(90):
                        exact=mp.mpf(precise['90']['delta']);diff=abs(mp.mpf(delta(p))-exact)
                        tolerance=limits['mp_absolute']+limits['mp_relative']*abs(exact)
                        digit_error=max(abs(mp.mpf(precise['60'][k])-mp.mpf(precise['90'][k])) for k in ('anchor','final','delta'))
                    mp_review[str(q)]=dict(error=float(diff),tolerance=float(tolerance),passed=bool(diff<=tolerance),
                        refinement=float(digit_error),refinement_passed=bool(digit_error<limits['mp_refinement']))
        r_cases={f'q{q}':dict(coordinates=records[f'q{q}_sd1']['coordinates']) for q in rules}
        dump(folder/'r_cases.json',r_cases)
        log=(folder/'r.log').open('x')
        process=subprocess.Popen(['Rscript',str(paths[-1]),str(folder/'input.json'),str(folder/'r_cases.json'),str(folder)],stdout=log,stderr=subprocess.STDOUT)
        r_jobs.append((name,process,log))
        cross_start=[];cross_order={}
        reference=load(DATA/name/'b20_q801_sd1.json')
        for q,(z,w) in rules.items():
            matched=[r for r in records.values() if r['q']==q]
            local=gh_likelihood_difference(problem,matched[0]['coordinates'],z,w)
            cross_start.append(dict(q=q,coordinate_range=float(np.max(np.ptp([r['coordinates'] for r in matched],axis=0))),
                sigma_range=float(np.ptp([r['sigma'] for r in matched])),
                nll_range=float(np.ptp([local(r['coordinates']) for r in matched])),
                eap_range=float(np.max(np.ptp([r['finite']['eap'] for r in matched],axis=0))),
                sd_range=float(np.max(np.ptp([r['finite']['sd'] for r in matched],axis=0)))))
            point=records[f'q{q}_sd1']
            cross_order[str(q)]=dict(coordinates=point['coordinates'],
                finite={str(other):finite_gh(problem,point['coordinates'],oz,ow) for other,(oz,ow) in rules.items()},
                coordinate_difference_to_reference=float(np.max(abs(point['coordinates']-reference['coordinates']))),
                sigma_difference_to_reference=point['sigma']-reference['sigma'],
                continuous_nll_excess_to_reference=point['continuous']['nll']-reference['continuous']['nll'],
                continuous_score_difference_to_reference=difference(point['continuous'],reference['continuous']))
        results[name]=dict(records=records,cross_start=cross_start,cross_order=cross_order,mp_review=mp_review)
        dump(folder/'python.json',results[name])
    independent={}
    for name,process,log in r_jobs:
        status=process.wait();log.close();cases={}
        dump(out/name/'r_execution.json',dict(exit_code=status))
        if status==0:
            for q in rules:
                r=load(out/name/f'q{q}_r.json');py=results[name]
                cases[str(q)]=dict(continuous=difference(py['records'][f'q{q}_sd1']['continuous'],r['continuous'][-1]),
                    finite={str(other):difference(py['cross_order'][str(q)]['finite'][str(other)],r['finite'][str(other)]) for other in rules})
        independent[name]=dict(exit_code=status,cases=cases)
    rows=[r for d in results.values() for r in d['records'].values()]
    executions=[r for v in rows for r in v['preliminary_results']+v['refinement_results']]
    r_errors=[e for d in independent.values() for r in d['cases'].values() for e in [r['continuous'],*r['finite'].values()]]
    arithmetic=dict(trajectories_retained=len(rows)==PROTOCOL['planned_trajectories'],
        continuous_refinement=all(max(r['continuous_refinement_error'].values())<limits['continuous_refinement'] for r in rows),
        continuous_error_bounds=all(r['continuous']['numeric_relative_mass_error_sum']<limits['continuous_error_bound'] and
            r['continuous']['tail_relative_mass_bound_sum']<limits['continuous_error_bound'] for r in rows),
        gradient_fd=all(r['gradient_fd_error'] is not None and r['gradient_fd_error']<limits['gradient_fd'] for r in rows),
        scalar_reconstruction=all(r['scalar_reconstruction_error']<limits['scalar_reconstruction'] for r in rows),
        mp_delta=all(v['passed'] for d in results.values() for v in d['mp_review'].values()),
        mp_refinement=all(v['refinement_passed'] for d in results.values() for v in d['mp_review'].values()),
        independent_r=bool(all(d['exit_code']==0 and len(d['cases'])==4 for d in independent.values()) and all(
            e['nll']<limits['r_nll'] and max(e['eap'],e['sd'])<limits['r_moments'] for e in r_errors)),
        sources_unchanged=all(sha(p)==h for p,h in hashes.items()),
        previous_summaries_unchanged=all(sha(p)==h for p,h in INPUT_HASHES.items()))
    dump(out/'summary.json',dict(classification=PROTOCOL['classification'],qualification_eligible=False,scientific_inference_ready=False,
        arithmetic_checks=arithmetic,arithmetic_checks_pass=all(arithmetic.values()),optimizer_calls=len(executions),
        all_optimizer_calls_successful=all(r['success'] for r in executions),
        all_refinement_checks_pass=all(r['refinement_checks_pass'] for r in rows),
        all_integration_targets_met=all(r['integration_target_met'] for r in rows),results=results,
        independent_r=independent,protocol_sha256=sha(out/'protocol.json'),
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}))
    print('Arithmetic:',arithmetic,'Refinement:',sum(r['refinement_checks_pass'] for r in rows),'/',len(rows),
          'Raw optimizer success:',sum(r['success'] for r in executions),'/',len(executions),flush=True)
    return 0 if all(arithmetic.values()) else 1


if __name__=='__main__':
    raise SystemExit(main())
