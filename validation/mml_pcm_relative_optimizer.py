#!/usr/bin/env python3
"""Frozen fixed-anchor optimizer comparison; preserve every termination outcome."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import platform

import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import minimize

from mml_pcm_continuous_refit import ROOT, dump, sha
from mml_pcm_independent_conditions import PCMIntegral, difference
from mml_pcm_external_review import fixed_grid
from mml_pcm_likelihood_difference import likelihood_difference
from mml_pcm_optimizer_resolution import mp_nll
from mfrm_app.mml_stationarity import (
    central_difference_gradient, information_diagnostics, projected_gradient,
)

GRID = ROOT/'validation/generated/mml_pcm_grid_range_density_20260914'
RATIO = ROOT/'validation/generated/mml_pcm_likelihood_difference_20260914'
INPUT_HASHES = {
    str(GRID/'summary.json'): '53e68b8425a6f40347c2b57a0a35a4e229e51a4ba5757970c722e08a0736656d',
    str(RATIO/'summary.json'): 'b6ca553543d2752ec4aaa9e4635fb9f368ada7d5b29cdffeb8d3d525b1bb5885',
}
PROTOCOL = dict(
    id='pcm_relative_optimizer_20260914_v1', classification='OBSERVED_DEVELOPMENT_ONLY',
    qualification_eligible=False, scientific_inference_ready=False,
    question='Does fixed-anchor likelihood-ratio arithmetic improve gradient attainment, and how do ftol and anchor distance affect it?',
    cells={'ordinary_1':'b20_q801','ordinary_2':'b20_q801',
           'wide_1':'b12_q241','wide_2':'b20_q801'},
    starts=['zero','retained'], arms=['raw','relative_zero','relative_retained'],
    anchor_policy='Zero or the selected endpoint in the frozen cell; unchanged for the entire primary/restart pair, independently of the chosen start',
    ftols={'zero':0.0,'historical':1e-15},
    optimizer=dict(gtol=1e-8,maxiter=250,maxls=50,maxcor=10,maxfun=15000),
    sigma_bounds=[0.05,10.0], passes=2, planned_calls=96,
    stopping='Same requested gtol in all arms. ftol=0 removes positive relative-reduction tolerance, but still allows CONV_F at zero/nonpositive computed reduction. Historical ftol is a scale-sensitivity control, not an equivalent stopping rule across offsets.',
    gradient='All arms use the unchanged analytical gradient of the unnormalized fixed-grid total NLL',
    review='Retain every query, accepted iterate, initial/returned/last queried coordinates, raw optimizer fields and fresh endpoint values. Never pick a preferred start after seeing results.',
    endpoint_checks='Both returns: raw/projected gradients. Final endpoint of every pair: local curvature, local-ratio and actual-objective FD at 1e-6, continuous-integral evaluation, score errors, restart change and cross-start distances.',
    mp_selection='For each cell, the retained anchor and the two final relative_retained ftol=0 endpoints; 60/90 digits, no outcome-based selection',
    limits=dict(scalar_reconstruction=1e-8,local_gradient_fd=1e-7,
                mp_delta_absolute=1e-18,mp_delta_relative=1e-11,mp_refinement=1e-45),
    integration_targets=dict(nll_per_person=1e-6,eap=1e-4,sd=1e-4),
    failure_policy='Keep all failures and caps; no extra retries or post-result settings changes. Execution/arithmetic flags are separate from optimizer success, requested-gradient attainment and inference.',
    scope='Four observed complete PCM cells; no new data, native TAM/ConQuest runs, GH adaptation, normalized-prior objective, missing data, weights or inference qualification',
)


def run_pair(evaluate, objective, start, bounds, options):
    point=np.array(start)
    runs=[]
    for _ in range(PROTOCOL['passes']):
        initial=point.copy()
        queries,accepted=[],[]
        def value_gradient(p):
            v=evaluate(p)
            f=objective(p,v)
            queries.append(dict(coordinates=p.copy(),objective=f,raw_nll=v['nll'],gradient=v['gradient']))
            return f,v['gradient']
        opt=minimize(value_gradient,initial,jac=True,method='L-BFGS-B',bounds=bounds,
                     options=options,callback=lambda p:accepted.append(p.copy()))
        point=opt.x.copy()
        v=evaluate(point)
        grad=projected_gradient(point,v['gradient'],bounds)
        runs.append(dict(initial=initial,returned=point,success=bool(opt.success),status=int(opt.status),
            message=str(opt.message),nit=int(opt.nit),nfev=int(opt.nfev),njev=int(opt.njev),
            reported_objective=float(opt.fun),returned_objective=objective(point,v),returned_value=v,
            gradient_supnorm=float(np.max(abs(v['gradient']))),projected_gradient_supnorm=float(np.max(abs(grad))),
            requested_gradient_pass=bool(np.max(abs(grad))<=options['gtol']),
            last_query=queries[-1],last_accepted=accepted[-1] if accepted else initial,
            queries=queries,accepted=accepted))
    return runs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    out=parser.parse_args().output_dir.resolve()
    load=lambda path:json.loads(path.read_text())
    for path,h in INPUT_HASHES.items():
        assert sha(path)==h,path
        for name,expected in load(Path(path))['artifact_sha256'].items():
            assert sha(Path(path).parent/name)==expected,name
    for path,h in load(RATIO/'input.json')['source_sha256'].items():
        assert sha(path)==h,path
    paths=[Path(__file__).resolve(),ROOT/'mfrm_app/mml_stationarity.py']+[
        Path(__file__).with_name('mml_pcm_'+name+'.py') for name in
        ['continuous_refit','independent_conditions','external_review','likelihood_difference','optimizer_resolution']]
    hashes={str(p):sha(p) for p in paths}
    out.mkdir(parents=True,exist_ok=False)
    # Freeze the whole comparison before any optimizer or new numerical evaluation.
    dump(out/'input.json',dict(PROTOCOL,input_sha256=INPUT_HASHES,source_sha256=hashes,
        environment=dict(python=platform.python_version(),numpy=np.__version__,
                         scipy=scipy.__version__,mpmath=mp.__version__)))
    bounds=[(None,None)]*23+[tuple(np.log(PROTOCOL['sigma_bounds']))]
    results,mp_reviews={},{}
    for name,cell_id in PROTOCOL['cells'].items():
        folder=out/name; folder.mkdir()
        original=load(GRID/name/'input.json')
        cell=load(GRID/name/(cell_id+'.json'))
        problem=PCMIntegral(original['data'])
        theta=np.linspace(-cell['config']['bound'],cell['config']['bound'],cell['config']['q'])
        retained=np.array(cell['fit']['coordinates'])
        starts=dict(zero=np.zeros(24),retained=retained)
        evaluate=lambda p:fixed_grid(problem,p,theta)
        zero_delta=likelihood_difference(problem,starts['zero'],theta)
        retained_delta=likelihood_difference(problem,retained,theta)
        objectives=dict(raw=lambda p,v:v['nll'],relative_zero=lambda p,v:zero_delta(p),
                        relative_retained=lambda p,v:retained_delta(p))
        offsets=dict(raw=0.,relative_zero=evaluate(starts['zero'])['nll'],
                     relative_retained=evaluate(retained)['nll'])
        records={}
        for tolerance,ftol in PROTOCOL['ftols'].items():
            for arm,objective in objectives.items():
                for start_name,start in starts.items():
                    key=f'{tolerance}_{arm}_{start_name}'
                    options=dict(PROTOCOL['optimizer'],ftol=ftol)
                    runs=run_pair(evaluate,objective,start,bounds,options)
                    dump(folder/(key+'_runs.json'),dict(options=options,arm=arm,start=start_name,
                        anchor=starts['zero'] if arm=='relative_zero' else retained if arm=='relative_retained' else None,
                        runs=runs))
                    p=runs[-1]['returned']; finite=runs[-1]['returned_value']
                    continuous=problem.evaluate(p,**original['integration'][1])
                    error=difference(finite,continuous)
                    local=likelihood_difference(problem,p,theta)
                    local_fd=central_difference_gradient(local,p,relative_step=1e-6)
                    actual_fd=central_difference_gradient(lambda a:objective(a,evaluate(a)),p,relative_step=1e-6)
                    def vg(a):
                        v=evaluate(a)
                        return v['nll'],v['gradient']
                    reconstruction=max(abs(r['returned_objective']+offsets[arm]-r['returned_value']['nll']) for r in runs)
                    record=dict(tolerance=tolerance,arm=arm,start=start_name,coordinates=p,
                        raw_nll=finite['nll'],continuous=continuous,integration_error=error,
                        integration_target_met=bool(error['nll']/len(problem.y)<1e-6 and error['eap']<1e-4 and error['sd']<1e-4),
                        sigma=float(np.exp(p[-1])),sigma_interior=bool(bounds[-1][0]<p[-1]<bounds[-1][1]),
                        gradient_supnorm=runs[-1]['gradient_supnorm'],requested_gradient_pass=runs[-1]['requested_gradient_pass'],
                        both_optimizer_calls_successful=all(r['success'] for r in runs),
                        both_requested_gradient_checks_pass=all(r['requested_gradient_pass'] for r in runs),
                        curvature=asdict(information_diagnostics(vg,p)),
                        local_gradient_fd_max=float(np.max(abs(local_fd-finite['gradient']))),
                        actual_objective_fd_max=float(np.max(abs(actual_fd-finite['gradient']))),
                        scalar_reconstruction_error=reconstruction,
                        last_reported_objective_error=abs(runs[-1]['reported_objective']-runs[-1]['returned_objective']),
                        restart_displacement=float(np.max(abs(runs[0]['returned']-p))),
                        restart_improvement=likelihood_difference(problem,p,theta)(runs[0]['returned']),
                        retained_to_final_nll_delta=retained_delta(p),
                        coordinate_distance_to_retained=float(np.max(abs(p-retained))),
                        optimizer_results=[{k:r[k] for k in ['success','status','message','nit','nfev','gradient_supnorm','requested_gradient_pass']} for r in runs],
                        scientific_inference_ready=False)
                    records[key]=record
                    dump(folder/(key+'.json'),record)
                    print(name,key,'g',record['gradient_supnorm'],'success',record['both_optimizer_calls_successful'],
                          'gtol',record['both_requested_gradient_checks_pass'],flush=True)
        distances=[]
        for tolerance in PROTOCOL['ftols']:
            for arm in objectives:
                a,b=[records[f'{tolerance}_{arm}_{start}']['coordinates'] for start in starts]
                distances.append(dict(tolerance=tolerance,arm=arm,
                    cross_start_coordinates=float(np.max(abs(a-b))),
                    cross_start_nll_delta=likelihood_difference(problem,a,theta)(b)))
        results[name]=dict(records=records,cross_start=distances)
        dump(folder/'comparison.json',results[name])
        # Exactly the predeclared endpoints; unsuccessful runs are not excluded.
        selected=dict(anchor=retained,**{start:records[f'zero_relative_retained_{start}']['coordinates'] for start in starts})
        precise={}
        for digits in (60,90):
            with mp.workdps(digits):
                values={k:mp_nll(problem,p,theta) for k,p in selected.items()}
                precise[str(digits)]=dict(nll={k:mp.nstr(v,digits) for k,v in values.items()},
                    delta={k:mp.nstr(v-values['anchor'],digits) for k,v in values.items()})
            dump(folder/f'mp_{digits}.json',precise[str(digits)])
            print(name,'MP',digits,'complete',flush=True)
        checks={}
        with mp.workdps(90):
            for start in starts:
                ref=mp.mpf(precise['90']['delta'][start])
                observed=records[f'zero_relative_retained_{start}']['retained_to_final_nll_delta']
                error=abs(mp.mpf(observed)-ref)
                limit=PROTOCOL['limits']['mp_delta_absolute']+PROTOCOL['limits']['mp_delta_relative']*abs(ref)
                refinement=max(abs(mp.mpf(precise['60']['nll'][k])-mp.mpf(precise['90']['nll'][k])) for k in ('anchor',start))
                checks[start]=dict(delta_error=float(error),delta_tolerance=float(limit),delta_pass=bool(error<=limit),
                    refinement=float(refinement),refinement_pass=bool(refinement<PROTOCOL['limits']['mp_refinement']))
        mp_reviews[name]=checks
        dump(folder/'mp_review.json',checks)
    rows=[r for c in results.values() for r in c['records'].values()]
    all_runs=[r for row in rows for r in row['optimizer_results']]
    groups={}
    for tolerance in PROTOCOL['ftols']:
        for arm in PROTOCOL['arms']:
            subset=[r for r in rows if r['tolerance']==tolerance and r['arm']==arm]
            runs=[r for row in subset for r in row['optimizer_results']]
            groups[f'{tolerance}_{arm}']=dict(calls=len(runs),successful_calls=sum(r['success'] for r in runs),
                requested_gradient_calls=sum(r['requested_gradient_pass'] for r in runs),
                final_requested_gradient_pairs=sum(r['requested_gradient_pass'] for r in subset),
                final_gradient_max=max(r['gradient_supnorm'] for r in subset),
                final_newton_correction_max=max(r['curvature']['newton_correction_supnorm'] for r in subset))
    checks=dict(planned_calls_retained=len(all_runs)==PROTOCOL['planned_calls'],
        finite_endpoints=all(np.all(np.isfinite(r['coordinates'])) and np.isfinite(r['raw_nll']) for r in rows),
        scalar_reconstruction=all(r['scalar_reconstruction_error']<PROTOCOL['limits']['scalar_reconstruction'] for r in rows),
        local_gradient_fd=all(r['local_gradient_fd_max']<PROTOCOL['limits']['local_gradient_fd'] for r in rows),
        mp_delta=all(v['delta_pass'] for c in mp_reviews.values() for v in c.values()),
        mp_refinement=all(v['refinement_pass'] for c in mp_reviews.values() for v in c.values()),
        inputs_unchanged=all(sha(p)==h for p,h in INPUT_HASHES.items()),
        sources_unchanged=all(sha(p)==h for p,h in hashes.items()),
        historical_failure_retained=load(GRID/'summary.json')['implementation_checks_pass'] is False)
    dump(out/'summary.json',dict(classification=PROTOCOL['classification'],scientific_inference_ready=False,
        qualification_eligible=False,execution_arithmetic_checks=checks,execution_arithmetic_checks_pass=all(checks.values()),
        all_optimizer_calls_successful=all(r['success'] for r in all_runs),
        all_requested_gradient_checks_pass=all(r['requested_gradient_pass'] for r in all_runs),
        historical_implementation_checks_pass=False,groups=groups,results=results,mp_reviews=mp_reviews,
        protocol_sha256=sha(out/'input.json'),
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}))
    print('Execution/arithmetic:',checks,'Optimizer success:',sum(r['success'] for r in all_runs),
          'Requested gradient:',sum(r['requested_gradient_pass'] for r in all_runs),'/',len(all_runs),flush=True)
    return 0 if all(checks.values()) else 1


if __name__=='__main__':
    raise SystemExit(main())
