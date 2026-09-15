#!/usr/bin/env python3
"""Audit initial, first, selected and last states of the direct-start experiment."""
import argparse
import json
from pathlib import Path
import re
import subprocess

import numpy as np
from scipy.special import roots_hermitenorm

from mml_pcm_conquest_direct_start import ORDERS, check_outputs
from mml_pcm_conquest_history_audit import coordinates, finite_gh
from mml_pcm_conquest_variance_replay import variance_updates, HISTORY, HISTORY_HASH
from mml_pcm_continuous_refit import PCMIntegral, PREVIOUS, dump, fd_check
from run_pcm_conquest_check import ROOT, sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    out=parser.parse_args().directory.resolve()
    spec=json.loads((out/'input.json').read_text())
    execution=json.loads((out/'execution.json').read_text())
    assert not execution['remaining_unattempted'] and len(execution['runs'])==8
    assert all(r['output_checks_pass'] and r['exit_code']==0 for r in execution['runs'].values())
    for record,field,base in ((spec,'source_sha256',ROOT),(spec,'input_sha256',out),(execution,'artifact_sha256',out)):
        for name,h in record[field].items():
            assert sha(base/name)==h,name
    assert sha(HISTORY/'summary.json')==HISTORY_HASH
    old=json.loads((HISTORY/'summary.json').read_text())
    original=json.loads((PREVIOUS/'input.json').read_text())
    problem=PCMIntegral(original['data'])
    reviewdir=out/'review';reviewdir.mkdir(exist_ok=False)
    r_script=ROOT/'validation/mml_pcm_quadrature_diagnosis.R'
    paths=[Path(__file__).resolve(),r_script]
    hashes=spec['source_sha256']|{str(p.relative_to(ROOT)):sha(p) for p in paths}
    evaluation={k:original[k] for k in ('classification','scientific_inference_ready','qualification_eligible','data','orders','integration')}
    evaluation.update(source_sha256=hashes,execution_sha256=sha(out/'execution.json'),
                      scope='Same-coordinate arithmetic checks; native stopping is recorded separately')
    dump(reviewdir/'evaluation_input.json',evaluation)
    results={};cases={}
    from mml_pcm_conquest_mc import rows
    for q in ORDERS:
        initial=spec['starts'][str(q)]
        start=np.array(initial['coordinates'])
        z,w=roots_hermitenorm(q); w/=np.sqrt(2*np.pi)
        evaluate=lambda par:finite_gh(problem,par,z,w)
        short=out/f'q{q}_iterations1'; full=out/f'q{q}_iterations2000'
        for directory in (short,full): check_outputs(directory)
        short_history=rows(short/'history.csv');history=rows(full/'history.csv')
        assert len(short_history)==1
        assert [int(r['Iteration']) for r in history]==list(range(1,len(history)+1))
        points=np.array([coordinates(r) for r in history])
        duplicate=np.array_equal(coordinates(short_history[0]),points[0]) and short_history[0]['LogLikelihood']==history[0]['LogLikelihood']
        text=(full/'review.txt').read_text()
        selected=int(re.search(r'The number of iterations:\s*(\d+)',text)[1])-1
        assert 0<=selected<len(history)
        parameters=rows(full/'parameters.csv')
        exported=np.array([float(r['Estimate']) for r in parameters])
        assert np.array_equal(np.r_[exported[5:8],exported[:5],exported[8:]],points[selected,:-1])
        assert float(rows(full/'covariance.csv')[0]['Covariance'])==float(history[selected]['wvar 1 1'])
        finite=[evaluate(p) for p in points]
        initial_again=evaluate(start)
        assert abs(initial_again['nll']-initial['finite']['nll'])<1e-8
        cases[f'q{q}_initial']=dict(q=q,coordinates=start,finite=initial_again,
            continuous=problem.evaluate(start,**spec['integration']))
        trace=[]
        for i,(p,v,row) in enumerate(zip(points,finite,history)):
            before=initial_again if i==0 else finite[i-1]
            predicted=variance_updates(before['eap'],before['sd'])
            actual=float(row['wvar 1 1'])
            trace.append(dict(iteration=i+1,sigma=float(np.exp(p[-1])),reported_nll=float(row['LogLikelihood'])/2,
                finite_nll=v['nll'],nll_excess_to_start=v['nll']-initial_again['nll'],
                previous_point_nll_error=float(row['LogLikelihood'])/2-before['nll'],
                centered_variance_error=actual-predicted['centered_mean'],
                fixed_zero_mean_variance_error=actual-predicted['fixed_zero_mean'],
                log_sd_gradient=float(v['gradient'][-1]),gradient_supnorm=float(np.max(abs(v['gradient'])))))
        dump(reviewdir/f'q{q}_trace.json',trace)
        checkpoints={}
        for role,index in dict(first=0,selected=selected,last=len(history)-1).items():
            key=f'q{q}_iteration{index+1}'
            if key not in cases:
                cases[key]=dict(q=q,coordinates=points[index],finite=finite[index],
                    continuous=problem.evaluate(points[index],**spec['integration']))
            checkpoints[role]=dict(case=key,iteration=index+1,
                coordinate_movement=float(np.max(abs(points[index]-start))),
                structural_movement=float(np.max(abs(points[index,:-1]-start[:-1]))),
                sigma=float(np.exp(points[index,-1])),sigma_movement=float(np.exp(points[index,-1])-np.exp(start[-1])),
                finite_nll_excess_to_start=finite[index]['nll']-initial_again['nll'],
                continuous_nll_change=cases[key]['continuous']['nll']-cases[f'q{q}_initial']['continuous']['nll'],
                gradient_supnorm=float(np.max(abs(finite[index]['gradient']))))
        fd=fd_check(evaluate,points[0]); assert max(c['max_difference'] for c in fd)<1e-5
        baseline=np.array(old['results'][f'cq_gh{q}']['endpoints']['last']['coordinates'])
        results[str(q)]=dict(initial=initial,short_first_matches_full_first=bool(duplicate),
            native_parameter_count=int(re.search(r'Total number of estimated parameters:\s*(-?\d+)',text)[1]),
            selected_iteration=selected+1,executed_rows=len(history),reached_iteration_cap=len(history)>=2000,
            termination=re.search(r'Iterations terminated[^\n]+',text)[0],checkpoints=checkpoints,
            first_update=trace[0],last_update=trace[-1],first_gradient_fd=fd,
            maximum_previous_nll_error=max(abs(t['previous_point_nll_error']) for t in trace),
            maximum_centered_variance_error=max(abs(t['centered_variance_error']) for t in trace),
            last_coordinate_distance_to_previous_native=float(np.max(abs(points[-1]-baseline))),
            last_finite_nll_difference_to_previous_native=finite[-1]['nll']-evaluate(baseline)['nll'])
        print(q,'first',checkpoints['first'],'iterations',len(history),flush=True)
    dump(reviewdir/'python.json',cases)
    subprocess.run(['Rscript',str(r_script),str(reviewdir/'evaluation_input.json'),str(reviewdir/'python.json'),str(reviewdir)],check=True)
    differences={}
    for name,case in cases.items():
        independent=json.loads((reviewdir/f'{name}_r.json').read_text())
        differences[name]={kind:{k:float(np.max(abs(np.asarray(case[kind][k])-reference[k])))
            for k in ('nll','eap','sd')} for kind,reference in
            (('finite',independent['finite'][str(case['q'])]),('continuous',independent['continuous'][1]))}
    arithmetic=all(d<1e-8 for c in differences.values() for kind in c.values() for d in kind.values())
    assert hashes=={p:sha(ROOT/p) for p in hashes}
    dump(reviewdir/'summary.json',dict(classification='OBSERVED_DEVELOPMENT_ONLY',scientific_inference_ready=False,
        qualification_eligible=False,independent_arithmetic_checks_pass=arithmetic,
        all_first_step_controls_match=all(r['short_first_matches_full_first'] for r in results.values()),
        all_native_parameter_counts_24=all(r['native_parameter_count']==24 for r in results.values()),
        results=results,cases=cases,independent_r_differences=differences,source_sha256=hashes,
        artifact_sha256={p.name:sha(p) for p in reviewdir.iterdir() if p.is_file()},
        limitations=['One observed PCM fixture and one direct stationary start per quadrature order',
            'Native exports are rounded; no bitwise equality to unrounded imported coordinates',
            'The one-step cap is intentional; complete output is not convergence qualification',
            'Conditional scores use deterministic integration here; native Monte Carlo scores and SEs are not qualified',
            'Finite-GH stationary starts are not necessarily continuous-integral optima']))
    print('Independent arithmetic:',arithmetic,'Saved',reviewdir)
    return 0 if arithmetic else 1


if __name__=='__main__':
    raise SystemExit(main())
