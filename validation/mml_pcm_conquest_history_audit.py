#!/usr/bin/env python3
"""Replay every retained ConQuest GH iteration without refitting ConQuest."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp, roots_hermitenorm

from mml_pcm_conquest_mc import SOURCE, SOURCE_HASH, rows
from mml_pcm_continuous_refit import PCMIntegral, PREVIOUS, SUMMARY_HASH, dump, fd_check
from run_pcm_conquest_check import ROOT, sha


def finite_gh(problem, par, z, w):
    par=np.asarray(par)
    assert par.shape==(24,) and np.all(np.isfinite(par)) and np.all(w>0)
    theta=np.exp(par[-1])*z
    logits=(problem.design@par[:-1])[:,:,None]+problem.k[None,:,None]*theta
    norm=logsumexp(logits,axis=1)
    prob=np.exp(logits-norm[:,None,:])
    joint=(theta[None,:]*problem.total[:,None]+(problem.observed@par[:-1])[:,None]
           -norm.sum(axis=0)+np.log(w))
    mass=logsumexp(joint,axis=1)
    post=np.exp(joint-mass[:,None])
    structural=-(problem.observed-post@np.einsum('ikq,ikd->qd',prob,problem.design)).sum(axis=0)
    # Differentiate moving GH nodes: d(theta)/d(log sigma) = theta.
    log_sd=-np.sum(post*theta*(problem.total[:,None]-(prob*problem.k[None,:,None]).sum(axis=(0,1))))
    eap=post@theta
    return dict(nll=float(-mass.sum()), gradient=np.r_[structural,log_sd], eap=eap,
                sd=np.sqrt(np.sum(post*(theta-eap[:,None])**2,axis=1)))


def coordinates(row):
    xsi=np.array([float(row[f'xsi {i}']) for i in range(1,24)])
    assert float(row['Dim 1 Var 1'])==0 and float(row['wvar 1 1'])>0
    return np.r_[xsi[5:8],xsi[:5],xsi[8:],.5*np.log(float(row['wvar 1 1']))]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    out=parser.parse_args().output_dir.resolve()
    assert sha(SOURCE/'review_summary.json')==SOURCE_HASH
    assert sha(PREVIOUS/'summary.json')==SUMMARY_HASH
    source=json.loads((SOURCE/'review_summary.json').read_text())
    previous=json.loads((PREVIOUS/'summary.json').read_text())
    for base,record in ((SOURCE,source),(PREVIOUS,previous)):
        for name,h in record['artifact_sha256'].items():
            assert sha(base/name)==h,name
    original=json.loads((PREVIOUS/'input.json').read_text())
    app_fits=json.loads((PREVIOUS/'python.json').read_text())
    out.mkdir(parents=True,exist_ok=False)
    hashes={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__).resolve(),
        ROOT/'validation/mml_pcm_continuous_refit.py',ROOT/'validation/mml_pcm_conquest_mc.py',
        ROOT/'validation/run_pcm_conquest_check.py')}
    dump(out/'input.json',dict(classification='OBSERVED_DEVELOPMENT_ONLY', scientific_inference_ready=False,
        qualification_eligible=False, source_summary_sha256=SOURCE_HASH, previous_summary_sha256=SUMMARY_HASH,
        source_sha256=hashes, orders=[31,61,121,181], independent_arithmetic_limit=1e-8,
        gradient_fd_limit=1e-5, integration=original['integration'][1],
        design='Every native history row: reported NLL, current and previous coordinate GH NLL; selected and last gradients',
        scope='One rounded observed fixture; no native refits, source-internal proof, SE/CI or coverage qualification'))
    problem=PCMIntegral(original['data'])
    results={}
    traces={}
    for q in (31,61,121,181):
        name=f'cq_gh{q}'
        native=source['fits'][name]
        history=rows(SOURCE/name/'history.csv')
        assert [int(r['Iteration']) for r in history]==list(range(1,len(history)+1))
        points=np.array([coordinates(r) for r in history])
        selected=native['selected_iteration']-1
        assert np.array_equal(points[selected],native['coordinates'])
        z,w=roots_hermitenorm(q); w=w/np.sqrt(2*np.pi)
        assert np.all(w>0) and abs(w.sum()-1)<1e-14
        evaluate=lambda par:finite_gh(problem,par,z,w)
        checks=fd_check(evaluate,points[selected])
        assert max(c['max_difference'] for c in checks)<1e-5
        values=[evaluate(p) for p in points]
        independent=json.loads((SOURCE/f'{name}_r.json').read_text())['finite'][str(q)]
        r_differences={k:float(np.max(abs(np.asarray(values[selected][k])-independent[k]))) for k in ('nll','eap','sd')}
        assert max(r_differences.values())<1e-8
        app=app_fits[f'q{q}_from_q31']
        app_value=evaluate(app['coordinates'])
        app_differences={k:float(np.max(abs(np.asarray(app_value[k])-app['finite'][str(q)][k]))) for k in ('nll','eap','sd')}
        assert max(app_differences.values())<1e-8
        reported=np.array([float(r['LogLikelihood'])/2 for r in history])
        current=np.array([v['nll'] for v in values])
        lagged=current[:-1]-reported[1:]
        same=current[1:]-reported[1:]
        trace=[dict(iteration=i+1,sigma=float(np.exp(points[i,-1])),reported_nll=float(reported[i]),
            current_coordinate_nll=float(current[i]),
            previous_coordinate_nll=None if i==0 else float(current[i-1]),
            gradient_supnorm=float(np.max(abs(v['gradient']))),log_sd_gradient=float(v['gradient'][-1]))
            for i,v in enumerate(values)]
        dump(out/f'{name}_trace.json',trace); traces[q]=trace
        endpoints={}
        for role,index in dict(selected=selected,last=len(history)-1).items():
            continuous=problem.evaluate(points[index],**original['integration'][1])
            endpoints[role]=dict(iteration=index+1,coordinates=points[index],finite=values[index],continuous=continuous,
                coordinate_distance_to_app=float(np.max(abs(points[index]-app['coordinates']))),
                finite_nll_excess_to_app=float(current[index]-app_value['nll']))
        best=int(np.argmin(current))
        results[name]=dict(selected_iteration=selected+1,executed_rows=len(history),termination=native['termination'],
            native_minimum_reported_iteration=int(np.argmin(reported))+1,
            minimum_recomputed_iteration=best+1,minimum_recomputed_nll=float(current[best]),
            selected_nll_excess_to_best_history=float(current[selected]-current[best]),
            same_row_error_max=float(np.max(abs(same))),same_row_error_median=float(np.median(abs(same))),
            previous_row_error_max=float(np.max(abs(lagged))),previous_row_error_median=float(np.median(abs(lagged))),
            comparison_rows=len(lagged),r_differences=r_differences,gradient_fd=checks,
            app_differences=app_differences,app_finite=app_value,app_sigma=app['sigma'],endpoints=endpoints,
            last_step_coordinate_change=float(np.max(abs(points[-1]-points[-2]))),
            tail20_sigma_range=float(np.ptp(np.exp(points[-20:,-1]))))
        print(name,{k:results[name][k] for k in ('selected_iteration','executed_rows','previous_row_error_max','same_row_error_max')},flush=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(9,6),layout='constrained')
    for ax,q in zip(axes.flat,(31,61,121,181)):
        trace=traces[q]; r=results[f'cq_gh{q}']; selected=r['selected_iteration']
        ax.plot([t['iteration'] for t in trace],[t['sigma'] for t in trace],color='#17617d',label='Native history')
        ax.plot(selected,trace[selected-1]['sigma'],'o',color='#b54b36',label='Exported checkpoint')
        ax.axhline(r['app_sigma'],color='#777777',linestyle='--',label='Direct finite-GH solution')
        ax.set(title=f'GH {q} points',xlabel='Iteration',ylabel='Population SD (logits)')
        ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=7)
    fig.suptitle('ConQuest history: stopping and exported checkpoints')
    fig.savefig(out/'gh_history.png',dpi=180); plt.close(fig)
    assert hashes=={p:sha(ROOT/p) for p in hashes}
    dump(out/'summary.json',dict(classification='OBSERVED_DEVELOPMENT_ONLY',scientific_inference_ready=False,
        qualification_eligible=False,independent_arithmetic_checks_pass=True,results=results,
        source_sha256=hashes,artifact_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()},
        limitations=['History coordinates and deviance are rounded; lagged differences include that error',
            'The first history row has no preceding saved coordinates and is excluded from lag comparisons',
            'One native start per GH order; no new ConQuest fitting or global-optimum claim',
            'Finite-objective stationarity does not establish continuous integration accuracy']))
    print('Saved',out)


if __name__=='__main__':
    main()
