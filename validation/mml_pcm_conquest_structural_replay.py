#!/usr/bin/env python3
"""Compare explicit structural M-step replays with saved ConQuest GH transitions."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp, roots_hermitenorm

from mml_pcm_conquest_history_audit import coordinates, finite_gh
from mml_pcm_conquest_variance_replay import variance_updates
from mml_pcm_conquest_mc import SOURCE, SOURCE_HASH, rows
from mml_pcm_continuous_refit import PCMIntegral, PREVIOUS, dump
from run_pcm_conquest_check import ROOT, sha

DIRECT=ROOT/'validation/generated/mml_pcm_conquest_direct_start_20260914'
DIRECT_HASH='3bd9c5d2982f18b6183295452bc0d315b140fb4cd0733192cb6921ff53c43338'
NATIVE_ORDER=[*range(3,8),*range(3),*range(8,23)]


def posterior(problem, x, theta, w):
    logits=(problem.design@x)[:,:,None]+problem.k[None,:,None]*theta
    normalizer=logsumexp(logits,axis=1)
    joint=theta[None,:]*problem.total[:,None]+(problem.observed@x)[:,None]-normalizer.sum(0)+np.log(w)
    h=np.exp(joint-logsumexp(joint,axis=1)[:,None])
    assert np.all(np.isfinite(h)) and np.max(abs(h.sum(1)-1))<1e-12
    return h


def conditional(problem, x, theta, counts):
    # Negative expected complete log likelihood, omitting a constant in x.
    logits=(problem.design@x)[:,:,None]+problem.k[None,:,None]*theta
    normalizer=logsumexp(logits,axis=1)
    p=np.exp(logits-normalizer[:,None,:])
    ed=np.einsum('ikq,ikd->iqd',p,problem.design)
    variance=np.einsum('ikq,ikd->iqd',p,problem.design**2)-ed**2
    observed=problem.observed.sum(0)
    value=float(counts@normalizer.sum(0)-observed@x)
    gradient=np.einsum('q,iqd->d',counts,ed)-observed
    diagonal=np.einsum('q,iqd->d',counts,variance)
    assert np.all(np.isfinite(gradient)) and np.all(diagonal>0)
    return value,gradient,diagonal


def replay(problem, x, theta, counts):
    x=x.copy();snapshots={}
    for sweep in range(1,6):
        for d in NATIVE_ORDER:
            for _ in range(10):
                _,gradient,diagonal=conditional(problem,x,theta,counts)
                change=float(np.clip(gradient[d]/diagonal[d],-1,1))
                x[d]-=change
                if abs(change)<1e-8:break
        if sweep in (1,5):snapshots[str(sweep)]=x.copy()
    return snapshots


def fd_check(problem,x,theta,counts):
    _,g,d=conditional(problem,x,theta,counts);checks=[]
    for step in (1e-4,3e-5):
        vg=[];hd=[]
        for j in range(23):
            offset=np.eye(23)[j]*step
            plus=conditional(problem,x+offset,theta,counts);minus=conditional(problem,x-offset,theta,counts)
            vg.append((plus[0]-minus[0])/(2*step));hd.append((plus[1][j]-minus[1][j])/(2*step))
        checks.append(dict(step=step,gradient_difference=float(np.max(abs(g-vg))),
                           hessian_diagonal_difference=float(np.max(abs(d-hd)))))
    assert max(max(c['gradient_difference'],c['hessian_diagonal_difference']) for c in checks)<1e-5
    return checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    out=parser.parse_args().output_dir.resolve()
    assert sorted(NATIVE_ORDER)==list(range(23))
    assert sha(DIRECT/'review/summary.json')==DIRECT_HASH and sha(SOURCE/'review_summary.json')==SOURCE_HASH
    direct=json.loads((DIRECT/'review/summary.json').read_text());cold=json.loads((SOURCE/'review_summary.json').read_text())
    for record,base in ((direct,DIRECT/'review'),(cold,SOURCE)):
        for name,h in record['artifact_sha256'].items():assert sha(base/name)==h,name
    problem=PCMIntegral(json.loads((PREVIOUS/'input.json').read_text())['data'])
    transitions=[]
    for q in (31,61,121,181):
        warm=direct['results'][str(q)];initial=np.array(warm['initial']['coordinates'])
        for origin,base,name,selected in (
            ('direct',DIRECT,f'q{q}_iterations2000',warm['selected_iteration']),
            ('default',SOURCE,f'cq_gh{q}',cold['fits'][f'cq_gh{q}']['selected_iteration'])):
            history=rows(base/name/'history.csv')
            indices=sorted(set(([0] if origin=='direct' else [])+[selected-1,len(history)-1]))
            for i in indices:
                before=initial if i==0 else coordinates(history[i-1])
                transitions.append(dict(id=f'{origin}_q{q}_iteration{i+1}',q=q,before=before,after=coordinates(history[i])))
    out.mkdir(parents=True,exist_ok=False)
    helpers=['mml_pcm_conquest_history_audit.py','mml_pcm_conquest_variance_replay.py','mml_pcm_conquest_mc.py',
             'mml_pcm_continuous_refit.py','run_pcm_conquest_check.py']
    paths=[Path(__file__).resolve(),*[ROOT/'validation'/p for p in helpers]]
    hashes={str(p.relative_to(ROOT)):sha(p) for p in paths}
    dump(out/'input.json',dict(classification='OBSERVED_DEVELOPMENT_ONLY',scientific_inference_ready=False,
        qualification_eligible=False,posthoc=True,transitions=transitions,source_sha256=hashes,
        source_summary_sha256=SOURCE_HASH,direct_summary_sha256=DIRECT_HASH,
        variants=['old_nodes_old_weights','new_nodes_old_weights','new_nodes_recomputed_weights'],
        sweeps=[1,5],coordinate_order=NATIVE_ORDER,max_newton_steps=10,step_cap=1,step_tolerance=1e-8,
        scope='Replay hypotheses from saved output; no claim that replay tolerance or pass count equals native internals'))
    results={}
    for case in transitions:
        q=case['q'];before=case['before'];after=case['after'];x=before[:-1]
        z,w=roots_hermitenorm(q);w/=np.sqrt(2*np.pi)
        theta=np.exp(before[-1])*z;h=posterior(problem,x,theta,w)
        old_finite=finite_gh(problem,before,z,w)
        assert np.max(abs(conditional(problem,x,theta,h.sum(0))[1]-old_finite['gradient'][:-1]))<1e-10
        predicted=variance_updates(old_finite['eap'],old_finite['sd'])['centered_mean']
        new_theta=np.sqrt(predicted)*z
        variants={}
        for label,t,post in (('old_nodes_old_weights',theta,h),('new_nodes_old_weights',new_theta,h),
                            ('new_nodes_recomputed_weights',new_theta,posterior(problem,x,new_theta,w))):
            counts=post.sum(0);initial_value=conditional(problem,x,t,counts)[0]
            native_value,native_g,_=conditional(problem,after[:-1],t,counts)
            estimates=replay(problem,x,t,counts)
            variants[label]=dict(native_conditional_objective_change=native_value-initial_value,
                native_conditional_gradient_supnorm=float(np.max(abs(native_g))),
                sweeps={n:dict(coordinates=p,max_difference=float(np.max(abs(p-after[:-1]))),
                    objective_change=conditional(problem,p,t,counts)[0]-initial_value) for n,p in estimates.items()})
        checks=fd_check(problem,x,new_theta,h.sum(0)) if case['id'].startswith('direct_') and case['id'].endswith('_iteration1') else []
        results[case['id']]=dict(q=q,predicted_variance=predicted,
            native_variance_error=float(np.exp(2*after[-1])-predicted),variants=variants,derivative_checks=checks)
        dump(out/f"{case['id']}.json",results[case['id']])
        print(case['id'],{k:[v['sweeps'][n]['max_difference'] for n in ('1','5')] for k,v in variants.items()},flush=True)
    assert hashes=={p:sha(ROOT/p) for p in hashes}
    dump(out/'summary.json',dict(classification='OBSERVED_DEVELOPMENT_ONLY',scientific_inference_ready=False,
        qualification_eligible=False,results=results,source_sha256=hashes,
        artifact_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()},
        limitations=['One observed PCM fixture and selected transitions from two starts per GH order',
            'Native coordinates rounded; inner stopping details remain unspecified by this replay',
            'Compare conditional objective changes within a variant, not its absolute values across variants',
            'No native engine changes or SE/CI/coverage qualification']))


if __name__=='__main__':main()
