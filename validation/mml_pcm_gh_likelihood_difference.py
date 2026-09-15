#!/usr/bin/env python3
"""Audit PCM likelihood differences with SD-dependent GH nodes; no fitting."""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import platform
import subprocess

import mpmath as mp
import numpy as np
import scipy
from scipy.special import logsumexp, roots_hermitenorm

from mml_pcm_continuous_refit import ROOT, dump, sha
from mml_pcm_independent_conditions import PCMIntegral, difference
from mml_pcm_conquest_history_audit import finite_gh
from mml_pcm_likelihood_difference import _log_average_exp, self_check

PREVIOUS = ROOT/'validation/generated/mml_pcm_two_stage_20260914'
PREVIOUS_HASH = 'c3558407a71b08d2ea555d860904f3750f0b23b41e1e537e7465bb5cd8018c15'
PROTOCOL = dict(
    id='pcm_gh_likelihood_difference_20260914_v1', classification='DEVELOPMENT_ONLY',
    qualification_eligible=False, scientific_inference_ready=False,
    question='Does the exact PCM likelihood ratio preserve moving-GH NLL differences and the log-SD derivative?',
    datasets=['ordinary_new_1', 'wide_new_1'], anchor='Retained b20_q801_sd1 endpoint, same coordinates at every GH order; no new fitting or best-order selection',
    orders=[31,61,121,181],
    targets='Identity, log-SD +1e-7 and -0.01, all coordinates +1e-7 or +0.25 times repeated [-1,0,1], and all-zero coordinates',
    measure='Fixed standard-normal z nodes and positive weights; theta=exp(log-SD)*z. No extra normal-density or Jacobian ratio.',
    mp_digits=[60,90], gradient_steps=[1e-4,1e-5,1e-6],
    r_selection='Anchor, log-SD -0.01, joint-large and zero for both datasets, all four orders; independent statmod nodes and scalar PCM probabilities',
    continuous_selection='Same retained anchor at all orders; loose/tight adaptive integration, without GH refitting',
    negative_controls='Freezing theta makes a pure SD move identically zero. Posterior density-moment log-SD score need not equal the finite moving-GH derivative.',
    limits=dict(delta_absolute=1e-18,delta_relative=1e-11,mp_refinement=1e-45,
        gradient_at_1e_minus6=1e-6,log_sd_mp=1e-9,composition=1e-9,
        raw_nll=1e-8,r_nll=1e-8,r_moments=1e-9,continuous_refinement=1e-8),
    integration_targets=dict(nll_per_person=1e-6,eap=1e-4,sd=1e-4),
    scope='Observed complete unit-weight PCM arithmetic only. No optimizer/stopping/API changes, native TAM/ConQuest fits, recovery or SE/CI/coverage claims.',
    failure_policy='Freeze sources and selection before evaluation; retain all failures, no retuned limits or replacement outputs',
)


def gh_likelihood_difference(problem, anchor, z, w):
    """Return NLL(point)-NLL(anchor) with the same finite standard-normal rule."""
    # ponytail: mean-zero, unit-slope complete PCM only; other models need their own ratio audit.
    anchor = np.asarray(anchor,dtype=float).copy()
    z,w = np.asarray(z,dtype=float).copy(),np.asarray(w,dtype=float).copy()
    def point_ok(p):
        if p.shape != (24,) or not np.all(np.isfinite(p)):
            raise ValueError('Expected 24 finite PCM coordinates')
        if not np.log(.05) <= p[-1] <= np.log(10):
            raise ValueError('This diagnostic is restricted to SD in [0.05, 10]')
    point_ok(anchor)
    if (z.ndim != 1 or len(z)<2 or w.shape != z.shape or
        not np.all(np.isfinite(z)) or not np.all(np.diff(z)>0) or
        not np.all(np.isfinite(w)) or not np.all(w>0) or
        abs(w.sum()-1)>1e-12 or abs(w@z)>1e-12 or abs(w@(z*z)-1)>1e-12):
        raise ValueError('Expected increasing standard-normal nodes with finite positive probability weights')
    theta = np.exp(anchor[-1])*z
    logits = (problem.design@anchor[:-1])[:,None,:]+theta[None,:,None]*problem.k
    normalizer = logsumexp(logits,axis=-1)
    log_category = logits-normalizer[:,:,None]
    log_posterior = problem.total[:,None]*theta-normalizer.sum(axis=0)+np.log(w)
    log_posterior -= logsumexp(log_posterior,axis=-1,keepdims=True)
    def delta(point):
        point = np.asarray(point,dtype=float)
        point_ok(point)
        change = point-anchor
        theta_change = theta*np.expm1(change[-1])
        item_change = (problem.design@change[:-1])[:,None,:]+theta_change[None,:,None]*problem.k
        normalizer_change = _log_average_exp(log_category,item_change).sum(axis=0)
        joint_change = ((problem.observed@change[:-1])[:,None]+
                        problem.total[:,None]*theta_change-normalizer_change)
        return -math.fsum(_log_average_exp(log_posterior,joint_change))
    return delta


def mp_gh(problem, point, z, w):
    """Independent scalar PCM expansion and moving-node log-SD derivative."""
    p = [mp.mpf(float(x)) for x in point]
    raters = p[:3]+[-mp.fsum(p[:3])]
    offsets = []
    for r in range(4):
        for c in range(5):
            steps = p[8+3*c:11+3*c]
            steps = steps+[-mp.fsum(steps)]
            offsets.append([-k*(raters[r]+p[3+c])-mp.fsum(steps[:k]) for k in range(5)])
    theta = [mp.exp(p[-1])*mp.mpf(float(v)) for v in z]
    bases,means = [],[]
    for t,weight in zip(theta,w):
        terms = [[mp.exp(o[k]+k*t) for k in range(5)] for o in offsets]
        denominators = [mp.fsum(a) for a in terms]
        bases.append(mp.log(mp.mpf(float(weight)))-mp.fsum(mp.log(a) for a in denominators))
        means.append(mp.fsum(mp.fsum(k*a[k] for k in range(5))/b for a,b in zip(terms,denominators)))
    observed = mp.fsum(offsets[i][int(y[i])] for y in problem.y for i in range(20))
    nll,score = -observed,mp.mpf(0)
    for total,count in Counter(map(int,problem.total)).items():
        mass = [mp.exp(total*t+b) for t,b in zip(theta,bases)]
        denominator = mp.fsum(mass)
        nll -= count*mp.log(denominator)
        score -= count*mp.fsum(a*t*(total-m) for a,t,m in zip(mass,theta,means))/denominator
    return dict(nll=nll,log_sd_gradient=score)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    out = parser.parse_args().output_dir.resolve()
    load = lambda p:json.loads(p.read_text())
    assert sha(PREVIOUS/'summary.json')==PREVIOUS_HASH
    for name,h in load(PREVIOUS/'summary.json')['artifact_sha256'].items():
        assert sha(PREVIOUS/name)==h,name
    for name,h in load(PREVIOUS/'protocol.json')['source_sha256'].items():
        assert sha(name)==h,name
    paths = [Path(__file__).resolve()]+[Path(__file__).with_name('mml_pcm_'+name) for name in
        ['continuous_refit.py','independent_conditions.py','conquest_history_audit.py',
         'likelihood_difference.py','independent_conditions.R']]
    hashes = {str(p):sha(p) for p in paths}
    out.mkdir(parents=True,exist_ok=False)
    dump(out/'protocol.json',dict(PROTOCOL,previous_summary_sha256=PREVIOUS_HASH,source_sha256=hashes,
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__)))
    self_check()
    limits = PROTOCOL['limits']; results = {}; r_jobs = []
    for name in PROTOCOL['datasets']:
        folder = out/name; folder.mkdir()
        inp = load(PREVIOUS/name/'input.json')
        inp['orders'] = PROTOCOL['orders']; dump(folder/'input.json',inp)
        problem = PCMIntegral(inp['data'])
        a = np.array(load(PREVIOUS/name/'b20_q801_sd1.json')['coordinates'])
        small=a.copy(); small[-1]+=1e-7
        shifted=a.copy(); shifted[-1]-=.01
        pattern=np.tile([-1.,0.,1.],8)
        targets=dict(identity=a,logsd_small=small,logsd_minus=shifted,
                     joint_small=a+1e-7*pattern,joint_large=a+.25*pattern,zero=np.zeros(24))
        dump(folder/'r_cases.json',{k:dict(coordinates=targets[k]) for k in ['identity','logsd_minus','joint_large','zero']})
        # R is an independent process; preserve its exit status and all emitted records.
        log = (folder/'r.log').open('x')
        process = subprocess.Popen(['Rscript',str(paths[-1]),str(folder/'input.json'),str(folder/'r_cases.json'),str(folder)],stdout=log,stderr=subprocess.STDOUT)
        r_jobs.append((name,process,log))
        continuous = [problem.evaluate(a,**s) for s in inp['integration']]
        dump(folder/'continuous.json',continuous)
        results[name]={}
        for q in PROTOCOL['orders']:
            z,w = roots_hermitenorm(q); w/=np.sqrt(2*np.pi)
            delta = gh_likelihood_difference(problem,a,z,w)
            # A zero quadrature weight must be rejected, never silently clipped.
            bad=w.copy(); bad[0]=0
            try:
                gh_likelihood_difference(problem,a,z,bad)
            except ValueError:
                pass
            else:
                raise AssertionError('Zero GH weight was accepted')
            baseline=finite_gh(problem,a,z,w)
            values={k:dict(coordinates=p,finite=finite_gh(problem,p,z,w),stable_delta=delta(p),
                naive_delta=finite_gh(problem,p,z,w)['nll']-baseline['nll'],
                reverse_delta=gh_likelihood_difference(problem,p,z,w)(a)) for k,p in targets.items()}
            gradients=[]
            for h in PROTOCOL['gradient_steps']:
                steps=h*np.maximum(1,abs(a))
                fd=np.array([(delta(a+d)-delta(a-d))/(2*s) for d,s in zip(np.diag(steps),steps)])
                gradients.append(dict(step=h,gradient=fd,max_difference=float(np.max(abs(fd-baseline['gradient'])))))
            composition=abs(delta(targets['joint_large'])-delta(shifted)-gh_likelihood_difference(problem,shifted,z,w)(targets['joint_large']))
            density_score=float(np.sum(1-(baseline['eap']**2+baseline['sd']**2)/np.exp(2*a[-1])))
            frozen=shifted.copy();frozen[-1]=a[-1]
            controls=dict(frozen_node_delta=delta(frozen),correct_sd_delta=delta(shifted),
                density_moment_score=density_score,moving_node_score=float(baseline['gradient'][-1]),
                density_score_difference=abs(density_score-baseline['gradient'][-1]))
            record=dict(anchor=a,nodes=z,weights=w,targets=values,gradient_fd=gradients,
                composition_error=composition,negative_controls=controls,
                integration_error=difference(baseline,continuous[-1]),
                continuous_refinement_error=difference(continuous[0],continuous[1]))
            dump(folder/f'q{q}_float64.json',record)
            precise={}
            for digits in PROTOCOL['mp_digits']:
                with mp.workdps(digits):
                    exact={k:mp_gh(problem,p,z,w) for k,p in targets.items()}
                    precise[str(digits)]={k:dict(nll=mp.nstr(v['nll'],digits),
                        delta=mp.nstr(v['nll']-exact['identity']['nll'],digits),
                        log_sd_gradient=mp.nstr(v['log_sd_gradient'],digits)) for k,v in exact.items()}
                dump(folder/f'q{q}_mp_{digits}.json',precise[str(digits)])
            with mp.workdps(90):
                for k,v in values.items():
                    ref=mp.mpf(precise['90'][k]['delta'])
                    tolerance=limits['delta_absolute']+limits['delta_relative']*abs(ref)
                    v.update(mp_error=float(abs(mp.mpf(v['stable_delta'])-ref)),mp_tolerance=float(tolerance),
                        mp_pass=bool(abs(mp.mpf(v['stable_delta'])-ref)<=tolerance),
                        antisymmetry_pass=bool(abs(v['stable_delta']+v['reverse_delta'])<=2*tolerance),
                        raw_nll_error=float(abs(mp.mpf(v['finite']['nll'])-mp.mpf(precise['90'][k]['nll']))),
                        log_sd_mp_error=float(abs(mp.mpf(v['finite']['gradient'][-1])-mp.mpf(precise['90'][k]['log_sd_gradient']))),
                        precision_refinement=max(float(abs(mp.mpf(precise['60'][k][field])-mp.mpf(precise['90'][k][field]))) for field in ('nll','delta','log_sd_gradient')))
            results[name][str(q)]=record;dump(folder/f'q{q}.json',record)
            print(name,q,'MP',all(v['mp_pass'] for v in values.values()),'FD',gradients[-1]['max_difference'],flush=True)
    r_results={}
    for name,process,log in r_jobs:
        status=process.wait();log.close()
        dump(out/name/'r_execution.json',dict(exit_code=status))
        assert status==0,(name,status)
        r_results[name]={}
        for key in ['identity','logsd_minus','joint_large','zero']:
            r=load(out/name/(key+'_r.json'))
            r_results[name][key]={str(q):difference(results[name][str(q)]['targets'][key]['finite'],r['finite'][str(q)]) for q in PROTOCOL['orders']}
        anchor_r=load(out/name/'identity_r.json')
        r_results[name]['continuous']=difference(load(out/name/'continuous.json')[-1],anchor_r['continuous'][-1])
    records=[r for d in results.values() for r in d.values()]
    comparisons=[v for r in records for v in r['targets'].values()]
    r_errors=[v for d in r_results.values() for key,row in d.items() for v in ([row] if key=='continuous' else row.values())]
    checks=dict(scalar_controls=True,zero_weights_rejected=True,
        identity=all(r['targets']['identity']['stable_delta']==0 for r in records),
        mp_precision=all(v['mp_pass'] for v in comparisons),
        antisymmetry=all(v['antisymmetry_pass'] for v in comparisons),
        mp_refinement=all(v['precision_refinement']<limits['mp_refinement'] for v in comparisons),
        raw_nll=all(v['raw_nll_error']<limits['raw_nll'] for v in comparisons),
        log_sd_gradient=all(v['log_sd_mp_error']<limits['log_sd_mp'] for v in comparisons),
        full_gradient_fd=all(r['gradient_fd'][-1]['max_difference']<limits['gradient_at_1e_minus6'] for r in records),
        composition=all(r['composition_error']<limits['composition'] for r in records),
        frozen_node_control=all(r['negative_controls']['frozen_node_delta']==0 and abs(r['negative_controls']['correct_sd_delta'])>1e-6 for r in records),
        density_score_control=any(r['negative_controls']['density_score_difference']>1e-3 for r in records),
        independent_r=all(e['nll']<limits['r_nll'] and e['eap']<limits['r_moments'] and e['sd']<limits['r_moments'] for e in r_errors),
        continuous_refinement=all(max(r['continuous_refinement_error'].values())<limits['continuous_refinement'] for r in records),
        sources_unchanged=all(sha(p)==h for p,h in hashes.items()),previous_unchanged=sha(PREVIOUS/'summary.json')==PREVIOUS_HASH)
    dump(out/'summary.json',dict(classification=PROTOCOL['classification'],qualification_eligible=False,
        scientific_inference_ready=False,arithmetic_checks_pass=all(checks.values()),checks=checks,
        results=results,independent_r=r_results,protocol_sha256=sha(out/'protocol.json'),
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}))
    print('Arithmetic checks:',checks,flush=True)
    return 0 if all(checks.values()) else 1


if __name__=='__main__':
    raise SystemExit(main())
