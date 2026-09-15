#!/usr/bin/env python3
"""Fixed-grid PCM likelihood ratios; development arithmetic, no optimizer changes."""
import argparse
import json
import math
from pathlib import Path
import platform

import mpmath as mp
import numpy as np
import scipy
from scipy.special import logsumexp

from mml_pcm_continuous_refit import ROOT, dump, sha
from mml_pcm_external_review import fixed_grid
from mml_pcm_independent_conditions import PCMIntegral
from mml_pcm_optimizer_resolution import mp_nll

PREVIOUS = ROOT/'validation/generated/mml_pcm_optimizer_resolution_20260914'
GRID = ROOT/'validation/generated/mml_pcm_grid_range_density_20260914'
INPUT_HASHES = {
    str(PREVIOUS/'summary.json'): '331d4d5e39cde5d9c008cfb48367aa7ae7a8015dc1d2751fb33b7f7214e290e7',
    str(GRID/'summary.json'): '53e68b8425a6f40347c2b57a0a35a4e229e51a4ba5757970c722e08a0736656d',
}
PROTOCOL = dict(
    id='pcm_likelihood_difference_20260914_v1', classification='OBSERVED_DEVELOPMENT_ONLY',
    qualification_eligible=False, scientific_inference_ready=False,
    question='Can an exact likelihood-ratio identity recover small finite-grid NLL differences without changing the target?',
    anchors=['ordinary_1 failed point', 'ordinary_1 failed point plus 0.01 in coordinate 0',
             'ordinary_1 all coordinates zero', 'wide_1 returned b12_q241 zero-start point'],
    grids='ordinary_1: +/-20, 801 points; wide_1: +/-12, 241 points; unchanged binary nodes',
    targets='Each anchor itself, coordinate 0 plus 0.01, log-SD minus 0.01, all coordinates plus 1e-7 or 0.25 times repeated [-1,0,1]; additionally the two retained Newton/line-query points',
    ratio='Unnormalized fixed-grid normal prior; mean zero, unit slopes, complete PCM crossings. Posterior and conditional log weights are normalized at the anchor.',
    evaluation='log1p(weighted expm1(change)) if every absolute change <= 0.5; otherwise logsumexp(log weights + change) minus logsumexp(log weights). Retain log weights to avoid discarding revived tail mass.',
    mp_digits=[60,90], gradient_steps=[1e-4,1e-5,1e-6,1e-7],
    checks=dict(delta_absolute=1e-18, delta_relative=1e-11, mp_refinement=1e-45,
                gradient_at_1e_minus6=1e-7, composition_absolute=1e-9),
    limits='Observed arithmetic only; no minimization, stopping-policy changes, continuous-integral qualification, missing-data or weighted-model claims',
    failure_policy='Keep every comparison and failed check; no numerical settings adjusted after results',
)


def _log_average_exp(log_weights, change):
    """Exact log expectation; use small-increment arithmetic where it is safe."""
    if not np.all(np.isfinite(change)):
        raise FloatingPointError('Nonfinite log ratio')
    log_weights = log_weights-logsumexp(log_weights,axis=-1,keepdims=True)
    if np.max(abs(change)) <= .5:
        return np.log1p(np.sum(np.exp(log_weights)*np.expm1(change),axis=-1))
    return logsumexp(log_weights+change,axis=-1)-logsumexp(log_weights,axis=-1)


def likelihood_difference(problem, anchor, theta):
    """Return NLL(point)-NLL(anchor), without subtracting rounded total NLLs."""
    # ponytail: fixed nodes and complete unit-weight PCM only; moving GH nodes need a separate identity audit.
    anchor = np.asarray(anchor,dtype=float).copy()
    theta = np.asarray(theta,dtype=float).copy()
    def point_ok(p):
        if p.shape != (24,) or not np.all(np.isfinite(p)):
            raise ValueError('Expected 24 finite PCM coordinates')
        if not np.log(.05) <= p[-1] <= np.log(10):
            raise ValueError('This diagnostic is restricted to SD in [0.05, 10]')
    point_ok(anchor)
    if (theta.ndim != 1 or len(theta)<2 or not np.all(np.isfinite(theta)) or
        not np.all(np.diff(theta)>0) or not np.allclose(np.diff(theta),theta[1]-theta[0],rtol=0,atol=1e-12)):
        raise ValueError('Expected finite, increasing, equally spaced nodes')
    logits = (problem.design@anchor[:-1])[:,None,:]+theta[None,:,None]*problem.k
    normalizer = logsumexp(logits,axis=-1)
    log_category = logits-normalizer[:,:,None]
    variance_ratio = (theta*np.exp(-anchor[-1]))**2
    # Person-specific observed offsets and grid constants cancel from posterior normalization.
    log_posterior = (problem.total[:,None]*theta-normalizer.sum(axis=0)-.5*variance_ratio)
    log_posterior -= logsumexp(log_posterior,axis=-1,keepdims=True)
    def delta(point):
        point = np.asarray(point,dtype=float)
        point_ok(point)
        change = point-anchor
        item_change = (problem.design@change[:-1])[:,None,:]
        normalizer_change = _log_average_exp(log_category,item_change).sum(axis=0)
        prior_change = -change[-1]-.5*variance_ratio*np.expm1(-2*change[-1])
        joint_change = (problem.observed@change[:-1])[:,None]-normalizer_change+prior_change
        return -math.fsum(_log_average_exp(log_posterior,joint_change))
    return delta


def self_check():
    # A log weight too small to exponentiate can still matter after reweighting.
    assert abs(_log_average_exp(np.array([0.,-1000.]),np.array([0.,1000.]))-np.log(2)) < 1e-14
    for c in (0.,1e-20,.5,np.nextafter(.5,1.),-.5,np.nextafter(-.5,-1.)):
        assert abs(_log_average_exp(np.log([.2,.8]),np.full(2,c))-c) < 1e-15


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    out = parser.parse_args().output_dir.resolve()
    load = lambda path: json.loads(path.read_text())
    for path,h in INPUT_HASHES.items():
        assert sha(path)==h,path
        for name,expected in load(Path(path))['artifact_sha256'].items():
            assert sha(Path(path).parent/name)==expected,name
    old_manifest = load(PREVIOUS/'input.json')
    for path,h in old_manifest['source_sha256'].items():
        assert sha(path)==h,path
    paths = [Path(__file__).resolve(),Path(__file__).with_name('mml_pcm_optimizer_resolution.py'),
             Path(__file__).with_name('mml_pcm_independent_conditions.py'),
             Path(__file__).with_name('mml_pcm_external_review.py'),
             Path(__file__).with_name('mml_pcm_continuous_refit.py')]
    hashes = {str(p):sha(p) for p in paths}
    out.mkdir(parents=True,exist_ok=False)
    dump(out/'input.json',dict(PROTOCOL,input_sha256=INPUT_HASHES,source_sha256=hashes,
        environment=dict(python=platform.python_version(),numpy=np.__version__,
                         scipy=scipy.__version__,mpmath=mp.__version__)))
    self_check()
    previous = load(PREVIOUS/'proposals.json')['coordinates']
    x = np.array(previous['initial']); displaced=x.copy(); displaced[0]+=.01
    wide = np.array(load(GRID/'wide_1/b12_q241.json')['fit']['starts']['zero']['coordinates'])
    conditions = dict(failed=('ordinary_1',x,np.linspace(-20,20,801)),
        displaced=('ordinary_1',displaced,np.linspace(-20,20,801)),
        zero=('ordinary_1',np.zeros(24),np.linspace(-20,20,801)),
        wide=('wide_1',wide,np.linspace(-12,12,241)))
    pattern = np.tile([-1.,0.,1.],8)
    results, mp_cache = {}, {}
    def precise(name,problem,point,theta,digits):
        key=(name,tuple(point),digits)
        if key not in mp_cache:
            mp_cache[key]=mp_nll(problem,point,theta)
        return mp_cache[key]
    for label,(name,a,t) in conditions.items():
        problem=PCMIntegral(load(GRID/name/'input.json')['data'])
        delta=likelihood_difference(problem,a,t)
        shifted=a.copy(); shifted[0]+=.01
        sd_shift=a.copy(); sd_shift[-1]-=.01
        targets=dict(identity=a,rater_plus=shifted,logsd_minus=sd_shift,
                     joint_small=a+1e-7*pattern,joint_large=a+.25*pattern)
        if label=='failed':
            targets.update(newton=np.array(previous['newton']),line_query=np.array(previous['line_query']))
        values={k:dict(coordinates=p,stable_delta=delta(p),
            naive_delta=fixed_grid(problem,p,t)['nll']-fixed_grid(problem,a,t)['nll'],
            reverse_delta=likelihood_difference(problem,p,t)(a)) for k,p in targets.items()}
        analytic=fixed_grid(problem,a,t)['gradient']
        gradients=[]
        for h in PROTOCOL['gradient_steps']:
            steps=h*np.maximum(1,abs(a))
            numerical=np.array([(delta(a+d)-delta(a-d))/(2*s) for d,s in zip(np.diag(steps),steps)])
            gradients.append(dict(step=h,gradient=numerical,
                                  max_difference=float(np.max(abs(numerical-analytic)))))
        via=likelihood_difference(problem,shifted,t)(targets['joint_large'])
        composition_error=abs(delta(targets['joint_large'])-delta(shifted)-via)
        dump(out/f'{label}_float64.json',dict(anchor=a,targets=values,gradient_fd=gradients,
                                           composition_error=composition_error))
        for digits in PROTOCOL['mp_digits']:
            with mp.workdps(digits):
                base=precise(name,problem,a,t,digits)
                for k,p in targets.items():
                    value=precise(name,problem,p,t,digits)
                    values[k][f'mp_nll_{digits}']=mp.nstr(value,digits)
                    values[k][f'mp_delta_{digits}']=mp.nstr(value-base,digits)
            dump(out/f'{label}_mp_{digits}.json',values)
            print(label,'precision',digits,'complete',flush=True)
        with mp.workdps(90):
            for v in values.values():
                ref=mp.mpf(v['mp_delta_90'])
                tolerance=mp.mpf(PROTOCOL['checks']['delta_absolute'])+PROTOCOL['checks']['delta_relative']*abs(ref)
                error=abs(mp.mpf(v['stable_delta'])-ref)
                v.update(absolute_error=float(error),tolerance=float(tolerance),
                    precision_pass=bool(error<=tolerance),
                    antisymmetry_pass=bool(abs(v['stable_delta']+v['reverse_delta'])<=2*tolerance),
                    mp_refinement=float(abs(mp.mpf(v['mp_nll_60'])-mp.mpf(v['mp_nll_90']))))
        results[label]=dict(targets=values,gradient_fd=gradients,composition_error=composition_error)
        dump(out/f'{label}.json',results[label])
    checks=dict(scalar_controls=True,identity=all(r['targets']['identity']['stable_delta']==0 for r in results.values()),
        mp_precision=all(v['precision_pass'] for r in results.values() for v in r['targets'].values()),
        antisymmetry=all(v['antisymmetry_pass'] for r in results.values() for v in r['targets'].values()),
        mp_refinement=all(v['mp_refinement']<PROTOCOL['checks']['mp_refinement'] for r in results.values() for v in r['targets'].values()),
        gradient=all(r['gradient_fd'][2]['max_difference']<PROTOCOL['checks']['gradient_at_1e_minus6'] for r in results.values()),
        composition=all(r['composition_error']<PROTOCOL['checks']['composition_absolute'] for r in results.values()),
        inputs_unchanged=all(sha(p)==h for p,h in INPUT_HASHES.items()),
        sources_unchanged=all(sha(p)==h for p,h in hashes.items()),
        historical_failure_retained=load(GRID/'summary.json')['implementation_checks_pass'] is False)
    dump(out/'summary.json',dict(classification=PROTOCOL['classification'],scientific_inference_ready=False,
        qualification_eligible=False,checks=checks,arithmetic_checks_pass=all(checks.values()),
        historical_implementation_checks_pass=False,results=results,protocol_sha256=sha(out/'input.json'),
        artifact_sha256={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}))
    print('Arithmetic checks:',checks,flush=True)
    return 0 if all(checks.values()) else 1


if __name__=='__main__':
    raise SystemExit(main())
