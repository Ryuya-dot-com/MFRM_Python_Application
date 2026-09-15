#!/usr/bin/env python3
"""Audit eight retained ConQuest regression/missingness histories without refitting."""
import argparse
from decimal import Decimal
import json
from pathlib import Path
import re
from types import SimpleNamespace

import numpy as np
from scipy.special import logsumexp

from mml_observed_pcm_refit import ObservedPCMReference
from mml_pcm_continuous_refit import ROOT, dump, sha
from mml_pcm_conquest_mc import rows
from mml_unit_weight_native import CELLS, app_parameters, native_parameters, grid
from mml_unit_weight_conquest_mc import PRIOR, PRIOR_SHA, verify_prior, verify_hashes

EPS, DISPLAY, ARITHMETIC = 5e-7, 2.5e-7, 1e-8


def app_point(u):
    return app_parameters(u[:8], [0., u[8]], u[9])


def history_point(row):
    assert float(row['Dim 1 Var 1']) == 0
    return np.array([float(row[f'xsi {i}']) for i in range(1, 9)] +
                    [float(row['Dim 2 Var 2']), float(row['wvar 1 1'])])


def evaluator(ref, bound):
    """Use native (eight xsi, slope, variance) coordinates for rounding and scores."""
    assert ref.n == 32 and np.all(ref.weight == 1)
    theta = np.linspace(-bound, bound, bound*20+1)
    base = app_point(np.r_[np.zeros(9), 1.])[:-1]
    basis = np.column_stack([app_point(np.r_[np.eye(8)[j], 0., 1.])[:-1]-base for j in range(8)])
    design = ref.design @ basis
    offset = ref.offset + ref.design @ base
    assert np.max(abs(offset)) < 1e-14
    observed = design[np.arange(len(ref.y)), ref.y]
    contrast = design-observed[:, None, :]
    offset_contrast = offset-offset[np.arange(len(ref.y)), ref.y, None]
    category = ref.k[None, :]-ref.y[:, None]
    # Jacobian from app structural coordinates to native (xsi, slope).
    zero_native = np.r_[native_parameters(np.zeros(10)), 0.]
    jacobian = np.column_stack([np.r_[native_parameters(p), p[0]]-zero_native for p in np.eye(10)[:9]])

    def evaluate(u):
        assert u.shape == (10,) and np.all(np.isfinite(u)) and u[-1] > 0
        mu = ref.x*u[8]; variance = u[-1]
        logits = (offset+design@u[:8])[:, None, :]+theta[None, :, None]*ref.k
        logp = logits-logsumexp(logits, axis=2, keepdims=True)
        prob = np.exp(logp)
        ll = np.zeros((ref.n, len(theta)))
        np.add.at(ll, ref.person, logp[np.arange(len(ref.y)), :, ref.y])
        distance = theta[None, :]-mu[:, None]
        logw = -.5*distance**2/variance-.5*np.log(2*np.pi*variance)+np.log(.1)
        mass = logsumexp(ll+logw, axis=1)
        prior_mass = logsumexp(logw, axis=1)
        post, prior = np.exp(ll+logw-mass[:, None]), np.exp(logw-prior_mass[:, None])
        expected = np.einsum('iqk,ikd->iqd', prob, design)
        structural = np.einsum('iq,iqd->d', post[ref.person], expected)-observed.sum(0)
        density_score = np.stack((ref.x[:, None]*distance/variance,
                                 -.5/variance+.5*distance**2/variance**2), axis=2)
        raw_gradient = np.r_[structural, -np.einsum('pq,pqd->d', post, density_score)]
        correction = np.r_[np.zeros(8), np.einsum('pq,pqd->d', prior, density_score)]
        gradient = raw_gradient+correction
        eap = post@theta
        transform = lambda g: np.r_[jacobian.T@g[:9], 2*variance*g[-1]]
        return dict(nll=float(-mass.sum()+prior_mass.sum()), raw_nll=float(-mass.sum()),
            log_prior_mass=prior_mass, gradient=gradient, raw_gradient=raw_gradient,
            app_gradient=transform(gradient), raw_app_gradient=transform(raw_gradient),
            eap=eap, sd=np.sqrt(np.sum(post*(theta[None, :]-eap[:, None])**2, axis=1)))

    def enclosure(u, epsilon):
        # Logit contrasts bound missing-data likelihoods one observed response at a time.
        assert u[-1] > epsilon >= 0
        center = (offset_contrast+contrast@u[:8])[:, :, None]+category[:, :, None]*theta
        radius = epsilon*abs(contrast).sum(2)[:, :, None]
        conditional = [-logsumexp(center+sign*radius, axis=1) for sign in (1, -1)]
        ll = [np.zeros((ref.n, len(theta))) for _ in range(2)]
        for a, b in zip(ll, conditional): np.add.at(a, ref.person, b)
        distance = abs(theta[None, :]-ref.x[:, None]*u[8])
        mean_radius = abs(ref.x[:, None])*epsilon
        # The factor 1/sqrt(2*pi*v) cancels when normalizing each person's prior.
        raw_lo = -.5*(distance+mean_radius)**2/(u[-1]-epsilon)
        raw_hi = -.5*np.maximum(distance-mean_radius, 0.)**2/(u[-1]+epsilon)
        logw_lo = raw_lo-logsumexp(raw_hi, axis=1, keepdims=True)
        logw_hi = raw_hi-logsumexp(raw_lo, axis=1, keepdims=True)
        return np.array([-logsumexp(ll[1]+logw_hi, axis=1).sum(),
                         -logsumexp(ll[0]+logw_lo, axis=1).sum()])
    return evaluate, enclosure


def controls(ref, bound, u, evaluate, enclosure):
    value = evaluate(u)
    old = grid(ref, app_point(u), bound)
    replay = max(abs(value['nll']-old['normalized_nll']), abs(value['raw_nll']-old['nll']),
                 *[float(np.max(abs(value[k]-old[k]))) for k in ('eap','sd','log_prior_mass')])
    zero = float(np.max(abs(enclosure(u, 0)-value['nll'])))
    assert max(replay, zero) < ARITHMETIC
    box = enclosure(u, EPS)
    for sign in (-1, 1):
        corner = u+sign*EPS*np.where(np.arange(10) % 2, 1, -1)
        assert box[0]-ARITHMETIC <= evaluate(corner)['nll'] <= box[1]+ARITHMETIC
    checks = []
    for h in (1e-4, 3e-5):
        raw, normalized = [], []
        for direction in np.eye(10):
            plus = grid(ref, app_point(u+h*direction), bound)
            minus = grid(ref, app_point(u-h*direction), bound)
            raw.append((plus['nll']-minus['nll'])/(2*h))
            normalized.append((plus['normalized_nll']-minus['normalized_nll'])/(2*h))
        checks.append(dict(step=h, raw_error=float(np.max(abs(value['raw_gradient']-raw))),
                           normalized_error=float(np.max(abs(value['gradient']-normalized)))))
    assert max(max(c['raw_error'], c['normalized_error']) for c in checks) < 1e-5
    return dict(replay_error=replay, zero_radius_error=zero, gradient_fd=checks)


def review(cell, bound, old_spec, out):
    folder = PRIOR/cell; native = folder/f'cq_b{bound}_refit'
    data = json.loads((folder/'input.json').read_text())
    ref = ObservedPCMReference(SimpleNamespace(idx=data['indices'],
        config=dict(n_person=32, population_model=dict(X=np.array(data['x'])[:, None]))))
    evaluate, enclosure = evaluator(ref, bound)
    p0 = np.array(old_spec['start']); initial = np.r_[native_parameters(p0), p0[0], np.exp(2*p0[-1])]
    assert np.max(abs(app_point(initial)-p0)) < 1e-14
    imported = np.r_[np.loadtxt(folder/'init_parameters.txt')[:, 1],
                     np.loadtxt(folder/'init_beta.txt')[-1], np.loadtxt(folder/'init_covariance.txt')[-1]]
    assert np.array_equal(initial, imported)
    history = rows(native/'history.csv')
    assert [int(r['Iteration']) for r in history] == list(range(1, len(history)+1))
    fields = ['LogLikelihood','Dim 2 Var 2','wvar 1 1', *[f'xsi {i}' for i in range(1,9)]]
    assert all(Decimal(row[k]).as_tuple().exponent >= -6 for row in history for k in fields)
    points = [history_point(r) for r in history]
    report = (native/'review.txt').read_text()
    selected = int(re.search(r'The number of iterations:\s*(\d+)', report)[1])-1
    assert 0 <= selected < len(history)
    beta = rows(native/'reg_coefficients.csv'); assert float(beta[0]['Estimate']) == 0
    returned = np.array([float(r['Estimate']) for r in rows(native/'parameters.csv')] +
                       [float(beta[1]['Estimate']), float(rows(native/'covariance.csv')[0]['Covariance'])])
    assert np.array_equal(points[selected], returned)
    assert int(re.search(r'Total number of estimated parameters:\s*(-?\d+)', report)[1]) == 10
    final_display = float(re.search(r'Final Deviance:\s*([\d.]+)', report)[1])/2
    assert abs(final_display-float(history[selected]['LogLikelihood'])/2) < 2.76e-6
    checked = {role:controls(ref, bound, p, evaluate, enclosure) for role, p in [('initial',initial),('returned',returned)]}
    before = evaluate(initial); before_box = enclosure(initial, 0); trace = []
    for i, (row, point) in enumerate(zip(history, points)):
        value, box = evaluate(point), enclosure(point, EPS)
        reported = float(row['LogLikelihood'])/2
        linear = DISPLAY+ARITHMETIC+(EPS*np.sum(abs(before['gradient'])) if i else 0.)
        compatible = lambda interval: bool(interval[0]-DISPLAY-ARITHMETIC <= reported <= interval[1]+DISPLAY+ARITHMETIC)
        trace.append(dict(iteration=i+1, reported_nll=reported, before_nll=before['nll'], after_nll=value['nll'],
            before_box=before_box, after_box=box, before_error=reported-before['nll'], after_error=reported-value['nll'],
            before_compatible=compatible(before_box), after_compatible=compatible(box), first_order_allowance=linear,
            within_first_order_allowance=bool(abs(reported-before['nll']) <= linear)))
        before, before_box = value, box
    end = evaluate(returned)
    independent = json.loads((folder/'r_grid.json').read_text())[f'cq_b{bound}_refit']
    r_error = max(abs(end['nll']-independent['normalized_nll']), abs(end['raw_nll']-independent['nll']),
        *[float(np.max(abs(end[k]-independent[k]))) for k in ('eap','sd','log_prior_mass')])
    assert r_error < ARITHMETIC
    previous = json.loads((folder/f'b{bound}_refit_review.json').read_text())['endpoints']['cq']
    assert np.array_equal(app_point(returned), previous['coordinates'])
    dump(out/f'{cell}_b{bound}_trace.json', trace)
    return dict(rows=len(trace), before_compatible_rows=sum(t['before_compatible'] for t in trace),
        same_row_rejected_rows=sum(not t['after_compatible'] for t in trace),
        first_order_compatible_rows=sum(t['within_first_order_allowance'] for t in trace),
        max_before_error=max(abs(t['before_error']) for t in trace), max_after_error=max(abs(t['after_error']) for t in trace),
        max_first_order_ratio=max(abs(t['before_error'])/t['first_order_allowance'] for t in trace),
        selected_iteration=selected+1, returned_equals_last=bool(np.array_equal(returned, points[-1])),
        initial_native=initial, returned_native=returned, returned_app=app_point(returned),
        final_display_nll=final_display, selected=trace[selected], returned=end,
        termination=re.search(r'Iterations terminated[^\n]+', report)[0],
        retained_continuous=previous['continuous'], retained_accuracy_pass=previous['accuracy_pass'],
        retained_total_score_error=previous['total_score_error'], retained_independent_r_error=r_error, controls=checked)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--audit', action='store_true')
    args = parser.parse_args(); out = args.output_dir.resolve()
    prior, old_spec = verify_prior()
    if args.audit:
        saved = json.loads((out/'summary.json').read_text()); spec = json.loads((out/'protocol.json').read_text())
        assert saved['protocol_sha256'] == sha(out/'protocol.json')
        verify_hashes(ROOT, spec['source_sha256']); verify_hashes(out, saved['artifact_sha256'])
        print('Integrity replay:', len(prior['artifact_sha256']), 'prior artifacts;', len(saved['artifact_sha256']), 'new artifacts')
        return
    sources = dict(old_spec['source_sha256'])
    for p in [Path(__file__).resolve(), ROOT/'validation/mml_unit_weight_conquest_mc.py']:
        sources[str(p.relative_to(ROOT))] = sha(p)
    out.mkdir(parents=True, exist_ok=False)
    dump(out/'protocol.json', dict(classification='POST_HOC_NUMERICAL_DIAGNOSIS', scientific_inference_ready=False,
        qualification_eligible=False, prior_summary_sha256=PRIOR_SHA, source_sha256=sources,
        cells=CELLS, bounds=[12,20], start=old_spec['start'],
        question='Does before-update NLL timing hold when the regression slope and variance change, with incomplete responses and native-coordinate rounding?',
        design='All eight retained refit histories; exact imported first point, then previous rounded native coordinates. Unit weights, mean beta*x, variance free, normalized finite prior per person.',
        rounding=dict(native_half_width=EPS, nll_display_half_width=DISPLAY, arithmetic_guard=ARITHMETIC),
        checks='Native affine mapping; existing grid evaluator, opposing box corners and zero-radius identity; two finite-difference steps for all 10 native coordinates at initial/returned points; saved independent R at returned points.',
        gradient_fd_limit=1e-5, method='Conservative category-contrast and per-person prior bounds; first-order score allowance separate. Float64 with guard, not directed-rounding intervals.',
        scope='Two paired response datasets with complete/missing designs; no new native runs, stopping or continuous-accuracy qualification, recovery, SE/CI/coverage or default changes.'))
    results = {}
    for cell in CELLS:
        for bound in (12,20):
            key = f'{cell}_b{bound}'; results[key] = review(cell, bound, old_spec, out)
            print(key, {k:results[key][k] for k in ('rows','before_compatible_rows','same_row_rejected_rows','max_before_error','max_first_order_ratio')}, flush=True)
    verify_prior(); verify_hashes(ROOT, sources)
    dump(out/'summary.json', dict(classification='POST_HOC_NUMERICAL_DIAGNOSIS', scientific_inference_ready=False,
        qualification_eligible=False, protocol_sha256=sha(out/'protocol.json'), results=results,
        all_before_compatible=all(r['before_compatible_rows']==r['rows'] for r in results.values()),
        arithmetic_controls_pass=True, artifact_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()}))
    print('Saved', out, flush=True)


if __name__ == '__main__':
    main()
