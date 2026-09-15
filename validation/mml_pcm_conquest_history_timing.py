#!/usr/bin/env python3
"""Post-hoc timing audit of retained GH and fixed-grid histories; no native fits."""
import argparse
from decimal import Decimal
import json
from pathlib import Path
import re

import numpy as np
from scipy.special import logsumexp, roots_hermitenorm

from mml_pcm_conquest_history_audit import coordinates, finite_gh
from mml_pcm_conquest_mc import SOURCE, SOURCE_HASH, rows
from mml_pcm_conquest_structural_replay import DIRECT, DIRECT_HASH
from mml_pcm_conquest_variance_replay import HISTORY, HISTORY_HASH
from mml_pcm_continuous_refit import PREVIOUS, ROOT, dump, sha
from mml_pcm_grid_conquest_review import normalized_grid
from mml_pcm_grid_check import nodes
from mml_pcm_independent_conditions import PCMIntegral

GENERATED = ROOT / 'validation/generated'
CHECKPOINT = GENERATED / 'mml_pcm_conquest_checkpoint_check_20260914'
GRID = GENERATED / 'mml_pcm_grid_conquest_20260914'
EPS = 5e-7  # Native xsi and variance are rounded to six decimal places.
DISPLAY = 2.5e-7  # Six-decimal deviance divided by two.
ARITHMETIC = 1e-8


def retained_inputs():
    manifest = {}
    records = [(SOURCE/'review_summary.json', SOURCE_HASH),
        (HISTORY/'summary.json', HISTORY_HASH),
        (DIRECT/'review/summary.json', DIRECT_HASH),
        (DIRECT/'input.json', '540ad2edb643026b47fe121f2235bd525e8acb9ad7138f51fdda5a4357e00547'),
        (DIRECT/'execution.json', '69584f5d1777a10d323a21d98c2899598dc800cc127de42fb3d729054f077c07'),
        (CHECKPOINT/'input.json', 'f5e9599254df533a39071a2c520495a9a707c38604031bacb6befa76a6fa0e71'),
        (CHECKPOINT/'summary.json', 'd4ad89fac3424ac67c4a2d7c74a4cc589f853ee866f68aa0ff6f5fedb2d1dba0'),
        (CHECKPOINT/'review/summary.json', '089b2841314309d600f7d91b00836371fff91be80ccd7cfb669852cf965c23ff'),
        (GRID/'input.json', '2133938d04c1515a3f2a06214958cefa6eb6e80c66d51f3731a7e3ccae54f298'),
        (GRID/'execution.json', '67f0d8aeb6a5ac1b6fae6688450ee0d514c9d6924fb603923aeedc6dd6cd4c18'),
        (GRID/'review/summary.json', '741baefb043bce158023b15b840fc69ba7461ddd7b62e25e76aa8ab75abbee26')]
    def check(path, expected):
        assert sha(path) == expected, path
        key = str(path.relative_to(ROOT))
        assert key not in manifest or manifest[key] == expected
        manifest[key] = expected
    for path, expected in records:
        check(path, expected)
        record = json.loads(path.read_text())
        for key, base in [('artifact_sha256', path.parent), ('input_sha256', path.parent), ('source_sha256', ROOT)]:
            for name, digest in record.get(key, {}).items():
                check(base/name, digest)
    # PREVIOUS/input.json is also bound by the original native review artifacts.
    data = PREVIOUS/'input.json'
    previous_summary = json.loads((PREVIOUS/'summary.json').read_text())
    check(data, previous_summary['artifact_sha256']['input.json'])
    return manifest


def enclosure(problem, point, rule, epsilon):
    """Conservative NLL range for a native xsi/variance rounding box.

    Bound every category logit contrast, then log-sum-exp, conditional likelihood
    products and the positive mixture. This is an algebraic interval enclosure
    evaluated in float64, with a separate arithmetic guard; not directed rounding.
    """
    variance = np.exp(2*point[-1])
    assert variance > epsilon >= 0
    if 'bound' in rule:
        theta = nodes(rule)
        theta_radius = np.zeros_like(theta)
        raw_lo = -.5*theta**2/(variance-epsilon)
        raw_hi = -.5*theta**2/(variance+epsilon)
        logw_lo, logw_hi = raw_lo-logsumexp(raw_hi), raw_hi-logsumexp(raw_lo)
    else:
        z, w = rule['z'], rule['w']
        theta = np.sqrt(variance)*z
        theta_radius = abs(z)*max(np.sqrt(variance+epsilon)-np.sqrt(variance),
                                  np.sqrt(variance)-np.sqrt(variance-epsilon))
        logw_lo = logw_hi = np.log(w)
    # Axes: item, observed category, alternative category, [coordinate or node].
    contrast = problem.design[:, None, :, :] - problem.design[:, :, None, :]
    category = problem.k[None, :] - problem.k[:, None]
    center = (contrast@point[:-1])[..., None] + category[None, ..., None]*theta
    radius = epsilon*abs(contrast).sum(-1)[..., None] + abs(category)[None, ..., None]*theta_radius
    logp_lo = -logsumexp(center+radius, axis=2)
    logp_hi = -logsumexp(center-radius, axis=2)
    index = np.arange(problem.y.shape[1])[None, :]
    joint_lo = logp_lo[index, problem.y].sum(1) + logw_lo
    joint_hi = logp_hi[index, problem.y].sum(1) + logw_hi
    return np.array([-logsumexp(joint_hi, axis=1).sum(), -logsumexp(joint_lo, axis=1).sum()])


def plan():
    data = json.loads((PREVIOUS/'input.json').read_text())['data']
    starts = json.loads((DIRECT/'input.json').read_text())['starts']
    cases = []
    for q in (31, 61, 121, 181):
        for group, folder, start in [('gh_default', SOURCE/f'cq_gh{q}', None),
                ('gh_direct', DIRECT/f'q{q}_iterations2000', starts[str(q)]['coordinates']),
                ('gh_one', DIRECT/f'q{q}_iterations1', starts[str(q)]['coordinates'])]:
            cases.append(dict(id=f'{group}_{q}', group=group, folder=folder, q=q, start=start, data=data))
    for name in ('full_keep_no', 'full_keep_yes', 'three_keep_no', 'five_keep_no'):
        cases.append(dict(id=name, group='gh_checkpoint', folder=CHECKPOINT/name, q=61,
                          start=starts['61']['coordinates'], data=data))
    for dataset in ('ordinary_1', 'ordinary_2', 'wide_1', 'wide_2'):
        data = json.loads((GRID/'review'/dataset/'input.json').read_text())['data']
        for bound in (12, 20):
            cases.append(dict(id=f'{dataset}_b{bound}', group='fixed_grid', dataset=dataset,
                folder=GRID/f'{dataset}_b{bound}', bound=bound, q=bound*20+1, spacing=.1,
                start=np.zeros(24), data=data))
    return cases


def audit_case(case, out):
    problem = PCMIntegral(case['data'])
    history = rows(case['folder']/'history.csv')
    assert [int(r['Iteration']) for r in history] == list(range(1, len(history)+1))
    for row in history:
        assert all(Decimal(row[k]).as_tuple().exponent >= -6
                   for k in ['LogLikelihood', 'wvar 1 1', *[f'xsi {i}' for i in range(1,24)]])
    points = [coordinates(r) for r in history]
    rule = {k:case[k] for k in ('bound','q','spacing')} if 'bound' in case else {}
    if rule:
        evaluate = lambda p: normalized_grid(problem, p, rule)
    else:
        z, w = roots_hermitenorm(case['q']); w /= np.sqrt(2*np.pi)
        assert np.all(w > 0) and abs(w.sum()-1) < 1e-14
        rule = dict(z=z, w=w)
        evaluate = lambda p: finite_gh(problem, p, z, w)
    report = (case['folder']/'review.txt').read_text()
    selected = int(re.search(r'The number of iterations:\s*(\d+)', report)[1])-1
    xsi = np.array([float(r['Estimate']) for r in rows(case['folder']/'parameters.csv')])
    assert np.array_equal(points[selected][:-1], np.r_[xsi[5:8], xsi[:5], xsi[8:]])
    assert float(history[selected]['wvar 1 1']) == float(rows(case['folder']/'covariance.csv')[0]['Covariance'])
    final_nll = float(re.search(r'Final Deviance:\s*([\d.]+)', report)[1])/2
    assert abs(final_nll-float(history[selected]['LogLikelihood'])/2) < 2.76e-6
    values = [evaluate(p) for p in points]
    boxes = [enclosure(problem, p, rule, EPS) for p in points]
    # Runnable arithmetic controls: zero-radius identity and opposing box corners.
    point = points[selected]
    zero_error = float(np.max(abs(enclosure(problem, point, rule, 0)-values[selected]['nll'])))
    assert zero_error < ARITHMETIC
    for sign in (-1, 1):
        perturbed = point.copy()
        perturbed[:-1] += sign*EPS*np.where(np.arange(23) % 2, 1, -1)
        perturbed[-1] = .5*np.log(np.exp(2*point[-1])+sign*EPS)
        assert boxes[selected][0]-ARITHMETIC <= evaluate(perturbed)['nll'] <= boxes[selected][1]+ARITHMETIC
    start = None if case['start'] is None else np.array(case['start'])
    before = None if start is None else evaluate(start)
    before_box = None if start is None else enclosure(problem, start, rule, 0)
    trace = []
    for i, (row, p, value, box) in enumerate(zip(history, points, values, boxes)):
        native = float(row['LogLikelihood'])/2
        after_ok = bool(box[0]-DISPLAY-ARITHMETIC <= native <= box[1]+DISPLAY+ARITHMETIC)
        record = dict(iteration=i+1, reported_nll=native, after_nll=value['nll'], after_box=box,
            same_row_error=native-value['nll'], after_compatible=after_ok)
        if before is not None:
            old = start if i == 0 else points[i-1]
            epsilon = 0 if i == 0 else EPS
            log_sd_radius = .5*np.log(np.exp(2*old[-1])/(np.exp(2*old[-1])-epsilon))
            linear = DISPLAY+ARITHMETIC+epsilon*np.sum(abs(before['gradient'][:-1]))+log_sd_radius*abs(before['gradient'][-1])
            record.update(before_nll=before['nll'], before_box=before_box,
                previous_row_error=native-before['nll'], first_order_allowance=linear,
                before_compatible=bool(before_box[0]-DISPLAY-ARITHMETIC <= native <= before_box[1]+DISPLAY+ARITHMETIC),
                within_first_order_allowance=bool(abs(native-before['nll']) <= linear))
        trace.append(record)
        before, before_box = value, box
    eligible = [r for r in trace if 'before_nll' in r]
    # Reuse the independent R values already retained at the returned coordinate.
    if case['group'] == 'gh_default':
        reference = json.loads((SOURCE/f'cq_gh{case["q"]}_r.json').read_text())['finite'][str(case['q'])]
    elif case['group'] in ('gh_direct', 'gh_one'):
        reference = json.loads((DIRECT/f'review/q{case["q"]}_iteration{selected+1}_r.json').read_text())['finite'][str(case['q'])]
    elif case['group'] == 'gh_checkpoint':
        reference = json.loads((CHECKPOINT/f'review/{case["id"]}_r.json').read_text())['finite'][str(case['q'])]
    else:
        reference = json.loads((GRID/f'review/{case["dataset"]}/b{case["bound"]}_returned_r.json').read_text())['finite'][f'b{case["bound"]}_q{case["q"]}']
        reference = dict(reference, nll=reference['normalized_prior_nll'])
    r_error = max(float(np.max(abs(np.asarray(values[selected][k])-reference[k]))) for k in ('nll','eap','sd'))
    assert r_error < ARITHMETIC
    dump(out/f'{case["id"]}.json', trace)
    return dict(group=case['group'], q=case['q'], executed_rows=len(history), compared_rows=len(eligible),
        before_compatible_rows=sum(r['before_compatible'] for r in eligible),
        same_row_rejected_rows=sum(not r['after_compatible'] for r in eligible),
        first_order_compatible_rows=sum(r['within_first_order_allowance'] for r in eligible),
        max_before_error=max(abs(r['previous_row_error']) for r in eligible),
        max_same_row_error=max(abs(r['same_row_error']) for r in eligible),
        max_first_order_ratio=max(abs(r['previous_row_error'])/r['first_order_allowance'] for r in eligible),
        max_before_box_width=max(float(np.ptp(r['before_box'])) for r in eligible),
        selected_iteration=selected+1, selected=trace[selected], final_report_nll=final_nll,
        returned_gradient_supnorm=float(np.max(abs(values[selected]['gradient']))),
        zero_radius_error=zero_error, retained_independent_r_error=r_error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--audit', action='store_true')
    args = parser.parse_args(); out = args.output_dir.resolve()
    manifest = retained_inputs()
    if args.audit:
        saved = json.loads((out/'summary.json').read_text())
        assert saved['retained_inputs'] == manifest
        for name, digest in saved['artifact_sha256'].items(): assert sha(out/name) == digest, name
        for name, digest in saved['source_sha256'].items(): assert sha(ROOT/name) == digest, name
        print('Integrity replay:', len(manifest), 'retained files;', len(saved['artifact_sha256']), 'new artifacts')
        return
    out.mkdir(parents=True, exist_ok=False)
    helpers = ('mml_pcm_conquest_history_audit.py','mml_pcm_grid_conquest_review.py',
        'mml_pcm_independent_conditions.py','mml_pcm_external_review.py','mml_pcm_grid_check.py',
        'mml_pcm_continuous_refit.py','mml_pcm_conquest_mc.py','mml_pcm_conquest_structural_replay.py',
        'mml_pcm_conquest_variance_replay.py','run_pcm_conquest_check.py')
    sources = {str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__).resolve(), *[ROOT/'validation'/n for n in helpers]]}
    cases = plan()
    dump(out/'protocol.json', dict(classification='POST_HOC_NUMERICAL_DIAGNOSIS', scientific_inference_ready=False,
        qualification_eligible=False, question='Are older displayed NLLs compatible with before or after coordinates once export rounding is propagated?',
        cases=[{k:v for k,v in c.items() if k not in ('data','folder')} for c in cases],
        retained_inputs=manifest, source_sha256=sources,
        limits=dict(native_coordinate_half_width=EPS, nll_display_half_width=DISPLAY, float64_guard=ARITHMETIC),
        method='All retained rows; unknown initial points excluded only from before/after comparisons. Category-contrast interval bounds, normalized finite prior; first-order sensitivity is a separate diagnostic.',
        scope='One 24-person GH fixture with duplicate controls; four 80-person grid datasets. No native execution, timing uniqueness near convergence, source-internal proof or inference qualification.'))
    results = {}
    for case in cases:
        results[case['id']] = audit_case(case, out)
        print(case['id'], {k:results[case['id']][k] for k in ('compared_rows','before_compatible_rows','same_row_rejected_rows','max_before_error','max_first_order_ratio')}, flush=True)
    assert retained_inputs() == manifest
    assert sources == {name:sha(ROOT/name) for name in sources}
    dump(out/'summary.json', dict(classification='POST_HOC_NUMERICAL_DIAGNOSIS', scientific_inference_ready=False,
        qualification_eligible=False, results=results, retained_inputs=manifest, source_sha256=sources,
        all_before_compatible=all(r['before_compatible_rows'] == r['compared_rows'] for r in results.values()),
        arithmetic_controls_pass=True, artifact_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()}))
    print('Saved', out, flush=True)


if __name__ == '__main__':
    main()
