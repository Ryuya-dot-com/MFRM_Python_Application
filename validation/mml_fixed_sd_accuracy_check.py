#!/usr/bin/env python3
"""Offline development check: same-point integration versus GH recalibration.

Reads a restricted PCM data.csv / analysis_snapshot.json pair. Never replaces a
fitted result or qualifies inference. Existing numerical references stay intact.
"""
import argparse
import io
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace
import zipfile

import numpy as np
import pandas as pd
from scipy.special import ndtr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import streamlit_app as app
from mml_fixed_sd_pcm_polish import Reference, PERM, differences
from mml_free_sd_stationarity_adapter import AppFreeSdProblem
from mml_pcm_continuous_refit import dump, sha
from mfrm_app.mml_stationarity import JointPolishOptions, audit_joint_gradient, polish_fixed_sd

PROTOCOL = dict(
    classification='OBSERVED_DEVELOPMENT_ONLY', inference_ready=False,
    question='Separate integration error at the original/returned point from recalibration movement on one continuous-integral scale.',
    orders=[31, 61, 121], starts=['original', 'zero'],
    continuous=[dict(bound=12., rel_tol=1e-10, abs_tol=1e-12),
                dict(bound=16., rel_tol=1e-12, abs_tol=1e-14)],
    optimizer=dict(maxiter=500, gtol=1e-9, ftol=0., maxls=50),
    gradient_audit=dict(offset=np.linspace(-.08, .09, 8).tolist(), relative_step=1e-5),
    engineering_limits=dict(finite_replay=1e-8, gradient_fd=1e-6, reference_refinement=1e-8,
        reference_moment_allowance=1e-8, reference_mass_allowance=1e-7,
        cross_start_coordinates=1e-7),
    failure_policy='Retain all six trajectories and every solver status; no retuning, replacement or best-start selection.',
    interpretation='The baseline is continuous evaluation of the ORIGINAL returned point, not a continuous-likelihood optimum. Signed per-person differences add; maxima do not. No SE/CI/coverage or public result import.',
)


def restore(data_bytes, snapshot_bytes):
    """Reconstruct observed rows, then reuse the actual current-fit preflight."""
    job = app._pcm_native.prepare_inputs(data_bytes, snapshot_bytes)
    saved = json.loads(snapshot_bytes, object_hook=lambda d:
        float(d['nonfinite']) if set(d) == {'nonfinite'} else d)
    cfg = saved['config']
    for spec in [cfg['theta_spec'], *cfg['facet_specs'].values()]:
        spec['anchors'] = np.asarray(spec['anchors'], dtype=float)
        spec['groups'] = np.asarray(spec['groups'], dtype=object)
    cfg['population_model']['X'] = np.asarray(cfg['population_model']['X'], dtype=float)
    data = pd.read_csv(io.BytesIO(data_bytes), keep_default_na=False,
                       dtype={f: str for f in ['Person', *cfg['facet_names']]})
    prep = app.prepare_mfrm_data(data, 'Person', cfg['facet_names'], 'Score', rating_min=0, rating_max=3)
    if prep['levels'] != saved['levels']:
        raise ValueError('Restored level order differs from the saved analysis')
    result = dict(config=cfg, prep=prep, params=saved['params'], opt=SimpleNamespace(
        x=np.asarray(saved['returned_free_coordinates']), fun=saved['optimizer_reported_fun'],
        success=saved['optimizer_reported_success'], message=saved['optimizer_message']))
    prepared = app.prepare_fixed_sd_pcm_check(result)
    prepared['reference_job'] = job
    return prepared


def reference_map(prepared):
    """Exact signed permutation: app coordinates -> retained reference basis."""
    cfg = prepared['config']; c = cfg['step_facet']
    r, t = [f for f in cfg['facet_names'] if f != c]
    native_rows = {c: [0, 1], r: [2], t: [3]}
    transform = np.zeros((8, 8)); position = 0
    for f in cfg['facet_names']:
        for row in native_rows[f]:
            transform[row, position] = -cfg['facet_signs'][f]
            position += 1
    transform[4:, 4:] = np.eye(4)
    return transform[PERM]


def moment_allowances(value, bound, sigma):
    """Conservative per-person tail bounds and propagated quadrature estimates.

    Unit-weight response likelihood is <= 1. The retained reference reports SUMS
    of relative mass bounds/estimates; each sum also bounds its largest component.
    Numerical quadrature error is an estimate, not a rigorous enclosure.
    """
    r = value['numeric_relative_mass_error_sum']
    tail_raw = value['tail_relative_mass_bound_sum']
    scores = np.asarray([value['eap'], value['sd']])
    if (not np.isfinite([r, tail_raw, bound, sigma]).all() or not 0 <= r < 1
            or tail_raw < 0 or not 0 < bound <= 16 or sigma <= 0
            or scores.size == 0 or not np.isfinite(scores).all() or np.any(scores[1] <= 0)):
        raise ValueError('Unusable relative numerical mass error')
    tail = tail_raw / (1-r)
    phi = np.exp(-bound*bound/2)/np.sqrt(2*np.pi)
    m0 = 2*ndtr(-bound)
    m1, m2 = 2*phi, 2*(bound*phi+ndtr(-bound))
    mean = float(np.max(np.abs(value['eap'])))
    second = float(np.max(np.asarray(value['sd'])**2+np.asarray(value['eap'])**2))
    num_mean = r*(sigma+mean)/(1-r)
    num_second = r*(sigma*sigma+second)/(1-r)
    mean_error = num_mean + tail*(sigma*m1/m0+mean+num_mean)
    second_error = num_second + tail*(sigma*sigma*m2/m0+second+num_second)
    variance_error = second_error+mean_error*(2*mean+mean_error)
    min_sd = float(np.min(scores[1]))
    sd_error = variance_error/(min_sd+np.sqrt(max(0., min_sd*min_sd-variance_error)))
    return dict(mean=mean_error, sd=float(sd_error),
                nll=r/(1-r)+tail, numerical_mass=r, tail_mass=tail)


def decompose(finite, continuous, baseline):
    """Keep signed vectors; only summarize magnitudes AFTER decomposition."""
    result = {}
    for key in ('nll', 'eap', 'sd'):
        integration = np.asarray(finite[key])-continuous[key]
        recalibration = np.asarray(continuous[key])-baseline[key]
        total = np.asarray(finite[key])-baseline[key]
        result[key] = dict(integration=integration, recalibration=recalibration, total=total,
            identity_residual=float(np.max(np.abs(total-integration-recalibration))),
            maximum_absolute={name: float(np.max(np.abs(v))) for name, v in
                              [('integration', integration), ('recalibration', recalibration), ('total', total)]})
    return result


def check(prepared, save):
    """Persist each result immediately; the caller owns a new output directory."""
    cfg, idx, sizes, point = (prepared[k] for k in ('config', 'idx', 'sizes', 'coordinates'))
    if len(idx['score_k']) > 1000:
        raise ValueError('Development check supports at most 1,000 response rows')
    initial = app.evaluate_fixed_sd_pcm_check(prepared)
    sigma = prepared['quad']['sd']; transform = reference_map(prepared)
    job = prepared['reference_job']
    if not np.allclose(transform@point, np.asarray(job['xsi'])[PERM], atol=1e-12, rtol=0):
        raise ValueError('Reference coordinate map disagrees with the exported calibration')
    reference = Reference(job)
    limits = PROTOCOL['engineering_limits']

    def continuous(p):
        values = [reference.continuous(transform@p, **s) for s in PROTOCOL['continuous']]
        for v in values:
            v['gradient'] = transform.T@v['gradient']
        refinement = differences(*values)
        allowances = [moment_allowances(v, s['bound'], sigma) for v,s in zip(values, PROTOCOL['continuous'])]
        passed = (max(refinement.values()) < limits['reference_refinement'] and
                  all(max(a['mean'], a['sd']) < limits['reference_moment_allowance'] and
                      max(a['numerical_mass'], a['tail_mass']) < limits['reference_mass_allowance'] for a in allowances))
        return dict(values=values, refinement=refinement, allowances=allowances, reference_checks_pass=bool(passed))

    original_continuous = continuous(point)
    save('original_continuous.json', original_continuous)
    baseline = original_continuous['values'][-1]
    audit_point = point + PROTOCOL['gradient_audit']['offset']
    def audit_vg(p):
        v = reference.continuous(transform@p, **PROTOCOL['continuous'][-1])
        return v['nll'], transform.T@v['gradient']
    audit = audit_joint_gradient(lambda p: audit_vg(p)[0], audit_vg, audit_point,
        relative_step=PROTOCOL['gradient_audit']['relative_step']).to_dict()
    audit.update(coordinates_checked=audit_point,
                 check_pass=audit['maximum_absolute_difference'] < limits['gradient_fd'])
    save('continuous_gradient_audit.json', audit)
    rows = []
    for q in PROTOCOL['orders']:
        quad = app.make_mml_quadrature(cfg, q)
        problem = AppFreeSdProblem(cfg, sizes, idx, tuple(point), sigma, q, float(len(idx['score_k'])),
                                   ((None, None),)*len(point), (.05, 10.))
        vg = lambda p: problem.value_gradient(p, sigma)
        log_sd = float(np.log(sigma))
        def factory(anchor):
            delta = problem.likelihood_difference(np.r_[anchor, log_sd])
            return lambda p: delta(np.r_[p, log_sd])
        def finite(p):
            actual, gradient = vg(p)
            score = app.compute_person_eap(idx, cfg, app.expand_params(p, sizes, cfg), quad)
            result = dict(nll=actual, gradient=gradient, eap=score.Estimate.to_numpy(), sd=score.SD.to_numpy())
            independent = reference.finite(transform@p, quad['nodes']/sigma, quad['weights'])
            independent['gradient'] = transform.T@independent['gradient']
            result['independent_differences'] = differences(result, independent)
            result['independent_check_pass'] = max(result['independent_differences'].values()) < limits['finite_replay']
            return result
        original = finite(point)
        save(f'q{q}_same_point.json', dict(finite=original, continuous=original_continuous,
             decomposition=decompose(original, baseline, baseline)))
        for label, start in [('original', point), ('zero', np.zeros_like(point))]:
            key = f'q{q}_{label}'
            run = polish_fixed_sd(start, vg, factory, options=JointPolishOptions(**PROTOCOL['optimizer']))
            save(key+'_optimizer.json', run)
            row = dict(key=key, run=run, evaluation_completed=False, evaluation_failure=None)
            if run['completed']:
                try:
                    candidate = np.asarray(run['candidate_coordinates'])
                    f = finite(candidate); row['finite'] = f
                    integral = continuous(candidate)
                    row.update(finite=f, continuous=integral, evaluation_completed=True,
                        decomposition=decompose(f, integral['values'][-1], baseline),
                        structural_change=candidate-point)
                except Exception as exc:
                    row['evaluation_failure'] = dict(type=type(exc).__name__, message=str(exc))
            rows.append(row); save(key+'.json', row)
            print(key, 'evaluated:', row['evaluation_completed'], 'gradient target:', run.get('gradient_tolerance_met'), flush=True)
    cross_start = {str(q): float(np.max(np.abs(np.asarray(rows[i]['run']['candidate_coordinates'])-
                                    rows[i+1]['run']['candidate_coordinates'])))
                   for i,q in zip(range(0,6,2), PROTOCOL['orders'])
                   if rows[i]['run']['completed'] and rows[i+1]['run']['completed']}
    return dict(initial=initial, coordinate_map=transform, person_order=job['person_labels'],
        original_continuous=original_continuous, continuous_gradient_audit=audit,
        trajectories=rows, cross_start_coordinates=cross_start,
        cross_start_checks={q: d < limits['cross_start_coordinates'] for q,d in cross_start.items()},
        expected_trajectories=6, evaluated_trajectories=sum(r['evaluation_completed'] for r in rows),
        inference_ready=False, continuous_refitted=False, original_result_replaced=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); args.output.mkdir(parents=True, exist_ok=False)
    sources = ['validation/mml_fixed_sd_accuracy_check.py', 'validation/mml_fixed_sd_pcm_polish.py',
        'validation/mml_free_sd_stationarity_adapter.py', 'validation/mml_pcm_continuous_refit.py',
        'streamlit_app.py', 'mfrm_app/mml_stationarity.py', 'mfrm_app/mml_engine_v2.py',
        'mfrm_app/legacy_compat/pcm_native_runner.py']
    inputs = {name:(args.input/name).read_bytes() for name in ('data.csv','analysis_snapshot.json')}
    for name, content in inputs.items(): (args.output/name).write_bytes(content)
    with zipfile.ZipFile(args.output/'sources.zip', 'x', zipfile.ZIP_DEFLATED) as z:
        for name in sources: z.write(ROOT/name, name)
    dump(args.output/'protocol.json', dict(PROTOCOL, source_sha256={n:sha(ROOT/n) for n in sources},
         input_sha256={n:sha(args.output/n) for n in inputs}))
    started = time.perf_counter()
    try:
        prepared = restore(inputs['data.csv'], inputs['analysis_snapshot.json'])
        dump(args.output/'runtime.json', prepared['implementation'])
        result = check(prepared, lambda name, obj: dump(args.output/name, obj))
        result['elapsed_seconds'] = time.perf_counter()-started
        dump(args.output/'summary.json', result)
        dump(args.output/'manifest.json', {p.name:sha(p) for p in sorted(args.output.iterdir()) if p.is_file()})
    except Exception as exc:
        dump(args.output/'failure.json', dict(type=type(exc).__name__, message=str(exc)))
        raise
    print('Saved development comparison:', args.output)


if __name__ == '__main__':
    main()
