#!/usr/bin/env python3
"""Localize the retained fixed-SD reference error without changing integration.

Offline, single-threaded instrumentation only. No optimization or result import.
"""
import argparse
import json
from pathlib import Path
import sys
import time
from unittest.mock import patch
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from validation import mml_fixed_sd_accuracy_check as prior

PROTOCOL = dict(
    classification='OBSERVED_DEVELOPMENT_ONLY', inference_ready=False,
    cases=['positive_rater', 'one_missing'], points='Original and all six retained endpoints in each case',
    integration='Exactly the two prior settings, same integrand, scaling, breakpoints and quad_vec return value; capture its raw integrals and shared max-norm error estimate.',
    comparison=['Prior sum of relative errors and global moment extremes',
                'Maximum relative error and global moment extremes',
                'Each relative error and its own moments'],
    tails='Retain the prior all-person tail-mass sum as an upper bound for EACH person; tail localization is not required for this comparison.',
    limits=dict(prior.PROTOCOL['engineering_limits'], saved_value_replay=1e-12),
    interpretation='A shared estimated absolute error divided by each scaled mass is not a measured individual error. Quadrature estimates are not rigorous enclosures. Total NLL keeps the SUM-based allowance. Preserve all prior flags.',
    failure_policy='Keep every point and both settings, including failures. No retuning or best-setting selection.',
)


def inspect_reference(reference, point, settings):
    """Observe the unchanged integral and reuse the existing moment propagation."""
    module = sys.modules[prior.Reference.__module__]
    original = module.quad_vec
    captured = []

    def record(*args, **kwargs):
        result = original(*args, **kwargs)
        captured.append((np.array(result[0], copy=True), float(result[1])))
        return result

    # ponytail: temporary module patch is only for this sequential offline audit;
    # use an explicit diagnostic return if later integrating with concurrent UI work.
    with patch.object(module, 'quad_vec', record):
        value = reference.continuous(point, **settings)
    if len(captured) != 1:
        raise ValueError('Expected exactly one vector integration call')
    integrals, error = captured[0]
    if (integrals.shape != (len(reference.y), 11) or not np.isfinite(integrals).all()
            or np.any(integrals[:, 0] <= 0) or not np.isfinite(error) or error < 0):
        raise ValueError('Invalid captured integration diagnostics')
    relative = error/integrals[:, 0]
    if float(relative.sum()) != value['numeric_relative_mass_error_sum']:
        raise ValueError('Captured relative errors do not reconstruct the retained sum')
    allowance = lambda v: prior.moment_allowances(v, settings['bound'], reference.sigma)
    pooled = allowance(value)
    max_relative = allowance(dict(value, numeric_relative_mass_error_sum=float(relative.max())))
    personal = [allowance(dict(value, numeric_relative_mass_error_sum=float(r),
                             eap=[value['eap'][i]], sd=[value['sd'][i]]))
                for i, r in enumerate(relative)]
    return dict(value=value, scaled_integrals=integrals, shared_absolute_error_estimate=error,
        relative_error_estimate=relative, pooled=pooled, maximum_relative_error=max_relative,
        per_person=personal, total_nll_allowance=pooled['nll'])


def audit_case(folder, save):
    saved = json.loads((folder/'summary.json').read_text())
    prepared = prior.restore((folder/'data.csv').read_bytes(), (folder/'analysis_snapshot.json').read_bytes())
    transform = prior.reference_map(prepared)
    reference = prior.Reference(prepared['reference_job'])
    points = [('original', prepared['coordinates'], saved['original_continuous'])]
    for row in saved['trajectories']:
        points.append((row['key'], row['run']['candidate_coordinates'], row.get('continuous')))
    records = []
    limits = PROTOCOL['limits']
    for label, point, old in points:
        row = dict(point=label, coordinates=point, prior_reference_checks_pass=None if old is None else old['reference_checks_pass'],
                   completed=False, failure=None, settings=[], localized_diagnostics_meet_thresholds=False)
        try:
            if point is None or old is None:
                raise ValueError('Retained endpoint or continuous evaluation unavailable')
            for i, settings in enumerate(prior.PROTOCOL['continuous']):
                result = inspect_reference(reference, transform@np.asarray(point), settings)
                result['value']['gradient'] = transform.T@result['value']['gradient']
                difference = prior.differences(result['value'], old['values'][i], keys=(
                    'nll', 'gradient', 'eap', 'sd', 'numeric_relative_mass_error_sum', 'tail_relative_mass_bound_sum'))
                result.update(settings=settings, saved_value_difference=difference)
                row['settings'].append(result)
                if (max(difference.values()) > limits['saved_value_replay']
                        or result['value']['quadrature_neval'] != old['values'][i]['quadrature_neval']):
                    raise ValueError('Reference calculation changed from the retained evaluation')
            refinement = prior.differences(*(r['value'] for r in row['settings']))
            meets = (max(refinement.values()) < limits['reference_refinement'] and
                all(max(a['mean'], a['sd']) < limits['reference_moment_allowance']
                    for r in row['settings'] for a in r['per_person']) and
                all(max(r['pooled']['numerical_mass'], r['pooled']['tail_mass']) < limits['reference_mass_allowance']
                    for r in row['settings']))
            row.update(completed=True, refinement=refinement, localized_diagnostics_meet_thresholds=bool(meets))
        except Exception as exc:
            row['failure'] = dict(type=type(exc).__name__, message=str(exc))
        records.append(row); save(label+'.json', row)
        print(folder.name, label, 'completed:', row['completed'], 'localized:', row['localized_diagnostics_meet_thresholds'], flush=True)
    return dict(person_order=prepared['reference_job']['person_labels'], expected_points=7, records=records,
                implementation=prepared['implementation'],
                prior_runtime_matches=prepared['implementation'] == json.loads((folder/'runtime.json').read_text()),
                inference_ready=False, prior_flags_replaced=False, refitted=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=False)
    qa = json.loads((args.previous/'qa.json').read_text())
    sources = set(); inputs = {}
    for case in PROTOCOL['cases']:
        folder = args.previous/(case+'_v2')
        if prior.sha(folder/'summary.json') != qa['cases'][case]['summary_sha256']:
            raise ValueError('Prior summary hash mismatch')
        if prior.sha(folder/'manifest.json') != qa['cases'][case]['manifest_sha256'][folder.name]:
            raise ValueError('Prior manifest hash mismatch')
        for name, digest in json.loads((folder/'manifest.json').read_text()).items():
            if prior.sha(folder/name) != digest:
                raise ValueError('Prior artifact hash mismatch: '+name)
        protocol = json.loads((folder/'protocol.json').read_text())
        for name, digest in protocol['source_sha256'].items():
            if prior.sha(ROOT/name) != digest:
                raise ValueError('Retained numerical source changed: '+name)
            sources.add(name)
        if protocol['continuous'] != prior.PROTOCOL['continuous']:
            raise ValueError('Prior integration settings changed')
        dest = args.output/'inputs'/case; dest.mkdir(parents=True)
        for name in ('summary.json', 'data.csv', 'analysis_snapshot.json', 'protocol.json', 'runtime.json'):
            (dest/name).write_bytes((folder/name).read_bytes())
        inputs[case] = {p.name:prior.sha(p) for p in dest.iterdir()}
    sources.update(['validation/mml_fixed_sd_reference_error_audit.py', 'tests/test_mml_fixed_sd_accuracy_check.py'])
    with zipfile.ZipFile(args.output/'sources.zip', 'x', zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(sources): archive.write(ROOT/name, name)
    prior.dump(args.output/'protocol.json', dict(PROTOCOL, input_sha256=inputs,
        source_sha256={n:prior.sha(ROOT/n) for n in sorted(sources)}, previous_qa_sha256=prior.sha(args.previous/'qa.json')))
    started = time.perf_counter(); results = {}
    for case in PROTOCOL['cases']:
        out = args.output/case; out.mkdir()
        try:
            results[case] = audit_case(args.output/'inputs'/case, lambda n,v: prior.dump(out/n,v))
        except Exception as exc:
            results[case] = dict(failure=dict(type=type(exc).__name__, message=str(exc)), inference_ready=False)
        prior.dump(out/'summary.json', results[case])
    prior.dump(args.output/'summary.json', dict(cases=results, elapsed_seconds=time.perf_counter()-started,
               inference_ready=False, prior_flags_replaced=False, refitted=False))
    prior.dump(args.output/'manifest.json', {str(p.relative_to(args.output)):prior.sha(p)
               for p in sorted(args.output.rglob('*')) if p.is_file()})
    print('Saved reference error audit:', args.output)


if __name__ == '__main__':
    main()
