"""Continuous GPCM joint-information accuracy at every retained study point.

No refitting or sampling claims. Frozen same-point finite-GH and independent R
continuous-moment artifacts are required; failed points remain in the report.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from validation.mml_gpcm_integral import GpcmIntegral
from validation.mml_joint_information_check import save
from validation.mml_rasch_quadrature_review import verify


SETTINGS = [dict(bound=12., rel_tol=1e-9, abs_tol=1e-11),
            dict(bound=16., rel_tol=1e-11, abs_tol=1e-13)]
LIMITS = dict(reference_nll=1e-8, reference_gradient=1e-8, reference_information=1e-7,
    reference_moments=1e-9, relative_mass_error=1e-7, tail_mass=1e-10,
    tail_gradient=1e-8, tail_information=1e-7, numeric_information_allowance=1e-5,
    r_nll=1e-8, r_moments=1e-9, r_sigma_score=1e-8,
    gh_nll_per_person=1e-6, gh_score_supnorm=1e-4, gh_moments=1e-4,
    gh_expanded_se_relative=.005)


def difference(a, b, key):
    return float(np.max(abs(np.asarray(a[key])-np.asarray(b[key]))))


def evaluate(case, point_record, r):
    reference = GpcmIntegral(pd.DataFrame(case['data']), case['config'])
    results = [reference.continuous(case['point'], **s, information=True) for s in SETTINGS]
    loose, tight = results
    refinement = {k: difference(loose, tight, k) for k in ['nll', 'gradient', 'information', 'eap', 'sd']}
    integration = {k: difference(point_record, tight, k) for k in ['nll', 'gradient', 'information', 'eap', 'sd']}
    r_error = {k: difference(tight, r, k) for k in ['nll', 'eap', 'sd', 'density_log_sigma_score']}
    jacobian = np.array(point_record['expanded_jacobian'])
    eigenvalues, covariance, expanded_se = {}, {}, {}
    for label, matrix in [('finite', np.array(point_record['information'])), ('continuous', tight['information'])]:
        eigenvalues[label] = np.linalg.eigvalsh(matrix)
        if eigenvalues[label].min() > 0:
            covariance[label] = np.linalg.inv(matrix)
            variance = np.einsum('ij,jk,ik->i', jacobian, covariance[label], jacobian)
            if np.any(variance < 0):
                raise FloatingPointError('Negative expanded variance; no flooring applied')
            expanded_se[label] = np.sqrt(variance)
    se_relative, worst = None, None
    if len(expanded_se) == 2:
        valid = expanded_se['continuous'] > 0
        if np.any((~valid) & (expanded_se['finite'] != 0)):
            raise ValueError('Fixed-coordinate variance changed')
        error = abs(expanded_se['finite'][valid]/expanded_se['continuous'][valid]-1)
        se_relative = float(error.max())
        worst = np.array(point_record['expanded_labels'])[valid][error.argmax()]
    checks = dict(
        reference_nll=refinement['nll'] < LIMITS['reference_nll'],
        reference_gradient=refinement['gradient'] < LIMITS['reference_gradient'],
        reference_information=refinement['information'] < LIMITS['reference_information'],
        reference_moments=max(refinement['eap'], refinement['sd']) < LIMITS['reference_moments'],
        reference_mass=tight['numeric_relative_mass_error_sum'] < LIMITS['relative_mass_error'],
        reference_tail_mass=tight['tail_relative_mass_bound_sum'] < LIMITS['tail_mass'],
        reference_tail_gradient=tight['tail_gradient_abs_bound'] < LIMITS['tail_gradient'],
        reference_tail_information=tight['tail_information_abs_bound'] < LIMITS['tail_information'],
        reference_information_allowance=tight['numeric_information_error_allowance'] < LIMITS['numeric_information_allowance'],
        r_nll=r_error['nll'] < LIMITS['r_nll'],
        r_moments=max(r_error['eap'], r_error['sd']) < LIMITS['r_moments'],
        r_sigma_score=r_error['density_log_sigma_score'] < LIMITS['r_sigma_score'],
        gh_nll=integration['nll']/reference.n < LIMITS['gh_nll_per_person'],
        gh_gradient=integration['gradient'] < LIMITS['gh_score_supnorm'],
        gh_moments=max(integration['eap'], integration['sd']) < LIMITS['gh_moments'],
        gh_expanded_se=se_relative is not None and se_relative < LIMITS['gh_expanded_se_relative'])
    return dict(id=case['id'], point=case['point'], continuous=results, refinement=refinement,
        integration=integration, R_error=r_error, eigenvalues=eigenvalues, covariance=covariance,
        expanded_se=expanded_se, expanded_se_relative_difference=se_relative,
        worst_expanded_se_coordinate=worst, checks=checks, all_checks_pass=all(checks.values()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('parent', type=Path)
    parser.add_argument('continuous_reference', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    verify(args.parent); verify(args.continuous_reference)
    parent_digest = hashlib.sha256((args.parent/'manifest.json').read_bytes()).hexdigest()
    reference_protocol = json.loads((args.continuous_reference/'protocol.json').read_text())
    if reference_protocol['parent_manifest_sha256'] != parent_digest:
        raise ValueError('R reference belongs to a different parent')
    parent_summary = json.loads((args.parent/'summary.json').read_text())
    cases = json.loads((args.parent/'reference_input.json').read_text())['cases']
    r_results = {c['id']: c['continuous'][-1] for c in json.loads((args.continuous_reference/'R_results.json').read_text())['cases']}
    out = args.output; out.mkdir(parents=True, exist_ok=False)
    save(out/'protocol.json', dict(classification='DEVELOPMENT_ONLY', scientific_inference_ready=False,
        parent_manifest_sha256=parent_digest,
        R_manifest_sha256=hashlib.sha256((args.continuous_reference/'manifest.json').read_bytes()).hexdigest(),
        parent_all_checks_pass=parent_summary['all_checks_pass'], planned_cases=parent_summary['planned'],
        settings=SETTINGS, limits=LIMITS,
        scope='Observed engineering points; continuous full joint information, including nonlinear slopes; no coverage claims',
        error_contract='Gaussian polynomial tail enclosures; adaptive numerical error is estimated, not rigorous'))
    sources = [Path(__file__), ROOT/'validation/mml_gpcm_integral.py', ROOT/'validation/mml_rasch_integral.py',
        ROOT/'validation/mml_rasch_quadrature_review.py', ROOT/'validation/mml_joint_information_check.py']
    with zipfile.ZipFile(out/'sources.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sources:
            archive.write(path, str(path.relative_to(ROOT)))
    save(out/'source_manifest.json', {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    save(out/'runtime.json', dict(python=sys.version, packages={d.metadata['Name']: d.version for d in importlib.metadata.distributions()}))
    records, failures = [], []
    for case in cases:
        try:
            point = json.loads((args.parent/(case['id']+'.json')).read_text())
            record = evaluate(case, point, r_results[case['id']])
            save(out/(case['id']+'.json'), record)
            records.append({k: record[k] for k in ['id', 'refinement', 'integration', 'R_error',
                'expanded_se_relative_difference', 'worst_expanded_se_coordinate', 'checks', 'all_checks_pass']})
            print(case['id'], 'complete; failed:', [k for k, v in record['checks'].items() if not v], flush=True)
        except Exception as exc:
            failures.append(dict(id=case['id'], type=type(exc).__name__, message=str(exc)))
            save(out/'failures.json', failures)
            print(case['id'], 'FAILED', repr(exc), flush=True)
    summary = dict(classification='DEVELOPMENT_ONLY', scientific_inference_ready=False,
        planned=parent_summary['planned'], completed=len(records), failures=failures, comparisons=records,
        all_checks_pass=len(records) == parent_summary['planned'] and not failures and all(r['all_checks_pass'] for r in records))
    save(out/'summary.json', summary)
    save(out/'manifest.json', {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file()})
    return 0 if summary['all_checks_pass'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
