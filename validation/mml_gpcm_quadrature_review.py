"""Consolidate every planned high-order GPCM development point and failure.

The highest-order summary does not erase lower-order failures or authorize
scientific inference. Independent starting points are never selected away.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from validation.mml_rasch_quadrature_review import verify


def review(fit, moments, information):
    planned = [f'GPCM_{design}_Q{q}_{start}' for design in ['regular', 'sparse', 'unbalanced']
               for q, start in [(121, 'from61'), (241, 'from61'), (321, 'from61'), (321, 'from31')]]
    fit_rows = {r['id']: r for r in fit['comparisons']}
    moment_rows = {r['id']: r for r in moments['comparisons']}
    info_rows = {r['id']: r for r in information['comparisons']}
    records = []
    for key in planned:
        a, b, c = fit_rows.get(key), moment_rows.get(key), info_rows.get(key)
        checks = dict(finite_reference=a is not None and a['all_checks_pass'],
            continuous_moments=b is not None and b['all_checks_pass'],
            continuous_information=c is not None and c['all_checks_pass'])
        q = int(key.split('_Q')[1].split('_')[0])
        if q == 321:
            distance = fit.get('cross_start_coordinate_difference', {}).get(key.split('_')[1])
            checks['cross_start'] = distance is not None and distance < 1e-6
        records.append(dict(id=key, q=q, checks=checks, all_checks_pass=all(checks.values()),
            expanded_se_relative_difference=c['expanded_se_relative_difference'] if c is not None else None))
    exact = all(len(s['comparisons']) == len(planned) and set(rows) == set(planned)
                for s, rows in [(fit, fit_rows), (moments, moment_rows), (information, info_rows)])
    complete = (exact and fit['completed'] == information['completed'] == len(planned)
                and fit['planned'] == information['planned'] == len(planned)
                and not fit['failures'] and not information['failures'] and moments['complete']
                and moments['closed_form']['passed'] and fit['R_returncode'] == 0)
    highest = [r for r in records if r['q'] == 321]
    return dict(classification='POST_RESULT_DEVELOPMENT_REVIEW', scientific_inference_ready=False,
        scope='Twelve retained high-order points; no recovery, coverage or fresh-cohort qualification',
        complete=complete, planned=len(planned), comparisons=records,
        highest_order_all_development_checks_pass=complete and len(highest) == 6 and all(r['all_checks_pass'] for r in highest),
        all_orders_all_checks_pass=complete and all(r['all_checks_pass'] for r in records))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['fit', 'moments', 'information', 'output']:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    sources = {k: getattr(args, k) for k in ['fit', 'moments', 'information']}
    for directory in sources.values():
        verify(directory)
    digests = {k: hashlib.sha256((d/'manifest.json').read_bytes()).hexdigest() for k, d in sources.items()}
    moment_protocol = json.loads((args.moments/'protocol.json').read_text())
    info_protocol = json.loads((args.information/'protocol.json').read_text())
    if (moment_protocol['parent_manifest_sha256'] != digests['fit']
        or info_protocol['parent_manifest_sha256'] != digests['fit']
        or info_protocol['R_manifest_sha256'] != digests['moments']):
        raise ValueError('The three evidence sets have different parent identities')
    contract = json.loads((args.fit/'contract.json').read_text())
    if contract.get('cross_start_tolerance') != 1e-6 or contract['planned_cases'] != 12:
        raise ValueError('Unexpected high-order study contract')
    result = review(*(json.loads((sources[k]/'summary.json').read_text()) for k in ['fit', 'moments', 'information']))
    result['parent_manifest_sha256'] = digests
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/'source.py').write_bytes(Path(__file__).read_bytes())
    (args.output/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    manifest = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(args.output.iterdir())}
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'comparisons'}, indent=2), flush=True)
    return 0 if result['highest_order_all_development_checks_pass'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
