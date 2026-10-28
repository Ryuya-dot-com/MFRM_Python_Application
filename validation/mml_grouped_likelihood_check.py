"""Retained-point grouped-likelihood/R checks, with separate diagnostic timing."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import time
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from validation.mml_grouped_likelihood import GroupedMmlLikelihood
from validation.mml_joint_information_check import save
from validation.mml_rasch_quadrature_review import verify
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data
from validation.mml_gpcm_joint_check import prepare_problem


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['rasch', 'rasch_reference', 'gpcm', 'output']:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    for p in [args.rasch, args.rasch_reference, args.gpcm]: verify(p)
    if (args.rasch/'reference_input.json').read_bytes() != (args.rasch_reference/'parent_reference_input.json').read_bytes():
        raise ValueError('RSM/PCM R reference input differs')
    out = args.output; out.mkdir(parents=True, exist_ok=False)
    limits = dict(native_nll=1e-8, native_gradient=1e-9, native_information=1e-4,
                  r_nll=1e-8, r_gradient=2e-5, r_information=1e-4)
    save(out/'protocol.json', dict(classification='DEVELOPMENT_ONLY', scientific_inference_ready=False,
        planned=18, no_refitting=True, limits=limits, timing_repeats=3,
        timing_scope='Shared Mac; concurrent R reference job. Diagnostic medians, not a general production speed claim.',
        parent_manifests={str(p): hashlib.sha256((p/'manifest.json').read_bytes()).hexdigest()
                          for p in [args.rasch, args.rasch_reference, args.gpcm]}))
    sources = [ROOT/'streamlit_app.py', *sorted((ROOT/'mfrm_app').glob('*.py')),
        *[ROOT/'validation'/name for name in ['mml_grouped_likelihood.py', 'mml_grouped_likelihood_check.py',
            'mml_rasch_integral.py', 'mml_gpcm_integral.py', 'mml_free_sd_stationarity_adapter.py',
            'mml_free_sd_quadrature_adapter.py', 'mml_joint_information_check.py', 'mml_rasch_quadrature_review.py', 'mml_gpcm_joint_check.py']]]
    with zipfile.ZipFile(out/'sources.zip', 'w', zipfile.ZIP_DEFLATED) as z:
        for p in sources: z.write(p, str(p.relative_to(ROOT)))
    save(out/'source_manifest.json', {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    save(out/'runtime.json', dict(python=sys.version, packages={d.metadata['Name']: d.version for d in importlib.metadata.distributions()}))
    records, failures = [], []
    for parent, oracle in [(args.rasch, args.rasch_reference), (args.gpcm, args.gpcm)]:
        cases = json.loads((parent/'reference_input.json').read_text())['cases']
        r = {c['id']: c for c in json.loads((oracle/'R_results.json').read_text())['cases']}
        for case in cases:
            try:
                cfg, data, point = case['config'], pd.DataFrame(case['data']), np.array(case['point'])
                q = len(case['nodes'])
                if cfg['model'] == 'GPCM':
                    _, native = prepare_problem(data, q)
                else:
                    _, native = prepare_app_free_sd_problem_from_data(data, person_col='Person', facet_cols=cfg['facet_names'],
                        score_col='Score', rating_min=cfg['rating_min'], rating_max=cfg['rating_max'], model=cfg['model'],
                        noncenter_facet=cfg['noncenter_facet'], step_facet=cfg['step_facet'],
                        positive_facets=tuple(f for f in cfg['facet_names'] if cfg['facet_signs'][f] == 1), quadrature_points=q)
                grouped = GroupedMmlLikelihood(data, cfg)
                result = grouped.evaluate(point, case['nodes'], case['weights'], information=True)
                saved = json.loads((parent/(case['id']+'.json')).read_text())
                errors = {prefix+'_'+key: float(np.max(abs(np.array(result[key])-np.array(other[key]))))
                          for prefix, other in [('native', saved), ('r', r[case['id']])]
                          for key in ['nll', 'gradient', 'information']}
                checks = {k: v < limits[k] for k, v in errors.items()}
                times = {}
                for label, calculate in [
                    ('native_value_gradient', lambda: native.joint_value_gradient(point)),
                    ('grouped_value_gradient', lambda: grouped.evaluate(point, case['nodes'], case['weights'])),
                    ('grouped_information', lambda: grouped.evaluate(point, case['nodes'], case['weights'], information=True))]:
                    values = []
                    for _ in range(3):
                        start = time.perf_counter(); calculate(); values.append(time.perf_counter()-start)
                    times[label] = dict(seconds=values, median=float(np.median(values)))
                record = dict(id=case['id'], point=point, grouped=result, errors=errors, checks=checks, timings=times,
                              all_checks_pass=all(checks.values()))
                save(out/(case['id']+'.json'), record)
                records.append({k: record[k] for k in ['id', 'errors', 'checks', 'timings', 'all_checks_pass']})
                print(case['id'], 'checks', all(checks.values()), flush=True)
            except Exception as exc:
                failures.append(dict(id=case['id'], type=type(exc).__name__, message=str(exc)))
                save(out/'failures.json', failures)
                print(case['id'], 'FAILED', repr(exc), flush=True)
    summary = dict(classification='DEVELOPMENT_ONLY', scientific_inference_ready=False, planned=18,
        completed=len(records), failures=failures, comparisons=records,
        all_checks_pass=len(records) == 18 and not failures and all(r['all_checks_pass'] for r in records))
    save(out/'summary.json', summary)
    save(out/'manifest.json', {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()})
    return 0 if summary['all_checks_pass'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
