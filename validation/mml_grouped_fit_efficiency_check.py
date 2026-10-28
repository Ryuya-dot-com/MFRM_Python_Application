"""Paired complete-fit development benchmark on retained development inputs.

No population-study inputs, changed solver settings, or selected retries.
Timing includes grouping, the existence audit, both mandatory solver stages,
joint observed information and covariance. It excludes problem preparation,
input/output, and application reporting. Shared-machine timing is diagnostic.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import shutil
import statistics
import time
import zipfile

import numpy as np
import pandas as pd
import scipy

from validation import mml_grouped_fit as solver
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data
from validation.mml_gpcm_joint_check import prepare_problem
from validation.mml_population_study import json_bytes


LIMITS = dict(point_abs=2e-5, nll_per_person=1e-8, covariance_scaled=1e-4,
              eap_abs=1e-5, posterior_sd_abs=1e-5)


def digest(path):
    raw = path.read_bytes()
    if len(raw) != path.stat().st_size:
        raise OSError(f'Incomplete local read: {path}')
    return hashlib.sha256(raw).hexdigest()


def timed_fit(engine, data, problem, q, *, fit_module=solver):
    phases = {key: dict(calls=0, seconds=0.) for key in
              ['construct', 'existence_audit', 'value_gradient', 'joint_information']}

    def record(key, function, *args, **kwargs):
        start = time.perf_counter()
        try:
            return function(*args, **kwargs)
        finally:
            phases[key]['seconds'] += time.perf_counter()-start
            phases[key]['calls'] += 1

    class Measured(engine):
        def __init__(self, *args, **kwargs):
            record('construct', super().__init__, *args, **kwargs)

        def structural_recession_audit(self, *args, **kwargs):
            return record('existence_audit', super().structural_recession_audit, *args, **kwargs)

        def evaluate(self, *args, **kwargs):
            key = 'joint_information' if kwargs.get('information', False) else 'value_gradient'
            return record(key, super().evaluate, *args, **kwargs)

    original = fit_module.GroupedMmlLikelihood
    fit_module.GroupedMmlLikelihood = Measured
    start = time.perf_counter()
    try:
        result = fit_module.fit_grouped_joint(data, problem, q=q)
    except Exception as error:
        result = dict(status='exception', numerical_pass=False,
                      error_type=type(error).__name__, error=str(error))
    finally:
        elapsed = time.perf_counter()-start
        fit_module.GroupedMmlLikelihood = original
    # Retain sufficient evidence to reconstruct the large posterior grid.
    if 'evaluated' in result:
        result['evaluated'].pop('posterior', None)
    return dict(seconds=elapsed, phases=phases, result=result)


def compare(old, new, persons):
    a, b = old['result'], new['result']
    checks = dict(both_numerical_pass=bool(a['numerical_pass'] and b['numerical_pass']),
                  same_status=a['status'] == b['status'], same_solver_checks=a.get('checks') == b.get('checks'))
    errors = {}
    if a.get('point') is not None and b.get('point') is not None:
        errors['point_abs'] = float(np.max(abs(a['point']-b['point'])))
        errors['nll_per_person'] = float(abs(a['evaluated']['nll']-b['evaluated']['nll'])/persons)
        for key, metric in [('eap', 'eap_abs'), ('sd', 'posterior_sd_abs')]:
            errors[metric] = float(np.max(abs(a['evaluated'][key]-b['evaluated'][key])))
    if a.get('covariance') is not None and b.get('covariance') is not None:
        scale = np.sqrt(np.outer(np.diag(a['covariance']), np.diag(a['covariance'])))
        errors['covariance_scaled'] = float(np.max(abs(a['covariance']-b['covariance'])/scale))
    checks.update({key: bool(key in errors and np.isfinite(errors[key]) and errors[key] <= limit)
                   for key, limit in LIMITS.items()})
    return dict(checks=checks, errors=errors, all_checks_pass=all(checks.values()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('parent', type=Path)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('retained_grouped_baseline', args.baseline)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    engines = [module.GroupedMmlLikelihood, solver.GroupedMmlLikelihood]
    cases = [f'{m}_{a}_N192_K4_Q641' for m in ['RSM', 'PCM', 'GPCM']
             for a in ['regular', 'sparse', 'unbalanced']]
    sources = json.loads((args.parent/'source_manifest.json').read_text())
    sources['validation/mml_grouped_fit_efficiency_check.py'] = digest(Path(__file__))
    sources['validation/mml_grouped_likelihood_baseline.py'] = digest(args.baseline)
    with zipfile.ZipFile(out/'sources.zip', 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for name, expected in sources.items():
            path = root/name
            if digest(path) != expected:
                raise ValueError(f'Source differs from retained development snapshot: {name}')
            archive.write(path, name)
    (out/'source_manifest.json').write_bytes(json_bytes(sources))
    inputs = {}
    for name in cases:
        metadata = json.loads((args.parent/(name+'.json')).read_text())
        if digest(args.parent/(name+'.csv')) != metadata['input_sha256']:
            raise ValueError(f'Retained input mismatch: {name}')
        for suffix in ['.csv', '.json']:
            source = args.parent/(name+suffix)
            inputs[name+suffix] = digest(source)
            shutil.copyfile(source, out/source.name)
    plan = dict(classification='DEVELOPMENT_COMPLETE_FIT_EFFICIENCY_ONLY',
        scientific_inference_ready=False, utc=datetime.now(timezone.utc).isoformat(),
        cases=cases, quadrature_order=641, repeats=3, comparison_limits=LIMITS,
        solver_options=solver.OPTIONS, solver_limits=solver.LIMITS, inputs_sha256=inputs,
        source_manifest_sha256=digest(out/'source_manifest.json'),
        parent_manifest_sha256=digest(args.parent/'source_manifest.json'),
        protocol='Serial paired alternating order; every fit starts cold and has the mandatory polish; every failure retained; no retries',
        scope=__doc__, omitted_result_fields=['evaluated.posterior'],
        runtime=dict(python=platform.python_version(), numpy=np.__version__,
                     scipy=scipy.__version__, pandas=pd.__version__, platform=platform.platform()))
    (out/'plan.json').write_bytes(json_bytes(plan))
    records = []
    for name in cases:
        data = pd.read_csv(out/(name+'.csv')); model = name.split('_')[0]
        if model == 'GPCM':
            _, problem = prepare_problem(data, 641)
        else:
            _, problem = prepare_app_free_sd_problem_from_data(data, person_col='Person',
                facet_cols=['Rater','Task','Criterion'], score_col='Score', rating_min=0, rating_max=3,
                model=model, noncenter_facet='Criterion', step_facet=None if model=='RSM' else 'Criterion',
                quadrature_points=641)
        pairs = []
        for repeat in range(plan['repeats']):
            pair = {}
            for j in ([0,1] if repeat%2==0 else [1,0]):
                label = ['baseline','candidate'][j]
                pair[label] = timed_fit(engines[j], data, problem, 641)
                (out/f'{name}_{repeat}_{label}.json').write_bytes(json_bytes(pair[label]))
                print(name, repeat, label, pair[label]['result']['status'], round(pair[label]['seconds'],3), flush=True)
            pair['comparison'] = compare(pair['baseline'], pair['candidate'], data['Person'].nunique())
            pairs.append(pair)
        medians = {label: statistics.median(p[label]['seconds'] for p in pairs) for label in ['baseline','candidate']}
        record = dict(id=name, median_seconds=medians, median_time_ratio=medians['baseline']/medians['candidate'],
                      comparisons=[p['comparison'] for p in pairs],
                      all_checks_pass=all(p['comparison']['all_checks_pass'] for p in pairs))
        records.append(record)
        (out/(name+'_comparison.json')).write_bytes(json_bytes(record))
    summary = dict(classification=plan['classification'], scientific_inference_ready=False,
                   records=records, all_checks_pass=all(r['all_checks_pass'] for r in records),
                   completed_utc=datetime.now(timezone.utc).isoformat())
    (out/'summary.json').write_bytes(json_bytes(summary))
    return 0 if summary['all_checks_pass'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
