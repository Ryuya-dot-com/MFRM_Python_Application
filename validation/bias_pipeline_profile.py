"""Profile a separate development cohort through the sealed G3 implementation.

Includes actual source/runtime checks, durable writes and selected independent
R references. Never writes to a formal study or confers statistical qualification.
Profiler overhead and concurrent formal workers preclude throughput forecasts.
"""
import argparse
from collections import defaultdict
from contextlib import contextmanager
import cProfile
from datetime import datetime, timezone
from functools import wraps
import hashlib
import json
from pathlib import Path
import pstats
import shutil
import sys
import time
import zipfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('study', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    root, study, out = (p.resolve() for p in (args.source, args.study, args.output))
    out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(root))
    from validation import bias_population_controller as controller
    from validation import bias_population_study as executor
    from validation import bias_population_references as references
    from validation.mml_grouped_likelihood import GroupedMmlLikelihood
    from validation.mml_population_study import AttemptLedger, json_bytes
    from validation.bias_population_protocol import planned_bias_runs

    write = lambda path, value: path.write_bytes(json_bytes(value))
    clock = defaultdict(lambda: dict(calls=0, inclusive=0., exclusive=0.))
    stack = []

    @contextmanager
    def phase(name):
        frame = [time.perf_counter(), 0.]
        stack.append(frame)
        try:
            yield
        finally:
            elapsed = time.perf_counter() - frame[0]
            stack.pop()
            entry = clock[name]
            entry['calls'] += 1
            entry['inclusive'] += elapsed
            entry['exclusive'] += elapsed - frame[1]
            if stack:
                stack[-1][1] += elapsed

    def wrap(owner, name, label):
        original = getattr(owner, name)
        @wraps(original)
        def timed(*a, **kw):
            with phase(label):
                return original(*a, **kw)
        setattr(owner, name, timed)

    with phase('startup_seal_capacity_and_R_runtime'):
        formal_plan, formal_identity = controller.verify_seal(study)
        capacity = controller.capacity_check(study, formal_plan, formal_identity)
        verify = controller.environment_verifier(study, formal_plan)
    selected = [s for s in formal_plan['scenarios']
                if (s['design']['persons'] == 192 and len(s['design']['steps'][0]) == 5)
                or (s['design']['persons'] == 48 and s['design']['assignment'] == 'regular' and s['null'])]
    assert len(selected) == 21
    plan = dict(formal_plan, classification='DEVELOPMENT_ONLY', qualification_execution_allowed=False,
                study_id='DEVELOPMENT_ONLY_g3_pipeline_profile_20261006', master_seed=2**117+2026100604,
                scenarios=[dict(s, replicates=1, independent_reference_replicates=[1] if
                    s['design']['persons'] == 192 and s['null'] else []) for s in selected])
    write(out/'plan.json', plan)
    identity = hashlib.sha256((out/'plan.json').read_bytes()).hexdigest()
    write(out/'measurement.json', dict(scope=__doc__, utc=datetime.now(timezone.utc).isoformat(),
        source_root=str(root), formal_T0_sha256=formal_identity, development_plan_sha256=identity,
        formal_study_access='read only', numerical_workers=1, capacity=capacity))
    with zipfile.ZipFile(out/'sources.zip', 'x', zipfile.ZIP_DEFLATED) as archive:
        sources = {}
        for path in sorted(root.rglob('*')):
            if path.is_file():
                name = str(path.relative_to(root))
                sources[name] = controller.digest(path)
                archive.write(path, name)
        archive.write(__file__, 'profiling_driver.py')
        sources['profiling_driver.py'] = controller.digest(__file__)
    write(out/'source_manifest.json', sources)

    def verifier():
        with phase('source_and_python_runtime_verification'):
            verify()
        # The verified formal-source binding is not a development cohort's seal.

    for name, label in [('generate_shifted_sample','generation'),
                        ('prepare_retained_problem','input_validation_and_preparation'),
                        ('prepare_bias_problem','bias_problem_preparation'),
                        ('screen_from_point','conditional_screen'),
                        ('fit_grouped_joint','fit_other'), ('compare_orders','quadrature_comparison')]:
        wrap(executor, name, label)
    wrap(GroupedMmlLikelihood, '__init__', 'grouped_construction')
    wrap(GroupedMmlLikelihood, 'structural_recession_audit', 'prefit_recession_audit')
    original_evaluate = GroupedMmlLikelihood.evaluate
    def evaluate(self, *a, **kw):
        with phase('likelihood_information' if kw.get('information') else 'likelihood_gradient'):
            return original_evaluate(self, *a, **kw)
    GroupedMmlLikelihood.evaluate = evaluate
    for name in ['append', 'get', 'verify']:
        wrap(AttemptLedger, name, 'ledger_'+name)
    wrap(references, '_r_family', 'independent_R_process')
    records = []
    profiler = cProfile.Profile()
    try:
        for run in planned_bias_runs(plan):
            scenario = next(s for s in plan['scenarios'] if s['id'] == run['scenario_id'])
            ledger = AttemptLedger(out/(scenario['id']+'.sqlite'), plan_sha256=identity)
            profiler.enable()
            try:
                with phase('development_pipeline_other'):
                    executor.execute_bias_run(ledger, run, scenario, plan, verify_environment=verifier)
                    chain = ledger.verify()
                reference = None
                if scenario['independent_reference_replicates']:
                    target = out/'references'/run['run_id']
                    with phase('reference_export'):
                        references.export_reference(target, ledger, scenario, 1, plan, identity)
                    with phase('reference_execution_other'):
                        reference = references.execute_reference_bundle(target, plan, identity, root,
                                                                        controller.RSCRIPT, verifier)
                records.append(dict(run_id=run['run_id'], terminal=ledger.get(run['run_id'],'terminal'),
                                    chain=chain, reference=reference))
            finally:
                profiler.disable()
                ledger.close()
            write(out/'timing_checkpoint.json', dict(clock=clock, completed=len(records)))
            print(run['run_id'], len(records), '/', len(selected), flush=True)
        verifier()
        write(out/'summary.json', dict(classification='DEVELOPMENT_ONLY', scientific_inference_ready=False,
            clock=clock, records=records, free_bytes=shutil.disk_usage(out).free,
            timing_interpretation='Exclusive phase times partition measured scope; inclusive times overlap. R is deliberately oversampled.'))
    finally:
        profiler.dump_stats(str(out/'pipeline.prof'))
        with (out/'profile.txt').open('w') as handle:
            pstats.Stats(profiler, stream=handle).strip_dirs().sort_stats('cumulative').print_stats(65)
        write(out/'completion_manifest.json', {str(p.relative_to(out)):controller.digest(p)
            for p in out.rglob('*') if p.is_file()})


if __name__ == '__main__':
    main()
