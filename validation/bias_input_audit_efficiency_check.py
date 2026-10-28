"""Paired G3 pipelines with only the input-audit function changed in memory.

Independent development seeds; actual per-phase source/runtime checks and
durable writes; no cProfile overhead. The declared function substitution is
development instrumentation, not a formal source binding. Independent R is
excluded here and measured separately by bias_pipeline_profile.py.
"""
import argparse
import ast
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import statistics
import sys
import time
import zipfile


def audit_function(path, app):
    tree = ast.parse(path.read_text())
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'build_response_data_audit')
    namespace = dict(vars(app))
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), namespace)
    return namespace[node.name]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('study', type=Path)
    parser.add_argument('candidate_app', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    root, formal, candidate, out = (p.resolve() for p in
                                    (args.source,args.study,args.candidate_app,args.output))
    out.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(root))
    import streamlit_app as app
    from validation import bias_population_controller as controller
    from validation.bias_population_protocol import planned_bias_runs
    from validation.bias_population_study import execute_bias_run
    from validation.mml_population_study import AttemptLedger,json_bytes

    write = lambda path,value: path.write_bytes(json_bytes(value))
    original = app.build_response_data_audit
    versions = dict(baseline=original,candidate=audit_function(candidate,app))
    original_plan, formal_identity = controller.verify_seal(formal)
    verify = controller.environment_verifier(formal,original_plan)
    def verifier():
        verify()  # Read-only check of the sealed base; no formal binding in dev ledger.
    plan = dict(original_plan,classification='DEVELOPMENT_ONLY',qualification_execution_allowed=False,
        study_id='DEVELOPMENT_ONLY_g3_input_audit_efficiency_20261006',master_seed=2**117+2026100605,
        scenarios=[dict(s,replicates=1,independent_reference_replicates=[]) for s in original_plan['scenarios']])
    write(out/'plan.json',plan)
    identity=controller.digest(out/'plan.json')
    write(out/'measurement.json',dict(scope=__doc__,utc=datetime.now(timezone.utc).isoformat(),
        numerical_workers=1,formal_T0_sha256=formal_identity,development_plan_sha256=identity,
        comparison='Exact phase inventory and SHA256 of every event payload; timestamps and ledger chains differ.',
        timing='One pair per design, alternating order; includes both Q fits, four screens, phase verification and ledger writes; excludes imports/startup and independent R. Other formal workers remain active.'))
    sources={}
    with zipfile.ZipFile(out/'sources.zip','x',zipfile.ZIP_DEFLATED) as z:
        for path in sorted(root.rglob('*')):
            if path.is_file():
                name=str(path.relative_to(root));sources[name]=controller.digest(path);z.write(path,name)
        for name,path in [('candidate_app.py',candidate),('efficiency_driver.py',Path(__file__))]:
            sources[name]=controller.digest(path);z.write(path,name)
    write(out/'source_manifest.json',sources)
    ledgers={label:AttemptLedger(out/(label+'.sqlite'),plan_sha256=identity) for label in versions}
    records=[]
    try:
        for index,run in enumerate(planned_bias_runs(plan)):
            scenario=next(s for s in plan['scenarios'] if s['id']==run['scenario_id'])
            timings={}
            for label in (list(versions) if index%2==0 else list(reversed(versions))):
                app.build_response_data_audit=versions[label]
                start=time.perf_counter()
                execute_bias_run(ledgers[label],run,scenario,plan,verify_environment=verifier)
                timings[label]=time.perf_counter()-start
            hashes=[dict(l.connection.execute('SELECT phase,raw_sha256 FROM event WHERE run_id=?',
                                              (run['run_id'],))) for l in ledgers.values()]
            equal=hashes[0]==hashes[1] and 'terminal' in hashes[0]
            record=dict(run_id=run['run_id'],seconds=timings,time_ratio=timings['baseline']/timings['candidate'],
                        exact_payload_agreement=equal,event_payload_sha256=dict(zip(versions,hashes)))
            records.append(record)
            write(out/(scenario['id']+'.json'),record)
            print(run['run_id'],equal,round(record['time_ratio'],3),flush=True)
        app.build_response_data_audit=original
        verifier()
        write(out/'summary.json',dict(classification='DEVELOPMENT_ONLY',scientific_inference_ready=False,
            all_checks_pass=all(r['exact_payload_agreement'] for r in records),records=records,
            median_ratio=statistics.median(r['time_ratio'] for r in records),
            total_seconds={v:sum(r['seconds'][v] for r in records) for v in versions},
            ledger_chains={v:l.verify() for v,l in ledgers.items()}))
    finally:
        app.build_response_data_audit=original
        for ledger in ledgers.values():ledger.close()
        write(out/'artifact_manifest.json',{str(p.relative_to(out)):controller.digest(p)
             for p in out.rglob('*') if p.is_file()})


if __name__=='__main__':
    main()
