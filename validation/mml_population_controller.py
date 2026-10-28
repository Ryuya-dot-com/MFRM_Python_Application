"""Freeze, verify and execute the prospectively specified population study.

Only SHA256+nonce requests leave the Mac. The sealed manifest includes the
full source bundle, protocol, ordered seeds, runtimes and timestamp root.
Actual data generation requires both verified RFC3161 receipt and an isolated
source copy. Worker processes are numerical jobs, not autonomous agents.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

from validation.mml_free_sd_prefit_transport import live_python_environment
from validation.mml_population_protocol import make_plan, planned_runs
from validation.mml_population_study import AttemptLedger, execute_registered_run, json_bytes
from validation.mml_population_analysis import summarize_scenario
from validation.mml_execution_amendment import expected_root, verified_sources


OPENSSL = Path('/opt/homebrew/bin/openssl')
RSCRIPT = Path('/usr/local/bin/Rscript')
TRUSTED_CA_SHA256 = '2151b61137ffa86bf664691ba67e7da0b19f98c758e3d228d5d8ebf27e044438'
THREAD_ENV = dict(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',MKL_NUM_THREADS='1')


def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_new(path,value):
    with Path(path).open('xb') as handle:
        handle.write(json_bytes(value)+b'\n');handle.flush();os.fsync(handle.fileno())


def python_runtime():
    if any(os.environ.get(k) != v for k,v in THREAD_ENV.items()):
        raise ValueError('Launch with the four BLAS/OpenMP thread environment variables fixed at one')
    return dict(environment=live_python_environment(),thread_environment={k:os.environ[k] for k in THREAD_ENV},
        executable=str(Path(sys.executable).resolve()),dont_write_bytecode=bool(sys.dont_write_bytecode))


def r_runtime():
    result=subprocess.run([str(RSCRIPT),'--vanilla',str(ROOT/'validation/mml_population_runtime.R')],
        capture_output=True,text=True,check=True)
    return dict(executable_sha256=digest(RSCRIPT),runtime=json.loads(result.stdout),stderr=result.stderr)


def source_paths():
    # Include complete local Python packages and the validation source layer;
    # generated datasets, historical result folders and venvs are excluded.
    return sorted({ROOT/'streamlit_app.py',ROOT/'CITATION.cff',*list((ROOT/'mfrm_app').rglob('*.py')),
        *list((ROOT/'validation').glob('*.py')),*list((ROOT/'validation').glob('*.R')),
        *list((ROOT/'locales').glob('*.json'))})


def freeze(out):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    plan=make_plan();write_new(out/'plan.json',plan)
    with (out/'runs.csv').open('x',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=['run_id','scenario_id','replicate','seed','cohort'],lineterminator='\n')
        writer.writeheader();writer.writerows(planned_runs(plan));handle.flush();os.fsync(handle.fileno())
    paths=source_paths()
    if any(p.is_symlink() for p in paths): raise ValueError('Source links are outside the frozen-copy contract')
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    write_new(out/'source_manifest.json',hashes)
    with zipfile.ZipFile(out/'sources.zip','x',zipfile.ZIP_DEFLATED) as archive:
        for p in paths:archive.write(p,str(p.relative_to(ROOT)))
    write_new(out/'python_runtime.json',python_runtime());write_new(out/'r_runtime.json',r_runtime())
    ca=ROOT/'validation/generated/g4_timestamp_probe_20261004/cacert.pem'
    if digest(ca)!=TRUSTED_CA_SHA256: raise ValueError('Timestamp root differs from the explicitly trusted probe root')
    shutil.copyfile(ca,out/'cacert.pem')
    manifest=dict(schema='joint-mml-population-seal-v1',study_id=plan['study_id'],
        created_utc=datetime.now(timezone.utc).isoformat(),scientific_inference_ready=False,
        files={p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file()},
        openssl_executable=str(OPENSSL),openssl_sha256=digest(OPENSSL),rscript_executable=str(RSCRIPT),
        trusted_ca_sha256=TRUSTED_CA_SHA256,protocol_identity='Exact bytes of this T0.json are timestamped')
    write_new(out/'T0.json',manifest)
    subprocess.run([str(OPENSSL),'ts','-query','-data',str(out/'T0.json'),'-sha256','-cert','-out',str(out/'query.tsq')],check=True)
    print(json.dumps(dict(study=str(out),T0_sha256=digest(out/'T0.json'),
        status='FROZEN_AWAITING_EXTERNAL_RECEIPT',datasets=plan['population_datasets'],negative_controls=plan['negative_controls'])),flush=True)


def verify_seal(out, *, receipt=True):
    out=Path(out).resolve();seal=json.loads((out/'T0.json').read_bytes())
    for name,sha in seal['files'].items():
        if Path(name).name!=name or digest(out/name)!=sha: raise ValueError('Frozen artifact changed: '+name)
    if digest(OPENSSL)!=seal['openssl_sha256'] or seal['trusted_ca_sha256']!=TRUSTED_CA_SHA256:
        raise ValueError('Receipt verifier or trust anchor changed')
    if receipt:
        for target in [('-queryfile',out/'query.tsq'),('-data',out/'T0.json')]:
            subprocess.run([str(OPENSSL),'ts','-verify',target[0],str(target[1]),'-in',str(out/'response.tsr'),
                '-CAfile',str(out/'cacert.pem')],capture_output=True,text=True,check=True)
    plan=json.loads((out/'plan.json').read_bytes())
    if json_bytes(plan)!=json_bytes(make_plan()): raise ValueError('Protocol differs from the frozen controller specification')
    with (out/'runs.csv').open(newline='',encoding='utf-8') as handle:
        rows=csv.DictReader(handle)
        for actual,expected in zip(rows,planned_runs(plan),strict=True):
            if actual!={k:str(v) for k,v in expected.items()}: raise ValueError('Ordered RunId/seed plan differs')
    return plan,digest(out/'T0.json')


def restore(out):
    plan,identity=verify_seal(out)
    destination=Path('/private/tmp')/('mfrm-population-'+identity[:16])
    hashes=json.loads((out/'source_manifest.json').read_bytes())
    if not destination.exists():
        staging=Path(tempfile.mkdtemp(prefix='mfrm-population-restoring-',dir='/private/tmp'))
        with zipfile.ZipFile(out/'sources.zip') as archive:
            if len(archive.namelist())!=len(hashes) or set(archive.namelist())!=set(hashes): raise ValueError('Source archive entries differ from manifest')
            for name in hashes:
                relative=Path(name)
                if relative.is_absolute() or '..' in relative.parts: raise ValueError('Nonlocal source archive path')
                target=staging/relative;target.parent.mkdir(parents=True,exist_ok=True)
                with target.open('xb') as handle:handle.write(archive.read(name))
                if digest(target)!=hashes[name]:raise ValueError('Restored source differs: '+name)
        for p in staging.rglob('*'):
            if p.is_file():p.chmod(0o444)
        staging.rename(destination)
    for name,sha in hashes.items():
        if digest(destination/name)!=sha: raise ValueError('Restored source differs: '+name)
    if {str(p.relative_to(destination)) for p in destination.rglob('*') if p.is_file()}!=set(hashes):
        raise ValueError('Unexpected files in isolated source root')
    return destination


def environment_verifier(out,plan,*,check_R=True):
    hashes=json.loads((out/'source_manifest.json').read_bytes())
    amended=verified_sources(out,ROOT,OPENSSL,TRUSTED_CA_SHA256)
    binding=None
    if amended is not None:hashes,binding=amended
    expected_python=json.loads((out/'python_runtime.json').read_bytes())
    if check_R and json_bytes(r_runtime())!=json_bytes(json.loads((out/'r_runtime.json').read_bytes())):
        raise ValueError('R/TAM runtime differs from sealed environment')
    def verify():
        if shutil.disk_usage(out).free < plan['execution']['disk_reserve_bytes']:
            raise RuntimeError('Disk reserve reached before starting additional work; no failed fit is retried')
        for name,sha in hashes.items():
            if digest(ROOT/name)!=sha: raise ValueError('Live frozen source changed: '+name)
        if json_bytes(python_runtime())!=json_bytes(expected_python): raise ValueError('Python source/runtime changed')
        if binding is not None:
            if digest(out/'implementation_activation.json')!=binding['activation_sha256']:
                raise ValueError('Execution activation changed during a run')
            return binding
    verify();return verify


def run_scenario(out,scenario_id):
    out=Path(out);plan,identity=verify_seal(out)
    scenario=next(s for s in plan['scenarios'] if s['id']==scenario_id)
    verifier=environment_verifier(out,plan)
    state=out/'attempts';state.mkdir(exist_ok=True)
    with (state/(scenario_id+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        ledger=AttemptLedger(state/(scenario_id+'.sqlite'),plan_sha256=identity)
        ledger.verify()
        try:
            for run in planned_runs(dict(plan,scenarios=[scenario])):
                if ledger.started(run['run_id']):
                    if ledger.get(run['run_id'],'terminal') is None:
                        ledger.append(run['run_id'],'terminal',dict(status='interrupted_attempt_retained',
                            truth_available=ledger.get(run['run_id'],'input_retained') is not None,
                            numerical_replay_prohibited=True,scientific_inference_ready=False))
                    continue
                execute_registered_run(ledger,run,scenario,plan,verify_environment=verifier)
                if run['replicate']==1 or run['replicate'] % 100==0:
                    print(scenario_id,run['replicate'],'/',scenario['replicates'],'attempt retained',flush=True)
            verifier()
            summary=summarize_scenario(ledger,scenario,plan)
            completed=out/'summaries';completed.mkdir(exist_ok=True)
            path=completed/(scenario_id+'.json')
            if path.exists():
                if path.read_bytes()!=json_bytes(summary)+b'\n': raise ValueError('Completed summary differs; no replacement')
            else:write_new(path,summary)
            return dict(scenario_id=scenario_id,status=summary['status'],ledger=ledger.verify())
        finally:ledger.close()


def run(out,workers):
    out=Path(out).resolve();plan,identity=verify_seal(out)
    if ROOT!=expected_root(out,identity,Path('/private/tmp')/('mfrm-population-'+identity[:16])):
        raise ValueError('Execute from the restored frozen source root before generating any data')
    if not 1 <= workers <= plan['execution']['maximum_workers']: raise ValueError('Worker count exceeds the registered limit')
    environment_verifier(out,plan)
    with (out/'controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        # Controls are run first, but positive cases are not relabelled or removed
        # if a control fails. All planned cases retain their original cohorts.
        scenarios=sorted(plan['scenarios'],key=lambda s:s['cohort']!='negative_control')
        with ProcessPoolExecutor(max_workers=workers) as pool:
            jobs={pool.submit(run_scenario,str(out),s['id']):s['id'] for s in scenarios}
            for job in as_completed(jobs):
                try: print(json.dumps(job.result()),flush=True)
                except Exception:
                    for pending in jobs:pending.cancel()
                    raise
        print('ALL_REGISTERED_SCENARIOS_TERMINAL; REFERENCES AND SCIENTIFIC REVIEW REMAIN',flush=True)


def export_references(out):
    from validation.mml_population_references import export_reference_bundle
    out=Path(out).resolve();plan,identity=verify_seal(out)
    if ROOT!=expected_root(out,identity,Path('/private/tmp')/('mfrm-population-'+identity[:16])):
        raise ValueError('Export references from the isolated frozen source root')
    verifier=environment_verifier(out,plan);summaries=[]
    for scenario in plan['scenarios']:
        selected=sorted(set(scenario['independent_R_replicates']+scenario['TAM_replicates']))
        if not selected:continue
        ledger_path=out/'attempts'/(scenario['id']+'.sqlite')
        if not ledger_path.exists():continue  # Unstarted scenarios stay pending.
        ledger=AttemptLedger.read_only(ledger_path,plan_sha256=identity)
        try:
            ledger.verify()
            for replicate in selected:
                run_id=f"{scenario['id']}__{replicate:04d}"
                if ledger.get(run_id,'terminal') is None:continue
                target=out/'references'/run_id
                if target.exists():
                    manifest=target/'manifest.json'
                    if not manifest.is_file():raise ValueError('Interrupted reference export retained; no replacement: '+run_id)
                    for name,sha in json.loads(manifest.read_bytes()).items():
                        if digest(target/name)!=sha:raise ValueError('Existing reference input changed: '+name)
                    if json.loads((target/'reference_plan.json').read_bytes())['study_sha256']!=identity:
                        raise ValueError('Reference export belongs to another study')
                    continue
                verifier();record=export_reference_bundle(target,ledger,scenario,replicate,plan,identity)
                summaries.append(record);print(json.dumps(record),flush=True)
        finally:ledger.close()
    return summaries


def run_references(out):
    from validation.mml_population_reference_execution import execute_bundle
    out=Path(out).resolve();plan,identity=verify_seal(out)
    if ROOT!=expected_root(out,identity,Path('/private/tmp')/('mfrm-population-'+identity[:16])):
        raise ValueError('Execute references from the isolated frozen source root')
    with (out/'controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        verifier=environment_verifier(out,plan)
        export_references(out)
        for scenario in plan['scenarios']:
            for replicate in sorted(set(scenario['independent_R_replicates']+scenario['TAM_replicates'])):
                target=out/'references'/f"{scenario['id']}__{replicate:04d}"
                if (target/'manifest.json').is_file():
                    result=execute_bundle(target,plan,identity,ROOT,RSCRIPT,verifier)
                    print(json.dumps(dict(run_id=result['run_id'],all_checks_pass=result['all_checks_pass'])),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=['freeze','verify','restore','run','export-references','run-references'])
    parser.add_argument('study',type=Path);parser.add_argument('--workers',type=int,default=2);args=parser.parse_args()
    if args.command=='freeze':freeze(args.study)
    elif args.command=='verify':print(verify_seal(args.study)[1])
    elif args.command=='restore':print(restore(args.study.resolve()))
    elif args.command=='export-references':export_references(args.study)
    elif args.command=='run-references':run_references(args.study)
    else:
        run(args.study,args.workers)
        run_references(args.study)


if __name__=='__main__':main()
