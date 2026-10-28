"""Seal, execute and collect the prespecified G3 conditional-screen study.

The controller never sends a timestamp request. Data generation requires its
verified receipt, isolated source, fixed runtimes, sufficient storage and
passing deterministic controls. Statistical qualification remains a separate
review even when every numerical job completes.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import csv
from datetime import datetime,timezone
import fcntl
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

from validation.bias_population_protocol import make_registered_bias_plan,planned_bias_runs
from validation.bias_population_study import execute_bias_run
from validation.bias_population_analysis import summarize_bias_scenario,assemble_evidence
from validation.bias_population_references import export_reference,execute_reference_bundle
from validation.bias_negative_controls import execute_controls
from validation.mml_population_study import AttemptLedger,json_bytes
from validation.mml_population_controller import (OPENSSL,RSCRIPT,TRUSTED_CA_SHA256,
    digest,write_new,python_runtime,r_runtime,environment_verifier)
from validation.mml_population_reference_execution import _once


def source_paths():
    # All app package/resources plus the actually imported validation graph.
    # Unrelated historical experiment drivers are not runtime dependencies.
    paths={ROOT/'streamlit_app.py',ROOT/'CITATION.cff',Path(__file__).resolve(),
        ROOT/'validation/bias_conditional_reference.R',ROOT/'validation/mml_population_runtime.R'}
    paths.update((ROOT/'mfrm_app').rglob('*.py'));paths.update((ROOT/'locales').glob('*.json'))
    for module in list(sys.modules.values()):
        filename=getattr(module,'__file__',None)
        if not filename:continue
        path=Path(filename).resolve()
        if path.is_relative_to(ROOT) and path.suffix=='.py' and path.relative_to(ROOT).parts[0]=='validation':
            paths.add(path)
    return sorted(paths)


def isolated_root(identity):return Path('/private/tmp')/('mfrm-bias-population-'+identity[:16])


def freeze(out):
    runtime=python_runtime()
    if runtime.get('dont_write_bytecode') is not True:raise ValueError('Freeze and execute with Python -B')
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False);plan=make_registered_bias_plan()
    write_new(out/'plan.json',plan)
    with (out/'runs.csv').open('x',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=['run_id','scenario_id','replicate','seed'],lineterminator='\n')
        writer.writeheader();writer.writerows(planned_bias_runs(plan))
    paths=source_paths()
    if any(p.is_symlink() for p in paths):raise ValueError('Source links are outside the frozen-copy contract')
    write_new(out/'source_manifest.json',{str(p.relative_to(ROOT)):digest(p) for p in paths})
    with zipfile.ZipFile(out/'sources.zip','x',zipfile.ZIP_DEFLATED) as archive:
        for p in paths:archive.write(p,str(p.relative_to(ROOT)))
    write_new(out/'python_runtime.json',runtime);write_new(out/'r_runtime.json',r_runtime())
    ca=ROOT/'validation/generated/g4_timestamp_probe_20261004/cacert.pem'
    if digest(ca)!=TRUSTED_CA_SHA256:raise ValueError('Timestamp trust root changed')
    (out/'cacert.pem').write_bytes(ca.read_bytes())
    write_new(out/'T0.json',dict(schema='bias-population-seal-v1',study_id=plan['study_id'],
        created_utc=datetime.now(timezone.utc).isoformat(),scientific_inference_ready=False,
        files={p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file()},
        openssl_sha256=digest(OPENSSL),trusted_ca_sha256=TRUSTED_CA_SHA256))
    subprocess.run([str(OPENSSL),'ts','-query','-data',str(out/'T0.json'),'-sha256','-cert','-out',str(out/'query.tsq')],check=True)
    return dict(study=str(out),T0_sha256=digest(out/'T0.json'),status='FROZEN_AWAITING_EXTERNAL_RECEIPT',
        planned_datasets=sum(s['replicates'] for s in plan['scenarios']))


def verify_seal(out,*,receipt=True):
    out=Path(out).resolve();seal=json.loads((out/'T0.json').read_bytes())
    if seal['schema']!='bias-population-seal-v1':raise ValueError('Different study seal schema')
    required={'plan.json','runs.csv','source_manifest.json','sources.zip','python_runtime.json','r_runtime.json','cacert.pem'}
    if set(seal['files'])!=required:raise ValueError('Incomplete frozen artifact inventory')
    for name,sha in seal['files'].items():
        if Path(name).name!=name or digest(out/name)!=sha:raise ValueError('Frozen artifact changed: '+name)
    if digest(OPENSSL)!=seal['openssl_sha256'] or seal['trusted_ca_sha256']!=TRUSTED_CA_SHA256:
        raise ValueError('Receipt verifier or trust anchor changed')
    if receipt:
        for flag,target in [('-queryfile','query.tsq'),('-data','T0.json')]:
            subprocess.run([str(OPENSSL),'ts','-verify',flag,str(out/target),'-in',str(out/'response.tsr'),
                '-CAfile',str(out/'cacert.pem')],capture_output=True,text=True,check=True)
    plan=json.loads((out/'plan.json').read_bytes())
    if seal['study_id']!=plan['study_id']:raise ValueError('Study identity differs from seal')
    if json.loads((out/'python_runtime.json').read_bytes()).get('dont_write_bytecode') is not True:
        raise ValueError('Frozen Python runtime allows source-tree bytecode writes')
    if json_bytes(plan)!=json_bytes(make_registered_bias_plan()):raise ValueError('Plan differs from the frozen specification')
    with (out/'runs.csv').open(newline='',encoding='utf-8') as handle:
        for actual,expected in zip(csv.DictReader(handle),planned_bias_runs(plan),strict=True):
            if actual!={k:str(v) for k,v in expected.items()}:raise ValueError('Ordered RunId/seed inventory changed')
    return plan,digest(out/'T0.json')


def restore(out):
    out=Path(out).resolve();_,identity=verify_seal(out);destination=isolated_root(identity)
    hashes=json.loads((out/'source_manifest.json').read_bytes())
    if not destination.exists():
        staging=Path(tempfile.mkdtemp(prefix='mfrm-bias-restoring-',dir='/private/tmp'))
        with zipfile.ZipFile(out/'sources.zip') as archive:
            if len(archive.namelist())!=len(hashes) or set(archive.namelist())!=set(hashes):raise ValueError('Source archive entries changed')
            for name,sha in hashes.items():
                relative=Path(name)
                if relative.is_absolute() or '..' in relative.parts:raise ValueError('Nonlocal source archive entry')
                target=staging/relative;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(archive.read(name))
                if digest(target)!=sha:raise ValueError('Restored source changed: '+name)
                target.chmod(0o444)
        staging.rename(destination)
    if {str(p.relative_to(destination)) for p in destination.rglob('*') if p.is_file()}!=set(hashes):
        raise ValueError('Unexpected files in isolated source root')
    for name,sha in hashes.items():
        if digest(destination/name)!=sha:raise ValueError('Restored source changed: '+name)
    return destination


def storage_requirement(out,plan,identity):
    """Remaining includes in-flight runs; a start alone does not free capacity."""
    remaining=0
    for scenario in plan['scenarios']:
        expected={f"{scenario['id']}__{i:04d}" for i in range(1,scenario['replicates']+1)}
        path=Path(out)/'attempts'/(scenario['id']+'.sqlite');terminal=set()
        if path.exists():
            ledger=AttemptLedger.read_only(path,plan_sha256=identity)
            try:
                existing={r[0] for r in ledger.connection.execute('SELECT DISTINCT run_id FROM event')}
                if existing-expected:raise ValueError('Unregistered attempt in storage inventory')
                terminal={r[0] for r in ledger.connection.execute("SELECT run_id FROM event WHERE phase='terminal'")}
            finally:ledger.close()
        remaining+=len(expected-terminal)
    settings=plan['execution']
    return dict(remaining_replicates=remaining,
        required_free_bytes=remaining*settings['planning_bytes_per_remaining_replicate']
            +settings['reference_storage_allowance_bytes']+settings['disk_reserve_bytes'],
        interpretation='Storage planning allowance, not a guarantee of final file size; the per-phase reserve remains enforced')


def capacity_check(out,plan,identity):
    result=storage_requirement(out,plan,identity);result['available_free_bytes']=shutil.disk_usage(out).free
    result['pass_']=result['available_free_bytes']>=result['required_free_bytes']
    if not result['pass_']:raise RuntimeError('Insufficient planned storage headroom: '+json.dumps(result))
    return result


def _write_completed(path,value):
    if path.exists():
        if path.read_bytes()!=json_bytes(value)+b'\n':raise ValueError('Completed artifact differs; no replacement: '+str(path))
    else:write_new(path,value)


def run_scenario(out,scenario_id):
    out=Path(out).resolve();plan,identity=verify_seal(out)
    if ROOT!=isolated_root(identity):raise ValueError('Scenario execution requires isolated frozen source')
    scenario=next(s for s in plan['scenarios'] if s['id']==scenario_id)
    capacity_check(out,plan,identity);verifier=environment_verifier(out,plan)
    if not (out/'controls/summary.json').exists():raise RuntimeError('Deterministic controls must complete before any scenario')
    controls=execute_controls(out/'controls',plan,identity,ROOT,RSCRIPT,verifier)
    if not controls['all_checks_pass']:raise RuntimeError('Deterministic controls failed; no scenario generation authorized')
    attempts=out/'attempts';attempts.mkdir(exist_ok=True)
    with (attempts/(scenario_id+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        ledger=AttemptLedger(attempts/(scenario_id+'.sqlite'),plan_sha256=identity)
        try:
            ledger.verify()
            for registered in planned_bias_runs(dict(plan,scenarios=[scenario])):
                key=registered['run_id']
                if ledger.started(key):
                    if ledger.get(key,'terminal') is None:
                        ledger.append(key,'terminal',dict(status='interrupted_attempt_retained',numerical_replay_prohibited=True,
                            truth_available=ledger.get(key,'input_retained') is not None,scientific_inference_ready=False))
                    continue
                execute_bias_run(ledger,registered,scenario,plan,verify_environment=verifier)
                if registered['replicate']==1 or registered['replicate']%100==0:
                    print(scenario_id,registered['replicate'],'/',scenario['replicates'],'attempt retained',flush=True)
            verifier();summary=summarize_bias_scenario(ledger,scenario,plan)
            summary.update(study_sha256=identity,source_ledger=ledger.verify())
            (out/'summaries').mkdir(exist_ok=True);_write_completed(out/'summaries'/(scenario_id+'.json'),summary)
            return dict(scenario_id=scenario_id,status=summary['status'],ledger=ledger.verify())
        finally:ledger.close()


def _references(out,plan,identity,verifier):
    results=[];batch=AttemptLedger(out/'reference_exports.sqlite',plan_sha256=identity)
    try:
        batch.verify()
        for scenario in plan['scenarios']:
            path=out/'attempts'/(scenario['id']+'.sqlite')
            if not path.exists():raise ValueError('Unstarted required scenario: '+scenario['id'])
            ledger=AttemptLedger.read_only(path,plan_sha256=identity)
            try:
                ledger.verify()
                if any(ledger.get(r['run_id'],'terminal') is None for r in planned_bias_runs(dict(plan,scenarios=[scenario]))):
                    raise ValueError('Reference batch requires every planned scenario attempt terminal')
                for replicate in scenario['independent_reference_replicates']:
                    key=f"{scenario['id']}__{replicate:04d}";target=out/'references'/key
                    required=len(plan['scope']['lanes'])*(len(scenario['design']['rater'])*len(scenario['design']['task'])+4)
                    def export():
                        export_reference(target,ledger,scenario,replicate,plan,identity)
                        return dict(manifest_sha256=digest(target/'manifest.json'))
                    exported=_once(batch,key,'export',export,verifier)
                    if exported['status']!='completed':
                        results.append(dict(run_id=key,status=exported['status'],required_slots=required,passed_slots=0,
                            all_checks_pass=False,stage='export'));continue
                    if digest(target/'manifest.json')!=exported['value']['manifest_sha256']:
                        raise ValueError('Retained reference export changed')
                    results.append(execute_reference_bundle(target,plan,identity,ROOT,RSCRIPT,verifier))
            finally:ledger.close()
        return results
    finally:batch.close()


def run(out,workers):
    out=Path(out).resolve();plan,identity=verify_seal(out)
    if ROOT!=isolated_root(identity):raise ValueError('Execute from the restored frozen source root before generating data')
    if type(workers) is not int or not 1<=workers<=plan['execution']['maximum_workers']:raise ValueError('Unregistered worker count')
    with (out/'controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        capacity_check(out,plan,identity);verifier=environment_verifier(out,plan)
        controls=execute_controls(out/'controls',plan,identity,ROOT,RSCRIPT,verifier)
        if not controls['all_checks_pass']:raise RuntimeError('Deterministic controls failed; no population generation authorized')
        with ProcessPoolExecutor(max_workers=workers) as pool:
            jobs=[pool.submit(run_scenario,str(out),s['id']) for s in plan['scenarios']]
            for job in as_completed(jobs):
                try:print(json.dumps(job.result()),flush=True)
                except Exception:
                    for pending in jobs:pending.cancel()
                    raise
        references=_references(out,plan,identity,verifier)
        summaries=[json.loads((out/'summaries'/(s['id']+'.json')).read_bytes()) for s in plan['scenarios']]
        packet=assemble_evidence(plan,summaries,references,controls)
        packet['study_sha256']=identity
        _write_completed(out/'evidence_packet.json',packet)
        return dict(status='all_planned_attempts_terminal',scientific_inference_ready=False,study_sha256=identity)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=['freeze','verify','restore','run','capacity'])
    parser.add_argument('study',type=Path);parser.add_argument('--workers',type=int,default=2);args=parser.parse_args()
    if args.command=='freeze':result=freeze(args.study)
    elif args.command=='verify':result=dict(T0_sha256=verify_seal(args.study)[1])
    elif args.command=='restore':result=dict(source_root=str(restore(args.study)))
    elif args.command=='capacity':
        plan,identity=verify_seal(args.study,receipt=False)
        result=storage_requirement(args.study,plan,identity);result['available_free_bytes']=shutil.disk_usage(args.study).free
    else:result=run(args.study,args.workers)
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
