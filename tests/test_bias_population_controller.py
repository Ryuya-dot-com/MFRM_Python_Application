"""No signed receipt or headroom: no generation. Started attempts stay closed."""
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from validation import bias_population_controller as controller
from validation.bias_population_protocol import make_registered_bias_plan,planned_bias_runs
from validation.mml_population_study import AttemptLedger,json_bytes


def fixture_plan():
    plan=make_registered_bias_plan();plan.update(classification='DEVELOPMENT_ONLY',
        study_id='DEVELOPMENT_ONLY_BIAS_CONTROLLER',master_seed=713580023,qualification_execution_allowed=False)
    scenario=plan['scenarios'][0];scenario.update(replicates=2,independent_reference_replicates=[1])
    plan['scenarios']=[scenario];plan['planned_population_datasets']=2
    plan['execution'].update(planning_bytes_per_remaining_replicate=16384,reference_storage_allowance_bytes=8192,disk_reserve_bytes=4096)
    return plan,scenario


def fixture_seal(tmp_path,monkeypatch):
    plan,_=fixture_plan();monkeypatch.setattr(controller,'make_registered_bias_plan',lambda:plan)
    monkeypatch.setattr(controller,'python_runtime',lambda:{'development_fixture':True,'dont_write_bytecode':True})
    monkeypatch.setattr(controller,'r_runtime',lambda:{'development_fixture':True})
    monkeypatch.setattr(controller,'source_paths',lambda:[controller.ROOT/'CITATION.cff',Path(controller.__file__).resolve()])
    out=tmp_path/'study';controller.freeze(out)
    return out,plan


def test_runtime_source_inventory_contains_plan_controller_and_version_asset():
    paths=controller.source_paths()
    for relative in ['CITATION.cff','streamlit_app.py','validation/bias_population_controller.py',
        'validation/bias_population_protocol.py','validation/bias_population_study.py',
        'validation/bias_population_references.py','validation/bias_conditional_reference.R']:
        assert controller.ROOT/relative in paths


def test_absent_or_probe_receipt_never_reaches_generation(tmp_path,monkeypatch):
    out,plan=fixture_seal(tmp_path,monkeypatch)
    assert json_bytes(controller.verify_seal(out,receipt=False)[0])==json_bytes(plan)
    with pytest.raises(subprocess.CalledProcessError):controller.run(out,1)
    probe=controller.ROOT/'validation/generated/g4_timestamp_probe_20261004/response.tsr'
    (out/'response.tsr').write_bytes(probe.read_bytes())
    with pytest.raises(subprocess.CalledProcessError):controller.run(out,1)
    assert not (out/'attempts').exists() and not (out/'controls').exists()


def test_seed_inventory_cannot_be_changed_by_updating_unsigned_local_manifest(tmp_path,monkeypatch):
    out,_=fixture_seal(tmp_path,monkeypatch)
    lines=(out/'runs.csv').read_text().splitlines();cells=lines[1].split(',');cells[-1]='123';lines[1]=','.join(cells)
    (out/'runs.csv').write_text('\n'.join(lines)+'\n')
    seal=json.loads((out/'T0.json').read_bytes());seal['files']['runs.csv']=controller.digest(out/'runs.csv')
    (out/'T0.json').write_bytes(json_bytes(seal))
    with pytest.raises(ValueError,match='Ordered RunId/seed'):controller.verify_seal(out,receipt=False)


def test_seal_cannot_omit_required_source_or_runtime_bindings(tmp_path,monkeypatch):
    out,_=fixture_seal(tmp_path,monkeypatch)
    seal=json.loads((out/'T0.json').read_bytes());del seal['files']['source_manifest.json']
    (out/'T0.json').write_bytes(json_bytes(seal))
    with pytest.raises(ValueError,match='Incomplete frozen artifact inventory'):controller.verify_seal(out,receipt=False)


def test_registration_does_not_change_the_statistical_design_or_tolerances():
    from validation.bias_population_protocol import make_bias_plan
    draft=make_bias_plan();registered=make_registered_bias_plan()
    for key in ['study_id','master_seed','scenarios','scope','estimator','screen','endpoints','precision',
        'draft_interpretation_limits','independent_reference','deterministic_controls']:
        assert json_bytes(draft[key])==json_bytes(registered[key])
    assert not draft['qualification_execution_allowed'] and registered['qualification_execution_allowed']
    assert not registered['scientific_inference_ready']


def test_restore_requires_exact_inventory_and_no_extra_files(tmp_path,monkeypatch):
    out,_=fixture_seal(tmp_path,monkeypatch)
    # A development-only local extraction exercise, not a receipt substitute
    # for any qualification cohort. No execution/generator is invoked here.
    original=controller.verify_seal
    monkeypatch.setattr(controller,'verify_seal',lambda p:original(p,receipt=False))
    isolated=tmp_path/'isolated';monkeypatch.setattr(controller,'isolated_root',lambda identity:isolated)
    assert controller.restore(out)==isolated
    (isolated/'unexpected.py').write_text('# extra')
    with pytest.raises(ValueError,match='Unexpected files'):controller.restore(out)


def test_capacity_counts_inflight_until_terminal_and_never_drops_unregistered_runs(tmp_path):
    plan,scenario=fixture_plan();runs=list(planned_bias_runs(plan));(tmp_path/'attempts').mkdir()
    base=controller.storage_requirement(tmp_path,plan,'fixture')
    ledger=AttemptLedger(tmp_path/'attempts'/(scenario['id']+'.sqlite'),plan_sha256='fixture')
    ledger.append(runs[0]['run_id'],'generation_started',runs[0])
    assert controller.storage_requirement(tmp_path,plan,'fixture')==base
    ledger.append(runs[0]['run_id'],'terminal',dict(status='generation_failed'))
    after=controller.storage_requirement(tmp_path,plan,'fixture')
    assert base['remaining_replicates']==2 and after['remaining_replicates']==1
    assert base['required_free_bytes']-after['required_free_bytes']==16384
    ledger.append('foreign','terminal',dict(status='generation_failed'))
    with pytest.raises(ValueError,match='Unregistered attempt'):controller.storage_requirement(tmp_path,plan,'fixture')
    ledger.close()


def test_insufficient_capacity_stops_before_control_or_generation(tmp_path,monkeypatch):
    plan,_=fixture_plan()
    monkeypatch.setattr(controller,'verify_seal',lambda out:(plan,'fixture'))
    monkeypatch.setattr(controller,'isolated_root',lambda identity:controller.ROOT)
    monkeypatch.setattr(controller.shutil,'disk_usage',lambda out:SimpleNamespace(free=4096))
    with pytest.raises(RuntimeError,match='Insufficient planned storage'):controller.run(tmp_path,1)
    assert not (tmp_path/'attempts').exists() and not (tmp_path/'controls').exists()


def test_failed_deterministic_control_stops_before_population_attempts(tmp_path,monkeypatch):
    plan,_=fixture_plan()
    monkeypatch.setattr(controller,'verify_seal',lambda out:(plan,'fixture'))
    monkeypatch.setattr(controller,'isolated_root',lambda identity:controller.ROOT)
    monkeypatch.setattr(controller,'capacity_check',lambda *args:None)
    monkeypatch.setattr(controller,'environment_verifier',lambda *args:lambda:None)
    monkeypatch.setattr(controller,'execute_controls',lambda *args:dict(all_checks_pass=False))
    with pytest.raises(RuntimeError,match='controls failed'):controller.run(tmp_path,1)
    assert not (tmp_path/'attempts').exists()


def test_worker_closes_interruption_and_executes_only_unstarted_run(tmp_path,monkeypatch):
    plan,scenario=fixture_plan();runs=list(planned_bias_runs(plan));attempts=tmp_path/'attempts';attempts.mkdir()
    ledger=AttemptLedger(attempts/(scenario['id']+'.sqlite'),plan_sha256='fixture')
    ledger.append(runs[0]['run_id'],'generation_started',runs[0]);ledger.close();called=[]
    monkeypatch.setattr(controller,'verify_seal',lambda out:(plan,'fixture'))
    monkeypatch.setattr(controller,'isolated_root',lambda identity:controller.ROOT)
    monkeypatch.setattr(controller,'capacity_check',lambda *args:None)
    monkeypatch.setattr(controller,'environment_verifier',lambda *args:lambda:None)
    (tmp_path/'controls').mkdir();(tmp_path/'controls/summary.json').write_text('{}')
    monkeypatch.setattr(controller,'execute_controls',lambda *args:dict(all_checks_pass=True))
    def execute(ledger,run,*args,**kwargs):
        called.append(run['run_id']);ledger.append(run['run_id'],'generation_started',run)
        ledger.append(run['run_id'],'terminal',dict(status='generation_failed'))
    monkeypatch.setattr(controller,'execute_bias_run',execute)
    result=controller.run_scenario(tmp_path,scenario['id'])
    assert result['status']=='complete' and called==[runs[1]['run_id']]
    summary=json.loads((tmp_path/'summaries'/(scenario['id']+'.json')).read_bytes())
    assert summary['execution_status_counts']=={'generation_failed':1,'interrupted_attempt_retained':1}
    assert summary['study_sha256']=='fixture' and summary['source_ledger']==result['ledger']
    assert controller.run_scenario(tmp_path,scenario['id'])==result and called==[runs[1]['run_id']]


def test_direct_worker_entry_cannot_skip_deterministic_controls(tmp_path,monkeypatch):
    plan,scenario=fixture_plan()
    monkeypatch.setattr(controller,'verify_seal',lambda out:(plan,'fixture'))
    monkeypatch.setattr(controller,'isolated_root',lambda identity:controller.ROOT)
    monkeypatch.setattr(controller,'capacity_check',lambda *args:None)
    monkeypatch.setattr(controller,'environment_verifier',lambda *args:lambda:None)
    with pytest.raises(RuntimeError,match='controls must complete'):controller.run_scenario(tmp_path,scenario['id'])
    assert not (tmp_path/'attempts').exists()


def test_interrupted_reference_export_is_retained_as_required_failure(tmp_path):
    plan,scenario=fixture_plan();attempts=tmp_path/'attempts';attempts.mkdir()
    ledger=AttemptLedger(attempts/(scenario['id']+'.sqlite'),plan_sha256='fixture')
    for run in planned_bias_runs(plan):
        ledger.append(run['run_id'],'generation_started',run);ledger.append(run['run_id'],'terminal',dict(status='generation_failed'))
    ledger.close();exports=AttemptLedger(tmp_path/'reference_exports.sqlite',plan_sha256='fixture')
    key=list(planned_bias_runs(plan))[0]['run_id'];exports.append(key,'export_started',dict(single_attempt=True));exports.close()
    result=controller._references(tmp_path,plan,'fixture',lambda:None)
    assert result==[dict(run_id=key,status='interrupted_reference_retained',required_slots=40,passed_slots=0,all_checks_pass=False,stage='export')]
