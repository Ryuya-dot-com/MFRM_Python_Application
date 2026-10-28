"""Independent family construction and failed-reference denominators."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import pytest

from validation import bias_population_references as reference
from validation import bias_population_study as study
from validation.bias_population_protocol import make_bias_plan
from validation.mfrm_recovery_design import replicate_seed
from validation.mml_population_protocol import design
from validation.mml_population_study import AttemptLedger,json_bytes

ROOT=Path(__file__).resolve().parents[1]
RSCRIPT=Path('/usr/local/bin/Rscript')


@pytest.fixture
def bundle(tmp_path,monkeypatch):
    plan=make_bias_plan();plan.update(classification='DEVELOPMENT_ONLY',study_id='DEVELOPMENT_ONLY_G3_REFERENCE_TEST',master_seed=91450031)
    scenario=plan['scenarios'][0];scenario.update(id='RSM_reference_test',design=asdict(design('RSM','regular',persons=24)),
        replicates=1,independent_reference_replicates=[1]);plan['scenarios']=[scenario]
    run=dict(run_id=scenario['id']+'__0001',scenario_id=scenario['id'],replicate=1,
        seed=replicate_seed(plan['master_seed'],plan['study_id']+':'+scenario['id'],1))
    def fail_fit(*args,**kwargs):raise ArithmeticError('Injected development calibration failure')
    monkeypatch.setattr(study,'fit_grouped_joint',fail_fit)
    ledger=AttemptLedger(tmp_path/'data.sqlite',plan_sha256='fixture')
    study.execute_bias_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    out=tmp_path/'reference';reference.export_reference(out,ledger,scenario,1,plan,'fixture');ledger.close()
    return out,plan


@pytest.mark.skipif(not RSCRIPT.exists(),reason='Independent R unavailable')
def test_reference_uses_all_required_slots_and_checks_cached_artifacts(bundle):
    out,plan=bundle
    result=reference.execute_reference_bundle(out,plan,'fixture',ROOT,RSCRIPT,lambda:None)
    assert result['required_slots']==40 and result['passed_slots']==20
    assert not result['all_checks_pass'] and not result['scientific_inference_ready']
    assert [r['status'] for r in result['families']]==['completed','completed','required_input_unavailable','required_input_unavailable']
    again=reference.execute_reference_bundle(out,plan,'fixture',ROOT,RSCRIPT,lambda:None)
    assert again==result
    artifact=out/'execution/known_calibration_known_person/result.json'
    artifact.write_text('{}')
    with pytest.raises(ValueError,match='Retained R artifact changed'):
        reference.execute_reference_bundle(out,plan,'fixture',ROOT,RSCRIPT,lambda:None)


def test_failed_R_is_never_rerun_or_removed_from_denominator(bundle,monkeypatch):
    out,plan=bundle;calls=[]
    def broken(*args):calls.append(True);raise RuntimeError('Injected R failure')
    monkeypatch.setattr(reference,'_r_family',broken)
    result=reference.execute_reference_bundle(out,plan,'fixture',ROOT,RSCRIPT,lambda:None)
    assert len(calls)==2 and result['required_slots']==40 and result['passed_slots']==0
    assert result['status_counts']=={'reference_failed':2,'required_input_unavailable':2}
    again=reference.execute_reference_bundle(out,plan,'fixture',ROOT,RSCRIPT,lambda:None)
    assert len(calls)==2 and result==again


def test_interrupted_reference_is_closed_without_replay(bundle,monkeypatch):
    out,plan=bundle;execution=out/'execution';execution.mkdir()
    ledger=AttemptLedger(execution/'events.sqlite',plan_sha256='fixture')
    selection=json.loads((out/'selection.json').read_bytes())
    for basis in plan['scope']['lanes'][:2]:
        ledger.append(selection['run_id']+'_'+basis,'R_started',dict(single_attempt=True))
    ledger.close()
    def forbidden(*args):raise AssertionError('An interrupted reference must never be re-executed')
    monkeypatch.setattr(reference,'_r_family',forbidden)
    result=reference.execute_reference_bundle(out,plan,'fixture',ROOT,RSCRIPT,lambda:None)
    assert result['status_counts']=={'interrupted_reference_retained':2,'required_input_unavailable':2}
    assert result['required_slots']==40 and result['passed_slots']==0


def test_reference_denominator_is_checked_against_plan_not_editable_selection(bundle):
    out,plan=bundle;path=out/'selection.json';selection=json.loads(path.read_bytes())
    selection['families'][0]['required_scalar_cells']=5;path.write_bytes(json_bytes(selection)+b'\n')
    manifest_path=out/'manifest.json';manifest=json.loads(manifest_path.read_bytes())
    manifest['selection.json']=hashlib.sha256(path.read_bytes()).hexdigest();manifest_path.write_bytes(json_bytes(manifest)+b'\n')
    with pytest.raises(ValueError,match='denominator changed'):
        reference.execute_reference_bundle(out,plan,'fixture',ROOT,RSCRIPT,lambda:None)


def test_original_response_transport_cannot_be_changed(bundle):
    out,plan=bundle
    path=out/'known_calibration_known_person/family.csv'
    path.write_bytes(path.read_bytes()+b'\n')
    with pytest.raises(ValueError,match='Reference input changed'):
        reference.execute_reference_bundle(out,plan,'fixture',ROOT,RSCRIPT,lambda:None)


def test_missing_generation_keeps_every_required_reference_slot(tmp_path):
    plan=make_bias_plan();plan.update(classification='DEVELOPMENT_ONLY',study_id='DEVELOPMENT_ONLY_REFERENCE_NO_DATA',master_seed=12350)
    scenario=plan['scenarios'][0];scenario.update(replicates=1,independent_reference_replicates=[1]);plan['scenarios']=[scenario]
    run_id=scenario['id']+'__0001';ledger=AttemptLedger(tmp_path/'events.sqlite',plan_sha256='fixture')
    ledger.append(run_id,'generation_started',dict(run_id=run_id,scenario_id=scenario['id'],replicate=1,
        seed=replicate_seed(plan['master_seed'],plan['study_id']+':'+scenario['id'],1)))
    ledger.append(run_id,'terminal',dict(status='generation_failed'))
    reference.export_reference(tmp_path/'reference',ledger,scenario,1,plan,'fixture');ledger.close()
    result=reference.execute_reference_bundle(tmp_path/'reference',plan,'fixture',ROOT,RSCRIPT,lambda:None)
    assert result['required_slots']==40 and result['passed_slots']==0
    assert result['status_counts']=={'required_input_unavailable':4}
