"""Retained development R results, never new population or reference runs.

The fixture copies g4_mml_final_review_20261006/development_reference (RSM N24).
Its two finite-R, two continuous and two TAM slots are implementation evidence.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile

import pytest

from validation import mml_population_final_review as review
from validation.mml_population_study import AttemptLedger, json_bytes


def write(path, value):
    path.write_bytes(json_bytes(value)+b'\n')


def hashes(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob('*') if p.is_file()}


@pytest.fixture
def retained(tmp_path, monkeypatch):
    fixture = Path(__file__).parent/'fixtures/mml_final_review_rsm_n24.zip'
    with zipfile.ZipFile(fixture) as archive:
        archive.extractall(tmp_path)
    bundle = tmp_path/'reference'
    plan = json.loads((tmp_path/'plan.json').read_bytes())
    identity = json.loads((bundle/'reference_plan.json').read_bytes())['study_sha256']
    from validation.mml_gpcm_integral import GpcmIntegral
    from validation.mml_rasch_integral import RaschIntegral
    def forbidden(*args, **kwargs):
        pytest.fail('Final review must not execute R or new continuous integration')
    monkeypatch.setattr(subprocess, 'run', forbidden)
    monkeypatch.setattr(RaschIntegral, 'continuous', forbidden)
    monkeypatch.setattr(GpcmIntegral, 'continuous', forbidden)
    return bundle, plan, identity


def check(retained):
    bundle, plan, identity = retained
    return review.review_reference(bundle, plan['scenarios'][0], 1, plan, identity)


def r_root(retained, lane):
    return next((retained[0]/'execution').glob('*_Q641_'+lane+'_R'))


def rewrite_file(root, name, value):
    """Keep the local manifest consistent, to exercise the event/input binding."""
    write(root/name, value)
    execution = json.loads((root/'execution.json').read_bytes())
    execution['files'][name] = hashlib.sha256((root/name).read_bytes()).hexdigest()
    write(root/'execution.json', execution)


def rewrite_event(retained, phase, change):
    """Invent a consistent event chain in a temporary copy, never formal files."""
    bundle, _, identity = retained
    path = bundle/'execution/events.sqlite'
    ledger = AttemptLedger.read_only(path, plan_sha256=identity)
    rows = [(key, phase, ledger.get(key, phase)) for key, phase in
            ledger.connection.execute('SELECT run_id,phase FROM event ORDER BY sequence')]
    ledger.close()
    path.unlink()
    ledger = AttemptLedger(path, plan_sha256=identity)
    changed = False
    for key, event_phase, value in rows:
        if not changed and event_phase == phase:
            value = change(value)
            changed = True
        if value is not None:
            ledger.append(key, event_phase, value)
    assert changed
    head = ledger.verify()
    ledger.close()
    summary = json.loads((bundle/'execution/summary.json').read_bytes())
    summary['ledger'] = head
    write(bundle/'execution/summary.json', summary)


def test_actual_retained_comparisons_are_read_only_and_no_references_rerun(retained):
    before = hashes(retained[0])
    result = check(retained)
    assert result['required_slots'] == result['passed_slots'] == 6
    assert all(s['execution_success'] is True for s in result['slots'])
    assert not result['scientific_inference_ready']
    assert hashes(retained[0]) == before


@pytest.mark.parametrize('lane', ['finite_R', 'continuous'])
@pytest.mark.parametrize('mutation', ['raw_result', 'rehashed_result', 'rehashed_input', 'missing_file_hash'])
def test_r_artifact_corruption_or_rebound_values_are_rejected(retained, lane, mutation):
    root = r_root(retained, lane)
    if mutation == 'raw_result':
        (root/'result.json').write_text('{}')
    elif mutation == 'rehashed_result':
        value = json.loads((root/'result.json').read_bytes())
        value['cases'][0]['id'] = 'OTHER_CASE'
        rewrite_file(root, 'result.json', value)
    elif mutation == 'rehashed_input':
        value = json.loads((root/'input.json').read_bytes())
        value['cases'][0]['point'][0] += 1
        rewrite_file(root, 'input.json', value)
    else:
        value = json.loads((root/'execution.json').read_bytes())
        value['files'].pop('stdout.txt')
        write(root/'execution.json', value)
    with pytest.raises(ValueError):
        check(retained)


@pytest.mark.parametrize('lane', ['finite_R', 'continuous'])
def test_nonzero_r_exit_cannot_be_credited(retained, lane):
    root = r_root(retained, lane)
    value = json.loads((root/'execution.json').read_bytes())
    value['returncode'] = 1
    write(root/'execution.json', value)
    result = check(retained)
    assert result['required_slots'] == 6 and result['passed_slots'] == 5
    assert not result['all_checks_pass']


@pytest.mark.parametrize('lane', ['finite_R', 'continuous'])
def test_completed_review_without_completed_r_result_retains_failed_slot(retained, lane):
    rewrite_event(retained, lane+'_R_finished', lambda value: None)
    result = check(retained)
    assert result['required_slots'] == 6 and result['passed_slots'] == 5
    assert not result['all_checks_pass']


@pytest.mark.parametrize('lane', ['finite_R', 'continuous'])
def test_plausible_passing_flags_do_not_replace_numeric_comparison(retained, lane):
    def change(value):
        field = 'errors' if lane == 'finite_R' else 'integration_error'
        value['value'][field]['nll_per_person' if lane == 'finite_R' else 'nll'] = 42.0
        return value
    rewrite_event(retained, lane+'_review_finished', change)
    with pytest.raises(ValueError, match='R numerical comparison'):
        check(retained)


def test_continuous_review_requires_saved_python_integral_values(retained):
    rewrite_event(retained, 'continuous_python_finished', lambda value: None)
    result = check(retained)
    assert result['required_slots'] == 6 and result['passed_slots'] == 5
    assert not result['all_checks_pass']


def test_artifact_symlink_cannot_escape_bundle(retained, tmp_path):
    root = r_root(retained, 'finite_R')
    outside = tmp_path/'outside_stdout.txt'
    outside.write_bytes((root/'stdout.txt').read_bytes())
    (root/'stdout.txt').unlink()
    (root/'stdout.txt').symlink_to(outside)
    with pytest.raises(ValueError, match='Retained reference artifact changed'):
        check(retained)
