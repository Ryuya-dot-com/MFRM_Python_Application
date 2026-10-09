"""The study cannot bypass its signed receipt or frozen input inventory."""
import json
import subprocess

import pytest

from validation import mml_population_controller as controller
from validation.mml_population_study import json_bytes
from tests.test_mml_population_study import fixture_plan


def fixture_seal(tmp_path,monkeypatch):
    plan,_,_=fixture_plan();monkeypatch.setattr(controller,'make_plan',lambda:plan)
    monkeypatch.setattr(controller,'python_runtime',lambda:{'fixture':True})
    monkeypatch.setattr(controller,'r_runtime',lambda:{'fixture':True})
    out=tmp_path/'study';controller.freeze(out)
    return out


def test_frozen_sources_include_import_time_citation_asset():
    assert controller.ROOT/'CITATION.cff' in controller.source_paths()


def test_absent_or_wrong_receipt_cannot_authorize_a_frozen_plan(tmp_path,monkeypatch):
    out=fixture_seal(tmp_path,monkeypatch)
    assert controller.verify_seal(out,receipt=False)[0]['study_id'].startswith('DEVELOPMENT_')
    with pytest.raises(subprocess.CalledProcessError):controller.verify_seal(out)
    probe=controller.ROOT/'validation/generated/g4_timestamp_probe_20261004/response.tsr'
    (out/'response.tsr').write_bytes(probe.read_bytes())
    with pytest.raises(subprocess.CalledProcessError):controller.verify_seal(out)
    assert not (out/'attempts').exists()


def test_source_bundle_change_fails_before_receipt_verification(tmp_path,monkeypatch):
    out=fixture_seal(tmp_path,monkeypatch)
    with (out/'sources.zip').open('ab') as handle:handle.write(b'changed')
    with pytest.raises(ValueError,match='Frozen artifact changed'):controller.verify_seal(out,receipt=False)


def test_ordered_seed_inventory_is_checked_semantically_even_with_updated_local_hash(tmp_path,monkeypatch):
    out=fixture_seal(tmp_path,monkeypatch)
    lines=(out/'runs.csv').read_text().splitlines();cells=lines[1].split(',');cells[3]='123';lines[1]=','.join(cells)
    (out/'runs.csv').write_text('\n'.join(lines)+'\n')
    manifest=json.loads((out/'T0.json').read_bytes());manifest['files']['runs.csv']=controller.digest(out/'runs.csv')
    (out/'T0.json').write_bytes(json_bytes(manifest))
    with pytest.raises(ValueError,match='Ordered RunId/seed'):controller.verify_seal(out,receipt=False)
