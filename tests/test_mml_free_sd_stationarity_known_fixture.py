from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tarfile

import pytest

from validation import mml_free_sd_stationarity_known_fixture as fixture


pytestmark = pytest.mark.retained_evidence

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "validation" / "mml_free_sd_stationarity_known_fixture_v2_20260811"
REGISTRY = ROOT / "validation" / "mml_free_sd_stationarity_fixture_registry_20260811.json"
SOURCES = ROOT / "validation" / "source_snapshots" / "mml_known_fixture_v2_20260811.tar.gz"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def historical_source_sha256(name: str) -> str:
    # The retained identity describes its original source, not today's optional API.
    # Missing or corrupted archived bytes fail; no fallback to updated hashes.
    with tarfile.open(SOURCES, "r:gz") as archive:
        source = archive.extractfile(name)
        assert source is not None, name
        return hashlib.sha256(source.read()).hexdigest()


def retained_payload() -> dict:
    if not ARTIFACT.is_dir():
        pytest.skip("append-only known fixture v2 artifact is not available")
    return json.loads((ARTIFACT / "replay.json").read_text(encoding="utf-8"))


def test_retained_v2_fixture_is_complete_hash_bound_and_nonqualifying() -> None:
    payload = retained_payload()
    assert sorted(path.name for path in ARTIFACT.iterdir()) == [
        "README.md",
        "identity.json",
        "replay.json",
    ]
    identity = json.loads((ARTIFACT / "identity.json").read_text(encoding="utf-8"))
    assert identity["schema_version"] == "mml-free-sd-stationarity-known-fixture-v2"
    assert payload["qualification_eligible"] is False
    assert payload["classification"] == fixture.CLASSIFICATION
    for name, expected in identity["artifact_sha256"].items():
        assert sha256_file(ARTIFACT / name) == expected
    for name, expected in identity["source_sha256"].items():
        assert historical_source_sha256(name) == expected
    assert fixture.validate_frozen_inputs(fixture.DEFAULT_STUDY) == identity[
        "frozen_evidence_sha256"
    ]


def test_historical_fixture_registry_sources_remain_hash_bound() -> None:
    registry = json.loads(REGISTRY.read_text(encoding="utf-8"))
    assert registry["scientific_thresholds_frozen"] is False
    assert registry["qualification_data_generated"] is False
    for item in registry["fixtures"]:
        assert item["qualification_eligible"] is False
        if "source_sha256" in item:
            assert historical_source_sha256(item["source"]) == item["source_sha256"]
        if "implementation_sha256" in item:
            assert (
                historical_source_sha256(item["implementation"])
                == item["implementation_sha256"]
            )
        if item.get("status") == "CURRENT_HARDENED_ENGINEERING_REPLAY":
            assert sha256_file(ARTIFACT / "identity.json") == item[
                "artifact_identity_sha256"
            ]


def test_file_atomic_publication_validates_before_identity_last(
    tmp_path: Path,
    monkeypatch,
) -> None:
    payload = retained_payload()
    monkeypatch.setattr(fixture, "replay", lambda _study: payload)
    output = tmp_path / "known_fixture_copy"
    fixture.publish(output, fixture.DEFAULT_STUDY)
    assert sorted(path.name for path in output.iterdir()) == [
        "README.md",
        "identity.json",
        "replay.json",
    ]
    identity = json.loads((output / "identity.json").read_text(encoding="utf-8"))
    for name, expected in identity["artifact_sha256"].items():
        assert sha256_file(output / name) == expected


def test_publication_failure_removes_owned_partial_directory(
    tmp_path: Path,
    monkeypatch,
) -> None:
    payload = retained_payload()
    monkeypatch.setattr(fixture, "replay", lambda _study: payload)
    real_replace = fixture.os.replace
    calls = 0

    def fail_second_replace(source, target):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected Dropbox file lock")
        return real_replace(source, target)

    monkeypatch.setattr(fixture.os, "replace", fail_second_replace)
    output = tmp_path / "known_fixture_failure"
    with pytest.raises(OSError, match="injected Dropbox file lock"):
        fixture.publish(output, fixture.DEFAULT_STUDY)
    assert not output.exists()
