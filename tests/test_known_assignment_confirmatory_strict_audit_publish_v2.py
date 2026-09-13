from __future__ import annotations

import os
from pathlib import Path

import pytest

from validation import known_assignment_confirmatory_strict_audit as v1
from validation import known_assignment_confirmatory_strict_audit_publish_v2 as publish_v2


def _frozen_v1() -> dict[str, str]:
    return v1._load_json(publish_v2.PLAN_PATH)["frozen_v1"]


def test_identity_uses_registered_frozen_v1_hashes_not_live_rehash() -> None:
    frozen = {
        "audit_plan_sha256": "1" * 64,
        "auditor_sha256": "2" * 64,
        "test_contract_sha256": "3" * 64,
        "execution_registration_sha256": "4" * 64,
    }
    identity = publish_v2._build_v1_identity(
        {"parent_hashes": {"study": "5" * 64}},
        b"assessment",
        b"markdown",
        frozen,
    )
    assert identity["audit_plan_sha256"] == frozen["audit_plan_sha256"]
    assert identity["auditor_sha256"] == frozen["auditor_sha256"]
    assert identity["test_contract_sha256"] == frozen["test_contract_sha256"]
    assert identity["audit_registration_sha256"] == frozen[
        "execution_registration_sha256"
    ]


@pytest.fixture(scope="module")
def audited_result() -> dict[str, object]:
    return v1.audit_study(v1.DEFAULT_STUDY_DIR)


@pytest.mark.retained_evidence
def test_file_atomic_transport_publishes_exact_v1_bundle(
    tmp_path: Path,
    audited_result: dict[str, object],
) -> None:
    output = tmp_path / "published"
    publish_v2.stage_and_publish(audited_result, output, _frozen_v1())
    assert {path.name for path in output.iterdir()} == set(publish_v2.FINAL_NAMES)
    v1.validate_published_bundle(output, audited_result)
    with pytest.raises(FileExistsError, match="overwrite"):
        publish_v2.stage_and_publish(audited_result, output, _frozen_v1())


@pytest.mark.retained_evidence
def test_file_atomic_transport_failure_removes_completion_and_owned_directory(
    tmp_path: Path,
    audited_result: dict[str, object],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output = tmp_path / "failed-publication"
    original_replace = os.replace
    calls = 0

    def fail_second_replace(source: object, destination: object) -> None:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise PermissionError("injected second-file replace failure")
        original_replace(source, destination)

    monkeypatch.setattr(publish_v2.os, "replace", fail_second_replace)
    with pytest.raises(PermissionError, match="injected second-file"):
        publish_v2.stage_and_publish(audited_result, output, _frozen_v1())
    assert not output.exists()


@pytest.mark.skipif(
    not publish_v2.REGISTRATION_PATH.is_file(),
    reason="Publish-v2 registration is created only after transport tests pass",
)
def test_publish_v2_registration_freezes_transport_and_v1() -> None:
    plan = publish_v2.validate_registration()
    assert plan["scope"] == "PUBLICATION_TRANSPORT_ONLY"
    assert plan["scientific_or_audit_logic_change"] is False
    assert plan["transport_contract"]["completion_marker"] == (
        "audit_identity.json is installed last; absence means incomplete publication"
    )
