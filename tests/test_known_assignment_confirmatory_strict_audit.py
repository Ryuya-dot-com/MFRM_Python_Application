from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from validation import known_assignment_confirmatory_strict_audit as audit


def test_json_loader_rejects_duplicate_keys_and_nonfinite_constants(tmp_path: Path) -> None:
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text('{"gate": true, "gate": false}\n', encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate JSON key"):
        audit._load_json(duplicate)

    nonfinite = tmp_path / "nonfinite.json"
    nonfinite.write_text('{"metric": NaN}\n', encoding="utf-8")
    with pytest.raises(ValueError, match="Nonfinite JSON constant"):
        audit._load_json(nonfinite)


@pytest.mark.parametrize("value", [1, 0, "true", "false", None])
def test_json_boolean_contract_rejects_truthy_substitutes(value: object) -> None:
    with pytest.raises(TypeError, match="JSON boolean"):
        audit._require_bool(value, "gate")


def test_csv_boolean_and_finite_contracts_are_fail_closed() -> None:
    assert audit._strict_bool_series(pd.Series(["True", "False"]), "gate").tolist() == [
        True,
        False,
    ]
    with pytest.raises(TypeError, match="noncanonical boolean"):
        audit._strict_bool_series(pd.Series(["TRUE"]), "gate")
    with pytest.raises(TypeError, match="noncanonical boolean"):
        audit._strict_bool_series(pd.Series([True]), "gate")
    with pytest.raises(ValueError, match="nonfinite numeric"):
        audit._require_finite(pd.DataFrame({"value": [np.nan]}), ("value",), "metric")


def test_csv_reader_preserves_boolean_lexemes_and_rejects_bad_headers(
    tmp_path: Path,
) -> None:
    lexical = tmp_path / "lexical.csv"
    lexical.write_text("gate\nTRUE\ntrue\nFALSE\n", encoding="utf-8")
    frame = audit._read_csv(lexical)
    assert frame["gate"].tolist() == ["TRUE", "true", "FALSE"]
    with pytest.raises(TypeError, match="noncanonical boolean"):
        audit._strict_bool_series(frame["gate"], "gate")

    duplicate = tmp_path / "duplicate-header.csv"
    duplicate.write_text("value,value\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="header"):
        audit._read_csv(duplicate)

    extra = tmp_path / "extra-header.csv"
    extra.write_text("value,extra\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="columns/order"):
        audit._read_csv(extra, ("value",))


@pytest.mark.parametrize("value", [True, 1.0, "1", None])
def test_json_integer_contract_rejects_bool_float_and_text(value: object) -> None:
    with pytest.raises(TypeError, match="JSON integer"):
        audit._require_json_int(value, "count")


def _valid_artifact_fixture(study: Path) -> dict[str, object]:
    artifact_root = study / "work" / "00001"
    artifact_root.mkdir(parents=True)
    for filename in audit.STANDARD_ATTEMPT_FILES:
        (artifact_root / filename).write_text(filename, encoding="utf-8")
    return {
        "artifact_root": "work/00001",
        "artifact_sha256": audit._file_hashes(artifact_root),
    }


def test_attempt_artifact_contract_rejects_empty_escape_and_unrecorded_file(
    tmp_path: Path,
) -> None:
    study = tmp_path / "study"
    completion = _valid_artifact_fixture(study)
    assert len(audit._validate_attempt_artifacts(study, 1, completion)) == 4

    empty = dict(completion, artifact_sha256={})
    with pytest.raises(ValueError, match="hash map is empty"):
        audit._validate_attempt_artifacts(study, 1, empty)

    escaping = dict(completion, artifact_root="../outside")
    with pytest.raises(ValueError, match="not canonical"):
        audit._validate_attempt_artifacts(study, 1, escaping)

    (study / "work" / "00001" / "unrecorded.txt").write_text("x", encoding="utf-8")
    with pytest.raises(ValueError, match="artifact set/hash mismatch"):
        audit._validate_attempt_artifacts(study, 1, completion)


def test_artifact_hashing_rejects_symlink_when_supported(tmp_path: Path) -> None:
    root = tmp_path / "root"
    root.mkdir()
    target = tmp_path / "target.txt"
    target.write_text("target", encoding="utf-8")
    link = root / "link.txt"
    try:
        os.symlink(target, link)
    except (OSError, NotImplementedError):
        pytest.skip("Symlink creation is unavailable in this Windows test context")
    with pytest.raises(ValueError, match="symlink"):
        audit._file_hashes(root)


def test_exact_group_key_contract_rejects_duplicate_key() -> None:
    frame = pd.DataFrame(
        {
            "RunId": ["run-1", "run-1"],
            "Coordinate": ["x", "x"],
        }
    )
    with pytest.raises(ValueError, match="duplicate keys"):
        audit._require_exact_group_keys(
            frame,
            "RunId",
            ("Coordinate",),
            {"run-1"},
            {"x"},
            "gradient",
        )


@pytest.mark.parametrize(
    "observed",
    [
        {},
        {"../escape.csv": "0" * 64},
        {**{name: "0" * 64 for name in audit.R_RAW_SOURCE_FILES}, "unknown.csv": "0" * 64},
    ],
)
def test_r_source_hash_map_contract_is_exact(observed: dict[str, str]) -> None:
    with pytest.raises(ValueError, match="keys differ"):
        audit._require_exact_keys(observed, audit.R_RAW_SOURCE_FILES, "R raw source hash map")


@pytest.mark.retained_evidence
def test_facets_pair_metrics_reject_fractional_tries_and_bad_replay() -> None:
    study = audit.DEFAULT_STUDY_DIR
    work = study / "work" / "00048"
    metrics = audit._load_json(work / "resilient_pair_metrics.json")
    identity = audit._load_json(study / "study_identity.json")
    kwargs = {
        "ordinal": 48,
        "expected_run_id": "gamma_neg_0p8__vector-05",
        "identity": identity,
        "work": work,
        "maximum_tries": 3,
        "replay_tolerance": 1e-12,
    }
    assert audit._validate_facets_pair_metrics(metrics, **kwargs) == 1

    fractional = dict(metrics, facets_tries=1.5)
    with pytest.raises(TypeError, match="JSON integer"):
        audit._validate_facets_pair_metrics(fractional, **kwargs)

    too_many = dict(metrics, facets_tries=4, facets_retry_count=3)
    with pytest.raises(ValueError, match="retry ledger"):
        audit._validate_facets_pair_metrics(too_many, **kwargs)

    nonfinite = dict(metrics, python_replay_main_max_abs_difference=float("nan"))
    with pytest.raises(ValueError, match="must be finite"):
        audit._validate_facets_pair_metrics(nonfinite, **kwargs)

    excessive = dict(metrics, python_replay_threshold_max_abs_difference=2e-12)
    with pytest.raises(ValueError, match="replay contract"):
        audit._validate_facets_pair_metrics(excessive, **kwargs)


def test_publish_rejects_unplanned_output_before_running_audit(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="planned output path"):
        audit.publish_audit(audit.DEFAULT_STUDY_DIR, tmp_path / "wrong-output")


def test_publish_write_failure_preserves_original_error_and_cleans_partial(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    study = tmp_path / "study"
    study.mkdir()
    output = tmp_path / "audit-output"
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(
        json.dumps({"immutability": {"output": str(output)}}), encoding="utf-8"
    )
    monkeypatch.setattr(audit, "PLAN_PATH", plan_path)
    monkeypatch.setattr(audit, "audit_study", lambda _: {"result": "unused"})

    def fail_write(*_args: object, **_kwargs: object) -> None:
        raise RuntimeError("injected write failure")

    monkeypatch.setattr(audit, "_write_json", fail_write)
    with pytest.raises(RuntimeError, match="injected write failure"):
        audit.publish_audit(study, output)
    assert not output.exists()
    assert not list(tmp_path.glob("audit-output.partial.*"))


@pytest.mark.skipif(
    not audit.REGISTRATION_PATH.is_file(),
    reason="Strict-audit execution registration is created only after unit tests pass",
)
@pytest.mark.retained_evidence
def test_registered_kac200_bundle_passes_strict_audit_without_publication(
    tmp_path: Path,
) -> None:
    result = audit.audit_study(audit.DEFAULT_STUDY_DIR)
    assert result["overall_audit_status"] == (
        "INTEGRITY_PASS_SCIENTIFIC_WORDING_QUALIFIED"
    )
    assert result["strict_evidence_integrity"]["status"] == "PASS"
    qualification = result["posthoc_scientific_qualification"]
    assert qualification["OptimizationStationarityQualification"] == "INCOMPLETE"
    assert qualification["exact_mle_stationarity_claim_authorized"] is False
    sensitivity = result["r_subset_descriptive_sensitivity"]
    assert sensitivity["role"] == "POSTHOC_DESCRIPTIVE_NUMERICAL_SENSITIVITY"
    assert sensitivity["new_confirmatory_test_constructed"] is False
    assert result["registered_endpoints"]["all_four_endpoint_statistics_recomputed"] is True
    assert result["registered_endpoints"]["secondary_holm_recomputed"] is True

    bundle = tmp_path / "published-bundle"
    bundle.mkdir()
    audit._write_json(bundle / "assessment.json", result)
    (bundle / "AUDIT_ADDENDUM.md").write_text(
        audit._markdown(result), encoding="utf-8"
    )
    identity = {
        "schema_version": f"{audit.SCHEMA_VERSION}_identity_v1",
        "parent_hashes": result["parent_hashes"],
        "audit_plan_sha256": audit.sha256_file(audit.PLAN_PATH),
        "audit_registration_sha256": audit.sha256_file(audit.REGISTRATION_PATH),
        "auditor_sha256": audit.sha256_file(Path(audit.__file__).resolve()),
        "test_contract_sha256": audit.sha256_file(audit.TEST_PATH),
        "artifact_sha256": {
            filename: audit.sha256_file(bundle / filename)
            for filename in ("assessment.json", "AUDIT_ADDENDUM.md")
        },
    }
    audit._write_json(bundle / "audit_identity.json", identity)
    audit.validate_published_bundle(bundle, result)

    (bundle / "AUDIT_ADDENDUM.md").write_text("tampered\n", encoding="utf-8")
    identity["artifact_sha256"]["AUDIT_ADDENDUM.md"] = audit.sha256_file(
        bundle / "AUDIT_ADDENDUM.md"
    )
    audit._write_json(bundle / "audit_identity.json", identity)
    with pytest.raises(ValueError, match="Markdown differs"):
        audit.validate_published_bundle(bundle, result)


def test_plan_does_not_authorize_posthoc_stability_or_exact_mle_claim() -> None:
    plan = audit._load_json(audit.PLAN_PATH)
    scientific = plan["scientific_qualification_contract"]
    sensitivity = plan["r_subset_sensitivity_contract"]
    assert scientific["optimization_stationarity_qualification"] == "INCOMPLETE"
    assert scientific["exact_mle_stationarity_claim_authorized"] is False
    assert sensitivity["role"] == "POSTHOC_DESCRIPTIVE_NUMERICAL_SENSITIVITY"
    assert sensitivity["new_confirmatory_test_constructed"] is False
    assert "stable" not in scientific["authorized_claim"].lower()
