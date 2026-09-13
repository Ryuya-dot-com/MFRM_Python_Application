#!/usr/bin/env python3
"""Dropbox-safe, file-atomic publisher for the frozen kac200 strict audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any, Iterable

from validation import known_assignment_confirmatory_strict_audit as v1


ROOT = Path(__file__).resolve().parents[1]
PLAN_PATH = (
    ROOT
    / "validation"
    / "known_assignment_confirmatory_strict_audit_publish_v2_plan_20260811.json"
)
REGISTRATION_PATH = (
    ROOT
    / "validation"
    / "known_assignment_confirmatory_strict_audit_publish_v2_execution_registration_20260811.json"
)
TEST_PATH = ROOT / "tests" / "test_known_assignment_confirmatory_strict_audit_publish_v2.py"
FINAL_NAMES = ("assessment.json", "AUDIT_ADDENDUM.md", "audit_identity.json")
REGISTRATION_KEYS = {
    "schema_version",
    "registered_date",
    "registered_before_v2_publication",
    "registered_after_two_v1_transport_failures",
    "scientific_or_audit_logic_change",
    "tests_passed_before_v2_publication",
    "test_command",
    "test_result",
    "plan_sha256",
    "publisher_sha256",
    "test_sha256",
    "frozen_v1",
    "output",
    "claim_boundary",
}


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")


def _durable_create(path: Path, payload: bytes) -> None:
    with path.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def validate_registration() -> dict[str, Any]:
    plan = v1._load_json(PLAN_PATH)
    registration = v1._load_json(REGISTRATION_PATH)
    if plan.get("schema_version") != "known_assignment_confirmatory_strict_audit_publish_v2_plan_v1":
        raise ValueError("Publish-v2 plan schema changed")
    if plan.get("scope") != "PUBLICATION_TRANSPORT_ONLY":
        raise ValueError("Publish-v2 plan scope changed")
    if v1._require_bool(
        plan.get("scientific_or_audit_logic_change"),
        "publish-v2 plan scientific_or_audit_logic_change",
    ):
        raise ValueError("Publish-v2 plan may not change audit/scientific logic")
    v1._require_exact_keys(registration, REGISTRATION_KEYS, "publish-v2 registration")
    if registration["schema_version"] != "known_assignment_confirmatory_strict_audit_publish_v2_registration_v1":
        raise ValueError("Publish-v2 registration schema changed")
    expected_hashes = {
        "plan_sha256": v1.sha256_file(PLAN_PATH),
        "publisher_sha256": v1.sha256_file(Path(__file__).resolve()),
        "test_sha256": v1.sha256_file(TEST_PATH),
    }
    for key, expected in expected_hashes.items():
        if registration[key] != expected:
            raise ValueError(f"Publish-v2 registration hash mismatch: {key}")
    for key in (
        "registered_before_v2_publication",
        "registered_after_two_v1_transport_failures",
        "tests_passed_before_v2_publication",
    ):
        if not v1._require_bool(registration[key], f"publish-v2 {key}"):
            raise ValueError(f"Publish-v2 registration gate is false: {key}")
    if v1._require_bool(
        registration["scientific_or_audit_logic_change"],
        "publish-v2 scientific_or_audit_logic_change",
    ):
        raise ValueError("Publish-v2 may not change audit/scientific logic")
    if registration["frozen_v1"] != plan["frozen_v1"]:
        raise ValueError("Publish-v2 frozen-v1 identity changed")
    if registration["output"] != plan["output"]:
        raise ValueError("Publish-v2 output changed")
    if registration["claim_boundary"] != plan["claim_boundary"]:
        raise ValueError("Publish-v2 claim boundary changed")
    frozen_paths = {
        "audit_plan_sha256": v1.PLAN_PATH,
        "auditor_sha256": Path(v1.__file__).resolve(),
        "test_contract_sha256": v1.TEST_PATH,
        "execution_registration_sha256": v1.REGISTRATION_PATH,
    }
    for key, path in frozen_paths.items():
        if v1.sha256_file(path) != plan["frozen_v1"][key]:
            raise ValueError(f"Frozen v1 changed before publish-v2: {key}")
    return plan


def _build_v1_identity(
    result: dict[str, Any],
    assessment: bytes,
    markdown: bytes,
    frozen_v1: dict[str, str],
) -> dict[str, Any]:
    return {
        "schema_version": f"{v1.SCHEMA_VERSION}_identity_v1",
        "parent_hashes": result["parent_hashes"],
        "audit_plan_sha256": frozen_v1["audit_plan_sha256"],
        "audit_registration_sha256": frozen_v1["execution_registration_sha256"],
        "auditor_sha256": frozen_v1["auditor_sha256"],
        "test_contract_sha256": frozen_v1["test_contract_sha256"],
        "artifact_sha256": {
            "assessment.json": _sha256_bytes(assessment),
            "AUDIT_ADDENDUM.md": _sha256_bytes(markdown),
        },
    }


def _cleanup_owned_output(output_dir: Path, names: Iterable[str], error: BaseException) -> None:
    cleanup_errors: list[str] = []
    for name in names:
        path = output_dir / name
        try:
            path.unlink(missing_ok=True)
        except OSError as cleanup_error:
            cleanup_errors.append(f"{path.name}: {cleanup_error}")
    try:
        output_dir.rmdir()
    except OSError as cleanup_error:
        cleanup_errors.append(f"{output_dir.name}: {cleanup_error}")
    if cleanup_errors and hasattr(error, "add_note"):
        error.add_note("Publish-v2 cleanup errors: " + " | ".join(cleanup_errors))


def stage_and_publish(
    result: dict[str, Any],
    output_dir: Path,
    frozen_v1: dict[str, str],
) -> None:
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite publish-v2 output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=False)
    nonce = f"{os.getpid()}.{time.time_ns()}"
    staged = {
        name: f".{name}.stage.{nonce}"
        for name in FINAL_NAMES
    }
    owned_names = [*staged.values(), *FINAL_NAMES]
    try:
        assessment = _json_bytes(result)
        markdown = v1._markdown(result).encode("utf-8")
        identity_value = _build_v1_identity(result, assessment, markdown, frozen_v1)
        identity = _json_bytes(identity_value)
        payloads = {
            "assessment.json": assessment,
            "AUDIT_ADDENDUM.md": markdown,
            "audit_identity.json": identity,
        }
        for final_name in FINAL_NAMES:
            _durable_create(output_dir / staged[final_name], payloads[final_name])
        if v1._load_json(output_dir / staged["assessment.json"]) != result:
            raise ValueError("Staged assessment differs from recomputation")
        if (output_dir / staged["AUDIT_ADDENDUM.md"]).read_bytes() != markdown:
            raise ValueError("Staged Markdown differs from recomputation")
        if v1._load_json(output_dir / staged["audit_identity.json"]) != identity_value:
            raise ValueError("Staged identity differs from recomputation")
        if {path.name for path in output_dir.iterdir()} != set(staged.values()):
            raise ValueError("Unexpected file entered the publish-v2 staging directory")
        for final_name in FINAL_NAMES:
            os.replace(output_dir / staged[final_name], output_dir / final_name)
        v1.validate_published_bundle(output_dir, result)
    except BaseException as error:
        # audit_identity.json is the completion marker and is invalidated first.
        cleanup_order = ["audit_identity.json", *owned_names]
        _cleanup_owned_output(output_dir, cleanup_order, error)
        raise


def publish(output_dir: Path = v1.DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    plan = validate_registration()
    output_dir = output_dir.resolve()
    expected_output = (ROOT / plan["output"]).resolve()
    if output_dir != expected_output:
        raise ValueError("Publish-v2 may write only to its registered output")
    study_dir = v1.DEFAULT_STUDY_DIR.resolve()
    if v1._inside(output_dir, study_dir):
        raise ValueError("Publish-v2 output must remain outside the frozen study")
    result = v1.audit_study(study_dir)
    plan_after_audit = validate_registration()
    if plan_after_audit != plan:
        raise ValueError("Publish-v2 plan changed during the audit")
    stage_and_publish(result, output_dir, plan["frozen_v1"])
    return result


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=v1.DEFAULT_OUTPUT_DIR)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    result = publish(parse_args(argv).output_dir)
    print(result["overall_audit_status"])


if __name__ == "__main__":
    main()
