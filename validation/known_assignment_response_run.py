"""Hash-bound resumable execution of known-assignment response screening."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pandas as pd

from validation.estimand_distribution_study import _run_one_attempt
from validation.facets_resilient_pair import (
    dependency_manifest,
    dependency_manifest_digest,
    fit_resilient_pair,
)
from validation.known_assignment_large_design_dp_shard import _sha256
from validation.known_assignment_response_prepare import ATTEMPT_TYPES, INPUT_FILES, STUDY_DIR


PLAN_PATH = ROOT / "validation" / "known_assignment_response_screening_plan_20260811.json"
REGISTRATION_PATH = ROOT / "validation" / "known_assignment_response_execution_registration_20260811.json"
INPUT_IDENTITY_PATH = STUDY_DIR / "input_identity.json"
EXECUTION_IDENTITY_PATH = STUDY_DIR / "execution_identity.json"
FACETS_EXE_DEFAULT = Path(r"C:\Facets\Facets.exe")


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}.{time.time_ns()}")
    _json_dump(temporary, value)
    os.replace(temporary, path)


def build_dependency_manifest() -> dict[str, str]:
    manifest = dependency_manifest()
    additions = {
        "validation/known_assignment_response_run.py": Path(__file__).resolve(),
        "validation/known_assignment_response_screening_plan_20260811.json": PLAN_PATH,
        "validation/known_assignment_response_execution_registration_20260811.json": REGISTRATION_PATH,
        "validation/known_assignment_response_prepare.py": ROOT / "validation" / "known_assignment_response_prepare.py",
        "validation/estimand_distribution_study.py": ROOT / "validation" / "estimand_distribution_study.py",
        "validation/facets_resilient_pair.py": ROOT / "validation" / "facets_resilient_pair.py",
    }
    for name, path in additions.items():
        if not path.is_file():
            raise FileNotFoundError(f"Execution dependency missing: {path}")
        manifest[name] = _sha256(path)
    return dict(sorted(manifest.items()))


def validate_registration() -> dict[str, Any]:
    value = json.loads(REGISTRATION_PATH.read_text(encoding="utf-8"))
    expected = {
        "plan_sha256": _sha256(PLAN_PATH),
        "runner_sha256": _sha256(Path(__file__).resolve()),
        "prepare_runner_sha256": _sha256(ROOT / "validation" / "known_assignment_response_prepare.py"),
        "estimand_runner_sha256": _sha256(ROOT / "validation" / "estimand_distribution_study.py"),
        "resilient_pair_sha256": _sha256(ROOT / "validation" / "facets_resilient_pair.py"),
    }
    for key, expected_hash in expected.items():
        if str(value.get(key, "")).lower() != expected_hash.lower():
            raise ValueError(f"Execution registration mismatch: {key}")
    if not bool(value.get("tests_passed_before_registration", False)):
        raise ValueError("Execution registration does not assert passing tests")
    return value


def validate_input_identity() -> dict[str, Any]:
    identity = json.loads(INPUT_IDENTITY_PATH.read_text(encoding="utf-8"))
    if identity["plan_sha256"] != _sha256(PLAN_PATH):
        raise ValueError("Prepared input plan hash mismatch")
    if not bool(identity.get("all_checks_pass", False)):
        raise ValueError("Prepared input audit did not pass")
    input_dir = STUDY_DIR / "retained_input"
    for name in INPUT_FILES:
        if _sha256(input_dir / name) != identity["retained_input_sha256"][name]:
            raise ValueError(f"Prepared input hash mismatch: {name}")
    return identity


def prepare_execution_identity(facets_exe: Path) -> dict[str, Any]:
    if EXECUTION_IDENTITY_PATH.exists():
        return validate_execution_identity(facets_exe)
    registration = validate_registration()
    input_identity = validate_input_identity()
    facets_exe = facets_exe.resolve()
    if not facets_exe.is_file():
        raise FileNotFoundError(f"FACETS executable missing: {facets_exe}")
    dependencies = build_dependency_manifest()
    retained_manifest_path = STUDY_DIR / "execution_dependency_manifest.json"
    _json_dump(retained_manifest_path, dependencies)
    identity = {
        "schema_version": "known_assignment_response_execution_identity_v1",
        "plan_sha256": _sha256(PLAN_PATH),
        "registration_sha256": _sha256(REGISTRATION_PATH),
        "runner_sha256": _sha256(Path(__file__).resolve()),
        "input_identity_sha256": _sha256(INPUT_IDENTITY_PATH),
        "dependency_manifest_sha256": dependency_manifest_digest(dependencies),
        "facets_executable": str(facets_exe),
        "facets_executable_sha256": _sha256(facets_exe),
        "python": sys.version,
        "platform": platform.platform(),
        "attempts": 120,
        "registration_tests": registration.get("test_command"),
        "input_checks": input_identity["checks"],
    }
    _json_dump(EXECUTION_IDENTITY_PATH, identity)
    return identity


def validate_execution_identity(facets_exe: Path | None = None) -> dict[str, Any]:
    identity = json.loads(EXECUTION_IDENTITY_PATH.read_text(encoding="utf-8"))
    validate_registration()
    validate_input_identity()
    expected = {
        "plan_sha256": _sha256(PLAN_PATH),
        "registration_sha256": _sha256(REGISTRATION_PATH),
        "runner_sha256": _sha256(Path(__file__).resolve()),
        "input_identity_sha256": _sha256(INPUT_IDENTITY_PATH),
    }
    for key, value in expected.items():
        if str(identity.get(key, "")).lower() != value.lower():
            raise ValueError(f"Execution identity mismatch: {key}")
    dependencies = build_dependency_manifest()
    retained = json.loads((STUDY_DIR / "execution_dependency_manifest.json").read_text(encoding="utf-8"))
    if retained != dependencies or identity["dependency_manifest_sha256"] != dependency_manifest_digest(dependencies):
        raise ValueError("Execution dependency manifest changed")
    recorded_facets = Path(identity["facets_executable"])
    if facets_exe is not None and facets_exe.resolve() != recorded_facets.resolve():
        raise ValueError("Requested FACETS executable differs from recorded executable")
    if _sha256(recorded_facets) != identity["facets_executable_sha256"]:
        raise ValueError("FACETS executable hash changed")
    return identity


def select_shard_attempts(attempts: pd.DataFrame, *, shard_index: int, shard_count: int) -> pd.DataFrame:
    if shard_count < 1 or shard_index < 0 or shard_index >= shard_count:
        raise ValueError("Require 0 <= shard-index < shard-count")
    ordinals = pd.to_numeric(attempts["AttemptOrdinal"], errors="raise").astype(int)
    if ordinals.duplicated().any():
        raise ValueError("Attempt ordinals must be unique")
    return attempts.loc[ordinals.mod(shard_count).eq(shard_index)].copy()


def _artifact_root(attempt: pd.Series) -> Path:
    return STUDY_DIR / "work" / f"{int(attempt['AttemptOrdinal']):05d}"


def _completion_path(attempt: pd.Series) -> Path:
    return STUDY_DIR / "attempts" / f"{int(attempt['AttemptOrdinal']):05d}" / "completion.json"


def _artifact_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): _sha256(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _validate_completion(completion: dict[str, Any], attempt: pd.Series, identity: dict[str, Any]) -> None:
    expected = {
        "attempt_id": str(attempt["AttemptId"]),
        "attempt_fingerprint": str(attempt["AttemptFingerprint"]),
        "run_input_sha256": str(attempt["RunInputSHA256"]),
        "runner_sha256": identity["runner_sha256"],
        "dependency_manifest_sha256": identity["dependency_manifest_sha256"],
        "facets_executable_sha256": identity["facets_executable_sha256"],
    }
    for key, value in expected.items():
        if str(completion.get(key, "")) != str(value):
            raise ValueError(f"Completion mismatch for {attempt['AttemptId']}: {key}")
    root = STUDY_DIR / str(completion["artifact_root"])
    for name, expected_hash in completion.get("artifact_sha256", {}).items():
        if not (root / name).is_file() or _sha256(root / name) != expected_hash:
            raise ValueError(f"Completion artifact mismatch: {root / name}")


def _write_standard(
    root: Path,
    runs: pd.DataFrame,
    recovery: pd.DataFrame,
    thresholds: pd.DataFrame,
    constraints: Any,
) -> None:
    runs.to_csv(root / "run_ledger.csv", index=False, lineterminator="\n")
    recovery.to_csv(root / "recovery.csv", index=False, lineterminator="\n")
    thresholds.to_csv(root / "thresholds.csv", index=False, lineterminator="\n")
    pd.DataFrame(constraints).to_csv(root / "constraints.csv", index=False, lineterminator="\n")


def run_shard(
    *,
    facets_exe: Path,
    shard_index: int,
    shard_count: int,
    resume: bool,
    timeout_seconds: float,
) -> dict[str, int]:
    identity = prepare_execution_identity(facets_exe)
    attempts = pd.read_csv(STUDY_DIR / "retained_input" / "attempt_manifest.csv")
    selected = select_shard_attempts(attempts, shard_index=shard_index, shard_count=shard_count)
    manifest = pd.read_csv(STUDY_DIR / "retained_input" / "manifest.csv").set_index("RunId", drop=False)
    ratings_all = pd.read_csv(STUDY_DIR / "retained_input" / "generated_ratings.csv")
    truth_all = pd.read_csv(STUDY_DIR / "retained_input" / "generated_facet_truth.csv")
    threshold_all = pd.read_csv(STUDY_DIR / "retained_input" / "generated_pcm_threshold_truth.csv")
    (STUDY_DIR / "attempts").mkdir(exist_ok=True)
    (STUDY_DIR / "locks").mkdir(exist_ok=True)
    summary = {"assigned": len(selected), "completed_now": 0, "skipped": 0, "failed": 0}
    for sequence, (_, attempt) in enumerate(selected.iterrows(), start=1):
        completion_path = _completion_path(attempt)
        if completion_path.is_file():
            completion = json.loads(completion_path.read_text(encoding="utf-8"))
            _validate_completion(completion, attempt, identity)
            if not resume:
                raise FileExistsError(f"Attempt already complete: {attempt['AttemptId']}")
            summary["skipped"] += 1
            print(f"[{sequence}/{len(selected)}] skip {attempt['AttemptId']}", flush=True)
            continue
        artifact_root = _artifact_root(attempt)
        if artifact_root.exists():
            raise FileExistsError(f"Unmarked artifact root requires adjudication: {artifact_root}")
        completion_path.parent.mkdir(parents=True, exist_ok=True)
        run_id = str(attempt["RunId"])
        manifest_row = manifest.loc[run_id]
        ratings = ratings_all.loc[ratings_all["RunId"].eq(run_id)].drop(columns="RunId")
        truth = truth_all.loc[truth_all["RunId"].eq(run_id)]
        thresholds = threshold_all.loc[threshold_all["RunId"].eq(run_id)]
        execution_completed = False
        statistical_ready = False
        calibration_ready: bool | None = None
        failure_reason = ""
        try:
            if str(attempt["AttemptType"]) == "RESILIENT_FACETS_PYTHON_JMLE_PCM":
                outcome = fit_resilient_pair(
                    manifest_row,
                    ratings,
                    truth,
                    thresholds,
                    facets_exe=Path(identity["facets_executable"]),
                    output_dir=artifact_root,
                    lock_path=STUDY_DIR / "locks" / "facets_global.lock",
                    timeout_seconds=timeout_seconds,
                )
                runs = pd.read_csv(artifact_root / "combined_run_ledger.csv")
                recovery = pd.read_csv(artifact_root / "combined_recovery.csv")
                threshold_rows = pd.read_csv(artifact_root / "combined_thresholds.csv")
                _write_standard(artifact_root, runs, recovery, threshold_rows, [])
                statistical_ready = bool(outcome["statistical_evidence_ready"])
                calibration_ready = bool(outcome["calibration_ready"])
            else:
                artifact_root.mkdir(parents=True)
                runs, recovery, threshold_rows, constraints = _run_one_attempt(
                    attempt,
                    manifest_row,
                    ratings,
                    truth,
                    thresholds,
                    facets_exe=Path(identity["facets_executable"]),
                    stage_dir=artifact_root,
                    timeout_seconds=timeout_seconds,
                )
                _write_standard(artifact_root, runs, recovery, threshold_rows, constraints)
                statistical_ready = bool(runs["IncludedInStudy"].fillna(False).astype(bool).all())
            execution_completed = True
        except Exception as exc:  # denominator is retained; no silent replacement
            failure_reason = f"{type(exc).__name__}: {exc}"
            artifact_root.mkdir(parents=True, exist_ok=True)
            _json_dump(artifact_root / "unhandled_failure.json", {"attempt_id": str(attempt["AttemptId"]), "failure_reason": failure_reason})
            summary["failed"] += 1
        completion = {
            "schema_version": "known_assignment_response_attempt_completion_v1",
            "attempt_id": str(attempt["AttemptId"]),
            "attempt_type": str(attempt["AttemptType"]),
            "attempt_fingerprint": str(attempt["AttemptFingerprint"]),
            "run_input_sha256": str(attempt["RunInputSHA256"]),
            "runner_sha256": identity["runner_sha256"],
            "dependency_manifest_sha256": identity["dependency_manifest_sha256"],
            "facets_executable_sha256": identity["facets_executable_sha256"],
            "execution_completed": execution_completed,
            "statistical_evidence_ready": statistical_ready,
            "facets_calibration_ready": calibration_ready,
            "failure_reason": failure_reason,
            "artifact_root": artifact_root.relative_to(STUDY_DIR).as_posix(),
            "artifact_sha256": _artifact_hashes(artifact_root),
        }
        _atomic_json(completion_path, completion)
        summary["completed_now"] += 1
        print(
            f"[{sequence}/{len(selected)}] complete {attempt['AttemptId']} "
            f"ready={statistical_ready} calibration={calibration_ready}",
            flush=True,
        )
    print(json.dumps(summary, sort_keys=True), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--facets-exe", type=Path, default=FACETS_EXE_DEFAULT)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--timeout-seconds", type=float, default=120.0)
    args = parser.parse_args()
    run_shard(
        facets_exe=args.facets_exe,
        shard_index=args.shard_index,
        shard_count=args.shard_count,
        resume=args.resume,
        timeout_seconds=args.timeout_seconds,
    )


if __name__ == "__main__":
    main()
