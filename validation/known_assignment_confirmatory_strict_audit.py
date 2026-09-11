#!/usr/bin/env python3
"""Append-only, fail-closed audit of the completed kac200 confirmation.

This auditor is intentionally separate from the hash-frozen v1 runner.  It
does not repair, replace, or reinterpret any v1 artifact.  It verifies the
retained evidence bundle and records a posthoc qualification of the wording
around Python optimizer stationarity.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import time
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
PLAN_PATH = ROOT / "validation" / "known_assignment_confirmatory_strict_audit_plan_20260811.json"
REGISTRATION_PATH = (
    ROOT
    / "validation"
    / "known_assignment_confirmatory_strict_audit_execution_registration_20260811.json"
)
TEST_PATH = ROOT / "tests" / "test_known_assignment_confirmatory_strict_audit.py"
DEFAULT_STUDY_DIR = ROOT / "validation" / "kac200_20260811"
DEFAULT_OUTPUT_DIR = ROOT / "validation" / "kac200_posthoc_strict_audit_20260811"

SCHEMA_VERSION = "kac200_posthoc_strict_audit_v1"
AUDIT_REGISTRATION_KEYS = {
    "schema_version",
    "registered_date",
    "registration_nature",
    "registered_before_end_to_end_audit_and_publication",
    "registered_after_confirmatory_endpoints_known",
    "original_v1_artifacts_modified",
    "scientific_acceptance_thresholds_selected_from_kac200",
    "posthoc_integrity_tolerances_labeled_posthoc",
    "tests_passed_before_audit",
    "strict_test_command_before_registration",
    "strict_test_result_before_registration",
    "existing_validation_test_command",
    "existing_validation_test_result",
    "audit_plan_sha256",
    "auditor_sha256",
    "test_contract_sha256",
    "target",
    "environment",
    "r_process_contract",
    "claim_boundary",
}
PAIR_ATTEMPT = "RESILIENT_FACETS_PYTHON_JMLE_PCM"
MML_FIXED_ATTEMPT = "PYTHON_MML_FIXED_SD08_Q31_PCM"
MML_FREE_ATTEMPT = "PYTHON_MML_FREE_SD_Q31_PCM"
CMLE_ATTEMPT = "PYTHON_EXACT_CMLE_PCM"
MML_FIXED_MODE = "PYTHON_MML_FIXED_SD08_Q31"
MML_FREE_MODE = "PYTHON_MML_FREE_SD_Q31"
FACETS_MODE = "FACETS_4_5_JMLE"
PRIMARY_ID = "KA1_FREE_MML_EXPECTED_WITHIN_DATASET_FIXED_RATER_RMSE_STRESS_CONTRAST"
ENDPOINT_ALTERNATIVES = {
    PRIMARY_ID: "greater",
    "KA2_FIXED_MML_EXPECTED_WITHIN_DATASET_FIXED_RATER_RMSE_STRESS_CONTRAST": "greater",
    "KA3_FREE_MML_SYMMETRIC_STRESS_POPULATION_SD_SHIFT": "less",
    "KA4_FREE_MML_DIRECTION_ALIGNED_RATER_ERROR_SLOPE": "greater",
}
CALIBRATION_TRIPLETS = (5, 22, 41, 57, 72, 99, 116, 128, 143, 161, 176, 200)
RATER_TRUTH = {"R01": -0.45, "R02": -0.15, "R03": 0.15, "R04": 0.45}
STANDARD_ATTEMPT_FILES = {
    "run_ledger.csv",
    "recovery.csv",
    "thresholds.csv",
    "constraints.csv",
}
ATTEMPT_MANIFEST_COLUMNS = (
    "AttemptOrdinal",
    "AttemptId",
    "RunId",
    "AttemptType",
    "PersonVector",
    "Replicate",
    "Gamma",
    "Design",
    "PersonDistribution",
    "RunInputSHA256",
    "DependencyManifestSHA256",
    "AttemptFingerprint",
)
R_BUNDLE_FILES = {
    "assessment.json",
    "crossfit_cross_evaluations.csv",
    "crossfit_gradients.csv",
    "crossfit_parameters.csv",
    "crossfit_runs.csv",
    "runtime_identity.csv",
    "r_process.json",
    "verification_identity.json",
}
R_RAW_SOURCE_FILES = {
    "crossfit_gradients.csv",
    "crossfit_parameters.csv",
    "crossfit_runs.csv",
    "runtime_identity.csv",
}
R_RETAINED_INPUT_FILES = {
    "attempt_manifest.csv",
    "generated_ratings.csv",
    "manifest.csv",
}
R_CSV_COLUMNS = {
    "crossfit_runs.csv": (
        "RunId", "PersonVector", "Gamma", "PythonRecordedLogLik",
        "PythonSolutionRLogLikQ31", "PythonSolutionRLogLikQ61",
        "PythonSolutionRGradientSupNormQ31", "PythonSigma",
        "PythonSigmaFixedPoint", "PythonSigmaFixedPointResidual",
        "RConvergenceCode", "RMessage", "RLogLikQ31",
        "RQ31SolutionRLogLikQ61", "RGradientSupNormQ31", "RSigma",
        "RFunctionEvaluations", "RGradientEvaluations", "RQ61ConvergenceCode",
        "RQ61Message", "ROptimizedLogLikQ61", "RGradientSupNormQ61",
        "RQ61Sigma", "RQ61FunctionEvaluations", "RQ61GradientEvaluations",
    ),
    "crossfit_parameters.csv": (
        "RunId", "Block", "Level", "PythonEstimate", "REstimateQ31",
        "REstimateQ61", "DifferenceRQ31MinusPython", "DifferenceRQ61MinusRQ31",
    ),
    "crossfit_gradients.csv": (
        "RunId", "Coordinate", "PythonSolutionRGradientQ31",
        "RQ31SolutionRGradientQ31", "RQ61SolutionRGradientQ61",
    ),
    "runtime_identity.csv": (
        "RVersion", "Platform", "BLAS", "LAPACK", "SelectedVectors", "Datasets",
    ),
    "crossfit_cross_evaluations.csv": (
        "RunId", "PythonSolutionPythonLogLikQ31", "PythonSolutionPythonLogLikQ61",
        "RQ31SolutionPythonLogLikQ31", "RQ31SolutionPythonLogLikQ61",
        "RQ61SolutionPythonLogLikQ61", "PersonVector", "Gamma",
        "PythonRecordedLogLik", "PythonSolutionRLogLikQ31",
        "PythonSolutionRLogLikQ61", "PythonSolutionRGradientSupNormQ31",
        "PythonSigma", "PythonSigmaFixedPoint", "PythonSigmaFixedPointResidual",
        "RConvergenceCode", "RMessage", "RLogLikQ31", "RQ31SolutionRLogLikQ61",
        "RGradientSupNormQ31", "RSigma", "RFunctionEvaluations",
        "RGradientEvaluations", "RQ61ConvergenceCode", "RQ61Message",
        "ROptimizedLogLikQ61", "RGradientSupNormQ61", "RQ61Sigma",
        "RQ61FunctionEvaluations", "RQ61GradientEvaluations",
        "AbsRAtPythonVsRecordedQ31", "AbsPythonAtPythonVsRecordedQ31",
        "AbsCrossLanguageAtPythonQ31", "AbsCrossLanguageAtRQ31",
        "AbsCrossLanguageAtRQ61", "RQ31ImprovementOverPythonQ31",
        "PythonQ61MinusQ31", "RQ61OptimizedMinusRQ31EvaluatedQ61",
    ),
}
FACETS_PAIR_METRIC_KEYS = {
    "calibration_ready",
    "confirmatory_evidence_replaced",
    "dependency_manifest_sha256",
    "direct_agreement_pass",
    "facets_executable_sha256",
    "facets_executable_size_bytes",
    "facets_final_error",
    "facets_reported_version",
    "facets_retry_count",
    "facets_tries",
    "lock_scope",
    "pair_fully_qualified",
    "python_replay_main_max_abs_difference",
    "python_replay_threshold_max_abs_difference",
    "python_replay_tolerance",
    "run_id",
    "schema_version",
    "statistical_evidence_ready",
}
COMPLETION_KEYS = {
    "schema_version",
    "attempt_id",
    "attempt_type",
    "attempt_fingerprint",
    "run_input_sha256",
    "dependency_manifest_sha256",
    "runner_sha256",
    "facets_executable_sha256",
    "execution_completed",
    "statistical_evidence_ready",
    "facets_calibration_ready",
    "failure_reason",
    "artifact_root",
    "artifact_sha256",
}
R_GATES = {
    "r_version_4_5_1",
    "datasets_complete",
    "r_q31_convergence_all",
    "r_q61_convergence_all",
    "r_at_python_matches_recorded_q31",
    "python_at_python_matches_recorded_q31",
    "cross_language_python_solution_q31",
    "cross_language_r_solution_q31",
    "cross_language_r_solution_q61",
    "r_q31_terminal_gradient",
    "r_q61_terminal_gradient",
    "python_sigma_fixed_point",
    "r_q31_vs_python_parameter",
    "r_q31_vs_python_sigma",
    "q61_vs_q31_parameter",
    "q61_vs_q31_sigma",
    "r_q31_loglik_improvement_bounded",
}
R_PARAMETER_KEYS = {
    *(('Rater', level) for level in ('R01', 'R02', 'R03', 'R04')),
    *(('Task', level) for level in ('T01', 'T02', 'T03')),
    *(('Criterion', level) for level in ('C01', 'C02')),
    *(('Step', f'{criterion}::{category}') for criterion in ('C01', 'C02') for category in (1, 2, 3)),
}
R_GRADIENT_KEYS = {
    "Rater::R01",
    "Rater::R02",
    "Rater::R03",
    "Task::T01",
    "Task::T02",
    "Criterion::C01",
    "Criterion::C02",
    "Step::C01::1",
    "Step::C01::2",
    "Step::C02::1",
    "Step::C02::2",
    "LogPopulationSD",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Duplicate JSON key: {key}")
        value[key] = item
    return value


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=_reject_duplicate_json_keys,
        parse_constant=lambda constant: (_ for _ in ()).throw(
            ValueError(f"Nonfinite JSON constant: {constant}")
        ),
    )
    if not isinstance(value, dict):
        raise TypeError(f"Expected a JSON object: {path}")
    return value


def _write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _require_exact_keys(value: dict[str, Any], expected: set[str], label: str) -> None:
    observed = set(value)
    if observed != expected:
        raise ValueError(
            f"{label} keys differ; missing={sorted(expected - observed)}, "
            f"extra={sorted(observed - expected)}"
        )


def _require_bool(value: Any, label: str) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"{label} must be a JSON boolean, got {type(value).__name__}")
    return value


def _require_json_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label} must be a JSON integer, got {type(value).__name__}")
    return value


def _require_json_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{label} must be a JSON number, got {type(value).__name__}")
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _strict_bool_series(series: pd.Series, label: str) -> pd.Series:
    def parse(value: Any) -> bool:
        if isinstance(value, str) and value in {"True", "False"}:
            return value == "True"
        raise TypeError(f"{label} contains a noncanonical boolean: {value!r}")

    return series.map(parse)


def _read_csv(
    path: Path,
    expected_columns: tuple[str, ...] | None = None,
) -> pd.DataFrame:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        try:
            header = next(csv.reader(handle))
        except StopIteration as error:
            raise ValueError(f"CSV is empty: {path}") from error
    if not header or any(not column for column in header) or len(header) != len(set(header)):
        raise ValueError(f"CSV header is empty, duplicated, or unnamed: {path}")
    if expected_columns is not None and tuple(header) != expected_columns:
        raise ValueError(f"CSV columns/order changed: {path}")
    return pd.read_csv(
        path,
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    )


def _require_finite(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    columns = tuple(columns)
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing numeric columns: {missing}")
    numeric = frame.loc[:, columns].apply(pd.to_numeric, errors="coerce")
    bad = ~np.isfinite(numeric.to_numpy(dtype=float))
    if bad.any():
        locations = np.argwhere(bad)[:10]
        detail = [
            f"row={int(row)}, column={columns[int(column)]}"
            for row, column in locations
        ]
        raise ValueError(f"{label} contains nonfinite numeric values: {detail}")


def _inside(path: Path, root: Path) -> bool:
    try:
        return os.path.commonpath((str(path.resolve()), str(root.resolve()))) == str(
            root.resolve()
        )
    except ValueError:
        return False


def _is_link_like(path: Path) -> bool:
    return path.is_symlink() or bool(
        hasattr(path, "is_junction") and path.is_junction()
    )


def _guard_study_roots(study_dir: Path) -> None:
    for name in ("retained_input", "attempts", "work", "aggregate", "r_mml"):
        path = study_dir / name
        if not path.is_dir() or _is_link_like(path) or not _inside(path, study_dir):
            raise ValueError(f"Study evidence root is missing, linked, or escaping: {name}")


def _file_hashes(root: Path) -> dict[str, str]:
    if _is_link_like(root):
        raise ValueError(f"Artifact root must not be a symlink or junction: {root}")
    hashes: dict[str, str] = {}
    for path in sorted(root.rglob("*"), key=lambda item: item.as_posix()):
        if _is_link_like(path):
            raise ValueError(f"Artifact tree contains a symlink or junction: {path}")
        if path.is_file():
            hashes[path.relative_to(root).as_posix()] = sha256_file(path)
    return hashes


def _canonical_digest(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _validate_attempt_artifacts(
    study_dir: Path,
    ordinal: int,
    completion: dict[str, Any],
) -> dict[str, str]:
    expected_relative = f"work/{ordinal:05d}"
    if completion["artifact_root"] != expected_relative:
        raise ValueError(f"Attempt {ordinal} artifact root is not canonical")
    artifact_root = (study_dir / completion["artifact_root"]).resolve()
    canonical_root = (study_dir / expected_relative).resolve()
    if artifact_root != canonical_root or not _inside(artifact_root, study_dir):
        raise ValueError(f"Attempt {ordinal} artifact root escapes the study")
    recorded_hashes = completion["artifact_sha256"]
    if not isinstance(recorded_hashes, dict) or not recorded_hashes:
        raise ValueError(f"Attempt {ordinal} artifact hash map is empty or invalid")
    actual_hashes = _file_hashes(artifact_root)
    if recorded_hashes != actual_hashes:
        missing = sorted(set(actual_hashes) - set(recorded_hashes))
        extra = sorted(set(recorded_hashes) - set(actual_hashes))
        changed = sorted(
            key
            for key in set(recorded_hashes) & set(actual_hashes)
            if recorded_hashes[key] != actual_hashes[key]
        )
        raise ValueError(
            f"Attempt {ordinal} artifact set/hash mismatch; "
            f"unrecorded={missing}, missing={extra}, changed={changed}"
        )
    if not STANDARD_ATTEMPT_FILES.issubset(actual_hashes):
        raise ValueError(f"Attempt {ordinal} lacks a standard artifact")
    return actual_hashes


def validate_audit_registration() -> tuple[dict[str, Any], dict[str, Any]]:
    plan = _load_json(PLAN_PATH)
    registration = _load_json(REGISTRATION_PATH)
    _require_exact_keys(
        registration, AUDIT_REGISTRATION_KEYS, "strict-audit execution registration"
    )
    expected = {
        "audit_plan_sha256": sha256_file(PLAN_PATH),
        "auditor_sha256": sha256_file(Path(__file__).resolve()),
        "test_contract_sha256": sha256_file(TEST_PATH),
    }
    for key, digest in expected.items():
        if str(registration.get(key, "")).lower() != digest:
            raise ValueError(f"Strict-audit registration mismatch: {key}")
    if registration.get("schema_version") != "known_assignment_confirmatory_strict_audit_registration_v1":
        raise ValueError("Strict-audit registration schema changed")
    if not _require_bool(registration.get("tests_passed_before_audit"), "tests_passed_before_audit"):
        raise ValueError("Strict-audit tests were not registered as passing")
    required_true = (
        "registered_before_end_to_end_audit_and_publication",
        "registered_after_confirmatory_endpoints_known",
        "posthoc_integrity_tolerances_labeled_posthoc",
    )
    if not all(
        _require_bool(registration[key], f"registration {key}") for key in required_true
    ):
        raise ValueError("Strict-audit timing/disclosure registration is false")
    required_false = (
        "original_v1_artifacts_modified",
        "scientific_acceptance_thresholds_selected_from_kac200",
    )
    if any(
        _require_bool(registration[key], f"registration {key}") for key in required_false
    ):
        raise ValueError("Strict-audit immutability/scientific-threshold disclosure changed")
    if registration["target"] != plan["target"]:
        raise ValueError("Strict-audit registration target differs from the plan")
    if registration["claim_boundary"] != plan["claim_boundary"]:
        raise ValueError("Strict-audit registration claim boundary differs from the plan")
    return plan, registration


def validate_frozen_identity(
    study_dir: Path,
    plan: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    target = plan["target"]
    identity_path = study_dir / "study_identity.json"
    if sha256_file(identity_path) != target["study_identity_sha256"]:
        raise ValueError("Target study identity hash changed")
    if sha256_file(study_dir / "aggregate" / "aggregate_identity.json") != target[
        "aggregate_identity_sha256"
    ]:
        raise ValueError("Target aggregate identity hash changed")
    if sha256_file(study_dir / "r_mml" / "verification_identity.json") != target[
        "r_verification_identity_sha256"
    ]:
        raise ValueError("Target R verification identity hash changed")
    if sha256_file(study_dir / "r_mml" / "assessment.json") != target[
        "r_assessment_sha256"
    ]:
        raise ValueError("Target R assessment hash changed")
    if sha256_file(study_dir / "r_mml" / "r_process.json") != target[
        "r_process_sha256"
    ]:
        raise ValueError("Target R process hash changed")
    if sha256_file(study_dir / "aggregate" / "assessment.json") != target[
        "aggregate_assessment_sha256"
    ]:
        raise ValueError("Target aggregate assessment hash changed")

    identity = _load_json(identity_path)
    dependency_manifest = _load_json(study_dir / "dependency_manifest.json")
    for relative, digest in dependency_manifest.items():
        path = ROOT / relative
        if not path.is_file() or sha256_file(path) != digest:
            raise ValueError(f"Frozen execution dependency changed: {relative}")
    if _canonical_digest(dependency_manifest) != identity["dependency_manifest_sha256"]:
        raise ValueError("Dependency manifest digest changed")
    for filename, digest in identity["retained_input_sha256"].items():
        path = study_dir / "retained_input" / filename
        if not path.is_file() or _is_link_like(path) or sha256_file(path) != digest:
            raise ValueError(f"Retained input changed: {filename}")
    retained_entries = list((study_dir / "retained_input").iterdir())
    if (
        {path.name for path in retained_entries} != set(identity["retained_input_sha256"])
        or any(not path.is_file() or _is_link_like(path) for path in retained_entries)
    ):
        raise ValueError("Retained-input file set changed or contains a link/non-file")
    for file_key, hash_key in (
        ("plan_file", "plan_sha256"),
        ("registration_file", "registration_sha256"),
    ):
        path = ROOT / identity[file_key]
        if not path.is_file() or sha256_file(path) != identity[hash_key]:
            raise ValueError(f"Frozen {file_key} changed")
    original_registration = _load_json(ROOT / identity["registration_file"])
    if sha256_file(ROOT / identity["registration_file"]) != target[
        "execution_registration_sha256"
    ]:
        raise ValueError("Target execution registration hash changed")
    for executable_key, hash_key in (
        ("facets_executable", "facets_executable_sha256"),
        ("rscript_executable", "rscript_executable_sha256"),
    ):
        path = Path(identity[executable_key])
        if not path.is_file() or sha256_file(path) != identity[hash_key]:
            raise ValueError(f"Frozen executable changed: {executable_key}")
        if identity[hash_key] != original_registration[hash_key]:
            raise ValueError(f"Study used an unregistered executable: {executable_key}")
    return identity, original_registration


def validate_completion_bundle(
    study_dir: Path,
    identity: dict[str, Any],
) -> dict[str, Any]:
    attempts = _read_csv(study_dir / "retained_input" / "attempt_manifest.csv")
    if tuple(attempts.columns) != ATTEMPT_MANIFEST_COLUMNS:
        raise ValueError("Attempt manifest columns/order changed")
    if len(attempts) != 1272:
        raise ValueError("Attempt denominator is not 1272")
    if attempts["AttemptId"].astype(str).duplicated().any():
        raise ValueError("AttemptId is not unique")
    ordinals = pd.to_numeric(attempts["AttemptOrdinal"], errors="raise").astype(int)
    if ordinals.duplicated().any():
        raise ValueError("AttemptOrdinal is not unique")
    counts = attempts["AttemptType"].astype(str).value_counts().to_dict()
    expected_counts = {
        MML_FIXED_ATTEMPT: 600,
        MML_FREE_ATTEMPT: 600,
        PAIR_ATTEMPT: 36,
        CMLE_ATTEMPT: 36,
    }
    if counts != expected_counts:
        raise ValueError(f"Attempt lane denominator changed: {counts}")

    marker_hashes: dict[str, str] = {}
    actual_attempt_directories: set[str] = set()
    facets_ready = 0
    for _, attempt in attempts.iterrows():
        ordinal = int(attempt["AttemptOrdinal"])
        directory = study_dir / "attempts" / f"{ordinal:05d}"
        marker_path = directory / "completion.json"
        if not marker_path.is_file():
            raise FileNotFoundError(f"Missing completion marker: {marker_path}")
        if _is_link_like(directory) or _is_link_like(marker_path):
            raise ValueError(f"Completion marker path contains a symlink: {marker_path}")
        if {path.name for path in directory.iterdir()} != {"completion.json"}:
            raise ValueError(f"Completion marker directory {ordinal} contains extra files")
        actual_attempt_directories.add(directory.name)
        completion = _load_json(marker_path)
        _require_exact_keys(completion, COMPLETION_KEYS, f"completion {ordinal}")
        expected_fields = {
            "schema_version": "known_assignment_confirmatory_v1",
            "attempt_id": str(attempt["AttemptId"]),
            "attempt_type": str(attempt["AttemptType"]),
            "attempt_fingerprint": str(attempt["AttemptFingerprint"]),
            "run_input_sha256": str(attempt["RunInputSHA256"]),
            "dependency_manifest_sha256": identity["dependency_manifest_sha256"],
            "runner_sha256": identity["runner_sha256"],
            "facets_executable_sha256": identity["facets_executable_sha256"],
            "failure_reason": "",
        }
        for key, expected in expected_fields.items():
            if completion[key] != expected:
                raise ValueError(f"Completion {ordinal} mismatch: {key}")
        if not _require_bool(completion["execution_completed"], f"completion {ordinal} execution"):
            raise ValueError(f"Attempt {ordinal} did not complete")
        if not _require_bool(
            completion["statistical_evidence_ready"], f"completion {ordinal} readiness"
        ):
            raise ValueError(f"Attempt {ordinal} is not evidence-ready")
        calibration = completion["facets_calibration_ready"]
        if str(attempt["AttemptType"]) == PAIR_ATTEMPT:
            if not _require_bool(calibration, f"completion {ordinal} FACETS calibration"):
                raise ValueError(f"FACETS pair {ordinal} did not calibrate")
            facets_ready += 1
        elif calibration is not None:
            raise TypeError(f"Non-FACETS attempt {ordinal} has a calibration boolean")

        _validate_attempt_artifacts(study_dir, ordinal, completion)
        marker_hashes[str(attempt["AttemptId"])] = sha256_file(marker_path)

    attempt_entries = list((study_dir / "attempts").iterdir())
    if any(_is_link_like(path) or not path.is_dir() for path in attempt_entries):
        raise ValueError("Attempts root contains a non-directory or symlink entry")
    observed_directories = {path.name for path in attempt_entries}
    if observed_directories != actual_attempt_directories:
        raise ValueError("Attempt directory set differs from the frozen manifest")
    work_entries = list((study_dir / "work").iterdir())
    if (
        {path.name for path in work_entries} != actual_attempt_directories
        or any(_is_link_like(path) or not path.is_dir() for path in work_entries)
    ):
        raise ValueError("Work directory set differs or contains a link/non-directory")
    return {
        "attempts": len(attempts),
        "facets_ready": facets_ready,
        "completion_marker_set_sha256": _canonical_digest(marker_hashes),
        "marker_hashes": marker_hashes,
    }


def _require_close(observed: Any, expected: float, label: str) -> None:
    observed_float = float(observed)
    if not np.isfinite(observed_float) or not np.isclose(
        observed_float, expected, rtol=1e-12, atol=1e-14
    ):
        raise ValueError(f"{label} differs: observed={observed_float}, expected={expected}")


def _recompute_endpoint(values: np.ndarray, alternative: str) -> dict[str, float]:
    if values.shape != (200,) or not np.isfinite(values).all():
        raise ValueError("Endpoint contrast vector must contain 200 finite values")
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=1))
    se = float(sd / np.sqrt(len(values)))
    statistic = float(mean / se)
    if alternative == "greater":
        raw_p = float(stats.t.sf(statistic, df=199))
    elif alternative == "less":
        raw_p = float(stats.t.cdf(statistic, df=199))
    else:
        raise ValueError(f"Unsupported endpoint alternative: {alternative}")
    half_width = float(stats.t.ppf(0.975, df=199) * se)
    return {
        "MeanContrast": mean,
        "MonteCarloSD": sd,
        "MonteCarloSE": se,
        "TStatistic": statistic,
        "DegreesOfFreedom": 199.0,
        "RawOneSidedP": raw_p,
        "Lower95": mean - half_width,
        "Upper95": mean + half_width,
        "TwoSided95HalfWidth": half_width,
    }


def _validate_primary_secondary(aggregate_dir: Path) -> dict[str, Any]:
    primary = _read_csv(aggregate_dir / "primary_result.csv")
    secondary = _read_csv(aggregate_dir / "secondary_results.csv")
    if len(primary) != 1 or str(primary.loc[0, "EndpointId"]) != PRIMARY_ID:
        raise ValueError("Primary endpoint identity changed")
    _require_finite(
        primary,
        (
            "FiniteTriplets",
            "RequiredTriplets",
            "MeanContrast",
            "MonteCarloSD",
            "MonteCarloSE",
            "TStatistic",
            "DegreesOfFreedom",
            "RawOneSidedP",
            "Lower95",
            "Upper95",
            "TwoSided95HalfWidth",
        ),
        "primary result",
    )
    for column in ("FullTripletGate", "DirectionPass", "DirectionConfirmed", "PrecisionQualified"):
        if not _strict_bool_series(primary[column], f"primary {column}").all():
            raise ValueError(f"Primary gate is false: {column}")
    if int(primary.loc[0, "FiniteTriplets"]) != 200 or int(
        primary.loc[0, "RequiredTriplets"]
    ) != 200:
        raise ValueError("Primary triplet denominator changed")
    if len(secondary) != 3 or set(secondary["EndpointId"].astype(str)) != (
        set(ENDPOINT_ALTERNATIVES) - {PRIMARY_ID}
    ):
        raise ValueError("Secondary family is not exactly three endpoints")
    _require_finite(
        secondary,
        (
            "FiniteTriplets",
            "RequiredTriplets",
            "RawOneSidedP",
            "HolmAdjustedP",
            "MeanContrast",
            "Lower95",
            "Upper95",
        ),
        "secondary results",
    )
    for column in ("FullTripletGate", "DirectionPass", "PrimaryGatePass", "DirectionConfirmed"):
        if not _strict_bool_series(secondary[column], f"secondary {column}").all():
            raise ValueError(f"Secondary gate is false: {column}")
    if not pd.to_numeric(secondary["FiniteTriplets"]).eq(200).all():
        raise ValueError("Secondary triplet denominator changed")
    if not pd.to_numeric(secondary["RequiredTriplets"]).eq(200).all():
        raise ValueError("Secondary required triplet denominator changed")
    if set(secondary["IntervalMultiplicityStatus"].astype(str)) != {
        "unadjusted_descriptive_95_percent"
    }:
        raise ValueError("Secondary interval multiplicity label changed")

    contrasts = _read_csv(aggregate_dir / "registered_endpoint_contrasts.csv")
    if (
        len(contrasts) != 800
        or contrasts.duplicated(["EndpointId", "PersonVector"]).any()
        or set(contrasts["EndpointId"].astype(str)) != set(ENDPOINT_ALTERNATIVES)
    ):
        raise ValueError("Registered endpoint contrast key set changed")
    expected_vectors = set(range(1, 201))
    recomputed: dict[str, dict[str, float]] = {}
    for endpoint_id, alternative in ENDPOINT_ALTERNATIVES.items():
        selected = contrasts.loc[
            contrasts["EndpointId"].astype(str).eq(endpoint_id)
        ].sort_values("PersonVector")
        if set(pd.to_numeric(selected["PersonVector"]).astype(int)) != expected_vectors:
            raise ValueError(f"Endpoint PersonVector set changed: {endpoint_id}")
        recomputed[endpoint_id] = _recompute_endpoint(
            pd.to_numeric(selected["Contrast"], errors="raise").to_numpy(dtype=float),
            alternative,
        )

    result_rows = pd.concat([primary, secondary], ignore_index=True)
    for _, row in result_rows.iterrows():
        endpoint_id = str(row["EndpointId"])
        alternative = ENDPOINT_ALTERNATIVES[endpoint_id]
        if str(row["Alternative"]) != alternative:
            raise ValueError(f"Endpoint alternative changed: {endpoint_id}")
        for column, expected in recomputed[endpoint_id].items():
            _require_close(row[column], expected, f"{endpoint_id}.{column}")

    raw_secondary = {
        str(row["EndpointId"]): float(row["RawOneSidedP"])
        for _, row in secondary.iterrows()
    }
    ordered = sorted(raw_secondary, key=raw_secondary.get)
    holm: dict[str, float] = {}
    running = 0.0
    family_size = len(ordered)
    for rank, endpoint_id in enumerate(ordered):
        candidate = min(1.0, (family_size - rank) * raw_secondary[endpoint_id])
        running = max(running, candidate)
        holm[endpoint_id] = running
    for _, row in secondary.iterrows():
        endpoint_id = str(row["EndpointId"])
        _require_close(
            row["HolmAdjustedP"], holm[endpoint_id], f"{endpoint_id}.HolmAdjustedP"
        )
    _require_close(
        primary.loc[0, "AdjustedP"],
        float(primary.loc[0, "RawOneSidedP"]),
        "primary.AdjustedP",
    )
    if not float(primary.loc[0, "TwoSided95HalfWidth"]) <= float(
        primary.loc[0, "PrecisionTarget"]
    ) == 0.015:
        raise ValueError("Primary precision target/gate changed")
    return {
        "primary_mean": float(primary.loc[0, "MeanContrast"]),
        "primary_lower95": float(primary.loc[0, "Lower95"]),
        "primary_upper95": float(primary.loc[0, "Upper95"]),
        "primary_half_width": float(primary.loc[0, "TwoSided95HalfWidth"]),
        "primary_direction_confirmed": True,
        "primary_precision_qualified": True,
        "secondary_family_size": 3,
        "secondary_all_holm_direction_confirmed": True,
        "secondary_intervals_multiplicity_status": sorted(
            secondary["IntervalMultiplicityStatus"].astype(str).unique().tolist()
        ),
        "all_four_endpoint_statistics_recomputed": True,
        "secondary_holm_recomputed": True,
    }


def _validate_mml_keys(aggregate_dir: Path) -> dict[str, Any]:
    recovery = _read_csv(aggregate_dir / "recovery.csv")
    selected = recovery.loc[
        recovery["EstimatorMode"].astype(str).isin((MML_FIXED_MODE, MML_FREE_MODE))
        & recovery["Facet"].astype(str).eq("Rater")
    ].copy()
    selected["IncludedStrict"] = _strict_bool_series(
        selected["IncludedInStudy"], "MML recovery IncludedInStudy"
    )
    if not selected["IncludedStrict"].all():
        raise ValueError("An MML Rater recovery row is excluded")
    _require_finite(
        selected,
        (
            "PersonVector",
            "Gamma",
            "Truth",
            "Estimate",
            "RawError",
            "EstimateAligned",
            "TruthAligned",
            "ErrorAligned",
        ),
        "MML Rater recovery",
    )
    numeric_recovery = {
        column: pd.to_numeric(selected[column], errors="raise").to_numpy(dtype=float)
        for column in (
            "Truth",
            "Estimate",
            "RawError",
            "EstimateAligned",
            "TruthAligned",
            "ErrorAligned",
        )
    }
    expected_truth = selected["Level"].astype(str).map(RATER_TRUTH).to_numpy(dtype=float)
    aligned_estimate = (
        pd.Series(numeric_recovery["Estimate"], index=selected.index)
        - pd.Series(numeric_recovery["Estimate"], index=selected.index)
        .groupby(selected["AttemptId"])
        .transform("mean")
    ).to_numpy(dtype=float)
    algebra_checks = {
        "Truth": (numeric_recovery["Truth"], expected_truth),
        "TruthAligned": (numeric_recovery["TruthAligned"], expected_truth),
        "RawError": (
            numeric_recovery["RawError"],
            numeric_recovery["Estimate"] - numeric_recovery["Truth"],
        ),
        "EstimateAligned": (numeric_recovery["EstimateAligned"], aligned_estimate),
        "ErrorAligned": (
            numeric_recovery["ErrorAligned"],
            numeric_recovery["EstimateAligned"] - numeric_recovery["TruthAligned"],
        ),
    }
    for label, (observed, expected) in algebra_checks.items():
        if not np.allclose(observed, expected, rtol=0, atol=1e-12):
            raise ValueError(f"MML Rater recovery algebra changed: {label}")
    groups = selected.groupby(["PersonVector", "Gamma", "EstimatorMode"], sort=False)
    if groups.ngroups != 1200:
        raise ValueError("MML Rater recovery does not contain 1200 cells")
    observed_cells = {
        (int(person_vector), float(gamma), str(mode))
        for person_vector, gamma, mode in groups.groups
    }
    expected_cells = {
        (person_vector, gamma, mode)
        for person_vector in range(1, 201)
        for gamma in (-0.8, 0.0, 0.8)
        for mode in (MML_FIXED_MODE, MML_FREE_MODE)
    }
    if observed_cells != expected_cells:
        raise ValueError("MML recovery cell keys differ from the registered direct product")
    expected_levels = set(RATER_TRUTH)
    for keys, group in groups:
        levels = group["Level"].astype(str).tolist()
        if len(levels) != 4 or set(levels) != expected_levels or len(set(levels)) != 4:
            raise ValueError(f"MML Rater key set changed for {keys}: {levels}")

    runs = _read_csv(aggregate_dir / "run_ledger.csv")
    mml_runs = runs.loc[runs["EstimatorMode"].astype(str).isin((MML_FIXED_MODE, MML_FREE_MODE))]
    if len(mml_runs) != 1200 or mml_runs["AttemptId"].astype(str).nunique() != 1200:
        raise ValueError("MML run ledger denominator/key changed")
    for column in ("Converged", "InferenceReady", "IncludedInStudy", "ConstraintPass"):
        if not _strict_bool_series(mml_runs[column], f"MML run {column}").all():
            raise ValueError(f"MML run gate is false: {column}")
    _require_finite(mml_runs, ("LogLik", "MaxAbsConstraintResidual"), "MML run ledger")

    manifest = _read_csv(aggregate_dir.parent / "retained_input" / "attempt_manifest.csv")
    manifest_mml = manifest.loc[
        manifest["AttemptType"].astype(str).isin((MML_FIXED_ATTEMPT, MML_FREE_ATTEMPT)),
        ["AttemptId", "RunId", "AttemptType", "PersonVector", "Gamma"],
    ].copy()
    manifest_mml["EstimatorMode"] = manifest_mml["AttemptType"].astype(str).map(
        {MML_FIXED_ATTEMPT: MML_FIXED_MODE, MML_FREE_ATTEMPT: MML_FREE_MODE}
    )
    expected_attempt_keys = set(
        map(
            tuple,
            manifest_mml[["AttemptId", "RunId", "EstimatorMode"]]
            .astype(str)
            .to_numpy(),
        )
    )
    run_attempt_keys = set(
        map(
            tuple,
            mml_runs[["AttemptId", "RunId", "EstimatorMode"]].astype(str).to_numpy(),
        )
    )
    if run_attempt_keys != expected_attempt_keys:
        raise ValueError("MML run ledger does not exactly join to the frozen attempt manifest")
    expected_recovery_keys = {
        (
            str(row.AttemptId),
            str(row.RunId),
            str(row.EstimatorMode),
            int(row.PersonVector),
            float(row.Gamma),
            level,
        )
        for row in manifest_mml.itertuples(index=False)
        for level in RATER_TRUTH
    }
    observed_recovery_keys = {
        (
            str(row.AttemptId),
            str(row.RunId),
            str(row.EstimatorMode),
            int(row.PersonVector),
            float(row.Gamma),
            str(row.Level),
        )
        for row in selected.itertuples(index=False)
    }
    if (
        len(selected) != len(expected_recovery_keys)
        or observed_recovery_keys != expected_recovery_keys
    ):
        raise ValueError("MML Rater recovery does not exactly join to attempt x Rater keys")

    constraints = _read_csv(aggregate_dir / "constraints.csv")
    constraints = constraints.loc[
        constraints["EstimatorMode"].astype(str).isin((MML_FIXED_MODE, MML_FREE_MODE))
    ].copy()
    if len(constraints) != 1200 or constraints["AttemptId"].astype(str).nunique() != 1200:
        raise ValueError("MML constraint denominator/key changed")
    if not _strict_bool_series(constraints["ConstraintPass"], "MML constraints").all():
        raise ValueError("An MML constraint failed")
    _require_finite(constraints, ("MaxAbsConstraintResidual",), "MML constraints")
    if (
        pd.to_numeric(constraints["MaxAbsConstraintResidual"], errors="raise").abs()
        > 1e-8
    ).any():
        raise ValueError("MML constraint residual exceeds the posthoc integrity tolerance")
    constraint_attempt_keys = set(
        map(
            tuple,
            constraints[["AttemptId", "RunId", "EstimatorMode"]]
            .astype(str)
            .to_numpy(),
        )
    )
    if constraint_attempt_keys != expected_attempt_keys:
        raise ValueError("MML constraints do not exactly join to the frozen attempt manifest")

    # Rebuild the registered primary endpoint only from the free-MML Rater
    # recovery rows.  This is a lineage check: FACETS display/fit fields cannot
    # silently enter KA1 merely because the aggregate summary still hashes.
    free = selected.loc[selected["EstimatorMode"].astype(str).eq(MML_FREE_MODE)].copy()
    free["PersonVector"] = pd.to_numeric(free["PersonVector"], errors="raise").astype(int)
    free["Gamma"] = pd.to_numeric(free["Gamma"], errors="raise").astype(float)
    free["SquaredAlignedError"] = np.square(
        pd.to_numeric(free["ErrorAligned"], errors="raise")
    )
    rmse = (
        free.groupby(["RunId", "PersonVector", "Gamma"], as_index=False)[
            "SquaredAlignedError"
        ]
        .mean()
        .rename(columns={"SquaredAlignedError": "MeanSquaredAlignedError"})
    )
    rmse["RMSE"] = np.sqrt(rmse["MeanSquaredAlignedError"])
    wide = rmse.pivot(index="PersonVector", columns="Gamma", values="RMSE")
    if len(wide) != 200 or set(wide.columns.astype(float)) != {-0.8, 0.0, 0.8}:
        raise ValueError("Free-MML recovery does not form 200 complete KA1 triplets")
    rebuilt = (
        (0.5 * (wide[-0.8] + wide[0.8]) - wide[0.0])
        .rename("Contrast")
        .reset_index()
        .sort_values("PersonVector")
    )
    registered = _read_csv(aggregate_dir / "registered_endpoint_contrasts.csv")
    registered = registered.loc[
        registered["EndpointId"].astype(str).eq(PRIMARY_ID),
        ["PersonVector", "Contrast"],
    ].copy()
    registered["PersonVector"] = pd.to_numeric(
        registered["PersonVector"], errors="raise"
    ).astype(int)
    registered = registered.sort_values("PersonVector")
    if len(registered) != 200 or not np.array_equal(
        registered["PersonVector"].to_numpy(dtype=int),
        rebuilt["PersonVector"].to_numpy(dtype=int),
    ) or not np.allclose(
        registered["Contrast"].to_numpy(dtype=float),
        rebuilt["Contrast"].to_numpy(dtype=float),
        rtol=0,
        atol=1e-12,
    ):
        raise ValueError("Registered KA1 contrasts do not reconstruct from free-MML Rater errors")
    return {
        "mml_attempts": 1200,
        "mml_rater_cells": 1200,
        "expected_cell_direct_product_exact": True,
        "attempt_manifest_join_exact": True,
        "rater_rows_per_cell": 4,
        "rater_key_set": sorted(expected_levels),
        "all_rater_errors_finite": True,
        "rater_recovery_algebra_recomputed": True,
        "all_constraints_pass": True,
        "ka1_lineage_reconstructed_from_free_mml_rater_errors": True,
        "facets_display_or_fit_values_in_ka1_lineage": False,
    }


def _validate_estimator_lane_counts(aggregate_dir: Path) -> dict[str, Any]:
    runs = _read_csv(aggregate_dir / "run_ledger.csv")
    expected = {
        MML_FREE_MODE: 600,
        MML_FIXED_MODE: 600,
        FACETS_MODE: 36,
        "PYTHON_JMLE": 36,
        "PYTHON_EXACT_CMLE": 36,
    }
    if runs["EstimatorMode"].astype(str).value_counts().to_dict() != expected:
        raise ValueError("Estimator-mode ledger denominator changed")
    result: dict[str, Any] = {}
    for mode, planned in expected.items():
        selected = runs.loc[runs["EstimatorMode"].astype(str).eq(mode)]
        lane = {"planned": planned, "attempted": len(selected)}
        for column, output_key in (
            ("FitReturned", "returned"),
            ("Converged", "converged"),
            ("InferenceReady", "inference_ready"),
            ("IncludedInStudy", "included_in_study"),
        ):
            lane[output_key] = int(
                _strict_bool_series(selected[column], f"{mode}.{column}").sum()
            )
        if any(value != planned for key, value in lane.items() if key != "planned"):
            raise ValueError(f"Estimator lane is incomplete: {mode}")
        result[mode] = lane
    return result


def validate_aggregate_bundle(
    study_dir: Path,
    identity: dict[str, Any],
    completion: dict[str, Any],
) -> dict[str, Any]:
    aggregate_dir = study_dir / "aggregate"
    aggregate_entries = list(aggregate_dir.iterdir())
    if any(_is_link_like(path) or not path.is_file() for path in aggregate_entries):
        raise ValueError("Aggregate contains a link or non-file entry")
    aggregate_identity = _load_json(aggregate_dir / "aggregate_identity.json")
    if aggregate_identity["study_identity_sha256"] != sha256_file(
        study_dir / "study_identity.json"
    ):
        raise ValueError("Aggregate points to a different study identity")
    if aggregate_identity["completion_marker_set_sha256"] != completion[
        "completion_marker_set_sha256"
    ]:
        raise ValueError("Aggregate completion-marker digest changed")
    r_assessment_path = study_dir / "r_mml" / "assessment.json"
    if aggregate_identity["r_mml_assessment_sha256"] != sha256_file(r_assessment_path):
        raise ValueError("Aggregate R assessment parent changed")
    recorded = aggregate_identity["artifact_sha256"]
    actual = {
        path.name: sha256_file(path)
        for path in aggregate_dir.iterdir()
        if path.is_file() and path.name != "aggregate_identity.json"
    }
    if recorded != actual:
        raise ValueError("Aggregate file set or hash changed")

    assessment = _load_json(aggregate_dir / "assessment.json")
    for key in (
        "ScientificEndpointComputed",
        "ScientificInferenceValidated",
        "FACETSComplementQualification",
        "FACETSValidatedComplementaryWorkbench",
        "primary_direction_confirmed",
        "primary_precision_qualified",
    ):
        if not _require_bool(assessment.get(key), f"aggregate assessment {key}"):
            raise ValueError(f"Original aggregate status is false: {key}")
    if _require_bool(
        assessment.get("ordinary_facets_fit_used_as_raw_input"),
        "ordinary_facets_fit_used_as_raw_input",
    ):
        raise ValueError("FACETS ordinary fit was used as raw input")
    for family in ("scientific_gates", "facets_gates", "cmle_descriptive_gates"):
        gates = assessment.get(family)
        if not isinstance(gates, dict) or not gates:
            raise ValueError(f"Missing aggregate gate family: {family}")
        if not all(_require_bool(value, f"{family}.{key}") for key, value in gates.items()):
            raise ValueError(f"Aggregate gate family contains a failure: {family}")
    return {
        "aggregate_file_hashes_exact": True,
        "original_status": {
            "OriginalScientificInferenceValidated": True,
            "OriginalFACETSComplementQualification": True,
            "OriginalFACETSValidatedComplementaryWorkbench": True,
            "meaning": "PASS under the prospectively frozen v1 gates",
            "overwritten_or_reinterpreted_as_original": False,
            "source": "aggregate/assessment.json",
            "source_sha256": sha256_file(aggregate_dir / "assessment.json"),
        },
        "endpoints": _validate_primary_secondary(aggregate_dir),
        "mml_keys": _validate_mml_keys(aggregate_dir),
        "estimator_lanes": _validate_estimator_lane_counts(aggregate_dir),
    }


def _expected_r_run_ids(attempts: pd.DataFrame) -> set[str]:
    vector = pd.to_numeric(attempts["PersonVector"], errors="raise").astype(int)
    selected = attempts.loc[
        attempts["AttemptType"].astype(str).eq(MML_FREE_ATTEMPT)
        & vector.isin(CALIBRATION_TRIPLETS)
    ]
    run_ids = set(selected["RunId"].astype(str))
    if len(selected) != 36 or len(run_ids) != 36:
        raise ValueError("Frozen R subset is not 36 unique datasets")
    return run_ids


def _require_exact_group_keys(
    frame: pd.DataFrame,
    group_column: str,
    key_columns: tuple[str, ...],
    expected_groups: set[str],
    expected_keys: set[Any],
    label: str,
) -> None:
    groups = set(frame[group_column].astype(str))
    if groups != expected_groups:
        raise ValueError(f"{label} RunId set changed")
    if frame.duplicated([group_column, *key_columns]).any():
        raise ValueError(f"{label} contains duplicate keys")
    for group_id, group in frame.groupby(group_column, sort=False):
        if len(key_columns) == 1:
            observed = set(group[key_columns[0]].astype(str))
        else:
            observed = set(map(tuple, group.loc[:, key_columns].astype(str).to_numpy()))
        if observed != expected_keys:
            raise ValueError(f"{label} key set changed for {group_id}")


def _r_subset_sensitivity(
    parameters: pd.DataFrame,
    runs: pd.DataFrame,
    aggregate_dir: Path,
) -> dict[str, Any]:
    rater = parameters.loc[parameters["Block"].astype(str).eq("Rater")].copy()
    metadata = runs[["RunId", "PersonVector", "Gamma"]].copy()
    metadata["PersonVector"] = pd.to_numeric(
        metadata["PersonVector"], errors="raise"
    ).astype(int)
    metadata["Gamma"] = pd.to_numeric(metadata["Gamma"], errors="raise").astype(float)
    solutions = {
        "Python": "PythonEstimate",
        "R_Q31": "REstimateQ31",
        "R_Q61": "REstimateQ61",
    }
    contrast_tables: dict[str, pd.DataFrame] = {}
    for solution, column in solutions.items():
        work = rater[["RunId", "Level", column]].copy()
        work["Truth"] = work["Level"].astype(str).map(RATER_TRUTH)
        _require_finite(work, (column, "Truth"), f"{solution} Rater parameters")
        work["RawError"] = pd.to_numeric(work[column]) - pd.to_numeric(work["Truth"])
        work["AlignedError"] = work["RawError"] - work.groupby("RunId")["RawError"].transform(
            "mean"
        )
        rmse = (
            work.groupby("RunId", as_index=False)["AlignedError"]
            .agg(RMSE=lambda values: float(np.sqrt(np.mean(np.square(values)))))
            .merge(metadata, on="RunId", validate="one_to_one")
        )
        wide = rmse.pivot(index="PersonVector", columns="Gamma", values="RMSE")
        if set(wide.columns.astype(float)) != {-0.8, 0.0, 0.8} or len(wide) != 12:
            raise ValueError(f"{solution} R subset is not 12 complete triplets")
        contrast = 0.5 * (wide[-0.8] + wide[0.8]) - wide[0.0]
        contrast_tables[solution] = contrast.rename("Contrast").reset_index()

    registered = _read_csv(aggregate_dir / "registered_endpoint_contrasts.csv")
    registered = registered.loc[
        registered["EndpointId"].astype(str).eq(PRIMARY_ID)
        & pd.to_numeric(registered["PersonVector"]).astype(int).isin(CALIBRATION_TRIPLETS),
        ["PersonVector", "Contrast"],
    ].copy()
    registered["PersonVector"] = pd.to_numeric(
        registered["PersonVector"], errors="raise"
    ).astype(int)
    registered = registered.sort_values("PersonVector")
    python = contrast_tables["Python"].sort_values("PersonVector")
    if not np.array_equal(
        registered["PersonVector"].to_numpy(dtype=int),
        python["PersonVector"].to_numpy(dtype=int),
    ) or not np.allclose(
        registered["Contrast"].to_numpy(dtype=float),
        python["Contrast"].to_numpy(dtype=float),
        rtol=0,
        atol=1e-12,
    ):
        raise ValueError("R-subset Python reconstruction differs from the registered endpoint")

    python_values = contrast_tables["Python"]["Contrast"].to_numpy(dtype=float)
    result: dict[str, Any] = {
        "role": "POSTHOC_DESCRIPTIVE_NUMERICAL_SENSITIVITY",
        "triplets": 12,
        "datasets": 36,
        "new_confirmatory_test_constructed": False,
        "full_600_dataset_r_validation_claim": False,
        "solutions": {},
    }
    for solution, table in contrast_tables.items():
        values = table["Contrast"].to_numpy(dtype=float)
        result["solutions"][solution] = {
            "mean_contrast": float(np.mean(values)),
            "positive_triplets": int(np.sum(values > 0)),
            "maximum_absolute_triplet_change_from_python": float(
                np.max(np.abs(values - python_values))
            ),
        }
    return result


def validate_r_bundle(
    study_dir: Path,
    original_registration: dict[str, Any],
    plan: dict[str, Any],
) -> dict[str, Any]:
    output_dir = study_dir / "r_mml"
    observed_files = {path.name for path in output_dir.iterdir() if path.is_file()}
    if observed_files != R_BUNDLE_FILES or any(
        _is_link_like(path) or not path.is_file() for path in output_dir.iterdir()
    ):
        raise ValueError("R bundle file set changed or contains a symlink/non-file")
    verification = _load_json(output_dir / "verification_identity.json")
    assessment = _load_json(output_dir / "assessment.json")
    process = _load_json(output_dir / "r_process.json")
    if verification.get("schema_version") != "known_assignment_mml_crossfit_v1_identity_v1":
        raise ValueError("R verification identity schema changed")
    if _require_json_int(
        verification.get("expected_datasets"), "R verification expected_datasets"
    ) != 36:
        raise ValueError("R verification denominator changed")
    expected_script_hashes = {
        "r_script_sha256": sha256_file(ROOT / "validation" / "known_assignment_mml_crossfit.R"),
        "python_verifier_sha256": sha256_file(
            ROOT / "validation" / "known_assignment_mml_crossfit.py"
        ),
    }
    for key, digest in expected_script_hashes.items():
        if verification.get(key) != digest or assessment.get(key) != digest:
            raise ValueError(f"R verifier source hash changed: {key}")
    if not isinstance(verification.get("source_sha256"), dict):
        raise TypeError("R source hash map must be an object")
    if not isinstance(verification.get("retained_input_sha256"), dict):
        raise TypeError("R retained-input hash map must be an object")
    _require_exact_keys(
        verification["source_sha256"], R_RAW_SOURCE_FILES, "R raw source hash map"
    )
    _require_exact_keys(
        verification["retained_input_sha256"],
        R_RETAINED_INPUT_FILES,
        "R retained-input hash map",
    )
    for filename, digest in verification["source_sha256"].items():
        path = output_dir / filename
        if Path(filename).name != filename or _is_link_like(path) or sha256_file(path) != digest:
            raise ValueError(f"Raw R output changed: {filename}")
    for filename, digest in verification["retained_input_sha256"].items():
        path = study_dir / "retained_input" / filename
        if Path(filename).name != filename or _is_link_like(path) or sha256_file(path) != digest:
            raise ValueError(f"R verifier retained input changed: {filename}")

    command = process.get("command")
    if not isinstance(command, list) or len(command) != 10:
        raise ValueError("R process command shape changed")
    if Path(command[0]).resolve() != Path(original_registration["rscript_executable"]).resolve():
        raise ValueError("R process used an unregistered Rscript")
    if Path(command[1]).resolve() != (ROOT / "validation" / "known_assignment_mml_crossfit.R"):
        raise ValueError("R process used a different script")
    options = dict(zip(command[2::2], command[3::2], strict=True))
    expected_options = {
        "--study": str(study_dir.resolve()),
        "--output": str(output_dir.resolve()),
        "--vectors": ",".join(map(str, CALIBRATION_TRIPLETS)),
        "--maxit": str(int(plan["strict_integrity_contract"]["independent_r"]["required_maxit"])),
    }
    if options != expected_options:
        raise ValueError(f"R process command options changed: {options}")
    if _require_json_int(process.get("returncode"), "R process returncode") != 0:
        raise ValueError("R process did not return zero")

    if assessment.get("schema_version") != "known_assignment_mml_crossfit_v1":
        raise ValueError("R assessment schema changed")
    if not _require_bool(assessment.get("pass"), "R assessment pass"):
        raise ValueError("R assessment is not PASS")
    if (
        _require_json_int(assessment.get("datasets"), "R assessment datasets") != 36
        or assessment.get("r_version") != "4.5.1"
    ):
        raise ValueError("R assessment runtime/denominator changed")
    if assessment.get("tolerance") != original_registration["r_mml_tolerance"]:
        raise ValueError("R assessment tolerance differs from registration")
    gates = assessment.get("gates")
    if not isinstance(gates, dict):
        raise TypeError("R assessment gates must be an object")
    _require_exact_keys(gates, R_GATES, "R assessment gates")
    if not all(_require_bool(value, f"R gate {key}") for key, value in gates.items()):
        raise ValueError("At least one registered R gate failed")
    metrics = assessment.get("metrics")
    if not isinstance(metrics, dict) or not metrics:
        raise ValueError("R assessment metrics are missing")
    metric_values = np.asarray(
        [
            _require_json_number(value, f"R assessment metric {key}")
            for key, value in metrics.items()
        ],
        dtype=float,
    )

    attempts = _read_csv(study_dir / "retained_input" / "attempt_manifest.csv")
    expected_run_ids = _expected_r_run_ids(attempts)
    runs = _read_csv(
        output_dir / "crossfit_runs.csv", R_CSV_COLUMNS["crossfit_runs.csv"]
    )
    parameters = _read_csv(
        output_dir / "crossfit_parameters.csv", R_CSV_COLUMNS["crossfit_parameters.csv"]
    )
    gradients = _read_csv(
        output_dir / "crossfit_gradients.csv", R_CSV_COLUMNS["crossfit_gradients.csv"]
    )
    runtime = _read_csv(
        output_dir / "runtime_identity.csv", R_CSV_COLUMNS["runtime_identity.csv"]
    )
    evaluations = _read_csv(
        output_dir / "crossfit_cross_evaluations.csv",
        R_CSV_COLUMNS["crossfit_cross_evaluations.csv"],
    )
    if len(runs) != 36 or set(runs["RunId"].astype(str)) != expected_run_ids or runs["RunId"].duplicated().any():
        raise ValueError("R cross-fit run keys changed")
    run_numeric = [
        column
        for column in runs.columns
        if column not in {"RunId", "RMessage", "RQ61Message"}
    ]
    _require_finite(runs, run_numeric, "R cross-fit runs")
    if not pd.to_numeric(runs["RConvergenceCode"]).eq(0).all() or not pd.to_numeric(
        runs["RQ61ConvergenceCode"]
    ).eq(0).all():
        raise ValueError("R optimizer convergence code changed")
    _require_exact_group_keys(
        parameters,
        "RunId",
        ("Block", "Level"),
        expected_run_ids,
        R_PARAMETER_KEYS,
        "R parameters",
    )
    _require_finite(
        parameters,
        (
            "PythonEstimate",
            "REstimateQ31",
            "REstimateQ61",
            "DifferenceRQ31MinusPython",
            "DifferenceRQ61MinusRQ31",
        ),
        "R parameters",
    )
    for column in (
        "PythonEstimate",
        "REstimateQ31",
        "REstimateQ61",
        "DifferenceRQ31MinusPython",
        "DifferenceRQ61MinusRQ31",
    ):
        parameters[column] = pd.to_numeric(parameters[column], errors="raise")
    if not np.allclose(
        pd.to_numeric(parameters["DifferenceRQ31MinusPython"]),
        pd.to_numeric(parameters["REstimateQ31"])
        - pd.to_numeric(parameters["PythonEstimate"]),
        rtol=0,
        atol=1e-12,
    ) or not np.allclose(
        pd.to_numeric(parameters["DifferenceRQ61MinusRQ31"]),
        pd.to_numeric(parameters["REstimateQ61"])
        - pd.to_numeric(parameters["REstimateQ31"]),
        rtol=0,
        atol=1e-12,
    ):
        raise ValueError("R parameter difference columns fail their defining algebra")
    maximum_parameter_constraint_residual = 0.0
    for estimate in ("PythonEstimate", "REstimateQ31", "REstimateQ61"):
        for block in ("Rater", "Task"):
            residual = (
                parameters.loc[parameters["Block"].astype(str).eq(block)]
                .groupby("RunId")[estimate]
                .sum()
                .abs()
                .max()
            )
            maximum_parameter_constraint_residual = max(
                maximum_parameter_constraint_residual, float(residual)
            )
        step = parameters.loc[parameters["Block"].astype(str).eq("Step")].copy()
        step["Criterion"] = step["Level"].astype(str).str.split("::", n=1).str[0]
        residual = step.groupby(["RunId", "Criterion"])[estimate].sum().abs().max()
        maximum_parameter_constraint_residual = max(
            maximum_parameter_constraint_residual, float(residual)
        )
    if maximum_parameter_constraint_residual > 1e-12:
        raise ValueError("R parameter identification constraints are not preserved")
    _require_exact_group_keys(
        gradients,
        "RunId",
        ("Coordinate",),
        expected_run_ids,
        R_GRADIENT_KEYS,
        "R gradients",
    )
    _require_finite(
        gradients,
        (
            "PythonSolutionRGradientQ31",
            "RQ31SolutionRGradientQ31",
            "RQ61SolutionRGradientQ61",
        ),
        "R gradients",
    )
    for coordinate_column, supnorm_column in (
        ("PythonSolutionRGradientQ31", "PythonSolutionRGradientSupNormQ31"),
        ("RQ31SolutionRGradientQ31", "RGradientSupNormQ31"),
        ("RQ61SolutionRGradientQ61", "RGradientSupNormQ61"),
    ):
        reconstructed = (
            gradients.assign(
                AbsoluteCoordinate=pd.to_numeric(
                    gradients[coordinate_column], errors="raise"
                ).abs()
            )
            .groupby("RunId")["AbsoluteCoordinate"]
            .max()
            .sort_index()
        )
        recorded_supnorm = (
            runs.set_index("RunId")[supnorm_column].astype(float).sort_index()
        )
        if not np.allclose(
            reconstructed.to_numpy(dtype=float),
            recorded_supnorm.to_numpy(dtype=float),
            rtol=0,
            atol=float(
                plan["strict_integrity_contract"]["independent_r"][
                    "gradient_supnorm_serialization_atol"
                ]
            ),
        ):
            raise ValueError(f"R gradient sup norm does not reconstruct: {supnorm_column}")
    if (
        len(evaluations) != 36
        or set(evaluations["RunId"].astype(str)) != expected_run_ids
        or evaluations["RunId"].astype(str).duplicated().any()
    ):
        raise ValueError("R cross-evaluation keys changed")
    evaluation_numeric = [
        column
        for column in evaluations.columns
        if column not in {"RunId", "RMessage", "RQ61Message"}
    ]
    _require_finite(evaluations, evaluation_numeric, "R cross-evaluations")
    minimum_improvement = float(
        plan["strict_integrity_contract"]["independent_r"][
            "minimum_matched_reoptimization_improvement_posthoc_integrity_tolerance"
        ]
    )
    if not (
        pd.to_numeric(evaluations["RQ31ImprovementOverPythonQ31"])
        >= minimum_improvement
    ).all() or not (
        pd.to_numeric(evaluations["RQ61OptimizedMinusRQ31EvaluatedQ61"])
        >= minimum_improvement
    ).all():
        raise ValueError("An R reoptimization worsened its matched likelihood")
    if len(runtime) != 1 or str(runtime.loc[0, "RVersion"]) != "4.5.1":
        raise ValueError("R runtime identity changed")
    if str(runtime.loc[0, "SelectedVectors"]) != ",".join(map(str, CALIBRATION_TRIPLETS)):
        raise ValueError("R runtime selected vectors changed")
    if int(runtime.loc[0, "Datasets"]) != 36:
        raise ValueError("R runtime dataset denominator changed")

    # Reconstruct the verifier output in a temporary directory.  This makes a
    # hand-written assessment.json insufficient to qualify the bundle.
    from validation import known_assignment_mml_crossfit as crossfit

    with tempfile.TemporaryDirectory(prefix="kac200_strict_r_recheck_") as temporary:
        temporary_path = Path(temporary)
        for filename in verification["source_sha256"]:
            shutil.copy2(output_dir / filename, temporary_path / filename)
        recomputed = crossfit.assess_crossfit(
            study_dir=study_dir,
            output_dir=temporary_path,
            expected_datasets=36,
            tolerance=dict(original_registration["r_mml_tolerance"]),
        )
        if recomputed != assessment:
            raise ValueError("R assessment cannot be reproduced from raw outputs")
        if sha256_file(temporary_path / "crossfit_cross_evaluations.csv") != sha256_file(
            output_dir / "crossfit_cross_evaluations.csv"
        ):
            raise ValueError("R cross-evaluations cannot be reproduced")
        if _load_json(temporary_path / "verification_identity.json") != verification:
            raise ValueError("R verification identity cannot be reproduced")

    python_gradient = pd.to_numeric(runs["PythonSolutionRGradientSupNormQ31"]).to_numpy(
        dtype=float
    )
    if "python_terminal_gradient_q31" in original_registration["r_mml_tolerance"]:
        raise ValueError("Unexpected retrospective Python terminal-gradient threshold")
    stationarity = {
        "python_terminal_score_was_recorded": True,
        "python_terminal_score_had_frozen_acceptance_threshold": False,
        "OptimizationStationarityQualification": "INCOMPLETE",
        "retrospective_stationarity_pass_fail_decision_made": False,
        "posthoc_threshold_selection_prohibited": True,
        "observed_values_role": "DESCRIPTIVE_DIAGNOSTIC_ONLY",
        "descriptive_cutpoints_selected_posthoc": True,
        "acceptance_or_rejection_semantics": False,
        "datasets": 36,
        "minimum": float(np.min(python_gradient)),
        "mean": float(np.mean(python_gradient)),
        "median": float(np.median(python_gradient)),
        "maximum": float(np.max(python_gradient)),
        "descriptive_count_above_0p01": int(np.sum(python_gradient > 0.01)),
        "descriptive_count_above_0p1": int(np.sum(python_gradient > 0.1)),
        "all_values_sha256": hashlib.sha256(python_gradient.tobytes()).hexdigest(),
        "maximum_r_q31_terminal_gradient": float(
            pd.to_numeric(runs["RGradientSupNormQ31"]).max()
        ),
        "maximum_r_q61_terminal_gradient": float(
            pd.to_numeric(runs["RGradientSupNormQ61"]).max()
        ),
    }
    return {
        "strict_r_bundle_integrity": True,
        "raw_output_hashes": {
            filename: sha256_file(output_dir / filename)
            for filename in sorted(verification["source_sha256"])
        },
        "cross_evaluations_sha256": sha256_file(
            output_dir / "crossfit_cross_evaluations.csv"
        ),
        "parameter_difference_algebra_exact_within_serialization_tolerance": True,
        "maximum_parameter_constraint_residual": maximum_parameter_constraint_residual,
        "gradient_supnorms_reconstructed": True,
        "matched_likelihood_reoptimizations_nonworsening": True,
        "assessment_recomputed": True,
        "stationarity_audit": stationarity,
        "r_subset_sensitivity": _r_subset_sensitivity(
            parameters, runs, study_dir / "aggregate"
        ),
    }


def _validate_facets_pair_metrics(
    metrics: dict[str, Any],
    *,
    ordinal: int,
    expected_run_id: str,
    identity: dict[str, Any],
    work: Path,
    maximum_tries: int,
    replay_tolerance: float,
) -> int:
    _require_exact_keys(metrics, FACETS_PAIR_METRIC_KEYS, f"FACETS pair metrics {ordinal}")
    for key in (
        "calibration_ready",
        "direct_agreement_pass",
        "pair_fully_qualified",
        "statistical_evidence_ready",
    ):
        if not _require_bool(metrics[key], f"FACETS pair {ordinal} {key}"):
            raise ValueError(f"FACETS pair {ordinal} failed {key}")
    if _require_bool(
        metrics["confirmatory_evidence_replaced"],
        f"FACETS pair {ordinal} confirmatory_evidence_replaced",
    ):
        raise ValueError(f"FACETS pair {ordinal} replaced scientific evidence")
    expected_pair_fields = {
        "schema_version": "mfrm-facets-resilient-pair-v1",
        "facets_executable_sha256": identity["facets_executable_sha256"],
        "facets_reported_version": "4.5.0",
        "facets_final_error": "",
        "run_id": expected_run_id,
    }
    for key, expected in expected_pair_fields.items():
        if metrics[key] != expected:
            raise ValueError(f"FACETS pair {ordinal} mismatch: {key}")
    tries = _require_json_int(metrics["facets_tries"], f"FACETS pair {ordinal} tries")
    retries = _require_json_int(
        metrics["facets_retry_count"], f"FACETS pair {ordinal} retry count"
    )
    if not 1 <= tries <= maximum_tries or retries != tries - 1:
        raise ValueError(f"FACETS pair {ordinal} retry ledger is inconsistent")
    executable_size = _require_json_int(
        metrics["facets_executable_size_bytes"],
        f"FACETS pair {ordinal} executable size",
    )
    if executable_size != Path(identity["facets_executable"]).stat().st_size:
        raise ValueError(f"FACETS pair {ordinal} executable size changed")
    replay_values = np.asarray(
        [
            _require_json_number(
                metrics["python_replay_main_max_abs_difference"],
                f"FACETS pair {ordinal} main replay difference",
            ),
            _require_json_number(
                metrics["python_replay_threshold_max_abs_difference"],
                f"FACETS pair {ordinal} threshold replay difference",
            ),
            _require_json_number(
                metrics["python_replay_tolerance"],
                f"FACETS pair {ordinal} replay tolerance",
            ),
        ],
        dtype=float,
    )
    if (
        not np.isfinite(replay_values).all()
        or replay_values[2] != replay_tolerance
        or replay_values[0] > replay_tolerance
        or replay_values[1] > replay_tolerance
    ):
        raise ValueError(f"FACETS pair {ordinal} Python replay contract changed")
    dependency_manifest = _load_json(work / "dependency_manifest.json")
    if _canonical_digest(dependency_manifest) != metrics["dependency_manifest_sha256"]:
        raise ValueError(f"FACETS pair {ordinal} dependency digest changed")
    return tries


def validate_facets_bundle(
    study_dir: Path,
    identity: dict[str, Any],
    original_registration: dict[str, Any],
    plan: dict[str, Any],
) -> dict[str, Any]:
    if identity["facets_executable_sha256"] != original_registration[
        "facets_executable_sha256"
    ]:
        raise ValueError("FACETS executable differs from the frozen registration")
    runs = _read_csv(study_dir / "aggregate" / "run_ledger.csv")
    facets = runs.loc[runs["EstimatorMode"].astype(str).eq(FACETS_MODE)].copy()
    if len(facets) != 36 or facets["AttemptId"].astype(str).nunique() != 36:
        raise ValueError("FACETS pair denominator/key changed")
    for column in ("Converged", "InferenceReady", "IncludedInStudy", "CalibrationReady"):
        if not _strict_bool_series(facets[column], f"FACETS {column}").all():
            raise ValueError(f"FACETS gate is false: {column}")
    if "DirectAgreementPass" in facets and not _strict_bool_series(
        facets["DirectAgreementPass"], "FACETS DirectAgreementPass"
    ).all():
        raise ValueError("A FACETS direct-agreement gate failed")
    _require_finite(
        facets,
        (
            "MainWeightedMAE",
            "MainMaxAbsDifference",
            "MinimumWithinFacetSpearman",
            "ThresholdWeightedMAE",
            "ThresholdMaxAbsDifference",
        ),
        "FACETS calibration",
    )
    assessment = _load_json(study_dir / "aggregate" / "assessment.json")
    calibration = assessment["facets_calibration"]
    attempts = _read_csv(study_dir / "retained_input" / "attempt_manifest.csv")
    pair_attempts = attempts.loc[
        attempts["AttemptType"].astype(str).eq(PAIR_ATTEMPT)
    ]
    if len(pair_attempts) != 36:
        raise ValueError("FACETS pair manifest denominator changed")
    expected_pair_attempt_ids = set(pair_attempts["AttemptId"].astype(str))
    expected_pair_run_ids = set(pair_attempts["RunId"].astype(str))
    if (
        set(facets["AttemptId"].astype(str)) != expected_pair_attempt_ids
        or set(facets["RunId"].astype(str)) != expected_pair_run_ids
    ):
        raise ValueError("FACETS ledger keys differ from the frozen calibration subset")
    outcomes = _read_csv(study_dir / "aggregate" / "attempt_outcomes.csv")
    pair_outcomes = outcomes.loc[
        outcomes["AttemptType"].astype(str).eq(PAIR_ATTEMPT)
    ]
    if (
        len(pair_outcomes) != 36
        or pair_outcomes["AttemptId"].astype(str).nunique() != 36
        or set(pair_outcomes["AttemptId"].astype(str)) != expected_pair_attempt_ids
        or set(pair_outcomes["RunId"].astype(str)) != expected_pair_run_ids
        or set(pair_outcomes["facets_reported_version"].astype(str)) != {"4.5.0"}
    ):
        raise ValueError("FACETS outcome version/denominator is not exactly 36 at 4.5.0")
    audit_facets_contract = plan["strict_integrity_contract"]["facets"]
    maximum_tries = int(audit_facets_contract["maximum_tries"])
    replay_tolerance = float(audit_facets_contract["python_replay_tolerance"])
    total_tries = 0
    for _, attempt in pair_attempts.iterrows():
        ordinal = int(attempt["AttemptOrdinal"])
        work = study_dir / "work" / f"{ordinal:05d}"
        metrics = _load_json(work / "resilient_pair_metrics.json")
        total_tries += _validate_facets_pair_metrics(
            metrics,
            ordinal=ordinal,
            expected_run_id=str(attempt["RunId"]),
            identity=identity,
            work=work,
            maximum_tries=maximum_tries,
            replay_tolerance=replay_tolerance,
        )
    recomputed_calibration = {
        "maximum_main_weighted_mae": float(pd.to_numeric(facets["MainWeightedMAE"]).max()),
        "maximum_main_absolute_difference": float(
            pd.to_numeric(facets["MainMaxAbsDifference"]).max()
        ),
        "maximum_threshold_weighted_mae": float(
            pd.to_numeric(facets["ThresholdWeightedMAE"]).max()
        ),
        "maximum_threshold_absolute_difference": float(
            pd.to_numeric(facets["ThresholdMaxAbsDifference"]).max()
        ),
        "minimum_within_facet_spearman": float(
            pd.to_numeric(facets["MinimumWithinFacetSpearman"]).min()
        ),
    }
    for key, expected in recomputed_calibration.items():
        _require_close(calibration[key], expected, f"FACETS aggregate {key}")
    original_plan = _load_json(ROOT / identity["plan_file"])
    facets_contract = original_plan["facets_contract"]
    registered_gates = {
        "maximum_main_weighted_mae": (
            recomputed_calibration["maximum_main_weighted_mae"],
            float(facets_contract["main_weighted_mae_max"]),
            "maximum",
        ),
        "maximum_main_absolute_difference": (
            recomputed_calibration["maximum_main_absolute_difference"],
            float(facets_contract["main_absolute_difference_max"]),
            "maximum",
        ),
        "maximum_threshold_weighted_mae": (
            recomputed_calibration["maximum_threshold_weighted_mae"],
            float(facets_contract["threshold_weighted_mae_max"]),
            "maximum",
        ),
        "maximum_threshold_absolute_difference": (
            recomputed_calibration["maximum_threshold_absolute_difference"],
            float(facets_contract["threshold_absolute_difference_max"]),
            "maximum",
        ),
        "minimum_within_facet_spearman": (
            recomputed_calibration["minimum_within_facet_spearman"],
            float(facets_contract["minimum_within_facet_spearman"]),
            "minimum",
        ),
    }
    for label, (observed, threshold, direction) in registered_gates.items():
        passed = observed <= threshold if direction == "maximum" else observed >= threshold
        if not passed:
            raise ValueError(f"FACETS registered gate failed: {label}")
    return {
        "facets_executable_sha256": identity["facets_executable_sha256"],
        "facets_reported_version": "4.5.0",
        "attempted_pairs": 36,
        "ready_pairs": 36,
        "failed_pairs": 0,
        "total_facets_tries": total_tries,
        "all_pair_metrics_revalidated": True,
        **recomputed_calibration,
        "aggregate_calibration_summary_recomputed": True,
        "prospectively_registered_calibration_thresholds_rechecked": True,
        "frozen_calibration_run_ids_exact": True,
        "ordinary_two_decimal_fit_used_as_raw_input": False,
    }


def audit_study(study_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    plan, registration = validate_audit_registration()
    expected_study = (ROOT / plan["target"]["study"]).resolve()
    if study_dir != expected_study:
        raise ValueError("Strict audit may inspect only its hash-frozen target study")
    _guard_study_roots(study_dir)
    identity, original_registration = validate_frozen_identity(study_dir, plan)
    completion = validate_completion_bundle(study_dir, identity)
    if completion["completion_marker_set_sha256"] != plan["target"][
        "completion_marker_set_sha256"
    ]:
        raise ValueError("Target completion-marker set digest changed")
    aggregate = validate_aggregate_bundle(study_dir, identity, completion)
    r_bundle = validate_r_bundle(study_dir, original_registration, plan)
    facets = validate_facets_bundle(study_dir, identity, original_registration, plan)
    scientific = plan["scientific_qualification_contract"]
    return {
        "schema_version": SCHEMA_VERSION,
        "audit_nature": plan["audit_nature"],
        "does_not_modify_frozen_v1": True,
        "overall_audit_status": "INTEGRITY_PASS_SCIENTIFIC_WORDING_QUALIFIED",
        "strict_evidence_integrity": {
            "status": "PASS",
            "completion_markers": completion["attempts"],
            "completion_marker_set_sha256": completion[
                "completion_marker_set_sha256"
            ],
            "all_completion_roots_canonical": True,
            "all_completion_file_sets_exact": True,
            "aggregate_file_hashes_exact": aggregate["aggregate_file_hashes_exact"],
            "r_bundle_recomputed": r_bundle["assessment_recomputed"],
            "mml_rater_key_sets_exact": True,
            "facets_registered_executable_and_version_exact": True,
        },
        "original_artifact_status": aggregate["original_status"],
        "posthoc_scientific_qualification": {
            "status": scientific["posthoc_status"],
            "original_result_revoked": False,
            "DirectionResultRetained": scientific["direction_result_retained"],
            "PrecisionResultRetained": scientific["precision_result_retained"],
            "OptimizationStationarityQualification": scientific[
                "optimization_stationarity_qualification"
            ],
            "exact_mle_stationarity_claim_authorized": scientific[
                "exact_mle_stationarity_claim_authorized"
            ],
            "authorized_claim": scientific["authorized_claim"],
            "practical_importance_claim_authorized": scientific[
                "practical_importance_claim_authorized"
            ],
            "estimator_superiority_claim_authorized": scientific[
                "estimator_superiority_claim_authorized"
            ],
        },
        "registered_endpoints": aggregate["endpoints"],
        "failure_denominator": {
            "attempt_units": {
                "planned": 1272,
                "execution_completed": completion["attempts"],
                "statistical_evidence_ready": completion["attempts"],
                "failed": 0,
            },
            "estimator_fit_lanes": aggregate["estimator_lanes"],
            "endpoint_triplets": {
                endpoint_id: {"required": 200, "finite": 200}
                for endpoint_id in ENDPOINT_ALTERNATIVES
            },
            "independent_r_datasets": {
                "planned": 36,
                "returned": 36,
                "converged": 36,
                "inference_ready": 36,
                "failed": 0,
            },
        },
        "mml_key_integrity": aggregate["mml_keys"],
        "stationarity_audit": r_bundle["stationarity_audit"],
        "r_subset_descriptive_sensitivity": r_bundle["r_subset_sensitivity"],
        "r_bundle_integrity": {
            key: value
            for key, value in r_bundle.items()
            if key not in {"stationarity_audit", "r_subset_sensitivity"}
        },
        "facets_calibration": facets,
        "validation_boundaries": {
            "FACETS": {
                "validated": [
                    "same-dataset FACETS JMLE versus Python JMLE numerical calibration",
                    "response/category encoding and main-facet/PCM-threshold alignment on 36 datasets",
                ],
                "not_validated": [
                    "Python MML likelihood or optimizer",
                    "free population-SD estimation",
                    "KA1 confirmatory endpoint",
                    "ordinary fit statistics as raw high-precision values",
                    "generalization beyond the registered DGP",
                ],
            },
            "R": {
                "validated": [
                    "independent base-R likelihood evaluation",
                    "local Q31 reoptimization from the Python solution",
                    "Q31 versus Q61 sensitivity on the frozen 36-dataset subset",
                ],
                "not_validated": [
                    "independent global-basin recovery",
                    "all 600 free-MML datasets",
                    "FACETS parity or the DGP itself",
                ],
            },
            "Python": {
                "produced": [
                    "KA1 through KA4 over 200 triplets",
                    "600 free-MML and 600 fixed-MML fits",
                ],
                "qualification_gap": [
                    "no prospectively frozen Python terminal-score acceptance criterion"
                ],
            },
        },
        "parent_hashes": {
            "study_identity_sha256": sha256_file(study_dir / "study_identity.json"),
            "aggregate_identity_sha256": sha256_file(
                study_dir / "aggregate" / "aggregate_identity.json"
            ),
            "r_verification_identity_sha256": sha256_file(
                study_dir / "r_mml" / "verification_identity.json"
            ),
            "r_assessment_sha256": sha256_file(study_dir / "r_mml" / "assessment.json"),
            "execution_registration_sha256": sha256_file(
                ROOT / identity["registration_file"]
            ),
            "aggregate_assessment_sha256": sha256_file(
                study_dir / "aggregate" / "assessment.json"
            ),
            "completion_marker_set_sha256": completion[
                "completion_marker_set_sha256"
            ],
            "r_process_sha256": sha256_file(study_dir / "r_mml" / "r_process.json"),
        },
        "claim_boundary": plan["claim_boundary"],
        "audit_registration_sha256": sha256_file(REGISTRATION_PATH),
        "auditor_sha256": registration["auditor_sha256"],
    }


def _markdown(result: dict[str, Any]) -> str:
    stationarity = result["stationarity_audit"]
    sensitivity = result["r_subset_descriptive_sensitivity"]["solutions"]
    return f"""# KAC200 append-only strict audit addendum

The hash-frozen v1 artifacts were not modified. Under the prospectively frozen
v1 gates, `ScientificInferenceValidated=true`. This posthoc audit does not
revoke the registered direction or precision result. It does qualify the
stronger optimizer wording: because no prospective acceptance threshold was
frozen for the Python-solution terminal Q31 score,
`OptimizationStationarityQualification=INCOMPLETE`.

No retrospective stationarity pass/fail threshold was selected from the
observed kac200 values.

## Integrity result

- Strict evidence integrity: **{result['strict_evidence_integrity']['status']}**
- Completion markers: {result['strict_evidence_integrity']['completion_markers']}/1272
- Completion roots and artifact file sets: exact
- Aggregate hashes and completion-marker digest: exact
- R raw bundle and assessment recomputation: exact
- FACETS executable/version: registered hash / 4.5.0
- MML Rater keys: R01--R04 exactly once in every one of 1200 cells

## Stationarity qualification

- Python Q31 terminal score datasets: {stationarity['datasets']}
- Minimum / median / maximum sup norm: {stationarity['minimum']:.12g} /
  {stationarity['median']:.12g} / {stationarity['maximum']:.12g}
- Counts above 0.01 / 0.1 (descriptive only):
  {stationarity['descriptive_count_above_0p01']} /
  {stationarity['descriptive_count_above_0p1']}
- Maximum R Q31 / Q61 terminal gradient:
  {stationarity['maximum_r_q31_terminal_gradient']:.12g} /
  {stationarity['maximum_r_q61_terminal_gradient']:.12g}

## Frozen R-subset numerical sensitivity

This is a posthoc descriptive sensitivity, not a new confirmatory test.

- Python KA1 mean: {sensitivity['Python']['mean_contrast']:.12g}
- R Q31 KA1 mean: {sensitivity['R_Q31']['mean_contrast']:.12g}
- R Q61 KA1 mean: {sensitivity['R_Q61']['mean_contrast']:.12g}
- Positive triplets (Python / R Q31 / R Q61):
  {sensitivity['Python']['positive_triplets']} /
  {sensitivity['R_Q31']['positive_triplets']} /
  {sensitivity['R_Q61']['positive_triplets']}

## Authorized wording

The registered algorithm produced a direction-confirmed and precision-qualified
RMSE stress contrast under the frozen DGP. On the frozen 12-triplet subset, the
observed Python, R Q31, and R Q61 contrasts had the reported means above and
12/12 positive signs; no stability pass/fail criterion was assigned.
Exact Python MLE stationarity is not qualified. FACETS calibrated the
same-estimand JMLE lane only; it did not validate MML or KA1. Practical
importance, estimator superiority, and broader-DGP robustness remain
unauthorized claims.
"""


def validate_published_bundle(output_dir: Path, expected_result: dict[str, Any]) -> None:
    output_dir = output_dir.resolve()
    expected_files = {"assessment.json", "AUDIT_ADDENDUM.md", "audit_identity.json"}
    entries = list(output_dir.iterdir())
    if (
        {path.name for path in entries} != expected_files
        or any(_is_link_like(path) or not path.is_file() for path in entries)
    ):
        raise ValueError("Published strict-audit file set changed or is not regular")
    assessment = _load_json(output_dir / "assessment.json")
    if assessment != expected_result:
        raise ValueError("Published strict-audit assessment differs from recomputation")
    if (output_dir / "AUDIT_ADDENDUM.md").read_text(encoding="utf-8") != _markdown(
        expected_result
    ):
        raise ValueError("Published strict-audit Markdown differs from recomputation")
    identity = _load_json(output_dir / "audit_identity.json")
    _require_exact_keys(
        identity,
        {
            "schema_version",
            "parent_hashes",
            "audit_plan_sha256",
            "audit_registration_sha256",
            "auditor_sha256",
            "test_contract_sha256",
            "artifact_sha256",
        },
        "published strict-audit identity",
    )
    expected_identity_fields = {
        "schema_version": f"{SCHEMA_VERSION}_identity_v1",
        "parent_hashes": expected_result["parent_hashes"],
        "audit_plan_sha256": sha256_file(PLAN_PATH),
        "audit_registration_sha256": sha256_file(REGISTRATION_PATH),
        "auditor_sha256": sha256_file(Path(__file__).resolve()),
        "test_contract_sha256": sha256_file(TEST_PATH),
    }
    for key, expected in expected_identity_fields.items():
        if identity[key] != expected:
            raise ValueError(f"Published strict-audit identity mismatch: {key}")
    expected_artifacts = {
        filename: sha256_file(output_dir / filename)
        for filename in ("assessment.json", "AUDIT_ADDENDUM.md")
    }
    if identity["artifact_sha256"] != expected_artifacts:
        raise ValueError("Published strict-audit artifact hashes changed")


def publish_audit(study_dir: Path, output_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    output_dir = output_dir.resolve()
    planned_output = (ROOT / _load_json(PLAN_PATH)["immutability"]["output"]).resolve()
    if output_dir != planned_output:
        raise ValueError("Strict audit may publish only to its planned output path")
    if _inside(output_dir, study_dir):
        raise ValueError("Strict audit output must remain outside the frozen study")
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite strict audit: {output_dir}")
    result = audit_study(study_dir)
    partial = output_dir.with_name(f"{output_dir.name}.partial.{os.getpid()}.{time.time_ns()}")
    if _inside(partial, study_dir):
        raise ValueError("Strict-audit partial output would enter the frozen study")
    partial.mkdir(parents=True, exist_ok=False)
    try:
        _write_json(partial / "assessment.json", result)
        (partial / "AUDIT_ADDENDUM.md").write_text(_markdown(result), encoding="utf-8")
        identity = {
            "schema_version": f"{SCHEMA_VERSION}_identity_v1",
            "parent_hashes": result["parent_hashes"],
            "audit_plan_sha256": sha256_file(PLAN_PATH),
            "audit_registration_sha256": sha256_file(REGISTRATION_PATH),
            "auditor_sha256": sha256_file(Path(__file__).resolve()),
            "test_contract_sha256": sha256_file(TEST_PATH),
            "artifact_sha256": {
                filename: sha256_file(partial / filename)
                for filename in ("assessment.json", "AUDIT_ADDENDUM.md")
            },
        }
        _write_json(partial / "audit_identity.json", identity)
        validate_published_bundle(partial, result)
        os.replace(partial, output_dir)
        validate_published_bundle(output_dir, result)
    except BaseException:
        if partial.exists():
            shutil.rmtree(partial, ignore_errors=True)
        raise
    return result


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-dir", type=Path, default=DEFAULT_STUDY_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--no-publish", action="store_true")
    action.add_argument("--verify-only", action="store_true")
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    if args.verify_only:
        result = audit_study(args.study_dir)
        validate_published_bundle(args.output_dir, result)
    elif args.no_publish:
        result = audit_study(args.study_dir)
    else:
        result = publish_audit(args.study_dir, args.output_dir)
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
