#!/usr/bin/env python3
"""Versioned resilient FACETS/Python JMLE pair component.

The independent Python artifact is completion-marked before FACETS is called.
FACETS calibration may use a narrowly bounded retry and a cross-process lock;
calibration failure never deletes or de-qualifies a ready Python artifact.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import sys
import threading
import time
from typing import Any, Callable, Iterator

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from validation.estimand_bridge_pilot import _normalize_recovery, _normalize_steps  # noqa: E402
from validation.estimand_distribution_study import (  # noqa: E402
    FACETS_MODE,
    FIT_MODEL,
    PYTHON_JMLE_MODE,
    THRESHOLD_CONDITION,
    _jmle_outputs,
)
from validation.facets_pcm_boundary_pilot import (  # noqa: E402
    _fit_facets_and_python,
    python_fit,
)
from validation.operating_characteristics_facets import sha256_file  # noqa: E402


SCHEMA_VERSION = "mfrm-facets-resilient-pair-v1"
PLAN_PATH = REPO_ROOT / "validation" / "facets_resilient_pair_plan_20260811.json"
PYTHON_REPLAY_TOLERANCE = 1e-12
DEFAULT_MAX_RETRIES = 2
PairFitCallable = Callable[..., tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]]

_THREAD_LOCKS: dict[str, threading.Lock] = {}
_THREAD_LOCKS_GUARD = threading.Lock()


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}.{time.time_ns()}")
    _json_dump(temporary, value)
    os.replace(temporary, path)


def _frame_digest(frame: pd.DataFrame) -> str:
    payload = frame.to_csv(index=False, lineterminator="\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def dependency_manifest(repo_root: Path = REPO_ROOT) -> dict[str, str]:
    """Hash direct validation code and the full local application Python surface."""

    repo_root = repo_root.resolve()
    paths = set((repo_root / "mfrm_app").rglob("*.py"))
    paths.update({
        repo_root / "streamlit_app.py",
        repo_root / "requirements.txt",
        repo_root / "requirements-dev.txt",
        repo_root / "validation" / "estimand_bridge_pilot.py",
        repo_root / "validation" / "estimand_distribution_study.py",
        repo_root / "validation" / "facets_pcm_boundary_pilot.py",
        repo_root / "validation" / "facets_pcm_known_truth_smoke.py",
        repo_root / "validation" / "facets_visible_launcher.py",
        repo_root / "validation" / "operating_characteristics_facets.py",
        repo_root / "validation" / "operating_characteristics_python_mml.py",
        Path(__file__).resolve(),
        PLAN_PATH.resolve(),
    })
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Dependency manifest files missing: {missing}")
    return {
        path.relative_to(repo_root).as_posix(): sha256_file(path)
        for path in sorted(paths, key=lambda item: item.as_posix().lower())
    }


def dependency_manifest_digest(manifest: dict[str, str]) -> str:
    return hashlib.sha256(
        json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _thread_lock(path: Path) -> threading.Lock:
    key = str(path.resolve()).lower()
    with _THREAD_LOCKS_GUARD:
        return _THREAD_LOCKS.setdefault(key, threading.Lock())


@contextmanager
def cross_process_file_lock(
    path: Path,
    *,
    timeout_seconds: float = 300.0,
    poll_seconds: float = 0.05,
) -> Iterator[None]:
    """Serialize FACETS calls across threads and processes on Windows/POSIX."""

    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    local_lock = _thread_lock(path)
    if not local_lock.acquire(timeout=timeout_seconds):
        raise TimeoutError(f"Timed out waiting for thread lock: {path}")
    handle = None
    acquired = False
    try:
        handle = path.open("a+b")
        handle.seek(0, os.SEEK_END)
        if handle.tell() == 0:
            handle.write(b"0")
            handle.flush()
        deadline = time.monotonic() + timeout_seconds
        while not acquired:
            try:
                handle.seek(0)
                if os.name == "nt":
                    import msvcrt  # pylint: disable=import-outside-toplevel

                    msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                else:
                    import fcntl  # type: ignore[import-not-found]  # pylint: disable=import-outside-toplevel

                    fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                acquired = True
            except (OSError, BlockingIOError):
                if time.monotonic() >= deadline:
                    raise TimeoutError(f"Timed out waiting for process lock: {path}")
                time.sleep(poll_seconds)
        yield
    finally:
        if handle is not None:
            if acquired:
                handle.seek(0)
                if os.name == "nt":
                    import msvcrt  # pylint: disable=import-outside-toplevel

                    msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
                else:
                    import fcntl  # type: ignore[import-not-found]  # pylint: disable=import-outside-toplevel

                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            handle.close()
        local_lock.release()


def retryable_missing_report_failure(error: BaseException, try_dir: Path) -> bool:
    """Allow retry only for FACETS exit-0 plus the named missing report."""

    message = f"{type(error).__name__}: {error}"
    match = re.search(r"FACETS (primary|auxiliary) failed .*exit=0", message)
    if match is None:
        return False
    expected = "report_u6.txt" if match.group(1) == "primary" else "report_u2.txt"
    return not any(try_dir.rglob(expected))


@dataclass
class RetryExecution:
    result: Any | None
    records: list[dict[str, Any]]
    final_error: str


def execute_with_bounded_retry(
    operation: Callable[[Path, int], Any],
    *,
    attempts_root: Path,
    maximum_retries: int = DEFAULT_MAX_RETRIES,
    lock_path: Path | None = None,
    lock_timeout_seconds: float = 300.0,
) -> RetryExecution:
    if maximum_retries < 0:
        raise ValueError("maximum_retries must be nonnegative")
    attempts_root.mkdir(parents=True, exist_ok=False)
    records: list[dict[str, Any]] = []
    final_error = ""
    for try_number in range(1, maximum_retries + 2):
        try_dir = attempts_root / f"pair_try_{try_number:02d}"
        try_dir.mkdir()
        started = time.perf_counter()
        try:
            if lock_path is None:
                result = operation(try_dir, try_number)
            else:
                with cross_process_file_lock(
                    lock_path, timeout_seconds=lock_timeout_seconds
                ):
                    result = operation(try_dir, try_number)
            record = {
                "Try": try_number,
                "Succeeded": True,
                "RetryEligible": False,
                "FailureReason": "",
                "ElapsedSeconds": time.perf_counter() - started,
            }
            records.append(record)
            _atomic_json(try_dir / "try_outcome.json", record)
            return RetryExecution(result=result, records=records, final_error="")
        except Exception as exc:  # retained operational denominator
            retryable = retryable_missing_report_failure(exc, try_dir)
            final_error = f"{type(exc).__name__}: {exc}"
            record = {
                "Try": try_number,
                "Succeeded": False,
                "RetryEligible": retryable,
                "FailureReason": final_error,
                "ElapsedSeconds": time.perf_counter() - started,
            }
            records.append(record)
            _atomic_json(try_dir / "try_outcome.json", record)
            if not retryable or try_number > maximum_retries:
                break
    return RetryExecution(result=None, records=records, final_error=final_error)


def _fit_python_first(
    manifest_row: pd.Series,
    ratings: pd.DataFrame,
    truth: pd.DataFrame,
    threshold_truth: pd.DataFrame,
    *,
    output_dir: Path,
    dependency_digest: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    result = python_fit(ratings, FIT_MODEL)
    summary = result["summary"].iloc[0]
    estimates = result["facets"]["others"].copy()
    estimates["SE"] = np.nan
    estimates["Status"] = np.nan
    thresholds = _normalize_steps(
        result,
        fit_model=FIT_MODEL,
        estimator_mode=PYTHON_JMLE_MODE,
        run_id=str(manifest_row["RunId"]),
        threshold_condition=THRESHOLD_CONDITION,
        threshold_truth=threshold_truth,
    )
    finite_main = pd.to_numeric(estimates["Estimate"], errors="coerce").notna().all()
    finite_thresholds = pd.to_numeric(thresholds["Estimate"], errors="coerce").notna().all()
    ready = bool(
        summary.get("Converged", False)
        and summary.get("InferenceReady", False)
        and finite_main
        and finite_thresholds
    )
    recovery = _normalize_recovery(
        estimates,
        truth,
        manifest_row,
        estimator_family="JMLE",
        estimator_mode=PYTHON_JMLE_MODE,
        fit_model=FIT_MODEL,
        included=ready,
    )
    recovery["PersonDistribution"] = manifest_row["PersonDistribution"]
    recovery["IncludedInStudy"] = recovery.pop("IncludedInBridge")
    thresholds["ConditionId"] = manifest_row["ConditionId"]
    thresholds["Design"] = manifest_row["Design"]
    thresholds["PersonDistribution"] = manifest_row["PersonDistribution"]
    thresholds["Replicate"] = int(manifest_row["Replicate"])
    thresholds["EstimatorFamily"] = "JMLE"
    thresholds["EstimandClass"] = "fixed_person_joint_likelihood"
    thresholds["Engine"] = "PythonApp"
    thresholds["IncludedInStudy"] = ready & thresholds["TruthError"].notna()
    thresholds["DisplayPrecisionContract"] = "native numeric estimate"
    run = pd.DataFrame([{
        "RunId": manifest_row["RunId"],
        "ConditionId": manifest_row["ConditionId"],
        "Design": manifest_row["Design"],
        "PersonDistribution": manifest_row["PersonDistribution"],
        "Replicate": int(manifest_row["Replicate"]),
        "EstimatorFamily": "JMLE",
        "EstimatorMode": PYTHON_JMLE_MODE,
        "EstimandClass": "fixed_person_joint_likelihood",
        "LikelihoodBasis": "joint",
        "Engine": "PythonApp",
        "Rows": int(len(ratings)),
        "FitReturned": True,
        "Converged": bool(summary.get("Converged", False)),
        "InferenceReady": bool(summary.get("InferenceReady", False)),
        "IncludedInStudy": ready,
        "StatisticalEvidenceReady": ready,
        "CalibrationReady": False,
        "PairFullyQualified": False,
        "ElapsedSeconds": time.perf_counter() - started,
        "FailureReason": "" if ready else "Independent Python JMLE readiness failed",
    }])
    artifacts = {
        "run_ledger.csv": run,
        "recovery.csv": recovery,
        "thresholds.csv": thresholds,
    }
    for filename, frame in artifacts.items():
        frame.to_csv(output_dir / filename, index=False, lineterminator="\n")
    completion = {
        "schema_version": SCHEMA_VERSION,
        "run_id": str(manifest_row["RunId"]),
        "dependency_manifest_sha256": dependency_digest,
        "statistical_evidence_ready": ready,
        "artifact_sha256": {
            filename: sha256_file(output_dir / filename) for filename in artifacts
        },
    }
    _atomic_json(output_dir / "python_completion.json", completion)
    return run, recovery, thresholds


def _python_replay_differences(
    independent_recovery: pd.DataFrame,
    independent_thresholds: pd.DataFrame,
    pair_recovery: pd.DataFrame,
    pair_thresholds: pd.DataFrame,
) -> tuple[float, float]:
    pair_python = pair_recovery[pair_recovery["EstimatorMode"].eq(PYTHON_JMLE_MODE)]
    main = independent_recovery[["Facet", "Level", "EstimateAligned"]].merge(
        pair_python[["Facet", "Level", "EstimateAligned"]],
        on=["Facet", "Level"], suffixes=("Independent", "Pair"), validate="one_to_one"
    )
    main_max = float(
        (main["EstimateAlignedIndependent"] - main["EstimateAlignedPair"]).abs().max()
    )
    pair_python_thresholds = pair_thresholds[
        pair_thresholds["EstimatorMode"].eq(PYTHON_JMLE_MODE)
    ]
    threshold = independent_thresholds[["StepFacetLevel", "Category", "Estimate"]].merge(
        pair_python_thresholds[["StepFacetLevel", "Category", "Estimate"]],
        on=["StepFacetLevel", "Category"],
        suffixes=("Independent", "Pair"),
        validate="one_to_one",
    )
    threshold_max = float(
        (threshold["EstimateIndependent"] - threshold["EstimatePair"]).abs().max()
    )
    return main_max, threshold_max


def fit_resilient_pair(
    manifest_row: pd.Series,
    ratings: pd.DataFrame,
    truth: pd.DataFrame,
    threshold_truth: pd.DataFrame,
    *,
    facets_exe: Path,
    output_dir: Path,
    lock_path: Path,
    timeout_seconds: float = 120.0,
    maximum_retries: int = DEFAULT_MAX_RETRIES,
    pair_fit: PairFitCallable = _fit_facets_and_python,
) -> dict[str, Any]:
    """Persist Python first, then obtain a separately gated FACETS calibration."""

    output_dir = output_dir.resolve()
    facets_exe = facets_exe.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Resilient pair output already exists: {output_dir}")
    if not facets_exe.is_file():
        raise FileNotFoundError(f"FACETS executable missing: {facets_exe}")
    output_dir.mkdir(parents=True)
    manifest = dependency_manifest()
    manifest_digest = dependency_manifest_digest(manifest)
    _json_dump(output_dir / "dependency_manifest.json", manifest)

    python_run, python_recovery, python_thresholds = _fit_python_first(
        manifest_row,
        ratings,
        truth,
        threshold_truth,
        output_dir=output_dir / "python_independent",
        dependency_digest=manifest_digest,
    )
    python_ready = bool(python_run.iloc[0]["StatisticalEvidenceReady"])

    def operation(try_dir: Path, _try_number: int) -> Any:
        work_root = try_dir / "facets_runs"
        work_root.mkdir()
        return pair_fit(
            manifest_row,
            ratings,
            truth,
            threshold_truth,
            fit_model=FIT_MODEL,
            facets_exe=facets_exe,
            work_root=work_root,
            timeout_seconds=timeout_seconds,
        )

    retry = execute_with_bounded_retry(
        operation,
        attempts_root=output_dir / "facets_attempts",
        maximum_retries=maximum_retries,
        lock_path=lock_path,
        lock_timeout_seconds=max(300.0, timeout_seconds * 3),
    )
    calibration_ready = False
    replay_main_max = np.nan
    replay_threshold_max = np.nan
    facets_run = pd.DataFrame([{
        "RunId": manifest_row["RunId"],
        "ConditionId": manifest_row["ConditionId"],
        "Design": manifest_row["Design"],
        "PersonDistribution": manifest_row["PersonDistribution"],
        "Replicate": int(manifest_row["Replicate"]),
        "EstimatorFamily": "JMLE",
        "EstimatorMode": FACETS_MODE,
        "EstimandClass": "fixed_person_joint_likelihood",
        "LikelihoodBasis": "joint",
        "Engine": "FACETS",
        "FitReturned": False,
        "Converged": False,
        "InferenceReady": False,
        "IncludedInStudy": False,
        "StatisticalEvidenceReady": False,
        "CalibrationReady": False,
        "PairFullyQualified": False,
        "FailureReason": retry.final_error,
    }])
    facets_recovery = pd.DataFrame()
    facets_thresholds = pd.DataFrame()
    pair_metrics: dict[str, Any] = {}
    if retry.result is not None:
        metrics, main, pair_threshold_pairs, _, _ = retry.result
        pair_runs, pair_recovery, pair_thresholds = _jmle_outputs(
            manifest_row, metrics, main, pair_threshold_pairs
        )
        replay_main_max, replay_threshold_max = _python_replay_differences(
            python_recovery, python_thresholds, pair_recovery, pair_thresholds
        )
        replay_pass = bool(
            replay_main_max <= PYTHON_REPLAY_TOLERANCE
            and replay_threshold_max <= PYTHON_REPLAY_TOLERANCE
        )
        calibration_ready = bool(metrics["DirectAgreementPass"] and replay_pass)
        facets_run = pair_runs[pair_runs["EstimatorMode"].eq(FACETS_MODE)].copy()
        facets_run["IncludedInStudy"] = calibration_ready
        facets_run["StatisticalEvidenceReady"] = False
        facets_run["CalibrationReady"] = calibration_ready
        facets_run["PairFullyQualified"] = python_ready and calibration_ready
        facets_run["FailureReason"] = "" if calibration_ready else "Calibration or Python replay gate failed"
        facets_recovery = pair_recovery[pair_recovery["EstimatorMode"].eq(FACETS_MODE)].copy()
        facets_recovery["IncludedInStudy"] = calibration_ready
        facets_thresholds = pair_thresholds[pair_thresholds["EstimatorMode"].eq(FACETS_MODE)].copy()
        facets_thresholds["IncludedInStudy"] = calibration_ready
        pair_metrics = metrics

    python_run["CalibrationReady"] = calibration_ready
    python_run["PairFullyQualified"] = python_ready and calibration_ready
    python_recovery["CalibrationReady"] = calibration_ready
    python_thresholds["CalibrationReady"] = calibration_ready
    combined_runs = pd.concat([facets_run, python_run], ignore_index=True, sort=False)
    combined_recovery = pd.concat(
        [facets_recovery, python_recovery], ignore_index=True, sort=False
    )
    combined_thresholds = pd.concat(
        [facets_thresholds, python_thresholds], ignore_index=True, sort=False
    )
    combined_runs.to_csv(output_dir / "combined_run_ledger.csv", index=False, lineterminator="\n")
    combined_recovery.to_csv(output_dir / "combined_recovery.csv", index=False, lineterminator="\n")
    combined_thresholds.to_csv(output_dir / "combined_thresholds.csv", index=False, lineterminator="\n")
    pd.DataFrame(retry.records).to_csv(
        output_dir / "facets_try_ledger.csv", index=False, lineterminator="\n"
    )
    outcome = {
        "schema_version": SCHEMA_VERSION,
        "run_id": str(manifest_row["RunId"]),
        "statistical_evidence_ready": python_ready,
        "calibration_ready": calibration_ready,
        "pair_fully_qualified": bool(python_ready and calibration_ready),
        "facets_tries": len(retry.records),
        "facets_retry_count": max(0, len(retry.records) - 1),
        "facets_final_error": retry.final_error,
        "python_replay_main_max_abs_difference": replay_main_max,
        "python_replay_threshold_max_abs_difference": replay_threshold_max,
        "python_replay_tolerance": PYTHON_REPLAY_TOLERANCE,
        "direct_agreement_pass": bool(pair_metrics.get("DirectAgreementPass", False)),
        "dependency_manifest_sha256": manifest_digest,
        "facets_executable_sha256": sha256_file(facets_exe),
        "facets_executable_size_bytes": facets_exe.stat().st_size,
        "facets_reported_version": pair_metrics.get("FACETSVersion"),
        "lock_scope": "legacy_pair_call_including_recomputed_python",
        "confirmatory_evidence_replaced": False,
    }
    _json_dump(output_dir / "resilient_pair_metrics.json", outcome)
    artifact_files = (
        "combined_run_ledger.csv",
        "combined_recovery.csv",
        "combined_thresholds.csv",
        "facets_try_ledger.csv",
        "resilient_pair_metrics.json",
        "dependency_manifest.json",
    )
    identity = {
        "schema_version": f"{SCHEMA_VERSION}-identity",
        "plan_sha256": sha256_file(PLAN_PATH),
        "component_sha256": sha256_file(Path(__file__).resolve()),
        "dependency_manifest_sha256": manifest_digest,
        "facets_executable_sha256": sha256_file(facets_exe),
        "python_version": sys.version,
        "platform": platform.platform(),
        "artifact_sha256": {
            filename: sha256_file(output_dir / filename) for filename in artifact_files
        },
    }
    _atomic_json(output_dir / "resilient_pair_identity.json", identity)
    return outcome


__all__ = [
    "DEFAULT_MAX_RETRIES",
    "PYTHON_REPLAY_TOLERANCE",
    "RetryExecution",
    "cross_process_file_lock",
    "dependency_manifest",
    "dependency_manifest_digest",
    "execute_with_bounded_retry",
    "fit_resilient_pair",
    "retryable_missing_report_failure",
]
