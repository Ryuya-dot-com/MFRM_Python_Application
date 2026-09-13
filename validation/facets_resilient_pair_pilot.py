#!/usr/bin/env python3
"""Run the registered resilient FACETS/Python pair fault-injection pilot."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import sys
import threading
import time
from typing import Any

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from validation.estimand_distribution_study import validate_study_identity  # noqa: E402
from validation.facets_pcm_boundary_pilot import _fit_facets_and_python  # noqa: E402
from validation.facets_resilient_pair import (  # noqa: E402
    PLAN_PATH,
    cross_process_file_lock,
    fit_resilient_pair,
)
from validation.operating_characteristics_facets import sha256_file, validate_bundle  # noqa: E402


TARGET_ORDINAL = 1024
ORIGINAL_MARKER_SHA256 = "6451a1075cf7bde124338ce12e6d44456b415cc8a41f7856bcc8aa3e455884fa"


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _load_target(study_dir: Path) -> tuple[pd.Series, pd.Series, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    validate_study_identity(study_dir)
    input_dir = study_dir / "retained_input"
    tables = validate_bundle(input_dir)
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    selected = attempts[
        pd.to_numeric(attempts["AttemptOrdinal"], errors="raise").astype(int).eq(TARGET_ORDINAL)
    ]
    if len(selected) != 1:
        raise ValueError("Registered target attempt is not unique")
    attempt = selected.iloc[0]
    run_id = str(attempt["RunId"])
    manifest = tables["manifest.csv"].set_index("RunId", drop=False)
    ratings = tables["generated_ratings.csv"]
    ratings = ratings[ratings["RunId"].astype(str).eq(run_id)].drop(columns="RunId")
    truth = tables["generated_facet_truth.csv"]
    truth = truth[truth["RunId"].astype(str).eq(run_id)]
    thresholds = pd.read_csv(input_dir / "generated_pcm_threshold_truth.csv")
    thresholds = thresholds[thresholds["RunId"].astype(str).eq(run_id)]
    return attempt, manifest.loc[run_id], ratings, truth, thresholds


def _lock_stress(lock_path: Path, output_path: Path) -> dict[str, Any]:
    active = 0
    maximum_active = 0
    guard = threading.Lock()
    rows: list[dict[str, Any]] = []

    def worker(index: int) -> None:
        nonlocal active, maximum_active
        requested = time.perf_counter()
        with cross_process_file_lock(lock_path, timeout_seconds=30):
            entered = time.perf_counter()
            with guard:
                active += 1
                maximum_active = max(maximum_active, active)
            time.sleep(0.10)
            with guard:
                active -= 1
            exited = time.perf_counter()
        rows.append({
            "Worker": index,
            "Requested": requested,
            "Entered": entered,
            "Exited": exited,
            "HoldSeconds": exited - entered,
        })

    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(worker, range(4)))
    pd.DataFrame(rows).sort_values("Entered").to_csv(
        output_path, index=False, lineterminator="\n"
    )
    return {
        "workers": 4,
        "maximum_concurrent_holders": maximum_active,
        "serialization_pass": maximum_active == 1,
    }


def run_pilot(*, study_dir: Path, facets_exe: Path, output_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    facets_exe = facets_exe.resolve()
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Pilot output already exists: {output_dir}")
    if not facets_exe.is_file():
        raise FileNotFoundError(f"FACETS executable missing: {facets_exe}")
    original_marker = study_dir / "attempts" / f"{TARGET_ORDINAL:05d}" / "completion.json"
    if sha256_file(original_marker) != ORIGINAL_MARKER_SHA256:
        raise ValueError("Frozen failed marker changed before resilient pilot")
    _, manifest_row, ratings, truth, thresholds = _load_target(study_dir)
    output_dir.mkdir(parents=True)
    lock_path = output_dir / "locks" / "facets_global.lock"

    first_call = {"count": 0}

    def fail_first_then_real(*args: Any, **kwargs: Any) -> Any:
        first_call["count"] += 1
        if first_call["count"] == 1:
            raise RuntimeError(
                f"FACETS primary failed for {manifest_row['RunId']}::PCM: exit=0"
            )
        return _fit_facets_and_python(*args, **kwargs)

    fi01 = fit_resilient_pair(
        manifest_row,
        ratings,
        truth,
        thresholds,
        facets_exe=facets_exe,
        output_dir=output_dir / "FI01_RETRYABLE_FIRST_FAILURE",
        lock_path=lock_path,
        pair_fit=fail_first_then_real,
    )

    def always_retryable(*_args: Any, **_kwargs: Any) -> Any:
        raise RuntimeError(
            f"FACETS primary failed for {manifest_row['RunId']}::PCM: exit=0"
        )

    fi02 = fit_resilient_pair(
        manifest_row,
        ratings,
        truth,
        thresholds,
        facets_exe=facets_exe,
        output_dir=output_dir / "FI02_RETRYABLE_EXHAUSTION",
        lock_path=lock_path,
        pair_fit=always_retryable,
    )

    def nonretryable(*_args: Any, **_kwargs: Any) -> Any:
        raise RuntimeError(
            f"FACETS primary failed for {manifest_row['RunId']}::PCM: exit=5"
        )

    fi03 = fit_resilient_pair(
        manifest_row,
        ratings,
        truth,
        thresholds,
        facets_exe=facets_exe,
        output_dir=output_dir / "FI03_NONRETRYABLE_FAILURE",
        lock_path=lock_path,
        pair_fit=nonretryable,
    )
    fi04 = _lock_stress(lock_path, output_dir / "FI04_lock_timeline.csv")

    scenario_rows = [
        {"Scenario": "FI01_RETRYABLE_FIRST_FAILURE", **fi01},
        {"Scenario": "FI02_RETRYABLE_EXHAUSTION", **fi02},
        {"Scenario": "FI03_NONRETRYABLE_FAILURE", **fi03},
        {"Scenario": "FI04_LOCK_SERIALIZATION", **fi04},
    ]
    pd.DataFrame(scenario_rows).to_csv(
        output_dir / "scenario_summary.csv", index=False, lineterminator="\n"
    )
    checks = {
        "fi01_python_ready": fi01["statistical_evidence_ready"] is True,
        "fi01_calibration_ready": fi01["calibration_ready"] is True,
        "fi01_two_tries": fi01["facets_tries"] == 2,
        "fi01_retry_count_one": fi01["facets_retry_count"] == 1,
        "fi01_python_replay": (
            fi01["python_replay_main_max_abs_difference"] <= 1e-12
            and fi01["python_replay_threshold_max_abs_difference"] <= 1e-12
        ),
        "fi02_python_survives": fi02["statistical_evidence_ready"] is True,
        "fi02_calibration_false": fi02["calibration_ready"] is False,
        "fi02_exactly_three_tries": fi02["facets_tries"] == 3,
        "fi03_python_survives": fi03["statistical_evidence_ready"] is True,
        "fi03_calibration_false": fi03["calibration_ready"] is False,
        "fi03_no_retry": fi03["facets_tries"] == 1,
        "fi04_serialization": fi04["serialization_pass"] is True,
        "original_marker_unchanged": sha256_file(original_marker) == ORIGINAL_MARKER_SHA256,
    }
    metrics = {
        "schema_version": "mfrm-facets-resilient-pair-pilot-v1",
        "target_run_id": str(manifest_row["RunId"]),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "confirmatory_evidence_replaced": False,
        "promotion_gate_pass": all(checks.values()),
        "claim_limit": "Operational remediation only; no confirmatory row replacement or estimator ranking.",
    }
    _json_dump(output_dir / "pilot_metrics.json", metrics)
    artifact_hashes = {
        path.relative_to(output_dir).as_posix(): sha256_file(path)
        for path in sorted(output_dir.rglob("*"))
        if path.is_file() and path.name != "pilot_identity.json"
    }
    identity = {
        "schema_version": "mfrm-facets-resilient-pair-pilot-identity-v1",
        "plan_sha256": sha256_file(PLAN_PATH),
        "pilot_script_sha256": sha256_file(Path(__file__).resolve()),
        "component_sha256": sha256_file(REPO_ROOT / "validation" / "facets_resilient_pair.py"),
        "study_identity_sha256": sha256_file(study_dir / "study_identity.json"),
        "original_completion_marker_sha256": sha256_file(original_marker),
        "facets_executable_sha256": sha256_file(facets_exe),
        "artifact_sha256": artifact_hashes,
        "artifact_hashes_sha256": hashlib.sha256(
            json.dumps(artifact_hashes, sort_keys=True).encode("utf-8")
        ).hexdigest(),
    }
    _json_dump(output_dir / "pilot_identity.json", identity)
    print(json.dumps(metrics, indent=2, sort_keys=True))
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-dir", type=Path, required=True)
    parser.add_argument("--facets-exe", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    run_pilot(
        study_dir=args.study_dir,
        facets_exe=args.facets_exe,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
