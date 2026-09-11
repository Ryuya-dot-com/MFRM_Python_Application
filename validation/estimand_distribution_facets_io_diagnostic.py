#!/usr/bin/env python3
"""Serially replay the frozen FACETS report-I/O failure outside the evidence set."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from validation.estimand_distribution_study import (
    _run_one_attempt,
    validate_study_identity,
)
from validation.operating_characteristics_facets import sha256_file, validate_bundle


PLAN_PATH = REPO_ROOT / "validation" / "estimand_distribution_facets_io_diagnostic_plan_20260811.json"


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run_diagnostic(*, study_dir: Path, facets_exe: Path, output_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Diagnostic output already exists: {output_dir}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    target = plan["target"]
    replay = plan["replay"]
    validate_study_identity(study_dir)
    if not facets_exe.is_file():
        raise FileNotFoundError(f"FACETS executable missing: {facets_exe}")
    input_dir = study_dir / "retained_input"
    tables = validate_bundle(input_dir)
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    selected = attempts[
        pd.to_numeric(attempts["AttemptOrdinal"], errors="raise").astype(int).eq(
            int(target["attempt_ordinal"])
        )
    ]
    if len(selected) != 1:
        raise ValueError("Registered diagnostic attempt ordinal is not unique")
    attempt = selected.iloc[0]
    if (
        str(attempt["AttemptId"]) != str(target["attempt_id"])
        or str(attempt["AttemptFingerprint"]) != str(target["attempt_fingerprint"])
    ):
        raise ValueError("Registered diagnostic target identity changed")
    original_marker = study_dir / "attempts" / f"{int(target['attempt_ordinal']):05d}" / "completion.json"
    if sha256_file(original_marker) != str(target["original_completion_marker_sha256"]):
        raise ValueError("Original failed completion marker changed before diagnostic")

    manifest = tables["manifest.csv"].set_index("RunId", drop=False)
    ratings_all = tables["generated_ratings.csv"]
    truth_all = tables["generated_facet_truth.csv"]
    thresholds_all = pd.read_csv(input_dir / "generated_pcm_threshold_truth.csv")
    run_id = str(attempt["RunId"])
    manifest_row = manifest.loc[run_id]
    ratings = ratings_all[ratings_all["RunId"].astype(str).eq(run_id)].drop(columns="RunId")
    truth = truth_all[truth_all["RunId"].astype(str).eq(run_id)]
    threshold_truth = thresholds_all[thresholds_all["RunId"].astype(str).eq(run_id)]

    output_dir.mkdir(parents=True)
    summary_rows: list[dict[str, Any]] = []
    for repeat in range(1, int(replay["serial_repeats"]) + 1):
        repeat_dir = output_dir / f"serial_{repeat:02d}"
        repeat_dir.mkdir()
        succeeded = True
        error = ""
        try:
            runs, recovery, thresholds, constraints = _run_one_attempt(
                attempt,
                manifest_row,
                ratings,
                truth,
                threshold_truth,
                facets_exe=facets_exe,
                stage_dir=repeat_dir,
                timeout_seconds=float(replay["timeout_seconds"]),
            )
        except Exception as exc:
            succeeded = False
            error = f"{type(exc).__name__}: {exc}"
            runs = recovery = thresholds = pd.DataFrame()
            constraints = []
        runs.to_csv(repeat_dir / "run_ledger.csv", index=False, lineterminator="\n")
        recovery.to_csv(repeat_dir / "recovery.csv", index=False, lineterminator="\n")
        thresholds.to_csv(repeat_dir / "thresholds.csv", index=False, lineterminator="\n")
        pd.DataFrame(constraints).to_csv(
            repeat_dir / "constraints.csv", index=False, lineterminator="\n"
        )
        facets_row = runs[runs.get("EstimatorMode", pd.Series(dtype=object)).eq("FACETS_4_5_JMLE")]
        summary_rows.append({
            "Repeat": repeat,
            "Succeeded": succeeded,
            "FailureReason": error,
            "FACETSRows": int(len(facets_row)),
            "DirectAgreementPass": bool(
                len(facets_row) == 1
                and facets_row["DirectAgreementPass"].fillna(False).astype(bool).all()
            ),
            "MainMaxAbsDifference": (
                float(facets_row.iloc[0]["MainMaxAbsDifference"])
                if len(facets_row) == 1 else None
            ),
            "ThresholdMaxAbsDifference": (
                float(facets_row.iloc[0]["ThresholdMaxAbsDifference"])
                if len(facets_row) == 1 else None
            ),
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(output_dir / "diagnostic_summary.csv", index=False, lineterminator="\n")
    marker_unchanged = sha256_file(original_marker) == str(target["original_completion_marker_sha256"])
    output_files = sorted(
        path for path in output_dir.rglob("*")
        if path.is_file() and path.name != "diagnostic_identity.json"
    )
    artifact_hashes = {
        path.relative_to(output_dir).as_posix(): sha256_file(path) for path in output_files
    }
    metrics = {
        "schema_version": "mfrm-estimand-distribution-facets-io-diagnostic-v1",
        "serial_repeats": int(len(summary)),
        "successful_repeats": int(summary["Succeeded"].sum()),
        "direct_agreement_repeats": int(summary["DirectAgreementPass"].sum()),
        "all_serial_repeats_succeeded": bool(summary["Succeeded"].all()),
        "original_completion_marker_unchanged": marker_unchanged,
        "confirmatory_evidence_replaced": False,
        "classification": (
            "transient operational/report-I-O failure supported; concurrency not proven"
            if summary["Succeeded"].all()
            else "serially reproducible FACETS failure"
        ),
    }
    _json_dump(output_dir / "diagnostic_metrics.json", metrics)
    identity = {
        "schema_version": "mfrm-estimand-distribution-facets-io-diagnostic-identity-v1",
        "plan_sha256": sha256_file(PLAN_PATH),
        "diagnostic_script_sha256": sha256_file(Path(__file__).resolve()),
        "study_identity_sha256": sha256_file(study_dir / "study_identity.json"),
        "original_completion_marker_sha256": sha256_file(original_marker),
        "artifact_sha256": artifact_hashes,
        "artifact_hashes_sha256": hashlib.sha256(
            json.dumps(artifact_hashes, sort_keys=True).encode("utf-8")
        ).hexdigest(),
    }
    _json_dump(output_dir / "diagnostic_identity.json", identity)
    print(json.dumps(metrics, indent=2, sort_keys=True))
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-dir", type=Path, required=True)
    parser.add_argument("--facets-exe", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    run_diagnostic(
        study_dir=args.study_dir,
        facets_exe=args.facets_exe.resolve(),
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
