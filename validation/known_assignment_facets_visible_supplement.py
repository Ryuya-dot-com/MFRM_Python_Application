"""Run separate visible-mode FACETS calibrations on the 12 frozen inputs."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pandas as pd

from validation.facets_resilient_pair import (
    dependency_manifest,
    dependency_manifest_digest,
    fit_resilient_pair,
)
from validation.operating_characteristics_facets import sha256_file
from validation.operating_characteristics_facets import (
    FACETS_WINDOWS_SAFE_PATH_CHARS,
    _facets_path_budget_violations,
    slugify,
)


PLAN_PATH = ROOT / "validation" / "known_assignment_facets_visible_supplement_plan_20260811.json"
REGISTRATION_PATH = ROOT / "validation" / "known_assignment_facets_visible_supplement_registration_v3_20260811.json"
PATH_AMENDMENT_PATH = ROOT / "validation" / "known_assignment_facets_visible_supplement_path_amendment_20260811.json"
SOURCE_STUDY = ROOT / "validation" / "known_assignment_multivector_preflight4_20260811"
OUTPUT_DIR = SOURCE_STUDY / "v"
ATTEMPT_TYPE = "RESILIENT_FACETS_PYTHON_JMLE_PCM"


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _artifact_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): sha256_file(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _validate_prerequisites() -> tuple[dict[str, Any], dict[str, Any]]:
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    registration = json.loads(REGISTRATION_PATH.read_text(encoding="utf-8"))
    expected = {
        "plan_sha256": sha256_file(PLAN_PATH),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "path_amendment_sha256": sha256_file(PATH_AMENDMENT_PATH),
    }
    for key, value in expected.items():
        if registration.get(key) != value:
            raise ValueError(f"Supplement registration mismatch: {key}")
    dependencies = dependency_manifest()
    if registration["dependency_manifest_sha256"] != dependency_manifest_digest(dependencies):
        raise ValueError("Supplement dependency manifest changed")
    for key, relative, digest_key in (
        ("source identity", plan["source_study_identity"], "source_study_identity_sha256"),
        ("preendpoint audit", plan["preendpoint_audit"], "preendpoint_audit_sha256"),
        ("helper qualification", plan["helper_qualification"], "helper_qualification_sha256"),
    ):
        if sha256_file(ROOT / relative).lower() != plan[digest_key].lower():
            raise ValueError(f"{key} hash mismatch")
    helper = json.loads((ROOT / plan["helper_qualification"]).read_text(encoding="utf-8"))
    if not bool(helper.get("pass")) or bool(helper.get("scientific_endpoint_read")):
        raise ValueError("Shared helper qualification is not a clean PASS")
    identity = json.loads((ROOT / plan["source_study_identity"]).read_text(encoding="utf-8"))
    for filename, expected_hash in identity["retained_input_sha256"].items():
        if sha256_file(SOURCE_STUDY / "retained_input" / filename) != expected_hash:
            raise ValueError(f"Frozen retained input changed: {filename}")
    facets_exe = Path(plan["facets_executable"])
    if sha256_file(facets_exe) != plan["facets_executable_sha256"]:
        raise ValueError("FACETS executable changed")
    return plan, registration


def run(*, resume: bool) -> dict[str, int]:
    plan, registration = _validate_prerequisites()
    if OUTPUT_DIR.exists() and not resume:
        raise FileExistsError(f"Refusing to overwrite supplemental calibration: {OUTPUT_DIR}")
    OUTPUT_DIR.mkdir(exist_ok=True)
    (OUTPUT_DIR / "markers").mkdir(exist_ok=True)
    (OUTPUT_DIR / "locks").mkdir(exist_ok=True)
    input_dir = SOURCE_STUDY / "retained_input"
    attempts = pd.read_csv(input_dir / "attempt_manifest.csv")
    attempts = attempts.loc[attempts["AttemptType"].astype(str).eq(ATTEMPT_TYPE)].copy()
    if len(attempts) != int(plan["execution"]["datasets"]):
        raise ValueError(f"Expected 12 JMLE attempts, found {len(attempts)}")
    manifest = pd.read_csv(input_dir / "manifest.csv").set_index("RunId", drop=False)
    ratings_all = pd.read_csv(input_dir / "generated_ratings.csv")
    truth_all = pd.read_csv(input_dir / "generated_facet_truth.csv")
    threshold_all = pd.read_csv(input_dir / "generated_pcm_threshold_truth.csv")
    os.environ["MFRM_FACETS_BATCH_MODE"] = "NO"
    summary = {"assigned": len(attempts), "completed_now": 0, "skipped": 0}
    for sequence, (_, attempt) in enumerate(attempts.iterrows(), start=1):
        ordinal = int(attempt["AttemptOrdinal"])
        marker_path = OUTPUT_DIR / "markers" / f"{ordinal:05d}.json"
        if marker_path.is_file():
            marker = json.loads(marker_path.read_text(encoding="utf-8"))
            root = OUTPUT_DIR / marker["artifact_root"]
            for relative, digest in marker["artifact_sha256"].items():
                if sha256_file(root / relative) != digest:
                    raise ValueError(f"Supplement marker artifact changed: {relative}")
            summary["skipped"] += 1
            print(f"[{sequence}/12] skip {attempt['RunId']}", flush=True)
            continue
        artifact_root = OUTPUT_DIR / f"{ordinal:04x}"
        if artifact_root.exists():
            raise FileExistsError(f"Unmarked supplemental artifact requires adjudication: {artifact_root}")
        run_id = str(attempt["RunId"])
        prospective_run_dir = (
            artifact_root.resolve()
            / "facets_attempts"
            / "pair_try_01"
            / "facets_runs"
            / f"{slugify(run_id)}__pcm"
        )
        prospective_paths = (
            prospective_run_dir / "analysis.txt",
            prospective_run_dir / "report_u6.txt",
            prospective_run_dir / "scores.4.txt",
            prospective_run_dir / "report_u2.txt",
            prospective_run_dir / "scores_u2.4.txt",
        )
        violations = _facets_path_budget_violations(prospective_paths)
        if violations:
            detail = "; ".join(f"{length}: {path}" for path, length in violations)
            raise OSError(
                "Supplemental FACETS path exceeds the registered "
                f"{FACETS_WINDOWS_SAFE_PATH_CHARS}-character budget: {detail}"
            )
        run_ratings_with_id = ratings_all.loc[ratings_all["RunId"].astype(str).eq(run_id)]
        run_truth = truth_all.loc[truth_all["RunId"].astype(str).eq(run_id)]
        run_thresholds = threshold_all.loc[threshold_all["RunId"].astype(str).eq(run_id)]
        input_hash = str(attempt["RunInputSHA256"])
        outcome = fit_resilient_pair(
            manifest.loc[run_id],
            run_ratings_with_id.drop(columns="RunId"),
            run_truth,
            run_thresholds,
            facets_exe=Path(plan["facets_executable"]),
            output_dir=artifact_root,
            lock_path=OUTPUT_DIR / "locks" / "facets_global.lock",
            timeout_seconds=float(plan["execution"]["timeout_seconds_per_facets_call"]),
            maximum_retries=int(plan["execution"]["maximum_retries_per_dataset"]),
        )
        marker = {
            "schema_version": "known_assignment_facets_visible_supplement_marker_v1",
            "run_id": run_id,
            "attempt_id": str(attempt["AttemptId"]),
            "attempt_fingerprint": str(attempt["AttemptFingerprint"]),
            "run_input_sha256": input_hash,
            "dependency_manifest_sha256": registration["dependency_manifest_sha256"],
            "runner_sha256": registration["runner_sha256"],
            "artifact_root": artifact_root.relative_to(OUTPUT_DIR).as_posix(),
            "artifact_sha256": _artifact_hashes(artifact_root),
            "calibration_ready": bool(outcome["calibration_ready"]),
            "direct_agreement_pass": bool(outcome["direct_agreement_pass"]),
            "facets_retry_count": int(outcome["facets_retry_count"]),
            "scientific_endpoint_read": False,
        }
        _json_dump(marker_path, marker)
        summary["completed_now"] += 1
        print(
            f"[{sequence}/12] {run_id}: calibration={marker['calibration_ready']}",
            flush=True,
        )
    return summary


def assess() -> dict[str, Any]:
    plan, registration = _validate_prerequisites()
    marker_paths = sorted((OUTPUT_DIR / "markers").glob("*.json"))
    if len(marker_paths) != int(plan["execution"]["datasets"]):
        raise ValueError(f"Expected 12 supplemental markers, found {len(marker_paths)}")
    rows = []
    marker_hashes = {}
    for marker_path in marker_paths:
        marker = json.loads(marker_path.read_text(encoding="utf-8"))
        root = OUTPUT_DIR / marker["artifact_root"]
        for relative, digest in marker["artifact_sha256"].items():
            if sha256_file(root / relative) != digest:
                raise ValueError(f"Supplement artifact mismatch: {relative}")
        metrics = json.loads((root / "resilient_pair_metrics.json").read_text(encoding="utf-8"))
        run_ledger = pd.read_csv(root / "combined_run_ledger.csv")
        facets_row = run_ledger.loc[run_ledger["EstimatorMode"].astype(str).eq("FACETS_4_5_JMLE")]
        if len(facets_row) != 1:
            raise ValueError(f"Expected one FACETS row for {marker['run_id']}")
        row = {
            **{key: marker[key] for key in (
                "run_id", "attempt_id", "run_input_sha256", "calibration_ready",
                "direct_agreement_pass", "facets_retry_count",
            )},
            "python_replay_main_max_abs_difference": metrics["python_replay_main_max_abs_difference"],
            "python_replay_threshold_max_abs_difference": metrics["python_replay_threshold_max_abs_difference"],
            "facets_reported_version": metrics["facets_reported_version"],
            "main_weighted_mae": facets_row.iloc[0].get("MainWeightedMAE"),
            "main_max_abs_difference": facets_row.iloc[0].get("MainMaxAbsDifference"),
            "threshold_weighted_mae": facets_row.iloc[0].get("ThresholdWeightedMAE"),
            "threshold_max_abs_difference": facets_row.iloc[0].get("ThresholdMaxAbsDifference"),
            "minimum_within_facet_spearman": facets_row.iloc[0].get("MinimumWithinFacetSpearman"),
        }
        rows.append(row)
        marker_hashes[marker_path.name] = sha256_file(marker_path)
    ledger = pd.DataFrame(rows)
    gates = {
        "completed_datasets_12": len(ledger) == 12,
        "calibration_ready_12": bool(ledger["calibration_ready"].all()),
        "direct_agreement_pass_12": bool(ledger["direct_agreement_pass"].all()),
        "zero_retries": int(ledger["facets_retry_count"].sum()) == 0,
        "python_replay_within_tolerance": bool(
            ledger[[
                "python_replay_main_max_abs_difference",
                "python_replay_threshold_max_abs_difference",
            ]].max().max() <= float(plan["operational_gate"]["python_replay_tolerance"])
        ),
        "all_input_hashes_match": True,
    }
    passed = bool(all(gates.values()))
    ledger.to_csv(OUTPUT_DIR / "calibration_ledger.csv", index=False, lineterminator="\n")
    assessment = {
        "schema_version": "known_assignment_facets_visible_supplement_assessment_v1",
        "pass": passed,
        "gates": gates,
        "datasets": len(ledger),
        "maximum_main_weighted_mae": float(pd.to_numeric(ledger["main_weighted_mae"]).max()),
        "maximum_main_absolute_difference": float(pd.to_numeric(ledger["main_max_abs_difference"]).max()),
        "maximum_threshold_weighted_mae": float(pd.to_numeric(ledger["threshold_weighted_mae"]).max()),
        "maximum_threshold_absolute_difference": float(pd.to_numeric(ledger["threshold_max_abs_difference"]).max()),
        "minimum_within_facet_spearman": float(pd.to_numeric(ledger["minimum_within_facet_spearman"]).min()),
        "marker_sha256": marker_hashes,
        "registration_sha256": sha256_file(REGISTRATION_PATH),
        "aggregate_join_authorized": passed,
        "scientific_endpoint_read": False,
        "aggregate_created": False,
        "claim_boundary": plan["claim_boundary"],
    }
    _json_dump(OUTPUT_DIR / "assessment.json", assessment)
    print(json.dumps(assessment, ensure_ascii=False))
    if not passed:
        raise SystemExit(2)
    return assessment


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("run", "assess"))
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.command == "run":
        print(json.dumps(run(resume=args.resume), ensure_ascii=False))
    else:
        assess()


if __name__ == "__main__":
    main()
