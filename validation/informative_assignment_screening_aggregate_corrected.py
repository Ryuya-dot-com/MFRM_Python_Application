#!/usr/bin/env python3
"""Result-blind metadata correction for the registered assignment screen."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Iterable

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from validation.estimand_bridge_pilot import CMLE_MODE  # noqa: E402
from validation.estimand_distribution_study import (  # noqa: E402
    FACETS_MODE,
    MML_FIXED_MODE,
    MML_FREE_MODE,
    PYTHON_JMLE_MODE,
    enrich_registered_metadata,
)
from validation.informative_assignment_screening import (  # noqa: E402
    SCIENTIFIC_MODES,
    screening_contrast_summary,
    validate_study_identity,
)
from validation.operating_characteristics_facets import sha256_file  # noqa: E402


SCHEMA_VERSION = "mfrm-informative-assignment-screening-aggregate-corrected-v1"
AMENDMENT_PATH = (
    REPO_ROOT
    / "validation"
    / "informative_assignment_screening_aggregation_amendment_20260811.json"
)
REGISTRATION_PATH = (
    REPO_ROOT
    / "validation"
    / "informative_assignment_screening_aggregation_execution_registration_20260811.json"
)


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _validate_registration() -> dict[str, Any]:
    registration = json.loads(REGISTRATION_PATH.read_text(encoding="utf-8"))
    expected = {
        "amendment_sha256": sha256_file(AMENDMENT_PATH),
        "aggregation_script_sha256": sha256_file(Path(__file__).resolve()),
    }
    for key, digest in expected.items():
        if str(registration.get(key, "")).lower() != digest.lower():
            raise ValueError(f"Aggregation registration mismatch: {key}")
    if not bool(registration.get("tests_passed_before_registration", False)):
        raise ValueError("Aggregation registration does not assert passing tests")
    return registration


def validate_original_aggregate(study_dir: Path) -> dict[str, Any]:
    aggregate = study_dir / "aggregate"
    identity_path = aggregate / "aggregate_identity.json"
    identity = json.loads(identity_path.read_text(encoding="utf-8"))
    for filename, expected in identity["artifact_sha256"].items():
        actual = sha256_file(aggregate / filename)
        if actual != expected:
            raise ValueError(f"Original aggregate artifact changed: {filename}")
    amendment = json.loads(AMENDMENT_PATH.read_text(encoding="utf-8"))
    expected_identity = amendment["trigger"]["original_aggregate_identity_sha256"]
    if sha256_file(identity_path) != expected_identity:
        raise ValueError("Original aggregate identity differs from registered amendment")
    return identity


def enrich_all_ledgers(
    runs: pd.DataFrame,
    recovery: pd.DataFrame,
    thresholds: pd.DataFrame,
    constraints: pd.DataFrame,
    manifest: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    frames = {
        "runs": runs,
        "recovery": recovery,
        "thresholds": thresholds,
        "constraints": constraints,
    }
    missing_before: dict[str, dict[str, int]] = {}
    enriched: dict[str, pd.DataFrame] = {}
    for name, frame in frames.items():
        missing_before[name] = {
            column: int(frame[column].isna().sum()) if column in frame else int(len(frame))
            for column in ("ConditionId", "Design", "PersonDistribution", "Replicate")
        }
        enriched[name] = enrich_registered_metadata(frame, manifest, label=name)
        if len(enriched[name]) and enriched[name]["RunId"].isna().any():
            raise ValueError(f"{name} contains an unknown RunId after enrichment")
    audit = {
        "rows": {name: int(len(frame)) for name, frame in enriched.items()},
        "missing_before": missing_before,
        "missing_after": {
            name: {
                column: int(frame[column].isna().sum())
                for column in ("ConditionId", "Design", "PersonDistribution", "Replicate")
            }
            for name, frame in enriched.items()
        },
        "all_registered_metadata_complete": all(
            frame[list(("ConditionId", "Design", "PersonDistribution", "Replicate"))]
            .notna()
            .all()
            .all()
            for frame in enriched.values()
            if len(frame)
        ),
    }
    return (
        enriched["runs"],
        enriched["recovery"],
        enriched["thresholds"],
        enriched["constraints"],
        audit,
    )


def build_corrected_outputs(
    runs: pd.DataFrame,
    recovery: pd.DataFrame,
    thresholds: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    eligible_recovery = recovery[
        recovery["EstimatorMode"].astype(str).isin(SCIENTIFIC_MODES)
        & recovery["IncludedInStudy"].fillna(False).astype(bool)
        & recovery["Facet"].astype(str).isin({"Rater", "Task", "Criterion"})
    ].copy()
    eligible_recovery["ErrorAligned"] = pd.to_numeric(
        eligible_recovery["ErrorAligned"], errors="coerce"
    )
    facet_loss = eligible_recovery.groupby(
        ["RunId", "Replicate", "Design", "EstimatorMode", "Facet"], as_index=False
    )["ErrorAligned"].agg(
        RMSE=lambda values: float(np.sqrt(np.mean(np.square(values)))),
        MAE=lambda values: float(np.mean(np.abs(values))),
    )
    facet_long = facet_loss.melt(
        id_vars=["RunId", "Replicate", "Design", "EstimatorMode", "Facet"],
        value_vars=["RMSE", "MAE"],
        var_name="Metric",
        value_name="Loss",
    )
    facet_long["RecoveryDomain"] = "Facet:" + facet_long["Facet"].astype(str)

    eligible_thresholds = thresholds[
        thresholds["EstimatorMode"].astype(str).isin(SCIENTIFIC_MODES)
        & thresholds["IncludedInStudy"].fillna(False).astype(bool)
    ].copy()
    eligible_thresholds["TruthError"] = pd.to_numeric(
        eligible_thresholds["TruthError"], errors="coerce"
    )
    threshold_loss = eligible_thresholds.groupby(
        ["RunId", "Replicate", "Design", "EstimatorMode"], as_index=False
    )["TruthError"].agg(
        ThresholdRows="count",
        RMSE=lambda values: float(np.sqrt(np.mean(np.square(values)))),
        MAE=lambda values: float(np.mean(np.abs(values))),
    )
    threshold_loss = threshold_loss[threshold_loss["ThresholdRows"].eq(6)]
    threshold_long = threshold_loss.melt(
        id_vars=["RunId", "Replicate", "Design", "EstimatorMode"],
        value_vars=["RMSE", "MAE"],
        var_name="Metric",
        value_name="Loss",
    )
    threshold_long["RecoveryDomain"] = "Threshold"
    run_loss = pd.concat(
        [
            facet_long.drop(columns="Facet"),
            threshold_long.drop(columns="ThresholdRows", errors="ignore"),
        ],
        ignore_index=True,
        sort=False,
    )
    planned = run_loss[run_loss["Design"].eq("planned_connected")].drop(
        columns=["RunId", "Design"]
    )
    aligned = run_loss[
        run_loss["Design"].eq("ability_severity_aligned_connected")
    ].drop(columns=["RunId", "Design"])
    contrasts = aligned.merge(
        planned,
        on=["Replicate", "EstimatorMode", "RecoveryDomain", "Metric"],
        suffixes=("Aligned", "Planned"),
        validate="one_to_one",
    )
    contrasts["ContrastAlignedMinusPlanned"] = (
        contrasts["LossAligned"] - contrasts["LossPlanned"]
    )
    screening = screening_contrast_summary(contrasts)

    rater_errors = eligible_recovery[eligible_recovery["Facet"].eq("Rater")]
    rater_planned = rater_errors[rater_errors["Design"].eq("planned_connected")][
        ["Replicate", "EstimatorMode", "Level", "ErrorAligned"]
    ]
    rater_aligned = rater_errors[
        rater_errors["Design"].eq("ability_severity_aligned_connected")
    ][["Replicate", "EstimatorMode", "Level", "ErrorAligned"]]
    rater_contrasts = rater_aligned.merge(
        rater_planned,
        on=["Replicate", "EstimatorMode", "Level"],
        suffixes=("Aligned", "Planned"),
        validate="one_to_one",
    )
    rater_contrasts["ErrorContrastAlignedMinusPlanned"] = (
        rater_contrasts["ErrorAlignedAligned"] - rater_contrasts["ErrorAlignedPlanned"]
    )
    rater_summary = rater_contrasts.groupby(
        ["EstimatorMode", "Level"], as_index=False
    )["ErrorContrastAlignedMinusPlanned"].agg(N="count", Mean="mean", SD="std")
    context = run_loss.groupby(
        ["EstimatorMode", "Design", "RecoveryDomain", "Metric"], as_index=False
    )["Loss"].agg(N="count", Mean="mean", SD="std")
    free_sd = runs[runs["EstimatorMode"].astype(str).eq(MML_FREE_MODE)].copy()
    free_sd_paired = free_sd.pivot(
        index="Replicate", columns="Design", values="EstimatedPopulationSD"
    ).reset_index()
    free_sd_paired["ContrastAlignedMinusPlanned"] = (
        free_sd_paired["ability_severity_aligned_connected"]
        - free_sd_paired["planned_connected"]
    )
    return {
        "run_loss.csv": run_loss,
        "paired_design_contrasts.csv": contrasts,
        "screening_contrast_summary.csv": screening,
        "rater_level_contrasts.csv": rater_contrasts,
        "rater_level_summary.csv": rater_summary,
        "design_context_summary.csv": context,
        "free_sd_paired.csv": free_sd_paired,
    }


def aggregate_corrected(study_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    validate_study_identity(study_dir)
    registration = _validate_registration()
    original_identity = validate_original_aggregate(study_dir)
    original_dir = study_dir / "aggregate"
    manifest = pd.read_csv(study_dir / "retained_input" / "manifest.csv")
    runs = pd.read_csv(original_dir / "run_ledger.csv")
    recovery = pd.read_csv(original_dir / "recovery.csv")
    thresholds = pd.read_csv(original_dir / "thresholds.csv")
    constraints = pd.read_csv(original_dir / "constraints.csv")
    outcomes = pd.read_csv(original_dir / "attempt_outcomes.csv")
    runs, recovery, thresholds, constraints, join_audit = enrich_all_ledgers(
        runs, recovery, thresholds, constraints, manifest
    )
    if not join_audit["all_registered_metadata_complete"]:
        raise ValueError("Registered metadata enrichment remained incomplete")
    calculated = build_corrected_outputs(runs, recovery, thresholds)
    screening = calculated["screening_contrast_summary.csv"]
    original_metrics = json.loads(
        (original_dir / "screening_metrics.json").read_text(encoding="utf-8")
    )
    original_operational = {
        key: value
        for key, value in original_metrics["gates"].items()
        if key != "screening_pairs_20_each"
    }
    gates = {
        **original_operational,
        "metadata_join_complete": join_audit["all_registered_metadata_complete"],
        "screening_cells_32": len(screening) == 32,
        "screening_pairs_20_each": len(screening) == 32
        and screening["FinitePairs"].eq(20).all(),
        "confirmatory_claims_withheld": not screening["ConfirmatoryClaimAllowed"].any(),
    }
    marker_digest = original_metrics["completion_marker_set_sha256"]
    metrics = {
        "schema_version": SCHEMA_VERSION,
        "phase": "screening",
        "datasets": 60,
        "attempts": 240,
        "gates": {key: bool(value) for key, value in gates.items()},
        "qualification_pass": all(bool(value) for value in gates.values()),
        "facets_calibration": original_metrics["facets_calibration"],
        "completion_marker_set_sha256": marker_digest,
        "metadata_correction_only": True,
        "original_aggregate_qualification_pass": False,
        "original_aggregate_preserved": True,
        "confirmatory_claims_allowed": False,
        "estimator_ranking_constructed": False,
        "cross_basis_likelihood_comparison": "prohibited",
        "claim_limit": original_metrics["claim_limit"],
    }
    output_dir = study_dir / "aggregate_corrected"
    if output_dir.exists():
        raise FileExistsError(f"Corrected aggregate already exists: {output_dir}")
    output_dir.mkdir()
    source_frames = {
        "attempt_outcomes.csv": outcomes,
        "run_ledger.csv": runs,
        "recovery.csv": recovery,
        "thresholds.csv": thresholds,
        "constraints.csv": constraints,
        **calculated,
    }
    for filename, frame in source_frames.items():
        frame.to_csv(output_dir / filename, index=False, lineterminator="\n")
    _json_dump(output_dir / "metadata_join_audit.json", join_audit)
    _json_dump(output_dir / "screening_metrics_corrected.json", metrics)
    artifacts = tuple(source_frames) + (
        "metadata_join_audit.json",
        "screening_metrics_corrected.json",
    )
    identity = {
        "schema_version": f"{SCHEMA_VERSION}-identity",
        "study_identity_sha256": sha256_file(study_dir / "study_identity.json"),
        "original_aggregate_identity_sha256": sha256_file(
            original_dir / "aggregate_identity.json"
        ),
        "original_aggregate_artifact_contract": original_identity["artifact_sha256"],
        "amendment_sha256": sha256_file(AMENDMENT_PATH),
        "aggregation_registration_sha256": sha256_file(REGISTRATION_PATH),
        "aggregation_script_sha256": sha256_file(Path(__file__).resolve()),
        "completion_marker_set_sha256": marker_digest,
        "artifact_sha256": {
            filename: sha256_file(output_dir / filename) for filename in artifacts
        },
        "registration": registration,
    }
    _json_dump(output_dir / "aggregate_identity.json", identity)
    return metrics


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    print(json.dumps(aggregate_corrected(args.study_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

