#!/usr/bin/env python3
"""Produce preregistered paired summaries for the 20-replicate shape screen."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Iterable

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from validation.estimand_distribution_study import (  # noqa: E402
    CMLE_MODE,
    DISTRIBUTIONS,
    FACETS_MODE,
    MML_FREE_MODE,
    MML_FIXED_MODE,
    PERSON_SD,
    PYTHON_JMLE_MODE,
)
from validation.operating_characteristics_facets import sha256_file  # noqa: E402


SCHEMA_VERSION = "mfrm-estimand-distribution-screening-analysis-v1"
PLAN_PATH = REPO_ROOT / "validation" / "estimand_distribution_screening_plan_20260811.json"
ANALYSIS_AMENDMENT_PATH = (
    REPO_ROOT / "validation" / "estimand_distribution_screening_analysis_amendment2_20260811.json"
)
ESTIMATOR_MODES = (
    FACETS_MODE,
    PYTHON_JMLE_MODE,
    MML_FIXED_MODE,
    MML_FREE_MODE,
    CMLE_MODE,
)
MAIN_FACETS = ("Rater", "Task", "Criterion")
NONNORMAL = tuple(value for value in DISTRIBUTIONS if value != "normal")


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _rmse(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    return float(np.sqrt(np.mean(np.square(numeric)))) if len(numeric) else np.nan


def _as_bool(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype(str).str.lower().eq("true")


def summarize_contrasts(
    pairs: pd.DataFrame,
    *,
    group_columns: list[str],
    metric_columns: list[str],
) -> pd.DataFrame:
    """Summarize registered paired differences without hypothesis tests."""

    rows: list[dict[str, Any]] = []
    for keys, group in pairs.groupby(group_columns, dropna=False, sort=False):
        key_values = keys if isinstance(keys, tuple) else (keys,)
        metadata = dict(zip(group_columns, key_values))
        for metric in metric_columns:
            difference_column = f"{metric}Contrast"
            values = pd.to_numeric(group[difference_column], errors="coerce").dropna()
            n = int(len(values))
            mean = float(values.mean()) if n else np.nan
            sd = float(values.std(ddof=1)) if n > 1 else np.nan
            se = float(sd / np.sqrt(n)) if n > 1 else np.nan
            rows.append({
                **metadata,
                "Metric": metric,
                "PairedReplicates": n,
                "MeanContrast": mean,
                "MonteCarloSD": sd,
                "MonteCarloSE": se,
                "ScreeningNormalApproxLower95": mean - 1.96 * se if np.isfinite(se) else np.nan,
                "ScreeningNormalApproxUpper95": mean + 1.96 * se if np.isfinite(se) else np.nan,
                "IntervalExcludesZero": bool(
                    np.isfinite(se) and (mean - 1.96 * se > 0 or mean + 1.96 * se < 0)
                ),
                "EvidenceTier": "20-replicate screening; not confirmatory",
            })
    return pd.DataFrame(rows)


def build_shape_pairs(
    run_metrics: pd.DataFrame,
    *,
    identity_columns: list[str],
    metric_columns: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    baseline = run_metrics[run_metrics["PersonDistribution"].eq("normal")].drop(
        columns=["RunId", "PersonDistribution"]
    )
    nonnormal = run_metrics[run_metrics["PersonDistribution"].isin(NONNORMAL)].copy()
    pairs = nonnormal.merge(
        baseline,
        on=identity_columns,
        how="inner",
        suffixes=("", "Normal"),
        validate="many_to_one",
    )
    for metric in metric_columns:
        pairs[f"{metric}Contrast"] = pairs[metric] - pairs[f"{metric}Normal"]
    group_columns = [
        column for column in ("EstimatorMode", "PersonDistribution", "Design", "Facet")
        if column in pairs
    ]
    summary = summarize_contrasts(
        pairs, group_columns=group_columns, metric_columns=metric_columns
    )
    summary["ContrastDefinition"] = "non-normal minus normal"
    return pairs, summary


def build_design_pairs(
    run_metrics: pd.DataFrame,
    *,
    identity_columns: list[str],
    metric_columns: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    complete = run_metrics[run_metrics["Design"].eq("complete")].drop(
        columns=["RunId", "Design"]
    )
    planned = run_metrics[run_metrics["Design"].eq("planned_connected")].copy()
    pairs = planned.merge(
        complete,
        on=identity_columns,
        how="inner",
        suffixes=("", "Complete"),
        validate="one_to_one",
    )
    for metric in metric_columns:
        pairs[f"{metric}Contrast"] = pairs[metric] - pairs[f"{metric}Complete"]
    group_columns = [
        column for column in ("EstimatorMode", "PersonDistribution", "Facet")
        if column in pairs
    ]
    summary = summarize_contrasts(
        pairs, group_columns=group_columns, metric_columns=metric_columns
    )
    summary["ContrastDefinition"] = "planned-connected minus complete"
    return pairs, summary


def run_analysis(study_dir: Path, output_dir: Path) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    aggregate_dir = study_dir / "aggregate"
    metrics = json.loads((aggregate_dir / "study_metrics.json").read_text(encoding="utf-8"))
    if not metrics.get("qualification_pass", False) or metrics.get("datasets") != 160:
        raise ValueError("Qualified 160-dataset screening aggregate is required")
    runs = pd.read_csv(aggregate_dir / "study_run_ledger.csv")
    recovery = pd.read_csv(aggregate_dir / "study_recovery.csv")
    thresholds = pd.read_csv(aggregate_dir / "study_thresholds.csv")
    attempts = pd.read_csv(aggregate_dir / "attempt_ledger.csv")
    if set(runs["EstimatorMode"].dropna().astype(str)) != set(ESTIMATOR_MODES):
        raise ValueError("Estimator mode set differs from the registered analysis")
    output_dir.mkdir(parents=True)

    eligible_recovery = recovery[
        _as_bool(recovery["IncludedInStudy"])
        & recovery["Facet"].astype(str).isin(MAIN_FACETS)
    ].copy()
    main_run = (
        eligible_recovery.groupby(
            ["RunId", "Replicate", "EstimatorMode", "PersonDistribution", "Design", "Facet"],
            dropna=False,
        )["ErrorAligned"]
        .agg(MeanError="mean", MAE=lambda x: float(x.abs().mean()), RMSE=_rmse)
        .reset_index()
    )
    main_shape_pairs, main_shape_summary = build_shape_pairs(
        main_run,
        identity_columns=["Replicate", "EstimatorMode", "Design", "Facet"],
        metric_columns=["MeanError", "MAE", "RMSE"],
    )
    main_design_pairs, main_design_summary = build_design_pairs(
        main_run,
        identity_columns=["Replicate", "EstimatorMode", "PersonDistribution", "Facet"],
        metric_columns=["MeanError", "MAE", "RMSE"],
    )

    main_element_condition = (
        eligible_recovery.groupby(
            ["EstimatorMode", "PersonDistribution", "Design", "Facet", "Level"],
            dropna=False,
        )["ErrorAligned"]
        .agg(
            Replicates="size",
            MeanError="mean",
            ErrorSD="std",
            MinimumError="min",
            MaximumError="max",
        )
        .reset_index()
    )
    main_element_condition["MonteCarloSE"] = (
        main_element_condition["ErrorSD"] / np.sqrt(main_element_condition["Replicates"])
    )
    main_element_condition["ScreeningNormalApproxLower95"] = (
        main_element_condition["MeanError"] - 1.96 * main_element_condition["MonteCarloSE"]
    )
    main_element_condition["ScreeningNormalApproxUpper95"] = (
        main_element_condition["MeanError"] + 1.96 * main_element_condition["MonteCarloSE"]
    )
    main_element_condition["IntervalExcludesZero"] = (
        main_element_condition["ScreeningNormalApproxLower95"].gt(0)
        | main_element_condition["ScreeningNormalApproxUpper95"].lt(0)
    )
    facet_bias_summary = (
        main_element_condition.groupby(
            ["EstimatorMode", "PersonDistribution", "Design", "Facet"], dropna=False
        )["MeanError"]
        .agg(
            Levels="size",
            RMSLevelBias=_rmse,
            MaxAbsLevelBias=lambda x: float(pd.to_numeric(x).abs().max()),
        )
        .reset_index()
    )
    normal_elements = eligible_recovery[
        eligible_recovery["PersonDistribution"].eq("normal")
    ].drop(columns=["RunId", "PersonDistribution"])
    nonnormal_elements = eligible_recovery[
        eligible_recovery["PersonDistribution"].isin(NONNORMAL)
    ].copy()
    main_element_shape_pairs = nonnormal_elements.merge(
        normal_elements,
        on=["Replicate", "EstimatorMode", "Design", "Facet", "Level"],
        suffixes=("", "Normal"),
        validate="many_to_one",
    )
    main_element_shape_pairs["ErrorAlignedContrast"] = (
        main_element_shape_pairs["ErrorAligned"]
        - main_element_shape_pairs["ErrorAlignedNormal"]
    )
    main_element_shape_summary = summarize_contrasts(
        main_element_shape_pairs,
        group_columns=["EstimatorMode", "PersonDistribution", "Design", "Facet", "Level"],
        metric_columns=["ErrorAligned"],
    )
    main_element_shape_summary["ContrastDefinition"] = "non-normal minus normal"
    complete_elements = eligible_recovery[
        eligible_recovery["Design"].eq("complete")
    ].drop(columns=["RunId", "Design"])
    planned_elements = eligible_recovery[
        eligible_recovery["Design"].eq("planned_connected")
    ].copy()
    main_element_design_pairs = planned_elements.merge(
        complete_elements,
        on=["Replicate", "EstimatorMode", "PersonDistribution", "Facet", "Level"],
        suffixes=("", "Complete"),
        validate="one_to_one",
    )
    main_element_design_pairs["ErrorAlignedContrast"] = (
        main_element_design_pairs["ErrorAligned"]
        - main_element_design_pairs["ErrorAlignedComplete"]
    )
    main_element_design_summary = summarize_contrasts(
        main_element_design_pairs,
        group_columns=["EstimatorMode", "PersonDistribution", "Facet", "Level"],
        metric_columns=["ErrorAligned"],
    )
    main_element_design_summary["ContrastDefinition"] = "planned-connected minus complete"

    eligible_thresholds = thresholds[_as_bool(thresholds["IncludedInStudy"])].copy()
    threshold_run = (
        eligible_thresholds.groupby(
            ["RunId", "Replicate", "EstimatorMode", "PersonDistribution", "Design"],
            dropna=False,
        )["TruthError"]
        .agg(MeanError="mean", MAE=lambda x: float(x.abs().mean()), RMSE=_rmse, Thresholds="size")
        .reset_index()
    )
    threshold_shape_pairs, threshold_shape_summary = build_shape_pairs(
        threshold_run,
        identity_columns=["Replicate", "EstimatorMode", "Design"],
        metric_columns=["MeanError", "MAE", "RMSE"],
    )
    threshold_design_pairs, threshold_design_summary = build_design_pairs(
        threshold_run,
        identity_columns=["Replicate", "EstimatorMode", "PersonDistribution"],
        metric_columns=["MeanError", "MAE", "RMSE"],
    )

    threshold_element_condition = (
        eligible_thresholds.groupby(
            [
                "EstimatorMode", "PersonDistribution", "Design",
                "StepFacetLevel", "Category",
            ],
            dropna=False,
        )["TruthError"]
        .agg(
            Replicates="size",
            MeanError="mean",
            ErrorSD="std",
            MinimumError="min",
            MaximumError="max",
        )
        .reset_index()
    )
    threshold_element_condition["MonteCarloSE"] = (
        threshold_element_condition["ErrorSD"]
        / np.sqrt(threshold_element_condition["Replicates"])
    )
    threshold_element_condition["ScreeningNormalApproxLower95"] = (
        threshold_element_condition["MeanError"]
        - 1.96 * threshold_element_condition["MonteCarloSE"]
    )
    threshold_element_condition["ScreeningNormalApproxUpper95"] = (
        threshold_element_condition["MeanError"]
        + 1.96 * threshold_element_condition["MonteCarloSE"]
    )
    threshold_element_condition["IntervalExcludesZero"] = (
        threshold_element_condition["ScreeningNormalApproxLower95"].gt(0)
        | threshold_element_condition["ScreeningNormalApproxUpper95"].lt(0)
    )
    threshold_bias_summary = (
        threshold_element_condition.groupby(
            ["EstimatorMode", "PersonDistribution", "Design"], dropna=False
        )["MeanError"]
        .agg(
            Thresholds="size",
            RMSThresholdBias=_rmse,
            MaxAbsThresholdBias=lambda x: float(pd.to_numeric(x).abs().max()),
        )
        .reset_index()
    )
    normal_thresholds = eligible_thresholds[
        eligible_thresholds["PersonDistribution"].eq("normal")
    ].drop(columns=["RunId", "PersonDistribution"])
    nonnormal_thresholds = eligible_thresholds[
        eligible_thresholds["PersonDistribution"].isin(NONNORMAL)
    ].copy()
    threshold_element_shape_pairs = nonnormal_thresholds.merge(
        normal_thresholds,
        on=[
            "Replicate", "EstimatorMode", "Design", "StepFacetLevel", "Category"
        ],
        suffixes=("", "Normal"),
        validate="many_to_one",
    )
    threshold_element_shape_pairs["TruthErrorContrast"] = (
        threshold_element_shape_pairs["TruthError"]
        - threshold_element_shape_pairs["TruthErrorNormal"]
    )
    threshold_element_shape_summary = summarize_contrasts(
        threshold_element_shape_pairs,
        group_columns=[
            "EstimatorMode", "PersonDistribution", "Design", "StepFacetLevel", "Category"
        ],
        metric_columns=["TruthError"],
    )
    threshold_element_shape_summary["ContrastDefinition"] = "non-normal minus normal"
    complete_thresholds = eligible_thresholds[
        eligible_thresholds["Design"].eq("complete")
    ].drop(columns=["RunId", "Design"])
    planned_thresholds = eligible_thresholds[
        eligible_thresholds["Design"].eq("planned_connected")
    ].copy()
    threshold_element_design_pairs = planned_thresholds.merge(
        complete_thresholds,
        on=[
            "Replicate", "EstimatorMode", "PersonDistribution", "StepFacetLevel", "Category"
        ],
        suffixes=("", "Complete"),
        validate="one_to_one",
    )
    threshold_element_design_pairs["TruthErrorContrast"] = (
        threshold_element_design_pairs["TruthError"]
        - threshold_element_design_pairs["TruthErrorComplete"]
    )
    threshold_element_design_summary = summarize_contrasts(
        threshold_element_design_pairs,
        group_columns=[
            "EstimatorMode", "PersonDistribution", "StepFacetLevel", "Category"
        ],
        metric_columns=["TruthError"],
    )
    threshold_element_design_summary["ContrastDefinition"] = "planned-connected minus complete"

    free_sd = runs[runs["EstimatorMode"].eq(MML_FREE_MODE)].copy()
    free_sd["EstimatedPopulationSD"] = pd.to_numeric(
        free_sd["EstimatedPopulationSD"], errors="coerce"
    )
    free_sd["SDError"] = free_sd["EstimatedPopulationSD"] - PERSON_SD
    free_sd_condition = (
        free_sd.groupby(["PersonDistribution", "Design"], dropna=False)
        .agg(
            Replicates=("RunId", "size"),
            MeanEstimatedSD=("EstimatedPopulationSD", "mean"),
            SDEstimateSD=("EstimatedPopulationSD", "std"),
            MeanSDError=("SDError", "mean"),
            RMSESDError=("SDError", _rmse),
            MinimumEstimatedSD=("EstimatedPopulationSD", "min"),
            MaximumEstimatedSD=("EstimatedPopulationSD", "max"),
        )
        .reset_index()
    )
    free_sd_condition["MonteCarloSEMeanSD"] = (
        free_sd_condition["SDEstimateSD"] / np.sqrt(free_sd_condition["Replicates"])
    )
    free_sd_shape_pairs, free_sd_shape_summary = build_shape_pairs(
        free_sd[["RunId", "Replicate", "PersonDistribution", "Design", "EstimatedPopulationSD"]],
        identity_columns=["Replicate", "Design"],
        metric_columns=["EstimatedPopulationSD"],
    )
    free_sd_design_pairs, free_sd_design_summary = build_design_pairs(
        free_sd[["RunId", "Replicate", "PersonDistribution", "Design", "EstimatedPopulationSD"]],
        identity_columns=["Replicate", "PersonDistribution"],
        metric_columns=["EstimatedPopulationSD"],
    )

    cmle = runs[runs["EstimatorMode"].eq(CMLE_MODE)].copy()
    cmle["PersonsExtreme"] = pd.to_numeric(cmle["PersonsExtreme"], errors="coerce")
    cmle["PersonsInformative"] = pd.to_numeric(cmle["PersonsInformative"], errors="coerce")
    cmle["PersonDenominator"] = cmle["PersonsExtreme"] + cmle["PersonsInformative"]
    cmle_extremes = (
        cmle.groupby(["PersonDistribution", "Design"], dropna=False)
        .agg(
            AttemptedRuns=("RunId", "size"),
            EligibleRuns=("IncludedInStudy", lambda x: int(_as_bool(x).sum())),
            RunsWithExtremePersons=("PersonsExtreme", lambda x: int(pd.to_numeric(x).gt(0).sum())),
            ExtremePersons=("PersonsExtreme", "sum"),
            InformativePersons=("PersonsInformative", "sum"),
            PersonDenominator=("PersonDenominator", "sum"),
            MaximumExtremePersonsPerRun=("PersonsExtreme", "max"),
        )
        .reset_index()
    )
    cmle_extremes["ExtremePersonRate"] = (
        cmle_extremes["ExtremePersons"] / cmle_extremes["PersonDenominator"]
    )

    facets = runs[runs["EstimatorMode"].eq(FACETS_MODE)].copy()
    parity = {
        "pairs": int(len(facets)),
        "eligible": int(_as_bool(facets["ComparisonEligible"]).sum()),
        "direct_agreement_passes": int(_as_bool(facets["DirectAgreementPass"]).sum()),
        "maximum_main_weighted_mae": float(pd.to_numeric(facets["MainWeightedMAE"]).max()),
        "maximum_main_absolute_difference": float(pd.to_numeric(facets["MainMaxAbsDifference"]).max()),
        "maximum_threshold_weighted_mae": float(pd.to_numeric(facets["ThresholdWeightedMAE"]).max()),
        "maximum_threshold_absolute_difference": float(pd.to_numeric(facets["ThresholdMaxAbsDifference"]).max()),
        "minimum_within_facet_spearman": float(pd.to_numeric(facets["MinimumWithinFacetSpearman"]).min()),
    }

    frames = {
        "main_run_metrics.csv": main_run,
        "main_shape_pairs.csv": main_shape_pairs,
        "main_shape_contrast_summary.csv": main_shape_summary,
        "main_design_pairs.csv": main_design_pairs,
        "main_design_contrast_summary.csv": main_design_summary,
        "main_element_condition_summary.csv": main_element_condition,
        "main_facet_bias_summary.csv": facet_bias_summary,
        "main_element_shape_pairs.csv": main_element_shape_pairs,
        "main_element_shape_contrast_summary.csv": main_element_shape_summary,
        "main_element_design_pairs.csv": main_element_design_pairs,
        "main_element_design_contrast_summary.csv": main_element_design_summary,
        "threshold_run_metrics.csv": threshold_run,
        "threshold_shape_pairs.csv": threshold_shape_pairs,
        "threshold_shape_contrast_summary.csv": threshold_shape_summary,
        "threshold_design_pairs.csv": threshold_design_pairs,
        "threshold_design_contrast_summary.csv": threshold_design_summary,
        "threshold_element_condition_summary.csv": threshold_element_condition,
        "threshold_bias_summary.csv": threshold_bias_summary,
        "threshold_element_shape_pairs.csv": threshold_element_shape_pairs,
        "threshold_element_shape_contrast_summary.csv": threshold_element_shape_summary,
        "threshold_element_design_pairs.csv": threshold_element_design_pairs,
        "threshold_element_design_contrast_summary.csv": threshold_element_design_summary,
        "free_sd_runs.csv": free_sd,
        "free_sd_condition_summary.csv": free_sd_condition,
        "free_sd_shape_pairs.csv": free_sd_shape_pairs,
        "free_sd_shape_contrast_summary.csv": free_sd_shape_summary,
        "free_sd_design_pairs.csv": free_sd_design_pairs,
        "free_sd_design_contrast_summary.csv": free_sd_design_summary,
        "cmle_extreme_summary.csv": cmle_extremes,
    }
    for filename, frame in frames.items():
        frame.to_csv(output_dir / filename, index=False, lineterminator="\n")

    expected = {
        "main_run_facet_rows": 2400,
        "threshold_run_rows": 800,
        "main_shape_pair_rows": 1800,
        "main_design_pair_rows": 1200,
        "threshold_shape_pair_rows": 600,
        "threshold_design_pair_rows": 400,
        "free_sd_shape_pair_rows": 120,
        "free_sd_design_pair_rows": 80,
        "cmle_condition_rows": 8,
        "main_element_condition_rows": 360,
        "main_element_shape_pair_rows": 5400,
        "main_element_shape_summary_rows": 270,
        "main_element_design_pair_rows": 3600,
        "main_element_design_summary_rows": 180,
        "threshold_element_condition_rows": 240,
        "threshold_element_shape_pair_rows": 3600,
        "threshold_element_shape_summary_rows": 180,
        "threshold_element_design_pair_rows": 2400,
        "threshold_element_design_summary_rows": 120,
    }
    observed = {
        "main_run_facet_rows": len(main_run),
        "threshold_run_rows": len(threshold_run),
        "main_shape_pair_rows": len(main_shape_pairs),
        "main_design_pair_rows": len(main_design_pairs),
        "threshold_shape_pair_rows": len(threshold_shape_pairs),
        "threshold_design_pair_rows": len(threshold_design_pairs),
        "free_sd_shape_pair_rows": len(free_sd_shape_pairs),
        "free_sd_design_pair_rows": len(free_sd_design_pairs),
        "cmle_condition_rows": len(cmle_extremes),
        "main_element_condition_rows": len(main_element_condition),
        "main_element_shape_pair_rows": len(main_element_shape_pairs),
        "main_element_shape_summary_rows": len(main_element_shape_summary),
        "main_element_design_pair_rows": len(main_element_design_pairs),
        "main_element_design_summary_rows": len(main_element_design_summary),
        "threshold_element_condition_rows": len(threshold_element_condition),
        "threshold_element_shape_pair_rows": len(threshold_element_shape_pairs),
        "threshold_element_shape_summary_rows": len(threshold_element_shape_summary),
        "threshold_element_design_pair_rows": len(threshold_element_design_pairs),
        "threshold_element_design_summary_rows": len(threshold_element_design_summary),
    }
    gates = {
        "source_aggregate_qualified": bool(metrics["qualification_pass"]),
        "all_attempts_retained": bool(
            len(attempts) == 640
            and attempts["Completed"].fillna(False).astype(bool).all()
            and attempts["AttemptSucceeded"].fillna(False).astype(bool).all()
        ),
        "registered_row_counts": observed == expected,
        "paired_replicates_20": bool(
            all(
                frame["PairedReplicates"].eq(20).all()
                for frame in (
                    main_shape_summary,
                    main_design_summary,
                    threshold_shape_summary,
                    threshold_design_summary,
                    free_sd_shape_summary,
                    free_sd_design_summary,
                    main_element_shape_summary,
                    main_element_design_summary,
                    threshold_element_shape_summary,
                    threshold_element_design_summary,
                )
            )
        ),
        "facets_python_parity": bool(
            parity["pairs"] == 160
            and parity["eligible"] == 160
            and parity["direct_agreement_passes"] == 160
        ),
        "cmle_denominator_complete": bool(
            cmle_extremes["PersonDenominator"].eq(20 * 80).all()
        ),
    }
    analysis_metrics = {
        "schema_version": SCHEMA_VERSION,
        "study": study_dir.name,
        "claim_limit": "20-replicate paired screening; no estimator ranking or confirmatory inference.",
        "expected_rows": expected,
        "observed_rows": observed,
        "facets_python_parity": parity,
        "cmle_extreme_person_total": int(cmle_extremes["ExtremePersons"].sum()),
        "cmle_person_denominator": int(cmle_extremes["PersonDenominator"].sum()),
        "free_sd_global_min": float(free_sd["EstimatedPopulationSD"].min()),
        "free_sd_global_max": float(free_sd["EstimatedPopulationSD"].max()),
        "gates": gates,
        "analysis_pass": bool(all(gates.values())),
        "cross_basis_likelihood_comparison": "prohibited",
    }
    _json_dump(output_dir / "screening_analysis_metrics.json", analysis_metrics)
    identity = {
        "schema_version": SCHEMA_VERSION,
        "plan_sha256": sha256_file(PLAN_PATH),
        "analysis_amendment_sha256": sha256_file(ANALYSIS_AMENDMENT_PATH),
        "analysis_script_sha256": sha256_file(Path(__file__).resolve()),
        "source_aggregate_identity_sha256": sha256_file(aggregate_dir / "aggregate_identity.json"),
        "source_metrics_sha256": sha256_file(aggregate_dir / "study_metrics.json"),
        "output_sha256": {
            filename: sha256_file(output_dir / filename)
            for filename in (*frames.keys(), "screening_analysis_metrics.json")
        },
    }
    _json_dump(output_dir / "screening_analysis_identity.json", identity)
    print(json.dumps(analysis_metrics, indent=2, sort_keys=True))
    return analysis_metrics


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    run_analysis(args.study_dir, args.output_dir)


if __name__ == "__main__":
    main()
