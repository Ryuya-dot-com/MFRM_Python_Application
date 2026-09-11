#!/usr/bin/env python3
"""Run the frozen confirmatory analysis for the Person-distribution study.

This module only compares registered replicate-level paired contrasts. It does
not compare likelihood values or rank estimators across distinct estimands.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PLAN_PATH = (
    REPO_ROOT / "validation" / "estimand_distribution_confirmatory_plan_20260811.json"
)
EXPECTED_ENDPOINT_IDS = (
    "H1_FREE_SD_RIGHT_SKEW_PLANNED_LT_NORMAL",
    "H2_FREE_SD_HEAVY_TAIL_PLANNED_LT_NORMAL",
    "H3_JMLE_MIXTURE_PLANNED_C02_STEP2_LT_NORMAL",
    "H4_FREE_MML_MIXTURE_PLANNED_C02_STEP2_LT_NORMAL",
    "H5_EXACT_CMLE_MIXTURE_PLANNED_C02_STEP2_LT_NORMAL",
    "H6_JMLE_NORMAL_THRESHOLD_RMSE_PLANNED_GT_COMPLETE",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _included(frame: pd.DataFrame) -> pd.Series:
    if "IncludedInStudy" not in frame:
        return pd.Series(False, index=frame.index)
    values = frame["IncludedInStudy"]
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def _paired_metric(
    frame: pd.DataFrame,
    *,
    endpoint_id: str,
    metric: str,
    fixed_filters: dict[str, Any],
    pair_column: str,
    left_value: str,
    right_value: str,
) -> pd.DataFrame:
    subset = frame[_included(frame)].copy()
    for column, value in fixed_filters.items():
        subset = subset[subset[column].astype(str).eq(str(value))]
    subset[metric] = pd.to_numeric(subset[metric], errors="coerce")
    subset["Replicate"] = pd.to_numeric(subset["Replicate"], errors="raise").astype(int)
    subset = subset[subset[pair_column].astype(str).isin({left_value, right_value})]
    keys = ["Replicate", pair_column]
    if subset.duplicated(keys).any():
        raise ValueError(f"{endpoint_id} has duplicate registered cells")
    wide = subset.pivot(index="Replicate", columns=pair_column, values=metric)
    for required in (left_value, right_value):
        if required not in wide:
            wide[required] = np.nan
    contrast = wide[left_value] - wide[right_value]
    return pd.DataFrame({
        "EndpointId": endpoint_id,
        "Replicate": contrast.index.astype(int),
        "Contrast": contrast.to_numpy(dtype=float),
    })


def compute_primary_contrasts(
    runs: pd.DataFrame,
    thresholds: pd.DataFrame,
) -> pd.DataFrame:
    """Construct exactly the six registered replicate-level contrast vectors."""

    parts = [
        _paired_metric(
            runs,
            endpoint_id=EXPECTED_ENDPOINT_IDS[0],
            metric="EstimatedPopulationSD",
            fixed_filters={
                "EstimatorMode": "PYTHON_MML_FREE_SD_Q31",
                "Design": "planned_connected",
            },
            pair_column="PersonDistribution",
            left_value="right_skew",
            right_value="normal",
        ),
        _paired_metric(
            runs,
            endpoint_id=EXPECTED_ENDPOINT_IDS[1],
            metric="EstimatedPopulationSD",
            fixed_filters={
                "EstimatorMode": "PYTHON_MML_FREE_SD_Q31",
                "Design": "planned_connected",
            },
            pair_column="PersonDistribution",
            left_value="heavy_tail_t3",
            right_value="normal",
        ),
    ]
    local_specs = (
        (EXPECTED_ENDPOINT_IDS[2], "PYTHON_JMLE"),
        (EXPECTED_ENDPOINT_IDS[3], "PYTHON_MML_FREE_SD_Q31"),
        (EXPECTED_ENDPOINT_IDS[4], "PYTHON_EXACT_CMLE"),
    )
    for endpoint_id, mode in local_specs:
        parts.append(_paired_metric(
            thresholds,
            endpoint_id=endpoint_id,
            metric="TruthError",
            fixed_filters={
                "EstimatorMode": mode,
                "Design": "planned_connected",
                "StepFacetLevel": "C02",
                "Category": 2,
            },
            pair_column="PersonDistribution",
            left_value="symmetric_mixture",
            right_value="normal",
        ))

    jmle_normal = thresholds[
        _included(thresholds)
        & thresholds["EstimatorMode"].astype(str).eq("PYTHON_JMLE")
        & thresholds["PersonDistribution"].astype(str).eq("normal")
        & thresholds["Design"].astype(str).isin({"complete", "planned_connected"})
    ].copy()
    jmle_normal["TruthError"] = pd.to_numeric(jmle_normal["TruthError"], errors="coerce")
    jmle_normal["Replicate"] = pd.to_numeric(
        jmle_normal["Replicate"], errors="raise"
    ).astype(int)
    threshold_keys = ["Replicate", "Design", "StepFacetLevel", "Category"]
    if jmle_normal.duplicated(threshold_keys).any():
        raise ValueError(f"{EXPECTED_ENDPOINT_IDS[5]} has duplicate threshold elements")
    element_counts = jmle_normal.groupby(["Replicate", "Design"])["TruthError"].size()
    if len(element_counts) and not element_counts.eq(6).all():
        raise ValueError(f"{EXPECTED_ENDPOINT_IDS[5]} requires exactly six threshold elements")
    rmse = (
        jmle_normal.assign(SquaredError=jmle_normal["TruthError"].pow(2))
        .groupby(["Replicate", "Design"])["SquaredError"]
        .mean()
        .pow(0.5)
        .rename("ThresholdRMSE")
        .reset_index()
    )
    parts.append(_paired_metric(
        rmse.assign(IncludedInStudy=True),
        endpoint_id=EXPECTED_ENDPOINT_IDS[5],
        metric="ThresholdRMSE",
        fixed_filters={},
        pair_column="Design",
        left_value="planned_connected",
        right_value="complete",
    ))
    contrasts = pd.concat(parts, ignore_index=True)
    observed_ids = tuple(dict.fromkeys(contrasts["EndpointId"].astype(str)))
    if observed_ids != EXPECTED_ENDPOINT_IDS:
        raise RuntimeError("Primary endpoint construction order changed")
    return contrasts


def holm_adjust(p_values: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(p_values), dtype=float)
    if len(values) == 0 or not np.isfinite(values).all():
        raise ValueError("Holm adjustment requires finite p values")
    order = np.argsort(values, kind="mergesort")
    adjusted_sorted = np.empty(len(values), dtype=float)
    running = 0.0
    for rank, index in enumerate(order):
        candidate = (len(values) - rank) * values[index]
        running = max(running, candidate)
        adjusted_sorted[rank] = min(1.0, running)
    adjusted = np.empty(len(values), dtype=float)
    adjusted[order] = adjusted_sorted
    return adjusted


def summarize_primary(
    contrasts: pd.DataFrame,
    plan: dict[str, Any],
) -> pd.DataFrame:
    endpoints = plan["primary_family"]["endpoints"]
    if tuple(endpoint["id"] for endpoint in endpoints) != EXPECTED_ENDPOINT_IDS:
        raise ValueError("Plan primary endpoint identities do not match frozen analysis")
    required_n = int(plan["primary_family"]["required_finite_paired_replicates_per_endpoint"])
    rows: list[dict[str, Any]] = []
    for endpoint in endpoints:
        values = pd.to_numeric(
            contrasts.loc[contrasts["EndpointId"].eq(endpoint["id"]), "Contrast"],
            errors="coerce",
        ).dropna().to_numpy(dtype=float)
        n = len(values)
        mean = float(np.mean(values)) if n else np.nan
        sd = float(np.std(values, ddof=1)) if n > 1 else np.nan
        se = sd / np.sqrt(n) if n > 1 else np.nan
        statistic = mean / se if np.isfinite(se) and se > 0 else np.nan
        df = n - 1
        alternative = str(endpoint["alternative"])
        if np.isfinite(statistic):
            raw_p = (
                float(student_t.cdf(statistic, df))
                if alternative == "less"
                else float(student_t.sf(statistic, df))
            )
            critical = float(student_t.ppf(0.975, df))
            half_width = critical * se
        else:
            raw_p = np.nan
            half_width = np.nan
        direction_pass = bool(
            np.isfinite(mean)
            and ((alternative == "less" and mean < 0) or (alternative == "greater" and mean > 0))
        )
        rows.append({
            "EndpointId": endpoint["id"],
            "Alternative": alternative,
            "FinitePairedReplicates": n,
            "RequiredPairedReplicates": required_n,
            "FullPairGate": n == required_n,
            "MeanContrast": mean,
            "MonteCarloSD": sd,
            "MonteCarloSE": se,
            "TStatistic": statistic,
            "DegreesOfFreedom": df,
            "RawOneSidedP": raw_p,
            "Lower95": mean - half_width if np.isfinite(half_width) else np.nan,
            "Upper95": mean + half_width if np.isfinite(half_width) else np.nan,
            "TwoSided95HalfWidth": half_width,
            "PrecisionTarget": float(endpoint["precision_target_two_sided_95_half_width"]),
            "PrecisionPass": bool(
                np.isfinite(half_width)
                and half_width <= float(endpoint["precision_target_two_sided_95_half_width"])
            ),
            "DirectionPass": direction_pass,
        })
    output = pd.DataFrame(rows)
    finite_p = output["RawOneSidedP"].notna().all()
    output["HolmAdjustedP"] = (
        holm_adjust(output["RawOneSidedP"]) if finite_p else np.nan
    )
    output["DirectionConfirmed"] = (
        output["FullPairGate"].astype(bool)
        & output["DirectionPass"].astype(bool)
        & pd.to_numeric(output["HolmAdjustedP"], errors="coerce").le(
            float(plan["primary_family"]["familywise_alpha"])
        )
    )
    return output


def _validate_evidence_bundle(
    study_dir: Path,
    aggregate_dir: Path,
    plan_path: Path,
    plan: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    identity_path = study_dir / "study_identity.json"
    aggregate_identity_path = aggregate_dir / "aggregate_identity.json"
    metrics_path = aggregate_dir / "study_metrics.json"
    identity = json.loads(identity_path.read_text(encoding="utf-8"))
    aggregate_identity = json.loads(aggregate_identity_path.read_text(encoding="utf-8"))
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    expected_start, expected_end = plan["evidence_status"]["confirmatory_replicates"]
    expected_n = int(plan["evidence_status"]["confirmatory_replicate_count"])
    plan_hash = sha256_file(plan_path)
    registered_analysis_hash = str(
        plan.get("analysis_implementation", {})
        .get("analysis_script", {})
        .get("sha256", "")
    )
    checks = {
        "analysis_script": registered_analysis_hash == sha256_file(Path(__file__).resolve()),
        "phase": identity.get("phase") == "confirmatory" and metrics.get("phase") == "confirmatory",
        "replicate_count": int(identity.get("replicates", -1)) == expected_n,
        "replicate_range": (
            int(identity.get("replicate_start", -1)) == int(expected_start)
            and int(identity.get("replicate_end", -1)) == int(expected_end)
        ),
        "study_plan": str(identity.get("plan_sha256")) == plan_hash,
        "aggregate_plan": str(aggregate_identity.get("plan_sha256")) == plan_hash,
        "study_identity": str(aggregate_identity.get("study_identity_sha256"))
        == sha256_file(identity_path),
        "metrics": str(aggregate_identity.get("metrics_sha256")) == sha256_file(metrics_path),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError(f"Confirmatory evidence identity failed: {failed}")
    return identity, aggregate_identity, metrics


def analyze_confirmatory(
    study_dir: Path,
    *,
    aggregate_name: str = "aggregate",
    output_name: str = "confirmatory_analysis",
    plan_path: Path = DEFAULT_PLAN_PATH,
) -> dict[str, Any]:
    study_dir = study_dir.resolve()
    aggregate_dir = study_dir / aggregate_name
    output_dir = study_dir / output_name
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    plan_path = plan_path.resolve()
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    identity, aggregate_identity, study_metrics = _validate_evidence_bundle(
        study_dir, aggregate_dir, plan_path, plan
    )
    runs = pd.read_csv(aggregate_dir / "study_run_ledger.csv")
    thresholds = pd.read_csv(aggregate_dir / "study_thresholds.csv")
    contrasts = compute_primary_contrasts(runs, thresholds)
    expected_start, expected_end = plan["evidence_status"]["confirmatory_replicates"]
    observed_replicates = set(pd.to_numeric(contrasts["Replicate"], errors="raise").astype(int))
    if not observed_replicates.issubset(set(range(int(expected_start), int(expected_end) + 1))):
        raise ValueError("Primary contrasts contain a screening or unregistered replicate")
    results = summarize_primary(contrasts, plan)
    output_dir.mkdir()
    contrasts.to_csv(output_dir / "primary_contrasts.csv", index=False, lineterminator="\n")
    results.to_csv(output_dir / "primary_results.csv", index=False, lineterminator="\n")
    metrics = {
        "schema_version": "mfrm-estimand-distribution-confirmatory-analysis-v1",
        "primary_endpoints": int(len(results)),
        "full_pair_endpoints": int(results["FullPairGate"].sum()),
        "direction_confirmed_endpoints": int(results["DirectionConfirmed"].sum()),
        "precision_qualified_endpoints": int(results["PrecisionPass"].sum()),
        "all_primary_pairs_complete": bool(results["FullPairGate"].all()),
        "workbench_qualification_pass": bool(study_metrics.get("qualification_pass", False)),
        "facets_python_parity_pass": bool(
            study_metrics.get("gates", {}).get("facets_python_parity", False)
        ),
        "familywise_alpha": float(plan["primary_family"]["familywise_alpha"]),
        "multiplicity_method": "Holm",
        "estimator_ranking_performed": False,
        "cross_basis_likelihood_comparison_performed": False,
        "claim_limit": "Only the six registered independent-confirmation endpoints support confirmatory language; precision and workbench qualification are separate gates."
    }
    _json_dump(output_dir / "confirmatory_metrics.json", metrics)
    output_hashes = {
        filename: sha256_file(output_dir / filename)
        for filename in (
            "primary_contrasts.csv",
            "primary_results.csv",
            "confirmatory_metrics.json",
        )
    }
    analysis_identity = {
        "schema_version": "mfrm-estimand-distribution-confirmatory-analysis-identity-v1",
        "plan_sha256": sha256_file(plan_path),
        "analysis_script_sha256": sha256_file(Path(__file__).resolve()),
        "study_identity_sha256": sha256_file(study_dir / "study_identity.json"),
        "aggregate_identity_sha256": sha256_file(aggregate_dir / "aggregate_identity.json"),
        "frozen_aggregate_plan_sha256": aggregate_identity["plan_sha256"],
        "output_sha256": output_hashes,
    }
    _json_dump(output_dir / "confirmatory_analysis_identity.json", analysis_identity)
    print(json.dumps(metrics, indent=2, sort_keys=True))
    return metrics


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-dir", type=Path, required=True)
    parser.add_argument("--aggregate-name", default="aggregate")
    parser.add_argument("--output-name", default="confirmatory_analysis")
    parser.add_argument("--plan-path", type=Path, default=DEFAULT_PLAN_PATH)
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    analyze_confirmatory(
        args.study_dir,
        aggregate_name=args.aggregate_name,
        output_name=args.output_name,
        plan_path=args.plan_path,
    )


if __name__ == "__main__":
    main()
