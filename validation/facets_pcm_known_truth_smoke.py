#!/usr/bin/env python3
"""Run the preregistered four-run FACETS/Python known-truth PCM smoke.

This is an implementation-qualification study, not an estimator-performance
study.  It retains one independently generated input bundle, sends identical
rows to FACETS 4.5 and the application JMLE, and evaluates only the locked
direct-agreement and direction gates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import spearmanr


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from validation.facets_pcm_probe import (  # noqa: E402
    AUXILIARY_DECIMALS,
    build_pcm_spec,
    map_pcm_scale_tables,
    normalize_python_pcm_steps,
    pcm_category_count_audit,
    pcm_threshold_pairs,
    python_pcm_fit,
)
from validation.operating_characteristics_facets import (  # noqa: E402
    FACET_COLUMNS,
    RECOVERY_FACETS,
    align_recovery,
    invoke_facets,
    parse_iteration_report,
    parse_score_file,
    parse_table8_categories,
    sha256_file,
    slugify,
    validate_bundle,
)


SCHEMA_VERSION = "mfrm-facets-pcm-known-truth-smoke-v1"
SEEDS = (731101, 731102)
N_PERSONS = 80
CATEGORIES = 4
RATER_TRUTH = {"R01": -0.45, "R02": -0.15, "R03": 0.15, "R04": 0.45}
TASK_TRUTH = {"T01": -0.35, "T02": 0.0, "T03": 0.35}
CRITERION_TRUTH = {"C01": -0.20, "C02": 0.20}
THRESHOLD_CONDITIONS = {
    "shared": {
        "C01": np.array([-1.10, 0.00, 1.10]),
        "C02": np.array([-1.10, 0.00, 1.10]),
    },
    "heterogeneous": {
        "C01": np.array([-1.60, 0.25, 1.35]),
        "C02": np.array([-0.45, -0.35, 0.80]),
    },
}
MAIN_MAE_MAX = 0.01
MAIN_ABS_MAX = 0.05
MIN_SPEARMAN = 0.999
THRESHOLD_MAE_MAX = 0.02
THRESHOLD_ABS_MAX = 0.05


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _frame_sha256(frame: pd.DataFrame, columns: list[str]) -> str:
    normalized = frame.loc[:, columns].copy()
    return _sha256_text(normalized.to_csv(index=False, lineterminator="\n"))


def pcm_probabilities(eta: float, thresholds: np.ndarray) -> np.ndarray:
    """Adjacent-category PCM probabilities under the locked sign convention."""

    thresholds = np.asarray(thresholds, dtype=float)
    cumulative = np.r_[0.0, np.cumsum(thresholds)]
    logits = np.arange(len(cumulative), dtype=float) * float(eta) - cumulative
    return np.exp(logits - logsumexp(logits))


def generate_latent_replicate(seed: int, *, n_persons: int = N_PERSONS) -> tuple[pd.DataFrame, dict[str, float]]:
    """Create the paired person measures and row uniforms for one replicate."""

    rng = np.random.default_rng(seed)
    person_values = rng.normal(loc=0.0, scale=0.80, size=n_persons)
    person_values = np.clip(person_values, -2.2, 2.2)
    person_values -= person_values.mean()
    person_truth = {f"P{index + 1:03d}": float(value) for index, value in enumerate(person_values)}
    rows: list[dict[str, Any]] = []
    for person, theta in person_truth.items():
        for rater, rater_value in RATER_TRUTH.items():
            for task, task_value in TASK_TRUTH.items():
                for criterion, criterion_value in CRITERION_TRUTH.items():
                    rows.append({
                        "Person": person,
                        "Rater": rater,
                        "Task": task,
                        "Criterion": criterion,
                        "Theta": theta,
                        "Eta": theta - rater_value - task_value - criterion_value,
                        "Uniform": float(rng.random()),
                    })
    return pd.DataFrame(rows), person_truth


def apply_threshold_condition(
    latent: pd.DataFrame,
    thresholds: dict[str, np.ndarray],
) -> pd.DataFrame:
    """Turn a paired latent table into scores without drawing new randomness."""

    rows = latent.copy()
    scores: list[int] = []
    for row in rows.itertuples(index=False):
        probabilities = pcm_probabilities(row.Eta, thresholds[str(row.Criterion)])
        score = int(np.searchsorted(np.cumsum(probabilities), float(row.Uniform), side="right"))
        scores.append(min(score, len(probabilities) - 1))
    rows["Score"] = scores
    return rows[["Person", "Rater", "Task", "Criterion", "Score"]]


def generate_bundle(
    *,
    seeds: Iterable[int] = SEEDS,
    n_persons: int = N_PERSONS,
) -> dict[str, pd.DataFrame]:
    manifests: list[dict[str, Any]] = []
    ratings_parts: list[pd.DataFrame] = []
    truth_parts: list[pd.DataFrame] = []
    threshold_parts: list[pd.DataFrame] = []
    pairing_parts: list[pd.DataFrame] = []
    for replicate, seed in enumerate(tuple(seeds), start=1):
        latent, person_truth = generate_latent_replicate(seed, n_persons=n_persons)
        latent_key = _frame_sha256(
            latent,
            ["Person", "Rater", "Task", "Criterion", "Theta", "Eta", "Uniform"],
        )
        for condition, thresholds in THRESHOLD_CONDITIONS.items():
            run_id = f"pcm_{condition}::rep-{replicate:05d}"
            manifests.append({
                "RunId": run_id,
                "ConditionId": f"pcm_{condition}",
                "Design": "complete_balanced_pcm_known_truth",
                "TruthBias": 0.0,
                "Replicate": replicate,
                "Seed": seed,
                "Categories": CATEGORIES,
                "ThresholdCondition": condition,
                "PairedLatentSHA256": latent_key,
            })
            run_ratings = apply_threshold_condition(latent, thresholds)
            run_ratings.insert(0, "RunId", run_id)
            ratings_parts.append(run_ratings)
            facet_rows = [
                {"RunId": run_id, "Facet": "Person", "Level": level, "Truth": value}
                for level, value in person_truth.items()
            ]
            for facet, values in (
                ("Rater", RATER_TRUTH),
                ("Task", TASK_TRUTH),
                ("Criterion", CRITERION_TRUTH),
            ):
                facet_rows.extend(
                    {"RunId": run_id, "Facet": facet, "Level": level, "Truth": value}
                    for level, value in values.items()
                )
            truth_parts.append(pd.DataFrame(facet_rows))
            for criterion, vector in thresholds.items():
                for category, value in enumerate(vector, start=1):
                    threshold_parts.append(pd.DataFrame([{
                        "RunId": run_id,
                        "ConditionId": f"pcm_{condition}",
                        "Replicate": replicate,
                        "Seed": seed,
                        "StepFacetLevel": criterion,
                        "Category": category,
                        "ThresholdTruth": float(value),
                    }]))
            pairing = latent[["Person", "Rater", "Task", "Criterion", "Uniform"]].copy()
            pairing.insert(0, "RunId", run_id)
            pairing_parts.append(pairing)
    return {
        "manifest.csv": pd.DataFrame(manifests),
        "generated_ratings.csv": pd.concat(ratings_parts, ignore_index=True),
        "generated_facet_truth.csv": pd.concat(truth_parts, ignore_index=True),
        "generated_anchors.csv": pd.DataFrame(columns=["RunId", "Facet", "Level", "Anchor"]),
        "generated_pcm_threshold_truth.csv": pd.concat(threshold_parts, ignore_index=True),
        "paired_uniforms.csv": pd.concat(pairing_parts, ignore_index=True),
    }


def write_bundle(bundle: dict[str, pd.DataFrame], input_dir: Path) -> None:
    input_dir.mkdir(parents=True, exist_ok=False)
    for filename, frame in bundle.items():
        frame.to_csv(input_dir / filename, index=False, lineterminator="\n")


def audit_parent_extension(input_dir: Path, parent_dir: Path) -> pd.DataFrame:
    """Verify byte-stable per-run content for every parent RunId."""

    child_manifest = pd.read_csv(input_dir / "manifest.csv")
    parent_manifest = pd.read_csv(parent_dir / "manifest.csv")
    child_ratings = pd.read_csv(input_dir / "generated_ratings.csv")
    parent_ratings = pd.read_csv(parent_dir / "generated_ratings.csv")
    child_thresholds = pd.read_csv(input_dir / "generated_pcm_threshold_truth.csv")
    parent_thresholds = pd.read_csv(parent_dir / "generated_pcm_threshold_truth.csv")
    rows = []
    for run_id in parent_manifest["RunId"].astype(str):
        child_manifest_row = child_manifest[child_manifest["RunId"].astype(str).eq(run_id)]
        parent_manifest_row = parent_manifest[parent_manifest["RunId"].astype(str).eq(run_id)]
        manifest_match = bool(
            len(child_manifest_row) == len(parent_manifest_row) == 1
            and child_manifest_row[parent_manifest.columns].reset_index(drop=True).equals(
                parent_manifest_row.reset_index(drop=True)
            )
        )
        rating_columns = ["RunId", "Person", "Rater", "Task", "Criterion", "Score"]
        child_rating_hash = _frame_sha256(
            child_ratings[child_ratings["RunId"].astype(str).eq(run_id)], rating_columns
        )
        parent_rating_hash = _frame_sha256(
            parent_ratings[parent_ratings["RunId"].astype(str).eq(run_id)], rating_columns
        )
        threshold_columns = [
            "RunId", "ConditionId", "Replicate", "Seed", "StepFacetLevel", "Category", "ThresholdTruth"
        ]
        child_threshold_hash = _frame_sha256(
            child_thresholds[child_thresholds["RunId"].astype(str).eq(run_id)], threshold_columns
        )
        parent_threshold_hash = _frame_sha256(
            parent_thresholds[parent_thresholds["RunId"].astype(str).eq(run_id)], threshold_columns
        )
        rows.append({
            "RunId": run_id,
            "ManifestMatch": manifest_match,
            "RatingsMatch": child_rating_hash == parent_rating_hash,
            "ThresholdTruthMatch": child_threshold_hash == parent_threshold_hash,
            "ChildRatingSHA256": child_rating_hash,
            "ParentRatingSHA256": parent_rating_hash,
            "ChildThresholdSHA256": child_threshold_hash,
            "ParentThresholdSHA256": parent_threshold_hash,
        })
    audit = pd.DataFrame(rows)
    audit["AllMatch"] = audit[["ManifestMatch", "RatingsMatch", "ThresholdTruthMatch"]].all(axis=1)
    if not audit["AllMatch"].all():
        raise ValueError("Expanded PCM pilot does not preserve every parent known-truth run")
    return audit


def expanded_pcm_design_audit(ratings: pd.DataFrame) -> dict[str, Any]:
    """Rank-check the constrained adjacent-category PCM design."""

    levels = {
        facet: list(dict.fromkeys(ratings[facet].astype(str)))
        for facet in ("Person", "Rater", "Task", "Criterion")
    }
    n_steps = CATEGORIES - 1
    columns: list[str] = [f"Person::{level}" for level in levels["Person"]]
    slices: dict[str, slice] = {"Person": slice(0, len(columns))}
    for facet in ("Rater", "Task", "Criterion"):
        start = len(columns)
        columns.extend(f"{facet}::{level}" for level in levels[facet][:-1])
        slices[facet] = slice(start, len(columns))
    threshold_slices: dict[str, slice] = {}
    for criterion in levels["Criterion"]:
        start = len(columns)
        columns.extend(f"Threshold::{criterion}::{step}" for step in range(1, n_steps))
        threshold_slices[criterion] = slice(start, len(columns))

    def centered_row(level: str, ordered: list[str]) -> np.ndarray:
        result = np.zeros(len(ordered) - 1, dtype=float)
        index = ordered.index(level)
        if index < len(ordered) - 1:
            result[index] = 1.0
        else:
            result[:] = -1.0
        return result

    matrix = np.zeros((len(ratings) * n_steps, len(columns)), dtype=float)
    out_row = 0
    for row in ratings.itertuples(index=False):
        person_index = levels["Person"].index(str(row.Person))
        for transition in range(1, n_steps + 1):
            matrix[out_row, person_index] = 1.0
            for facet in ("Rater", "Task", "Criterion"):
                code = centered_row(str(getattr(row, facet)), levels[facet])
                matrix[out_row, slices[facet]] = -code
            threshold_slice = threshold_slices[str(row.Criterion)]
            if transition < n_steps:
                matrix[out_row, threshold_slice.start + transition - 1] = -1.0
            else:
                matrix[out_row, threshold_slice] = 1.0
            out_row += 1
    rank = int(np.linalg.matrix_rank(matrix))
    return {
        "Rows": int(matrix.shape[0]),
        "Columns": int(matrix.shape[1]),
        "Rank": rank,
        "Nullity": int(matrix.shape[1] - rank),
        "ConstraintContract": (
            "Person noncentered; Rater/Task/Criterion and each Criterion threshold vector sum-to-zero"
        ),
    }


def _observation_weights(ratings: pd.DataFrame) -> pd.DataFrame:
    parts = []
    for facet in RECOVERY_FACETS:
        part = ratings.groupby(facet, observed=False).size().rename("ObservationWeight").reset_index()
        part = part.rename(columns={facet: "Level"})
        part.insert(0, "Facet", facet)
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def _weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    numeric_weights = pd.to_numeric(weights, errors="coerce").to_numpy(dtype=float)
    finite = np.isfinite(numeric) & np.isfinite(numeric_weights) & (numeric_weights > 0)
    return float(np.average(numeric[finite], weights=numeric_weights[finite])) if finite.any() else float("nan")


def fit_run(
    manifest_row: pd.Series,
    *,
    ratings: pd.DataFrame,
    truth: pd.DataFrame,
    threshold_truth: pd.DataFrame,
    facets_exe: Path,
    run_root: Path,
    timeout_seconds: float,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    run_id = str(manifest_row["RunId"])
    run_dir = run_root / slugify(run_id)
    run_dir.mkdir(parents=True, exist_ok=False)
    spec_path = run_dir / "analysis.txt"
    report_path = run_dir / "report_u6.txt"
    aux_report_path = run_dir / "report_u2.txt"
    score_base = run_dir / "scores.txt"
    aux_score_base = run_dir / "scores_u2.txt"
    spec, level_maps = build_pcm_spec(manifest_row, ratings, pd.DataFrame(), score_base=score_base)
    spec_path.write_text(spec, encoding="utf-8", newline="\n")
    started = time.perf_counter()
    primary = invoke_facets(
        facets_exe, spec_path, report_path, timeout_seconds=timeout_seconds
    )
    (run_dir / "stdout_u6.txt").write_text(primary.stdout or "", encoding="utf-8")
    (run_dir / "stderr_u6.txt").write_text(primary.stderr or "", encoding="utf-8")
    if primary.returncode != 0 or not report_path.is_file():
        raise RuntimeError(f"FACETS primary pass failed for {run_id}: exit={primary.returncode}")
    auxiliary = invoke_facets(
        facets_exe,
        spec_path,
        aux_report_path,
        timeout_seconds=timeout_seconds,
        extra_specs=(
            f"Umean=0,1,{AUXILIARY_DECIMALS}",
            f"Scorefile={aux_score_base}",
        ),
    )
    (run_dir / "stdout_u2.txt").write_text(auxiliary.stdout or "", encoding="utf-8")
    (run_dir / "stderr_u2.txt").write_text(auxiliary.stderr or "", encoding="utf-8")
    if auxiliary.returncode != 0 or not aux_report_path.is_file():
        raise RuntimeError(f"FACETS auxiliary pass failed for {run_id}: exit={auxiliary.returncode}")

    report_info = parse_iteration_report(report_path)
    scores = pd.concat([
        parse_score_file(
            run_dir / f"scores.{facet_number}.txt",
            facet_number=facet_number,
            facet_name=facet_name,
        )
        for facet_number, facet_name in FACET_COLUMNS
    ], ignore_index=True)
    facets_recovery = align_recovery(scores, truth, pd.DataFrame(), manifest_row)
    primary_categories = map_pcm_scale_tables(
        parse_table8_categories(report_path), criterion_map=level_maps["Criterion"]
    )
    auxiliary_categories = map_pcm_scale_tables(
        parse_table8_categories(aux_report_path), criterion_map=level_maps["Criterion"]
    )
    count_audit = pcm_category_count_audit(auxiliary_categories, ratings)

    python_result = python_pcm_fit(ratings, pd.DataFrame(), categories=CATEGORIES)
    python_summary = python_result["summary"].iloc[0]
    python_estimates = python_result["facets"]["others"].copy()
    python_estimates["SE"] = np.nan
    python_estimates["Status"] = np.nan
    python_recovery = align_recovery(python_estimates, truth, pd.DataFrame(), manifest_row)
    keys = ["RunId", "Facet", "Level"]
    recovery_pairs = facets_recovery.merge(
        python_recovery[keys + ["EstimateAligned", "ErrorAligned"]].rename(columns={
            "EstimateAligned": "PythonEstimateAligned",
            "ErrorAligned": "PythonTruthErrorAligned",
        }),
        on=keys,
        validate="one_to_one",
    )
    recovery_pairs = recovery_pairs.rename(columns={
        "EstimateAligned": "FACETSEstimateAligned",
        "ErrorAligned": "FACETSTruthErrorAligned",
    })
    recovery_pairs["EstimateDifference"] = (
        recovery_pairs["FACETSEstimateAligned"] - recovery_pairs["PythonEstimateAligned"]
    )
    recovery_pairs["AbsoluteEstimateDifference"] = recovery_pairs["EstimateDifference"].abs()
    recovery_pairs = recovery_pairs.merge(
        _observation_weights(ratings), on=["Facet", "Level"], validate="one_to_one"
    )

    threshold_pairs = pcm_threshold_pairs(
        primary_categories, normalize_python_pcm_steps(python_result)
    ).merge(
        threshold_truth[["StepFacetLevel", "Category", "ThresholdTruth"]],
        on=["StepFacetLevel", "Category"],
        validate="one_to_one",
    )
    threshold_pairs["FACETSTruthError"] = (
        threshold_pairs["ThresholdMeasureDisplayed"] - threshold_pairs["ThresholdTruth"]
    )
    threshold_pairs["PythonTruthError"] = (
        threshold_pairs["PythonThreshold"] - threshold_pairs["ThresholdTruth"]
    )
    criterion_weights = ratings.groupby("Criterion", observed=False).size().to_dict()
    threshold_pairs["ObservationWeight"] = threshold_pairs["StepFacetLevel"].map(criterion_weights)

    design = expanded_pcm_design_audit(ratings)
    common_eligible = bool(
        report_info.get("Converged", False)
        and python_summary.get("Converged", False)
        and python_summary.get("InferenceReady", False)
        and design["Nullity"] == 0
        and count_audit["CountMatches"].all()
    )
    recovery_pairs["ComparisonEligible"] = common_eligible
    threshold_pairs["ComparisonEligible"] = common_eligible
    correlations: dict[str, float] = {}
    for facet, group in recovery_pairs.groupby("Facet", sort=False):
        value = spearmanr(
            group["FACETSEstimateAligned"], group["PythonEstimateAligned"]
        ).statistic
        correlations[str(facet)] = float(value) if np.isfinite(value) else float("nan")
    direct_main = recovery_pairs[recovery_pairs["ComparisonEligible"]]
    direct_threshold = threshold_pairs[threshold_pairs["ComparisonEligible"]]
    run_metrics = {
        "RunId": run_id,
        "ConditionId": manifest_row["ConditionId"],
        "Replicate": int(manifest_row["Replicate"]),
        "Seed": int(manifest_row["Seed"]),
        "Rows": int(len(ratings)),
        "RatingRowSHA256": _frame_sha256(
            ratings, ["Person", "Rater", "Task", "Criterion", "Score"]
        ),
        "FACETSVersion": report_info.get("FacetsVersion"),
        "FACETSConverged": bool(report_info.get("Converged", False)),
        "PythonConverged": bool(python_summary.get("Converged", False)),
        "PythonInferenceReady": bool(python_summary.get("InferenceReady", False)),
        "ExpandedPCMDesignNullity": int(design["Nullity"]),
        "Table8ScaleTables": int(primary_categories["TableNumber"].nunique()),
        "CategoryCountRows": int(len(count_audit)),
        "CategoryCountMatches": int(count_audit["CountMatches"].sum()),
        "ComparisonEligible": common_eligible,
        "DirectMainParameters": int(len(direct_main)),
        "MainWeightedMAE": _weighted_mean(
            direct_main["AbsoluteEstimateDifference"], direct_main["ObservationWeight"]
        ),
        "MainMaxAbsDifference": float(direct_main["AbsoluteEstimateDifference"].max()) if len(direct_main) else np.nan,
        "MinimumWithinFacetSpearman": float(min(correlations.values())) if correlations else np.nan,
        "WithinFacetSpearman": json.dumps(correlations, sort_keys=True),
        "DirectThresholds": int(len(direct_threshold)),
        "ThresholdWeightedMAE": _weighted_mean(
            direct_threshold["AbsoluteThresholdDifference"], direct_threshold["ObservationWeight"]
        ),
        "ThresholdMaxAbsDifference": float(direct_threshold["AbsoluteThresholdDifference"].max()) if len(direct_threshold) else np.nan,
        "ElapsedSeconds": time.perf_counter() - started,
        "SpecificationSHA256": sha256_file(spec_path),
        "PrimaryReportSHA256": sha256_file(report_path),
        "AuxiliaryReportSHA256": sha256_file(aux_report_path),
    }
    run_metrics["DirectAgreementPass"] = bool(
        common_eligible
        and run_metrics["MainWeightedMAE"] <= MAIN_MAE_MAX
        and run_metrics["MainMaxAbsDifference"] <= MAIN_ABS_MAX
        and run_metrics["MinimumWithinFacetSpearman"] >= MIN_SPEARMAN
        and run_metrics["ThresholdWeightedMAE"] <= THRESHOLD_MAE_MAX
        and run_metrics["ThresholdMaxAbsDifference"] <= THRESHOLD_ABS_MAX
    )
    primary_categories.to_csv(run_dir / "table8_u6.csv", index=False)
    auxiliary_categories.to_csv(run_dir / "table8_u2.csv", index=False)
    count_audit.to_csv(run_dir / "category_count_audit.csv", index=False)
    recovery_pairs.to_csv(run_dir / "facets_python_main_pairs.csv", index=False)
    threshold_pairs.to_csv(run_dir / "facets_python_threshold_pairs.csv", index=False)
    (run_dir / "expanded_pcm_design_audit.json").write_text(
        json.dumps(design, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return run_metrics, recovery_pairs, threshold_pairs, count_audit, primary_categories


def _threshold_distance(group: pd.DataFrame, estimate_column: str) -> float:
    pivot = group.pivot(index="Category", columns="StepFacetLevel", values=estimate_column)
    if set(pivot.columns) != set(CRITERION_TRUTH):
        return float("nan")
    difference = pivot["C01"] - pivot["C02"]
    return float(np.sqrt(np.mean(np.square(difference))))


def run_smoke(args: argparse.Namespace) -> None:
    output_dir = args.output_dir.resolve()
    facets_exe = args.facets_exe.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    if not facets_exe.is_file():
        raise FileNotFoundError(f"FACETS executable is missing: {facets_exe}")
    output_dir.mkdir(parents=True)
    input_dir = output_dir / "retained_input"
    if args.replicates < 1:
        raise ValueError("--replicates must be positive")
    seeds = tuple(int(args.seed_start) + index for index in range(int(args.replicates)))
    write_bundle(generate_bundle(seeds=seeds, n_persons=int(args.persons)), input_dir)
    parent_audit = pd.DataFrame()
    if args.parent_dir is not None:
        parent_audit = audit_parent_extension(
            input_dir, args.parent_dir.resolve() / "retained_input"
        )
        parent_audit.to_csv(output_dir / "pcm_parent_extension_audit.csv", index=False)
    tables = validate_bundle(input_dir)
    manifest = tables["manifest.csv"]
    ratings_all = tables["generated_ratings.csv"]
    truth_all = tables["generated_facet_truth.csv"]
    anchors_all = tables["generated_anchors.csv"]
    threshold_truth_all = pd.read_csv(input_dir / "generated_pcm_threshold_truth.csv")
    if not anchors_all.empty:
        raise ValueError("Known-truth PCM smoke unexpectedly contains anchors")
    work_root = output_dir / "facets_runs"
    work_root.mkdir()
    run_rows: list[dict[str, Any]] = []
    main_parts: list[pd.DataFrame] = []
    threshold_parts: list[pd.DataFrame] = []
    count_parts: list[pd.DataFrame] = []
    category_parts: list[pd.DataFrame] = []
    for _, manifest_row in manifest.iterrows():
        run_id = str(manifest_row["RunId"])
        ratings = ratings_all[ratings_all["RunId"].astype(str).eq(run_id)].drop(columns="RunId")
        truth = truth_all[truth_all["RunId"].astype(str).eq(run_id)]
        threshold_truth = threshold_truth_all[
            threshold_truth_all["RunId"].astype(str).eq(run_id)
        ]
        try:
            run, main, threshold, counts, categories = fit_run(
                manifest_row,
                ratings=ratings,
                truth=truth,
                threshold_truth=threshold_truth,
                facets_exe=facets_exe,
                run_root=work_root,
                timeout_seconds=args.timeout_seconds,
            )
            run["Attempted"] = True
            run["FailureReason"] = ""
            run_rows.append(run)
            main_parts.append(main)
            threshold_parts.append(threshold.assign(
                RunId=run_id,
                Replicate=manifest_row["Replicate"],
                ConditionId=manifest_row["ConditionId"],
            ))
            count_parts.append(counts.assign(RunId=run_id))
            category_parts.append(categories.assign(RunId=run_id))
        except Exception as exc:  # retain every attempted run; continue fail-closed
            run_rows.append({
                "RunId": run_id,
                "ConditionId": manifest_row["ConditionId"],
                "Replicate": int(manifest_row["Replicate"]),
                "Seed": int(manifest_row["Seed"]),
                "Rows": int(len(ratings)),
                "Attempted": True,
                "FACETSConverged": False,
                "PythonConverged": False,
                "PythonInferenceReady": False,
                "ExpandedPCMDesignNullity": np.nan,
                "Table8ScaleTables": 0,
                "CategoryCountRows": 0,
                "CategoryCountMatches": 0,
                "ComparisonEligible": False,
                "DirectAgreementPass": False,
                "FailureReason": f"{type(exc).__name__}: {exc}",
            })

    runs = pd.DataFrame(run_rows)
    main_pairs = pd.concat(main_parts, ignore_index=True) if main_parts else pd.DataFrame()
    threshold_pairs = pd.concat(threshold_parts, ignore_index=True) if threshold_parts else pd.DataFrame()
    counts = pd.concat(count_parts, ignore_index=True) if count_parts else pd.DataFrame()
    categories = pd.concat(category_parts, ignore_index=True) if category_parts else pd.DataFrame()
    eligible_main = (
        main_pairs[main_pairs["ComparisonEligible"]] if not main_pairs.empty else main_pairs
    )
    eligible_thresholds = (
        threshold_pairs[threshold_pairs["ComparisonEligible"]]
        if not threshold_pairs.empty else threshold_pairs
    )
    direction_rows = []
    for replicate in sorted(manifest["Replicate"].unique()):
        replicate_rows = (
            threshold_pairs[threshold_pairs["Replicate"].eq(replicate)]
            if not threshold_pairs.empty else threshold_pairs
        )
        row: dict[str, Any] = {"Replicate": int(replicate)}
        for engine_column, engine in (
            ("ThresholdMeasureDisplayed", "FACETS"),
            ("PythonThreshold", "Python"),
        ):
            distances = {}
            if not replicate_rows.empty:
                for condition_id, condition_rows in replicate_rows.groupby("ConditionId"):
                    distances[str(condition_id)] = _threshold_distance(condition_rows, engine_column)
            shared = distances.get("pcm_shared", np.nan)
            heterogeneous = distances.get("pcm_heterogeneous", np.nan)
            row[f"{engine}SharedRMSDistance"] = shared
            row[f"{engine}HeterogeneousRMSDistance"] = heterogeneous
            row[f"{engine}DirectionPass"] = bool(
                np.isfinite(shared) and np.isfinite(heterogeneous) and heterogeneous > shared
            )
        direction_rows.append(row)
    directions = pd.DataFrame(direction_rows)

    direct_passes = int(runs.get("DirectAgreementPass", pd.Series(dtype=bool)).fillna(False).sum())
    direct_pass_rate = direct_passes / len(runs) if len(runs) else float("nan")
    comparison_eligible_rate = float(runs["ComparisonEligible"].fillna(False).mean())
    expected_count_rows = len(runs) * len(CRITERION_TRUTH) * CATEGORIES
    count_matches = int(counts["CountMatches"].sum()) if not counts.empty else 0
    all_count_rows_present = bool(len(counts) == expected_count_rows)
    direction_passes = int(
        directions[["FACETSDirectionPass", "PythonDirectionPass"]].sum().sum()
    )
    direction_checks = int(len(directions) * 2)
    direction_pass_rate = direction_passes / direction_checks if direction_checks else float("nan")
    aggregate = {
        "schema_version": SCHEMA_VERSION,
        "study_depth": "known_truth_smoke" if int(args.replicates) == 2 else "pipeline_pilot",
        "claim_limit": (
            f"{len(runs)}-run known-truth implementation qualification; "
            "no bias, coverage, false-positive, or power claim."
        ),
        "runs": int(len(runs)),
        "facets_version": sorted(runs["FACETSVersion"].dropna().astype(str).unique().tolist()),
        "all_facets_converged": bool(runs["FACETSConverged"].all()),
        "all_python_converged": bool(runs["PythonConverged"].all()),
        "all_python_inference_ready": bool(runs["PythonInferenceReady"].all()),
        "all_expanded_pcm_nullity_zero": bool(
            runs.loc[runs["ComparisonEligible"].fillna(False), "ExpandedPCMDesignNullity"].eq(0).all()
        ),
        "all_scale_tables_mapped": bool(runs["Table8ScaleTables"].eq(len(CRITERION_TRUTH)).all()),
        "category_count_rows": int(len(counts)),
        "expected_category_count_rows": int(expected_count_rows),
        "category_count_matches": count_matches,
        "all_category_counts_match": bool(
            all_count_rows_present and not counts.empty and counts["CountMatches"].all()
        ),
        "all_runs_comparison_eligible": bool(runs["ComparisonEligible"].all()),
        "comparison_eligible_rate": comparison_eligible_rate,
        "direct_main_parameters": int(len(eligible_main)),
        "main_weighted_mae_logits": (
            _weighted_mean(eligible_main["AbsoluteEstimateDifference"], eligible_main["ObservationWeight"])
            if not eligible_main.empty else None
        ),
        "main_max_abs_difference_logits": (
            float(eligible_main["AbsoluteEstimateDifference"].max()) if not eligible_main.empty else None
        ),
        "minimum_within_run_facet_spearman": (
            float(runs["MinimumWithinFacetSpearman"].min())
            if "MinimumWithinFacetSpearman" in runs and runs["MinimumWithinFacetSpearman"].notna().any()
            else None
        ),
        "direct_thresholds": int(len(eligible_thresholds)),
        "threshold_weighted_mae_logits": (
            _weighted_mean(eligible_thresholds["AbsoluteThresholdDifference"], eligible_thresholds["ObservationWeight"])
            if not eligible_thresholds.empty else None
        ),
        "threshold_max_abs_difference_logits": (
            float(eligible_thresholds["AbsoluteThresholdDifference"].max())
            if not eligible_thresholds.empty else None
        ),
        "paired_direction_checks": direction_checks,
        "paired_direction_passes": direction_passes,
        "paired_direction_pass_rate": direction_pass_rate,
        "per_run_direct_agreement_passes": direct_passes,
        "per_run_direct_agreement_pass_rate": direct_pass_rate,
        "parent_runs_audited": int(len(parent_audit)),
        "parent_runs_preserved": int(parent_audit["AllMatch"].sum()) if len(parent_audit) else 0,
    }
    gates = {
        "return_convergence_and_mapping": bool(
            aggregate["all_facets_converged"]
            and aggregate["all_python_converged"]
            and aggregate["all_python_inference_ready"]
            and aggregate["all_expanded_pcm_nullity_zero"]
            and aggregate["all_scale_tables_mapped"]
            and aggregate["all_category_counts_match"]
            and aggregate["all_runs_comparison_eligible"]
        ),
        "main_weighted_mae": bool(
            aggregate["main_weighted_mae_logits"] is not None
            and aggregate["main_weighted_mae_logits"] <= MAIN_MAE_MAX
        ),
        "main_max_abs_difference": bool(
            aggregate["main_max_abs_difference_logits"] is not None
            and aggregate["main_max_abs_difference_logits"] <= MAIN_ABS_MAX
        ),
        "within_facet_spearman": bool(
            aggregate["minimum_within_run_facet_spearman"] is not None
            and aggregate["minimum_within_run_facet_spearman"] >= MIN_SPEARMAN
        ),
        "threshold_weighted_mae": bool(
            aggregate["threshold_weighted_mae_logits"] is not None
            and aggregate["threshold_weighted_mae_logits"] <= THRESHOLD_MAE_MAX
        ),
        "threshold_max_abs_difference": bool(
            aggregate["threshold_max_abs_difference_logits"] is not None
            and aggregate["threshold_max_abs_difference_logits"] <= THRESHOLD_ABS_MAX
        ),
        "paired_heterogeneity_direction": bool(
            aggregate["paired_direction_passes"] == aggregate["paired_direction_checks"]
        ),
    }
    aggregate["locked_tolerances"] = {
        "main_weighted_mae_max": MAIN_MAE_MAX,
        "main_absolute_difference_max": MAIN_ABS_MAX,
        "minimum_spearman": MIN_SPEARMAN,
        "threshold_weighted_mae_max": THRESHOLD_MAE_MAX,
        "threshold_absolute_difference_max": THRESHOLD_ABS_MAX,
    }
    aggregate["gates"] = gates
    if int(args.replicates) == 2:
        qualification_gates = gates
    else:
        qualification_gates = {
            "parent_extension": bool(
                len(parent_audit) > 0 and parent_audit["AllMatch"].all()
            ),
            "input_and_counts": aggregate["all_category_counts_match"],
            "identification": aggregate["all_expanded_pcm_nullity_zero"],
            "engine_accounting_95pct": comparison_eligible_rate >= 0.95,
            "per_run_direct_agreement_95pct": direct_pass_rate >= 0.95,
            "aggregate_direct_tolerances": bool(
                gates["main_weighted_mae"]
                and gates["main_max_abs_difference"]
                and gates["within_facet_spearman"]
                and gates["threshold_weighted_mae"]
                and gates["threshold_max_abs_difference"]
            ),
            "paired_direction_90pct": direction_pass_rate >= 0.90,
        }
    aggregate["qualification_gates"] = qualification_gates
    aggregate["qualification_pass"] = bool(all(qualification_gates.values()))

    runs.to_csv(output_dir / "pcm_known_truth_run_metrics.csv", index=False)
    main_pairs.to_csv(output_dir / "pcm_known_truth_main_pairs.csv", index=False)
    threshold_pairs.to_csv(output_dir / "pcm_known_truth_threshold_pairs.csv", index=False)
    counts.to_csv(output_dir / "pcm_known_truth_count_audit.csv", index=False)
    categories.to_csv(output_dir / "pcm_known_truth_table8_u6.csv", index=False)
    directions.to_csv(output_dir / "pcm_known_truth_direction_checks.csv", index=False)
    (output_dir / "pcm_known_truth_metrics.json").write_text(
        json.dumps(aggregate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    identity = {
        "schema_version": SCHEMA_VERSION,
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "probe_adapter_sha256": sha256_file(REPO_ROOT / "validation" / "facets_pcm_probe.py"),
        "app_sha256": sha256_file(REPO_ROOT / "streamlit_app.py"),
        "facets_executable": str(facets_exe),
        "facets_executable_sha256": sha256_file(facets_exe),
        "python_version": sys.version,
        "platform": platform.platform(),
        "retained_input_sha256": {
            path.name: sha256_file(path) for path in sorted(input_dir.glob("*.csv"))
        },
        "input_contract": "Each engine consumed the same per-RunId rows read back from retained_input/generated_ratings.csv.",
    }
    (output_dir / "pcm_known_truth_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--facets-exe", type=Path, default=Path(r"C:\Facets\Facets.exe"))
    parser.add_argument("--timeout-seconds", type=float, default=120.0)
    parser.add_argument("--replicates", type=int, default=2)
    parser.add_argument("--seed-start", type=int, default=SEEDS[0])
    parser.add_argument("--persons", type=int, default=N_PERSONS)
    parser.add_argument("--parent-dir", type=Path)
    return parser.parse_args(argv)


def main() -> None:
    run_smoke(parse_args())


if __name__ == "__main__":
    main()
