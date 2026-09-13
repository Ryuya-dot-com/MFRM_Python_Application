#!/usr/bin/env python3
"""Run the preregistered PCM missingness/misspecification boundary pilot.

The pilot fits both PCM and RSM to the same known-truth PCM rows in FACETS
4.5 and the application JMLE.  Engine agreement, model misspecification, and
structural nonidentification are kept as separate result layers.
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
from scipy.stats import spearmanr


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from validation.facets_pcm_known_truth_smoke import (  # noqa: E402
    AUXILIARY_DECIMALS,
    CATEGORIES,
    CRITERION_TRUTH,
    N_PERSONS,
    RATER_TRUTH,
    TASK_TRUTH,
    THRESHOLD_CONDITIONS,
    _frame_sha256,
    _observation_weights,
    _weighted_mean,
    apply_threshold_condition,
    generate_latent_replicate,
)
from validation.facets_pcm_probe import (  # noqa: E402
    build_pcm_spec,
    map_pcm_scale_tables,
    normalize_python_pcm_steps,
    pcm_category_count_audit,
    pcm_threshold_pairs,
)
from validation.operating_characteristics_facets import (  # noqa: E402
    FACET_COLUMNS,
    RECOVERY_FACETS,
    align_recovery,
    build_facets_spec,
    invoke_facets,
    parse_iteration_report,
    parse_score_file,
    parse_table8_categories,
    sha256_file,
    slugify,
    validate_bundle,
)
from validation.operating_characteristics_python_table8 import (  # noqa: E402
    build_jmle_kwargs,
)


SCHEMA_VERSION = "mfrm-facets-pcm-boundary-pilot-v1"
SEEDS = (842201, 842202)
DESIGNS = ("complete", "planned_connected", "disconnected_negative_control")
FIT_MODELS = ("PCM", "RSM")
MAIN_MAE_MAX = 0.01
MAIN_ABS_MAX = 0.05
MIN_SPEARMAN = 0.999
THRESHOLD_MAE_MAX = 0.02
THRESHOLD_ABS_MAX = 0.05


def apply_observation_design(frame: pd.DataFrame, design: str) -> pd.DataFrame:
    """Apply a deterministic complete, ring-connected, or disconnected mask."""

    if design == "complete":
        output = frame.copy()
    elif design == "planned_connected":
        rater_levels = list(RATER_TRUTH)
        person_index = frame["Person"].astype(str).str.extract(r"(\d+)$", expand=False).astype(int) - 1
        rater_index = frame["Rater"].map({level: index for index, level in enumerate(rater_levels)})
        keep = rater_index.eq(person_index.mod(len(rater_levels))) | rater_index.eq(
            person_index.add(1).mod(len(rater_levels))
        )
        output = frame.loc[keep].copy()
    elif design == "disconnected_negative_control":
        person_index = frame["Person"].astype(str).str.extract(r"(\d+)$", expand=False).astype(int)
        first = person_index.le(N_PERSONS // 2) & frame["Rater"].isin(["R01", "R02"])
        second = person_index.gt(N_PERSONS // 2) & frame["Rater"].isin(["R03", "R04"])
        output = frame.loc[first | second].copy()
    else:
        raise ValueError(f"Unknown observation design: {design}")
    return output.reset_index(drop=True)


def person_rater_components(ratings: pd.DataFrame) -> list[dict[str, Any]]:
    """Return connected components in the bipartite Person-Rater graph."""

    adjacency: dict[str, set[str]] = {}
    for person, rater in ratings[["Person", "Rater"]].astype(str).itertuples(index=False):
        p_node = f"Person::{person}"
        r_node = f"Rater::{rater}"
        adjacency.setdefault(p_node, set()).add(r_node)
        adjacency.setdefault(r_node, set()).add(p_node)
    seen: set[str] = set()
    components = []
    for start in sorted(adjacency):
        if start in seen:
            continue
        stack = [start]
        nodes: set[str] = set()
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            nodes.add(node)
            stack.extend(adjacency[node] - seen)
        components.append({
            "Persons": sum(node.startswith("Person::") for node in nodes),
            "Raters": sum(node.startswith("Rater::") for node in nodes),
            "Nodes": sorted(nodes),
        })
    return components


def constrained_adjacent_design_audit(ratings: pd.DataFrame, model: str) -> dict[str, Any]:
    """Audit the full constrained adjacent-category design for RSM or PCM."""

    model = str(model).upper()
    if model not in FIT_MODELS:
        raise ValueError(f"Unsupported audit model: {model}")
    levels = {
        facet: list(dict.fromkeys(ratings[facet].astype(str)))
        for facet in ("Person", "Rater", "Task", "Criterion")
    }
    n_steps = CATEGORIES - 1
    columns = [f"Person::{level}" for level in levels["Person"]]
    slices: dict[str, slice] = {"Person": slice(0, len(columns))}
    for facet in ("Rater", "Task", "Criterion"):
        start = len(columns)
        columns.extend(f"{facet}::{level}" for level in levels[facet][:-1])
        slices[facet] = slice(start, len(columns))
    threshold_slices: dict[str, slice] = {}
    scale_levels = levels["Criterion"] if model == "PCM" else ["__COMMON__"]
    for scale in scale_levels:
        start = len(columns)
        columns.extend(f"Threshold::{scale}::{step}" for step in range(1, n_steps))
        threshold_slices[scale] = slice(start, len(columns))

    def centered_code(level: str, ordered: list[str]) -> np.ndarray:
        code = np.zeros(len(ordered) - 1, dtype=float)
        index = ordered.index(level)
        if index < len(ordered) - 1:
            code[index] = 1.0
        else:
            code[:] = -1.0
        return code

    matrix = np.zeros((len(ratings) * n_steps, len(columns)), dtype=float)
    out_row = 0
    for row in ratings.itertuples(index=False):
        person_index = levels["Person"].index(str(row.Person))
        scale = str(row.Criterion) if model == "PCM" else "__COMMON__"
        for transition in range(1, n_steps + 1):
            matrix[out_row, person_index] = 1.0
            for facet in ("Rater", "Task", "Criterion"):
                matrix[out_row, slices[facet]] = -centered_code(
                    str(getattr(row, facet)), levels[facet]
                )
            threshold_slice = threshold_slices[scale]
            if transition < n_steps:
                matrix[out_row, threshold_slice.start + transition - 1] = -1.0
            else:
                matrix[out_row, threshold_slice] = 1.0
            out_row += 1
    rank = int(np.linalg.matrix_rank(matrix))
    components = person_rater_components(ratings)
    return {
        "Model": model,
        "Rows": int(matrix.shape[0]),
        "Columns": int(matrix.shape[1]),
        "Rank": rank,
        "Nullity": int(matrix.shape[1] - rank),
        "PersonRaterComponents": int(len(components)),
        "ComponentSizes": [
            {"Persons": component["Persons"], "Raters": component["Raters"]}
            for component in components
        ],
        "ConstraintContract": (
            "Person noncentered; Rater/Task/Criterion sum-to-zero; "
            + ("each Criterion threshold vector" if model == "PCM" else "one common threshold vector")
            + " sum-to-zero"
        ),
    }


def generate_boundary_bundle() -> dict[str, pd.DataFrame]:
    manifests = []
    ratings_parts = []
    truth_parts = []
    threshold_parts = []
    for replicate, seed in enumerate(SEEDS, start=1):
        latent, person_truth = generate_latent_replicate(seed)
        latent_hash = _frame_sha256(
            latent,
            ["Person", "Rater", "Task", "Criterion", "Theta", "Eta", "Uniform"],
        )
        for condition, thresholds in THRESHOLD_CONDITIONS.items():
            complete = apply_threshold_condition(latent, thresholds)
            for design in DESIGNS:
                run_id = f"{design}__{condition}::rep-{replicate:05d}"
                run_ratings = apply_observation_design(complete, design)
                expected_nullity = 1 if design == "disconnected_negative_control" else 0
                manifests.append({
                    "RunId": run_id,
                    "ConditionId": f"{design}__{condition}",
                    "Design": design,
                    "TruthBias": 0.0,
                    "Replicate": replicate,
                    "Seed": seed,
                    "Categories": CATEGORIES,
                    "ThresholdCondition": condition,
                    "ExpectedStructuralNullity": expected_nullity,
                    "PairedLatentSHA256": latent_hash,
                    "ExpectedRows": int(len(run_ratings)),
                })
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
                        threshold_parts.append({
                            "RunId": run_id,
                            "ThresholdCondition": condition,
                            "StepFacetLevel": criterion,
                            "Category": category,
                            "ThresholdTruth": float(value),
                        })
    return {
        "manifest.csv": pd.DataFrame(manifests),
        "generated_ratings.csv": pd.concat(ratings_parts, ignore_index=True),
        "generated_facet_truth.csv": pd.concat(truth_parts, ignore_index=True),
        "generated_anchors.csv": pd.DataFrame(columns=["RunId", "Facet", "Level", "Anchor"]),
        "generated_pcm_threshold_truth.csv": pd.DataFrame(threshold_parts),
    }


def write_bundle(bundle: dict[str, pd.DataFrame], input_dir: Path) -> None:
    input_dir.mkdir(parents=True, exist_ok=False)
    for filename, frame in bundle.items():
        frame.to_csv(input_dir / filename, index=False, lineterminator="\n")


def build_rsm_spec(
    manifest_row: pd.Series,
    ratings: pd.DataFrame,
    *,
    score_base: Path,
) -> tuple[str, dict[str, dict[str, int]]]:
    """Build a single-model RSM spec without the historical bias statement."""

    base, maps = build_facets_spec(
        manifest_row, ratings, pd.DataFrame(), score_base=score_base
    )
    source = "Models=\n?,?,?,?,R3\n?,?B,?B,?,R3\n*"
    replacement = "Models=\n?,?,?,?,R3\n*"
    if base.count(source) != 1:
        raise ValueError("Registered RSM Models block changed; conversion fails closed")
    output = base.replace(source, replacement)
    output = output.replace("Title=FACETS replay ", "Title=FACETS RSM boundary pilot ", 1)
    return output, maps


def python_fit(ratings: pd.DataFrame, model: str) -> dict[str, Any]:
    import streamlit_app as app  # pylint: disable=import-outside-toplevel

    kwargs = build_jmle_kwargs(CATEGORIES)
    kwargs.update({"model": model, "maxit": 400})
    if model == "PCM":
        kwargs["step_facet"] = "Criterion"
    return app.mfrm_estimate(
        data=ratings[["Person", "Rater", "Task", "Criterion", "Score"]].copy(),
        **kwargs,
    )


def rsm_count_audit(categories: pd.DataFrame, ratings: pd.DataFrame) -> pd.DataFrame:
    if categories["TableNumber"].nunique() != 1:
        raise ValueError("RSM Table 8 must contain exactly one modeled scale")
    expected = ratings.groupby("Score", observed=False).size().rename("GeneratedCount").reset_index()
    expected = expected.rename(columns={"Score": "Category"})
    audit = categories.merge(
        expected, on="Category", how="outer", validate="one_to_one", indicator=True
    )
    audit["CountMatches"] = (
        audit["_merge"].eq("both")
        & pd.to_numeric(audit["TotalCount"], errors="coerce").eq(
            pd.to_numeric(audit["GeneratedCount"], errors="coerce")
        )
    )
    return audit


def normalize_rsm_steps(result: dict[str, Any]) -> pd.DataFrame:
    source = result.get("steps", pd.DataFrame()).copy()
    if source.empty or not {"Step", "Estimate"}.issubset(source.columns):
        raise ValueError("Python RSM steps are unavailable")
    output = pd.DataFrame({
        "Category": source["Step"].astype(str).str.extract(r"(\d+)$", expand=False).astype(int),
        "Step": source["Step"].astype(str),
        "PythonThreshold": pd.to_numeric(source["Estimate"], errors="coerce"),
    })
    if not np.isfinite(output["PythonThreshold"]).all():
        raise ValueError("Python RSM thresholds contain non-finite values")
    return output


def rsm_threshold_pairs(categories: pd.DataFrame, python_steps: pd.DataFrame) -> pd.DataFrame:
    facets = categories.loc[
        pd.to_numeric(categories["Category"], errors="coerce").gt(0),
        [
            "TableNumber", "Model", "Category", "ThresholdMeasureToken",
            "ThresholdMeasureDisplayed", "ThresholdMeasureDisplayDecimals",
            "ThresholdMeasureRawLowerBound", "ThresholdMeasureRawUpperBound",
        ],
    ].copy()
    pairs = facets.merge(
        python_steps, on="Category", how="outer", validate="one_to_one", indicator=True
    )
    if not pairs["_merge"].eq("both").all():
        raise ValueError("FACETS and Python RSM threshold keys do not match")
    pairs = pairs.drop(columns="_merge")
    pairs["StepFacetLevel"] = "__COMMON__"
    pairs["ThresholdDifference"] = (
        pairs["ThresholdMeasureDisplayed"] - pairs["PythonThreshold"]
    )
    pairs["AbsoluteThresholdDifference"] = pairs["ThresholdDifference"].abs()
    return pairs


def _fit_facets_and_python(
    manifest_row: pd.Series,
    ratings: pd.DataFrame,
    truth: pd.DataFrame,
    threshold_truth: pd.DataFrame,
    *,
    fit_model: str,
    facets_exe: Path,
    work_root: Path,
    timeout_seconds: float,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    run_id = str(manifest_row["RunId"])
    fit_id = f"{run_id}::{fit_model}"
    run_dir = work_root / f"{slugify(run_id)}__{fit_model.lower()}"
    run_dir.mkdir(parents=True, exist_ok=False)
    design = constrained_adjacent_design_audit(ratings, fit_model)

    started = time.perf_counter()
    python_result = python_fit(ratings, fit_model)
    python_summary = python_result["summary"].iloc[0]
    python_estimates = python_result["facets"]["others"].copy()
    python_estimates["SE"] = np.nan
    python_estimates["Status"] = np.nan
    python_recovery = align_recovery(
        python_estimates, truth, pd.DataFrame(), manifest_row
    )

    spec_path = run_dir / "analysis.txt"
    report_path = run_dir / "report_u6.txt"
    aux_report_path = run_dir / "report_u2.txt"
    score_base = run_dir / "scores.txt"
    aux_score_base = run_dir / "scores_u2.txt"
    if fit_model == "PCM":
        spec, level_maps = build_pcm_spec(
            manifest_row, ratings, pd.DataFrame(), score_base=score_base
        )
    else:
        spec, level_maps = build_rsm_spec(
            manifest_row, ratings, score_base=score_base
        )
    spec_path.write_text(spec, encoding="utf-8", newline="\n")
    primary = invoke_facets(
        facets_exe, spec_path, report_path, timeout_seconds=timeout_seconds
    )
    (run_dir / "stdout_u6.txt").write_text(primary.stdout or "", encoding="utf-8")
    (run_dir / "stderr_u6.txt").write_text(primary.stderr or "", encoding="utf-8")
    if primary.returncode != 0 or not report_path.is_file():
        raise RuntimeError(f"FACETS primary failed for {fit_id}: exit={primary.returncode}")
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
        raise RuntimeError(f"FACETS auxiliary failed for {fit_id}: exit={auxiliary.returncode}")

    report_info = parse_iteration_report(report_path)
    scores = pd.concat([
        parse_score_file(
            run_dir / f"scores.{number}.txt", facet_number=number, facet_name=name
        )
        for number, name in FACET_COLUMNS
    ], ignore_index=True)
    facets_recovery = align_recovery(scores, truth, pd.DataFrame(), manifest_row)
    primary_categories = parse_table8_categories(report_path)
    auxiliary_categories = parse_table8_categories(aux_report_path)
    if fit_model == "PCM":
        primary_categories = map_pcm_scale_tables(
            primary_categories, criterion_map=level_maps["Criterion"]
        )
        auxiliary_categories = map_pcm_scale_tables(
            auxiliary_categories, criterion_map=level_maps["Criterion"]
        )
        count_audit = pcm_category_count_audit(auxiliary_categories, ratings)
        threshold_pairs = pcm_threshold_pairs(
            primary_categories, normalize_python_pcm_steps(python_result)
        )
        threshold_pairs = threshold_pairs.merge(
            threshold_truth[["StepFacetLevel", "Category", "ThresholdTruth"]],
            on=["StepFacetLevel", "Category"],
            validate="one_to_one",
        )
        threshold_pairs["TruthTargetStatus"] = "data_generating_pcm_threshold"
        criterion_weights = ratings.groupby("Criterion", observed=False).size().to_dict()
        threshold_pairs["ObservationWeight"] = threshold_pairs["StepFacetLevel"].map(
            criterion_weights
        )
    else:
        primary_categories = primary_categories.assign(StepFacetLevel="__COMMON__")
        auxiliary_categories = auxiliary_categories.assign(StepFacetLevel="__COMMON__")
        count_audit = rsm_count_audit(auxiliary_categories, ratings)
        threshold_pairs = rsm_threshold_pairs(
            primary_categories, normalize_rsm_steps(python_result)
        )
        if str(manifest_row["ThresholdCondition"]) == "shared":
            common_truth = (
                threshold_truth.groupby("Category", as_index=False)["ThresholdTruth"].mean()
            )
            threshold_pairs = threshold_pairs.merge(
                common_truth, on="Category", validate="one_to_one"
            )
            threshold_pairs["TruthTargetStatus"] = "data_generating_common_threshold"
        else:
            threshold_pairs["ThresholdTruth"] = np.nan
            threshold_pairs["TruthTargetStatus"] = "undefined_under_pcm_to_rsm_misspecification"
        threshold_pairs["ObservationWeight"] = len(ratings)

    keys = ["RunId", "Facet", "Level"]
    recovery_pairs = facets_recovery.merge(
        python_recovery[keys + ["EstimateAligned", "ErrorAligned"]].rename(columns={
            "EstimateAligned": "PythonEstimateAligned",
            "ErrorAligned": "PythonTruthErrorAligned",
        }),
        on=keys,
        validate="one_to_one",
    ).rename(columns={
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
    common_eligible = bool(
        design["Nullity"] == 0
        and report_info.get("Converged", False)
        and python_summary.get("Converged", False)
        and python_summary.get("InferenceReady", False)
        and count_audit["CountMatches"].all()
    )
    recovery_pairs["ComparisonEligible"] = common_eligible
    threshold_pairs["ComparisonEligible"] = common_eligible
    threshold_pairs["FACETSTruthError"] = (
        threshold_pairs["ThresholdMeasureDisplayed"] - threshold_pairs["ThresholdTruth"]
    )
    threshold_pairs["PythonTruthError"] = (
        threshold_pairs["PythonThreshold"] - threshold_pairs["ThresholdTruth"]
    )
    correlations = {}
    for facet, group in recovery_pairs.groupby("Facet", sort=False):
        statistic = spearmanr(
            group["FACETSEstimateAligned"], group["PythonEstimateAligned"]
        ).statistic
        correlations[str(facet)] = float(statistic) if np.isfinite(statistic) else np.nan
    direct_main = recovery_pairs[recovery_pairs["ComparisonEligible"]]
    direct_threshold = threshold_pairs[threshold_pairs["ComparisonEligible"]]
    main_mae = _weighted_mean(
        direct_main["AbsoluteEstimateDifference"], direct_main["ObservationWeight"]
    ) if len(direct_main) else np.nan
    main_max = float(direct_main["AbsoluteEstimateDifference"].max()) if len(direct_main) else np.nan
    threshold_mae = _weighted_mean(
        direct_threshold["AbsoluteThresholdDifference"], direct_threshold["ObservationWeight"]
    ) if len(direct_threshold) else np.nan
    threshold_max = float(direct_threshold["AbsoluteThresholdDifference"].max()) if len(direct_threshold) else np.nan
    min_spearman = float(min(correlations.values())) if correlations else np.nan
    metrics = {
        "FitId": fit_id,
        "RunId": run_id,
        "FitModel": fit_model,
        "Design": manifest_row["Design"],
        "ThresholdCondition": manifest_row["ThresholdCondition"],
        "Replicate": int(manifest_row["Replicate"]),
        "Seed": int(manifest_row["Seed"]),
        "Rows": int(len(ratings)),
        "IndependentNullity": int(design["Nullity"]),
        "PersonRaterComponents": int(design["PersonRaterComponents"]),
        "FACETSVersion": report_info.get("FacetsVersion"),
        "FACETSConverged": bool(report_info.get("Converged", False)),
        "FACETSSubsetConnected": bool(report_info.get("SubsetConnected", False)),
        "PythonConverged": bool(python_summary.get("Converged", False)),
        "PythonInferenceReady": bool(python_summary.get("InferenceReady", False)),
        "PythonEtaStructuralNullity": python_summary.get("EtaStructuralNullity"),
        "ComparisonEligible": common_eligible,
        "Table8ScaleTables": int(primary_categories["TableNumber"].nunique()),
        "CategoryCountRows": int(len(count_audit)),
        "CategoryCountMatches": int(count_audit["CountMatches"].sum()),
        "MainWeightedMAE": main_mae,
        "MainMaxAbsDifference": main_max,
        "MinimumWithinFacetSpearman": min_spearman,
        "ThresholdWeightedMAE": threshold_mae,
        "ThresholdMaxAbsDifference": threshold_max,
        "PythonLogLik": float(python_summary["LogLik"]),
        "PythonLogLikPerObs": float(python_summary["LogLikPerObs"]),
        "PythonAIC": float(python_summary["AIC"]),
        "PythonBIC": float(python_summary["BIC"]),
        "ElapsedSeconds": time.perf_counter() - started,
        "FailureReason": "",
    }
    metrics["DirectAgreementPass"] = bool(
        common_eligible
        and main_mae <= MAIN_MAE_MAX
        and main_max <= MAIN_ABS_MAX
        and min_spearman >= MIN_SPEARMAN
        and threshold_mae <= THRESHOLD_MAE_MAX
        and threshold_max <= THRESHOLD_ABS_MAX
    )
    for frame in (recovery_pairs, threshold_pairs, count_audit, primary_categories):
        frame.insert(0, "FitModel", fit_model)
        frame.insert(0, "FitId", fit_id)
    recovery_pairs.to_csv(run_dir / "main_pairs.csv", index=False)
    threshold_pairs.to_csv(run_dir / "threshold_pairs.csv", index=False)
    count_audit.to_csv(run_dir / "category_count_audit.csv", index=False)
    primary_categories.to_csv(run_dir / "table8_u6.csv", index=False)
    (run_dir / "design_audit.json").write_text(
        json.dumps(design, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return metrics, recovery_pairs, threshold_pairs, count_audit, primary_categories


def build_model_selection(fits: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "RunId", "Design", "ThresholdCondition", "Replicate", "Seed",
        "IndependentNullity", "ComparisonEligible", "FitModel", "PythonLogLik",
        "PythonLogLikPerObs", "PythonAIC", "PythonBIC",
    ]
    source = fits[columns].copy()
    pcm = source[source["FitModel"].eq("PCM")].drop(columns="FitModel")
    rsm = source[source["FitModel"].eq("RSM")].drop(columns="FitModel")
    pairs = pcm.merge(
        rsm,
        on=["RunId", "Design", "ThresholdCondition", "Replicate", "Seed"],
        suffixes=("_PCM", "_RSM"),
        validate="one_to_one",
    )
    pairs["ModelComparisonEligible"] = (
        pairs["ComparisonEligible_PCM"]
        & pairs["ComparisonEligible_RSM"]
        & pairs["IndependentNullity_PCM"].eq(0)
        & pairs["IndependentNullity_RSM"].eq(0)
    )
    pairs["PCMMinusRSMLogLik"] = pairs["PythonLogLik_PCM"] - pairs["PythonLogLik_RSM"]
    pairs["PCMMinusRSMLogLikPerObs"] = (
        pairs["PythonLogLikPerObs_PCM"] - pairs["PythonLogLikPerObs_RSM"]
    )
    pairs["RSMMinusPCMAIC"] = pairs["PythonAIC_RSM"] - pairs["PythonAIC_PCM"]
    pairs["RSMMinusPCMBIC"] = pairs["PythonBIC_RSM"] - pairs["PythonBIC_PCM"]
    return pairs


def build_direction_checks(model_selection: pd.DataFrame) -> pd.DataFrame:
    eligible = model_selection[model_selection["ModelComparisonEligible"]]
    rows = []
    for (design, replicate), group in eligible.groupby(["Design", "Replicate"]):
        lookup = group.set_index("ThresholdCondition")
        shared = lookup.loc["shared", "PCMMinusRSMLogLikPerObs"] if "shared" in lookup.index else np.nan
        heterogeneous = (
            lookup.loc["heterogeneous", "PCMMinusRSMLogLikPerObs"]
            if "heterogeneous" in lookup.index else np.nan
        )
        rows.append({
            "Design": design,
            "Replicate": int(replicate),
            "SharedPCMMinusRSMLogLikPerObs": shared,
            "HeterogeneousPCMMinusRSMLogLikPerObs": heterogeneous,
            "DirectionPass": bool(
                np.isfinite(shared) and np.isfinite(heterogeneous) and heterogeneous > shared
            ),
        })
    return pd.DataFrame(rows)


def run_pilot(args: argparse.Namespace) -> None:
    output_dir = args.output_dir.resolve()
    facets_exe = args.facets_exe.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    if not facets_exe.is_file():
        raise FileNotFoundError(f"FACETS executable is missing: {facets_exe}")
    output_dir.mkdir(parents=True)
    input_dir = output_dir / "retained_input"
    write_bundle(generate_boundary_bundle(), input_dir)
    tables = validate_bundle(input_dir)
    manifest = tables["manifest.csv"]
    ratings_all = tables["generated_ratings.csv"]
    truth_all = tables["generated_facet_truth.csv"]
    threshold_truth_all = pd.read_csv(input_dir / "generated_pcm_threshold_truth.csv")
    work_root = output_dir / "fit_runs"
    work_root.mkdir()
    fit_rows = []
    main_parts = []
    threshold_parts = []
    count_parts = []
    table8_parts = []
    for _, manifest_row in manifest.iterrows():
        run_id = str(manifest_row["RunId"])
        ratings = ratings_all[ratings_all["RunId"].astype(str).eq(run_id)].drop(columns="RunId")
        truth = truth_all[truth_all["RunId"].astype(str).eq(run_id)]
        threshold_truth = threshold_truth_all[
            threshold_truth_all["RunId"].astype(str).eq(run_id)
        ]
        for fit_model in FIT_MODELS:
            design = constrained_adjacent_design_audit(ratings, fit_model)
            try:
                metrics, main, thresholds, counts, table8 = _fit_facets_and_python(
                    manifest_row,
                    ratings,
                    truth,
                    threshold_truth,
                    fit_model=fit_model,
                    facets_exe=facets_exe,
                    work_root=work_root,
                    timeout_seconds=args.timeout_seconds,
                )
                fit_rows.append(metrics)
                main_parts.append(main)
                threshold_parts.append(thresholds)
                count_parts.append(counts)
                table8_parts.append(table8)
            except Exception as exc:  # retain all attempted fit contracts
                fit_rows.append({
                    "FitId": f"{run_id}::{fit_model}",
                    "RunId": run_id,
                    "FitModel": fit_model,
                    "Design": manifest_row["Design"],
                    "ThresholdCondition": manifest_row["ThresholdCondition"],
                    "Replicate": int(manifest_row["Replicate"]),
                    "Seed": int(manifest_row["Seed"]),
                    "Rows": int(len(ratings)),
                    "IndependentNullity": int(design["Nullity"]),
                    "PersonRaterComponents": int(design["PersonRaterComponents"]),
                    "FACETSConverged": False,
                    "FACETSSubsetConnected": False,
                    "PythonConverged": False,
                    "PythonInferenceReady": False,
                    "ComparisonEligible": False,
                    "CategoryCountRows": 0,
                    "CategoryCountMatches": 0,
                    "DirectAgreementPass": False,
                    "FailureReason": f"{type(exc).__name__}: {exc}",
                })
    fits = pd.DataFrame(fit_rows)
    main_pairs = pd.concat(main_parts, ignore_index=True) if main_parts else pd.DataFrame()
    threshold_pairs = pd.concat(threshold_parts, ignore_index=True) if threshold_parts else pd.DataFrame()
    counts = pd.concat(count_parts, ignore_index=True) if count_parts else pd.DataFrame()
    table8 = pd.concat(table8_parts, ignore_index=True) if table8_parts else pd.DataFrame()
    model_selection = build_model_selection(fits)
    directions = build_direction_checks(model_selection)

    rank_full = fits[fits["Design"].isin(["complete", "planned_connected"])]
    disconnected = fits[fits["Design"].eq("disconnected_negative_control")]
    expected_count_rows = int(sum(
        (len(CRITERION_TRUTH) if model == "PCM" else 1) * CATEGORIES
        for _ in range(len(manifest)) for model in FIT_MODELS
    ))
    count_gate = bool(
        len(counts) == expected_count_rows
        and not counts.empty
        and counts["CountMatches"].all()
    )
    rank_gate = bool(rank_full["IndependentNullity"].eq(0).all())
    negative_gate = bool(
        disconnected["IndependentNullity"].ge(1).all()
        and ~disconnected["PythonInferenceReady"].fillna(False).any()
        and ~disconnected["ComparisonEligible"].fillna(False).any()
    )
    engine_gate = bool(
        rank_full["ComparisonEligible"].all()
        and rank_full["DirectAgreementPass"].all()
    )
    direction_gate = bool(len(directions) == 4 and directions["DirectionPass"].all())
    aggregate = {
        "schema_version": SCHEMA_VERSION,
        "claim_limit": "Two-replicate boundary pilot; no estimator-performance claim.",
        "generated_runs": int(len(manifest)),
        "attempted_fit_contracts": int(len(fits)),
        "successful_fit_contracts": int(fits["FailureReason"].fillna("").eq("").sum()),
        "rank_full_fit_contracts": int(len(rank_full)),
        "rank_full_comparison_eligible": int(rank_full["ComparisonEligible"].sum()),
        "rank_full_direct_agreement_passes": int(rank_full["DirectAgreementPass"].sum()),
        "disconnected_fit_contracts": int(len(disconnected)),
        "disconnected_with_independent_nullity": int(disconnected["IndependentNullity"].ge(1).sum()),
        "disconnected_python_inference_withheld": int((~disconnected["PythonInferenceReady"].fillna(False)).sum()),
        "facets_subset_connected_in_disconnected": int(disconnected["FACETSSubsetConnected"].fillna(False).sum()),
        "table8_count_rows": int(len(counts)),
        "expected_table8_count_rows": expected_count_rows,
        "table8_count_matches": int(counts["CountMatches"].sum()) if len(counts) else 0,
        "direction_checks": int(len(directions)),
        "direction_passes": int(directions["DirectionPass"].sum()) if len(directions) else 0,
        "gates": {
            "input_counts": count_gate,
            "rank_full_designs": rank_gate,
            "disconnected_negative_control": negative_gate,
            "rank_full_engine_agreement": engine_gate,
            "misspecification_direction": direction_gate,
        },
    }
    aggregate["qualification_pass"] = bool(all(aggregate["gates"].values()))
    if not main_pairs.empty:
        eligible_main = main_pairs[main_pairs["ComparisonEligible"]]
        aggregate["eligible_main_parameters"] = int(len(eligible_main))
        aggregate["main_weighted_mae_logits"] = _weighted_mean(
            eligible_main["AbsoluteEstimateDifference"], eligible_main["ObservationWeight"]
        )
        aggregate["main_max_abs_difference_logits"] = float(
            eligible_main["AbsoluteEstimateDifference"].max()
        )
    if not threshold_pairs.empty:
        eligible_threshold = threshold_pairs[threshold_pairs["ComparisonEligible"]]
        aggregate["eligible_thresholds"] = int(len(eligible_threshold))
        aggregate["threshold_weighted_mae_logits"] = _weighted_mean(
            eligible_threshold["AbsoluteThresholdDifference"],
            eligible_threshold["ObservationWeight"],
        )
        aggregate["threshold_max_abs_difference_logits"] = float(
            eligible_threshold["AbsoluteThresholdDifference"].max()
        )

    fits.to_csv(output_dir / "boundary_fit_ledger.csv", index=False)
    main_pairs.to_csv(output_dir / "boundary_main_pairs.csv", index=False)
    threshold_pairs.to_csv(output_dir / "boundary_threshold_pairs.csv", index=False)
    counts.to_csv(output_dir / "boundary_count_audit.csv", index=False)
    table8.to_csv(output_dir / "boundary_table8_u6.csv", index=False)
    model_selection.to_csv(output_dir / "boundary_model_selection.csv", index=False)
    directions.to_csv(output_dir / "boundary_direction_checks.csv", index=False)
    (output_dir / "boundary_metrics.json").write_text(
        json.dumps(aggregate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    identity = {
        "schema_version": SCHEMA_VERSION,
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "plan_sha256": sha256_file(REPO_ROOT / "validation" / "facets_pcm_boundary_pilot_plan_20260811.json"),
        "app_sha256": sha256_file(REPO_ROOT / "streamlit_app.py"),
        "facets_executable": str(facets_exe),
        "facets_executable_sha256": sha256_file(facets_exe),
        "python_version": sys.version,
        "platform": platform.platform(),
        "retained_input_sha256": {
            path.name: sha256_file(path) for path in sorted(input_dir.glob("*.csv"))
        },
    }
    (output_dir / "boundary_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--facets-exe", type=Path, default=Path(r"C:\Facets\Facets.exe"))
    parser.add_argument("--timeout-seconds", type=float, default=120.0)
    return parser.parse_args(argv)


def main() -> None:
    run_pilot(parse_args())


if __name__ == "__main__":
    main()
