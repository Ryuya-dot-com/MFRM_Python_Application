#!/usr/bin/env python3
"""Run a reproducible Python operating-characteristics pilot.

This runner deliberately separates three profiles:

``smoke``
    Two replicates per condition.  Verifies orchestration and failure
    accounting only; never use the rates as performance evidence.
``pilot``
    Twenty replicates per condition.  Useful for debugging the design matrix
    and estimating runtime, still below the default study-depth floor.
``study``
    One hundred replicates per condition.  A starting point for Monte Carlo
    screening, not a universal validation threshold.

The current increment fits the standalone Python JMLE engine.  The emitted
manifest and result tables are engine-neutral so mfrmr, TAM, immer, and sirt
can be appended without changing the aggregation contract.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
import platform
from pathlib import Path
import sys
import time
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import streamlit_app as app  # noqa: E402
from mfrm_app import decision_stability as ds  # noqa: E402
from mfrm_app import operating_characteristics as oc  # noqa: E402


PROFILE_REPLICATES = {"smoke": 2, "pilot": 20, "study": 100}
DEFAULT_OUTPUT = REPO / "validation" / "operating_characteristics_20260809"
DEFAULT_PRECISION_PLAN = REPO / "validation" / "operating_characteristics_precision_plan_20260809.json"
FOCAL_RATER = "R01"
FOCAL_TASK = "T01"
GENERATED_BUNDLE_FILES = {
    "generated_ratings.csv": [
        "SchemaVersion", "RunId", "ConditionId", "Design", "TruthBias",
        "Replicate", "Seed", "Person", "Rater", "Task", "Criterion", "Score",
    ],
    "generated_facet_truth.csv": [
        "SchemaVersion", "RunId", "ConditionId", "Design", "TruthBias",
        "Replicate", "Seed", "Facet", "Level", "Truth",
    ],
    "generated_anchors.csv": [
        "SchemaVersion", "RunId", "ConditionId", "Design", "TruthBias",
        "Replicate", "Seed", "Facet", "Level", "Anchor",
    ],
}
DERIVED_CROSS_ENGINE_OUTPUT_FILES = (
    "PILOT20_ASSESSMENT.md",
    "pilot_integrity_checks.csv",
    "pilot_budget_assessment.csv",
    "pilot_gradient_readiness.csv",
    "pilot_anchor_sensitivity.csv",
    "pilot_first_read_summary.csv",
    "pilot_eligibility_funnel.png",
    "pilot_gradient_readiness.png",
    "pilot_float_boundary_mismatches.png",
    "pilot_anchor_contamination.png",
    "BRIDGE_RESULTS.md",
    "bridge_engine_availability_r.csv",
    "bridge_import_summary_r.csv",
    "bridge_validation_r.csv",
    "MFRMR_RESULTS.md",
    "mfrmr_adapter_identity.csv",
    "mfrmr_bias_decisions.csv",
    "mfrmr_parameter_recovery.csv",
    "mfrmr_runs.csv",
    "cross_engine_runs_oc.csv",
    "cross_engine_failure_accounting_oc.csv",
    "cross_engine_failure_reasons_oc.csv",
    "cross_engine_bias_decisions_oc.csv",
    "cross_engine_binary_oc.csv",
    "cross_engine_parameter_recovery_oc.csv",
    "cross_engine_estimation_oc.csv",
    "cross_engine_first_read_summary.csv",
    "python_mfrmr_parameter_pairs.csv",
    "python_mfrmr_parameter_agreement.csv",
    "python_mfrmr_bias_pairs.csv",
    "python_mfrmr_bias_agreement.csv",
    "python_mfrmr_readiness_contrast.csv",
    "python_mfrmr_parameter_agreement.png",
    "python_mfrmr_bias_agreement.png",
    "python_mfrmr_readiness.png",
    "python_cmle_runs.csv",
    "python_cmle_coefficients.csv",
    "python_cmle_parameter_recovery.csv",
    "python_cmle_adapter_identity.csv",
    "immer_cmle_runs.csv",
    "immer_cmle_coefficients.csv",
    "immer_cmle_adapter_identity.csv",
    "python_immer_cmle_pairs.csv",
    "python_immer_cmle_agreement.csv",
    "python_immer_cmle_readiness.csv",
    "python_cmle_omitted_bias_sensitivity.csv",
    "python_cmle_omitted_bias_summary.csv",
    "python_immer_cmle_parameter_agreement.png",
    "python_immer_cmle_readiness.png",
    "CMLE_IMMER_RESULTS.md",
    "sirt_runs.csv",
    "sirt_rater_recovery.csv",
    "sirt_item_recovery.csv",
    "sirt_person_estimates.csv",
    "sirt_adapter_identity.csv",
    "python_sirt_rater_pairs.csv",
    "python_sirt_rater_agreement.csv",
    "sirt_quadrature_pairs.csv",
    "sirt_quadrature_summary.csv",
    "sirt_anchor_sensitivity.csv",
    "sirt_omitted_bias_sensitivity.csv",
    "sirt_sparse_identification.csv",
    "sirt_first_read_summary.csv",
    "python_sirt_rater_agreement.png",
    "sirt_quadrature_sensitivity.png",
    "SIRT_RESULTS.md",
    "tam_runs.csv",
    "tam_facet_recovery.csv",
    "tam_step_estimates.csv",
    "tam_person_estimates.csv",
    "tam_surface_estimates.csv",
    "tam_adapter_identity.csv",
    "tam_quadrature_pairs.csv",
    "tam_quadrature_summary.csv",
    "tam_quadrature_run_pairs.csv",
    "python_tam_facet_pairs.csv",
    "python_tam_facet_agreement.csv",
    "tam_sirt_rater_pairs.csv",
    "tam_sirt_rater_agreement.csv",
    "tam_anchor_sensitivity.csv",
    "tam_omitted_bias_sensitivity.csv",
    "tam_variance_boundary_audit.csv",
    "tam_progress_threshold_summary.csv",
    "tam_sparse_identification.csv",
    "tam_first_read_summary.csv",
    "python_tam_facet_agreement.png",
    "tam_quadrature_sensitivity.png",
    "tam_anchor_constraint_transmission.png",
    "TAM_RESULTS.md",
)


@dataclass(frozen=True)
class GeneratedDesign:
    data: pd.DataFrame
    facet_truth: pd.DataFrame
    anchors: pd.DataFrame
    realized_missing_rate: float
    category_counts: pd.DataFrame


def sha256_file(path: Path) -> str:
    """Hash retained bridge bytes so R can validate the exact same inputs."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _identified_generated_frame(
    frame: pd.DataFrame,
    manifest_row: pd.Series,
    columns: list[str],
) -> pd.DataFrame:
    identity = {
        "SchemaVersion": oc.SCHEMA_VERSION,
        "RunId": str(manifest_row["RunId"]),
        "ConditionId": str(manifest_row["ConditionId"]),
        "Design": str(manifest_row["Design"]),
        "TruthBias": float(manifest_row["TruthBias"]),
        "Replicate": int(manifest_row["Replicate"]),
        "Seed": int(manifest_row["Seed"]),
    }
    output = frame.copy()
    for column, value in reversed(list(identity.items())):
        output.insert(0, column, value)
    return output.reindex(columns=columns)


def _append_generated_bundle(
    output: Path,
    manifest_row: pd.Series,
    generated: GeneratedDesign,
    *,
    first: bool,
) -> dict[str, object]:
    """Append one run to atomic staging CSVs and return content identity."""

    frames = {
        "generated_ratings.csv": generated.data,
        "generated_facet_truth.csv": generated.facet_truth,
        "generated_anchors.csv": generated.anchors,
    }
    identified: dict[str, pd.DataFrame] = {}
    for filename, frame in frames.items():
        prepared = _identified_generated_frame(
            frame,
            manifest_row,
            GENERATED_BUNDLE_FILES[filename],
        )
        identified[filename] = prepared
        prepared.to_csv(
            output / f".{filename}.partial",
            mode="w" if first else "a",
            header=first,
            index=False,
        )

    ratings_fingerprint = oc.frame_fingerprint(generated.data.reset_index(drop=True))
    truth_fingerprint = oc.frame_fingerprint(generated.facet_truth.reset_index(drop=True))
    anchors_for_hash = generated.anchors.reindex(columns=["Facet", "Level", "Anchor"])
    anchor_fingerprint = oc.frame_fingerprint(anchors_for_hash.reset_index(drop=True))
    data_id = oc.canonical_condition_id(
        {
            "SeedGroup": manifest_row["SeedGroup"],
            "TruthBias": manifest_row["TruthBias"],
            "Replicate": manifest_row["Replicate"],
            "Seed": manifest_row["Seed"],
            "RatingsFingerprint": ratings_fingerprint,
            "FacetTruthFingerprint": truth_fingerprint,
        },
        prefix="generated-data",
    )
    fit_input_id = oc.canonical_condition_id(
        {"DataId": data_id, "AnchorFingerprint": anchor_fingerprint},
        prefix="fit-input",
    )
    return {
        "SchemaVersion": oc.SCHEMA_VERSION,
        "RunId": str(manifest_row["RunId"]),
        "ConditionId": str(manifest_row["ConditionId"]),
        "DataId": data_id,
        "FitInputId": fit_input_id,
        "SeedGroup": str(manifest_row["SeedGroup"]),
        "TruthBias": float(manifest_row["TruthBias"]),
        "Replicate": int(manifest_row["Replicate"]),
        "Seed": int(manifest_row["Seed"]),
        "RatingsRows": int(len(identified["generated_ratings.csv"])),
        "FacetTruthRows": int(len(identified["generated_facet_truth.csv"])),
        "AnchorRows": int(len(identified["generated_anchors.csv"])),
        "RatingsFingerprint": ratings_fingerprint,
        "FacetTruthFingerprint": truth_fingerprint,
        "AnchorFingerprint": anchor_fingerprint,
    }


def _finalize_generated_bundle(
    output: Path,
    identities: pd.DataFrame,
) -> pd.DataFrame:
    """Atomically publish staged bridge inputs and return their file hashes."""

    identity_name = "generated_data_identity.csv"
    identity_partial = output / f".{identity_name}.partial"
    identities.to_csv(identity_partial, index=False)
    identity_partial.replace(output / identity_name)

    inventory_rows: list[dict[str, object]] = []
    row_counts = {
        "generated_ratings.csv": int(identities["RatingsRows"].sum()),
        "generated_facet_truth.csv": int(identities["FacetTruthRows"].sum()),
        "generated_anchors.csv": int(identities["AnchorRows"].sum()),
        identity_name: int(len(identities)),
    }
    for filename in GENERATED_BUNDLE_FILES:
        partial = output / f".{filename}.partial"
        partial.replace(output / filename)
    for filename, rows in row_counts.items():
        path = output / filename
        inventory_rows.append({
            "SchemaVersion": oc.SCHEMA_VERSION,
            "File": filename,
            "Rows": rows,
            "SHA256": sha256_file(path),
            "Boundary": "Byte-level identity for cross-engine import; do not regenerate in R.",
        })
    inventory = pd.DataFrame(inventory_rows)
    inventory.to_csv(output / "generated_bundle_files.csv", index=False)
    return inventory


def validate_generated_bundle(output: Path) -> pd.DataFrame:
    """Fail closed when bridge files, hashes, run identities, or counts differ."""

    inventory = pd.read_csv(output / "generated_bundle_files.csv")
    identities = pd.read_csv(output / "generated_data_identity.csv")
    manifest = pd.read_csv(output / "manifest.csv")
    checks: list[dict[str, object]] = []

    def record(check: str, passed: bool, evidence: str) -> None:
        checks.append({"Check": check, "Passed": bool(passed), "Evidence": evidence})

    hashes_match = True
    for row in inventory.itertuples(index=False):
        path = output / str(row.File)
        actual = sha256_file(path) if path.exists() else "missing"
        hashes_match &= actual == str(row.SHA256)
    record("bundle_file_sha256", hashes_match, f"validated {len(inventory)} retained file hashes")

    manifest_runs = set(manifest["RunId"].astype(str))
    identity_runs = set(identities["RunId"].astype(str))
    record(
        "manifest_run_identity",
        manifest_runs == identity_runs and identities["RunId"].is_unique,
        f"manifest={len(manifest_runs)}; identity={len(identity_runs)}",
    )
    for filename, count_column in [
        ("generated_ratings.csv", "RatingsRows"),
        ("generated_facet_truth.csv", "FacetTruthRows"),
        ("generated_anchors.csv", "AnchorRows"),
    ]:
        frame = pd.read_csv(output / filename)
        actual = frame.groupby("RunId").size().reindex(identities["RunId"], fill_value=0).to_numpy()
        expected = pd.to_numeric(identities[count_column], errors="coerce").fillna(-1).to_numpy()
        record(
            f"{filename}_run_counts",
            bool(np.array_equal(actual, expected)),
            f"rows={len(frame)}; expected={int(expected.sum())}",
        )
    result = pd.DataFrame(checks)
    if not bool(result["Passed"].all()):
        failed = result.loc[~result["Passed"], "Check"].astype(str).tolist()
        raise ValueError(f"generated bundle validation failed: {failed}")
    return result


def study_conditions() -> list[dict[str, object]]:
    """Return paired null/alternative conditions across four design stresses."""

    designs = [
        {
            "Design": "balanced_small",
            "Persons": 30,
            "Raters": 4,
            "RatersPerPerson": 4,
            "MissingRate": 0.0,
            "AnchorShare": 0.0,
            "AnchorDrift": 0.0,
        },
        {
            "Design": "balanced_large_anchors",
            "SeedGroup": "large_anchor_sensitivity",
            "Persons": 80,
            "Raters": 4,
            "RatersPerPerson": 4,
            "MissingRate": 0.0,
            "AnchorShare": 0.5,
            "AnchorDrift": 0.0,
        },
        {
            "Design": "sparse_missing",
            "Persons": 18,
            "Raters": 8,
            "RatersPerPerson": 1,
            "MissingRate": 0.35,
            "AnchorShare": 0.0,
            "AnchorDrift": 0.0,
        },
        {
            "Design": "anchor_drift",
            "SeedGroup": "large_anchor_sensitivity",
            "Persons": 80,
            "Raters": 4,
            "RatersPerPerson": 4,
            "MissingRate": 0.0,
            "AnchorShare": 0.5,
            "AnchorDrift": 0.25,
        },
    ]
    rows: list[dict[str, object]] = []
    for design in designs:
        for truth_bias in (0.0, 0.6):
            condition_id = f"{design['Design']}__bias_{str(truth_bias).replace('.', 'p')}"
            rows.append({
                "ConditionId": condition_id,
                **design,
                "SeedGroup": str(design.get("SeedGroup", design["Design"])),
                "Tasks": 3,
                "Criteria": 2,
                "Categories": 4,
                "TruthBias": truth_bias,
                "TruthPositive": bool(abs(truth_bias) > 0),
                "BiasCell": f"{FOCAL_RATER} x {FOCAL_TASK}",
            })
    return rows


def _centred_normal(rng: np.random.Generator, count: int, sd: float) -> np.ndarray:
    values = rng.normal(0.0, float(sd), int(count))
    return values - float(np.mean(values))


def generate_condition_data(condition: pd.Series) -> GeneratedDesign:
    """Generate a paired RSM dataset with an optional focal local interaction."""

    rng = np.random.default_rng(int(condition["Seed"]))
    n_person = int(condition["Persons"])
    n_rater = int(condition["Raters"])
    raters_per_person = int(condition["RatersPerPerson"])
    n_task = int(condition["Tasks"])
    n_criterion = int(condition["Criteria"])
    n_cat = int(condition["Categories"])
    missing_rate = float(condition["MissingRate"])
    truth_bias = float(condition["TruthBias"])

    persons = [f"P{index:03d}" for index in range(1, n_person + 1)]
    raters = [f"R{index:02d}" for index in range(1, n_rater + 1)]
    tasks = [f"T{index:02d}" for index in range(1, n_task + 1)]
    criteria = [f"C{index:02d}" for index in range(1, n_criterion + 1)]
    theta = _centred_normal(rng, n_person, 1.0)
    rater = _centred_normal(rng, n_rater, 0.40)
    task = _centred_normal(rng, n_task, 0.30)
    criterion = _centred_normal(rng, n_criterion, 0.25)
    adjacent_steps = np.linspace(-1.2, 1.2, n_cat - 1)
    cumulative_steps = np.concatenate([[0.0], np.cumsum(adjacent_steps)])
    categories = np.arange(n_cat, dtype=int)

    design_rows: list[tuple[int, int, int, int]] = []
    for person_index in range(n_person):
        selected_raters = [
            (person_index + offset) % n_rater
            for offset in range(min(raters_per_person, n_rater))
        ]
        for rater_index in selected_raters:
            for task_index in range(n_task):
                for criterion_index in range(n_criterion):
                    design_rows.append((person_index, rater_index, task_index, criterion_index))

    uniforms = rng.random(len(design_rows))
    scores: list[int] = []
    for uniform, (person_index, rater_index, task_index, criterion_index) in zip(uniforms, design_rows):
        local_shift = (
            truth_bias
            if raters[rater_index] == FOCAL_RATER and tasks[task_index] == FOCAL_TASK
            else 0.0
        )
        eta = (
            theta[person_index]
            - rater[rater_index]
            - task[task_index]
            - criterion[criterion_index]
            + local_shift
        )
        log_weights = categories.astype(float) * eta - cumulative_steps
        log_weights -= float(np.max(log_weights))
        probabilities = np.exp(log_weights)
        probabilities /= float(probabilities.sum())
        scores.append(int(np.searchsorted(np.cumsum(probabilities), uniform, side="right")))

    data = pd.DataFrame({
        "Person": [persons[p] for p, _, _, _ in design_rows],
        "Rater": [raters[r] for _, r, _, _ in design_rows],
        "Task": [tasks[t] for _, _, t, _ in design_rows],
        "Criterion": [criteria[c] for _, _, _, c in design_rows],
        "Score": scores,
    })
    before_missing = int(len(data))
    if missing_rate > 0 and not data.empty:
        keep = rng.random(len(data)) >= missing_rate
        data = data.loc[keep].reset_index(drop=True)
    realized_missing = 1.0 - len(data) / max(before_missing, 1)

    truth_parts = [
        pd.DataFrame({"Facet": "Rater", "Level": raters, "Truth": rater}),
        pd.DataFrame({"Facet": "Task", "Level": tasks, "Truth": task}),
        pd.DataFrame({"Facet": "Criterion", "Level": criteria, "Truth": criterion}),
    ]
    truth = pd.concat(truth_parts, ignore_index=True)

    anchor_share = float(condition["AnchorShare"])
    n_anchor = int(math.ceil(anchor_share * n_rater)) if anchor_share > 0 else 0
    anchor_levels = raters[:n_anchor]
    anchors = truth.loc[
        (truth["Facet"] == "Rater") & truth["Level"].isin(anchor_levels),
        ["Facet", "Level", "Truth"],
    ].rename(columns={"Truth": "Anchor"})
    if not anchors.empty:
        anchors["Anchor"] = anchors["Anchor"] + float(condition["AnchorDrift"])

    counts = (
        data["Score"].value_counts().reindex(range(n_cat), fill_value=0)
        .rename_axis("Score").reset_index(name="Count")
    )
    return GeneratedDesign(
        data=data,
        facet_truth=truth,
        anchors=anchors,
        realized_missing_rate=realized_missing,
        category_counts=counts,
    )


def _summary_value(summary: pd.DataFrame, column: str, default: object = np.nan) -> object:
    if isinstance(summary, pd.DataFrame) and not summary.empty and column in summary.columns:
        return summary.iloc[0][column]
    return default


def _parameter_rows(
    manifest_row: pd.Series,
    diagnostics: dict,
    generated: GeneratedDesign,
    *,
    include: bool,
) -> pd.DataFrame:
    measures = diagnostics.get("measures", pd.DataFrame()) if isinstance(diagnostics, dict) else pd.DataFrame()
    if not isinstance(measures, pd.DataFrame) or measures.empty:
        return pd.DataFrame()
    level_column = "Level" if "Level" in measures.columns else "Element" if "Element" in measures.columns else None
    if level_column is None or not {"Facet", "Estimate"}.issubset(measures.columns):
        return pd.DataFrame()
    estimated = measures.loc[
        measures["Facet"].astype(str).isin(["Rater", "Task", "Criterion"]),
        [column for column in ["Facet", level_column, "Estimate", "SE"] if column in measures.columns],
    ].copy()
    estimated = estimated.rename(columns={level_column: "Level"})
    estimated["Facet"] = estimated["Facet"].astype(str)
    estimated["Level"] = estimated["Level"].astype(str)
    merged = generated.facet_truth.merge(estimated, on=["Facet", "Level"], how="left")
    merged["Estimate"] = pd.to_numeric(merged["Estimate"], errors="coerce")
    merged["Truth"] = pd.to_numeric(merged["Truth"], errors="coerce")
    merged["SE"] = pd.to_numeric(merged.get("SE", np.nan), errors="coerce")
    merged["RawError"] = merged["Estimate"] - merged["Truth"]
    merged["EstimateAligned"] = np.nan
    merged["TruthAligned"] = np.nan
    merged["ErrorAligned"] = np.nan
    merged["ComparisonScale"] = "mean_aligned_location"
    anchored_facets = set(generated.anchors["Facet"].astype(str)) if not generated.anchors.empty else set()
    for facet, indices in merged.groupby("Facet", sort=False).groups.items():
        subset = merged.loc[indices]
        finite = subset["Estimate"].notna() & subset["Truth"].notna()
        if str(facet) in anchored_facets:
            # Fixed anchors identify the absolute facet scale. Re-centering
            # would erase anchor contamination and can manufacture error even
            # when the supplied anchors equal truth.
            merged.loc[indices, "EstimateAligned"] = subset["Estimate"]
            merged.loc[indices, "TruthAligned"] = subset["Truth"]
            merged.loc[indices, "ErrorAligned"] = subset["RawError"]
            merged.loc[indices, "ComparisonScale"] = "anchor_identified_absolute"
        else:
            shift = float(subset.loc[finite, "RawError"].mean()) if finite.any() else np.nan
            merged.loc[indices, "EstimateAligned"] = subset["Estimate"] - shift
            merged.loc[indices, "TruthAligned"] = subset["Truth"]
            merged.loc[indices, "ErrorAligned"] = merged.loc[indices, "EstimateAligned"] - merged.loc[indices, "TruthAligned"]
    anchor_keys = set(map(tuple, generated.anchors[["Facet", "Level"]].astype(str).to_numpy())) if not generated.anchors.empty else set()
    merged["Anchored"] = [
        (str(facet), str(level)) in anchor_keys
        for facet, level in zip(merged["Facet"], merged["Level"])
    ]
    merged["IncludedInSummary"] = bool(include) & merged["ErrorAligned"].notna()
    identity_columns = [
        "RunId", "ConditionId", "Design", "TruthBias", "Replicate", "Seed",
    ]
    for column in reversed(identity_columns):
        merged.insert(0, column, manifest_row.get(column))
    merged.insert(2, "Engine", "PythonApp")
    merged.insert(3, "Estimator", "JMLE")
    merged["ParameterType"] = merged["Facet"]
    merged["CoverageMethod"] = np.where(
        merged["ComparisonScale"].eq("anchor_identified_absolute"),
        "anchor-identified absolute conditional-Wald diagnostic",
        "mean-aligned conditional-Wald diagnostic",
    )
    return merged


def fit_manifest_row(
    manifest_row: pd.Series,
    *,
    generated: GeneratedDesign | None = None,
) -> tuple[
    dict[str, object],
    pd.DataFrame,
    dict[str, object],
    pd.DataFrame,
    pd.DataFrame,
]:
    """Generate, fit, and evaluate one manifest row without dropping failures."""

    started = time.perf_counter()
    base = {column: manifest_row[column] for column in manifest_row.index}
    run: dict[str, object] = {
        **base,
        "Engine": "PythonApp",
        "Estimator": "JMLE",
        "FitReturned": False,
        "Converged": False,
        "InferenceReady": False,
        "BiasAvailable": False,
        "FocalCellSparse": True,
        "AnalysisEligible": False,
        "FailureStage": "generate",
        "FailureReason": "",
        "Warnings": "",
        "Rows": 0,
        "RealizedMissingRate": np.nan,
        "ZeroCategories": np.nan,
        "ZeroPairCells": np.nan,
        "LowPairCells": np.nan,
        "FitBoundaryStatistics": np.nan,
        "FitDisplayDecisionMismatches": np.nan,
        "BiasBoundaryStatistics": np.nan,
        "BiasDisplayDecisionMismatches": np.nan,
        "TerminalGradientSupNorm": np.nan,
        "GradientReviewTolerance": 1e-4,
        "GradientReady1e4": False,
        "GradientReady1e5": False,
        "GradientReady1e6": False,
        "StrictInferenceReady1e4": False,
        "StrictAnalysisEligible1e4": False,
        "ElapsedSeconds": np.nan,
    }
    bias_decision: dict[str, object] = {
        "RunId": manifest_row["RunId"],
        "ConditionId": manifest_row["ConditionId"],
        "Design": manifest_row["Design"],
        "TruthBias": float(manifest_row["TruthBias"]),
        "TruthPositive": bool(manifest_row["TruthPositive"]),
        "Replicate": int(manifest_row["Replicate"]),
        "Seed": int(manifest_row["Seed"]),
        "Engine": "PythonApp",
        "Estimator": "JMLE",
        "AnalysisEligible": False,
        "StrictAnalysisEligible1e4": False,
        "BiasEstimate": np.nan,
        "BiasSE": np.nan,
        "p_holm": np.nan,
        "AbsBias": np.nan,
        "SparseCell": True,
        "DecisionHolmRaw": pd.NA,
        "DecisionHolmDisplayed": pd.NA,
        "DecisionPracticalRaw": pd.NA,
        "DecisionPracticalDisplayed": pd.NA,
        "DecisionStrongRaw": pd.NA,
        "DecisionStrongDisplayed": pd.NA,
        "DecisionAnyNonSparseFlag": pd.NA,
    }
    diagnostics: dict = {}
    fit_audit = pd.DataFrame()
    bias_audit = pd.DataFrame()
    caught_messages: list[str] = []
    try:
        if generated is None:
            generated = generate_condition_data(manifest_row)
        run["FailureStage"] = "fit"
        run["Rows"] = int(len(generated.data))
        run["RealizedMissingRate"] = generated.realized_missing_rate
        run["ZeroCategories"] = int((generated.category_counts["Count"] == 0).sum())
        sparse_bundle = {
            "data": generated.data,
            "category_counts": generated.category_counts,
            "meta": {
                "facet_names": ["Rater", "Task", "Criterion"],
                "facet_level_counts": [
                    int(manifest_row["Raters"]),
                    int(manifest_row["Tasks"]),
                    int(manifest_row["Criteria"]),
                ],
                "full_rows": int(
                    int(manifest_row["Persons"])
                    * int(manifest_row["RatersPerPerson"])
                    * int(manifest_row["Tasks"])
                    * int(manifest_row["Criteria"])
                ),
                "missing_rate": float(manifest_row["MissingRate"]),
            },
        }
        pair_cells = app.build_custom_simulation_sparse_pair_cells(sparse_bundle)
        run["ZeroPairCells"] = int(pd.to_numeric(pair_cells.get("ZeroCells", 0), errors="coerce").fillna(0).sum()) if not pair_cells.empty else 0
        run["LowPairCells"] = int(pd.to_numeric(pair_cells.get("LowCountCells", 0), errors="coerce").fillna(0).sum()) if not pair_cells.empty else 0

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = app.mfrm_estimate(
                generated.data,
                person_col="Person",
                facet_cols=["Rater", "Task", "Criterion"],
                score_col="Score",
                rating_min=0,
                rating_max=int(manifest_row["Categories"]) - 1,
                model="RSM",
                method="JMLE",
                noncenter_facet="Person",
                anchor_df=generated.anchors if not generated.anchors.empty else None,
                anchor_policy="warn",
                min_common_anchors=2,
                maxit=160,
                reltol=1e-6,
                keep_original=True,
            )
            diagnostics = app.mfrm_diagnostics(
                result,
                compute_pca=False,
                compute_marginal=False,
            )
            caught_messages.extend(str(item.message) for item in caught)
        run["FitReturned"] = True
        summary = result.get("summary", pd.DataFrame())
        run["Converged"] = bool(_summary_value(summary, "Converged", False))
        run["InferenceReady"] = bool(
            _summary_value(summary, "InferenceReady", run["Converged"])
        )
        run["Iterations"] = _summary_value(summary, "Iterations", np.nan)
        run["GradientNorm"] = _summary_value(summary, "GradientNorm", np.nan)
        run["LogLik"] = _summary_value(summary, "LogLik", np.nan)
        optimizer = result.get("opt") if isinstance(result, dict) else None
        optimizer_gradient = getattr(optimizer, "jac", None)
        if optimizer_gradient is not None:
            gradient_values = np.asarray(optimizer_gradient, dtype=float).reshape(-1)
            finite_gradient = gradient_values[np.isfinite(gradient_values)]
            if finite_gradient.size:
                run["TerminalGradientSupNorm"] = float(np.max(np.abs(finite_gradient)))
        gradient_sup = float(run["TerminalGradientSupNorm"])
        run["GradientReady1e4"] = bool(
            np.isfinite(gradient_sup) and gradient_sup <= 1e-4
        )
        run["GradientReady1e5"] = bool(
            np.isfinite(gradient_sup) and gradient_sup <= 1e-5
        )
        run["GradientReady1e6"] = bool(
            np.isfinite(gradient_sup) and gradient_sup <= 1e-6
        )

        fit_source = diagnostics.get("fit", diagnostics.get("measures", pd.DataFrame()))
        fit_audit = ds.audit_fit_decision_stability(fit_source)
        fit_counts = ds.summarize_boundary_audit(fit_audit)
        run["FitBoundaryStatistics"] = int(
            fit_counts["numerical_boundary"] + fit_counts["display_rounding_boundary"]
        )
        run["FitDisplayDecisionMismatches"] = int(fit_counts["display_decision_mismatch"])

        run["FailureStage"] = "bias"
        bias_bundle = app.estimate_bias_interaction(
            result,
            diagnostics,
            "Rater",
            "Task",
            omit_extreme=False,
        )
        if isinstance(bias_bundle, dict) and "table" in bias_bundle:
            dff = app.build_dff_bias_screening_table(
                bias_bundle,
                alpha=0.05,
                min_n=5,
                practical_logit=0.50,
            )
            focal = dff.loc[
                (dff["FacetA_Level"].astype(str) == FOCAL_RATER)
                & (dff["FacetB_Level"].astype(str) == FOCAL_TASK)
            ]
            if not focal.empty:
                focal_row = focal.iloc[0]
                p_holm = float(pd.to_numeric(pd.Series([focal_row.get("p_holm")]), errors="coerce").iloc[0])
                abs_bias = float(pd.to_numeric(pd.Series([focal_row.get("AbsBias")]), errors="coerce").iloc[0])
                sparse = bool(focal_row.get("SparseCell", True))
                bias_decision.update({
                    "BiasEstimate": float(focal_row.get("BiasSize", np.nan)),
                    "BiasSE": float(focal_row.get("SE", np.nan)),
                    "p_holm": p_holm,
                    "AbsBias": abs_bias,
                    "SparseCell": sparse,
                    "DecisionHolmRaw": bool(np.isfinite(p_holm) and p_holm < 0.05),
                    "DecisionHolmDisplayed": bool(np.isfinite(p_holm) and round(p_holm, 4) < 0.05),
                    "DecisionPracticalRaw": bool(np.isfinite(abs_bias) and abs_bias >= 0.50),
                    "DecisionPracticalDisplayed": bool(np.isfinite(abs_bias) and round(abs_bias, 4) >= 0.50),
                    "DecisionAnyNonSparseFlag": bool(focal_row.get("Flag", False) and not sparse),
                })
                bias_decision["DecisionStrongRaw"] = bool(
                    bias_decision["DecisionHolmRaw"]
                    and bias_decision["DecisionPracticalRaw"]
                    and not sparse
                )
                bias_decision["DecisionStrongDisplayed"] = bool(
                    bias_decision["DecisionHolmDisplayed"]
                    and bias_decision["DecisionPracticalDisplayed"]
                    and not sparse
                )
                run["BiasAvailable"] = bool(np.isfinite(p_holm) and np.isfinite(abs_bias))
                run["FocalCellSparse"] = sparse
                boundary_audit = ds.audit_bias_decision_stability(dff)
                bias_audit = boundary_audit
                boundary_counts = ds.summarize_boundary_audit(bias_audit)
                run["BiasBoundaryStatistics"] = int(
                    boundary_counts["numerical_boundary"]
                    + boundary_counts["display_rounding_boundary"]
                )
                run["BiasDisplayDecisionMismatches"] = int(
                    boundary_counts["display_decision_mismatch"]
                )
            else:
                run["FailureReason"] = "focal Rater x Task cell missing from DFF table"
        else:
            run["FailureReason"] = str(
                bias_bundle.get("_skip_reason", "bias table unavailable")
                if isinstance(bias_bundle, dict) else
                "bias result unavailable"
            )

        eligible = bool(
            run["FitReturned"]
            and run["Converged"]
            and run["BiasAvailable"]
            and not run["FocalCellSparse"]
        )
        run["AnalysisEligible"] = eligible
        bias_decision["AnalysisEligible"] = eligible
        run["StrictInferenceReady1e4"] = bool(
            run["FitReturned"]
            and run["Converged"]
            and run["GradientReady1e4"]
        )
        run["StrictAnalysisEligible1e4"] = bool(
            eligible and run["StrictInferenceReady1e4"]
        )
        bias_decision["StrictAnalysisEligible1e4"] = run["StrictAnalysisEligible1e4"]
        if eligible:
            run["FailureStage"] = ""
            run["FailureReason"] = ""
        elif not run["FailureReason"]:
            if not run["Converged"]:
                run["FailureStage"] = "convergence"
                run["FailureReason"] = "fit returned without convergence"
            elif run["FocalCellSparse"]:
                run["FailureStage"] = "bias"
                run["FailureReason"] = "focal bias cell below min_n"
            else:
                run["FailureReason"] = "analysis eligibility not satisfied"
    except Exception as exc:  # retain every attempted run
        run["FailureReason"] = f"{type(exc).__name__}: {exc}"[:1000]
    finally:
        run["Warnings"] = " | ".join(caught_messages)[:2000]
        run["ElapsedSeconds"] = time.perf_counter() - started

    parameter_rows = (
        _parameter_rows(
            manifest_row,
            diagnostics,
            generated,
            include=bool(run["FitReturned"] and run["Converged"]),
        )
        if generated is not None else
        pd.DataFrame()
    )
    identity = {
        "RunId": manifest_row["RunId"],
        "ConditionId": manifest_row["ConditionId"],
        "Design": manifest_row["Design"],
        "TruthBias": float(manifest_row["TruthBias"]),
        "Replicate": int(manifest_row["Replicate"]),
        "Seed": int(manifest_row["Seed"]),
        "Engine": "PythonApp",
        "Estimator": "JMLE",
    }
    for audit_type, audit in (("fit", fit_audit), ("bias", bias_audit)):
        if audit.empty:
            continue
        audit.insert(0, "AuditType", audit_type)
        for column, value in reversed(identity.items()):
            audit.insert(0, column, value)
    return run, parameter_rows, bias_decision, fit_audit, bias_audit


def _plot_decision_rates(operating: pd.DataFrame, output: Path) -> None:
    selected = operating.loc[operating["DecisionRule"] == "DecisionStrongRaw"].copy()
    if selected.empty:
        return
    selected["Label"] = selected["Design"].astype(str) + "\nbias=" + selected["TruthBias"].astype(str)
    x = np.arange(len(selected))
    rate = pd.to_numeric(selected["DecisionRate"], errors="coerce").to_numpy(dtype=float)
    eligible = pd.to_numeric(selected["EligibleDecisions"], errors="coerce").fillna(0).to_numpy(dtype=int)
    lower = pd.to_numeric(selected["DecisionRateWilsonLower95"], errors="coerce").to_numpy(dtype=float)
    upper = pd.to_numeric(selected["DecisionRateWilsonUpper95"], errors="coerce").to_numpy(dtype=float)
    yerr = np.vstack([rate - lower, upper - rate])
    fig, ax = plt.subplots(figsize=(11, 5.5))
    colours = ["#d95f02" if bool(value) else "#1b9e77" for value in selected["TruthPositive"]]
    ax.bar(x, rate, color=colours, alpha=0.85)
    finite = np.isfinite(rate) & np.isfinite(yerr).all(axis=0)
    if finite.any():
        ax.errorbar(x[finite], rate[finite], yerr=yerr[:, finite], fmt="none", ecolor="#222222", capsize=4)
    for index in np.flatnonzero(eligible <= 0):
        ax.text(
            index,
            0.50,
            "No eligible\ndecisions",
            ha="center",
            va="center",
            color="#666666",
            fontsize=9,
            fontweight="bold",
        )
    ax.set_xticks(x, selected["Label"], rotation=30, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Eligible-decision rate")
    ax.set_title("Strong bias-screen decision rate (Wilson 95% interval; smoke only)")
    ax.text(
        0.01, -0.31,
        "Green=null condition; orange=positive local-bias condition. Pilot rates are not validation claims.",
        transform=ax.transAxes,
    )
    fig.tight_layout()
    fig.savefig(output / "decision_rates.png", dpi=180)
    plt.close(fig)


def _plot_run_accounting(accounting: pd.DataFrame, output: Path) -> None:
    if accounting.empty:
        return
    work = accounting.copy()
    work["Label"] = work["ConditionId"].astype(str)
    x = np.arange(len(work))
    width = 0.36
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.bar(x - width / 2, work["ConvergenceRate"], width, label="Converged", color="#4C78A8")
    ax.bar(x + width / 2, work["AnalysisEligibleRate"], width, label="Bias decision eligible", color="#F58518")
    ax.set_xticks(x, work["Label"], rotation=30, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Share of attempted runs")
    ax.set_title("Run accounting keeps non-convergence and unavailable decisions visible", pad=34)
    for index, value in enumerate(pd.to_numeric(work["AnalysisEligible"], errors="coerce").fillna(0)):
        if int(value) == 0:
            ax.text(index + width / 2, 0.025, "0 eligible", rotation=90, ha="center", va="bottom", color="#A64B00")
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2)
    fig.tight_layout()
    fig.savefig(output / "run_accounting.png", dpi=180)
    plt.close(fig)


def _plot_estimation(estimation: pd.DataFrame, output: Path) -> None:
    if estimation.empty:
        return
    work = estimation.copy()
    design_labels = {
        "balanced_small": "balanced small",
        "balanced_large_anchors": "correct anchors",
        "sparse_missing": "sparse + missing",
        "anchor_drift": "anchor drift",
    }
    work["Label"] = (
        work["Design"].astype(str).map(design_labels).fillna(work["Design"].astype(str))
        + "\n" + work["ParameterType"].astype(str)
        + "; bias=" + work["TruthBias"].astype(str)
    )
    x = np.arange(len(work))
    fig, axes = plt.subplots(2, 1, figsize=(15, 8), sharex=True)
    axes[0].bar(x, work["RMSE"], color="#4C78A8")
    axes[0].set_ylabel("RMSE (logits; scale in CSV)")
    axes[0].set_title("Parameter recovery among included returned fits (descriptive smoke only)")
    axes[1].bar(x, work["Coverage"], color="#59A14F")
    axes[1].axhline(0.95, color="#B22222", linestyle="--", linewidth=1.2, label="Nominal .95")
    axes[1].set_ylim(0, 1.05)
    axes[1].set_ylabel("Conditional-Wald coverage")
    axes[1].legend(frameon=False)
    axes[1].set_xticks(x, work["Label"], rotation=52, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "estimation_rmse_coverage.png", dpi=180)
    plt.close(fig)


def write_results(
    output: Path,
    *,
    profile: str,
    manifest: pd.DataFrame,
    runs: pd.DataFrame,
    operating: pd.DataFrame,
    estimation: pd.DataFrame,
    accounting: pd.DataFrame,
    sensitivity: pd.DataFrame,
    first_read: pd.DataFrame,
    precision_plan_sha256: str,
    parent_manifest_path: Path | None,
) -> str:
    attempts = int(len(runs))
    returned = int(runs["FitReturned"].astype(bool).sum()) if attempts else 0
    converged = int(runs["Converged"].astype(bool).sum()) if attempts else 0
    eligible = int(runs["AnalysisEligible"].astype(bool).sum()) if attempts else 0
    fit_boundary = int(pd.to_numeric(runs["FitBoundaryStatistics"], errors="coerce").fillna(0).sum()) if attempts else 0
    fit_mismatch = int(pd.to_numeric(runs["FitDisplayDecisionMismatches"], errors="coerce").fillna(0).sum()) if attempts else 0
    bias_boundary = int(pd.to_numeric(runs["BiasBoundaryStatistics"], errors="coerce").fillna(0).sum()) if attempts else 0
    bias_mismatch = int(pd.to_numeric(runs["BiasDisplayDecisionMismatches"], errors="coerce").fillna(0).sum()) if attempts else 0
    max_zero = int(pd.to_numeric(runs["ZeroPairCells"], errors="coerce").fillna(0).max()) if attempts else 0
    max_low = int(pd.to_numeric(runs["LowPairCells"], errors="coerce").fillna(0).max()) if attempts else 0
    changed = int(pd.to_numeric(sensitivity.get("ChangedConclusions", 0), errors="coerce").fillna(0).sum()) if not sensitivity.empty else 0
    elapsed = pd.to_numeric(runs.get("ElapsedSeconds", np.nan), errors="coerce")
    elapsed = elapsed.loc[elapsed.notna() & np.isfinite(elapsed) & (elapsed >= 0)]
    total_elapsed = float(elapsed.sum()) if len(elapsed) else np.nan
    median_elapsed = float(elapsed.median()) if len(elapsed) else np.nan
    p95_elapsed = float(elapsed.quantile(0.95)) if len(elapsed) else np.nan
    evidence_tiers = sorted(set(operating.get("EvidenceTier", pd.Series(dtype=str)).astype(str)))
    first_read_lines = "\n".join(
        f"- **{row['Check']} — {row['Status']}:** {row['Evidence']}"
        for _, row in first_read.sort_values("Priority").iterrows()
    )
    try:
        output_display = output.relative_to(REPO)
    except ValueError:
        output_display = output
    parent_argument = ""
    if parent_manifest_path is not None:
        try:
            parent_display = parent_manifest_path.relative_to(REPO)
        except ValueError:
            parent_display = parent_manifest_path
        parent_argument = f" --parent-manifest {parent_display}"

    def rater_rmse(design: str, truth_bias: float) -> float:
        if estimation.empty:
            return np.nan
        selected = estimation.loc[
            estimation["Design"].astype(str).eq(design)
            & estimation["ParameterType"].astype(str).eq("Rater")
            & pd.to_numeric(estimation["TruthBias"], errors="coerce").eq(float(truth_bias))
        ]
        return float(selected.iloc[0]["RMSE"]) if not selected.empty else np.nan

    clean_anchor_rmse = rater_rmse("balanced_large_anchors", 0.0)
    drift_anchor_rmse = rater_rmse("anchor_drift", 0.0)
    sparse_rater_rmse = rater_rmse("sparse_missing", 0.0)
    diagnostic_lines = "\n".join([
        f"- Clean-anchor Rater RMSE (null condition): {clean_anchor_rmse:.4f} logits.",
        f"- +0.25-logit contaminated-anchor Rater RMSE (paired null condition): {drift_anchor_rmse:.4f} logits.",
        f"- Sparse-design Rater RMSE (null condition): {sparse_rater_rmse:.4f} logits.",
    ])
    text = f"""# Python operating-characteristics pilot

Schema: `{oc.SCHEMA_VERSION}`  
Profile: `{profile}`  
Platform: `{platform.platform()}`  
Python: `{platform.python_version()}`  
Manifest SHA-256: `{oc.frame_fingerprint(manifest)}`
Registered precision-plan SHA-256: `{precision_plan_sha256}`

## Outcome

- {len(manifest)} attempted condition-replicate fits across {manifest['ConditionId'].nunique()} conditions.
- {returned}/{attempts} fits returned; {converged}/{attempts} converged; {eligible}/{attempts} produced an eligible focal bias decision.
- Recorded fit time was {total_elapsed:.1f} seconds in total (median {median_elapsed:.3f}, 95th percentile {p95_elapsed:.3f} seconds/run); wall-clock orchestration also includes bundle and plot work.
- The sparsest retained design had up to {max_zero} empty and {max_low} low-count facet-pair cells.
- Fit boundary auditing found {fit_boundary} near-threshold statistic(s), including {fit_mismatch} raw/display classification mismatch(es).
- Bias boundary auditing found {bias_boundary} near-threshold statistic(s), including {bias_mismatch} raw/display decision mismatch(es), across full bias-cell families.
- Conclusion sensitivity found {changed} focal raw-versus-displayed decision change(s).
- Evidence tier(s): {', '.join(evidence_tiers) or 'none'}.

## First read

{first_read_lines}

## Descriptive stress signals

{diagnostic_lines}

These values are diagnostic pilot signals from paired generated conditions,
not stable estimates of estimator performance. Their purpose is to verify that
the negative controls are visible before spending a larger simulation budget.

## Critical interpretation

This is an orchestration and failure-accounting pilot. The `{profile}` profile
does not estimate stable false-positive rates, power, coverage, or sample-size
requirements. Failed fits, non-converged fits, sparse focal cells, and
unavailable decisions remain explicit and are never counted as negative bias
decisions.

The local-bias generator uses common random numbers for each paired null and
0.60-logit alternative. Python JMLE is the only engine in this increment.
Cross-engine operating-characteristic claims require appending mfrmr, TAM,
immer, and sirt results under the same manifest and decision schema.

Anchor values are hard constraints. The drift condition deliberately supplies
anchors shifted by +0.25 logits; its results study sensitivity to contaminated
linking inputs, not empirical anchor invariance. Recovery uses absolute errors
for anchor-identified facets and mean-aligned errors for unidentified facet
locations. Conditional-Wald coverage includes only returned finite positive
SEs and does not propagate alignment or model-selection uncertainty.

## Retained evidence

- `manifest.csv` and `conditions.csv`
- `registered_precision_plan.json` and `study_identity.json`
- `manifest_extension_audit.json` when a parent manifest is supplied
- `generated_ratings.csv`, `generated_facet_truth.csv`, and `generated_anchors.csv`
- `generated_data_identity.csv`, `generated_bundle_files.csv`, and `bridge_validation_python.csv`
- `runs.csv`, `failure_accounting.csv`, and `failure_reasons.csv`
- `parameter_recovery.csv` and `estimation_operating_characteristics.csv`
- `bias_decisions.csv` and `binary_operating_characteristics.csv`
- `conclusion_sensitivity.csv`
- `fit_decision_stability.csv` and `bias_decision_stability.csv`
- `first_read_summary.csv`
- `decision_rates.png`, `run_accounting.png`, and `estimation_rmse_coverage.png`

Reproduce the retained profile with:

```bash
MPLCONFIGDIR=/tmp/mfrm_matplotlib_cache \\
  python3 validation/operating_characteristics_pilot.py --profile {profile} \\
  --output {output_display}{parent_argument}
```
"""
    (output / "RESULTS.md").write_text(text, encoding="utf-8")
    return text


def run_study(
    *,
    profile: str,
    replicates: int,
    base_seed: int,
    output: Path,
    max_runs: int,
    precision_plan_path: Path = DEFAULT_PRECISION_PLAN,
    parent_manifest_path: Path | None = None,
) -> str:
    precision_plan_path = precision_plan_path.resolve()
    precision_plan = oc.load_precision_plan(precision_plan_path)
    registered_replicates = int(
        precision_plan["profiles"][profile]["replicates_per_condition"]
    )
    conditions = study_conditions()
    manifest = oc.build_replicate_manifest(
        conditions,
        replicates=replicates,
        base_seed=base_seed,
        max_runs=max_runs,
    )
    manifest_extension_audit: dict[str, object] | None = None
    if parent_manifest_path is not None:
        parent_manifest_path = parent_manifest_path.resolve()
        parent_manifest = pd.read_csv(parent_manifest_path)
        manifest_extension_audit = oc.audit_manifest_extension(parent_manifest, manifest)
    output.mkdir(parents=True, exist_ok=True)
    precision_plan_bytes = precision_plan_path.read_bytes()
    (output / "registered_precision_plan.json").write_bytes(precision_plan_bytes)
    precision_plan_sha256 = hashlib.sha256(precision_plan_bytes).hexdigest()
    try:
        precision_plan_source = str(precision_plan_path.relative_to(REPO))
    except ValueError:
        precision_plan_source = precision_plan_path.name
    # A regenerated Python bundle invalidates every derived R-side validation
    # artifact until the byte-hash bridge is run again.
    for filename in DERIVED_CROSS_ENGINE_OUTPUT_FILES:
        (output / filename).unlink(missing_ok=True)
    pd.DataFrame(conditions).to_csv(output / "conditions.csv", index=False)
    manifest.to_csv(output / "manifest.csv", index=False)
    if manifest_extension_audit is not None:
        (output / "manifest_extension_audit.json").write_text(
            json.dumps(manifest_extension_audit, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    study_identity = {
        "schema_version": oc.SCHEMA_VERSION,
        "profile": profile,
        "replicates": int(replicates),
        "registered_profile_replicates": registered_replicates,
        "profile_replicate_override": bool(int(replicates) != registered_replicates),
        "base_seed": int(base_seed),
        "manifest_sha256": oc.frame_fingerprint(manifest),
        "precision_plan_schema_version": precision_plan["schema_version"],
        "precision_plan_source": precision_plan_source,
        "precision_plan_sha256": precision_plan_sha256,
        "parent_manifest_sha256": (
            manifest_extension_audit["PriorManifestSHA256"]
            if manifest_extension_audit is not None else
            None
        ),
        "manifest_extension_audit_passed": (
            bool(manifest_extension_audit["Passed"])
            if manifest_extension_audit is not None else
            None
        ),
        "engine": "PythonApp",
        "estimator": "JMLE",
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "streamlit_app_sha256": sha256_file(REPO / "streamlit_app.py"),
        "operating_characteristics_module_sha256": sha256_file(
            REPO / "mfrm_app" / "operating_characteristics.py"
        ),
        "decision_stability_module_sha256": sha256_file(
            REPO / "mfrm_app" / "decision_stability.py"
        ),
        "numerical_controls": {
            "maxit": 160,
            "reltol": 1e-6,
            "gradient_sup_norm_sensitivity_thresholds": [1e-4, 1e-5, 1e-6],
            "strict_gradient_sensitivity_is_post_pilot_audit_not_registered_primary": True,
        },
    }

    run_rows: list[dict[str, object]] = []
    parameter_parts: list[pd.DataFrame] = []
    bias_rows: list[dict[str, object]] = []
    fit_audit_parts: list[pd.DataFrame] = []
    bias_audit_parts: list[pd.DataFrame] = []
    generated_identity_rows: list[dict[str, object]] = []
    for run_index, (_, manifest_row) in enumerate(manifest.iterrows()):
        generated = generate_condition_data(manifest_row)
        generated_identity_rows.append(
            _append_generated_bundle(
                output,
                manifest_row,
                generated,
                first=run_index == 0,
            )
        )
        run, parameters, bias, fit_audit, bias_audit = fit_manifest_row(
            manifest_row,
            generated=generated,
        )
        run_rows.append(run)
        if not parameters.empty:
            parameter_parts.append(parameters)
        bias_rows.append(bias)
        if not fit_audit.empty:
            fit_audit_parts.append(fit_audit)
        if not bias_audit.empty:
            bias_audit_parts.append(bias_audit)

    generated_identities = pd.DataFrame(generated_identity_rows)
    bundle_inventory = _finalize_generated_bundle(output, generated_identities)
    bridge_validation = validate_generated_bundle(output)
    bridge_validation.to_csv(output / "bridge_validation_python.csv", index=False)
    study_identity["generated_bundle_inventory_sha256"] = oc.frame_fingerprint(bundle_inventory)
    study_identity["generated_bundle_files"] = {
        str(row.File): str(row.SHA256)
        for row in bundle_inventory.itertuples(index=False)
    }
    (output / "study_identity.json").write_text(
        json.dumps(study_identity, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    runs = pd.DataFrame(run_rows)
    parameters = pd.concat(parameter_parts, ignore_index=True) if parameter_parts else pd.DataFrame()
    bias_decisions = pd.DataFrame(bias_rows)
    fit_decision_stability = (
        pd.concat(fit_audit_parts, ignore_index=True)
        if fit_audit_parts else
        pd.DataFrame()
    )
    bias_decision_stability = (
        pd.concat(bias_audit_parts, ignore_index=True)
        if bias_audit_parts else
        pd.DataFrame()
    )
    runs.to_csv(output / "runs.csv", index=False)
    parameters.to_csv(output / "parameter_recovery.csv", index=False)
    bias_decisions.to_csv(output / "bias_decisions.csv", index=False)
    fit_decision_stability.to_csv(output / "fit_decision_stability.csv", index=False)
    bias_decision_stability.to_csv(output / "bias_decision_stability.csv", index=False)

    operating = oc.summarize_binary_operating_characteristics(
        bias_decisions,
        decision_columns=[
            "DecisionHolmRaw",
            "DecisionPracticalRaw",
            "DecisionStrongRaw",
            "DecisionAnyNonSparseFlag",
        ],
        group_columns=["ConditionId", "Design", "TruthBias", "TruthPositive", "Engine", "Estimator"],
    )
    estimation = oc.summarize_estimation_operating_characteristics(
        parameters,
        group_columns=[
            "ConditionId", "Design", "TruthBias", "Engine", "Estimator",
            "ParameterType", "ComparisonScale",
        ],
    ) if not parameters.empty else pd.DataFrame()
    accounting, failure_reasons = oc.summarize_run_accounting(
        runs,
        group_columns=["ConditionId", "Design", "TruthBias", "Engine", "Estimator"],
    )
    sensitivity = oc.audit_conclusion_sensitivity(
        bias_decisions,
        raw_column="DecisionStrongRaw",
        comparison_columns=["DecisionStrongDisplayed"],
        group_columns=["ConditionId", "Design", "TruthBias", "Engine", "Estimator"],
    )
    first_read = oc.build_operating_characteristics_first_read(
        accounting,
        operating,
        estimation,
        profile=profile,
    )
    operating.to_csv(output / "binary_operating_characteristics.csv", index=False)
    estimation.to_csv(output / "estimation_operating_characteristics.csv", index=False)
    accounting.to_csv(output / "failure_accounting.csv", index=False)
    failure_reasons.to_csv(output / "failure_reasons.csv", index=False)
    sensitivity.to_csv(output / "conclusion_sensitivity.csv", index=False)
    first_read.to_csv(output / "first_read_summary.csv", index=False)

    _plot_decision_rates(operating, output)
    _plot_run_accounting(accounting, output)
    _plot_estimation(estimation, output)
    return write_results(
        output,
        profile=profile,
        manifest=manifest,
        runs=runs,
        operating=operating,
        estimation=estimation,
        accounting=accounting,
        sensitivity=sensitivity,
        first_read=first_read,
        precision_plan_sha256=precision_plan_sha256,
        parent_manifest_path=parent_manifest_path,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=sorted(PROFILE_REPLICATES), default="smoke")
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument("--base-seed", type=int, default=20260809)
    parser.add_argument("--max-runs", type=int, default=1_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--precision-plan", type=Path, default=DEFAULT_PRECISION_PLAN)
    parser.add_argument(
        "--parent-manifest",
        type=Path,
        default=None,
        help="Optional smaller manifest whose rows and seeds must be preserved.",
    )
    args = parser.parse_args()
    replicates = int(args.replicates or PROFILE_REPLICATES[args.profile])
    text = run_study(
        profile=args.profile,
        replicates=replicates,
        base_seed=int(args.base_seed),
        output=args.output.resolve(),
        max_runs=int(args.max_runs),
        precision_plan_path=args.precision_plan,
        parent_manifest_path=args.parent_manifest,
    )
    try:
        print(text)
    except UnicodeEncodeError:
        # Japanese Windows commonly exposes a CP932 console. The retained
        # UTF-8 report is already written; console rendering must not turn a
        # completed simulation into a false process failure.
        encoding = getattr(sys.stdout, "encoding", None) or "utf-8"
        safe_text = text.encode(encoding, errors="backslashreplace").decode(encoding)
        print(safe_text)


if __name__ == "__main__":
    main()
