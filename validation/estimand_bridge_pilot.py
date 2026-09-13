#!/usr/bin/env python3
"""Bridge FACETS/Python JMLE, native MML, and exact CMLE estimands.

The script consumes the retained rank-full rows from the PCM boundary pilot.
It does not regenerate data and never treats cross-estimator differences as
direct engine-parity errors.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import sys
import time
from typing import Any, Iterable
import warnings

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from mfrm_app.cmle import audit_cmle_eligibility, fit_cmle  # noqa: E402
from validation.facets_pcm_boundary_pilot import (  # noqa: E402
    FIT_MODELS,
    normalize_rsm_steps,
)
from validation.facets_pcm_probe import normalize_python_pcm_steps  # noqa: E402
from validation.operating_characteristics_facets import (  # noqa: E402
    RECOVERY_FACETS,
    align_recovery,
    sha256_file,
    validate_bundle,
)
from validation.operating_characteristics_python_mml import (  # noqa: E402
    MODE_SPECS,
    build_mml_kwargs,
)


SCHEMA_VERSION = "mfrm-estimand-bridge-pilot-v1"
SELECTED_DESIGNS = ("complete", "planned_connected")
MML_MODES = ("PYTHON_MML_FIXED_SD1_Q31", "PYTHON_MML_FREE_SD_Q31")
CMLE_MODE = "PYTHON_EXACT_CMLE"
PARENT_FACETS_MODE = "FACETS_4_5_JMLE_FROM_PARENT"
PARENT_PYTHON_MODE = "PYTHON_JMLE_FROM_PARENT"
CONSTRAINT_TOLERANCE = 1e-8


def validate_parent_identity(parent_dir: Path) -> dict[str, str]:
    """Re-hash retained inputs against the parent boundary identity."""

    identity_path = parent_dir / "boundary_identity.json"
    if not identity_path.is_file():
        raise FileNotFoundError(f"Parent identity is missing: {identity_path}")
    identity = json.loads(identity_path.read_text(encoding="utf-8"))
    expected = identity.get("retained_input_sha256", {})
    if not expected:
        raise ValueError("Parent identity does not contain retained input hashes")
    observed = {}
    for filename, digest in expected.items():
        path = parent_dir / "retained_input" / filename
        if not path.is_file():
            raise FileNotFoundError(f"Parent retained input is missing: {path}")
        actual = sha256_file(path)
        if actual.lower() != str(digest).lower():
            raise ValueError(
                f"Parent retained input hash mismatch for {filename}: {actual} != {digest}"
            )
        observed[filename] = actual
    return observed


def selected_manifest(parent_dir: Path) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    tables = validate_bundle(parent_dir / "retained_input")
    manifest = tables["manifest.csv"].copy()
    selected = manifest[manifest["Design"].astype(str).isin(SELECTED_DESIGNS)].copy()
    if len(selected) != 8:
        raise ValueError(f"Expected 8 rank-full parent runs, found {len(selected)}")
    expected_cells = selected.groupby(["Design", "ThresholdCondition", "Replicate"]).size()
    if len(expected_cells) != 8 or not expected_cells.eq(1).all():
        raise ValueError("Selected manifest is not the locked 2x2x2 bridge design")
    if selected["ExpectedStructuralNullity"].astype(int).ne(0).any():
        raise ValueError("Selected bridge manifest includes a structurally nonidentified run")
    return tables, selected


def _constraint_audit(
    facets: pd.DataFrame,
    steps: pd.DataFrame,
    *,
    estimator_mode: str,
    fit_model: str,
) -> dict[str, Any]:
    """Check only the sum-to-zero constraints declared by each estimator."""

    facet_sums = facets.groupby("Facet")["Estimate"].sum().to_dict()
    if estimator_mode in MML_MODES:
        # MML fixes the latent population mean and leaves Criterion as the
        # explicit location-bearing observed facet.
        constrained_facet_sums = {
            facet: value for facet, value in facet_sums.items()
            if facet in {"Rater", "Task"}
        }
    else:
        constrained_facet_sums = facet_sums
    if fit_model == "PCM":
        step_sums = steps.groupby("StepFacetLevel")["Estimate"].sum().to_dict()
    else:
        step_sums = {"__COMMON__": float(steps["Estimate"].sum())}
    residuals = [abs(float(value)) for value in constrained_facet_sums.values()]
    residuals.extend(abs(float(value)) for value in step_sums.values())
    maximum = max(residuals) if residuals else np.nan
    return {
        "ConstrainedFacetSums": constrained_facet_sums,
        "StepSums": step_sums,
        "MaxAbsConstraintResidual": maximum,
        "ConstraintPass": bool(np.isfinite(maximum) and maximum <= CONSTRAINT_TOLERANCE),
    }


def _normalize_steps(
    result: dict[str, Any],
    *,
    fit_model: str,
    estimator_mode: str,
    run_id: str,
    threshold_condition: str,
    threshold_truth: pd.DataFrame,
) -> pd.DataFrame:
    if estimator_mode == CMLE_MODE:
        source = result["steps"].copy()
        output = pd.DataFrame({
            "StepFacetLevel": (
                source["Level"].astype(str)
                if fit_model == "PCM" else pd.Series("__COMMON__", index=source.index)
            ),
            "Category": pd.to_numeric(source["Step"], errors="raise").astype(int),
            "Step": source["Step"],
            "Estimate": pd.to_numeric(source["Estimate"], errors="coerce"),
            "SE": pd.to_numeric(source["SE"], errors="coerce"),
        })
    elif fit_model == "PCM":
        output = normalize_python_pcm_steps(result).rename(
            columns={"PythonThreshold": "Estimate"}
        )
        output["SE"] = np.nan
    else:
        output = normalize_rsm_steps(result).rename(
            columns={"PythonThreshold": "Estimate"}
        )
        output["StepFacetLevel"] = "__COMMON__"
        output["SE"] = np.nan
    output.insert(0, "RunId", run_id)
    output.insert(1, "FitModel", fit_model)
    output.insert(2, "EstimatorMode", estimator_mode)
    output["ThresholdCondition"] = threshold_condition
    if fit_model == "PCM":
        output = output.merge(
            threshold_truth[["StepFacetLevel", "Category", "ThresholdTruth"]],
            on=["StepFacetLevel", "Category"],
            how="left",
            validate="one_to_one",
        )
        output["TruthTargetStatus"] = "data_generating_pcm_threshold"
    elif threshold_condition == "shared":
        common = threshold_truth.groupby("Category", as_index=False)["ThresholdTruth"].mean()
        output = output.merge(common, on="Category", how="left", validate="one_to_one")
        output["TruthTargetStatus"] = "data_generating_common_threshold"
    else:
        output["ThresholdTruth"] = np.nan
        output["TruthTargetStatus"] = "undefined_under_pcm_to_rsm_misspecification"
    output["TruthError"] = output["Estimate"] - output["ThresholdTruth"]
    return output


def _normalize_recovery(
    estimates: pd.DataFrame,
    truth: pd.DataFrame,
    manifest_row: pd.Series,
    *,
    estimator_family: str,
    estimator_mode: str,
    fit_model: str,
    included: bool,
) -> pd.DataFrame:
    source = estimates[[
        column for column in ("Facet", "Level", "Estimate", "SE", "Status")
        if column in estimates.columns
    ]].copy()
    if "SE" not in source:
        source["SE"] = np.nan
    if "Status" not in source:
        source["Status"] = np.nan
    recovery = align_recovery(source, truth, pd.DataFrame(), manifest_row)
    recovery["EstimatorFamily"] = estimator_family
    recovery["EstimatorMode"] = estimator_mode
    recovery["FitModel"] = fit_model
    recovery["ThresholdCondition"] = manifest_row["ThresholdCondition"]
    recovery["Engine"] = "PythonApp"
    recovery["Estimator"] = estimator_family
    recovery["Mode"] = estimator_mode
    recovery["IncludedInBridge"] = bool(included) & recovery["ErrorAligned"].notna()
    recovery["EstimandClass"] = (
        "fixed_person_joint_likelihood" if estimator_family == "JMLE"
        else "normal_person_population_marginal_likelihood" if estimator_family == "MML"
        else "person_total_conditional_likelihood"
    )
    return recovery


def import_parent_jmle(
    parent_dir: Path,
    selected_ids: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ledger = pd.read_csv(parent_dir / "boundary_fit_ledger.csv")
    ledger = ledger[ledger["RunId"].astype(str).isin(selected_ids)].copy()
    if len(ledger) != 16 or not ledger["ComparisonEligible"].all():
        raise ValueError("All 16 selected parent JMLE contracts must remain eligible")
    main = pd.read_csv(parent_dir / "boundary_main_pairs.csv")
    main = main[main["RunId"].astype(str).isin(selected_ids)].copy()
    metadata = ledger[["RunId", "ThresholdCondition"]].drop_duplicates("RunId")
    main = main.merge(metadata, on="RunId", how="left", validate="many_to_one")
    rows = []
    for mode, estimate_column, error_column, engine in (
        (PARENT_FACETS_MODE, "FACETSEstimateAligned", "FACETSTruthErrorAligned", "FACETS"),
        (PARENT_PYTHON_MODE, "PythonEstimateAligned", "PythonTruthErrorAligned", "PythonApp"),
    ):
        part = pd.DataFrame({
            "RunId": main["RunId"],
            "ConditionId": main["ConditionId"],
            "Design": main["Design"],
            "Replicate": main["Replicate"],
            "Seed": main["Seed"],
            "FitModel": main["FitModel"],
            "ThresholdCondition": main["ThresholdCondition"],
            "EstimatorFamily": "JMLE",
            "EstimatorMode": mode,
            "EstimandClass": "fixed_person_joint_likelihood",
            "Engine": engine,
            "Facet": main["Facet"],
            "Level": main["Level"],
            "Truth": main["Truth"],
            "EstimateAligned": main[estimate_column],
            "TruthAligned": main["TruthAligned"],
            "ErrorAligned": main[error_column],
            "SE": main["SE"] if mode == PARENT_FACETS_MODE else np.nan,
            "IncludedInBridge": main["ComparisonEligible"],
        })
        rows.append(part)
    parent_recovery = pd.concat(rows, ignore_index=True)

    thresholds = pd.read_csv(parent_dir / "boundary_threshold_pairs.csv")
    thresholds["RunId"] = thresholds["FitId"].astype(str).str.rsplit("::", n=1).str[0]
    thresholds = thresholds[thresholds["RunId"].isin(selected_ids)].copy()
    threshold_metadata = ledger[
        ["RunId", "Design", "ThresholdCondition", "Replicate", "Seed"]
    ].drop_duplicates("RunId")
    thresholds = thresholds.merge(
        threshold_metadata, on="RunId", how="left", validate="many_to_one"
    )
    threshold_rows = []
    for mode, estimate_column, error_column, engine in (
        (PARENT_FACETS_MODE, "ThresholdMeasureDisplayed", "FACETSTruthError", "FACETS"),
        (PARENT_PYTHON_MODE, "PythonThreshold", "PythonTruthError", "PythonApp"),
    ):
        part = pd.DataFrame({
            "RunId": thresholds["RunId"],
            "FitModel": thresholds["FitModel"],
            "Design": thresholds["Design"],
            "ThresholdCondition": thresholds["ThresholdCondition"],
            "Replicate": thresholds["Replicate"],
            "Seed": thresholds["Seed"],
            "EstimatorFamily": "JMLE",
            "EstimatorMode": mode,
            "EstimandClass": "fixed_person_joint_likelihood",
            "Engine": engine,
            "StepFacetLevel": thresholds["StepFacetLevel"],
            "Category": thresholds["Category"],
            "Estimate": thresholds[estimate_column],
            "SE": np.nan,
            "ThresholdTruth": thresholds["ThresholdTruth"],
            "TruthTargetStatus": thresholds["TruthTargetStatus"],
            "TruthError": thresholds[error_column],
            "IncludedInBridge": thresholds["ComparisonEligible"],
        })
        threshold_rows.append(part)
    parent_thresholds = pd.concat(threshold_rows, ignore_index=True)
    return ledger, parent_recovery, parent_thresholds


def fit_mml(
    app: Any,
    ratings: pd.DataFrame,
    truth: pd.DataFrame,
    threshold_truth: pd.DataFrame,
    manifest_row: pd.Series,
    *,
    fit_model: str,
    mode: str,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    started = time.perf_counter()
    kwargs = build_mml_kwargs(mode, int(manifest_row["Categories"]))
    kwargs["model"] = fit_model
    if fit_model == "PCM":
        kwargs["step_facet"] = "Criterion"
    warning_messages = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = app.mfrm_estimate(
            data=ratings[["Person", "Rater", "Task", "Criterion", "Score"]].copy(),
            **kwargs,
        )
        warning_messages = [str(item.message) for item in caught]
    summary = result["summary"].iloc[0]
    returned = True
    converged = bool(summary.get("Converged", False))
    ready = bool(summary.get("InferenceReady", False))
    facets = result["facets"]["others"].copy()
    finite_main = pd.to_numeric(facets["Estimate"], errors="coerce").notna().all()
    thresholds = _normalize_steps(
        result,
        fit_model=fit_model,
        estimator_mode=mode,
        run_id=str(manifest_row["RunId"]),
        threshold_condition=str(manifest_row["ThresholdCondition"]),
        threshold_truth=threshold_truth,
    )
    finite_steps = pd.to_numeric(thresholds["Estimate"], errors="coerce").notna().all()
    constraint = _constraint_audit(
        facets, thresholds.rename(columns={"Estimate": "Estimate"}),
        estimator_mode=mode, fit_model=fit_model,
    )
    included = bool(returned and converged and ready and finite_main and finite_steps and constraint["ConstraintPass"])
    recovery = _normalize_recovery(
        facets,
        truth,
        manifest_row,
        estimator_family="MML",
        estimator_mode=mode,
        fit_model=fit_model,
        included=included,
    )
    run = {
        "RunId": manifest_row["RunId"],
        "Design": manifest_row["Design"],
        "ThresholdCondition": manifest_row["ThresholdCondition"],
        "Replicate": int(manifest_row["Replicate"]),
        "Rows": int(len(ratings)),
        "FitModel": fit_model,
        "EstimatorFamily": "MML",
        "EstimatorMode": mode,
        "LikelihoodBasis": "normal_population_marginal",
        "FitReturned": returned,
        "Converged": converged,
        "InferenceReady": ready,
        "IncludedInBridge": included,
        "FiniteMainEstimates": bool(finite_main),
        "FiniteThresholdEstimates": bool(finite_steps),
        "MaxAbsConstraintResidual": constraint["MaxAbsConstraintResidual"],
        "ConstraintPass": constraint["ConstraintPass"],
        "PopulationSDMode": "estimated" if MODE_SPECS[mode]["estimate_population_sd"] else "fixed",
        "PopulationSDInput": MODE_SPECS[mode]["population_prior_sd"],
        "EstimatedPopulationSD": summary.get("EstimatedPopulationSD"),
        "LogLik": summary.get("LogLik"),
        "LogLikPerObs": summary.get("LogLikPerObs"),
        "AIC": summary.get("AIC"),
        "BIC": summary.get("BIC"),
        "GradientNorm": summary.get("GradientNorm"),
        "Warnings": " | ".join(warning_messages),
        "FailureReason": "" if included else "MML operational readiness/constraint gate failed",
        "ElapsedSeconds": time.perf_counter() - started,
    }
    return run, recovery, thresholds, constraint


def fit_exact_cmle(
    ratings: pd.DataFrame,
    truth: pd.DataFrame,
    threshold_truth: pd.DataFrame,
    manifest_row: pd.Series,
    *,
    fit_model: str,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    started = time.perf_counter()
    kwargs = {
        "person_col": "Person",
        "facet_cols": ["Rater", "Task", "Criterion"],
        "score_col": "Score",
        "rating_min": 0,
        "rating_max": int(manifest_row["Categories"]) - 1,
        "model": fit_model,
        "step_facet": "Criterion" if fit_model == "PCM" else None,
    }
    audit = audit_cmle_eligibility(ratings, **kwargs)
    audit_eligible = bool(audit.get("eligible", False))
    if not audit_eligible:
        issues = audit.get("issues_table", pd.DataFrame())
        raise ValueError("CMLE prefit audit failed: " + ";".join(issues.get("Code", [])))
    result = fit_cmle(
        ratings,
        **kwargs,
        maxiter=800,
        gtol=1e-8,
        newton_polish_maxiter=12,
        finite_mle_gate=True,
    )
    summary = result["summary"].iloc[0]
    facets = result["facets"]["others"].copy()
    thresholds = _normalize_steps(
        result,
        fit_model=fit_model,
        estimator_mode=CMLE_MODE,
        run_id=str(manifest_row["RunId"]),
        threshold_condition=str(manifest_row["ThresholdCondition"]),
        threshold_truth=threshold_truth,
    )
    finite_main = (
        pd.to_numeric(facets["Estimate"], errors="coerce").notna().all()
        and pd.to_numeric(facets["SE"], errors="coerce").notna().all()
    )
    finite_steps = (
        pd.to_numeric(thresholds["Estimate"], errors="coerce").notna().all()
        and pd.to_numeric(thresholds["SE"], errors="coerce").notna().all()
    )
    constraint = _constraint_audit(
        facets, thresholds, estimator_mode=CMLE_MODE, fit_model=fit_model
    )
    converged = bool(summary["Converged"])
    ready = bool(summary["InferenceReady"])
    finite_mle = bool(summary["FiniteMLEExistenceQualified"])
    included = bool(
        audit_eligible and converged and ready and finite_mle
        and finite_main and finite_steps and constraint["ConstraintPass"]
    )
    recovery = _normalize_recovery(
        facets,
        truth,
        manifest_row,
        estimator_family="CMLE",
        estimator_mode=CMLE_MODE,
        fit_model=fit_model,
        included=included,
    )
    rows = len(ratings)
    run = {
        "RunId": manifest_row["RunId"],
        "Design": manifest_row["Design"],
        "ThresholdCondition": manifest_row["ThresholdCondition"],
        "Replicate": int(manifest_row["Replicate"]),
        "Rows": rows,
        "FitModel": fit_model,
        "EstimatorFamily": "CMLE",
        "EstimatorMode": CMLE_MODE,
        "LikelihoodBasis": "exact_person_total_conditional",
        "FitReturned": True,
        "ConditionalDesignEligible": audit_eligible,
        "ConditionalRank": audit.get("conditional_rank"),
        "ConditionalNullity": audit.get("conditional_nullity"),
        "PersonsInformative": audit.get("persons_informative"),
        "PersonsExtreme": audit.get("persons_extreme"),
        "Converged": converged,
        "InferenceReady": ready,
        "FiniteMLEExistenceQualified": finite_mle,
        "IncludedInBridge": included,
        "FiniteMainEstimates": bool(finite_main),
        "FiniteThresholdEstimates": bool(finite_steps),
        "MaxAbsConstraintResidual": constraint["MaxAbsConstraintResidual"],
        "ConstraintPass": constraint["ConstraintPass"],
        "LogLik": summary["ConditionalLogLik"],
        "LogLikPerObs": float(summary["ConditionalLogLik"]) / rows,
        "AIC": summary["ConditionalAIC"],
        "BIC": np.nan,
        "GradientNorm": summary["GradientNorm"],
        "FailureReason": "" if included else str(summary["ReadinessReasons"]),
        "ElapsedSeconds": time.perf_counter() - started,
    }
    return run, recovery, thresholds, constraint


def build_model_direction_checks(
    parent_ledger: pd.DataFrame,
    new_runs: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    parent = parent_ledger.copy()
    parent = parent.assign(
        EstimatorMode=PARENT_PYTHON_MODE,
        LikelihoodBasis="fixed_person_joint",
        IncludedInBridge=parent["ComparisonEligible"],
        LogLik=parent["PythonLogLik"],
        LogLikPerObs=parent["PythonLogLikPerObs"],
        AIC=parent["PythonAIC"],
        BIC=parent["PythonBIC"],
    )
    columns = [
        "RunId", "Design", "ThresholdCondition", "Replicate", "Rows",
        "FitModel", "EstimatorMode", "LikelihoodBasis", "IncludedInBridge",
        "LogLik", "LogLikPerObs", "AIC", "BIC",
    ]
    combined = pd.concat([parent[columns], new_runs[columns]], ignore_index=True)
    pair_rows = []
    for (mode, run_id), group in combined.groupby(["EstimatorMode", "RunId"]):
        if set(group["FitModel"]) != {"PCM", "RSM"} or not group["IncludedInBridge"].all():
            continue
        lookup = group.set_index("FitModel")
        pair_rows.append({
            "EstimatorMode": mode,
            "LikelihoodBasis": lookup.loc["PCM", "LikelihoodBasis"],
            "RunId": run_id,
            "Design": lookup.loc["PCM", "Design"],
            "ThresholdCondition": lookup.loc["PCM", "ThresholdCondition"],
            "Replicate": int(lookup.loc["PCM", "Replicate"]),
            "PCMMinusRSMLogLik": lookup.loc["PCM", "LogLik"] - lookup.loc["RSM", "LogLik"],
            "PCMMinusRSMLogLikPerObs": (
                lookup.loc["PCM", "LogLikPerObs"] - lookup.loc["RSM", "LogLikPerObs"]
            ),
            "RSMMinusPCMAIC": lookup.loc["RSM", "AIC"] - lookup.loc["PCM", "AIC"],
            "RSMMinusPCMBIC": (
                lookup.loc["RSM", "BIC"] - lookup.loc["PCM", "BIC"]
                if np.isfinite(lookup.loc["RSM", "BIC"]) and np.isfinite(lookup.loc["PCM", "BIC"])
                else np.nan
            ),
        })
    pairs = pd.DataFrame(pair_rows)
    checks = []
    for (mode, design, replicate), group in pairs.groupby(
        ["EstimatorMode", "Design", "Replicate"]
    ):
        lookup = group.set_index("ThresholdCondition")
        shared = lookup.loc["shared", "PCMMinusRSMLogLikPerObs"] if "shared" in lookup.index else np.nan
        hetero = lookup.loc["heterogeneous", "PCMMinusRSMLogLikPerObs"] if "heterogeneous" in lookup.index else np.nan
        checks.append({
            "EstimatorMode": mode,
            "LikelihoodBasis": group["LikelihoodBasis"].iloc[0],
            "Design": design,
            "Replicate": int(replicate),
            "SharedGainPerObs": shared,
            "HeterogeneousGainPerObs": hetero,
            "DirectionPass": bool(np.isfinite(shared) and np.isfinite(hetero) and hetero > shared),
        })
    return pairs, pd.DataFrame(checks)


def build_estimator_shifts(recovery: pd.DataFrame) -> pd.DataFrame:
    base = recovery[
        recovery["EstimatorMode"].eq(PARENT_PYTHON_MODE)
    ][["RunId", "FitModel", "Facet", "Level", "EstimateAligned"]].rename(
        columns={"EstimateAligned": "PythonJMLEEstimateAligned"}
    )
    output = recovery.merge(
        base, on=["RunId", "FitModel", "Facet", "Level"], how="left", validate="many_to_one"
    )
    output["ShiftFromPythonJMLE"] = (
        output["EstimateAligned"] - output["PythonJMLEEstimateAligned"]
    )
    output["ComparisonClass"] = np.where(
        output["EstimatorMode"].eq(PARENT_FACETS_MODE),
        "same_estimand_engine_comparison",
        np.where(
            output["EstimatorMode"].eq(PARENT_PYTHON_MODE),
            "identity_reference",
            "different_estimand_descriptive",
        ),
    )
    return output


def _rmse(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    return float(np.sqrt(np.mean(np.square(numeric)))) if len(numeric) else np.nan


def run_bridge(args: argparse.Namespace) -> None:
    parent_dir = args.parent_dir.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    input_hashes = validate_parent_identity(parent_dir)
    tables, manifest = selected_manifest(parent_dir)
    selected_ids = set(manifest["RunId"].astype(str))
    parent_ledger, parent_recovery, parent_thresholds = import_parent_jmle(
        parent_dir, selected_ids
    )
    output_dir.mkdir(parents=True)
    ratings_all = tables["generated_ratings.csv"]
    truth_all = tables["generated_facet_truth.csv"]
    threshold_truth_all = pd.read_csv(
        parent_dir / "retained_input" / "generated_pcm_threshold_truth.csv"
    )
    import streamlit_app as app  # pylint: disable=import-outside-toplevel

    run_rows = []
    recovery_parts = [parent_recovery]
    threshold_parts = [parent_thresholds]
    constraint_rows = []
    attempted = len(manifest) * len(FIT_MODELS) * (len(MML_MODES) + 1)
    attempt_number = 0
    for _, manifest_row in manifest.iterrows():
        run_id = str(manifest_row["RunId"])
        ratings = ratings_all[ratings_all["RunId"].astype(str).eq(run_id)].copy()
        truth = truth_all[truth_all["RunId"].astype(str).eq(run_id)].copy()
        threshold_truth = threshold_truth_all[
            threshold_truth_all["RunId"].astype(str).eq(run_id)
        ].copy()
        for fit_model in FIT_MODELS:
            for mode in MML_MODES:
                attempt_number += 1
                try:
                    run, recovery, thresholds, constraint = fit_mml(
                        app,
                        ratings,
                        truth,
                        threshold_truth,
                        manifest_row,
                        fit_model=fit_model,
                        mode=mode,
                    )
                    run_rows.append(run)
                    recovery_parts.append(recovery)
                    threshold_parts.append(thresholds.assign(
                        EstimatorFamily="MML",
                        EstimandClass="normal_person_population_marginal_likelihood",
                        Design=manifest_row["Design"],
                        Replicate=manifest_row["Replicate"],
                        Seed=manifest_row["Seed"],
                        IncludedInBridge=run["IncludedInBridge"],
                    ))
                    constraint_rows.append({
                        "RunId": run_id,
                        "FitModel": fit_model,
                        "EstimatorMode": mode,
                        **constraint,
                    })
                except Exception as exc:  # preserve full attempt ledger
                    run_rows.append({
                        "RunId": run_id, "Design": manifest_row["Design"],
                        "ThresholdCondition": manifest_row["ThresholdCondition"],
                        "Replicate": int(manifest_row["Replicate"]), "Rows": len(ratings),
                        "FitModel": fit_model, "EstimatorFamily": "MML",
                        "EstimatorMode": mode, "LikelihoodBasis": "normal_population_marginal",
                        "FitReturned": False, "Converged": False, "InferenceReady": False,
                        "IncludedInBridge": False,
                        "FailureReason": f"{type(exc).__name__}: {exc}",
                    })
                print(
                    f"[{attempt_number:02d}/{attempted:02d}] {run_id} {fit_model} {mode}: "
                    f"included={run_rows[-1].get('IncludedInBridge', False)}",
                    flush=True,
                )
            attempt_number += 1
            try:
                run, recovery, thresholds, constraint = fit_exact_cmle(
                    ratings,
                    truth,
                    threshold_truth,
                    manifest_row,
                    fit_model=fit_model,
                )
                run_rows.append(run)
                recovery_parts.append(recovery)
                threshold_parts.append(thresholds.assign(
                    EstimatorFamily="CMLE",
                    EstimandClass="person_total_conditional_likelihood",
                    Design=manifest_row["Design"],
                    Replicate=manifest_row["Replicate"],
                    Seed=manifest_row["Seed"],
                    IncludedInBridge=run["IncludedInBridge"],
                ))
                constraint_rows.append({
                    "RunId": run_id,
                    "FitModel": fit_model,
                    "EstimatorMode": CMLE_MODE,
                    **constraint,
                })
            except Exception as exc:
                run_rows.append({
                    "RunId": run_id, "Design": manifest_row["Design"],
                    "ThresholdCondition": manifest_row["ThresholdCondition"],
                    "Replicate": int(manifest_row["Replicate"]), "Rows": len(ratings),
                    "FitModel": fit_model, "EstimatorFamily": "CMLE",
                    "EstimatorMode": CMLE_MODE, "LikelihoodBasis": "exact_person_total_conditional",
                    "FitReturned": False, "Converged": False, "InferenceReady": False,
                    "IncludedInBridge": False,
                    "FailureReason": f"{type(exc).__name__}: {exc}",
                })
            print(
                f"[{attempt_number:02d}/{attempted:02d}] {run_id} {fit_model} {CMLE_MODE}: "
                f"included={run_rows[-1].get('IncludedInBridge', False)}",
                flush=True,
            )

    new_runs = pd.DataFrame(run_rows)
    recovery = pd.concat(recovery_parts, ignore_index=True)
    thresholds = pd.concat(threshold_parts, ignore_index=True)
    constraints = pd.DataFrame(constraint_rows)
    model_pairs, direction_checks = build_model_direction_checks(parent_ledger, new_runs)
    shifts = build_estimator_shifts(recovery)
    included_recovery = recovery[recovery["IncludedInBridge"]]
    recovery_summary = (
        included_recovery.groupby(
            ["EstimatorMode", "FitModel", "Design", "ThresholdCondition", "Facet"],
            dropna=False,
        )
        .agg(
            Parameters=("ErrorAligned", "size"),
            MeanError=("ErrorAligned", "mean"),
            RMSE=("ErrorAligned", _rmse),
            MAE=("ErrorAligned", lambda values: float(pd.to_numeric(values).abs().mean())),
        )
        .reset_index()
    )
    included_thresholds = thresholds[
        thresholds["IncludedInBridge"] & thresholds["ThresholdTruth"].notna()
    ]
    threshold_summary = (
        included_thresholds.groupby(
            ["EstimatorMode", "FitModel", "ThresholdCondition"], dropna=False
        )
        .agg(
            Thresholds=("TruthError", "size"),
            MeanError=("TruthError", "mean"),
            RMSE=("TruthError", _rmse),
            MAE=("TruthError", lambda values: float(pd.to_numeric(values).abs().mean())),
        )
        .reset_index()
    )

    mml_runs = new_runs[new_runs["EstimatorFamily"].eq("MML")]
    cmle_runs = new_runs[new_runs["EstimatorFamily"].eq("CMLE")]
    free_sd = mml_runs[mml_runs["EstimatorMode"].eq("PYTHON_MML_FREE_SD_Q31")]
    free_sd_valid = pd.to_numeric(free_sd["EstimatedPopulationSD"], errors="coerce")
    gates = {
        "input_identity": True,
        "parent_jmle_eligible": bool(len(parent_ledger) == 16 and parent_ledger["ComparisonEligible"].all()),
        "mml_32_operational": bool(
            len(mml_runs) == 32
            and mml_runs["FitReturned"].all()
            and mml_runs["Converged"].all()
            and mml_runs["InferenceReady"].all()
            and mml_runs["IncludedInBridge"].all()
            and free_sd_valid.notna().all()
            and free_sd_valid.gt(0).all()
        ),
        "cmle_16_operational": bool(
            len(cmle_runs) == 16
            and cmle_runs["FitReturned"].all()
            and cmle_runs["ConditionalDesignEligible"].all()
            and cmle_runs["ConditionalNullity"].eq(0).all()
            and cmle_runs["Converged"].all()
            and cmle_runs["InferenceReady"].all()
            and cmle_runs["FiniteMLEExistenceQualified"].all()
            and cmle_runs["IncludedInBridge"].all()
        ),
        "constraints": bool(
            len(constraints) == 48
            and constraints["ConstraintPass"].all()
            and pd.to_numeric(constraints["MaxAbsConstraintResidual"], errors="coerce")
            .le(CONSTRAINT_TOLERANCE).all()
        ),
        "model_direction": bool(
            len(direction_checks) == 16 and direction_checks["DirectionPass"].all()
        ),
    }
    metrics = {
        "schema_version": SCHEMA_VERSION,
        "claim_limit": "Two-replicate estimand bridge; no estimator ranking or operating-characteristic claim.",
        "selected_parent_runs": int(len(manifest)),
        "reused_parent_jmle_contracts": int(len(parent_ledger) * 2),
        "new_fit_attempts": int(len(new_runs)),
        "mml_attempts": int(len(mml_runs)),
        "mml_included": int(mml_runs["IncludedInBridge"].sum()),
        "cmle_attempts": int(len(cmle_runs)),
        "cmle_included": int(cmle_runs["IncludedInBridge"].sum()),
        "recovery_rows": int(len(recovery)),
        "threshold_rows": int(len(thresholds)),
        "direction_checks": int(len(direction_checks)),
        "direction_passes": int(direction_checks["DirectionPass"].sum()),
        "free_sd_min": float(free_sd_valid.min()),
        "free_sd_max": float(free_sd_valid.max()),
        "max_abs_constraint_residual": float(
            pd.to_numeric(constraints["MaxAbsConstraintResidual"], errors="coerce").max()
        ),
        "gates": gates,
        "qualification_pass": bool(all(gates.values())),
        "likelihood_comparison_contract": {
            "joint": PARENT_PYTHON_MODE,
            "marginal": list(MML_MODES),
            "conditional": CMLE_MODE,
            "cross_basis_numeric_comparison": "prohibited",
        },
    }
    new_runs.to_csv(output_dir / "estimand_bridge_run_ledger.csv", index=False)
    recovery.to_csv(output_dir / "estimand_bridge_recovery.csv", index=False)
    thresholds.to_csv(output_dir / "estimand_bridge_thresholds.csv", index=False)
    constraints.to_csv(output_dir / "estimand_bridge_constraints.csv", index=False)
    shifts.to_csv(output_dir / "estimand_bridge_shifts_from_jmle.csv", index=False)
    recovery_summary.to_csv(output_dir / "estimand_bridge_recovery_summary.csv", index=False)
    threshold_summary.to_csv(output_dir / "estimand_bridge_threshold_summary.csv", index=False)
    model_pairs.to_csv(output_dir / "estimand_bridge_model_pairs.csv", index=False)
    direction_checks.to_csv(output_dir / "estimand_bridge_direction_checks.csv", index=False)
    (output_dir / "estimand_bridge_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    identity = {
        "schema_version": SCHEMA_VERSION,
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "plan_sha256": sha256_file(REPO_ROOT / "validation" / "estimand_bridge_pilot_plan_20260811.json"),
        "app_sha256": sha256_file(REPO_ROOT / "streamlit_app.py"),
        "cmle_sha256": sha256_file(REPO_ROOT / "mfrm_app" / "cmle.py"),
        "parent_dir": str(parent_dir),
        "parent_metrics_sha256": sha256_file(parent_dir / "boundary_metrics.json"),
        "parent_identity_sha256": sha256_file(parent_dir / "boundary_identity.json"),
        "retained_input_sha256": input_hashes,
        "python_version": sys.version,
        "platform": platform.platform(),
    }
    (output_dir / "estimand_bridge_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main() -> None:
    run_bridge(parse_args())


if __name__ == "__main__":
    main()
