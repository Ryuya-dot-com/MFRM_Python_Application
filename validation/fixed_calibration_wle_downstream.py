#!/usr/bin/env python3
"""Replay rank-full Stage-B2 fits and audit fixed-calibration WLE sensitivity."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import streamlit_app as app  # noqa: E402
from mfrm_app.person_scoring import score_fixed_calibration_persons  # noqa: E402
from validation import operating_characteristics_pilot as pilot  # noqa: E402
from validation import operating_characteristics_strict_jmle as stage_a  # noqa: E402
from validation import operating_characteristics_strict_jmle_a2 as stage_a2  # noqa: E402


DEFAULT_INPUT = ROOT / "validation/operating_characteristics_pilot20_20260809"
DEFAULT_B2 = ROOT / "validation/operating_characteristics_strict_jmle_b2_20260809"
DEFAULT_IDENT = ROOT / "validation/operating_characteristics_identifiability_20260809"
DEFAULT_PLAN = ROOT / "validation/fixed_calibration_wle_downstream_plan_20260809.json"
DEFAULT_STAGE_A_PLAN = ROOT / "validation/operating_characteristics_strict_jmle_plan_20260809.json"
DEFAULT_OUTPUT = ROOT / "validation/fixed_calibration_wle_downstream_20260809"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_and_validate_plan(path: Path, b2_dir: Path, ident_dir: Path) -> dict:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != "mfrm-fixed-calibration-wle-downstream-plan-v1":
        raise ValueError("Unexpected downstream WLE plan schema.")
    expected = plan["input_identity"]
    paths = {
        "stage_b2_identity_sha256": b2_dir / "strict_jmle_b2_identity.json",
        "stage_b2_runs_sha256": b2_dir / "strict_jmle_b2_runs.csv",
        "identifiability_runs_sha256": ident_dir / "jmle_identifiability_runs.csv",
        "wle_parity_decision_sha256": ROOT / "validation/fixed_calibration_wle_20260809/parity_decision.json",
        "person_scoring_sha256": ROOT / "mfrm_app/person_scoring.py",
        "streamlit_app_sha256": ROOT / "streamlit_app.py",
        "stage_b2_adapter_sha256": ROOT / "validation/operating_characteristics_strict_jmle_b2.py",
        "stage_b2_plan_sha256": ROOT / "validation/operating_characteristics_strict_jmle_b2_plan_20260809.json",
    }
    failures = [key for key, file in paths.items() if sha256_file(file) != str(expected[key])]
    if failures:
        raise ValueError(f"Downstream WLE input identity failed: {failures}")
    parity = json.loads(paths["wle_parity_decision_sha256"].read_text(encoding="utf-8"))
    if parity.get("overall_passed") is not True:
        raise ValueError("Fixed-calibration Python-TAM WLE parity is not qualified.")
    return plan


def _patch_result(result: dict, params: dict, optimizer, objective: float, person_scoring: pd.DataFrame) -> dict:
    out = copy.deepcopy(result)
    out["params"] = params
    out["opt"] = optimizer
    if isinstance(out.get("summary"), pd.DataFrame) and not out["summary"].empty:
        out["summary"].loc[:, "LogLik"] = -float(objective)
        out["summary"].loc[:, "Converged"] = bool(getattr(optimizer, "success", False))

    person_table = out["facets"]["person"].copy()
    scoring = person_scoring.copy()
    scoring["Person"] = scoring["Person"].astype(str)
    person_table["Person"] = person_table["Person"].astype(str)
    estimate_map = scoring.set_index("Person")["Estimate"]
    se_map = scoring.set_index("Person")["StandardError"]
    person_table["Estimate"] = person_table["Person"].map(estimate_map).astype(float)
    if "SD" in person_table.columns:
        person_table["SD"] = person_table["Person"].map(se_map).astype(float)
    for column in (
        "ExtremeScorePattern",
        "ExtremeScoreDirection",
        "EstimateRole",
    ):
        if column in scoring.columns:
            person_table[column] = person_table["Person"].map(scoring.set_index("Person")[column])
    is_wle = person_scoring["Estimator"].astype(str).eq("WLE_fixed_calibration").all()
    if is_wle:
        person_table["FiniteJMLEEstimate"] = False
        person_table["PersonInferenceReady"] = False
        person_table["ReportableEstimate"] = np.nan
    else:
        interior = ~person_table["ExtremeScorePattern"].fillna(False).astype(bool)
        person_table["FiniteJMLEEstimate"] = interior
        person_table["PersonInferenceReady"] = interior
        person_table["ReportableEstimate"] = person_table["Estimate"].where(interior)
    out["facets"]["person"] = person_table

    facet_table = out["facets"]["others"].copy()
    for facet, values in params["facets"].items():
        levels = [str(value) for value in out["config"]["facet_levels"][facet]]
        lookup = dict(zip(levels, np.asarray(values, dtype=float)))
        selected = facet_table["Facet"].astype(str).eq(str(facet))
        facet_table.loc[selected, "Estimate"] = (
            facet_table.loc[selected, "Level"].astype(str).map(lookup).to_numpy()
        )
    out["facets"]["others"] = facet_table
    return out


def _person_score_from_polished(result: dict, params: dict, idx: dict) -> pd.DataFrame:
    if result["config"]["model"] != "RSM":
        raise ValueError("The frozen Stage-B2 bundle is expected to contain RSM fits only.")
    k = np.arange(int(result["config"]["n_cat"]), dtype=float)
    base_eta = app.compute_base_eta(idx, params, result["config"])
    step_cumulative = np.concatenate([[0.0], np.cumsum(params["steps"])])
    intercepts = base_eta[:, None] * k[None, :] - step_cumulative[None, :]
    slopes = np.tile(k, (len(base_eta), 1))
    person_levels = [str(value) for value in result["prep"]["levels"]["Person"]]
    persons = np.asarray(person_levels, dtype=object)[idx["person"]]
    return score_fixed_calibration_persons(
        persons,
        idx["score_k"],
        intercepts,
        slopes,
        row_weights=idx.get("weight"),
        person_levels=person_levels,
    )


def _fit_zone(value: float) -> str:
    if not np.isfinite(value):
        return "unavailable"
    if value < 0.5:
        return "below_0.5"
    if value <= 1.5:
        return "0.5_to_1.5"
    if value < 2.0:
        return "above_1.5_below_2.0"
    return "at_or_above_2.0"


def _fit_comparison(run_meta: dict, diagnostics_jmle: dict, diagnostics_wle: dict) -> pd.DataFrame:
    keys = ["Facet", "Level"]
    columns = keys + ["Infit", "Outfit"]
    left = diagnostics_jmle["measures"][columns].rename(
        columns={"Infit": "InfitJMLE", "Outfit": "OutfitJMLE"}
    )
    right = diagnostics_wle["measures"][columns].rename(
        columns={"Infit": "InfitWLE", "Outfit": "OutfitWLE"}
    )
    merged = left.merge(right, on=keys, how="outer", validate="one_to_one")
    rows: list[dict[str, object]] = []
    for record in merged.itertuples(index=False):
        for metric in ("Infit", "Outfit"):
            jmle = float(getattr(record, f"{metric}JMLE"))
            wle = float(getattr(record, f"{metric}WLE"))
            jmle_display = float(format(jmle, ".3g")) if np.isfinite(jmle) else np.nan
            wle_display = float(format(wle, ".3g")) if np.isfinite(wle) else np.nan
            rows.append(
                {
                    **run_meta,
                    "Facet": str(record.Facet),
                    "Level": str(record.Level),
                    "Metric": metric,
                    "JMLE": jmle,
                    "WLE": wle,
                    "WLEMinusJMLE": wle - jmle,
                    "ZoneJMLE": _fit_zone(jmle),
                    "ZoneWLE": _fit_zone(wle),
                    "RawEstimatorDecisionSwitch": _fit_zone(jmle) != _fit_zone(wle),
                    "DisplayJMLE": jmle_display,
                    "DisplayWLE": wle_display,
                    "RawDisplayMismatchJMLE": _fit_zone(jmle) != _fit_zone(jmle_display),
                    "RawDisplayMismatchWLE": _fit_zone(wle) != _fit_zone(wle_display),
                }
            )
    return pd.DataFrame(rows)


def _focal_bias(result: dict, diagnostics: dict) -> dict[str, object]:
    bundle = app.estimate_bias_interaction(
        result,
        diagnostics,
        "Rater",
        "Task",
        omit_extreme=False,
    )
    if not isinstance(bundle, dict) or "table" not in bundle:
        return {"Available": False, "Reason": str(bundle.get("_skip_reason", "unavailable"))}
    dff = app.build_dff_bias_screening_table(
        bundle,
        alpha=0.05,
        min_n=5,
        practical_logit=0.50,
    )
    focal = dff.loc[
        dff["FacetA_Level"].astype(str).eq(pilot.FOCAL_RATER)
        & dff["FacetB_Level"].astype(str).eq(pilot.FOCAL_TASK)
    ]
    if focal.empty:
        return {"Available": False, "Reason": "focal cell absent"}
    row = focal.iloc[0]
    p_holm = float(row["p_holm"])
    abs_bias = float(row["AbsBias"])
    sparse = bool(row["SparseCell"])
    holm = bool(np.isfinite(p_holm) and p_holm < 0.05)
    practical = bool(np.isfinite(abs_bias) and abs_bias >= 0.50)
    return {
        "Available": bool(np.isfinite(p_holm) and np.isfinite(abs_bias)),
        "Reason": "",
        "BiasEstimate": float(row["BiasSize"]),
        "BiasSE": float(row["SE"]),
        "PHolm": p_holm,
        "AbsBias": abs_bias,
        "SparseCell": sparse,
        "DecisionHolmRaw": holm,
        "DecisionPracticalRaw": practical,
        "DecisionStrongRaw": bool(holm and practical and not sparse),
        "DecisionHolmDisplayed": bool(np.isfinite(p_holm) and round(p_holm, 4) < 0.05),
        "DecisionPracticalDisplayed": bool(np.isfinite(abs_bias) and round(abs_bias, 4) >= 0.50),
    }


def _plot_persons(persons: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3))
    normal = ~persons["ExtremeScorePattern"].astype(bool)
    axes[0].scatter(
        persons.loc[normal, "JMLE"], persons.loc[normal, "WLE"],
        s=12, alpha=0.35, color="#1f77b4", label="Interior",
    )
    extreme = ~normal
    axes[0].scatter(
        persons.loc[extreme, "JMLE"], persons.loc[extreme, "WLE"],
        s=38, alpha=0.9, color="#d62728", label="Exact extreme",
    )
    finite = persons[["JMLE", "WLE"]].to_numpy(dtype=float)
    low, high = float(np.nanmin(finite)), float(np.nanmax(finite))
    axes[0].plot([low, high], [low, high], color="black", linewidth=1, linestyle="--")
    axes[0].set_xlabel("Polished JMLE Person coordinate")
    axes[0].set_ylabel("Fixed-calibration Warm WLE")
    axes[0].legend(frameon=False)
    axes[0].set_title("Person-coordinate sensitivity")

    interior_persons = persons.loc[normal]
    central_limit = max(0.1, float(interior_persons["WLEMinusJMLE"].abs().quantile(0.995)))
    bins = np.linspace(-central_limit, central_limit, 41)
    for design, frame in interior_persons.groupby("Design", sort=False):
        axes[1].hist(
            frame["WLEMinusJMLE"], bins=bins, alpha=0.55, label=str(design)
        )
    axes[1].axvline(0.0, color="black", linewidth=1, linestyle="--")
    axes[1].set_xlabel("WLE minus JMLE, interior (logits)")
    axes[1].set_ylabel("Person count")
    axes[1].legend(frameon=False)
    axes[1].set_title("Interior shift distribution")

    extreme_frame = persons.loc[extreme].copy()
    if not extreme_frame.empty:
        extreme_frame["Label"] = (
            extreme_frame["Design"].astype(str).str.replace("balanced_large_anchors", "large_anchor")
            + "\n"
            + extreme_frame["TruthBias"].map(lambda value: f"bias={value:g}")
        )
        x = np.arange(len(extreme_frame))
        axes[2].stem(
            x,
            extreme_frame["WLEMinusJMLE"],
            linefmt="#d62728",
            markerfmt="o",
            basefmt="black",
        )
        axes[2].set_xticks(x, extreme_frame["Label"], rotation=25, ha="right")
    axes[2].set_ylabel("WLE minus JMLE (logits)")
    axes[2].set_title("Exact-extreme shifts")
    fig.tight_layout()
    fig.savefig(output / "wle_person_sensitivity.png", dpi=180)
    plt.close(fig)


def _plot_decisions(fit: pd.DataFrame, bias: pd.DataFrame, output: Path) -> None:
    fit_counts = (
        fit.groupby("ConditionId", sort=False)["RawEstimatorDecisionSwitch"]
        .sum()
        .astype(int)
    )
    bias_counts = (
        bias.groupby("ConditionId", sort=False)["AnyRawBiasDecisionSwitch"]
        .sum()
        .astype(int)
    )
    labels = list(dict.fromkeys([*fit_counts.index.tolist(), *bias_counts.index.tolist()]))
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.bar(x - 0.18, [fit_counts.get(label, 0) for label in labels], width=0.36, label="Fit-zone switches")
    ax.bar(x + 0.18, [bias_counts.get(label, 0) for label in labels], width=0.36, label="Any focal bias-rule switches")
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.set_ylabel("Raw decision switches")
    ax.set_title("Downstream decision sensitivity by frozen condition")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "wle_downstream_decision_switches.png", dpi=180)
    plt.close(fig)


def run(input_dir: Path, b2_dir: Path, ident_dir: Path, stage_a_plan_path: Path, plan_path: Path, output: Path) -> None:
    started = time.perf_counter()
    input_dir = input_dir.resolve()
    b2_dir = b2_dir.resolve()
    ident_dir = ident_dir.resolve()
    output = output.resolve()
    plan = load_and_validate_plan(plan_path.resolve(), b2_dir, ident_dir)
    strict_plan = stage_a.load_plan(stage_a_plan_path.resolve())
    manifest = pd.read_csv(input_dir / "manifest.csv")
    ratings = pd.read_csv(input_dir / "generated_ratings.csv")
    truth = pd.read_csv(input_dir / "generated_facet_truth.csv")
    anchors = pd.read_csv(input_dir / "generated_anchors.csv")
    retained = pd.read_csv(b2_dir / "strict_jmle_b2_runs.csv").set_index("RunId")
    ident = pd.read_csv(ident_dir / "jmle_identifiability_runs.csv").set_index("RunId")
    included = manifest["RunId"].map(ident["StructurallyIdentified"]).fillna(False).astype(bool)
    if int(included.sum()) != int(plan["analysis_population"]["included_rank_full_runs"]):
        raise ValueError("The rank-full run count differs from the frozen plan.")

    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_wle_downstream_plan.json").write_bytes(plan_path.read_bytes())
    person_parts: list[pd.DataFrame] = []
    fit_parts: list[pd.DataFrame] = []
    bias_rows: list[dict[str, object]] = []
    run_rows: list[dict[str, object]] = []
    options = dict(plan["replay_controls"]["polish_options"])
    for manifest_row in manifest.loc[included].itertuples(index=False):
        row = pd.Series(manifest_row._asdict())
        run_id = str(row["RunId"])
        generated = stage_a.generated_design_for_run(row, ratings, truth, anchors)
        base, _, _ = stage_a2.fit_stage_a_baseline(
            row,
            generated,
            strict_plan["numerical_controls"],
        )
        start, idx, config, sizes = stage_a2.optimizer_context(base)
        optimizer = minimize(
            app.mfrm_loglik_jmle_value_grad,
            start,
            args=(idx, config, sizes),
            jac=True,
            method="L-BFGS-B",
            bounds=app.build_optimizer_bounds(sizes, config),
            options=options,
        )
        objective, gradient = app.mfrm_loglik_jmle_value_grad(
            optimizer.x, idx, config, sizes
        )
        gradient_sup = float(np.max(np.abs(gradient)))
        retained_row = retained.loc[run_id]
        objective_difference = float(objective - retained_row["Objective"])
        gradient_difference = float(
            gradient_sup - retained_row["TerminalGradientSupNorm"]
        )
        replay_passed = bool(
            abs(objective_difference)
            <= float(plan["replay_controls"]["objective_max_abs_difference_from_retained_b2"])
            and abs(gradient_difference)
            <= float(plan["replay_controls"]["gradient_sup_norm_max_abs_difference_from_retained_b2"])
        )
        params_jmle = app.expand_params(optimizer.x, sizes, config)
        scored = _person_score_from_polished(base, params_jmle, idx)
        params_wle = copy.deepcopy(params_jmle)
        params_wle["theta"] = scored["Estimate"].to_numpy(dtype=float)
        calibration_change = max(
            [0.0]
            + [
                float(np.max(np.abs(np.asarray(params_wle["facets"][facet]) - np.asarray(params_jmle["facets"][facet]))))
                for facet in config["facet_names"]
            ]
            + [float(np.max(np.abs(np.asarray(params_wle["steps"]) - np.asarray(params_jmle["steps"]))))]
        )

        jmle_scoring = scored.copy()
        jmle_scoring["Estimator"] = "JMLE_polished"
        jmle_scoring["Estimate"] = params_jmle["theta"]
        jmle_scoring["StandardError"] = np.nan
        jmle_scoring["EstimateRole"] = "polished_jmle_person_coordinate"
        result_jmle = _patch_result(base, params_jmle, optimizer, objective, jmle_scoring)
        result_wle = _patch_result(base, params_wle, optimizer, objective, scored)
        diagnostics_jmle = app.mfrm_diagnostics(result_jmle, compute_pca=False, compute_marginal=False)
        diagnostics_wle = app.mfrm_diagnostics(result_wle, compute_pca=False, compute_marginal=False)
        meta = {
            "RunId": run_id,
            "ConditionId": str(row["ConditionId"]),
            "Design": str(row["Design"]),
            "TruthBias": float(row["TruthBias"]),
            "Replicate": int(row["Replicate"]),
        }
        persons = pd.DataFrame(
            {
                **{key: value for key, value in meta.items()},
                "Person": scored["Person"].astype(str),
                "JMLE": params_jmle["theta"],
                "WLE": scored["Estimate"].to_numpy(dtype=float),
                "WLEStandardError": scored["StandardError"].to_numpy(dtype=float),
                "WLEMinusJMLE": scored["Estimate"].to_numpy(dtype=float) - params_jmle["theta"],
                "ExtremeScorePattern": scored["ExtremeScorePattern"].astype(bool),
                "ExtremeScoreDirection": scored["ExtremeScoreDirection"].astype(str),
                "WLEStatus": scored["Status"].astype(str),
                "WLEAdjustedScoreResidual": scored["AdjustedScoreResidual"].to_numpy(dtype=float),
            }
        )
        person_parts.append(persons)
        fit_parts.append(_fit_comparison(meta, diagnostics_jmle, diagnostics_wle))
        bias_jmle = _focal_bias(result_jmle, diagnostics_jmle)
        bias_wle = _focal_bias(result_wle, diagnostics_wle)
        bias_row: dict[str, object] = {**meta}
        for prefix, values in (("JMLE", bias_jmle), ("WLE", bias_wle)):
            for key, value in values.items():
                bias_row[f"{prefix}{key}"] = value
        bias_row["BiasEstimateWLEMinusJMLE"] = (
            float(bias_wle.get("BiasEstimate", np.nan))
            - float(bias_jmle.get("BiasEstimate", np.nan))
        )
        bias_row["DecisionStrongSwitch"] = (
            bias_jmle.get("DecisionStrongRaw") != bias_wle.get("DecisionStrongRaw")
        )
        bias_row["DecisionHolmSwitch"] = (
            bias_jmle.get("DecisionHolmRaw") != bias_wle.get("DecisionHolmRaw")
        )
        bias_row["DecisionPracticalSwitch"] = (
            bias_jmle.get("DecisionPracticalRaw") != bias_wle.get("DecisionPracticalRaw")
        )
        bias_row["AnyRawBiasDecisionSwitch"] = bool(
            bias_row["DecisionHolmSwitch"]
            or bias_row["DecisionPracticalSwitch"]
            or bias_row["DecisionStrongSwitch"]
        )
        bias_row["RawDisplayMismatchJMLE"] = (
            bias_jmle.get("DecisionStrongRaw")
            != bool(
                bias_jmle.get("DecisionHolmDisplayed", False)
                and bias_jmle.get("DecisionPracticalDisplayed", False)
                and not bias_jmle.get("SparseCell", True)
            )
        )
        bias_row["RawDisplayMismatchWLE"] = (
            bias_wle.get("DecisionStrongRaw")
            != bool(
                bias_wle.get("DecisionHolmDisplayed", False)
                and bias_wle.get("DecisionPracticalDisplayed", False)
                and not bias_wle.get("SparseCell", True)
            )
        )
        bias_rows.append(bias_row)
        run_rows.append(
            {
                **meta,
                "ReplayOptimizerSuccess": bool(optimizer.success),
                "ReplayObjectiveDifference": objective_difference,
                "ReplayGradientSupNormDifference": gradient_difference,
                "ReplayPassed": replay_passed,
                "Persons": len(scored),
                "ExtremePersons": int(scored["ExtremeScorePattern"].sum()),
                "WLEStatusOK": int(scored["Status"].eq("ok").sum()),
                "WLEMaxAdjustedScoreResidual": float(scored["AdjustedScoreResidual"].max()),
                "CalibrationMaxAbsChange": calibration_change,
            }
        )

    persons = pd.concat(person_parts, ignore_index=True)
    fit = pd.concat(fit_parts, ignore_index=True)
    bias = pd.DataFrame(bias_rows)
    runs = pd.DataFrame(run_rows)
    withheld = manifest.loc[~included, ["RunId", "ConditionId", "Design", "TruthBias", "Replicate"]].copy()
    withheld["StructuralNullity"] = withheld["RunId"].map(ident["StructuralNullity"])
    withheld["Status"] = "withheld_structurally_rank_deficient"
    withheld["Reason"] = "Fixed-calibration WLE cannot repair unidentified facet calibration."

    summary_rows = []
    for design, frame in persons.groupby("Design", sort=False):
        run_ids = frame["RunId"].unique()
        fit_frame = fit.loc[fit["RunId"].isin(run_ids)]
        bias_frame = bias.loc[bias["RunId"].isin(run_ids)]
        summary_rows.append(
            {
                "Design": design,
                "Runs": len(run_ids),
                "Persons": len(frame),
                "ExtremePersons": int(frame["ExtremeScorePattern"].sum()),
                "MedianAbsPersonShift": float(frame["WLEMinusJMLE"].abs().median()),
                "P95AbsPersonShift": float(frame["WLEMinusJMLE"].abs().quantile(0.95)),
                "MaxAbsPersonShift": float(frame["WLEMinusJMLE"].abs().max()),
                "FitComparisons": len(fit_frame),
                "FitZoneSwitches": int(fit_frame["RawEstimatorDecisionSwitch"].sum()),
                "BiasHolmDecisionSwitches": int(bias_frame["DecisionHolmSwitch"].sum()),
                "BiasPracticalDecisionSwitches": int(bias_frame["DecisionPracticalSwitch"].sum()),
                "BiasStrongDecisionSwitches": int(bias_frame["DecisionStrongSwitch"].sum()),
                "MaxAbsFocalBiasShift": float(bias_frame["BiasEstimateWLEMinusJMLE"].abs().max()),
            }
        )
    summary = pd.DataFrame(summary_rows)
    fit_transitions = (
        fit.groupby(
            ["Facet", "Metric", "ZoneJMLE", "ZoneWLE"],
            dropna=False,
            sort=False,
        )
        .size()
        .reset_index(name="Comparisons")
    )
    fit_transitions["RawEstimatorDecisionSwitch"] = (
        fit_transitions["ZoneJMLE"] != fit_transitions["ZoneWLE"]
    )
    float_display = pd.DataFrame(
        [
            {
                "Surface": "JMLE fit",
                "Format": ".3g",
                "RawDisplayDecisionMismatches": int(fit["RawDisplayMismatchJMLE"].sum()),
            },
            {
                "Surface": "WLE fit",
                "Format": ".3g",
                "RawDisplayDecisionMismatches": int(fit["RawDisplayMismatchWLE"].sum()),
            },
            {
                "Surface": "JMLE focal bias",
                "Format": "4 decimals",
                "RawDisplayDecisionMismatches": int(bias["RawDisplayMismatchJMLE"].sum()),
            },
            {
                "Surface": "WLE focal bias",
                "Format": "4 decimals",
                "RawDisplayDecisionMismatches": int(bias["RawDisplayMismatchWLE"].sum()),
            },
        ]
    )
    bias_transitions = pd.DataFrame(
        [
            {
                "DecisionRule": "Holm p < 0.05",
                "JMLEPositive": int(bias["JMLEDecisionHolmRaw"].sum()),
                "WLEPositive": int(bias["WLEDecisionHolmRaw"].sum()),
                "Switches": int(bias["DecisionHolmSwitch"].sum()),
            },
            {
                "DecisionRule": "Abs bias >= 0.50",
                "JMLEPositive": int(bias["JMLEDecisionPracticalRaw"].sum()),
                "WLEPositive": int(bias["WLEDecisionPracticalRaw"].sum()),
                "Switches": int(bias["DecisionPracticalSwitch"].sum()),
            },
            {
                "DecisionRule": "Strong = Holm and practical and non-sparse",
                "JMLEPositive": int(bias["JMLEDecisionStrongRaw"].sum()),
                "WLEPositive": int(bias["WLEDecisionStrongRaw"].sum()),
                "Switches": int(bias["DecisionStrongSwitch"].sum()),
            },
        ]
    )
    all_gates_passed = bool(
        runs["ReplayPassed"].all()
        and runs["ReplayOptimizerSuccess"].all()
        and persons["WLEStatus"].eq("ok").all()
        and persons["WLEAdjustedScoreResidual"].le(1e-8).all()
        and np.isfinite(persons.loc[persons["ExtremeScorePattern"], "WLE"]).all()
        and runs["CalibrationMaxAbsChange"].eq(0.0).all()
        and len(withheld) == 40
    )

    runs.to_csv(output / "wle_downstream_runs.csv", index=False, float_format="%.17g")
    persons.to_csv(output / "wle_person_pairs.csv", index=False, float_format="%.17g")
    fit.to_csv(output / "wle_fit_decision_sensitivity.csv", index=False, float_format="%.17g")
    bias.to_csv(output / "wle_bias_decision_sensitivity.csv", index=False, float_format="%.17g")
    withheld.to_csv(output / "wle_structurally_withheld_runs.csv", index=False, float_format="%.17g")
    summary.to_csv(output / "wle_downstream_summary.csv", index=False, float_format="%.17g")
    fit_transitions.to_csv(output / "wle_fit_transition_summary.csv", index=False, float_format="%.17g")
    float_display.to_csv(output / "wle_float_display_summary.csv", index=False, float_format="%.17g")
    bias_transitions.to_csv(output / "wle_bias_transition_summary.csv", index=False, float_format="%.17g")
    _plot_persons(persons, output)
    _plot_decisions(fit, bias, output)

    decision = {
        "schema_version": "mfrm-fixed-calibration-wle-downstream-result-v1",
        "execution_gates_passed": all_gates_passed,
        "application_integration_authorized": False,
        "included_rank_full_runs": len(runs),
        "withheld_rank_deficient_runs": len(withheld),
        "persons": len(persons),
        "extreme_persons": int(persons["ExtremeScorePattern"].sum()),
        "maximum_absolute_person_shift": float(persons["WLEMinusJMLE"].abs().max()),
        "fit_zone_switches": int(fit["RawEstimatorDecisionSwitch"].sum()),
        "fit_raw_display_mismatches": int(
            fit["RawDisplayMismatchJMLE"].sum() + fit["RawDisplayMismatchWLE"].sum()
        ),
        "focal_bias_strong_decision_switches": int(bias["DecisionStrongSwitch"].sum()),
        "focal_bias_holm_decision_switches": int(bias["DecisionHolmSwitch"].sum()),
        "focal_bias_practical_decision_switches": int(bias["DecisionPracticalSwitch"].sum()),
        "focal_bias_any_rule_switches": int(bias["AnyRawBiasDecisionSwitch"].sum()),
        "bias_raw_display_mismatches": int(
            bias["RawDisplayMismatchJMLE"].sum() + bias["RawDisplayMismatchWLE"].sum()
        ),
        "maximum_absolute_focal_bias_shift": float(
            bias["BiasEstimateWLEMinusJMLE"].abs().max()
        ),
        "decision_source": "unrounded retained values",
        "elapsed_seconds": time.perf_counter() - started,
    }
    (output / "wle_downstream_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    report = f"""# Fixed-calibration Warm WLE downstream sensitivity

## Decision

Execution gates: **{'PASS' if all_gates_passed else 'FAIL'}**. Application integration remains **WITHHELD** by the frozen plan. The comparison holds polished facets and steps fixed and changes only the Person coordinate from polished JMLE to Warm WLE.

## Outcome

- Rank-full runs replayed: {len(runs)}/120; structurally rank-deficient sparse runs withheld: {len(withheld)}/40.
- Persons compared: {len(persons)}; exact extreme patterns: {int(persons['ExtremeScorePattern'].sum())}.
- Maximum absolute Person-coordinate shift: {persons['WLEMinusJMLE'].abs().max():.9g} logits.
- Fit-zone switches on raw values: {int(fit['RawEstimatorDecisionSwitch'].sum())}/{len(fit)} element-metric comparisons.
- Fit raw/display mismatches at `.3g`: {int(fit['RawDisplayMismatchJMLE'].sum() + fit['RawDisplayMismatchWLE'].sum())} across both estimator surfaces.
- Focal Holm decision switches: {int(bias['DecisionHolmSwitch'].sum())}/{len(bias)} runs.
- Focal practical-size decision switches: {int(bias['DecisionPracticalSwitch'].sum())}/{len(bias)} runs.
- Focal strong-bias decision switches: {int(bias['DecisionStrongSwitch'].sum())}/{len(bias)} runs.
- Maximum absolute focal bias-estimate shift: {bias['BiasEstimateWLEMinusJMLE'].abs().max():.9g} logits.
- Bias raw/display mismatches at four decimals: {int(bias['RawDisplayMismatchJMLE'].sum() + bias['RawDisplayMismatchWLE'].sum())} across both estimator surfaces.

## By design

{summary.to_string(index=False)}

## Interpretation boundary

Warm WLE removes the infinite fixed-calibration Person-score boundary, but it does not repair structurally unidentified facet calibration and does not propagate facet-calibration uncertainty. Any fit or bias decision switch is therefore evidence that Person-score convention is part of the downstream estimand, not evidence that WLE is automatically preferable. Decisions above use unrounded retained values; display-rounding mismatches are reported separately. No Streamlit estimator option or automatic fallback is enabled.
"""
    (output / "FIXED_CALIBRATION_WLE_DOWNSTREAM_RESULTS.md").write_text(report, encoding="utf-8")
    if not all_gates_passed:
        raise SystemExit("One or more frozen downstream WLE execution gates failed.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--b2", type=Path, default=DEFAULT_B2)
    parser.add_argument("--identifiability", type=Path, default=DEFAULT_IDENT)
    parser.add_argument("--stage-a-plan", type=Path, default=DEFAULT_STAGE_A_PLAN)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.input, args.b2, args.identifiability, args.stage_a_plan, args.plan, args.output)


if __name__ == "__main__":
    main()
