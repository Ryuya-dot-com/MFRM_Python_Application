#!/usr/bin/env python3
"""Run prospective Stage-A strict Python JMLE numerical qualification."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
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
from validation import operating_characteristics_pilot as pilot  # noqa: E402


PLAN_SCHEMA = "mfrm-strict-jmle-qualification-plan-v1"
DEFAULT_INPUT = REPO / "validation" / "operating_characteristics_20260809"
DEFAULT_OUTPUT = REPO / "validation" / "operating_characteristics_strict_jmle_smoke_20260809"
DEFAULT_PLAN = REPO / "validation" / "operating_characteristics_strict_jmle_plan_20260809.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_plan(path: Path) -> dict:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != PLAN_SCHEMA:
        raise ValueError(f"strict JMLE plan schema must be {PLAN_SCHEMA}")
    stage = plan.get("stage_a", {})
    controls = plan.get("numerical_controls", {})
    verification = plan.get("gradient_verification", {})
    gates = plan.get("stage_a_gates", {})
    expected = {
        "input_replicates_per_condition": 2,
        "conditions": 8,
        "attempted_fits": 16,
    }
    for key, value in expected.items():
        if int(stage.get(key, -1)) != value:
            raise ValueError(f"strict JMLE stage_a {key} must remain {value}")
    if int(controls.get("maxit", 0)) != 500:
        raise ValueError("strict JMLE maxit must remain 500")
    if float(controls.get("reltol", np.nan)) != 1e-9:
        raise ValueError("strict JMLE reltol must remain 1e-9")
    if float(controls.get("terminal_gradient_sup_norm_threshold", np.nan)) != 1e-4:
        raise ValueError("strict JMLE terminal-gradient threshold must remain 1e-4")
    if int(verification.get("finite_difference_coordinates_per_fit", 0)) != 8:
        raise ValueError("strict JMLE finite-difference coordinate count must remain 8")
    if float(verification.get("finite_difference_max_abs_error_tolerance", np.nan)) != 1e-5:
        raise ValueError("strict JMLE finite-difference tolerance must remain 1e-5")
    if not bool(gates.get("all_gates_must_pass_for_stage_b", False)):
        raise ValueError("strict JMLE Stage B must require every Stage-A gate")
    if plan.get("scope_boundaries", {}).get("public_surface_status") != "Withheld":
        raise ValueError("strict JMLE public surface must remain Withheld")
    return plan


def validate_input(input_dir: Path, plan: dict) -> pd.DataFrame:
    stage = plan["stage_a"]
    expected_files = {
        "manifest.csv": stage["input_manifest_file_sha256"],
        "generated_ratings.csv": stage["generated_ratings_sha256"],
        "generated_facet_truth.csv": stage["generated_facet_truth_sha256"],
        "generated_anchors.csv": stage["generated_anchors_sha256"],
    }
    rows: list[dict[str, object]] = []
    for filename, expected in expected_files.items():
        actual = sha256_file(input_dir / filename)
        rows.append({
            "Check": f"input_sha256::{filename}",
            "Passed": actual == str(expected),
            "Evidence": f"actual={actual}; expected={expected}",
        })
    manifest = pd.read_csv(input_dir / "manifest.csv")
    actual_frame = oc.frame_fingerprint(manifest)
    rows.append({
        "Check": "manifest_frame_sha256",
        "Passed": actual_frame == stage["input_manifest_frame_sha256"],
        "Evidence": f"actual={actual_frame}; expected={stage['input_manifest_frame_sha256']}",
    })
    rows.append({
        "Check": "manifest_shape",
        "Passed": (
            len(manifest) == int(stage["attempted_fits"])
            and manifest["ConditionId"].nunique() == int(stage["conditions"])
            and manifest.groupby("ConditionId")["Replicate"].nunique().eq(
                int(stage["input_replicates_per_condition"])
            ).all()
        ),
        "Evidence": (
            f"runs={len(manifest)}; conditions={manifest['ConditionId'].nunique()}; "
            f"replicate_range={manifest['Replicate'].min()}-{manifest['Replicate'].max()}"
        ),
    })
    bridge = pilot.validate_generated_bundle(input_dir)
    rows.append({
        "Check": "python_generated_bundle_validation",
        "Passed": bool(bridge["Passed"].all()),
        "Evidence": f"passed={int(bridge['Passed'].sum())}/{len(bridge)}",
    })
    checks = pd.DataFrame(rows)
    if not checks["Passed"].all():
        failed = checks.loc[~checks["Passed"], "Check"].astype(str).tolist()
        raise ValueError(f"strict JMLE input validation failed: {failed}")
    return checks


def generated_design_for_run(
    manifest_row: pd.Series,
    ratings: pd.DataFrame,
    truth: pd.DataFrame,
    anchors: pd.DataFrame,
) -> pilot.GeneratedDesign:
    run_id = str(manifest_row["RunId"])
    data = ratings.loc[
        ratings["RunId"].astype(str).eq(run_id),
        ["Person", "Rater", "Task", "Criterion", "Score"],
    ].copy()
    facet_truth = truth.loc[
        truth["RunId"].astype(str).eq(run_id),
        ["Facet", "Level", "Truth"],
    ].copy()
    run_anchors = anchors.loc[
        anchors["RunId"].astype(str).eq(run_id),
        ["Facet", "Level", "Anchor"],
    ].copy()
    full_rows = (
        int(manifest_row["Persons"])
        * int(manifest_row["RatersPerPerson"])
        * int(manifest_row["Tasks"])
        * int(manifest_row["Criteria"])
    )
    n_categories = int(manifest_row["Categories"])
    category_counts = (
        data["Score"].value_counts().reindex(range(n_categories), fill_value=0)
        .rename_axis("Score").reset_index(name="Count")
    )
    return pilot.GeneratedDesign(
        data=data,
        facet_truth=facet_truth,
        anchors=run_anchors,
        realized_missing_rate=1.0 - len(data) / max(full_rows, 1),
        category_counts=category_counts,
    )


def finite_difference_check(
    parameters: np.ndarray,
    analytical_gradient: np.ndarray,
    idx: dict,
    config: dict,
    sizes: dict,
    coordinate_count: int,
) -> pd.DataFrame:
    dimension = int(parameters.size)
    if dimension == 0:
        return pd.DataFrame()
    count = min(max(1, int(coordinate_count)), dimension)
    coordinates = np.unique(np.linspace(0, dimension - 1, count, dtype=int))
    epsilon_step = np.finfo(float).eps ** (1.0 / 3.0)
    rows: list[dict[str, object]] = []
    for coordinate in coordinates:
        step = epsilon_step * (1.0 + abs(float(parameters[coordinate])))
        plus = parameters.copy()
        minus = parameters.copy()
        plus[coordinate] += step
        minus[coordinate] -= step
        value_plus = app.mfrm_loglik_jmle_value_grad(plus, idx, config, sizes)[0]
        value_minus = app.mfrm_loglik_jmle_value_grad(minus, idx, config, sizes)[0]
        finite_difference = (float(value_plus) - float(value_minus)) / (2.0 * step)
        analytical = float(analytical_gradient[coordinate])
        rows.append({
            "Coordinate": int(coordinate),
            "Step": step,
            "AnalyticalGradient": analytical,
            "FiniteDifferenceGradient": finite_difference,
            "AbsDifference": abs(analytical - finite_difference),
        })
    return pd.DataFrame(rows)


def anchor_check(result: dict, anchors: pd.DataFrame, tolerance: float) -> tuple[bool, float]:
    if anchors.empty:
        return True, 0.0
    estimates = result["facets"]["others"][["Facet", "Level", "Estimate"]].copy()
    merged = anchors.merge(estimates, on=["Facet", "Level"], how="left", validate="one_to_one")
    deviation = (
        pd.to_numeric(merged["Estimate"], errors="coerce")
        - pd.to_numeric(merged["Anchor"], errors="coerce")
    ).abs()
    maximum = float(deviation.max()) if len(deviation) else np.nan
    passed = bool(len(deviation) == len(anchors) and deviation.notna().all() and (deviation <= tolerance).all())
    return passed, maximum


def focal_bias(result: dict, diagnostics: dict) -> dict[str, object]:
    output: dict[str, object] = {
        "BiasAvailable": False,
        "FocalCellSparse": True,
        "BiasEstimate": np.nan,
        "BiasSE": np.nan,
        "p_holm": np.nan,
        "AbsBias": np.nan,
        "DecisionStrongRaw": pd.NA,
    }
    bundle = app.estimate_bias_interaction(result, diagnostics, "Rater", "Task", omit_extreme=False)
    if not isinstance(bundle, dict) or "table" not in bundle:
        return output
    table = app.build_dff_bias_screening_table(
        bundle,
        alpha=0.05,
        min_n=5,
        practical_logit=0.50,
    )
    focal = table.loc[
        table["FacetA_Level"].astype(str).eq(pilot.FOCAL_RATER)
        & table["FacetB_Level"].astype(str).eq(pilot.FOCAL_TASK)
    ]
    if focal.empty:
        return output
    row = focal.iloc[0]
    p_holm = float(pd.to_numeric(pd.Series([row.get("p_holm")]), errors="coerce").iloc[0])
    abs_bias = float(pd.to_numeric(pd.Series([row.get("AbsBias")]), errors="coerce").iloc[0])
    sparse = bool(row.get("SparseCell", True))
    output.update({
        "BiasAvailable": bool(np.isfinite(p_holm) and np.isfinite(abs_bias)),
        "FocalCellSparse": sparse,
        "BiasEstimate": float(row.get("BiasSize", np.nan)),
        "BiasSE": float(row.get("SE", np.nan)),
        "p_holm": p_holm,
        "AbsBias": abs_bias,
        "DecisionStrongRaw": bool(
            np.isfinite(p_holm)
            and p_holm < 0.05
            and np.isfinite(abs_bias)
            and abs_bias >= 0.50
            and not sparse
        ),
    })
    return output


def fit_run(
    manifest_row: pd.Series,
    generated: pilot.GeneratedDesign,
    matched_row: pd.Series,
    plan: dict,
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    controls = plan["numerical_controls"]
    verification = plan["gradient_verification"]
    run: dict[str, object] = {
        "RunId": str(manifest_row["RunId"]),
        "ConditionId": str(manifest_row["ConditionId"]),
        "Design": str(manifest_row["Design"]),
        "TruthBias": float(manifest_row["TruthBias"]),
        "Replicate": int(manifest_row["Replicate"]),
        "Seed": int(manifest_row["Seed"]),
        "Engine": "PythonApp",
        "Estimator": "JMLE",
        "Mode": plan["stage_a"]["label"],
        "RequestedMaxIt": int(controls["maxit"]),
        "RequestedRelTol": float(controls["reltol"]),
        "Rows": int(len(generated.data)),
        "RequestedAnchorRows": int(len(generated.anchors)),
        "FitReturned": False,
        "OptimizerSuccess": False,
        "OptimizerStatus": np.nan,
        "OptimizerMessage": "",
        "NumericalQualificationPassed": False,
        "AnalysisEligible": False,
        "FailureStage": "fit",
        "FailureReason": "",
        "Iterations": np.nan,
        "FunctionEvaluations": np.nan,
        "GradientEvaluations": np.nan,
        "LogLik": np.nan,
        "MatchedLogLik": float(matched_row.get("LogLik", np.nan)),
        "LogLikDeltaStrictMinusMatched": np.nan,
        "RecomputedObjective": np.nan,
        "TerminalGradientL2Norm": np.nan,
        "TerminalGradientSupNorm": np.nan,
        "GradientSupNormThreshold": float(controls["terminal_gradient_sup_norm_threshold"]),
        "GradientGatePassed": False,
        "OptimizerJacRecomputedMaxAbsDifference": np.nan,
        "OptimizerJacGatePassed": False,
        "FiniteDifferenceCoordinates": 0,
        "FiniteDifferenceMaxAbsError": np.nan,
        "FiniteDifferenceGatePassed": False,
        "AnchorMaxAbsDeviation": np.nan,
        "AnchorGatePassed": False,
        "MatchedLogLikGatePassed": False,
        "BiasAvailable": False,
        "FocalCellSparse": True,
        "ElapsedSeconds": np.nan,
        "Warnings": "",
    }
    parameter_rows = pd.DataFrame()
    fit_audit = pd.DataFrame()
    fd_rows = pd.DataFrame()
    started = time.perf_counter()
    caught_messages: list[str] = []
    diagnostics: dict = {}
    try:
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
                maxit=int(controls["maxit"]),
                reltol=float(controls["reltol"]),
                keep_original=True,
            )
            diagnostics = app.mfrm_diagnostics(result, compute_pca=False, compute_marginal=False)
            caught_messages.extend(str(item.message) for item in caught)
        run["FitReturned"] = True
        summary = result["summary"].iloc[0]
        optimizer = result["opt"]
        parameters = np.asarray(optimizer.x, dtype=float)
        config = result["config"]
        sizes = app.build_param_sizes(config)
        idx = app.build_indices(
            result["prep"],
            step_facet=config.get("step_facet"),
            slope_facet=config.get("slope_facet"),
        )
        objective, gradient = app.mfrm_loglik_jmle_value_grad(parameters, idx, config, sizes)
        optimizer_jac = np.asarray(getattr(optimizer, "jac", np.full_like(gradient, np.nan)), dtype=float)
        gradient_sup = float(np.max(np.abs(gradient))) if gradient.size else 0.0
        gradient_l2 = float(np.linalg.norm(gradient))
        jac_difference = (
            float(np.max(np.abs(optimizer_jac - gradient)))
            if optimizer_jac.shape == gradient.shape and np.isfinite(optimizer_jac).all() else
            np.nan
        )
        fd_rows = finite_difference_check(
            parameters,
            gradient,
            idx,
            config,
            sizes,
            int(verification["finite_difference_coordinates_per_fit"]),
        )
        fd_max = float(fd_rows["AbsDifference"].max()) if not fd_rows.empty else np.nan
        anchors_passed, anchor_max = anchor_check(
            result,
            generated.anchors,
            float(controls["hard_anchor_tolerance"]),
        )
        loglik = float(summary["LogLik"])
        matched_loglik = float(run["MatchedLogLik"])
        loglik_delta = loglik - matched_loglik
        optimizer_success = bool(getattr(optimizer, "success", False))
        gradient_passed = bool(np.isfinite(gradient_sup) and gradient_sup <= float(controls["terminal_gradient_sup_norm_threshold"]))
        jac_passed = bool(np.isfinite(jac_difference) and jac_difference <= float(verification["optimizer_jac_recomputed_max_abs_tolerance"]))
        fd_passed = bool(np.isfinite(fd_max) and fd_max <= float(verification["finite_difference_max_abs_error_tolerance"]))
        loglik_passed = bool(
            np.isfinite(loglik_delta)
            and loglik_delta >= -float(controls["matched_loglik_degradation_tolerance"])
        )
        qualification = bool(
            optimizer_success
            and gradient_passed
            and jac_passed
            and fd_passed
            and anchors_passed
            and loglik_passed
        )
        bias = focal_bias(result, diagnostics)
        analysis_eligible = bool(
            qualification and bias["BiasAvailable"] and not bias["FocalCellSparse"]
        )
        run.update({
            "OptimizerSuccess": optimizer_success,
            "OptimizerStatus": int(getattr(optimizer, "status", -1)),
            "OptimizerMessage": str(getattr(optimizer, "message", "")),
            "NumericalQualificationPassed": qualification,
            "AnalysisEligible": analysis_eligible,
            "Iterations": int(getattr(optimizer, "nit", 0)),
            "FunctionEvaluations": int(getattr(optimizer, "nfev", 0)),
            "GradientEvaluations": int(getattr(optimizer, "njev", 0)),
            "LogLik": loglik,
            "LogLikDeltaStrictMinusMatched": loglik_delta,
            "RecomputedObjective": float(objective),
            "TerminalGradientL2Norm": gradient_l2,
            "TerminalGradientSupNorm": gradient_sup,
            "GradientGatePassed": gradient_passed,
            "OptimizerJacRecomputedMaxAbsDifference": jac_difference,
            "OptimizerJacGatePassed": jac_passed,
            "FiniteDifferenceCoordinates": int(len(fd_rows)),
            "FiniteDifferenceMaxAbsError": fd_max,
            "FiniteDifferenceGatePassed": fd_passed,
            "AnchorMaxAbsDeviation": anchor_max,
            "AnchorGatePassed": anchors_passed,
            "MatchedLogLikGatePassed": loglik_passed,
            **bias,
        })
        failures = []
        for passed, reason in (
            (optimizer_success, "optimizer_success_false"),
            (gradient_passed, "gradient_sup_norm_above_threshold"),
            (jac_passed, "optimizer_jac_recomputation_mismatch"),
            (fd_passed, "finite_difference_gradient_mismatch"),
            (anchors_passed, "hard_anchor_constraint_mismatch"),
            (loglik_passed, "strict_loglik_worse_than_matched"),
        ):
            if not passed:
                failures.append(reason)
        if failures:
            run["FailureStage"] = "numerical_qualification"
            run["FailureReason"] = ";".join(failures)
        elif not bias["BiasAvailable"] or bias["FocalCellSparse"]:
            run["FailureStage"] = "bias_scope"
            run["FailureReason"] = (
                "numerically qualified; focal bias cell sparse or unavailable"
            )
        else:
            run["FailureStage"] = ""
            run["FailureReason"] = ""
        parameter_rows = pilot._parameter_rows(
            manifest_row,
            diagnostics,
            generated,
            include=qualification,
        )
        if not parameter_rows.empty:
            parameter_rows.insert(4, "Mode", plan["stage_a"]["label"])
            parameter_rows["NumericalQualificationPassed"] = qualification
        fit_source = diagnostics.get("fit", diagnostics.get("measures", pd.DataFrame()))
        fit_audit = ds.audit_fit_decision_stability(fit_source)
    except Exception as exc:
        run["FailureReason"] = f"{type(exc).__name__}: {exc}"[:1000]
    finally:
        run["Warnings"] = " | ".join(caught_messages)[:2000]
        run["ElapsedSeconds"] = time.perf_counter() - started

    identity = {
        "RunId": run["RunId"],
        "ConditionId": run["ConditionId"],
        "Design": run["Design"],
        "TruthBias": run["TruthBias"],
        "Replicate": run["Replicate"],
        "Mode": run["Mode"],
    }
    for frame in (fit_audit, fd_rows):
        if frame.empty:
            continue
        for column, value in reversed(identity.items()):
            frame.insert(0, column, value)
    return run, parameter_rows, fit_audit, fd_rows


def stage_gates(runs: pd.DataFrame, input_checks: pd.DataFrame, plan: dict) -> pd.DataFrame:
    expected = int(plan["stage_a"]["attempted_fits"])
    anchor_runs = runs.loc[pd.to_numeric(runs["RequestedAnchorRows"], errors="coerce").fillna(0).gt(0)]
    definitions = [
        ("input_identity", len(input_checks), int(input_checks["Passed"].sum())),
        ("fit_returned", expected, int(runs["FitReturned"].astype(bool).sum())),
        ("optimizer_success", expected, int(runs["OptimizerSuccess"].astype(bool).sum())),
        ("finite_recomputed_gradient", expected, int(pd.to_numeric(runs["TerminalGradientSupNorm"], errors="coerce").notna().sum())),
        ("terminal_gradient_sup_norm", expected, int(runs["GradientGatePassed"].astype(bool).sum())),
        ("optimizer_jac_recomputed", expected, int(runs["OptimizerJacGatePassed"].astype(bool).sum())),
        ("finite_difference_gradient", expected, int(runs["FiniteDifferenceGatePassed"].astype(bool).sum())),
        ("hard_anchor_constraint", len(anchor_runs), int(anchor_runs["AnchorGatePassed"].astype(bool).sum())),
        ("matched_loglik_non_degradation", expected, int(runs["MatchedLogLikGatePassed"].astype(bool).sum())),
        ("numerical_qualification", expected, int(runs["NumericalQualificationPassed"].astype(bool).sum())),
    ]
    rows = []
    for gate, required, passed in definitions:
        rows.append({
            "Gate": gate,
            "Required": int(required),
            "Passed": int(passed),
            "GatePassed": int(passed) == int(required),
        })
    all_passed = all(bool(row["GatePassed"]) for row in rows)
    rows.append({
        "Gate": "stage_b_authorization",
        "Required": len(rows),
        "Passed": sum(bool(row["GatePassed"]) for row in rows),
        "GatePassed": all_passed,
    })
    return pd.DataFrame(rows)


def plot_gradients(runs: pd.DataFrame, output: Path, threshold: float) -> None:
    labels = [
        f"{row.Design.replace('_', ' ')}\nbias={row.TruthBias:g} r{int(row.Replicate)}"
        for row in runs.itertuples(index=False)
    ]
    values = pd.to_numeric(runs["TerminalGradientSupNorm"], errors="coerce")
    colors = np.where(runs["GradientGatePassed"].astype(bool), "#59A14F", "#E15759")
    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.scatter(np.arange(len(runs)), values, color=colors, s=48)
    ax.axhline(threshold, color="#222222", linestyle="--", linewidth=1.2, label=f"registered threshold {threshold:.0e}")
    ax.set_yscale("log")
    ax.set_xticks(np.arange(len(runs)), labels, rotation=42, ha="right")
    ax.set_ylabel("Recomputed terminal gradient sup norm")
    ax.set_title("Prospective strict Python JMLE Stage-A numerical qualification")
    ax.legend(frameon=False)
    ax.grid(axis="y", which="both", alpha=0.22)
    fig.tight_layout()
    fig.savefig(output / "strict_jmle_gradient_qualification.png", dpi=180)
    plt.close(fig)


def plot_loglik(runs: pd.DataFrame, output: Path) -> None:
    delta = pd.to_numeric(runs["LogLikDeltaStrictMinusMatched"], errors="coerce")
    fig, ax = plt.subplots(figsize=(10.5, 6))
    ax.bar(np.arange(len(runs)), delta, color=np.where(delta >= 0, "#59A14F", "#E15759"))
    ax.axhline(0, color="#222222", linewidth=1)
    ax.set_xlabel("Frozen smoke RunId order")
    ax.set_ylabel("Strict minus matched-control log likelihood")
    ax.set_title("Strict optimization must not degrade the retained objective")
    ax.grid(axis="y", alpha=0.22)
    fig.tight_layout()
    fig.savefig(output / "strict_jmle_loglik_change.png", dpi=180)
    plt.close(fig)


def write_report(output: Path, runs: pd.DataFrame, gates: pd.DataFrame, plan: dict) -> None:
    authorized = bool(gates.loc[gates["Gate"].eq("stage_b_authorization"), "GatePassed"].iloc[0])
    qualified = int(runs["NumericalQualificationPassed"].astype(bool).sum())
    analysis_eligible = int(runs["AnalysisEligible"].astype(bool).sum())
    gradient = pd.to_numeric(runs["TerminalGradientSupNorm"], errors="coerce")
    fd = pd.to_numeric(runs["FiniteDifferenceMaxAbsError"], errors="coerce")
    jac = pd.to_numeric(runs["OptimizerJacRecomputedMaxAbsDifference"], errors="coerce")
    delta = pd.to_numeric(runs["LogLikDeltaStrictMinusMatched"], errors="coerce")
    failed_lines = "\n".join(
        f"- `{row.Gate}`: {int(row.Passed)}/{int(row.Required)}"
        for row in gates.loc[~gates["GatePassed"]].itertuples(index=False)
    ) or "- None."
    decision = (
        "Stage B is authorized under the frozen Stage-A gates."
        if authorized else
        "Stage B is not authorized; diagnose the failed numerical gates without changing this plan post hoc."
    )
    report = f"""# Strict Python JMLE Stage-A qualification

## Decision

{decision}

## Outcome

- {len(runs)}/{plan['stage_a']['attempted_fits']} frozen smoke RunIds attempted.
- {int(runs['FitReturned'].astype(bool).sum())}/{len(runs)} fits returned; {int(runs['OptimizerSuccess'].astype(bool).sum())}/{len(runs)} optimizer success flags.
- {qualified}/{len(runs)} passed every numerical-qualification gate.
- {analysis_eligible}/{len(runs)} also had an available non-sparse focal bias decision; this separate count does not affect numerical qualification.
- Recomputed terminal-gradient sup norm range: {gradient.min():.6g} to {gradient.max():.6g} (registered maximum `{plan['numerical_controls']['terminal_gradient_sup_norm_threshold']}`).
- Maximum optimizer-jac/recomputed-gradient difference: {jac.max():.6g}.
- Maximum centered finite-difference discrepancy: {fd.max():.6g}.
- Strict-minus-matched log-likelihood range: {delta.min():.6g} to {delta.max():.6g}.
- Recorded fitting and verification time: {runs['ElapsedSeconds'].sum():.2f} seconds.

## Failed gates

{failed_lines}

## Interpretation boundary

This is a numerical-qualification result on frozen synthetic smoke inputs. It
does not estimate power, false-positive rate, coverage, anchor robustness, or
package superiority. Sparse focal bias cells remain unavailable rather than
becoming negative decisions. Stage B may run only if the recorded authorization
gate passes; public application status remains `Withheld`.

## Retained artifacts

- `strict_jmle_runs.csv`
- `strict_jmle_parameter_recovery.csv`
- `strict_jmle_fit_decision_stability.csv`
- `strict_jmle_finite_difference_gradient.csv`
- `strict_jmle_stage_gates.csv`
- `strict_jmle_input_checks.csv`
- `strict_jmle_first_read_summary.csv`
- `strict_jmle_identity.json`
- `strict_jmle_gradient_qualification.png`
- `strict_jmle_loglik_change.png`
"""
    (output / "STRICT_JMLE_RESULTS.md").write_text(report, encoding="utf-8")


def run(input_dir: Path, output: Path, plan_path: Path) -> None:
    input_dir = input_dir.resolve()
    output = output.resolve()
    plan_path = plan_path.resolve()
    plan = load_plan(plan_path)
    input_checks = validate_input(input_dir, plan)
    manifest = pd.read_csv(input_dir / "manifest.csv")
    ratings = pd.read_csv(input_dir / "generated_ratings.csv")
    truth = pd.read_csv(input_dir / "generated_facet_truth.csv")
    anchors = pd.read_csv(input_dir / "generated_anchors.csv")
    matched = pd.read_csv(input_dir / "runs.csv").set_index("RunId", drop=False)
    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_strict_jmle_plan.json").write_bytes(plan_path.read_bytes())

    run_rows: list[dict[str, object]] = []
    parameter_parts: list[pd.DataFrame] = []
    audit_parts: list[pd.DataFrame] = []
    fd_parts: list[pd.DataFrame] = []
    for _, manifest_row in manifest.iterrows():
        run_id = str(manifest_row["RunId"])
        generated = generated_design_for_run(manifest_row, ratings, truth, anchors)
        run_row, parameters, fit_audit, fd_rows = fit_run(
            manifest_row,
            generated,
            matched.loc[run_id],
            plan,
        )
        run_rows.append(run_row)
        if not parameters.empty:
            parameter_parts.append(parameters)
        if not fit_audit.empty:
            audit_parts.append(fit_audit)
        if not fd_rows.empty:
            fd_parts.append(fd_rows)

    runs = pd.DataFrame(run_rows)
    parameters = pd.concat(parameter_parts, ignore_index=True) if parameter_parts else pd.DataFrame()
    audits = pd.concat(audit_parts, ignore_index=True) if audit_parts else pd.DataFrame()
    finite_difference = pd.concat(fd_parts, ignore_index=True) if fd_parts else pd.DataFrame()
    gates = stage_gates(runs, input_checks, plan)
    authorized = bool(gates.loc[gates["Gate"].eq("stage_b_authorization"), "GatePassed"].iloc[0])
    first_read = pd.DataFrame([
        {
            "Priority": 1,
            "Check": "Strict JMLE Stage A",
            "Status": "Pass" if authorized else "Blocked",
            "Evidence": f"{int(runs['NumericalQualificationPassed'].sum())}/{len(runs)} passed all numerical gates.",
            "NextAction": (
                "Run Stage B with unchanged controls in a separate directory."
                if authorized else
                "Diagnose failed gates; do not alter the registered Stage-A thresholds after seeing results."
            ),
            "PublicSurfaceEnabled": False,
        },
        {
            "Priority": 2,
            "Check": "Bias scope",
            "Status": "Separate",
            "Evidence": f"{int(runs['AnalysisEligible'].sum())}/{len(runs)} were additionally eligible for the focal bias decision.",
            "NextAction": "Keep sparse-cell availability separate from optimizer qualification.",
            "PublicSurfaceEnabled": False,
        },
        {
            "Priority": 3,
            "Check": "Public application surface",
            "Status": "Withheld",
            "Evidence": "Repository-only numerical qualification; no performance or package-ranking claim.",
            "NextAction": "Retain disabled until the operating-characteristics promotion gates pass.",
            "PublicSurfaceEnabled": False,
        },
    ])
    runs.to_csv(output / "strict_jmle_runs.csv", index=False)
    parameters.to_csv(output / "strict_jmle_parameter_recovery.csv", index=False)
    audits.to_csv(output / "strict_jmle_fit_decision_stability.csv", index=False)
    finite_difference.to_csv(output / "strict_jmle_finite_difference_gradient.csv", index=False)
    gates.to_csv(output / "strict_jmle_stage_gates.csv", index=False)
    input_checks.to_csv(output / "strict_jmle_input_checks.csv", index=False)
    first_read.to_csv(output / "strict_jmle_first_read_summary.csv", index=False)
    plot_gradients(
        runs,
        output,
        float(plan["numerical_controls"]["terminal_gradient_sup_norm_threshold"]),
    )
    plot_loglik(runs, output)
    write_report(output, runs, gates, plan)
    artifact_names = [
        "strict_jmle_runs.csv",
        "strict_jmle_parameter_recovery.csv",
        "strict_jmle_fit_decision_stability.csv",
        "strict_jmle_finite_difference_gradient.csv",
        "strict_jmle_stage_gates.csv",
        "strict_jmle_input_checks.csv",
        "strict_jmle_first_read_summary.csv",
        "strict_jmle_gradient_qualification.png",
        "strict_jmle_loglik_change.png",
        "STRICT_JMLE_RESULTS.md",
    ]
    identity = {
        "schema_version": PLAN_SCHEMA,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "stage_b_authorized": authorized,
        "strict_plan_sha256": sha256_file(plan_path),
        "adapter_script_sha256": sha256_file(Path(__file__).resolve()),
        "streamlit_app_sha256": sha256_file(REPO / "streamlit_app.py"),
        "pilot_runner_sha256": sha256_file(REPO / "validation" / "operating_characteristics_pilot.py"),
        "input_manifest_frame_sha256": oc.frame_fingerprint(manifest),
        "artifacts": {name: sha256_file(output / name) for name in artifact_names},
    }
    (output / "strict_jmle_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    args = parser.parse_args()
    run(args.input, args.output, args.plan)


if __name__ == "__main__":
    main()
