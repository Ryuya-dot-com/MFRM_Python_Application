#!/usr/bin/env python3
"""Run prospective Stage-A2 gradient-directed Python JMLE polish qualification."""

from __future__ import annotations

import argparse
import json
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
from scipy.optimize import minimize


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import streamlit_app as app  # noqa: E402
from validation import operating_characteristics_pilot as pilot  # noqa: E402
from validation import operating_characteristics_strict_jmle as stage_a  # noqa: E402


PLAN_SCHEMA = "mfrm-strict-jmle-a2-polish-plan-v1"
DEFAULT_INPUT = REPO / "validation" / "operating_characteristics_20260809"
DEFAULT_STAGE_A = REPO / "validation" / "operating_characteristics_strict_jmle_smoke_20260809"
DEFAULT_STAGE_A_PLAN = REPO / "validation" / "operating_characteristics_strict_jmle_plan_20260809.json"
DEFAULT_PLAN = REPO / "validation" / "operating_characteristics_strict_jmle_a2_plan_20260809.json"
DEFAULT_OUTPUT = REPO / "validation" / "operating_characteristics_strict_jmle_a2_smoke_20260809"


def load_plan(path: Path) -> dict:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != PLAN_SCHEMA:
        raise ValueError(f"Stage-A2 plan schema must be {PLAN_SCHEMA}")
    qualification = plan.get("qualification", {})
    if int(qualification.get("attempted_fits_per_candidate", 0)) != 16:
        raise ValueError("Stage-A2 must retain 16 attempted fits per candidate")
    if float(qualification.get("terminal_gradient_sup_norm_threshold", np.nan)) != 1e-4:
        raise ValueError("Stage-A2 terminal-gradient threshold must remain 1e-4")
    if float(qualification.get("optimizer_jac_recomputed_max_abs_tolerance", np.nan)) != 1e-12:
        raise ValueError("Stage-A2 optimizer-jac tolerance must remain 1e-12")
    if float(qualification.get("finite_difference_max_abs_error_tolerance", np.nan)) != 1e-5:
        raise ValueError("Stage-A2 finite-difference tolerance must remain 1e-5")
    if int(qualification.get("finite_difference_coordinates_per_fit", 0)) != 8:
        raise ValueError("Stage-A2 finite-difference coordinate count must remain 8")
    candidates = plan.get("candidates", [])
    expected = {
        "lbfgsb_precision": (
            1,
            "L-BFGS-B",
            {"maxiter": 1500, "gtol": 1e-7, "ftol": 1e-15, "maxls": 50, "maxcor": 20},
        ),
        "bfgs_gradient": (
            2,
            "BFGS",
            {"maxiter": 1000, "gtol": 1e-7, "xrtol": 0.0, "c1": 1e-4, "c2": 0.9},
        ),
    }
    if [candidate.get("candidate_id") for candidate in candidates] != list(expected):
        raise ValueError("Stage-A2 candidate order or identity changed")
    for candidate in candidates:
        priority, method, options = expected[str(candidate["candidate_id"])]
        if int(candidate.get("priority", 0)) != priority:
            raise ValueError(f"Stage-A2 priority changed for {candidate['candidate_id']}")
        if candidate.get("method") != method or candidate.get("options") != options:
            raise ValueError(f"Stage-A2 controls changed for {candidate['candidate_id']}")
    boundaries = plan.get("scope_boundaries", {})
    if boundaries.get("public_surface_status") != "Withheld":
        raise ValueError("Stage-A2 public surface must remain Withheld")
    if bool(boundaries.get("application_core_change_authorized_by_this_plan", True)):
        raise ValueError("Stage-A2 must not authorize an application-core change")
    return plan


def validate_input_identity(
    input_dir: Path,
    stage_a_dir: Path,
    stage_a_plan_path: Path,
    plan: dict,
) -> pd.DataFrame:
    expected = plan["input_identity"]
    paths = {
        "stage_a_plan_sha256": stage_a_plan_path,
        "stage_a_adapter_sha256": Path(stage_a.__file__).resolve(),
        "streamlit_app_sha256": REPO / "streamlit_app.py",
        "stage_a_runs_sha256": stage_a_dir / "strict_jmle_runs.csv",
        "stage_a_identity_sha256": stage_a_dir / "strict_jmle_identity.json",
        "manifest_sha256": input_dir / "manifest.csv",
        "generated_ratings_sha256": input_dir / "generated_ratings.csv",
        "generated_facet_truth_sha256": input_dir / "generated_facet_truth.csv",
        "generated_anchors_sha256": input_dir / "generated_anchors.csv",
    }
    rows: list[dict[str, object]] = []
    for key, path in paths.items():
        actual = stage_a.sha256_file(path)
        target = str(expected[key])
        rows.append({
            "Check": key,
            "Passed": actual == target,
            "Evidence": f"actual={actual}; expected={target}",
        })
    prior_identity = json.loads((stage_a_dir / "strict_jmle_identity.json").read_text(encoding="utf-8"))
    rows.append({
        "Check": "original_stage_b_remains_unauthorized",
        "Passed": prior_identity.get("stage_b_authorized") is False,
        "Evidence": f"stage_b_authorized={prior_identity.get('stage_b_authorized')}",
    })
    prior_runs = pd.read_csv(stage_a_dir / "strict_jmle_runs.csv")
    common_message = str(plan["frozen_diagnosis"]["common_optimizer_message"])
    diagnosis_passed = bool(
        len(prior_runs) == 16
        and prior_runs["OptimizerSuccess"].astype(bool).sum() == 16
        and prior_runs["GradientGatePassed"].astype(bool).sum() == 0
        and prior_runs["OptimizerMessage"].astype(str).eq(common_message).all()
    )
    rows.append({
        "Check": "frozen_stage_a_diagnosis",
        "Passed": diagnosis_passed,
        "Evidence": (
            f"runs={len(prior_runs)}; optimizer_success={int(prior_runs['OptimizerSuccess'].astype(bool).sum())}; "
            f"gradient_pass={int(prior_runs['GradientGatePassed'].astype(bool).sum())}; "
            f"message_match={int(prior_runs['OptimizerMessage'].astype(str).eq(common_message).sum())}"
        ),
    })
    checks = pd.DataFrame(rows)
    if not checks["Passed"].all():
        failed = checks.loc[~checks["Passed"], "Check"].astype(str).tolist()
        raise ValueError(f"Stage-A2 input identity failed: {failed}")
    return checks


def fit_stage_a_baseline(
    manifest_row: pd.Series,
    generated: pilot.GeneratedDesign,
    controls: dict,
) -> tuple[dict, float, str]:
    started = time.perf_counter()
    caught_messages: list[str] = []
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
        caught_messages.extend(str(item.message) for item in caught)
    return result, time.perf_counter() - started, " | ".join(caught_messages)[:2000]


def optimizer_context(result: dict) -> tuple[np.ndarray, dict, dict, dict]:
    optimizer = result["opt"]
    config = result["config"]
    sizes = app.build_param_sizes(config)
    idx = app.build_indices(
        result["prep"],
        step_facet=config.get("step_facet"),
        slope_facet=config.get("slope_facet"),
    )
    return np.asarray(optimizer.x, dtype=float), idx, config, sizes


def baseline_audit(
    manifest_row: pd.Series,
    result: dict,
    elapsed: float,
    warning_text: str,
    stored: pd.Series,
    plan: dict,
) -> tuple[dict[str, object], np.ndarray, dict, dict, dict, float, np.ndarray]:
    parameters, idx, config, sizes = optimizer_context(result)
    objective, gradient = app.mfrm_loglik_jmle_value_grad(parameters, idx, config, sizes)
    gradient_sup = float(np.max(np.abs(gradient))) if gradient.size else 0.0
    optimizer = result["opt"]
    loglik = float(result["summary"].iloc[0]["LogLik"])
    reproduction = plan["baseline_reproduction"]
    loglik_difference = abs(loglik - float(stored["LogLik"]))
    gradient_difference = abs(gradient_sup - float(stored["TerminalGradientSupNorm"]))
    message = str(getattr(optimizer, "message", ""))
    message_match = message == str(stored["OptimizerMessage"])
    reproduced = bool(
        loglik_difference <= float(reproduction["loglik_max_abs_difference_tolerance"])
        and gradient_difference <= float(reproduction["terminal_gradient_sup_norm_max_abs_difference_tolerance"])
        and message_match
    )
    row = {
        "RunId": str(manifest_row["RunId"]),
        "ConditionId": str(manifest_row["ConditionId"]),
        "Design": str(manifest_row["Design"]),
        "TruthBias": float(manifest_row["TruthBias"]),
        "Replicate": int(manifest_row["Replicate"]),
        "Rows": int(len(result["prep"]["data"])),
        "OptimizerDimension": int(parameters.size),
        "BaselineOptimizerSuccess": bool(getattr(optimizer, "success", False)),
        "BaselineOptimizerStatus": int(getattr(optimizer, "status", -1)),
        "BaselineOptimizerMessage": message,
        "BaselineIterations": int(getattr(optimizer, "nit", 0)),
        "BaselineFunctionEvaluations": int(getattr(optimizer, "nfev", 0)),
        "BaselineObjective": float(objective),
        "BaselineLogLik": loglik,
        "StoredStageALogLik": float(stored["LogLik"]),
        "BaselineLogLikMaxAbsDifference": loglik_difference,
        "BaselineTerminalGradientSupNorm": gradient_sup,
        "StoredStageATerminalGradientSupNorm": float(stored["TerminalGradientSupNorm"]),
        "BaselineGradientMaxAbsDifference": gradient_difference,
        "BaselineOptimizerMessageMatched": message_match,
        "BaselineReproduced": reproduced,
        "BaselineElapsedSeconds": elapsed,
        "Warnings": warning_text,
    }
    return row, parameters, idx, config, sizes, float(objective), np.asarray(gradient, dtype=float)


def polished_anchor_check(
    parameters: np.ndarray,
    config: dict,
    sizes: dict,
    anchors: pd.DataFrame,
    tolerance: float,
) -> tuple[bool, float]:
    if anchors.empty:
        return True, 0.0
    expanded = app.expand_params(parameters, sizes, config)
    deviations: list[float] = []
    for anchor in anchors.itertuples(index=False):
        facet = str(anchor.Facet)
        level = str(anchor.Level)
        levels = [str(value) for value in config["facet_levels"][facet]]
        if level not in levels:
            return False, np.nan
        estimate = float(expanded["facets"][facet][levels.index(level)])
        deviations.append(abs(estimate - float(anchor.Anchor)))
    maximum = max(deviations) if deviations else 0.0
    return bool(np.isfinite(maximum) and maximum <= tolerance), float(maximum)


def polish_one(
    baseline_row: dict[str, object],
    start: np.ndarray,
    idx: dict,
    config: dict,
    sizes: dict,
    baseline_objective: float,
    anchors: pd.DataFrame,
    candidate: dict,
    plan: dict,
) -> tuple[dict[str, object], pd.DataFrame]:
    qualification = plan["qualification"]
    row: dict[str, object] = {
        "RunId": baseline_row["RunId"],
        "ConditionId": baseline_row["ConditionId"],
        "Design": baseline_row["Design"],
        "TruthBias": baseline_row["TruthBias"],
        "Replicate": baseline_row["Replicate"],
        "Rows": baseline_row["Rows"],
        "OptimizerDimension": baseline_row["OptimizerDimension"],
        "CandidateId": str(candidate["candidate_id"]),
        "CandidatePriority": int(candidate["priority"]),
        "Method": str(candidate["method"]),
        "BaselineReproduced": bool(baseline_row["BaselineReproduced"]),
        "BaselineObjective": baseline_objective,
        "BaselineTerminalGradientSupNorm": baseline_row["BaselineTerminalGradientSupNorm"],
        "FitReturned": False,
        "OptimizerSuccess": False,
        "OptimizerStatus": np.nan,
        "OptimizerMessage": "",
        "Iterations": np.nan,
        "FunctionEvaluations": np.nan,
        "GradientEvaluations": np.nan,
        "Objective": np.nan,
        "LogLikGainFromBaseline": np.nan,
        "ObjectiveNonDegradationGatePassed": False,
        "TerminalGradientL2Norm": np.nan,
        "TerminalGradientSupNorm": np.nan,
        "TerminalGradientSupNormPerObservation": np.nan,
        "GradientGatePassed": False,
        "OptimizerJacRecomputedMaxAbsDifference": np.nan,
        "OptimizerJacGatePassed": False,
        "FiniteDifferenceCoordinates": 0,
        "FiniteDifferenceMaxAbsError": np.nan,
        "FiniteDifferenceGatePassed": False,
        "AnchorMaxAbsDeviation": np.nan,
        "AnchorGatePassed": False,
        "FiniteObjectiveAndGradientGatePassed": False,
        "MaxAbsParameterMovement": np.nan,
        "L2ParameterMovement": np.nan,
        "NumericalQualificationPassed": False,
        "FailureReason": "",
        "ElapsedSeconds": np.nan,
        "Warnings": "",
    }
    fd_rows = pd.DataFrame()
    started = time.perf_counter()
    caught_messages: list[str] = []
    try:
        minimize_kwargs: dict[str, object] = {
            "fun": app.mfrm_loglik_jmle_value_grad,
            "x0": np.array(start, dtype=float, copy=True),
            "args": (idx, config, sizes),
            "jac": True,
            "method": str(candidate["method"]),
            "options": dict(candidate["options"]),
        }
        if str(candidate["method"]) == "L-BFGS-B":
            minimize_kwargs["bounds"] = app.build_optimizer_bounds(sizes, config)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            optimizer = minimize(**minimize_kwargs)
            caught_messages.extend(str(item.message) for item in caught)
        parameters = np.asarray(optimizer.x, dtype=float)
        objective, gradient = app.mfrm_loglik_jmle_value_grad(parameters, idx, config, sizes)
        gradient = np.asarray(gradient, dtype=float)
        optimizer_jac = np.asarray(getattr(optimizer, "jac", np.full_like(gradient, np.nan)), dtype=float)
        gradient_sup = float(np.max(np.abs(gradient))) if gradient.size else 0.0
        gradient_l2 = float(np.linalg.norm(gradient))
        movement = parameters - np.asarray(start, dtype=float)
        jac_difference = (
            float(np.max(np.abs(optimizer_jac - gradient)))
            if optimizer_jac.shape == gradient.shape and np.isfinite(optimizer_jac).all()
            else np.nan
        )
        fd_rows = stage_a.finite_difference_check(
            parameters,
            gradient,
            idx,
            config,
            sizes,
            int(qualification["finite_difference_coordinates_per_fit"]),
        )
        fd_max = float(fd_rows["AbsDifference"].max()) if not fd_rows.empty else np.nan
        anchor_passed, anchor_max = polished_anchor_check(
            parameters,
            config,
            sizes,
            anchors,
            float(qualification["hard_anchor_tolerance"]),
        )
        finite_passed = bool(np.isfinite(objective) and np.isfinite(gradient).all())
        gradient_passed = bool(
            finite_passed
            and gradient_sup <= float(qualification["terminal_gradient_sup_norm_threshold"])
        )
        jac_passed = bool(
            np.isfinite(jac_difference)
            and jac_difference <= float(qualification["optimizer_jac_recomputed_max_abs_tolerance"])
        )
        fd_passed = bool(
            np.isfinite(fd_max)
            and fd_max <= float(qualification["finite_difference_max_abs_error_tolerance"])
        )
        objective_passed = bool(
            np.isfinite(objective)
            and objective - baseline_objective
            <= float(qualification["objective_degradation_tolerance_from_stage_a"])
        )
        optimizer_success = bool(getattr(optimizer, "success", False))
        qualified = bool(
            baseline_row["BaselineReproduced"]
            and optimizer_success
            and finite_passed
            and gradient_passed
            and jac_passed
            and fd_passed
            and anchor_passed
            and objective_passed
        )
        row.update({
            "FitReturned": True,
            "OptimizerSuccess": optimizer_success,
            "OptimizerStatus": int(getattr(optimizer, "status", -1)),
            "OptimizerMessage": str(getattr(optimizer, "message", "")),
            "Iterations": int(getattr(optimizer, "nit", 0)),
            "FunctionEvaluations": int(getattr(optimizer, "nfev", 0)),
            "GradientEvaluations": int(getattr(optimizer, "njev", 0)),
            "Objective": float(objective),
            "LogLikGainFromBaseline": float(baseline_objective - objective),
            "ObjectiveNonDegradationGatePassed": objective_passed,
            "TerminalGradientL2Norm": gradient_l2,
            "TerminalGradientSupNorm": gradient_sup,
            "TerminalGradientSupNormPerObservation": gradient_sup / max(int(baseline_row["Rows"]), 1),
            "GradientGatePassed": gradient_passed,
            "OptimizerJacRecomputedMaxAbsDifference": jac_difference,
            "OptimizerJacGatePassed": jac_passed,
            "FiniteDifferenceCoordinates": int(len(fd_rows)),
            "FiniteDifferenceMaxAbsError": fd_max,
            "FiniteDifferenceGatePassed": fd_passed,
            "AnchorMaxAbsDeviation": anchor_max,
            "AnchorGatePassed": anchor_passed,
            "FiniteObjectiveAndGradientGatePassed": finite_passed,
            "MaxAbsParameterMovement": float(np.max(np.abs(movement))) if movement.size else 0.0,
            "L2ParameterMovement": float(np.linalg.norm(movement)),
            "NumericalQualificationPassed": qualified,
        })
        failed = []
        for passed, reason in (
            (bool(baseline_row["BaselineReproduced"]), "baseline_not_reproduced"),
            (optimizer_success, "optimizer_success_false"),
            (finite_passed, "nonfinite_objective_or_gradient"),
            (gradient_passed, "gradient_sup_norm_above_threshold"),
            (jac_passed, "optimizer_jac_recomputation_mismatch"),
            (fd_passed, "finite_difference_gradient_mismatch"),
            (anchor_passed, "hard_anchor_constraint_mismatch"),
            (objective_passed, "objective_degraded_from_stage_a"),
        ):
            if not passed:
                failed.append(reason)
        row["FailureReason"] = ";".join(failed)
    except Exception as exc:
        row["FailureReason"] = f"{type(exc).__name__}: {exc}"[:1000]
    finally:
        row["ElapsedSeconds"] = time.perf_counter() - started
        row["Warnings"] = " | ".join(caught_messages)[:2000]

    if not fd_rows.empty:
        for column, value in reversed({
            "RunId": row["RunId"],
            "ConditionId": row["ConditionId"],
            "Design": row["Design"],
            "TruthBias": row["TruthBias"],
            "Replicate": row["Replicate"],
            "CandidateId": row["CandidateId"],
        }.items()):
            fd_rows.insert(0, column, value)
    return row, fd_rows


def candidate_gates(
    runs: pd.DataFrame,
    baseline: pd.DataFrame,
    input_checks: pd.DataFrame,
    plan: dict,
) -> pd.DataFrame:
    expected = int(plan["qualification"]["attempted_fits_per_candidate"])
    rows: list[dict[str, object]] = []
    for candidate in plan["candidates"]:
        candidate_id = str(candidate["candidate_id"])
        frame = runs.loc[runs["CandidateId"].eq(candidate_id)]
        anchor_runs = frame.loc[
            frame["RunId"].isin(
                baseline.loc[baseline["RequestedAnchorRows"].gt(0), "RunId"]
            )
        ]
        definitions = [
            ("input_identity", len(input_checks), int(input_checks["Passed"].sum())),
            ("baseline_reproduction", expected, int(frame["BaselineReproduced"].astype(bool).sum())),
            ("fit_returned", expected, int(frame["FitReturned"].astype(bool).sum())),
            ("optimizer_success", expected, int(frame["OptimizerSuccess"].astype(bool).sum())),
            ("finite_objective_and_gradient", expected, int(frame["FiniteObjectiveAndGradientGatePassed"].astype(bool).sum())),
            ("terminal_gradient_sup_norm", expected, int(frame["GradientGatePassed"].astype(bool).sum())),
            ("optimizer_jac_recomputed", expected, int(frame["OptimizerJacGatePassed"].astype(bool).sum())),
            ("finite_difference_gradient", expected, int(frame["FiniteDifferenceGatePassed"].astype(bool).sum())),
            ("hard_anchor_constraint", len(anchor_runs), int(anchor_runs["AnchorGatePassed"].astype(bool).sum())),
            ("objective_non_degradation", expected, int(frame["ObjectiveNonDegradationGatePassed"].astype(bool).sum())),
            ("numerical_qualification", expected, int(frame["NumericalQualificationPassed"].astype(bool).sum())),
        ]
        gate_results = []
        for gate, required, passed in definitions:
            gate_passed = int(passed) == int(required)
            gate_results.append(gate_passed)
            rows.append({
                "CandidateId": candidate_id,
                "CandidatePriority": int(candidate["priority"]),
                "Gate": gate,
                "Required": int(required),
                "Passed": int(passed),
                "GatePassed": gate_passed,
            })
        rows.append({
            "CandidateId": candidate_id,
            "CandidatePriority": int(candidate["priority"]),
            "Gate": "candidate_qualified",
            "Required": len(gate_results),
            "Passed": sum(gate_results),
            "GatePassed": all(gate_results),
        })
    return pd.DataFrame(rows)


def selection_table(gates: pd.DataFrame, plan: dict) -> pd.DataFrame:
    qualified: list[str] = []
    for candidate in sorted(plan["candidates"], key=lambda value: int(value["priority"])):
        candidate_id = str(candidate["candidate_id"])
        passed = bool(
            gates.loc[
                gates["CandidateId"].eq(candidate_id) & gates["Gate"].eq("candidate_qualified"),
                "GatePassed",
            ].iloc[0]
        )
        if passed:
            qualified.append(candidate_id)
    selected = qualified[0] if qualified else ""
    return pd.DataFrame([
        {
            "Decision": "OriginalStageB",
            "Authorized": False,
            "SelectedCandidate": "",
            "Evidence": "Stage A remains failed and cannot be retroactively changed.",
            "PublicSurfaceEnabled": False,
        },
        {
            "Decision": "RevisedStageB2",
            "Authorized": bool(selected),
            "SelectedCandidate": selected,
            "Evidence": (
                f"Selected lowest-priority qualifying candidate: {selected}."
                if selected else
                "No A2 candidate passed every frozen gate on all 16 RunIds."
            ),
            "PublicSurfaceEnabled": False,
        },
        {
            "Decision": "ApplicationCoreChange",
            "Authorized": False,
            "SelectedCandidate": selected,
            "Evidence": "A separate implementation and regression gate is required even if Stage B2 is authorized.",
            "PublicSurfaceEnabled": False,
        },
    ])


def plot_gradient_comparison(runs: pd.DataFrame, output: Path, threshold: float) -> None:
    order = runs["RunId"].drop_duplicates().tolist()
    baseline = runs.drop_duplicates("RunId").set_index("RunId")["BaselineTerminalGradientSupNorm"].reindex(order)
    positions = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.scatter(positions, baseline, color="#4E79A7", marker="x", s=50, label="Stage A baseline")
    offsets = {"lbfgsb_precision": -0.14, "bfgs_gradient": 0.14}
    colors = {"lbfgsb_precision": "#59A14F", "bfgs_gradient": "#F28E2B"}
    for candidate_id, frame in runs.groupby("CandidateId", sort=False):
        values = frame.set_index("RunId")["TerminalGradientSupNorm"].reindex(order)
        ax.scatter(
            positions + offsets.get(str(candidate_id), 0.0),
            values,
            s=42,
            color=colors.get(str(candidate_id), "#777777"),
            label=str(candidate_id),
        )
    ax.axhline(threshold, color="#222222", linestyle="--", linewidth=1.2, label=f"frozen threshold {threshold:.0e}")
    ax.set_yscale("log")
    ax.set_xticks(positions, [value.split("::")[0].replace("__", "\n") for value in order], rotation=40, ha="right")
    ax.set_ylabel("Recomputed terminal gradient sup norm")
    ax.set_title("Stage-A2 gradient-directed polish qualification")
    ax.grid(axis="y", which="both", alpha=0.22)
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(output / "strict_jmle_a2_gradient_comparison.png", dpi=180)
    plt.close(fig)


def plot_effort(runs: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.5))
    order = [str(candidate["candidate_id"]) for candidate in runs.attrs["plan_candidates"]]
    data = [pd.to_numeric(runs.loc[runs["CandidateId"].eq(value), "FunctionEvaluations"], errors="coerce").dropna() for value in order]
    axes[0].boxplot(data, tick_labels=order, showmeans=True)
    axes[0].set_ylabel("Polish function/gradient evaluations")
    axes[0].set_title("Additional optimizer effort")
    gain = [pd.to_numeric(runs.loc[runs["CandidateId"].eq(value), "LogLikGainFromBaseline"], errors="coerce").dropna() for value in order]
    axes[1].boxplot(gain, tick_labels=order, showmeans=True)
    axes[1].axhline(0, color="#222222", linewidth=1)
    axes[1].set_ylabel("Log-likelihood gain from Stage A")
    axes[1].set_title("Objective improvement is small but audited")
    for axis in axes:
        axis.tick_params(axis="x", rotation=15)
        axis.grid(axis="y", alpha=0.22)
    fig.tight_layout()
    fig.savefig(output / "strict_jmle_a2_effort_and_gain.png", dpi=180)
    plt.close(fig)


def write_report(
    output: Path,
    baseline: pd.DataFrame,
    runs: pd.DataFrame,
    gates: pd.DataFrame,
    selection: pd.DataFrame,
    plan: dict,
) -> None:
    selected = str(selection.loc[selection["Decision"].eq("RevisedStageB2"), "SelectedCandidate"].iloc[0])
    authorized = bool(selection.loc[selection["Decision"].eq("RevisedStageB2"), "Authorized"].iloc[0])
    sections: list[str] = []
    for candidate in plan["candidates"]:
        candidate_id = str(candidate["candidate_id"])
        frame = runs.loc[runs["CandidateId"].eq(candidate_id)]
        failed = gates.loc[gates["CandidateId"].eq(candidate_id) & ~gates["GatePassed"]]
        failed_text = ", ".join(
            f"{row.Gate}={int(row.Passed)}/{int(row.Required)}"
            for row in failed.itertuples(index=False)
        ) or "none"
        messages = "; ".join(
            f"{count}x {message}"
            for message, count in frame["OptimizerMessage"].value_counts().items()
        )
        sections.append(
            f"### `{candidate_id}`\n\n"
            f"- {int(frame['NumericalQualificationPassed'].astype(bool).sum())}/{len(frame)} runs passed every gate.\n"
            f"- Terminal-gradient sup norm range: {pd.to_numeric(frame['TerminalGradientSupNorm'], errors='coerce').min():.6g} to {pd.to_numeric(frame['TerminalGradientSupNorm'], errors='coerce').max():.6g}.\n"
            f"- Log-likelihood gain range: {pd.to_numeric(frame['LogLikGainFromBaseline'], errors='coerce').min():.6g} to {pd.to_numeric(frame['LogLikGainFromBaseline'], errors='coerce').max():.6g}.\n"
            f"- Maximum absolute free-coordinate movement: {pd.to_numeric(frame['MaxAbsParameterMovement'], errors='coerce').max():.6g} (descriptive only; not a post-hoc gate).\n"
            f"- Maximum finite-difference error: {pd.to_numeric(frame['FiniteDifferenceMaxAbsError'], errors='coerce').max():.6g}.\n"
            f"- Failed gates: {failed_text}.\n"
            f"- Termination messages: {messages}."
        )
    decision = (
        f"Revised Stage B2 is authorized with `{selected}` under the frozen A2 rule."
        if authorized else
        "Revised Stage B2 is not authorized because no candidate passed every frozen A2 gate."
    )
    report = f"""# Strict Python JMLE Stage-A2 polish qualification

## Decision

{decision} Original Stage B remains unauthorized, application-core integration
is not authorized by this experiment, and the public surface remains withheld.

## Baseline reproduction

- {int(baseline['BaselineReproduced'].astype(bool).sum())}/{len(baseline)} Stage-A terminal states reproduced within the frozen tolerances.
- {int(baseline['BaselineOptimizerMessageMatched'].astype(bool).sum())}/{len(baseline)} optimizer messages matched.
- Maximum reproduced log-likelihood difference: {pd.to_numeric(baseline['BaselineLogLikMaxAbsDifference'], errors='coerce').max():.6g}.
- Maximum reproduced terminal-gradient difference: {pd.to_numeric(baseline['BaselineGradientMaxAbsDifference'], errors='coerce').max():.6g}.

## Candidate results

{chr(10).join(sections)}

## Interpretation boundary

This is a numerical solver qualification on 16 frozen synthetic smoke inputs.
It is not an operating-characteristics result and does not establish bias power,
false-positive control, coverage, sparse-data robustness, equivalence with R
packages, or package superiority. A qualifying candidate may advance only to a
separate Stage B2 on the already frozen 160 RunIds. Statistical decisions still
require raw-value auditing near display thresholds.

## Retained artifacts

- `strict_jmle_a2_baseline.csv`
- `strict_jmle_a2_runs.csv`
- `strict_jmle_a2_finite_difference_gradient.csv`
- `strict_jmle_a2_candidate_gates.csv`
- `strict_jmle_a2_selection.csv`
- `strict_jmle_a2_input_checks.csv`
- `strict_jmle_a2_first_read_summary.csv`
- `strict_jmle_a2_identity.json`
- `strict_jmle_a2_gradient_comparison.png`
- `strict_jmle_a2_effort_and_gain.png`
"""
    (output / "STRICT_JMLE_A2_RESULTS.md").write_text(report, encoding="utf-8")


def run(
    input_dir: Path,
    stage_a_dir: Path,
    stage_a_plan_path: Path,
    output: Path,
    plan_path: Path,
) -> None:
    input_dir = input_dir.resolve()
    stage_a_dir = stage_a_dir.resolve()
    stage_a_plan_path = stage_a_plan_path.resolve()
    output = output.resolve()
    plan_path = plan_path.resolve()
    plan = load_plan(plan_path)
    input_checks = validate_input_identity(input_dir, stage_a_dir, stage_a_plan_path, plan)
    strict_plan = stage_a.load_plan(stage_a_plan_path)
    controls = strict_plan["numerical_controls"]
    manifest = pd.read_csv(input_dir / "manifest.csv")
    ratings = pd.read_csv(input_dir / "generated_ratings.csv")
    truth = pd.read_csv(input_dir / "generated_facet_truth.csv")
    anchors = pd.read_csv(input_dir / "generated_anchors.csv")
    stored = pd.read_csv(stage_a_dir / "strict_jmle_runs.csv").set_index("RunId", drop=False)
    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_strict_jmle_a2_plan.json").write_bytes(plan_path.read_bytes())

    baseline_rows: list[dict[str, object]] = []
    run_rows: list[dict[str, object]] = []
    fd_parts: list[pd.DataFrame] = []
    for _, manifest_row in manifest.iterrows():
        run_id = str(manifest_row["RunId"])
        generated = stage_a.generated_design_for_run(manifest_row, ratings, truth, anchors)
        base_result, base_elapsed, base_warnings = fit_stage_a_baseline(manifest_row, generated, controls)
        baseline_row, start, idx, config, sizes, objective, _ = baseline_audit(
            manifest_row,
            base_result,
            base_elapsed,
            base_warnings,
            stored.loc[run_id],
            plan,
        )
        baseline_row["RequestedAnchorRows"] = int(len(generated.anchors))
        baseline_rows.append(baseline_row)
        for candidate in plan["candidates"]:
            candidate_row, fd_rows = polish_one(
                baseline_row,
                start,
                idx,
                config,
                sizes,
                objective,
                generated.anchors,
                candidate,
                plan,
            )
            run_rows.append(candidate_row)
            if not fd_rows.empty:
                fd_parts.append(fd_rows)

    baseline = pd.DataFrame(baseline_rows)
    runs = pd.DataFrame(run_rows)
    finite_difference = pd.concat(fd_parts, ignore_index=True) if fd_parts else pd.DataFrame()
    gates = candidate_gates(runs, baseline, input_checks, plan)
    selection = selection_table(gates, plan)
    selected = str(selection.loc[selection["Decision"].eq("RevisedStageB2"), "SelectedCandidate"].iloc[0])
    authorized = bool(selection.loc[selection["Decision"].eq("RevisedStageB2"), "Authorized"].iloc[0])
    first_read = pd.DataFrame([
        {
            "Priority": 1,
            "Check": "Stage-A2 solver qualification",
            "Status": "Pass" if authorized else "Blocked",
            "Evidence": f"selected={selected or 'none'}; original Stage B remains unauthorized",
            "NextAction": (
                "Run distinct Stage B2 with the unchanged selected candidate."
                if authorized else
                "Stop; diagnose without changing the frozen A2 thresholds or options."
            ),
            "PublicSurfaceEnabled": False,
        },
        {
            "Priority": 2,
            "Check": "Application integration",
            "Status": "Withheld",
            "Evidence": "A2 evaluates solver behavior only and authorizes no core change.",
            "NextAction": "Require a separate implementation, regression, and user-facing evidence gate.",
            "PublicSurfaceEnabled": False,
        },
        {
            "Priority": 3,
            "Check": "Statistical evidence",
            "Status": "Not assessed",
            "Evidence": "16 smoke RunIds cannot estimate operating characteristics.",
            "NextAction": "Preserve sparse availability and raw/display threshold audits in Stage B2.",
            "PublicSurfaceEnabled": False,
        },
    ])
    baseline.to_csv(output / "strict_jmle_a2_baseline.csv", index=False)
    runs.to_csv(output / "strict_jmle_a2_runs.csv", index=False)
    finite_difference.to_csv(output / "strict_jmle_a2_finite_difference_gradient.csv", index=False)
    gates.to_csv(output / "strict_jmle_a2_candidate_gates.csv", index=False)
    selection.to_csv(output / "strict_jmle_a2_selection.csv", index=False)
    input_checks.to_csv(output / "strict_jmle_a2_input_checks.csv", index=False)
    first_read.to_csv(output / "strict_jmle_a2_first_read_summary.csv", index=False)
    plot_gradient_comparison(
        runs,
        output,
        float(plan["qualification"]["terminal_gradient_sup_norm_threshold"]),
    )
    runs.attrs["plan_candidates"] = plan["candidates"]
    plot_effort(runs, output)
    write_report(output, baseline, runs, gates, selection, plan)

    artifact_names = [
        "strict_jmle_a2_baseline.csv",
        "strict_jmle_a2_runs.csv",
        "strict_jmle_a2_finite_difference_gradient.csv",
        "strict_jmle_a2_candidate_gates.csv",
        "strict_jmle_a2_selection.csv",
        "strict_jmle_a2_input_checks.csv",
        "strict_jmle_a2_first_read_summary.csv",
        "strict_jmle_a2_gradient_comparison.png",
        "strict_jmle_a2_effort_and_gain.png",
        "STRICT_JMLE_A2_RESULTS.md",
    ]
    identity = {
        "schema_version": PLAN_SCHEMA,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "original_stage_b_authorized": False,
        "revised_stage_b2_authorized": authorized,
        "selected_candidate": selected,
        "a2_plan_sha256": stage_a.sha256_file(plan_path),
        "a2_adapter_sha256": stage_a.sha256_file(Path(__file__).resolve()),
        "streamlit_app_sha256": stage_a.sha256_file(REPO / "streamlit_app.py"),
        "artifacts": {name: stage_a.sha256_file(output / name) for name in artifact_names},
    }
    (output / "strict_jmle_a2_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--stage-a", type=Path, default=DEFAULT_STAGE_A)
    parser.add_argument("--stage-a-plan", type=Path, default=DEFAULT_STAGE_A_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    args = parser.parse_args()
    run(args.input, args.stage_a, args.stage_a_plan, args.output, args.plan)


if __name__ == "__main__":
    main()
