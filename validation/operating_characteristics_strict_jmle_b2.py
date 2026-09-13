#!/usr/bin/env python3
"""Run frozen 160-fit Stage-B2 precision JMLE numerical qualification."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import streamlit_app as app  # noqa: E402
from validation import operating_characteristics_pilot as pilot  # noqa: E402
from validation import operating_characteristics_strict_jmle as stage_a  # noqa: E402
from validation import operating_characteristics_strict_jmle_a2 as stage_a2  # noqa: E402


PLAN_SCHEMA = "mfrm-strict-jmle-b2-numerical-plan-v1"
DEFAULT_INPUT = REPO / "validation" / "operating_characteristics_pilot20_20260809"
DEFAULT_STAGE_A_PLAN = REPO / "validation" / "operating_characteristics_strict_jmle_plan_20260809.json"
DEFAULT_A2 = REPO / "validation" / "operating_characteristics_strict_jmle_a2_smoke_20260809"
DEFAULT_PLAN = REPO / "validation" / "operating_characteristics_strict_jmle_b2_plan_20260809.json"
DEFAULT_OUTPUT = REPO / "validation" / "operating_characteristics_strict_jmle_b2_20260809"


def load_plan(path: Path) -> dict:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != PLAN_SCHEMA:
        raise ValueError(f"Stage-B2 plan schema must be {PLAN_SCHEMA}")
    baseline = plan.get("strict_baseline", {})
    polish = plan.get("selected_polish", {})
    qualification = plan.get("qualification", {})
    floating = plan.get("floating_point_audit", {})
    if (int(baseline.get("maxit", 0)), float(baseline.get("reltol", np.nan))) != (500, 1e-9):
        raise ValueError("Stage-B2 strict baseline controls changed")
    if polish.get("candidate_id") != "lbfgsb_precision" or polish.get("method") != "L-BFGS-B":
        raise ValueError("Stage-B2 selected candidate changed")
    expected_options = {"maxiter": 1500, "gtol": 1e-7, "ftol": 1e-15, "maxls": 50, "maxcor": 20}
    if polish.get("options") != expected_options:
        raise ValueError("Stage-B2 polish options changed")
    if int(qualification.get("attempted_fits", 0)) != 160:
        raise ValueError("Stage-B2 must retain 160 attempted fits")
    if float(qualification.get("terminal_gradient_sup_norm_threshold", np.nan)) != 1e-4:
        raise ValueError("Stage-B2 terminal-gradient threshold must remain 1e-4")
    if int(qualification.get("anchor_run_count", 0)) != 80:
        raise ValueError("Stage-B2 anchor-run count must remain 80")
    if int(qualification.get("finite_difference_coordinates_per_fit", 0)) != 8:
        raise ValueError("Stage-B2 finite-difference coordinate count must remain 8")
    if floating.get("decision_source") != "unrounded raw terminal gradient sup norm":
        raise ValueError("Stage-B2 decisions must use unrounded raw gradients")
    if [float(value) for value in floating.get("sensitivity_thresholds", [])] != [1e-4, 1e-5, 1e-6]:
        raise ValueError("Stage-B2 sensitivity thresholds changed")
    boundaries = plan.get("scope_boundaries", {})
    if boundaries.get("public_surface_status") != "Withheld":
        raise ValueError("Stage-B2 public surface must remain Withheld")
    if bool(boundaries.get("application_core_change_authorized_by_this_plan", True)):
        raise ValueError("Stage-B2 cannot authorize an application-core change")
    return plan


def validate_input_identity(
    input_dir: Path,
    a2_dir: Path,
    stage_a_plan_path: Path,
    plan: dict,
) -> pd.DataFrame:
    expected = plan["input_identity"]
    paths = {
        "streamlit_app_sha256": REPO / "streamlit_app.py",
        "stage_a_plan_sha256": stage_a_plan_path,
        "stage_a_adapter_sha256": Path(stage_a.__file__).resolve(),
        "stage_a2_plan_sha256": REPO / "validation" / "operating_characteristics_strict_jmle_a2_plan_20260809.json",
        "stage_a2_adapter_sha256": Path(stage_a2.__file__).resolve(),
        "stage_a2_selection_sha256": a2_dir / "strict_jmle_a2_selection.csv",
        "stage_a2_identity_sha256": a2_dir / "strict_jmle_a2_identity.json",
        "pilot_manifest_file_sha256": input_dir / "manifest.csv",
        "pilot_generated_ratings_sha256": input_dir / "generated_ratings.csv",
        "pilot_generated_facet_truth_sha256": input_dir / "generated_facet_truth.csv",
        "pilot_generated_anchors_sha256": input_dir / "generated_anchors.csv",
        "pilot_runs_sha256": input_dir / "runs.csv",
        "pilot_study_identity_sha256": input_dir / "study_identity.json",
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

    a2_selection = pd.read_csv(a2_dir / "strict_jmle_a2_selection.csv")
    revised = a2_selection.loc[a2_selection["Decision"].eq("RevisedStageB2")].iloc[0]
    selection_passed = bool(revised["Authorized"]) and str(revised["SelectedCandidate"]) == "lbfgsb_precision"
    rows.append({
        "Check": "stage_a2_authorization_chain",
        "Passed": selection_passed,
        "Evidence": f"authorized={revised['Authorized']}; selected={revised['SelectedCandidate']}",
    })
    a2_identity = json.loads((a2_dir / "strict_jmle_a2_identity.json").read_text(encoding="utf-8"))
    identity_passed = bool(
        a2_identity.get("original_stage_b_authorized") is False
        and a2_identity.get("revised_stage_b2_authorized") is True
        and a2_identity.get("selected_candidate") == "lbfgsb_precision"
    )
    rows.append({
        "Check": "stage_a2_identity_decision",
        "Passed": identity_passed,
        "Evidence": (
            f"original={a2_identity.get('original_stage_b_authorized')}; "
            f"revised={a2_identity.get('revised_stage_b2_authorized')}; "
            f"selected={a2_identity.get('selected_candidate')}"
        ),
    })
    manifest = pd.read_csv(input_dir / "manifest.csv")
    manifest_passed = bool(
        len(manifest) == 160
        and manifest["ConditionId"].nunique() == 8
        and manifest.groupby("ConditionId")["Replicate"].nunique().eq(20).all()
    )
    rows.append({
        "Check": "pilot_manifest_shape",
        "Passed": manifest_passed,
        "Evidence": (
            f"runs={len(manifest)}; conditions={manifest['ConditionId'].nunique()}; "
            f"replicate_range={manifest['Replicate'].min()}-{manifest['Replicate'].max()}"
        ),
    })
    generated_checks = pilot.validate_generated_bundle(input_dir)
    rows.append({
        "Check": "pilot_generated_bundle_validation",
        "Passed": bool(generated_checks["Passed"].all()),
        "Evidence": f"passed={int(generated_checks['Passed'].sum())}/{len(generated_checks)}",
    })
    checks = pd.DataFrame(rows)
    if not checks["Passed"].all():
        failed = checks.loc[~checks["Passed"], "Check"].astype(str).tolist()
        raise ValueError(f"Stage-B2 input identity failed: {failed}")
    return checks


def strict_baseline_audit(
    manifest_row: pd.Series,
    result: dict,
    elapsed: float,
    warning_text: str,
    pilot_row: pd.Series,
) -> tuple[dict[str, object], np.ndarray, dict, dict, dict, float]:
    parameters, idx, config, sizes = stage_a2.optimizer_context(result)
    objective, gradient = app.mfrm_loglik_jmle_value_grad(parameters, idx, config, sizes)
    gradient = np.asarray(gradient, dtype=float)
    gradient_sup = float(np.max(np.abs(gradient))) if gradient.size else 0.0
    optimizer = result["opt"]
    loglik = float(result["summary"].iloc[0]["LogLik"])
    ready = bool(
        np.isfinite(parameters).all()
        and np.isfinite(objective)
        and np.isfinite(gradient).all()
    )
    row = {
        "RunId": str(manifest_row["RunId"]),
        "ConditionId": str(manifest_row["ConditionId"]),
        "Design": str(manifest_row["Design"]),
        "TruthBias": float(manifest_row["TruthBias"]),
        "Replicate": int(manifest_row["Replicate"]),
        "Rows": int(len(result["prep"]["data"])),
        "OptimizerDimension": int(parameters.size),
        "RequestedAnchorRows": 0,
        "StrictBaselineReady": ready,
        "StrictBaselineOptimizerSuccess": bool(getattr(optimizer, "success", False)),
        "StrictBaselineOptimizerStatus": int(getattr(optimizer, "status", -1)),
        "StrictBaselineOptimizerMessage": str(getattr(optimizer, "message", "")),
        "StrictBaselineIterations": int(getattr(optimizer, "nit", 0)),
        "StrictBaselineFunctionEvaluations": int(getattr(optimizer, "nfev", 0)),
        "StrictBaselineObjective": float(objective),
        "StrictBaselineLogLik": loglik,
        "StrictBaselineTerminalGradientSupNorm": gradient_sup,
        "StrictBaselineElapsedSeconds": elapsed,
        "PilotConverged": bool(pilot_row["Converged"]),
        "PilotLogLik": float(pilot_row["LogLik"]),
        "PilotTerminalGradientSupNorm": float(pilot_row["TerminalGradientSupNorm"]),
        "StrictMinusPilotLogLik": loglik - float(pilot_row["LogLik"]),
        "StrictMinusPilotGradientSupNorm": gradient_sup - float(pilot_row["TerminalGradientSupNorm"]),
        "Warnings": warning_text,
    }
    return row, parameters, idx, config, sizes, float(objective)


def add_floating_point_audit(row: dict[str, object], plan: dict) -> dict[str, object]:
    floating = plan["floating_point_audit"]
    threshold = float(plan["qualification"]["terminal_gradient_sup_norm_threshold"])
    raw = float(row["TerminalGradientSupNorm"])
    display_text = format(raw, str(floating["display_format"]))
    display_value = float(display_text)
    raw_decision = bool(raw <= threshold)
    display_decision = bool(display_value <= threshold)
    row.update({
        "GradientDecisionSource": "raw_unrounded",
        "DisplayedTerminalGradientSupNorm": display_text,
        "DisplayRoundTripTerminalGradientSupNorm": display_value,
        "RawGradientGatePassed": raw_decision,
        "DisplayRoundTripGradientGatePassed": display_decision,
        "RawDisplayDecisionMismatch": raw_decision != display_decision,
        "DistanceRawMinusThreshold": raw - threshold,
        "NearGradientThreshold": abs(raw - threshold) <= float(floating["near_boundary_absolute_band"]),
    })
    for value in floating["sensitivity_thresholds"]:
        label = f"GradientReady{float(value):.0e}".replace("-0", "-").replace("+", "")
        row[label] = bool(raw <= float(value))
    return row


def qualification_gates(
    runs: pd.DataFrame,
    baseline: pd.DataFrame,
    input_checks: pd.DataFrame,
    plan: dict,
) -> pd.DataFrame:
    expected = int(plan["qualification"]["attempted_fits"])
    anchor_runs = runs.loc[
        runs["RunId"].isin(baseline.loc[baseline["RequestedAnchorRows"].gt(0), "RunId"])
    ]
    definitions = [
        ("input_identity", len(input_checks), int(input_checks["Passed"].sum())),
        ("strict_baseline_ready", expected, int(runs["StrictBaselineReady"].astype(bool).sum())),
        ("polish_fit_returned", expected, int(runs["FitReturned"].astype(bool).sum())),
        ("polish_optimizer_success", expected, int(runs["OptimizerSuccess"].astype(bool).sum())),
        ("finite_objective_and_gradient", expected, int(runs["FiniteObjectiveAndGradientGatePassed"].astype(bool).sum())),
        ("terminal_gradient_sup_norm", expected, int(runs["GradientGatePassed"].astype(bool).sum())),
        ("optimizer_jac_recomputed", expected, int(runs["OptimizerJacGatePassed"].astype(bool).sum())),
        ("finite_difference_gradient", expected, int(runs["FiniteDifferenceGatePassed"].astype(bool).sum())),
        ("hard_anchor_constraint", len(anchor_runs), int(anchor_runs["AnchorGatePassed"].astype(bool).sum())),
        ("objective_non_degradation", expected, int(runs["ObjectiveNonDegradationGatePassed"].astype(bool).sum())),
        ("raw_value_decision_source", expected, int(runs["GradientDecisionSource"].eq("raw_unrounded").sum())),
        ("numerical_qualification", expected, int(runs["NumericalQualificationPassed"].astype(bool).sum())),
    ]
    rows = []
    for gate, required, passed in definitions:
        rows.append({
            "Gate": gate,
            "Required": int(required),
            "Passed": int(passed),
            "GatePassed": int(required) == int(passed),
        })
    authorized = all(bool(row["GatePassed"]) for row in rows)
    rows.append({
        "Gate": "integration_experiment_authorization",
        "Required": len(rows),
        "Passed": sum(bool(row["GatePassed"]) for row in rows),
        "GatePassed": authorized,
    })
    return pd.DataFrame(rows)


def condition_summary(runs: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for keys, frame in runs.groupby(["ConditionId", "Design", "TruthBias"], sort=False):
        gradient = pd.to_numeric(frame["TerminalGradientSupNorm"], errors="coerce")
        movement = pd.to_numeric(frame["MaxAbsParameterMovement"], errors="coerce")
        gain = pd.to_numeric(frame["LogLikGainFromBaseline"], errors="coerce")
        rows.append({
            "ConditionId": keys[0],
            "Design": keys[1],
            "TruthBias": keys[2],
            "Attempted": len(frame),
            "Qualified": int(frame["NumericalQualificationPassed"].astype(bool).sum()),
            "QualificationRate": float(frame["NumericalQualificationPassed"].astype(bool).mean()),
            "OptimizerSuccess": int(frame["OptimizerSuccess"].astype(bool).sum()),
            "GradientMedian": float(gradient.median()),
            "GradientP95": float(gradient.quantile(0.95)),
            "GradientMax": float(gradient.max()),
            "MaxAbsParameterMovementMedian": float(movement.median()),
            "MaxAbsParameterMovementP95": float(movement.quantile(0.95)),
            "MaxAbsParameterMovementMax": float(movement.max()),
            "LogLikGainMin": float(gain.min()),
            "LogLikGainMedian": float(gain.median()),
            "LogLikGainMax": float(gain.max()),
            "NearGradientThreshold": int(frame["NearGradientThreshold"].astype(bool).sum()),
            "RawDisplayDecisionMismatches": int(frame["RawDisplayDecisionMismatch"].astype(bool).sum()),
        })
    return pd.DataFrame(rows)


def sensitivity_summary(runs: pd.DataFrame, plan: dict) -> pd.DataFrame:
    frames = [("ALL", runs)]
    frames.extend((str(condition), frame) for condition, frame in runs.groupby("ConditionId", sort=False))
    rows = []
    for condition, frame in frames:
        values = pd.to_numeric(frame["TerminalGradientSupNorm"], errors="coerce")
        for threshold in plan["floating_point_audit"]["sensitivity_thresholds"]:
            threshold = float(threshold)
            passed = values.le(threshold) & values.notna()
            display_values = values.map(lambda value: float(format(float(value), ".3g")) if np.isfinite(value) else np.nan)
            display_passed = display_values.le(threshold) & display_values.notna()
            rows.append({
                "ConditionId": condition,
                "Threshold": threshold,
                "Attempted": len(frame),
                "RawPass": int(passed.sum()),
                "RawPassRate": float(passed.mean()),
                "DisplayRoundTripPass": int(display_passed.sum()),
                "RawDisplayDecisionMismatches": int((passed != display_passed).sum()),
            })
    return pd.DataFrame(rows)


def decision_table(gates: pd.DataFrame) -> pd.DataFrame:
    authorized = bool(
        gates.loc[gates["Gate"].eq("integration_experiment_authorization"), "GatePassed"].iloc[0]
    )
    return pd.DataFrame([
        {
            "Decision": "OriginalStageB",
            "Authorized": False,
            "Evidence": "Stage A remains failed; no retroactive authorization.",
            "PublicSurfaceEnabled": False,
        },
        {
            "Decision": "StageB2NumericalQualification",
            "Authorized": authorized,
            "Evidence": "All frozen Stage-B2 numerical gates passed." if authorized else "At least one frozen Stage-B2 numerical gate failed.",
            "PublicSurfaceEnabled": False,
        },
        {
            "Decision": "IntegrationExperiment",
            "Authorized": authorized,
            "Evidence": "May enter a separate implementation/regression experiment." if authorized else "Integration experiment remains blocked.",
            "PublicSurfaceEnabled": False,
        },
        {
            "Decision": "ApplicationCoreChange",
            "Authorized": False,
            "Evidence": "Requires separate implementation evidence and review.",
            "PublicSurfaceEnabled": False,
        },
        {
            "Decision": "PublicSurface",
            "Authorized": False,
            "Evidence": "Operating characteristics and polished bias decisions were not re-estimated.",
            "PublicSurfaceEnabled": False,
        },
    ])


def plot_gradient_by_condition(runs: pd.DataFrame, output: Path, threshold: float) -> None:
    order = runs["ConditionId"].drop_duplicates().tolist()
    values = [
        pd.to_numeric(runs.loc[runs["ConditionId"].eq(condition), "TerminalGradientSupNorm"], errors="coerce").dropna()
        for condition in order
    ]
    fig, ax = plt.subplots(figsize=(13, 6.5))
    box = ax.boxplot(values, tick_labels=[value.replace("__", "\n") for value in order], showmeans=True, patch_artist=True)
    for patch in box["boxes"]:
        patch.set_facecolor("#A0CBE8")
        patch.set_alpha(0.65)
    ax.axhline(threshold, color="#E15759", linestyle="--", linewidth=1.3, label=f"frozen raw threshold {threshold:.0e}")
    ax.set_yscale("log")
    ax.set_ylabel("Polished terminal gradient sup norm")
    ax.set_title("Stage-B2 precision JMLE across 160 frozen RunIds")
    ax.tick_params(axis="x", rotation=35)
    ax.grid(axis="y", which="both", alpha=0.22)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "strict_jmle_b2_gradient_by_condition.png", dpi=180)
    plt.close(fig)


def plot_baseline_and_movement(runs: pd.DataFrame, output: Path) -> None:
    colors = {
        "balanced_small": "#4E79A7",
        "balanced_large_anchors": "#59A14F",
        "sparse_missing": "#F28E2B",
        "anchor_drift": "#E15759",
    }
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.7))
    for design, frame in runs.groupby("Design", sort=False):
        axes[0].scatter(
            frame["BaselineTerminalGradientSupNorm"],
            frame["TerminalGradientSupNorm"],
            s=28,
            alpha=0.75,
            color=colors.get(str(design), "#777777"),
            label=str(design),
        )
    limits = [1e-8, max(float(runs["BaselineTerminalGradientSupNorm"].max()), 1e-1)]
    axes[0].plot(limits, limits, color="#222222", linestyle=":", linewidth=1)
    axes[0].axhline(1e-4, color="#E15759", linestyle="--", linewidth=1)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlim(limits)
    axes[0].set_xlabel("Strict baseline gradient sup norm")
    axes[0].set_ylabel("Polished gradient sup norm")
    axes[0].set_title("Gradient reduction")
    axes[0].legend(frameon=False, fontsize=8)
    order = runs["Design"].drop_duplicates().tolist()
    movement = [
        pd.to_numeric(runs.loc[runs["Design"].eq(design), "MaxAbsParameterMovement"], errors="coerce").dropna()
        for design in order
    ]
    axes[1].boxplot(movement, tick_labels=[value.replace("_", "\n") for value in order], showmeans=True)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Maximum absolute free-coordinate movement")
    axes[1].set_title("Polish movement (descriptive, not a gate)")
    axes[1].tick_params(axis="x", rotation=20)
    for axis in axes:
        axis.grid(axis="y", which="both", alpha=0.22)
    fig.tight_layout()
    fig.savefig(output / "strict_jmle_b2_gradient_and_movement.png", dpi=180)
    plt.close(fig)


def write_report(
    output: Path,
    baseline: pd.DataFrame,
    runs: pd.DataFrame,
    gates: pd.DataFrame,
    decisions: pd.DataFrame,
    plan: dict,
) -> None:
    authorized = bool(decisions.loc[decisions["Decision"].eq("StageB2NumericalQualification"), "Authorized"].iloc[0])
    gradient = pd.to_numeric(runs["TerminalGradientSupNorm"], errors="coerce")
    movement = pd.to_numeric(runs["MaxAbsParameterMovement"], errors="coerce")
    gain = pd.to_numeric(runs["LogLikGainFromBaseline"], errors="coerce")
    failed = gates.loc[~gates["GatePassed"]]
    failed_text = "\n".join(
        f"- `{row.Gate}`: {int(row.Passed)}/{int(row.Required)}"
        for row in failed.itertuples(index=False)
    ) or "- None."
    baseline_messages = "; ".join(
        f"{count}x {message}"
        for message, count in baseline["StrictBaselineOptimizerMessage"].value_counts().items()
    )
    polish_messages = "; ".join(
        f"{count}x {message}"
        for message, count in runs["OptimizerMessage"].value_counts().items()
    )
    decision = (
        "Stage-B2 numerical qualification passed; a separate integration experiment may begin."
        if authorized else
        "Stage-B2 numerical qualification failed; application integration remains blocked."
    )
    report = f"""# Strict Python JMLE Stage-B2 numerical qualification

## Decision

{decision} This does not authorize an application-core change or a public
performance claim. Original Stage B remains unauthorized.

## Outcome

- {len(runs)}/160 frozen pilot RunIds attempted.
- {int(baseline['StrictBaselineReady'].astype(bool).sum())}/{len(baseline)} strict baselines returned finite objectives and gradients.
- {int(runs['OptimizerSuccess'].astype(bool).sum())}/{len(runs)} precision-polish optimizer success flags.
- {int(runs['NumericalQualificationPassed'].astype(bool).sum())}/{len(runs)} passed every frozen numerical gate.
- Raw polished terminal-gradient sup norm range: {gradient.min():.6g} to {gradient.max():.6g}; frozen gate `1e-4`.
- Maximum absolute free-coordinate movement: {movement.max():.6g} (descriptive only).
- Log-likelihood gain range: {gain.min():.6g} to {gain.max():.6g}.
- Near-boundary raw gradients (within `5e-7` of `1e-4`): {int(runs['NearGradientThreshold'].astype(bool).sum())}.
- Raw/display-roundtrip decision mismatches at 3 significant digits: {int(runs['RawDisplayDecisionMismatch'].astype(bool).sum())}; decisions use raw values.
- Strict-baseline time: {baseline['StrictBaselineElapsedSeconds'].sum():.2f} seconds; polish/verification time: {runs['ElapsedSeconds'].sum():.2f} seconds.

## Termination messages

- Strict baseline: {baseline_messages}.
- Precision polish: {polish_messages}.

## Failed gates

{failed_text}

## Interpretation boundary

This stage qualifies numerical termination on the frozen 160-run pilot bundle.
It does not recompute bias decisions, power, false-positive rate, coverage, or
cross-package equivalence after polishing. The retained pilot used different
baseline iteration controls, so its likelihood and gradient are comparison
fields rather than exact-reproduction gates. Sparse focal cells remain
unavailable rather than negative. The next step, if authorized, is a separate
core-integration experiment that reconstructs all downstream tables from the
polished terminal vector and reruns the statistical and UI decision audits.

## Retained artifacts

- `strict_jmle_b2_baseline.csv`
- `strict_jmle_b2_runs.csv`
- `strict_jmle_b2_finite_difference_gradient.csv`
- `strict_jmle_b2_condition_summary.csv`
- `strict_jmle_b2_gradient_sensitivity.csv`
- `strict_jmle_b2_gates.csv`
- `strict_jmle_b2_decisions.csv`
- `strict_jmle_b2_input_checks.csv`
- `strict_jmle_b2_first_read_summary.csv`
- `strict_jmle_b2_identity.json`
- `strict_jmle_b2_gradient_by_condition.png`
- `strict_jmle_b2_gradient_and_movement.png`
"""
    (output / "STRICT_JMLE_B2_RESULTS.md").write_text(report, encoding="utf-8")


def run(input_dir: Path, a2_dir: Path, stage_a_plan_path: Path, output: Path, plan_path: Path) -> None:
    input_dir = input_dir.resolve()
    a2_dir = a2_dir.resolve()
    stage_a_plan_path = stage_a_plan_path.resolve()
    output = output.resolve()
    plan_path = plan_path.resolve()
    plan = load_plan(plan_path)
    input_checks = validate_input_identity(input_dir, a2_dir, stage_a_plan_path, plan)
    strict_plan = stage_a.load_plan(stage_a_plan_path)
    manifest = pd.read_csv(input_dir / "manifest.csv")
    ratings = pd.read_csv(input_dir / "generated_ratings.csv")
    truth = pd.read_csv(input_dir / "generated_facet_truth.csv")
    anchors = pd.read_csv(input_dir / "generated_anchors.csv")
    pilot_runs = pd.read_csv(input_dir / "runs.csv").set_index("RunId", drop=False)
    candidate = {
        "candidate_id": plan["selected_polish"]["candidate_id"],
        "priority": 1,
        "method": plan["selected_polish"]["method"],
        "options": plan["selected_polish"]["options"],
    }
    polish_plan = {
        "qualification": {
            "finite_difference_coordinates_per_fit": plan["qualification"]["finite_difference_coordinates_per_fit"],
            "hard_anchor_tolerance": plan["qualification"]["hard_anchor_tolerance"],
            "terminal_gradient_sup_norm_threshold": plan["qualification"]["terminal_gradient_sup_norm_threshold"],
            "optimizer_jac_recomputed_max_abs_tolerance": plan["qualification"]["optimizer_jac_recomputed_max_abs_tolerance"],
            "finite_difference_max_abs_error_tolerance": plan["qualification"]["finite_difference_max_abs_error_tolerance"],
            "objective_degradation_tolerance_from_stage_a": plan["qualification"]["objective_degradation_tolerance_from_strict_baseline"],
        }
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_strict_jmle_b2_plan.json").write_bytes(plan_path.read_bytes())
    baseline_rows: list[dict[str, object]] = []
    run_rows: list[dict[str, object]] = []
    fd_parts: list[pd.DataFrame] = []
    started = time.perf_counter()
    for _, manifest_row in manifest.iterrows():
        run_id = str(manifest_row["RunId"])
        generated = stage_a.generated_design_for_run(manifest_row, ratings, truth, anchors)
        base_result, base_elapsed, warning_text = stage_a2.fit_stage_a_baseline(
            manifest_row,
            generated,
            strict_plan["numerical_controls"],
        )
        baseline_row, start, idx, config, sizes, objective = strict_baseline_audit(
            manifest_row,
            base_result,
            base_elapsed,
            warning_text,
            pilot_runs.loc[run_id],
        )
        baseline_row["RequestedAnchorRows"] = int(len(generated.anchors))
        baseline_rows.append(baseline_row)
        proxy = {
            "RunId": baseline_row["RunId"],
            "ConditionId": baseline_row["ConditionId"],
            "Design": baseline_row["Design"],
            "TruthBias": baseline_row["TruthBias"],
            "Replicate": baseline_row["Replicate"],
            "Rows": baseline_row["Rows"],
            "OptimizerDimension": baseline_row["OptimizerDimension"],
            "BaselineReproduced": baseline_row["StrictBaselineReady"],
            "BaselineTerminalGradientSupNorm": baseline_row["StrictBaselineTerminalGradientSupNorm"],
        }
        run_row, fd_rows = stage_a2.polish_one(
            proxy,
            start,
            idx,
            config,
            sizes,
            objective,
            generated.anchors,
            candidate,
            polish_plan,
        )
        run_row["StrictBaselineReady"] = bool(run_row.pop("BaselineReproduced"))
        run_row["FailureReason"] = str(run_row["FailureReason"]).replace(
            "baseline_not_reproduced", "strict_baseline_not_ready"
        )
        run_row["PilotConverged"] = baseline_row["PilotConverged"]
        run_row["PilotTerminalGradientSupNorm"] = baseline_row["PilotTerminalGradientSupNorm"]
        run_row["BaselineTerminalGradientSupNorm"] = baseline_row["StrictBaselineTerminalGradientSupNorm"]
        run_row = add_floating_point_audit(run_row, plan)
        run_rows.append(run_row)
        if not fd_rows.empty:
            fd_parts.append(fd_rows)

    baseline = pd.DataFrame(baseline_rows)
    runs = pd.DataFrame(run_rows)
    finite_difference = pd.concat(fd_parts, ignore_index=True) if fd_parts else pd.DataFrame()
    conditions = condition_summary(runs)
    sensitivity = sensitivity_summary(runs, plan)
    gates = qualification_gates(runs, baseline, input_checks, plan)
    decisions = decision_table(gates)
    authorized = bool(decisions.loc[decisions["Decision"].eq("StageB2NumericalQualification"), "Authorized"].iloc[0])
    first_read = pd.DataFrame([
        {
            "Priority": 1,
            "Check": "Stage-B2 numerical qualification",
            "Status": "Pass" if authorized else "Blocked",
            "Evidence": f"{int(runs['NumericalQualificationPassed'].astype(bool).sum())}/{len(runs)} passed every frozen gate.",
            "NextAction": (
                "Begin a separate application-integration experiment."
                if authorized else
                "Diagnose failures without changing the frozen plan."
            ),
            "PublicSurfaceEnabled": False,
        },
        {
            "Priority": 2,
            "Check": "Floating-point boundary",
            "Status": "Audited",
            "Evidence": (
                f"near={int(runs['NearGradientThreshold'].astype(bool).sum())}; "
                f"raw/display mismatches={int(runs['RawDisplayDecisionMismatch'].astype(bool).sum())}; decisions use raw values"
            ),
            "NextAction": "Expose raw-value basis and boundary warning in any later UI integration.",
            "PublicSurfaceEnabled": False,
        },
        {
            "Priority": 3,
            "Check": "Bias and sparse scope",
            "Status": "Not recomputed",
            "Evidence": "Numerical runner does not reconstruct downstream bias tables after polish.",
            "NextAction": "Recompute all downstream outputs only after core integration is separately authorized.",
            "PublicSurfaceEnabled": False,
        },
        {
            "Priority": 4,
            "Check": "Public application surface",
            "Status": "Withheld",
            "Evidence": "No public performance, equivalence, or package-ranking claim is authorized.",
            "NextAction": "Retain repository-only evidence.",
            "PublicSurfaceEnabled": False,
        },
    ])
    baseline.to_csv(output / "strict_jmle_b2_baseline.csv", index=False)
    runs.to_csv(output / "strict_jmle_b2_runs.csv", index=False)
    finite_difference.to_csv(output / "strict_jmle_b2_finite_difference_gradient.csv", index=False)
    conditions.to_csv(output / "strict_jmle_b2_condition_summary.csv", index=False)
    sensitivity.to_csv(output / "strict_jmle_b2_gradient_sensitivity.csv", index=False)
    gates.to_csv(output / "strict_jmle_b2_gates.csv", index=False)
    decisions.to_csv(output / "strict_jmle_b2_decisions.csv", index=False)
    input_checks.to_csv(output / "strict_jmle_b2_input_checks.csv", index=False)
    first_read.to_csv(output / "strict_jmle_b2_first_read_summary.csv", index=False)
    plot_gradient_by_condition(
        runs,
        output,
        float(plan["qualification"]["terminal_gradient_sup_norm_threshold"]),
    )
    plot_baseline_and_movement(runs, output)
    write_report(output, baseline, runs, gates, decisions, plan)
    artifact_names = [
        "strict_jmle_b2_baseline.csv",
        "strict_jmle_b2_runs.csv",
        "strict_jmle_b2_finite_difference_gradient.csv",
        "strict_jmle_b2_condition_summary.csv",
        "strict_jmle_b2_gradient_sensitivity.csv",
        "strict_jmle_b2_gates.csv",
        "strict_jmle_b2_decisions.csv",
        "strict_jmle_b2_input_checks.csv",
        "strict_jmle_b2_first_read_summary.csv",
        "strict_jmle_b2_gradient_by_condition.png",
        "strict_jmle_b2_gradient_and_movement.png",
        "STRICT_JMLE_B2_RESULTS.md",
    ]
    identity = {
        "schema_version": PLAN_SCHEMA,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "wall_seconds": time.perf_counter() - started,
        "original_stage_b_authorized": False,
        "stage_b2_numerically_qualified": authorized,
        "integration_experiment_authorized": authorized,
        "application_core_change_authorized": False,
        "public_surface_enabled": False,
        "b2_plan_sha256": stage_a.sha256_file(plan_path),
        "b2_adapter_sha256": stage_a.sha256_file(Path(__file__).resolve()),
        "streamlit_app_sha256": stage_a.sha256_file(REPO / "streamlit_app.py"),
        "artifacts": {name: stage_a.sha256_file(output / name) for name in artifact_names},
    }
    (output / "strict_jmle_b2_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--stage-a2", type=Path, default=DEFAULT_A2)
    parser.add_argument("--stage-a-plan", type=Path, default=DEFAULT_STAGE_A_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    args = parser.parse_args()
    run(args.input, args.stage_a2, args.stage_a_plan, args.output, args.plan)


if __name__ == "__main__":
    main()
