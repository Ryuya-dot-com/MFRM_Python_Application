#!/usr/bin/env python3
"""Assess the frozen 20-replicate Python pilot without promoting its outcomes.

Budget calculations use only runtime and attempted/returned/converged/eligible
accounting allowed by the registered precision plan. Decision rates, recovery,
and coverage are reported only as explicitly unstable pilot diagnostics.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from mfrm_app import operating_characteristics as oc  # noqa: E402


DEFAULT_INPUT = REPO / "validation" / "operating_characteristics_pilot20_20260809"
DEFAULT_PLAN = REPO / "validation" / "operating_characteristics_precision_plan_20260809.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series.dtype):
        return series.fillna(False).astype(bool)
    return series.astype("string").str.strip().str.lower().map({
        "true": True,
        "1": True,
        "yes": True,
        "false": False,
        "0": False,
        "no": False,
    }).fillna(False).astype(bool)


def build_budget_assessment(runs: pd.DataFrame, plan: dict) -> pd.DataFrame:
    """Return condition budgets from permitted pilot accounting fields only."""

    required = {
        "ConditionId", "Design", "TruthBias", "RunId", "FitReturned",
        "Converged", "InferenceReady", "AnalysisEligible",
        "StrictAnalysisEligible1e4", "ElapsedSeconds",
    }
    missing = required - set(runs.columns)
    if missing:
        raise ValueError(f"pilot runs missing required columns: {sorted(missing)}")
    target = max(
        int(entry["minimum_eligible_per_condition"])
        for entry in plan["confirmatory_precision"].values()
    )
    cap = int(plan["budget_rule"]["maximum_attempts_per_condition"])
    rows: list[dict[str, object]] = []
    groups = ["ConditionId", "Design", "TruthBias"]
    for keys, part in runs.groupby(groups, sort=False, dropna=False):
        identity = dict(zip(groups, keys))
        attempts = int(len(part))
        returned = int(_bool(part["FitReturned"]).sum())
        converged = int(_bool(part["Converged"]).sum())
        eligible = int(_bool(part["AnalysisEligible"]).sum())
        strict_eligible = int(_bool(part["StrictAnalysisEligible1e4"]).sum())
        current = oc.plan_confirmatory_attempts(
            eligible,
            attempts,
            target_eligible=target,
            maximum_attempts=cap,
        )
        strict = oc.plan_confirmatory_attempts(
            strict_eligible,
            attempts,
            target_eligible=target,
            maximum_attempts=cap,
        )
        elapsed = pd.to_numeric(part["ElapsedSeconds"], errors="coerce")
        elapsed = elapsed.loc[elapsed.notna() & np.isfinite(elapsed) & (elapsed >= 0)]
        mean_seconds = float(elapsed.mean()) if len(elapsed) else np.nan
        current_attempts = current["RequiredFixedAttempts"]
        projected_hours = (
            mean_seconds * int(current_attempts) / 3600.0
            if current_attempts is not None and np.isfinite(mean_seconds) else
            np.nan
        )
        if strict_eligible == 0:
            disposition = "blocked_pending_strict_jmle_numerical_qualification"
        elif strict["Status"] != "fixed_attempt_budget_available":
            disposition = str(strict["Status"])
        else:
            disposition = "eligible_for_fixed_manifest_review_not_promotion"
        rows.append({
            **identity,
            "PilotAttempts": attempts,
            "FitsReturned": returned,
            "ConvergedOptimizerFlag": converged,
            "RegisteredAnalysisEligible": eligible,
            "RegisteredEligibilityRate": current["ObservedEligibilityRate"],
            "RegisteredEligibilityWilsonLower95": current["EligibilityWilsonLower95"],
            "RegisteredEligibilityWilsonUpper95": current["EligibilityWilsonUpper95"],
            "TargetEligible": target,
            "RequiredFixedAttemptsFromRegisteredEligibility": current_attempts,
            "RegisteredBudgetStatus": current["Status"],
            "StrictGradientEligible1e4": strict_eligible,
            "StrictEligibilityWilsonLower95": strict["EligibilityWilsonLower95"],
            "RequiredFixedAttemptsStrict1e4": strict["RequiredFixedAttempts"],
            "StrictBudgetStatus": strict["Status"],
            "MeanFitSeconds": mean_seconds,
            "MedianFitSeconds": float(elapsed.median()) if len(elapsed) else np.nan,
            "P95FitSeconds": float(elapsed.quantile(0.95)) if len(elapsed) else np.nan,
            "ProjectedSerialFitHoursFromRegisteredEligibility": projected_hours,
            "PlanningDisposition": disposition,
            "PlanningBoundary": (
                "Current-attempt calculations use registered app eligibility only. "
                "The post-pilot 1e-4 sup-norm audit is a blocker sensitivity, not a "
                "retroactively promoted primary rule. Pilot decision outcomes are not used."
            ),
        })
    return pd.DataFrame(rows)


def build_gradient_readiness(runs: pd.DataFrame) -> pd.DataFrame:
    """Summarize terminal-gradient sensitivity without selecting a threshold."""

    thresholds = (1e-4, 1e-5, 1e-6)
    rows: list[dict[str, object]] = []
    group_columns = ["ConditionId", "Design", "TruthBias"]
    grouped = list(runs.groupby(group_columns, sort=False, dropna=False))
    grouped.append((
        ("ALL", "all_conditions", np.nan),
        runs,
    ))
    for keys, part in grouped:
        gradient = pd.to_numeric(part["TerminalGradientSupNorm"], errors="coerce")
        finite = gradient.notna() & np.isfinite(gradient)
        converged = _bool(part["Converged"])
        for threshold in thresholds:
            ready = finite & converged & (gradient <= threshold)
            rows.append({
                "ConditionId": keys[0],
                "Design": keys[1],
                "TruthBias": keys[2],
                "Threshold": threshold,
                "Attempts": int(len(part)),
                "FiniteGradient": int(finite.sum()),
                "OptimizerFlagConverged": int(converged.sum()),
                "GradientReady": int(ready.sum()),
                "GradientReadyRate": float(ready.mean()) if len(part) else np.nan,
                "MinTerminalGradientSupNorm": float(gradient.loc[finite].min()) if finite.any() else np.nan,
                "MedianTerminalGradientSupNorm": float(gradient.loc[finite].median()) if finite.any() else np.nan,
                "MaxTerminalGradientSupNorm": float(gradient.loc[finite].max()) if finite.any() else np.nan,
                "Boundary": (
                    "Sensitivity audit only. Optimizer success and a finite returned fit "
                    "do not by themselves establish numerical stationarity."
                ),
            })
    return pd.DataFrame(rows)


def build_anchor_sensitivity(parameters: pd.DataFrame) -> pd.DataFrame:
    """Compare paired clean and contaminated hard-anchor Rater estimates."""

    selected = parameters.loc[
        parameters["Design"].astype(str).isin(["balanced_large_anchors", "anchor_drift"])
        & parameters["Facet"].astype(str).eq("Rater")
    ].copy()
    keys = ["TruthBias", "Replicate", "Seed", "Facet", "Level", "Truth"]
    clean = selected.loc[selected["Design"].eq("balanced_large_anchors")]
    drift = selected.loc[selected["Design"].eq("anchor_drift")]
    paired = clean.merge(
        drift,
        on=keys,
        how="inner",
        suffixes=("Clean", "Drift"),
        validate="one_to_one",
    )
    paired["EstimateShiftDriftMinusClean"] = (
        pd.to_numeric(paired["EstimateDrift"], errors="coerce")
        - pd.to_numeric(paired["EstimateClean"], errors="coerce")
    )
    paired["Anchored"] = _bool(paired["AnchoredClean"])
    paired["InputAnchorShift"] = np.where(paired["Anchored"], 0.25, 0.0)
    paired["ShiftErrorAgainstInput"] = (
        paired["EstimateShiftDriftMinusClean"] - paired["InputAnchorShift"]
    )
    return paired[[
        "TruthBias", "Replicate", "Seed", "Facet", "Level", "Truth",
        "Anchored", "EstimateClean", "EstimateDrift",
        "EstimateShiftDriftMinusClean", "InputAnchorShift",
        "ShiftErrorAgainstInput",
    ]]


def integrity_checks(input_dir: Path, plan_path: Path) -> pd.DataFrame:
    """Audit retained hashes, row identities, and extension evidence."""

    identity = json.loads((input_dir / "study_identity.json").read_text(encoding="utf-8"))
    extension = json.loads((input_dir / "manifest_extension_audit.json").read_text(encoding="utf-8"))
    manifest = pd.read_csv(input_dir / "manifest.csv")
    runs = pd.read_csv(input_dir / "runs.csv")
    bridge = pd.read_csv(input_dir / "bridge_validation_python.csv")
    checks: list[dict[str, object]] = []

    def add(check: str, passed: bool, evidence: str) -> None:
        checks.append({"Check": check, "Passed": bool(passed), "Evidence": evidence})

    source_plan_hash = sha256_file(plan_path)
    retained_plan_hash = sha256_file(input_dir / "registered_precision_plan.json")
    oc.load_precision_plan(plan_path)
    oc.load_precision_plan(input_dir / "registered_precision_plan.json")
    add(
        "precision_plan_hash",
        source_plan_hash == retained_plan_hash == identity.get("precision_plan_sha256"),
        f"source={source_plan_hash}; retained={retained_plan_hash}",
    )
    manifest_hash = oc.frame_fingerprint(manifest)
    add(
        "manifest_hash",
        manifest_hash == identity.get("manifest_sha256"),
        f"computed={manifest_hash}; recorded={identity.get('manifest_sha256')}",
    )
    add(
        "run_manifest_identity",
        set(manifest["RunId"].astype(str)) == set(runs["RunId"].astype(str))
        and len(manifest) == len(runs) == 160,
        f"manifest={len(manifest)}; runs={len(runs)}; unique_runs={runs['RunId'].nunique()}",
    )
    counts = manifest.groupby("ConditionId")["Replicate"].nunique()
    add(
        "registered_pilot_depth",
        len(counts) == 8 and counts.eq(20).all() and not bool(identity.get("profile_replicate_override")),
        f"conditions={len(counts)}; min/max replicates={counts.min()}/{counts.max()}",
    )
    add(
        "manifest_extension",
        bool(extension.get("Passed"))
        and extension.get("PriorManifestSHA256") == extension.get("RetainedRowsSHA256")
        and extension.get("ExpandedManifestSHA256") == manifest_hash,
        f"prior={extension.get('PriorRuns')}; added={extension.get('AddedRuns')}; expanded={extension.get('ExpandedRuns')}",
    )
    add(
        "generated_bundle_validation",
        not bridge.empty and _bool(bridge["Passed"]).all(),
        f"passed={int(_bool(bridge['Passed']).sum())}/{len(bridge)}",
    )
    source_files = {
        "runner_sha256": Path(__file__).with_name("operating_characteristics_pilot.py"),
        "streamlit_app_sha256": REPO / "streamlit_app.py",
        "operating_characteristics_module_sha256": REPO / "mfrm_app" / "operating_characteristics.py",
        "decision_stability_module_sha256": REPO / "mfrm_app" / "decision_stability.py",
    }
    for key, path in source_files.items():
        computed = sha256_file(path)
        add(key, computed == identity.get(key), f"computed={computed}; recorded={identity.get(key)}")
    return pd.DataFrame(checks)


def build_first_read(
    checks: pd.DataFrame,
    runs: pd.DataFrame,
    budget: pd.DataFrame,
    fit_audit: pd.DataFrame,
    anchor: pd.DataFrame,
) -> pd.DataFrame:
    attempts = int(len(runs))
    returned = int(_bool(runs["FitReturned"]).sum())
    converged = int(_bool(runs["Converged"]).sum())
    eligible = int(_bool(runs["AnalysisEligible"]).sum())
    strict = int(_bool(runs["StrictAnalysisEligible1e4"]).sum())
    mismatch = int((~_bool(fit_audit["DisplayDecisionConsistent"])).sum())
    sparse = budget.loc[budget["Design"].eq("sparse_missing")]
    sparse_eligible = " and ".join(
        f"{int(value)}/{int(attempts)}"
        for value, attempts in zip(
            sparse["RegisteredAnalysisEligible"],
            sparse["PilotAttempts"],
        )
    )
    anchored = anchor.loc[anchor["Anchored"]]
    unanchored = anchor.loc[~anchor["Anchored"]]
    rows = [
        (1, "Artifact integrity", "Pass" if checks["Passed"].all() else "Fail",
         f"{int(checks['Passed'].sum())}/{len(checks)} identity and bundle checks passed.",
         "Resolve any failed hash or identity check before interpretation."),
        (2, "Evidence scope", "Pilot only",
         "20 replicates per condition; pilot outcomes are excluded from confirmatory performance claims.",
         "Use runtime and failure/eligibility accounting only for the next budget decision."),
        (3, "Numerical readiness", "Blocked",
         f"{returned}/{attempts} returned and {converged}/{attempts} optimizer flags converged, but {strict}/{attempts} passed the post-pilot 1e-4 sup-norm eligibility sensitivity.",
         "Qualify a stricter Python JMLE numerical mode before fixing a confirmatory manifest."),
        (4, "Sparse design", "Redesign",
         f"The two sparse conditions produced {sparse_eligible} eligible focal decisions; their registered-budget inflation exceeds the 2,500-attempt cap.",
         "Change overlap/cell information or the estimand; do not recode sparse cells as negative decisions."),
        (5, "Registered budget", "Not yet fixable",
         "Non-sparse app-eligibility calculations imply roughly 597-655 attempts/condition, but strict numerical eligibility is zero and overrides scheduling.",
         "Rerun a separately labelled strict numerical pilot, then recompute from eligibility only."),
        (6, "Floating-point/display thresholds", "Review",
         f"{mismatch} fit classifications differ between unrounded values and 3-decimal display; the stored raw classification remains authoritative.",
         "Show the raw value and boundary label in any future one-click result surface."),
        (7, "Anchor contamination", "Sensitivity visible",
         f"Anchored Raters shifted {anchored['EstimateShiftDriftMinusClean'].mean():.6f} logits on average for a +0.25 input; unanchored Raters averaged {unanchored['EstimateShiftDriftMinusClean'].mean():.6f}.",
         "Describe this as hard-anchor transmission, not robustness or an anchor-share recommendation."),
        (8, "Bias operating characteristics", "Too imprecise",
         "Decision-rate Wilson intervals remain wide, especially with only four eligible sparse decisions.",
         "Do not quote pilot false-positive rate or power as performance evidence."),
        (9, "Cross-engine depth", "Pending",
         "mfrmr/TAM/immer/sirt adapters remain smoke-depth; no 20-replicate cross-engine performance run was executed.",
         "After Python numerical qualification, extend only matched estimands and keep MML sensitivities separate."),
        (10, "Public application surface", "Withheld",
         "Repository-only evidence; no automatic sample-size, anchor-share, or package ranking is enabled.",
         "Keep the future button disabled until all protocol promotion gates pass."),
    ]
    return pd.DataFrame([
        {
            "Priority": priority,
            "Check": check,
            "Status": status,
            "Evidence": evidence,
            "NextAction": action,
            "PublicSurfaceEnabled": False,
        }
        for priority, check, status, evidence, action in rows
    ])


def plot_eligibility_funnel(budget: pd.DataFrame, output: Path) -> None:
    labels = [
        f"{row.Design.replace('_', ' ')}\nbias={row.TruthBias:g}"
        for row in budget.itertuples(index=False)
    ]
    x = np.arange(len(budget))
    width = 0.2
    fig, ax = plt.subplots(figsize=(13, 6.8))
    series = [
        ("PilotAttempts", "attempted", "#4C78A8"),
        ("ConvergedOptimizerFlag", "optimizer flag converged", "#59A14F"),
        ("RegisteredAnalysisEligible", "registered app eligible", "#F28E2B"),
        ("StrictGradientEligible1e4", "also sup-norm <= 1e-4", "#E15759"),
    ]
    for offset, (column, label, color) in enumerate(series):
        ax.bar(x + (offset - 1.5) * width, budget[column], width, label=label, color=color)
    ax.set_xticks(x, labels, rotation=28, ha="right")
    ax.set_ylabel("Runs out of 20")
    ax.set_ylim(0, 22)
    ax.set_title("Pilot readiness funnel: optimizer success is not numerical qualification")
    ax.legend(ncol=2, frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "pilot_eligibility_funnel.png", dpi=180)
    plt.close(fig)


def plot_gradients(runs: pd.DataFrame, output: Path) -> None:
    designs = list(dict.fromkeys(runs["Design"].astype(str)))
    fig, ax = plt.subplots(figsize=(10.5, 6.8))
    for index, design in enumerate(designs):
        part = runs.loc[runs["Design"].astype(str).eq(design)].copy()
        values = pd.to_numeric(part["TerminalGradientSupNorm"], errors="coerce")
        finite = values.notna() & np.isfinite(values) & (values > 0)
        offsets = np.linspace(-0.24, 0.24, int(finite.sum())) if finite.any() else np.array([])
        colors = np.where(_bool(part.loc[finite, "Converged"]), "#4C78A8", "#E15759")
        ax.scatter(index + offsets, values.loc[finite], c=colors, alpha=0.72, s=24, edgecolor="none")
    for threshold, style in ((1e-4, "-"), (1e-5, "--"), (1e-6, ":")):
        ax.axhline(threshold, color="#222222", linestyle=style, linewidth=1.1, label=f"{threshold:.0e}")
    ax.set_yscale("log")
    ax.set_xticks(range(len(designs)), [value.replace("_", " ") for value in designs], rotation=18, ha="right")
    ax.set_ylabel("Terminal gradient sup norm (log scale)")
    ax.set_title("All 160 matched-control fits remain above the strict sensitivity bands")
    ax.legend(title="sensitivity threshold", frameon=False)
    ax.grid(axis="y", which="both", alpha=0.22)
    fig.tight_layout()
    fig.savefig(output / "pilot_gradient_readiness.png", dpi=180)
    plt.close(fig)


def plot_float_mismatches(fit_audit: pd.DataFrame, output: Path) -> None:
    mismatches = fit_audit.loc[~_bool(fit_audit["DisplayDecisionConsistent"])].copy()
    if mismatches.empty:
        return
    mismatches["SignedDisplayHalfUnits"] = (
        (pd.to_numeric(mismatches["RawValue"], errors="coerce")
         - pd.to_numeric(mismatches["NearestThreshold"], errors="coerce"))
        / pd.to_numeric(mismatches["DisplayHalfUnit"], errors="coerce")
    )
    labels = [
        f"{row.RunId.split('::')[0]} rep {int(row.Replicate)}\n{row.Level} {row.Statistic}: {row.RawDecision} -> displayed {row.DisplayDecision}"
        for row in mismatches.itertuples(index=False)
    ]
    y = np.arange(len(mismatches))
    fig, ax = plt.subplots(figsize=(11.5, 6.5))
    values = mismatches["SignedDisplayHalfUnits"].to_numpy(dtype=float)
    ax.barh(y, values, color=np.where(values < 0, "#4C78A8", "#E15759"))
    ax.axvline(0, color="#222222", linewidth=1.2)
    ax.axvline(-1, color="#777777", linestyle="--", linewidth=1)
    ax.axvline(1, color="#777777", linestyle="--", linewidth=1)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Signed distance from threshold / 3-decimal display half-unit")
    ax.set_title("Five fit labels would change if the displayed value drove classification")
    ax.grid(axis="x", alpha=0.22)
    fig.tight_layout()
    fig.savefig(output / "pilot_float_boundary_mismatches.png", dpi=180)
    plt.close(fig)


def plot_anchor_sensitivity(anchor: pd.DataFrame, output: Path) -> None:
    levels = list(dict.fromkeys(anchor["Level"].astype(str)))
    fig, ax = plt.subplots(figsize=(9.5, 6.4))
    for index, level in enumerate(levels):
        values = pd.to_numeric(
            anchor.loc[anchor["Level"].astype(str).eq(level), "EstimateShiftDriftMinusClean"],
            errors="coerce",
        ).dropna()
        offsets = np.linspace(-0.18, 0.18, len(values))
        color = "#E15759" if bool(anchor.loc[anchor["Level"].eq(level), "Anchored"].iloc[0]) else "#4C78A8"
        ax.scatter(index + offsets, values, color=color, alpha=0.55, s=20)
        ax.plot([index - 0.25, index + 0.25], [values.mean(), values.mean()], color="#222222", linewidth=2)
    ax.axhline(0.25, color="#E15759", linestyle="--", linewidth=1, label="injected anchor shift (+0.25)")
    ax.axhline(0.0, color="#4C78A8", linestyle=":", linewidth=1, label="no input shift")
    ax.set_xticks(range(len(levels)), levels)
    ax.set_ylabel("Estimate shift: contaminated minus clean anchors (logits)")
    ax.set_title("Hard-anchor contamination is transmitted exactly to fixed Raters")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.22)
    fig.tight_layout()
    fig.savefig(output / "pilot_anchor_contamination.png", dpi=180)
    plt.close(fig)


def write_report(
    output: Path,
    checks: pd.DataFrame,
    runs: pd.DataFrame,
    budget: pd.DataFrame,
    gradient: pd.DataFrame,
    fit_audit: pd.DataFrame,
    anchor: pd.DataFrame,
    binary: pd.DataFrame,
    first_read: pd.DataFrame,
) -> None:
    attempts = len(runs)
    current_eligible = int(_bool(runs["AnalysisEligible"]).sum())
    strict_eligible = int(_bool(runs["StrictAnalysisEligible1e4"]).sum())
    all_gradient = gradient.loc[
        gradient["ConditionId"].eq("ALL") & gradient["Threshold"].eq(1e-4)
    ].iloc[0]
    mismatches = fit_audit.loc[~_bool(fit_audit["DisplayDecisionConsistent"])]
    sparse = budget.loc[budget["Design"].eq("sparse_missing")]
    ordinary = budget.loc[~budget["Design"].eq("sparse_missing")]
    current_attempts = pd.to_numeric(
        ordinary["RequiredFixedAttemptsFromRegisteredEligibility"], errors="coerce"
    ).dropna()
    anchored = anchor.loc[anchor["Anchored"], "EstimateShiftDriftMinusClean"]
    unanchored = anchor.loc[~anchor["Anchored"], "EstimateShiftDriftMinusClean"]
    strong = binary.loc[binary["DecisionRule"].eq("DecisionStrongRaw")]
    widest = float(
        (pd.to_numeric(strong["DecisionRateWilsonUpper95"], errors="coerce")
         - pd.to_numeric(strong["DecisionRateWilsonLower95"], errors="coerce")).max()
    )
    first_lines = "\n".join(
        f"- **{row.Check} — {row.Status}:** {row.Evidence}"
        for row in first_read.sort_values("Priority").itertuples(index=False)
    )
    report = f"""# Python 20-replicate pilot assessment

## Decision

Do **not** advance the current matched-control Python JMLE lane directly to the
100-replicate screening profile or a confirmatory run. The artifact is intact
and the orchestration scales, but numerical readiness and sparse-cell
eligibility fail before Monte Carlo precision becomes the limiting issue.

## First read

{first_lines}

## What the pilot established

- All {attempts} fits returned in {runs['ElapsedSeconds'].sum():.1f} recorded fit-seconds; {int(_bool(runs['Converged']).sum())}/{attempts} carried the optimizer convergence flag and {current_eligible}/{attempts} met the registered app eligibility rule.
- All {int(checks['Passed'].sum())}/{len(checks)} retained hash, manifest-extension, RunId, profile-depth, and generated-bundle checks passed.
- The original 16 smoke rows are unchanged inside the 160-row pilot manifest; 144 consecutive replicate rows were added.
- Non-sparse eligibility alone would imply {int(current_attempts.min())}-{int(current_attempts.max())} fixed attempts per condition to target 500 eligible results using the registered Wilson-lower-bound inflation rule.

## Why confirmatory scheduling is blocked

The app's `Converged` flag records optimizer termination, not a separately
qualified stationarity rule. The reconstructed terminal gradient sup norm
ranged from {all_gradient.MinTerminalGradientSupNorm:.6g} to
{all_gradient.MaxTerminalGradientSupNorm:.6g}. Consequently,
{strict_eligible}/{attempts} runs meet the post-pilot `1e-4` sup-norm plus
analysis-eligibility sensitivity (and none meet `1e-5` or `1e-6`). This
threshold is not retroactively substituted for the registered primary rule;
it is a blocker demonstrating that a strict Python numerical mode must be
qualified and prospectively frozen first.

The sparse null and alternative conditions each yielded only
{int(sparse.iloc[0].RegisteredAnalysisEligible)}/20 and
{int(sparse.iloc[1].RegisteredAnalysisEligible)}/20 eligible focal decisions.
Their Wilson-lower-bound attempt inflation exceeds the registered 2,500-run
cap. More brute-force replicates would estimate a mostly unavailable estimand;
the design or focal-cell estimand needs revision.

## Floating-point decision audit

The row-level fit audit contains {len(fit_audit):,} statistic evaluations,
{int(fit_audit['BoundaryStatus'].eq('display_rounding_boundary').sum())}
display-boundary rows, and {len(mismatches)} raw/display classification
mismatches. Four raw Infit/Outfit values just below `0.50` would display as
`0.500` and appear acceptable; one raw value `1.500475` would display as
`1.500` and hide a noisy classification. Raw finite values remain
authoritative. The focal bias screen had no raw/display decision mismatch in
this pilot, but that absence is not a general stability claim.

## Anchors and bias-screen scope

Because clean and contaminated anchor runs reuse identical rating bytes, the
paired contrast isolates the supplied hard-anchor shift. Anchored Raters moved
{anchored.mean():.6f} logits on average (range {anchored.min():.6f} to
{anchored.max():.6f}) under the injected `+0.25`; unanchored Raters averaged
{unanchored.mean():.6f} (range {unanchored.min():.6f} to
{unanchored.max():.6f}). This is constraint transmission, not robustness.

The pilot's strongest-decision Wilson intervals have a maximum full width of
{widest:.3f}. Decision rates, false-positive rates, power, RMSE, and coverage
remain debugging diagnostics and were not used to calculate the budget.

## Next registered action

1. Define a separately labelled strict Python JMLE mode with prospective
   max-iteration, relative-tolerance, terminal sup-norm, and failure rules.
2. Run a small numerical-qualification extension without pooling its outcomes
   with this pilot; require materially improved stationarity and retain all
   failures.
3. Redesign the sparse focal cell or explicitly declare that bias detection is
   unavailable there; do not solve it by counting sparse runs as negatives.
4. Only then freeze a fixed confirmatory manifest. The frozen v1 precision
   plan's continuous-bias section also conflicts with its forbidden pilot
   outcome list, so continuous-outcome precision needs a prospective amendment
   or an independent variance-planning source before use.
5. Keep the public one-click simulation/cross-engine surface `Withheld`.

## Retained assessment artifacts

- `pilot_integrity_checks.csv`
- `pilot_budget_assessment.csv`
- `pilot_gradient_readiness.csv`
- `pilot_anchor_sensitivity.csv`
- `pilot_first_read_summary.csv`
- `pilot_eligibility_funnel.png`
- `pilot_gradient_readiness.png`
- `pilot_float_boundary_mismatches.png`
- `pilot_anchor_contamination.png`

The four charts are decision aids for repository review, not public validation
figures or automated design recommendations.
"""
    (output / "PILOT20_ASSESSMENT.md").write_text(report, encoding="utf-8")


def run(input_dir: Path, plan_path: Path) -> None:
    input_dir = input_dir.resolve()
    plan_path = plan_path.resolve()
    plan = oc.load_precision_plan(plan_path)
    runs = pd.read_csv(input_dir / "runs.csv")
    parameters = pd.read_csv(input_dir / "parameter_recovery.csv")
    fit_audit = pd.read_csv(input_dir / "fit_decision_stability.csv")
    binary = pd.read_csv(input_dir / "binary_operating_characteristics.csv")
    checks = integrity_checks(input_dir, plan_path)
    if not checks["Passed"].all():
        failed = checks.loc[~checks["Passed"], "Check"].astype(str).tolist()
        raise ValueError(f"pilot integrity checks failed: {failed}")
    budget = build_budget_assessment(runs, plan)
    gradient = build_gradient_readiness(runs)
    anchor = build_anchor_sensitivity(parameters)
    first_read = build_first_read(checks, runs, budget, fit_audit, anchor)
    checks.to_csv(input_dir / "pilot_integrity_checks.csv", index=False)
    budget.to_csv(input_dir / "pilot_budget_assessment.csv", index=False)
    gradient.to_csv(input_dir / "pilot_gradient_readiness.csv", index=False)
    anchor.to_csv(input_dir / "pilot_anchor_sensitivity.csv", index=False)
    first_read.to_csv(input_dir / "pilot_first_read_summary.csv", index=False)
    plot_eligibility_funnel(budget, input_dir)
    plot_gradients(runs, input_dir)
    plot_float_mismatches(fit_audit, input_dir)
    plot_anchor_sensitivity(anchor, input_dir)
    write_report(
        input_dir,
        checks,
        runs,
        budget,
        gradient,
        fit_audit,
        anchor,
        binary,
        first_read,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--precision-plan", type=Path, default=DEFAULT_PLAN)
    args = parser.parse_args()
    run(args.input, args.precision_plan)


if __name__ == "__main__":
    main()
