#!/usr/bin/env python3
"""Normalize sirt MML sensitivity evidence against the Python OC baseline."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


Q30 = "SIRT_MML_Q30_SENSITIVITY"
Q61 = "SIRT_MML_Q61_PRIMARY"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _quadrature_block(
    frame: pd.DataFrame,
    *,
    parameter_type: str,
    parameter_column: str,
    estimate_column: str,
    inclusion_column: str,
) -> pd.DataFrame:
    keys = ["RunId", parameter_column]
    q30 = frame.loc[frame["Mode"].eq(Q30), keys + [estimate_column, inclusion_column]].rename(
        columns={
            parameter_column: "Parameter",
            estimate_column: "EstimateQ30",
            inclusion_column: "IncludedQ30",
        }
    )
    q61 = frame.loc[frame["Mode"].eq(Q61), keys + [estimate_column, inclusion_column]].rename(
        columns={
            parameter_column: "Parameter",
            estimate_column: "EstimateQ61",
            inclusion_column: "IncludedQ61",
        }
    )
    out = q30.merge(q61, on=["RunId", "Parameter"], how="outer", validate="one_to_one", indicator=True)
    out["ParameterType"] = parameter_type
    out["BothPresent"] = out["_merge"].eq("both")
    out["BothIncluded"] = (
        out["BothPresent"]
        & out["IncludedQ30"].fillna(False).astype(bool)
        & out["IncludedQ61"].fillna(False).astype(bool)
    )
    out["DifferenceQ61MinusQ30"] = out["EstimateQ61"] - out["EstimateQ30"]
    out["AbsDifference"] = out["DifferenceQ61MinusQ30"].abs()
    return out.drop(columns="_merge")


def build_quadrature_pairs(
    rater: pd.DataFrame,
    item: pd.DataFrame,
    person: pd.DataFrame,
    runs: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare q30/q61 without pooling rater, item, and Person estimands."""
    blocks = [
        _quadrature_block(
            rater,
            parameter_type="Rater severity",
            parameter_column="Level",
            estimate_column="EstimateAligned",
            inclusion_column="IncludedInSummary",
        ),
        _quadrature_block(
            item,
            parameter_type="Virtual-item centered location",
            parameter_column="Item",
            estimate_column="EstimateAligned",
            inclusion_column="IncludedInSummary",
        ),
        _quadrature_block(
            item,
            parameter_type="Virtual-item raw location",
            parameter_column="Item",
            estimate_column="Estimate",
            inclusion_column="IncludedInSummary",
        ),
        _quadrature_block(
            person,
            parameter_type="Person EAP",
            parameter_column="Person",
            estimate_column="Estimate",
            inclusion_column="IncludedInSummary",
        ),
    ]
    metadata = runs.loc[runs["Mode"].eq(Q61), [
        "RunId", "ConditionId", "Design", "TruthBias", "Replicate", "Seed"
    ]]
    pairs = pd.concat(blocks, ignore_index=True).merge(
        metadata, on="RunId", how="left", validate="many_to_one"
    )

    rows: list[dict[str, object]] = []
    for (parameter_type, design), group in pairs.groupby(
        ["ParameterType", "Design"], sort=False
    ):
        for scope, mask in [
            ("all_jointly_returned", group["BothPresent"]),
            ("jointly_analysis_eligible", group["BothIncluded"]),
        ]:
            scoped = group.loc[mask]
            run_max = scoped.groupby("RunId")["AbsDifference"].max()
            rows.append(
                {
                    "ParameterType": parameter_type,
                    "Design": design,
                    "Scope": scope,
                    "RunIds": int(scoped["RunId"].nunique()),
                    "Pairs": int(len(scoped)),
                    "MeanAbsDifference": float(scoped["AbsDifference"].mean()) if len(scoped) else np.nan,
                    "MaxAbsDifference": float(scoped["AbsDifference"].max()) if len(scoped) else np.nan,
                    "RunsMaxDifferenceLE1e1": int((run_max <= 1e-1).sum()),
                    "RunsMaxDifferenceLE1e2": int((run_max <= 1e-2).sum()),
                    "RunsMaxDifferenceLE1e3": int((run_max <= 1e-3).sum()),
                    "ThresholdDenominator": int(len(run_max)),
                }
            )
    run_pairs = runs.loc[runs["Mode"].eq(Q30)].merge(
        runs.loc[runs["Mode"].eq(Q61)],
        on="RunId",
        suffixes=("Q30", "Q61"),
        validate="one_to_one",
    )
    for metric in ["PopulationSD", "LogLik", "EAPReliability"]:
        run_pairs[f"{metric}DifferenceQ61MinusQ30"] = (
            run_pairs[f"{metric}Q61"] - run_pairs[f"{metric}Q30"]
        )
        run_pairs[f"Abs{metric}Difference"] = run_pairs[
            f"{metric}DifferenceQ61MinusQ30"
        ].abs()
    run_summary = []
    for design, group in run_pairs.groupby("DesignQ61", sort=False):
        eligible = group["AnalysisEligibleQ30"].astype(bool) & group["AnalysisEligibleQ61"].astype(bool)
        for scope, mask in [
            ("all_jointly_returned", pd.Series(True, index=group.index)),
            ("jointly_analysis_eligible", eligible),
        ]:
            scoped = group.loc[mask]
            run_summary.append(
                {
                    "ParameterType": "Run-level distribution/loglik",
                    "Design": design,
                    "Scope": scope,
                    "RunIds": int(len(scoped)),
                    "Pairs": int(len(scoped)),
                    "MeanAbsDifference": float(scoped["AbsPopulationSDDifference"].mean()) if len(scoped) else np.nan,
                    "MaxAbsDifference": float(scoped["AbsPopulationSDDifference"].max()) if len(scoped) else np.nan,
                    "RunsMaxDifferenceLE1e1": int((scoped["AbsPopulationSDDifference"] <= 1e-1).sum()),
                    "RunsMaxDifferenceLE1e2": int((scoped["AbsPopulationSDDifference"] <= 1e-2).sum()),
                    "RunsMaxDifferenceLE1e3": int((scoped["AbsPopulationSDDifference"] <= 1e-3).sum()),
                    "ThresholdDenominator": int(len(scoped)),
                    "MaxAbsLogLikDifference": float(scoped["AbsLogLikDifference"].max()) if len(scoped) else np.nan,
                    "MaxAbsEAPReliabilityDifference": float(scoped["AbsEAPReliabilityDifference"].max()) if len(scoped) else np.nan,
                }
            )
    summary = pd.concat([pd.DataFrame(rows), pd.DataFrame(run_summary)], ignore_index=True)
    return pairs, summary


def build_python_sirt_rater_agreement(
    python_recovery: pd.DataFrame,
    sirt_recovery: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare rater recovery as cross-estimator sensitivity, never parity."""
    python = python_recovery.loc[python_recovery["Facet"].eq("Rater")].copy()
    sirt = sirt_recovery.loc[sirt_recovery["Mode"].eq(Q61)].copy()
    keys = ["RunId", "Level", "ComparisonScale"]
    left = python[keys + [
        "ConditionId", "Design", "TruthBias", "Replicate", "Seed", "TruthAligned",
        "EstimateAligned", "SE", "ErrorAligned", "IncludedInSummary"
    ]].rename(
        columns={
            "EstimateAligned": "EstimatePythonJMLE",
            "SE": "SEPythonJMLE",
            "ErrorAligned": "ErrorPythonJMLE",
            "IncludedInSummary": "IncludedPythonJMLE",
        }
    )
    right = sirt[keys + [
        "EstimateAligned", "SE", "ErrorAligned", "IncludedInSummary", "Anchored"
    ]].rename(
        columns={
            "EstimateAligned": "EstimateSirtMML",
            "SE": "SESirtMML",
            "ErrorAligned": "ErrorSirtMML",
            "IncludedInSummary": "IncludedSirtMML",
        }
    )
    pairs = left.merge(right, on=keys, how="outer", validate="one_to_one", indicator=True)
    pairs["BothPresent"] = pairs["_merge"].eq("both")
    pairs["BothIncluded"] = (
        pairs["BothPresent"]
        & pairs["IncludedPythonJMLE"].fillna(False).astype(bool)
        & pairs["IncludedSirtMML"].fillna(False).astype(bool)
    )
    pairs["DifferencePythonJMLEMinusSirtMML"] = (
        pairs["EstimatePythonJMLE"] - pairs["EstimateSirtMML"]
    )
    pairs["AbsDifference"] = pairs["DifferencePythonJMLEMinusSirtMML"].abs()
    pairs["InterpretationBoundary"] = (
        "cross-estimator sensitivity: Python additive RSM JMLE versus sirt virtual-item PCM MML"
    )
    pairs = pairs.drop(columns="_merge")

    rows: list[dict[str, object]] = []
    for design, group in pairs.groupby("Design", sort=False):
        included = group.loc[group["BothIncluded"]]
        rank_rows = []
        for run_id, run_group in included.groupby("RunId", sort=False):
            if len(run_group) < 2:
                continue
            rank_rows.append(
                {
                    "RunId": run_id,
                    "Spearman": run_group["EstimatePythonJMLE"].corr(
                        run_group["EstimateSirtMML"], method="spearman"
                    ),
                    "ExactOrder": run_group.sort_values("EstimatePythonJMLE")["Level"].tolist()
                    == run_group.sort_values("EstimateSirtMML")["Level"].tolist(),
                }
            )
        rank = pd.DataFrame(rank_rows)
        rows.append(
            {
                "Design": design,
                "RunIds": int(group["RunId"].nunique()),
                "BothIncludedRunIds": int(included["RunId"].nunique()),
                "Pairs": int(len(group)),
                "BothIncludedPairs": int(len(included)),
                "MeanAbsDifference": float(included["AbsDifference"].mean()) if len(included) else np.nan,
                "MaxAbsDifference": float(included["AbsDifference"].max()) if len(included) else np.nan,
                "PythonJMLERMSE": float(np.sqrt(np.mean(np.square(included["ErrorPythonJMLE"])))) if len(included) else np.nan,
                "SirtMMLRMSE": float(np.sqrt(np.mean(np.square(included["ErrorSirtMML"])))) if len(included) else np.nan,
                "MeanRunSpearman": float(rank["Spearman"].mean()) if len(rank) else np.nan,
                "ExactOrderRunIds": int(rank["ExactOrder"].sum()) if len(rank) else 0,
                "OrderRunIdDenominator": int(len(rank)),
                "InterpretationBoundary": (
                    "different estimator, item-threshold structure, and Person treatment; not parity"
                ),
            }
        )
    return pairs, pd.DataFrame(rows)


def build_anchor_sensitivity(
    sirt_recovery: pd.DataFrame,
    sirt_runs: pd.DataFrame,
) -> pd.DataFrame:
    clean = sirt_recovery.loc[sirt_recovery["Design"].eq("balanced_large_anchors")]
    drift = sirt_recovery.loc[sirt_recovery["Design"].eq("anchor_drift")]
    keys = ["Mode", "TruthBias", "Replicate", "Seed", "Level"]
    pairs = clean.merge(drift, on=keys, suffixes=("Clean", "Drift"), validate="one_to_one")
    pairs["EstimateShiftDriftMinusClean"] = pairs["EstimateDrift"] - pairs["EstimateClean"]
    pairs["ShiftErrorFromInjected0p25"] = pairs["EstimateShiftDriftMinusClean"] - 0.25
    run_clean = sirt_runs.loc[sirt_runs["Design"].eq("balanced_large_anchors")]
    run_drift = sirt_runs.loc[sirt_runs["Design"].eq("anchor_drift")]
    run_pairs = run_clean.merge(
        run_drift,
        on=["Mode", "TruthBias", "Replicate", "Seed"],
        suffixes=("Clean", "Drift"),
        validate="one_to_one",
    )
    run_pairs["DevianceShiftDriftMinusClean"] = run_pairs["DevianceDrift"] - run_pairs["DevianceClean"]
    run_pairs["PopulationSDShiftDriftMinusClean"] = run_pairs["PopulationSDDrift"] - run_pairs["PopulationSDClean"]
    return pairs.merge(
        run_pairs[[
            "Mode", "TruthBias", "Replicate", "Seed",
            "DevianceShiftDriftMinusClean", "PopulationSDShiftDriftMinusClean"
        ]],
        on=["Mode", "TruthBias", "Replicate", "Seed"],
        validate="many_to_one",
    )


def build_omitted_bias_sensitivity(sirt_recovery: pd.DataFrame) -> pd.DataFrame:
    primary = sirt_recovery.loc[sirt_recovery["Mode"].eq(Q61)]
    null = primary.loc[primary["TruthBias"].eq(0)].rename(
        columns={
            "RunId": "RunIdNull", "EstimateAligned": "EstimateNull",
            "IncludedInSummary": "IncludedNull"
        }
    )
    alternative = primary.loc[primary["TruthBias"].ne(0)].rename(
        columns={
            "RunId": "RunIdAlternative", "TruthBias": "TruthBiasAlternative",
            "EstimateAligned": "EstimateAlternative",
            "IncludedInSummary": "IncludedAlternative"
        }
    )
    keys = ["Design", "Replicate", "Seed", "Level"]
    pairs = null[keys + ["RunIdNull", "EstimateNull", "IncludedNull"]].merge(
        alternative[keys + [
            "RunIdAlternative", "TruthBiasAlternative", "EstimateAlternative", "IncludedAlternative"
        ]],
        on=keys,
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    pairs["PairEligible"] = (
        pairs["_merge"].eq("both")
        & pairs["IncludedNull"].fillna(False).astype(bool)
        & pairs["IncludedAlternative"].fillna(False).astype(bool)
    )
    pairs["EstimateChangeAlternativeMinusNull"] = pairs["EstimateAlternative"] - pairs["EstimateNull"]
    pairs["AbsEstimateChange"] = pairs["EstimateChangeAlternativeMinusNull"].abs()
    pairs["InterpretationBoundary"] = (
        "omitted Rater x Task interaction leakage into sirt rater main effects; not bias detection"
    )
    return pairs.drop(columns="_merge")


def build_sparse_identification(
    sirt_runs: pd.DataFrame,
    python_runs: pd.DataFrame,
    python_cmle_runs: pd.DataFrame,
    mfrmr_runs: pd.DataFrame,
) -> pd.DataFrame:
    sirt = sirt_runs.loc[
        sirt_runs["Mode"].eq(Q61) & sirt_runs["Design"].eq("sparse_missing")
    ].copy()
    python = python_runs.loc[python_runs["Design"].eq("sparse_missing"), [
        "RunId", "FitReturned", "Converged", "InferenceReady"
    ]].rename(columns=lambda name: name if name == "RunId" else f"PythonJMLE{name}")
    cmle = python_cmle_runs.loc[python_cmle_runs["Design"].eq("sparse_missing"), [
        "RunId", "ConditionalDesignEligible", "ConditionalRank", "ConditionalNullity"
    ]].rename(columns=lambda name: name if name == "RunId" else f"PythonCMLE{name}")
    mfrmr = mfrmr_runs.loc[
        mfrmr_runs["Design"].eq("sparse_missing") & mfrmr_runs["Mode"].eq("MFRMR_JML_STRICT"),
        ["RunId", "FitReturned", "Converged", "InferenceReady", "FailureReason"],
    ].rename(columns=lambda name: name if name == "RunId" else f"Mfrmr{name}")
    out = sirt.merge(python, on="RunId", validate="one_to_one").merge(
        cmle, on="RunId", validate="one_to_one"
    ).merge(mfrmr, on="RunId", validate="one_to_one")
    out["IdentificationBoundary"] = (
        "sirt MML can bridge disconnected raters through the estimated common Person distribution; "
        "this is not observed-design connectivity or CMLE/JMLE identification"
    )
    return out


def build_first_read_summary(
    runs: pd.DataFrame,
    quadrature_pairs: pd.DataFrame,
    rater_agreement: pd.DataFrame,
    anchor_sensitivity: pd.DataFrame,
    sparse: pd.DataFrame,
) -> pd.DataFrame:
    """Create a future UI projection without promoting the smoke evidence."""
    primary = runs.loc[runs["Mode"].eq(Q61)]
    eligible_q = quadrature_pairs.loc[quadrature_pairs["BothIncluded"]]
    qmax = eligible_q.groupby("ParameterType")["AbsDifference"].max()
    balanced = rater_agreement.loc[
        rater_agreement["Design"].eq("balanced_small")
    ].iloc[0]
    return pd.DataFrame(
        [
            {
                "Order": 1,
                "Section": "Evidence scope",
                "Status": "Pilot only",
                "Headline": "Two-replicate smoke sensitivity",
                "Evidence": "Failure accounting and estimand boundaries are tested; operating characteristics are not estimated.",
                "UserAction": "Do not use these rates as performance or sample-size evidence.",
            },
            {
                "Order": 2,
                "Section": "Run accounting",
                "Status": "Review",
                "Headline": f"{int(runs['FitReturned'].sum())}/{len(runs)} returned; {int(runs['AnalysisEligible'].sum())}/{len(runs)} eligible",
                "Evidence": "One sparse +0.60 RunId reached the iteration cap in both quadrature modes.",
                "UserAction": "Inspect returned-but-excluded rows before interpreting numerical contrasts.",
            },
            {
                "Order": 3,
                "Section": "Estimator boundary",
                "Status": "Sensitivity only",
                "Headline": "sirt MML is not Python JMLE/CMLE parity",
                "Evidence": f"Balanced-small Rater MAE={balanced['MeanAbsDifference']:.6g}; maximum={balanced['MaxAbsDifference']:.6g} logits.",
                "UserAction": "Interpret closeness as design-specific sensitivity, not estimator identity.",
            },
            {
                "Order": 4,
                "Section": "Quadrature",
                "Status": "Review",
                "Headline": "Raw and identified scales differ",
                "Evidence": (
                    f"Eligible maxima: Rater={qmax['Rater severity']:.6g}; "
                    f"item centered={qmax['Virtual-item centered location']:.6g}; "
                    f"item raw={qmax['Virtual-item raw location']:.6g}; "
                    f"Person EAP={qmax['Person EAP']:.6g}."
                ),
                "UserAction": "Use unrounded values and the declared comparison scale for threshold decisions.",
            },
            {
                "Order": 5,
                "Section": "Sparse identification",
                "Status": "Caution",
                "Headline": f"Primary sirt eligible for {int(sparse['AnalysisEligible'].sum())}/{len(sparse)} sparse RunIds",
                "Evidence": "The common estimated Person distribution can bridge Raters rejected by observed-design rank gates.",
                "UserAction": "Treat this as assumption-based MML identification, not observed connectivity.",
            },
            {
                "Order": 6,
                "Section": "Anchor contamination",
                "Status": "Caution",
                "Headline": "Injected anchor shift is transmitted",
                "Evidence": f"Maximum deviation from the injected +0.25 shift={anchor_sensitivity['ShiftErrorFromInjected0p25'].abs().max():.6g} logits.",
                "UserAction": "Audit anchor provenance; do not interpret the result as anchor robustness.",
            },
            {
                "Order": 7,
                "Section": "Public application surface",
                "Status": "Withheld",
                "Headline": "No sirt or cross-engine result button is enabled",
                "Evidence": f"Primary sirt eligible runs={int(primary['AnalysisEligible'].sum())}/{len(primary)} under a smoke-only profile.",
                "UserAction": "Complete TAM, study-depth, platform, and independent UX review gates first.",
            },
        ]
    ).assign(PublicSurfaceEnabled=False)


def plot_rater_agreement(pairs: pd.DataFrame, path: Path) -> None:
    data = pairs.loc[pairs["BothPresent"]].copy()
    colors = {
        "balanced_small": "#2563eb",
        "balanced_large_anchors": "#059669",
        "anchor_drift": "#d97706",
        "sparse_missing": "#dc2626",
    }
    fig, ax = plt.subplots(figsize=(8.2, 7.0))
    for design, group in data.groupby("Design", sort=False):
        included = group.loc[group["BothIncluded"]]
        excluded = group.loc[~group["BothIncluded"]]
        ax.scatter(
            included["EstimatePythonJMLE"], included["EstimateSirtMML"],
            s=42, alpha=0.72, color=colors.get(design, "#6b7280"), label=design
        )
        if len(excluded):
            ax.scatter(
                excluded["EstimatePythonJMLE"], excluded["EstimateSirtMML"],
                marker="x", s=58, linewidth=1.6, alpha=0.9,
                color=colors.get(design, "#6b7280"),
            )
    limits = [
        float(np.nanmin(data[["EstimatePythonJMLE", "EstimateSirtMML"]].to_numpy())) - 0.15,
        float(np.nanmax(data[["EstimatePythonJMLE", "EstimateSirtMML"]].to_numpy())) + 0.15,
    ]
    ax.plot(limits, limits, "--", color="#111827", linewidth=1.2, label="Identity reference")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_xlabel("Python additive RSM JMLE rater estimate (logits)")
    ax.set_ylabel("sirt virtual-item PCM MML rater estimate (logits)")
    ax.set_title("Rater estimates: cross-estimator sensitivity, not parity")
    ax.grid(alpha=0.2)
    handles, labels = ax.get_legend_handles_labels()
    handles.extend([
        Line2D([0], [0], marker="o", linestyle="none", color="#4b5563", label="jointly analysis eligible"),
        Line2D([0], [0], marker="x", linestyle="none", color="#4b5563", label="returned but excluded"),
    ])
    labels.extend(["jointly analysis eligible", "returned but excluded"])
    ax.legend(handles, labels, frameon=False, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_quadrature_sensitivity(pairs: pd.DataFrame, path: Path) -> None:
    data = pairs.loc[pairs["BothPresent"]].copy()
    order = [
        "Rater severity",
        "Virtual-item centered location",
        "Virtual-item raw location",
        "Person EAP",
    ]
    fig, ax = plt.subplots(figsize=(11.2, 6.3))
    rng = np.random.default_rng(20260809)
    for index, parameter_type in enumerate(order):
        block = data.loc[data["ParameterType"].eq(parameter_type)].dropna(subset=["AbsDifference"])
        included = block.loc[block["BothIncluded"]]
        excluded = block.loc[~block["BothIncluded"]]
        jitter = rng.uniform(-0.14, 0.14, size=len(included))
        ax.scatter(
            np.full(len(included), index) + jitter,
            included["AbsDifference"], s=22, alpha=0.55,
        )
        if len(excluded):
            jitter = rng.uniform(-0.14, 0.14, size=len(excluded))
            ax.scatter(
                np.full(len(excluded), index) + jitter,
                excluded["AbsDifference"], marker="x", s=30,
                color="#dc2626", alpha=0.8,
            )
        if len(included):
            ax.scatter(index, included["AbsDifference"].median(), marker="D", s=70, color="#111827", zorder=4)
    for value, style, label in [(1e-1, ":", "0.1"), (1e-2, "--", "0.01"), (1e-3, ":", "0.001")]:
        ax.axhline(value, linestyle=style, linewidth=1.1, label=f"{label} logit reference")
    ax.set_yscale("log")
    display = ["Rater severity", "Virtual-item\ncentered", "Virtual-item\nraw", "Person EAP"]
    ax.set_xticks(range(len(order)), display)
    ax.set_ylabel("Absolute q61 - q30 estimate difference")
    ax.set_title("sirt quadrature sensitivity differs by output type")
    ax.grid(alpha=0.2, which="both", axis="y")
    handles, labels = ax.get_legend_handles_labels()
    handles.extend([
        Line2D([0], [0], marker="D", linestyle="none", color="#111827", label="eligible median"),
        Line2D([0], [0], marker="x", linestyle="none", color="#dc2626", label="returned but excluded"),
    ])
    labels.extend(["eligible median", "returned but excluded"])
    ax.legend(handles, labels, frameon=False, fontsize=8.5, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_results(
    output_dir: Path,
    runs: pd.DataFrame,
    quadrature_pairs: pd.DataFrame,
    quadrature_summary: pd.DataFrame,
    rater_agreement: pd.DataFrame,
    anchor_sensitivity: pd.DataFrame,
    omitted_bias: pd.DataFrame,
    sparse: pd.DataFrame,
) -> None:
    primary = runs.loc[runs["Mode"].eq(Q61)]
    qmax = quadrature_pairs.loc[quadrature_pairs["BothIncluded"]].groupby(
        "ParameterType"
    )["AbsDifference"].max()
    population_sd_max = float(
        runs.loc[runs["Mode"].eq(Q30), ["RunId", "PopulationSD", "AnalysisEligible"]].merge(
            runs.loc[runs["Mode"].eq(Q61), ["RunId", "PopulationSD", "AnalysisEligible"]],
            on="RunId", suffixes=("Q30", "Q61"), validate="one_to_one"
        ).loc[
            lambda frame: frame["AnalysisEligibleQ30"] & frame["AnalysisEligibleQ61"]
        ].assign(
            Difference=lambda frame: (frame["PopulationSDQ61"] - frame["PopulationSDQ30"]).abs()
        )["Difference"].max()
    )
    balanced = rater_agreement.loc[rater_agreement["Design"].eq("balanced_small")].iloc[0]
    sparse_agreement = rater_agreement.loc[rater_agreement["Design"].eq("sparse_missing")].iloc[0]
    anchor_shift_max_error = float(anchor_sensitivity["ShiftErrorFromInjected0p25"].abs().max())
    eligible_omitted = omitted_bias.loc[
        omitted_bias["PairEligible"] & omitted_bias["Design"].eq("balanced_small")
    ]
    text = f"""# sirt rm.facets operating-characteristics smoke sensitivity

## First read

- `sirt::rm.facets()` is an MML estimator, not a JMLE/CMLE parity target. This adapter maps Task x Criterion to six virtual PCM items, estimates a common normal Person distribution, and retains Rater severity. Python additive-RSM JMLE comparisons are therefore labelled cross-estimator sensitivity.
- All {int(runs['FitReturned'].sum())}/{len(runs)} q30/q61 attempts returned. {int(runs['AnalysisEligible'].sum())}/{len(runs)} were analysis eligible. One sparse +0.60 RunId reached the 1,200-iteration cap in both modes even though its relative deviance change was small; the unreturned internal global-parameter criterion remained unresolved.
- Among jointly analysis-eligible runs, quadrature sensitivity was output-dependent. The maximum absolute q61-q30 difference was `{qmax['Rater severity']:.4g}` logits for Rater severity, `{qmax['Virtual-item centered location']:.4g}` for mean-centered virtual-item location, `{qmax['Virtual-item raw location']:.4g}` for raw virtual-item location, `{qmax['Person EAP']:.4g}` for EAP, and `{population_sd_max:.4g}` for population SD. Returned-but-excluded values are retained under a separate scope; raw and centered item differences must not be conflated.
- In balanced-small data, Python JMLE versus primary sirt MML Rater estimates differed by `{balanced['MeanAbsDifference']:.4g}` logits on average and `{balanced['MaxAbsDifference']:.4g}` at most; all {int(balanced['ExactOrderRunIds'])}/{int(balanced['OrderRunIdDenominator'])} RunIds preserved the exact Rater ordering. This close result does not establish estimator identity.
- Sparse data exposed the estimand difference: primary sirt MML was ready for {int(sparse['AnalysisEligible'].sum())}/{len(sparse)} RunIds, while exact CMLE and strict mfrmr rejected all four structurally and Python JMLE returned all four. Among jointly included Rater rows, Python/sirt mean absolute difference was `{sparse_agreement['MeanAbsDifference']:.4g}` and the maximum was `{sparse_agreement['MaxAbsDifference']:.4g}` logits. sirt's common population distribution can bridge disconnected Raters; this is an assumption-based MML identification route, not observed connectivity.
- The +0.25 contaminated-anchor input shifted sirt Rater estimates by essentially +0.25: the maximum deviation from the injected shift was `{anchor_shift_max_error:.3g}` logits, while deviance and population SD were nearly unchanged. This demonstrates transmission of anchor contamination, not robustness.
- The local Rater x Task interaction is not estimated by this sirt specification. In paired balanced-small runs it moved Rater main effects by `{eligible_omitted['AbsEstimateChange'].mean():.4g}` logits on average and `{eligible_omitted['AbsEstimateChange'].max():.4g}` at most; these are leakage diagnostics, not false-positive or power estimates.
- Information criteria are withheld. In sirt 4.2-133, `rm_facets_ic` subtracts the numeric `b.rater.center` code from the Rater parameter count; unanchored center mode 2 reports `RR-2`, which is not treated here as a comparable free-parameter count.

## Interpretation boundary

This remains a two-replicate smoke study. It does not establish package superiority, stable recovery, coverage, a sample-size rule, or equivalence between MML and JMLE. Fixed-anchor SE values are not interpreted as uncertainty, `delta.item` covariance is unavailable, and the public Streamlit surface remains withheld.

## Retained evidence

- `sirt_runs.csv`, `sirt_rater_recovery.csv`, `sirt_item_recovery.csv`, `sirt_person_estimates.csv`
- `python_sirt_rater_pairs.csv`, `python_sirt_rater_agreement.csv`
- `sirt_quadrature_pairs.csv`, `sirt_quadrature_summary.csv`
- `sirt_anchor_sensitivity.csv`, `sirt_omitted_bias_sensitivity.csv`, `sirt_sparse_identification.csv`
- `sirt_first_read_summary.csv` is the ordered future-UI projection; its public-surface row remains `Withheld`.
- `python_sirt_rater_agreement.png`, `sirt_quadrature_sensitivity.png`
- `sirt_adapter_identity.csv` records bundle, function, source-body, and adapter-script hashes.
"""
    (output_dir / "SIRT_RESULTS.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    input_dir = args.input.resolve()
    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    runs = pd.read_csv(input_dir / "sirt_runs.csv")
    rater = pd.read_csv(input_dir / "sirt_rater_recovery.csv")
    item = pd.read_csv(input_dir / "sirt_item_recovery.csv")
    person = pd.read_csv(input_dir / "sirt_person_estimates.csv")
    identity = pd.read_csv(input_dir / "sirt_adapter_identity.csv").iloc[0]
    bundle_hash = sha256_file(input_dir / "generated_bundle_files.csv")
    if str(identity["BundleInventorySHA256"]) != bundle_hash:
        raise RuntimeError("sirt adapter is stale for the current generated bundle.")
    if sha256_file(Path(__file__).with_name("operating_characteristics_sirt.R")) != str(identity["AdapterScriptSHA256"]):
        raise RuntimeError("sirt adapter script changed after its retained run.")

    python_recovery = pd.read_csv(input_dir / "parameter_recovery.csv")
    python_runs = pd.read_csv(input_dir / "runs.csv")
    python_cmle_runs = pd.read_csv(input_dir / "python_cmle_runs.csv")
    mfrmr_runs = pd.read_csv(input_dir / "mfrmr_runs.csv")
    quadrature_pairs, quadrature_summary = build_quadrature_pairs(rater, item, person, runs)
    rater_pairs, rater_agreement = build_python_sirt_rater_agreement(
        python_recovery, rater
    )
    anchor_sensitivity = build_anchor_sensitivity(rater, runs)
    omitted_bias = build_omitted_bias_sensitivity(rater)
    sparse = build_sparse_identification(
        runs, python_runs, python_cmle_runs, mfrmr_runs
    )
    first_read = build_first_read_summary(
        runs, quadrature_pairs, rater_agreement, anchor_sensitivity, sparse
    )

    quadrature_pairs.to_csv(output_dir / "sirt_quadrature_pairs.csv", index=False)
    quadrature_summary.to_csv(output_dir / "sirt_quadrature_summary.csv", index=False)
    rater_pairs.to_csv(output_dir / "python_sirt_rater_pairs.csv", index=False)
    rater_agreement.to_csv(output_dir / "python_sirt_rater_agreement.csv", index=False)
    anchor_sensitivity.to_csv(output_dir / "sirt_anchor_sensitivity.csv", index=False)
    omitted_bias.to_csv(output_dir / "sirt_omitted_bias_sensitivity.csv", index=False)
    sparse.to_csv(output_dir / "sirt_sparse_identification.csv", index=False)
    first_read.to_csv(output_dir / "sirt_first_read_summary.csv", index=False)
    plot_rater_agreement(rater_pairs, output_dir / "python_sirt_rater_agreement.png")
    plot_quadrature_sensitivity(quadrature_pairs, output_dir / "sirt_quadrature_sensitivity.png")
    write_results(
        output_dir,
        runs,
        quadrature_pairs,
        quadrature_summary,
        rater_agreement,
        anchor_sensitivity,
        omitted_bias,
        sparse,
    )


if __name__ == "__main__":
    main()
