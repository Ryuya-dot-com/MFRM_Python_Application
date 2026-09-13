#!/usr/bin/env python3
"""Normalize Python and mfrmr operating-characteristics smoke evidence."""

from __future__ import annotations

import argparse
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


DEFAULT_INPUT = REPO / "validation" / "operating_characteristics_20260809"


def _bool(series: pd.Series) -> pd.Series:
    return series.astype("string").str.lower().map({"true": True, "false": False}).astype("boolean")


def build_parameter_agreement(
    python_parameters: pd.DataFrame,
    mfrmr_parameters: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pair aligned facet estimates and summarize numerical agreement by mode."""

    keys = ["RunId", "ConditionId", "Design", "TruthBias", "Facet", "Level", "ComparisonScale"]
    py = python_parameters[keys + ["EstimateAligned", "SE", "IncludedInSummary"]].copy()
    mf = mfrmr_parameters[keys + ["Mode", "EstimateAligned", "SE", "IncludedInSummary"]].copy()
    paired = py.merge(mf, on=keys, how="inner", suffixes=("Python", "Mfrmr"), validate="one_to_many")
    paired["EstimateDifferenceMfrmrMinusPython"] = (
        pd.to_numeric(paired["EstimateAlignedMfrmr"], errors="coerce")
        - pd.to_numeric(paired["EstimateAlignedPython"], errors="coerce")
    )
    paired["SEDifferenceMfrmrMinusPython"] = (
        pd.to_numeric(paired["SEMfrmr"], errors="coerce")
        - pd.to_numeric(paired["SEPython"], errors="coerce")
    )
    paired["FinitePair"] = (
        np.isfinite(pd.to_numeric(paired["EstimateAlignedPython"], errors="coerce"))
        & np.isfinite(pd.to_numeric(paired["EstimateAlignedMfrmr"], errors="coerce"))
    )
    paired["BothIncluded"] = (
        _bool(paired["IncludedInSummaryPython"]).fillna(False)
        & _bool(paired["IncludedInSummaryMfrmr"]).fillna(False)
        & paired["FinitePair"]
    )

    rows: list[dict[str, object]] = []
    groups = ["Mode", "ConditionId", "Design", "TruthBias", "Facet", "ComparisonScale"]
    for group_key, part in paired.groupby(groups, dropna=False, sort=False):
        identity = dict(zip(groups, group_key if isinstance(group_key, tuple) else (group_key,)))
        finite = part.loc[part["FinitePair"]]
        included = part.loc[part["BothIncluded"]]
        difference = pd.to_numeric(finite["EstimateDifferenceMfrmrMinusPython"], errors="coerce")
        included_difference = pd.to_numeric(
            included["EstimateDifferenceMfrmrMinusPython"], errors="coerce"
        )
        correlation = (
            float(finite["EstimateAlignedPython"].corr(finite["EstimateAlignedMfrmr"]))
            if len(finite) > 1 else np.nan
        )
        rows.append({
            "SchemaVersion": oc.SCHEMA_VERSION,
            **identity,
            "RowsPaired": int(len(part)),
            "FinitePairs": int(len(finite)),
            "BothIncludedPairs": int(len(included)),
            "MeanDifferenceMfrmrMinusPython": float(difference.mean()) if len(difference) else np.nan,
            "MAEDifference": float(difference.abs().mean()) if len(difference) else np.nan,
            "RMSEDifference": float(np.sqrt(np.mean(np.square(difference)))) if len(difference) else np.nan,
            "MaxAbsDifference": float(difference.abs().max()) if len(difference) else np.nan,
            "Correlation": correlation,
            "Within0p001": int((difference.abs() <= 0.001).sum()),
            "Within0p01": int((difference.abs() <= 0.01).sum()),
            "Within0p05": int((difference.abs() <= 0.05).sum()),
            "IncludedMAEDifference": (
                float(included_difference.abs().mean()) if len(included_difference) else np.nan
            ),
            "Boundary": (
                "Aligned finite pairs describe numerical agreement. BothIncludedPairs additionally require "
                "each engine's inference/recovery inclusion contract; matched and strict modes are not pooled."
            ),
        })
    return paired, pd.DataFrame(rows)


def build_bias_agreement(
    python_bias: pd.DataFrame,
    mfrmr_bias: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pair focal bias screens without treating unavailable rows as agreement."""

    keys = ["RunId", "ConditionId", "Design", "TruthBias", "TruthPositive", "Replicate", "Seed"]
    columns = [
        "AnalysisEligible", "BiasEstimate", "BiasSE", "p_holm", "AbsBias",
        "SparseCell", "DecisionStrongRaw", "DecisionStrongDisplayed",
    ]
    py = python_bias[keys + columns].copy()
    mf = mfrmr_bias[keys + ["Mode", *columns]].copy()
    paired = py.merge(mf, on=keys, how="inner", suffixes=("Python", "Mfrmr"), validate="one_to_many")
    for column in ("BiasEstimate", "BiasSE", "p_holm", "AbsBias"):
        paired[f"{column}DifferenceMfrmrMinusPython"] = (
            pd.to_numeric(paired[f"{column}Mfrmr"], errors="coerce")
            - pd.to_numeric(paired[f"{column}Python"], errors="coerce")
        )
    paired["FiniteBiasPair"] = (
        np.isfinite(pd.to_numeric(paired["BiasEstimatePython"], errors="coerce"))
        & np.isfinite(pd.to_numeric(paired["BiasEstimateMfrmr"], errors="coerce"))
    )
    paired["BothEligible"] = (
        _bool(paired["AnalysisEligiblePython"]).fillna(False)
        & _bool(paired["AnalysisEligibleMfrmr"]).fillna(False)
    )
    python_decision = _bool(paired["DecisionStrongRawPython"])
    mfrmr_decision = _bool(paired["DecisionStrongRawMfrmr"])
    paired["DecisionComparable"] = paired["BothEligible"] & python_decision.notna() & mfrmr_decision.notna()
    paired["StrongDecisionAgreement"] = pd.NA
    comparable = paired["DecisionComparable"]
    paired.loc[comparable, "StrongDecisionAgreement"] = (
        python_decision.loc[comparable].astype(bool).to_numpy()
        == mfrmr_decision.loc[comparable].astype(bool).to_numpy()
    )

    rows: list[dict[str, object]] = []
    for mode, part in paired.groupby("Mode", dropna=False, sort=False):
        finite = part.loc[part["FiniteBiasPair"]]
        eligible = part.loc[part["BothEligible"]]
        comparable_part = part.loc[part["DecisionComparable"]]
        difference = pd.to_numeric(finite["BiasEstimateDifferenceMfrmrMinusPython"], errors="coerce")
        agreement = _bool(comparable_part["StrongDecisionAgreement"]).fillna(False)
        rows.append({
            "SchemaVersion": oc.SCHEMA_VERSION,
            "Mode": mode,
            "Pairs": int(len(part)),
            "FiniteBiasPairs": int(len(finite)),
            "BothEligiblePairs": int(len(eligible)),
            "DecisionComparablePairs": int(len(comparable_part)),
            "StrongDecisionAgreements": int(agreement.sum()),
            "StrongDecisionDisagreements": int(len(agreement) - agreement.sum()),
            "BiasMAEDifference": float(difference.abs().mean()) if len(difference) else np.nan,
            "BiasMaxAbsDifference": float(difference.abs().max()) if len(difference) else np.nan,
            "SEMAEDifference": float(
                pd.to_numeric(finite["BiasSEDifferenceMfrmrMinusPython"], errors="coerce").abs().mean()
            ) if len(finite) else np.nan,
            "Boundary": (
                "Decision agreement is counted only when both engines mark the focal analysis eligible. "
                "Finite but inference-unready screens remain descriptive pairs."
            ),
        })
    return paired, pd.DataFrame(rows)


def build_readiness_contrast(python_runs: pd.DataFrame, mfrmr_runs: pd.DataFrame) -> pd.DataFrame:
    """Retain asymmetric fit-return and readiness outcomes by numerical mode."""

    py = python_runs[["RunId", "FitReturned", "Converged", "InferenceReady", "AnalysisEligible"]].copy()
    mf = mfrmr_runs[[
        "RunId", "Mode", "FitReturned", "Converged", "InferenceReady", "AnalysisEligible",
        "FailureStage", "FailureReason", "GradientNorm",
    ]].copy()
    paired = py.merge(mf, on="RunId", how="inner", suffixes=("Python", "Mfrmr"), validate="one_to_many")
    for column in (
        "FitReturnedPython", "ConvergedPython", "InferenceReadyPython", "AnalysisEligiblePython",
        "FitReturnedMfrmr", "ConvergedMfrmr", "InferenceReadyMfrmr", "AnalysisEligibleMfrmr",
    ):
        paired[column] = _bool(paired[column]).fillna(False)
    paired["Contrast"] = np.select(
        [
            paired["InferenceReadyPython"] & ~paired["FitReturnedMfrmr"],
            paired["InferenceReadyPython"] & paired["FitReturnedMfrmr"] & ~paired["InferenceReadyMfrmr"],
            paired["InferenceReadyPython"] & paired["InferenceReadyMfrmr"],
        ],
        [
            "Python ready; mfrmr structurally rejected/failed",
            "Python ready; mfrmr returned but readiness withheld",
            "both inference ready",
        ],
        default="other readiness combination",
    )
    return (
        paired.groupby(["Mode", "Contrast"], dropna=False, sort=False)
        .size().reset_index(name="Runs")
    )


def _plot_parameter_agreement(pairs: pd.DataFrame, output: Path) -> None:
    work = pairs.loc[
        pairs["Mode"].eq("MFRMR_JML_STRICT") & pairs["BothIncluded"]
    ].copy()
    if work.empty:
        return
    x = pd.to_numeric(work["EstimateAlignedPython"], errors="coerce")
    y = pd.to_numeric(work["EstimateAlignedMfrmr"], errors="coerce")
    finite = np.isfinite(x) & np.isfinite(y)
    work, x, y = work.loc[finite], x.loc[finite], y.loc[finite]
    palette = {
        "balanced_small": "#4C78A8",
        "balanced_large_anchors": "#59A14F",
        "anchor_drift": "#E45756",
    }
    fig, ax = plt.subplots(figsize=(7.5, 7))
    for design, part in work.groupby("Design", sort=False):
        ax.scatter(
            part["EstimateAlignedPython"],
            part["EstimateAlignedMfrmr"],
            label=str(design).replace("_", " "),
            color=palette.get(str(design), "#777777"),
            alpha=0.78,
            edgecolor="white",
            linewidth=0.4,
        )
    lower = float(min(x.min(), y.min()))
    upper = float(max(x.max(), y.max()))
    pad = max((upper - lower) * 0.05, 0.02)
    ax.plot([lower - pad, upper + pad], [lower - pad, upper + pad], "--", color="#222222", linewidth=1)
    ax.set_xlim(lower - pad, upper + pad)
    ax.set_ylim(lower - pad, upper + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Python aligned facet estimate (logits)")
    ax.set_ylabel("mfrmr strict aligned facet estimate (logits)")
    ax.set_title("Python–mfrmr facet agreement among jointly included estimates")
    ax.legend(frameon=False)
    ax.text(
        0.02,
        0.98,
        "Smoke only; sparse rejected fits are not plotted.",
        transform=ax.transAxes,
        va="top",
        color="#555555",
    )
    fig.tight_layout()
    fig.savefig(output / "python_mfrmr_parameter_agreement.png", dpi=180)
    plt.close(fig)


def _plot_bias_agreement(pairs: pd.DataFrame, output: Path) -> None:
    work = pairs.loc[
        pairs["Mode"].eq("MFRMR_JML_STRICT") & pairs["BothEligible"] & pairs["FiniteBiasPair"]
    ].copy()
    if work.empty:
        return
    x = pd.to_numeric(work["BiasEstimatePython"], errors="coerce")
    y = pd.to_numeric(work["BiasEstimateMfrmr"], errors="coerce")
    colours = ["#E45756" if bool(value) else "#59A14F" for value in work["TruthPositive"]]
    fig, ax = plt.subplots(figsize=(7.5, 7))
    ax.scatter(x, y, color=colours, s=58, alpha=0.85, edgecolor="white", linewidth=0.5)
    lower = float(min(x.min(), y.min()))
    upper = float(max(x.max(), y.max()))
    pad = max((upper - lower) * 0.08, 0.03)
    ax.plot([lower - pad, upper + pad], [lower - pad, upper + pad], "--", color="#222222", linewidth=1)
    ax.set_xlim(lower - pad, upper + pad)
    ax.set_ylim(lower - pad, upper + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Python focal bias estimate (logits)")
    ax.set_ylabel("mfrmr strict focal bias estimate (logits)")
    ax.set_title("Focal conditional-bias agreement among jointly eligible screens")
    ax.text(
        0.02,
        0.98,
        "Green=null; red=+0.60-logit generated local bias. Smoke only.",
        transform=ax.transAxes,
        va="top",
        color="#555555",
    )
    fig.tight_layout()
    fig.savefig(output / "python_mfrmr_bias_agreement.png", dpi=180)
    plt.close(fig)


def _plot_readiness(python_runs: pd.DataFrame, mfrmr_runs: pd.DataFrame, output: Path) -> None:
    frames = []
    python = python_runs.copy()
    python["DisplayMode"] = "Python app\nJML"
    frames.append(python)
    for mode, label in [
        ("MFRMR_JML_MATCHED_CONTROL", "mfrmr\nmatched control"),
        ("MFRMR_JML_STRICT", "mfrmr\nstrict"),
    ]:
        part = mfrmr_runs.loc[mfrmr_runs["Mode"].eq(mode)].copy()
        part["DisplayMode"] = label
        frames.append(part)
    work = pd.concat(frames, ignore_index=True, sort=False)
    labels = ["Python app\nJML", "mfrmr\nmatched control", "mfrmr\nstrict"]
    metrics = ["FitReturned", "Converged", "InferenceReady", "AnalysisEligible"]
    metric_labels = ["Fit returned", "Converged", "Inference ready", "Bias eligible"]
    colours = ["#4C78A8", "#72B7B2", "#F58518", "#E45756"]
    x = np.arange(len(labels))
    width = 0.18
    fig, ax = plt.subplots(figsize=(9, 5.6))
    for index, (metric, metric_label, colour) in enumerate(zip(metrics, metric_labels, colours)):
        rates = []
        for label in labels:
            values = _bool(work.loc[work["DisplayMode"].eq(label), metric]).fillna(False)
            rates.append(float(values.mean()) if len(values) else np.nan)
        ax.bar(x + (index - 1.5) * width, rates, width, label=metric_label, color=colour)
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 1.18)
    ax.set_ylabel("Share of 16 attempted RunIds")
    ax.set_title("Readiness policies diverge before numerical agreement is assessed", pad=14)
    ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 0.99))
    ax.text(
        0.01,
        -0.18,
        "mfrmr rejects four rank-deficient sparse designs; matched-control returned fits fail its gradient review.",
        transform=ax.transAxes,
        color="#555555",
    )
    fig.tight_layout()
    fig.savefig(output / "python_mfrmr_readiness.png", dpi=180)
    plt.close(fig)


def summarize(input_dir: Path) -> str:
    python_runs = pd.read_csv(input_dir / "runs.csv")
    python_parameters = pd.read_csv(input_dir / "parameter_recovery.csv")
    python_bias = pd.read_csv(input_dir / "bias_decisions.csv")
    mfrmr_runs = pd.read_csv(input_dir / "mfrmr_runs.csv")
    mfrmr_parameters = pd.read_csv(input_dir / "mfrmr_parameter_recovery.csv")
    mfrmr_bias = pd.read_csv(input_dir / "mfrmr_bias_decisions.csv")

    python_runs["Mode"] = "PYTHON_APP_JML"
    python_parameters["Mode"] = "PYTHON_APP_JML"
    python_bias["Mode"] = "PYTHON_APP_JML"
    combined_runs = pd.concat([python_runs, mfrmr_runs], ignore_index=True, sort=False)
    combined_parameters = pd.concat([python_parameters, mfrmr_parameters], ignore_index=True, sort=False)
    combined_bias = pd.concat([python_bias, mfrmr_bias], ignore_index=True, sort=False)

    accounting, failure_reasons = oc.summarize_run_accounting(
        combined_runs,
        group_columns=["ConditionId", "Design", "TruthBias", "Engine", "Estimator", "Mode"],
    )
    binary = oc.summarize_binary_operating_characteristics(
        combined_bias,
        decision_columns=[
            "DecisionHolmRaw", "DecisionPracticalRaw", "DecisionStrongRaw",
            "DecisionAnyNonSparseFlag",
        ],
        group_columns=[
            "ConditionId", "Design", "TruthBias", "TruthPositive", "Engine", "Estimator", "Mode",
        ],
    )
    estimation = oc.summarize_estimation_operating_characteristics(
        combined_parameters,
        group_columns=[
            "ConditionId", "Design", "TruthBias", "Engine", "Estimator", "Mode",
            "ParameterType", "ComparisonScale",
        ],
    )
    first_read = oc.build_operating_characteristics_first_read(
        accounting,
        binary,
        estimation,
        profile="cross-engine-smoke",
    )
    parameter_pairs, parameter_agreement = build_parameter_agreement(
        python_parameters,
        mfrmr_parameters,
    )
    bias_pairs, bias_agreement = build_bias_agreement(python_bias, mfrmr_bias)
    readiness = build_readiness_contrast(python_runs, mfrmr_runs)

    outputs = {
        "cross_engine_runs_oc.csv": combined_runs,
        "cross_engine_failure_accounting_oc.csv": accounting,
        "cross_engine_failure_reasons_oc.csv": failure_reasons,
        "cross_engine_bias_decisions_oc.csv": combined_bias,
        "cross_engine_binary_oc.csv": binary,
        "cross_engine_parameter_recovery_oc.csv": combined_parameters,
        "cross_engine_estimation_oc.csv": estimation,
        "cross_engine_first_read_summary.csv": first_read,
        "python_mfrmr_parameter_pairs.csv": parameter_pairs,
        "python_mfrmr_parameter_agreement.csv": parameter_agreement,
        "python_mfrmr_bias_pairs.csv": bias_pairs,
        "python_mfrmr_bias_agreement.csv": bias_agreement,
        "python_mfrmr_readiness_contrast.csv": readiness,
    }
    for filename, frame in outputs.items():
        frame.to_csv(input_dir / filename, index=False)

    _plot_parameter_agreement(parameter_pairs, input_dir)
    _plot_bias_agreement(bias_pairs, input_dir)
    _plot_readiness(python_runs, mfrmr_runs, input_dir)

    strict_parameter = parameter_agreement.loc[
        parameter_agreement["Mode"].eq("MFRMR_JML_STRICT")
    ]
    strict_bias = bias_agreement.loc[bias_agreement["Mode"].eq("MFRMR_JML_STRICT")].iloc[0]
    strict_runs = mfrmr_runs.loc[mfrmr_runs["Mode"].eq("MFRMR_JML_STRICT")]
    matched_runs = mfrmr_runs.loc[mfrmr_runs["Mode"].eq("MFRMR_JML_MATCHED_CONTROL")]
    strict_finite = strict_parameter.loc[strict_parameter["BothIncludedPairs"] > 0]
    parameter_max = float(strict_finite["MaxAbsDifference"].max())
    parameter_mae = float(
        np.average(
            strict_finite["MAEDifference"],
            weights=strict_finite["FinitePairs"],
        )
    )
    strict_ready = int(_bool(strict_runs["InferenceReady"]).fillna(False).sum())
    matched_ready = int(_bool(matched_runs["InferenceReady"]).fillna(False).sum())
    strict_rejected = int((~_bool(strict_runs["FitReturned"]).fillna(False)).sum())
    matched_gradient_min = float(pd.to_numeric(matched_runs["GradientNorm"], errors="coerce").min())
    matched_gradient_max = float(pd.to_numeric(matched_runs["GradientNorm"], errors="coerce").max())
    strict_gradient_min = float(pd.to_numeric(strict_runs["GradientNorm"], errors="coerce").min())
    strict_gradient_max = float(pd.to_numeric(strict_runs["GradientNorm"], errors="coerce").max())
    text = f"""# Python–mfrmr operating-characteristics smoke comparison

## Outcome

- mfrmr attempted {len(mfrmr_runs)} fits: 16 RunIds under two separately retained numerical modes.
- Strict mfrmr mode: {strict_ready}/16 inference-ready; {strict_rejected}/16 structurally rejected before optimization.
- Python-control mfrmr mode: {matched_ready}/16 inference-ready; returned-fit terminal gradients ranged from {matched_gradient_min:.6g} to {matched_gradient_max:.6g}.
- Strict returned-fit terminal gradients ranged from {strict_gradient_min:.6g} to {strict_gradient_max:.6g}.
- Among strict, jointly included facet estimates, weighted mean absolute Python–mfrmr difference was {parameter_mae:.6g} logits; maximum grouped absolute difference was {parameter_max:.6g} logits.
- Strict focal bias pairs: {int(strict_bias['BothEligiblePairs'])} jointly eligible; {int(strict_bias['StrongDecisionAgreements'])} strong-decision agreements and {int(strict_bias['StrongDecisionDisagreements'])} disagreements.
- Strict focal bias mean absolute difference was {float(strict_bias['BiasMAEDifference']):.6g} logits.

## Critical interpretation

This smoke run shows near numerical agreement when both engines accept the same
identified fit and mfrmr reaches its strict terminal-gradient contract. It also
shows that optimizer code 0 is insufficient: the looser matched-control mode
returned 12 fits but withheld inference readiness for every one.

The sparse design is a deliberate negative control. Python returned numerical
fits and later withheld the focal bias decision because the cell was sparse;
mfrmr rejected all four sparse RunIds before optimization because the
estimator-specific constrained design had rank 23 of 30 (nullity 7). This is an
important readiness-policy difference, not evidence that either returned
estimate is automatically correct.

Bias probabilities are conditional plug-in screening quantities in both
implementations. Two replicates per condition cannot estimate false-positive
rate, power, coverage, or package superiority. Matched-control and strict
mfrmr modes must not be pooled.

## Retained evidence

- `mfrmr_runs.csv`, `mfrmr_parameter_recovery.csv`, and `mfrmr_bias_decisions.csv`
- `mfrmr_adapter_identity.csv`
- `python_mfrmr_parameter_pairs.csv` and `python_mfrmr_parameter_agreement.csv`
- `python_mfrmr_bias_pairs.csv` and `python_mfrmr_bias_agreement.csv`
- `python_mfrmr_readiness_contrast.csv`
- `python_mfrmr_parameter_agreement.png`, `python_mfrmr_bias_agreement.png`, and `python_mfrmr_readiness.png`
- cross-engine accounting, binary, estimation, and first-read CSVs
"""
    (input_dir / "MFRMR_RESULTS.md").write_text(text, encoding="utf-8")
    return text


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()
    print(summarize(args.input.resolve()))


if __name__ == "__main__":
    main()
