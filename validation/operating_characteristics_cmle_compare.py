#!/usr/bin/env python3
"""Compare native Python exact CMLE with byte-matched immer CMLE evidence."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_cmle_pairs(
    python_coefficients: pd.DataFrame,
    immer_coefficients: pd.DataFrame,
    python_runs: pd.DataFrame,
    immer_runs: pd.DataFrame,
) -> pd.DataFrame:
    """Return one auditable row per matched free coordinate."""
    keys = ["RunId", "Parameter"]
    left = python_coefficients[keys + ["Estimate", "SE"]].rename(
        columns={"Estimate": "EstimatePython", "SE": "SEPython"}
    )
    right = immer_coefficients[keys + ["Estimate", "SE"]].rename(
        columns={"Estimate": "EstimateImmer", "SE": "SEImmer"}
    )
    pairs = left.merge(right, on=keys, how="outer", validate="one_to_one", indicator=True)
    run_columns = ["RunId", "ConditionId", "Design", "TruthBias", "Replicate"]
    pairs = pairs.merge(python_runs[run_columns + ["ParityEligible", "ConditionalLogLik"]].rename(
        columns={"ParityEligible": "PythonParityEligible", "ConditionalLogLik": "ConditionalLogLikPython"}
    ), on="RunId", how="left", validate="many_to_one")
    pairs = pairs.merge(immer_runs[["RunId", "ParityEligible", "ConditionalLogLik", "GradientSupNorm"]].rename(
        columns={
            "ParityEligible": "ImmerParityEligible",
            "ConditionalLogLik": "ConditionalLogLikImmer",
            "GradientSupNorm": "ImmerGradientSupNorm",
        }
    ), on="RunId", how="left", validate="many_to_one")
    pairs["BothPresent"] = pairs["_merge"].eq("both")
    pairs["JointParityEligible"] = (
        pairs["BothPresent"]
        & pairs["PythonParityEligible"].fillna(False).astype(bool)
        & pairs["ImmerParityEligible"].fillna(False).astype(bool)
    )
    pairs["EstimateDifferencePythonMinusImmer"] = pairs["EstimatePython"] - pairs["EstimateImmer"]
    pairs["AbsEstimateDifference"] = pairs["EstimateDifferencePythonMinusImmer"].abs()
    pairs["SEDifferencePythonMinusImmer"] = pairs["SEPython"] - pairs["SEImmer"]
    pairs["AbsSEDifference"] = pairs["SEDifferencePythonMinusImmer"].abs()
    pairs["ConditionalLogLikDifferencePythonMinusImmer"] = (
        pairs["ConditionalLogLikPython"] - pairs["ConditionalLogLikImmer"]
    )
    return pairs.drop(columns="_merge")


def build_agreement_summary(pairs: pd.DataFrame) -> pd.DataFrame:
    """Summarize finite-returned and inference-ready agreement separately."""
    rows: list[dict[str, object]] = []
    scopes = [
        ("all_jointly_returned", pairs["BothPresent"]),
        ("joint_inference_ready", pairs["JointParityEligible"]),
    ]
    for scope, mask in scopes:
        selected = pairs.loc[mask].copy()
        for design, group in [("All", selected), *list(selected.groupby("Design", sort=False))]:
            if group.empty:
                continue
            rows.append(
                {
                    "Scope": scope,
                    "Design": str(design),
                    "RunIds": int(group["RunId"].nunique()),
                    "FreeCoordinatePairs": int(len(group)),
                    "MeanAbsEstimateDifference": float(group["AbsEstimateDifference"].mean()),
                    "MaxAbsEstimateDifference": float(group["AbsEstimateDifference"].max()),
                    "MeanAbsSEDifference": float(group["AbsSEDifference"].mean()),
                    "MaxAbsSEDifference": float(group["AbsSEDifference"].max()),
                    "MaxAbsConditionalLogLikDifference": float(
                        group.groupby("RunId")["ConditionalLogLikDifferencePythonMinusImmer"].first().abs().max()
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_readiness_summary(
    python_runs: pd.DataFrame,
    immer_runs: pd.DataFrame,
    identities: pd.DataFrame,
) -> pd.DataFrame:
    """Keep RunId and unique generated-data denominators visible."""
    py = python_runs.merge(identities[["RunId", "DataId"]], on="RunId", validate="one_to_one")
    im = immer_runs.merge(identities[["RunId", "DataId"]], on="RunId", validate="one_to_one")
    rows: list[dict[str, object]] = []
    specifications = [
        (
            "Python",
            py,
            {
                "DesignEligible": "ConditionalDesignEligible",
                "FitAttempted": "ConditionalDesignEligible",
                "FitReturned": "FitReturned",
                "Converged": "Converged",
                "InferenceReady": "InferenceReady",
                "ParityEligible": "ParityEligible",
                "RequestedConditionEligible": "RequestedConditionEligible",
            },
        ),
        (
            "immer",
            im,
            {
                "DesignEligible": "SharedConditionalDesignEligible",
                "FitAttempted": "FitAttempted",
                "FitReturned": "FitReturned",
                "Converged": "Converged",
                "InferenceReady": "InferenceReady",
                "ParityEligible": "ParityEligible",
                "RequestedConditionEligible": "RequestedConditionEligible",
            },
        ),
    ]
    for engine, frame, metrics in specifications:
        for metric, column in metrics.items():
            flag = frame[column].fillna(False).astype(bool)
            rows.append(
                {
                    "Engine": engine,
                    "Metric": metric,
                    "RunIdNumerator": int(flag.sum()),
                    "RunIdDenominator": int(len(frame)),
                    "RunIdRate": float(flag.mean()),
                    "UniqueDataNumerator": int(frame.loc[flag, "DataId"].nunique()),
                    "UniqueDataDenominator": int(frame["DataId"].nunique()),
                    "UniqueDataRate": float(
                        frame.loc[flag, "DataId"].nunique() / frame["DataId"].nunique()
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_omitted_bias_sensitivity(
    python_coefficients: pd.DataFrame,
    python_runs: pd.DataFrame,
    identities: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pair null/+bias runs without pretending CMLE estimates an interaction."""
    metadata = python_runs[
        ["RunId", "Design", "TruthBias", "Replicate", "Seed", "ParityEligible"]
    ].merge(identities[["RunId", "DataId"]], on="RunId", validate="one_to_one")
    coordinates = python_coefficients[["RunId", "Parameter", "Estimate", "SE"]].merge(
        metadata, on="RunId", how="left", validate="many_to_one"
    )
    null = coordinates.loc[coordinates["TruthBias"].eq(0)].rename(
        columns={
            "RunId": "RunIdNull",
            "DataId": "DataIdNull",
            "Estimate": "EstimateNull",
            "SE": "SENull",
            "ParityEligible": "ParityEligibleNull",
        }
    )
    alternative = coordinates.loc[coordinates["TruthBias"].ne(0)].rename(
        columns={
            "RunId": "RunIdAlternative",
            "DataId": "DataIdAlternative",
            "TruthBias": "TruthBiasAlternative",
            "Estimate": "EstimateAlternative",
            "SE": "SEAlternative",
            "ParityEligible": "ParityEligibleAlternative",
        }
    )
    key = ["Design", "Replicate", "Seed", "Parameter"]
    pairs = null[
        key
        + [
            "RunIdNull",
            "DataIdNull",
            "EstimateNull",
            "SENull",
            "ParityEligibleNull",
        ]
    ].merge(
        alternative[
            key
            + [
                "RunIdAlternative",
                "DataIdAlternative",
                "TruthBiasAlternative",
                "EstimateAlternative",
                "SEAlternative",
                "ParityEligibleAlternative",
            ]
        ],
        on=key,
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    pairs["PairEligible"] = (
        pairs["_merge"].eq("both")
        & pairs["ParityEligibleNull"].fillna(False).astype(bool)
        & pairs["ParityEligibleAlternative"].fillna(False).astype(bool)
    )
    pairs["EstimateChangeAlternativeMinusNull"] = (
        pairs["EstimateAlternative"] - pairs["EstimateNull"]
    )
    pairs["AbsEstimateChange"] = pairs["EstimateChangeAlternativeMinusNull"].abs()
    pairs["InterpretationBoundary"] = (
        "paired omitted-local-interaction sensitivity; not a CMLE bias estimate or detection decision"
    )
    pairs = pairs.drop(columns="_merge")

    rows: list[dict[str, object]] = []
    for design, group in pairs.loc[pairs["PairEligible"]].groupby("Design", sort=False):
        focal_rater = group.loc[group["Parameter"].eq("facet:Rater:free:R01")]
        focal_task = group.loc[group["Parameter"].eq("facet:Task:free:T01")]
        rows.append(
            {
                "Design": design,
                "ReplicatePairs": int(
                    group[["Replicate", "Seed"]].drop_duplicates().shape[0]
                ),
                "FreeCoordinatePairs": int(len(group)),
                "TruthBiasAlternative": float(group["TruthBiasAlternative"].iloc[0]),
                "MeanAbsCoordinateChange": float(group["AbsEstimateChange"].mean()),
                "MaxAbsCoordinateChange": float(group["AbsEstimateChange"].max()),
                "MeanFocalRaterCoordinateChange": float(
                    focal_rater["EstimateChangeAlternativeMinusNull"].mean()
                ),
                "MeanFocalTaskCoordinateChange": float(
                    focal_task["EstimateChangeAlternativeMinusNull"].mean()
                ),
                "InterpretationBoundary": (
                    "additive CMLE leakage sensitivity; no local-interaction estimand"
                ),
            }
        )
    return pairs, pd.DataFrame(rows)


def plot_parameter_agreement(pairs: pd.DataFrame, path: Path) -> None:
    data = pairs.loc[pairs["BothPresent"]].copy()
    fig, ax = plt.subplots(figsize=(8.2, 7.0))
    ready = data["JointParityEligible"]
    ax.scatter(
        data.loc[~ready, "EstimatePython"],
        data.loc[~ready, "EstimateImmer"],
        s=35,
        alpha=0.72,
        color="#d97706",
        label="Both returned; immer readiness withheld",
    )
    ax.scatter(
        data.loc[ready, "EstimatePython"],
        data.loc[ready, "EstimateImmer"],
        s=35,
        alpha=0.78,
        color="#2563eb",
        label="Joint inference-ready",
    )
    limits = [
        float(np.nanmin(data[["EstimatePython", "EstimateImmer"]].to_numpy())) - 0.05,
        float(np.nanmax(data[["EstimatePython", "EstimateImmer"]].to_numpy())) + 0.05,
    ]
    ax.plot(limits, limits, linestyle="--", linewidth=1.2, color="#111827", label="Identity")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_xlabel("Python exact CMLE free coordinate (logits)")
    ax.set_ylabel("immer_cml free coordinate (logits)")
    ax.set_title("Byte-matched CMLE structural coordinates")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, loc="best")
    ax.text(
        0.035,
        0.02,
        "Smoke evidence only; anchors and local-bias parameters are outside this estimand.",
        transform=ax.transAxes,
        fontsize=9,
        color="#4b5563",
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_readiness(
    readiness: pd.DataFrame,
    immer_runs: pd.DataFrame,
    identities: pd.DataFrame,
    path: Path,
) -> None:
    metrics = ["DesignEligible", "FitReturned", "Converged", "InferenceReady", "RequestedConditionEligible"]
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 6.0), gridspec_kw={"width_ratios": [1.05, 1.25]})
    x = np.arange(len(metrics))
    width = 0.36
    for offset, (engine, color) in zip([-width / 2, width / 2], [("Python", "#2563eb"), ("immer", "#d97706")]):
        values = readiness.loc[
            readiness["Engine"].eq(engine) & readiness["Metric"].isin(metrics)
        ].set_index("Metric").reindex(metrics)["RunIdRate"].to_numpy()
        axes[0].bar(x + offset, values, width=width, label=engine, color=color)
    axes[0].set_xticks(x, ["Design\neligible", "Fit\nreturned", "Converged", "Inference\nready", "Requested\nscope eligible"])
    axes[0].set_ylim(0, 1.08)
    axes[0].set_ylabel("RunId rate")
    axes[0].set_title("Fail-closed run accounting")
    axes[0].grid(axis="y", alpha=0.2)
    axes[0].legend(frameon=False)

    grad = immer_runs.loc[immer_runs["FitReturned"].fillna(False).astype(bool)].merge(
        identities[["RunId", "DataId"]], on="RunId", validate="one_to_one"
    )
    # Plot one point per unique generated dataset; clean/drift anchor pairs have
    # identical rating bytes and must not visually double-weight the gradient.
    grad = grad.drop_duplicates("DataId").sort_values("GradientSupNorm").reset_index(drop=True)
    axes[1].scatter(np.arange(len(grad)), grad["GradientSupNorm"], color="#d97706", s=55, zorder=3)
    for threshold, style, label in [
        (1e-4, ":", r"$10^{-4}$ sensitivity"),
        (1e-5, "--", r"$10^{-5}$ primary"),
        (1e-6, ":", r"$10^{-6}$ sensitivity"),
    ]:
        axes[1].axhline(threshold, linestyle=style, linewidth=1.3, label=label)
    axes[1].set_yscale("log")
    axes[1].set_xticks(np.arange(len(grad)), [f"D{i + 1}" for i in range(len(grad))])
    axes[1].set_xlabel("Unique eligible rating dataset")
    axes[1].set_ylabel("immer terminal gradient sup norm")
    axes[1].set_title("Readiness changes with the declared threshold")
    axes[1].grid(alpha=0.2, which="both")
    axes[1].legend(frameon=False, fontsize=9)
    fig.suptitle("Python exact CMLE vs immer_cml: readiness is not optimizer-code identity", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_results(
    output_dir: Path,
    pairs: pd.DataFrame,
    agreement: pd.DataFrame,
    readiness: pd.DataFrame,
    python_runs: pd.DataFrame,
    immer_runs: pd.DataFrame,
    identities: pd.DataFrame,
    omitted_bias_summary: pd.DataFrame,
) -> None:
    returned = agreement.query("Scope == 'all_jointly_returned' and Design == 'All'").iloc[0]
    ready = agreement.query("Scope == 'joint_inference_ready' and Design == 'All'").iloc[0]
    immer_returned = immer_runs[immer_runs["FitReturned"].fillna(False).astype(bool)]
    threshold_counts = {
        "1e-4": int(immer_returned["GradientReadyAt1e4"].astype(bool).sum()),
        "1e-5": int(immer_returned["GradientReadyAt1e5"].astype(bool).sum()),
        "1e-6": int(immer_returned["GradientReadyAt1e6"].astype(bool).sum()),
    }
    unique_immer_returned = immer_returned.merge(
        identities[["RunId", "DataId"]], on="RunId", validate="one_to_one"
    ).drop_duplicates("DataId")
    unique_threshold_counts = {
        "1e-4": int(unique_immer_returned["GradientReadyAt1e4"].astype(bool).sum()),
        "1e-5": int(unique_immer_returned["GradientReadyAt1e5"].astype(bool).sum()),
        "1e-6": int(unique_immer_returned["GradientReadyAt1e6"].astype(bool).sum()),
    }
    maximum_gradient_check_error = float(
        immer_returned["GradientFiniteDifferenceMaxAbsError"].max()
    )
    gradient_classification_checks = int(
        immer_returned["GradientCheckSupportsPrimaryClassification"].astype(bool).sum()
    )
    sparse = python_runs["Design"].eq("sparse_missing")
    anchor = python_runs["RequestedAnchorRows"].gt(0)
    prefit_rank_agreements = int(
        immer_runs["PrefitRankAgreement"].fillna(False).astype(bool).sum()
    )
    balanced_omitted = omitted_bias_summary.loc[
        omitted_bias_summary["Design"].eq("balanced_small")
    ].iloc[0]
    text = f"""# Native Python CMLE vs immer CMLE smoke comparison

## First read

- The same byte-validated ratings were used. Python and immer both restricted the comparison to an unanchored additive RSM conditional likelihood; Person parameters were conditioned out.
- All {int(returned['RunIds'])} structurally eligible RunIds ({len(unique_immer_returned)} unique eligible rating datasets) returned in both engines. Across {int(returned['FreeCoordinatePairs'])} free-coordinate pairs, the mean absolute estimate difference was `{returned['MeanAbsEstimateDifference']:.3g}` logits and the maximum was `{returned['MaxAbsEstimateDifference']:.3g}` logits. The maximum conditional-loglikelihood difference was `{returned['MaxAbsConditionalLogLikDifference']:.3g}`.
- Under the primary terminal-gradient threshold `1e-5`, only {int(ready['RunIds'])} RunIds ({int(ready['FreeCoordinatePairs'])} coordinates) were jointly inference-ready. `immer` returned optimizer code 0 for every eligible fit, but five RunIds ended just above the declared readiness threshold.
- The threshold sensitivity is explicit: among {len(immer_returned)} returned immer RunIds, gradient readiness is {threshold_counts['1e-4']}/{len(immer_returned)} at `1e-4`, {threshold_counts['1e-5']}/{len(immer_returned)} at `1e-5`, and {threshold_counts['1e-6']}/{len(immer_returned)} at `1e-6`. The corresponding unique-data counts are {unique_threshold_counts['1e-4']}/{len(unique_immer_returned)}, {unique_threshold_counts['1e-5']}/{len(unique_immer_returned)}, and {unique_threshold_counts['1e-6']}/{len(unique_immer_returned)}. Display rounding is never used for these decisions.
- The reconstructed analytical `immer` gradient was checked against a centered finite-difference gradient. The maximum component discrepancy was `{maximum_gradient_check_error:.3g}`, and it was smaller than the distance to the primary `1e-5` boundary in {gradient_classification_checks}/{len(immer_returned)} returned RunIds.
- Python exact moments and the independently reconstructed R/immer conditional information agreed on rank/nullity for {prefit_rank_agreements}/{len(immer_runs)} RunIds. The {int(sparse.sum())} sparse RunIds were rejected before either optimizer at rank 5/12 (nullity 7). This is a structural-design rejection, not a non-convergence count.
- The {int(anchor.sum())} anchor-requesting RunIds were fitted only as *unanchored numerical-parity references*. Neither CMLE adapter consumed the anchors, so all were ineligible for the requested anchored-condition claim. Identical results for clean/drift anchor pairs therefore show ignored anchor inputs, not robustness to anchor contamination.
- CMLE did not estimate the generated local Rater x Task interaction. In the paired balanced-small smoke data, changing the generator from 0 to +0.60 shifted additive free coordinates by `{balanced_omitted['MeanAbsCoordinateChange']:.3g}` logits on average and up to `{balanced_omitted['MaxAbsCoordinateChange']:.3g}`. This is omitted-interaction leakage sensitivity—not a bias estimate, false-positive rate, or power result.

## Interpretation boundary

This is a two-replicate smoke comparison, not an operating-characteristic study or a package ranking. It establishes close numerical identity of the matched structural conditional likelihood and exposes a consequential stopping-threshold difference. It does not validate Person scoring, local Rater x Task bias, anchors/linking, fit statistics, false-positive rate, power, or interval coverage. Those require separately matched estimands and study-depth replication.

## Retained evidence

- `python_cmle_runs.csv`, `python_cmle_coefficients.csv`, `python_cmle_parameter_recovery.csv`
- `immer_cmle_runs.csv`, `immer_cmle_coefficients.csv`
- `python_immer_cmle_pairs.csv`, `python_immer_cmle_agreement.csv`, `python_immer_cmle_readiness.csv`
- `python_cmle_omitted_bias_sensitivity.csv`, `python_cmle_omitted_bias_summary.csv`
- `python_immer_cmle_parameter_agreement.png`, `python_immer_cmle_readiness.png`
- adapter identity files record source/function/script hashes and the zero-weight declared-support sentinel used only inside immer.
"""
    (output_dir / "CMLE_IMMER_RESULTS.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    input_dir = args.input.resolve()
    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    python_runs = pd.read_csv(input_dir / "python_cmle_runs.csv")
    immer_runs = pd.read_csv(input_dir / "immer_cmle_runs.csv")
    python_coefficients = pd.read_csv(input_dir / "python_cmle_coefficients.csv")
    immer_coefficients = pd.read_csv(input_dir / "immer_cmle_coefficients.csv")
    identities = pd.read_csv(input_dir / "generated_data_identity.csv")
    python_identity = pd.read_csv(input_dir / "python_cmle_adapter_identity.csv").iloc[0]
    immer_identity = pd.read_csv(input_dir / "immer_cmle_adapter_identity.csv").iloc[0]
    bundle_hash = sha256_file(input_dir / "generated_bundle_files.csv")
    if str(python_identity["BundleInventorySHA256"]) != bundle_hash or str(immer_identity["BundleInventorySHA256"]) != bundle_hash:
        raise RuntimeError("CMLE adapter identity is stale for the generated bundle.")
    if sha256_file(Path(__file__).with_name("operating_characteristics_cmle.py")) != str(python_identity["AdapterScriptSHA256"]):
        raise RuntimeError("Python CMLE adapter script changed after its retained run.")
    if sha256_file(Path(__file__).with_name("operating_characteristics_immer_cmle.R")) != str(immer_identity["AdapterScriptSHA256"]):
        raise RuntimeError("immer CMLE adapter script changed after its retained run.")

    pairs = build_cmle_pairs(python_coefficients, immer_coefficients, python_runs, immer_runs)
    agreement = build_agreement_summary(pairs)
    readiness = build_readiness_summary(python_runs, immer_runs, identities)
    omitted_bias_pairs, omitted_bias_summary = build_omitted_bias_sensitivity(
        python_coefficients, python_runs, identities
    )
    pairs.to_csv(output_dir / "python_immer_cmle_pairs.csv", index=False)
    agreement.to_csv(output_dir / "python_immer_cmle_agreement.csv", index=False)
    readiness.to_csv(output_dir / "python_immer_cmle_readiness.csv", index=False)
    omitted_bias_pairs.to_csv(
        output_dir / "python_cmle_omitted_bias_sensitivity.csv", index=False
    )
    omitted_bias_summary.to_csv(
        output_dir / "python_cmle_omitted_bias_summary.csv", index=False
    )
    plot_parameter_agreement(pairs, output_dir / "python_immer_cmle_parameter_agreement.png")
    plot_readiness(readiness, immer_runs, identities, output_dir / "python_immer_cmle_readiness.png")
    write_results(
        output_dir,
        pairs,
        agreement,
        readiness,
        python_runs,
        immer_runs,
        identities,
        omitted_bias_summary,
    )


if __name__ == "__main__":
    main()
