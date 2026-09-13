#!/usr/bin/env python3
"""Run the registered bias and effective-threshold addendum."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.fit_threshold_pilot_addendum import (  # noqa: E402
    canonical_measure_source_pairs,
    recovery_by_fit_flag,
    replicate_wle_recovery,
    summarize_measure_source_pairs,
    summarize_wle_recovery,
    threshold_dimension_audit,
    wle_person_recovery,
)


DEFAULT_ADDENDUM = ROOT / "validation/cmle_wle_fit_known_truth_pilot_addendum_20260810.json"
PARENT = ROOT / "validation/cmle_wle_fit_known_truth_pilot_20260810"
DEFAULT_OUTPUT = ROOT / "validation/cmle_wle_fit_known_truth_pilot_addendum_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, lineterminator="\n", float_format="%.17g")


def validate_addendum(path: Path) -> dict[str, object]:
    addendum = json.loads(path.read_text(encoding="utf-8"))
    if addendum.get("schema_version") != "mfrm-cmle-wle-fit-known-truth-pilot-addendum-v1":
        raise ValueError("Unexpected known-truth pilot addendum schema.")
    identity = addendum["immutable_parent_evidence"]
    sources = {
        "registered_plan_sha256": ROOT / "validation/cmle_wle_fit_known_truth_pilot_plan_20260810.json",
        "pilot_decision_sha256": PARENT / "known_truth_pilot_decision.json",
        "person_fit_sha256": PARENT / "known_truth_person_fit.csv",
        "responses_sha256": PARENT / "known_truth_responses.csv",
        "canonical_operating_summary_sha256": PARENT / "known_truth_canonical_operating_summary.csv",
        "threshold_operating_summary_sha256": PARENT / "known_truth_threshold_operating_summary.csv",
        "pilot_runner_sha256": ROOT / "validation/cmle_wle_fit_known_truth_pilot.py",
        "known_truth_core_sha256": ROOT / "mfrm_app/fit_threshold_operating_characteristics.py",
    }
    mismatches = [
        key
        for key, source in sources.items()
        if not source.exists() or sha256_file(source) != str(identity[key])
    ]
    if mismatches:
        raise ValueError(f"Known-truth addendum parent identity failed: {mismatches}")
    if addendum["frozen_gates"]["performance_value_success_gate"] is not None:
        raise ValueError("The post-result addendum cannot add a performance gate.")
    return addendum


def plot_recovery(summary: pd.DataFrame, output: Path) -> None:
    frame = summary.loc[
        ~summary["WLEExactExtremePattern"].astype(bool)
        & summary["Metric"].isin(["MeanBias", "RootMeanSquaredError"])
    ].copy()
    frame["Label"] = frame["ConditionId"] + " | " + frame["TruthGroup"]
    labels = sorted(frame["Label"].unique())
    colors = {"MeanBias": "#2166ac", "RootMeanSquaredError": "#b2182b"}
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 7.0), sharey=True)
    for ax, metric, title in zip(
        axes,
        ("MeanBias", "RootMeanSquaredError"),
        ("Mean WLE bias", "Mean replicate RMSE"),
    ):
        metric_frame = frame.loc[frame["Metric"].eq(metric)].set_index("Label").reindex(labels)
        y = np.arange(len(labels))
        ax.errorbar(
            metric_frame["MeanReplicateMetric"],
            y,
            xerr=metric_frame["ReplicateMetricMCSE"].fillna(0.0),
            fmt="o",
            color=colors[metric],
            capsize=3,
        )
        if metric == "MeanBias":
            ax.axvline(0.0, color="black", linestyle="--", linewidth=1)
        ax.set_title(title)
        ax.set_xlabel("Logits")
        ax.grid(axis="x", alpha=0.25)
        ax.set_yticks(y, labels, fontsize=8)
        ax.invert_yaxis()
    fig.suptitle("Known-truth WLE recovery (non-extreme response patterns)\nError bars are Monte Carlo SE, not confidence intervals")
    fig.tight_layout()
    fig.savefig(output / "known_truth_wle_recovery.png", dpi=180)
    plt.close(fig)


def plot_measure_pairs(summary: pd.DataFrame, output: Path) -> None:
    frame = summary.copy()
    frame["Label"] = frame["ConditionId"] + " | " + frame["TruthGroup"]
    frame = frame.sort_values("Label").reset_index(drop=True)
    columns = ["BothUnflagged", "WLEOnly", "GeneratingThetaOnly", "BothFlagged"]
    colors = ["#d9d9d9", "#2166ac", "#fdae61", "#b2182b"]
    labels = ["both unflagged", "WLE only", "generating theta only", "both flagged"]
    fig, ax = plt.subplots(figsize=(11.0, 7.2))
    y = np.arange(len(frame))
    left = np.zeros(len(frame), dtype=float)
    for column, color, label in zip(columns, colors, labels):
        share = frame[column] / frame["Persons"]
        ax.barh(y, share, left=left, color=color, label=label)
        left += share.to_numpy(dtype=float)
    ax.set_yticks(y, frame["Label"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("Person share")
    ax.set_title("Canonical either-upper flag by Person-measure source\nSensitivity comparison; neither source is selected post hoc")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(output / "known_truth_measure_source_flags.png", dpi=180)
    plt.close(fig)


def corrected_threshold_plot(parent_projection: pd.DataFrame, output: Path) -> None:
    frame = parent_projection.sort_values("Label").reset_index(drop=True)
    y = np.arange(len(frame))
    fig, ax = plt.subplots(figsize=(10.8, 7.0))
    ax.hlines(y, frame["Minimum"], frame["Maximum"], color="#777777", linewidth=4)
    ax.scatter(frame["Minimum"], y, color="#2166ac", s=24)
    ax.scatter(frame["Maximum"], y, color="#2166ac", s=24)
    ax.scatter(frame["Canonical"], y, color="#b2182b", s=38, zorder=3)
    ax.set_yticks(y, frame["Label"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Mean replicate-specific either-upper flag rate")
    ax.set_title("Either-upper range: five effective acceptable-boundary inputs\n125 nominal triplets collapse because lower/noisy boundaries do not enter this rule")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "known_truth_threshold_ranges_corrected.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--addendum", type=Path, default=DEFAULT_ADDENDUM)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    addendum = validate_addendum(args.addendum)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_known_truth_pilot_addendum.json").write_text(
        args.addendum.read_text(encoding="utf-8"), encoding="utf-8"
    )

    persons = pd.read_csv(PARENT / "known_truth_person_fit.csv")
    rates = pd.read_csv(PARENT / "known_truth_threshold_replicate_surface.csv")
    parent_projection = pd.read_csv(PARENT / "known_truth_threshold_range_projection.csv")
    recovery = wle_person_recovery(persons)
    replicate_recovery = replicate_wle_recovery(recovery)
    recovery_summary = summarize_wle_recovery(replicate_recovery)
    flag_recovery = recovery_by_fit_flag(recovery)
    pairs = canonical_measure_source_pairs(persons)
    pair_summary = summarize_measure_source_pairs(pairs)
    dimension_audit = threshold_dimension_audit(rates)

    outputs = {
        "known_truth_wle_person_recovery.csv": recovery,
        "known_truth_wle_replicate_recovery.csv": replicate_recovery,
        "known_truth_wle_recovery_summary.csv": recovery_summary,
        "known_truth_wle_recovery_by_fit_flag.csv": flag_recovery,
        "known_truth_measure_source_person_pairs.csv": pairs,
        "known_truth_measure_source_flag_summary.csv": pair_summary,
        "known_truth_threshold_dimension_audit.csv": dimension_audit,
    }
    for filename, frame in outputs.items():
        write_csv(frame, output / filename)
    plot_recovery(recovery_summary, output)
    plot_measure_pairs(pair_summary, output)
    corrected_threshold_plot(parent_projection, output)

    identity_mismatches = int((~pairs["BothPresent"]).sum())
    dimension_mismatches = int(dimension_audit["IrrelevantDimensionMismatchGroups"].sum())
    contract_passed = bool(
        identity_mismatches
        == int(addendum["frozen_gates"]["person_identity_mismatches_allowed"])
        and dimension_mismatches
        == int(
            addendum["frozen_gates"][
                "threshold_irrelevant_dimension_mismatches_allowed"
            ]
        )
        and dimension_audit["Passed"].all()
    )
    nonextreme = recovery_summary.loc[
        ~recovery_summary["WLEExactExtremePattern"].astype(bool)
    ]
    bias = nonextreme.loc[nonextreme["Metric"].eq("MeanBias")]
    rmse = nonextreme.loc[nonextreme["Metric"].eq("RootMeanSquaredError")]
    totals = pair_summary[
        ["BothUnflagged", "WLEOnly", "GeneratingThetaOnly", "BothFlagged", "FlagDisagreements"]
    ].sum()
    decision = {
        "schema_version": "mfrm-cmle-wle-fit-known-truth-pilot-addendum-result-v1",
        "analysis_executed": True,
        "contract_passed": contract_passed,
        "overall_status": (
            "post_result_bias_and_dimension_audit_complete_performance_withheld"
            if contract_passed
            else "post_result_addendum_contract_failed"
        ),
        "wle_person_rows": int(len(recovery)),
        "wle_exact_extreme_person_rows": int(recovery["WLEExactExtremePattern"].sum()),
        "nonextreme_mean_bias_range": {
            "minimum": float(bias["MeanReplicateMetric"].min()),
            "maximum": float(bias["MeanReplicateMetric"].max()),
        },
        "nonextreme_mean_replicate_rmse_range": {
            "minimum": float(rmse["MeanReplicateMetric"].min()),
            "maximum": float(rmse["MeanReplicateMetric"].max()),
        },
        "canonical_measure_source_flag_cells": {
            key: int(totals[key])
            for key in ("BothUnflagged", "WLEOnly", "GeneratingThetaOnly", "BothFlagged")
        },
        "canonical_measure_source_flag_disagreements": int(totals["FlagDisagreements"]),
        "person_identity_mismatches": identity_mismatches,
        "threshold_irrelevant_dimension_mismatches": dimension_mismatches,
        "threshold_effective_input_configurations": {
            row.Rule: int(row.EffectiveInputConfigurations)
            for row in dimension_audit.itertuples(index=False)
        },
        "threshold_observed_distinct_flag_count_vectors": {
            row.Rule: int(row.ObservedDistinctFlagCountVectors)
            for row in dimension_audit.itertuples(index=False)
        },
        "post_result_descriptive": True,
        "confirmatory_recovery_claim": False,
        "causal_fit_flag_error_claim": False,
        "automatic_threshold_or_measure_source_selection": False,
        "anchor_or_connectivity_evaluated": False,
        "public_ui_integration_authorized": False,
        "registered_addendum_sha256": sha256_file(args.addendum),
        "addendum_core_sha256": sha256_file(
            ROOT / "mfrm_app/fit_threshold_pilot_addendum.py"
        ),
        "runner_source_sha256": sha256_file(Path(__file__)),
    }
    (output / "known_truth_pilot_addendum_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    dimension_lines = [
        f"- {row.Rule}: {row.NominalThresholdTriplets} nominal triplets, "
        f"{row.EffectiveInputConfigurations} effective input configurations, "
        f"{row.ObservedDistinctFlagCountVectors} observed count vectors."
        for row in dimension_audit.itertuples(index=False)
    ]
    report = f"""# Known-truth Person-fit pilot: bias and threshold-dimension addendum

## Decision

**{decision['overall_status']}.** The addendum was registered after the parent pilot was inspected. It preserves the parent evidence and adds descriptive WLE recovery, Person-measure-source flag sensitivity, and an audit of duplicated threshold dimensions. No performance value is a success gate.

## WLE recovery

- Non-extreme mean replicate bias ranged from {decision['nonextreme_mean_bias_range']['minimum']:.3f} to {decision['nonextreme_mean_bias_range']['maximum']:.3f} logits across condition/truth strata.
- Non-extreme mean replicate RMSE ranged from {decision['nonextreme_mean_replicate_rmse_range']['minimum']:.3f} to {decision['nonextreme_mean_replicate_rmse_range']['maximum']:.3f} logits.
- {decision['wle_exact_extreme_person_rows']} exact-extreme WLE Person rows are retained separately rather than hidden in the non-extreme plot.
- Error is WLE minus generating theta on the absolute generating scale; no post-hoc recentering was applied.

## Person-measure-source sensitivity

Under the same canonical raw either-upper rule, the paired cells were: both unflagged {int(totals['BothUnflagged'])}, WLE only {int(totals['WLEOnly'])}, generating theta only {int(totals['GeneratingThetaOnly'])}, and both flagged {int(totals['BothFlagged'])}. The {decision['canonical_measure_source_flag_disagreements']} disagreements show that estimating Person location from the same responses changes residual fit conclusions. This is sensitivity evidence; generating theta does not replace fit-sample WLE as the primary operational estimand.

## Effective threshold dimensions

{chr(10).join(dimension_lines)}

All irrelevant-dimension invariance checks passed with {dimension_mismatches} mismatches. In particular, the parent either-upper range has five effective acceptable_upper inputs, not 125 independent pieces of information. The corrected figure makes that projection explicit.

## Boundary of use

Recovery and flag/error associations are post-result descriptive pilot outputs. They do not validate the thresholds, establish causality, supply confirmatory bias/RMSE, or resolve calibration, anchor, connectivity, and cross-engine uncertainty. Public UI integration remains withheld.
"""
    (output / "CMLE_WLE_FIT_KNOWN_TRUTH_PILOT_ADDENDUM.md").write_text(
        report, encoding="utf-8"
    )
    if not contract_passed:
        raise SystemExit("Known-truth pilot addendum contract failed.")


if __name__ == "__main__":
    main()
