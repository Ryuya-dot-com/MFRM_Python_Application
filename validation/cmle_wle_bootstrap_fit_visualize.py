#!/usr/bin/env python3
"""Create a readable boundary-focused view of frozen bootstrap fit evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE = (
    ROOT / "validation/cmle_wle_bootstrap_fit_extension_corrected_20260810"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    args = parser.parse_args()
    evidence = args.evidence.resolve()
    draws_path = evidence / "bootstrap_fit_person_draws.csv"
    decision_path = evidence / "bootstrap_fit_extension_decision.json"
    baseline_path = evidence / "bootstrap_fit_baseline_persons.csv"
    draws = pd.read_csv(draws_path)
    baseline = pd.read_csv(baseline_path).drop_duplicates(["Model", "Person"])
    decision = json.loads(decision_path.read_text(encoding="utf-8"))
    ready = draws.loc[draws["PersonFitReady"].fillna(False).astype(bool)].copy()

    boundary_rows: list[dict[str, object]] = []
    mismatch_rows: list[dict[str, object]] = []
    for statistic in ("Infit", "Outfit"):
        status_column = f"{statistic}BoundaryStatus"
        consistent_column = f"{statistic}DisplayDecisionConsistent"
        counts = ready[status_column].value_counts()
        boundary_rows.append(
            {
                "Statistic": statistic,
                "Stable": int(counts.get("stable", 0)),
                "DisplayRoundingBoundary": int(
                    counts.get("display_rounding_boundary", 0)
                ),
                "NumericalBoundary": int(counts.get("numerical_boundary", 0)),
                "RawDisplayMismatches": int(
                    (~ready[consistent_column].fillna(False).astype(bool)).sum()
                ),
            }
        )
        for (model, lane), frame in ready.groupby(["Model", "Lane"], sort=False):
            mismatch_rows.append(
                {
                    "Model": str(model),
                    "Lane": str(lane),
                    "Statistic": statistic,
                    "RawDisplayMismatches": int(
                        (~frame[consistent_column].fillna(False).astype(bool)).sum()
                    ),
                }
            )
    boundary = pd.DataFrame(boundary_rows)
    mismatch = pd.DataFrame(mismatch_rows)
    detail = mismatch.merge(boundary, on="Statistic", validate="many_to_one")
    detail.to_csv(
        evidence / "bootstrap_fit_boundary_mismatch_detail.csv",
        index=False,
        lineterminator="\n",
    )

    direction_parts: list[pd.DataFrame] = []
    for statistic in ("Infit", "Outfit"):
        part = (
            ready.groupby(
                [
                    "Model",
                    "Lane",
                    f"Baseline{statistic}Class",
                    f"{statistic}Class",
                ],
                as_index=False,
            )
            .size()
            .rename(
                columns={
                    f"Baseline{statistic}Class": "BaselineClass",
                    f"{statistic}Class": "ReplicateClass",
                    "size": "PersonReplicates",
                }
            )
        )
        part.insert(2, "Statistic", statistic)
        part["ClassChanged"] = part["BaselineClass"].ne(part["ReplicateClass"])
        direction_parts.append(part)
    directions = pd.concat(direction_parts, ignore_index=True)
    directions.to_csv(
        evidence / "bootstrap_fit_transition_direction_summary.csv",
        index=False,
        lineterminator="\n",
    )

    baseline_class_parts: list[pd.DataFrame] = []
    for statistic in ("Infit", "Outfit"):
        part = (
            baseline.groupby(["Model", f"{statistic}Class"], as_index=False)
            .size()
            .rename(
                columns={
                    f"{statistic}Class": "BaselineClass",
                    "size": "Persons",
                }
            )
        )
        part.insert(1, "Statistic", statistic)
        baseline_class_parts.append(part)
    baseline_classes = pd.concat(baseline_class_parts, ignore_index=True)
    baseline_classes.to_csv(
        evidence / "bootstrap_fit_baseline_class_summary.csv",
        index=False,
        lineterminator="\n",
    )

    person_summary = pd.read_csv(evidence / "bootstrap_fit_person_summary.csv")
    person_summary["AnyTransitionShare"] = (
        person_summary["AnyFitZoneTransitions"]
        / person_summary["SuccessfulFitReplicates"]
    )
    person_concentration = person_summary.sort_values(
        ["AnyTransitionShare", "RawDisplayMismatches"], ascending=[False, False]
    )
    person_concentration.to_csv(
        evidence / "bootstrap_fit_person_transition_concentration.csv",
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )
    extreme_summary = (
        person_summary.groupby(
            ["Model", "Lane", "BaselineExtremeScorePattern"], as_index=False
        )[["AnyFitZoneTransitions", "SuccessfulFitReplicates"]]
        .sum()
    )
    extreme_summary["AnyTransitionShare"] = (
        extreme_summary["AnyFitZoneTransitions"]
        / extreme_summary["SuccessfulFitReplicates"]
    )
    extreme_summary.to_csv(
        evidence / "bootstrap_fit_baseline_extreme_transition_summary.csv",
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw={"width_ratios": [0.9, 1.35]})
    x = np.arange(len(boundary))
    width = 0.34
    display_bars = axes[0].bar(
        x - width / 2,
        boundary["DisplayRoundingBoundary"],
        width,
        color="#f0ad00",
        label="display-rounding boundary",
    )
    numerical_bars = axes[0].bar(
        x + width / 2,
        boundary["NumericalBoundary"],
        width,
        color="#c9302c",
        label="numerical boundary",
    )
    for bars in (display_bars, numerical_bars):
        for bar in bars:
            axes[0].text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.35,
                str(int(bar.get_height())),
                ha="center",
                va="bottom",
                fontsize=10,
            )
    for index, row in boundary.iterrows():
        axes[0].text(
            index,
            -1.8,
            f"stable: {int(row['Stable']):,}",
            ha="center",
            va="top",
            fontsize=9,
            color="#287a2c",
        )
    axes[0].set_xticks(x, boundary["Statistic"])
    axes[0].set_ylabel("Boundary Person-replicates")
    axes[0].set_ylim(-2.2, max(18.0, float(boundary[["DisplayRoundingBoundary", "NumericalBoundary"]].to_numpy().max()) + 3.5))
    axes[0].set_title("Boundary values enlarged")
    axes[0].legend(frameon=False, fontsize=9, loc="upper left")
    axes[0].grid(axis="y", alpha=0.25)

    lane_short = {
        "fixed_score_conditional_pattern": "fixed-score",
        "joint_plugin_parametric": "joint plug-in",
    }
    model_lane = list(dict.fromkeys(zip(mismatch["Model"], mismatch["Lane"])))
    row_labels = [f"{model}\n{lane_short[lane]}" for model, lane in model_lane]
    heat = np.zeros((len(model_lane), 2), dtype=int)
    for row_index, (model, lane) in enumerate(model_lane):
        for column_index, statistic in enumerate(("Infit", "Outfit")):
            selected = mismatch.loc[
                mismatch["Model"].eq(model)
                & mismatch["Lane"].eq(lane)
                & mismatch["Statistic"].eq(statistic),
                "RawDisplayMismatches",
            ]
            heat[row_index, column_index] = int(selected.iloc[0])
    image = axes[1].imshow(heat, cmap="YlOrRd", vmin=0, vmax=max(3, int(heat.max())))
    for row_index in range(heat.shape[0]):
        for column_index in range(heat.shape[1]):
            axes[1].text(
                column_index,
                row_index,
                str(int(heat[row_index, column_index])),
                ha="center",
                va="center",
                fontsize=12,
                color="black",
            )
    axes[1].set_xticks([0, 1], ["Infit", "Outfit"])
    axes[1].set_yticks(np.arange(len(row_labels)), row_labels)
    axes[1].set_title("Raw/display classification mismatches")
    fig.colorbar(image, ax=axes[1], shrink=0.78, label="mismatch count")
    fig.suptitle(
        "CMLE-WLE floating-point decision audit\n"
        "24 display-boundary statistics; 9 classification mismatches; decisions used raw MnSq",
        fontsize=15,
    )
    fig.tight_layout()
    output_path = evidence / "bootstrap_fit_boundary_mismatch_detail.png"
    fig.savefig(output_path, dpi=180)
    plt.close(fig)

    rate_frame = pd.read_csv(evidence / "bootstrap_fit_transition_rates.csv")
    fixed_rate_min = float(
        rate_frame.loc[
            rate_frame["Lane"].eq("fixed_score_conditional_pattern"),
            "TransitionRate",
        ].min()
    )
    fixed_rate_max = float(
        rate_frame.loc[
            rate_frame["Lane"].eq("fixed_score_conditional_pattern"),
            "TransitionRate",
        ].max()
    )
    joint_rate_min = float(
        rate_frame.loc[
            rate_frame["Lane"].eq("joint_plugin_parametric"),
            "TransitionRate",
        ].min()
    )
    joint_rate_max = float(
        rate_frame.loc[
            rate_frame["Lane"].eq("joint_plugin_parametric"),
            "TransitionRate",
        ].max()
    )
    baseline_extreme = extreme_summary.loc[
        extreme_summary["BaselineExtremeScorePattern"].astype(bool)
    ]
    fixed_extreme = baseline_extreme.loc[
        baseline_extreme["Lane"].eq("fixed_score_conditional_pattern")
    ]
    joint_extreme = baseline_extreme.loc[
        baseline_extreme["Lane"].eq("joint_plugin_parametric")
    ]
    critical_review = f"""# Critical decomposition of CMLE-WLE fit transitions

This addendum is descriptive. It does not recompute the fitted statistics, validate the 0.50/1.50/2.00 thresholds, or turn dependent Person-replicates into independent binomial trials.

## What the aggregate count hides

- Per-statistic transition shares were {fixed_rate_min:.3%}--{fixed_rate_max:.3%} in the fixed-score lane and {joint_rate_min:.3%}--{joint_rate_max:.3%} in the joint plug-in lane. The lanes answer different questions and must not be pooled.
- At baseline, each model had 14 acceptable and two overfit Persons for both Infit and Outfit. The direction table shows that most changes leave the broad acceptable band toward overfit or noisy classifications; the total is not a model-misfit prevalence estimate.
- Baseline exact-extreme Persons had {int(fixed_extreme['AnyFitZoneTransitions'].sum())}/{int(fixed_extreme['SuccessfulFitReplicates'].sum())} any-statistic transitions in the fixed-score lane. Their response pattern is degenerate when the total is fixed. In the joint lane they had {int(joint_extreme['AnyFitZoneTransitions'].sum())}/{int(joint_extreme['SuccessfulFitReplicates'].sum())} transitions ({joint_extreme['AnyFitZoneTransitions'].sum() / joint_extreme['SuccessfulFitReplicates'].sum():.3%}) because totals and extreme status may change.
- The highest Person-level any-transition share was {person_concentration['AnyTransitionShare'].max():.3%}. This concentration and the small five-to-six-observation Person records make a universal threshold-stability claim inappropriate.
- Twenty-four statistics were close enough to a threshold for three-decimal display to hide the boundary, and nine displays implied a different class. Counts and plots use raw retained MnSq, not the displayed value.

## Next evidence gate

Vary observations per Person, missingness/connectivity, categories, local dependence, anchor proportion/contamination, and model misspecification under known truth. Report class-transition behavior and false-positive/power properties by raw distance to threshold. Until then, show MnSq as descriptive sensitivity with the raw value and boundary status; keep ZSTD, p-values, public traffic-light claims, and automatic decisions unavailable.
"""
    (evidence / "CMLE_WLE_BOOTSTRAP_FIT_CRITICAL_REVIEW.md").write_text(
        critical_review, encoding="utf-8"
    )

    manifest = {
        "schema_version": "mfrm-cmle-wle-bootstrap-fit-visualization-v1",
        "descriptive_only": True,
        "statistical_recomputation": False,
        "source_person_draws_sha256": sha256_file(draws_path),
        "source_decision_sha256": sha256_file(decision_path),
        "source_contract_passed": bool(decision["contract_passed"]),
        "boundary_statistics": int(
            boundary["DisplayRoundingBoundary"].sum()
            + boundary["NumericalBoundary"].sum()
        ),
        "raw_display_mismatches": int(boundary["RawDisplayMismatches"].sum()),
        "detail_csv_sha256": sha256_file(
            evidence / "bootstrap_fit_boundary_mismatch_detail.csv"
        ),
        "detail_png_sha256": sha256_file(output_path),
        "transition_direction_summary_sha256": sha256_file(
            evidence / "bootstrap_fit_transition_direction_summary.csv"
        ),
        "baseline_class_summary_sha256": sha256_file(
            evidence / "bootstrap_fit_baseline_class_summary.csv"
        ),
        "person_transition_concentration_sha256": sha256_file(
            evidence / "bootstrap_fit_person_transition_concentration.csv"
        ),
        "baseline_extreme_transition_summary_sha256": sha256_file(
            evidence / "bootstrap_fit_baseline_extreme_transition_summary.csv"
        ),
        "critical_review_sha256": sha256_file(
            evidence / "CMLE_WLE_BOOTSTRAP_FIT_CRITICAL_REVIEW.md"
        ),
        "generator_sha256": sha256_file(Path(__file__)),
    }
    (evidence / "bootstrap_fit_visualization_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
