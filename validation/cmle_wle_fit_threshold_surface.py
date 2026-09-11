#!/usr/bin/env python3
"""Map post-result CMLE-WLE MnSq threshold and display sensitivity."""

from __future__ import annotations

import argparse
from itertools import product
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

from mfrm_app.fit_threshold_sensitivity import (  # noqa: E402
    CANONICAL_FIT_THRESHOLDS,
    classify_fit_mnsq_array,
    evaluate_display_precision_sensitivity,
    evaluate_fit_threshold_surface,
)


DEFAULT_PLAN = ROOT / "validation/cmle_wle_fit_threshold_surface_plan_20260810.json"
SOURCE = ROOT / "validation/cmle_wle_bootstrap_fit_extension_corrected_20260810"
DEFAULT_OUTPUT = ROOT / "validation/cmle_wle_fit_threshold_surface_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_plan(plan_path: Path) -> dict[str, object]:
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != "mfrm-cmle-wle-fit-threshold-surface-plan-v1":
        raise ValueError("Unexpected threshold-surface plan schema.")
    paths = {
        "decision_stability_sha256": ROOT / "mfrm_app/decision_stability.py",
        "cmle_wle_person_fit_core_sha256": ROOT / "mfrm_app/cmle_wle_fit.py",
        "bootstrap_person_draws_sha256": SOURCE / "bootstrap_fit_person_draws.csv",
        "bootstrap_baseline_persons_sha256": SOURCE / "bootstrap_fit_baseline_persons.csv",
        "bootstrap_fit_extension_decision_sha256": SOURCE / "bootstrap_fit_extension_decision.json",
        "bootstrap_fit_critical_review_sha256": SOURCE / "CMLE_WLE_BOOTSTRAP_FIT_CRITICAL_REVIEW.md",
    }
    failed = [
        key
        for key, path in paths.items()
        if not path.exists() or sha256_file(path) != str(plan["input_identity"][key])
    ]
    if failed:
        raise ValueError(f"Threshold-surface input identity failed: {failed}")
    if bool(plan["frozen_gates"]["automatic_threshold_selection"]):
        raise ValueError("Threshold-surface plan cannot select a threshold automatically.")
    return plan


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(
        path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )


def threshold_triplets(plan: dict[str, object]) -> list[tuple[float, float, float]]:
    grid = plan["threshold_grid"]
    return [
        tuple(float(value) for value in values)
        for values in product(
            grid["overfit_upper"],
            grid["acceptable_upper"],
            grid["noisy_upper"],
        )
    ]


def canonical_reproduction(draws: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    rows: list[dict[str, object]] = []
    counts: dict[str, int] = {}
    any_recomputed = np.zeros(len(draws), dtype=bool)
    any_stored = draws["AnyFitZoneTransition"].astype(bool).to_numpy()
    for statistic in ("Infit", "Outfit"):
        current = classify_fit_mnsq_array(
            draws[statistic].to_numpy(dtype=float), CANONICAL_FIT_THRESHOLDS
        )
        baseline = classify_fit_mnsq_array(
            draws[f"Baseline{statistic}"].to_numpy(dtype=float),
            CANONICAL_FIT_THRESHOLDS,
        )
        changed = current != baseline
        stored_changed = draws[f"{statistic}ZoneChanged"].astype(bool).to_numpy()
        current_mismatch = current != draws[f"{statistic}Class"].astype(str).to_numpy()
        baseline_mismatch = (
            baseline
            != draws[f"Baseline{statistic}Class"].astype(str).to_numpy()
        )
        transition_mismatch = changed != stored_changed
        any_recomputed |= changed
        counts[f"{statistic}Transitions"] = int(changed.sum())
        rows.append(
            {
                "Statistic": statistic,
                "Rows": int(len(draws)),
                "StoredCurrentClassMismatches": int(current_mismatch.sum()),
                "StoredBaselineClassMismatches": int(baseline_mismatch.sum()),
                "StoredTransitionIndicatorMismatches": int(
                    transition_mismatch.sum()
                ),
                "RecomputedClassTransitions": int(changed.sum()),
                "Passed": bool(
                    not current_mismatch.any()
                    and not baseline_mismatch.any()
                    and not transition_mismatch.any()
                ),
            }
        )
    any_mismatch = any_recomputed != any_stored
    counts["AnyStatisticTransitions"] = int(any_recomputed.sum())
    rows.append(
        {
            "Statistic": "Any",
            "Rows": int(len(draws)),
            "StoredCurrentClassMismatches": np.nan,
            "StoredBaselineClassMismatches": np.nan,
            "StoredTransitionIndicatorMismatches": int(any_mismatch.sum()),
            "RecomputedClassTransitions": int(any_recomputed.sum()),
            "Passed": bool(not any_mismatch.any()),
        }
    )
    return pd.DataFrame(rows), counts


def aggregate_surface(surface: pd.DataFrame) -> pd.DataFrame:
    group_columns = [
        "Statistic",
        "OverfitUpper",
        "AcceptableUpper",
        "NoisyUpper",
        "CanonicalThresholds",
    ]
    aggregate = (
        surface.groupby(group_columns, as_index=False)[
            ["PersonReplicates", "ClassTransitions"]
        ]
        .sum()
    )
    aggregate["TransitionShare"] = (
        aggregate["ClassTransitions"] / aggregate["PersonReplicates"]
    )
    aggregate["ClassificationInput"] = "finite_unrounded_mnsq"
    aggregate["DependentPersonReplicates"] = True
    return aggregate


def one_at_a_time(surface: pd.DataFrame) -> pd.DataFrame:
    lower0, acceptable0, noisy0 = CANONICAL_FIT_THRESHOLDS
    parts: list[pd.DataFrame] = []
    specifications = [
        (
            "overfit_upper",
            "OverfitUpper",
            surface["AcceptableUpper"].eq(acceptable0)
            & surface["NoisyUpper"].eq(noisy0),
            lower0,
        ),
        (
            "acceptable_upper",
            "AcceptableUpper",
            surface["OverfitUpper"].eq(lower0)
            & surface["NoisyUpper"].eq(noisy0),
            acceptable0,
        ),
        (
            "noisy_upper",
            "NoisyUpper",
            surface["OverfitUpper"].eq(lower0)
            & surface["AcceptableUpper"].eq(acceptable0),
            noisy0,
        ),
    ]
    for name, column, mask, canonical in specifications:
        part = surface.loc[mask].copy()
        part["VaryingThreshold"] = name
        part["ThresholdValue"] = part[column]
        part["CanonicalThresholdValue"] = canonical
        part["OffsetFromCanonical"] = part[column] - canonical
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def surface_range(surface: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for identity, frame in surface.groupby(
        ["Model", "Lane", "Statistic"], sort=False
    ):
        canonical = frame.loc[frame["CanonicalThresholds"]]
        minimum = frame["ClassTransitions"].min()
        maximum = frame["ClassTransitions"].max()
        rows.append(
            {
                "Model": identity[0],
                "Lane": identity[1],
                "Statistic": identity[2],
                "PersonReplicates": int(frame["PersonReplicates"].iloc[0]),
                "CanonicalTransitions": int(canonical["ClassTransitions"].iloc[0]),
                "MinimumTransitions": int(minimum),
                "MaximumTransitions": int(maximum),
                "TransitionCountRange": int(maximum - minimum),
                "CanonicalShare": float(canonical["TransitionShare"].iloc[0]),
                "MinimumShare": float(frame["TransitionShare"].min()),
                "MaximumShare": float(frame["TransitionShare"].max()),
                "ThresholdTriplets": int(len(frame)),
                "ThresholdValidityClaim": False,
            }
        )
    return pd.DataFrame(rows)


def proximity_profile(
    draws: pd.DataFrame,
    baseline: pd.DataFrame,
    bands: list[float],
) -> pd.DataFrame:
    baseline_identity = baseline[
        ["Model", "Lane", "Person", "ExtremeScorePattern"]
    ].rename(columns={"ExtremeScorePattern": "BaselineExtremeScorePattern"})
    work = draws.merge(
        baseline_identity,
        on=["Model", "Lane", "Person"],
        how="left",
        validate="many_to_one",
    )
    baseline_unique = baseline.copy()
    baseline_unique["BaselineExtremeScorePattern"] = baseline_unique[
        "ExtremeScorePattern"
    ].astype(bool)
    threshold_specs = [
        ("overfit_upper", 0.5),
        ("acceptable_upper", 1.5),
        ("noisy_upper", 2.0),
    ]
    rows: list[dict[str, object]] = []
    for source_name, source in (
        ("Replicate", work),
        ("Baseline", baseline_unique),
    ):
        for statistic in ("Infit", "Outfit"):
            for (model, lane, extreme), frame in source.groupby(
                ["Model", "Lane", "BaselineExtremeScorePattern"], sort=False
            ):
                values = frame[statistic].to_numpy(dtype=float)
                for threshold_name, threshold in threshold_specs:
                    distance = np.abs(values - threshold)
                    for band in bands:
                        count = int(np.sum(distance <= band))
                        rows.append(
                            {
                                "Source": source_name,
                                "Model": model,
                                "Lane": lane,
                                "BaselineExtremeScorePattern": bool(extreme),
                                "Statistic": statistic,
                                "Threshold": threshold_name,
                                "ThresholdValue": threshold,
                                "AbsoluteDistanceBand": float(band),
                                "Rows": int(len(frame)),
                                "WithinBand": count,
                                "WithinBandShare": float(count / len(frame)),
                                "UncertaintyInterval": False,
                            }
                        )
    return pd.DataFrame(rows)


def person_concentration(
    draws: pd.DataFrame,
    triplets: list[tuple[float, float, float]],
) -> pd.DataFrame:
    surface = evaluate_fit_threshold_surface(
        draws,
        threshold_triplets=triplets,
        group_columns=("Model", "Lane", "Person"),
    )
    rows: list[dict[str, object]] = []
    for identity, frame in surface.groupby(
        ["Model", "Lane", "Person", "Statistic"], sort=False
    ):
        canonical = frame.loc[frame["CanonicalThresholds"]]
        rows.append(
            {
                "Model": identity[0],
                "Lane": identity[1],
                "Person": identity[2],
                "Statistic": identity[3],
                "PersonReplicates": int(frame["PersonReplicates"].iloc[0]),
                "CanonicalTransitionShare": float(
                    canonical["TransitionShare"].iloc[0]
                ),
                "MinimumTransitionShare": float(frame["TransitionShare"].min()),
                "MaximumTransitionShare": float(frame["TransitionShare"].max()),
                "TransitionShareRange": float(
                    frame["TransitionShare"].max()
                    - frame["TransitionShare"].min()
                ),
                "ThresholdTriplets": int(len(frame)),
            }
        )
    return pd.DataFrame(rows)


def first_read(
    *,
    canonical_passed: bool,
    precision: pd.DataFrame,
    ranges: pd.DataFrame,
) -> pd.DataFrame:
    precision_three = precision.loc[precision["DisplayDecimals"].eq(3)]
    precision_disagreements = int(
        precision_three.loc[
            precision_three["Statistic"].isin(["Infit", "Outfit"]),
            "ReplicateClassMismatches",
        ].sum()
    )
    maximum_range = int(ranges["TransitionCountRange"].max())
    return pd.DataFrame(
        [
            {
                "CardOrder": 1,
                "CardId": "canonical_reproduction",
                "Status": "research_ready" if canonical_passed else "blocked",
                "DisplayValue": "raw canonical classes reproduced" if canonical_passed else "canonical mismatch",
                "Interpretation": "The sensitivity map starts from the retained unrounded decisions.",
                "NextAction": "Stop if stored and recomputed canonical decisions differ.",
            },
            {
                "CardOrder": 2,
                "CardId": "display_precision",
                "Status": "caution" if precision_disagreements else "research_ready",
                "DisplayValue": f"{precision_disagreements} three-decimal class mismatches",
                "Interpretation": "Rounded values are counterfactual display decisions, not stored decisions.",
                "NextAction": "Show raw MnSq and boundary status for affected rows.",
            },
            {
                "CardOrder": 3,
                "CardId": "threshold_dependence",
                "Status": "caution",
                "DisplayValue": f"maximum group/stat transition range {maximum_range}",
                "Interpretation": "Classification sensitivity depends on the chosen threshold triplet.",
                "NextAction": "Inspect one-at-a-time and full-factorial surfaces without selecting a favorable threshold.",
            },
            {
                "CardOrder": 4,
                "CardId": "threshold_validity",
                "Status": "withheld",
                "DisplayValue": "false-positive and power unknown",
                "Interpretation": "A post-result surface cannot validate or optimize thresholds.",
                "NextAction": "Run a prospectively registered known-truth repeated simulation.",
            },
            {
                "CardOrder": 5,
                "CardId": "public_ui",
                "Status": "withheld",
                "DisplayValue": "no adjustable-threshold decision control",
                "Interpretation": "Research sensitivity controls must not silently become operational cutoffs.",
                "NextAction": "Complete validity and comprehension gates before UI integration.",
            },
        ]
    )


def plot_outputs(
    aggregate: pd.DataFrame,
    oat: pd.DataFrame,
    ranges: pd.DataFrame,
    proximity: pd.DataFrame,
    output: Path,
) -> None:
    oat_aggregate = (
        oat.groupby(
            ["VaryingThreshold", "ThresholdValue", "Statistic"], as_index=False
        )[["PersonReplicates", "ClassTransitions"]]
        .sum()
    )
    oat_aggregate["TransitionShare"] = (
        oat_aggregate["ClassTransitions"] / oat_aggregate["PersonReplicates"]
    )
    write_csv(oat_aggregate, output / "threshold_one_at_a_time_aggregate.csv")
    names = ["overfit_upper", "acceptable_upper", "noisy_upper"]
    titles = ["Overfit boundary", "Acceptable/noisy boundary", "Noisy/distorting boundary"]
    canonical = dict(zip(names, CANONICAL_FIT_THRESHOLDS))
    colors = {"Infit": "#1f77b4", "Outfit": "#d95f02", "Any": "#6a3d9a"}
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
    for ax, name, title in zip(axes, names, titles):
        frame = oat_aggregate.loc[oat_aggregate["VaryingThreshold"].eq(name)]
        for statistic in ("Infit", "Outfit", "Any"):
            line = frame.loc[frame["Statistic"].eq(statistic)].sort_values(
                "ThresholdValue"
            )
            ax.plot(
                line["ThresholdValue"],
                line["TransitionShare"],
                marker="o",
                color=colors[statistic],
                label=statistic,
            )
        ax.axvline(canonical[name], color="black", linestyle="--", linewidth=1)
        ax.set_title(title)
        ax.set_xlabel("Threshold value")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("Raw class-transition share")
    axes[-1].legend(frameon=False)
    fig.suptitle("One-at-a-time CMLE-WLE threshold sensitivity\nOther thresholds fixed at 0.50 / 1.50 / 2.00")
    fig.tight_layout()
    fig.savefig(output / "threshold_one_at_a_time_profiles.png", dpi=180)
    plt.close(fig)

    range_plot = ranges.copy()
    range_plot["Label"] = (
        range_plot["Model"].astype(str)
        + " | "
        + range_plot["Lane"].map(
            {
                "fixed_score_conditional_pattern": "fixed-score",
                "joint_plugin_parametric": "joint plug-in",
            }
        )
        + " | "
        + range_plot["Statistic"].astype(str)
    )
    range_plot = range_plot.sort_values(
        ["Statistic", "Model", "Lane"], ascending=[True, False, True]
    ).reset_index(drop=True)
    y = np.arange(len(range_plot))
    fig, ax = plt.subplots(figsize=(10.5, 7.0))
    ax.hlines(
        y,
        range_plot["MinimumShare"],
        range_plot["MaximumShare"],
        color="#888888",
        linewidth=4,
    )
    ax.scatter(
        range_plot["CanonicalShare"], y, color="#b2182b", s=42, zorder=3, label="canonical"
    )
    ax.scatter(
        range_plot["MinimumShare"], y, color="#2166ac", s=25, zorder=3, label="grid minimum/maximum"
    )
    ax.scatter(range_plot["MaximumShare"], y, color="#2166ac", s=25, zorder=3)
    ax.set_yticks(y, range_plot["Label"], fontsize=8)
    ax.set_xlabel("Person-replicate transition share")
    ax.set_title("Full 125-triplet sensitivity range\nDescriptive only; no threshold is selected")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "threshold_full_factorial_ranges.png", dpi=180)
    plt.close(fig)

    prox = proximity.loc[proximity["Source"].eq("Replicate")].copy()
    prox = (
        prox.groupby(
            ["Statistic", "Threshold", "ThresholdValue", "AbsoluteDistanceBand"],
            as_index=False,
        )[["Rows", "WithinBand"]]
        .sum()
    )
    prox["WithinBandShare"] = prox["WithinBand"] / prox["Rows"]
    write_csv(prox, output / "threshold_proximity_aggregate.csv")
    fig, ax = plt.subplots(figsize=(9.5, 5.3))
    line_styles = {"Infit": "-", "Outfit": "--"}
    threshold_colors = {
        "overfit_upper": "#1b9e77",
        "acceptable_upper": "#d95f02",
        "noisy_upper": "#7570b3",
    }
    for (statistic, threshold), frame in prox.groupby(["Statistic", "Threshold"]):
        ax.plot(
            frame["AbsoluteDistanceBand"],
            frame["WithinBandShare"],
            marker="o",
            linestyle=line_styles[statistic],
            color=threshold_colors[threshold],
            label=f"{statistic} | {threshold}",
        )
    ax.set_xscale("log")
    ax.set_xlabel("Absolute raw distance from canonical threshold")
    ax.set_ylabel("Replicate-statistic share within band")
    ax.set_title("Empirical mass near canonical MnSq thresholds")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(output / "threshold_raw_proximity.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    plan = validate_plan(args.plan)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_threshold_surface_plan.json").write_text(
        args.plan.read_text(encoding="utf-8"), encoding="utf-8"
    )

    draws = pd.read_csv(SOURCE / "bootstrap_fit_person_draws.csv")
    baseline = pd.read_csv(SOURCE / "bootstrap_fit_baseline_persons.csv")
    if len(draws) != int(plan["frozen_gates"]["input_rows_must_equal"]):
        raise ValueError("Threshold-surface input row count differs from the plan.")
    ready = draws["PersonFitReady"].fillna(False).astype(bool)
    if not ready.all() or not np.isfinite(
        draws[["Infit", "Outfit", "BaselineInfit", "BaselineOutfit"]].to_numpy(
            dtype=float
        )
    ).all():
        raise ValueError("Threshold-surface inputs contain unavailable fit statistics.")

    triplets = threshold_triplets(plan)
    surface = evaluate_fit_threshold_surface(
        draws, threshold_triplets=triplets
    )
    shuffled = evaluate_fit_threshold_surface(
        draws.sample(frac=1.0, random_state=20260810).reset_index(drop=True),
        threshold_triplets=triplets,
    )
    sort_columns = [
        "Model",
        "Lane",
        "Statistic",
        "OverfitUpper",
        "AcceptableUpper",
        "NoisyUpper",
    ]
    row_order_invariant = surface.sort_values(sort_columns).reset_index(
        drop=True
    ).equals(shuffled.sort_values(sort_columns).reset_index(drop=True))
    aggregate = aggregate_surface(surface)
    oat = one_at_a_time(surface)
    ranges = surface_range(surface)
    precision = evaluate_display_precision_sensitivity(
        draws,
        decimals=plan["display_precision_grid"]["decimals"],
    )
    precision_aggregate = (
        precision.groupby(["DisplayDecimals", "Statistic"], as_index=False)[
            [
                "PersonReplicates",
                "RawClassTransitions",
                "RoundedClassTransitions",
                "RoundedMinusRawTransitions",
                "TransitionIndicatorDisagreements",
                "ReplicateClassMismatches",
                "BaselineClassMismatches",
            ]
        ]
        .sum(min_count=1)
    )
    proximity = proximity_profile(
        draws,
        baseline,
        [float(value) for value in plan["proximity_bands"]["absolute_raw_distance"]],
    )
    concentration = person_concentration(draws, triplets)
    reproduction, reproduced_counts = canonical_reproduction(draws)
    expected_counts = plan["canonical_contract"]["expected_input_identity_counts"]
    canonical_count_mismatches = sum(
        int(reproduced_counts[key] != int(expected_counts[key]))
        for key in expected_counts
    )
    canonical_class_mismatches = int(
        reproduction[
            [
                "StoredCurrentClassMismatches",
                "StoredBaselineClassMismatches",
                "StoredTransitionIndicatorMismatches",
            ]
        ]
        .fillna(0)
        .to_numpy()
        .sum()
    )
    unique_triplets = int(
        surface[["OverfitUpper", "AcceptableUpper", "NoisyUpper"]]
        .drop_duplicates()
        .shape[0]
    )
    canonical_passed = bool(
        reproduction["Passed"].all()
        and canonical_count_mismatches == 0
        and canonical_class_mismatches == 0
    )
    contract_passed = bool(
        canonical_passed
        and unique_triplets
        == int(plan["frozen_gates"]["threshold_triplets_must_equal"])
        and row_order_invariant
    )
    first_read_frame = first_read(
        canonical_passed=canonical_passed,
        precision=precision_aggregate,
        ranges=ranges,
    )

    outputs = {
        "threshold_canonical_reproduction.csv": reproduction,
        "threshold_full_surface.csv": surface,
        "threshold_full_surface_aggregate.csv": aggregate,
        "threshold_one_at_a_time.csv": oat,
        "threshold_full_factorial_range.csv": ranges,
        "threshold_display_precision.csv": precision,
        "threshold_display_precision_aggregate.csv": precision_aggregate,
        "threshold_raw_proximity.csv": proximity,
        "threshold_person_concentration.csv": concentration,
        "threshold_first_read_projection.csv": first_read_frame,
    }
    for filename, frame in outputs.items():
        write_csv(frame, output / filename)
    plot_outputs(aggregate, oat, ranges, proximity, output)

    aggregate_ranges = (
        aggregate.groupby("Statistic")["ClassTransitions"]
        .agg(["min", "max"])
        .reset_index()
    )
    precision_three = precision_aggregate.loc[
        precision_aggregate["DisplayDecimals"].eq(3)
    ]
    decision = {
        "schema_version": "mfrm-cmle-wle-fit-threshold-surface-result-v1",
        "analysis_executed": True,
        "contract_passed": contract_passed,
        "overall_status": "post_result_descriptive_surface_complete_validity_withheld"
        if contract_passed
        else "threshold_surface_contract_failed",
        "input_person_replicates": int(len(draws)),
        "threshold_triplets": unique_triplets,
        "canonical_class_or_transition_mismatches": canonical_class_mismatches,
        "canonical_count_mismatches": canonical_count_mismatches,
        "row_order_invariant": row_order_invariant,
        "canonical_transition_counts": reproduced_counts,
        "full_factorial_aggregate_transition_ranges": {
            row.Statistic: {"minimum": int(row.min), "maximum": int(row.max)}
            for row in aggregate_ranges.itertuples(index=False)
        },
        "maximum_group_statistic_transition_count_range": int(
            ranges["TransitionCountRange"].max()
        ),
        "maximum_person_transition_share_range": float(
            concentration["TransitionShareRange"].max()
        ),
        "three_decimal_replicate_class_mismatches": int(
            precision_three.loc[
                precision_three["Statistic"].isin(["Infit", "Outfit"]),
                "ReplicateClassMismatches",
            ].sum()
        ),
        "three_decimal_any_transition_indicator_disagreements": int(
            precision_three.loc[
                precision_three["Statistic"].eq("Any"),
                "TransitionIndicatorDisagreements",
            ].sum()
        ),
        "classification_uses_unrounded_values": True,
        "automatic_threshold_selection": False,
        "threshold_optimality_validated": False,
        "false_positive_or_power_evaluated": False,
        "known_truth_simulation": False,
        "independent_binomial_interpretation": False,
        "ZSTD_or_p_value_authorized": False,
        "streamlit_integration_authorized": False,
        "registered_plan_sha256": sha256_file(args.plan),
        "sensitivity_core_sha256": sha256_file(
            ROOT / "mfrm_app/fit_threshold_sensitivity.py"
        ),
        "runner_source_sha256": sha256_file(Path(__file__)),
    }
    (output / "threshold_surface_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    range_text = ", ".join(
        f"{row.Statistic} {int(row.min)}--{int(row.max)}"
        for row in aggregate_ranges.itertuples(index=False)
    )
    report = f"""# CMLE-WLE MnSq threshold and display sensitivity

## Decision

**{decision['overall_status']}.** The retained raw canonical classifications and transition indicators were reproduced with {canonical_class_mismatches} mismatches, and all {unique_triplets} frozen threshold triplets were evaluated. This is a post-result descriptive surface. It does not validate, optimize, or recommend a threshold.

## Reproduction and display precision

- Canonical raw transitions: Infit {reproduced_counts['InfitTransitions']}, Outfit {reproduced_counts['OutfitTransitions']}, either statistic {reproduced_counts['AnyStatisticTransitions']}.
- Three-decimal replicate class mismatches: {decision['three_decimal_replicate_class_mismatches']} of {2 * len(draws)} statistic rows.
- Three-decimal any-transition indicator disagreements: {decision['three_decimal_any_transition_indicator_disagreements']} of {len(draws)} Person-replicates.
- All stored decisions remain based on finite unrounded MnSq; rounded classifications are counterfactual sensitivity output only.

## Threshold dependence

- Across the full factorial grid, aggregate transition-count ranges were: {range_text}.
- The maximum within-group/statistic transition-count range was {decision['maximum_group_statistic_transition_count_range']} Person-replicates.
- The largest Person-level transition-share range across the grid was {decision['maximum_person_transition_share_range']:.3f}.
- Person-replicates are repeated and dependent; these descriptive shares do not have a simple independent-binomial interpretation.

## Boundary of use

The surface answers how conclusions move when a user changes a rule after fit. It cannot show which rule is correct. No threshold was selected from the observed minimum or maximum. False-positive rate, power, known-truth recovery, local-dependence robustness, anchor contamination, ZSTD/p-values, and public UI controls remain unavailable. The next gate is a prospectively registered known-truth repeated simulation varying observations per Person, missingness/connectivity, categories, local dependence, model misspecification, and anchors.
"""
    (output / "CMLE_WLE_FIT_THRESHOLD_SURFACE_RESULTS.md").write_text(
        report, encoding="utf-8"
    )
    if not contract_passed:
        raise SystemExit("CMLE-WLE threshold-surface gates failed.")


if __name__ == "__main__":
    main()
