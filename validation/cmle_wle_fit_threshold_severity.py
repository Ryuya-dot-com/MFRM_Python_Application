#!/usr/bin/env python3
"""Decompose the post-result MnSq threshold surface into nominal classes."""

from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.fit_threshold_sensitivity import CANONICAL_FIT_THRESHOLDS  # noqa: E402
from mfrm_app.fit_threshold_severity import (  # noqa: E402
    FIT_CLASSES,
    evaluate_fit_class_transition_surface,
)


AMENDMENT = ROOT / "validation/cmle_wle_fit_threshold_surface_severity_amendment_20260810.json"
PARENT_PLAN = ROOT / "validation/cmle_wle_fit_threshold_surface_plan_20260810.json"
PARENT_RUNNER = ROOT / "validation/cmle_wle_fit_threshold_surface.py"
PARENT_OUTPUT = ROOT / "validation/cmle_wle_fit_threshold_surface_20260810"
BOOTSTRAP_OUTPUT = ROOT / "validation/cmle_wle_bootstrap_fit_extension_corrected_20260810"
DEFAULT_OUTPUT = ROOT / "validation/cmle_wle_fit_threshold_surface_severity_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, lineterminator="\n", float_format="%.17g")


def validate_amendment(path: Path) -> tuple[dict[str, object], dict[str, object]]:
    amendment = json.loads(path.read_text(encoding="utf-8"))
    if amendment.get("schema_version") != "mfrm-cmle-wle-fit-threshold-severity-amendment-v1":
        raise ValueError("Unexpected threshold-severity amendment schema.")
    identities = amendment["immutable_parent_evidence"]
    actual_paths = {
        "threshold_surface_plan_sha256": PARENT_PLAN,
        "threshold_surface_decision_sha256": PARENT_OUTPUT / "threshold_surface_decision.json",
        "threshold_surface_runner_sha256": PARENT_RUNNER,
        "threshold_full_surface_sha256": PARENT_OUTPUT / "threshold_full_surface.csv",
        "threshold_display_precision_sha256": PARENT_OUTPUT / "threshold_display_precision.csv",
        "bootstrap_person_draws_sha256": BOOTSTRAP_OUTPUT / "bootstrap_fit_person_draws.csv",
        "bootstrap_baseline_persons_sha256": BOOTSTRAP_OUTPUT / "bootstrap_fit_baseline_persons.csv",
    }
    mismatches = [
        key
        for key, source in actual_paths.items()
        if not source.exists() or sha256_file(source) != identities[key]
    ]
    parent_decision = json.loads(
        (PARENT_OUTPUT / "threshold_surface_decision.json").read_text(encoding="utf-8")
    )
    if parent_decision.get("sensitivity_core_sha256") != identities["threshold_sensitivity_core_sha256"]:
        mismatches.append("threshold_sensitivity_core_sha256")
    if mismatches:
        raise ValueError(f"Immutable parent evidence mismatch: {sorted(mismatches)}")
    if not amendment.get("parent_evidence_must_not_be_overwritten"):
        raise ValueError("The amendment must protect the parent evidence.")
    return amendment, parent_decision


def threshold_triplets() -> list[tuple[float, float, float]]:
    plan = json.loads(PARENT_PLAN.read_text(encoding="utf-8"))
    grid = plan["threshold_grid"]
    return [
        tuple(float(value) for value in values)
        for values in product(
            grid["overfit_upper"],
            grid["acceptable_upper"],
            grid["noisy_upper"],
        )
    ]


MATRIX_KEYS = [
    "Model",
    "Lane",
    "Statistic",
    "OverfitUpper",
    "AcceptableUpper",
    "NoisyUpper",
]


def matrix_checks(matrices: pd.DataFrame, parent: pd.DataFrame) -> pd.DataFrame:
    totals = (
        matrices.groupby(MATRIX_KEYS, as_index=False)
        .agg(MatrixCellTotal=("CellCount", "sum"), PersonReplicates=("PersonReplicates", "first"))
    )
    off_diagonal = (
        matrices.loc[matrices["ClassChanged"]]
        .groupby(MATRIX_KEYS, as_index=False)["CellCount"]
        .sum()
        .rename(columns={"CellCount": "MatrixOffDiagonal"})
    )
    expected = parent.loc[parent["Statistic"].ne("Any"), MATRIX_KEYS + ["ClassTransitions"]]
    checks = totals.merge(off_diagonal, on=MATRIX_KEYS, validate="one_to_one")
    checks = checks.merge(expected, on=MATRIX_KEYS, validate="one_to_one")
    checks["MatrixTotalMatches"] = checks["MatrixCellTotal"].eq(checks["PersonReplicates"])
    checks["OffDiagonalMatchesParent"] = checks["MatrixOffDiagonal"].eq(checks["ClassTransitions"])
    return checks


def canonical_storage_check(matrices: pd.DataFrame) -> pd.DataFrame:
    canonical = matrices.loc[matrices["CanonicalThresholds"]].copy()
    keys = ["Model", "Lane", "Statistic", "BaselineClass", "ReplicateClass"]
    canonical = canonical[keys + ["CellCount"]].rename(columns={"CellCount": "RecomputedCellCount"})
    stored = pd.read_csv(BOOTSTRAP_OUTPUT / "bootstrap_fit_transition_matrix.csv")
    stored = stored.rename(columns={"PersonReplicates": "StoredCellCount"})
    compared = canonical.merge(stored, on=keys, how="outer", validate="one_to_one")
    compared[["RecomputedCellCount", "StoredCellCount"]] = compared[
        ["RecomputedCellCount", "StoredCellCount"]
    ].fillna(0).astype(int)
    compared["CountDifference"] = compared["RecomputedCellCount"] - compared["StoredCellCount"]
    compared["Passed"] = compared["CountDifference"].eq(0)
    return compared


def noisy_upper_profile(matrices: pd.DataFrame) -> pd.DataFrame:
    lower, acceptable, _ = CANONICAL_FIT_THRESHOLDS
    selected = matrices.loc[
        matrices["OverfitUpper"].eq(lower)
        & matrices["AcceptableUpper"].eq(acceptable)
    ].copy()
    rows: list[dict[str, object]] = []
    for identity, frame in selected.groupby(
        ["Model", "Lane", "Statistic", "NoisyUpper"], sort=False
    ):
        denominator = int(frame["PersonReplicates"].iloc[0])
        for source, column in (("Baseline", "BaselineClass"), ("Replicate", "ReplicateClass")):
            counts = frame.groupby(column)["CellCount"].sum().to_dict()
            for fit_class in FIT_CLASSES:
                count = int(counts.get(fit_class, 0))
                rows.append(
                    {
                        "Model": identity[0],
                        "Lane": identity[1],
                        "Statistic": identity[2],
                        "NoisyUpper": float(identity[3]),
                        "ClassSource": source,
                        "FitClass": fit_class,
                        "PersonReplicates": denominator,
                        "ClassCount": count,
                        "ClassShare": float(count / denominator),
                        "ClassificationInput": "finite_unrounded_mnsq",
                    }
                )
    return pd.DataFrame(rows)


def aggregate_profile(profile: pd.DataFrame) -> pd.DataFrame:
    aggregate = (
        profile.groupby(
            ["Statistic", "NoisyUpper", "ClassSource", "FitClass"],
            as_index=False,
        )[["PersonReplicates", "ClassCount"]]
        .sum()
    )
    aggregate["ClassShare"] = aggregate["ClassCount"] / aggregate["PersonReplicates"]
    return aggregate


def plot_profile(aggregate: pd.DataFrame, output: Path) -> None:
    frame = aggregate.loc[
        aggregate["ClassSource"].eq("Replicate")
        & aggregate["FitClass"].isin(["noisy", "distorting"])
    ]
    colors = {"noisy": "#d95f02", "distorting": "#7b3294"}
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.7), sharey=True)
    for ax, statistic in zip(axes, ("Infit", "Outfit")):
        statistic_frame = frame.loc[frame["Statistic"].eq(statistic)]
        for fit_class in ("noisy", "distorting"):
            line = statistic_frame.loc[statistic_frame["FitClass"].eq(fit_class)].sort_values("NoisyUpper")
            ax.plot(
                line["NoisyUpper"],
                line["ClassShare"],
                marker="o",
                color=colors[fit_class],
                label=fit_class,
            )
        ax.axvline(2.0, color="black", linestyle="--", linewidth=1)
        ax.set_title(statistic)
        ax.set_xlabel("Noisy/distorting boundary")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("Replicate-statistic class share")
    axes[-1].legend(frameon=False)
    fig.suptitle("Flat binary transitions can conceal noisy/distorting relabelling\nRaw MnSq; other boundaries fixed at 0.50 and 1.50")
    fig.tight_layout()
    fig.savefig(output / "threshold_noisy_upper_class_redistribution.png", dpi=180)
    plt.close(fig)


def endpoint_counts(aggregate: pd.DataFrame) -> dict[str, dict[str, dict[str, int]]]:
    replicate = aggregate.loc[
        aggregate["ClassSource"].eq("Replicate")
        & aggregate["FitClass"].isin(["noisy", "distorting"])
        & aggregate["NoisyUpper"].isin([1.9, 2.1])
    ]
    result: dict[str, dict[str, dict[str, int]]] = {}
    for statistic in ("Infit", "Outfit"):
        result[statistic] = {}
        for fit_class in ("noisy", "distorting"):
            rows = replicate.loc[
                replicate["Statistic"].eq(statistic)
                & replicate["FitClass"].eq(fit_class)
            ].set_index("NoisyUpper")
            result[statistic][fit_class] = {
                "at_1_9": int(rows.loc[1.9, "ClassCount"]),
                "at_2_1": int(rows.loc[2.1, "ClassCount"]),
                "change_2_1_minus_1_9": int(
                    rows.loc[2.1, "ClassCount"] - rows.loc[1.9, "ClassCount"]
                ),
            }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--amendment", type=Path, default=AMENDMENT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    amendment, parent_decision = validate_amendment(args.amendment)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_threshold_severity_amendment.json").write_text(
        args.amendment.read_text(encoding="utf-8"), encoding="utf-8"
    )

    draws = pd.read_csv(BOOTSTRAP_OUTPUT / "bootstrap_fit_person_draws.csv")
    parent_surface = pd.read_csv(PARENT_OUTPUT / "threshold_full_surface.csv")
    triplets = threshold_triplets()
    matrices = evaluate_fit_class_transition_surface(draws, threshold_triplets=triplets)
    shuffled = evaluate_fit_class_transition_surface(
        draws.sample(frac=1.0, random_state=20260810).reset_index(drop=True),
        threshold_triplets=triplets,
    )
    sort_columns = MATRIX_KEYS + ["BaselineClass", "ReplicateClass"]
    row_order_invariant = matrices.sort_values(sort_columns).reset_index(drop=True).equals(
        shuffled.sort_values(sort_columns).reset_index(drop=True)
    )
    checks = matrix_checks(matrices, parent_surface)
    canonical_check = canonical_storage_check(matrices)
    profile = noisy_upper_profile(matrices)
    profile_aggregate = aggregate_profile(profile)

    combined = profile_aggregate.loc[
        profile_aggregate["FitClass"].isin(["noisy", "distorting"])
    ].groupby(
        ["Statistic", "NoisyUpper", "ClassSource"], as_index=False
    )["ClassCount"].sum()
    combined_invariant = bool(
        combined.groupby(["Statistic", "ClassSource"])["ClassCount"].nunique().eq(1).all()
    )
    parent_oat = parent_surface.loc[
        parent_surface["OverfitUpper"].eq(0.5)
        & parent_surface["AcceptableUpper"].eq(1.5)
    ]
    parent_oat_aggregate = (
        parent_oat.groupby(["Statistic", "NoisyUpper"], as_index=False)[
            "ClassTransitions"
        ].sum()
    )
    binary_ranges = (
        parent_oat_aggregate.groupby("Statistic")["ClassTransitions"]
        .agg(["min", "max"])
        .reset_index()
    )
    endpoints = endpoint_counts(profile_aggregate)

    contract_passed = bool(
        len(triplets) == 125
        and len(matrices) == 4 * 2 * 125 * 16
        and checks["MatrixTotalMatches"].all()
        and checks["OffDiagonalMatchesParent"].all()
        and canonical_check["Passed"].all()
        and row_order_invariant
        and combined_invariant
    )
    write_csv(matrices, output / "threshold_class_transition_full_surface.csv")
    write_csv(checks, output / "threshold_class_transition_checks.csv")
    write_csv(canonical_check, output / "threshold_canonical_stored_matrix_check.csv")
    write_csv(profile, output / "threshold_noisy_upper_class_profile.csv")
    write_csv(profile_aggregate, output / "threshold_noisy_upper_class_profile_aggregate.csv")
    write_csv(binary_ranges, output / "threshold_noisy_upper_binary_transition_ranges.csv")
    plot_profile(profile_aggregate, output)

    decision = {
        "schema_version": "mfrm-cmle-wle-fit-threshold-severity-result-v1",
        "analysis_executed": True,
        "contract_passed": contract_passed,
        "overall_status": (
            "post_result_class_decomposition_complete_validity_withheld"
            if contract_passed
            else "class_decomposition_contract_failed"
        ),
        "input_person_replicates": int(len(draws)),
        "threshold_triplets": int(len(triplets)),
        "complete_matrix_rows": int(len(matrices)),
        "matrix_total_mismatches": int((~checks["MatrixTotalMatches"]).sum()),
        "parent_off_diagonal_mismatches": int((~checks["OffDiagonalMatchesParent"]).sum()),
        "canonical_stored_cell_mismatches": int((~canonical_check["Passed"]).sum()),
        "row_order_invariant": row_order_invariant,
        "noisy_plus_distorting_count_invariant": combined_invariant,
        "binary_transition_ranges_across_noisy_upper": {
            row.Statistic: {"minimum": int(row.min), "maximum": int(row.max)}
            for row in binary_ranges.itertuples(index=False)
        },
        "replicate_class_endpoint_counts": endpoints,
        "flat_binary_profile_implies_severity_invariance": False,
        "classification_uses_unrounded_values": True,
        "automatic_threshold_selection": False,
        "threshold_optimality_validated": False,
        "known_truth_simulation": False,
        "public_ui_integration_authorized": False,
        "parent_overall_status": parent_decision["overall_status"],
        "registered_amendment_sha256": sha256_file(args.amendment),
        "parent_sensitivity_core_sha256": sha256_file(ROOT / "mfrm_app/fit_threshold_sensitivity.py"),
        "severity_addendum_core_sha256": sha256_file(ROOT / "mfrm_app/fit_threshold_severity.py"),
        "runner_source_sha256": sha256_file(Path(__file__)),
    }
    (output / "threshold_severity_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    endpoint_lines = []
    for statistic in ("Infit", "Outfit"):
        values = endpoints[statistic]
        endpoint_lines.append(
            f"- {statistic}: noisy {values['noisy']['at_1_9']} -> {values['noisy']['at_2_1']} "
            f"({values['noisy']['change_2_1_minus_1_9']:+d}); distorting "
            f"{values['distorting']['at_1_9']} -> {values['distorting']['at_2_1']} "
            f"({values['distorting']['change_2_1_minus_1_9']:+d})."
        )
    binary_lines = [
        f"- {row.Statistic}: {int(row.min)}--{int(row.max)} binary transitions."
        for row in binary_ranges.itertuples(index=False)
    ]
    report = f"""# CMLE-WLE threshold class-decomposition addendum

## Decision

**{decision['overall_status']}.** All 16 nominal baseline-by-replicate cells were emitted for each Model x Lane x Statistic x 125 threshold triplets. Matrix totals, parent off-diagonal transition counts, and the stored canonical matrix were reproduced with zero mismatches. This addendum was registered only after the flat noisy-boundary binary profile was observed, so it is diagnostic post-result evidence, not confirmation.

## What the binary transition indicator concealed

With overfit and acceptable boundaries fixed at 0.50 and 1.50, changing the noisy/distorting boundary from 1.90 to 2.10 reallocated replicate labels as follows:

{chr(10).join(endpoint_lines)}

The noisy-plus-distorting total stayed constant, as required by the fixed 1.50 boundary. The original changed/not-changed indicator had these aggregate ranges across the same noisy_upper values:

{chr(10).join(binary_lines)}

Therefore a flat or nearly flat binary profile is not evidence that the severity labels are stable. It is a projection artifact: noisy and distorting are different nominal outcomes but can both remain off-diagonal relative to the baseline class.

## Boundary of use

All classifications use retained finite unrounded MnSq. No ordinal distance was imposed on the four labels, no threshold was selected, and no false-positive, power, known-truth, or public-UI claim is authorized. The parent threshold surface remains immutable. A prospectively registered known-truth simulation is still required before operational threshold controls can be considered.
"""
    (output / "CMLE_WLE_FIT_THRESHOLD_SEVERITY_ADDENDUM.md").write_text(
        report, encoding="utf-8"
    )
    if not contract_passed:
        raise SystemExit("CMLE-WLE threshold class-decomposition gates failed.")


if __name__ == "__main__":
    main()
