#!/usr/bin/env python3
"""Run the prospectively registered known-truth Person MnSq pilot."""

from __future__ import annotations

import argparse
from itertools import product
import hashlib
import json
from pathlib import Path
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.fit_threshold_operating_characteristics import (  # noqa: E402
    PERSON_MEASURE_TRUE,
    PERSON_MEASURE_WLE,
    replicate_fit_rule_rates,
    rounded_rule_disagreements,
    score_known_truth_person_fit,
    simulate_known_truth_person_fit_sample,
    summarize_replicate_fit_rates,
)
from mfrm_app.fit_threshold_sensitivity import (  # noqa: E402
    CANONICAL_FIT_THRESHOLDS,
    classify_fit_mnsq_array,
)
from mfrm_app.operating_characteristics import (  # noqa: E402
    build_replicate_manifest,
    frame_fingerprint,
)


DEFAULT_PLAN = ROOT / "validation/cmle_wle_fit_known_truth_pilot_plan_20260810.json"
DEFAULT_OUTPUT = ROOT / "validation/cmle_wle_fit_known_truth_pilot_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, lineterminator="\n", float_format="%.17g")


def validate_plan(path: Path) -> dict[str, object]:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != "mfrm-cmle-wle-fit-known-truth-pilot-plan-v1":
        raise ValueError("Unexpected known-truth pilot plan schema.")
    identity = plan["input_identity"]
    sources = {
        "fixed_calibration_person_scoring_sha256": ROOT / "mfrm_app/person_scoring.py",
        "cmle_wle_person_fit_core_sha256": ROOT / "mfrm_app/cmle_wle_fit.py",
        "decision_stability_sha256": ROOT / "mfrm_app/decision_stability.py",
        "threshold_sensitivity_core_sha256": ROOT / "mfrm_app/fit_threshold_sensitivity.py",
        "operating_characteristics_core_sha256": ROOT / "mfrm_app/operating_characteristics.py",
        "cmle_core_sha256": ROOT / "mfrm_app/cmle.py",
        "operating_characteristics_protocol_sha256": ROOT / "validation/OPERATING_CHARACTERISTICS_PROTOCOL_20260809.md",
        "threshold_surface_decision_sha256": ROOT / "validation/cmle_wle_fit_threshold_surface_20260810/threshold_surface_decision.json",
        "threshold_class_decomposition_decision_sha256": ROOT / "validation/cmle_wle_fit_threshold_surface_severity_20260810/threshold_severity_decision.json",
        "python_sirt_person_fit_parity_decision_sha256": ROOT / "validation/cmle_wle_person_fit_20260810/person_fit_parity_decision.json",
    }
    mismatches = [
        key
        for key, source in sources.items()
        if not source.exists() or sha256_file(source) != str(identity[key])
    ]
    if mismatches:
        raise ValueError(f"Known-truth pilot input identity failed: {mismatches}")
    gates = plan["frozen_pilot_gates"]
    if gates["performance_value_success_gate"] is not None:
        raise ValueError("A debugging pilot must not have a performance success gate.")
    if bool(gates["automatic_threshold_selection"]):
        raise ValueError("A debugging pilot cannot select a threshold automatically.")
    return plan


def threshold_triplets(plan: dict[str, object]) -> list[tuple[float, float, float]]:
    surface = plan["threshold_surface"]
    return [
        tuple(float(value) for value in values)
        for values in product(
            surface["overfit_upper"],
            surface["acceptable_upper"],
            surface["noisy_upper"],
        )
    ]


def run_one(
    row: pd.Series,
    plan: dict[str, object],
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame]:
    constants = plan["simulation_constants"]
    thresholds = constants["category_thresholds"][str(int(row["Categories"]))]
    deviations = (
        plan["threshold_heterogeneity_contract"][
            "criterion_step_deviations_for_four_categories"
        ]
        if str(row["Mechanism"])
        == "affected_person_pcm_threshold_heterogeneity"
        else None
    )
    started = time.perf_counter()
    sample = simulate_known_truth_person_fit_sample(
        row.to_dict(),
        seed=int(row["Seed"]),
        persons=int(constants["persons_per_replicate"]),
        rater_effects=constants["rater_effects"],
        criterion_effects=constants["criterion_effects"],
        category_thresholds=thresholds,
        threshold_deviations=deviations,
    )
    persons = score_known_truth_person_fit(sample)
    responses = sample["responses"].copy()
    identity = {
        "ConditionId": str(row["ConditionId"]),
        "Replicate": int(row["Replicate"]),
        "RunId": str(row["RunId"]),
        "Seed": int(row["Seed"]),
        "SeedGroup": str(row["SeedGroup"]),
        "Categories": int(row["Categories"]),
        "ObservationsPerPerson": int(row["ObservationsPerPerson"]),
    }
    for column, value in reversed(list(identity.items())):
        responses.insert(0, column, value)
        persons.insert(0, column, value)
    audit = sample["audit"]
    ledger = {
        **identity,
        "Attempted": True,
        "Returned": True,
        "FailureStage": "none",
        "FailureReason": "",
        "ResponseRows": int(len(responses)),
        "Persons": int(constants["persons_per_replicate"]),
        "AffectedPersons": int(audit["AffectedPersons"]),
        "MechanismAppliedRows": int(audit["MechanismAppliedRows"]),
        "ModifiedRows": int(audit["ModifiedRows"]),
        "CategorySupportViolations": int(audit["CategorySupportViolations"]),
        "WLEPersonRows": int(persons["PersonMeasureSource"].eq(PERSON_MEASURE_WLE).sum()),
        "TrueThetaPersonRows": int(persons["PersonMeasureSource"].eq(PERSON_MEASURE_TRUE).sum()),
        "PersonFitReadyRows": int(persons["PersonFitReady"].sum()),
        "ElapsedSeconds": float(time.perf_counter() - started),
    }
    return ledger, responses, persons


def failure_ledger(row: pd.Series, exc: Exception, elapsed: float) -> dict[str, object]:
    return {
        "ConditionId": str(row["ConditionId"]),
        "Replicate": int(row["Replicate"]),
        "RunId": str(row["RunId"]),
        "Seed": int(row["Seed"]),
        "SeedGroup": str(row["SeedGroup"]),
        "Categories": int(row["Categories"]),
        "ObservationsPerPerson": int(row["ObservationsPerPerson"]),
        "Attempted": True,
        "Returned": False,
        "FailureStage": "generation_or_person_scoring",
        "FailureReason": f"{type(exc).__name__}: {str(exc)[:300]}",
        "ResponseRows": 0,
        "Persons": 0,
        "AffectedPersons": 0,
        "MechanismAppliedRows": 0,
        "ModifiedRows": 0,
        "CategorySupportViolations": 0,
        "WLEPersonRows": 0,
        "TrueThetaPersonRows": 0,
        "PersonFitReadyRows": 0,
        "ElapsedSeconds": float(elapsed),
    }


def replay_checks(
    manifest: pd.DataFrame,
    plan: dict[str, object],
    responses: pd.DataFrame,
    persons: pd.DataFrame,
) -> pd.DataFrame:
    selected_ids = [
        "rsm4_random_dense::rep-00001",
        "rsm4_random_sparse::rep-00001",
    ]
    rows: list[dict[str, object]] = []
    for run_id in selected_ids:
        manifest_row = manifest.loc[manifest["RunId"].eq(run_id)].iloc[0]
        _, replay_responses, replay_persons = run_one(manifest_row, plan)
        original_responses = responses.loc[responses["RunId"].eq(run_id)].reset_index(drop=True)
        original_persons = persons.loc[persons["RunId"].eq(run_id)].reset_index(drop=True)
        original_response_hash = frame_fingerprint(original_responses)
        replay_response_hash = frame_fingerprint(replay_responses)
        original_person_hash = frame_fingerprint(original_persons)
        replay_person_hash = frame_fingerprint(replay_persons)
        rows.append(
            {
                "RunId": run_id,
                "OriginalResponseFingerprint": original_response_hash,
                "ReplayResponseFingerprint": replay_response_hash,
                "OriginalPersonFingerprint": original_person_hash,
                "ReplayPersonFingerprint": replay_person_hash,
                "ResponseReplayExact": original_response_hash == replay_response_hash,
                "PersonReplayExact": original_person_hash == replay_person_hash,
                "Passed": original_response_hash == replay_response_hash
                and original_person_hash == replay_person_hash,
            }
        )
    return pd.DataFrame(rows)


def canonical_recomputation(persons: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for statistic in ("Infit", "Outfit"):
        recomputed = classify_fit_mnsq_array(
            persons[statistic].to_numpy(dtype=float), CANONICAL_FIT_THRESHOLDS
        )
        stored = persons[f"{statistic}Class"].astype(str).to_numpy()
        rows.append(
            {
                "Statistic": statistic,
                "PersonRows": int(len(persons)),
                "ClassMismatches": int(np.sum(recomputed != stored)),
                "Passed": bool(np.array_equal(recomputed, stored)),
            }
        )
    return pd.DataFrame(rows)


def truth_identity_mismatches(persons: pd.DataFrame) -> int:
    keys = ["ConditionId", "Replicate", "Person"]
    checks = persons.groupby(keys).agg(
        Sources=("PersonMeasureSource", "nunique"),
        TruthGroups=("TruthGroup", "nunique"),
        AffectedValues=("Affected", "nunique"),
        ThetaValues=("TrueTheta", "nunique"),
    )
    return int(
        (
            checks["Sources"].ne(2)
            | checks["TruthGroups"].ne(1)
            | checks["AffectedValues"].ne(1)
            | checks["ThetaValues"].ne(1)
        ).sum()
    )


def availability_summary(persons: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for identity, frame in persons.groupby(
        ["ConditionId", "PersonMeasureSource", "TruthGroup"], sort=False
    ):
        ready = frame["PersonFitReady"].fillna(False).astype(bool)
        extreme = np.fromiter(
            (
                False if pd.isna(value) else bool(value)
                for value in frame["WLEExtremeScorePattern"]
            ),
            dtype=bool,
            count=len(frame),
        )
        rows.append(
            {
                "ConditionId": identity[0],
                "PersonMeasureSource": identity[1],
                "TruthGroup": identity[2],
                "PersonRows": int(len(frame)),
                "PersonFitReady": int(ready.sum()),
                "PersonFitUnavailable": int((~ready).sum()),
                "WLEExactExtremePatterns": int(extreme.sum())
                if identity[1] == PERSON_MEASURE_WLE
                else np.nan,
                "MeanObservations": float(frame["NObservations"].mean()),
                "MedianInfit": float(frame.loc[ready, "Infit"].median())
                if ready.any()
                else np.nan,
                "MedianOutfit": float(frame.loc[ready, "Outfit"].median())
                if ready.any()
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


def first_read(
    *,
    contract_passed: bool,
    ledger: pd.DataFrame,
    canonical_summary: pd.DataFrame,
    rounding: pd.DataFrame,
) -> pd.DataFrame:
    wle = canonical_summary.loc[
        canonical_summary["PersonMeasureSource"].eq(PERSON_MEASURE_WLE)
        & canonical_summary["Rule"].eq("either_upper")
    ]
    max_clean = float(
        wle.loc[wle["TruthGroup"].eq("clean"), "MeanReplicateRate"].max()
    )
    affected = wle.loc[wle["TruthGroup"].eq("affected"), "MeanReplicateRate"]
    affected_range = (
        f"{affected.min():.3f}--{affected.max():.3f}"
        if len(affected)
        else "unavailable"
    )
    three = rounding.loc[
        rounding["DisplayDecimals"].eq(3)
        & rounding["Rule"].eq("either_upper")
        & rounding["PersonMeasureSource"].eq(PERSON_MEASURE_WLE)
    ]
    disagreements = int(three["RuleDisagreements"].sum())
    return pd.DataFrame(
        [
            {
                "CardOrder": 1,
                "CardId": "pilot_contract",
                "Status": "research_ready" if contract_passed else "blocked",
                "DisplayValue": f"{int(ledger['Returned'].sum())}/{len(ledger)} runs returned",
                "Interpretation": "This checks the frozen pilot machinery, not diagnostic validity.",
                "NextAction": "Stop and inspect the ledger if any contract gate fails.",
            },
            {
                "CardOrder": 2,
                "CardId": "clean_person_flags",
                "Status": "caution",
                "DisplayValue": f"maximum pilot clean-person rate {max_clean:.3f}",
                "Interpretation": "Rates are means of replicate-specific rates under known truth.",
                "NextAction": "Inspect observation-count and category strata; do not select a threshold.",
            },
            {
                "CardOrder": 3,
                "CardId": "affected_person_detection",
                "Status": "caution",
                "DisplayValue": f"pilot affected-person range {affected_range}",
                "Interpretation": "Mechanism-specific detection is not confirmatory power.",
                "NextAction": "Retain low and high detection mechanisms without post-hoc relabelling.",
            },
            {
                "CardOrder": 4,
                "CardId": "rounding",
                "Status": "caution" if disagreements else "research_ready",
                "DisplayValue": f"{disagreements} three-decimal either-upper disagreements",
                "Interpretation": "Raw decisions remain authoritative.",
                "NextAction": "Expose raw values and boundary status in any future research UI.",
            },
            {
                "CardOrder": 5,
                "CardId": "anchors_and_connectivity",
                "Status": "withheld",
                "DisplayValue": "not estimated in fixed-calibration Phase A",
                "Interpretation": "Python CMLE lacks the required native hard-anchor contract.",
                "NextAction": "Validate anchors and conditional connectivity before a refit extension.",
            },
            {
                "CardOrder": 6,
                "CardId": "public_ui",
                "Status": "withheld",
                "DisplayValue": "no traffic-light or threshold control",
                "Interpretation": "Twenty replicates per condition are design debugging only.",
                "NextAction": "Register precision and comprehension gates before public integration.",
            },
        ]
    )


def plot_outputs(
    canonical: pd.DataFrame,
    threshold_summary: pd.DataFrame,
    rounding: pd.DataFrame,
    output: Path,
) -> None:
    primary = canonical.loc[
        canonical["PersonMeasureSource"].eq(PERSON_MEASURE_WLE)
        & canonical["Rule"].eq("either_upper")
    ].copy()
    primary["Label"] = primary["ConditionId"] + " | " + primary["TruthGroup"]
    primary = primary.sort_values(["ConditionId", "TruthGroup"]).reset_index(drop=True)
    colors = primary["TruthGroup"].map({"clean": "#2166ac", "affected": "#b2182b"})
    fig, ax = plt.subplots(figsize=(10.8, 6.8))
    y = np.arange(len(primary))
    ax.barh(
        y,
        primary["MeanReplicateRate"],
        xerr=primary["ReplicateRateMCSE"].fillna(0.0),
        color=colors,
        alpha=0.88,
        capsize=3,
    )
    ax.set_yticks(y, primary["Label"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Mean replicate-specific either-upper flag rate")
    ax.set_title("Known-truth MnSq pilot: clean and affected Persons\nError bars are Monte Carlo SE, not confidence intervals")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "known_truth_canonical_rates.png", dpi=180)
    plt.close(fig)

    rounded = rounding.loc[
        rounding["PersonMeasureSource"].eq(PERSON_MEASURE_WLE)
    ].groupby(["DisplayDecimals", "Rule"], as_index=False)["RuleDisagreements"].sum()
    fig, ax = plt.subplots(figsize=(9.2, 5.1))
    for rule, frame in rounded.groupby("Rule"):
        ax.plot(
            frame["DisplayDecimals"],
            frame["RuleDisagreements"],
            marker="o",
            label=rule,
        )
    ax.set_xlabel("Displayed decimals")
    ax.set_ylabel("Raw-versus-rounded rule disagreements")
    ax.set_title("Known-truth pilot display-precision sensitivity\nRaw decisions remain authoritative")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "known_truth_display_precision.png", dpi=180)
    plt.close(fig)

    surface = threshold_summary.loc[
        threshold_summary["PersonMeasureSource"].eq(PERSON_MEASURE_WLE)
        & threshold_summary["Rule"].eq("either_upper")
    ].copy()
    rows: list[dict[str, object]] = []
    for identity, frame in surface.groupby(["ConditionId", "TruthGroup"]):
        canonical_row = frame.loc[frame["CanonicalThresholds"]]
        rows.append(
            {
                "Label": f"{identity[0]} | {identity[1]}",
                "Minimum": float(frame["MeanReplicateRate"].min()),
                "Maximum": float(frame["MeanReplicateRate"].max()),
                "Canonical": float(canonical_row["MeanReplicateRate"].iloc[0]),
            }
        )
    ranges = pd.DataFrame(rows).sort_values("Label").reset_index(drop=True)
    write_csv(ranges, output / "known_truth_threshold_range_projection.csv")
    fig, ax = plt.subplots(figsize=(10.8, 7.0))
    y = np.arange(len(ranges))
    ax.hlines(y, ranges["Minimum"], ranges["Maximum"], color="#777777", linewidth=4)
    ax.scatter(ranges["Minimum"], y, color="#2166ac", s=24)
    ax.scatter(ranges["Maximum"], y, color="#2166ac", s=24)
    ax.scatter(ranges["Canonical"], y, color="#b2182b", s=38, zorder=3)
    ax.set_yticks(y, ranges["Label"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Mean replicate-specific either-upper flag rate")
    ax.set_title("Frozen 125-triplet range\nDescriptive pilot only; no threshold selected")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "known_truth_threshold_ranges.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    plan = validate_plan(args.plan)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_known_truth_pilot_plan.json").write_text(
        args.plan.read_text(encoding="utf-8"), encoding="utf-8"
    )

    constants = plan["simulation_constants"]
    conditions = pd.DataFrame(plan["conditions"])
    manifest = build_replicate_manifest(
        plan["conditions"],
        replicates=int(constants["replicates_per_condition"]),
        base_seed=int(constants["base_seed"]),
    )
    write_csv(conditions, output / "known_truth_conditions.csv")
    write_csv(manifest, output / "known_truth_replicate_manifest.csv")

    ledger_rows: list[dict[str, object]] = []
    response_parts: list[pd.DataFrame] = []
    person_parts: list[pd.DataFrame] = []
    for row in manifest.itertuples(index=False):
        series = pd.Series(row._asdict())
        started = time.perf_counter()
        try:
            ledger, responses, persons = run_one(series, plan)
            ledger_rows.append(ledger)
            response_parts.append(responses)
            person_parts.append(persons)
        except Exception as exc:  # retain every attempted run
            ledger_rows.append(failure_ledger(series, exc, time.perf_counter() - started))
    ledger = pd.DataFrame(ledger_rows)
    responses = pd.concat(response_parts, ignore_index=True) if response_parts else pd.DataFrame()
    persons = pd.concat(person_parts, ignore_index=True) if person_parts else pd.DataFrame()
    write_csv(ledger, output / "known_truth_attempt_ledger.csv")
    write_csv(responses, output / "known_truth_responses.csv")
    write_csv(persons, output / "known_truth_person_fit.csv")

    if persons.empty:
        raise SystemExit("Known-truth pilot returned no Person-fit rows.")
    replay = replay_checks(manifest, plan, responses, persons)
    reproduction = canonical_recomputation(persons)
    identity_mismatches = truth_identity_mismatches(persons)
    triplets = threshold_triplets(plan)
    canonical_rates = replicate_fit_rule_rates(
        persons, threshold_triplets=[CANONICAL_FIT_THRESHOLDS]
    )
    canonical_summary = summarize_replicate_fit_rates(canonical_rates)
    wle_persons = persons.loc[
        persons["PersonMeasureSource"].eq(PERSON_MEASURE_WLE)
    ].copy()
    threshold_rates = replicate_fit_rule_rates(
        wle_persons, threshold_triplets=triplets
    )
    shuffled_rates = replicate_fit_rule_rates(
        wle_persons.sample(frac=1.0, random_state=20260810).reset_index(drop=True),
        threshold_triplets=triplets,
    )
    sort_columns = [
        "ConditionId",
        "Replicate",
        "PersonMeasureSource",
        "TruthGroup",
        "OverfitUpper",
        "AcceptableUpper",
        "NoisyUpper",
        "Rule",
    ]
    row_order_invariant = threshold_rates.sort_values(sort_columns).reset_index(drop=True).equals(
        shuffled_rates.sort_values(sort_columns).reset_index(drop=True)
    )
    threshold_summary = summarize_replicate_fit_rates(threshold_rates)
    rounding = rounded_rule_disagreements(
        persons,
        decimals=plan["display_precision"]["decimals"],
    )
    availability = availability_summary(persons)

    output_frames = {
        "known_truth_replay_checks.csv": replay,
        "known_truth_canonical_reproduction.csv": reproduction,
        "known_truth_canonical_replicate_rates.csv": canonical_rates,
        "known_truth_canonical_operating_summary.csv": canonical_summary,
        "known_truth_threshold_replicate_surface.csv": threshold_rates,
        "known_truth_threshold_operating_summary.csv": threshold_summary,
        "known_truth_display_precision.csv": rounding,
        "known_truth_availability_extremes.csv": availability,
    }
    for filename, frame in output_frames.items():
        write_csv(frame, output / filename)

    seed_coupling_passed = bool(
        manifest.groupby(["SeedGroup", "Replicate"])["Seed"].nunique().eq(1).all()
    )
    source_counts = persons.groupby("PersonMeasureSource").size().to_dict()
    gates = plan["frozen_pilot_gates"]
    unique_triplets = int(
        threshold_rates[["OverfitUpper", "AcceptableUpper", "NoisyUpper"]]
        .drop_duplicates()
        .shape[0]
    )
    contract_passed = bool(
        len(manifest) == int(gates["manifest_rows_must_equal"])
        and len(ledger) == int(gates["attempt_ledger_rows_must_equal"])
        and int(source_counts.get(PERSON_MEASURE_WLE, 0))
        == int(gates["generated_person_rows_per_measure_source_must_equal"])
        and int(source_counts.get(PERSON_MEASURE_TRUE, 0))
        == int(gates["generated_person_rows_per_measure_source_must_equal"])
        and unique_triplets == int(gates["threshold_triplets_must_equal"])
        and seed_coupling_passed
        and replay["Passed"].all()
        and int(ledger["CategorySupportViolations"].sum())
        == int(gates["category_support_violations_allowed"])
        and identity_mismatches
        == int(gates["person_identity_or_truth_mismatches_allowed"])
        and int(reproduction["ClassMismatches"].sum())
        == int(gates["raw_canonical_recomputation_mismatches_allowed"])
        and row_order_invariant
    )
    first_read_frame = first_read(
        contract_passed=contract_passed,
        ledger=ledger,
        canonical_summary=canonical_summary,
        rounding=rounding,
    )
    write_csv(first_read_frame, output / "known_truth_first_read_projection.csv")
    plot_outputs(canonical_summary, threshold_summary, rounding, output)

    primary = canonical_summary.loc[
        canonical_summary["PersonMeasureSource"].eq(PERSON_MEASURE_WLE)
        & canonical_summary["Rule"].eq("either_upper")
    ]
    clean = primary.loc[primary["TruthGroup"].eq("clean")]
    affected = primary.loc[primary["TruthGroup"].eq("affected")]
    three_decimal = rounding.loc[
        rounding["PersonMeasureSource"].eq(PERSON_MEASURE_WLE)
        & rounding["DisplayDecimals"].eq(3)
        & rounding["Rule"].eq("either_upper")
    ]
    decision = {
        "schema_version": "mfrm-cmle-wle-fit-known-truth-pilot-result-v1",
        "analysis_executed": True,
        "contract_passed": contract_passed,
        "overall_status": (
            "pilot_contract_complete_performance_withheld"
            if contract_passed
            else "pilot_contract_failed"
        ),
        "attempted_runs": int(len(ledger)),
        "returned_runs": int(ledger["Returned"].sum()),
        "failed_runs": int((~ledger["Returned"]).sum()),
        "response_rows": int(len(responses)),
        "person_rows_by_measure_source": {
            str(key): int(value) for key, value in source_counts.items()
        },
        "threshold_triplets": unique_triplets,
        "seed_coupling_passed": seed_coupling_passed,
        "same_seed_replay_passed": bool(replay["Passed"].all()),
        "category_support_violations": int(ledger["CategorySupportViolations"].sum()),
        "person_truth_identity_mismatches": identity_mismatches,
        "canonical_raw_class_mismatches": int(reproduction["ClassMismatches"].sum()),
        "threshold_surface_row_order_invariant": row_order_invariant,
        "pilot_clean_person_either_upper_rate_range": {
            "minimum": float(clean["MeanReplicateRate"].min()),
            "maximum": float(clean["MeanReplicateRate"].max()),
        },
        "pilot_affected_person_either_upper_detection_range": {
            "minimum": float(affected["MeanReplicateRate"].min()),
            "maximum": float(affected["MeanReplicateRate"].max()),
        },
        "three_decimal_either_upper_rule_disagreements": int(
            three_decimal["RuleDisagreements"].sum()
        ),
        "independent_monte_carlo_unit": "replicate",
        "independent_binomial_person_interpretation": False,
        "performance_value_success_gate_applied": False,
        "threshold_validated": False,
        "confirmatory_false_positive_or_power": False,
        "anchor_or_connectivity_evaluated": False,
        "cross_engine_repeated_fit_evaluated": False,
        "automatic_threshold_selection": False,
        "public_ui_integration_authorized": False,
        "registered_plan_sha256": sha256_file(args.plan),
        "known_truth_core_sha256": sha256_file(
            ROOT / "mfrm_app/fit_threshold_operating_characteristics.py"
        ),
        "runner_source_sha256": sha256_file(Path(__file__)),
        "response_evidence_sha256": sha256_file(output / "known_truth_responses.csv"),
    }
    (output / "known_truth_pilot_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    clean_rows = [
        f"- {row.ConditionId}: {row.MeanReplicateRate:.3f} "
        f"(MCSE {row.ReplicateRateMCSE:.3f})."
        for row in clean.itertuples(index=False)
    ]
    affected_rows = [
        f"- {row.ConditionId}: {row.MeanReplicateRate:.3f} "
        f"(MCSE {row.ReplicateRateMCSE:.3f})."
        for row in affected.itertuples(index=False)
    ]
    report = f"""# Known-truth CMLE-WLE Person MnSq pilot

## Decision

**{decision['overall_status']}.** The frozen machinery completed {decision['returned_runs']}/{decision['attempted_runs']} runs and retained {decision['response_rows']:,} response rows. This is a 20-replicate-per-condition design pilot. No observed false-positive or detection value was used as a success gate, and none is a confirmatory estimate.

## Canonical raw either-upper rule

The rule is `Infit > 1.50 or Outfit > 1.50`, evaluated from finite unrounded MnSq after fit-sample WLE scoring. Clean-Person pilot rates were:

{chr(10).join(clean_rows)}

Affected-Person mechanism-specific detection rates were:

{chr(10).join(affected_rows)}

These are means of replicate-specific rates. Persons within a replicate are clustered; the displayed MCSE is the between-replicate SD divided by sqrt(20), not a Person-level Wilson interval.

## Floating point and threshold dependence

- Canonical stored raw classes were reproduced with {decision['canonical_raw_class_mismatches']} mismatches.
- Three-decimal counterfactual display classification changed the either-upper rule for {decision['three_decimal_either_upper_rule_disagreements']} WLE Person rows across the full pilot.
- All 125 previously frozen threshold triplets were evaluated without selecting a favorable value. The retained range figure shows how both clean and affected rates move.

## Scope boundary

Phase A fixes the known working structural calibration. It therefore isolates Person scoring and fit-rule behavior but does not evaluate estimated-CMLE calibration error, sparse graph identification, anchor proportion/contamination, or repeated cross-engine fits. Python CMLE still lacks the native hard-anchor contract required for an honest anchor stress test; no post-hoc recentering substitute was used. The next extension must be registered before adding that capability or fitting the retained generated bytes with mfrmr, TAM, immer, and sirt.

Public traffic lights, p-values/ZSTD, threshold recommendations, sample-size claims, automatic estimator switches, and Streamlit integration remain withheld.
"""
    (output / "CMLE_WLE_FIT_KNOWN_TRUTH_PILOT_RESULTS.md").write_text(
        report, encoding="utf-8"
    )
    if not contract_passed:
        raise SystemExit("Known-truth Person-fit pilot contract failed.")


if __name__ == "__main__":
    main()
