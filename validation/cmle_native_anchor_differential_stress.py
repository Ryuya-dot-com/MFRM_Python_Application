#!/usr/bin/env python3
"""Run the prospectively registered native-CMLE differential-anchor stress."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
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

from mfrm_app.cmle import fit_cmle  # noqa: E402
from mfrm_app.cmle_person_scoring import score_cmle_persons_wle  # noqa: E402
from mfrm_app.cmle_wle_fit import compute_cmle_wle_person_fit  # noqa: E402
from mfrm_app.fit_threshold_operating_characteristics import fit_rule_flags  # noqa: E402
from mfrm_app.operating_characteristics import frame_fingerprint  # noqa: E402


DEFAULT_PLAN = ROOT / "validation/cmle_native_anchor_differential_stress_plan_20260810.json"
DEFAULT_OUTPUT = ROOT / "validation/cmle_native_anchor_differential_stress_20260810"


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
    if plan.get("study_id") != "cmle_native_anchor_differential_stress_v1":
        raise ValueError("Unexpected differential-anchor stress plan.")
    identities = dict(plan["input_identity"])
    identities.update(
        {
            value["path"]: value["sha256"]
            for value in plan["parent_evidence"].values()
        }
    )
    mismatches = [
        relative
        for relative, digest in identities.items()
        if not (ROOT / relative).exists() or sha256_file(ROOT / relative) != digest
    ]
    if mismatches:
        raise ValueError(f"Differential-anchor input identity failed: {mismatches}")
    if plan["integrity_gates"]["performance_value_success_gate"] is not None:
        raise ValueError("Debugging stress cannot have a performance success gate.")
    return plan


def select_responses(plan: dict[str, object]) -> pd.DataFrame:
    source_path = ROOT / "validation/cmle_wle_fit_known_truth_pilot_20260810/known_truth_responses.csv"
    source = pd.read_csv(source_path, float_precision="round_trip")
    contract = plan["retained_response_selection"]
    selected = source.loc[
        source["ConditionId"].isin(contract["condition_ids"])
        & source["Replicate"].isin(contract["replicates"])
    ].reset_index(drop=True)
    if len(selected) != int(contract["response_rows"]):
        raise ValueError("Selected response-row count differs from plan.")
    if selected["RunId"].nunique() != int(contract["run_ids"]):
        raise ValueError("Selected RunId count differs from plan.")
    if selected[["RunId", "Person"]].drop_duplicates().shape[0] != int(
        contract["persons"]
    ):
        raise ValueError("Selected Person count differs from plan.")
    return selected


def expected_kparams(plan: dict[str, object], anchors: int) -> int:
    return int(plan["fit_contract"]["expected_kparams_by_anchor_count"][str(anchors)])


def truth_for(plan: dict[str, object], row: pd.Series) -> float:
    truth = plan["truth"]
    if row["ParameterType"] == "Facet":
        return float(truth[str(row["Facet"])][str(row["Level"])])
    return float(truth["Step"][str(int(row["Step"]))])


def run_one(
    run_frame: pd.DataFrame,
    scenario: dict[str, object],
    plan: dict[str, object],
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame]:
    condition = str(run_frame["ConditionId"].iloc[0])
    replicate = int(run_frame["Replicate"].iloc[0])
    run_id = str(run_frame["RunId"].iloc[0])
    name = str(scenario["scenario"])
    anchors = pd.DataFrame(scenario["anchors"])
    hard_anchors = None if anchors.empty else anchors
    analysis = run_frame[["Person", "Rater", "Criterion", "ObservedCategory"]].copy()
    started = time.perf_counter()
    fit = fit_cmle(
        analysis,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="ObservedCategory",
        rating_min=0,
        rating_max=3,
        model="RSM",
        hard_anchors=hard_anchors,
        gtol=float(plan["fit_contract"]["gtol"]),
        maxiter=int(plan["fit_contract"]["maxiter"]),
        newton_polish_maxiter=int(plan["fit_contract"]["newton_polish_maxiter"]),
    )
    summary = fit["summary"].iloc[0]
    fitted_anchors = fit["facets"]["others"].loc[
        fit["facets"]["others"]["Anchored"]
    ]
    expected_anchors = {
        (str(row.Facet), str(row.Level)): float(row.Value)
        for row in anchors.itertuples(index=False)
    }
    anchor_exact = len(fitted_anchors) == len(expected_anchors) and all(
        float(row.Estimate) == expected_anchors[(str(row.Facet), str(row.Level))]
        and float(row.SE) == 0.0
        for row in fitted_anchors.itertuples(index=False)
    )
    expected_k = expected_kparams(plan, len(anchors))
    ledger = {
        "ConditionId": condition,
        "Replicate": replicate,
        "RunId": run_id,
        "ScenarioOrder": int(scenario["order"]),
        "Scenario": name,
        "ScenarioFamily": str(scenario["family"]),
        "AnchorShare": float(scenario["anchor_share"]),
        "Anchors": int(len(anchors)),
        "InputFingerprint": frame_fingerprint(analysis),
        "Attempted": True,
        "Returned": True,
        "FailureStage": "none",
        "FailureReason": "",
        "ResponseRows": int(len(analysis)),
        "Persons": int(analysis["Person"].nunique()),
        "AnchorExact": bool(anchor_exact),
        "ExpectedKParams": expected_k,
        "KParams": int(summary["KParams"]),
        "ExpectedKParamsExact": int(summary["KParams"]) == expected_k,
        "Eligible": bool(fit["audit"]["eligible"]),
        "Converged": bool(summary["Converged"]),
        "InferenceReady": bool(summary["InferenceReady"]),
        "ConditionalLogLik": float(summary["ConditionalLogLik"]),
        "ConditionalAIC": float(summary["ConditionalAIC"]),
        "InformationRank": int(summary["InformationRank"]),
        "InformationNullity": int(summary["InformationNullity"]),
        "InformationConditionNumber": float(summary["InformationConditionNumber"]),
        "GradientSupNorm": float(summary["GradientSupNorm"]),
        "ElapsedSeconds": float(time.perf_counter() - started),
    }

    structural = pd.concat(
        [fit["facets"]["others"], fit["steps"]], ignore_index=True
    )
    for column, value in reversed(
        list(
            {
                "ConditionId": condition,
                "Replicate": replicate,
                "RunId": run_id,
                "ScenarioOrder": int(scenario["order"]),
                "Scenario": name,
                "ScenarioFamily": str(scenario["family"]),
                "AnchorShare": float(scenario["anchor_share"]),
            }.items()
        )
    ):
        structural.insert(0, column, value)
    structural["Truth"] = structural.apply(lambda row: truth_for(plan, row), axis=1)
    structural["Error"] = structural["Estimate"] - structural["Truth"]
    structural["AbsoluteError"] = structural["Error"].abs()
    structural["SquaredError"] = structural["Error"] ** 2

    if not bool(summary["InferenceReady"]):
        return ledger, structural, pd.DataFrame()
    wle = score_cmle_persons_wle(fit)[
        ["Person", "Estimate", "StandardError", "ExtremeScorePattern", "Status"]
    ].rename(
        columns={
            "Estimate": "WLEEstimate",
            "StandardError": "WLEStandardError",
            "ExtremeScorePattern": "WLEExactExtreme",
            "Status": "WLEStatus",
        }
    )
    person_fit = compute_cmle_wle_person_fit(fit)["persons"]
    persons = wle.merge(
        person_fit[["Person", "Infit", "Outfit", "PersonFitReady"]],
        on="Person",
        how="left",
        validate="one_to_one",
    )
    truth = (
        run_frame.groupby("Person", sort=False)
        .agg(
            TrueTheta=("TrueTheta", "first"),
            Affected=("Affected", "first"),
            TruthGroup=("TruthGroup", "first"),
        )
        .reset_index()
    )
    persons = persons.merge(truth, on="Person", how="left", validate="one_to_one")
    for column, value in reversed(
        list(
            {
                "ConditionId": condition,
                "Replicate": replicate,
                "RunId": run_id,
                "ScenarioOrder": int(scenario["order"]),
                "Scenario": name,
                "ScenarioFamily": str(scenario["family"]),
                "AnchorShare": float(scenario["anchor_share"]),
            }.items()
        )
    ):
        persons.insert(0, column, value)
    persons["WLEError"] = persons["WLEEstimate"] - persons["TrueTheta"]
    raw = fit_rule_flags(persons["Infit"], persons["Outfit"])
    persons["FitEligibleRaw"] = raw["eligible"]
    persons["EitherUpperRaw"] = raw["either_upper"]
    persons["EitherDistortingRaw"] = raw["either_distorting"]
    persons["EitherOverfitRaw"] = raw["either_overfit"]
    for decimals in plan["raw_fit_contract"]["display_counterfactual_decimals"]:
        rounded = fit_rule_flags(
            np.round(persons["Infit"].to_numpy(dtype=float), int(decimals)),
            np.round(persons["Outfit"].to_numpy(dtype=float), int(decimals)),
        )["either_upper"]
        persons[f"EitherUpperRounded{int(decimals)}"] = rounded
        persons[f"RawRounded{int(decimals)}Disagreement"] = (
            rounded != raw["either_upper"]
        )
    persons["InfitDistanceTo1.5"] = np.abs(persons["Infit"] - 1.5)
    persons["OutfitDistanceTo1.5"] = np.abs(persons["Outfit"] - 1.5)
    persons["EitherDistanceTo1.5"] = persons[
        ["InfitDistanceTo1.5", "OutfitDistanceTo1.5"]
    ].min(axis=1)
    persons["DecisionsUseUnroundedValues"] = True
    return ledger, structural, persons


def failure_row(
    run_frame: pd.DataFrame,
    scenario: dict[str, object],
    plan: dict[str, object],
    exc: Exception,
    elapsed: float,
) -> dict[str, object]:
    anchors = len(scenario["anchors"])
    analysis = run_frame[["Person", "Rater", "Criterion", "ObservedCategory"]]
    return {
        "ConditionId": str(run_frame["ConditionId"].iloc[0]),
        "Replicate": int(run_frame["Replicate"].iloc[0]),
        "RunId": str(run_frame["RunId"].iloc[0]),
        "ScenarioOrder": int(scenario["order"]),
        "Scenario": str(scenario["scenario"]),
        "ScenarioFamily": str(scenario["family"]),
        "AnchorShare": float(scenario["anchor_share"]),
        "Anchors": anchors,
        "InputFingerprint": frame_fingerprint(analysis),
        "Attempted": True,
        "Returned": False,
        "FailureStage": "fit_or_downstream",
        "FailureReason": f"{type(exc).__name__}: {str(exc)[:500]}",
        "ResponseRows": int(len(analysis)),
        "Persons": int(analysis["Person"].nunique()),
        "AnchorExact": False,
        "ExpectedKParams": expected_kparams(plan, anchors),
        "KParams": np.nan,
        "ExpectedKParamsExact": False,
        "Eligible": False,
        "Converged": False,
        "InferenceReady": False,
        "ConditionalLogLik": np.nan,
        "ConditionalAIC": np.nan,
        "InformationRank": np.nan,
        "InformationNullity": np.nan,
        "InformationConditionNumber": np.nan,
        "GradientSupNorm": np.nan,
        "ElapsedSeconds": float(elapsed),
    }


def replicate_person_summary(persons: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group = [
        "ConditionId",
        "Replicate",
        "RunId",
        "ScenarioOrder",
        "Scenario",
        "ScenarioFamily",
        "AnchorShare",
        "TruthGroup",
    ]
    for identity, frame in persons.groupby(group, sort=False):
        error = frame["WLEError"].to_numpy(dtype=float)
        rows.append(
            {
                **dict(zip(group, identity)),
                "Persons": int(len(frame)),
                "WLEBias": float(np.mean(error)),
                "WLEMAE": float(np.mean(np.abs(error))),
                "WLERMSE": float(np.sqrt(np.mean(error**2))),
                "WLEExactExtremeRate": float(frame["WLEExactExtreme"].mean()),
                "EitherUpperRawRate": float(frame["EitherUpperRaw"].mean()),
                "EitherDistortingRawRate": float(frame["EitherDistortingRaw"].mean()),
                "EitherOverfitRawRate": float(frame["EitherOverfitRaw"].mean()),
                "Rounded3Disagreements": int(frame["RawRounded3Disagreement"].sum()),
                "Rounded6Disagreements": int(frame["RawRounded6Disagreement"].sum()),
                "Within0.0005Of1.5": int((frame["EitherDistanceTo1.5"] <= 0.0005).sum()),
                "Within0.001Of1.5": int((frame["EitherDistanceTo1.5"] <= 0.001).sum()),
            }
        )
    return pd.DataFrame(rows)


def aggregate_replicates(replicates: pd.DataFrame) -> pd.DataFrame:
    identities = [
        "ConditionId",
        "ScenarioOrder",
        "Scenario",
        "ScenarioFamily",
        "AnchorShare",
        "TruthGroup",
    ]
    metrics = [
        "WLEBias",
        "WLEMAE",
        "WLERMSE",
        "WLEExactExtremeRate",
        "EitherUpperRawRate",
        "EitherDistortingRawRate",
        "EitherOverfitRawRate",
    ]
    rows = []
    for identity, frame in replicates.groupby(identities, sort=False):
        row = {
            **dict(zip(identities, identity)),
            "Replicates": int(frame["Replicate"].nunique()),
            "Persons": int(frame["Persons"].sum()),
            "Rounded3Disagreements": int(frame["Rounded3Disagreements"].sum()),
            "Rounded6Disagreements": int(frame["Rounded6Disagreements"].sum()),
            "Within0.0005Of1.5": int(frame["Within0.0005Of1.5"].sum()),
            "Within0.001Of1.5": int(frame["Within0.001Of1.5"].sum()),
        }
        n = int(row["Replicates"])
        for metric in metrics:
            values = frame[metric].to_numpy(dtype=float)
            row[f"Mean{metric}"] = float(np.mean(values))
            row[f"SD{metric}"] = float(np.std(values, ddof=1)) if n > 1 else np.nan
            row[f"MCSE{metric}"] = (
                float(np.std(values, ddof=1) / np.sqrt(n)) if n > 1 else np.nan
            )
        rows.append(row)
    return pd.DataFrame(rows)


def structural_summaries(structural: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    group = [
        "ConditionId",
        "Replicate",
        "RunId",
        "ScenarioOrder",
        "Scenario",
        "ScenarioFamily",
        "AnchorShare",
        "ParameterType",
        "Facet",
        "Anchored",
    ]
    rows = []
    for identity, frame in structural.groupby(group, sort=False, dropna=False):
        error = frame["Error"].to_numpy(dtype=float)
        rows.append(
            {
                **dict(zip(group, identity)),
                "Estimates": int(len(frame)),
                "Bias": float(np.mean(error)),
                "MAE": float(np.mean(np.abs(error))),
                "RMSE": float(np.sqrt(np.mean(error**2))),
                "MeanSE": float(frame["SE"].mean()),
            }
        )
    replicate = pd.DataFrame(rows)
    aggregate_group = [
        "ConditionId",
        "ScenarioOrder",
        "Scenario",
        "ScenarioFamily",
        "AnchorShare",
        "ParameterType",
        "Facet",
        "Anchored",
    ]
    aggregate_rows = []
    for identity, frame in replicate.groupby(aggregate_group, sort=False, dropna=False):
        n = int(frame["Replicate"].nunique())
        row = {
            **dict(zip(aggregate_group, identity)),
            "Replicates": n,
            "Estimates": int(frame["Estimates"].sum()),
        }
        for metric in ["Bias", "MAE", "RMSE", "MeanSE"]:
            values = frame[metric].to_numpy(dtype=float)
            row[f"Mean{metric}"] = float(np.mean(values))
            row[f"SD{metric}"] = float(np.std(values, ddof=1)) if n > 1 else np.nan
            row[f"MCSE{metric}"] = (
                float(np.std(values, ddof=1) / np.sqrt(n)) if n > 1 else np.nan
            )
        aggregate_rows.append(row)
    return replicate, pd.DataFrame(aggregate_rows)


def transition_evidence(persons: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    key = ["ConditionId", "Replicate", "RunId", "Person"]
    reference_scenarios = {
        "unanchored": "a0_unanchored",
        "three_correct": "a3_correct_R01_R03_R05",
    }
    parts = []
    for label, scenario in reference_scenarios.items():
        reference = persons.loc[
            persons["Scenario"].eq(scenario), [*key, "EitherUpperRaw", "WLEEstimate"]
        ].rename(
            columns={
                "EitherUpperRaw": "ReferenceEitherUpperRaw",
                "WLEEstimate": "ReferenceWLEEstimate",
            }
        )
        compared = persons.merge(reference, on=key, how="left", validate="many_to_one")
        compared["Reference"] = label
        compared["RawFlagChanged"] = (
            compared["EitherUpperRaw"] != compared["ReferenceEitherUpperRaw"]
        )
        compared["WLEShift"] = compared["WLEEstimate"] - compared["ReferenceWLEEstimate"]
        parts.append(compared)
    row_evidence = pd.concat(parts, ignore_index=True)
    group = [
        "ConditionId",
        "ScenarioOrder",
        "Scenario",
        "ScenarioFamily",
        "TruthGroup",
        "Reference",
    ]
    rows = []
    for identity, frame in row_evidence.groupby(group, sort=False):
        before = frame["ReferenceEitherUpperRaw"].astype(bool)
        after = frame["EitherUpperRaw"].astype(bool)
        shift = frame["WLEShift"].to_numpy(dtype=float)
        rows.append(
            {
                **dict(zip(group, identity)),
                "Pairs": int(len(frame)),
                "BothUnflagged": int((~before & ~after).sum()),
                "ScenarioOnlyFlagged": int((~before & after).sum()),
                "ReferenceOnlyFlagged": int((before & ~after).sum()),
                "BothFlagged": int((before & after).sum()),
                "RawFlagChanged": int((before != after).sum()),
                "RawFlagChangeRate": float((before != after).mean()),
                "MeanWLEShift": float(np.mean(shift)),
                "MaxAbsWLEShift": float(np.max(np.abs(shift))),
            }
        )
    return row_evidence, pd.DataFrame(rows)


def save_figures(
    person_aggregate: pd.DataFrame,
    structural_aggregate: pd.DataFrame,
    transitions: pd.DataFrame,
    output: Path,
) -> None:
    scenarios = (
        person_aggregate[["ScenarioOrder", "Scenario"]]
        .drop_duplicates()
        .sort_values("ScenarioOrder")
    )
    order = scenarios["Scenario"].tolist()
    labels = [
        "0",
        "1 correct",
        "2 correct",
        "3 correct",
        "5 correct",
        "3 common",
        "3 differential",
        "3 focal",
    ]
    conditions = person_aggregate["ConditionId"].drop_duplicates().tolist()
    figure, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
    for axis, condition in zip(axes.flat, conditions):
        subset = person_aggregate.loc[
            person_aggregate["ConditionId"].eq(condition)
            & person_aggregate["TruthGroup"].eq("clean")
        ].set_index("Scenario").reindex(order)
        axis.bar(labels, subset["MeanWLERMSE"], color="#5C8374")
        axis.set_title(condition)
        axis.set_ylabel("Mean replicate WLE RMSE")
        axis.tick_params(axis="x", rotation=27)
        axis.grid(axis="y", alpha=0.25)
    figure.suptitle("Differential-anchor stress: clean-person WLE RMSE (10 replicates)")
    figure.savefig(output / "differential_anchor_wle_rmse.png", dpi=180)
    plt.close(figure)

    free_rater = structural_aggregate.loc[
        structural_aggregate["ParameterType"].eq("Facet")
        & structural_aggregate["Facet"].eq("Rater")
        & ~structural_aggregate["Anchored"].astype(bool)
    ]
    figure, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
    for axis, condition in zip(axes.flat, conditions):
        subset = free_rater.loc[free_rater["ConditionId"].eq(condition)].set_index(
            "Scenario"
        ).reindex(order)
        axis.bar(labels, subset["MeanRMSE"], color="#526D82")
        axis.set_title(condition)
        axis.set_ylabel("Mean replicate free-Rater RMSE")
        axis.tick_params(axis="x", rotation=27)
        axis.grid(axis="y", alpha=0.25)
    figure.suptitle("Differential-anchor stress: free-Rater recovery (10 replicates)")
    figure.savefig(output / "differential_anchor_rater_rmse.png", dpi=180)
    plt.close(figure)

    contrast_names = [
        "a3_common_plus_0.25",
        "a3_differential_minus0.25_0_plus0.25",
        "a3_focal_R03_plus_0.25",
    ]
    contrast_labels = ["common +.25", "differential", "focal R03 +.25"]
    subset = transitions.loc[
        transitions["Reference"].eq("three_correct")
        & transitions["Scenario"].isin(contrast_names)
    ].copy()
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    for axis, condition in zip(axes.flat, conditions):
        frame = subset.loc[subset["ConditionId"].eq(condition)]
        truth_groups = frame["TruthGroup"].drop_duplicates().tolist()
        width = 0.36 if len(truth_groups) > 1 else 0.6
        x = np.arange(len(contrast_names), dtype=float)
        for index, truth_group in enumerate(truth_groups):
            values = (
                frame.loc[frame["TruthGroup"].eq(truth_group)]
                .set_index("Scenario")
                .reindex(contrast_names)["RawFlagChangeRate"]
            )
            offset = (index - (len(truth_groups) - 1) / 2) * width
            axis.bar(x + offset, values, width=width, label=truth_group)
        axis.set_xticks(x, contrast_labels, rotation=20)
        axis.set_title(condition)
        axis.set_ylabel("Raw either-upper change rate vs 3 correct")
        axis.grid(axis="y", alpha=0.25)
        if len(truth_groups) > 1:
            axis.legend()
    figure.suptitle("Fit-decision sensitivity to anchor contamination (descriptive)")
    figure.savefig(output / "differential_anchor_fit_transitions.png", dpi=180)
    plt.close(figure)


def markdown_table(frame: pd.DataFrame) -> str:
    def cell(value: object) -> str:
        if pd.isna(value):
            return "NA"
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.6f}"
        return str(value).replace("|", "\\|")

    headers = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend(
        "| " + " | ".join(cell(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    plan = validate_plan(arguments.plan)
    output = arguments.output
    output.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(arguments.plan, output / "registered_differential_anchor_plan.json")

    selected = select_responses(plan)
    selected_path = output / "differential_anchor_retained_responses.csv"
    selected_fingerprint = frame_fingerprint(selected)
    write_csv(selected, selected_path)
    reloaded = pd.read_csv(selected_path, float_precision="round_trip")
    reloaded_fingerprint = frame_fingerprint(reloaded)

    groups = list(selected.groupby(["ConditionId", "Replicate"], sort=False))
    scenarios = plan["anchor_scenarios"]
    total = len(groups) * len(scenarios)
    ledgers: list[dict[str, object]] = []
    structural_parts: list[pd.DataFrame] = []
    person_parts: list[pd.DataFrame] = []
    attempt = 0
    for _, run_frame in groups:
        run_frame = run_frame.reset_index(drop=True)
        for scenario in scenarios:
            attempt += 1
            started = time.perf_counter()
            try:
                ledger, structural, persons = run_one(run_frame, scenario, plan)
                structural_parts.append(structural)
                if not persons.empty:
                    person_parts.append(persons)
            except Exception as exc:  # complete attempt accounting is required
                ledger = failure_row(
                    run_frame, scenario, plan, exc, time.perf_counter() - started
                )
            ledgers.append(ledger)
            if attempt % 8 == 0 or not ledger["Returned"]:
                write_csv(pd.DataFrame(ledgers), output / "_working_attempt_ledger.csv")
            print(
                f"[{attempt:03d}/{total:03d}] {ledger['RunId']} {ledger['Scenario']} "
                f"returned={ledger['Returned']} ready={ledger['InferenceReady']} "
                f"elapsed={ledger['ElapsedSeconds']:.2f}s",
                flush=True,
            )

    ledger = pd.DataFrame(ledgers)
    structural = pd.concat(structural_parts, ignore_index=True) if structural_parts else pd.DataFrame()
    persons = pd.concat(person_parts, ignore_index=True) if person_parts else pd.DataFrame()
    person_replicates = replicate_person_summary(persons)
    person_aggregate = aggregate_replicates(person_replicates)
    structural_replicates, structural_aggregate = structural_summaries(structural)
    transition_rows, transition_summary = transition_evidence(persons)

    raw_recomputed = fit_rule_flags(persons["Infit"], persons["Outfit"])[
        "either_upper"
    ]
    raw_passed = bool(
        np.array_equal(raw_recomputed, persons["EitherUpperRaw"].to_numpy(dtype=bool))
    )
    same_byte_passed = bool(
        ledger.groupby("RunId")["InputFingerprint"].nunique().eq(1).all()
    )
    attempts_passed = len(ledger) == int(plan["fit_contract"]["expected_attempts"]) and bool(
        ledger["Attempted"].all()
    )
    returned = ledger["Returned"].astype(bool)
    anchor_passed = bool(ledger.loc[returned, "AnchorExact"].all())
    kparams_passed = bool(ledger.loc[returned, "ExpectedKParamsExact"].all())
    source_passed = bool(
        len(selected) == int(plan["retained_response_selection"]["response_rows"])
        and selected_fingerprint == reloaded_fingerprint
    )
    common = transition_rows.loc[
        transition_rows["Reference"].eq("three_correct")
        & transition_rows["Scenario"].eq("a3_common_plus_0.25")
    ]
    common_shift_error = float(np.max(np.abs(common["WLEShift"] - 0.25)))
    common_fit_mismatches = int(common["RawFlagChanged"].sum())
    common_passed = bool(
        common_shift_error
        <= float(
            plan["integrity_gates"]["common_shift_negative_control"][
                "absolute_tolerance"
            ]
        )
        and common_fit_mismatches
        == int(
            plan["integrity_gates"]["common_shift_negative_control"][
                "raw_fit_mismatches"
            ]
        )
    )
    decision = {
        "schema_version": "mfrm-cmle-native-anchor-differential-stress-decision-v1",
        "study_id": plan["study_id"],
        "contract_passed": bool(
            attempts_passed
            and source_passed
            and same_byte_passed
            and anchor_passed
            and kparams_passed
            and raw_passed
            and common_passed
        ),
        "attempt_accounting_passed": attempts_passed,
        "retained_source_identity_passed": source_passed,
        "same_byte_pairing_passed": same_byte_passed,
        "anchor_exact_passed": anchor_passed,
        "expected_kparams_passed": kparams_passed,
        "raw_fit_recomputation_passed": raw_passed,
        "common_shift_negative_control_passed": common_passed,
        "common_shift_error_max": common_shift_error,
        "common_shift_raw_fit_mismatches": common_fit_mismatches,
        "attempts": int(len(ledger)),
        "returned": int(ledger["Returned"].sum()),
        "inference_ready": int(ledger["InferenceReady"].sum()),
        "run_ids": int(selected["RunId"].nunique()),
        "response_rows": int(len(selected)),
        "person_rows": int(len(persons)),
        "structural_rows": int(len(structural)),
        "rounded3_disagreements": int(persons["RawRounded3Disagreement"].sum()),
        "rounded6_disagreements": int(persons["RawRounded6Disagreement"].sum()),
        "within_0.0005_of_1.5": int((persons["EitherDistanceTo1.5"] <= 0.0005).sum()),
        "within_0.001_of_1.5": int((persons["EitherDistanceTo1.5"] <= 0.001).sum()),
        "performance_value_success_gate": None,
        "performance_claims_ready": False,
        "universal_anchor_share_ready": False,
        "cross_engine_claims_ready": False,
        "public_ui_ready": False,
        "decisions_use_unrounded_values": True,
        "selected_source_fingerprint": selected_fingerprint,
        "reloaded_source_fingerprint": reloaded_fingerprint,
    }

    outputs = {
        "differential_anchor_attempt_ledger.csv": ledger,
        "differential_anchor_structural_recovery.csv": structural,
        "differential_anchor_structural_replicates.csv": structural_replicates,
        "differential_anchor_structural_aggregate.csv": structural_aggregate,
        "differential_anchor_person_results.csv": persons,
        "differential_anchor_person_replicates.csv": person_replicates,
        "differential_anchor_person_aggregate.csv": person_aggregate,
        "differential_anchor_transition_rows.csv": transition_rows,
        "differential_anchor_transition_summary.csv": transition_summary,
    }
    for filename, frame in outputs.items():
        write_csv(frame, output / filename)
    save_figures(person_aggregate, structural_aggregate, transition_summary, output)
    decision["retained_response_sha256"] = sha256_file(selected_path)
    decision["output_sha256"] = {
        filename: sha256_file(output / filename) for filename in outputs
    }
    decision_path = output / "differential_anchor_decision.json"
    decision_path.write_text(
        json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    key_person = person_aggregate.loc[
        person_aggregate["TruthGroup"].eq("clean"),
        [
            "ConditionId",
            "Scenario",
            "Replicates",
            "MeanWLEBias",
            "MeanWLERMSE",
            "MeanEitherUpperRawRate",
            "MCSEEitherUpperRawRate",
            "Rounded3Disagreements",
        ],
    ]
    key_transition = transition_summary.loc[
        transition_summary["Reference"].eq("three_correct")
        & transition_summary["Scenario"].isin(
            [
                "a3_common_plus_0.25",
                "a3_differential_minus0.25_0_plus0.25",
                "a3_focal_R03_plus_0.25",
            ]
        ),
        [
            "ConditionId",
            "Scenario",
            "TruthGroup",
            "Pairs",
            "RawFlagChanged",
            "RawFlagChangeRate",
            "MeanWLEShift",
        ],
    ]
    report = "\n".join(
        [
            "# Native CMLE differential-anchor stress",
            "",
            "> Ten-replicate debugging evidence. This does not establish a universal anchor share, contamination detection power, cross-engine equivalence, or public UI readiness.",
            "",
            "## Contract",
            "",
            f"- Contract passed: `{decision['contract_passed']}`",
            f"- Returned/inference-ready: `{decision['returned']}/{decision['attempts']}` and `{decision['inference_ready']}/{decision['attempts']}`",
            f"- Common +0.25 shift maximum error: `{decision['common_shift_error_max']:.3e}`; raw fit mismatches: `{decision['common_shift_raw_fit_mismatches']}`",
            f"- Raw versus 3-decimal counterfactual disagreements: `{decision['rounded3_disagreements']}`; 6-decimal: `{decision['rounded6_disagreements']}`",
            "",
            "## Clean-person recovery and raw fit rates",
            "",
            markdown_table(key_person),
            "",
            "## Contamination contrasts versus three correct anchors",
            "",
            markdown_table(key_transition),
            "",
            "## Boundary",
            "",
            "Anchor share is inseparable here from which deterministic Rater levels are anchored. Common-mode contamination is an origin negative control and cannot be detected by fit. Differential/focal patterns can alter relative calibration and fit, but unchanged fit is not evidence of unbiased measures. All decisions use raw values; rounded results are counterfactual audits only.",
            "",
        ]
    )
    (output / "CMLE_NATIVE_ANCHOR_DIFFERENTIAL_STRESS.md").write_text(
        report, encoding="utf-8"
    )
    working = output / "_working_attempt_ledger.csv"
    if working.exists():
        working.unlink()
    print(json.dumps(decision, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
