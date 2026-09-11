#!/usr/bin/env python3
"""Run the prospectively amended native hard-anchor Phase-B smoke study."""

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


DEFAULT_AMENDMENT = ROOT / "validation/cmle_native_hard_anchor_smoke_amendment_20260810.json"
DEFAULT_OUTPUT = ROOT / "validation/cmle_native_hard_anchor_smoke_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, lineterminator="\n", float_format="%.17g")


def validate_amendment(path: Path) -> dict[str, object]:
    amendment = json.loads(path.read_text(encoding="utf-8"))
    if amendment.get("study_id") != "cmle_native_hard_anchor_phase_b_smoke_v1":
        raise ValueError("Unexpected hard-anchor smoke amendment.")
    expected = {
        amendment["parent_plan"]["sha256"]: ROOT / amendment["parent_plan"]["path"],
        amendment["retained_response_source"]["sha256"]: ROOT
        / amendment["retained_response_source"]["path"],
        amendment["retained_response_source"]["generator_plan_sha256"]: ROOT
        / amendment["retained_response_source"]["generator_plan_path"],
    }
    expected.update(
        {
            digest: ROOT / relative
            for relative, digest in amendment["implemented_input_identity"].items()
            if relative != "test_gate"
        }
    )
    mismatches = [
        str(source.relative_to(ROOT))
        for digest, source in expected.items()
        if not source.exists() or sha256_file(source) != digest
    ]
    if mismatches:
        raise ValueError(f"Hard-anchor smoke input identity failed: {mismatches}")
    return amendment


def selected_source(amendment: dict[str, object]) -> pd.DataFrame:
    source_contract = amendment["retained_response_source"]
    source = pd.read_csv(
        ROOT / source_contract["path"], float_precision="round_trip"
    )
    selection = source_contract["selection"]
    selected = source.loc[
        source["ConditionId"].isin(selection["condition_ids"])
        & source["Replicate"].isin(selection["replicates"])
    ].copy()
    selected.reset_index(drop=True, inplace=True)
    if len(selected) != int(selection["expected_response_rows"]):
        raise ValueError("Selected response-row count differs from the amendment.")
    if selected["RunId"].nunique() != int(selection["run_count"]):
        raise ValueError("Selected RunId count differs from the amendment.")
    return selected


def structural_truth(amendment: dict[str, object], row: pd.Series) -> float:
    truth = amendment["truth"]
    if row["ParameterType"] == "Facet":
        block = (
            truth["rater_effects"]
            if row["Facet"] == "Rater"
            else truth["criterion_effects"]
        )
        return float(block[str(row["Level"])])
    return float(truth["step_effects"][str(int(row["Step"]))])


def run_fit(
    run_frame: pd.DataFrame,
    *,
    scenario: dict[str, object],
    amendment: dict[str, object],
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame]:
    run_id = str(run_frame["RunId"].iloc[0])
    condition = str(run_frame["ConditionId"].iloc[0])
    replicate = int(run_frame["Replicate"].iloc[0])
    anchors = pd.DataFrame(scenario["anchors"])
    hard_anchors = None if anchors.empty else anchors
    scenario_name = str(scenario["scenario"])
    analysis = run_frame[["Person", "Rater", "Criterion", "ObservedCategory"]].copy()
    input_fingerprint = frame_fingerprint(analysis)
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
        gtol=float(amendment["fit_contract"]["gtol"]),
        maxiter=int(amendment["fit_contract"]["maxiter"]),
        newton_polish_maxiter=int(
            amendment["fit_contract"]["newton_polish_maxiter"]
        ),
    )
    summary = fit["summary"].iloc[0]
    expected_parameters = 10 if len(anchors) <= 1 else 8
    fitted_anchor_rows = fit["facets"]["others"].loc[
        fit["facets"]["others"]["Anchored"]
    ]
    expected_anchor_values = {
        (str(row.Facet), str(row.Level)): float(row.Value)
        for row in anchors.itertuples(index=False)
    }
    anchor_exact = len(fitted_anchor_rows) == len(expected_anchor_values) and all(
        float(row.Estimate) == expected_anchor_values[(str(row.Facet), str(row.Level))]
        and float(row.SE) == 0.0
        for row in fitted_anchor_rows.itertuples(index=False)
    )
    ledger = {
        "ConditionId": condition,
        "Replicate": replicate,
        "RunId": run_id,
        "Scenario": scenario_name,
        "InputFingerprint": input_fingerprint,
        "Attempted": True,
        "Returned": True,
        "FailureStage": "none",
        "FailureReason": "",
        "ResponseRows": int(len(analysis)),
        "Persons": int(analysis["Person"].nunique()),
        "Anchors": int(len(anchors)),
        "ExpectedKParams": expected_parameters,
        "KParams": int(summary["KParams"]),
        "ExpectedKParamsExact": int(summary["KParams"]) == expected_parameters,
        "AnchorExact": bool(anchor_exact),
        "Eligible": bool(fit["audit"]["eligible"]),
        "Converged": bool(summary["Converged"]),
        "InferenceReady": bool(summary["InferenceReady"]),
        "InformationRank": int(summary["InformationRank"]),
        "InformationNullity": int(summary["InformationNullity"]),
        "InformationConditionNumber": float(summary["InformationConditionNumber"]),
        "GradientSupNorm": float(summary["GradientSupNorm"]),
        "ElapsedSeconds": float(time.perf_counter() - started),
    }

    expanded = pd.concat(
        [fit["facets"]["others"], fit["steps"]], ignore_index=True
    )
    expanded.insert(0, "Scenario", scenario_name)
    expanded.insert(0, "RunId", run_id)
    expanded.insert(0, "Replicate", replicate)
    expanded.insert(0, "ConditionId", condition)
    expanded["Truth"] = expanded.apply(
        lambda row: structural_truth(amendment, row), axis=1
    )
    expanded["Error"] = expanded["Estimate"] - expanded["Truth"]
    expanded["AbsoluteError"] = expanded["Error"].abs()
    expanded["SquaredError"] = expanded["Error"] ** 2

    if not bool(summary["InferenceReady"]):
        return ledger, expanded, pd.DataFrame()
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
    fit_persons = compute_cmle_wle_person_fit(fit)["persons"]
    person = wle.merge(
        fit_persons[["Person", "Infit", "Outfit", "PersonFitReady"]],
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
    person = person.merge(truth, on="Person", how="left", validate="one_to_one")
    person.insert(0, "Scenario", scenario_name)
    person.insert(0, "RunId", run_id)
    person.insert(0, "Replicate", replicate)
    person.insert(0, "ConditionId", condition)
    person["WLEError"] = person["WLEEstimate"] - person["TrueTheta"]
    flags = fit_rule_flags(person["Infit"], person["Outfit"])
    person["EitherUpperRaw"] = flags["either_upper"]
    person["FitEligibleRaw"] = flags["eligible"]
    person["DecisionsUseUnroundedValues"] = True
    return ledger, expanded, person


def failure_ledger(
    run_frame: pd.DataFrame,
    scenario: dict[str, object],
    exc: Exception,
    elapsed: float,
) -> dict[str, object]:
    anchors = scenario["anchors"]
    return {
        "ConditionId": str(run_frame["ConditionId"].iloc[0]),
        "Replicate": int(run_frame["Replicate"].iloc[0]),
        "RunId": str(run_frame["RunId"].iloc[0]),
        "Scenario": str(scenario["scenario"]),
        "InputFingerprint": frame_fingerprint(
            run_frame[["Person", "Rater", "Criterion", "ObservedCategory"]]
        ),
        "Attempted": True,
        "Returned": False,
        "FailureStage": "fit_or_downstream",
        "FailureReason": f"{type(exc).__name__}: {str(exc)[:500]}",
        "ResponseRows": int(len(run_frame)),
        "Persons": int(run_frame["Person"].nunique()),
        "Anchors": int(len(anchors)),
        "ExpectedKParams": 10 if len(anchors) <= 1 else 8,
        "KParams": np.nan,
        "ExpectedKParamsExact": False,
        "AnchorExact": False,
        "Eligible": False,
        "Converged": False,
        "InferenceReady": False,
        "InformationRank": np.nan,
        "InformationNullity": np.nan,
        "InformationConditionNumber": np.nan,
        "GradientSupNorm": np.nan,
        "ElapsedSeconds": float(elapsed),
    }


def person_summary(persons: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for identity, frame in persons.groupby(
        ["ConditionId", "Scenario", "TruthGroup"], sort=False
    ):
        error = frame["WLEError"].to_numpy(dtype=float)
        rows.append(
            {
                "ConditionId": identity[0],
                "Scenario": identity[1],
                "TruthGroup": identity[2],
                "Runs": int(frame["RunId"].nunique()),
                "Persons": int(len(frame)),
                "WLEBias": float(np.mean(error)),
                "WLEMAE": float(np.mean(np.abs(error))),
                "WLERMSE": float(np.sqrt(np.mean(error**2))),
                "WLEExactExtremes": int(frame["WLEExactExtreme"].sum()),
                "WLEExactExtremeRate": float(frame["WLEExactExtreme"].mean()),
                "EitherUpperRawCount": int(frame["EitherUpperRaw"].sum()),
                "EitherUpperRawRate": float(frame["EitherUpperRaw"].mean()),
                "MeanInfit": float(frame["Infit"].mean()),
                "MeanOutfit": float(frame["Outfit"].mean()),
            }
        )
    return pd.DataFrame(rows)


def structural_summary(structural: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for identity, frame in structural.groupby(
        ["ConditionId", "Scenario", "ParameterType", "Facet", "Anchored"],
        sort=False,
        dropna=False,
    ):
        error = frame["Error"].to_numpy(dtype=float)
        rows.append(
            {
                "ConditionId": identity[0],
                "Scenario": identity[1],
                "ParameterType": identity[2],
                "Facet": identity[3],
                "Anchored": bool(identity[4]),
                "Runs": int(frame["RunId"].nunique()),
                "Estimates": int(len(frame)),
                "Bias": float(np.mean(error)),
                "MAE": float(np.mean(np.abs(error))),
                "RMSE": float(np.sqrt(np.mean(error**2))),
                "MeanSE": float(frame["SE"].mean()),
            }
        )
    return pd.DataFrame(rows)


def fit_transitions(persons: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    key = ["ConditionId", "Replicate", "RunId", "Person"]
    baseline = persons.loc[
        persons["Scenario"].eq("unanchored_sum_zero"),
        [*key, "EitherUpperRaw"],
    ].rename(columns={"EitherUpperRaw": "UnanchoredEitherUpperRaw"})
    anchored = persons.loc[
        ~persons["Scenario"].eq("unanchored_sum_zero")
    ].merge(baseline, on=key, how="left", validate="many_to_one")
    anchored["FlagChanged"] = (
        anchored["EitherUpperRaw"] != anchored["UnanchoredEitherUpperRaw"]
    )
    rows = []
    for identity, frame in anchored.groupby(
        ["ConditionId", "Scenario", "TruthGroup"], sort=False
    ):
        before = frame["UnanchoredEitherUpperRaw"].astype(bool)
        after = frame["EitherUpperRaw"].astype(bool)
        rows.append(
            {
                "ConditionId": identity[0],
                "Scenario": identity[1],
                "TruthGroup": identity[2],
                "Pairs": int(len(frame)),
                "BothUnflagged": int((~before & ~after).sum()),
                "AnchoredOnlyFlagged": int((~before & after).sum()),
                "UnanchoredOnlyFlagged": int((before & ~after).sum()),
                "BothFlagged": int((before & after).sum()),
                "Changed": int((before != after).sum()),
                "ChangeRate": float((before != after).mean()),
            }
        )
    return anchored, pd.DataFrame(rows)


def save_figures(
    person: pd.DataFrame,
    structural: pd.DataFrame,
    output: Path,
) -> None:
    scenario_order = [
        "unanchored_sum_zero",
        "one_of_six_raters_correct",
        "three_of_six_raters_correct",
        "three_of_six_raters_contaminated_plus_0.25",
    ]
    labels = ["0 anchors", "1 correct", "3 correct", "3 contaminated"]
    conditions = list(dict.fromkeys(person["ConditionId"].tolist()))
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    for axis, condition in zip(axes.flat, conditions):
        subset = person.loc[
            person["ConditionId"].eq(condition) & person["TruthGroup"].eq("clean")
        ]
        values = [
            np.sqrt(np.mean(subset.loc[subset["Scenario"].eq(scenario), "WLEError"] ** 2))
            for scenario in scenario_order
        ]
        axis.bar(labels, values, color=["#526D82", "#5C8374", "#84A98C", "#C1666B"])
        axis.set_title(condition)
        axis.set_ylabel("Clean-person WLE RMSE")
        axis.tick_params(axis="x", rotation=18)
        axis.grid(axis="y", alpha=0.25)
    figure.suptitle("Phase-B smoke: WLE recovery (2 replicates; descriptive only)")
    figure.savefig(output / "hard_anchor_wle_recovery.png", dpi=180)
    plt.close(figure)

    rater = structural.loc[
        structural["ParameterType"].eq("Facet")
        & structural["Facet"].eq("Rater")
        & ~structural["Anchored"].astype(bool)
    ]
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    for axis, condition in zip(axes.flat, conditions):
        subset = rater.loc[rater["ConditionId"].eq(condition)]
        values = [
            np.sqrt(np.mean(subset.loc[subset["Scenario"].eq(scenario), "Error"] ** 2))
            for scenario in scenario_order
        ]
        axis.bar(labels, values, color=["#526D82", "#5C8374", "#84A98C", "#C1666B"])
        axis.set_title(condition)
        axis.set_ylabel("Unanchored-rater RMSE")
        axis.tick_params(axis="x", rotation=18)
        axis.grid(axis="y", alpha=0.25)
    figure.suptitle("Phase-B smoke: free-rater recovery (2 replicates; descriptive only)")
    figure.savefig(output / "hard_anchor_rater_recovery.png", dpi=180)
    plt.close(figure)


def report_markdown(
    decision: dict[str, object],
    ledger: pd.DataFrame,
    person_summary_frame: pd.DataFrame,
    structural_summary_frame: pd.DataFrame,
    transitions: pd.DataFrame,
) -> str:
    def markdown_table(frame: pd.DataFrame) -> str:
        def cell(value: object) -> str:
            if pd.isna(value):
                return "NA"
            if isinstance(value, (float, np.floating)):
                return f"{float(value):.6f}"
            return str(value).replace("|", "\\|").replace("\n", " ")

        headers = [str(column).replace("|", "\\|") for column in frame.columns]
        lines = [
            "| " + " | ".join(headers) + " |",
            "| " + " | ".join("---" for _ in headers) + " |",
        ]
        lines.extend(
            "| " + " | ".join(cell(value) for value in row) + " |"
            for row in frame.itertuples(index=False, name=None)
        )
        return "\n".join(lines)

    key_person = person_summary_frame.loc[
        person_summary_frame["TruthGroup"].eq("clean"),
        ["ConditionId", "Scenario", "Persons", "WLEBias", "WLERMSE", "WLEExactExtremeRate", "EitherUpperRawRate"],
    ]
    key_rater = structural_summary_frame.loc[
        structural_summary_frame["ParameterType"].eq("Facet")
        & structural_summary_frame["Facet"].eq("Rater")
        & ~structural_summary_frame["Anchored"],
        ["ConditionId", "Scenario", "Estimates", "Bias", "RMSE", "MeanSE"],
    ]
    return "\n".join(
        [
            "# Native CMLE hard-anchor Phase-B smoke",
            "",
            "> Engineering smoke only: two replicates per condition. These values do not establish robustness, power, coverage, cross-engine equivalence, or UI release readiness.",
            "",
            "## Contract result",
            "",
            f"- Contract passed: `{decision['contract_passed']}`",
            f"- Attempts returned: `{int(ledger['Returned'].sum())}/{len(ledger)}`",
            f"- Inference-ready fits: `{int(ledger['InferenceReady'].sum())}/{len(ledger)}`",
            f"- Fixed anchors exact with SE zero: `{decision['anchor_exact_passed']}`",
            f"- Same-byte input identity: `{decision['same_byte_pairing_passed']}`",
            f"- Raw fit-flag recomputation: `{decision['raw_flag_recomputation_passed']}`",
            "",
            "## Clean-person WLE diagnostics",
            "",
            markdown_table(key_person),
            "",
            "## Free-rater recovery",
            "",
            markdown_table(key_rater),
            "",
            "## Same-byte raw either-upper transitions",
            "",
            markdown_table(transitions),
            "",
            "## Interpretation boundary",
            "",
            "The smoke isolates whether native hard anchors enter the exact conditional likelihood and propagate coherently. The contaminated anchors are deliberately wrong by +0.25. Any observed movement is descriptive direction-of-effect evidence, not a calibrated robustness limit. Raw Infit and Outfit values remain authoritative; rounded display values were not used.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--amendment", type=Path, default=DEFAULT_AMENDMENT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()

    amendment = validate_amendment(arguments.amendment)
    output = arguments.output
    output.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(arguments.amendment, output / "registered_hard_anchor_smoke_amendment.json")
    parent_plan = ROOT / amendment["parent_plan"]["path"]
    shutil.copyfile(parent_plan, output / "registered_hard_anchor_plan.json")

    selected = selected_source(amendment)
    selected_path = output / "hard_anchor_smoke_retained_responses.csv"
    selected_fingerprint = frame_fingerprint(selected)
    write_csv(selected, selected_path)
    reloaded = pd.read_csv(selected_path, float_precision="round_trip")
    reloaded_fingerprint = frame_fingerprint(reloaded)

    ledger_rows: list[dict[str, object]] = []
    structural_parts: list[pd.DataFrame] = []
    person_parts: list[pd.DataFrame] = []
    grouped = list(selected.groupby(["ConditionId", "Replicate"], sort=False))
    total = len(grouped) * len(amendment["frozen_anchor_scenarios"])
    attempt = 0
    for _, run_frame in grouped:
        run_frame = run_frame.reset_index(drop=True)
        for scenario in amendment["frozen_anchor_scenarios"]:
            attempt += 1
            started = time.perf_counter()
            try:
                ledger, structural, person = run_fit(
                    run_frame, scenario=scenario, amendment=amendment
                )
                structural_parts.append(structural)
                if not person.empty:
                    person_parts.append(person)
            except Exception as exc:  # complete failure accounting is part of the contract
                ledger = failure_ledger(
                    run_frame, scenario, exc, time.perf_counter() - started
                )
            ledger_rows.append(ledger)
            print(
                f"[{attempt:02d}/{total:02d}] {ledger['RunId']} {ledger['Scenario']} "
                f"returned={ledger['Returned']} ready={ledger['InferenceReady']} "
                f"elapsed={ledger['ElapsedSeconds']:.2f}s",
                flush=True,
            )

    ledger = pd.DataFrame(ledger_rows)
    structural = pd.concat(structural_parts, ignore_index=True) if structural_parts else pd.DataFrame()
    persons = pd.concat(person_parts, ignore_index=True) if person_parts else pd.DataFrame()
    person_summary_frame = person_summary(persons)
    structural_summary_frame = structural_summary(structural)
    transition_rows, transition_summary = fit_transitions(persons)

    recomputed = fit_rule_flags(persons["Infit"], persons["Outfit"])["either_upper"]
    raw_flag_passed = bool(
        np.array_equal(recomputed, persons["EitherUpperRaw"].to_numpy(dtype=bool))
    )
    same_byte_counts = ledger.groupby("RunId")["InputFingerprint"].nunique()
    same_byte_passed = bool(same_byte_counts.eq(1).all())
    attempts_passed = len(ledger) == 32 and bool(ledger["Attempted"].all())
    anchor_passed = bool(
        ledger.loc[ledger["Returned"] & ledger["Anchors"].gt(0), "AnchorExact"].all()
    )
    kparams_passed = bool(
        ledger.loc[ledger["Returned"], "ExpectedKParamsExact"].all()
    )
    source_passed = bool(
        len(selected) == 24000
        and selected["RunId"].nunique() == 8
        and selected_fingerprint == reloaded_fingerprint
    )
    decision = {
        "schema_version": "mfrm-cmle-native-hard-anchor-smoke-decision-v1",
        "study_id": amendment["study_id"],
        "contract_passed": bool(
            attempts_passed
            and anchor_passed
            and kparams_passed
            and source_passed
            and same_byte_passed
            and raw_flag_passed
        ),
        "attempt_accounting_passed": attempts_passed,
        "anchor_exact_passed": anchor_passed,
        "expected_free_parameter_count_passed": kparams_passed,
        "retained_source_identity_passed": source_passed,
        "same_byte_pairing_passed": same_byte_passed,
        "raw_flag_recomputation_passed": raw_flag_passed,
        "decisions_use_unrounded_values": True,
        "attempts": int(len(ledger)),
        "returned": int(ledger["Returned"].sum()),
        "inference_ready": int(ledger["InferenceReady"].sum()),
        "response_rows": int(len(selected)),
        "run_ids": int(selected["RunId"].nunique()),
        "person_result_rows": int(len(persons)),
        "structural_result_rows": int(len(structural)),
        "selected_source_fingerprint": selected_fingerprint,
        "reloaded_source_fingerprint": reloaded_fingerprint,
        "performance_success_gate": None,
        "performance_claims_ready": False,
        "cross_engine_claims_ready": False,
        "ui_release_ready": False,
        "interpretation": "two-replicate engineering smoke only",
    }

    outputs = {
        "hard_anchor_smoke_attempt_ledger.csv": ledger,
        "hard_anchor_smoke_structural_recovery.csv": structural,
        "hard_anchor_smoke_structural_summary.csv": structural_summary_frame,
        "hard_anchor_smoke_person_results.csv": persons,
        "hard_anchor_smoke_person_summary.csv": person_summary_frame,
        "hard_anchor_smoke_person_transitions.csv": transition_rows,
        "hard_anchor_smoke_transition_summary.csv": transition_summary,
    }
    for filename, frame in outputs.items():
        write_csv(frame, output / filename)
    save_figures(persons, structural, output)
    decision["output_sha256"] = {
        filename: sha256_file(output / filename) for filename in outputs
    }
    decision["retained_response_file_sha256"] = sha256_file(selected_path)
    decision_path = output / "hard_anchor_smoke_decision.json"
    decision_path.write_text(
        json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    report = report_markdown(
        decision, ledger, person_summary_frame, structural_summary_frame, transition_summary
    )
    (output / "CMLE_NATIVE_HARD_ANCHOR_SMOKE.md").write_text(report, encoding="utf-8")
    print(json.dumps(decision, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
