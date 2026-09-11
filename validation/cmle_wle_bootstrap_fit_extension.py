#!/usr/bin/env python3
"""Run the prospectively registered CMLE-WLE bootstrap Person-fit extension."""

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

from mfrm_app.cmle import fit_cmle  # noqa: E402
from mfrm_app.cmle_wle_bootstrap import (  # noqa: E402
    FIXED_SCORE_LANE,
    JOINT_PLUGIN_LANE,
    run_cmle_wle_bootstrap,
)


DEFAULT_PLAN = ROOT / "validation/cmle_wle_bootstrap_fit_extension_plan_20260810.json"
ORIGINAL_DIR = ROOT / "validation/cmle_wle_bootstrap_pilot_20260810"
FAILED_DIR = ROOT / "validation/cmle_wle_bootstrap_fit_extension_20260810"
DEFAULT_AMENDMENT = (
    ROOT / "validation/cmle_wle_bootstrap_fit_extension_amendment_20260810.json"
)
DEFAULT_OUTPUT = ROOT / "validation/cmle_wle_bootstrap_fit_extension_corrected_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def directory_inventory(directory: Path) -> dict[str, str]:
    return {
        str(path.relative_to(directory)): sha256_file(path)
        for path in sorted(directory.rglob("*"))
        if path.is_file()
    }


def validate_plan(plan_path: Path) -> dict[str, object]:
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != "mfrm-cmle-wle-bootstrap-fit-extension-plan-v1":
        raise ValueError("Unexpected bootstrap fit-extension plan schema.")
    original_decision_path = ORIGINAL_DIR / "bootstrap_pilot_decision.json"
    original_decision = json.loads(original_decision_path.read_text(encoding="utf-8"))
    identity = plan["input_identity"]
    checks = {
        "bootstrap_core_pre_extension_sha256": (
            str(original_decision["bootstrap_core_sha256"])
        ),
        "person_fit_core_sha256": sha256_file(ROOT / "mfrm_app/cmle_wle_fit.py"),
        "person_fit_registered_plan_sha256": sha256_file(
            ROOT / "validation/cmle_wle_person_fit_plan_20260810.json"
        ),
        "person_fit_parity_decision_sha256": sha256_file(
            ROOT
            / "validation/cmle_wle_person_fit_20260810/person_fit_parity_decision.json"
        ),
        "original_bootstrap_plan_sha256": sha256_file(
            ROOT / "validation/cmle_wle_bootstrap_plan_20260810.json"
        ),
        "original_bootstrap_pilot_decision_sha256": sha256_file(
            original_decision_path
        ),
        "original_bootstrap_fixture_sha256": sha256_file(
            ORIGINAL_DIR / "bootstrap_fixture.csv"
        ),
        "decision_stability_sha256": sha256_file(
            ROOT / "mfrm_app/decision_stability.py"
        ),
    }
    failed = [key for key, actual in checks.items() if str(identity[key]) != actual]
    if failed:
        raise ValueError(f"Bootstrap fit-extension input identity failed: {failed}")
    if plan["frozen_gates"].get(
        "minimum_person_fit_available_share_for_public_promotion"
    ) is not None:
        raise ValueError("The descriptive extension must not add a public promotion gate.")
    return plan


def validate_amendment(amendment_path: Path, plan_path: Path) -> dict[str, object]:
    amendment = json.loads(amendment_path.read_text(encoding="utf-8"))
    if amendment.get("schema_version") != "mfrm-cmle-wle-bootstrap-fit-extension-amendment-v1":
        raise ValueError("Unexpected bootstrap fit-extension amendment schema.")
    failed_identity = amendment["failed_run_identity"]
    failed_decision_path = FAILED_DIR / "bootstrap_fit_extension_decision.json"
    failed_decision = json.loads(failed_decision_path.read_text(encoding="utf-8"))
    checks = {
        "original_plan_sha256": sha256_file(plan_path),
        "decision_sha256": sha256_file(failed_decision_path),
        "runner_sha256": str(failed_decision["runner_source_sha256"]),
        "extension_person_draws_sha256": sha256_file(
            FAILED_DIR / "bootstrap_fit_person_draws.csv"
        ),
        "original_person_draws_sha256": sha256_file(
            ORIGINAL_DIR / "bootstrap_person_draws.csv"
        ),
    }
    failed = []
    if str(amendment["original_plan_sha256"]) != checks["original_plan_sha256"]:
        failed.append("original_plan_sha256")
    for key in (
        "decision_sha256",
        "runner_sha256",
        "extension_person_draws_sha256",
        "original_person_draws_sha256",
    ):
        if str(failed_identity[key]) != checks[key]:
            failed.append(key)
    if failed:
        raise ValueError(f"Bootstrap fit-extension amendment identity failed: {failed}")
    correction = amendment["prospective_correction"]
    if bool(correction["exact_zero_gate_changed"]) or bool(
        correction["tolerance_changed"]
    ):
        raise ValueError("Corrected rerun must retain the original exact-zero gate.")
    return amendment


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(
        path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )


def wilson_interval(successes: int, attempts: int) -> tuple[float, float]:
    if attempts <= 0:
        return np.nan, np.nan
    z = 1.959963984540054
    proportion = successes / attempts
    denominator = 1.0 + z * z / attempts
    center = (proportion + z * z / (2.0 * attempts)) / denominator
    half = z * np.sqrt(
        proportion * (1.0 - proportion) / attempts
        + z * z / (4.0 * attempts**2)
    ) / denominator
    return float(center - half), float(center + half)


def fit_fixture(frame: pd.DataFrame, model: str) -> dict[str, object]:
    return fit_cmle(
        frame,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model=model,
        step_facet="Criterion" if model == "PCM" else None,
        gtol=1e-8,
    )


def summarize_person_fit(
    results: list[dict[str, object]], attempts: int
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for result in results:
        model = str(result["summary"].iloc[0]["Model"])
        lane = str(result["summary"].iloc[0]["Lane"])
        baseline = result["baseline_persons"].assign(
            Person=lambda frame: frame["Person"].astype(str)
        )
        draws = result["person_draws"].assign(
            Person=lambda frame: frame["Person"].astype(str)
        )
        for baseline_row in baseline.itertuples(index=False):
            person = str(baseline_row.Person)
            person_draws = draws.loc[draws["Person"].eq(person)]
            ready = person_draws.loc[
                person_draws["PersonFitReady"].fillna(False).astype(bool)
            ]
            rows.append(
                {
                    "Model": model,
                    "Lane": lane,
                    "Person": person,
                    "BaselineExtremeScorePattern": bool(
                        baseline_row.ExtremeScorePattern
                    ),
                    "BaselineInfit": float(baseline_row.Infit),
                    "BaselineOutfit": float(baseline_row.Outfit),
                    "BaselineInfitClass": str(baseline_row.InfitClass),
                    "BaselineOutfitClass": str(baseline_row.OutfitClass),
                    "BaselineRawDisplayMismatch": bool(
                        not baseline_row.InfitDisplayDecisionConsistent
                        or not baseline_row.OutfitDisplayDecisionConsistent
                    ),
                    "Attempts": int(attempts),
                    "SuccessfulFitReplicates": int(len(ready)),
                    "UnavailableFitReplicates": int(attempts - len(ready)),
                    "InfitZoneTransitions": int(
                        ready["InfitZoneChanged"].astype(bool).sum()
                    ),
                    "OutfitZoneTransitions": int(
                        ready["OutfitZoneChanged"].astype(bool).sum()
                    ),
                    "AnyFitZoneTransitions": int(
                        ready["AnyFitZoneTransition"].astype(bool).sum()
                    ),
                    "RawDisplayMismatches": int(
                        ready["RawDisplayMismatch"].astype(bool).sum()
                    ),
                    "InfitBoundaryReplicates": int(
                        ready["InfitBoundaryStatus"].ne("stable").sum()
                    ),
                    "OutfitBoundaryReplicates": int(
                        ready["OutfitBoundaryStatus"].ne("stable").sum()
                    ),
                    "InfitBootstrapMedian": float(ready["Infit"].median())
                    if len(ready)
                    else np.nan,
                    "OutfitBootstrapMedian": float(ready["Outfit"].median())
                    if len(ready)
                    else np.nan,
                    "InfitBootstrapMinimum": float(ready["Infit"].min())
                    if len(ready)
                    else np.nan,
                    "InfitBootstrapMaximum": float(ready["Infit"].max())
                    if len(ready)
                    else np.nan,
                    "OutfitBootstrapMinimum": float(ready["Outfit"].min())
                    if len(ready)
                    else np.nan,
                    "OutfitBootstrapMaximum": float(ready["Outfit"].max())
                    if len(ready)
                    else np.nan,
                    "TransitionEvidenceDescriptiveOnly": True,
                    "ThresholdOptimalityValidated": False,
                }
            )
    return pd.DataFrame(rows)


def transition_matrix(person_draws: pd.DataFrame) -> pd.DataFrame:
    ready = person_draws.loc[
        person_draws["PersonFitReady"].fillna(False).astype(bool)
    ]
    parts: list[pd.DataFrame] = []
    for statistic in ("Infit", "Outfit"):
        baseline_column = f"Baseline{statistic}Class"
        replicate_column = f"{statistic}Class"
        part = (
            ready.groupby(
                ["Model", "Lane", baseline_column, replicate_column],
                as_index=False,
                dropna=False,
            )
            .size()
            .rename(
                columns={
                    baseline_column: "BaselineClass",
                    replicate_column: "ReplicateClass",
                    "size": "PersonReplicates",
                }
            )
        )
        part.insert(2, "Statistic", statistic)
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def wle_identity_check(
    persisted_extension_draws: Path,
) -> tuple[pd.DataFrame, dict[str, object]]:
    original = pd.read_csv(ORIGINAL_DIR / "bootstrap_person_draws.csv")
    person_draws = pd.read_csv(persisted_extension_draws)
    keys = ["Lane", "Model", "Replicate", "Person"]
    numeric = [
        "Estimate",
        "StandardError",
        "AdjustedScoreResidual",
        "BaselineEstimate",
        "BootstrapMinusBaseline",
    ]
    pairs = original[keys + numeric].merge(
        person_draws[keys + numeric],
        on=keys,
        how="outer",
        suffixes=("Original", "Extension"),
        validate="one_to_one",
        indicator=True,
    )
    maximum = 0.0
    for column in numeric:
        difference = (
            pairs[f"{column}Original"] - pairs[f"{column}Extension"]
        ).abs()
        pairs[f"Abs{column}Difference"] = difference
        finite = difference[np.isfinite(difference)]
        if len(finite):
            maximum = max(maximum, float(finite.max()))
    summary = {
        "OriginalRows": int(len(original)),
        "ExtensionRows": int(len(person_draws)),
        "BothPresentRows": int(pairs["_merge"].eq("both").sum()),
        "IdentityMismatches": int(pairs["_merge"].ne("both").sum()),
        "MaximumAbsoluteNumericDifference": maximum,
    }
    return pairs.drop(columns="_merge"), summary


def first_read_projection(
    model_summary: pd.DataFrame,
    person_summary: pd.DataFrame,
    person_draws: pd.DataFrame,
) -> pd.DataFrame:
    fit_share = float(model_summary["PersonFitAvailableShare"].min())
    any_transitions = int(person_summary["AnyFitZoneTransitions"].sum())
    raw_display = int(person_summary["RawDisplayMismatches"].sum())
    boundary = int(
        person_summary[
            ["InfitBoundaryReplicates", "OutfitBoundaryReplicates"]
        ].to_numpy().sum()
    )
    return pd.DataFrame(
        [
            {
                "CardOrder": 1,
                "CardId": "person_fit_parity",
                "Status": "research_ready",
                "DisplayValue": "Python-sirt parity passed",
                "Interpretation": "Untrimmed fixed-theta Person MnSq translation is validated for RSM/PCM.",
                "NextAction": "Keep TAM trimmed Outfit labelled as a different default estimand.",
            },
            {
                "CardOrder": 2,
                "CardId": "bootstrap_fit_readiness",
                "Status": "research_ready" if fit_share == 1.0 else "caution",
                "DisplayValue": f"minimum fit availability {100.0 * fit_share:.1f}%",
                "Interpretation": "Denominator is every attempted replicate, not successful fits only.",
                "NextAction": "Inspect ledger failures and Wilson bounds.",
            },
            {
                "CardOrder": 3,
                "CardId": "fit_zone_transitions",
                "Status": "caution",
                "DisplayValue": f"{any_transitions} Person-replicate transitions",
                "Interpretation": "Transitions are descriptive sensitivity under each bootstrap lane.",
                "NextAction": "Inspect Infit and Outfit separately; do not pool the two lanes.",
            },
            {
                "CardOrder": 4,
                "CardId": "floating_display_stability",
                "Status": "research_ready" if raw_display == 0 and boundary == 0 else "caution",
                "DisplayValue": f"{boundary} boundary values; {raw_display} raw/display mismatches",
                "Interpretation": "All transition decisions used unrounded MnSq.",
                "NextAction": "Show raw value and boundary status whenever caution is present.",
            },
            {
                "CardOrder": 5,
                "CardId": "fit_reference_distribution",
                "Status": "withheld",
                "DisplayValue": "ZSTD and p-values unavailable",
                "Interpretation": "No validated CMLE-WLE finite-sample reference distribution exists here.",
                "NextAction": "Run a prospectively registered repeated-truth calibration study.",
            },
            {
                "CardOrder": 6,
                "CardId": "interval_claim",
                "Status": "withheld",
                "DisplayValue": "fit intervals and total SE withheld",
                "Interpretation": "This extension does not qualify coverage or variance composition.",
                "NextAction": "Keep conditional WLE SE and bootstrap sensitivity distinct.",
            },
            {
                "CardOrder": 7,
                "CardId": "public_ui",
                "Status": "withheld",
                "DisplayValue": "Streamlit integration not authorized",
                "Interpretation": "Research parity and an implementation pilot are not public validation.",
                "NextAction": "Complete repeated-truth, sparse/anchor, threshold, and usability gates.",
            },
        ]
    )


def plot_fit_extension(
    person_draws: pd.DataFrame,
    person_summary: pd.DataFrame,
    output: Path,
) -> None:
    lane_label = {
        FIXED_SCORE_LANE: "fixed-score",
        JOINT_PLUGIN_LANE: "joint plug-in",
    }
    ready = person_draws.loc[
        person_draws["PersonFitReady"].fillna(False).astype(bool)
    ].copy()
    rates: list[dict[str, object]] = []
    for (model, lane), frame in ready.groupby(["Model", "Lane"], sort=False):
        for statistic in ("Infit", "Outfit"):
            rates.append(
                {
                    "Model": model,
                    "Lane": lane,
                    "Statistic": statistic,
                    "TransitionRate": float(frame[f"{statistic}ZoneChanged"].mean()),
                    "Transitions": int(frame[f"{statistic}ZoneChanged"].sum()),
                    "PersonReplicates": int(len(frame)),
                }
            )
    rate_frame = pd.DataFrame(rates)
    write_csv(rate_frame, output / "bootstrap_fit_transition_rates.csv")

    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    labels = [
        f"{row.Model}\n{lane_label[row.Lane]}\n{row.Statistic}"
        for row in rate_frame.itertuples()
    ]
    colors = ["#1f77b4" if row.Statistic == "Infit" else "#d95f02" for row in rate_frame.itertuples()]
    bars = ax.bar(np.arange(len(rate_frame)), rate_frame["TransitionRate"], color=colors)
    for bar, row in zip(bars, rate_frame.itertuples()):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.005,
            f"{row.Transitions}/{row.PersonReplicates}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax.set_xticks(np.arange(len(labels)), labels, fontsize=8)
    ax.set_ylabel("Raw-class transition share")
    ax.set_title("CMLE-WLE Person fit-zone sensitivity\nInfit and Outfit shown separately; lanes are not pooled")
    ax.set_ylim(0, max(0.05, float(rate_frame["TransitionRate"].max()) * 1.22))
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "bootstrap_fit_zone_transition_rates.png", dpi=180)
    plt.close(fig)

    boundary_rows: list[dict[str, object]] = []
    for statistic in ("Infit", "Outfit"):
        counts = ready[f"{statistic}BoundaryStatus"].value_counts()
        for status in ("stable", "display_rounding_boundary", "numerical_boundary"):
            boundary_rows.append(
                {
                    "Statistic": statistic,
                    "BoundaryStatus": status,
                    "PersonReplicates": int(counts.get(status, 0)),
                }
            )
    boundary_frame = pd.DataFrame(boundary_rows)
    write_csv(boundary_frame, output / "bootstrap_fit_boundary_status.csv")
    pivot = boundary_frame.pivot(
        index="Statistic", columns="BoundaryStatus", values="PersonReplicates"
    ).fillna(0)
    pivot = pivot.reindex(
        columns=["stable", "display_rounding_boundary", "numerical_boundary"],
        fill_value=0,
    )
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    colors = ["#4daf4a", "#ffbf00", "#d73027"]
    bottom = np.zeros(len(pivot))
    for column, color in zip(pivot.columns, colors):
        values = pivot[column].to_numpy(dtype=float)
        ax.bar(pivot.index, values, bottom=bottom, label=column, color=color)
        bottom += values
    mismatch_count = int(ready["RawDisplayMismatch"].astype(bool).sum())
    ax.set_ylabel("Successful Person-replicates")
    ax.set_title(
        "Floating-point and display-boundary audit\n"
        f"Raw/display classification mismatches: {mismatch_count}"
    )
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "bootstrap_fit_boundary_stability.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--amendment", type=Path, default=DEFAULT_AMENDMENT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    plan = validate_plan(args.plan)
    amendment = validate_amendment(args.amendment, args.plan)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    original_before = directory_inventory(ORIGINAL_DIR)
    failed_before = directory_inventory(FAILED_DIR)
    (output / "original_pilot_inventory_before.json").write_text(
        json.dumps(original_before, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (output / "registered_fit_extension_plan.json").write_text(
        args.plan.read_text(encoding="utf-8"), encoding="utf-8"
    )
    (output / "registered_fit_extension_amendment.json").write_text(
        args.amendment.read_text(encoding="utf-8"), encoding="utf-8"
    )

    matrix = plan["run_matrix"]
    fixture = pd.read_csv(ORIGINAL_DIR / "bootstrap_fixture.csv")
    replicates = int(matrix["replicates_per_model_lane"])
    seed = int(matrix["primary_seed"])
    results: list[dict[str, object]] = []
    baseline_rows: list[dict[str, object]] = []
    for model in matrix["models"]:
        fit = fit_fixture(fixture, str(model))
        fit_row = fit["summary"].iloc[0]
        baseline_rows.append(
            {
                "Model": str(model),
                "Rows": int(len(fit["design"].data)),
                "Persons": int(fit_row["PersonsTotal"]),
                "ExactExtremes": int(fit_row["PersonsExtreme"]),
                "MissingnessPatterns": int(fit_row["Patterns"]),
                "CMLEInferenceReady": bool(fit_row["InferenceReady"]),
            }
        )
        for lane in matrix["lanes"]:
            def progress(row: dict[str, object], *, model_name=str(model), lane_name=str(lane)) -> None:
                replicate = int(row["Replicate"])
                if replicate == 1 or replicate % 50 == 0 or replicate == replicates:
                    print(
                        f"{model_name} {lane_name}: {replicate}/{replicates} "
                        f"fit={bool(row['PersonFitAvailable'])}",
                        flush=True,
                    )

            results.append(
                run_cmle_wle_bootstrap(
                    fit,
                    lane=str(lane),
                    n_replicates=replicates,
                    seed=seed,
                    gtol=1e-8,
                    progress_callback=progress,
                )
            )

    baseline = pd.DataFrame(baseline_rows)
    ledger = pd.concat([result["ledger"] for result in results], ignore_index=True)
    generator = pd.concat(
        [result["generator_audit"] for result in results], ignore_index=True
    )
    model_summary = pd.concat(
        [result["summary"] for result in results], ignore_index=True
    )
    person_draws = pd.concat(
        [result["person_draws"] for result in results], ignore_index=True
    )
    baseline_persons = pd.concat(
        [
            result["baseline_persons"].assign(
                Model=str(result["summary"].iloc[0]["Model"]),
                Lane=str(result["summary"].iloc[0]["Lane"]),
            )
            for result in results
        ],
        ignore_index=True,
    )
    person_summary = summarize_person_fit(results, replicates)
    transitions = transition_matrix(person_draws)
    for index, row in model_summary.iterrows():
        lower, upper = wilson_interval(
            int(row["PersonFitAvailableReplicates"]),
            int(row["AttemptedReplicates"]),
        )
        model_summary.loc[index, "PersonFitAvailableWilson95Lower"] = lower
        model_summary.loc[index, "PersonFitAvailableWilson95Upper"] = upper
    failure_summary = (
        ledger.assign(FailureStage=ledger["FailureStage"].replace("", "none"))
        .groupby(["Model", "Lane", "FailureStage"], as_index=False)
        .size()
    )
    # The original pilot evidence is persisted at %.17g.  Persist and read the
    # extension under the identical contract before enforcing the frozen exact-
    # zero gate; comparing a parsed CSV against an in-memory pre-serialization
    # float is an asymmetric representation check.
    persisted_person_draws = output / "bootstrap_fit_person_draws.csv"
    write_csv(person_draws, persisted_person_draws)
    wle_pairs, wle_identity = wle_identity_check(persisted_person_draws)
    first_read = first_read_projection(model_summary, person_summary, person_draws)

    artifacts = {
        "bootstrap_fit_baseline_models.csv": baseline,
        "bootstrap_fit_baseline_persons.csv": baseline_persons,
        "bootstrap_fit_replicate_ledger.csv": ledger,
        "bootstrap_fit_generator_audit.csv": generator,
        "bootstrap_fit_model_summary.csv": model_summary,
        "bootstrap_fit_failure_summary.csv": failure_summary,
        "bootstrap_fit_person_draws.csv": person_draws,
        "bootstrap_fit_person_summary.csv": person_summary,
        "bootstrap_fit_transition_matrix.csv": transitions,
        "bootstrap_fit_wle_identity_pairs.csv": wle_pairs,
        "bootstrap_fit_first_read_projection.csv": first_read,
    }
    for filename, frame in artifacts.items():
        write_csv(frame, output / filename)
    plot_fit_extension(person_draws, person_summary, output)

    original_after = directory_inventory(ORIGINAL_DIR)
    failed_after = directory_inventory(FAILED_DIR)
    (output / "original_pilot_inventory_after.json").write_text(
        json.dumps(original_after, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    changed_original_files = sorted(
        set(original_before) ^ set(original_after)
        | {
            path
            for path in set(original_before) & set(original_after)
            if original_before[path] != original_after[path]
        }
    )
    changed_failed_files = sorted(
        set(failed_before) ^ set(failed_after)
        | {
            path
            for path in set(failed_before) & set(failed_after)
            if failed_before[path] != failed_after[path]
        }
    )
    fixed_generator = generator.loc[generator["Lane"].eq(FIXED_SCORE_LANE)]
    attempted_expected = int(plan["frozen_gates"]["attempted_replicates_must_equal"])
    wle_difference_gate = float(
        plan["frozen_gates"][
            "maximum_difference_from_original_wle_draws_for_same_seed"
        ]
    )
    contract_passed = bool(
        len(ledger) == attempted_expected
        and ledger["PersonFitAvailable"].all()
        and fixed_generator["ChangedPersonTotals"].eq(0).all()
        and wle_identity["IdentityMismatches"] == 0
        and wle_identity["MaximumAbsoluteNumericDifference"] <= wle_difference_gate
        and not changed_original_files
        and not changed_failed_files
    )
    ready_draws = person_draws.loc[
        person_draws["PersonFitReady"].fillna(False).astype(bool)
    ]
    infit_transitions = int(ready_draws["InfitZoneChanged"].astype(bool).sum())
    outfit_transitions = int(ready_draws["OutfitZoneChanged"].astype(bool).sum())
    any_transitions = int(ready_draws["AnyFitZoneTransition"].astype(bool).sum())
    raw_display = int(ready_draws["RawDisplayMismatch"].astype(bool).sum())
    display_boundaries = int(
        ready_draws["InfitBoundaryStatus"].eq("display_rounding_boundary").sum()
        + ready_draws["OutfitBoundaryStatus"].eq("display_rounding_boundary").sum()
    )
    numerical_boundaries = int(
        ready_draws["InfitBoundaryStatus"].eq("numerical_boundary").sum()
        + ready_draws["OutfitBoundaryStatus"].eq("numerical_boundary").sum()
    )
    decision = {
        "schema_version": "mfrm-cmle-wle-bootstrap-fit-extension-result-v1",
        "extension_executed": True,
        "contract_passed": contract_passed,
        "overall_status": "descriptive_fit_extension_complete_public_promotion_withheld"
        if contract_passed
        else "fit_extension_contract_failed",
        "attempted_replicates": int(len(ledger)),
        "person_replicates": int(len(person_draws)),
        "minimum_person_fit_available_share": float(
            model_summary["PersonFitAvailableShare"].min()
        ),
        "minimum_person_fit_available_wilson95_lower": float(
            model_summary["PersonFitAvailableWilson95Lower"].min()
        ),
        "infit_zone_transitions": infit_transitions,
        "outfit_zone_transitions": outfit_transitions,
        "any_fit_zone_transitions": any_transitions,
        "display_rounding_boundary_statistics": display_boundaries,
        "numerical_boundary_statistics": numerical_boundaries,
        "raw_display_decision_mismatches": raw_display,
        "maximum_absolute_wle_difference_from_original_pilot": float(
            wle_identity["MaximumAbsoluteNumericDifference"]
        ),
        "wle_identity_mismatches_from_original_pilot": int(
            wle_identity["IdentityMismatches"]
        ),
        "fixed_score_person_total_mismatches": int(
            fixed_generator["ChangedPersonTotals"].sum()
        ),
        "original_pilot_files_changed": changed_original_files,
        "failed_extension_files_changed": changed_failed_files,
        "wle_identity_comparison": "symmetric %.17g persisted CSV readback",
        "fit_statistics_use_unrounded_values": True,
        "outlier_trimming_applied": False,
        "fit_transition_evidence_descriptive_only": True,
        "threshold_optimality_validated": False,
        "coverage_qualified": False,
        "ZSTD_authorized": False,
        "p_value_authorized": False,
        "confidence_interval_authorized": False,
        "total_inferential_se_authorized": False,
        "streamlit_integration_authorized": False,
        "promotion_threshold_registered": False,
        "registered_extension_plan_sha256": sha256_file(args.plan),
        "registered_extension_amendment_sha256": sha256_file(args.amendment),
        "original_bootstrap_core_sha256": plan["input_identity"][
            "bootstrap_core_pre_extension_sha256"
        ],
        "extended_bootstrap_core_sha256": sha256_file(
            ROOT / "mfrm_app/cmle_wle_bootstrap.py"
        ),
        "person_fit_core_sha256": sha256_file(ROOT / "mfrm_app/cmle_wle_fit.py"),
        "runner_source_sha256": sha256_file(Path(__file__)),
    }
    (output / "bootstrap_fit_extension_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    report = f"""# CMLE-WLE bootstrap Person-fit extension

## Decision

**{decision['overall_status']}.** The prospectively registered fit extension used a new evidence directory and left the original pilot inventory unchanged. It completes descriptive untrimmed Person Infit/Outfit transition output; it does not authorize ZSTD, p-values, confidence intervals, a total SE, threshold optimality, or Streamlit integration.

## Attempt and identity accounting

- Attempted replicates: {len(ledger)}; Person-replicates: {len(person_draws)}.
- Minimum Person-fit availability: {model_summary['PersonFitAvailableShare'].min():.3f}; minimum Wilson 95% lower bound: {model_summary['PersonFitAvailableWilson95Lower'].min():.3f}.
- Same-seed WLE maximum difference from the original pilot: {wle_identity['MaximumAbsoluteNumericDifference']:.12g}; identity mismatches: {wle_identity['IdentityMismatches']}.
- Fixed-score Person-total mismatches: {int(fixed_generator['ChangedPersonTotals'].sum())}.
- Original pilot files changed during the extension: {len(changed_original_files)}.
- Failed first-run files changed during the corrected rerun: {len(changed_failed_files)}.
- WLE identity comparison used symmetric `%.17g` persisted CSV readback under the unchanged exact-zero gate.

## Fit and floating-point sensitivity

- Infit class transitions: {infit_transitions} of {len(ready_draws)} successful Person-replicates.
- Outfit class transitions: {outfit_transitions} of {len(ready_draws)} successful Person-replicates.
- At least one of Infit/Outfit changed class: {any_transitions} Person-replicates.
- Display-rounding boundary statistics: {display_boundaries}; numerical-boundary statistics: {numerical_boundaries}.
- Raw/display decision mismatches: {raw_display}. All transition counts used finite unrounded MnSq; displayed three-decimal values never drove classification.

## Critical boundary

These transitions are bootstrap-lane-specific sensitivity summaries. The fixed-score lane preserves every Person total; the joint plug-in lane does not. Their transition counts must not be pooled as one sampling distribution. TAM's default trimmed Outfit remains a different estimand from this untrimmed sirt-parity target.

ZSTD and p-values remain unavailable rather than zero. Repeated-truth coverage, sparse-connectivity and anchor-contamination stress, model misspecification/local dependence, and user-facing comprehension testing remain required before public integration.
"""
    (output / "CMLE_WLE_BOOTSTRAP_FIT_EXTENSION_RESULTS.md").write_text(
        report, encoding="utf-8"
    )
    if not contract_passed:
        raise SystemExit("CMLE-WLE bootstrap fit-extension gates failed.")


if __name__ == "__main__":
    main()
