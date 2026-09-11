#!/usr/bin/env python3
"""Run the frozen CMLE-WLE two-lane bootstrap implementation pilot."""

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
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.special import logsumexp


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle import fit_cmle  # noqa: E402
from mfrm_app.cmle_person_scoring import score_cmle_persons_wle  # noqa: E402
from mfrm_app.cmle_wle_bootstrap import (  # noqa: E402
    FIXED_SCORE_LANE,
    JOINT_PLUGIN_LANE,
    _row_category_kernels,
    _suffix_log_coefficients,
    generate_cmle_wle_bootstrap_sample,
    run_cmle_wle_bootstrap,
)


DEFAULT_PLAN = ROOT / "validation/cmle_wle_bootstrap_plan_20260810.json"
DEFAULT_OUTPUT = ROOT / "validation/cmle_wle_bootstrap_pilot_20260810"
REPRODUCIBILITY_REPLICATES = 10
JOINT_MONTE_CARLO_Z_TOLERANCE = 5.0


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_plan(plan_path: Path) -> dict[str, object]:
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != "mfrm-cmle-wle-bootstrap-plan-v1":
        raise ValueError("Unexpected CMLE-WLE bootstrap plan schema.")
    paths = {
        "cmle_core_sha256": ROOT / "mfrm_app/cmle.py",
        "fixed_calibration_person_scoring_sha256": ROOT / "mfrm_app/person_scoring.py",
        "cmle_wle_bridge_sha256": ROOT / "mfrm_app/cmle_person_scoring.py",
        "cmle_wle_uncertainty_sha256": ROOT / "mfrm_app/cmle_wle_uncertainty.py",
        "cmle_wle_bridge_decision_sha256": ROOT / "validation/cmle_wle_bridge_20260810/cmle_wle_bridge_decision.json",
        "wle_parity_decision_sha256": ROOT / "validation/fixed_calibration_wle_20260809/parity_decision.json",
        "calibration_sensitivity_decision_sha256": ROOT / "validation/cmle_wle_calibration_sensitivity_20260810/calibration_sensitivity_decision.json",
        "calibration_sensitivity_documentation_amendment_sha256": ROOT / "validation/cmle_wle_calibration_sensitivity_documentation_amendment_20260810.json",
        "cmle_phase0_design_sha256": ROOT / "docs/cmle_phase0_design.md",
        "bootstrap_design_sha256": ROOT / "docs/cmle_wle_bootstrap_design.md",
    }
    failed = [
        key
        for key, path in paths.items()
        if sha256_file(path) != str(plan["input_identity"][key])
    ]
    if failed:
        raise ValueError(f"CMLE-WLE bootstrap input identity failed: {failed}")
    if plan["implementation_pilot"].get("minimum_successful_share_for_promotion") is not None:
        raise ValueError("Implementation pilot must not contain a post-result promotion threshold.")
    return plan


def fixture() -> pd.DataFrame:
    scores = {
        "P01": [0, 1, 1, 2, 2, 0],
        "P02": [1, 2, 0, 1, 2, 1],
        "P03": [2, 0, 2, 1, 0, 1],
        "P04": [1, 0, 1, 2, 0, 2],
        "P05": [0, 2, 1, 0, 2, 1],
        "P06": [2, 1, 0, 2, 1, 0],
        "P07": [0, 1, 2, 0, 1, 2],
        "P08": [2, 2, 1, 0, 0, 1],
        "P09": [1, 1, 2, 2, 0, 0],
        "P10": [0, 2, 2, 1, 1, 0],
        "P11": [2, 0, 1, 1, 2, 0],
        "P12": [1, 2, 1, 0, 2, 0],
        "P13": [0, 0, 2, 1, 1, 2],
        "P14": [2, 1, 0, 0, 1, 2],
        "P15": [0, 0, 0, 0, 0, 0],
        "P16": [2, 2, 2, 2, 2, 2],
    }
    units = [
        ("R1", "C1"),
        ("R1", "C2"),
        ("R2", "C1"),
        ("R2", "C2"),
        ("R3", "C1"),
        ("R3", "C2"),
    ]
    omitted = {
        "P02": {5},
        "P03": {1},
        "P04": {2},
        "P05": {0, 5},
        "P06": {3},
        "P07": {4},
        "P08": {1, 4},
        "P09": {0},
        "P10": {2, 5},
        "P11": {3},
        "P12": {1},
        "P13": {4},
    }
    rows: list[tuple[str, str, str, int]] = []
    for person, values in scores.items():
        for index, ((rater, criterion), score) in enumerate(zip(units, values)):
            if index not in omitted.get(person, set()):
                rows.append((person, rater, criterion, score))
    return pd.DataFrame(rows, columns=["Person", "Rater", "Criterion", "Score"])


def fit_fixture(model: str) -> dict[str, object]:
    return fit_cmle(
        fixture(),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model=model,
        step_facet="Criterion" if model == "PCM" else None,
        gtol=1e-8,
    )


def recurrence_checks(tolerance: float) -> pd.DataFrame:
    cases = {
        "RSM_kernel": np.array(
            [[0.0, 0.2, -0.1], [0.0, -0.3, 0.5], [0.0, 0.7, -0.2]],
            dtype=float,
        ),
        "PCM_kernel": np.array(
            [[0.0, -0.5, 0.4], [0.0, 0.8, -0.1], [0.0, -0.2, 0.6]],
            dtype=float,
        ),
    }
    rows = []
    target = 3
    for case, kernels in cases.items():
        suffix = _suffix_log_coefficients(kernels, target)
        patterns = [
            values
            for values in product(range(kernels.shape[1]), repeat=kernels.shape[0])
            if sum(values) == target
        ]
        brute = logsumexp(
            [
                sum(kernels[row, value] for row, value in enumerate(values))
                for values in patterns
            ]
        )
        difference = abs(float(suffix[0, target]) - float(brute))
        rows.append(
            {
                "Case": case,
                "Target": target,
                "Patterns": len(patterns),
                "RecurrenceLogNormalizer": float(suffix[0, target]),
                "EnumeratedLogNormalizer": float(brute),
                "AbsoluteDifference": difference,
                "Tolerance": tolerance,
                "Passed": difference <= tolerance,
            }
        )
    return pd.DataFrame(rows)


def _child_seeds(seed: int, replicates: int) -> list[int]:
    return [
        int(child.generate_state(1, dtype=np.uint64)[0])
        for child in np.random.SeedSequence(seed).spawn(replicates)
    ]


def joint_frequency_check(
    fit: dict[str, object],
    *,
    replicates: int,
    seed: int,
) -> pd.DataFrame:
    design, kernels = _row_category_kernels(fit)
    wle = score_cmle_persons_wle(fit).set_index("Person")
    persons = design.data[design.person_col].astype(str).to_numpy(dtype=object)
    theta = pd.Series(persons).map(wle["Estimate"]).to_numpy(dtype=float)
    logits = kernels + theta[:, None] * np.arange(design.n_categories)[None, :]
    logits -= logsumexp(logits, axis=1, keepdims=True)
    probabilities = np.exp(logits)
    observed = np.zeros(design.n_categories, dtype=int)
    for child_seed in _child_seeds(seed, replicates):
        generated = generate_cmle_wle_bootstrap_sample(
            fit, lane=JOINT_PLUGIN_LANE, seed=child_seed
        )["data"]
        categories = (
            generated[design.score_col].to_numpy(dtype=int) - design.rating_min
        )
        observed += np.bincount(categories, minlength=design.n_categories)
    expected = replicates * probabilities.sum(axis=0)
    variance = replicates * np.sum(probabilities * (1.0 - probabilities), axis=0)
    standardized = np.divide(
        observed - expected,
        np.sqrt(variance),
        out=np.zeros_like(expected, dtype=float),
        where=variance > 0,
    )
    return pd.DataFrame(
        {
            "Model": design.model,
            "Category": np.arange(design.n_categories) + design.rating_min,
            "ObservedCount": observed,
            "ExpectedCount": expected,
            "MonteCarloSD": np.sqrt(variance),
            "StandardizedResidual": standardized,
            "AbsoluteStandardizedResidual": np.abs(standardized),
            "Tolerance": JOINT_MONTE_CARLO_Z_TOLERANCE,
            "Passed": np.abs(standardized) <= JOINT_MONTE_CARLO_Z_TOLERANCE,
        }
    )


def reproducibility_check(
    primary: dict[str, object],
    repeated: dict[str, object],
) -> dict[str, object]:
    count = int(repeated["summary"].iloc[0]["AttemptedReplicates"])
    ledger_columns = [
        "Lane",
        "Model",
        "Replicate",
        "SeedIdentity",
        "CMLEConverged",
        "CMLEInferenceReady",
        "Rank",
        "Nullity",
        "WLEAvailable",
        "FailureStage",
        "FailureReason",
    ]
    left_ledger = primary["ledger"].iloc[:count][ledger_columns].reset_index(drop=True)
    right_ledger = repeated["ledger"][ledger_columns].reset_index(drop=True)
    status_mismatches = int(
        np.sum(left_ledger.astype(str).to_numpy() != right_ledger.astype(str).to_numpy())
    )
    keys = ["Lane", "Model", "Replicate", "Person"]
    numeric = [
        "Estimate",
        "StandardError",
        "AdjustedScoreResidual",
        "BaselineEstimate",
        "BootstrapMinusBaseline",
    ]
    left_persons = primary["person_draws"].loc[
        primary["person_draws"]["Replicate"].le(count), keys + numeric
    ]
    right_persons = repeated["person_draws"][keys + numeric]
    merged = left_persons.merge(
        right_persons,
        on=keys,
        how="outer",
        suffixes=("Primary", "Repeated"),
        indicator=True,
    )
    identity_mismatches = int(merged["_merge"].ne("both").sum())
    differences = []
    for column in numeric:
        differences.extend(
            np.abs(
                merged[f"{column}Primary"].to_numpy(dtype=float)
                - merged[f"{column}Repeated"].to_numpy(dtype=float)
            ).tolist()
        )
    finite_differences = np.asarray(differences, dtype=float)
    finite_differences = finite_differences[np.isfinite(finite_differences)]
    maximum_difference = (
        float(np.max(finite_differences)) if len(finite_differences) else 0.0
    )
    return {
        "Model": str(primary["summary"].iloc[0]["Model"]),
        "Lane": str(primary["summary"].iloc[0]["Lane"]),
        "ReplayedReplicates": count,
        "LedgerCellMismatches": status_mismatches,
        "PersonIdentityMismatches": identity_mismatches,
        "MaximumAbsNumericDifference": maximum_difference,
    }


def wilson_interval(successes: int, attempts: int) -> tuple[float, float]:
    if attempts <= 0:
        return np.nan, np.nan
    z = 1.959963984540054
    proportion = successes / attempts
    denominator = 1.0 + z * z / attempts
    center = (proportion + z * z / (2.0 * attempts)) / denominator
    half = z * np.sqrt(
        proportion * (1.0 - proportion) / attempts + z * z / (4.0 * attempts**2)
    ) / denominator
    return float(center - half), float(center + half)


def summarize_persons(
    results: list[dict[str, object]],
    replicates: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for result in results:
        model = str(result["summary"].iloc[0]["Model"])
        lane = str(result["summary"].iloc[0]["Lane"])
        baseline = result["baseline_persons"].copy()
        baseline["Person"] = baseline["Person"].astype(str)
        draws = result["person_draws"].copy()
        if not draws.empty:
            draws["Person"] = draws["Person"].astype(str)
        for baseline_row in baseline.itertuples(index=False):
            person = str(baseline_row.Person)
            frame = draws.loc[draws["Person"].eq(person)] if not draws.empty else draws
            estimates = (
                frame["Estimate"].to_numpy(dtype=float)
                if not frame.empty else np.asarray([], dtype=float)
            )
            successful = int(np.isfinite(estimates).sum())
            values = estimates[np.isfinite(estimates)]
            baseline_estimate = float(baseline_row.Estimate)
            baseline_extreme = bool(baseline_row.ExtremeScorePattern)
            extreme_transitions = (
                int(np.sum(frame["ExtremeScorePattern"].astype(bool) != baseline_extreme))
                if not frame.empty else 0
            )
            sign_transitions = (
                int(np.sum(np.sign(values) != np.sign(baseline_estimate)))
                if len(values) else 0
            )
            rows.append(
                {
                    "Model": model,
                    "Lane": lane,
                    "Person": person,
                    "BaselineEstimate": baseline_estimate,
                    "BaselineConditionalWLESE": float(baseline_row.StandardError),
                    "BaselineExtremeStatus": baseline_extreme,
                    "Attempts": replicates,
                    "SuccessfulScores": successful,
                    "BootstrapMean": float(np.mean(values)) if len(values) else np.nan,
                    "BootstrapSD": (
                        float(np.std(values, ddof=1)) if len(values) >= 2 else np.nan
                    ),
                    "BootstrapQ025": (
                        float(np.quantile(values, 0.025)) if len(values) else np.nan
                    ),
                    "BootstrapQ500": (
                        float(np.quantile(values, 0.5)) if len(values) else np.nan
                    ),
                    "BootstrapQ975": (
                        float(np.quantile(values, 0.975)) if len(values) else np.nan
                    ),
                    "BootstrapMeanMinusBaseline": (
                        float(np.mean(values) - baseline_estimate) if len(values) else np.nan
                    ),
                    "BootstrapSDToBaselineConditionalWLESERatio": (
                        float(np.std(values, ddof=1) / float(baseline_row.StandardError))
                        if len(values) >= 2 and float(baseline_row.StandardError) > 0
                        else np.nan
                    ),
                    "ExtremeStateTransitions": extreme_transitions,
                    "SignTransitions": sign_transitions,
                    "BaselineNearZeroForSign": abs(baseline_estimate) < 0.05,
                    "SignTransitionsDescriptiveOnly": True,
                    "FitZoneTransitions": np.nan,
                    "FitZoneTransitionsAvailable": False,
                    "RawDisplayMismatches": np.nan,
                    "RawDisplayMismatchesAvailable": False,
                    "BootstrapQuantilesAreConfidenceInterval": False,
                    "CoverageQualified": False,
                }
            )
    return pd.DataFrame(rows)


def summarize_extreme_transitions(
    results: list[dict[str, object]],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for result in results:
        model = str(result["summary"].iloc[0]["Model"])
        lane = str(result["summary"].iloc[0]["Lane"])
        baseline = (
            result["baseline_persons"]
            .assign(Person=lambda frame: frame["Person"].astype(str))
            .set_index("Person")["ExtremeScorePattern"]
            .astype(bool)
        )
        draws = result["person_draws"].copy()
        draws["Person"] = draws["Person"].astype(str)
        draws["BaselineExtreme"] = draws["Person"].map(baseline).astype(bool)
        draws["CurrentExtreme"] = draws["ExtremeScorePattern"].astype(bool)
        transition = np.select(
            [
                draws["BaselineExtreme"] & ~draws["CurrentExtreme"],
                ~draws["BaselineExtreme"] & draws["CurrentExtreme"],
                draws["BaselineExtreme"] & draws["CurrentExtreme"],
            ],
            ["extreme_to_interior", "interior_to_extreme", "unchanged_extreme"],
            default="unchanged_interior",
        )
        counts = pd.Series(transition).value_counts()
        for state in (
            "extreme_to_interior",
            "interior_to_extreme",
            "unchanged_extreme",
            "unchanged_interior",
        ):
            rows.append(
                {
                    "Model": model,
                    "Lane": lane,
                    "Transition": state,
                    "PersonReplicates": int(counts.get(state, 0)),
                }
            )
    return pd.DataFrame(rows)


def first_read_projection(
    *,
    sampler_passed: bool,
    model_summary: pd.DataFrame,
    person_summary: pd.DataFrame,
) -> pd.DataFrame:
    minimum_ready = float(model_summary["InferenceReadyShare"].min())
    maximum_sd = float(person_summary["BootstrapSD"].max())
    transitions = int(person_summary["ExtremeStateTransitions"].sum())
    return pd.DataFrame(
        [
            {
                "CardOrder": 1,
                "CardId": "bootstrap_generator",
                "Status": "research_ready" if sampler_passed else "blocked",
                "DisplayValue": "sampler contracts passed" if sampler_passed else "sampler blocked",
                "Interpretation": "Enumeration, total preservation, seed, and frequency gates.",
                "NextAction": "Inspect lane-specific generator evidence.",
            },
            {
                "CardOrder": 2,
                "CardId": "cmle_refit_readiness",
                "Status": "research_ready" if minimum_ready == 1.0 else "caution",
                "DisplayValue": f"minimum readiness {100.0 * minimum_ready:.1f}%",
                "Interpretation": "Readiness uses all attempted replicates as denominator.",
                "NextAction": "Inspect rank/nullity and failure-stage tables.",
            },
            {
                "CardOrder": 3,
                "CardId": "person_sensitivity",
                "Status": "caution",
                "DisplayValue": f"maximum bootstrap SD {maximum_sd:.3f} logits",
                "Interpretation": "The two lanes answer different sensitivity questions.",
                "NextAction": "Inspect Person-level distributions without pooling lanes.",
            },
            {
                "CardOrder": 4,
                "CardId": "exact_extremes",
                "Status": "caution",
                "DisplayValue": f"{transitions} extreme-state transitions",
                "Interpretation": "Fixed-score and joint plug-in extreme behavior are not equivalent.",
                "NextAction": "Show baseline and replicate extreme status together.",
            },
            {
                "CardOrder": 5,
                "CardId": "fit_thresholds",
                "Status": "withheld",
                "DisplayValue": "CMLE-WLE fit-zone transitions not qualified",
                "Interpretation": "No unvalidated fit statistic was synthesized for this pilot.",
                "NextAction": "Prospectively validate CMLE-WLE residual fit before thresholds.",
            },
            {
                "CardOrder": 6,
                "CardId": "interval_claim",
                "Status": "withheld",
                "DisplayValue": "confidence interval withheld",
                "Interpretation": "No repeated-truth coverage evidence exists.",
                "NextAction": "Register a post-pilot ADEMP coverage study.",
            },
            {
                "CardOrder": 7,
                "CardId": "public_ui",
                "Status": "withheld",
                "DisplayValue": "Streamlit integration not authorized",
                "Interpretation": "A completed implementation pilot is not public validation.",
                "NextAction": "Complete output, coverage, new-Person, anchor, and UI gates.",
            },
        ]
    )


def plot_results(
    model_summary: pd.DataFrame,
    person_summary: pd.DataFrame,
    extreme_transition_summary: pd.DataFrame,
    output: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    lane_short = {
        FIXED_SCORE_LANE: "fixed-score",
        JOINT_PLUGIN_LANE: "joint plug-in",
    }
    labels = [f"{row.Model}\n{lane_short[row.Lane]}" for row in model_summary.itertuples()]
    x = np.arange(len(labels))
    axes[0].bar(x - 0.18, model_summary["InferenceReadyShare"], width=0.36, label="CMLE ready")
    axes[0].bar(x + 0.18, model_summary["WLEAvailableShare"], width=0.36, label="WLE available")
    axes[0].set_ylim(0, 1.05)
    axes[0].set_ylabel("Share of all attempted replicates")
    axes[0].set_xticks(x, labels, rotation=20, ha="right", fontsize=8)
    axes[0].set_title(
        "Bootstrap refit accounting\n"
        f"200/200 each; minimum Wilson 95% lower {model_summary['InferenceReadyWilson95Lower'].min():.3f}"
    )
    axes[0].legend(frameon=False)

    lane_colors = {FIXED_SCORE_LANE: "#1f77b4", JOINT_PLUGIN_LANE: "#d95f02"}
    for (model, lane), frame in person_summary.groupby(["Model", "Lane"], sort=False):
        extreme = frame["BaselineExtremeStatus"].astype(bool)
        axes[1].scatter(
            frame.loc[~extreme, "BaselineEstimate"],
            frame.loc[~extreme, "BootstrapSD"],
            s=28,
            alpha=0.7,
            color=lane_colors[lane],
            marker="o" if model == "RSM" else "s",
        )
        axes[1].scatter(
            frame.loc[extreme, "BaselineEstimate"],
            frame.loc[extreme, "BootstrapSD"],
            s=75,
            color=lane_colors[lane],
            marker="X",
        )
    axes[1].set_xlabel("Baseline fixed-calibration WLE estimate")
    axes[1].set_ylabel("Bootstrap SD (successful refits)")
    axes[1].set_title("Person sensitivity; lanes not pooled")
    axes[1].legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=lane_colors[FIXED_SCORE_LANE], label="Fixed-score lane"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=lane_colors[JOINT_PLUGIN_LANE], label="Joint plug-in lane"),
            Line2D([0], [0], marker="o", color="gray", linestyle="none", label="RSM interior"),
            Line2D([0], [0], marker="s", color="gray", linestyle="none", label="PCM interior"),
            Line2D([0], [0], marker="X", color="gray", linestyle="none", label="Baseline exact extreme"),
        ],
        frameon=False,
        fontsize=7,
        ncol=2,
    )
    fig.tight_layout()
    fig.savefig(output / "bootstrap_readiness_and_person_sd.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    conditions = extreme_transition_summary[["Model", "Lane"]].drop_duplicates()
    conditions["Label"] = [
        f"{row.Model}\n{lane_short[row.Lane]}" for row in conditions.itertuples()
    ]
    directions = (
        extreme_transition_summary.loc[
            extreme_transition_summary["Transition"].isin(
                ["interior_to_extreme", "extreme_to_interior"]
            )
        ]
        .pivot(index=["Model", "Lane"], columns="Transition", values="PersonReplicates")
        .fillna(0)
        .reindex(pd.MultiIndex.from_frame(conditions[["Model", "Lane"]]))
    )
    x = np.arange(len(conditions))
    first = ax.bar(
        x,
        directions["interior_to_extreme"],
        color="#fdae6b",
        label="interior → extreme",
    )
    second = ax.bar(
        x,
        directions["extreme_to_interior"],
        bottom=directions["interior_to_extreme"],
        color="#d95f02",
        label="extreme → interior",
    )
    ax.bar_label(
        first,
        labels=[f"{int(value)}" if value > 0 else "" for value in directions["interior_to_extreme"]],
        label_type="center",
        fontsize=8,
    )
    ax.bar_label(
        second,
        labels=[f"{int(value)}" if value > 0 else "" for value in directions["extreme_to_interior"]],
        label_type="center",
        fontsize=8,
    )
    totals = directions.sum(axis=1).to_numpy(dtype=float)
    for position, total in zip(x, totals):
        if total > 0:
            ax.text(position, total + max(3.0, totals.max() * 0.015), f"{int(total)}", ha="center", va="bottom")
    ax.set_xticks(x, conditions["Label"], fontsize=8)
    ax.set_ylabel("Person-replicate state transitions")
    ax.set_title("Direction of exact-extreme status transitions")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "bootstrap_extreme_state_transitions.png", dpi=180)
    plt.close(fig)


def run(plan_path: Path, output: Path) -> None:
    plan_path = plan_path.resolve()
    output = output.resolve()
    plan = validate_plan(plan_path)
    matrix = plan["implementation_pilot"]
    models = [str(value) for value in matrix["models"]]
    replicates = int(matrix["replicates_per_model_lane"])
    seed = int(matrix["primary_seed"])
    alternate_seed = int(matrix["alternate_seed"])
    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_bootstrap_plan.json").write_bytes(plan_path.read_bytes())
    fixture().to_csv(output / "bootstrap_fixture.csv", index=False)

    recurrence = recurrence_checks(
        float(matrix["conditional_probability_enumeration_tolerance"])
    )
    results: list[dict[str, object]] = []
    baseline_rows: list[dict[str, object]] = []
    reproducibility_rows: list[dict[str, object]] = []
    different_seed_rows: list[dict[str, object]] = []
    joint_frequency_parts: list[pd.DataFrame] = []
    for model in models:
        fit = fit_fixture(model)
        fit_row = fit["summary"].iloc[0]
        baseline_rows.append(
            {
                "Model": model,
                "Rows": len(fit["design"].data),
                "Persons": int(fit_row["PersonsTotal"]),
                "ExactExtremes": int(fit_row["PersonsExtreme"]),
                "MissingnessPatterns": int(fit_row["Patterns"]),
                "Converged": bool(fit_row["Converged"]),
                "InferenceReady": bool(fit_row["InferenceReady"]),
                "InformationRank": int(fit_row["InformationRank"]),
                "InformationNullity": int(fit_row["InformationNullity"]),
                "GradientSupNorm": float(fit_row["GradientSupNorm"]),
            }
        )
        joint_frequency_parts.append(
            joint_frequency_check(fit, replicates=replicates, seed=seed)
        )
        for lane in (FIXED_SCORE_LANE, JOINT_PLUGIN_LANE):
            result = run_cmle_wle_bootstrap(
                fit,
                lane=lane,
                n_replicates=replicates,
                seed=seed,
                gtol=1e-8,
            )
            results.append(result)
            repeated = run_cmle_wle_bootstrap(
                fit,
                lane=lane,
                n_replicates=REPRODUCIBILITY_REPLICATES,
                seed=seed,
                gtol=1e-8,
            )
            reproducibility_rows.append(reproducibility_check(result, repeated))
            first_seed = generate_cmle_wle_bootstrap_sample(
                fit, lane=lane, seed=seed
            )["data"][fit["design"].score_col].to_numpy(dtype=int)
            second_seed = generate_cmle_wle_bootstrap_sample(
                fit, lane=lane, seed=alternate_seed
            )["data"][fit["design"].score_col].to_numpy(dtype=int)
            different_seed_rows.append(
                {
                    "Model": model,
                    "Lane": lane,
                    "DifferentSeedChangedResponses": int(np.sum(first_seed != second_seed)),
                    "Passed": bool(np.any(first_seed != second_seed)),
                }
            )

    baseline = pd.DataFrame(baseline_rows)
    model_summary = pd.concat([result["summary"] for result in results], ignore_index=True)
    ledger = pd.concat([result["ledger"] for result in results], ignore_index=True)
    generator = pd.concat([result["generator_audit"] for result in results], ignore_index=True)
    person_draws = pd.concat(
        [result["person_draws"] for result in results], ignore_index=True
    )
    person_summary = summarize_persons(results, replicates)
    extreme_transition_summary = summarize_extreme_transitions(results)
    reproducibility = pd.DataFrame(reproducibility_rows)
    different_seed = pd.DataFrame(different_seed_rows)
    joint_frequency = pd.concat(joint_frequency_parts, ignore_index=True)

    for index, row in model_summary.iterrows():
        attempts = int(row["AttemptedReplicates"])
        ready = int(row["CMLEInferenceReadyReplicates"])
        lower, upper = wilson_interval(ready, attempts)
        model_summary.loc[index, "InferenceReadyWilson95Lower"] = lower
        model_summary.loc[index, "InferenceReadyWilson95Upper"] = upper
        available = int(row["WLEAvailableReplicates"])
        lower, upper = wilson_interval(available, attempts)
        model_summary.loc[index, "WLEAvailableWilson95Lower"] = lower
        model_summary.loc[index, "WLEAvailableWilson95Upper"] = upper

    failure_summary = ledger.assign(
        FailureStage=ledger["FailureStage"].replace("", "none")
    ).groupby(["Model", "Lane", "FailureStage"], as_index=False).size()
    rank_summary = (
        ledger.groupby(["Model", "Lane", "Rank", "Nullity"], dropna=False, as_index=False)
        .size()
    )

    fixed_generator = generator.loc[generator["Lane"].eq(FIXED_SCORE_LANE)]
    same_seed_passed = bool(
        reproducibility["LedgerCellMismatches"].eq(0).all()
        and reproducibility["PersonIdentityMismatches"].eq(0).all()
        and reproducibility["MaximumAbsNumericDifference"].le(
            float(matrix["same_seed_max_abs_numeric_difference"])
        ).all()
    )
    sampler_passed = bool(
        baseline["InferenceReady"].all()
        and recurrence["Passed"].all()
        and len(ledger) == len(models) * 2 * replicates
        and len(generator) == len(models) * 2 * replicates
        and fixed_generator["ChangedPersonTotals"].eq(
            int(matrix["fixed_score_total_mismatches_allowed"])
        ).all()
        and ledger["Rank"].notna().all()
        and ledger["Nullity"].notna().all()
        and same_seed_passed
        and different_seed["Passed"].all()
        and joint_frequency["Passed"].all()
    )
    output_contract_complete = bool(
        person_summary["FitZoneTransitionsAvailable"].all()
        and person_summary["RawDisplayMismatchesAvailable"].all()
    )
    first_read = first_read_projection(
        sampler_passed=sampler_passed,
        model_summary=model_summary,
        person_summary=person_summary,
    )

    artifacts = {
        "bootstrap_baseline_fits.csv": baseline,
        "bootstrap_replicate_ledger.csv": ledger,
        "bootstrap_generator_audit.csv": generator,
        "bootstrap_model_summary.csv": model_summary,
        "bootstrap_failure_summary.csv": failure_summary,
        "bootstrap_rank_summary.csv": rank_summary,
        "bootstrap_person_draws.csv": person_draws,
        "bootstrap_person_summary.csv": person_summary,
        "bootstrap_extreme_transition_summary.csv": extreme_transition_summary,
        "bootstrap_recurrence_check.csv": recurrence,
        "bootstrap_joint_frequency_check.csv": joint_frequency,
        "bootstrap_reproducibility.csv": reproducibility,
        "bootstrap_different_seed_check.csv": different_seed,
        "cmle_wle_bootstrap_first_read_projection.csv": first_read,
    }
    for filename, frame in artifacts.items():
        frame.to_csv(output / filename, index=False, float_format="%.17g")
    plot_results(model_summary, person_summary, extreme_transition_summary, output)

    maximum_sd = float(person_summary["BootstrapSD"].max())
    median_by_lane = (
        person_summary.groupby("Lane")["BootstrapSD"].median().to_dict()
    )
    median_ratio_by_lane = (
        person_summary.groupby("Lane")["BootstrapSDToBaselineConditionalWLESERatio"]
        .median()
        .to_dict()
    )
    total_transitions = int(person_summary["ExtremeStateTransitions"].sum())
    extreme_to_interior = int(
        extreme_transition_summary.loc[
            extreme_transition_summary["Transition"].eq("extreme_to_interior"),
            "PersonReplicates",
        ].sum()
    )
    interior_to_extreme = int(
        extreme_transition_summary.loc[
            extreme_transition_summary["Transition"].eq("interior_to_extreme"),
            "PersonReplicates",
        ].sum()
    )
    decision = {
        "schema_version": "mfrm-cmle-wle-bootstrap-pilot-result-v1",
        "pilot_executed": True,
        "sampler_contract_passed": sampler_passed,
        "registered_output_contract_complete": output_contract_complete,
        "overall_status": (
            "pilot_complete_promotion_withheld"
            if sampler_passed else "pilot_complete_sampler_contract_failed"
        ),
        "models": models,
        "lanes": [FIXED_SCORE_LANE, JOINT_PLUGIN_LANE],
        "replicates_per_model_lane": replicates,
        "attempted_replicates": len(ledger),
        "minimum_inference_ready_share": float(model_summary["InferenceReadyShare"].min()),
        "minimum_wle_available_share": float(model_summary["WLEAvailableShare"].min()),
        "maximum_person_bootstrap_sd": maximum_sd,
        "median_person_bootstrap_sd_by_lane": median_by_lane,
        "median_bootstrap_sd_to_baseline_conditional_wle_se_ratio_by_lane": median_ratio_by_lane,
        "maximum_abs_bootstrap_mean_minus_baseline": float(
            person_summary["BootstrapMeanMinusBaseline"].abs().max()
        ),
        "extreme_state_transitions": total_transitions,
        "extreme_to_interior_transitions": extreme_to_interior,
        "interior_to_extreme_transitions": interior_to_extreme,
        "fit_zone_transitions_evaluated": False,
        "raw_display_mismatches_evaluated": False,
        "coverage_qualified": False,
        "confidence_interval_authorized": False,
        "total_inferential_se_authorized": False,
        "streamlit_integration_authorized": False,
        "promotion_threshold_registered": False,
        "bootstrap_core_sha256": sha256_file(ROOT / "mfrm_app/cmle_wle_bootstrap.py"),
        "runner_source_sha256": sha256_file(Path(__file__).resolve()),
    }
    (output / "bootstrap_pilot_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    report = f"""# CMLE-WLE two-lane bootstrap implementation pilot

## Decision

**{decision['overall_status']}.** The frozen 200-replicate RSM/PCM implementation matrix was executed without adding a post-result success threshold. Sampler-contract status is **{'PASS' if sampler_passed else 'FAIL'}**. Promotion remains withheld because fit-zone/raw-display outputs and repeated-truth interval coverage are not qualified.

## Attempt accounting

- Attempted replicates: {len(ledger)} ({replicates} per model/lane).
- Minimum CMLE inference-ready share: {model_summary['InferenceReadyShare'].min():.3f}.
- Minimum WLE-available share: {model_summary['WLEAvailableShare'].min():.3f}.
- Fixed-score Person-total mismatches: {int(fixed_generator['ChangedPersonTotals'].sum())}.
- Maximum joint-generator absolute standardized category residual: {joint_frequency['AbsoluteStandardizedResidual'].max():.3f} (frozen tolerance {JOINT_MONTE_CARLO_Z_TOLERANCE:.1f}).
- Same-seed maximum numeric difference: {reproducibility['MaximumAbsNumericDifference'].max():.3g}; status/identity mismatches: {int(reproducibility['LedgerCellMismatches'].sum() + reproducibility['PersonIdentityMismatches'].sum())}.

## Person sensitivity

- Median fixed-score conditional-pattern bootstrap SD: {median_by_lane[FIXED_SCORE_LANE]:.6g} logits.
- Median joint plug-in bootstrap SD: {median_by_lane[JOINT_PLUGIN_LANE]:.6g} logits.
- Maximum Person bootstrap SD: {maximum_sd:.6g} logits.
- Median bootstrap-SD / baseline-conditional-WLE-SE ratio: {median_ratio_by_lane[FIXED_SCORE_LANE]:.6g} fixed-score and {median_ratio_by_lane[JOINT_PLUGIN_LANE]:.6g} joint plug-in.
- Maximum absolute bootstrap-mean minus baseline displacement: {person_summary['BootstrapMeanMinusBaseline'].abs().max():.6g} logits.
- Extreme-state transitions: {total_transitions} Person-replicates ({extreme_to_interior} extreme to interior; {interior_to_extreme} interior to extreme).

These lane dispersions are not interchangeable and their difference is not a variance decomposition. Percentile columns are not confidence intervals.

## Withheld outputs

This pilot does not synthesize an unvalidated CMLE-WLE Person-fit statistic. Fit-zone transitions and raw/display mismatches are exported as unavailable, not zero. Repeated-truth coverage, new-Person/unseen-unit scoring, anchors, total inferential SE, and Streamlit integration remain withheld.
"""
    (output / "CMLE_WLE_BOOTSTRAP_PILOT_RESULTS.md").write_text(
        report, encoding="utf-8"
    )
    if not sampler_passed:
        raise SystemExit("CMLE-WLE bootstrap sampler contract failed.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.plan, args.output)


if __name__ == "__main__":
    main()
