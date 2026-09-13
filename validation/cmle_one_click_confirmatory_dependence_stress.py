#!/usr/bin/env python3
"""Run the registered CMLE confirmatory dependence/MNAR stress surface."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

_MPL_CONFIG = Path(tempfile.gettempdir()) / "mfrm_app_matplotlib"
_MPL_CONFIG.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CONFIG))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle_one_click_confirmatory_dependence import (  # noqa: E402
    DOMAINS,
    LANGUAGES,
    dependence_stress_human_gate_status,
    derive_condition_seed,
    simulate_dependence_condition,
    simulate_registered_dependence_surface,
)


PLAN = ROOT / "validation/cmle_one_click_confirmatory_dependence_stress_plan_20260810.json"
OUTPUT = ROOT / "validation/cmle_one_click_confirmatory_dependence_stress_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dataframe_sha256(frame: pd.DataFrame) -> str:
    buffer = io.StringIO()
    frame.to_csv(buffer, index=False, float_format="%.17g")
    return hashlib.sha256(buffer.getvalue().encode("utf-8")).hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, float_format="%.17g")


def validate_registration() -> dict[str, object]:
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    if plan.get("study_id") != "cmle_one_click_confirmatory_dependence_stress_v1":
        raise ValueError("Unexpected dependence stress study identity.")
    mismatches = [
        relative
        for relative, expected in plan["parent_identity"].items()
        if not (ROOT / relative).is_file()
        or sha256_file(ROOT / relative) != expected
    ]
    if mismatches:
        raise ValueError(f"Dependence stress parent identity failed: {mismatches}")
    fixed = plan["fixed_design"]
    if tuple(fixed["languages"]) != LANGUAGES:
        raise ValueError("Registered language order differs from implementation.")
    if tuple(fixed["domains"]) != DOMAINS:
        raise ValueError("Registered domain order differs from implementation.")
    if fixed["candidate_planned_slots_per_language"] != [100, 200, 300, 500]:
        raise ValueError("Unexpected planned-slot surface.")
    if int(fixed["monte_carlo_replicates_per_scenario_n"]) != 2000:
        raise ValueError("Unexpected Monte Carlo replicate count.")
    if [item["scenario"] for item in plan["scenarios"]] != [
        "independent_safe",
        "participant_dependence_safe",
        "moderator_cluster_safe",
        "pooled_domain_hotspot",
        "blocked_mechanism_hotspot",
        "mnar_danger_underretained",
        "mnar_danger_overretained",
        "combined_adversarial",
    ]:
        raise ValueError("Registered scenario order changed.")
    return plan


def build_registration_audit(plan: dict[str, object]) -> pd.DataFrame:
    rows = []
    for relative, expected in plan["parent_identity"].items():
        path = ROOT / relative
        actual = sha256_file(path) if path.is_file() else ""
        rows.append(
            {
                "Artifact": relative,
                "ExpectedSHA256": expected,
                "ActualSHA256": actual,
                "Passed": bool(actual == expected),
            }
        )
    return pd.DataFrame(rows)


def build_seed_audit(
    summary: pd.DataFrame, *, base_seed: int, scenarios: list[dict[str, object]]
) -> pd.DataFrame:
    rows = []
    order = {item["scenario"]: index for index, item in enumerate(scenarios)}
    for _, row in summary.iterrows():
        expected = derive_condition_seed(
            base_seed, order[str(row["Scenario"])], int(row["PlannedSlotsPerLanguage"])
        )
        rows.append(
            {
                "Scenario": row["Scenario"],
                "PlannedSlotsPerLanguage": row["PlannedSlotsPerLanguage"],
                "ExpectedSeed": expected,
                "ActualSeed": int(row["Seed"]),
                "Passed": bool(expected == int(row["Seed"])),
            }
        )
    return pd.DataFrame(rows)


def build_marginal_calibration_audit(cells: pd.DataFrame) -> pd.DataFrame:
    result = cells[
        [
            "Scenario",
            "PlannedSlotsPerLanguage",
            "Language",
            "Domain",
            "TruthDangerousProbabilityRaw",
            "MeanCompleteDangerousProportionRaw",
            "CompleteDangerousProportionMonteCarloSERaw",
        ]
    ].copy()
    result["AbsoluteDifferenceRaw"] = (
        result["MeanCompleteDangerousProportionRaw"]
        - result["TruthDangerousProbabilityRaw"]
    ).abs()
    result["RegisteredToleranceRaw"] = (
        6.0 * result["CompleteDangerousProportionMonteCarloSERaw"] + 0.002
    )
    result["Passed"] = result["AbsoluteDifferenceRaw"].le(
        result["RegisteredToleranceRaw"]
    )
    return result


def build_baseline_binomial_audit(
    cells: pd.DataFrame, replicates: int
) -> pd.DataFrame:
    result = cells.loc[
        cells["Scenario"].eq("independent_safe"),
        [
            "Scenario",
            "PlannedSlotsPerLanguage",
            "Language",
            "Domain",
            "PrimaryPassRateRaw",
            "PrimaryPassMonteCarloSERaw",
            "ExactBinomialPrimaryPassProbabilityRaw",
        ],
    ].copy()
    exact = result["ExactBinomialPrimaryPassProbabilityRaw"]
    result["AbsoluteDifferenceRaw"] = (
        result["PrimaryPassRateRaw"] - exact
    ).abs()
    result["RegisteredToleranceRaw"] = (
        6.0
        * np.sqrt(
            result["PrimaryPassMonteCarloSERaw"] ** 2
            + exact * (1.0 - exact) / replicates
        )
        + 0.01
    )
    result["Passed"] = exact.notna() & result["AbsoluteDifferenceRaw"].le(
        result["RegisteredToleranceRaw"]
    )
    return result


def build_retention_calibration_audit(cells: pd.DataFrame) -> pd.DataFrame:
    result = cells[
        [
            "Scenario",
            "PlannedSlotsPerLanguage",
            "Language",
            "Domain",
            "ExpectedValidFractionRaw",
            "MeanValidFractionRaw",
            "ValidFractionMonteCarloSERaw",
            "ExpectedObservedDangerousProbabilityRaw",
            "MeanObservedDangerousProportionRaw",
            "ObservedDangerousProportionMonteCarloSERaw",
        ]
    ].copy()
    result["ValidFractionAbsoluteDifferenceRaw"] = (
        result["MeanValidFractionRaw"] - result["ExpectedValidFractionRaw"]
    ).abs()
    result["ValidFractionToleranceRaw"] = (
        6.0 * result["ValidFractionMonteCarloSERaw"] + 0.002
    )
    result["ObservedDangerAbsoluteDifferenceRaw"] = (
        result["MeanObservedDangerousProportionRaw"]
        - result["ExpectedObservedDangerousProbabilityRaw"]
    ).abs()
    result["ObservedDangerToleranceRaw"] = (
        6.0 * result["ObservedDangerousProportionMonteCarloSERaw"] + 0.002
    )
    result["Passed"] = result["ValidFractionAbsoluteDifferenceRaw"].le(
        result["ValidFractionToleranceRaw"]
    ) & result["ObservedDangerAbsoluteDifferenceRaw"].le(
        result["ObservedDangerToleranceRaw"]
    )
    return result


def build_logic_audit(
    summary: pd.DataFrame, replicate: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    grouped = replicate.groupby(
        ["Scenario", "PlannedSlotsPerLanguage"], sort=False
    )
    for (scenario, planned_n), frame in grouped:
        source = summary.loc[
            summary["Scenario"].eq(scenario)
            & summary["PlannedSlotsPerLanguage"].eq(planned_n)
        ].iloc[0]
        familywise_violations = int(
            (frame["All12FamilywisePass"] & ~frame["All12PrimaryPass"]).sum()
        )
        pooled_flag_violations = (
            0
            if bool(source["AnyPooledDomainTruthAtOrAboveThreshold"])
            else int(frame["FalseReassuringPooledTruth"].sum())
        )
        mechanism_flag_violations = (
            0
            if bool(source["AnyBlockedMechanismTruthAtOrAboveThreshold"])
            else int(frame["FalseReassuringMechanismTruth"].sum())
        )
        summary_mismatches = int(
            int(frame["All12PrimaryPass"].sum())
            != int(source["All12PrimaryPassSuccesses"])
        ) + int(
            int(frame["All12FamilywisePass"].sum())
            != int(source["All12FamilywisePassSuccesses"])
        )
        rows.append(
            {
                "Scenario": scenario,
                "PlannedSlotsPerLanguage": planned_n,
                "ReplicateRows": len(frame),
                "FamilywiseWithoutPrimaryViolations": familywise_violations,
                "FalsePooledWithoutTruthRiskViolations": pooled_flag_violations,
                "FalseMechanismWithoutTruthRiskViolations": mechanism_flag_violations,
                "SummaryCountMismatches": summary_mismatches,
                "Passed": bool(
                    familywise_violations == 0
                    and pooled_flag_violations == 0
                    and mechanism_flag_violations == 0
                    and summary_mismatches == 0
                ),
            }
        )
    return pd.DataFrame(rows)


def build_mechanism_audit(mechanism: pd.DataFrame) -> pd.DataFrame:
    rows = []
    groups = mechanism.groupby(
        ["Scenario", "PlannedSlotsPerLanguage", "Language", "Domain"],
        sort=False,
    )
    for keys, frame in groups:
        planned_n = int(keys[1])
        spread = int(
            frame["ExposurePerReplicate"].max()
            - frame["ExposurePerReplicate"].min()
        )
        total = int(frame["ExposurePerReplicate"].sum())
        rows.append(
            {
                "Scenario": keys[0],
                "PlannedSlotsPerLanguage": planned_n,
                "Language": keys[2],
                "Domain": keys[3],
                "MechanismCount": frame["BlockedMechanism"].nunique(),
                "TotalExposurePerReplicate": total,
                "ExposureSpread": spread,
                "ActivatedMechanismSpecificGateCount": int(
                    frame["MechanismSpecificGateActivated"].sum()
                ),
                "Passed": bool(
                    frame["BlockedMechanism"].nunique() == 3
                    and total == planned_n
                    and spread <= 1
                    and not frame["MechanismSpecificGateActivated"].any()
                ),
            }
        )
    return pd.DataFrame(rows)


def build_precision_audit(
    summary: pd.DataFrame, cells: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    summary_prefixes = (
        "All12PrimaryPass",
        "All12FamilywisePass",
        "FalseReassuringPooledTruth",
        "FalseReassuringMechanismTruth",
    )
    cell_prefixes = ("PrimaryPass", "FamilywisePass")
    for table, frame, prefixes in (
        ("scenario_summary", summary, summary_prefixes),
        ("cell_summary", cells, cell_prefixes),
    ):
        for prefix in prefixes:
            rate = frame[f"{prefix}RateRaw"].to_numpy(dtype=float)
            replicates = frame[f"{prefix}Replicates"].to_numpy(dtype=float)
            mcse = frame[f"{prefix}MonteCarloSERaw"].to_numpy(dtype=float)
            lower = frame[f"{prefix}WilsonLower95Raw"].to_numpy(dtype=float)
            upper = frame[f"{prefix}WilsonUpper95Raw"].to_numpy(dtype=float)
            expected_mcse = np.sqrt(rate * (1.0 - rate) / replicates)
            passed = bool(
                np.all(replicates == 2000)
                and np.allclose(mcse, expected_mcse, rtol=0.0, atol=2e-15)
                and np.all(lower <= rate)
                and np.all(rate <= upper)
            )
            rows.append(
                {
                    "Table": table,
                    "EndpointPrefix": prefix,
                    "RowsChecked": len(frame),
                    "MinimumReplicates": int(replicates.min()),
                    "MaximumMCSE": float(mcse.max()),
                    "Passed": passed,
                }
            )
    return pd.DataFrame(rows)


def build_determinism_audit(
    first_scenario: dict[str, object], *, base_seed: int
) -> pd.DataFrame:
    seed = derive_condition_seed(base_seed, 0, 100)
    first = simulate_dependence_condition(
        first_scenario, 100, replicates=64, seed=seed
    )
    second = simulate_dependence_condition(
        first_scenario, 100, replicates=64, seed=seed
    )
    rows = []
    for name in first:
        first_hash = dataframe_sha256(first[name])
        second_hash = dataframe_sha256(second[name])
        rows.append(
            {
                "Table": name,
                "Seed": seed,
                "FirstSHA256": first_hash,
                "SecondSHA256": second_hash,
                "Passed": bool(first_hash == second_hash),
            }
        )
    return pd.DataFrame(rows)


def render_figures(summary: pd.DataFrame, output: Path) -> None:
    labels = {
        "independent_safe": "independent safe",
        "participant_dependence_safe": "participant dependence",
        "moderator_cluster_safe": "moderator cluster",
        "pooled_domain_hotspot": "domain hotspot",
        "blocked_mechanism_hotspot": "mechanism hotspot",
        "mnar_danger_underretained": "MNAR danger under-retained",
        "mnar_danger_overretained": "MNAR danger over-retained",
        "combined_adversarial": "combined adversarial",
    }
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    for scenario, frame in summary.groupby("Scenario", sort=False):
        ax.plot(
            frame["PlannedSlotsPerLanguage"],
            frame["All12PrimaryPassRateRaw"],
            marker="o",
            linewidth=1.5,
            label=labels.get(scenario, scenario),
        )
    ax.set(
        xlabel="Planned slots per language (not selected n)",
        ylabel="Empirical probability all 12 pooled primary cells pass",
        title="Dependence, heterogeneity, and MNAR change joint gate behavior",
        ylim=(-0.02, 1.02),
    )
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, ncol=2, loc="best")
    fig.tight_layout()
    fig.savefig(output / "joint_primary_pass_stress.png", dpi=180)
    plt.close(fig)

    n500 = summary.loc[summary["PlannedSlotsPerLanguage"].eq(500)].copy()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
    y_positions = np.arange(len(n500), dtype=float)
    complete_values = n500["AggregateCompleteDangerousProportionRaw"].to_numpy()
    observed_values = n500["AggregateObservedDangerousProportionRaw"].to_numpy()
    limit = max(
        0.12,
        float(
            n500[
                [
                    "AggregateCompleteDangerousProportionRaw",
                    "AggregateObservedDangerousProportionRaw",
                ]
            ].to_numpy().max()
        )
        + 0.01,
    )
    for y_value, complete_value, observed_value in zip(
        y_positions, complete_values, observed_values
    ):
        axes[0].plot(
            [complete_value, observed_value],
            [y_value, y_value],
            color="#9ca3af",
            linewidth=1.4,
            zorder=1,
        )
    axes[0].scatter(
        complete_values,
        y_positions,
        s=45,
        label="complete-data truth",
        zorder=2,
    )
    axes[0].scatter(
        observed_values,
        y_positions,
        s=50,
        marker="x",
        linewidth=1.8,
        label="observed valid sample",
        zorder=3,
    )
    axes[0].set_yticks(
        y_positions,
        [labels.get(value, value) for value in n500["Scenario"]],
        fontsize=8,
    )
    axes[0].invert_yaxis()
    axes[0].set(
        xlabel="Dangerous-error proportion",
        ylabel="",
        title="Outcome-dependent invalidity distorts observed risk (n=500)",
        xlim=(0, limit),
    )
    axes[0].grid(alpha=0.2)
    axes[0].legend(fontsize=8, loc="lower right")

    risk = summary.loc[
        summary["AnyPooledDomainTruthAtOrAboveThreshold"]
        | summary["AnyBlockedMechanismTruthAtOrAboveThreshold"]
    ]
    for scenario, frame in risk.groupby("Scenario", sort=False):
        false_rate = np.maximum(
            frame["FalseReassuringPooledTruthRateRaw"],
            frame["FalseReassuringMechanismTruthRateRaw"],
        )
        axes[1].plot(
            frame["PlannedSlotsPerLanguage"],
            false_rate,
            marker="o",
            label=labels.get(scenario, scenario),
        )
    axes[1].set(
        xlabel="Planned slots per language (not selected n)",
        ylabel="False-reassuring all-12 pooled-gate rate",
        title="Pooling and MNAR can preserve false reassurance",
        ylim=(-0.02, 1.02),
    )
    axes[1].grid(alpha=0.2)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "mnar_mechanism_false_reassurance.png", dpi=180)
    plt.close(fig)


def run_tests(output: Path) -> dict[str, object]:
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-q",
        "tests/test_cmle_one_click_confirmatory_dependence.py",
        "tests/test_cmle_one_click_confirmatory_planning.py",
        "tests/test_cmle_one_click_confirmatory_gate.py",
        "tests/test_cmle_one_click_cognitive_interview.py",
        "tests/test_cmle_one_click_comprehension.py",
        "tests/test_cmle_one_click.py",
        "tests/test_cmle_one_click_archive.py",
        "tests/test_decision_stability.py",
        "tests/test_threshold_decision_integration.py",
    ]
    completed = subprocess.run(
        command, cwd=ROOT, text=True, capture_output=True, check=False
    )
    (output / "selected_tests_stdout.txt").write_text(
        completed.stdout, encoding="utf-8"
    )
    (output / "selected_tests_stderr.txt").write_text(
        completed.stderr, encoding="utf-8"
    )
    return {
        "passed": completed.returncode == 0,
        "returncode": completed.returncode,
        "command": command,
    }


def _row(summary: pd.DataFrame, scenario: str, planned_n: int) -> pd.Series:
    return summary.loc[
        summary["Scenario"].eq(scenario)
        & summary["PlannedSlotsPerLanguage"].eq(planned_n)
    ].iloc[0]


def write_review(output: Path, results: dict[str, object]) -> None:
    text = f"""# Confirmatory dependence and MNAR stress critical review

## Decision

The no-human-data dependence-stress contract **{'passed' if results['contract_passed'] else 'failed'}**. It does not select a minimum valid n, a planned recruitment total, or a public claim.

## Dependence changes the joint gate

At 300 planned slots per language with marginal dangerous-error probability 0.05 and no missingness, the all-12 primary pass rate was `{results['independent_n300_joint']:.6f}` under independent outcomes, `{results['participant_n300_joint']:.6f}` with participant latent dependence, and `{results['cluster_n300_joint']:.6f}` with participant plus moderator/site-like clustering. The corresponding empirical joint-minus-product values are retained in the scenario table. These are conditional Gaussian-threshold results; latent rho is not an observed binary ICC.

## Outcome-dependent invalidity can reverse the conclusion

In the registered MNAR danger-under-retained scenario, complete-data truth is 0.10 in every cell, but the aggregate observed valid-sample rate at n=500 was `{results['mnar_n500_observed']:.6f}` because retention was 0.95 for safe outcomes and 0.55 for dangerous outcomes. The pooled all-12 gate was falsely reassuring in `{results['mnar_n500_false']:.6f}` of replicates (95% Monte Carlo Wilson interval `{results['mnar_n500_false_lower']:.6f}` to `{results['mnar_n500_false_upper']:.6f}`). A Wilson interval on selected valid records does not repair this MNAR bias.

## Pooled blocked mechanisms can hide a hotspot

The blocked-mechanism hotspot assigns truth probabilities 0.02, 0.02, and 0.16 to the three rotated failure mechanisms, yielding a pooled blocked-domain truth below 0.10. At n=500 the ordinary pooled 12-cell gate was nevertheless falsely reassuring about the unacceptable mechanism in `{results['mechanism_n500_false']:.6f}` of replicates. No mechanism-specific gate was retroactively activated; the result shows that a future protocol must decide prospectively whether pooled evidence is scientifically acceptable.

The combined adversarial scenario adds participant dependence, clusters of 10, and lower retention of dangerous responses. Its n=500 false-reassuring mechanism rate was `{results['combined_n500_false']:.6f}`. More planned observations do not eliminate selection bias or a misspecified pooled estimand.

## Limits

- All outcomes are synthetic and conditional on eight constructed scenarios; they do not estimate future user behavior.
- Identical English/Japanese generating parameters are a software check, not language equivalence.
- Consecutive fixed-size clusters and Gaussian latent correlations are stress devices, not validated moderator/site models.
- The MNAR mechanisms are not identifiable from observed records without external information.
- Monte Carlo intervals quantify simulation error only, not design, model, or recruitment uncertainty.

Human participants remain zero. No sample size is selected, no confirmatory result exists, and the public UI stays withheld.
"""
    (output / "CONFIRMATORY_DEPENDENCE_STRESS_CRITICAL_REVIEW.md").write_text(
        text, encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    plan = validate_registration()
    if args.output.exists():
        raise FileExistsError(f"Evidence output exists: {args.output}")
    args.output.mkdir(parents=True)

    fixed = plan["fixed_design"]
    scenarios = plan["scenarios"]
    planned_n_values = fixed["candidate_planned_slots_per_language"]
    replicates = int(fixed["monte_carlo_replicates_per_scenario_n"])
    base_seed = int(fixed["base_seed"])
    surface = simulate_registered_dependence_surface(
        scenarios,
        planned_n_values,
        replicates=replicates,
        base_seed=base_seed,
    )
    summary = surface["scenario_summary"]
    cells = surface["cell_summary"]
    mechanism = surface["mechanism_summary"]
    replicate = surface["replicate_summary"]
    human_gate = dependence_stress_human_gate_status()

    registration_audit = build_registration_audit(plan)
    seed_audit = build_seed_audit(
        summary, base_seed=base_seed, scenarios=scenarios
    )
    marginal_audit = build_marginal_calibration_audit(cells)
    baseline_audit = build_baseline_binomial_audit(cells, replicates)
    retention_audit = build_retention_calibration_audit(cells)
    logic_audit = build_logic_audit(summary, replicate)
    mechanism_audit = build_mechanism_audit(mechanism)
    precision_audit = build_precision_audit(summary, cells)
    determinism_audit = build_determinism_audit(
        scenarios[0], base_seed=base_seed
    )

    outputs = {
        "scenario_joint_gate_summary.csv": summary,
        "language_domain_cell_summary.csv": cells,
        "blocked_mechanism_summary.csv": mechanism,
        "replicate_joint_gate_audit.csv": replicate,
        "registration_identity_audit.csv": registration_audit,
        "seed_derivation_audit.csv": seed_audit,
        "marginal_probability_calibration_audit.csv": marginal_audit,
        "independent_exact_binomial_calibration_audit.csv": baseline_audit,
        "retention_mnar_calibration_audit.csv": retention_audit,
        "joint_gate_logic_audit.csv": logic_audit,
        "blocked_mechanism_allocation_audit.csv": mechanism_audit,
        "monte_carlo_precision_audit.csv": precision_audit,
        "determinism_audit.csv": determinism_audit,
        "human_gate_status.csv": human_gate,
    }
    for name, frame in outputs.items():
        write_csv(frame, args.output / name)
    render_figures(summary, args.output)
    tests = run_tests(args.output)

    expected_conditions = len(scenarios) * len(planned_n_values)
    human = human_gate.iloc[0]
    gates = {
        "identity_passed": bool(registration_audit["Passed"].all()),
        "registration_passed": bool(
            len(summary) == expected_conditions
            and len(cells) == expected_conditions * len(LANGUAGES) * len(DOMAINS)
            and len(mechanism)
            == expected_conditions * len(LANGUAGES) * 2 * 3
            and len(replicate) == expected_conditions * replicates
            and seed_audit["Passed"].all()
        ),
        "determinism_passed": bool(determinism_audit["Passed"].all()),
        "marginal_calibration_passed": bool(marginal_audit["Passed"].all()),
        "baseline_binomial_calibration_passed": bool(
            len(baseline_audit) == len(planned_n_values) * len(LANGUAGES) * len(DOMAINS)
            and baseline_audit["Passed"].all()
        ),
        "retention_calibration_passed": bool(retention_audit["Passed"].all()),
        "joint_logic_passed": bool(
            logic_audit["Passed"].all()
            and logic_audit["ReplicateRows"].eq(replicates).all()
        ),
        "mechanism_allocation_passed": bool(mechanism_audit["Passed"].all()),
        "monte_carlo_precision_passed": bool(precision_audit["Passed"].all()),
        "human_and_selection_withheld_passed": bool(
            int(human["HumanParticipants"]) == 0
            and not bool(human["MinimumValidPerCellRegistered"])
            and not bool(human["PlannedRecruitmentNSelected"])
            and not bool(human["ConfirmatoryResultAvailable"])
            and not bool(human["LanguageEquivalenceEstablished"])
            and not bool(human["ClusterModelValidated"])
            and not bool(human["MNARModelIdentified"])
            and not bool(human["PublicSurfaceEnabled"])
            and not summary["MinimumValidNRegistered"].any()
            and not summary["PlannedRecruitmentNSelected"].any()
            and not summary["PublicSurfaceEnabled"].any()
        ),
        "tests_passed": bool(tests["passed"]),
    }
    contract_passed = all(gates.values())

    independent_n300 = _row(summary, "independent_safe", 300)
    participant_n300 = _row(summary, "participant_dependence_safe", 300)
    cluster_n300 = _row(summary, "moderator_cluster_safe", 300)
    mnar_n500 = _row(summary, "mnar_danger_underretained", 500)
    mechanism_n500 = _row(summary, "blocked_mechanism_hotspot", 500)
    combined_n500 = _row(summary, "combined_adversarial", 500)
    pooled_hotspot_n500 = _row(summary, "pooled_domain_hotspot", 500)
    results = {
        **gates,
        "contract_passed": contract_passed,
        "scenario_count": len(scenarios),
        "planned_n_count": len(planned_n_values),
        "condition_count": expected_conditions,
        "replicates_per_condition": replicates,
        "total_synthetic_replicates": len(replicate),
        "cell_summary_rows": len(cells),
        "mechanism_summary_rows": len(mechanism),
        "independent_n300_joint": float(independent_n300["All12PrimaryPassRateRaw"]),
        "participant_n300_joint": float(participant_n300["All12PrimaryPassRateRaw"]),
        "cluster_n300_joint": float(cluster_n300["All12PrimaryPassRateRaw"]),
        "mnar_n500_complete": float(
            mnar_n500["AggregateCompleteDangerousProportionRaw"]
        ),
        "mnar_n500_observed": float(
            mnar_n500["AggregateObservedDangerousProportionRaw"]
        ),
        "mnar_n500_false": float(
            mnar_n500["FalseReassuringPooledTruthRateRaw"]
        ),
        "mnar_n500_false_lower": float(
            mnar_n500["FalseReassuringPooledTruthWilsonLower95Raw"]
        ),
        "mnar_n500_false_upper": float(
            mnar_n500["FalseReassuringPooledTruthWilsonUpper95Raw"]
        ),
        "mechanism_n500_false": float(
            mechanism_n500["FalseReassuringMechanismTruthRateRaw"]
        ),
        "combined_n500_false": float(
            combined_n500["FalseReassuringMechanismTruthRateRaw"]
        ),
        "pooled_hotspot_n500_false": float(
            pooled_hotspot_n500["FalseReassuringPooledTruthRateRaw"]
        ),
        "maximum_reported_monte_carlo_se": float(
            precision_audit["MaximumMCSE"].max()
        ),
        "human_participants": 0,
        "minimum_valid_per_cell_registered": False,
        "planned_recruitment_n_selected": False,
        "confirmatory_result_available": False,
        "public_surface_enabled": False,
        "selected_tests": tests,
    }
    write_review(args.output, results)
    output_files = sorted(
        path.relative_to(args.output).as_posix()
        for path in args.output.rglob("*")
        if path.is_file() and path.name != "decision.json"
    )
    decision = json_safe(
        {
            "study_id": plan["study_id"],
            "plan_sha256": sha256_file(PLAN),
            "contract_passed": contract_passed,
            "contract_interpretation": "synthetic_dependence_mnar_stress_no_n_selected",
            "implementation_sha256": {
                "mfrm_app/cmle_one_click_confirmatory_dependence.py": sha256_file(
                    ROOT / "mfrm_app/cmle_one_click_confirmatory_dependence.py"
                ),
                "tests/test_cmle_one_click_confirmatory_dependence.py": sha256_file(
                    ROOT / "tests/test_cmle_one_click_confirmatory_dependence.py"
                ),
                "validation/cmle_one_click_confirmatory_dependence_stress.py": sha256_file(
                    Path(__file__)
                ),
            },
            "results": results,
            "output_sha256": {
                name: sha256_file(args.output / name) for name in output_files
            },
            "interpretation": {
                "simulation": "conditional synthetic stress, not user evidence",
                "dependence": "latent Gaussian rho is not observed binary ICC",
                "mnar": "observed-data Wilson gate does not correct informative invalidity",
                "mechanism": "pooled pass may conceal an unacceptable blocked mechanism",
                "sample_size": "not selected",
                "human_data": "none",
                "public_ui": "withheld",
            },
        }
    )
    (args.output / "decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(decision, indent=2, ensure_ascii=False, allow_nan=False))
    if not contract_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
