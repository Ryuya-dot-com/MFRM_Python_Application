#!/usr/bin/env python3
"""Run no-human-data sample-size sensitivity for the directional gate."""

from __future__ import annotations

import argparse
import hashlib
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

from mfrm_app.cmle_one_click_confirmatory_gate import (  # noqa: E402
    maximum_passing_errors,
)
from mfrm_app.cmle_one_click_confirmatory_planning import (  # noqa: E402
    FAMILYWISE_CONFIDENCE,
    PRIMARY_CONFIDENCE,
    binomial_cdf_direct,
    binomial_cdf_logspace,
    build_attrition_assurance_table,
    build_blocked_mechanism_exposure_table,
    build_cell_probability_surface,
    build_cluster_sensitivity_table,
    build_minimum_n_sensitivity,
    cell_gate_pass_probability,
    cluster_design_effect,
    confirmatory_planning_human_gate_status,
    maximum_passing_errors_fast,
)


PLAN = ROOT / "validation/cmle_one_click_confirmatory_sample_size_sensitivity_plan_20260810.json"
AMENDMENT = ROOT / "validation/cmle_one_click_confirmatory_sample_size_sensitivity_amendment_20260810.json"
OUTPUT = ROOT / "validation/cmle_one_click_confirmatory_sample_size_sensitivity_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def validate_registration() -> tuple[dict[str, object], dict[str, object]]:
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    amendment = json.loads(AMENDMENT.read_text(encoding="utf-8"))
    if plan.get("study_id") != "cmle_one_click_confirmatory_sample_size_sensitivity_v1":
        raise ValueError("Unexpected confirmatory planning study.")
    if amendment.get("study_id") != plan["study_id"]:
        raise ValueError("Amendment study identity mismatch.")
    if amendment.get("parent_plan_sha256") != sha256_file(PLAN):
        raise ValueError("Amendment parent plan identity mismatch.")
    mismatches = [
        relative
        for relative, expected in plan["parent_identity"].items()
        if not (ROOT / relative).is_file()
        or sha256_file(ROOT / relative) != expected
    ]
    if mismatches:
        raise ValueError(f"Confirmatory planning parent identity failed: {mismatches}")
    return plan, amendment


def probability_formula_audit() -> pd.DataFrame:
    rows = []
    for n in (5, 10, 25, 50, 100):
        for probability in (0.01, 0.05, 0.10, 0.50):
            for maximum in (0, 1, min(3, n), n):
                direct = binomial_cdf_direct(maximum, n, probability)
                logspace = binomial_cdf_logspace(maximum, n, probability)
                rows.append(
                    {
                        "Trials": n,
                        "TrueErrorProbability": probability,
                        "MaximumErrors": maximum,
                        "DirectSumProbability": direct,
                        "LogspaceProbability": logspace,
                        "AbsoluteDifference": abs(direct - logspace),
                        "Passed": math.isclose(
                            direct, logspace, rel_tol=2e-13, abs_tol=2e-15
                        ),
                    }
                )
    return pd.DataFrame(rows)


def passing_count_audit(candidate_n: list[int]) -> pd.DataFrame:
    rows = []
    for n in candidate_n:
        for scheme, confidence in (
            ("primary_per_cell", PRIMARY_CONFIDENCE),
            ("bonferroni_12_cell_sensitivity", FAMILYWISE_CONFIDENCE),
        ):
            exhaustive = maximum_passing_errors(n, confidence=confidence)
            fast = maximum_passing_errors_fast(n, confidence=confidence)
            rows.append(
                {
                    "ValidEligibleN": n,
                    "ConfidenceScheme": scheme,
                    "ExhaustiveMaximumPassingErrors": exhaustive,
                    "BinarySearchMaximumPassingErrors": fast,
                    "Passed": exhaustive == fast,
                }
            )
    return pd.DataFrame(rows)


def crossing_contract_audit(
    cell_minimum: pd.DataFrame, joint_minimum: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    for table_name, frame, target_column in (
        ("cell", cell_minimum, "CellPassProbabilityTarget"),
        ("joint_union_bound", joint_minimum, "JointUnionBoundTarget"),
    ):
        for index, row in frame.iterrows():
            target = float(row[target_column])
            first_ok = bool(
                not row["FirstCrossingAvailable"]
                or (
                    row["ProbabilityAtFirstCrossing"] >= target
                    and row["ProbabilityAtFirstCrossingPredecessor"] < target
                )
            )
            sustained_ok = bool(
                not row["SustainedThroughSearchAvailable"]
                or row["MinimumProbabilityFromSustainedThroughSearchMax"] >= target
            )
            rows.append(
                {
                    "Table": table_name,
                    "SourceRow": index,
                    "Target": target,
                    "FirstCrossingContractPassed": first_ok,
                    "SustainedThroughSearchContractPassed": sustained_ok,
                    "DownwardAdjacentStepCount": row["DownwardAdjacentStepCount"],
                    "MaximumDownwardAdjacentStep": row["MaximumDownwardAdjacentStep"],
                    "Passed": first_ok and sustained_ok,
                }
            )
    return pd.DataFrame(rows)


def attrition_contract_audit(attrition: pd.DataFrame) -> pd.DataFrame:
    result = attrition.copy()
    result["MinimumMeetsTarget"] = result["AssuranceAtMinimumRaw"] >= result[
        "AssuranceTarget"
    ]
    result["PredecessorFailsTarget"] = result[
        "AssuranceAtPredecessorRaw"
    ] < result["AssuranceTarget"]
    result["Passed"] = result["MinimumMeetsTarget"] & result[
        "PredecessorFailsTarget"
    ]
    return result


def cluster_contract_audit(cluster: pd.DataFrame) -> pd.DataFrame:
    result = cluster.copy()
    result["RecomputedDesignEffect"] = 1.0 + (
        result["AverageClusterSizeAssumed"] - 1.0
    ) * result["ICCAssumed"]
    result["RecomputedNominalN"] = np.ceil(
        result["EffectiveTargetN"] * result["RecomputedDesignEffect"]
    ).astype(int)
    result["Passed"] = np.isclose(
        result["DesignEffect"], result["RecomputedDesignEffect"], rtol=0, atol=1e-15
    ) & result["HeuristicNominalValidN"].eq(result["RecomputedNominalN"])
    return result


def mechanism_contract_audit(mechanism: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (candidate, language), group in mechanism.groupby(
        ["CandidateSlotsPerLanguage", "Language"], sort=False
    ):
        rows.append(
            {
                "CandidateSlotsPerLanguage": candidate,
                "Language": language,
                "MechanismCount": group["BlockedMechanism"].nunique(),
                "TotalExposure": int(group["MechanismExposureN"].sum()),
                "MinimumExposure": int(group["MechanismExposureN"].min()),
                "MaximumExposure": int(group["MechanismExposureN"].max()),
                "ExposureSpread": int(
                    group["MechanismExposureN"].max()
                    - group["MechanismExposureN"].min()
                ),
                "Passed": bool(
                    group["BlockedMechanism"].nunique() == 3
                    and int(group["MechanismExposureN"].sum()) == int(candidate)
                    and int(
                        group["MechanismExposureN"].max()
                        - group["MechanismExposureN"].min()
                    )
                    <= 1
                    and not group["MechanismSpecificPrimaryGateActivated"].any()
                ),
            }
        )
    return pd.DataFrame(rows)


def build_sawtooth_plot_data() -> pd.DataFrame:
    rows = []
    for probability in (0.025, 0.05, 0.075):
        for n in range(20, 601):
            result = cell_gate_pass_probability(n, probability)
            rows.append(
                {
                    "TrueDangerousErrorProbability": probability,
                    "ValidEligibleN": n,
                    "MaximumPassingErrors": result["MaximumPassingErrors"],
                    "CellPassProbabilityRaw": result["CellPassProbabilityRaw"],
                }
            )
    return pd.DataFrame(rows)


def render_figures(
    sawtooth: pd.DataFrame,
    attrition: pd.DataFrame,
    cluster: pd.DataFrame,
    output: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(10, 5.8))
    for probability, group in sawtooth.groupby("TrueDangerousErrorProbability"):
        ax.plot(
            group["ValidEligibleN"],
            group["CellPassProbabilityRaw"],
            linewidth=1.3,
            label=f"assumed p={probability:.3f}",
        )
    for target in (0.80, 0.90, 0.95):
        ax.axhline(target, color="#6b7280", linewidth=0.8, linestyle="--")
    ax.set(
        xlabel="Valid eligible n per cell",
        ylabel="Exact probability the primary Wilson cell passes",
        title="Discrete Wilson pass probability is sawtooth, not monotone",
        xlim=(20, 600),
        ylim=(0, 1.02),
    )
    ax.grid(alpha=0.2)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(output / "wilson_pass_probability_sawtooth.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    attrition_plot = attrition.loc[attrition["AssuranceTarget"].eq(0.95)]
    for retention, group in attrition_plot.groupby("RetentionProbabilityAssumed"):
        axes[0].plot(
            group["TargetValidN"],
            group["MinimumPlannedSlots"],
            marker="o",
            label=f"retention={retention:.2f}",
        )
    axes[0].plot([0, 520], [0, 520], color="#6b7280", linestyle="--", linewidth=0.8)
    axes[0].set(
        xlabel="Target valid n",
        ylabel="Minimum planned slots (95% assurance)",
        title="Attrition arithmetic sensitivity",
        xlim=(0, 520),
    )
    axes[0].grid(alpha=0.2)
    axes[0].legend()

    cluster_plot = cluster.loc[
        cluster["EffectiveTargetN"].eq(100)
        & cluster["AverageClusterSizeAssumed"].gt(1)
    ]
    for size, group in cluster_plot.groupby("AverageClusterSizeAssumed"):
        axes[1].plot(
            group["ICCAssumed"],
            group["HeuristicNominalValidN"],
            marker="o",
            label=f"mean cluster={size}",
        )
    axes[1].set(
        xlabel="Assumed ICC",
        ylabel="Heuristic nominal valid n for effective n=100",
        title="Cluster design-effect heuristic (not exact gate correction)",
    )
    axes[1].grid(alpha=0.2)
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(output / "attrition_cluster_sensitivity.png", dpi=180)
    plt.close(fig)


def run_tests(output: Path) -> dict[str, object]:
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-q",
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
    (output / "selected_tests_stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output / "selected_tests_stderr.txt").write_text(completed.stderr, encoding="utf-8")
    return {"passed": completed.returncode == 0, "returncode": completed.returncode, "command": command}


def write_review(output: Path, results: dict[str, object]) -> None:
    text = f"""# Confirmatory sample-size sensitivity critical review

## Decision

The no-human-data sensitivity contract **{'passed' if results['contract_passed'] else 'failed'}**. It maps assumptions; it does not select a minimum valid n or planned recruitment n.

## Discreteness changes the planning question

The probability of passing a Wilson cell is sawtooth because the maximum allowed integer error count changes in steps. Under assumed true dangerous-error probability `p=0.05`, a 90% per-cell pass target first crosses at `n={results['p05_cell90_first_crossing']}`, but it is not continuously at or above 90% through the registered search range until `n={results['p05_cell90_sustained']}`. For the dependence-robust 12-cell union-bound target of 90%, the corresponding values are `n={results['p05_joint90_first_crossing']}` and `n={results['p05_joint90_sustained']}`. A first crossing is therefore not a defensible stand-alone sample-size choice.

## Attrition and clustering are separate assumptions

For a target of 200 valid slots, assumed homogeneous 90% retention, and 95% assurance, the arithmetic minimum planned slots is {results['attrition_target200_retention90_assurance95']}. This does not show attrition is random. The cluster table uses `1+(m-1)ICC` only as a heuristic; for effective n=100, mean cluster size 10, and ICC 0.05, it gives nominal n=145. That is not an exact correction for a binary, discontinuous Wilson decision.

## Blocked-state heterogeneity can be hidden

The combined blocked-primary cell rotates boundary, structural-design, and invalid-input mechanisms. With 75 slots per language, each receives 25 exposures. Smaller totals cannot give all three mechanisms 25 observations. The current mechanism table is a warning only; it does not silently add a new primary gate. A future protocol must decide whether pooled false-green evidence is acceptable or mechanism-specific protection is required.

## Remaining risks

- Exact-binomial probabilities assume a common fixed p within each cell and do not model participant heterogeneity, learning, moderator effects, device/accessibility conditions, or recruitment selection.
- The independence projection is optimistic when cell pass events are dependent. The union-bound lower assurance is conservative but cannot repair biased marginal estimates.
- First-crossing and sustained-through-5000 results are both conditional on the assumed p and search range. Sustained through 5000 is not a proof beyond 5000.
- Attrition assurance assumes independent homogeneous retention; invalidity may be informative.
- Separate English/Japanese planning does not establish language equivalence.

The full surface is retained to prevent post-result assumption cherry-picking. Human participants remain zero, no n is selected, and the public UI stays withheld.
"""
    (output / "CONFIRMATORY_SAMPLE_SIZE_SENSITIVITY_CRITICAL_REVIEW.md").write_text(
        text, encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    plan, amendment = validate_registration()
    if args.output.exists():
        raise FileExistsError(f"Evidence output exists: {args.output}")
    args.output.mkdir(parents=True)

    probabilities = plan["exact_binomial_pass_probability"][
        "true_error_probabilities"
    ]
    candidate_n = plan["exact_binomial_pass_probability"]["candidate_valid_n"]
    targets = plan["exact_binomial_pass_probability"]["cell_probability_targets"]
    search_max = plan["exact_binomial_pass_probability"]["search_max_n"]
    surface = build_cell_probability_surface(probabilities, candidate_n)
    cell_minimum, joint_minimum = build_minimum_n_sensitivity(
        probabilities, targets, search_max_n=search_max
    )
    attrition = build_attrition_assurance_table(
        plan["attrition_assurance"]["target_valid_n"],
        plan["attrition_assurance"]["retention_probabilities"],
        plan["attrition_assurance"]["assurance_targets"],
    )
    cluster = build_cluster_sensitivity_table(
        plan["cluster_sensitivity"]["effective_targets"],
        plan["cluster_sensitivity"]["average_cluster_sizes"],
        plan["cluster_sensitivity"]["intraclass_correlations"],
    )
    mechanism = build_blocked_mechanism_exposure_table(
        plan["blocked_mechanism_heterogeneity"]["candidate_slots_per_language"]
    )
    formula_audit = probability_formula_audit()
    count_audit = passing_count_audit(candidate_n)
    crossing_audit = crossing_contract_audit(cell_minimum, joint_minimum)
    attrition_audit = attrition_contract_audit(attrition)
    cluster_audit = cluster_contract_audit(cluster)
    mechanism_audit = mechanism_contract_audit(mechanism)
    sawtooth = build_sawtooth_plot_data()
    human_gate = confirmatory_planning_human_gate_status()

    write_csv(surface, args.output / "cell_pass_probability_surface.csv")
    write_csv(cell_minimum, args.output / "cell_probability_target_crossings.csv")
    write_csv(joint_minimum, args.output / "joint_union_bound_target_crossings.csv")
    write_csv(attrition, args.output / "attrition_assurance_sensitivity.csv")
    write_csv(cluster, args.output / "cluster_design_effect_sensitivity.csv")
    write_csv(mechanism, args.output / "blocked_mechanism_exposure_sensitivity.csv")
    write_csv(formula_audit, args.output / "binomial_probability_formula_audit.csv")
    write_csv(count_audit, args.output / "maximum_passing_error_count_audit.csv")
    write_csv(crossing_audit, args.output / "crossing_and_sustained_contract_audit.csv")
    write_csv(attrition_audit, args.output / "attrition_minimality_audit.csv")
    write_csv(cluster_audit, args.output / "cluster_arithmetic_audit.csv")
    write_csv(mechanism_audit, args.output / "blocked_mechanism_balance_audit.csv")
    write_csv(sawtooth, args.output / "wilson_sawtooth_plot_data.csv")
    write_csv(human_gate, args.output / "human_gate_status.csv")
    render_figures(sawtooth, attrition, cluster, args.output)
    tests = run_tests(args.output)

    p05_cell90 = cell_minimum.loc[
        cell_minimum["TrueDangerousErrorProbability"].eq(0.05)
        & cell_minimum["ConfidenceScheme"].eq("primary_per_cell")
        & cell_minimum["CellPassProbabilityTarget"].eq(0.90)
    ].iloc[0]
    p05_joint90 = joint_minimum.loc[
        joint_minimum["TrueDangerousErrorProbability"].eq(0.05)
        & joint_minimum["JointUnionBoundTarget"].eq(0.90)
    ].iloc[0]
    attrition_example = attrition.loc[
        attrition["TargetValidN"].eq(200)
        & attrition["RetentionProbabilityAssumed"].eq(0.90)
        & attrition["AssuranceTarget"].eq(0.95)
    ].iloc[0]
    cluster_example = cluster.loc[
        cluster["EffectiveTargetN"].eq(100)
        & cluster["AverageClusterSizeAssumed"].eq(10)
        & cluster["ICCAssumed"].eq(0.05)
    ].iloc[0]
    mechanism75 = mechanism.loc[
        mechanism["CandidateSlotsPerLanguage"].eq(75)
    ]

    gates = {
        "identity_passed": True,
        "registration_and_amendment_passed": bool(
            amendment["parent_plan_sha256"] == sha256_file(PLAN)
            and amendment["study_id"] == plan["study_id"]
        ),
        "binomial_probability_passed": bool(
            len(formula_audit) == 80
            and formula_audit["Passed"].all()
            and count_audit["Passed"].all()
        ),
        "surface_passed": bool(
            len(surface) == len(probabilities) * len(candidate_n) * 2
            and surface["CellPassProbabilityRaw"].between(0, 1).all()
            and not surface["SampleSizeSelected"].any()
        ),
        "crossing_and_joint_passed": bool(
            len(cell_minimum) == len(probabilities) * len(targets) * 2
            and len(joint_minimum) == len(probabilities) * len(targets)
            and crossing_audit["Passed"].all()
            and crossing_audit["DownwardAdjacentStepCount"].gt(0).any()
            and not cell_minimum["SampleSizeSelected"].any()
            and not joint_minimum["SampleSizeSelected"].any()
        ),
        "attrition_passed": bool(
            len(attrition) == 72
            and attrition_audit["Passed"].all()
            and not attrition["PlannedRecruitmentNSelected"].any()
        ),
        "cluster_passed": bool(
            len(cluster) == 96
            and cluster_audit["Passed"].all()
            and cluster["Interpretation"].str.startswith("heuristic_").all()
            and not cluster["SampleSizeSelected"].any()
        ),
        "mechanism_passed": bool(
            len(mechanism) == 48
            and mechanism_audit["Passed"].all()
            and mechanism75["MechanismExposureN"].eq(25).all()
        ),
        "human_and_selection_withheld_passed": bool(
            len(human_gate) == 1
            and int(human_gate.iloc[0]["HumanParticipants"]) == 0
            and not bool(human_gate.iloc[0]["MinimumValidPerCellRegistered"])
            and not bool(human_gate.iloc[0]["PlannedRecruitmentNSelected"])
            and not bool(human_gate.iloc[0]["PilotResultAvailable"])
            and not bool(human_gate.iloc[0]["ConfirmatoryResultAvailable"])
            and not bool(human_gate.iloc[0]["PublicSurfaceEnabled"])
        ),
        "tests_passed": bool(tests["passed"]),
    }
    contract_passed = all(gates.values())
    results = {
        **gates,
        "contract_passed": contract_passed,
        "surface_rows": len(surface),
        "cell_target_rows": len(cell_minimum),
        "joint_target_rows": len(joint_minimum),
        "p05_cell90_first_crossing": int(p05_cell90["FirstCrossingN"]),
        "p05_cell90_sustained": int(p05_cell90["SustainedThroughSearchFromN"]),
        "p05_joint90_first_crossing": int(p05_joint90["FirstCrossingN"]),
        "p05_joint90_sustained": int(p05_joint90["SustainedThroughSearchFromN"]),
        "p05_cell_maximum_downward_step": float(
            p05_cell90["MaximumDownwardAdjacentStep"]
        ),
        "attrition_target200_retention90_assurance95": int(
            attrition_example["MinimumPlannedSlots"]
        ),
        "cluster_effective100_m10_icc05_nominal": int(
            cluster_example["HeuristicNominalValidN"]
        ),
        "mechanism_exposure_at_75_slots": 25,
        "human_participants": 0,
        "minimum_valid_per_cell_registered": False,
        "planned_recruitment_n_selected": False,
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
            "amendment_sha256": sha256_file(AMENDMENT),
            "contract_passed": contract_passed,
            "contract_interpretation": "no_human_data_sensitivity_surface_no_n_selected",
            "implementation_sha256": {
                "mfrm_app/cmle_one_click_confirmatory_planning.py": sha256_file(
                    ROOT / "mfrm_app/cmle_one_click_confirmatory_planning.py"
                ),
                "tests/test_cmle_one_click_confirmatory_planning.py": sha256_file(
                    ROOT / "tests/test_cmle_one_click_confirmatory_planning.py"
                ),
                "validation/cmle_one_click_confirmatory_sample_size_sensitivity.py": sha256_file(Path(__file__)),
            },
            "results": results,
            "output_sha256": {
                name: sha256_file(args.output / name) for name in output_files
            },
            "interpretation": {
                "first_crossing": "not a monotone or selected minimum n",
                "sustained": "only through registered search maximum 5000",
                "attrition": "independent homogeneous retention sensitivity only",
                "cluster": "heuristic design effect, not exact Wilson correction",
                "mechanism": "pooled blocked cells may mask heterogeneity",
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
