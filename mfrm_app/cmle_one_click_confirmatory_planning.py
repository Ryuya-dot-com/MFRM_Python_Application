"""No-human-data sample-size sensitivity for the directional CMLE gate."""

from __future__ import annotations

from collections.abc import Sequence
import math

import numpy as np
import pandas as pd
from scipy.special import gammaln, logsumexp

from mfrm_app.cmle_one_click_comprehension import HUMAN_STUDY_STATUS
from mfrm_app.cmle_one_click_confirmatory_gate import (
    BLOCKED_PRIMARY_CASES,
    build_confirmatory_assignment_schedule,
    wilson_upper_bound,
)


CONFIRMATORY_PLANNING_SCHEMA_VERSION = "cmle_confirmatory_sample_size_sensitivity_v1"
PRIMARY_CONFIDENCE = 0.95
FAMILYWISE_CONFIDENCE = 1.0 - 0.05 / 12.0
DECISION_THRESHOLD = 0.10


def binomial_cdf_logspace(max_errors: int, trials: int, probability: float) -> float:
    """Return P[X <= max_errors] for a Binomial(trials, probability)."""

    if not isinstance(max_errors, int) or not isinstance(trials, int):
        raise ValueError("max_errors and trials must be integers.")
    if trials < 0 or not 0.0 <= probability <= 1.0:
        raise ValueError("Require trials >= 0 and probability in [0, 1].")
    if max_errors < 0:
        return 0.0
    if max_errors >= trials:
        return 1.0
    if probability == 0.0:
        return 1.0
    if probability == 1.0:
        return 0.0
    values = np.arange(max_errors + 1, dtype=float)
    n = float(trials)
    log_terms = (
        gammaln(n + 1.0)
        - gammaln(values + 1.0)
        - gammaln(n - values + 1.0)
        + values * math.log(probability)
        + (n - values) * math.log1p(-probability)
    )
    return float(min(1.0, math.exp(float(logsumexp(log_terms)))))


def binomial_cdf_direct(max_errors: int, trials: int, probability: float) -> float:
    """Small-n independent direct-sum audit implementation."""

    if max_errors < 0:
        return 0.0
    upper = min(max_errors, trials)
    return float(
        sum(
            math.comb(trials, errors)
            * probability**errors
            * (1.0 - probability) ** (trials - errors)
            for errors in range(upper + 1)
        )
    )


def maximum_passing_errors_fast(
    trials: int,
    *,
    confidence: float = PRIMARY_CONFIDENCE,
    threshold: float = DECISION_THRESHOLD,
) -> int:
    """Binary-search the largest count with raw Wilson upper < threshold."""

    if not isinstance(trials, int) or trials <= 0:
        raise ValueError("trials must be a positive integer.")
    if wilson_upper_bound(0, trials, confidence=confidence) >= threshold:
        return -1
    low, high = 0, trials
    while low < high:
        middle = (low + high + 1) // 2
        if wilson_upper_bound(middle, trials, confidence=confidence) < threshold:
            low = middle
        else:
            high = middle - 1
    return low


def cell_gate_pass_probability(
    trials: int,
    true_error_probability: float,
    *,
    confidence: float = PRIMARY_CONFIDENCE,
    threshold: float = DECISION_THRESHOLD,
) -> dict[str, float | int]:
    """Return the exact-binomial probability that one Wilson cell passes."""

    maximum = maximum_passing_errors_fast(
        trials, confidence=confidence, threshold=threshold
    )
    probability = binomial_cdf_logspace(maximum, trials, true_error_probability)
    return {
        "MaximumPassingErrors": maximum,
        "CellPassProbabilityRaw": probability,
    }


def joint_pass_projections(cell_pass_probability: float, *, cell_count: int = 12) -> dict[str, float]:
    """Separate independence projection from a dependence-robust union bound."""

    if not 0.0 <= cell_pass_probability <= 1.0 or cell_count <= 0:
        raise ValueError("Invalid cell probability or cell_count.")
    return {
        "IndependenceProjectionRaw": cell_pass_probability**cell_count,
        "UnionBoundLowerRaw": max(
            0.0, 1.0 - cell_count * (1.0 - cell_pass_probability)
        ),
    }


def build_cell_probability_surface(
    true_error_probabilities: Sequence[float] = (0.0, 0.01, 0.025, 0.05, 0.075, 0.10),
    candidate_valid_n: Sequence[int] = (24, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500),
) -> pd.DataFrame:
    """Return the full registered p x n x confidence-scheme surface."""

    schemes = (
        ("primary_per_cell", PRIMARY_CONFIDENCE),
        ("bonferroni_12_cell_sensitivity", FAMILYWISE_CONFIDENCE),
    )
    rows = []
    for probability in true_error_probabilities:
        for n in candidate_valid_n:
            for scheme, confidence in schemes:
                result = cell_gate_pass_probability(
                    int(n), float(probability), confidence=confidence
                )
                joint = joint_pass_projections(
                    float(result["CellPassProbabilityRaw"])
                )
                rows.append(
                    {
                        "TrueDangerousErrorProbability": float(probability),
                        "ValidEligibleN": int(n),
                        "ConfidenceScheme": scheme,
                        "OneSidedConfidence": confidence,
                        "MaximumPassingErrors": result["MaximumPassingErrors"],
                        "CellPassProbabilityRaw": result["CellPassProbabilityRaw"],
                        **joint,
                        "DecisionThresholdRaw": DECISION_THRESHOLD,
                        "SampleSizeSelected": False,
                    }
                )
    return pd.DataFrame(rows)


def _probability_sequence(
    true_error_probability: float,
    confidence: float,
    search_max_n: int,
) -> list[float]:
    return [
        float(
            cell_gate_pass_probability(
                n, true_error_probability, confidence=confidence
            )["CellPassProbabilityRaw"]
        )
        for n in range(1, search_max_n + 1)
    ]


def _crossing_summary(values: Sequence[float], target: float) -> dict[str, object]:
    first_matches = [index + 1 for index, value in enumerate(values) if value >= target]
    first = first_matches[0] if first_matches else None
    suffix_minimum = [0.0] * len(values)
    running = 1.0
    for index in range(len(values) - 1, -1, -1):
        running = min(running, float(values[index]))
        suffix_minimum[index] = running
    sustained_matches = [
        index + 1 for index, value in enumerate(suffix_minimum) if value >= target
    ]
    sustained = sustained_matches[0] if sustained_matches else None
    downward = [
        max(0.0, float(left) - float(right))
        for left, right in zip(values, values[1:])
    ]
    return {
        "FirstCrossingN": first if first is not None else np.nan,
        "FirstCrossingAvailable": first is not None,
        "ProbabilityAtFirstCrossing": values[first - 1] if first is not None else np.nan,
        "ProbabilityAtFirstCrossingPredecessor": (
            values[first - 2] if first is not None and first > 1 else 0.0
        ),
        "SustainedThroughSearchFromN": sustained if sustained is not None else np.nan,
        "SustainedThroughSearchAvailable": sustained is not None,
        "MinimumProbabilityFromSustainedThroughSearchMax": (
            suffix_minimum[sustained - 1] if sustained is not None else np.nan
        ),
        "DownwardAdjacentStepCount": sum(value > 0.0 for value in downward),
        "MaximumDownwardAdjacentStep": max(downward, default=0.0),
    }


def build_minimum_n_sensitivity(
    true_error_probabilities: Sequence[float] = (0.0, 0.01, 0.025, 0.05, 0.075, 0.10),
    probability_targets: Sequence[float] = (0.80, 0.90, 0.95),
    *,
    search_max_n: int = 5000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Search cell and union-bound targets without selecting a design row."""

    cell_rows = []
    joint_rows = []
    schemes = (
        ("primary_per_cell", PRIMARY_CONFIDENCE),
        ("bonferroni_12_cell_sensitivity", FAMILYWISE_CONFIDENCE),
    )
    primary_cache: dict[float, list[float]] = {}
    for probability in true_error_probabilities:
        for scheme, confidence in schemes:
            values = _probability_sequence(
                float(probability), confidence, search_max_n
            )
            if scheme == "primary_per_cell":
                primary_cache[float(probability)] = values
            for target in probability_targets:
                crossing = _crossing_summary(values, float(target))
                cell_rows.append(
                    {
                        "TrueDangerousErrorProbability": float(probability),
                        "ConfidenceScheme": scheme,
                        "CellPassProbabilityTarget": float(target),
                        **crossing,
                        "SearchMaxN": search_max_n,
                        "SampleSizeSelected": False,
                    }
                )
    for probability, values in primary_cache.items():
        union_values = [
            joint_pass_projections(value)["UnionBoundLowerRaw"] for value in values
        ]
        independence_values = [
            joint_pass_projections(value)["IndependenceProjectionRaw"]
            for value in values
        ]
        for target in probability_targets:
            crossing = _crossing_summary(union_values, float(target))
            first = crossing["FirstCrossingN"]
            first_index = int(first) - 1 if crossing["FirstCrossingAvailable"] else None
            joint_rows.append(
                {
                    "TrueDangerousErrorProbability": probability,
                    "JointUnionBoundTarget": float(target),
                    **crossing,
                    "IndependenceProjectionAtFirstCrossing": (
                        independence_values[first_index]
                        if first_index is not None
                        else np.nan
                    ),
                    "SearchMaxN": search_max_n,
                    "SampleSizeSelected": False,
                }
            )
    return pd.DataFrame(cell_rows), pd.DataFrame(joint_rows)


def minimum_planned_slots_for_valid_target(
    target_valid_n: int,
    retention_probability: float,
    assurance_target: float,
    *,
    search_extra: int = 10000,
) -> dict[str, float | int]:
    """Find minimum planned N with binomial valid-n assurance."""

    if target_valid_n <= 0:
        raise ValueError("target_valid_n must be positive.")
    if not 0.0 < retention_probability <= 1.0:
        raise ValueError("retention_probability must be in (0, 1].")
    if not 0.0 < assurance_target < 1.0:
        raise ValueError("assurance_target must be in (0, 1).")
    for planned in range(target_valid_n, target_valid_n + search_extra + 1):
        assurance = 1.0 - binomial_cdf_logspace(
            target_valid_n - 1, planned, retention_probability
        )
        if assurance >= assurance_target:
            predecessor = (
                1.0
                - binomial_cdf_logspace(
                    target_valid_n - 1, planned - 1, retention_probability
                )
                if planned > target_valid_n
                else 0.0
            )
            return {
                "MinimumPlannedSlots": planned,
                "AssuranceAtMinimumRaw": assurance,
                "AssuranceAtPredecessorRaw": predecessor,
            }
    raise ValueError("Attrition assurance search limit was insufficient.")


def build_attrition_assurance_table(
    target_valid_n: Sequence[int] = (25, 50, 75, 100, 150, 200, 300, 500),
    retention_probabilities: Sequence[float] = (0.95, 0.90, 0.80),
    assurance_targets: Sequence[float] = (0.90, 0.95, 0.99),
) -> pd.DataFrame:
    rows = []
    for target in target_valid_n:
        for retention in retention_probabilities:
            for assurance in assurance_targets:
                result = minimum_planned_slots_for_valid_target(
                    int(target), float(retention), float(assurance)
                )
                rows.append(
                    {
                        "TargetValidN": int(target),
                        "RetentionProbabilityAssumed": float(retention),
                        "AssuranceTarget": float(assurance),
                        **result,
                        "Assumption": "independent_homogeneous_retention_sensitivity_only",
                        "PlannedRecruitmentNSelected": False,
                    }
                )
    return pd.DataFrame(rows)


def cluster_design_effect(average_cluster_size: float, intraclass_correlation: float) -> float:
    """Return the standard approximate equal-size cluster design effect."""

    if average_cluster_size < 1.0 or not 0.0 <= intraclass_correlation < 1.0:
        raise ValueError("Require average_cluster_size >= 1 and ICC in [0, 1).")
    return 1.0 + (average_cluster_size - 1.0) * intraclass_correlation


def build_cluster_sensitivity_table(
    effective_targets: Sequence[int] = (25, 50, 75, 100, 150, 200),
    average_cluster_sizes: Sequence[int] = (1, 5, 10, 20),
    intraclass_correlations: Sequence[float] = (0.0, 0.01, 0.05, 0.10),
) -> pd.DataFrame:
    rows = []
    for target in effective_targets:
        for size in average_cluster_sizes:
            for correlation in intraclass_correlations:
                effect = cluster_design_effect(float(size), float(correlation))
                rows.append(
                    {
                        "EffectiveTargetN": int(target),
                        "AverageClusterSizeAssumed": int(size),
                        "ICCAssumed": float(correlation),
                        "DesignEffect": effect,
                        "HeuristicNominalValidN": math.ceil(int(target) * effect),
                        "Interpretation": "heuristic_not_exact_wilson_or_binary_cluster_correction",
                        "SampleSizeSelected": False,
                    }
                )
    return pd.DataFrame(rows)


def build_blocked_mechanism_exposure_table(
    candidate_slots_per_language: Sequence[int] = (25, 50, 75, 100, 150, 200, 300, 500),
) -> pd.DataFrame:
    rows = []
    for slots in candidate_slots_per_language:
        schedule = build_confirmatory_assignment_schedule(int(slots))
        blocked = schedule.loc[schedule["CaseRole"].eq("blocked_primary")]
        counts = (
            blocked.groupby(["Language", "CaseId"]).size().rename("MechanismExposureN")
        )
        for (language, case_id), count in counts.items():
            rows.append(
                {
                    "CandidateSlotsPerLanguage": int(slots),
                    "Language": language,
                    "BlockedMechanism": case_id,
                    "MechanismExposureN": int(count),
                    "AllBlockedPrimaryExposureN": int(slots),
                    "Reaches25": int(count) >= 25,
                    "Reaches50": int(count) >= 50,
                    "Reaches75": int(count) >= 75,
                    "MechanismSpecificPrimaryGateActivated": False,
                    "Risk": "combined_blocked_cell_can_mask_mechanism_heterogeneity",
                }
            )
    return pd.DataFrame(rows)


def confirmatory_planning_human_gate_status() -> pd.DataFrame:
    """Return the deliberately unselected, no-human-data planning state."""

    return pd.DataFrame(
        [
            {
                "PlanningSchemaVersion": CONFIRMATORY_PLANNING_SCHEMA_VERSION,
                "HumanStudyStatus": HUMAN_STUDY_STATUS,
                "HumanParticipants": 0,
                "MinimumValidPerCellRegistered": False,
                "PlannedRecruitmentNSelected": False,
                "PilotResultAvailable": False,
                "ConfirmatoryResultAvailable": False,
                "LanguageEquivalenceAvailable": False,
                "PublicSurfaceEnabled": False,
                "NextGate": "choose_and_register_assumptions_before_any_confirmatory_response",
            }
        ]
    )


__all__ = [
    "CONFIRMATORY_PLANNING_SCHEMA_VERSION",
    "DECISION_THRESHOLD",
    "FAMILYWISE_CONFIDENCE",
    "PRIMARY_CONFIDENCE",
    "binomial_cdf_direct",
    "binomial_cdf_logspace",
    "build_attrition_assurance_table",
    "build_blocked_mechanism_exposure_table",
    "build_cell_probability_surface",
    "build_cluster_sensitivity_table",
    "build_minimum_n_sensitivity",
    "cell_gate_pass_probability",
    "cluster_design_effect",
    "confirmatory_planning_human_gate_status",
    "joint_pass_projections",
    "maximum_passing_errors_fast",
    "minimum_planned_slots_for_valid_target",
]
