"""Dependence, heterogeneity, and MNAR stress for the private CMLE gate.

This module simulates synthetic confirmatory-slot outcomes only. It does not
select a sample size, analyse human responses, or enable a public interface.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import math
from statistics import NormalDist

import numpy as np
import pandas as pd
from scipy.special import ndtri

from mfrm_app.cmle_one_click_confirmatory_gate import (
    BLOCKED_PRIMARY_CASES,
    DIRECTIONAL_RULES,
)
from mfrm_app.cmle_one_click_confirmatory_planning import (
    DECISION_THRESHOLD,
    FAMILYWISE_CONFIDENCE,
    PRIMARY_CONFIDENCE,
    cell_gate_pass_probability,
)


CONFIRMATORY_DEPENDENCE_SCHEMA_VERSION = "cmle_confirmatory_dependence_stress_v1"
LANGUAGES = ("en", "ja")
DOMAINS = tuple(str(rule["Domain"]) for rule in DIRECTIONAL_RULES)


def validate_dependence_scenario(scenario: Mapping[str, object]) -> None:
    """Validate one prospectively registered stress scenario."""

    required = {
        "scenario",
        "domain_probabilities",
        "blocked_mechanism_probabilities",
        "rho_participant",
        "rho_cluster",
        "cluster_size",
        "retention_safe",
        "retention_danger",
    }
    missing = required.difference(scenario)
    if missing:
        raise ValueError(f"Scenario lacks required fields: {sorted(missing)}")
    probabilities = tuple(float(value) for value in scenario["domain_probabilities"])
    if len(probabilities) != len(DOMAINS) or any(
        not 0.0 < value < 1.0 for value in probabilities
    ):
        raise ValueError("domain_probabilities must contain six values in (0, 1).")
    mechanism = scenario["blocked_mechanism_probabilities"]
    if mechanism is not None:
        mechanism_values = tuple(float(value) for value in mechanism)
        if len(mechanism_values) != len(BLOCKED_PRIMARY_CASES) or any(
            not 0.0 < value < 1.0 for value in mechanism_values
        ):
            raise ValueError(
                "blocked_mechanism_probabilities must be null or three values in (0, 1)."
            )
    rho_participant = float(scenario["rho_participant"])
    rho_cluster = float(scenario["rho_cluster"])
    if rho_participant < 0.0 or rho_cluster < 0.0:
        raise ValueError("Latent correlation components cannot be negative.")
    if rho_participant + rho_cluster >= 1.0:
        raise ValueError("Latent correlation components must sum to less than one.")
    cluster_size = scenario["cluster_size"]
    if not isinstance(cluster_size, int) or cluster_size <= 0:
        raise ValueError("cluster_size must be a positive integer.")
    for name in ("retention_safe", "retention_danger"):
        value = float(scenario[name])
        if not 0.0 < value <= 1.0:
            raise ValueError(f"{name} must be in (0, 1].")


def derive_condition_seed(base_seed: int, scenario_index: int, planned_n: int) -> int:
    """Return the registered deterministic seed for one scenario-by-n cell."""

    if base_seed < 0 or scenario_index < 0 or planned_n <= 0:
        raise ValueError("Invalid seed components.")
    return int(base_seed + scenario_index * 10_000 + planned_n)


def build_slot_truth_probabilities(
    scenario: Mapping[str, object], planned_n: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return slot-by-domain truth probabilities and rotated mechanism labels."""

    validate_dependence_scenario(scenario)
    if not isinstance(planned_n, int) or planned_n <= 0:
        raise ValueError("planned_n must be a positive integer.")
    base = np.asarray(scenario["domain_probabilities"], dtype=float)
    probabilities = np.broadcast_to(base, (planned_n, len(DOMAINS))).copy()
    indices = np.arange(planned_n, dtype=int) % len(BLOCKED_PRIMARY_CASES)
    mechanism_labels = np.asarray(BLOCKED_PRIMARY_CASES, dtype=object)[indices]
    mechanism = scenario["blocked_mechanism_probabilities"]
    if mechanism is not None:
        mechanism_values = np.asarray(mechanism, dtype=float)
        probabilities[:, 0] = mechanism_values[indices]
        probabilities[:, 1] = mechanism_values[indices]
    return probabilities, mechanism_labels


def analytic_valid_fraction(
    dangerous_probability: np.ndarray | float,
    retention_safe: float,
    retention_danger: float,
) -> np.ndarray:
    """Return P(valid) under outcome-dependent retention."""

    probability = np.asarray(dangerous_probability, dtype=float)
    return probability * retention_danger + (1.0 - probability) * retention_safe


def analytic_observed_danger_probability(
    dangerous_probability: np.ndarray | float,
    retention_safe: float,
    retention_danger: float,
) -> np.ndarray:
    """Return P(danger | valid) under the registered retention mechanism."""

    probability = np.asarray(dangerous_probability, dtype=float)
    denominator = analytic_valid_fraction(
        probability, retention_safe, retention_danger
    )
    return probability * retention_danger / denominator


def _wilson_upper_array(
    errors: np.ndarray,
    trials: np.ndarray,
    *,
    confidence: float,
) -> np.ndarray:
    """Vectorised one-sided Wilson upper bound with zero-valid fail closed."""

    errors_array = np.asarray(errors, dtype=float)
    trials_array = np.asarray(trials, dtype=float)
    valid = trials_array > 0.0
    safe_trials = np.where(valid, trials_array, 1.0)
    proportion = errors_array / safe_trials
    z = NormalDist().inv_cdf(confidence)
    z_squared = z * z
    denominator = 1.0 + z_squared / safe_trials
    center = proportion + z_squared / (2.0 * safe_trials)
    margin = z * np.sqrt(
        proportion * (1.0 - proportion) / safe_trials
        + z_squared / (4.0 * safe_trials * safe_trials)
    )
    upper = (center + margin) / denominator
    return np.where(valid, np.minimum(1.0, upper), 1.0)


def monte_carlo_binary_summary(values: np.ndarray) -> dict[str, float | int]:
    """Summarise a Boolean Monte Carlo endpoint with MCSE and Wilson interval."""

    array = np.asarray(values, dtype=bool).reshape(-1)
    replicates = int(array.size)
    if replicates <= 0:
        raise ValueError("At least one replicate is required.")
    successes = int(array.sum())
    rate = successes / replicates
    mcse = math.sqrt(rate * (1.0 - rate) / replicates)
    z = NormalDist().inv_cdf(0.975)
    z_squared = z * z
    denominator = 1.0 + z_squared / replicates
    center = rate + z_squared / (2.0 * replicates)
    margin = z * math.sqrt(
        rate * (1.0 - rate) / replicates
        + z_squared / (4.0 * replicates * replicates)
    )
    raw_lower = (center - margin) / denominator
    raw_upper = (center + margin) / denominator
    return {
        "Successes": successes,
        "Replicates": replicates,
        "RateRaw": rate,
        "MonteCarloSERaw": mcse,
        # The exact endpoint is zero/one when the observed rate is zero/one,
        # but subtraction can leave a tiny value on the wrong side. Preserve
        # the mathematical support and interval-containment contract.
        "WilsonLower95Raw": max(0.0, min(rate, raw_lower)),
        "WilsonUpper95Raw": min(1.0, max(rate, raw_upper)),
    }


def _rate_fields(prefix: str, values: np.ndarray) -> dict[str, float | int]:
    summary = monte_carlo_binary_summary(values)
    return {
        f"{prefix}Successes": summary["Successes"],
        f"{prefix}Replicates": summary["Replicates"],
        f"{prefix}RateRaw": summary["RateRaw"],
        f"{prefix}MonteCarloSERaw": summary["MonteCarloSERaw"],
        f"{prefix}WilsonLower95Raw": summary["WilsonLower95Raw"],
        f"{prefix}WilsonUpper95Raw": summary["WilsonUpper95Raw"],
    }


def _mean_mcse(values: np.ndarray) -> tuple[float, float]:
    array = np.asarray(values, dtype=float).reshape(-1)
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return math.nan, math.nan
    mean = float(finite.mean())
    mcse = float(finite.std(ddof=1) / math.sqrt(finite.size)) if finite.size > 1 else 0.0
    return mean, mcse


def simulate_dependence_condition(
    scenario: Mapping[str, object],
    planned_n: int,
    *,
    replicates: int,
    seed: int,
) -> dict[str, pd.DataFrame]:
    """Simulate one registered scenario-by-planned-n condition."""

    validate_dependence_scenario(scenario)
    if not isinstance(replicates, int) or replicates <= 0:
        raise ValueError("replicates must be a positive integer.")
    probabilities, mechanism_labels = build_slot_truth_probabilities(
        scenario, planned_n
    )
    rho_participant = float(scenario["rho_participant"])
    rho_cluster = float(scenario["rho_cluster"])
    cluster_size = int(scenario["cluster_size"])
    retention_safe = float(scenario["retention_safe"])
    retention_danger = float(scenario["retention_danger"])
    rng = np.random.default_rng(seed)

    participant = rng.standard_normal((replicates, len(LANGUAGES), planned_n, 1))
    group_indices = np.arange(planned_n, dtype=int) // cluster_size
    group_count = int(group_indices.max()) + 1
    cluster_base = rng.standard_normal(
        (replicates, len(LANGUAGES), group_count, 1)
    )
    cluster = cluster_base[:, :, group_indices, :]
    residual = rng.standard_normal(
        (replicates, len(LANGUAGES), planned_n, len(DOMAINS))
    )
    latent = (
        math.sqrt(rho_participant) * participant
        + math.sqrt(rho_cluster) * cluster
        + math.sqrt(1.0 - rho_participant - rho_cluster) * residual
    )
    thresholds = ndtri(probabilities)[None, None, :, :]
    dangerous = latent < thresholds
    retention_probability = np.where(
        dangerous, retention_danger, retention_safe
    )
    valid = rng.random(dangerous.shape) < retention_probability

    complete_errors = dangerous.sum(axis=2, dtype=np.int64)
    valid_counts = valid.sum(axis=2, dtype=np.int64)
    valid_errors = np.logical_and(dangerous, valid).sum(axis=2, dtype=np.int64)
    primary_upper = _wilson_upper_array(
        valid_errors, valid_counts, confidence=PRIMARY_CONFIDENCE
    )
    familywise_upper = _wilson_upper_array(
        valid_errors, valid_counts, confidence=FAMILYWISE_CONFIDENCE
    )
    primary_pass = primary_upper < DECISION_THRESHOLD
    familywise_pass = familywise_upper < DECISION_THRESHOLD
    primary_joint = primary_pass.all(axis=(1, 2))
    familywise_joint = familywise_pass.all(axis=(1, 2))

    truth_by_domain = probabilities.mean(axis=0)
    pooled_truth_risk = bool(np.any(truth_by_domain >= DECISION_THRESHOLD))
    mechanism_values = scenario["blocked_mechanism_probabilities"]
    mechanism_truth_risk = bool(
        mechanism_values is not None
        and np.any(np.asarray(mechanism_values, dtype=float) >= DECISION_THRESHOLD)
    )
    false_pooled = primary_joint & pooled_truth_risk
    false_mechanism = primary_joint & mechanism_truth_risk

    cell_rows: list[dict[str, object]] = []
    for language_index, language in enumerate(LANGUAGES):
        for domain_index, domain in enumerate(DOMAINS):
            p_slot = probabilities[:, domain_index]
            expected_valid_slot = analytic_valid_fraction(
                p_slot, retention_safe, retention_danger
            )
            expected_observed = float(
                (p_slot * retention_danger).sum() / expected_valid_slot.sum()
            )
            complete_rep = complete_errors[:, language_index, domain_index] / planned_n
            observed_rep = np.divide(
                valid_errors[:, language_index, domain_index],
                valid_counts[:, language_index, domain_index],
                out=np.full(replicates, np.nan, dtype=float),
                where=valid_counts[:, language_index, domain_index] > 0,
            )
            valid_fraction_rep = valid_counts[:, language_index, domain_index] / planned_n
            complete_mean, complete_mcse = _mean_mcse(complete_rep)
            observed_mean, observed_mcse = _mean_mcse(observed_rep)
            valid_mean, valid_mcse = _mean_mcse(valid_fraction_rep)
            primary_fields = _rate_fields(
                "PrimaryPass", primary_pass[:, language_index, domain_index]
            )
            familywise_fields = _rate_fields(
                "FamilywisePass", familywise_pass[:, language_index, domain_index]
            )
            exact_probability = math.nan
            if (
                str(scenario["scenario"]) == "independent_safe"
                and np.allclose(p_slot, p_slot[0], rtol=0.0, atol=0.0)
                and retention_safe == 1.0
                and retention_danger == 1.0
            ):
                exact_probability = float(
                    cell_gate_pass_probability(planned_n, float(p_slot[0]))[
                        "CellPassProbabilityRaw"
                    ]
                )
            cell_rows.append(
                {
                    "SchemaVersion": CONFIRMATORY_DEPENDENCE_SCHEMA_VERSION,
                    "Scenario": str(scenario["scenario"]),
                    "Seed": seed,
                    "PlannedSlotsPerLanguage": planned_n,
                    "Language": language,
                    "Domain": domain,
                    "TruthDangerousProbabilityRaw": float(p_slot.mean()),
                    "ExpectedValidFractionRaw": float(expected_valid_slot.mean()),
                    "ExpectedObservedDangerousProbabilityRaw": expected_observed,
                    "MeanCompleteDangerousProportionRaw": complete_mean,
                    "CompleteDangerousProportionMonteCarloSERaw": complete_mcse,
                    "MeanObservedDangerousProportionRaw": observed_mean,
                    "ObservedDangerousProportionMonteCarloSERaw": observed_mcse,
                    "MeanValidFractionRaw": valid_mean,
                    "ValidFractionMonteCarloSERaw": valid_mcse,
                    "MeanValidEligibleNRaw": float(
                        valid_counts[:, language_index, domain_index].mean()
                    ),
                    "MinimumValidEligibleN": int(
                        valid_counts[:, language_index, domain_index].min()
                    ),
                    "MaximumValidEligibleN": int(
                        valid_counts[:, language_index, domain_index].max()
                    ),
                    "ExactBinomialPrimaryPassProbabilityRaw": exact_probability,
                    **primary_fields,
                    **familywise_fields,
                    "MinimumValidNRegistered": False,
                    "SampleSizeSelected": False,
                }
            )
    cell_summary = pd.DataFrame(cell_rows)

    marginal_primary = cell_summary["PrimaryPassRateRaw"].to_numpy(dtype=float)
    independence_projection = float(np.prod(marginal_primary))
    union_lower = float(max(0.0, 1.0 - np.sum(1.0 - marginal_primary)))
    complete_overall = float(complete_errors.sum() / dangerous.size)
    observed_overall = float(valid_errors.sum() / valid_counts.sum())
    valid_overall = float(valid_counts.sum() / dangerous.size)
    scenario_summary = pd.DataFrame(
        [
            {
                "SchemaVersion": CONFIRMATORY_DEPENDENCE_SCHEMA_VERSION,
                "Scenario": str(scenario["scenario"]),
                "Seed": seed,
                "PlannedSlotsPerLanguage": planned_n,
                "RhoParticipantLatent": rho_participant,
                "RhoClusterLatent": rho_cluster,
                "ClusterSize": cluster_size,
                "RetentionSafe": retention_safe,
                "RetentionDanger": retention_danger,
                "AnyPooledDomainTruthAtOrAboveThreshold": pooled_truth_risk,
                "AnyBlockedMechanismTruthAtOrAboveThreshold": mechanism_truth_risk,
                "AggregateCompleteDangerousProportionRaw": complete_overall,
                "AggregateObservedDangerousProportionRaw": observed_overall,
                "AggregateValidFractionRaw": valid_overall,
                "MarginalPrimaryPassProductRaw": independence_projection,
                "MarginalPrimaryPassUnionBoundLowerRaw": union_lower,
                "MeanMarginalPrimaryPassRateRaw": float(marginal_primary.mean()),
                "MinimumMarginalPrimaryPassRateRaw": float(marginal_primary.min()),
                **_rate_fields("All12PrimaryPass", primary_joint),
                **_rate_fields("All12FamilywisePass", familywise_joint),
                **_rate_fields("FalseReassuringPooledTruth", false_pooled),
                **_rate_fields("FalseReassuringMechanismTruth", false_mechanism),
                "PrimaryJointMinusMarginalProductRaw": float(
                    primary_joint.mean() - independence_projection
                ),
                "MinimumValidNRegistered": False,
                "PlannedRecruitmentNSelected": False,
                "HumanParticipants": 0,
                "ConfirmatoryResultAvailable": False,
                "PublicSurfaceEnabled": False,
            }
        ]
    )

    mechanism_rows: list[dict[str, object]] = []
    mechanism_probability = (
        np.asarray(mechanism_values, dtype=float)
        if mechanism_values is not None
        else np.asarray(scenario["domain_probabilities"][:2], dtype=float)
    )
    for language_index, language in enumerate(LANGUAGES):
        for domain_index, domain in enumerate(DOMAINS[:2]):
            for mechanism_index, mechanism in enumerate(BLOCKED_PRIMARY_CASES):
                slot_mask = mechanism_labels == mechanism
                exposure = int(slot_mask.sum())
                complete = dangerous[:, language_index, slot_mask, domain_index]
                retained = valid[:, language_index, slot_mask, domain_index]
                observed_errors = np.logical_and(complete, retained).sum()
                observed_valid = retained.sum()
                truth_probability = (
                    float(mechanism_probability[mechanism_index])
                    if mechanism_values is not None
                    else float(probabilities[slot_mask, domain_index].mean())
                )
                expected_observed = float(
                    analytic_observed_danger_probability(
                        truth_probability, retention_safe, retention_danger
                    )
                )
                mechanism_rows.append(
                    {
                        "SchemaVersion": CONFIRMATORY_DEPENDENCE_SCHEMA_VERSION,
                        "Scenario": str(scenario["scenario"]),
                        "Seed": seed,
                        "PlannedSlotsPerLanguage": planned_n,
                        "Language": language,
                        "Domain": domain,
                        "BlockedMechanism": mechanism,
                        "ExposurePerReplicate": exposure,
                        "TruthDangerousProbabilityRaw": truth_probability,
                        "ExpectedObservedDangerousProbabilityRaw": expected_observed,
                        "AggregateCompleteDangerousProportionRaw": float(
                            complete.mean()
                        ),
                        "AggregateObservedDangerousProportionRaw": (
                            float(observed_errors / observed_valid)
                            if observed_valid > 0
                            else math.nan
                        ),
                        "MeanValidExposureRaw": float(
                            retained.sum(axis=1).mean()
                        ),
                        "TruthAtOrAboveDecisionThreshold": bool(
                            truth_probability >= DECISION_THRESHOLD
                        ),
                        "MechanismSpecificGateActivated": False,
                        "SampleSizeSelected": False,
                    }
                )
    mechanism_summary = pd.DataFrame(mechanism_rows)

    replicate_summary = pd.DataFrame(
        {
            "SchemaVersion": CONFIRMATORY_DEPENDENCE_SCHEMA_VERSION,
            "Scenario": str(scenario["scenario"]),
            "Seed": seed,
            "PlannedSlotsPerLanguage": planned_n,
            "Replicate": np.arange(1, replicates + 1, dtype=int),
            "All12PrimaryPass": primary_joint,
            "All12FamilywisePass": familywise_joint,
            "FalseReassuringPooledTruth": false_pooled,
            "FalseReassuringMechanismTruth": false_mechanism,
            "MinimumValidEligibleN": valid_counts.min(axis=(1, 2)),
            "MaximumValidEligibleN": valid_counts.max(axis=(1, 2)),
            "TotalCompleteDangerousErrors": complete_errors.sum(axis=(1, 2)),
            "TotalValidDangerousErrors": valid_errors.sum(axis=(1, 2)),
            "TotalValidEligible": valid_counts.sum(axis=(1, 2)),
        }
    )
    return {
        "scenario_summary": scenario_summary,
        "cell_summary": cell_summary,
        "mechanism_summary": mechanism_summary,
        "replicate_summary": replicate_summary,
    }


def simulate_registered_dependence_surface(
    scenarios: Sequence[Mapping[str, object]],
    planned_n_values: Sequence[int],
    *,
    replicates: int,
    base_seed: int,
) -> dict[str, pd.DataFrame]:
    """Run the complete prospectively registered synthetic stress surface."""

    if not scenarios or not planned_n_values:
        raise ValueError("Scenarios and planned_n_values cannot be empty.")
    if len({str(scenario["scenario"]) for scenario in scenarios}) != len(scenarios):
        raise ValueError("Scenario names must be unique.")
    collected: dict[str, list[pd.DataFrame]] = {
        "scenario_summary": [],
        "cell_summary": [],
        "mechanism_summary": [],
        "replicate_summary": [],
    }
    for scenario_index, scenario in enumerate(scenarios):
        validate_dependence_scenario(scenario)
        for planned_n in planned_n_values:
            seed = derive_condition_seed(base_seed, scenario_index, int(planned_n))
            result = simulate_dependence_condition(
                scenario,
                int(planned_n),
                replicates=replicates,
                seed=seed,
            )
            for name, frame in result.items():
                collected[name].append(frame)
    return {
        name: pd.concat(frames, ignore_index=True)
        for name, frames in collected.items()
    }


def dependence_stress_human_gate_status() -> pd.DataFrame:
    """Return the fail-closed human/public status for this synthetic stress."""

    return pd.DataFrame(
        [
            {
                "SchemaVersion": CONFIRMATORY_DEPENDENCE_SCHEMA_VERSION,
                "HumanParticipants": 0,
                "HumanStudyStatus": "not_started_no_human_data",
                "MinimumValidPerCellRegistered": False,
                "PlannedRecruitmentNSelected": False,
                "ConfirmatoryResultAvailable": False,
                "LanguageEquivalenceEstablished": False,
                "ClusterModelValidated": False,
                "MNARModelIdentified": False,
                "PublicSurfaceEnabled": False,
            }
        ]
    )
