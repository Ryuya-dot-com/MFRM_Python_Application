from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from mfrm_app.cmle_one_click_confirmatory_dependence import (
    DOMAINS,
    analytic_observed_danger_probability,
    analytic_valid_fraction,
    build_slot_truth_probabilities,
    dependence_stress_human_gate_status,
    derive_condition_seed,
    monte_carlo_binary_summary,
    simulate_dependence_condition,
    simulate_registered_dependence_surface,
    validate_dependence_scenario,
)


def scenario(
    name: str = "test",
    *,
    probabilities: tuple[float, ...] = (0.05,) * 6,
    mechanism: tuple[float, ...] | None = None,
    rho_participant: float = 0.0,
    rho_cluster: float = 0.0,
    cluster_size: int = 1,
    retention_safe: float = 1.0,
    retention_danger: float = 1.0,
) -> dict[str, object]:
    return {
        "scenario": name,
        "domain_probabilities": list(probabilities),
        "blocked_mechanism_probabilities": (
            list(mechanism) if mechanism is not None else None
        ),
        "rho_participant": rho_participant,
        "rho_cluster": rho_cluster,
        "cluster_size": cluster_size,
        "retention_safe": retention_safe,
        "retention_danger": retention_danger,
    }


def test_scenario_validation_rejects_invalid_parameters() -> None:
    validate_dependence_scenario(scenario())
    with pytest.raises(ValueError, match="six values"):
        validate_dependence_scenario(scenario(probabilities=(0.05,) * 5))
    with pytest.raises(ValueError, match="sum to less than one"):
        validate_dependence_scenario(
            scenario(rho_participant=0.7, rho_cluster=0.3)
        )
    invalid = scenario()
    invalid["retention_danger"] = 0.0
    with pytest.raises(ValueError, match="retention_danger"):
        validate_dependence_scenario(invalid)


def test_truth_matrix_rotates_and_balances_blocked_mechanisms() -> None:
    probabilities, mechanisms = build_slot_truth_probabilities(
        scenario(mechanism=(0.02, 0.04, 0.16)), 10
    )
    assert probabilities.shape == (10, len(DOMAINS))
    assert probabilities[:3, 0].tolist() == [0.02, 0.04, 0.16]
    assert probabilities[:3, 1].tolist() == [0.02, 0.04, 0.16]
    assert np.all(probabilities[:, 2:] == 0.05)
    counts = pd.Series(mechanisms).value_counts()
    assert int(counts.max() - counts.min()) <= 1


def test_outcome_dependent_retention_formula() -> None:
    probability = np.asarray([0.05, 0.10])
    valid = analytic_valid_fraction(probability, 0.95, 0.55)
    observed = analytic_observed_danger_probability(probability, 0.95, 0.55)
    assert np.allclose(valid, [0.93, 0.91])
    assert np.allclose(observed, probability * 0.55 / valid)
    assert np.all(observed < probability)


def test_seed_derivation_is_stable_and_condition_specific() -> None:
    assert derive_condition_seed(202608101, 0, 100) == 202608201
    assert derive_condition_seed(202608101, 2, 100) == 202628201
    assert derive_condition_seed(202608101, 2, 200) != derive_condition_seed(
        202608101, 2, 100
    )


def test_monte_carlo_summary_includes_precision() -> None:
    result = monte_carlo_binary_summary(np.asarray([True, False, True, True]))
    assert result["Successes"] == 3
    assert result["Replicates"] == 4
    assert result["RateRaw"] == 0.75
    assert math.isclose(result["MonteCarloSERaw"], math.sqrt(0.75 * 0.25 / 4))
    assert result["WilsonLower95Raw"] < 0.75 < result["WilsonUpper95Raw"]


def test_same_seed_returns_identical_tables() -> None:
    condition = scenario(
        rho_participant=0.2,
        rho_cluster=0.1,
        cluster_size=5,
        retention_safe=0.9,
        retention_danger=0.7,
    )
    first = simulate_dependence_condition(condition, 30, replicates=40, seed=1234)
    second = simulate_dependence_condition(condition, 30, replicates=40, seed=1234)
    for name in first:
        pd.testing.assert_frame_equal(first[name], second[name], check_exact=True)


def test_independent_baseline_matches_exact_cell_pass_probability() -> None:
    condition = scenario(name="independent_safe")
    result = simulate_dependence_condition(
        condition, 100, replicates=3000, seed=202608201
    )
    cells = result["cell_summary"]
    assert cells["ExactBinomialPrimaryPassProbabilityRaw"].notna().all()
    difference = (
        cells["PrimaryPassRateRaw"]
        - cells["ExactBinomialPrimaryPassProbabilityRaw"]
    ).abs()
    tolerance = (
        6.0
        * np.sqrt(
            cells["PrimaryPassMonteCarloSERaw"] ** 2
            + cells["ExactBinomialPrimaryPassProbabilityRaw"]
            * (1.0 - cells["ExactBinomialPrimaryPassProbabilityRaw"])
            / 3000.0
        )
        + 0.01
    )
    assert (difference <= tolerance).all()


def test_mnar_underretention_matches_analytic_distortion() -> None:
    condition = scenario(
        name="mnar_danger_underretained",
        probabilities=(0.10,) * 6,
        rho_participant=0.15,
        retention_safe=0.95,
        retention_danger=0.55,
    )
    result = simulate_dependence_condition(
        condition, 200, replicates=2000, seed=202668301
    )
    cells = result["cell_summary"]
    assert (
        cells["ExpectedObservedDangerousProbabilityRaw"]
        < cells["TruthDangerousProbabilityRaw"]
    ).all()
    difference = (
        cells["MeanObservedDangerousProportionRaw"]
        - cells["ExpectedObservedDangerousProbabilityRaw"]
    ).abs()
    tolerance = 6.0 * cells["ObservedDangerousProportionMonteCarloSERaw"] + 0.01
    assert (difference <= tolerance).all()


def test_mechanism_hotspot_is_not_relabelled_as_pooled_truth_risk() -> None:
    condition = scenario(
        name="blocked_mechanism_hotspot",
        mechanism=(0.02, 0.02, 0.16),
        rho_participant=0.15,
    )
    result = simulate_dependence_condition(
        condition, 300, replicates=300, seed=202648401
    )
    summary = result["scenario_summary"].iloc[0]
    assert not bool(summary["AnyPooledDomainTruthAtOrAboveThreshold"])
    assert bool(summary["AnyBlockedMechanismTruthAtOrAboveThreshold"])
    assert summary["FalseReassuringPooledTruthRateRaw"] == 0.0
    assert (
        summary["FalseReassuringMechanismTruthRateRaw"]
        == summary["All12PrimaryPassRateRaw"]
    )
    mechanism = result["mechanism_summary"]
    exposure_range = mechanism.groupby(["Language", "Domain"])[
        "ExposurePerReplicate"
    ].agg(lambda values: int(values.max() - values.min()))
    assert (exposure_range <= 1).all()
    assert mechanism["TruthAtOrAboveDecisionThreshold"].any()
    assert not mechanism["MechanismSpecificGateActivated"].any()


def test_familywise_joint_pass_is_subset_of_primary_joint_pass() -> None:
    result = simulate_dependence_condition(
        scenario(name="participant", rho_participant=0.35),
        200,
        replicates=500,
        seed=999,
    )
    replicate = result["replicate_summary"]
    assert (
        ~replicate["All12FamilywisePass"] | replicate["All12PrimaryPass"]
    ).all()
    summary = result["scenario_summary"].iloc[0]
    assert summary["All12FamilywisePassRateRaw"] <= summary["All12PrimaryPassRateRaw"]


def test_registered_surface_shapes_and_unique_seeds() -> None:
    scenarios = [scenario("first"), scenario("second", rho_participant=0.2)]
    result = simulate_registered_dependence_surface(
        scenarios, [25, 30], replicates=20, base_seed=1000
    )
    assert len(result["scenario_summary"]) == 4
    assert len(result["cell_summary"]) == 4 * 2 * 6
    assert len(result["mechanism_summary"]) == 4 * 2 * 2 * 3
    assert len(result["replicate_summary"]) == 4 * 20
    assert result["scenario_summary"]["Seed"].nunique() == 4


def test_human_gate_remains_closed() -> None:
    status = dependence_stress_human_gate_status().iloc[0]
    assert status["HumanParticipants"] == 0
    assert status["HumanStudyStatus"] == "not_started_no_human_data"
    assert not bool(status["MinimumValidPerCellRegistered"])
    assert not bool(status["PlannedRecruitmentNSelected"])
    assert not bool(status["ConfirmatoryResultAvailable"])
    assert not bool(status["LanguageEquivalenceEstablished"])
    assert not bool(status["ClusterModelValidated"])
    assert not bool(status["MNARModelIdentified"])
    assert not bool(status["PublicSurfaceEnabled"])
