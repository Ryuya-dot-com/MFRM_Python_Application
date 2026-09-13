"""Contracts for the known-truth Person MnSq pilot core."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mfrm_app.fit_threshold_operating_characteristics import (
    FIT_RULES,
    PERSON_MEASURE_TRUE,
    PERSON_MEASURE_WLE,
    fit_rule_flags,
    replicate_fit_rule_rates,
    rounded_rule_disagreements,
    rsm_category_intercepts,
    score_known_truth_person_fit,
    simulate_known_truth_person_fit_sample,
    summarize_replicate_fit_rates,
)


def _condition(
    mechanism: str = "clean",
    *,
    observations: int = 4,
) -> dict[str, object]:
    affected = 0.0 if mechanism == "clean" else 0.25
    probability = 0.0 if mechanism == "clean" else 0.5
    return {
        "ConditionId": f"test_{mechanism}",
        "Categories": 3,
        "ObservationsPerPerson": observations,
        "Mechanism": mechanism,
        "AffectedShare": affected,
        "MechanismProbability": probability,
    }


def _sample(condition: dict[str, object], seed: int = 314) -> dict[str, object]:
    return simulate_known_truth_person_fit_sample(
        condition,
        seed=seed,
        persons=12,
        rater_effects=[-0.25, 0.25],
        criterion_effects=[-0.15, 0.15],
        category_thresholds=[-0.6, 0.6],
    )


def test_rsm_category_intercepts_follow_cumulative_threshold_contract():
    observed = rsm_category_intercepts([0.0, 0.5], [-0.6, 0.6])
    expected = np.array([[0.0, 0.6, 0.0], [0.0, 0.1, -1.0]])
    assert np.allclose(observed, expected)
    with pytest.raises(ValueError, match="increasing"):
        rsm_category_intercepts([0.0], [0.2, -0.2])


def test_common_seed_preserves_baseline_and_random_changes_only_affected_rows():
    clean = _sample(_condition("clean"))
    random = _sample(_condition("uniform_random_replacement"))
    clean_rows = clean["responses"]
    random_rows = random["responses"]

    assert clean_rows["TrueTheta"].equals(random_rows["TrueTheta"])
    assert clean_rows["BaselineCategory"].equals(random_rows["BaselineCategory"])
    unaffected = ~random_rows["Affected"]
    assert random_rows.loc[unaffected, "ObservedCategory"].equals(
        random_rows.loc[unaffected, "BaselineCategory"]
    )
    assert not random_rows.loc[~random_rows["Affected"], "MechanismApplied"].any()
    assert random["audit"]["AffectedPersons"] == 3


def test_sparse_administration_has_distinct_frozen_count_per_person():
    sample = _sample(_condition("clean", observations=2))
    rows = sample["responses"]
    summary = rows.groupby("Person")["UnitIndex"].agg(["size", "nunique"])
    assert summary["size"].eq(2).all()
    assert summary["nunique"].eq(2).all()


def test_threshold_heterogeneity_changes_only_affected_true_kernels():
    condition = _condition("affected_person_pcm_threshold_heterogeneity")
    sample = simulate_known_truth_person_fit_sample(
        condition,
        seed=19,
        persons=12,
        rater_effects=[-0.25, 0.25],
        criterion_effects=[-0.15, 0.15],
        category_thresholds=[-0.6, 0.6],
        threshold_deviations={
            "C01": [-0.2, 0.2],
            "C02": [0.2, -0.2],
        },
    )
    rows = sample["responses"]
    differences = np.any(
        sample["working_intercepts"] != sample["true_intercepts"], axis=1
    )
    assert not differences[~rows["Affected"].to_numpy()].any()
    assert differences[rows["Affected"].to_numpy()].all()


def test_scoring_returns_wle_and_true_theta_fit_for_every_person():
    persons = score_known_truth_person_fit(_sample(_condition("clean")))
    assert len(persons) == 24
    assert set(persons["PersonMeasureSource"]) == {
        PERSON_MEASURE_WLE,
        PERSON_MEASURE_TRUE,
    }
    assert persons["PersonFitReady"].all()
    assert np.isfinite(persons[["Infit", "Outfit", "PersonMeasure"]]).all().all()
    assert persons["DecisionInput"].eq("finite_unrounded_mnsq").all()


def test_fit_rule_endpoints_are_raw_and_exact():
    values = np.array(
        [
            0.5,
            np.nextafter(0.5, -np.inf),
            1.5,
            np.nextafter(1.5, np.inf),
            2.0,
            np.nextafter(2.0, np.inf),
        ]
    )
    flags = fit_rule_flags(values, np.ones(len(values)))
    assert flags["either_upper"].tolist() == [False, False, False, True, True, True]
    assert flags["either_distorting"].tolist() == [False, False, False, False, False, True]
    assert flags["either_overfit"].tolist() == [False, True, False, False, False, False]
    assert flags["either_nonacceptable"].tolist() == [False, True, False, True, True, True]


def test_replicate_rates_keep_truth_groups_and_rules_separate():
    frame = pd.DataFrame(
        {
            "ConditionId": ["C"] * 4,
            "Replicate": [1] * 4,
            "PersonMeasureSource": [PERSON_MEASURE_WLE] * 4,
            "TruthGroup": ["clean", "clean", "affected", "affected"],
            "Infit": [1.0, 1.6, 1.0, 2.1],
            "Outfit": [1.0] * 4,
        }
    )
    rates = replicate_fit_rule_rates(
        frame, threshold_triplets=[(0.5, 1.5, 2.0)]
    )
    assert len(rates) == 2 * len(FIT_RULES)
    upper = rates.loc[rates["Rule"].eq("either_upper")]
    assert upper["PersonsEligible"].eq(2).all()
    assert upper["Flagged"].eq(1).all()
    assert upper["Rate"].eq(0.5).all()
    assert rates["PersonRowsClustered"].all()


def test_operating_summary_uses_replicate_mean_not_naive_pooled_rate():
    frame = pd.DataFrame(
        {
            "ConditionId": ["C", "C"],
            "PersonMeasureSource": [PERSON_MEASURE_WLE] * 2,
            "TruthGroup": ["clean"] * 2,
            "OverfitUpper": [0.5] * 2,
            "AcceptableUpper": [1.5] * 2,
            "NoisyUpper": [2.0] * 2,
            "CanonicalThresholds": [True] * 2,
            "Rule": ["either_upper"] * 2,
            "Replicate": [1, 2],
            "PersonsAttempted": [2, 18],
            "PersonsEligible": [2, 18],
            "PersonsUnavailable": [0, 0],
            "Flagged": [0, 9],
            "Rate": [0.0, 0.5],
        }
    )
    summary = summarize_replicate_fit_rates(frame).iloc[0]
    assert summary["MeanReplicateRate"] == pytest.approx(0.25)
    assert summary["ReplicateRateMCSE"] == pytest.approx(0.25)
    assert summary["PooledPersonRateDescriptive"] == pytest.approx(0.45)
    assert not bool(summary["IndependentBinomialWilsonIntervalAuthorized"])


def test_rounding_audit_never_replaces_raw_rule():
    frame = pd.DataFrame(
        {
            "ConditionId": ["C"],
            "Replicate": [1],
            "PersonMeasureSource": [PERSON_MEASURE_WLE],
            "TruthGroup": ["affected"],
            "Infit": [1.5004],
            "Outfit": [1.0],
        }
    )
    result = rounded_rule_disagreements(frame, decimals=[3, 4])
    upper = result.loc[result["Rule"].eq("either_upper")].set_index(
        "DisplayDecimals"
    )
    assert int(upper.loc[3, "RawFlagged"]) == 1
    assert int(upper.loc[3, "RoundedFlagged"]) == 0
    assert int(upper.loc[3, "RuleDisagreements"]) == 1
    assert int(upper.loc[4, "RuleDisagreements"]) == 0
    assert result["RawDecisionRetained"].all()
