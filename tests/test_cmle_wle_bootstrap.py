"""Contracts for the prospective exact-CMLE plus Warm-WLE bootstrap core."""

from __future__ import annotations

from collections import Counter
from itertools import product

import numpy as np
import pandas as pd
import pytest
from scipy.special import logsumexp

from mfrm_app.cmle import fit_cmle
from mfrm_app.cmle_wle_bootstrap import (
    FIXED_SCORE_LANE,
    JOINT_PLUGIN_LANE,
    _enumerated_fixed_score_distribution,
    _sample_fixed_score_pattern,
    _suffix_log_coefficients,
    generate_cmle_wle_bootstrap_sample,
    run_cmle_wle_bootstrap,
)


def _fit(model: str = "RSM") -> dict:
    scores = {
        "P1": [0, 0, 1, 2],
        "P2": [1, 0, 2, 2],
        "P3": [2, 1, 1, 0],
        "P4": [1, 2, 0, 2],
        "P5": [0, 2, 1, 0],
        "P6": [2, 0, 2, 1],
        "P7": [0, 0, 0, 0],
        "P8": [2, 2, 2, 2],
    }
    units = [("R1", "C1"), ("R1", "C2"), ("R2", "C1"), ("R2", "C2")]
    frame = pd.DataFrame(
        [
            (person, rater, criterion, score)
            for person, values in scores.items()
            for (rater, criterion), score in zip(units, values)
        ],
        columns=["Person", "Rater", "Criterion", "Score"],
    )
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


def test_conditional_recurrence_matches_brute_force_and_sampler_distribution():
    kernels = np.array(
        [[0.0, 0.2, -0.1], [0.0, -0.3, 0.5], [0.0, 0.7, -0.2]],
        dtype=float,
    )
    target = 3
    suffix = _suffix_log_coefficients(kernels, target)
    patterns = [
        values
        for values in product(range(3), repeat=3)
        if sum(values) == target
    ]
    brute_log_weights = [
        sum(kernels[row, value] for row, value in enumerate(values))
        for values in patterns
    ]
    assert abs(suffix[0, target] - logsumexp(brute_log_weights)) < 1e-12

    exact = _enumerated_fixed_score_distribution(kernels, target).set_index("Pattern")
    rng = np.random.default_rng(20260810)
    count = Counter(
        tuple(_sample_fixed_score_pattern(kernels, target, rng))
        for _ in range(20_000)
    )
    aligned = exact[["Probability"]].copy()
    aligned["Empirical"] = [count[pattern] / 20_000 for pattern in aligned.index]
    assert np.max(np.abs(aligned["Probability"] - aligned["Empirical"])) < 0.015


@pytest.mark.parametrize("model", ["RSM", "PCM"])
def test_fixed_score_sample_is_reproducible_and_preserves_every_person_total(model):
    fit = _fit(model)
    first = generate_cmle_wle_bootstrap_sample(
        fit, lane=FIXED_SCORE_LANE, seed=19
    )
    second = generate_cmle_wle_bootstrap_sample(
        fit, lane=FIXED_SCORE_LANE, seed=19
    )
    design = fit["design"]
    score_col = design.score_col
    assert np.array_equal(first["data"][score_col], second["data"][score_col])
    original_totals = design.data.groupby(design.person_col)[score_col].sum()
    sampled_totals = first["data"].groupby(design.person_col)[score_col].sum()
    pd.testing.assert_series_equal(original_totals, sampled_totals)
    assert int(first["audit"].iloc[0]["ChangedPersonTotals"]) == 0
    sampled = first["data"].set_index(design.person_col)
    assert sampled.loc["P7", score_col].eq(0).all()
    assert sampled.loc["P8", score_col].eq(2).all()


def test_joint_plugin_sample_is_reproducible_but_does_not_fix_totals():
    fit = _fit("PCM")
    first = generate_cmle_wle_bootstrap_sample(
        fit, lane=JOINT_PLUGIN_LANE, seed=20260810
    )
    second = generate_cmle_wle_bootstrap_sample(
        fit, lane=JOINT_PLUGIN_LANE, seed=20260810
    )
    different = generate_cmle_wle_bootstrap_sample(
        fit, lane=JOINT_PLUGIN_LANE, seed=20260811
    )
    score_col = fit["design"].score_col
    assert np.array_equal(first["data"][score_col], second["data"][score_col])
    assert not np.array_equal(first["data"][score_col], different["data"][score_col])
    assert not bool(first["audit"].iloc[0]["TotalPreservationRequired"])
    assert int(first["audit"].iloc[0]["ChangedPersonTotals"]) > 0


def test_bootstrap_runner_retains_every_attempt_and_is_seed_reproducible():
    fit = _fit("RSM")
    first = run_cmle_wle_bootstrap(
        fit,
        lane=FIXED_SCORE_LANE,
        n_replicates=4,
        seed=20260810,
    )
    second = run_cmle_wle_bootstrap(
        fit,
        lane=FIXED_SCORE_LANE,
        n_replicates=4,
        seed=20260810,
    )
    assert len(first["ledger"]) == 4
    assert len(first["generator_audit"]) == 4
    columns = [
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
    pd.testing.assert_frame_equal(first["ledger"][columns], second["ledger"][columns])
    assert first["ledger"]["Rank"].notna().all()
    assert first["ledger"]["Nullity"].notna().all()
    assert (
        first["ledger"]["PersonFitAvailable"]
        <= first["ledger"]["WLEAvailable"]
    ).all()
    assert first["ledger"].loc[
        first["ledger"]["PersonFitAvailable"], "PersonsFitReady"
    ].eq(8).all()
    assert first["ledger"].loc[
        ~first["ledger"]["PersonFitAvailable"], "PersonsFitReady"
    ].eq(0).all()
    assert first["ledger"]["PersonsFitTotal"].eq(8).all()
    required_fit_columns = {
        "Infit",
        "Outfit",
        "InfitClass",
        "OutfitClass",
        "BaselineInfit",
        "BaselineOutfit",
        "BaselineInfitClass",
        "BaselineOutfitClass",
        "InfitZoneTransition",
        "OutfitZoneTransition",
        "AnyFitZoneTransition",
        "RawDisplayMismatch",
    }
    assert required_fit_columns.issubset(first["person_draws"].columns)
    assert first["person_draws"]["PersonFitReady"].all()
    assert np.isfinite(first["person_draws"][["Infit", "Outfit"]]).all().all()
    assert first["baseline_persons"]["PersonFitReady"].all()
    summary = first["summary"].iloc[0]
    assert int(summary["AttemptedReplicates"]) == 4
    assert int(summary["PersonFitAvailableReplicates"]) == int(
        first["ledger"]["PersonFitAvailable"].sum()
    )
    assert float(summary["PersonFitAvailableShare"]) == float(
        first["ledger"]["PersonFitAvailable"].mean()
    )
    assert not bool(summary["IntervalCoverageQualified"])
    assert not bool(summary["TotalInferentialSEQualified"])
    assert not bool(summary["FitZSTDQualified"])
    assert not bool(summary["FitPValueQualified"])
    assert not bool(summary["StreamlitIntegrationAuthorized"])


def test_bootstrap_fail_closed_seed_lane_and_resource_guards():
    fit = _fit()
    with pytest.raises(ValueError, match="lane"):
        generate_cmle_wle_bootstrap_sample(fit, lane="automatic", seed=1)
    with pytest.raises(ValueError, match="boolean"):
        generate_cmle_wle_bootstrap_sample(fit, lane=FIXED_SCORE_LANE, seed=True)
    with pytest.raises(ValueError, match="maximum_replicates"):
        run_cmle_wle_bootstrap(
            fit,
            lane=FIXED_SCORE_LANE,
            n_replicates=3,
            seed=1,
            maximum_replicates=2,
        )
    with pytest.raises(ValueError, match="row-replicate work"):
        run_cmle_wle_bootstrap(
            fit,
            lane=FIXED_SCORE_LANE,
            n_replicates=3,
            seed=1,
            maximum_row_replicate_work=10,
        )

    not_ready = _fit()
    not_ready["summary"] = not_ready["summary"].copy()
    not_ready["summary"].loc[0, "InferenceReady"] = False
    with pytest.raises(ValueError, match="inference-ready"):
        generate_cmle_wle_bootstrap_sample(
            not_ready, lane=FIXED_SCORE_LANE, seed=1
        )
