from __future__ import annotations

import numpy as np
import pandas as pd
import pandas.testing as pdt

from mfrm_app.assignment_mechanism import (
    compare_chain_to_exact_oracle,
    enumerate_degree_conditioned_assignments,
    materialize_exchangeable_assignment_design,
    sample_degree_conditioned_assignment_chain,
    sample_degree_conditioned_assignment_dense_dp,
    sample_degree_conditioned_assignment_dp,
    sample_degree_conditioned_assignment_heatbath_chain,
)


PEOPLE = [f"p{index}" for index in range(6)]
RATERS = [f"r{index}" for index in range(3)]
PERSON_DEGREES = {person: 2 for person in PEOPLE}
RATER_DEGREES = {rater: 4 for rater in RATERS}
PERSON_COORDINATES = {person: float(index) for index, person in enumerate(PEOPLE)}
RATER_COORDINATES = {rater: float(index) for index, rater in enumerate(RATERS)}
INITIAL_EDGES = {
    ("p0", "r0"), ("p0", "r1"),
    ("p1", "r0"), ("p1", "r2"),
    ("p2", "r0"), ("p2", "r1"),
    ("p3", "r1"), ("p3", "r2"),
    ("p4", "r0"), ("p4", "r2"),
    ("p5", "r1"), ("p5", "r2"),
}


def _exchangeable_rows() -> pd.DataFrame:
    rows = []
    for person, rater in sorted(INITIAL_EDGES):
        for task in ("t1", "t2"):
            for criterion in ("c1", "c2"):
                rows.append(
                    {
                        "Person": person,
                        "Rater": rater,
                        "Task": task,
                        "Criterion": criterion,
                        "Score": (int(person[1:]) + int(rater[1:])) % 4,
                        "Weight": 1.0,
                    }
                )
    return pd.DataFrame(rows)


def _oracle(gamma: float):
    return enumerate_degree_conditioned_assignments(
        person_degrees=PERSON_DEGREES,
        rater_degrees=RATER_DEGREES,
        person_coordinates=PERSON_COORDINATES,
        rater_coordinates=RATER_COORDINATES,
        gamma=gamma,
        require_connected=True,
    )


def test_small_fixed_degree_state_space_is_exact_and_uniform_at_gamma_zero():
    oracle = _oracle(0.0)
    states = oracle["states"]
    summary = oracle["summary"].iloc[0]

    assert oracle["available"] is True
    assert len(states) == 90
    assert states["Connected"].all()
    assert np.isclose(states["Probability"].sum(), 1.0)
    assert np.allclose(states["Probability"], 1.0 / 90.0)
    assert summary["States"] == 90
    assert np.isclose(summary["ProbabilitySum"], 1.0)


def test_positive_gamma_increases_exact_ability_severity_alignment():
    neutral = _oracle(0.0)["summary"].iloc[0]
    aligned = _oracle(0.8)["summary"].iloc[0]
    anti = _oracle(-0.8)["summary"].iloc[0]

    assert aligned["ExpectedStatistic"] > neutral["ExpectedStatistic"]
    assert neutral["ExpectedStatistic"] > anti["ExpectedStatistic"]
    assert aligned["ExpectedAssignmentCorrelation"] > neutral["ExpectedAssignmentCorrelation"]
    assert neutral["ExpectedAssignmentCorrelation"] > anti["ExpectedAssignmentCorrelation"]


def test_metropolis_chain_preserves_margins_connectivity_and_is_reproducible():
    kwargs = dict(
        initial_edges=INITIAL_EDGES,
        person_coordinates=PERSON_COORDINATES,
        rater_coordinates=RATER_COORDINATES,
        gamma=0.8,
        n_samples=2_000,
        burnin=1_000,
        thin=5,
        seed=260811,
        require_connected=True,
    )
    first = sample_degree_conditioned_assignment_chain(**kwargs)
    second = sample_degree_conditioned_assignment_chain(**kwargs)

    assert first["available"] is True
    assert first["diagnostics"].iloc[0]["DegreeMarginsPreserved"]
    assert first["diagnostics"].iloc[0]["FinalOverlapComponents"] == 1
    assert first["diagnostics"].iloc[0]["MaximumStatisticUpdateResidual"] < 1e-10
    assert first["diagnostics"].iloc[0]["UniqueSampledStates"] > 30
    assert first["samples"].equals(second["samples"])
    assert first["final_edges"] == second["final_edges"]
    assert {person for person, _ in first["final_edges"]} == set(PEOPLE)
    assert {rater for _, rater in first["final_edges"]} == set(RATERS)


def test_seeded_chain_recovers_exact_state_distribution_and_moment():
    gamma = 0.8
    oracle = _oracle(gamma)
    chain = sample_degree_conditioned_assignment_chain(
        initial_edges=INITIAL_EDGES,
        person_coordinates=PERSON_COORDINATES,
        rater_coordinates=RATER_COORDINATES,
        gamma=gamma,
        n_samples=30_000,
        burnin=5_000,
        thin=3,
        seed=260812,
        require_connected=True,
    )
    comparison = compare_chain_to_exact_oracle(chain, oracle)
    summary = comparison["summary"].iloc[0]

    assert summary["ChainProbabilityOutsideOracle"] == 0.0
    assert summary["SampledStates"] >= 75
    assert summary["TotalVariationDistance"] < 0.06
    assert abs(summary["MeanStatisticError"]) < 0.06


def test_heatbath_chain_recovers_exact_oracle_with_better_block_movement():
    gamma = -0.8
    oracle = _oracle(gamma)
    chain = sample_degree_conditioned_assignment_heatbath_chain(
        initial_edges=INITIAL_EDGES,
        person_coordinates=PERSON_COORDINATES,
        rater_coordinates=RATER_COORDINATES,
        gamma=gamma,
        n_samples=12_000,
        burnin=2_000,
        thin=2,
        seed=260814,
        require_connected=True,
    )
    comparison = compare_chain_to_exact_oracle(chain, oracle)
    summary = comparison["summary"].iloc[0]
    diagnostics = chain["diagnostics"].iloc[0]

    assert chain["available"] is True
    assert diagnostics["DegreeMarginsPreserved"]
    assert diagnostics["FinalOverlapComponents"] == 1
    assert diagnostics["MaximumStatisticUpdateResidual"] < 1e-10
    assert diagnostics["MovementRateAllUpdates"] > 0.05
    assert summary["ChainProbabilityOutsideOracle"] == 0.0
    assert summary["SampledStates"] >= 75
    assert summary["TotalVariationDistance"] < 0.07
    assert abs(summary["MeanStatisticError"]) < 0.07


def test_dp_sampler_produces_independent_exact_oracle_draws():
    gamma = 0.8
    oracle = _oracle(gamma)
    sampled = sample_degree_conditioned_assignment_dp(
        person_degrees=PERSON_DEGREES,
        rater_degrees=RATER_DEGREES,
        person_coordinates=PERSON_COORDINATES,
        rater_coordinates=RATER_COORDINATES,
        gamma=gamma,
        n_samples=12_000,
        seed=260815,
        require_connected=True,
    )
    comparison = compare_chain_to_exact_oracle(sampled, oracle)
    summary = comparison["summary"].iloc[0]
    diagnostics = sampled["diagnostics"].iloc[0]

    assert sampled["available"] is True
    assert diagnostics["ExactIndependentSamples"]
    assert diagnostics["ConnectedRejections"] == 0
    assert diagnostics["DPStatesEvaluated"] < diagnostics["DPStateCap"]
    assert diagnostics["DegreeMarginsPreserved"]
    assert diagnostics["FinalOverlapComponents"] == 1
    assert diagnostics["MaximumStatisticUpdateResidual"] < 1e-10
    assert summary["ChainProbabilityOutsideOracle"] == 0.0
    assert summary["SampledStates"] >= 75
    assert summary["TotalVariationDistance"] < 0.06
    assert abs(summary["MeanStatisticError"]) < 0.06


def test_dense_dp_sampler_matches_recursive_partition_and_exact_oracle():
    gamma = -0.8
    oracle = _oracle(gamma)
    recursive = sample_degree_conditioned_assignment_dp(
        person_degrees=PERSON_DEGREES,
        rater_degrees=RATER_DEGREES,
        person_coordinates=PERSON_COORDINATES,
        rater_coordinates=RATER_COORDINATES,
        gamma=gamma,
        n_samples=1,
        seed=260816,
        require_connected=True,
    )
    dense = sample_degree_conditioned_assignment_dense_dp(
        person_degrees=PERSON_DEGREES,
        rater_degrees=RATER_DEGREES,
        person_coordinates=PERSON_COORDINATES,
        rater_coordinates=RATER_COORDINATES,
        gamma=gamma,
        n_samples=12_000,
        seed=260817,
        require_connected=True,
    )
    comparison = compare_chain_to_exact_oracle(dense, oracle)
    summary = comparison["summary"].iloc[0]
    diagnostics = dense["diagnostics"].iloc[0]

    assert dense["available"] is True
    assert np.isclose(
        diagnostics["DPLogNormalizerUnconditioned"],
        recursive["diagnostics"].iloc[0]["DPLogNormalizerUnconditioned"],
        atol=1e-12,
    )
    assert diagnostics["MaximumConditionalLogResidual"] < 1e-10
    assert diagnostics["MaximumStatisticUpdateResidual"] < 1e-10
    assert diagnostics["ConnectedRejections"] == 0
    assert summary["ChainProbabilityOutsideOracle"] == 0.0
    assert summary["SampledStates"] >= 75
    assert summary["TotalVariationDistance"] < 0.06
    assert abs(summary["MeanStatisticError"]) < 0.06


def test_dense_dp_fails_closed_when_registered_cell_cap_is_too_small():
    with np.testing.assert_raises(ValueError):
        sample_degree_conditioned_assignment_dense_dp(
            person_degrees=PERSON_DEGREES,
            rater_degrees=RATER_DEGREES,
            person_coordinates=PERSON_COORDINATES,
            rater_coordinates=RATER_COORDINATES,
            gamma=0.0,
            n_samples=1,
            max_dense_dp_cells=10,
        )


def test_invalid_degree_contract_and_disconnected_start_fail_closed():
    with np.testing.assert_raises(ValueError):
        enumerate_degree_conditioned_assignments(
            person_degrees={"p0": 1, "p1": 1},
            rater_degrees={"r0": 1, "r1": 2},
            person_coordinates={"p0": 0.0, "p1": 1.0},
            rater_coordinates={"r0": 0.0, "r1": 1.0},
            gamma=0.0,
        )

    disconnected = {("p0", "r0"), ("p1", "r0"), ("p2", "r1"), ("p3", "r1")}
    with np.testing.assert_raises(ValueError):
        sample_degree_conditioned_assignment_chain(
            initial_edges=disconnected,
            person_coordinates={f"p{i}": float(i) for i in range(4)},
            rater_coordinates={"r0": 0.0, "r1": 1.0},
            gamma=0.0,
            n_samples=10,
            require_connected=True,
        )


def test_sampled_graph_materialization_is_score_free_and_preserves_all_exposures():
    target = set(INITIAL_EDGES)
    target.remove(("p0", "r0"))
    target.remove(("p5", "r2"))
    target.add(("p0", "r2"))
    target.add(("p5", "r0"))
    source = _exchangeable_rows()
    bundle = materialize_exchangeable_assignment_design(
        source,
        target,
        facet_cols=["Rater", "Task", "Criterion"],
        context_cols=["Task", "Criterion"],
        rater_mapping_confirmed=True,
    )

    assert bundle["available"] is True
    assert bundle["invariants"]["Passed"].all()
    assert bundle["changed_edges"] == 2
    assert "Score" not in bundle["design"].columns
    assert bundle["design"]["_SourceRow"].tolist() == list(range(len(source)))
    materialized = set(
        bundle["design"][["Person", "Rater"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    )
    assert materialized == target
    pdt.assert_series_equal(
        source.groupby(["Rater", "Task", "Criterion"]).size().sort_index(),
        bundle["design"].groupby(["Rater", "Task", "Criterion"]).size().sort_index(),
    )


def test_materialization_does_not_depend_on_scores_and_fails_closed_on_bad_margin():
    target = set(INITIAL_EDGES)
    target.remove(("p0", "r0"))
    target.remove(("p5", "r2"))
    target.add(("p0", "r2"))
    target.add(("p5", "r0"))
    source = _exchangeable_rows()
    changed_scores = source.copy()
    changed_scores["Score"] = np.arange(len(changed_scores)) * 17
    kwargs = dict(
        target_edges=target,
        facet_cols=["Rater", "Task", "Criterion"],
        context_cols=["Task", "Criterion"],
        rater_mapping_confirmed=True,
    )
    first = materialize_exchangeable_assignment_design(source, **kwargs)
    second = materialize_exchangeable_assignment_design(changed_scores, **kwargs)

    pdt.assert_frame_equal(first["design"], second["design"])
    pdt.assert_frame_equal(first["assignment_map"], second["assignment_map"])

    invalid = set(target)
    invalid.remove(("p0", "r1"))
    failed = materialize_exchangeable_assignment_design(source, invalid, **{k: v for k, v in kwargs.items() if k != "target_edges"})
    assert failed["available"] is False
    assert not failed["invariants"]["Passed"].all()
    assert failed["design"].empty
