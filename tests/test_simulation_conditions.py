"""Scientific-condition and common-random-number contracts."""

from __future__ import annotations

from dataclasses import replace
import json
import math
import os
from pathlib import Path
import pickle
import subprocess
import sys

import pytest

import mfrm_app.simulation.conditions as condition_module
import mfrm_app.simulation as simulation_module

from mfrm_app.simulation import (
    ESTIMATOR_JMLE,
    ESTIMATOR_MML,
    FORMAL_MONTE_CARLO_MIN_REPLICATES,
    MAX_CENTRAL_TENDENCY_LOG_SCALE_MEAN,
    MAX_ESTIMATOR_ITERATIONS,
    MAX_GENERATING_STANDARD_DEVIATION,
    MAX_MONTE_CARLO_REPLICATES,
    PAIR_COMPARABILITY_REQUIREMENT_V1,
    RANDOMIZATION_SHARED_RECORD_SCOPE_TWO_ARM,
    RANDOM_STREAM_CRITERION_DIFFICULTY,
    RANDOM_STREAM_KEY_SCHEMAS,
    RANDOM_STREAM_PERSON_THETA,
    RANDOM_STREAM_RATER_CENTRAL_TENDENCY,
    RANDOM_STREAM_RATER_PERSON_LOCAL_BIAS,
    RANDOM_STREAM_RATER_SEVERITY,
    RANDOM_STREAM_RESPONSE,
    ClassificationRuleV1,
    EstimatorSpecV1,
    EvaluationPolicyV1,
    MonteCarloSpecV1,
    RandomizationSpecV1,
    RaterEffectConditionV1,
    RsmTruthSpecV1,
    SimulationConditionValidationError,
    SimulationScenarioSpecV1,
    classification_rule_fingerprint,
    estimator_spec_fingerprint,
    evaluation_policy_fingerprint,
    jmle_estimator_spec,
    keyed_standard_normal,
    keyed_uniform01,
    mml_estimator_spec,
    monte_carlo_spec_fingerprint,
    normalize_classification_rule,
    normalize_monte_carlo_spec,
    normalize_rsm_truth_spec,
    normalize_simulation_scenario_spec,
    randomization_spec_fingerprint,
    rater_effect_condition_fingerprint,
    rsm_truth_spec_fingerprint,
    simulation_scenario_spec_fingerprint,
    symmetric_adjacent_thresholds,
)


def _truth(*, effects: RaterEffectConditionV1 | None = None) -> RsmTruthSpecV1:
    return RsmTruthSpecV1(
        rating_min=0,
        n_categories=5,
        adjacent_thresholds=symmetric_adjacent_thresholds(5, span=2.0),
        rater_effects=effects or RaterEffectConditionV1(),
    )


def _scenario(method: str = ESTIMATOR_JMLE) -> SimulationScenarioSpecV1:
    estimator = (
        jmle_estimator_spec(max_iterations=100, relative_tolerance=1e-5)
        if method == ESTIMATOR_JMLE
        else mml_estimator_spec(
            max_iterations=100,
            relative_tolerance=1e-5,
            quadrature_nodes=15,
            population_prior_sd=1.0,
        )
    )
    return SimulationScenarioSpecV1(
        truth=_truth(),
        estimator=estimator,
        monte_carlo=MonteCarloSpecV1(
            requested_replicates=100,
            randomization=RandomizationSpecV1(master_seed=20260722),
        ),
        classification_rule=ClassificationRuleV1(cut_scores=(0.0,)),
    )


def test_scientific_scenario_round_trip_is_strict_json_and_deterministic():
    scenario = _scenario()
    payload = scenario.to_dict()
    restored = normalize_simulation_scenario_spec(
        json.loads(json.dumps(payload, allow_nan=False))
    )

    assert restored == scenario
    assert simulation_scenario_spec_fingerprint(restored) == (
        simulation_scenario_spec_fingerprint(scenario)
    )
    assert rsm_truth_spec_fingerprint(restored.truth) == rsm_truth_spec_fingerprint(
        scenario.truth
    )
    assert estimator_spec_fingerprint(restored.estimator) == (
        estimator_spec_fingerprint(scenario.estimator)
    )
    assert monte_carlo_spec_fingerprint(restored.monte_carlo) == (
        monte_carlo_spec_fingerprint(scenario.monte_carlo)
    )


def test_public_condition_fingerprints_are_frozen_across_internal_refactors():
    scenario = _scenario()
    expected = {
        "rater_effect": "9129519ff80194ed",
        "truth": "1d6ed76dc8cd979f",
        "estimator": "18d2e24456241522",
        "randomization": "14f16020e87f5428",
        "monte_carlo": "3641c809280f4af1",
        "classification": "15824596bc8b82f0",
        "evaluation": "720596a979f27eea",
        "scenario": "4e1a451b2700c8a6",
    }
    actual = {
        "rater_effect": rater_effect_condition_fingerprint(
            scenario.truth.rater_effects
        ),
        "truth": rsm_truth_spec_fingerprint(scenario.truth),
        "estimator": estimator_spec_fingerprint(scenario.estimator),
        "randomization": randomization_spec_fingerprint(
            scenario.monte_carlo.randomization
        ),
        "monte_carlo": monte_carlo_spec_fingerprint(scenario.monte_carlo),
        "classification": classification_rule_fingerprint(
            scenario.classification_rule
        ),
        "evaluation": evaluation_policy_fingerprint(scenario.evaluation_policy),
        "scenario": simulation_scenario_spec_fingerprint(scenario),
    }

    assert actual == expected


def test_public_condition_class_paths_and_pickle_round_trip_remain_compatible():
    public_types = (
        RaterEffectConditionV1,
        RsmTruthSpecV1,
        EstimatorSpecV1,
        RandomizationSpecV1,
        MonteCarloSpecV1,
        ClassificationRuleV1,
        EvaluationPolicyV1,
        SimulationScenarioSpecV1,
        SimulationConditionValidationError,
    )
    assert all(
        public_type.__module__ == "mfrm_app.simulation.conditions"
        for public_type in public_types
    )
    scenario = _scenario()
    assert pickle.loads(pickle.dumps(scenario)) == scenario


def test_public_condition_function_paths_remain_facade_owned():
    public_functions = (
        classification_rule_fingerprint,
        estimator_spec_fingerprint,
        evaluation_policy_fingerprint,
        jmle_estimator_spec,
        keyed_standard_normal,
        keyed_uniform01,
        mml_estimator_spec,
        monte_carlo_spec_fingerprint,
        normalize_classification_rule,
        normalize_monte_carlo_spec,
        normalize_rsm_truth_spec,
        normalize_simulation_scenario_spec,
        randomization_spec_fingerprint,
        rater_effect_condition_fingerprint,
        rsm_truth_spec_fingerprint,
        simulation_scenario_spec_fingerprint,
        symmetric_adjacent_thresholds,
    )
    assert all(
        public_function.__module__ == "mfrm_app.simulation.conditions"
        for public_function in public_functions
    )
    assert len(condition_module.__all__) == len(set(condition_module.__all__))
    assert all(hasattr(condition_module, name) for name in condition_module.__all__)


def test_saved_scenario_rejects_missing_and_unknown_fields():
    payload = _scenario().to_dict()
    missing = dict(payload)
    missing.pop("truth")
    with pytest.raises(SimulationConditionValidationError, match="missing fields"):
        normalize_simulation_scenario_spec(missing)

    unknown = dict(payload, display_label="fast")
    with pytest.raises(SimulationConditionValidationError, match="unknown fields"):
        normalize_simulation_scenario_spec(unknown)


@pytest.mark.parametrize("bad", (-0.1, float("nan"), float("inf"), True))
def test_rater_effect_standard_deviations_reject_invalid_values(bad):
    with pytest.raises(SimulationConditionValidationError):
        RaterEffectConditionV1(rater_severity_sd=bad)


def test_rater_misspecification_stress_is_explicit():
    base = _truth()
    stressed = _truth(effects=RaterEffectConditionV1(
        rater_severity_sd=0.5,
        central_tendency_log_scale_mean=0.2,
        central_tendency_log_scale_sd=0.1,
        rater_person_local_bias_sd=0.3,
    ))

    assert base.is_misspecification_stress is False
    assert stressed.is_misspecification_stress is True
    assert rsm_truth_spec_fingerprint(base) != rsm_truth_spec_fingerprint(stressed)


def test_truth_distribution_centering_and_probability_machine_are_explicit():
    truth = _truth()

    assert truth.person_distribution == "normal"
    assert truth.person_centering == "population_mean_no_sample_centering"
    assert truth.criterion_difficulty_distribution == "normal"
    assert (
        truth.criterion_difficulty_sd_basis
        == "pre_centering_superpopulation_draw_sd"
    )
    assert truth.criterion_difficulty_constraint == "sample_mean_zero"
    assert truth.response_probability_normalization == "logsumexp_stable"


def test_signed_zero_does_not_split_semantically_equal_scientific_identity():
    positive = RaterEffectConditionV1(
        rater_severity_sd=0.0,
        central_tendency_log_scale_mean=0.0,
        central_tendency_log_scale_sd=0.0,
        rater_person_local_bias_sd=0.0,
    )
    negative = RaterEffectConditionV1(
        rater_severity_sd=-0.0,
        central_tendency_log_scale_mean=-0.0,
        central_tendency_log_scale_sd=-0.0,
        rater_person_local_bias_sd=-0.0,
    )
    assert positive == negative
    assert positive.to_dict() == negative.to_dict()
    assert rater_effect_condition_fingerprint(positive) == (
        rater_effect_condition_fingerprint(negative)
    )
    assert rsm_truth_spec_fingerprint(_truth(effects=positive)) == (
        rsm_truth_spec_fingerprint(_truth(effects=negative))
    )


def test_truth_threshold_contract_rejects_length_order_location_and_binary_ct():
    with pytest.raises(SimulationConditionValidationError, match="length"):
        RsmTruthSpecV1(
            rating_min=0,
            n_categories=4,
            adjacent_thresholds=(-1.0, 1.0),
        )
    with pytest.raises(SimulationConditionValidationError, match="strictly increasing"):
        RsmTruthSpecV1(
            rating_min=0,
            n_categories=4,
            adjacent_thresholds=(-1.0, 1.0, 0.0),
        )
    with pytest.raises(SimulationConditionValidationError, match="mean zero"):
        RsmTruthSpecV1(
            rating_min=0,
            n_categories=4,
            adjacent_thresholds=(-0.5, 0.5, 1.0),
        )
    with pytest.raises(SimulationConditionValidationError, match="three"):
        RsmTruthSpecV1(
            rating_min=0,
            n_categories=2,
            adjacent_thresholds=(0.0,),
            rater_effects=RaterEffectConditionV1(
                central_tendency_log_scale_mean=0.2
            ),
        )


def test_symmetric_threshold_builder_is_centered_and_strictly_increasing():
    for n_categories in (2, 3, 4, 5, 10):
        thresholds = symmetric_adjacent_thresholds(n_categories, span=2.4)
        assert len(thresholds) == n_categories - 1
        assert math.fsum(thresholds) / len(thresholds) == pytest.approx(0.0)
        assert all(right > left for left, right in zip(thresholds, thresholds[1:]))

    with pytest.raises(SimulationConditionValidationError, match="too small"):
        symmetric_adjacent_thresholds(4, span=5e-324)

    for n_categories in range(2, 101):
        thresholds = symmetric_adjacent_thresholds(n_categories, span=200.0)
        assert RsmTruthSpecV1(
            rating_min=0,
            n_categories=n_categories,
            adjacent_thresholds=thresholds,
        ).adjacent_thresholds == thresholds


def test_generating_scale_and_estimator_iteration_bounds_fail_closed():
    with pytest.raises(SimulationConditionValidationError, match="<="):
        RaterEffectConditionV1(
            central_tendency_log_scale_mean=(
                MAX_CENTRAL_TENDENCY_LOG_SCALE_MEAN + 1.0
            )
        )
    with pytest.raises(SimulationConditionValidationError, match="<="):
        RaterEffectConditionV1(
            rater_person_local_bias_sd=MAX_GENERATING_STANDARD_DEVIATION + 1.0
        )
    with pytest.raises(SimulationConditionValidationError, match="<="):
        replace(_truth(), person_sd=MAX_GENERATING_STANDARD_DEVIATION + 1.0)
    with pytest.raises(SimulationConditionValidationError, match="<="):
        jmle_estimator_spec(max_iterations=MAX_ESTIMATOR_ITERATIONS + 1)


def test_estimator_method_specific_fields_fail_closed():
    jmle = jmle_estimator_spec()
    assert jmle.method == ESTIMATOR_JMLE
    assert jmle.reported_person_uncertainty == (
        "conditional_element_information_se_other_parameters_fixed"
    )
    assert jmle.threshold_constraint == "sum_to_zero"
    with pytest.raises(SimulationConditionValidationError, match="must be None"):
        replace(jmle, mml_engine="EM")
    with pytest.raises(SimulationConditionValidationError, match="threshold_constraint"):
        replace(jmle, threshold_constraint="first_zero")

    mml = mml_estimator_spec(quadrature_nodes=15)
    assert mml.method == ESTIMATOR_MML
    with pytest.raises(SimulationConditionValidationError, match="odd"):
        replace(mml, quadrature_nodes=14)
    with pytest.raises(SimulationConditionValidationError, match="fallback"):
        replace(mml, fallback="auto")
    with pytest.raises(SimulationConditionValidationError, match="fixed population SD"):
        replace(mml, estimate_population_sd=True)

    with pytest.raises(SimulationConditionValidationError, match="method"):
        EstimatorSpecV1(method="Bayes")


def test_mml_truth_prior_must_match_in_v1():
    with pytest.raises(SimulationConditionValidationError, match="prior SD"):
        SimulationScenarioSpecV1(
            truth=_truth(),
            estimator=mml_estimator_spec(population_prior_sd=1.5),
            monte_carlo=MonteCarloSpecV1(requested_replicates=20),
        )


def test_formal_monte_carlo_and_failure_policy_have_nonnegotiable_denominators():
    with pytest.raises(SimulationConditionValidationError, match=">= 20"):
        MonteCarloSpecV1(
            requested_replicates=FORMAL_MONTE_CARLO_MIN_REPLICATES - 1
        )
    policy = EvaluationPolicyV1()
    randomization = RandomizationSpecV1()
    assert (
        randomization.shared_record_scope
        == RANDOMIZATION_SHARED_RECORD_SCOPE_TWO_ARM
    )
    assert policy.pair_inclusion == "both_arms_success_same_replicate"
    assert (
        policy.pair_comparability_requirement
        == PAIR_COMPARABILITY_REQUIREMENT_V1
    )
    assert policy.failure_denominator == "all_requested_replicates"
    assert policy.partial_run_policy == "no_final_aggregate"
    with pytest.raises(SimulationConditionValidationError, match="retry_policy"):
        replace(policy, retry_policy="retry_with_more_iterations")
    with pytest.raises(SimulationConditionValidationError, match="<="):
        replace(policy, minimum_complete_pairs=MAX_MONTE_CARLO_REPLICATES + 1)
    with pytest.raises(SimulationConditionValidationError, match="<="):
        replace(policy, minimum_complete_pairs=10**5_000)


def test_minimum_complete_pairs_cannot_exceed_requested_replicates():
    with pytest.raises(SimulationConditionValidationError, match="cannot exceed"):
        SimulationScenarioSpecV1(
            truth=_truth(),
            estimator=jmle_estimator_spec(),
            monte_carlo=MonteCarloSpecV1(requested_replicates=20),
            evaluation_policy=EvaluationPolicyV1(minimum_complete_pairs=21),
        )


def test_keyed_uniform_has_frozen_golden_and_open_interval():
    randomization = RandomizationSpecV1(master_seed=20260722)
    value = keyed_uniform01(
        randomization,
        replicate_index=1,
        stream=RANDOM_STREAM_RESPONSE,
        key_parts=(
            "study",
            "P000001",
            "P000001-A001",
            "R000001",
            "C001",
        ),
    )
    assert value.hex() == "0x1.6cf543c6f8317p-2"
    assert 0.0 < value < 1.0


def test_keyed_standard_normal_has_frozen_inputs_and_tight_numeric_golden():
    randomization = RandomizationSpecV1(master_seed=20260722)
    key_parts = ("R000001",)
    first = keyed_uniform01(
        randomization,
        replicate_index=1,
        stream=RANDOM_STREAM_RATER_SEVERITY,
        key_parts=key_parts,
        lane=0,
    )
    second = keyed_uniform01(
        randomization,
        replicate_index=1,
        stream=RANDOM_STREAM_RATER_SEVERITY,
        key_parts=key_parts,
        lane=1,
    )
    value = keyed_standard_normal(
        randomization,
        replicate_index=1,
        stream=RANDOM_STREAM_RATER_SEVERITY,
        key_parts=key_parts,
    )
    assert first.hex() == "0x1.bacb3181dcd6dp-1"
    assert second.hex() == "0x1.c634da8283e0dp-1"
    assert value == pytest.approx(0.4089700664468539, rel=0.0, abs=5e-16)


@pytest.mark.parametrize("digest_byte", (b"\x00", b"\xff"))
def test_keyed_uniform_excludes_endpoints_for_extreme_digests(
    monkeypatch,
    digest_byte,
):
    class _Digest:
        def digest(self):
            return digest_byte * 32

    monkeypatch.setattr(condition_module.hashlib, "sha256", lambda _: _Digest())
    value = keyed_uniform01(
        RandomizationSpecV1(master_seed=0),
        replicate_index=1,
        stream=RANDOM_STREAM_RATER_SEVERITY,
        key_parts=("R000001",),
    )
    assert 0.0 < value < 1.0


def test_keyed_randomness_is_order_batch_worker_and_replicate_total_invariant():
    randomization = RandomizationSpecV1(master_seed=20260722)
    keys = tuple(("study", f"P{index:06d}") for index in range(1, 6))
    forward = {
        key: keyed_uniform01(
            randomization,
            replicate_index=7,
            stream=RANDOM_STREAM_PERSON_THETA,
            key_parts=key,
        )
        for key in keys
    }
    reverse = {
        key: keyed_uniform01(
            randomization,
            replicate_index=7,
            stream=RANDOM_STREAM_PERSON_THETA,
            key_parts=key,
        )
        for key in reversed(keys)
    }
    assert forward == reverse

    short = MonteCarloSpecV1(requested_replicates=20, randomization=randomization)
    long = MonteCarloSpecV1(requested_replicates=500, randomization=randomization)
    assert randomization_spec_fingerprint(short.randomization) == (
        randomization_spec_fingerprint(long.randomization)
    )
    assert monte_carlo_spec_fingerprint(short) != monte_carlo_spec_fingerprint(long)
    assert forward[keys[0]] == keyed_uniform01(
        long.randomization,
        replicate_index=7,
        stream=RANDOM_STREAM_PERSON_THETA,
        key_parts=keys[0],
    )


def test_keyed_randomness_changes_with_owned_identity_fields():
    randomization = RandomizationSpecV1(master_seed=20260722)
    key = ("R000001",)
    baseline = keyed_uniform01(
        randomization,
        replicate_index=1,
        stream=RANDOM_STREAM_RATER_SEVERITY,
        key_parts=key,
    )
    assert baseline != keyed_uniform01(
        replace(randomization, master_seed=20260723),
        replicate_index=1,
        stream=RANDOM_STREAM_RATER_SEVERITY,
        key_parts=key,
    )
    assert baseline != keyed_uniform01(
        randomization,
        replicate_index=2,
        stream=RANDOM_STREAM_RATER_SEVERITY,
        key_parts=key,
    )
    assert math.isfinite(keyed_standard_normal(
        randomization,
        replicate_index=1,
        stream=RANDOM_STREAM_RATER_SEVERITY,
        key_parts=key,
    ))


def test_random_stream_key_arity_prevents_ambiguous_response_identity():
    with pytest.raises(SimulationConditionValidationError, match="exactly 5"):
        keyed_uniform01(
            RandomizationSpecV1(),
            replicate_index=1,
            stream=RANDOM_STREAM_RESPONSE,
            key_parts=("P000001", "R000001"),
        )


def test_random_stream_semantic_key_schemas_are_public_and_fail_closed():
    assert RANDOM_STREAM_KEY_SCHEMAS == (
        (RANDOM_STREAM_PERSON_THETA, ("person_role", "person_id")),
        (RANDOM_STREAM_RATER_SEVERITY, ("rater_id",)),
        (RANDOM_STREAM_RATER_CENTRAL_TENDENCY, ("rater_id",)),
        (
            RANDOM_STREAM_RATER_PERSON_LOCAL_BIAS,
            ("person_role", "person_id", "rater_id"),
        ),
        (RANDOM_STREAM_CRITERION_DIFFICULTY, ("criterion_id",)),
        (
            RANDOM_STREAM_RESPONSE,
            (
                "person_role",
                "person_id",
                "artifact_id",
                "rater_id",
                "criterion_id",
            ),
        ),
    )

    invalid_calls = (
        (RANDOM_STREAM_PERSON_THETA, ("P000001", "study")),
        (
            RANDOM_STREAM_RATER_PERSON_LOCAL_BIAS,
            ("study", "R000001", "P000001"),
        ),
        (
            RANDOM_STREAM_RESPONSE,
            ("study", "P000001", "P000001-A001", "C001", "R000001"),
        ),
        (RANDOM_STREAM_PERSON_THETA, ("study", "P000000")),
        (RANDOM_STREAM_PERSON_THETA, ("study", "P0000001")),
        (RANDOM_STREAM_RATER_SEVERITY, ("R000000",)),
        (RANDOM_STREAM_RATER_SEVERITY, ("R" + "1" * 5_000,)),
        (RANDOM_STREAM_CRITERION_DIFFICULTY, ("C000",)),
        (
            RANDOM_STREAM_RESPONSE,
            ("study", "P000001", "P000001-A000", "R000001", "C001"),
        ),
        (
            RANDOM_STREAM_RESPONSE,
            (
                "common_anchor",
                "CA000001",
                "CA000001-A002",
                "R000001",
                "C001",
            ),
        ),
    )
    for stream, key_parts in invalid_calls:
        with pytest.raises(SimulationConditionValidationError):
            keyed_uniform01(
                RandomizationSpecV1(),
                replicate_index=1,
                stream=stream,
                key_parts=key_parts,
            )

    with pytest.raises(SimulationConditionValidationError, match="UTF-8"):
        keyed_uniform01(
            RandomizationSpecV1(),
            replicate_index=1,
            stream=RANDOM_STREAM_RATER_SEVERITY,
            key_parts=("R00000\ud800",),
        )
    with pytest.raises(SimulationConditionValidationError, match="<="):
        keyed_uniform01(
            RandomizationSpecV1(),
            replicate_index=MAX_MONTE_CARLO_REPLICATES + 1,
            stream=RANDOM_STREAM_RATER_SEVERITY,
            key_parts=("R000001",),
        )


@pytest.mark.parametrize(
    ("normalizer", "payload", "field_name"),
    (
        (normalize_rsm_truth_spec, _truth().to_dict(), "adjacent_thresholds"),
        (
            normalize_monte_carlo_spec,
            MonteCarloSpecV1(requested_replicates=20).to_dict(),
            "quantile_probs",
        ),
        (
            normalize_classification_rule,
            ClassificationRuleV1(cut_scores=(0.0,)).to_dict(),
            "cut_scores",
        ),
    ),
)
def test_saved_array_fields_reject_non_arrays_with_domain_error(
    normalizer,
    payload,
    field_name,
):
    malformed = dict(payload)
    malformed[field_name] = None
    with pytest.raises(SimulationConditionValidationError, match=field_name):
        normalizer(malformed)


def test_response_stream_is_uniform_lane_zero_only():
    key = (
        "study",
        "P000001",
        "P000001-A001",
        "R000001",
        "C001",
    )
    with pytest.raises(SimulationConditionValidationError, match="lane 0"):
        keyed_uniform01(
            RandomizationSpecV1(),
            replicate_index=1,
            stream=RANDOM_STREAM_RESPONSE,
            key_parts=key,
            lane=1,
        )
    with pytest.raises(SimulationConditionValidationError, match="uniform-only"):
        keyed_standard_normal(
            RandomizationSpecV1(),
            replicate_index=1,
            stream=RANDOM_STREAM_RESPONSE,
            key_parts=key,
        )


def test_classification_and_failure_policy_change_metric_identity_not_truth():
    first = _scenario()
    second = replace(
        first,
        classification_rule=ClassificationRuleV1(cut_scores=(-0.5, 0.5)),
    )
    assert simulation_scenario_spec_fingerprint(first) != (
        simulation_scenario_spec_fingerprint(second)
    )
    assert rsm_truth_spec_fingerprint(first.truth) == rsm_truth_spec_fingerprint(
        second.truth
    )
    assert randomization_spec_fingerprint(first.monte_carlo.randomization) == (
        randomization_spec_fingerprint(second.monte_carlo.randomization)
    )
    assert classification_rule_fingerprint(first.classification_rule) != (
        classification_rule_fingerprint(second.classification_rule)
    )
    assert evaluation_policy_fingerprint(first.evaluation_policy) == (
        evaluation_policy_fingerprint(second.evaluation_policy)
    )


def test_fingerprint_is_stable_in_a_fresh_process():
    expected = simulation_scenario_spec_fingerprint(_scenario())
    project_root = Path(__file__).resolve().parents[1]
    code = "\n".join((
        "from mfrm_app.simulation import *",
        "truth = RsmTruthSpecV1(rating_min=0, n_categories=5, adjacent_thresholds=symmetric_adjacent_thresholds(5, span=2.0))",
        "scenario = SimulationScenarioSpecV1(truth=truth, estimator=jmle_estimator_spec(max_iterations=100, relative_tolerance=1e-5), monte_carlo=MonteCarloSpecV1(requested_replicates=100, randomization=RandomizationSpecV1(master_seed=20260722)), classification_rule=ClassificationRuleV1(cut_scores=(0.0,)))",
        "print(simulation_scenario_spec_fingerprint(scenario))",
    ))
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=project_root,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    assert completed.stdout.strip() == expected


def test_conditions_import_remains_free_of_heavy_and_ui_dependencies():
    project_root = Path(__file__).resolve().parents[1]
    code = "\n".join((
        "import sys",
        "import mfrm_app.simulation.conditions",
        "for name in ('streamlit', 'pandas', 'numpy', 'networkx', 'scipy', 'plotly'):",
        "    assert name not in sys.modules, name",
    ))
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    subprocess.run(
        [sys.executable, "-c", code],
        cwd=project_root,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )


def test_simulation_package_public_exports_are_unique_and_resolvable():
    assert len(simulation_module.__all__) == len(set(simulation_module.__all__))
    assert all(hasattr(simulation_module, name) for name in simulation_module.__all__)
    for name in (
        "BASE_JMLE_BLOCKED",
        "BASE_MML_PRIOR_CONDITIONAL",
        "EXECUTION_NOT_STARTED",
        "PAIR_COMPARABILITY_REQUIREMENT_V1",
        "RANDOM_STREAM_KEY_SCHEMAS",
        "RANDOMIZATION_SHARED_RECORD_SCOPE_TWO_ARM",
        "SIMULATION_SCENARIO_SPEC_VERSION",
    ):
        assert name in simulation_module.__all__
