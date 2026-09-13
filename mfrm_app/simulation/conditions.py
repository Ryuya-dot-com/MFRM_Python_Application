"""Public scientific-condition contracts for prospective MFRM simulation.

The Streamlit application treats research design as a first-class, portable
protocol rather than as UI state around an estimator.  This module is the
stable, standard-library-only facade for that protocol: data-generating truth,
estimator settings, keyed randomization, Monte Carlo policy, decision rules,
and failure policy remain separately versioned behind one import surface.

Internal modules may be reorganized without changing saved JSON identities,
public import paths, or pickle paths.  Response generation and fitting remain
outside this facade so importing a design protocol cannot start computation.
"""

from __future__ import annotations

# ``hashlib`` remains reachable here for the historical test/monkeypatch
# surface.  The private randomization module imports the same module object.
import hashlib

from ._condition_estimator import (
    ESTIMATOR_JMLE,
    ESTIMATOR_METHOD_CHOICES,
    ESTIMATOR_MML,
    ESTIMATOR_SPEC_VERSION,
    EstimatorSpecV1,
    estimator_spec_fingerprint,
    jmle_estimator_spec,
    mml_estimator_spec,
    normalize_estimator_spec,
)
from ._condition_evaluation import (
    CLASSIFICATION_RULE_VERSION,
    EVALUATION_FAILURE_CODES,
    EVALUATION_POLICY_VERSION,
    MONTE_CARLO_SPEC_VERSION,
    PAIR_COMPARABILITY_REQUIREMENT_V1,
    ClassificationRuleV1,
    EvaluationPolicyV1,
    MonteCarloSpecV1,
    classification_rule_fingerprint,
    evaluation_policy_fingerprint,
    monte_carlo_spec_fingerprint,
    normalize_classification_rule,
    normalize_evaluation_policy,
    normalize_monte_carlo_spec,
)
from ._condition_randomization import (
    RANDOMIZATION_ALGORITHM_V1,
    RANDOMIZATION_SHARED_RECORD_SCOPE_TWO_ARM,
    RANDOMIZATION_SPEC_VERSION,
    RANDOM_STREAM_CHOICES,
    RANDOM_STREAM_CRITERION_DIFFICULTY,
    RANDOM_STREAM_KEY_SCHEMAS,
    RANDOM_STREAM_PERSON_THETA,
    RANDOM_STREAM_RATER_CENTRAL_TENDENCY,
    RANDOM_STREAM_RATER_PERSON_LOCAL_BIAS,
    RANDOM_STREAM_RATER_SEVERITY,
    RANDOM_STREAM_RESPONSE,
    RandomizationSpecV1,
    keyed_standard_normal,
    keyed_uniform01,
    normalize_randomization_spec,
    randomization_spec_fingerprint,
)
from ._condition_scenario import (
    SIMULATION_SCENARIO_SPEC_VERSION,
    SimulationScenarioSpecV1,
    normalize_simulation_scenario_spec,
    simulation_scenario_spec_fingerprint,
)
from ._condition_shared import (
    FORMAL_MONTE_CARLO_MIN_REPLICATES,
    MAX_ABS_ADJACENT_THRESHOLD,
    MAX_CENTRAL_TENDENCY_LOG_SCALE_MEAN,
    MAX_ESTIMATOR_ITERATIONS,
    MAX_GENERATING_STANDARD_DEVIATION,
    MAX_MONTE_CARLO_REPLICATES,
    MAX_RATING_MIN_ABS,
    MAX_SIMULATION_ARTIFACT_INDEX,
    MAX_SIMULATION_CRITERION_INDEX,
    MAX_SIMULATION_PERSON_INDEX,
    MAX_SIMULATION_RATER_INDEX,
    MAX_THRESHOLD_SPAN,
    MODEL_RSM,
    SimulationConditionValidationError,
)
from ._condition_truth import (
    RATER_EFFECT_CONDITION_VERSION,
    RSM_TRUTH_SPEC_VERSION,
    RaterEffectConditionV1,
    RsmTruthSpecV1,
    normalize_rater_effect_condition,
    normalize_rsm_truth_spec,
    rater_effect_condition_fingerprint,
    rsm_truth_spec_fingerprint,
    symmetric_adjacent_thresholds,
)


# Preserve the defining public path used by historical pickles, documentation,
# and downstream introspection even though implementations now live in focused
# private modules.  Function paths are preserved as well for a uniform facade.
_PUBLIC_MODULE = __name__
_PUBLIC_TYPES = (
    ClassificationRuleV1,
    EstimatorSpecV1,
    EvaluationPolicyV1,
    MonteCarloSpecV1,
    RandomizationSpecV1,
    RaterEffectConditionV1,
    RsmTruthSpecV1,
    SimulationConditionValidationError,
    SimulationScenarioSpecV1,
)
_PUBLIC_FUNCTIONS = (
    classification_rule_fingerprint,
    estimator_spec_fingerprint,
    evaluation_policy_fingerprint,
    jmle_estimator_spec,
    keyed_standard_normal,
    keyed_uniform01,
    mml_estimator_spec,
    monte_carlo_spec_fingerprint,
    normalize_classification_rule,
    normalize_estimator_spec,
    normalize_evaluation_policy,
    normalize_monte_carlo_spec,
    normalize_randomization_spec,
    normalize_rater_effect_condition,
    normalize_rsm_truth_spec,
    normalize_simulation_scenario_spec,
    randomization_spec_fingerprint,
    rater_effect_condition_fingerprint,
    rsm_truth_spec_fingerprint,
    simulation_scenario_spec_fingerprint,
    symmetric_adjacent_thresholds,
)
for _public_object in (*_PUBLIC_TYPES, *_PUBLIC_FUNCTIONS):
    _public_object.__module__ = _PUBLIC_MODULE
del _public_object


__all__ = [
    "CLASSIFICATION_RULE_VERSION",
    "ESTIMATOR_JMLE",
    "ESTIMATOR_METHOD_CHOICES",
    "ESTIMATOR_MML",
    "ESTIMATOR_SPEC_VERSION",
    "EVALUATION_FAILURE_CODES",
    "EVALUATION_POLICY_VERSION",
    "FORMAL_MONTE_CARLO_MIN_REPLICATES",
    "MAX_ABS_ADJACENT_THRESHOLD",
    "MAX_CENTRAL_TENDENCY_LOG_SCALE_MEAN",
    "MAX_ESTIMATOR_ITERATIONS",
    "MAX_GENERATING_STANDARD_DEVIATION",
    "MAX_MONTE_CARLO_REPLICATES",
    "MAX_RATING_MIN_ABS",
    "MAX_SIMULATION_ARTIFACT_INDEX",
    "MAX_SIMULATION_CRITERION_INDEX",
    "MAX_SIMULATION_PERSON_INDEX",
    "MAX_SIMULATION_RATER_INDEX",
    "MAX_THRESHOLD_SPAN",
    "MODEL_RSM",
    "MONTE_CARLO_SPEC_VERSION",
    "RANDOMIZATION_ALGORITHM_V1",
    "PAIR_COMPARABILITY_REQUIREMENT_V1",
    "RANDOMIZATION_SHARED_RECORD_SCOPE_TWO_ARM",
    "RANDOMIZATION_SPEC_VERSION",
    "RANDOM_STREAM_CHOICES",
    "RANDOM_STREAM_KEY_SCHEMAS",
    "RANDOM_STREAM_CRITERION_DIFFICULTY",
    "RANDOM_STREAM_PERSON_THETA",
    "RANDOM_STREAM_RATER_CENTRAL_TENDENCY",
    "RANDOM_STREAM_RATER_PERSON_LOCAL_BIAS",
    "RANDOM_STREAM_RATER_SEVERITY",
    "RANDOM_STREAM_RESPONSE",
    "RATER_EFFECT_CONDITION_VERSION",
    "RSM_TRUTH_SPEC_VERSION",
    "SIMULATION_SCENARIO_SPEC_VERSION",
    "ClassificationRuleV1",
    "EstimatorSpecV1",
    "EvaluationPolicyV1",
    "MonteCarloSpecV1",
    "RandomizationSpecV1",
    "RaterEffectConditionV1",
    "RsmTruthSpecV1",
    "SimulationConditionValidationError",
    "SimulationScenarioSpecV1",
    "classification_rule_fingerprint",
    "estimator_spec_fingerprint",
    "evaluation_policy_fingerprint",
    "jmle_estimator_spec",
    "keyed_standard_normal",
    "keyed_uniform01",
    "mml_estimator_spec",
    "monte_carlo_spec_fingerprint",
    "normalize_classification_rule",
    "normalize_estimator_spec",
    "normalize_evaluation_policy",
    "normalize_monte_carlo_spec",
    "normalize_randomization_spec",
    "normalize_rater_effect_condition",
    "normalize_rsm_truth_spec",
    "normalize_simulation_scenario_spec",
    "randomization_spec_fingerprint",
    "rater_effect_condition_fingerprint",
    "rsm_truth_spec_fingerprint",
    "simulation_scenario_spec_fingerprint",
    "symmetric_adjacent_thresholds",
]
