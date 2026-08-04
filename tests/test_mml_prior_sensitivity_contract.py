from __future__ import annotations

import copy

import pytest

from mfrm_app import evidence
from mfrm_app import mml_prior_sensitivity as contract


THRESHOLDS = {
    "max_abs_measure_shift": 0.25,
    "measure_rmse": 0.10,
    "rank_correlation": 0.98,
    "population_coefficient_shift": 0.10,
}


def _identities():
    baseline = evidence.build_analysis_identity(
        input_data_fingerprint="data_sha256",
        resolved_settings={"method": "MML", "population_prior_sd": 1.0},
        analysis_type="mfrm.mml.fit",
        engine_version="test",
    )
    variants = {}
    for name, sd in (("lower", 0.75), ("upper", 1.25)):
        variants[name] = evidence.build_analysis_identity(
            input_data_fingerprint=baseline.input_data_fingerprint,
            resolved_settings={"method": "MML", "population_prior_sd": sd},
            analysis_type=baseline.analysis_type,
            engine_version=baseline.engine_version,
            variant_id=name,
            baseline=baseline,
        )
    return baseline, variants


def _outcome(**metric_overrides):
    metrics = {
        "max_abs_measure_shift": 0.10,
        "measure_rmse": 0.05,
        "rank_correlation": 0.99,
    }
    metrics.update(metric_overrides)
    return {
        "run_ok": True,
        "converged": True,
        "comparable": True,
        "metrics": metrics,
        "summary": "comparison completed",
    }


def _bundle(outcomes, *, latent_regression=False):
    baseline, variants = _identities()
    plan = contract.build_sensitivity_plan(
        baseline,
        variants,
        thresholds=THRESHOLDS,
        latent_regression=latent_regression,
    )
    bundle = contract.build_contract_bundle(
        baseline=baseline,
        variants=variants,
        plan=plan,
        outcomes=outcomes,
        baseline_converged=True,
        latent_regression=latent_regression,
    )
    return baseline, variants, plan, bundle


def test_mml_contract_is_stable_at_or_inside_exact_thresholds_and_round_trips():
    outcomes = {
        "lower": _outcome(max_abs_measure_shift=0.25, measure_rmse=0.10),
        "upper": _outcome(rank_correlation=0.98),
    }
    baseline, variants, plan, bundle = _bundle(outcomes)

    assert bundle["schema_version"] == contract.MML_PRIOR_SENSITIVITY_BUNDLE_SCHEMA_VERSION
    assert bundle["sensitivity_decision"]["stability_state"] == "STABLE"
    assert "population_coefficient_shift" not in {
        criterion["metric"] for criterion in plan.decision_rule["criteria"]
    }
    restored_plan = evidence.SensitivityPlan.from_payload(bundle["sensitivity_plan"])
    restored_records = tuple(
        evidence.SensitivityRecord.from_payload(payload)
        for payload in bundle["sensitivity_records"]
    )
    restored = evidence.SensitivityDecision.from_payload(
        bundle["sensitivity_decision"],
        plan=restored_plan,
        records=restored_records,
    )
    assert restored.stability_state is evidence.StabilityState.STABLE
    assert restored_plan.baseline_analysis_id == baseline.analysis_id
    assert set(restored_plan.required_variant_ids) == set(variants)


@pytest.mark.parametrize(
    "changed_metric",
    [
        {"max_abs_measure_shift": 0.250001},
        {"measure_rmse": 0.100001},
        {"rank_correlation": 0.979999},
    ],
)
def test_mml_contract_marks_any_threshold_crossing_sensitive(changed_metric):
    _, _, _, bundle = _bundle(
        {"lower": _outcome(**changed_metric), "upper": _outcome()}
    )
    assert bundle["sensitivity_decision"]["stability_state"] == "SENSITIVE"


def test_mml_contract_fails_closed_for_partial_and_total_run_failure():
    failed = {
        "run_ok": False,
        "converged": False,
        "comparable": False,
        "metrics": {},
        "summary": "refit failed",
    }
    _, _, _, partial = _bundle({"lower": _outcome(), "upper": failed})
    _, _, _, total = _bundle({"lower": failed, "upper": failed})

    assert partial["sensitivity_decision"]["stability_state"] == "CONDITIONALLY_STABLE"
    assert total["sensitivity_decision"]["stability_state"] == "NOT_ASSESSED"
    assert all(
        record["computation_state"] == "HOLD"
        for record in total["sensitivity_records"]
    )


def test_mml_contract_requires_population_metric_only_for_latent_regression():
    outcomes = {"lower": _outcome(), "upper": _outcome()}
    _, _, plan, bundle = _bundle(outcomes, latent_regression=True)

    assert "population_coefficient_shift" in {
        criterion["metric"] for criterion in plan.decision_rule["criteria"]
    }
    assert bundle["sensitivity_decision"]["stability_state"] == "NOT_ASSESSED"
    assert bundle["sensitivity_decision"]["reason_code"] == (
        "sensitivity.rule_not_evaluable"
    )


@pytest.mark.parametrize("field", ["run_ok", "converged", "comparable"])
def test_mml_contract_rejects_string_boolean_outcome_flags(field):
    outcomes = {"lower": _outcome(), "upper": _outcome()}
    outcomes["lower"][field] = "false"

    with pytest.raises(evidence.ContractValidationError, match="must be a boolean"):
        _bundle(outcomes)


def test_saved_mml_contract_rejects_a_mutated_decision_with_same_baseline():
    _, _, _, bundle = _bundle(
        {"lower": _outcome(), "upper": _outcome()}
    )
    corrupted = copy.deepcopy(bundle)
    corrupted["sensitivity_decision"]["summary"] = "tampered summary"

    with pytest.raises(
        evidence.ContractValidationError,
        match="SensitivityDecisionID|does not reproduce",
    ):
        contract.validate_contract_bundle(corrupted)
