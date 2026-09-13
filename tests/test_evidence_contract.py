from __future__ import annotations

import json
from pathlib import Path

import pytest

from mfrm_app import evidence


def _baseline(**settings) -> evidence.AnalysisIdentity:
    return evidence.build_analysis_identity(
        input_data_fingerprint="data_0123456789abcdef",
        resolved_settings={"model": "RSM", "method": "JMLE", **settings},
        analysis_type="mfrm.fit",
        engine_version="0.3.0-beta",
        seed=20260724,
    )


def _variant(
    baseline: evidence.AnalysisIdentity,
    variant_id: str,
    **settings,
) -> evidence.AnalysisIdentity:
    return evidence.build_analysis_identity(
        input_data_fingerprint=baseline.input_data_fingerprint,
        resolved_settings={"model": "RSM", "method": "JMLE", **settings},
        analysis_type="mfrm.fit",
        engine_version=baseline.engine_version,
        variant_id=variant_id,
        baseline=baseline,
        seed=baseline.seed,
    )


def _plan(
    baseline: evidence.AnalysisIdentity,
    *variant_ids: str,
    conclusion_id: str = "rater_ordering",
) -> evidence.SensitivityPlan:
    variants = [
        _variant(baseline, variant_id, variant=variant_id)
        for variant_id in variant_ids
    ]
    return evidence.make_sensitivity_plan(
        baseline,
        conclusion_id=conclusion_id,
        required_variants=variants,
        decision_rule=evidence.make_sensitivity_rule(
            criteria=(
                {"metric": "rank_correlation", "operator": "lt", "value": 0.98},
                {"metric": "ordering_changed", "operator": "eq", "value": True},
            )
        ),
    )


def _available_evidence(
    identity: evidence.AnalysisIdentity,
    *,
    key: str = "fit.convergence",
    required: bool = True,
    stability: evidence.StabilityState = evidence.StabilityState.NOT_ASSESSED,
) -> evidence.EvidenceRecord:
    return evidence.make_evidence_record(
        identity,
        evidence_key=key,
        domain="model_fit",
        question="Did the Python estimator converge?",
        computation_state=evidence.ComputationState.AVAILABLE,
        stability_state=stability,
        required=required,
        summary="The optimizer reported convergence.",
        observed={"converged": True},
        interpretation_boundary="Convergence alone does not establish model fit.",
        recommended_action="Inspect residual and uncertainty evidence next.",
        source_artifacts=("convergence.csv",),
        sensitivity_ids=("sen_link",) if stability is not evidence.StabilityState.NOT_ASSESSED else (),
    )


def _sensitivity(
    baseline: evidence.AnalysisIdentity,
    variant_id: str,
    *,
    changed: bool | None,
    state: evidence.ComputationState = evidence.ComputationState.AVAILABLE,
    reason: evidence.ReasonCode | str | None = None,
    required: bool = True,
) -> evidence.SensitivityRecord:
    variant = _variant(baseline, variant_id, variant=variant_id)
    computed = state in {
        evidence.ComputationState.AVAILABLE,
        evidence.ComputationState.CAUTION,
    }
    return evidence.make_sensitivity_record(
        baseline,
        variant,
        conclusion_id="rater_ordering",
        required=required,
        credible=True,
        computation_state=state,
        comparable=computed,
        summary=f"Sensitivity result for {variant_id}.",
        reason_code=reason,
        metrics=(
            {
                "rank_correlation": 0.90 if changed else 0.99,
                "ordering_changed": bool(changed),
            }
            if computed
            else {}
        ),
        evidence_ids=(_available_evidence(baseline).evidence_id,) if computed else (),
        source_artifacts=("sensitivity.csv",) if computed else (),
    )


def test_public_state_vocabularies_are_exact_and_separate():
    assert [state.value for state in evidence.ComputationState] == [
        "AVAILABLE",
        "CAUTION",
        "HOLD",
        "NOT_ASSESSABLE",
    ]
    assert [state.value for state in evidence.StabilityState] == [
        "STABLE",
        "CONDITIONALLY_STABLE",
        "SENSITIVE",
        "NOT_ASSESSED",
    ]


def test_analysis_identity_is_order_independent_and_sensitive_to_lineage_inputs():
    first = _baseline(alpha=1, nested={"b": 2, "a": 1})
    reordered = evidence.build_analysis_identity(
        input_data_fingerprint="data_0123456789abcdef",
        resolved_settings={
            "nested": {"a": 1, "b": 2},
            "alpha": 1,
            "method": "JMLE",
            "model": "RSM",
        },
        analysis_type="mfrm.fit",
        engine_version="0.3.0-beta",
        seed=20260724,
    )
    changed_setting = _baseline(alpha=2, nested={"b": 2, "a": 1})
    changed_data = evidence.build_analysis_identity(
        input_data_fingerprint="different_data",
        resolved_settings={
            "model": "RSM",
            "method": "JMLE",
            "alpha": 1,
            "nested": {"a": 1, "b": 2},
        },
        analysis_type="mfrm.fit",
        engine_version="0.3.0-beta",
        seed=20260724,
    )
    variant = _variant(first, "leave_one_rater_out.r1", alpha=1)

    assert first.analysis_id == reordered.analysis_id
    assert len(
        {first.analysis_id, changed_setting.analysis_id, changed_data.analysis_id, variant.analysis_id}
    ) == 4
    assert first.baseline_analysis_id == first.analysis_id
    assert variant.baseline_analysis_id == first.analysis_id
    assert evidence.AnalysisIdentity.from_payload(first.to_payload()) == first


def test_variant_identity_requires_an_explicit_baseline():
    with pytest.raises(evidence.ContractValidationError, match="must name its baseline"):
        evidence.build_analysis_identity(
            input_data_fingerprint="data",
            resolved_settings={},
            analysis_type="mfrm.fit",
            engine_version="test",
            variant_id="alternative",
        )


def test_analysis_identity_payload_does_not_coerce_boolean_seed_to_integer():
    payload = _baseline().to_payload()
    payload["seed"] = True
    with pytest.raises(evidence.ContractValidationError, match="seed must be"):
        evidence.AnalysisIdentity.from_payload(payload)


def test_contract_json_rejects_nonstring_keys_and_nonfinite_values():
    assert evidence.canonical_json({"b": 2, "a": 1}) == '{"a":1,"b":2}'
    with pytest.raises(evidence.ContractValidationError, match="NaN or infinity"):
        evidence.payload_fingerprint({"bad": float("nan")})
    with pytest.raises(evidence.ContractValidationError, match="keys must be strings"):
        evidence.canonical_json({1: "integer"})


def test_evidence_round_trip_and_frame_use_fixed_export_contract():
    identity = _baseline()
    record = _available_evidence(identity)
    payload = record.to_payload()
    restored = evidence.EvidenceRecord.from_payload(json.loads(json.dumps(payload)))
    frame = evidence.evidence_records_to_frame([restored], identity=identity)

    assert restored == record
    assert tuple(frame.columns) == evidence.EVIDENCE_ROW_COLUMNS
    assert frame.loc[0, "ComputationState"] == "AVAILABLE"
    assert frame.loc[0, "StabilityState"] == "NOT_ASSESSED"
    assert frame.loc[0, "SourceArtifactsJSON"] == '["convergence.csv"]'


def test_evidence_content_is_deeply_immutable_and_payload_is_plain_json():
    identity = _baseline()
    source = {"nested": {"values": [1, 2]}, "path": Path("results/fit.csv")}
    record = evidence.make_evidence_record(
        identity,
        evidence_key="fit.details",
        domain="model_fit",
        question="What did the fit produce?",
        computation_state=evidence.ComputationState.AVAILABLE,
        summary="Fit details were archived.",
        observed=source,
        interpretation_boundary="The archive does not establish validity.",
        recommended_action="Review the archived diagnostics.",
        source_artifacts=("fit.csv",),
    )
    source["nested"]["values"].append(3)

    assert record.observed["nested"]["values"] == (1, 2)
    assert record.to_payload()["observed"] == {
        "nested": {"values": [1, 2]},
        "path": "results/fit.csv",
    }
    json.dumps(record.to_payload())
    with pytest.raises(TypeError):
        record.observed["new"] = True
    with pytest.raises(TypeError):
        record.observed["nested"]["values"][0] = 9


def test_unavailable_and_caution_evidence_require_reason_codes():
    identity = _baseline()
    common = dict(
        identity=identity,
        evidence_key="dimensionality.residual_pca",
        domain="dimensionality",
        question="Was a residual PCA screen assessable?",
        summary="The prerequisite residual matrix was unavailable.",
        interpretation_boundary="Absence of PCA is not evidence of unidimensionality.",
        recommended_action="Collect enough linked observations and rerun.",
    )
    with pytest.raises(evidence.ContractValidationError, match="requires a stable ReasonCode"):
        evidence.make_evidence_record(
            **common,
            computation_state=evidence.ComputationState.NOT_ASSESSABLE,
        )
    held = evidence.make_evidence_record(
        **common,
        computation_state=evidence.ComputationState.HOLD,
        reason_code=evidence.ReasonCode.PREREQUISITE_MISSING,
    )
    caution = evidence.make_evidence_record(
        **common,
        computation_state=evidence.ComputationState.CAUTION,
        reason_code="pca.too_few_columns",
        source_artifacts=("pca_stability_audit.csv",),
    )

    assert held.reason_code == "evidence.prerequisite_missing"
    assert caution.reason_code == "pca.too_few_columns"
    assert evidence.EvidenceRecord.from_payload(caution.to_payload()).reason_code == caution.reason_code


def test_assessed_evidence_stability_requires_a_sensitivity_link():
    with pytest.raises(evidence.ContractValidationError, match="requires at least one SensitivityID"):
        evidence.make_evidence_record(
            _baseline(),
            evidence_key="conclusion.ordering",
            domain="sensitivity",
            question="Was the ordering stable?",
            computation_state=evidence.ComputationState.AVAILABLE,
            stability_state=evidence.StabilityState.STABLE,
            summary="The ordering was stable.",
            interpretation_boundary="Only the planned variants were assessed.",
            recommended_action="Report the sensitivity plan.",
            source_artifacts=("sensitivity.csv",),
        )


def test_evidence_ledger_rejects_identity_mismatch_and_duplicate_ids():
    first = _baseline(alpha=1)
    second = _baseline(alpha=2)
    first_record = _available_evidence(first)
    second_record = _available_evidence(second)

    with pytest.raises(evidence.ContractValidationError, match="identity mismatch"):
        evidence.validate_evidence_records([first_record, second_record], identity=first)
    with pytest.raises(evidence.ContractValidationError, match="must be unique"):
        evidence.validate_evidence_records([first_record, first_record])


def test_legacy_status_and_disposition_normalization_stay_separate():
    assert evidence.normalize_legacy_status(
        "OK", vocabulary=evidence.LegacyVocabulary.RESOURCE_PREFLIGHT
    ) is evidence.ComputationState.AVAILABLE
    assert evidence.normalize_legacy_status(
        "Missing", vocabulary=evidence.LegacyVocabulary.FINAL_READINESS
    ) is evidence.ComputationState.NOT_ASSESSABLE
    assert evidence.normalize_legacy_disposition(
        "Do not claim",
        vocabulary=evidence.LegacyDispositionVocabulary.CLAIM_GATE,
    ) is evidence.DecisionDisposition.WITHHOLD
    with pytest.raises(evidence.ContractValidationError, match="supply explicit_state"):
        evidence.normalize_legacy_status(
            "Review", vocabulary=evidence.LegacyVocabulary.FINAL_READINESS
        )
    assert evidence.normalize_legacy_status(
        "Review",
        vocabulary=evidence.LegacyVocabulary.FINAL_READINESS,
        explicit_state=evidence.ComputationState.HOLD,
    ) is evidence.ComputationState.HOLD
    with pytest.raises(evidence.ContractValidationError, match="conflicts"):
        evidence.normalize_legacy_status(
            "OK",
            vocabulary=evidence.LegacyVocabulary.RESOURCE_PREFLIGHT,
            explicit_state=evidence.ComputationState.HOLD,
        )


def test_sensitivity_plan_is_order_independent_versioned_and_rule_sensitive():
    baseline = _baseline()
    first = _plan(baseline, "prior_sd.low", "prior_sd.high")
    reordered = _plan(baseline, "prior_sd.high", "prior_sd.low")
    changed_rule = evidence.make_sensitivity_plan(
        baseline,
        conclusion_id="rater_ordering",
        required_variants=(
            _variant(baseline, "prior_sd.low", variant="prior_sd.low"),
            _variant(baseline, "prior_sd.high", variant="prior_sd.high"),
        ),
        decision_rule=evidence.make_sensitivity_rule(
            criteria=(
                {"metric": "rank_correlation", "operator": "lt", "value": 0.95},
            )
        ),
    )

    assert first == reordered
    assert first.plan_id != changed_rule.plan_id
    assert evidence.SensitivityPlan.from_payload(json.loads(json.dumps(first.to_payload()))) == first


def test_sensitivity_plan_fixes_the_exact_variant_analysis_identity():
    baseline = _baseline()
    plan = _plan(baseline, "prior_sd.low")
    off_plan_variant = _variant(
        baseline,
        "prior_sd.low",
        variant="prior_sd.low",
        prior_scale=99,
    )
    off_plan_record = evidence.make_sensitivity_record(
        baseline,
        off_plan_variant,
        conclusion_id="rater_ordering",
        required=True,
        credible=True,
        computation_state=evidence.ComputationState.AVAILABLE,
        comparable=True,
        summary="This reused the label with a different configuration.",
        metrics={"rank_correlation": 0.99, "ordering_changed": False},
        evidence_ids=(_available_evidence(baseline).evidence_id,),
        source_artifacts=("off_plan.csv",),
    )

    with pytest.raises(evidence.ContractValidationError, match="prespecified AnalysisID"):
        evidence.synthesize_sensitivity_decision([off_plan_record], plan=plan)


def test_sensitivity_plan_embeds_and_validates_its_identity_ledger():
    baseline = _baseline()
    plan = _plan(baseline, "prior_sd.low")
    foreign_baseline = _baseline(alpha=999)
    foreign_variant = _variant(
        foreign_baseline,
        "prior_sd.low",
        variant="prior_sd.low",
    )
    payload = plan.to_payload()
    payload["required_variant_analysis_ids"]["prior_sd.low"] = (
        foreign_variant.analysis_id
    )
    payload["required_variant_identity_payloads"]["prior_sd.low"] = (
        foreign_variant.to_payload()
    )
    payload["plan_id"] = "spl_" + evidence.payload_fingerprint(
        {
            key: payload[key]
            for key in (
                "schema_version",
                "baseline_analysis_id",
                "conclusion_id",
                "baseline_identity_payload",
                "required_variant_analysis_ids",
                "required_variant_identity_payloads",
                "decision_rule",
            )
        },
        length=24,
    )

    with pytest.raises(evidence.ContractValidationError, match="plan baseline"):
        evidence.SensitivityPlan.from_payload(payload)


def test_executable_sensitivity_rule_controls_the_conclusion():
    baseline = _baseline()
    variant = _variant(baseline, "prior_sd.low", variant="prior_sd.low")
    record = _sensitivity(baseline, "prior_sd.low", changed=False)
    stable_plan = _plan(baseline, "prior_sd.low")
    sensitive_plan = evidence.make_sensitivity_plan(
        baseline,
        conclusion_id="rater_ordering",
        required_variants=(variant,),
        decision_rule=evidence.make_sensitivity_rule(
            criteria=(
                {"metric": "rank_correlation", "operator": "lt", "value": 1.0},
            )
        ),
    )

    stable = evidence.synthesize_sensitivity_decision([record], plan=stable_plan)
    sensitive = evidence.synthesize_sensitivity_decision([record], plan=sensitive_plan)

    assert stable.stability_state is evidence.StabilityState.STABLE
    assert sensitive.stability_state is evidence.StabilityState.SENSITIVE
    assert sensitive.changed_variant_ids == ("prior_sd.low",)


def test_sensitivity_rule_uses_fail_closed_three_valued_logic():
    baseline = _baseline()
    variant = _variant(baseline, "prior_sd.low", variant="prior_sd.low")
    plan = _plan(baseline, "prior_sd.low")

    def comparison(rank_correlation: float) -> evidence.SensitivityRecord:
        return evidence.make_sensitivity_record(
            baseline,
            variant,
            conclusion_id="rater_ordering",
            required=True,
            credible=True,
            computation_state=evidence.ComputationState.AVAILABLE,
            comparable=True,
            summary="Only one of the two planned rule metrics was archived.",
            metrics={"rank_correlation": rank_correlation},
            evidence_ids=(_available_evidence(baseline).evidence_id,),
            source_artifacts=("partial_metrics.csv",),
        )

    triggered = evidence.synthesize_sensitivity_decision(
        [comparison(0.90)], plan=plan
    )
    indeterminate = evidence.synthesize_sensitivity_decision(
        [comparison(0.99)], plan=plan
    )

    assert triggered.stability_state is evidence.StabilityState.SENSITIVE
    assert indeterminate.stability_state is evidence.StabilityState.NOT_ASSESSED
    assert indeterminate.reason_code == "sensitivity.rule_not_evaluable"


def test_sensitivity_rule_rejects_unknown_shape_and_boolean_thresholds():
    with pytest.raises(evidence.ContractValidationError, match="exactly"):
        evidence.make_sensitivity_rule(
            criteria=(
                {
                    "metric": "rank_correlation",
                    "operator": "lt",
                    "value": 0.98,
                    "unexpected": True,
                },
            )
        )
    with pytest.raises(evidence.ContractValidationError, match="numeric value"):
        evidence.make_sensitivity_rule(
            criteria=(
                {"metric": "rank_correlation", "operator": "lt", "value": True},
            )
        )


@pytest.mark.parametrize("wrong_type", ["true", {"value": True}, [True], 0, None])
def test_boolean_rule_type_mismatch_is_indeterminate(wrong_type):
    baseline = _baseline()
    variant = _variant(baseline, "prior_sd.low", variant="prior_sd.low")
    plan = evidence.make_sensitivity_plan(
        baseline,
        conclusion_id="rater_ordering",
        required_variants=(variant,),
        decision_rule=evidence.make_sensitivity_rule(
            criteria=(
                {"metric": "ordering_changed", "operator": "eq", "value": True},
            )
        ),
    )
    record = evidence.make_sensitivity_record(
        baseline,
        variant,
        conclusion_id="rater_ordering",
        required=True,
        credible=True,
        computation_state=evidence.ComputationState.AVAILABLE,
        comparable=True,
        summary="The rule metric has the wrong JSON type.",
        metrics={"ordering_changed": wrong_type},
        evidence_ids=(_available_evidence(baseline).evidence_id,),
        source_artifacts=("typed_metrics.json",),
    )

    decision = evidence.synthesize_sensitivity_decision([record], plan=plan)

    assert decision.stability_state is evidence.StabilityState.NOT_ASSESSED
    assert decision.reason_code == "sensitivity.rule_not_evaluable"


def test_numeric_rule_evaluator_is_total_for_large_json_integers():
    baseline = _baseline()
    variant = _variant(baseline, "prior_sd.low", variant="prior_sd.low")
    plan = evidence.make_sensitivity_plan(
        baseline,
        conclusion_id="rater_ordering",
        required_variants=(variant,),
        decision_rule=evidence.make_sensitivity_rule(
            criteria=(
                {"metric": "rank_correlation", "operator": "lt", "value": 0.98},
            )
        ),
    )
    record = evidence.make_sensitivity_record(
        baseline,
        variant,
        conclusion_id="rater_ordering",
        required=True,
        credible=True,
        computation_state=evidence.ComputationState.AVAILABLE,
        comparable=True,
        summary="A very large JSON integer remains evaluable without float coercion.",
        metrics={"rank_correlation": 10**1000},
        evidence_ids=(_available_evidence(baseline).evidence_id,),
        source_artifacts=("large_integer.json",),
    )

    decision = evidence.synthesize_sensitivity_decision([record], plan=plan)

    assert decision.stability_state is evidence.StabilityState.STABLE


def test_oversized_integers_fail_as_contract_errors_at_every_entry_point():
    too_large = 10**5000
    with pytest.raises(evidence.ContractValidationError, match="4096 bits"):
        evidence.payload_fingerprint({"value": too_large})
    with pytest.raises(evidence.ContractValidationError, match="4096 bits"):
        evidence.make_sensitivity_rule(
            criteria=(
                {"metric": "rank_correlation", "operator": "lt", "value": too_large},
            )
        )

    baseline = _baseline()
    variant = _variant(baseline, "prior_sd.low", variant="prior_sd.low")
    with pytest.raises(evidence.ContractValidationError, match="4096 bits"):
        evidence.make_sensitivity_record(
            baseline,
            variant,
            conclusion_id="rater_ordering",
            required=True,
            credible=True,
            computation_state=evidence.ComputationState.AVAILABLE,
            comparable=True,
            summary="The metric exceeds the canonical JSON integer boundary.",
            metrics={"rank_correlation": too_large},
            evidence_ids=(_available_evidence(baseline).evidence_id,),
            source_artifacts=("oversized_integer.json",),
        )


def test_sensitivity_record_metrics_are_immutable_and_json_serializable():
    baseline = _baseline()
    record = _sensitivity(baseline, "prior_sd.low", changed=False)
    json.dumps(record.to_payload())
    with pytest.raises(TypeError):
        record.metrics["rank_correlation"] = 0.5


def test_sensitivity_decision_is_stable_only_when_every_planned_variant_completes():
    baseline = _baseline()
    plan = _plan(baseline, "prior_sd.low", "prior_sd.high")
    records = [
        _sensitivity(baseline, "prior_sd.low", changed=False),
        _sensitivity(baseline, "prior_sd.high", changed=False),
    ]
    forward = evidence.synthesize_sensitivity_decision(records, plan=plan)
    reverse = evidence.synthesize_sensitivity_decision(reversed(records), plan=plan)

    assert forward == reverse
    assert forward.stability_state is evidence.StabilityState.STABLE
    assert forward.reason_code == "sensitivity.stable"
    assert evidence.SensitivityDecision.from_payload(
        json.loads(json.dumps(forward.to_payload())),
        plan=plan,
        records=records,
    ) == forward


def test_changed_planned_variant_makes_conclusion_sensitive_even_if_another_fails():
    baseline = _baseline()
    plan = _plan(baseline, "prior_sd.low", "prior_sd.high")
    records = [
        _sensitivity(baseline, "prior_sd.low", changed=True),
        _sensitivity(
            baseline,
            "prior_sd.high",
            changed=None,
            state=evidence.ComputationState.HOLD,
            reason=evidence.ReasonCode.RUN_FAILED,
        ),
    ]
    decision = evidence.synthesize_sensitivity_decision(records, plan=plan)

    assert decision.stability_state is evidence.StabilityState.SENSITIVE
    assert decision.reason_code == "sensitivity.conclusion_changed"
    assert decision.failed_variant_ids == ("prior_sd.high",)


def test_missing_limited_or_noncredible_sensitivity_cannot_be_called_stable():
    baseline = _baseline()
    plan = _plan(baseline, "prior_sd.low", "prior_sd.high")
    completed = _sensitivity(baseline, "prior_sd.low", changed=False)
    incomplete = evidence.synthesize_sensitivity_decision([completed], plan=plan)
    cautious = _sensitivity(
        baseline,
        "prior_sd.high",
        changed=False,
        state=evidence.ComputationState.CAUTION,
        reason=evidence.ReasonCode.SENSITIVITY_THRESHOLD_NEAR,
    )
    limited = evidence.synthesize_sensitivity_decision([completed, cautious], plan=plan)
    noncredible_variant = _variant(baseline, "prior_sd.high", variant="prior_sd.high")
    noncredible = evidence.make_sensitivity_record(
        baseline,
        noncredible_variant,
        conclusion_id="rater_ordering",
        required=True,
        credible=False,
        computation_state=evidence.ComputationState.NOT_ASSESSABLE,
        comparable=False,
        summary="The variant was not a credible comparison.",
        reason_code=evidence.ReasonCode.NONCOMPARABLE,
    )
    credibility_limited = evidence.synthesize_sensitivity_decision(
        [completed, noncredible], plan=plan
    )

    assert incomplete.stability_state is evidence.StabilityState.CONDITIONALLY_STABLE
    assert incomplete.missing_variant_ids == ("prior_sd.high",)
    assert limited.stability_state is evidence.StabilityState.CONDITIONALLY_STABLE
    assert limited.caution_variant_ids == ("prior_sd.high",)
    assert credibility_limited.stability_state is evidence.StabilityState.CONDITIONALLY_STABLE
    assert credibility_limited.noncredible_variant_ids == ("prior_sd.high",)


def test_empty_planned_sensitivity_is_traceable_and_not_assessed():
    baseline = _baseline()
    plan = _plan(baseline, "prior_sd.low")
    decision = evidence.synthesize_sensitivity_decision([], plan=plan)

    assert decision.baseline_analysis_id == baseline.analysis_id
    assert decision.stability_state is evidence.StabilityState.NOT_ASSESSED
    assert decision.reason_code == "sensitivity.missing_required"
    assert decision.missing_variant_ids == ("prior_sd.low",)


def test_sensitivity_decision_payload_rejects_incomplete_variant_partition():
    baseline = _baseline()
    plan = _plan(baseline, "prior_sd.low")
    decision = evidence.synthesize_sensitivity_decision([], plan=plan)
    payload = decision.to_payload()
    payload["missing_variant_ids"] = []

    with pytest.raises(evidence.ContractValidationError, match="partition"):
        evidence.SensitivityDecision.from_payload(payload, plan=plan, records=())


def test_sensitivity_decision_rejects_duplicate_variant_links():
    baseline = _baseline()
    plan = _plan(baseline, "prior_sd.low", "prior_sd.high")
    records = (
        _sensitivity(baseline, "prior_sd.low", changed=False),
        _sensitivity(baseline, "prior_sd.high", changed=False),
    )
    payload = evidence.synthesize_sensitivity_decision(records, plan=plan).to_payload()
    payload["variant_sensitivity_ids"]["prior_sd.high"] = payload[
        "variant_sensitivity_ids"
    ]["prior_sd.low"]

    with pytest.raises(evidence.ContractValidationError, match="distinct SensitivityID"):
        evidence.SensitivityDecision.from_payload(payload, plan=plan, records=records)


def test_sensitivity_decision_rejects_rehashed_but_unreproducible_outcome():
    baseline = _baseline()
    plan = _plan(baseline, "prior_sd.low")
    records = (_sensitivity(baseline, "prior_sd.low", changed=False),)
    payload = evidence.synthesize_sensitivity_decision(records, plan=plan).to_payload()
    payload["stability_state"] = "SENSITIVE"
    payload["reason_code"] = "sensitivity.conclusion_changed"
    payload["changed_variant_ids"] = ["prior_sd.low"]
    fingerprint_fields = (
        "schema_version",
        "reason_registry_version",
        "plan_id",
        "baseline_analysis_id",
        "conclusion_id",
        "stability_state",
        "reason_code",
        "summary",
        "required_variant_analysis_ids",
        "sensitivity_ids",
        "variant_sensitivity_ids",
        "record_fingerprints",
        "changed_variant_ids",
        "missing_variant_ids",
        "failed_variant_ids",
        "noncredible_variant_ids",
        "caution_variant_ids",
    )
    payload["decision_id"] = "sdc_" + evidence.payload_fingerprint(
        {key: payload[key] for key in fingerprint_fields}, length=24
    )

    with pytest.raises(evidence.ContractValidationError, match="does not reproduce"):
        evidence.SensitivityDecision.from_payload(payload, plan=plan, records=records)


def test_sensitivity_decision_restore_requires_bundle_context():
    baseline = _baseline()
    plan = _plan(baseline, "prior_sd.low")
    record = _sensitivity(baseline, "prior_sd.low", changed=False)
    payload = evidence.synthesize_sensitivity_decision([record], plan=plan).to_payload()

    with pytest.raises(evidence.ContractValidationError, match="requires linked plan"):
        evidence.SensitivityDecision.from_payload(payload)


def test_decision_record_links_matching_evidence_plan_and_sensitivity():
    baseline = _baseline()
    fit = _available_evidence(baseline)
    plan = _plan(baseline, "prior_sd.low")
    sensitivity_record = _sensitivity(baseline, "prior_sd.low", changed=False)
    sensitivity = evidence.synthesize_sensitivity_decision(
        [sensitivity_record], plan=plan
    )
    decision = evidence.make_decision_record(
        baseline,
        decision_key="report.rater_ordering",
        conclusion_id="rater_ordering",
        title="Rater ordering",
        evidence=[fit],
        sensitivity_plan=plan,
        sensitivity_records=(sensitivity_record,),
        rationale="The required evidence was available and the tested conclusion was stable.",
        interpretation_boundary="This does not classify individual raters as good or bad.",
        recommended_action="Report the tested sensitivity range.",
    )
    restored = evidence.DecisionRecord.from_payload(
        json.loads(json.dumps(decision.to_payload())),
        identity=baseline,
        evidence=(fit,),
        sensitivity_plan=plan,
        sensitivity_records=(sensitivity_record,),
    )

    assert restored == decision
    assert decision.disposition is evidence.DecisionDisposition.USE
    assert decision.evidence_ids == (fit.evidence_id,)
    assert decision.sensitivity_ids == (sensitivity_record.sensitivity_id,)
    assert decision.sensitivity_decision_id == sensitivity.decision_id
    assert decision.to_payload()["reason_codes"] == ["sensitivity.stable"]


def test_decision_rejects_sensitivity_for_a_different_conclusion():
    baseline = _baseline()
    fit = _available_evidence(baseline)
    plan = _plan(baseline, "prior_sd.low")
    sensitivity = evidence.synthesize_sensitivity_decision(
        [_sensitivity(baseline, "prior_sd.low", changed=False)], plan=plan
    )
    with pytest.raises(evidence.ContractValidationError, match="ConclusionID"):
        evidence.make_decision_record(
            baseline,
            decision_key="report.dimensionality",
            conclusion_id="dimensionality",
            title="Dimensionality",
            evidence=[fit],
            sensitivity_plan=plan,
            sensitivity_records=(_sensitivity(baseline, "prior_sd.low", changed=False),),
            rationale="Wrong sensitivity conclusion.",
            interpretation_boundary="Conclusions must match.",
            recommended_action="Use the correct plan.",
        )


def test_decision_can_archive_a_planned_but_not_assessed_sensitivity_attempt():
    baseline = _baseline()
    fit = _available_evidence(baseline)
    plan = _plan(baseline, "prior_sd.low")
    sensitivity = evidence.synthesize_sensitivity_decision([], plan=plan)
    decision = evidence.make_decision_record(
        baseline,
        decision_key="report.rater_ordering",
        conclusion_id="rater_ordering",
        title="Rater ordering",
        evidence=[fit],
        sensitivity_plan=plan,
        sensitivity_records=(),
        rationale="The sensitivity plan was recorded but its variant did not run.",
        interpretation_boundary="Conclusion stability was not assessed.",
        recommended_action="Run the missing planned variant.",
    )

    assert decision.stability_state is evidence.StabilityState.NOT_ASSESSED
    assert decision.sensitivity_decision_id == sensitivity.decision_id
    assert decision.sensitivity_plan_id == plan.plan_id
    assert decision.disposition is evidence.DecisionDisposition.USE_WITH_CAVEAT


def test_decision_without_sensitivity_records_explicit_not_run_reason():
    baseline = _baseline()
    decision = evidence.make_decision_record(
        baseline,
        decision_key="report.fit",
        conclusion_id="fit_reporting",
        title="Fit reporting",
        evidence=(_available_evidence(baseline),),
        rationale="Fit evidence was computed.",
        interpretation_boundary="No sensitivity conclusion was attempted.",
        recommended_action="Run the prespecified sensitivity plan.",
    )

    assert decision.stability_state is evidence.StabilityState.NOT_ASSESSED
    assert "sensitivity.not_run" in decision.reason_codes


def test_decision_restore_requires_and_rechecks_linked_evidence():
    baseline = _baseline()
    fit = _available_evidence(baseline)
    decision = evidence.make_decision_record(
        baseline,
        decision_key="report.fit",
        conclusion_id="fit_reporting",
        title="Fit reporting",
        evidence=(fit,),
        rationale="Fit evidence was computed.",
        interpretation_boundary="Sensitivity was not assessed.",
        recommended_action="Run the sensitivity plan.",
    )
    with pytest.raises(evidence.ContractValidationError, match="requires identity"):
        evidence.DecisionRecord.from_payload(decision.to_payload())

    altered_payload = fit.to_payload()
    altered_payload["summary"] = "The same logical EvidenceID now carries altered content."
    altered_fit = evidence.EvidenceRecord.from_payload(altered_payload)
    with pytest.raises(evidence.ContractValidationError, match="does not reproduce"):
        evidence.DecisionRecord.from_payload(
            decision.to_payload(),
            identity=baseline,
            evidence=(altered_fit,),
        )


def test_payload_restoration_rejects_unknown_keys_registry_and_collection_coercion():
    baseline = _baseline()
    fit = _available_evidence(baseline)
    decision = evidence.make_decision_record(
        baseline,
        decision_key="report.fit",
        conclusion_id="fit_reporting",
        title="Fit reporting",
        evidence=(fit,),
        rationale="Fit evidence was computed.",
        interpretation_boundary="Sensitivity was not assessed.",
        recommended_action="Run the sensitivity plan.",
    )

    unknown = decision.to_payload()
    unknown["future_override"] = "USE"
    with pytest.raises(evidence.ContractValidationError, match="shape mismatch"):
        evidence.DecisionRecord.from_payload(
            unknown,
            identity=baseline,
            evidence=(fit,),
        )

    wrong_registry = decision.to_payload()
    wrong_registry["reason_registry_version"] = "future_registry"
    with pytest.raises(evidence.ContractVersionError, match="reason registry"):
        evidence.DecisionRecord.from_payload(
            wrong_registry,
            identity=baseline,
            evidence=(fit,),
        )

    mapping_instead_of_array = decision.to_payload()
    mapping_instead_of_array["reason_codes"] = {"sensitivity.not_run": True}
    with pytest.raises(evidence.ContractValidationError, match="JSON array"):
        evidence.DecisionRecord.from_payload(
            mapping_instead_of_array,
            identity=baseline,
            evidence=(fit,),
        )

    string_instead_of_array = fit.to_payload()
    string_instead_of_array["source_artifacts"] = "xy"
    with pytest.raises(evidence.ContractValidationError, match="JSON array"):
        evidence.EvidenceRecord.from_payload(string_instead_of_array)


def test_decision_bundle_rejects_unresolved_sensitivity_evidence_link():
    baseline = _baseline()
    fit = _available_evidence(baseline)
    plan = _plan(baseline, "prior_sd.low")
    variant = _variant(baseline, "prior_sd.low", variant="prior_sd.low")
    orphaned = evidence.make_sensitivity_record(
        baseline,
        variant,
        conclusion_id="rater_ordering",
        required=True,
        credible=True,
        computation_state=evidence.ComputationState.AVAILABLE,
        comparable=True,
        summary="The comparison names evidence outside this decision bundle.",
        metrics={"rank_correlation": 0.99, "ordering_changed": False},
        evidence_ids=("evi_not_in_bundle",),
        source_artifacts=("sensitivity.csv",),
    )

    with pytest.raises(evidence.ContractValidationError, match="outside the decision bundle"):
        evidence.make_decision_record(
            baseline,
            decision_key="report.rater_ordering",
            conclusion_id="rater_ordering",
            title="Rater ordering",
            evidence=(fit,),
            sensitivity_plan=plan,
            sensitivity_records=(orphaned,),
            rationale="The link cannot be resolved.",
            interpretation_boundary="Every sensitivity result needs evidence provenance.",
            recommended_action="Include the source evidence record.",
        )


def test_boundary_only_intent_is_explicit_and_cannot_be_inferred_from_disposition():
    baseline = _baseline()
    decision = evidence.make_decision_record(
        baseline,
        decision_key="report.boundary",
        conclusion_id="report_boundary",
        title="Interpretation boundary",
        evidence=(_available_evidence(baseline),),
        boundary_only=True,
        rationale="Only a boundary statement is supported.",
        interpretation_boundary="No substantive stability claim is made.",
        recommended_action="Collect sensitivity evidence before a stronger claim.",
    )
    payload = decision.to_payload()
    payload["boundary_only"] = False

    with pytest.raises(evidence.ContractValidationError, match="inconsistent"):
        evidence.DecisionRecord.from_payload(
            payload,
            identity=baseline,
            evidence=(_available_evidence(baseline),),
        )


def test_required_evidence_fails_closed_while_optional_unavailable_evidence_does_not():
    baseline = _baseline()
    required_fit = _available_evidence(baseline)
    optional_missing = evidence.make_evidence_record(
        baseline,
        evidence_key="optional.pca",
        domain="dimensionality",
        question="Was optional PCA available?",
        computation_state=evidence.ComputationState.NOT_ASSESSABLE,
        required=False,
        reason_code=evidence.ReasonCode.NOT_APPLICABLE,
        summary="Optional PCA was not assessable.",
        interpretation_boundary="This optional screen does not block fit reporting.",
        recommended_action="Collect more overlap if this screen is needed.",
    )
    optional_decision = evidence.make_decision_record(
        baseline,
        decision_key="report.fit",
        conclusion_id="fit_reporting",
        title="Fit reporting",
        evidence=[required_fit, optional_missing],
        rationale="Required fit evidence was available.",
        interpretation_boundary="Optional PCA remains unavailable.",
        recommended_action="Report fit without a dimensionality claim.",
    )
    required_missing = evidence.make_evidence_record(
        baseline,
        evidence_key="required.pca",
        domain="dimensionality",
        question="Was required PCA available?",
        computation_state=evidence.ComputationState.NOT_ASSESSABLE,
        required=True,
        reason_code=evidence.ReasonCode.PREREQUISITE_MISSING,
        summary="Required PCA was not assessable.",
        interpretation_boundary="The planned dimensionality conclusion is unavailable.",
        recommended_action="Collect more linked observations.",
    )
    required_decision = evidence.make_decision_record(
        baseline,
        decision_key="report.dimensionality",
        conclusion_id="dimensionality",
        title="Dimensionality",
        evidence=[required_fit, required_missing],
        rationale="A required screen was unavailable.",
        interpretation_boundary="No dimensionality conclusion is supported.",
        recommended_action="Resolve the missing prerequisite.",
    )

    assert optional_decision.computation_state is evidence.ComputationState.AVAILABLE
    assert required_decision.computation_state is evidence.ComputationState.NOT_ASSESSABLE
    assert required_decision.disposition is evidence.DecisionDisposition.NOT_EVALUATED


def test_hold_evidence_cannot_be_masked_by_boundary_only_routing():
    baseline = _baseline()
    held = evidence.make_evidence_record(
        baseline,
        evidence_key="fit.convergence",
        domain="model_fit",
        question="Did the estimator converge?",
        computation_state=evidence.ComputationState.HOLD,
        reason_code=evidence.ReasonCode.RUN_FAILED,
        summary="The optimizer did not converge.",
        interpretation_boundary="Unconverged measures are not final evidence.",
        recommended_action="Resolve convergence before interpretation.",
    )
    decision = evidence.make_decision_record(
        baseline,
        decision_key="report.final_interpretation",
        conclusion_id="final_interpretation",
        title="Final interpretation",
        evidence=[held],
        boundary_only=True,
        rationale="A prerequisite failed.",
        interpretation_boundary="No final parameter interpretation is supported.",
        recommended_action="Refit after resolving the numerical failure.",
    )

    assert decision.computation_state is evidence.ComputationState.HOLD
    assert decision.disposition is evidence.DecisionDisposition.WITHHOLD


def test_evidence_stability_does_not_replace_a_validated_sensitivity_decision():
    baseline = _baseline()
    sensitive_evidence = _available_evidence(
        baseline,
        stability=evidence.StabilityState.SENSITIVE,
    )
    decision = evidence.make_decision_record(
        baseline,
        decision_key="report.sensitive_result",
        conclusion_id="sensitive_result",
        title="Sensitive result",
        evidence=[sensitive_evidence],
        rationale="No aggregate sensitivity decision was supplied.",
        interpretation_boundary="The evidence annotation alone is not an aggregate decision.",
        recommended_action="Attach the matching sensitivity plan and decision.",
    )

    assert decision.stability_state is evidence.StabilityState.NOT_ASSESSED
    assert decision.disposition is evidence.DecisionDisposition.USE_WITH_CAVEAT


def test_direct_decision_payload_rejects_inconsistent_disposition():
    baseline = _baseline()
    decision = evidence.make_decision_record(
        baseline,
        decision_key="report.fit",
        conclusion_id="fit_reporting",
        title="Fit reporting",
        evidence=[_available_evidence(baseline)],
        rationale="Fit evidence is available.",
        interpretation_boundary="Sensitivity was not assessed.",
        recommended_action="Add sensitivity evidence before a stronger conclusion.",
    )
    payload = decision.to_payload()
    payload["disposition"] = "USE"

    with pytest.raises(evidence.ContractValidationError, match="inconsistent"):
        evidence.DecisionRecord.from_payload(
            payload,
            identity=baseline,
            evidence=(_available_evidence(baseline),),
        )
