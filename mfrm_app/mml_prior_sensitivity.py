"""Pure contract adapter for fixed-prior-SD MML sensitivity results."""

from __future__ import annotations

import math
from typing import Mapping

from . import evidence


MML_PRIOR_SENSITIVITY_BUNDLE_SCHEMA_VERSION = (
    "mfrm_mml_prior_sd_contract_bundle_v1"
)
MML_PRIOR_SENSITIVITY_CONCLUSION_ID = "mml.fixed_prior_sd.measure_stability"


def validate_contract_bundle(bundle: Mapping[str, object]) -> dict[str, object]:
    """Restore and cross-check every linked object in a saved MML contract."""
    expected_keys = {
        "schema_version",
        "analysis_identities",
        "evidence_records",
        "sensitivity_plan",
        "sensitivity_records",
        "sensitivity_decision",
    }
    if not isinstance(bundle, Mapping) or set(bundle) != expected_keys:
        raise evidence.ContractValidationError(
            "MML prior-SD contract bundle has an invalid payload shape"
        )
    if bundle.get("schema_version") != MML_PRIOR_SENSITIVITY_BUNDLE_SCHEMA_VERSION:
        raise evidence.ContractVersionError(
            f"Unsupported MML prior-SD bundle version: {bundle.get('schema_version')!r}"
        )
    identities = bundle.get("analysis_identities")
    if not isinstance(identities, Mapping) or set(identities) != {"baseline", "variants"}:
        raise evidence.ContractValidationError(
            "MML prior-SD analysis identities must contain baseline and variants"
        )
    baseline_payload = identities.get("baseline")
    variant_payloads = identities.get("variants")
    if not isinstance(baseline_payload, Mapping) or not isinstance(
        variant_payloads, Mapping
    ):
        raise evidence.ContractValidationError(
            "MML prior-SD identity payloads must be mappings"
        )
    baseline = evidence.AnalysisIdentity.from_payload(baseline_payload)
    variants = {
        str(variant_id): evidence.AnalysisIdentity.from_payload(payload)
        for variant_id, payload in variant_payloads.items()
        if isinstance(payload, Mapping)
    }
    if len(variants) != len(variant_payloads):
        raise evidence.ContractValidationError(
            "Every MML prior-SD variant identity must be a mapping"
        )

    raw_evidence_records = bundle.get("evidence_records")
    raw_sensitivity_records = bundle.get("sensitivity_records")
    if not isinstance(raw_evidence_records, list) or not isinstance(
        raw_sensitivity_records, list
    ):
        raise evidence.ContractValidationError(
            "MML prior-SD evidence and sensitivity records must be arrays"
        )
    evidence_records = tuple(
        evidence.EvidenceRecord.from_payload(payload)
        for payload in raw_evidence_records
    )
    evidence.validate_evidence_records(evidence_records, identity=baseline)
    plan_payload = bundle.get("sensitivity_plan")
    if not isinstance(plan_payload, Mapping):
        raise evidence.ContractValidationError(
            "MML prior-SD sensitivity plan must be a mapping"
        )
    plan = evidence.SensitivityPlan.from_payload(plan_payload)
    if plan.baseline_analysis_id != baseline.analysis_id:
        raise evidence.ContractValidationError(
            "MML prior-SD plan baseline does not match the bundled identity"
        )
    if set(variants) != set(plan.required_variant_ids):
        raise evidence.ContractValidationError(
            "MML prior-SD bundled variants do not match the plan"
        )
    for variant_id, variant in variants.items():
        if variant.to_payload() != plan.required_variant_identity_payloads[variant_id]:
            raise evidence.ContractValidationError(
                f"MML prior-SD variant {variant_id!r} does not match the plan"
            )
    sensitivity_records = tuple(
        evidence.SensitivityRecord.from_payload(payload)
        for payload in raw_sensitivity_records
    )
    decision_payload = bundle.get("sensitivity_decision")
    if not isinstance(decision_payload, Mapping):
        raise evidence.ContractValidationError(
            "MML prior-SD sensitivity decision must be a mapping"
        )
    decision = evidence.SensitivityDecision.from_payload(
        decision_payload,
        plan=plan,
        records=sensitivity_records,
    )
    evidence_ids = {record.evidence_id for record in evidence_records}
    unresolved = sorted({
        evidence_id
        for record in sensitivity_records
        for evidence_id in record.evidence_ids
        if evidence_id not in evidence_ids
    })
    if unresolved:
        raise evidence.ContractValidationError(
            f"MML prior-SD sensitivity records reference unknown EvidenceIDs: {unresolved!r}"
        )
    return {
        "baseline": baseline,
        "variants": variants,
        "evidence_records": evidence_records,
        "sensitivity_plan": plan,
        "sensitivity_records": sensitivity_records,
        "sensitivity_decision": decision,
    }


def build_decision_rule(
    thresholds: Mapping[str, float],
    *,
    latent_regression: bool,
) -> Mapping[str, object]:
    """Return the executable form of the existing MML review thresholds."""
    criteria: list[dict[str, object]] = [
        {
            "metric": "max_abs_measure_shift",
            "operator": "gt",
            "value": float(thresholds["max_abs_measure_shift"]),
        },
        {
            "metric": "measure_rmse",
            "operator": "gt",
            "value": float(thresholds["measure_rmse"]),
        },
        {
            "metric": "rank_correlation",
            "operator": "lt",
            "value": float(thresholds["rank_correlation"]),
        },
    ]
    if latent_regression:
        criteria.append(
            {
                "metric": "population_coefficient_shift",
                "operator": "gt",
                "value": float(thresholds["population_coefficient_shift"]),
            }
        )
    return evidence.make_sensitivity_rule(criteria=criteria, combine="any")


def build_sensitivity_plan(
    baseline: evidence.AnalysisIdentity,
    variants: Mapping[str, evidence.AnalysisIdentity],
    *,
    thresholds: Mapping[str, float],
    latent_regression: bool,
) -> evidence.SensitivityPlan:
    """Prespecify the exact non-baseline variants and decision rule."""
    if not variants:
        raise evidence.ContractValidationError(
            "MML prior-SD sensitivity requires at least one non-baseline variant"
        )
    return evidence.make_sensitivity_plan(
        baseline,
        conclusion_id=MML_PRIOR_SENSITIVITY_CONCLUSION_ID,
        required_variants=tuple(variants[key] for key in sorted(variants)),
        decision_rule=build_decision_rule(
            thresholds,
            latent_regression=latent_regression,
        ),
    )


def _required_rule_metrics(plan: evidence.SensitivityPlan) -> tuple[str, ...]:
    return tuple(
        str(criterion["metric"])
        for criterion in plan.decision_rule["criteria"]
    )


def _finite_rule_metrics(raw: object) -> dict[str, int | float | bool | str]:
    if not isinstance(raw, Mapping):
        return {}
    normalized: dict[str, int | float | bool | str] = {}
    for key, value in raw.items():
        if isinstance(value, bool):
            normalized[str(key)] = value
        elif isinstance(value, int):
            normalized[str(key)] = value
        elif isinstance(value, float) and math.isfinite(value):
            normalized[str(key)] = value
        elif isinstance(value, str):
            normalized[str(key)] = value
    return normalized


def _required_outcome_boolean(
    outcome: Mapping[str, object],
    key: str,
    *,
    variant_id: str,
) -> bool:
    """Read one outcome flag without truthiness coercion.

    Contract payloads often pass through JSON/CSV boundaries where the strings
    ``"true"`` and ``"false"`` are both truthy in Python.  Accepting those
    values could therefore turn a failed comparison into credible evidence.
    """
    value = outcome.get(key)
    if not isinstance(value, bool):
        raise evidence.ContractValidationError(
            f"MML sensitivity outcome {variant_id!r} field {key!r} must be a boolean"
        )
    return value


def build_contract_bundle(
    *,
    baseline: evidence.AnalysisIdentity,
    variants: Mapping[str, evidence.AnalysisIdentity],
    plan: evidence.SensitivityPlan,
    outcomes: Mapping[str, Mapping[str, object]],
    baseline_converged: bool,
    latent_regression: bool,
) -> dict[str, object]:
    """Create records and a fail-closed decision from completed refits."""
    if set(variants) != set(plan.required_variant_ids):
        raise evidence.ContractValidationError(
            "MML sensitivity variants do not match the prespecified plan"
        )
    if set(outcomes) != set(variants):
        raise evidence.ContractValidationError(
            "MML sensitivity outcomes do not match the prespecified variants"
        )
    if not isinstance(baseline_converged, bool):
        raise evidence.ContractValidationError("baseline_converged must be a boolean")

    required_metrics = set(_required_rule_metrics(plan))
    artifacts = [
        "mml_prior_sd_sensitivity_summary.csv",
        "mml_prior_sd_sensitivity_measure_deltas.csv",
        "mml_prior_sd_sensitivity_settings.json",
    ]
    if latent_regression:
        artifacts.append("mml_prior_sd_sensitivity_population_deltas.csv")
    artifact_tuple = tuple(artifacts)

    classified: dict[str, dict[str, object]] = {}
    for variant_id in sorted(variants):
        outcome = outcomes[variant_id]
        if not isinstance(outcome, Mapping):
            raise evidence.ContractValidationError(
                f"MML sensitivity outcome {variant_id!r} must be a mapping"
            )
        run_ok = _required_outcome_boolean(
            outcome,
            "run_ok",
            variant_id=variant_id,
        )
        converged = _required_outcome_boolean(
            outcome,
            "converged",
            variant_id=variant_id,
        )
        comparable = _required_outcome_boolean(
            outcome,
            "comparable",
            variant_id=variant_id,
        )
        metrics = _finite_rule_metrics(outcome.get("metrics", {}))
        missing_metrics = sorted(required_metrics.difference(metrics))
        if not run_ok:
            state = evidence.ComputationState.HOLD
            reason = evidence.ReasonCode.RUN_FAILED
            comparable = False
            credible = False
        elif not comparable:
            state = evidence.ComputationState.NOT_ASSESSABLE
            reason = evidence.ReasonCode.NONCOMPARABLE
            credible = False
        elif missing_metrics:
            state = evidence.ComputationState.CAUTION
            reason = evidence.ReasonCode.SENSITIVITY_RULE_NOT_EVALUABLE
            credible = False
        elif not baseline_converged or not converged:
            state = evidence.ComputationState.CAUTION
            reason = evidence.ReasonCode.LIMITED_EVIDENCE
            credible = False
        else:
            state = evidence.ComputationState.AVAILABLE
            reason = None
            credible = True
        classified[variant_id] = {
            "state": state,
            "reason": reason,
            "credible": credible,
            "comparable": comparable,
            "metrics": metrics,
            "missing_metrics": missing_metrics,
            "summary": str(
                outcome.get("summary")
                or f"Fixed-prior-SD comparison for {variant_id}."
            ),
        }

    computed = [
        item
        for item in classified.values()
        if item["state"]
        in {evidence.ComputationState.AVAILABLE, evidence.ComputationState.CAUTION}
    ]
    if not computed:
        aggregate_state = evidence.ComputationState.NOT_ASSESSABLE
        aggregate_reason = evidence.ReasonCode.SENSITIVITY_MISSING_REQUIRED
        aggregate_sources: tuple[str, ...] = ()
    elif any(item["state"] is evidence.ComputationState.CAUTION for item in computed) or len(computed) < len(classified):
        aggregate_state = evidence.ComputationState.CAUTION
        aggregate_reason = evidence.ReasonCode.LIMITED_EVIDENCE
        aggregate_sources = artifact_tuple
    else:
        aggregate_state = evidence.ComputationState.AVAILABLE
        aggregate_reason = None
        aggregate_sources = artifact_tuple

    aggregate_evidence = evidence.make_evidence_record(
        baseline,
        evidence_key="mml.prior_sd.comparison_artifacts",
        domain="sensitivity",
        question=(
            "Are fitted measures stable across the prespecified fixed population "
            "prior SD variants?"
        ),
        computation_state=aggregate_state,
        reason_code=aggregate_reason,
        summary=(
            f"{len(computed)} of {len(classified)} required fixed-prior-SD "
            "comparison(s) produced comparable evidence."
        ),
        observed={
            "completed_comparisons": len(computed),
            "required_comparisons": len(classified),
            "baseline_converged": baseline_converged,
        },
        interpretation_boundary=(
            "This screen varies a fixed quadrature population SD. It does not "
            "estimate latent variance or validate the full MML model."
        ),
        recommended_action=(
            "Inspect the prespecified variant ledger and report any sensitivity, "
            "failed run, or non-convergence."
        ),
        source_artifacts=aggregate_sources,
    )

    records: list[evidence.SensitivityRecord] = []
    for variant_id in sorted(variants):
        item = classified[variant_id]
        state = item["state"]
        is_computed = state in {
            evidence.ComputationState.AVAILABLE,
            evidence.ComputationState.CAUTION,
        }
        records.append(
            evidence.make_sensitivity_record(
                baseline,
                variants[variant_id],
                conclusion_id=MML_PRIOR_SENSITIVITY_CONCLUSION_ID,
                required=True,
                credible=bool(item["credible"]),
                computation_state=state,
                comparable=bool(item["comparable"]) if is_computed else False,
                reason_code=item["reason"],
                summary=str(item["summary"]),
                metrics=item["metrics"],
                evidence_ids=(aggregate_evidence.evidence_id,) if is_computed else (),
                source_artifacts=artifact_tuple if is_computed else (),
            )
        )

    decision = evidence.synthesize_sensitivity_decision(records, plan=plan)
    bundle = {
        "schema_version": MML_PRIOR_SENSITIVITY_BUNDLE_SCHEMA_VERSION,
        "analysis_identities": {
            "baseline": baseline.to_payload(),
            "variants": {
                key: variants[key].to_payload() for key in sorted(variants)
            },
        },
        "evidence_records": [aggregate_evidence.to_payload()],
        "sensitivity_plan": plan.to_payload(),
        "sensitivity_records": [record.to_payload() for record in records],
        "sensitivity_decision": decision.to_payload(),
    }
    validate_contract_bundle(bundle)
    return bundle
