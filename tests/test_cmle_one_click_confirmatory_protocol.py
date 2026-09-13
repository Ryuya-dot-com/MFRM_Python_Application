from __future__ import annotations

import io
import math
import zipfile

import pandas as pd

from mfrm_app.cmle_one_click_confirmatory_gate import (
    BLOCKED_PRIMARY_CASES,
    INSTRUMENT_VERSION,
)
from mfrm_app.cmle_one_click_confirmatory_protocol import (
    LEDGER_COLUMNS,
    REGISTERED_INVALIDITY_CODES,
    build_attempt_denominator_data_dictionary,
    build_attempt_denominator_ledger_template,
    build_ethics_governance_evidence_template,
    build_invalidity_reason_codebook,
    build_partial_identification_sensitivity,
    build_private_protocol_preflight_bundle,
    build_protocol_decision_register,
    compute_protocol_content_identity,
    evaluate_protocol_readiness,
    invalidity_partial_identification,
    protocol_preflight_human_gate_status,
    validate_attempt_denominator_ledger,
    validate_private_protocol_preflight_bundle,
)


def valid_ledger_row(slot: str = "CF-EN-0001") -> dict[str, object]:
    return {
        "StudyRecordId": "REC-0001",
        "ConfirmatorySlotId": slot,
        "InstrumentVersion": INSTRUMENT_VERSION,
        "InstrumentContentSHA256": "a" * 64,
        "Language": "en",
        "SiteCode": "SITE-01",
        "ModeratorCode": "MOD-01",
        "CollectionWave": "WAVE-01",
        "BlockedMechanism": BLOCKED_PRIMARY_CASES[0],
        "ConsentDisposition": "consented",
        "EligibilityDisposition": "eligible",
        "SessionStarted": True,
        "ResponseLockCompleted": True,
        "CompletionDisposition": "complete",
        "InvalidityReasonCode": "none",
        "PrimaryAnalysisDisposition": "included",
        "ExclusionDecisionStage": "not_applicable",
        "AccessibilityStratumCode": "ACCESS-01",
        "DeviceStratumCode": "DEVICE-01",
        "ProtocolDeviationCode": "none",
        "PIIReviewed": True,
        "AuditReviewed": True,
    }


def protocol_components() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, str]:
    register = build_protocol_decision_register()
    codebook = build_invalidity_reason_codebook()
    dictionary = build_attempt_denominator_data_dictionary()
    sensitivity = build_partial_identification_sensitivity()
    identity = compute_protocol_content_identity(
        register, codebook, dictionary, sensitivity
    )
    return register, codebook, dictionary, sensitivity, identity


def test_default_register_is_deliberately_blocked() -> None:
    register = build_protocol_decision_register()
    assert len(register) == 18
    assert register["DecisionId"].is_unique
    fixed = register.loc[~register["Blocking"]]
    blocking = register.loc[register["Blocking"]]
    assert len(fixed) == 5
    assert len(blocking) == 13
    assert fixed["Status"].eq("fixed_by_parent_contract").all()
    assert blocking["Status"].eq("blocked_unresolved").all()
    assert blocking["EvidenceReference"].eq("").all()
    assert not register["RepositoryMaySelfApprove"].any()


def test_repository_template_cannot_self_authorize() -> None:
    register, _, _, _, identity = protocol_components()
    ethics = build_ethics_governance_evidence_template(identity)
    result = evaluate_protocol_readiness(
        register, ethics, expected_protocol_content_sha256=identity
    )
    assert not result["RecruitmentReady"]
    assert result["Status"] == "blocked_pre_recruitment"
    assert "unresolved_blocking_decisions" in result["BlockingCodes"]
    assert "repository_template_is_not_external_approval" in result["BlockingCodes"]
    assert "ethics_approval_unavailable" in result["BlockingCodes"]
    assert "recruitment_not_authorized" in result["BlockingCodes"]


def test_synthetic_external_evidence_path_requires_every_condition() -> None:
    register, _, _, _, identity = protocol_components()
    resolved = register.copy()
    mask = resolved["Blocking"]
    resolved.loc[mask, "Status"] = "resolved_prospectively"
    resolved.loc[mask, "EvidenceReference"] = "SYNTHETIC-EXTERNAL-EVIDENCE"
    ethics = build_ethics_governance_evidence_template(identity)
    ethics.loc[0, [
        "ReviewBodyReference",
        "EthicsApprovalReference",
        "RecruitmentAuthorizationReference",
        "PrivacyReviewReference",
        "DataRetentionPolicyReference",
        "IncidentResponseReference",
    ]] = "SYNTHETIC-REFERENCE"
    ethics.loc[0, "EvidenceSource"] = "synthetic_external_fixture"
    ethics.loc[0, "EthicsApprovalAvailable"] = True
    ethics.loc[0, "RecruitmentAuthorized"] = True
    ethics.loc[0, "ScopeIdentityVerified"] = True
    result = evaluate_protocol_readiness(
        resolved, ethics, expected_protocol_content_sha256=identity
    )
    assert result["RecruitmentReady"]
    ethics.loc[0, "ProtocolContentSHA256"] = "b" * 64
    mismatch = evaluate_protocol_readiness(
        resolved, ethics, expected_protocol_content_sha256=identity
    )
    assert not mismatch["RecruitmentReady"]
    assert "ethics_scope_identity_mismatch" in mismatch["BlockingCodes"]


def test_zero_row_ledger_template_has_exact_schema() -> None:
    template = build_attempt_denominator_ledger_template()
    assert len(template) == 0
    assert tuple(template.columns) == LEDGER_COLUMNS
    assert validate_attempt_denominator_ledger(template)["valid"]
    missing = validate_attempt_denominator_ledger(
        template, expected_slot_ids=["CF-EN-0001"]
    )
    assert not missing["valid"]
    assert "missing_expected_slots" in missing["failure_codes"]


def test_valid_populated_ledger_and_frozen_coverage_pass() -> None:
    ledger = pd.DataFrame([valid_ledger_row()], columns=LEDGER_COLUMNS)
    result = validate_attempt_denominator_ledger(
        ledger,
        expected_slot_ids=["CF-EN-0001"],
        instrument_version=INSTRUMENT_VERSION,
        instrument_content_sha256="a" * 64,
    )
    assert result["valid"]
    assert result["failure_codes"] == ()


def test_direct_identifier_and_free_text_columns_fail() -> None:
    ledger = pd.DataFrame([valid_ledger_row()], columns=LEDGER_COLUMNS)
    ledger["ParticipantEmail"] = "forbidden@example.invalid"
    result = validate_attempt_denominator_ledger(ledger)
    assert not result["valid"]
    assert "banned_identifier_or_free_text_column" in result["failure_codes"]
    assert "invalid_exact_column_schema" in result["failure_codes"]
    assert result["banned_columns"] == ("ParticipantEmail",)


def test_missing_duplicate_and_unexpected_slots_fail() -> None:
    first = valid_ledger_row("CF-EN-0001")
    second = valid_ledger_row("CF-EN-9999")
    second["StudyRecordId"] = "REC-0002"
    ledger = pd.DataFrame([first, second], columns=LEDGER_COLUMNS)
    result = validate_attempt_denominator_ledger(
        ledger, expected_slot_ids=["CF-EN-0001", "CF-EN-0002"]
    )
    assert not result["valid"]
    assert "missing_expected_slots" in result["failure_codes"]
    assert "unexpected_slots" in result["failure_codes"]
    duplicate = ledger.copy()
    duplicate.loc[1, "ConfirmatorySlotId"] = "CF-EN-0001"
    duplicate_result = validate_attempt_denominator_ledger(duplicate)
    assert "duplicate_confirmatory_slot_id" in duplicate_result["failure_codes"]


def test_unregistered_reason_and_post_outcome_exclusion_fail() -> None:
    row = valid_ledger_row()
    row["InvalidityReasonCode"] = "researcher_discretion"
    row["PrimaryAnalysisDisposition"] = "excluded"
    row["ExclusionDecisionStage"] = "post_outcome"
    ledger = pd.DataFrame([row], columns=LEDGER_COLUMNS)
    result = validate_attempt_denominator_ledger(ledger)
    assert not result["valid"]
    assert "unregistered_invalidity_reason" in result["failure_codes"]
    assert "post_outcome_or_unknown_exclusion_stage" in result["failure_codes"]


def test_invalidity_codebook_is_closed_and_retains_attempts() -> None:
    codebook = build_invalidity_reason_codebook()
    assert tuple(codebook["InvalidityReasonCode"]) == REGISTERED_INVALIDITY_CODES
    assert codebook["CountsInAttemptedDenominator"].all()
    assert codebook["CountsAsValidEligible"].sum() == 1
    assert not codebook["FreeTextSubstitutionAllowed"].any()
    amendment = codebook.loc[
        codebook["InvalidityReasonCode"].eq("other_requires_protocol_amendment")
    ].iloc[0]
    assert bool(amendment["RequiresProtocolAmendment"])


def test_partial_identification_reproduces_integer_formulas() -> None:
    result = invalidity_partial_identification(5, 100, 10)
    assert result["ObservedValidDangerousProportionRaw"] == 0.05
    assert result["AllInvalidSafeLowerRiskRaw"] == 5 / 110
    assert result["AllInvalidDangerousUpperRiskRaw"] == 15 / 110
    assert result["AllInvalidDangerousWilsonUpper95Raw"] > result[
        "ObservedValidWilsonUpper95Raw"
    ]
    assert not result["RobustToInvalidOutcomes"]
    assert not result["DecisionUsesDisplayedRounding"]
    assert not result["SampleSizeSelected"]


def test_partial_identification_surface_is_complete_and_monotone() -> None:
    surface = build_partial_identification_sensitivity()
    assert len(surface) == 4 * 3 * 6
    assert not surface["SampleSizeSelected"].any()
    assert not surface["MinimumValidNRegistered"].any()
    assert not surface["PlannedRecruitmentNSelected"].any()
    assert surface["HumanParticipants"].eq(0).all()
    groups = surface.groupby(
        ["ValidEligibleN", "ObservedDangerousErrorRateRegistered"], sort=False
    )
    for _, frame in groups:
        ordered = frame.sort_values("InvalidEligibleN")
        assert ordered["AllInvalidDangerousUpperRiskRaw"].is_monotonic_increasing
        assert ordered["AllInvalidDangerousWilsonUpper95Raw"].is_monotonic_increasing


def test_robustness_is_stricter_than_observed_valid_pass() -> None:
    surface = build_partial_identification_sensitivity()
    assert (
        ~surface["RobustToInvalidOutcomes"] | surface["ObservedValidWilsonPass"]
    ).all()
    assert (
        surface["RobustToInvalidOutcomes"]
        == surface["AllInvalidDangerousWilsonUpper95Raw"].lt(0.10)
    ).all()
    observed_only = surface.loc[
        surface["ObservedValidWilsonPass"] & ~surface["RobustToInvalidOutcomes"]
    ]
    assert len(observed_only) > 0


def test_protocol_content_identity_changes_with_scientific_content() -> None:
    register, codebook, dictionary, sensitivity, identity = protocol_components()
    changed = register.copy()
    changed.loc[changed["DecisionId"].eq("SCI-01"), "Decision"] += " changed"
    changed_identity = compute_protocol_content_identity(
        changed, codebook, dictionary, sensitivity
    )
    assert changed_identity != identity


def test_private_bundle_is_deterministic_and_valid() -> None:
    first = build_private_protocol_preflight_bundle()
    second = build_private_protocol_preflight_bundle()
    assert first["bundle_bytes"] == second["bundle_bytes"]
    validation = validate_private_protocol_preflight_bundle(first["bundle_bytes"])
    assert validation["valid"]
    assert first["readiness"]["Status"] == "blocked_pre_recruitment"
    assert not first["readiness"]["RecruitmentReady"]
    assert len(first["ledger_template"]) == 0
    assert first["bundle_manifest"]["HumanRows"].eq(0).all()


def test_bundle_member_tampering_fails_hash_validation() -> None:
    result = build_private_protocol_preflight_bundle()
    source = zipfile.ZipFile(io.BytesIO(result["bundle_bytes"]), "r")
    payloads = {name: source.read(name) for name in source.namelist()}
    source.close()
    payloads["START_HERE.md"] = payloads["START_HERE.md"].replace(
        b"BLOCKED", b"READY  ", 1
    )
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(payloads):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            info.create_system = 3
            archive.writestr(info, payloads[name])
    validation = validate_private_protocol_preflight_bundle(buffer.getvalue())
    assert not validation["valid"]
    assert "manifest_hash_mismatch" in validation["failure_codes"]
    assert "first_read_not_fail_closed" in validation["failure_codes"]


def test_human_gate_remains_closed_for_default_bundle() -> None:
    bundle = build_private_protocol_preflight_bundle()
    status = protocol_preflight_human_gate_status(bundle["readiness"]).iloc[0]
    assert status["HumanParticipants"] == 0
    assert status["HumanStudyStatus"] == "not_started_no_human_data"
    assert not bool(status["EthicsApprovalAvailable"])
    assert not bool(status["RecruitmentAuthorized"])
    assert not bool(status["RecruitmentReady"])
    assert not bool(status["MinimumValidPerCellRegistered"])
    assert not bool(status["PlannedRecruitmentNSelected"])
    assert not bool(status["ConfirmatoryResultAvailable"])
    assert not bool(status["PublicSurfaceEnabled"])
