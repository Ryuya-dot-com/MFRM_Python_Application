"""Fail-closed private protocol preflight for the CMLE confirmatory study.

The software contract is successful when unresolved scientific and external
governance decisions remain visibly blocked. Nothing in this module grants
ethics approval, recruitment authority, or public-interface eligibility.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
import hashlib
import io
import math
from pathlib import Path
import re
import zipfile

import pandas as pd

from mfrm_app.cmle_one_click_confirmatory_gate import (
    BLOCKED_PRIMARY_CASES,
    INSTRUMENT_VERSION,
    wilson_upper_bound,
)


PROTOCOL_PREFLIGHT_SCHEMA_VERSION = "cmle_confirmatory_protocol_preflight_v1"
DECISION_THRESHOLD = 0.10
PRIMARY_CONFIDENCE = 0.95
LEDGER_COLUMNS = (
    "StudyRecordId",
    "ConfirmatorySlotId",
    "InstrumentVersion",
    "InstrumentContentSHA256",
    "Language",
    "SiteCode",
    "ModeratorCode",
    "CollectionWave",
    "BlockedMechanism",
    "ConsentDisposition",
    "EligibilityDisposition",
    "SessionStarted",
    "ResponseLockCompleted",
    "CompletionDisposition",
    "InvalidityReasonCode",
    "PrimaryAnalysisDisposition",
    "ExclusionDecisionStage",
    "AccessibilityStratumCode",
    "DeviceStratumCode",
    "ProtocolDeviationCode",
    "PIIReviewed",
    "AuditReviewed",
)
REGISTERED_INVALIDITY_CODES = (
    "none",
    "declined_before_start",
    "withdrawn_before_response_lock",
    "withdrawn_after_response_lock",
    "ineligible_preregistered_rule",
    "technical_failure_before_response_lock",
    "technical_failure_after_response_lock",
    "assignment_or_version_mismatch",
    "duplicate_study_record",
    "unknown_response_code",
    "incomplete_primary_items",
    "moderator_interruption",
    "accessibility_barrier",
    "protocol_deviation_blinded_review",
    "other_requires_protocol_amendment",
)
PRIVATE_BUNDLE_MEMBERS = (
    "BUNDLE_MANIFEST.csv",
    "START_HERE.md",
    "attempt_denominator_data_dictionary.csv",
    "attempt_denominator_ledger_template.csv",
    "ethics_governance_evidence_template.csv",
    "invalidity_reason_codebook.csv",
    "partial_identification_sensitivity.csv",
    "pre_recruitment_readiness_summary.csv",
    "protocol_decision_register.csv",
)
_BANNED_COLUMN_PATTERNS = (
    "name",
    "email",
    "phone",
    "address",
    "dateofbirth",
    "ipaddress",
    "contact",
    "freetext",
    "notes",
    "comment",
)


def _canonical_csv_bytes(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(index=False, float_format="%.17g").encode("utf-8")


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _truth(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def build_protocol_decision_register() -> pd.DataFrame:
    """Return fixed and prospectively unresolved decisions in display order."""

    fixed = (
        (
            "FIX-01",
            "instrument",
            "Frozen content identity and version separation",
            "frozen_by_parent_contract",
        ),
        (
            "FIX-02",
            "estimand",
            "Six directional domains remain separate in English and Japanese",
            "frozen_by_parent_contract",
        ),
        (
            "FIX-03",
            "decision",
            "Integer counts and unrounded binary64 Wilson bounds drive decisions",
            "frozen_by_parent_contract",
        ),
        (
            "FIX-04",
            "denominator",
            "Attempted, invalid, valid, and eligible denominators remain separate",
            "frozen_by_parent_contract",
        ),
        (
            "FIX-05",
            "release",
            "Public UI requires a separately governed release decision",
            "frozen_by_parent_contract",
        ),
    )
    unresolved = (
        ("SCI-01", "estimand", "Target-population estimand under invalid or withdrawn outcomes", "scientific_owner"),
        ("SCI-02", "sample_size", "Minimum valid eligible n per cell", "scientific_owner"),
        ("SCI-03", "recruitment", "Planned recruitment and retention allowance", "scientific_owner"),
        ("SCI-04", "multiplicity", "Role of the 12-cell familywise sensitivity", "statistical_owner"),
        ("SCI-05", "missingness", "Invalidity estimand and MNAR sensitivity policy", "statistical_owner"),
        ("SCI-06", "clustering", "Cluster unit, allocation, and cluster-aware analysis", "statistical_owner"),
        ("SCI-07", "mechanism", "Mechanism-specific protection or justified pooling", "scientific_owner"),
        ("SCI-08", "language", "English/Japanese equivalence objective and method", "language_method_owner"),
        ("SCI-09", "accessibility", "Accessibility, device, and proficiency strata and minima", "accessibility_owner"),
        ("OPS-01", "operations", "Stopping, deviation, and blinded exclusion rules", "study_operations_owner"),
        ("ETH-01", "ethics", "Independent ethics review matching frozen protocol identity", "authorized_institutional_role"),
        ("ETH-02", "consent", "Consent and recruitment authorization matching approval scope", "authorized_institutional_role"),
        ("GOV-01", "governance", "Privacy review, data retention, and incident response", "data_governance_owner"),
    )
    rows: list[dict[str, object]] = []
    order = 0
    for decision_id, category, decision, evidence in fixed:
        order += 1
        rows.append(
            {
                "SchemaVersion": PROTOCOL_PREFLIGHT_SCHEMA_VERSION,
                "DisplayOrder": order,
                "DecisionId": decision_id,
                "Category": category,
                "Decision": decision,
                "Status": "fixed_by_parent_contract",
                "Blocking": False,
                "OwnerRole": "repository_contract",
                "EvidenceRequired": evidence,
                "EvidenceReference": evidence,
                "RepositoryMaySelfApprove": False,
            }
        )
    for decision_id, category, decision, owner in unresolved:
        order += 1
        rows.append(
            {
                "SchemaVersion": PROTOCOL_PREFLIGHT_SCHEMA_VERSION,
                "DisplayOrder": order,
                "DecisionId": decision_id,
                "Category": category,
                "Decision": decision,
                "Status": "blocked_unresolved",
                "Blocking": True,
                "OwnerRole": owner,
                "EvidenceRequired": "external_prospective_resolution",
                "EvidenceReference": "",
                "RepositoryMaySelfApprove": False,
            }
        )
    return pd.DataFrame(rows)


def build_invalidity_reason_codebook() -> pd.DataFrame:
    """Return the closed categorical invalidity/deviation vocabulary."""

    descriptions = {
        "none": ("no_registered_invalidity", "none", False, False),
        "declined_before_start": ("consent declined before session start", "before_lock", False, False),
        "withdrawn_before_response_lock": ("withdrawal before response lock", "before_lock", False, True),
        "withdrawn_after_response_lock": ("withdrawal after response lock", "after_lock", False, True),
        "ineligible_preregistered_rule": ("failed a prospectively registered eligibility rule", "before_lock", False, False),
        "technical_failure_before_response_lock": ("technical failure before response lock", "before_lock", False, False),
        "technical_failure_after_response_lock": ("technical failure after response lock", "after_lock", False, True),
        "assignment_or_version_mismatch": ("assignment or frozen version mismatch", "either", False, False),
        "duplicate_study_record": ("duplicate study-scoped record", "either", False, False),
        "unknown_response_code": ("response code outside frozen vocabulary", "after_lock", False, True),
        "incomplete_primary_items": ("one or more primary items incomplete", "after_lock", False, True),
        "moderator_interruption": ("moderator intervention compromised response", "either", False, True),
        "accessibility_barrier": ("registered task could not be completed accessibly", "either", False, True),
        "protocol_deviation_blinded_review": ("deviation excluded by documented blinded review", "either", False, False),
        "other_requires_protocol_amendment": ("unregistered reason requiring dated amendment", "either", True, True),
    }
    rows = []
    for order, code in enumerate(REGISTERED_INVALIDITY_CODES, start=1):
        description, timing, amendment, outcome_related = descriptions[code]
        rows.append(
            {
                "SchemaVersion": PROTOCOL_PREFLIGHT_SCHEMA_VERSION,
                "DisplayOrder": order,
                "InvalidityReasonCode": code,
                "Description": description,
                "OutcomeTiming": timing,
                "CountsInAttemptedDenominator": True,
                "CountsAsValidEligible": code == "none",
                "RequiresProtocolAmendment": amendment,
                "MayBeOutcomeRelated": outcome_related,
                "FreeTextSubstitutionAllowed": False,
            }
        )
    return pd.DataFrame(rows)


def build_attempt_denominator_ledger_template() -> pd.DataFrame:
    """Return the exact zero-row denominator-ledger schema."""

    return pd.DataFrame(columns=list(LEDGER_COLUMNS))


def build_attempt_denominator_data_dictionary() -> pd.DataFrame:
    """Describe every ledger field without permitting unrestricted text."""

    allowed: dict[str, str] = {
        "StudyRecordId": "non-identifying study-scoped code",
        "ConfirmatorySlotId": "frozen assignment-manifest slot code",
        "InstrumentVersion": INSTRUMENT_VERSION,
        "InstrumentContentSHA256": "64 lowercase hexadecimal characters",
        "Language": "en|ja",
        "SiteCode": "prospectively registered categorical code",
        "ModeratorCode": "prospectively registered categorical code",
        "CollectionWave": "prospectively registered categorical code",
        "BlockedMechanism": "|".join(BLOCKED_PRIMARY_CASES),
        "ConsentDisposition": "consented|declined|withdrawn",
        "EligibilityDisposition": "eligible|ineligible|pending",
        "SessionStarted": "true|false",
        "ResponseLockCompleted": "true|false",
        "CompletionDisposition": "complete|partial|not_started",
        "InvalidityReasonCode": "closed invalidity_reason_codebook.csv vocabulary",
        "PrimaryAnalysisDisposition": "included|excluded|pending",
        "ExclusionDecisionStage": "not_applicable|preregistered_automatic|blinded_review|pending",
        "AccessibilityStratumCode": "prospectively registered categorical code",
        "DeviceStratumCode": "prospectively registered categorical code",
        "ProtocolDeviationCode": "none|prospectively registered categorical code",
        "PIIReviewed": "true|false attestation",
        "AuditReviewed": "true|false attestation",
    }
    rows = []
    for order, column in enumerate(LEDGER_COLUMNS, start=1):
        rows.append(
            {
                "SchemaVersion": PROTOCOL_PREFLIGHT_SCHEMA_VERSION,
                "ColumnOrder": order,
                "Column": column,
                "RequiredForPopulatedRow": True,
                "DataType": (
                    "boolean"
                    if column
                    in {
                        "SessionStarted",
                        "ResponseLockCompleted",
                        "PIIReviewed",
                        "AuditReviewed",
                    }
                    else "categorical_code"
                ),
                "AllowedOrDefinition": allowed[column],
                "DirectIdentifierAllowed": False,
                "UnrestrictedFreeTextAllowed": False,
            }
        )
    return pd.DataFrame(rows)


def _banned_columns(columns: Iterable[object]) -> tuple[str, ...]:
    banned = []
    allowed_name = {re.sub(r"[^a-z0-9]", "", value.lower()) for value in LEDGER_COLUMNS}
    for column in columns:
        normalized = re.sub(r"[^a-z0-9]", "", str(column).lower())
        if normalized in allowed_name:
            continue
        if any(pattern in normalized for pattern in _BANNED_COLUMN_PATTERNS):
            banned.append(str(column))
    return tuple(banned)


def validate_attempt_denominator_ledger(
    ledger: pd.DataFrame,
    *,
    expected_slot_ids: Sequence[str] | None = None,
    instrument_version: str | None = None,
    instrument_content_sha256: str | None = None,
) -> dict[str, object]:
    """Fail closed on schema, PII-like columns, coverage, and outcome exclusion."""

    if not isinstance(ledger, pd.DataFrame):
        return {
            "valid": False,
            "failure_codes": ("ledger_not_dataframe",),
            "row_count": 0,
        }
    failures: list[str] = []
    banned = _banned_columns(ledger.columns)
    if banned:
        failures.append("banned_identifier_or_free_text_column")
    if tuple(ledger.columns) != LEDGER_COLUMNS:
        failures.append("invalid_exact_column_schema")
        return {
            "valid": False,
            "failure_codes": tuple(dict.fromkeys(failures)),
            "banned_columns": banned,
            "row_count": len(ledger),
        }
    expected = tuple(str(value) for value in (expected_slot_ids or ()))
    if ledger.empty:
        if expected:
            failures.append("missing_expected_slots")
        return {
            "valid": not failures,
            "failure_codes": tuple(dict.fromkeys(failures)),
            "banned_columns": banned,
            "row_count": 0,
        }
    text = ledger.astype(str)
    if text["StudyRecordId"].str.strip().eq("").any():
        failures.append("blank_study_record_id")
    if text["StudyRecordId"].duplicated().any():
        failures.append("duplicate_study_record_id")
    if text["ConfirmatorySlotId"].duplicated().any():
        failures.append("duplicate_confirmatory_slot_id")
    if expected:
        actual = set(text["ConfirmatorySlotId"])
        expected_set = set(expected)
        if expected_set.difference(actual):
            failures.append("missing_expected_slots")
        if actual.difference(expected_set):
            failures.append("unexpected_slots")
        if len(expected) != len(expected_set):
            failures.append("duplicate_expected_slots")
    if instrument_version is not None and not text["InstrumentVersion"].eq(
        str(instrument_version)
    ).all():
        failures.append("instrument_version_mismatch")
    if instrument_content_sha256 is not None and not text[
        "InstrumentContentSHA256"
    ].eq(str(instrument_content_sha256)).all():
        failures.append("instrument_content_identity_mismatch")
    if not text["Language"].isin({"en", "ja"}).all():
        failures.append("unknown_language")
    if not text["BlockedMechanism"].isin(set(BLOCKED_PRIMARY_CASES)).all():
        failures.append("unknown_blocked_mechanism")
    if not text["ConsentDisposition"].isin(
        {"consented", "declined", "withdrawn"}
    ).all():
        failures.append("unknown_consent_disposition")
    if not text["EligibilityDisposition"].isin(
        {"eligible", "ineligible", "pending"}
    ).all():
        failures.append("unknown_eligibility_disposition")
    if not text["CompletionDisposition"].isin(
        {"complete", "partial", "not_started"}
    ).all():
        failures.append("unknown_completion_disposition")
    if not text["PrimaryAnalysisDisposition"].isin(
        {"included", "excluded", "pending"}
    ).all():
        failures.append("unknown_analysis_disposition")
    if not text["ExclusionDecisionStage"].isin(
        {"not_applicable", "preregistered_automatic", "blinded_review", "pending"}
    ).all():
        failures.append("post_outcome_or_unknown_exclusion_stage")
    if not text["InvalidityReasonCode"].isin(set(REGISTERED_INVALIDITY_CODES)).all():
        failures.append("unregistered_invalidity_reason")
    invalid_without_reason = (
        text["PrimaryAnalysisDisposition"].eq("excluded")
        | ~text["CompletionDisposition"].eq("complete")
        | text["ConsentDisposition"].isin({"declined", "withdrawn"})
    ) & text["InvalidityReasonCode"].eq("none")
    if invalid_without_reason.any():
        failures.append("invalid_or_withdrawn_without_reason")
    included_with_reason = text["PrimaryAnalysisDisposition"].eq(
        "included"
    ) & ~text["InvalidityReasonCode"].eq("none")
    if included_with_reason.any():
        failures.append("included_record_has_invalidity_reason")
    for column in (
        "SessionStarted",
        "ResponseLockCompleted",
        "PIIReviewed",
        "AuditReviewed",
    ):
        if not text[column].str.strip().str.lower().isin(
            {"true", "false", "1", "0", "yes", "no"}
        ).all():
            failures.append(f"invalid_boolean_{column}")
    if not ledger["PIIReviewed"].map(_truth).all():
        failures.append("pii_review_incomplete")
    if not ledger["AuditReviewed"].map(_truth).all():
        failures.append("audit_review_incomplete")
    return {
        "valid": not failures,
        "failure_codes": tuple(dict.fromkeys(failures)),
        "banned_columns": banned,
        "row_count": len(ledger),
    }


def invalidity_partial_identification(
    dangerous_errors: int,
    valid_eligible_n: int,
    invalid_eligible_n: int,
) -> dict[str, object]:
    """Return observed and worst-case invalid-outcome bounds and Wilson decisions."""

    values = (dangerous_errors, valid_eligible_n, invalid_eligible_n)
    if any(not isinstance(value, int) for value in values):
        raise ValueError("Counts must be integers.")
    if valid_eligible_n <= 0 or invalid_eligible_n < 0:
        raise ValueError("Require positive valid n and nonnegative invalid n.")
    if dangerous_errors < 0 or dangerous_errors > valid_eligible_n:
        raise ValueError("Dangerous errors must lie within valid eligible n.")
    total = valid_eligible_n + invalid_eligible_n
    observed = dangerous_errors / valid_eligible_n
    lower = dangerous_errors / total
    upper = (dangerous_errors + invalid_eligible_n) / total
    observed_wilson = wilson_upper_bound(
        dangerous_errors, valid_eligible_n, confidence=PRIMARY_CONFIDENCE
    )
    all_safe_wilson = wilson_upper_bound(
        dangerous_errors, total, confidence=PRIMARY_CONFIDENCE
    )
    worst_wilson = wilson_upper_bound(
        dangerous_errors + invalid_eligible_n,
        total,
        confidence=PRIMARY_CONFIDENCE,
    )
    return {
        "DangerousErrors": dangerous_errors,
        "ValidEligibleN": valid_eligible_n,
        "InvalidEligibleN": invalid_eligible_n,
        "TotalEligibleIfObserved": total,
        "ObservedValidDangerousProportionRaw": observed,
        "AllInvalidSafeLowerRiskRaw": lower,
        "AllInvalidDangerousUpperRiskRaw": upper,
        "ObservedValidWilsonUpper95Raw": observed_wilson,
        "AllInvalidSafeWilsonUpper95Raw": all_safe_wilson,
        "AllInvalidDangerousWilsonUpper95Raw": worst_wilson,
        "ObservedValidWilsonPass": observed_wilson < DECISION_THRESHOLD,
        "RobustToInvalidOutcomes": worst_wilson < DECISION_THRESHOLD,
        "DecisionThresholdRaw": DECISION_THRESHOLD,
        "DecisionUsesDisplayedRounding": False,
        "SampleSizeSelected": False,
    }


def build_partial_identification_sensitivity(
    valid_n_values: Sequence[int] = (100, 200, 300, 500),
    observed_error_rates: Sequence[float] = (0.0, 0.025, 0.05),
    invalid_to_valid_ratios: Sequence[float] = (0.0, 0.01, 0.025, 0.05, 0.10, 0.20),
) -> pd.DataFrame:
    """Build the registered no-selection invalid-outcome sensitivity surface."""

    rows = []
    for valid_n in valid_n_values:
        if not isinstance(valid_n, int) or valid_n <= 0:
            raise ValueError("valid_n_values must contain positive integers.")
        for rate in observed_error_rates:
            if not 0.0 <= float(rate) <= 1.0:
                raise ValueError("observed_error_rates must lie in [0, 1].")
            errors = int(math.floor(valid_n * float(rate) + 0.5))
            for ratio in invalid_to_valid_ratios:
                if float(ratio) < 0.0:
                    raise ValueError("invalid_to_valid_ratios cannot be negative.")
                invalid_n = int(math.ceil(valid_n * float(ratio)))
                result = invalidity_partial_identification(
                    errors, valid_n, invalid_n
                )
                rows.append(
                    {
                        "SchemaVersion": PROTOCOL_PREFLIGHT_SCHEMA_VERSION,
                        "ObservedDangerousErrorRateRegistered": float(rate),
                        "InvalidToValidRatioRegistered": float(ratio),
                        "ActualInvalidToValidRatioRaw": invalid_n / valid_n,
                        **result,
                        "MinimumValidNRegistered": False,
                        "PlannedRecruitmentNSelected": False,
                        "HumanParticipants": 0,
                    }
                )
    return pd.DataFrame(rows)


def compute_protocol_content_identity(
    decision_register: pd.DataFrame,
    codebook: pd.DataFrame,
    data_dictionary: pd.DataFrame,
    sensitivity: pd.DataFrame,
) -> str:
    """Bind the scientific decision/template content without claiming approval."""

    digest = hashlib.sha256()
    for label, frame in (
        ("decision_register", decision_register),
        ("invalidity_codebook", codebook),
        ("ledger_dictionary", data_dictionary),
        ("partial_identification", sensitivity),
    ):
        label_bytes = label.encode("utf-8")
        payload = _canonical_csv_bytes(frame)
        digest.update(len(label_bytes).to_bytes(8, "big"))
        digest.update(label_bytes)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return digest.hexdigest()


def build_ethics_governance_evidence_template(
    protocol_content_sha256: str,
) -> pd.DataFrame:
    """Return an explicitly unapproved external-evidence template."""

    return pd.DataFrame(
        [
            {
                "SchemaVersion": PROTOCOL_PREFLIGHT_SCHEMA_VERSION,
                "ProtocolContentSHA256": protocol_content_sha256,
                "ReviewBodyReference": "",
                "EthicsApprovalReference": "",
                "ApprovalScopeReference": "",
                "ApprovalDate": "",
                "ExpiryDate": "",
                "ConsentVersion": "",
                "RecruitmentAuthorizationReference": "",
                "PrivacyReviewReference": "",
                "DataRetentionPolicyReference": "",
                "IncidentResponseReference": "",
                "EvidenceSource": "repository_generated_blank_template",
                "EthicsApprovalAvailable": False,
                "RecruitmentAuthorized": False,
                "ScopeIdentityVerified": False,
                "RepositoryGeneratedEvidenceMayAuthorize": False,
            }
        ]
    )


def evaluate_protocol_readiness(
    decision_register: pd.DataFrame,
    ethics_evidence: pd.DataFrame,
    *,
    expected_protocol_content_sha256: str,
) -> dict[str, object]:
    """Require external evidence and complete prospective decisions."""

    blockers: list[str] = []
    required_register = {
        "DecisionId",
        "Status",
        "Blocking",
        "EvidenceReference",
        "RepositoryMaySelfApprove",
    }
    if not isinstance(decision_register, pd.DataFrame) or not required_register.issubset(
        decision_register.columns
    ):
        return {
            "RecruitmentReady": False,
            "Status": "blocked_invalid_decision_register",
            "BlockingCodes": ("invalid_decision_register",),
        }
    if decision_register["DecisionId"].astype(str).duplicated().any():
        blockers.append("duplicate_decision_id")
    blocking = decision_register.loc[decision_register["Blocking"].map(_truth)]
    unresolved = blocking.loc[~blocking["Status"].astype(str).eq("resolved_prospectively")]
    if len(unresolved):
        blockers.append("unresolved_blocking_decisions")
    blank_evidence = blocking["EvidenceReference"].astype(str).str.strip().eq("")
    if blank_evidence.any():
        blockers.append("missing_blocking_decision_evidence")
    if decision_register["RepositoryMaySelfApprove"].map(_truth).any():
        blockers.append("repository_self_approval_forbidden")
    required_ethics = {
        "ProtocolContentSHA256",
        "ReviewBodyReference",
        "EthicsApprovalReference",
        "RecruitmentAuthorizationReference",
        "PrivacyReviewReference",
        "DataRetentionPolicyReference",
        "IncidentResponseReference",
        "EvidenceSource",
        "EthicsApprovalAvailable",
        "RecruitmentAuthorized",
        "ScopeIdentityVerified",
        "RepositoryGeneratedEvidenceMayAuthorize",
    }
    if not isinstance(ethics_evidence, pd.DataFrame) or len(ethics_evidence) != 1 or not required_ethics.issubset(
        ethics_evidence.columns
    ):
        blockers.append("invalid_ethics_evidence_schema")
    else:
        row = ethics_evidence.iloc[0]
        if str(row["ProtocolContentSHA256"]) != expected_protocol_content_sha256:
            blockers.append("ethics_scope_identity_mismatch")
        references = (
            "ReviewBodyReference",
            "EthicsApprovalReference",
            "RecruitmentAuthorizationReference",
            "PrivacyReviewReference",
            "DataRetentionPolicyReference",
            "IncidentResponseReference",
        )
        if any(not str(row[column]).strip() for column in references):
            blockers.append("missing_external_governance_references")
        if str(row["EvidenceSource"]).startswith("repository_generated"):
            blockers.append("repository_template_is_not_external_approval")
        if not _truth(row["EthicsApprovalAvailable"]):
            blockers.append("ethics_approval_unavailable")
        if not _truth(row["RecruitmentAuthorized"]):
            blockers.append("recruitment_not_authorized")
        if not _truth(row["ScopeIdentityVerified"]):
            blockers.append("approval_scope_not_verified")
        if _truth(row["RepositoryGeneratedEvidenceMayAuthorize"]):
            blockers.append("repository_generated_authorization_forbidden")
    unique = tuple(dict.fromkeys(blockers))
    ready = not unique
    return {
        "RecruitmentReady": ready,
        "Status": "ready_external_authority_verified" if ready else "blocked_pre_recruitment",
        "BlockingCodes": unique,
        "BlockingDecisionCount": int(len(unresolved)),
        "MissingBlockingEvidenceCount": int(blank_evidence.sum()),
        "ProtocolContentSHA256": expected_protocol_content_sha256,
    }


def build_pre_recruitment_readiness_summary(
    readiness: Mapping[str, object],
) -> pd.DataFrame:
    """Return the ordered private first-read cards."""

    cards = (
        (1, "overall", "blocked", "Recruitment remains blocked; external decisions and approvals are absent.", True),
        (2, "instrument", "fixed", "Instrument identity and directional gate are frozen by parent evidence.", False),
        (3, "scientific_design", "blocked", "Estimand, n, multiplicity, clustering, mechanisms, and strata remain unresolved.", True),
        (4, "denominator", "template_ready", "Every scheduled slot has a no-silent-drop categorical ledger schema.", False),
        (5, "invalidity", "sensitivity_ready", "Worst-case invalid-outcome bounds are available but do not identify MNAR.", False),
        (6, "ethics_governance", "blocked", "Repository templates are not ethics, privacy, consent, or recruitment approval.", True),
        (7, "public_release", "withheld", "No public Streamlit route or evidence button is enabled.", True),
    )
    return pd.DataFrame(
        [
            {
                "SchemaVersion": PROTOCOL_PREFLIGHT_SCHEMA_VERSION,
                "CardOrder": order,
                "CardId": card_id,
                "Status": status,
                "FirstRead": first_read,
                "Blocking": blocking,
                "RecruitmentReady": bool(readiness["RecruitmentReady"]),
                "PublicSurfaceEnabled": False,
            }
            for order, card_id, status, first_read, blocking in cards
        ]
    )


def protocol_preflight_human_gate_status(
    readiness: Mapping[str, object],
) -> pd.DataFrame:
    """Return the fail-closed human, ethics, selection, and public status."""

    return pd.DataFrame(
        [
            {
                "SchemaVersion": PROTOCOL_PREFLIGHT_SCHEMA_VERSION,
                "HumanParticipants": 0,
                "HumanStudyStatus": "not_started_no_human_data",
                "EthicsApprovalAvailable": False,
                "RecruitmentAuthorized": False,
                "RecruitmentReady": bool(readiness["RecruitmentReady"]),
                "MinimumValidPerCellRegistered": False,
                "PlannedRecruitmentNSelected": False,
                "ConfirmatoryResultAvailable": False,
                "LanguageEquivalenceEstablished": False,
                "PublicSurfaceEnabled": False,
            }
        ]
    )


def _start_here_text(
    readiness: Mapping[str, object], protocol_content_sha256: str
) -> str:
    blockers = ", ".join(str(value) for value in readiness["BlockingCodes"])
    return f"""# BLOCKED — private CMLE confirmatory protocol preflight

RecruitmentReady=false. PublicSurfaceEnabled=false. HumanParticipants=0.

This deterministic private bundle is a pre-recruitment audit aid. It is not ethics approval, consent authorization, privacy certification, recruitment authority, a selected sample size, or a confirmatory result.

ProtocolContentSHA256: `{protocol_content_sha256}`

Blocking codes: {blockers}

## 日本語

募集準備は未完了です。この非公開パッケージは事前点検用であり、倫理承認、
同意・募集の許可、個人情報保護の認証、標本数の決定、確認的研究結果を
意味しません。未決事項と外部証跡を、結果を見る前に解決してください。
"""


def _deterministic_zip(payloads: Mapping[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(
        buffer, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as archive:
        for name in sorted(payloads):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            info.create_system = 3
            archive.writestr(info, payloads[name], compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    return buffer.getvalue()


def build_private_protocol_preflight_bundle() -> dict[str, object]:
    """Return all no-human-data artifacts and deterministic private ZIP bytes."""

    register = build_protocol_decision_register()
    codebook = build_invalidity_reason_codebook()
    ledger = build_attempt_denominator_ledger_template()
    dictionary = build_attempt_denominator_data_dictionary()
    sensitivity = build_partial_identification_sensitivity()
    protocol_sha = compute_protocol_content_identity(
        register, codebook, dictionary, sensitivity
    )
    ethics = build_ethics_governance_evidence_template(protocol_sha)
    readiness = evaluate_protocol_readiness(
        register, ethics, expected_protocol_content_sha256=protocol_sha
    )
    preflight = build_pre_recruitment_readiness_summary(readiness)
    payloads: dict[str, bytes] = {
        "START_HERE.md": _start_here_text(readiness, protocol_sha).encode("utf-8"),
        "protocol_decision_register.csv": _canonical_csv_bytes(register),
        "pre_recruitment_readiness_summary.csv": _canonical_csv_bytes(preflight),
        "attempt_denominator_ledger_template.csv": _canonical_csv_bytes(ledger),
        "attempt_denominator_data_dictionary.csv": _canonical_csv_bytes(dictionary),
        "invalidity_reason_codebook.csv": _canonical_csv_bytes(codebook),
        "ethics_governance_evidence_template.csv": _canonical_csv_bytes(ethics),
        "partial_identification_sensitivity.csv": _canonical_csv_bytes(sensitivity),
    }
    manifest = pd.DataFrame(
        [
            {
                "Member": name,
                "ByteCount": len(payload),
                "SHA256": _sha256_bytes(payload),
                "HumanRows": 0,
                "ContainsDirectIdentifiers": False,
                "PublicArtifact": False,
            }
            for name, payload in sorted(payloads.items())
        ]
    )
    payloads["BUNDLE_MANIFEST.csv"] = _canonical_csv_bytes(manifest)
    bundle_bytes = _deterministic_zip(payloads)
    return {
        "decision_register": register,
        "invalidity_reason_codebook": codebook,
        "ledger_template": ledger,
        "ledger_data_dictionary": dictionary,
        "partial_identification_sensitivity": sensitivity,
        "ethics_template": ethics,
        "readiness": readiness,
        "preflight_summary": preflight,
        "protocol_content_sha256": protocol_sha,
        "bundle_manifest": manifest,
        "bundle_bytes": bundle_bytes,
    }


def validate_private_protocol_preflight_bundle(bundle_bytes: bytes) -> dict[str, object]:
    """Validate exact members, deterministic metadata, hashes, and blocked first read."""

    failures: list[str] = []
    try:
        with zipfile.ZipFile(io.BytesIO(bundle_bytes), "r") as archive:
            infos = archive.infolist()
            names = tuple(info.filename for info in infos)
            if names != tuple(sorted(PRIVATE_BUNDLE_MEMBERS)):
                failures.append("member_order_or_membership_mismatch")
            if len(names) != len(set(names)):
                failures.append("duplicate_members")
            for info in infos:
                if info.date_time != (1980, 1, 1, 0, 0, 0):
                    failures.append("non_deterministic_timestamp")
                if (info.external_attr >> 16) & 0o777 != 0o644:
                    failures.append("non_deterministic_permissions")
            if "BUNDLE_MANIFEST.csv" not in names:
                failures.append("missing_manifest")
                manifest = pd.DataFrame()
            else:
                manifest = pd.read_csv(
                    io.BytesIO(archive.read("BUNDLE_MANIFEST.csv")), dtype=str
                )
            required = {
                "Member",
                "ByteCount",
                "SHA256",
                "HumanRows",
                "ContainsDirectIdentifiers",
                "PublicArtifact",
            }
            if not required.issubset(manifest.columns):
                failures.append("invalid_manifest_schema")
            else:
                expected_manifest_members = set(names).difference(
                    {"BUNDLE_MANIFEST.csv"}
                )
                if set(manifest["Member"]) != expected_manifest_members:
                    failures.append("manifest_member_mismatch")
                for _, row in manifest.iterrows():
                    name = str(row["Member"])
                    if name not in names:
                        continue
                    payload = archive.read(name)
                    if int(row["ByteCount"]) != len(payload):
                        failures.append("manifest_byte_count_mismatch")
                    if str(row["SHA256"]) != _sha256_bytes(payload):
                        failures.append("manifest_hash_mismatch")
                    if int(row["HumanRows"]) != 0:
                        failures.append("human_rows_present")
                    if _truth(row["ContainsDirectIdentifiers"]):
                        failures.append("direct_identifiers_present")
                    if _truth(row["PublicArtifact"]):
                        failures.append("public_artifact_present")
            if "START_HERE.md" in names:
                start = archive.read("START_HERE.md").decode("utf-8")
                required_text = (
                    "BLOCKED",
                    "RecruitmentReady=false",
                    "PublicSurfaceEnabled=false",
                    "not ethics approval",
                )
                if any(value not in start for value in required_text):
                    failures.append("first_read_not_fail_closed")
            if "attempt_denominator_ledger_template.csv" in names:
                ledger = pd.read_csv(
                    io.BytesIO(archive.read("attempt_denominator_ledger_template.csv"))
                )
                if len(ledger) != 0 or tuple(ledger.columns) != LEDGER_COLUMNS:
                    failures.append("ledger_template_not_zero_row_exact_schema")
    except (zipfile.BadZipFile, KeyError, UnicodeDecodeError, ValueError) as error:
        failures.append(f"invalid_zip:{type(error).__name__}")
    unique = tuple(dict.fromkeys(failures))
    return {
        "valid": not unique,
        "failure_codes": unique,
        "bundle_sha256": _sha256_bytes(bundle_bytes),
    }
