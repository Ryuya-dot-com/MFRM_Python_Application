"""Versioned directional confirmatory gates for the private CMLE instrument.

This module supplies identity, assignment, and statistical mechanics only. It
does not choose a confirmatory sample size, collect human data, or qualify a
public user interface.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist

import pandas as pd

from mfrm_app.cmle_one_click_cognitive_interview import INSTRUMENT_VERSION
from mfrm_app.cmle_one_click_comprehension import HUMAN_STUDY_STATUS, SUPPORTED_LANGUAGES


CONFIRMATORY_GATE_SCHEMA_VERSION = "cmle_directional_confirmatory_gate_v1"
ANCHOR_PRIMARY_CASE = "pcm_interior_differential_hard_anchors"
BLOCKED_PRIMARY_CASES = (
    "binary_conditional_support_boundary",
    "rsm_structurally_disconnected",
    "invalid_missing_score_column",
)
CONFIRMATORY_RESPONSE_COLUMNS = (
    "InstrumentVersion",
    "ConfirmatorySlotId",
    "Language",
    "CaseId",
    "ItemId",
    "ResponseCode",
)
DIRECTIONAL_RULES = (
    {
        "ItemId": "calibration_readiness",
        "Domain": "false_green_calibration",
        "PrimaryCaseRole": "blocked_primary",
        "EligibleCorrectOption": "not_ready",
        "DangerousResponseCodes": ("ready",),
    },
    {
        "ItemId": "person_scoring_performed",
        "Domain": "invented_person_scores",
        "PrimaryCaseRole": "blocked_primary",
        "EligibleCorrectOption": "not_performed",
        "DangerousResponseCodes": ("performed",),
    },
    {
        "ItemId": "rounding_vignette",
        "Domain": "rounded_mnsq_reclassification",
        "PrimaryCaseRole": "anchored_primary",
        "EligibleCorrectOption": "raw_noisy",
        "DangerousResponseCodes": ("display_acceptable",),
    },
    {
        "ItemId": "calibration_uncertainty",
        "Domain": "ignored_calibration_uncertainty",
        "PrimaryCaseRole": "anchored_primary",
        "EligibleCorrectOption": "not_propagated",
        "DangerousResponseCodes": ("propagated",),
    },
    {
        "ItemId": "anchor_validity",
        "Domain": "anchor_validity_overclaim",
        "PrimaryCaseRole": "anchored_primary",
        "EligibleCorrectOption": "does_not_validate",
        "DangerousResponseCodes": ("validates",),
    },
    {
        "ItemId": "public_sharing",
        "Domain": "private_archive_public_sharing",
        "PrimaryCaseRole": "anchored_primary",
        "EligibleCorrectOption": "not_permitted",
        "DangerousResponseCodes": ("permitted",),
    },
)


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _canonical_json_bytes(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, ensure_ascii=False, separators=(",", ":")
    ).encode("utf-8")


def _composite_identity(components: Sequence[tuple[str, bytes]]) -> str:
    digest = hashlib.sha256()
    for label, payload in components:
        label_bytes = label.encode("utf-8")
        digest.update(len(label_bytes).to_bytes(8, "big"))
        digest.update(label_bytes)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return digest.hexdigest()


def build_preview_identity_manifest(preview_dir: Path) -> pd.DataFrame:
    """Hash the frozen preview bytes in deterministic relative-path order."""

    path = Path(preview_dir)
    files = sorted(item for item in path.rglob("*.html") if item.is_file())
    if len(files) != 10:
        raise ValueError("Exactly ten frozen HTML previews are required.")
    return pd.DataFrame(
        [
            {
                "RelativePath": item.relative_to(path).as_posix(),
                "ByteCount": item.stat().st_size,
                "SHA256": _sha256_bytes(item.read_bytes()),
            }
            for item in files
        ]
    )


def compute_instrument_content_identity(
    participant_task_bytes: bytes,
    researcher_key_bytes: bytes,
    preview_manifest: pd.DataFrame,
    *,
    directional_rules: Sequence[Mapping[str, object]] = DIRECTIONAL_RULES,
) -> dict[str, str]:
    """Return component and composite identities for one frozen instrument."""

    required = {"RelativePath", "ByteCount", "SHA256"}
    if not isinstance(preview_manifest, pd.DataFrame) or not required.issubset(
        preview_manifest.columns
    ):
        raise ValueError("preview_manifest lacks required identity fields.")
    preview_records = preview_manifest.sort_values("RelativePath", kind="stable")[
        ["RelativePath", "ByteCount", "SHA256"]
    ].to_dict("records")
    rule_records = [
        {
            "ItemId": str(rule["ItemId"]),
            "Domain": str(rule["Domain"]),
            "PrimaryCaseRole": str(rule["PrimaryCaseRole"]),
            "EligibleCorrectOption": str(rule["EligibleCorrectOption"]),
            "DangerousResponseCodes": sorted(
                str(value) for value in rule["DangerousResponseCodes"]
            ),
        }
        for rule in directional_rules
    ]
    preview_bytes = _canonical_json_bytes(preview_records)
    rules_bytes = _canonical_json_bytes(rule_records)
    components = [
        ("participant_tasks", bytes(participant_task_bytes)),
        ("researcher_key", bytes(researcher_key_bytes)),
        ("preview_manifest", preview_bytes),
        ("directional_rules", rules_bytes),
    ]
    return {
        "ParticipantTasksSHA256": _sha256_bytes(bytes(participant_task_bytes)),
        "ResearcherKeySHA256": _sha256_bytes(bytes(researcher_key_bytes)),
        "PreviewManifestSHA256": _sha256_bytes(preview_bytes),
        "DirectionalRulesSHA256": _sha256_bytes(rules_bytes),
        "InstrumentContentSHA256": _composite_identity(components),
    }


def build_instrument_version_record(
    participant_task_path: Path,
    researcher_key_path: Path,
    preview_dir: Path,
    *,
    instrument_version: str = INSTRUMENT_VERSION,
    frozen_date: str = "2026-08-10",
    parent_version: str = "",
    change_class: str = "initial_frozen_instrument",
    change_summary: str = "Initial frozen bilingual instrument; no human responses.",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build a one-row no-human-data instrument-version ledger."""

    preview_manifest = build_preview_identity_manifest(Path(preview_dir))
    identity = compute_instrument_content_identity(
        Path(participant_task_path).read_bytes(),
        Path(researcher_key_path).read_bytes(),
        preview_manifest,
    )
    ledger = pd.DataFrame(
        [
            {
                "GateSchemaVersion": CONFIRMATORY_GATE_SCHEMA_VERSION,
                "InstrumentVersion": str(instrument_version),
                "FrozenDate": str(frozen_date),
                "ParentVersion": str(parent_version),
                "ChangeClass": str(change_class),
                "ChangeSummary": str(change_summary),
                **identity,
                "FreezeStatus": "frozen_before_human_responses",
                "HumanResponsesBeforeFreeze": 0,
                "PoolingWithDifferentContentAllowed": False,
                "PublicSurfaceEnabled": False,
            }
        ]
    )
    return ledger, preview_manifest


def validate_instrument_version_ledger(ledger: pd.DataFrame) -> dict[str, object]:
    """Reject ambiguous versions, post-response freezes, and pooling permission."""

    required = {
        "InstrumentVersion",
        "InstrumentContentSHA256",
        "FreezeStatus",
        "HumanResponsesBeforeFreeze",
        "PoolingWithDifferentContentAllowed",
        "PublicSurfaceEnabled",
    }
    failures = []
    if not isinstance(ledger, pd.DataFrame) or not required.issubset(ledger.columns):
        return {"valid": False, "failure_codes": ("invalid_ledger_schema",)}
    if ledger["InstrumentVersion"].astype(str).duplicated().any():
        failures.append("duplicate_version_label")
    if ledger["InstrumentContentSHA256"].astype(str).duplicated().any():
        failures.append("duplicate_content_identity")
    if not ledger["FreezeStatus"].astype(str).eq("frozen_before_human_responses").all():
        failures.append("invalid_freeze_status")
    if not pd.to_numeric(ledger["HumanResponsesBeforeFreeze"], errors="coerce").eq(0).all():
        failures.append("responses_before_freeze")
    def truth(value: object) -> bool:
        if isinstance(value, bool):
            return value
        return str(value).strip().lower() in {"true", "1", "yes"}

    if ledger["PoolingWithDifferentContentAllowed"].map(truth).any():
        failures.append("cross_version_pooling_allowed")
    if ledger["PublicSurfaceEnabled"].map(truth).any():
        failures.append("public_surface_enabled")
    return {"valid": not failures, "failure_codes": tuple(dict.fromkeys(failures))}


def build_confirmatory_assignment_schedule(
    slots_per_language: int,
    *,
    instrument_version: str = INSTRUMENT_VERSION,
    languages: Sequence[str] = SUPPORTED_LANGUAGES,
) -> pd.DataFrame:
    """Assign one anchored and one rotating blocked primary case per slot."""

    if not isinstance(slots_per_language, int) or slots_per_language <= 0:
        raise ValueError("slots_per_language must be a positive integer.")
    if list(languages) != list(SUPPORTED_LANGUAGES):
        raise ValueError(f"languages must preserve {SUPPORTED_LANGUAGES}.")
    rows = []
    for language in languages:
        for slot_number in range(1, slots_per_language + 1):
            blocked = BLOCKED_PRIMARY_CASES[(slot_number - 1) % len(BLOCKED_PRIMARY_CASES)]
            roles = (
                (("anchored_primary", ANCHOR_PRIMARY_CASE), ("blocked_primary", blocked))
                if slot_number % 2 == 1
                else (("blocked_primary", blocked), ("anchored_primary", ANCHOR_PRIMARY_CASE))
            )
            for case_order, (role, case_id) in enumerate(roles, start=1):
                rows.append(
                    {
                        "GateSchemaVersion": CONFIRMATORY_GATE_SCHEMA_VERSION,
                        "InstrumentVersion": str(instrument_version),
                        "ConfirmatorySlotId": f"CF-{language.upper()}-{slot_number:04d}",
                        "SlotNumber": slot_number,
                        "Language": language,
                        "CaseOrder": case_order,
                        "CaseRole": role,
                        "CaseId": case_id,
                    }
                )
    return pd.DataFrame(rows)


def wilson_upper_bound(
    errors: int,
    trials: int,
    *,
    confidence: float = 0.95,
) -> float:
    """Compute a one-sided Wilson upper confidence bound."""

    if not isinstance(errors, int) or not isinstance(trials, int):
        raise ValueError("errors and trials must be integers.")
    if trials <= 0 or errors < 0 or errors > trials:
        raise ValueError("Require 0 <= errors <= trials and trials > 0.")
    if not 0.5 < confidence < 1.0:
        raise ValueError("confidence must lie strictly between 0.5 and 1.")
    z = NormalDist().inv_cdf(confidence)
    proportion = errors / trials
    z2 = z * z
    denominator = 1.0 + z2 / trials
    center = (proportion + z2 / (2.0 * trials)) / denominator
    half_width = (
        z
        / denominator
        * math.sqrt(
            proportion * (1.0 - proportion) / trials
            + z2 / (4.0 * trials * trials)
        )
    )
    return min(1.0, center + half_width)


def maximum_passing_errors(
    trials: int,
    *,
    threshold: float = 0.10,
    confidence: float = 0.95,
) -> int:
    """Return the largest integer error count whose raw bound is < threshold."""

    if not 0.0 < threshold < 1.0:
        raise ValueError("threshold must lie strictly between zero and one.")
    passing = [
        errors
        for errors in range(trials + 1)
        if wilson_upper_bound(errors, trials, confidence=confidence) < threshold
    ]
    return max(passing, default=-1)


def build_wilson_planning_table(
    candidate_n: Sequence[int] = (24, 25, 30, 40, 50, 75, 100, 150, 200),
    *,
    alpha: float = 0.05,
    cell_count: int = 12,
    threshold: float = 0.10,
) -> pd.DataFrame:
    """Build arithmetic planning values, not a sample-size recommendation."""

    primary_confidence = 1.0 - alpha
    familywise_confidence = 1.0 - alpha / cell_count
    rows = []
    for value in candidate_n:
        n = int(value)
        if n <= 0 or n != value:
            raise ValueError("candidate_n values must be positive integers.")
        rows.append(
            {
                "ValidEligibleSlots": n,
                "ThresholdRaw": threshold,
                "PrimaryConfidence": primary_confidence,
                "PrimaryZeroErrorUpper": wilson_upper_bound(0, n, confidence=primary_confidence),
                "PrimaryMaximumPassingErrors": maximum_passing_errors(
                    n, threshold=threshold, confidence=primary_confidence
                ),
                "FamilywiseSensitivityConfidence": familywise_confidence,
                "FamilywiseZeroErrorUpper": wilson_upper_bound(
                    0, n, confidence=familywise_confidence
                ),
                "FamilywiseMaximumPassingErrors": maximum_passing_errors(
                    n, threshold=threshold, confidence=familywise_confidence
                ),
                "Interpretation": "arithmetic_only_not_sample_size_recommendation",
            }
        )
    return pd.DataFrame(rows)


def _option_codes(options_json: object) -> set[str]:
    parsed = json.loads(str(options_json))
    if not isinstance(parsed, list):
        raise ValueError("OptionsJSON must be a list.")
    return {str(value["code"]) for value in parsed}


def evaluate_directional_confirmatory_gate(
    responses: pd.DataFrame,
    researcher_key: pd.DataFrame,
    assignment: pd.DataFrame,
    *,
    registered_instrument_version: str,
    minimum_valid_per_cell: int,
    threshold: float = 0.10,
    alpha: float = 0.05,
) -> dict[str, object]:
    """Evaluate 12 directional cells without pooling language or case roles."""

    if not isinstance(minimum_valid_per_cell, int) or minimum_valid_per_cell <= 0:
        raise ValueError("minimum_valid_per_cell must be prospectively set to a positive integer.")
    if set(responses.columns) != set(CONFIRMATORY_RESPONSE_COLUMNS):
        raise ValueError("responses must use the exact confirmatory response schema.")
    key_required = {"CaseId", "Language", "ItemId", "CorrectOption", "OptionsJSON"}
    if not key_required.issubset(researcher_key.columns):
        raise ValueError("researcher_key lacks required fields.")
    assignment_required = {
        "InstrumentVersion",
        "ConfirmatorySlotId",
        "Language",
        "CaseRole",
        "CaseId",
    }
    if not assignment_required.issubset(assignment.columns):
        raise ValueError("assignment lacks required fields.")
    if assignment["InstrumentVersion"].astype(str).nunique() != 1 or not assignment[
        "InstrumentVersion"
    ].astype(str).eq(str(registered_instrument_version)).all():
        raise ValueError("assignment does not reference exactly the registered instrument version.")
    if assignment.duplicated(["ConfirmatorySlotId", "CaseRole"]).any():
        raise ValueError("assignment contains duplicate case roles within a slot.")

    key = researcher_key[["CaseId", "Language", "ItemId", "CorrectOption", "OptionsJSON"]].copy()
    expected = assignment.merge(
        key, on=["CaseId", "Language"], how="left", validate="many_to_many"
    )
    if expected[["ItemId", "CorrectOption", "OptionsJSON"]].isna().any().any():
        raise ValueError("assignment references a case without a complete key.")
    response_work = responses.copy()
    for column in CONFIRMATORY_RESPONSE_COLUMNS:
        response_work[column] = response_work[column].astype(str)

    planned_slots = assignment["ConfirmatorySlotId"].astype(str).drop_duplicates().tolist()
    observed_slots = response_work["ConfirmatorySlotId"].drop_duplicates().tolist()
    all_slots = planned_slots + [slot for slot in observed_slots if slot not in planned_slots]
    slot_rows = []
    response_audit_rows = []
    expected_identity_columns = ["ConfirmatorySlotId", "CaseId", "ItemId"]
    for slot_id in all_slots:
        expected_slot = expected.loc[
            expected["ConfirmatorySlotId"].astype(str).eq(slot_id)
        ].copy()
        observed_slot = response_work.loc[
            response_work["ConfirmatorySlotId"].eq(slot_id)
        ].copy()
        reasons = []
        if expected_slot.empty:
            reasons.append("unregistered_slot")
        if observed_slot.empty:
            state = "not_attempted"
            slot_rows.append(
                {
                    "ConfirmatorySlotId": slot_id,
                    "Language": expected_slot.iloc[0]["Language"] if not expected_slot.empty else "",
                    "ExpectedRows": len(expected_slot),
                    "ObservedRows": 0,
                    "Attempted": False,
                    "ValidSlot": False,
                    "InvalidReasons": "",
                    "SlotState": state,
                }
            )
            continue
        if observed_slot.duplicated(expected_identity_columns).any():
            reasons.append("duplicate_response")
        expected_ids = set(map(tuple, expected_slot[expected_identity_columns].astype(str).to_numpy()))
        observed_ids = set(map(tuple, observed_slot[expected_identity_columns].astype(str).to_numpy()))
        if expected_ids - observed_ids:
            reasons.append("missing_response")
        if observed_ids - expected_ids:
            reasons.append("unregistered_response_identity")
        if not observed_slot["InstrumentVersion"].eq(str(registered_instrument_version)).all():
            reasons.append("mixed_or_unregistered_instrument_version")
        if not expected_slot.empty:
            expected_language = str(expected_slot.iloc[0]["Language"])
            if not observed_slot["Language"].eq(expected_language).all():
                reasons.append("language_assignment_mismatch")
        joined = expected_slot.merge(
            observed_slot.drop_duplicates(expected_identity_columns, keep="first"),
            on=expected_identity_columns,
            how="left",
            suffixes=("Expected", "Observed"),
        )
        unknown_option = False
        if not joined.empty:
            for _, row in joined.iterrows():
                response_code = str(row.get("ResponseCode", ""))
                valid_option = response_code in _option_codes(row["OptionsJSON"])
                unknown_option = unknown_option or not valid_option
                response_audit_rows.append(
                    {
                        "ConfirmatorySlotId": slot_id,
                        "Language": row["LanguageExpected"],
                        "CaseId": row["CaseId"],
                        "CaseRole": row["CaseRole"],
                        "ItemId": row["ItemId"],
                        "CorrectOption": row["CorrectOption"],
                        "ResponseCode": response_code,
                        "ValidOption": valid_option,
                        "Correct": valid_option and response_code == str(row["CorrectOption"]),
                    }
                )
        if unknown_option:
            reasons.append("unknown_response_option")
        reasons = list(dict.fromkeys(reasons))
        valid = not reasons and len(observed_slot) == len(expected_slot)
        if len(observed_slot) != len(expected_slot) and not {
            "missing_response",
            "duplicate_response",
            "unregistered_response_identity",
        }.intersection(reasons):
            reasons.append("response_row_count")
            valid = False
        slot_rows.append(
            {
                "ConfirmatorySlotId": slot_id,
                "Language": expected_slot.iloc[0]["Language"] if not expected_slot.empty else observed_slot.iloc[0]["Language"],
                "ExpectedRows": len(expected_slot),
                "ObservedRows": len(observed_slot),
                "Attempted": True,
                "ValidSlot": valid,
                "InvalidReasons": ";".join(reasons),
                "SlotState": "valid" if valid else "invalid",
            }
        )

    slot_audit = pd.DataFrame(slot_rows)
    response_audit = pd.DataFrame(
        response_audit_rows,
        columns=[
            "ConfirmatorySlotId",
            "Language",
            "CaseId",
            "CaseRole",
            "ItemId",
            "CorrectOption",
            "ResponseCode",
            "ValidOption",
            "Correct",
        ],
    )
    valid_slots = set(slot_audit.loc[slot_audit["ValidSlot"], "ConfirmatorySlotId"].astype(str))
    response_audit["ValidSlot"] = response_audit["ConfirmatorySlotId"].astype(str).isin(valid_slots)
    response_audit["DirectionalPrimary"] = False
    response_audit["DangerousError"] = False
    response_audit["DirectionalDomain"] = ""
    cell_rows = []
    primary_confidence = 1.0 - alpha
    familywise_confidence = 1.0 - alpha / (len(DIRECTIONAL_RULES) * len(SUPPORTED_LANGUAGES))
    for language in SUPPORTED_LANGUAGES:
        language_slots = slot_audit.loc[slot_audit["Language"].astype(str).eq(language)]
        for rule in DIRECTIONAL_RULES:
            mask = (
                response_audit["Language"].astype(str).eq(language)
                & response_audit["CaseRole"].astype(str).eq(str(rule["PrimaryCaseRole"]))
                & response_audit["ItemId"].astype(str).eq(str(rule["ItemId"]))
            )
            target = response_audit.loc[mask].copy()
            if not target.empty and not target["CorrectOption"].astype(str).eq(
                str(rule["EligibleCorrectOption"])
            ).all():
                raise ValueError(f"Directional rule/key mismatch for {rule['Domain']} in {language}.")
            response_audit.loc[mask, "DirectionalPrimary"] = True
            response_audit.loc[mask, "DirectionalDomain"] = str(rule["Domain"])
            dangerous_codes = {str(value) for value in rule["DangerousResponseCodes"]}
            dangerous_mask = mask & response_audit["ResponseCode"].astype(str).isin(dangerous_codes)
            response_audit.loc[dangerous_mask, "DangerousError"] = True
            valid_target = target.loc[target["ValidSlot"]]
            n = len(valid_target)
            dangerous = int(valid_target["ResponseCode"].astype(str).isin(dangerous_codes).sum())
            any_errors = int((~valid_target["Correct"]).sum())
            primary_upper = (
                wilson_upper_bound(dangerous, n, confidence=primary_confidence)
                if n
                else math.nan
            )
            familywise_upper = (
                wilson_upper_bound(dangerous, n, confidence=familywise_confidence)
                if n
                else math.nan
            )
            cell_rows.append(
                {
                    "Language": language,
                    "ItemId": rule["ItemId"],
                    "Domain": rule["Domain"],
                    "PrimaryCaseRole": rule["PrimaryCaseRole"],
                    "EligibleCorrectOption": rule["EligibleCorrectOption"],
                    "DangerousResponseCodesJSON": json.dumps(
                        list(rule["DangerousResponseCodes"]), separators=(",", ":")
                    ),
                    "PlannedSlots": len(language_slots),
                    "AttemptedSlots": int(language_slots["Attempted"].sum()),
                    "InvalidSlots": int((language_slots["Attempted"] & ~language_slots["ValidSlot"]).sum()),
                    "ValidEligibleSlots": n,
                    "AnyComprehensionErrors": any_errors,
                    "DangerousErrors": dangerous,
                    "DangerousErrorProportionRaw": dangerous / n if n else math.nan,
                    "PrimaryWilsonUpperRaw": primary_upper,
                    "PrimaryWilsonUpperDisplay6": round(primary_upper, 6) if n else math.nan,
                    "PrimaryPass": bool(
                        n >= minimum_valid_per_cell and primary_upper < threshold
                    ) if n else False,
                    "FamilywiseSensitivityUpperRaw": familywise_upper,
                    "FamilywiseSensitivityPass": bool(
                        n >= minimum_valid_per_cell and familywise_upper < threshold
                    ) if n else False,
                    "DecisionThresholdRaw": threshold,
                    "DecisionUsesRawBound": True,
                }
            )
    cell_summary = pd.DataFrame(cell_rows)
    overall_pass = bool(
        len(cell_summary) == 12
        and cell_summary["PrimaryPass"].all()
    )
    study_summary = pd.DataFrame(
        [
            {
                "GateSchemaVersion": CONFIRMATORY_GATE_SCHEMA_VERSION,
                "InstrumentVersion": registered_instrument_version,
                "PlannedSlots": int(slot_audit["ConfirmatorySlotId"].isin(planned_slots).sum()),
                "AttemptedSlots": int(slot_audit["Attempted"].sum()),
                "ValidSlots": int(slot_audit["ValidSlot"].sum()),
                "InvalidAttemptedSlots": int((slot_audit["Attempted"] & ~slot_audit["ValidSlot"]).sum()),
                "MinimumValidPerCell": minimum_valid_per_cell,
                "PrimaryCellCount": len(cell_summary),
                "PrimaryPassingCells": int(cell_summary["PrimaryPass"].sum()),
                "FamilywiseSensitivityPassingCells": int(
                    cell_summary["FamilywiseSensitivityPass"].sum()
                ),
                "OverallPrimaryGatePassed": overall_pass,
                "DecisionUsesUnroundedWilsonBound": True,
                "PublicSurfaceEnabled": False,
            }
        ]
    )
    return {
        "slot_audit": slot_audit,
        "response_audit": response_audit,
        "cell_summary": cell_summary,
        "study_summary": study_summary,
    }


def confirmatory_human_gate_status() -> pd.DataFrame:
    """Return the no-human-data state before pilot or confirmatory registration."""

    return pd.DataFrame(
        [
            {
                "GateSchemaVersion": CONFIRMATORY_GATE_SCHEMA_VERSION,
                "InstrumentVersion": INSTRUMENT_VERSION,
                "HumanStudyStatus": HUMAN_STUDY_STATUS,
                "HumanParticipants": 0,
                "PilotResultAvailable": False,
                "ConfirmatoryResultAvailable": False,
                "MinimumValidPerCellRegistered": False,
                "LanguageEquivalenceAvailable": False,
                "PublicSurfaceEnabled": False,
                "NextGate": "approved_cognitive_interviews_then_new_frozen_pilot_registration",
            }
        ]
    )


__all__ = [
    "ANCHOR_PRIMARY_CASE",
    "BLOCKED_PRIMARY_CASES",
    "CONFIRMATORY_GATE_SCHEMA_VERSION",
    "CONFIRMATORY_RESPONSE_COLUMNS",
    "DIRECTIONAL_RULES",
    "build_confirmatory_assignment_schedule",
    "build_instrument_version_record",
    "build_preview_identity_manifest",
    "build_wilson_planning_table",
    "compute_instrument_content_identity",
    "confirmatory_human_gate_status",
    "evaluate_directional_confirmatory_gate",
    "maximum_passing_errors",
    "validate_instrument_version_ledger",
    "wilson_upper_bound",
]
