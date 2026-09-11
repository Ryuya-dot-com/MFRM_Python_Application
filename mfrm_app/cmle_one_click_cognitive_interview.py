"""Private cognitive-interview operations for the bilingual CMLE instrument.

The functions in this module prepare schedule and record templates. They do
not recruit people, grant ethics approval, or create human-study evidence.
"""

from __future__ import annotations

from collections.abc import Sequence
import json
import math
import re

import pandas as pd

from mfrm_app.cmle_one_click_comprehension import (
    HUMAN_STUDY_STATUS,
    ITEM_IDS,
    SUPPORTED_LANGUAGES,
)


COGNITIVE_INTERVIEW_SCHEMA_VERSION = "cmle_cognitive_interview_operations_v1"
INSTRUMENT_VERSION = "cmle_one_click_comprehension_v1_frozen_20260810"
EXPERIENCE_STRATA = ("novice_target", "experienced_target")
ISSUE_CODES = (
    "wording",
    "status_discrimination",
    "numeracy_rounding",
    "estimator_role",
    "uncertainty_scope",
    "anchor_interpretation",
    "privacy_sharing",
    "next_action",
    "navigation",
    "visual_access",
    "assistive_technology",
    "translation_cultural",
)
ISSUE_SEVERITIES = ("none", "observation", "minor", "major", "critical")
RESEARCHER_ONLY_FIELDS = {
    "CorrectOption",
    "Critical",
    "MisconceptionCode",
    "Rationale",
    "ModeratorProbe",
    "ProbeTiming",
    "SeverityGuidance",
}
FORBIDDEN_IDENTITY_FIELDS = {
    "name",
    "email",
    "phone",
    "address",
    "organization",
    "employeeid",
    "studentid",
    "ipaddress",
    "deviceid",
    "contact",
}
RECORD_COLUMNS = (
    "OperationsSchemaVersion",
    "InstrumentVersion",
    "SessionSlotId",
    "Language",
    "ExperienceStratum",
    "PairId",
    "CaseOrder",
    "CaseId",
    "ItemOrder",
    "ItemId",
    "ResponseCode",
    "ResponseLocked",
    "ProbeStartedAfterLock",
    "ResponseTimeSeconds",
    "Confidence",
    "IssueCodesJSON",
    "IssueSeverity",
    "ParaphrasedObservation",
    "PIIReviewed",
    "Withdrawn",
    "CompletionState",
)

_PROBES = {
    "calibration_readiness": (
        "What on the screen led you to decide whether calibration was ready?",
        "較正がreadyかどうかを、画面のどこから判断しましたか？",
    ),
    "optimization_attempted": (
        "What tells you whether optimization ran or stopped earlier?",
        "最適化が実行されたか、その前に停止したかを何から判断しましたか？",
    ),
    "person_scoring_performed": (
        "Did you see an actual Person result or only a downstream status? Explain.",
        "実際のPerson結果を見たのか、下流状態だけを見たのか説明してください。",
    ),
    "fit_decision_input": (
        "Which number would you use if the displayed and raw values differ, and why?",
        "表示値と生値が異なるとき、どちらを使いますか。理由も説明してください。",
    ),
    "rounding_vignette": (
        "Please explain the classification without rounding 1.5004 in your reasoning.",
        "1.5004を丸めずに、分類理由を説明してください。",
    ),
    "calibration_uncertainty": (
        "What source of uncertainty is absent from the displayed conditional SE?",
        "表示された条件付きSEに含まれない不確実性は何ですか？",
    ),
    "anchor_validity": (
        "What can exact enforcement establish, and what can it not establish?",
        "exactな固定で示せることと、示せないことは何ですか？",
    ),
    "public_sharing": (
        "What would still be required before this result could be shared publicly?",
        "この結果を公開共有する前に、何が必要ですか？",
    ),
    "next_action": (
        "Why is that action safer than silently switching to another estimator?",
        "別の推定量へ黙って切り替えるより、なぜその行動が安全ですか？",
    ),
    "person_estimator_role": (
        "How would you describe the Person estimate without calling it a JMLE Person MLE?",
        "JMLE Person MLEと呼ばずに、このPerson推定値をどう説明しますか？",
    ),
}


def _normalized_column(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def _is_blank(value: object) -> bool:
    if value is None or value is pd.NA:
        return True
    try:
        if bool(pd.isna(value)):
            return True
    except (TypeError, ValueError):
        pass
    return str(value).strip() == ""


def _is_true(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if _is_blank(value):
        return False
    return str(value).strip().lower() in {"true", "1", "yes"}


def build_cognitive_interview_schedule(
    case_ids: Sequence[str],
    *,
    languages: Sequence[str] = SUPPORTED_LANGUAGES,
    instrument_version: str = INSTRUMENT_VERSION,
) -> pd.DataFrame:
    """Build a balanced 10-slot complete-pair schedule per language."""

    cases = [str(value) for value in case_ids]
    if len(cases) != 5 or len(set(cases)) != 5:
        raise ValueError("Exactly five unique CaseIds are required.")
    requested_languages = [str(value) for value in languages]
    if requested_languages != list(SUPPORTED_LANGUAGES):
        raise ValueError(f"languages must preserve {SUPPORTED_LANGUAGES}.")
    novice_edges = {
        frozenset((0, 1)),
        frozenset((1, 2)),
        frozenset((2, 3)),
        frozenset((3, 4)),
        frozenset((4, 0)),
    }
    pairs = [(left, right) for left in range(5) for right in range(left + 1, 5)]
    rows = []
    for language in requested_languages:
        for session_order, (left, right) in enumerate(pairs, start=1):
            difference = (right - left) % 5
            first, second = (left, right) if difference in {1, 2} else (right, left)
            stratum = (
                "novice_target"
                if frozenset((left, right)) in novice_edges
                else "experienced_target"
            )
            slot_id = f"CI-{language.upper()}-{session_order:02d}"
            pair_id = f"PAIR-{session_order:02d}"
            for case_order, index in enumerate((first, second), start=1):
                rows.append(
                    {
                        "OperationsSchemaVersion": COGNITIVE_INTERVIEW_SCHEMA_VERSION,
                        "InstrumentVersion": str(instrument_version),
                        "SessionSlotId": slot_id,
                        "SessionOrder": session_order,
                        "Language": language,
                        "ExperienceStratum": stratum,
                        "PairId": pair_id,
                        "CaseOrder": case_order,
                        "CaseId": cases[index],
                    }
                )
    schedule = pd.DataFrame(rows)
    if len(schedule) != 40:
        raise RuntimeError("Balanced schedule did not produce 40 assignment rows.")
    return schedule


def build_participant_session_tasks(
    participant_tasks: pd.DataFrame,
    schedule: pd.DataFrame,
) -> pd.DataFrame:
    """Expand frozen participant tasks into 20 private session packets."""

    required_tasks = {
        "SchemaVersion",
        "CaseId",
        "Language",
        "Order",
        "ItemId",
        "ScenarioHeadline",
        "Question",
        "OptionsJSON",
        "Required",
    }
    required_schedule = {
        "OperationsSchemaVersion",
        "InstrumentVersion",
        "SessionSlotId",
        "SessionOrder",
        "Language",
        "ExperienceStratum",
        "PairId",
        "CaseOrder",
        "CaseId",
    }
    if not isinstance(participant_tasks, pd.DataFrame) or not required_tasks.issubset(
        participant_tasks.columns
    ):
        raise ValueError("participant_tasks lacks required frozen fields.")
    leakage = RESEARCHER_ONLY_FIELDS.intersection(participant_tasks.columns)
    if leakage:
        raise ValueError(f"Researcher-only fields leaked into participant tasks: {sorted(leakage)}")
    if not isinstance(schedule, pd.DataFrame) or not required_schedule.issubset(
        schedule.columns
    ):
        raise ValueError("schedule lacks required assignment fields.")
    merged = schedule.merge(
        participant_tasks,
        on=["CaseId", "Language"],
        how="left",
        validate="many_to_many",
    )
    if merged["ItemId"].isna().any():
        raise ValueError("At least one scheduled case lacks participant tasks.")
    merged["TaskPosition"] = (merged["CaseOrder"].astype(int) - 1) * len(
        ITEM_IDS
    ) + merged["Order"].astype(int)
    merged["PreviewFile"] = merged.apply(
        lambda row: f"previews/{row['CaseId']}_{row['Language']}.html", axis=1
    )
    columns = [
        "OperationsSchemaVersion",
        "InstrumentVersion",
        "SessionSlotId",
        "SessionOrder",
        "Language",
        "ExperienceStratum",
        "PairId",
        "CaseOrder",
        "CaseId",
        "PreviewFile",
        "TaskPosition",
        "SchemaVersion",
        "Order",
        "ItemId",
        "ScenarioHeadline",
        "Question",
        "OptionsJSON",
        "Required",
    ]
    result = merged[columns].sort_values(
        ["Language", "SessionOrder", "TaskPosition"], kind="stable"
    ).reset_index(drop=True)
    if RESEARCHER_ONLY_FIELDS.intersection(result.columns):
        raise RuntimeError("Researcher-only fields leaked into a session packet.")
    return result


def build_moderator_probe_guide(researcher_key: pd.DataFrame) -> pd.DataFrame:
    """Build a private answer-bearing guide used only after response lock."""

    required = {
        "CaseId",
        "Language",
        "Order",
        "ItemId",
        "CorrectOption",
        "Critical",
        "MisconceptionCode",
        "Rationale",
    }
    if not isinstance(researcher_key, pd.DataFrame) or not required.issubset(
        researcher_key.columns
    ):
        raise ValueError("researcher_key lacks required fields.")
    result = researcher_key.copy()
    result["ModeratorProbe"] = result.apply(
        lambda row: _PROBES[str(row["ItemId"])][0 if row["Language"] == "en" else 1],
        axis=1,
    )
    result["ProbeTiming"] = "after_response_locked_only"
    result["SeverityGuidance"] = result["Critical"].map(
        {
            True: "critical_if_dangerous_misconception_path_is_plausible",
            False: "code_observed_barrier_without_inventing_criticality",
        }
    )
    result["PrivacyReminder"] = (
        "Paraphrase; exclude direct identifiers; set PIIReviewed before retention."
    )
    return result


def build_blank_cognitive_record_template(
    participant_session_tasks: pd.DataFrame,
) -> pd.DataFrame:
    """Return a no-human-values operational record template."""

    required = {
        "OperationsSchemaVersion",
        "InstrumentVersion",
        "SessionSlotId",
        "Language",
        "ExperienceStratum",
        "PairId",
        "CaseOrder",
        "CaseId",
        "Order",
        "ItemId",
    }
    if not isinstance(participant_session_tasks, pd.DataFrame) or not required.issubset(
        participant_session_tasks.columns
    ):
        raise ValueError("participant_session_tasks lacks required fields.")
    if RESEARCHER_ONLY_FIELDS.intersection(participant_session_tasks.columns):
        raise ValueError("Participant session tasks contain researcher-only fields.")
    result = participant_session_tasks[
        [
            "OperationsSchemaVersion",
            "InstrumentVersion",
            "SessionSlotId",
            "Language",
            "ExperienceStratum",
            "PairId",
            "CaseOrder",
            "CaseId",
            "Order",
            "ItemId",
        ]
    ].rename(columns={"Order": "ItemOrder"})
    result = result.assign(
        ResponseCode="",
        ResponseLocked=False,
        ProbeStartedAfterLock=False,
        ResponseTimeSeconds="",
        Confidence="",
        IssueCodesJSON="",
        IssueSeverity="none",
        ParaphrasedObservation="",
        PIIReviewed=False,
        Withdrawn=False,
        CompletionState="not_started",
    )
    return result[list(RECORD_COLUMNS)].reset_index(drop=True)


def _options_by_identity(researcher_key: pd.DataFrame) -> dict[tuple[str, str, str], set[str]]:
    options = {}
    for _, row in researcher_key.iterrows():
        parsed = json.loads(str(row["OptionsJSON"]))
        options[(str(row["CaseId"]), str(row["Language"]), str(row["ItemId"]))] = {
            str(item["code"]) for item in parsed
        }
    return options


def validate_cognitive_interview_records(
    records: pd.DataFrame,
    researcher_key: pd.DataFrame,
    schedule: pd.DataFrame,
    *,
    template_only: bool = False,
) -> dict[str, object]:
    """Fail closed on schema, assignment, response, timing, and privacy errors."""

    failures: list[dict[str, str]] = []

    def fail(code: str, location: object, detail: str) -> None:
        failures.append({"Code": code, "Location": str(location), "Detail": detail})

    if not isinstance(records, pd.DataFrame):
        raise ValueError("records must be a DataFrame.")
    observed_columns = list(records.columns)
    missing_columns = [column for column in RECORD_COLUMNS if column not in observed_columns]
    extra_columns = [column for column in observed_columns if column not in RECORD_COLUMNS]
    for column in missing_columns:
        fail("missing_required_column", column, "Required record column is absent.")
    for column in extra_columns:
        normalized = _normalized_column(column)
        if normalized in FORBIDDEN_IDENTITY_FIELDS:
            fail("forbidden_pii_column", column, "Direct/contact identifier field is forbidden.")
        fail("unexpected_column", column, "Record schema is an exact allowlist.")
    if missing_columns:
        return _validation_result(failures, template_only=template_only)

    assignment_columns = [
        "InstrumentVersion",
        "SessionSlotId",
        "Language",
        "ExperienceStratum",
        "PairId",
        "CaseOrder",
        "CaseId",
    ]
    expected = schedule[assignment_columns].merge(
        researcher_key[["CaseId", "Language", "Order", "ItemId"]].rename(
            columns={"Order": "ItemOrder"}
        ),
        on=["CaseId", "Language"],
        validate="many_to_many",
    )
    identity = ["SessionSlotId", "CaseOrder", "ItemId"]
    if len(records) != len(expected):
        fail("record_row_count", "records", f"Expected {len(expected)} rows; observed {len(records)}.")
    duplicate = records.duplicated(identity, keep=False)
    if duplicate.any():
        fail("duplicate_item", "records", f"Duplicate identity rows: {int(duplicate.sum())}.")
    actual_identity = set(map(tuple, records[identity].astype(str).to_numpy()))
    expected_identity = set(map(tuple, expected[identity].astype(str).to_numpy()))
    if expected_identity - actual_identity:
        fail("missing_item", "records", f"Missing identities: {len(expected_identity - actual_identity)}.")
    if actual_identity - expected_identity:
        fail("unregistered_assignment", "records", f"Unexpected identities: {len(actual_identity - expected_identity)}.")

    expected_lookup = expected.copy()
    for column in identity:
        expected_lookup[column] = expected_lookup[column].astype(str)
    expected_map = expected_lookup.set_index(identity)
    for index, row in records.iterrows():
        row_identity = (
            str(row["SessionSlotId"]),
            str(row["CaseOrder"]),
            str(row["ItemId"]),
        )
        if row_identity in expected_map.index:
            expected_row = expected_map.loc[row_identity]
            if isinstance(expected_row, pd.DataFrame):
                fail("duplicate_expected_assignment", index, "Schedule/key identity is not unique.")
            else:
                comparison_columns = [
                    column for column in assignment_columns if column not in identity
                ] + ["ItemOrder"]
                for column in comparison_columns:
                    if str(row[column]) != str(expected_row[column]):
                        fail("assignment_mismatch", index, f"{column} differs from registered assignment.")

    if template_only:
        human_value_columns = [
            "ResponseCode",
            "ResponseTimeSeconds",
            "Confidence",
            "IssueCodesJSON",
            "ParaphrasedObservation",
        ]
        for column in human_value_columns:
            if records[column].map(lambda value: not _is_blank(value)).any():
                fail("template_contains_human_value", column, "Blank template contains a response/observation value.")
        if not records["CompletionState"].astype(str).eq("not_started").all():
            fail("template_completion_state", "CompletionState", "Template must remain not_started.")
        if records["ResponseLocked"].map(_is_true).any() or records[
            "ProbeStartedAfterLock"
        ].map(_is_true).any():
            fail("template_interaction_state", "records", "Blank template cannot contain locked/probed states.")
        return _validation_result(failures, template_only=True)

    option_map = _options_by_identity(researcher_key)
    for index, row in records.iterrows():
        state = str(row["CompletionState"]).strip()
        withdrawn = _is_true(row["Withdrawn"])
        if state not in {"completed", "withdrawn"}:
            fail("invalid_completion_state", index, "Records must be completed or withdrawn.")
            continue
        if state == "withdrawn":
            if not withdrawn:
                fail("withdrawal_state_mismatch", index, "Withdrawn state requires Withdrawn=true.")
            for column in ("ResponseCode", "IssueCodesJSON", "ParaphrasedObservation"):
                if not _is_blank(row[column]):
                    fail("withdrawn_content_retained", index, f"{column} must be blank in this template contract.")
            continue
        if withdrawn:
            fail("withdrawal_state_mismatch", index, "Completed state cannot set Withdrawn=true.")
        response = str(row["ResponseCode"]).strip()
        key = (str(row["CaseId"]), str(row["Language"]), str(row["ItemId"]))
        if response not in option_map.get(key, set()):
            fail("unknown_response_option", index, "ResponseCode is not registered for this item.")
        if not _is_true(row["ResponseLocked"]):
            fail("response_not_locked", index, "Completed response must be locked before probing.")
        if not _is_true(row["ProbeStartedAfterLock"]):
            fail("probe_before_response_lock", index, "Probe timing was not confirmed after response lock.")
        try:
            seconds = float(row["ResponseTimeSeconds"])
            if not math.isfinite(seconds) or seconds < 0:
                raise ValueError
        except (TypeError, ValueError):
            fail("invalid_response_time", index, "Response time must be finite and nonnegative.")
        try:
            confidence = float(row["Confidence"])
            if not math.isfinite(confidence) or confidence not in {1, 2, 3, 4, 5}:
                raise ValueError
        except (TypeError, ValueError):
            fail("invalid_confidence", index, "Confidence must be an integer from 1 to 5.")
        issues: list[object] = []
        try:
            issues = json.loads(str(row["IssueCodesJSON"]))
            if not isinstance(issues, list):
                raise ValueError
            unknown = [str(code) for code in issues if str(code) not in ISSUE_CODES]
            if unknown:
                fail("unknown_issue_code", index, f"Unknown issue codes: {unknown}.")
        except (json.JSONDecodeError, ValueError):
            fail("invalid_issue_codes_json", index, "IssueCodesJSON must be a list of registered codes.")
        severity = str(row["IssueSeverity"]).strip()
        if severity not in ISSUE_SEVERITIES:
            fail("invalid_issue_severity", index, "IssueSeverity is not registered.")
        elif (not issues and severity != "none") or (issues and severity == "none"):
            fail(
                "issue_severity_mismatch",
                index,
                "IssueSeverity must be none exactly when the issue-code list is empty.",
            )
        if not _is_blank(row["ParaphrasedObservation"]) and not _is_true(row["PIIReviewed"]):
            fail("unreviewed_note", index, "Non-empty paraphrase requires PIIReviewed=true.")
    return _validation_result(failures, template_only=False)


def _validation_result(
    failures: list[dict[str, str]], *, template_only: bool
) -> dict[str, object]:
    failure_frame = pd.DataFrame(failures, columns=["Code", "Location", "Detail"])
    codes = tuple(dict.fromkeys(failure_frame["Code"].tolist())) if not failure_frame.empty else ()
    return {
        "valid": not failures,
        "status": (
            "valid_template"
            if template_only and not failures
            else "valid_record"
            if not template_only and not failures
            else "rejected"
        ),
        "failure_codes": codes,
        "failures": failure_frame,
    }


def cognitive_interview_human_gate_status() -> pd.DataFrame:
    """Return the deliberately non-promotable state before recruitment."""

    return pd.DataFrame(
        [
            {
                "OperationsSchemaVersion": COGNITIVE_INTERVIEW_SCHEMA_VERSION,
                "InstrumentVersion": INSTRUMENT_VERSION,
                "HumanStudyStatus": HUMAN_STUDY_STATUS,
                "PlannedSessionSlots": 20,
                "HumanParticipants": 0,
                "RecruitmentAuthorizedByRepository": False,
                "ComprehensionRateAvailable": False,
                "LanguageEquivalenceAvailable": False,
                "PublicSurfaceEnabled": False,
                "NextGate": "applicable_ethics_approval_then_bilingual_cognitive_interviews",
            }
        ]
    )


__all__ = [
    "COGNITIVE_INTERVIEW_SCHEMA_VERSION",
    "EXPERIENCE_STRATA",
    "FORBIDDEN_IDENTITY_FIELDS",
    "INSTRUMENT_VERSION",
    "ISSUE_CODES",
    "ISSUE_SEVERITIES",
    "RECORD_COLUMNS",
    "RESEARCHER_ONLY_FIELDS",
    "build_blank_cognitive_record_template",
    "build_cognitive_interview_schedule",
    "build_moderator_probe_guide",
    "build_participant_session_tasks",
    "cognitive_interview_human_gate_status",
    "validate_cognitive_interview_records",
]
