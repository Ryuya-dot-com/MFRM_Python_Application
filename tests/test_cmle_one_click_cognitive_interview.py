"""Contracts for private CMLE cognitive-interview operations."""

from __future__ import annotations

import copy
import json
from tempfile import TemporaryDirectory

import pandas as pd

from mfrm_app.cmle_one_click_cognitive_interview import (
    ISSUE_CODES,
    RESEARCHER_ONLY_FIELDS,
    build_blank_cognitive_record_template,
    build_cognitive_interview_schedule,
    build_moderator_probe_guide,
    build_participant_session_tasks,
    cognitive_interview_human_gate_status,
    validate_cognitive_interview_records,
)
from mfrm_app.cmle_one_click_comprehension import HUMAN_STUDY_STATUS
from validation.cmle_one_click_comprehension_readiness import run_cases
from pathlib import Path


def _materials():
    with TemporaryDirectory() as temporary:
        _, _, participant, key, _, _, _ = run_cases(Path(temporary))
    cases = participant["CaseId"].drop_duplicates().tolist()
    schedule = build_cognitive_interview_schedule(cases)
    tasks = build_participant_session_tasks(participant, schedule)
    template = build_blank_cognitive_record_template(tasks)
    return participant, key, schedule, tasks, template


def _complete_synthetic(template: pd.DataFrame, key: pd.DataFrame) -> pd.DataFrame:
    result = template.copy()
    options = {
        (str(row["CaseId"]), str(row["Language"]), str(row["ItemId"])): json.loads(
            row["OptionsJSON"]
        )[0]["code"]
        for _, row in key.iterrows()
    }
    result["ResponseCode"] = result.apply(
        lambda row: options[(str(row["CaseId"]), str(row["Language"]), str(row["ItemId"]))],
        axis=1,
    )
    result["ResponseLocked"] = True
    result["ProbeStartedAfterLock"] = True
    result["ResponseTimeSeconds"] = 12.5
    result["Confidence"] = 3
    result["IssueCodesJSON"] = "[]"
    result["IssueSeverity"] = "none"
    result["PIIReviewed"] = False
    result["CompletionState"] = "completed"
    return result


def test_schedule_is_pair_order_stratum_and_language_balanced() -> None:
    _, _, schedule, _, _ = _materials()
    assert len(schedule) == 40
    assert schedule["SessionSlotId"].nunique() == 20
    assert schedule.groupby("SessionSlotId").size().eq(2).all()
    for _, language in schedule.groupby("Language"):
        assert language["SessionSlotId"].nunique() == 10
        unordered = language.groupby("SessionSlotId")["CaseId"].apply(
            lambda values: tuple(sorted(values))
        )
        assert unordered.nunique() == 10
        exposure = language.groupby("CaseId").size()
        assert exposure.eq(4).all()
        position = language.groupby(["CaseId", "CaseOrder"]).size().unstack(fill_value=0)
        assert position.eq(2).all().all()
        strata = language.groupby(["CaseId", "ExperienceStratum"]).size().unstack(fill_value=0)
        assert strata.eq(2).all().all()
    parity = schedule.pivot_table(
        index=["SessionOrder", "CaseOrder"], columns="Language", values="CaseId", aggfunc="first"
    )
    assert parity["en"].equals(parity["ja"])


def test_participant_packets_are_complete_and_answer_free() -> None:
    _, _, _, tasks, _ = _materials()
    assert len(tasks) == 400
    assert tasks.groupby("SessionSlotId").size().eq(20).all()
    assert tasks.groupby(["SessionSlotId", "CaseOrder"]).size().eq(10).all()
    assert tasks.groupby("SessionSlotId")["TaskPosition"].apply(list).apply(
        lambda values: values == list(range(1, 21))
    ).all()
    assert not RESEARCHER_ONLY_FIELDS.intersection(tasks.columns)
    assert tasks["PreviewFile"].str.endswith(".html").all()


def test_answer_leakage_fails_before_packet_construction() -> None:
    participant, _, schedule, _, _ = _materials()
    contaminated = participant.assign(CorrectOption="leaked")
    try:
        build_participant_session_tasks(contaminated, schedule)
    except ValueError as exc:
        assert "Researcher-only" in str(exc)
    else:
        raise AssertionError("Answer-bearing participant packet was accepted.")


def test_moderator_guide_is_separate_complete_and_post_response() -> None:
    _, key, _, tasks, _ = _materials()
    guide = build_moderator_probe_guide(key)
    assert len(guide) == 100
    assert guide["ModeratorProbe"].str.len().gt(0).all()
    assert guide["ProbeTiming"].eq("after_response_locked_only").all()
    assert guide["PrivacyReminder"].str.contains("PIIReviewed", regex=False).all()
    assert "CorrectOption" in guide and "CorrectOption" not in tasks


def test_blank_template_and_complete_synthetic_record_validate() -> None:
    _, key, schedule, _, template = _materials()
    blank = validate_cognitive_interview_records(
        template, key, schedule, template_only=True
    )
    assert blank["valid"] and blank["status"] == "valid_template"
    complete = validate_cognitive_interview_records(
        _complete_synthetic(template, key), key, schedule
    )
    assert complete["valid"] and complete["status"] == "valid_record"


def test_registered_privacy_assignment_and_response_mutations_fail_closed() -> None:
    _, key, schedule, _, template = _materials()
    base = _complete_synthetic(template, key)
    mutations = []

    pii = base.assign(Email="person@example.invalid")
    mutations.append((pii, {"forbidden_pii_column", "unexpected_column"}))

    assignment = base.copy()
    assignment.loc[0, "CaseId"] = "not_assigned_here"
    mutations.append((assignment, {"assignment_mismatch"}))

    duplicate = pd.concat([base, base.iloc[[0]]], ignore_index=True)
    mutations.append((duplicate, {"record_row_count", "duplicate_item"}))

    missing = base.iloc[:-1].copy()
    mutations.append((missing, {"record_row_count", "missing_item"}))

    response = base.copy()
    response.loc[0, "ResponseCode"] = "unknown_option"
    mutations.append((response, {"unknown_response_option"}))

    issue = base.copy()
    issue.loc[0, "IssueCodesJSON"] = json.dumps(["not_registered"])
    mutations.append((issue, {"unknown_issue_code"}))

    note = base.copy()
    note.loc[0, "ParaphrasedObservation"] = "Paraphrased usability observation."
    note.loc[0, "PIIReviewed"] = False
    mutations.append((note, {"unreviewed_note"}))

    probe = base.copy()
    probe.loc[0, "ProbeStartedAfterLock"] = False
    mutations.append((probe, {"probe_before_response_lock"}))

    numeric = base.copy()
    numeric.loc[0, "ResponseTimeSeconds"] = -1
    numeric.loc[1, "Confidence"] = 8
    mutations.append((numeric, {"invalid_response_time", "invalid_confidence"}))

    for records, expected_codes in mutations:
        result = validate_cognitive_interview_records(records, key, schedule)
        assert not result["valid"]
        assert expected_codes.issubset(set(result["failure_codes"]))


def test_controlled_issue_codes_and_human_gate_stay_private() -> None:
    assert len(ISSUE_CODES) == 12 and len(set(ISSUE_CODES)) == 12
    status = cognitive_interview_human_gate_status().iloc[0]
    assert status["HumanStudyStatus"] == HUMAN_STUDY_STATUS
    assert status["HumanParticipants"] == 0
    assert not bool(status["RecruitmentAuthorizedByRepository"])
    assert not bool(status["ComprehensionRateAvailable"])
    assert not bool(status["LanguageEquivalenceAvailable"])
    assert not bool(status["PublicSurfaceEnabled"])
