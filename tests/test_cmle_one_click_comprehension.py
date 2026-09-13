"""Contracts for the non-public bilingual CMLE comprehension preparation."""

from __future__ import annotations

import copy
import json

import pandas as pd

from mfrm_app.cmle_one_click import run_cmle_one_click_analysis
from mfrm_app.cmle_one_click_comprehension import (
    HUMAN_STUDY_STATUS,
    ITEM_IDS,
    build_comprehension_materials,
    comprehension_human_gate_status,
    render_cmle_one_click_preview_html,
    score_comprehension_responses,
)


def _interior_frame() -> pd.DataFrame:
    scores = {
        "P1": [0, 1, 1, 2],
        "P2": [1, 0, 2, 1],
        "P3": [2, 1, 1, 0],
        "P4": [1, 2, 0, 1],
        "P5": [0, 2, 1, 1],
        "P6": [2, 0, 1, 1],
    }
    units = [("R1", "C1"), ("R1", "C2"), ("R2", "C1"), ("R2", "C2")]
    return pd.DataFrame(
        [
            (person, rater, criterion, score)
            for person, values in scores.items()
            for (rater, criterion), score in zip(units, values, strict=True)
        ],
        columns=["Person", "Rater", "Criterion", "Score"],
    )


def _boundary_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            (f"P{index:03d}", rater, score)
            for index in range(20)
            for rater, score in (("R1", 1), ("R2", 0))
        ],
        columns=["Person", "Rater", "Score"],
    )


def _run(frame: pd.DataFrame, *, binary=False):
    return run_cmle_one_click_analysis(
        frame,
        person_col="Person",
        facet_cols=["Rater"] if binary else ["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=1 if binary else 2,
        model="RSM",
        gtol=1e-8,
        maxiter=800,
        display_decimals=3,
    )


def _responses_from_key(key: pd.DataFrame, participant: str, *, wrong=False):
    rows = []
    for _, row in key.iterrows():
        codes = [item["code"] for item in json.loads(row["OptionsJSON"])]
        response = row["CorrectOption"]
        if wrong:
            response = next(code for code in codes if code != response)
        rows.append(
            {
                "ParticipantId": participant,
                "CaseId": row["CaseId"],
                "Language": row["Language"],
                "ItemId": row["ItemId"],
                "ResponseCode": response,
                "ResponseTimeSeconds": 3.0,
            }
        )
    return pd.DataFrame(rows)


def test_preview_has_six_textual_cards_no_script_and_escapes_content() -> None:
    result = _run(_interior_frame())
    for language in ("en", "ja"):
        preview = render_cmle_one_click_preview_html(
            result, language=language, case_id="ready"
        )
        assert f'<html lang="{language}">' in preview
        assert preview.count('<section class="result-card') == 6
        assert "data-public-surface=\"false\"" in preview
        assert "<script" not in preview.lower()
        assert "http://" not in preview and "https://" not in preview
        assert "status-text" in preview
    injected = copy.deepcopy(result)
    injected["cards"].loc[0, "HeadlineEn"] = "<script>alert(1)</script>"
    escaped = render_cmle_one_click_preview_html(
        injected, language="en", case_id="escape"
    )
    assert "<script>alert" not in escaped
    assert "&lt;script&gt;alert" in escaped


def test_participant_packet_has_no_answer_fields_and_language_key_parity() -> None:
    result = _run(_interior_frame())
    materials = [
        build_comprehension_materials(result, case_id="ready", language=language)
        for language in ("en", "ja")
    ]
    participant = pd.concat([item["participant_tasks"] for item in materials])
    key = pd.concat([item["researcher_key"] for item in materials])
    assert len(participant) == 20 and len(key) == 20
    assert participant.groupby("Language")["ItemId"].apply(list).apply(
        lambda values: values == list(ITEM_IDS)
    ).all()
    assert not {
        "CorrectOption",
        "Critical",
        "MisconceptionCode",
        "Rationale",
    }.intersection(participant.columns)
    pivot = key.pivot(index="ItemId", columns="Language", values="CorrectOption")
    assert pivot["en"].equals(pivot["ja"])
    critical = key.drop_duplicates("ItemId").set_index("ItemId")["Critical"]
    assert int(critical.sum()) == 6


def test_state_specific_answers_distinguish_ready_and_boundary() -> None:
    ready = build_comprehension_materials(
        _run(_interior_frame()), case_id="ready", language="en"
    )["researcher_key"].set_index("ItemId")
    boundary = build_comprehension_materials(
        _run(_boundary_frame(), binary=True), case_id="boundary", language="en"
    )["researcher_key"].set_index("ItemId")
    assert ready.loc["calibration_readiness", "CorrectOption"] == "ready"
    assert ready.loc["optimization_attempted", "CorrectOption"] == "attempted"
    assert ready.loc["person_scoring_performed", "CorrectOption"] == "performed"
    assert ready.loc["person_estimator_role", "CorrectOption"] == "fixed_calibration_wle"
    assert boundary.loc["calibration_readiness", "CorrectOption"] == "not_ready"
    assert boundary.loc["optimization_attempted", "CorrectOption"] == "not_attempted"
    assert boundary.loc["person_scoring_performed", "CorrectOption"] == "not_performed"
    assert boundary.loc["next_action", "CorrectOption"] == "review_boundary"
    assert boundary.loc["person_estimator_role", "CorrectOption"] == "not_computed"


def test_scoring_perfect_wrong_and_false_green_packets() -> None:
    key = build_comprehension_materials(
        _run(_boundary_frame(), binary=True), case_id="boundary", language="ja"
    )["researcher_key"]
    perfect = _responses_from_key(key, "perfect")
    wrong = _responses_from_key(key, "wrong", wrong=True)
    false_green = _responses_from_key(key, "false_green")
    false_green.loc[
        false_green["ItemId"].eq("calibration_readiness"), "ResponseCode"
    ] = "ready"
    scored = score_comprehension_responses(
        key, pd.concat([perfect, wrong, false_green], ignore_index=True)
    )["participant_summary"].set_index("ParticipantId")
    assert bool(scored.loc["perfect", "ParticipantCaseReady"])
    assert not bool(scored.loc["wrong", "ParticipantCaseReady"])
    assert scored.loc["wrong", "CriticalErrors"] == 6
    assert not bool(scored.loc["false_green", "ParticipantCaseReady"])
    assert scored.loc["false_green", "CriticalErrors"] == 1


def test_scoring_invalid_packets_fail_closed() -> None:
    key = build_comprehension_materials(
        _run(_interior_frame()), case_id="ready", language="en"
    )["researcher_key"]
    base = _responses_from_key(key, "base")
    packets = []
    duplicate = pd.concat([base.assign(ParticipantId="duplicate"), base.iloc[[0]].assign(ParticipantId="duplicate")])
    packets.append(duplicate)
    packets.append(base.iloc[:-1].assign(ParticipantId="missing"))
    unknown = base.assign(ParticipantId="unknown")
    unknown.loc[unknown.index[0], "ResponseCode"] = "not_an_option"
    packets.append(unknown)
    invalid_time = base.assign(ParticipantId="time")
    invalid_time.loc[invalid_time.index[0], "ResponseTimeSeconds"] = -1
    packets.append(invalid_time)
    scored = score_comprehension_responses(
        key, pd.concat(packets, ignore_index=True)
    )["participant_summary"].set_index("ParticipantId")
    assert not scored["ValidPacket"].any()
    assert not scored["ParticipantCaseReady"].any()
    assert "duplicate_response" in scored.loc["duplicate", "InvalidReasons"]
    assert "missing_response" in scored.loc["missing", "InvalidReasons"]
    assert "unknown_option" in scored.loc["unknown", "InvalidReasons"]
    assert "invalid_response_time" in scored.loc["time", "InvalidReasons"]


def test_human_gate_remains_not_started_and_public_false() -> None:
    status = comprehension_human_gate_status().iloc[0]
    assert status["HumanStudyStatus"] == HUMAN_STUDY_STATUS
    assert status["HumanParticipants"] == 0
    assert not bool(status["ComprehensionRateAvailable"])
    assert not bool(status["LanguageEquivalenceAvailable"])
    assert not bool(status["PublicSurfaceEnabled"])
