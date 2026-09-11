"""Non-public bilingual preview and comprehension-study preparation for CMLE.

Machine checks in this module validate instrument structure and scoring logic;
they do not constitute human usability, accessibility, translation-equivalence,
or comprehension evidence.
"""

from __future__ import annotations

from collections.abc import Mapping
import html
import json

import numpy as np
import pandas as pd

from mfrm_app.cmle_one_click import CMLE_ONE_CLICK_CARDS


CMLE_COMPREHENSION_SCHEMA_VERSION = "cmle_one_click_comprehension_v1"
HUMAN_STUDY_STATUS = "not_started_no_human_data"
SUPPORTED_LANGUAGES = ("en", "ja")
ITEM_IDS = (
    "calibration_readiness",
    "optimization_attempted",
    "person_scoring_performed",
    "fit_decision_input",
    "rounding_vignette",
    "calibration_uncertainty",
    "anchor_validity",
    "public_sharing",
    "next_action",
    "person_estimator_role",
)
_CRITICAL_ITEMS = {
    "calibration_readiness": "false_green_calibration",
    "person_scoring_performed": "invented_person_scores",
    "rounding_vignette": "rounded_mnsq_reclassification",
    "calibration_uncertainty": "ignored_calibration_uncertainty",
    "anchor_validity": "anchor_validity_overclaim",
    "public_sharing": "private_archive_public_sharing",
}

_STATUS_LABELS = {
    "en": {
        "ready": "Ready",
        "available": "Available",
        "blocked": "Stopped",
        "not_run": "Not run",
        "review": "Review",
        "caution": "Caution",
        "withheld": "Withheld",
    },
    "ja": {
        "ready": "利用可能",
        "available": "利用可能",
        "blocked": "停止",
        "not_run": "未実行",
        "review": "要確認",
        "caution": "注意",
        "withheld": "保留",
    },
}

_OPTIONS = {
    "calibration_readiness": (
        ("ready", "Yes, calibration is ready", "はい、較正はreadyです"),
        ("not_ready", "No, calibration is not ready", "いいえ、較正はreadyではありません"),
    ),
    "optimization_attempted": (
        ("attempted", "Optimization was attempted", "最適化を実行しました"),
        ("not_attempted", "Optimization was not attempted", "最適化を実行していません"),
    ),
    "person_scoring_performed": (
        ("performed", "Person scoring was performed", "Person scoringを実行しました"),
        ("not_performed", "Person scoring was not performed", "Person scoringを実行していません"),
    ),
    "fit_decision_input": (
        ("raw_unrounded", "Finite unrounded MnSq", "有限な丸め前MnSq"),
        ("displayed", "Displayed rounded MnSq", "表示用に丸めたMnSq"),
    ),
    "rounding_vignette": (
        ("raw_noisy", "Noisy, because raw Infit is 1.5004", "生のInfitが1.5004なのでnoisy"),
        ("display_acceptable", "Acceptable, because display is 1.500", "表示が1.500なのでacceptable"),
        ("cannot_decide", "Cannot decide", "判断できない"),
    ),
    "calibration_uncertainty": (
        ("not_propagated", "No, fitted CMLE calibration is treated as fixed", "いいえ、推定CMLE較正を固定値として扱います"),
        ("propagated", "Yes, it is included in the WLE SE", "はい、WLE SEに含まれます"),
    ),
    "anchor_validity": (
        ("does_not_validate", "No, exact enforcement does not validate an anchor", "いいえ、exactな固定はアンカーを妥当化しません"),
        ("validates", "Yes, exact reproduction proves the anchor is unbiased", "はい、exactな再現はアンカーが不偏である証拠です"),
        ("not_applicable_no_anchor", "No hard anchor is active in this case", "この条件ではhard anchorを使用していません"),
    ),
    "public_sharing": (
        ("not_permitted", "No, public sharing is withheld", "いいえ、公開共有は保留されています"),
        ("permitted", "Yes, the archive is ready for public sharing", "はい、アーカイブを公開共有できます"),
    ),
    "next_action": (
        ("correct_input", "Correct the input and rerun", "入力を修正して再実行する"),
        ("repair_design", "Review overlap/connectivity and design", "重なり・連結性・デザインを見直す"),
        ("review_boundary", "Review the boundary; choose any different estimator explicitly", "境界を確認し、別推定量なら明示的に選ぶ"),
        ("review_cautions", "Review raw fit, anchors, extremes, and withheld uncertainty", "生のfit、アンカー、extreme、保留中の不確実性を確認する"),
        ("repair_downstream", "Inspect the downstream scoring error", "下流scoringのエラーを確認する"),
    ),
    "person_estimator_role": (
        ("fixed_calibration_wle", "Warm WLE with fitted CMLE calibration held fixed", "推定CMLE較正を固定したWarm WLE"),
        ("jmle_person_mle", "A finite JMLE Person MLE", "有限なJMLE Person MLE"),
        ("not_computed", "No Person estimate was computed", "Person推定値は計算されていない"),
    ),
}

_QUESTIONS = {
    "calibration_readiness": (
        "Is the structural calibration ready for interpretation in this result?",
        "この結果では、構造較正を解釈できるready状態ですか？",
    ),
    "optimization_attempted": (
        "Was structural optimization attempted?",
        "構造最適化は実行されましたか？",
    ),
    "person_scoring_performed": (
        "Were fixed-calibration Person scores computed?",
        "固定較正Person得点は計算されましたか？",
    ),
    "fit_decision_input": (
        "Which value is allowed to drive an Infit/Outfit classification?",
        "Infit／Outfitの分類に使用してよい値はどれですか？",
    ),
    "rounding_vignette": (
        "Constructed vignette: raw Infit is 1.5004 and displayed Infit is 1.500. Which classification drives the decision?",
        "構成例：生のInfitは1.5004、表示Infitは1.500です。意思決定に使う分類はどれですか？",
    ),
    "calibration_uncertainty": (
        "Does the displayed conditional WLE SE propagate fitted-CMLE calibration uncertainty?",
        "表示された条件付きWLE SEには、推定CMLE較正の不確実性が伝播されていますか？",
    ),
    "anchor_validity": (
        "What does exact hard-anchor reproduction establish in this case?",
        "この条件でhard anchorをexactに再現することは何を示しますか？",
    ),
    "public_sharing": (
        "May this research/private result be treated as approved for public sharing?",
        "このresearch／private結果は、公開共有が承認済みと扱えますか？",
    ),
    "next_action": (
        "What is the most appropriate next action?",
        "最も適切な次の行動はどれですか？",
    ),
    "person_estimator_role": (
        "What is the role of the Person estimate in this result?",
        "この結果のPerson推定値はどのような位置づけですか？",
    ),
}

_RATIONALES = {
    "calibration_readiness": (
        "Readiness comes from the machine workflow status, not from the presence of finite-looking output.",
        "readinessは有限に見える出力の有無ではなく、機械workflowの状態で決まります。",
    ),
    "optimization_attempted": (
        "Input, design, and finite-existence gates can stop before optimization.",
        "入力、デザイン、有限存在の判定により、最適化前に停止する場合があります。",
    ),
    "person_scoring_performed": (
        "Person scoring is run only after calibration readiness.",
        "Person scoringは較正がreadyの場合にのみ実行されます。",
    ),
    "fit_decision_input": (
        "Displayed rounding is presentation only; finite raw MnSq drives classification.",
        "表示用丸めは提示のみで、有限な生MnSqが分類を決めます。",
    ),
    "rounding_vignette": (
        "Raw 1.5004 is above 1.50 and is noisy; displayed 1.500 must not overwrite it.",
        "生値1.5004は1.50を超えるためnoisyで、表示1.500で上書きしてはいけません。",
    ),
    "calibration_uncertainty": (
        "The conditional WLE SE treats fitted CMLE calibration as fixed.",
        "条件付きWLE SEは推定CMLE較正を固定値として扱います。",
    ),
    "anchor_validity": (
        "Exact implementation can enforce a supplied anchor but cannot prove its provenance or lack of bias.",
        "exact実装は指定アンカーを固定できますが、その根拠や不偏性は証明できません。",
    ),
    "public_sharing": (
        "The research preview and private archive are not de-identified or approved public artifacts.",
        "research previewとprivate archiveは、非識別化済み・公開承認済みの成果物ではありません。",
    ),
    "next_action": (
        "The next action follows the first blocking or caution state and never silently changes estimator.",
        "次の行動は最初の停止・注意状態に従い、推定量を黙って変更しません。",
    ),
    "person_estimator_role": (
        "Ready Person scores are Warm WLE conditional on fitted CMLE calibration, not JMLE Person MLEs.",
        "readyなPerson得点は推定CMLE較正を固定したWarm WLEで、JMLE Person MLEではありません。",
    ),
}


def _validate_result(result: Mapping[str, object]) -> tuple[pd.Series, pd.DataFrame, pd.Series]:
    if not isinstance(result, Mapping):
        raise ValueError("result must be a CMLE one-click result mapping.")
    summary = result.get("summary")
    cards = result.get("cards")
    calibration = result.get("calibration")
    if not isinstance(summary, pd.DataFrame) or len(summary) != 1:
        raise ValueError("One-click summary must contain one row.")
    if not isinstance(cards, pd.DataFrame) or cards["Card"].tolist() != list(
        CMLE_ONE_CLICK_CARDS
    ):
        raise ValueError("One-click cards do not match the six-card contract.")
    if not isinstance(calibration, Mapping):
        raise ValueError("Retained calibration workflow is missing.")
    calibration_summary = calibration.get("summary")
    if not isinstance(calibration_summary, pd.DataFrame) or len(calibration_summary) != 1:
        raise ValueError("Calibration workflow summary must contain one row.")
    return summary.iloc[0], cards, calibration_summary.iloc[0]


def _correct_answers(result: Mapping[str, object]) -> dict[str, str]:
    summary, _, calibration = _validate_result(result)
    terminal = str(summary["TerminalStatus"])
    calibration_ready = bool(summary["CalibrationReady"])
    person_ready = bool(summary["PersonScoringReady"])
    anchored = int(summary["AnchoredMainEffectLevelCount"])
    if terminal == "input_invalid":
        next_action = "correct_input"
    elif terminal == "design_not_identified":
        next_action = "repair_design"
    elif terminal in {"finite_mle_boundary", "finite_mle_unavailable"}:
        next_action = "review_boundary"
    elif terminal == "downstream_person_scoring_not_ready":
        next_action = "repair_downstream"
    else:
        next_action = "review_cautions"
    return {
        "calibration_readiness": "ready" if calibration_ready else "not_ready",
        "optimization_attempted": (
            "attempted" if bool(calibration["OptimizationAttempted"]) else "not_attempted"
        ),
        "person_scoring_performed": "performed" if person_ready else "not_performed",
        "fit_decision_input": "raw_unrounded",
        "rounding_vignette": "raw_noisy",
        "calibration_uncertainty": "not_propagated",
        "anchor_validity": (
            "does_not_validate" if anchored else "not_applicable_no_anchor"
        ),
        "public_sharing": "not_permitted",
        "next_action": next_action,
        "person_estimator_role": (
            "fixed_calibration_wle" if person_ready else "not_computed"
        ),
    }


def build_comprehension_materials(
    result: Mapping[str, object],
    *,
    case_id: str,
    language: str,
) -> dict[str, pd.DataFrame]:
    """Return participant tasks and a separately retained researcher key."""

    if language not in SUPPORTED_LANGUAGES:
        raise ValueError(f"language must be one of {SUPPORTED_LANGUAGES}.")
    summary, _, _ = _validate_result(result)
    answers = _correct_answers(result)
    lang_index = 0 if language == "en" else 1
    headline = str(summary["HeadlineEn" if language == "en" else "HeadlineJa"])
    task_rows = []
    key_rows = []
    for order, item_id in enumerate(ITEM_IDS, start=1):
        options = [
            {
                "code": row[0],
                "label": row[1] if language == "en" else row[2],
            }
            for row in _OPTIONS[item_id]
        ]
        question = _QUESTIONS[item_id][lang_index]
        common = {
            "SchemaVersion": CMLE_COMPREHENSION_SCHEMA_VERSION,
            "CaseId": str(case_id),
            "Language": language,
            "Order": order,
            "ItemId": item_id,
            "ScenarioHeadline": headline,
            "Question": question,
            "OptionsJSON": json.dumps(options, ensure_ascii=False, separators=(",", ":")),
            "Required": True,
        }
        task_rows.append(common)
        rationale = _RATIONALES[item_id][lang_index]
        key_rows.append(
            {
                **common,
                "CorrectOption": answers[item_id],
                "Critical": item_id in _CRITICAL_ITEMS,
                "MisconceptionCode": _CRITICAL_ITEMS.get(item_id, ""),
                "Rationale": rationale,
            }
        )
    participant = pd.DataFrame(task_rows)
    researcher = pd.DataFrame(key_rows)
    forbidden = {"CorrectOption", "Critical", "MisconceptionCode", "Rationale"}
    if forbidden.intersection(participant.columns):
        raise RuntimeError("Researcher answer fields leaked into participant tasks.")
    return {"participant_tasks": participant, "researcher_key": researcher}


def render_cmle_one_click_preview_html(
    result: Mapping[str, object],
    *,
    language: str,
    case_id: str,
) -> str:
    """Render a static, escaped, non-public six-card HTML preview."""

    if language not in SUPPORTED_LANGUAGES:
        raise ValueError(f"language must be one of {SUPPORTED_LANGUAGES}.")
    summary, cards, _ = _validate_result(result)
    is_ja = language == "ja"
    title = "Exact CMLE 結果プレビュー" if is_ja else "Exact CMLE result preview"
    warning = (
        "研究用・privateプレビューです。公開UI、公開共有、理解度合格を意味しません。"
        if is_ja
        else "Research/private preview only. This is not a public UI, sharing approval, or comprehension pass."
    )
    terminal_label = "最終状態" if is_ja else "Terminal status"
    summary_headline = str(summary["HeadlineJa" if is_ja else "HeadlineEn"])
    counts = (
        f"Person: {int(summary['PersonCount'])}、exact extreme: {int(summary['ExtremePersonCount'])}、"
        f"生値／表示値不一致: {int(summary['RawDisplayMismatchCount'])}、anchor: "
        f"{int(summary['AnchoredMainEffectLevelCount'])}"
        if is_ja
        else (
            f"Persons: {int(summary['PersonCount'])}; exact extremes: "
            f"{int(summary['ExtremePersonCount'])}; raw/display mismatches: "
            f"{int(summary['RawDisplayMismatchCount'])}; anchors: "
            f"{int(summary['AnchoredMainEffectLevelCount'])}"
        )
    )
    sections = []
    for _, row in cards.iterrows():
        status = str(row["Status"])
        status_text = _STATUS_LABELS[language].get(status, status)
        headline = str(row["HeadlineJa" if is_ja else "HeadlineEn"])
        detail = str(row["DetailJa" if is_ja else "DetailEn"])
        next_action = str(row["NextActionJa" if is_ja else "NextActionEn"])
        next_label = "次の行動" if is_ja else "Next action"
        section_id = f"card-{int(row['Order'])}"
        sections.append(
            f'<section class="result-card status-{html.escape(status)}" '
            f'data-status="{html.escape(status)}" aria-labelledby="{section_id}">'
            f'<div class="status-text">{html.escape(status_text)}</div>'
            f'<h2 id="{section_id}">{int(row["Order"])}. {html.escape(headline)}</h2>'
            f'<p>{html.escape(detail)}</p>'
            f'<p><strong>{html.escape(next_label)}:</strong> {html.escape(next_action)}</p>'
            "</section>"
        )
    return f"""<!doctype html>
<html lang="{language}">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>
body{{font-family:system-ui,sans-serif;line-height:1.55;max-width:960px;margin:0 auto;padding:24px;color:#172033;background:#f6f7f9}}
.research-warning{{border:3px solid #6b4f00;background:#fff4c2;padding:16px;font-weight:700}}
.summary{{background:#fff;border:1px solid #aab2c0;padding:16px;margin:16px 0}}
.result-card{{background:#fff;border:2px solid #677184;border-left-width:10px;padding:16px;margin:14px 0}}
.status-ready,.status-available{{border-left-color:#176b45}} .status-blocked{{border-left-color:#a61b1b}}
.status-caution,.status-review{{border-left-color:#956400}} .status-not_run,.status-withheld{{border-left-color:#596273}}
.status-text{{font-weight:800;text-transform:uppercase;letter-spacing:.04em}}
h1,h2{{line-height:1.25}} code{{overflow-wrap:anywhere}}
</style>
</head>
<body data-public-surface="false" data-case-id="{html.escape(str(case_id))}">
<main aria-label="{html.escape(title)}">
<div class="research-warning" role="note">{html.escape(warning)}</div>
<h1>{html.escape(title)}</h1>
<div class="summary">
<p><strong>{html.escape(terminal_label)}:</strong> <code>{html.escape(str(summary['TerminalStatus']))}</code></p>
<p>{html.escape(summary_headline)}</p>
<p>{html.escape(counts)}</p>
</div>
{''.join(sections)}
</main>
</body>
</html>
"""


def _option_codes(options_json: str) -> set[str]:
    try:
        values = json.loads(options_json)
    except json.JSONDecodeError as exc:
        raise ValueError("Researcher key OptionsJSON is invalid.") from exc
    if not isinstance(values, list):
        raise ValueError("Researcher key OptionsJSON must be a list.")
    return {str(value["code"]) for value in values}


def score_comprehension_responses(
    researcher_key: pd.DataFrame,
    responses: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    """Score participant-case packets while retaining invalid denominators."""

    key_required = {
        "CaseId",
        "Language",
        "ItemId",
        "CorrectOption",
        "Critical",
        "MisconceptionCode",
        "OptionsJSON",
    }
    response_required = {
        "ParticipantId",
        "CaseId",
        "Language",
        "ItemId",
        "ResponseCode",
    }
    if not isinstance(researcher_key, pd.DataFrame) or not key_required.issubset(
        researcher_key.columns
    ):
        raise ValueError("researcher_key lacks required scoring fields.")
    if not isinstance(responses, pd.DataFrame) or not response_required.issubset(
        responses.columns
    ):
        raise ValueError("responses lacks required packet fields.")
    key_identity = ["CaseId", "Language", "ItemId"]
    if researcher_key.duplicated(key_identity).any():
        raise ValueError("researcher_key contains duplicate item identities.")
    work = responses.copy()
    for column in response_required:
        work[column] = work[column].astype(str)
    packet_ids = ["ParticipantId", "CaseId", "Language"]
    packet_rows = []
    item_rows = []
    for packet, observed in work.groupby(packet_ids, sort=False, dropna=False):
        participant_id, case_id, language = packet
        expected = researcher_key.loc[
            researcher_key["CaseId"].astype(str).eq(case_id)
            & researcher_key["Language"].astype(str).eq(language)
        ].copy()
        reasons = []
        if expected.empty:
            reasons.append("unknown_case_or_language")
        duplicate = observed["ItemId"].duplicated(keep=False)
        if duplicate.any():
            reasons.append("duplicate_response")
        expected_items = set(expected["ItemId"].astype(str))
        observed_items = set(observed["ItemId"].astype(str))
        if expected_items - observed_items:
            reasons.append("missing_response")
        if observed_items - expected_items:
            reasons.append("unknown_item")
        invalid_time = False
        if "ResponseTimeSeconds" in observed.columns:
            time_values = pd.to_numeric(
                observed["ResponseTimeSeconds"], errors="coerce"
            ).to_numpy(dtype=float)
            invalid_time = bool(
                np.any(~np.isfinite(time_values)) or np.any(time_values < 0)
            )
            if invalid_time:
                reasons.append("invalid_response_time")
        joined = expected.merge(
            observed.drop_duplicates("ItemId", keep="first"),
            on=["CaseId", "Language", "ItemId"],
            how="left",
            validate="one_to_one" if not duplicate.any() else "one_to_many",
        ) if not expected.empty else pd.DataFrame()
        unknown_option = False
        if not joined.empty:
            for _, row in joined.iterrows():
                response_code = str(row.get("ResponseCode", ""))
                valid_option = response_code in _option_codes(str(row["OptionsJSON"]))
                if not valid_option:
                    unknown_option = True
                correct = valid_option and response_code == str(row["CorrectOption"])
                item_rows.append(
                    {
                        "ParticipantId": participant_id,
                        "CaseId": case_id,
                        "Language": language,
                        "ItemId": row["ItemId"],
                        "ResponseCode": response_code,
                        "CorrectOption": row["CorrectOption"],
                        "Critical": bool(row["Critical"]),
                        "MisconceptionCode": row["MisconceptionCode"],
                        "ValidOption": valid_option,
                        "Correct": correct,
                    }
                )
        if unknown_option:
            reasons.append("unknown_option")
        reasons = list(dict.fromkeys(reasons))
        packet_items = [
            row
            for row in item_rows
            if row["ParticipantId"] == participant_id
            and row["CaseId"] == case_id
            and row["Language"] == language
        ]
        valid = not reasons and len(packet_items) == len(expected)
        critical_errors = sum(
            bool(row["Critical"]) and not bool(row["Correct"])
            for row in packet_items
        )
        total_correct = sum(bool(row["Correct"]) for row in packet_items)
        ready = bool(valid and critical_errors == 0)
        packet_rows.append(
            {
                "ParticipantId": participant_id,
                "CaseId": case_id,
                "Language": language,
                "ExpectedItems": len(expected),
                "ObservedRows": len(observed),
                "ValidPacket": valid,
                "InvalidReasons": ";".join(reasons),
                "TotalCorrect": total_correct,
                "TotalItems": len(expected),
                "CriticalErrors": critical_errors,
                "ParticipantCaseReady": ready,
                "Disposition": (
                    "ready"
                    if ready
                    else "invalid"
                    if not valid
                    else "critical_error"
                ),
            }
        )
    item_audit = pd.DataFrame(item_rows)
    participant_summary = pd.DataFrame(packet_rows)
    study_summary = pd.DataFrame(
        [
            {
                "HumanStudyStatus": HUMAN_STUDY_STATUS,
                "AttemptedPackets": len(participant_summary),
                "ValidPackets": int(
                    participant_summary.get("ValidPacket", pd.Series(dtype=bool)).sum()
                ),
                "InvalidPackets": int(
                    (~participant_summary.get("ValidPacket", pd.Series(dtype=bool))).sum()
                )
                if not participant_summary.empty
                else 0,
                "ReadyPackets": int(
                    participant_summary.get(
                        "ParticipantCaseReady", pd.Series(dtype=bool)
                    ).sum()
                ),
                "CriticalErrorPackets": int(
                    participant_summary.get(
                        "CriticalErrors", pd.Series(dtype=int)
                    ).gt(0).sum()
                ),
                "PublicSurfaceEnabled": False,
                "Interpretation": "synthetic_scoring_logic_only_no_human_outcomes",
            }
        ]
    )
    return {
        "item_audit": item_audit,
        "participant_summary": participant_summary,
        "study_summary": study_summary,
    }


def comprehension_human_gate_status() -> pd.DataFrame:
    """Return the non-promotable human-study state before recruitment."""

    return pd.DataFrame(
        [
            {
                "SchemaVersion": CMLE_COMPREHENSION_SCHEMA_VERSION,
                "HumanStudyStatus": HUMAN_STUDY_STATUS,
                "HumanParticipants": 0,
                "ComprehensionRateAvailable": False,
                "LanguageEquivalenceAvailable": False,
                "PublicSurfaceEnabled": False,
                "NextGate": "cognitive_interviews_then_frozen_human_pilot",
            }
        ]
    )


__all__ = [
    "CMLE_COMPREHENSION_SCHEMA_VERSION",
    "HUMAN_STUDY_STATUS",
    "ITEM_IDS",
    "SUPPORTED_LANGUAGES",
    "build_comprehension_materials",
    "comprehension_human_gate_status",
    "render_cmle_one_click_preview_html",
    "score_comprehension_responses",
]
