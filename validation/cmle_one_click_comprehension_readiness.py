#!/usr/bin/env python3
"""Validate non-public bilingual CMLE preview/comprehension readiness.

This runner uses analytical fixtures and synthetic response packets only.  It
cannot produce human comprehension, accessibility, or language-equivalence
evidence.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle_one_click import (  # noqa: E402
    CMLE_ONE_CLICK_CARDS,
    run_cmle_one_click_analysis,
)
from mfrm_app.cmle_one_click_comprehension import (  # noqa: E402
    HUMAN_STUDY_STATUS,
    ITEM_IDS,
    SUPPORTED_LANGUAGES,
    build_comprehension_materials,
    comprehension_human_gate_status,
    render_cmle_one_click_preview_html,
    score_comprehension_responses,
)


PLAN = ROOT / "validation/cmle_one_click_comprehension_readiness_plan_20260810.json"
OUTPUT = ROOT / "validation/cmle_one_click_comprehension_readiness_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, float_format="%.17g")


def validate_plan() -> dict[str, object]:
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    if plan.get("study_id") != "cmle_one_click_comprehension_readiness_v1":
        raise ValueError("Unexpected comprehension-readiness plan.")
    mismatches = [
        relative
        for relative, expected in plan["parent_identity"].items()
        if not (ROOT / relative).is_file()
        or sha256_file(ROOT / relative) != expected
    ]
    if mismatches:
        raise ValueError(f"Comprehension parent identity failed: {mismatches}")
    return plan


def interior_frame(*, extremes: bool) -> pd.DataFrame:
    scores = {
        "P1": [0, 1, 1, 2],
        "P2": [1, 0, 2, 1],
        "P3": [2, 1, 1, 0],
        "P4": [1, 2, 0, 1],
        "P5": [0, 2, 1, 1],
        "P6": [2, 0, 1, 1],
    }
    if extremes:
        scores.update({"P7": [0, 0, 0, 0], "P8": [2, 2, 2, 2]})
    units = [("R1", "C1"), ("R1", "C2"), ("R2", "C1"), ("R2", "C2")]
    return pd.DataFrame(
        [
            (person, rater, criterion, score)
            for person, values in scores.items()
            for (rater, criterion), score in zip(units, values, strict=True)
        ],
        columns=["Person", "Rater", "Criterion", "Score"],
    )


def boundary_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            (f"P{index:03d}", rater, score)
            for index in range(20)
            for rater, score in (("R1", 1), ("R2", 0))
        ],
        columns=["Person", "Rater", "Score"],
    )


def structural_frame() -> pd.DataFrame:
    patterns = ([0, 1, 2], [1, 2, 0], [2, 0, 1], [0, 2, 1])
    rows = []
    for person_index in range(12):
        rater = "R1" if person_index < 6 else "R2"
        rows.extend(
            (f"P{person_index}", rater, criterion, score)
            for criterion, score in zip(
                ["C1", "C2", "C3"], patterns[person_index % 4], strict=True
            )
        )
    return pd.DataFrame(rows, columns=["Person", "Rater", "Criterion", "Score"])


def case_inputs():
    anchors = pd.DataFrame(
        [
            {"ParameterType": "Facet", "Facet": "Rater", "Level": "R1", "Value": 0.25},
            {"ParameterType": "Facet", "Facet": "Criterion", "Level": "C1", "Value": -0.2},
        ]
    )
    common = {
        "person_col": "Person",
        "score_col": "Score",
        "rating_min": 0,
        "gtol": 1e-8,
        "maxiter": 800,
        "display_decimals": 3,
    }
    return [
        (
            "rsm_interior_fit_sample_extremes",
            interior_frame(extremes=True),
            {**common, "facet_cols": ["Rater", "Criterion"], "rating_max": 2, "model": "RSM"},
        ),
        (
            "pcm_interior_differential_hard_anchors",
            interior_frame(extremes=False),
            {
                **common,
                "facet_cols": ["Rater", "Criterion"],
                "rating_max": 2,
                "model": "PCM",
                "step_facet": "Criterion",
                "hard_anchors": anchors,
            },
        ),
        (
            "binary_conditional_support_boundary",
            boundary_frame(),
            {**common, "facet_cols": ["Rater"], "rating_max": 1, "model": "RSM"},
        ),
        (
            "rsm_structurally_disconnected",
            structural_frame(),
            {**common, "facet_cols": ["Rater", "Criterion"], "rating_max": 2, "model": "RSM"},
        ),
        (
            "invalid_missing_score_column",
            interior_frame(extremes=False).drop(columns=["Score"]),
            {**common, "facet_cols": ["Rater", "Criterion"], "rating_max": 2, "model": "RSM"},
        ),
    ]


class PreviewParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.sections = 0
        self.status_texts = 0
        self.h2_ids: list[str] = []
        self.aria_labels: list[str] = []
        self.scripts = 0
        self.external_links = 0
        self.lang = ""
        self.public_surface = ""
        self.warning_notes = 0

    def handle_starttag(self, tag: str, attrs) -> None:
        values = dict(attrs)
        classes = str(values.get("class", "")).split()
        if tag == "html":
            self.lang = str(values.get("lang", ""))
        if tag == "body":
            self.public_surface = str(values.get("data-public-surface", ""))
        if tag == "section" and "result-card" in classes:
            self.sections += 1
            self.aria_labels.append(str(values.get("aria-labelledby", "")))
        if tag == "div" and "status-text" in classes:
            self.status_texts += 1
        if tag == "div" and values.get("role") == "note":
            self.warning_notes += 1
        if tag == "h2":
            self.h2_ids.append(str(values.get("id", "")))
        if tag == "script":
            self.scripts += 1
        for attribute in ("href", "src"):
            target = str(values.get(attribute, ""))
            if target.startswith(("http://", "https://", "//")):
                self.external_links += 1


def audit_preview(html_text: str, *, case_id: str, language: str) -> dict[str, object]:
    parser = PreviewParser()
    parser.feed(html_text)
    return {
        "CaseId": case_id,
        "Language": language,
        "SectionCount": parser.sections,
        "HeadingCount": len(parser.h2_ids),
        "StatusTextCount": parser.status_texts,
        "WarningNoteCount": parser.warning_notes,
        "LanguageMatch": parser.lang == language,
        "AriaHeadingLinksMatch": parser.aria_labels == parser.h2_ids,
        "UniqueHeadingIds": len(set(parser.h2_ids)) == len(parser.h2_ids),
        "ScriptCount": parser.scripts,
        "ExternalResourceCount": parser.external_links,
        "PublicSurfaceFalse": parser.public_surface == "false",
        "AuditPassed": bool(
            parser.sections == 6
            and len(parser.h2_ids) == 6
            and parser.status_texts == 6
            and parser.warning_notes == 1
            and parser.lang == language
            and parser.aria_labels == parser.h2_ids
            and len(set(parser.h2_ids)) == 6
            and parser.scripts == 0
            and parser.external_links == 0
            and parser.public_surface == "false"
        ),
    }


def run_cases(output: Path):
    expected_cases = {
        "rsm_interior_fit_sample_extremes": {
            "TerminalStatus": "analysis_ready_with_cautions",
            "Optimization": "attempted",
            "Calibration": "ready",
            "Person": "performed",
            "Anchor": "not_applicable_no_anchor",
            "Next": "review_cautions",
            "Estimator": "fixed_calibration_wle",
        },
        "pcm_interior_differential_hard_anchors": {
            "TerminalStatus": "analysis_ready_with_cautions",
            "Optimization": "attempted",
            "Calibration": "ready",
            "Person": "performed",
            "Anchor": "does_not_validate",
            "Next": "review_cautions",
            "Estimator": "fixed_calibration_wle",
        },
        "binary_conditional_support_boundary": {
            "TerminalStatus": "finite_mle_boundary",
            "Optimization": "not_attempted",
            "Calibration": "not_ready",
            "Person": "not_performed",
            "Anchor": "not_applicable_no_anchor",
            "Next": "review_boundary",
            "Estimator": "not_computed",
        },
        "rsm_structurally_disconnected": {
            "TerminalStatus": "design_not_identified",
            "Optimization": "not_attempted",
            "Calibration": "not_ready",
            "Person": "not_performed",
            "Anchor": "not_applicable_no_anchor",
            "Next": "repair_design",
            "Estimator": "not_computed",
        },
        "invalid_missing_score_column": {
            "TerminalStatus": "input_invalid",
            "Optimization": "not_attempted",
            "Calibration": "not_ready",
            "Person": "not_performed",
            "Anchor": "not_applicable_no_anchor",
            "Next": "correct_input",
            "Estimator": "not_computed",
        },
    }
    previews = output / "previews"
    previews.mkdir()
    ledger_rows = []
    audit_rows = []
    participant_frames = []
    key_frames = []
    retained_results = {}
    for case_id, frame, kwargs in case_inputs():
        result = run_cmle_one_click_analysis(frame, **kwargs)
        retained_results[case_id] = result
        summary = result["summary"].iloc[0]
        expected = expected_cases[case_id]
        ledger_rows.append(
            {
                "CaseId": case_id,
                "ExpectedTerminalStatus": expected["TerminalStatus"],
                "ObservedTerminalStatus": summary["TerminalStatus"],
                "TerminalStatusMatch": summary["TerminalStatus"] == expected["TerminalStatus"],
                "CalibrationReady": summary["CalibrationReady"],
                "PersonScoringReady": summary["PersonScoringReady"],
                "CardCount": len(result["cards"]),
                "CardOrderMatch": result["cards"]["Card"].tolist() == list(CMLE_ONE_CLICK_CARDS),
                "AnchoredMainEffectLevelCount": summary["AnchoredMainEffectLevelCount"],
                "PublicSurfaceEnabled": summary["PublicSurfaceEnabled"],
            }
        )
        for language in SUPPORTED_LANGUAGES:
            html_text = render_cmle_one_click_preview_html(
                result, language=language, case_id=case_id
            )
            (previews / f"{case_id}_{language}.html").write_text(
                html_text, encoding="utf-8"
            )
            audit_rows.append(audit_preview(html_text, case_id=case_id, language=language))
            materials = build_comprehension_materials(
                result, case_id=case_id, language=language
            )
            participant_frames.append(materials["participant_tasks"])
            key_frames.append(materials["researcher_key"])
    ledger = pd.DataFrame(ledger_rows)
    audit = pd.DataFrame(audit_rows)
    participant = pd.concat(participant_frames, ignore_index=True)
    key = pd.concat(key_frames, ignore_index=True)

    expected_answer_rows = []
    for case_id, expected in expected_cases.items():
        expected_answers = {
            "calibration_readiness": expected["Calibration"],
            "optimization_attempted": expected["Optimization"],
            "person_scoring_performed": expected["Person"],
            "fit_decision_input": "raw_unrounded",
            "rounding_vignette": "raw_noisy",
            "calibration_uncertainty": "not_propagated",
            "anchor_validity": expected["Anchor"],
            "public_sharing": "not_permitted",
            "next_action": expected["Next"],
            "person_estimator_role": expected["Estimator"],
        }
        for language in SUPPORTED_LANGUAGES:
            expected_answer_rows.extend(
                {
                    "CaseId": case_id,
                    "Language": language,
                    "ItemId": item_id,
                    "RegisteredCorrectOption": answer,
                }
                for item_id, answer in expected_answers.items()
            )
    expected_answers = pd.DataFrame(expected_answer_rows)
    answer_audit = key.merge(
        expected_answers,
        on=["CaseId", "Language", "ItemId"],
        validate="one_to_one",
    )
    answer_audit["AnswerMatch"] = answer_audit["CorrectOption"].eq(
        answer_audit["RegisteredCorrectOption"]
    )

    parity_rows = []
    for (case_id, item_id), group in key.groupby(["CaseId", "ItemId"], sort=False):
        by_language = group.set_index("Language")
        option_codes = {
            language: [row["code"] for row in json.loads(by_language.loc[language, "OptionsJSON"])]
            for language in SUPPORTED_LANGUAGES
        }
        parity_rows.append(
            {
                "CaseId": case_id,
                "ItemId": item_id,
                "LanguageRows": len(group),
                "OptionCodeParity": option_codes["en"] == option_codes["ja"],
                "CorrectOptionParity": by_language.loc["en", "CorrectOption"] == by_language.loc["ja", "CorrectOption"],
                "CriticalityParity": bool(by_language.loc["en", "Critical"]) == bool(by_language.loc["ja", "Critical"]),
                "MisconceptionParity": by_language.loc["en", "MisconceptionCode"] == by_language.loc["ja", "MisconceptionCode"],
            }
        )
    parity = pd.DataFrame(parity_rows)
    parity["ParityPassed"] = parity[
        ["OptionCodeParity", "CorrectOptionParity", "CriticalityParity", "MisconceptionParity"]
    ].all(axis=1) & parity["LanguageRows"].eq(2)

    injected = copy.deepcopy(retained_results["rsm_interior_fit_sample_extremes"])
    injected["cards"].loc[0, "HeadlineEn"] = '<script src="https://invalid.example/x.js">alert(1)</script>'
    escaped = render_cmle_one_click_preview_html(
        injected, language="en", case_id="escape_injection"
    )
    escape_audit = pd.DataFrame(
        [
            {
                "RawScriptAbsent": "<script" not in escaped.lower(),
                "ExternalURLIsEscapedText": "https://invalid.example" in escaped,
                "EscapedScriptPresent": "&lt;script" in escaped,
                "EscapePassed": "<script" not in escaped.lower() and "&lt;script" in escaped,
            }
        ]
    )

    write_csv(ledger, output / "scenario_ledger.csv")
    write_csv(audit, output / "html_preview_audit.csv")
    write_csv(participant, output / "participant_task_packet.csv")
    write_csv(key, output / "researcher_scoring_key.csv")
    write_csv(answer_audit, output / "registered_answer_audit.csv")
    write_csv(parity, output / "language_structure_parity.csv")
    write_csv(escape_audit, output / "html_escape_audit.csv")
    return ledger, audit, participant, key, answer_audit, parity, escape_audit


def response_rows(key: pd.DataFrame, participant_id: str, mode: str) -> pd.DataFrame:
    rows = []
    for _, item in key.iterrows():
        codes = [entry["code"] for entry in json.loads(item["OptionsJSON"])]
        response = str(item["CorrectOption"])
        if mode == "all_wrong":
            response = next(code for code in codes if code != response)
        elif mode == "false_green" and item["ItemId"] == "calibration_readiness":
            response = "ready"
        rows.append(
            {
                "ParticipantId": participant_id,
                "CaseId": item["CaseId"],
                "Language": item["Language"],
                "ItemId": item["ItemId"],
                "ResponseCode": response,
                "ResponseTimeSeconds": 3.0,
                "SyntheticMode": mode,
            }
        )
    return pd.DataFrame(rows)


def synthetic_scoring(key: pd.DataFrame, output: Path):
    packets = []
    blocked_cases = {
        "binary_conditional_support_boundary",
        "rsm_structurally_disconnected",
        "invalid_missing_score_column",
    }
    for (case_id, language), case_key in key.groupby(["CaseId", "Language"], sort=False):
        packets.append(response_rows(case_key, f"perfect::{case_id}::{language}", "perfect"))
        packets.append(response_rows(case_key, f"wrong::{case_id}::{language}", "all_wrong"))
        if case_id in blocked_cases:
            packets.append(response_rows(case_key, f"false_green::{case_id}::{language}", "false_green"))

    base_key = key.loc[
        key["CaseId"].eq("rsm_interior_fit_sample_extremes")
        & key["Language"].eq("en")
    ]
    base = response_rows(base_key, "invalid_base", "perfect")
    duplicate = pd.concat(
        [
            base.assign(ParticipantId="invalid::duplicate", SyntheticMode="duplicate"),
            base.iloc[[0]].assign(ParticipantId="invalid::duplicate", SyntheticMode="duplicate"),
        ],
        ignore_index=True,
    )
    missing = base.iloc[:-1].assign(ParticipantId="invalid::missing", SyntheticMode="missing")
    unknown = base.assign(ParticipantId="invalid::unknown_option", SyntheticMode="unknown_option")
    unknown.loc[unknown.index[0], "ResponseCode"] = "not_a_registered_option"
    bad_time = base.assign(ParticipantId="invalid::time", SyntheticMode="invalid_time")
    bad_time.loc[bad_time.index[0], "ResponseTimeSeconds"] = -1.0
    packets.extend([duplicate, missing, unknown, bad_time])
    responses = pd.concat(packets, ignore_index=True)
    scored = score_comprehension_responses(key, responses)
    write_csv(responses, output / "synthetic_responses.csv")
    write_csv(scored["item_audit"], output / "synthetic_item_audit.csv")
    write_csv(scored["participant_summary"], output / "synthetic_packet_summary.csv")
    write_csv(scored["study_summary"], output / "synthetic_scoring_summary.csv")
    return responses, scored


def run_tests(output: Path) -> dict[str, object]:
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-q",
        "tests/test_cmle_one_click_comprehension.py",
        "tests/test_cmle_one_click.py",
        "tests/test_cmle_one_click_archive.py",
        "tests/test_decision_stability.py",
        "tests/test_threshold_decision_integration.py",
    ]
    completed = subprocess.run(
        command, cwd=ROOT, text=True, capture_output=True, check=False
    )
    (output / "selected_tests_stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output / "selected_tests_stderr.txt").write_text(completed.stderr, encoding="utf-8")
    return {
        "passed": completed.returncode == 0,
        "returncode": completed.returncode,
        "command": command,
    }


def write_protocol(output: Path) -> None:
    text = """# Prospective human comprehension protocol (not yet started)

## Current state

No human responses have been collected. The HTML files are research/private instrument previews, not a public application surface. Synthetic packets validate scoring code only.

## Phase 1: cognitive interviews

Recruit 6–10 intended users for each language form across novice and experienced strata. Use think-aloud plus retrospective probing. Record wording failures and revise the instrument; this phase cannot promote the interface.

## Phase 2: pilot after instrument freeze

After freezing revisions, collect at least 24 complete participants per language. Retain attempted, incomplete, invalid, and valid denominators. Report every item error, critical misconception, completion time, and confidence. Pilot results alone cannot promote the interface.

## Fresh confirmatory gate

Register the final instrument again before confirmatory data. For every critical item within each language, require the one-sided 95% Wilson upper bound for the dangerous-error proportion to be below 0.10. Total scores and pooled languages cannot compensate for any failing item-language cell.

## Ethics, privacy, and language claims

Use synthetic analytical scenarios only; collect no operational rating records. Minimize participant identifiers and establish consent, access, and retention controls before recruitment. Do not claim Japanese–English equivalence from separate convenience samples; that requires a separately justified randomized bilingual design.
"""
    (output / "PROSPECTIVE_HUMAN_STUDY_PROTOCOL.md").write_text(text, encoding="utf-8")


def write_review(output: Path, results: dict[str, object]) -> None:
    text = f"""# Bilingual CMLE comprehension-readiness critical review

## Decision

The automated instrument-readiness contract **{'passed' if results['contract_passed'] else 'failed'}**. Human comprehension remains **not assessed**: zero participants, no comprehension rate, no language-equivalence estimate, and no public surface.

## What the machine evidence establishes

The retained five analytical states generated 10 escaped, standalone English/Japanese previews with six ordered cards each. The participant packet and separately stored key each contain 100 rows. English and Japanese share stable item and option codes, correct answers, criticality, and misconception domains. The 1.5004-to-1.500 vignette is scored from raw MnSq in every form.

Thirty synthetic participant-case packets exercised the scorer: 10 perfect packets passed; 10 all-wrong and 6 single false-green packets failed; duplicate, missing, unknown-option, and invalid-time packets were invalidated rather than silently omitted.

## What it does not establish

- Structural bilingual parity is not translation equivalence or cultural validity.
- Static HTML checks are not assistive-technology, keyboard, visual, or cognitive accessibility evidence.
- Synthetic packets are not observations of user comprehension or task completion.
- Six-card consistency does not show that terminology such as CMLE, WLE, MnSq, hard anchor, or conditional SE is understood.
- The previews and private archives are neither de-identified nor approved for public sharing.

## Promotion rule

Keep the Streamlit/public UI withheld. First conduct bilingual cognitive interviews, revise and freeze the instrument, then pilot it. A later confirmatory gate must pass each of six dangerous-misconception items separately within each language; no average score may conceal false-green calibration, invented Person scores, rounded-MnSq reclassification, ignored calibration uncertainty, anchor-validity overclaim, or inappropriate public sharing.
"""
    (output / "HTML_PREVIEW_INSTRUMENT_CRITICAL_REVIEW.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    plan = validate_plan()
    if args.output.exists():
        raise FileExistsError(f"Evidence output exists: {args.output}")
    args.output.mkdir(parents=True)

    ledger, html_audit, participant, key, answer_audit, parity, escape_audit = run_cases(args.output)
    responses, scored = synthetic_scoring(key, args.output)
    human_gate = comprehension_human_gate_status()
    write_csv(human_gate, args.output / "human_gate_status.csv")
    write_protocol(args.output)
    tests = run_tests(args.output)

    packet_summary = scored["participant_summary"]
    participant_forbidden = {"CorrectOption", "Critical", "MisconceptionCode", "Rationale"}
    critical_by_form = key.groupby(["CaseId", "Language"])["Critical"].sum()
    rounding_rows = key.loc[key["ItemId"].eq("rounding_vignette")]
    gates = {
        "identity_passed": True,
        "scenario_state_passed": bool(
            len(ledger) == 5
            and ledger["TerminalStatusMatch"].all()
            and ledger["CardOrderMatch"].all()
            and ledger["CardCount"].eq(6).all()
            and not ledger["PublicSurfaceEnabled"].any()
        ),
        "html_preview_passed": bool(
            len(html_audit) == 10
            and html_audit["AuditPassed"].all()
            and bool(escape_audit.iloc[0]["EscapePassed"])
        ),
        "participant_key_separation_passed": bool(
            len(participant) == 100
            and len(key) == 100
            and not participant_forbidden.intersection(participant.columns)
            and not participant.duplicated(["CaseId", "Language", "ItemId"]).any()
            and not key.duplicated(["CaseId", "Language", "ItemId"]).any()
        ),
        "language_parity_passed": bool(len(parity) == 50 and parity["ParityPassed"].all()),
        "registered_answers_passed": bool(
            len(answer_audit) == 100 and answer_audit["AnswerMatch"].all()
        ),
        "critical_domain_passed": bool(
            critical_by_form.eq(6).all()
            and key.loc[key["Critical"], "MisconceptionCode"].nunique() == 6
        ),
        "rounding_contract_passed": bool(
            len(rounding_rows) == 10
            and rounding_rows["CorrectOption"].eq("raw_noisy").all()
            and rounding_rows["Critical"].all()
            and rounding_rows["MisconceptionCode"].eq("rounded_mnsq_reclassification").all()
        ),
        "synthetic_scoring_passed": bool(
            len(packet_summary) == 30
            and packet_summary["Disposition"].value_counts().to_dict()
            == {"critical_error": 16, "ready": 10, "invalid": 4}
            and packet_summary.loc[
                packet_summary["ParticipantId"].str.startswith("perfect::"),
                "ParticipantCaseReady",
            ].all()
            and not packet_summary.loc[
                ~packet_summary["ParticipantId"].str.startswith("perfect::"),
                "ParticipantCaseReady",
            ].any()
            and int(packet_summary.loc[
                packet_summary["ParticipantId"].str.startswith("wrong::"),
                "CriticalErrors",
            ].min()) == 6
        ),
        "human_and_public_withheld_passed": bool(
            len(human_gate) == 1
            and human_gate.iloc[0]["HumanStudyStatus"] == HUMAN_STUDY_STATUS
            and int(human_gate.iloc[0]["HumanParticipants"]) == 0
            and not bool(human_gate.iloc[0]["ComprehensionRateAvailable"])
            and not bool(human_gate.iloc[0]["LanguageEquivalenceAvailable"])
            and not bool(human_gate.iloc[0]["PublicSurfaceEnabled"])
        ),
        "tests_passed": bool(tests["passed"]),
    }
    contract_passed = all(gates.values())
    results = {
        **gates,
        "contract_passed": contract_passed,
        "preview_count": len(html_audit),
        "participant_task_rows": len(participant),
        "researcher_key_rows": len(key),
        "synthetic_packet_count": len(packet_summary),
        "synthetic_dispositions": packet_summary["Disposition"].value_counts().to_dict(),
        "human_participants": 0,
        "human_comprehension_passed": False,
        "human_comprehension_status": "not_assessed",
        "language_equivalence_status": "not_assessed",
        "public_surface_enabled": False,
        "selected_tests": tests,
    }
    write_review(args.output, results)
    output_files = sorted(
        path.relative_to(args.output).as_posix()
        for path in args.output.rglob("*")
        if path.is_file() and path.name != "decision.json"
    )
    decision = json_safe(
        {
            "study_id": plan["study_id"],
            "plan_sha256": sha256_file(PLAN),
            "contract_passed": contract_passed,
            "contract_interpretation": "automated_instrument_readiness_only",
            "implementation_sha256": {
                "mfrm_app/cmle_one_click_comprehension.py": sha256_file(
                    ROOT / "mfrm_app/cmle_one_click_comprehension.py"
                ),
                "tests/test_cmle_one_click_comprehension.py": sha256_file(
                    ROOT / "tests/test_cmle_one_click_comprehension.py"
                ),
                "validation/cmle_one_click_comprehension_readiness.py": sha256_file(Path(__file__)),
            },
            "results": results,
            "output_sha256": {
                name: sha256_file(args.output / name) for name in output_files
            },
            "interpretation": {
                "human_comprehension": "not assessed; zero participants",
                "language_equivalence": "not assessed",
                "accessibility": "static structure only; assistive-technology testing not performed",
                "decision_input": "finite unrounded MnSq",
                "public_ui": "withheld",
            },
        }
    )
    (args.output / "decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(decision, indent=2, ensure_ascii=False, allow_nan=False))
    if not contract_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
