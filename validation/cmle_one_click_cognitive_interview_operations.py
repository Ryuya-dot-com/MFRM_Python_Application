#!/usr/bin/env python3
"""Build and validate a no-human-data CMLE cognitive-interview kit."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle_one_click_cognitive_interview import (  # noqa: E402
    EXPERIENCE_STRATA,
    ISSUE_CODES,
    RESEARCHER_ONLY_FIELDS,
    build_blank_cognitive_record_template,
    build_cognitive_interview_schedule,
    build_moderator_probe_guide,
    build_participant_session_tasks,
    cognitive_interview_human_gate_status,
    validate_cognitive_interview_records,
)


PLAN = ROOT / "validation/cmle_one_click_cognitive_interview_operations_plan_20260810.json"
SOURCE = ROOT / "validation/cmle_one_click_comprehension_readiness_20260810"
OUTPUT = ROOT / "validation/cmle_one_click_cognitive_interview_operations_20260810"


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
    if plan.get("study_id") != "cmle_one_click_cognitive_interview_operations_v1":
        raise ValueError("Unexpected cognitive-interview operations plan.")
    mismatches = [
        relative
        for relative, expected in plan["parent_identity"].items()
        if not (ROOT / relative).is_file()
        or sha256_file(ROOT / relative) != expected
    ]
    if mismatches:
        raise ValueError(f"Cognitive-interview parent identity failed: {mismatches}")
    return plan


def load_source_materials() -> tuple[pd.DataFrame, pd.DataFrame]:
    participant = pd.read_csv(SOURCE / "participant_task_packet.csv")
    key = pd.read_csv(SOURCE / "researcher_scoring_key.csv")
    return participant, key


def schedule_audits(schedule: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    case_rows = []
    pair_rows = []
    for language, language_frame in schedule.groupby("Language", sort=False):
        for case_id, case_frame in language_frame.groupby("CaseId", sort=False):
            case_rows.append(
                {
                    "Language": language,
                    "CaseId": case_id,
                    "Exposures": len(case_frame),
                    "FirstPosition": int(case_frame["CaseOrder"].eq(1).sum()),
                    "SecondPosition": int(case_frame["CaseOrder"].eq(2).sum()),
                    "NoviceTarget": int(case_frame["ExperienceStratum"].eq("novice_target").sum()),
                    "ExperiencedTarget": int(case_frame["ExperienceStratum"].eq("experienced_target").sum()),
                }
            )
        for slot_id, slot in language_frame.groupby("SessionSlotId", sort=False):
            ordered = slot.sort_values("CaseOrder")
            pair_rows.append(
                {
                    "Language": language,
                    "SessionSlotId": slot_id,
                    "PairId": ordered.iloc[0]["PairId"],
                    "ExperienceStratum": ordered.iloc[0]["ExperienceStratum"],
                    "FirstCase": ordered.iloc[0]["CaseId"],
                    "SecondCase": ordered.iloc[1]["CaseId"],
                    "UnorderedPair": "|".join(sorted(ordered["CaseId"].tolist())),
                    "AssignmentRows": len(ordered),
                }
            )
    return pd.DataFrame(case_rows), pd.DataFrame(pair_rows)


def write_session_packets(tasks: pd.DataFrame, output: Path) -> pd.DataFrame:
    packet_dir = output / "participant_session_packets"
    packet_dir.mkdir()
    manifest_rows = []
    for slot_id, packet in tasks.groupby("SessionSlotId", sort=False):
        path = packet_dir / f"{slot_id}.csv"
        write_csv(packet, path)
        manifest_rows.append(
            {
                "SessionSlotId": slot_id,
                "Language": packet.iloc[0]["Language"],
                "ExperienceStratum": packet.iloc[0]["ExperienceStratum"],
                "TaskRows": len(packet),
                "CaseCount": packet["CaseId"].nunique(),
                "ResearcherFieldCount": len(RESEARCHER_ONLY_FIELDS.intersection(packet.columns)),
                "SHA256": sha256_file(path),
                "RelativePath": path.relative_to(output).as_posix(),
            }
        )
    return pd.DataFrame(manifest_rows)


def complete_synthetic_record(template: pd.DataFrame, key: pd.DataFrame) -> pd.DataFrame:
    result = template.copy()
    option_map = {
        (str(row["CaseId"]), str(row["Language"]), str(row["ItemId"])): json.loads(
            str(row["OptionsJSON"])
        )[0]["code"]
        for _, row in key.iterrows()
    }
    result["ResponseCode"] = result.apply(
        lambda row: option_map[(str(row["CaseId"]), str(row["Language"]), str(row["ItemId"]))],
        axis=1,
    )
    result["ResponseLocked"] = True
    result["ProbeStartedAfterLock"] = True
    result["ResponseTimeSeconds"] = 10.25
    result["Confidence"] = 3
    result["IssueCodesJSON"] = "[]"
    result["IssueSeverity"] = "none"
    result["ParaphrasedObservation"] = ""
    result["PIIReviewed"] = False
    result["Withdrawn"] = False
    result["CompletionState"] = "completed"
    return result


def run_synthetic_validator(
    participant: pd.DataFrame,
    key: pd.DataFrame,
    schedule: pd.DataFrame,
    template: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    complete = complete_synthetic_record(template, key)
    cases: list[tuple[str, pd.DataFrame, bool, set[str], bool]] = [
        ("empty_template", template, True, set(), True),
        ("complete_valid_synthetic", complete, False, set(), True),
    ]
    pii = complete.assign(Email="synthetic@example.invalid")
    cases.append(("unexpected_pii_column", pii, False, {"forbidden_pii_column", "unexpected_column"}, False))
    assignment = complete.copy()
    assignment.loc[0, "CaseId"] = "unregistered_case"
    cases.append(("unregistered_case_assignment", assignment, False, {"assignment_mismatch"}, False))
    session = complete.copy()
    session.loc[0, "SessionSlotId"] = "CI-UNREGISTERED-99"
    cases.append(("unregistered_session", session, False, {"missing_item", "unregistered_assignment"}, False))
    duplicate = pd.concat([complete, complete.iloc[[0]]], ignore_index=True)
    cases.append(("duplicate_item", duplicate, False, {"record_row_count", "duplicate_item"}, False))
    missing = complete.iloc[:-1].copy()
    cases.append(("missing_item", missing, False, {"record_row_count", "missing_item"}, False))
    response = complete.copy()
    response.loc[0, "ResponseCode"] = "not_a_registered_option"
    cases.append(("unknown_response_option", response, False, {"unknown_response_option"}, False))
    issue = complete.copy()
    issue.loc[0, "IssueCodesJSON"] = '["not_a_registered_issue"]'
    issue.loc[0, "IssueSeverity"] = "major"
    cases.append(("unknown_issue_code", issue, False, {"unknown_issue_code"}, False))
    severity = complete.copy()
    severity.loc[0, "IssueCodesJSON"] = f'["{ISSUE_CODES[0]}"]'
    severity.loc[0, "IssueSeverity"] = "none"
    cases.append(("issue_severity_mismatch", severity, False, {"issue_severity_mismatch"}, False))
    note = complete.copy()
    note.loc[0, "ParaphrasedObservation"] = "Synthetic paraphrased observation."
    note.loc[0, "PIIReviewed"] = False
    cases.append(("unreviewed_nonempty_note", note, False, {"unreviewed_note"}, False))
    probe = complete.copy()
    probe.loc[0, "ProbeStartedAfterLock"] = False
    cases.append(("probe_before_response", probe, False, {"probe_before_response_lock"}, False))
    numeric = complete.copy()
    numeric.loc[0, "ResponseTimeSeconds"] = -1
    numeric.loc[1, "Confidence"] = 6
    cases.append(("invalid_time_and_confidence", numeric, False, {"invalid_response_time", "invalid_confidence"}, False))

    ledger_rows = []
    failure_frames = []
    for case_id, records, template_only, expected_codes, expected_valid in cases:
        result = validate_cognitive_interview_records(
            records, key, schedule, template_only=template_only
        )
        observed_codes = set(result["failure_codes"])
        contract_match = bool(result["valid"] == expected_valid and expected_codes.issubset(observed_codes))
        ledger_rows.append(
            {
                "SyntheticCaseId": case_id,
                "ExpectedValid": expected_valid,
                "ObservedValid": result["valid"],
                "ExpectedCodes": ";".join(sorted(expected_codes)),
                "ObservedCodes": ";".join(result["failure_codes"]),
                "TemplateOnly": template_only,
                "ContractMatch": contract_match,
                "Interpretation": "validator_logic_only_not_human_outcome",
            }
        )
        if not result["failures"].empty:
            failures = result["failures"].copy()
            failures.insert(0, "SyntheticCaseId", case_id)
            failure_frames.append(failures)

    try:
        build_participant_session_tasks(
            participant.assign(CorrectOption="forbidden_answer_leakage"), schedule
        )
        leakage_rejected = False
    except ValueError:
        leakage_rejected = True
    ledger_rows.append(
        {
            "SyntheticCaseId": "answer_key_field_in_participant_packet",
            "ExpectedValid": False,
            "ObservedValid": not leakage_rejected,
            "ExpectedCodes": "answer_key_leakage",
            "ObservedCodes": "answer_key_leakage" if leakage_rejected else "",
            "TemplateOnly": False,
            "ContractMatch": leakage_rejected,
            "Interpretation": "validator_logic_only_not_human_outcome",
        }
    )
    failure_audit = (
        pd.concat(failure_frames, ignore_index=True)
        if failure_frames
        else pd.DataFrame(columns=["SyntheticCaseId", "Code", "Location", "Detail"])
    )
    return pd.DataFrame(ledger_rows), failure_audit


def write_operations_documents(output: Path) -> None:
    private_readme = """# Private cognitive-interview operations kit

This folder contains templates only. It does not authorize recruitment, satisfy institutional ethics review, provide final consent language, or contain human responses.

## Before any session

1. Obtain all applicable ethics, privacy, consent, accessibility, and retention approvals outside this repository.
2. Keep names, email addresses, contact rosters, consent records, and recordings outside the analytical kit with separately governed access.
3. Assign a private `SessionSlotId`; do not replace it with a person identifier.
4. Confirm the frozen `InstrumentVersion`. Any wording or layout edit requires a new version and must not be pooled silently.

## During a session

Present the assigned preview and question, collect and lock the unaided response, and only then use the moderator probe. Do not teach CMLE, WLE, MnSq, anchor validity, uncertainty scope, or public-sharing status before response lock. Stop on withdrawal or burden.

## After a session

Use controlled issue codes and paraphrased observations. A non-empty observation requires documented PII review. Preserve attempted, incomplete, and withdrawn audit states under the approved protocol. Formative interviews may revise the instrument but may not produce a public comprehension rate or promote the UI.
"""
    (output / "PRIVATE_OPERATIONS_README.md").write_text(private_readme, encoding="utf-8")

    scripts = {
        "SESSION_SCRIPT_EN.md": """# Standardized session script — English

Private formative-research template; adapt only after applicable ethics review.

1. Confirm consent and remind the participant that they may pause or withdraw.
2. Explain that the product is being evaluated, not the participant; do not explain statistical terms.
3. Show the assigned scenario. Ask the participant to think aloud and answer using the displayed options.
4. Lock the response before clarification. Say: “Thank you. I will now ask how you interpreted the screen.”
5. Use only the registered retrospective probe. Do not reveal the answer until all unaided tasks for that scenario are complete.
6. Ask about visual, keyboard, zoom, screen-reader, terminology, and translation barriers as applicable; moderator observation is not accessibility certification.
7. Debrief, remind the participant how withdrawal is handled, and follow the approved retention procedure.
""",
        "SESSION_SCRIPT_JA.md": """# 標準セッションスクリプト — 日本語

privateな形成的研究用templateです。適用される倫理審査後にのみ調整してください。

1. 同意を確認し、いつでも中断・撤回できることを伝えます。
2. 評価対象は製品であり参加者ではないと説明します。統計用語は説明しません。
3. 割り当てたscenarioを提示し、think-aloudと表示選択肢による回答を依頼します。
4. 説明や確認をする前に回答をlockします。「ありがとうございます。ここから画面をどう解釈したか質問します」と伝えます。
5. 登録済みretrospective probeだけを使用します。そのscenarioの非誘導taskが終わるまで正解を示しません。
6. 必要に応じて視覚、keyboard、zoom、screen reader、用語、翻訳上の障壁を確認します。moderator観察はaccessibility認証ではありません。
7. debriefを行い、撤回時の扱いと承認済みretention手順を再確認します。
""",
    }
    for name, text in scripts.items():
        (output / name).write_text(text, encoding="utf-8")

    dictionary = """# Record data dictionary

| Field | Meaning | Constraint |
|---|---|---|
| `SessionSlotId` | Balanced schedule slot, not a person identifier | Registered values only |
| `ExperienceStratum` | Recruitment target stratum | `novice_target` or `experienced_target`; not an achieved proficiency measure |
| `ResponseCode` | Locked option code | Must occur in frozen item options |
| `ResponseLocked` | Unassisted answer locked | Required for completed rows |
| `ProbeStartedAfterLock` | Retrospective probe timing confirmation | Required for completed rows |
| `ResponseTimeSeconds` | Item response duration | Finite and nonnegative; not a speed norm |
| `Confidence` | Self-report | Integer 1–5; not ability |
| `IssueCodesJSON` | Controlled qualitative issue codes | JSON list from registered taxonomy |
| `IssueSeverity` | Formative issue severity | `none`, `observation`, `minor`, `major`, or `critical` |
| `ParaphrasedObservation` | PII-minimized paraphrase | Never a raw quotation or contact field |
| `PIIReviewed` | Manual review attestation | Required for any non-empty paraphrase; not an automated privacy guarantee |
| `Withdrawn` | Withdrawal marker | Content handling follows approved consent/protocol |
| `CompletionState` | Row state | `not_started`, `completed`, or `withdrawn` as applicable |
"""
    (output / "RECORD_DATA_DICTIONARY.md").write_text(dictionary, encoding="utf-8")


def run_tests(output: Path) -> dict[str, object]:
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-q",
        "tests/test_cmle_one_click_cognitive_interview.py",
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
    return {"passed": completed.returncode == 0, "returncode": completed.returncode, "command": command}


def write_critical_review(output: Path, results: dict[str, object]) -> None:
    text = f"""# CMLE cognitive-interview operations critical review

## Decision

The private no-human-data operations contract **{'passed' if results['contract_passed'] else 'failed'}**. This means the schedule, packet separation, blank record schema, and synthetic validator are ready for governed formative use. It is not ethics approval, recruitment authorization, a usability result, or a public-interface gate.

## Design strengths

Within each language, all ten unordered pairs of five analytical states occur once. Each state appears four times, twice in each position and twice in each target experience stratum. English and Japanese share pair/order/stratum codes. Twenty participant packet files contain 20 tasks each and no answer, rationale, criticality, misconception, or moderator-probe fields. The private moderator guide requires response lock before probing.

## Fail-closed data handling

The 400-row blank record template contains no responses, timings, confidence values, issue codes, notes, or human identifiers. Exact-column validation rejects unexpected and direct/contact identifier fields. It also rejects assignment drift, missing/duplicate items, unknown options or issue codes, inconsistent severity, an unreviewed note, a probe not confirmed after response lock, and invalid time/confidence values. These checks reduce accidental leakage; a `PIIReviewed` flag is a manual attestation, not automated de-identification.

## Residual threats

- Balanced slots do not balance who is recruited, language proficiency, moderator behavior, attrition, learning, or device/accessibility conditions.
- Two ten-item scenarios may create fatigue or learning; timing and burden must be examined during formative interviews rather than silently shortening the instrument.
- `novice_target` and `experienced_target` are recruitment intentions, not measured expertise or evidence of comparable groups.
- Qualitative coding can be inconsistent across moderators; training, double-coding, and an adjudication log remain necessary if inference is attempted.
- Moderator observation cannot replace browser, keyboard, screen-reader, zoom, contrast, or other assistive-technology testing.
- Formative interviews must not be converted into a comprehension prevalence estimate, a Japanese–English equivalence claim, or a public promotion decision.

## Next gate

After applicable institutional approval, conduct the sessions privately and log every instrument revision. Then freeze and prospectively register a new pilot instrument; do not pool responses across revisions. Human participants remain zero and `PublicSurfaceEnabled=false` in this evidence.
"""
    (output / "COGNITIVE_INTERVIEW_OPERATIONS_CRITICAL_REVIEW.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    plan = validate_plan()
    if args.output.exists():
        raise FileExistsError(f"Evidence output exists: {args.output}")
    args.output.mkdir(parents=True)

    participant, key = load_source_materials()
    cases = participant["CaseId"].drop_duplicates().tolist()
    schedule = build_cognitive_interview_schedule(cases)
    tasks = build_participant_session_tasks(participant, schedule)
    guide = build_moderator_probe_guide(key)
    template = build_blank_cognitive_record_template(tasks)
    case_audit, pair_audit = schedule_audits(schedule)
    packet_manifest = write_session_packets(tasks, args.output)
    validator_ledger, validator_failures = run_synthetic_validator(
        participant, key, schedule, template
    )
    human_gate = cognitive_interview_human_gate_status()

    write_csv(schedule, args.output / "balanced_session_schedule.csv")
    write_csv(case_audit, args.output / "case_balance_audit.csv")
    write_csv(pair_audit, args.output / "pair_assignment_audit.csv")
    write_csv(tasks, args.output / "all_participant_session_tasks.csv")
    write_csv(guide, args.output / "PRIVATE_MODERATOR_GUIDE.csv")
    write_csv(template, args.output / "blank_cognitive_interview_record_template.csv")
    write_csv(packet_manifest, args.output / "participant_packet_manifest.csv")
    write_csv(validator_ledger, args.output / "synthetic_validator_ledger.csv")
    write_csv(validator_failures, args.output / "synthetic_validator_failure_details.csv")
    write_csv(human_gate, args.output / "human_gate_status.csv")
    write_operations_documents(args.output)
    tests = run_tests(args.output)

    blank_validation = validate_cognitive_interview_records(
        template, key, schedule, template_only=True
    )
    case_balance = bool(
        len(case_audit) == 10
        and case_audit[["Exposures", "FirstPosition", "SecondPosition", "NoviceTarget", "ExperiencedTarget"]]
        .eq([4, 2, 2, 2, 2])
        .all()
        .all()
    )
    pair_balance = bool(
        len(pair_audit) == 20
        and pair_audit["AssignmentRows"].eq(2).all()
        and pair_audit.groupby("Language")["UnorderedPair"].nunique().eq(10).all()
        and pair_audit.groupby(["Language", "ExperienceStratum"]).size().eq(5).all()
    )
    parity = pair_audit.pivot(
        index="PairId", columns="Language", values=["ExperienceStratum", "FirstCase", "SecondCase"]
    )
    language_parity = all(
        parity[(field, "en")].equals(parity[(field, "ja")])
        for field in ("ExperienceStratum", "FirstCase", "SecondCase")
    )
    blank_human_columns = [
        "ResponseCode",
        "ResponseTimeSeconds",
        "Confidence",
        "IssueCodesJSON",
        "ParaphrasedObservation",
    ]
    gates = {
        "identity_passed": True,
        "schedule_balance_passed": bool(
            len(schedule) == 40
            and schedule["SessionSlotId"].nunique() == 20
            and case_balance
            and pair_balance
            and language_parity
        ),
        "participant_packet_passed": bool(
            len(tasks) == 400
            and tasks.groupby("SessionSlotId").size().eq(20).all()
            and not RESEARCHER_ONLY_FIELDS.intersection(tasks.columns)
            and len(packet_manifest) == 20
            and packet_manifest["TaskRows"].eq(20).all()
            and packet_manifest["CaseCount"].eq(2).all()
            and packet_manifest["ResearcherFieldCount"].eq(0).all()
        ),
        "moderator_guide_passed": bool(
            len(guide) == 100
            and guide["ProbeTiming"].eq("after_response_locked_only").all()
            and guide["ModeratorProbe"].str.len().gt(0).all()
            and {"CorrectOption", "Critical", "MisconceptionCode", "Rationale"}.issubset(guide.columns)
        ),
        "blank_record_template_passed": bool(
            len(template) == 400
            and blank_validation["valid"]
            and all(
                template[column].map(lambda value: str(value).strip() == "").all()
                for column in blank_human_columns
            )
            and template["CompletionState"].eq("not_started").all()
        ),
        "synthetic_validator_passed": bool(
            len(validator_ledger) == 14 and validator_ledger["ContractMatch"].all()
        ),
        "human_and_public_withheld_passed": bool(
            len(human_gate) == 1
            and int(human_gate.iloc[0]["HumanParticipants"]) == 0
            and not bool(human_gate.iloc[0]["RecruitmentAuthorizedByRepository"])
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
        "schedule_rows": len(schedule),
        "planned_session_slots": schedule["SessionSlotId"].nunique(),
        "participant_task_rows": len(tasks),
        "private_moderator_guide_rows": len(guide),
        "blank_record_rows": len(template),
        "synthetic_validator_cases": len(validator_ledger),
        "human_participants": 0,
        "human_comprehension_status": "not_assessed",
        "recruitment_authorized_by_repository": False,
        "public_surface_enabled": False,
        "selected_tests": tests,
    }
    write_critical_review(args.output, results)
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
            "contract_interpretation": "private_no_human_data_operations_readiness_only",
            "implementation_sha256": {
                "mfrm_app/cmle_one_click_cognitive_interview.py": sha256_file(
                    ROOT / "mfrm_app/cmle_one_click_cognitive_interview.py"
                ),
                "tests/test_cmle_one_click_cognitive_interview.py": sha256_file(
                    ROOT / "tests/test_cmle_one_click_cognitive_interview.py"
                ),
                "validation/cmle_one_click_cognitive_interview_operations.py": sha256_file(Path(__file__)),
            },
            "results": results,
            "output_sha256": {
                name: sha256_file(args.output / name) for name in output_files
            },
            "interpretation": {
                "ethics_approval": "not provided by this repository",
                "human_data": "none",
                "formative_comprehension": "not assessed",
                "language_equivalence": "not assessed",
                "accessibility": "not assessed by real users or assistive technologies",
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
