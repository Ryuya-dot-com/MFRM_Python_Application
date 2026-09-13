"""Fail-closed browser accessibility acceptance contract.

This module defines what must be exercised once an interactive browser is
available.  It intentionally does not infer browser acceptance from Streamlit
AppTest, CSS inspection, or the existence of a runbook.  Every required
case/check pair remains ``NOT_RUN`` until portable evidence is attached.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
import hashlib
import json
from pathlib import Path, PurePosixPath
from typing import Iterable, Mapping


BROWSER_ACCEPTANCE_SCHEMA_VERSION = "mfrm_browser_accessibility_acceptance_v1"
BROWSER_EVIDENCE_SCHEMA_VERSION = "mfrm_browser_accessibility_evidence_v1"


class BrowserAcceptanceError(ValueError):
    """Raised when the acceptance catalog or evidence is not trustworthy."""


class BrowserAcceptanceStatus(str, Enum):
    NOT_RUN = "NOT_RUN"
    PASS = "PASS"
    FAIL = "FAIL"
    BLOCKED = "BLOCKED"
    NOT_APPLICABLE = "NOT_APPLICABLE"


@dataclass(frozen=True, slots=True)
class BrowserProfile:
    profile_id: str
    css_width: int
    css_height: int
    zoom_percent: int
    input_mode: str
    reduced_motion: bool = False
    screen_reader: bool = False


@dataclass(frozen=True, slots=True)
class BrowserCheck:
    check_id: str
    criterion: str
    evidence_types: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class BrowserTask:
    task_id: str
    task: str
    profile_ids: tuple[str, ...]
    required_check_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class BrowserAcceptanceCase:
    case_id: str
    task_id: str
    locale: str
    profile_id: str
    required_check_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class BrowserEvidenceRecord:
    case_id: str
    check_id: str
    status: BrowserAcceptanceStatus = BrowserAcceptanceStatus.NOT_RUN
    evidence: Mapping[str, str] = field(default_factory=dict)
    observed: str = ""
    tested_at_utc: str = ""
    browser_build: str = ""
    reviewer: str = ""
    schema_version: str = BROWSER_EVIDENCE_SCHEMA_VERSION

    def to_row(self) -> dict[str, str]:
        return {
            "SchemaVersion": self.schema_version,
            "CatalogFingerprintSHA256": _catalog_fingerprint(),
            "CaseID": self.case_id,
            "CheckID": self.check_id,
            "Status": self.status.value,
            "EvidenceJSON": json.dumps(
                dict(sorted(self.evidence.items())),
                ensure_ascii=False,
                separators=(",", ":"),
            ),
            "Observed": self.observed,
            "TestedAtUTC": self.tested_at_utc,
            "BrowserBuild": self.browser_build,
            "Reviewer": self.reviewer,
        }

    @classmethod
    def from_row(cls, row: Mapping[str, object]) -> "BrowserEvidenceRecord":
        required = {
            "SchemaVersion",
            "CatalogFingerprintSHA256",
            "CaseID",
            "CheckID",
            "Status",
            "EvidenceJSON",
            "Observed",
            "TestedAtUTC",
            "BrowserBuild",
            "Reviewer",
        }
        if set(row) != required:
            raise BrowserAcceptanceError("Browser evidence row shape mismatch")
        if str(row["CatalogFingerprintSHA256"]) != _catalog_fingerprint():
            raise BrowserAcceptanceError("Browser evidence catalog fingerprint mismatch")
        try:
            evidence = json.loads(str(row["EvidenceJSON"]))
        except json.JSONDecodeError as exc:
            raise BrowserAcceptanceError("EvidenceJSON must be valid JSON") from exc
        if not isinstance(evidence, dict):
            raise BrowserAcceptanceError("EvidenceJSON must encode an object")
        try:
            status = BrowserAcceptanceStatus(str(row["Status"]))
        except ValueError as exc:
            raise BrowserAcceptanceError("Unknown browser evidence status") from exc
        return cls(
            schema_version=str(row["SchemaVersion"]),
            case_id=str(row["CaseID"]),
            check_id=str(row["CheckID"]),
            status=status,
            evidence={str(key): str(value) for key, value in evidence.items()},
            observed=str(row["Observed"]),
            tested_at_utc=str(row["TestedAtUTC"]),
            browser_build=str(row["BrowserBuild"]),
            reviewer=str(row["Reviewer"]),
        )


BROWSER_PROFILES: tuple[BrowserProfile, ...] = (
    BrowserProfile("desktop_keyboard", 1440, 900, 100, "keyboard"),
    BrowserProfile("narrow_keyboard", 320, 568, 100, "keyboard"),
    BrowserProfile("narrow_touch", 320, 568, 100, "touch"),
    BrowserProfile("zoom_200_keyboard", 640, 450, 200, "keyboard"),
    BrowserProfile("zoom_400_keyboard", 320, 225, 400, "keyboard"),
    BrowserProfile(
        "reduced_motion_keyboard", 1440, 900, 100, "keyboard", reduced_motion=True
    ),
    BrowserProfile(
        "screen_reader_keyboard", 1440, 900, 100, "keyboard", screen_reader=True
    ),
)


BROWSER_CHECKS: tuple[BrowserCheck, ...] = (
    BrowserCheck(
        "one_primary_action",
        "Exactly one visually primary next action is exposed for the current phase.",
        ("dom_snapshot", "screenshot"),
    ),
    BrowserCheck(
        "accessible_names",
        "Every required control has a unique, meaningful accessible name.",
        ("accessibility_tree", "dom_snapshot"),
    ),
    BrowserCheck(
        "keyboard_completion",
        "The task completes with Tab, Shift+Tab, Enter, Space, and arrow keys only.",
        ("keyboard_trace", "dom_snapshot"),
    ),
    BrowserCheck(
        "visible_focus",
        "Keyboard focus is never lost and its indicator remains visibly distinct.",
        ("keyboard_trace", "screenshot"),
    ),
    BrowserCheck(
        "focus_after_rerun",
        "A rerun does not move focus to an unrelated control or hide status context.",
        ("keyboard_trace", "dom_snapshot"),
    ),
    BrowserCheck(
        "heading_reading_order",
        "Heading levels and accessibility-tree order match the visual task sequence.",
        ("accessibility_tree", "dom_snapshot"),
    ),
    BrowserCheck(
        "no_hidden_required_control",
        "No required control is clipped, collapsed without signposting, or hover-only.",
        ("dom_snapshot", "screenshot"),
    ),
    BrowserCheck(
        "no_horizontal_page_scroll",
        "The page reflows without page-level horizontal scrolling.",
        ("dom_snapshot", "screenshot"),
    ),
    BrowserCheck(
        "sticky_not_obscure",
        "Sticky result navigation does not cover warnings, content, or actions.",
        ("dom_snapshot", "screenshot"),
    ),
    BrowserCheck(
        "non_color_status",
        "Status and severity remain understandable without color perception.",
        ("accessibility_tree", "screenshot"),
    ),
    BrowserCheck(
        "reduced_motion",
        "Reduced-motion preference removes non-essential animation and smooth scrolling.",
        ("dom_snapshot", "moderator_note"),
    ),
    BrowserCheck(
        "touch_target_size",
        "Required touch targets meet the supported WCAG 2.2 AA target-size scope.",
        ("dom_snapshot", "screenshot"),
    ),
    BrowserCheck(
        "screen_reader_announcement",
        "Task changes and run completion are announced once with sufficient context.",
        ("accessibility_tree", "moderator_note"),
    ),
    BrowserCheck(
        "long_japanese_reflow",
        "Japanese labels remain readable without clipping or ambiguous truncation.",
        ("dom_snapshot", "screenshot"),
    ),
    BrowserCheck(
        "privacy_origin_message",
        "Synthetic and user-data sources expose the correct non-color privacy severity.",
        ("accessibility_tree", "screenshot"),
    ),
    BrowserCheck(
        "raw_rows_collapsed",
        "Raw response rows are not exposed before an explicit disclosure action.",
        ("dom_snapshot", "screenshot"),
    ),
    BrowserCheck(
        "analysis_identity_stable",
        "Presentation-only actions preserve the current AnalysisID and do not refit.",
        ("app_identity_record", "moderator_note"),
    ),
    BrowserCheck(
        "weight_explanation",
        "Guided mode explains equal weighting and never offers Weight as a facet.",
        ("dom_snapshot", "screenshot"),
    ),
    BrowserCheck(
        "download_keyboard",
        "The result archive can be reached and initiated without pointer input.",
        ("keyboard_trace", "download_artifact"),
    ),
    BrowserCheck(
        "help_return_focus",
        "Contextual Help returns to the registered source without stealing moved focus.",
        ("keyboard_trace", "dom_snapshot", "screenshot"),
    ),
)


_KEYBOARD_PROFILE_CHECKS = (
    "keyboard_completion",
    "visible_focus",
    "focus_after_rerun",
    "heading_reading_order",
)
_NARROW_PROFILE_CHECKS = ("no_horizontal_page_scroll", "sticky_not_obscure")


BROWSER_TASKS: tuple[BrowserTask, ...] = (
    BrowserTask(
        "landing_sample",
        "Choose the sample-learning route from a new session.",
        ("desktop_keyboard", "narrow_keyboard", "zoom_400_keyboard"),
        ("one_primary_action", "accessible_names", "privacy_origin_message"),
    ),
    BrowserTask(
        "guided_sample_complete",
        "Complete the five-step sample route and save its non-scientific checkpoint.",
        (
            "desktop_keyboard",
            "narrow_keyboard",
            "reduced_motion_keyboard",
            "screen_reader_keyboard",
        ),
        (
            "one_primary_action",
            "accessible_names",
            "no_hidden_required_control",
            "non_color_status",
            "analysis_identity_stable",
        ),
    ),
    BrowserTask(
        "own_data_paste",
        "Paste synthetic fixture rows, confirm mapping/readiness, and run Guided defaults.",
        ("desktop_keyboard", "narrow_keyboard", "zoom_200_keyboard"),
        (
            "one_primary_action",
            "accessible_names",
            "privacy_origin_message",
            "raw_rows_collapsed",
        ),
    ),
    BrowserTask(
        "guided_weight_explanation",
        "Load a synthetic Weight column and verify the Guided equal-weight explanation.",
        ("desktop_keyboard", "narrow_keyboard"),
        ("accessible_names", "weight_explanation", "no_hidden_required_control"),
    ),
    BrowserTask(
        "advanced_anchor_setup",
        "Switch to Advanced, configure a synthetic anchor table, and return to setup.",
        ("desktop_keyboard", "narrow_keyboard", "zoom_200_keyboard"),
        ("accessible_names", "no_hidden_required_control", "analysis_identity_stable"),
    ),
    BrowserTask(
        "sparse_warning_recovery",
        "Run the sparse sample and follow its highest-priority reversible action.",
        ("desktop_keyboard", "narrow_keyboard"),
        ("one_primary_action", "non_color_status", "accessible_names"),
    ),
    BrowserTask(
        "result_navigation_download",
        "Navigate the fitted result and obtain the ZIP archive.",
        (
            "desktop_keyboard",
            "narrow_keyboard",
            "narrow_touch",
            "zoom_400_keyboard",
            "screen_reader_keyboard",
        ),
        (
            "one_primary_action",
            "accessible_names",
            "sticky_not_obscure",
            "download_keyboard",
            "non_color_status",
        ),
    ),
    BrowserTask(
        "contextual_help_return",
        "Open contextual Help and return to the exact registered source heading.",
        (
            "desktop_keyboard",
            "narrow_keyboard",
            "reduced_motion_keyboard",
            "screen_reader_keyboard",
        ),
        ("accessible_names", "help_return_focus", "analysis_identity_stable"),
    ),
)


def _catalog_maps() -> tuple[dict[str, BrowserProfile], dict[str, BrowserCheck], dict[str, BrowserTask]]:
    return (
        {row.profile_id: row for row in BROWSER_PROFILES},
        {row.check_id: row for row in BROWSER_CHECKS},
        {row.task_id: row for row in BROWSER_TASKS},
    )


def validate_browser_acceptance_catalog() -> None:
    profiles, checks, tasks = _catalog_maps()
    if len(profiles) != len(BROWSER_PROFILES):
        raise BrowserAcceptanceError("Duplicate browser profile ID")
    if len(checks) != len(BROWSER_CHECKS):
        raise BrowserAcceptanceError("Duplicate browser check ID")
    if len(tasks) != len(BROWSER_TASKS):
        raise BrowserAcceptanceError("Duplicate browser task ID")
    for profile in BROWSER_PROFILES:
        if profile.css_width < 320 or profile.css_height < 200:
            raise BrowserAcceptanceError("Browser profiles must cover at least 320 CSS px")
        if profile.zoom_percent not in {100, 200, 400}:
            raise BrowserAcceptanceError("Unsupported browser zoom profile")
        if profile.input_mode not in {"keyboard", "touch"}:
            raise BrowserAcceptanceError("Unknown browser input mode")
    for check in BROWSER_CHECKS:
        if not check.evidence_types or len(set(check.evidence_types)) != len(check.evidence_types):
            raise BrowserAcceptanceError("Checks require unique evidence types")
    for task in BROWSER_TASKS:
        unknown_profiles = set(task.profile_ids).difference(profiles)
        unknown_checks = set(task.required_check_ids).difference(checks)
        if unknown_profiles or unknown_checks:
            raise BrowserAcceptanceError(
                f"Task {task.task_id} has unknown profiles/checks: "
                f"{sorted(unknown_profiles)!r}/{sorted(unknown_checks)!r}"
            )


def _profile_check_ids(profile: BrowserProfile) -> tuple[str, ...]:
    checks: list[str] = []
    if profile.input_mode == "keyboard":
        checks.extend(_KEYBOARD_PROFILE_CHECKS)
    if profile.css_width <= 640 or profile.zoom_percent > 100:
        checks.extend(_NARROW_PROFILE_CHECKS)
    if profile.input_mode == "touch":
        checks.append("touch_target_size")
    if profile.reduced_motion:
        checks.append("reduced_motion")
    if profile.screen_reader:
        checks.append("screen_reader_announcement")
    return tuple(checks)


def build_browser_acceptance_cases() -> tuple[BrowserAcceptanceCase, ...]:
    validate_browser_acceptance_catalog()
    profiles, checks, _ = _catalog_maps()
    order = {check.check_id: index for index, check in enumerate(BROWSER_CHECKS)}
    cases: list[BrowserAcceptanceCase] = []
    for task in BROWSER_TASKS:
        for locale in ("en", "ja"):
            for profile_id in task.profile_ids:
                required = set(task.required_check_ids)
                required.update(_profile_check_ids(profiles[profile_id]))
                if locale == "ja":
                    required.add("long_japanese_reflow")
                required_ids = tuple(sorted(required, key=order.__getitem__))
                if not set(required_ids).issubset(checks):
                    raise BrowserAcceptanceError("Acceptance case references an unknown check")
                cases.append(BrowserAcceptanceCase(
                    case_id=f"{task.task_id}__{locale}__{profile_id}",
                    task_id=task.task_id,
                    locale=locale,
                    profile_id=profile_id,
                    required_check_ids=required_ids,
                ))
    case_ids = [case.case_id for case in cases]
    if len(case_ids) != len(set(case_ids)):
        raise BrowserAcceptanceError("Duplicate browser acceptance case ID")
    return tuple(cases)


def blank_browser_evidence_records() -> tuple[BrowserEvidenceRecord, ...]:
    return tuple(
        BrowserEvidenceRecord(case_id=case.case_id, check_id=check_id)
        for case in build_browser_acceptance_cases()
        for check_id in case.required_check_ids
    )


def _validate_portable_evidence_ref(value: str) -> None:
    ref = str(value).strip()
    path = PurePosixPath(ref)
    if (
        not ref
        or ref.startswith("/")
        or "://" in ref
        or "\\" in ref
        or ".." in path.parts
    ):
        raise BrowserAcceptanceError("Evidence references must be portable relative paths")


def validate_browser_evidence_record(record: BrowserEvidenceRecord) -> None:
    if record.schema_version != BROWSER_EVIDENCE_SCHEMA_VERSION:
        raise BrowserAcceptanceError("Unsupported browser evidence schema version")
    cases = {case.case_id: case for case in build_browser_acceptance_cases()}
    checks = {check.check_id: check for check in BROWSER_CHECKS}
    case = cases.get(record.case_id)
    check = checks.get(record.check_id)
    if case is None or check is None or record.check_id not in case.required_check_ids:
        raise BrowserAcceptanceError("Evidence does not match a required case/check pair")
    evidence = dict(record.evidence)
    allowed_evidence = set(check.evidence_types).union({"moderator_note"})
    if not set(evidence).issubset(allowed_evidence):
        raise BrowserAcceptanceError("Evidence type is not declared for this check")
    for ref in evidence.values():
        _validate_portable_evidence_ref(ref)
    if record.status is BrowserAcceptanceStatus.PASS:
        if set(evidence) != set(check.evidence_types):
            raise BrowserAcceptanceError("PASS requires every declared evidence type")
    elif record.status is BrowserAcceptanceStatus.FAIL:
        if not evidence:
            raise BrowserAcceptanceError("FAIL requires retained failure evidence")
    elif record.status is BrowserAcceptanceStatus.BLOCKED:
        if "moderator_note" not in evidence:
            raise BrowserAcceptanceError("BLOCKED requires a retained moderator note")
    elif record.status is BrowserAcceptanceStatus.NOT_RUN:
        if evidence:
            raise BrowserAcceptanceError("NOT_RUN cannot carry browser evidence")
        if any(
            str(value).strip()
            for value in (
                record.observed,
                record.tested_at_utc,
                record.browser_build,
                record.reviewer,
            )
        ):
            raise BrowserAcceptanceError("NOT_RUN cannot carry execution metadata")
    else:
        raise BrowserAcceptanceError("Required checks cannot be NOT_APPLICABLE")
    if record.status is not BrowserAcceptanceStatus.NOT_RUN:
        if not all(
            str(value).strip()
            for value in (
                record.observed,
                record.tested_at_utc,
                record.browser_build,
                record.reviewer,
            )
        ):
            raise BrowserAcceptanceError("Executed evidence requires audit metadata")
        try:
            tested_at = datetime.fromisoformat(record.tested_at_utc.replace("Z", "+00:00"))
        except ValueError as exc:
            raise BrowserAcceptanceError("TestedAtUTC must be an ISO-8601 UTC timestamp") from exc
        if (
            not record.tested_at_utc.endswith("Z")
            or tested_at.tzinfo is None
            or tested_at.utcoffset() != timezone.utc.utcoffset(tested_at)
        ):
            raise BrowserAcceptanceError("TestedAtUTC must use the UTC Z suffix")


def summarize_browser_acceptance(
    records: Iterable[BrowserEvidenceRecord],
) -> dict[str, object]:
    expected = {
        (case.case_id, check_id)
        for case in build_browser_acceptance_cases()
        for check_id in case.required_check_ids
    }
    indexed: dict[tuple[str, str], BrowserEvidenceRecord] = {}
    for record in records:
        validate_browser_evidence_record(record)
        key = (record.case_id, record.check_id)
        if key in indexed:
            raise BrowserAcceptanceError("Duplicate browser evidence record")
        indexed[key] = record
    extra = set(indexed).difference(expected)
    if extra:
        raise BrowserAcceptanceError("Unexpected browser evidence record")
    counts = {status.value: 0 for status in BrowserAcceptanceStatus}
    for record in indexed.values():
        counts[record.status.value] += 1
    missing = len(expected.difference(indexed))
    accepted = (
        missing == 0
        and len(indexed) == len(expected)
        and counts[BrowserAcceptanceStatus.PASS.value] == len(expected)
    )
    return {
        "schema_version": BROWSER_ACCEPTANCE_SCHEMA_VERSION,
        "decision": "ACCEPTED" if accepted else "NOT_ACCEPTED",
        "case_count": len(build_browser_acceptance_cases()),
        "required_record_count": len(expected),
        "received_record_count": len(indexed),
        "missing_record_count": missing,
        "status_counts": counts,
    }


def _catalog_fingerprint() -> str:
    payload = {
        "schema_version": BROWSER_ACCEPTANCE_SCHEMA_VERSION,
        "profiles": [asdict(profile) for profile in BROWSER_PROFILES],
        "checks": [asdict(check) for check in BROWSER_CHECKS],
        "tasks": [asdict(task) for task in BROWSER_TASKS],
        "cases": [asdict(case) for case in build_browser_acceptance_cases()],
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def write_browser_acceptance_bundle(output_dir: str | Path) -> dict[str, object]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    profiles = {profile.profile_id: profile for profile in BROWSER_PROFILES}
    tasks = {task.task_id: task for task in BROWSER_TASKS}
    cases = build_browser_acceptance_cases()
    matrix_fields = [
        "SchemaVersion", "CaseID", "TaskID", "Task", "Locale", "ProfileID",
        "CSSWidth", "CSSHeight", "ZoomPercent", "InputMode", "ReducedMotion",
        "ScreenReader", "RequiredChecks",
    ]
    with (output / "browser_acceptance_matrix.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=matrix_fields)
        writer.writeheader()
        for case in cases:
            profile = profiles[case.profile_id]
            writer.writerow({
                "SchemaVersion": BROWSER_ACCEPTANCE_SCHEMA_VERSION,
                "CaseID": case.case_id,
                "TaskID": case.task_id,
                "Task": tasks[case.task_id].task,
                "Locale": case.locale,
                "ProfileID": case.profile_id,
                "CSSWidth": profile.css_width,
                "CSSHeight": profile.css_height,
                "ZoomPercent": profile.zoom_percent,
                "InputMode": profile.input_mode,
                "ReducedMotion": str(profile.reduced_motion).lower(),
                "ScreenReader": str(profile.screen_reader).lower(),
                "RequiredChecks": "|".join(case.required_check_ids),
            })
    check_fields = ["SchemaVersion", "CheckID", "Criterion", "EvidenceTypes"]
    with (output / "browser_acceptance_checks.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=check_fields)
        writer.writeheader()
        for check in BROWSER_CHECKS:
            writer.writerow({
                "SchemaVersion": BROWSER_ACCEPTANCE_SCHEMA_VERSION,
                "CheckID": check.check_id,
                "Criterion": check.criterion,
                "EvidenceTypes": "|".join(check.evidence_types),
            })
    evidence_fields = list(blank_browser_evidence_records()[0].to_row())
    with (output / "browser_evidence_template.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=evidence_fields)
        writer.writeheader()
        for record in blank_browser_evidence_records():
            writer.writerow(record.to_row())
    summary = summarize_browser_acceptance(blank_browser_evidence_records())
    manifest = {
        **summary,
        "catalog_fingerprint_sha256": _catalog_fingerprint(),
        "files": [
            "browser_acceptance_matrix.csv",
            "browser_acceptance_checks.csv",
            "browser_evidence_template.csv",
        ],
        "boundary": (
            "Templates and static tests are not browser acceptance evidence. "
            "Decision remains NOT_ACCEPTED until every required record is PASS."
        ),
    }
    (output / "browser_acceptance_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def read_browser_evidence_csv(path: str | Path) -> tuple[BrowserEvidenceRecord, ...]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return tuple(
            BrowserEvidenceRecord.from_row(row)
            for row in csv.DictReader(handle)
        )


def write_browser_acceptance_decision(
    evidence_csv: str | Path,
    output_json: str | Path,
) -> dict[str, object]:
    summary = summarize_browser_acceptance(read_browser_evidence_csv(evidence_csv))
    decision = {
        **summary,
        "catalog_fingerprint_sha256": _catalog_fingerprint(),
        "evidence_file": Path(evidence_csv).name,
        "boundary": (
            "ACCEPTED requires every version-matched required record to be PASS "
            "with complete portable evidence and audit metadata."
        ),
    }
    output = Path(output_json)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(decision, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return decision


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Manage browser accessibility acceptance evidence")
    commands = parser.add_subparsers(dest="command", required=True)
    export_parser = commands.add_parser("export", help="Export a fail-closed blank bundle")
    export_parser.add_argument(
        "--output", required=True, help="Directory for the versioned acceptance bundle"
    )
    evaluate_parser = commands.add_parser("evaluate", help="Evaluate completed evidence")
    evaluate_parser.add_argument("--evidence", required=True, help="Completed evidence CSV")
    evaluate_parser.add_argument("--output", required=True, help="Decision JSON path")
    args = parser.parse_args(argv)
    if args.command == "export":
        result = write_browser_acceptance_bundle(args.output)
    else:
        result = write_browser_acceptance_decision(args.evidence, args.output)
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    if args.command == "export":
        return 0
    return 0 if result["decision"] == "ACCEPTED" else 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
