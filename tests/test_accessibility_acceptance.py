from __future__ import annotations

import csv
from dataclasses import replace
import inspect
import json
from pathlib import Path

import pandas as pd
import pytest

from mfrm_app import accessibility_acceptance as acceptance
import streamlit_app as app


def test_browser_acceptance_catalog_is_complete_and_stable() -> None:
    acceptance.validate_browser_acceptance_catalog()
    cases = acceptance.build_browser_acceptance_cases()

    assert len(cases) == 52
    assert len({case.case_id for case in cases}) == len(cases)
    assert {case.locale for case in cases} == {"en", "ja"}
    assert {case.task_id for case in cases} == {
        task.task_id for task in acceptance.BROWSER_TASKS
    }
    assert {case.profile_id for case in cases} == {
        profile.profile_id for profile in acceptance.BROWSER_PROFILES
    }


def test_profile_and_locale_requirements_are_added_to_every_case() -> None:
    profiles = {
        profile.profile_id: profile for profile in acceptance.BROWSER_PROFILES
    }
    for case in acceptance.build_browser_acceptance_cases():
        required = set(case.required_check_ids)
        profile = profiles[case.profile_id]
        if profile.input_mode == "keyboard":
            assert {"keyboard_completion", "visible_focus"}.issubset(required)
        if profile.css_width <= 640 or profile.zoom_percent > 100:
            assert {
                "no_horizontal_page_scroll",
                "sticky_not_obscure",
            }.issubset(required)
        if profile.input_mode == "touch":
            assert "touch_target_size" in required
        if profile.reduced_motion:
            assert "reduced_motion" in required
        if profile.screen_reader:
            assert "screen_reader_announcement" in required
        if case.locale == "ja":
            assert "long_japanese_reflow" in required


def test_blank_browser_evidence_is_explicitly_not_accepted() -> None:
    records = acceptance.blank_browser_evidence_records()
    summary = acceptance.summarize_browser_acceptance(records)

    assert summary["decision"] == "NOT_ACCEPTED"
    assert summary["missing_record_count"] == 0
    assert summary["status_counts"]["NOT_RUN"] == len(records)
    assert summary["status_counts"]["PASS"] == 0


def _passing_records() -> tuple[acceptance.BrowserEvidenceRecord, ...]:
    checks = {
        check.check_id: check for check in acceptance.BROWSER_CHECKS
    }
    records = []
    for blank in acceptance.blank_browser_evidence_records():
        evidence = {
            evidence_type: (
                f"evidence/{blank.case_id}/{blank.check_id}_{evidence_type}.txt"
            )
            for evidence_type in checks[blank.check_id].evidence_types
        }
        records.append(replace(
            blank,
            status=acceptance.BrowserAcceptanceStatus.PASS,
            evidence=evidence,
            observed="Criterion observed in the rendered browser task.",
            tested_at_utc="2026-08-10T12:00:00Z",
            browser_build="Chromium test-build",
            reviewer="accessibility-reviewer",
        ))
    return tuple(records)


def test_only_complete_evidence_can_produce_accepted_decision() -> None:
    passing = _passing_records()
    summary = acceptance.summarize_browser_acceptance(passing)
    assert summary["decision"] == "ACCEPTED"
    assert summary["status_counts"]["PASS"] == len(passing)

    incomplete = acceptance.summarize_browser_acceptance(passing[:-1])
    assert incomplete["decision"] == "NOT_ACCEPTED"
    assert incomplete["missing_record_count"] == 1


def test_pass_and_failure_evidence_are_validated_fail_closed() -> None:
    blank = acceptance.blank_browser_evidence_records()[0]
    with pytest.raises(acceptance.BrowserAcceptanceError, match="every declared"):
        acceptance.validate_browser_evidence_record(replace(
            blank,
            status=acceptance.BrowserAcceptanceStatus.PASS,
            observed="Observed",
            tested_at_utc="2026-08-10T12:00:00Z",
            browser_build="Browser",
            reviewer="Reviewer",
        ))

    with pytest.raises(acceptance.BrowserAcceptanceError, match="portable"):
        acceptance.validate_browser_evidence_record(replace(
            _passing_records()[0],
            evidence={
                key: "/private/evidence.txt"
                for key in _passing_records()[0].evidence
            },
        ))


def test_evidence_csv_row_roundtrip_preserves_contract() -> None:
    original = _passing_records()[0]
    row = original.to_row()
    restored = acceptance.BrowserEvidenceRecord.from_row(row)
    assert restored == original
    acceptance.validate_browser_evidence_record(restored)

    stale = dict(row)
    stale["CatalogFingerprintSHA256"] = "0" * 64
    with pytest.raises(acceptance.BrowserAcceptanceError, match="fingerprint"):
        acceptance.BrowserEvidenceRecord.from_row(stale)


def test_browser_acceptance_bundle_is_versioned_and_fail_closed(tmp_path: Path) -> None:
    manifest = acceptance.write_browser_acceptance_bundle(tmp_path)

    assert manifest["decision"] == "NOT_ACCEPTED"
    assert manifest["case_count"] == 52
    assert len(manifest["catalog_fingerprint_sha256"]) == 64
    assert (tmp_path / "browser_acceptance_matrix.csv").exists()
    assert (tmp_path / "browser_acceptance_checks.csv").exists()
    assert (tmp_path / "browser_evidence_template.csv").exists()
    stored = json.loads(
        (tmp_path / "browser_acceptance_manifest.json").read_text(encoding="utf-8")
    )
    assert stored == manifest

    with (tmp_path / "browser_acceptance_matrix.csv").open(
        encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 52
    assert {row["Locale"] for row in rows} == {"en", "ja"}

    decision = acceptance.write_browser_acceptance_decision(
        tmp_path / "browser_evidence_template.csv",
        tmp_path / "browser_acceptance_decision.json",
    )
    assert decision["decision"] == "NOT_ACCEPTED"
    assert decision["status_counts"]["NOT_RUN"] == decision["required_record_count"]


def test_completed_evidence_csv_evaluates_to_accepted(tmp_path: Path) -> None:
    records = _passing_records()
    evidence_path = tmp_path / "browser_evidence_completed.csv"
    with evidence_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0].to_row()))
        writer.writeheader()
        for record in records:
            writer.writerow(record.to_row())

    decision = acceptance.write_browser_acceptance_decision(
        evidence_path,
        tmp_path / "decision.json",
    )
    assert decision["decision"] == "ACCEPTED"
    assert decision["status_counts"]["PASS"] == len(records)


def test_compact_html_tables_have_accessible_name_and_header_scope() -> None:
    frame = pd.DataFrame({"A": ["<private>"], "B": [2]}, index=["row-1"])
    html = app._accessible_dataframe_html(
        frame,
        include_index=True,
        accessible_label='Review "table"',
    )

    assert 'aria-label="Review &quot;table&quot;"' in html
    assert '<th scope="col">A</th>' in html
    assert '<th scope="row">row-1</th>' in html
    assert "&lt;private&gt;" in html
    assert "<private>" not in html


def test_static_accessibility_css_covers_interactive_roles_and_preferences() -> None:
    source = inspect.getsource(app._inject_desktop_readability_css)

    for selector in (
        "a:focus-visible",
        "summary:focus-visible",
        '[role="button"]:focus-visible',
        '[role="combobox"]:focus-visible',
        '[role="switch"]:focus-visible',
        '[tabindex]:not([tabindex="-1"]):focus-visible',
    ):
        assert selector in source
    assert "var(--text-color, CanvasText)" in source
    assert "var(--background-color, Canvas)" in source
    assert "@media (forced-colors: active)" in source
    assert "@media (prefers-reduced-motion: reduce)" in source
    assert "@media (pointer: coarse)" in source
    assert "min-height: 2.75rem" in source


def test_app_does_not_hide_required_widget_labels_or_add_custom_hotkeys() -> None:
    source = Path(app.__file__).read_text(encoding="utf-8")
    assert 'label_visibility="collapsed"' not in source
    assert "label_visibility='collapsed'" not in source
    assert "keydown" not in source.casefold()
    assert "keyup" not in source.casefold()
    assert "accesskey=" not in source.casefold()
