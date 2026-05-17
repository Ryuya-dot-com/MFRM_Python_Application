"""Locale parity tests for the bilingual i18n scaffold.

These tests guard the ``locales/{en,ja}.json`` translation dictionaries:

1. Both locales expose the exact same set of leaf keys (no key drift).
2. The same ``{placeholder}`` set appears in every paired translation
   (so ``t("...", release_label=...)`` cannot break in one language).
3. No translated string is empty, which would silently render as a blank
   label in the UI.

Future locales added to ``SUPPORTED_LANGS`` are picked up automatically
via the ``LOCALE_FILES`` glob.
"""

from __future__ import annotations

import json
import string
from pathlib import Path

import pytest


LOCALES_DIR = Path(__file__).resolve().parents[1] / "locales"
LOCALE_FILES = sorted(LOCALES_DIR.glob("*.json"))


def _flatten_leaves(node, prefix: str = "") -> dict:
    """Walk a nested locale dict and return ``{dotted.key: text}`` pairs.

    Keys whose top-level segment starts with ``_`` (e.g. ``_metadata``) are
    skipped so they cannot accidentally collide with user-visible keys.
    """
    out = {}
    if not isinstance(node, dict):
        return out
    for k, v in node.items():
        if not prefix and k.startswith("_"):
            continue
        nk = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            out.update(_flatten_leaves(v, nk))
        elif isinstance(v, str):
            out[nk] = v
    return out


def _placeholders(text: str) -> set:
    """Return the set of ``{name}`` placeholders in a format string."""
    return {field for _, field, _, _ in string.Formatter().parse(text) if field}


def _load_locale(path: Path) -> dict:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def locales() -> dict:
    assert LOCALE_FILES, f"No locale JSON files found under {LOCALES_DIR}"
    return {path.stem: _load_locale(path) for path in LOCALE_FILES}


def test_locale_files_present():
    assert (LOCALES_DIR / "en.json").exists(), "en.json (default locale) is required"
    assert (LOCALES_DIR / "ja.json").exists(), "ja.json (Japanese locale) is required"


def test_locale_metadata_codes_match_filenames(locales):
    for lang, data in locales.items():
        meta = data.get("_metadata", {})
        code = meta.get("code")
        assert code == lang, (
            f"_metadata.code mismatch for {lang}.json: got {code!r}, expected {lang!r}"
        )


def test_locale_key_parity(locales):
    en_keys = set(_flatten_leaves(locales["en"]).keys())
    for lang, data in locales.items():
        if lang == "en":
            continue
        keys = set(_flatten_leaves(data).keys())
        only_in_en = en_keys - keys
        only_in_other = keys - en_keys
        assert not only_in_en and not only_in_other, (
            f"Key drift between en and {lang}: "
            f"only in en={sorted(only_in_en)}, only in {lang}={sorted(only_in_other)}"
        )


def test_locale_placeholder_parity(locales):
    en_leaves = _flatten_leaves(locales["en"])
    for lang, data in locales.items():
        if lang == "en":
            continue
        other_leaves = _flatten_leaves(data)
        for key, en_text in en_leaves.items():
            other_text = other_leaves.get(key, "")
            en_ph = _placeholders(en_text)
            other_ph = _placeholders(other_text)
            assert en_ph == other_ph, (
                f"Placeholder mismatch at {key}: en={sorted(en_ph)}, "
                f"{lang}={sorted(other_ph)}"
            )


def test_locale_no_empty_strings(locales):
    for lang, data in locales.items():
        for key, text in _flatten_leaves(data).items():
            assert text.strip(), f"Empty translation at {lang}:{key}"


def test_japanese_safety_and_guidance_copy_uses_task_centered_terms(locales):
    ja = _flatten_leaves(locales["ja"])
    checked_keys = {
        "sidebar.view_density_essential",
        "sidebar.view_density_full",
        "sidebar.view_density_help",
        "onboarding.banner_steps_label",
        "dimensionality.warning_pca_failed_template",
        "dimensionality.dimtest_caption",
        "dimensionality.dimtest_reject_warning",
        "data_source.upload_size_blocked_template",
        "data_source.paste_parse_error",
        "fit_details.guide_expander",
        "downloads.public_export_mode_label",
        "downloads.public_export_mode_help",
        "downloads.public_export_mode_info",
        "downloads.private_export_mode_warning",
        "resource_preflight.blocked_error",
        "resource_preflight.review_warning",
        "resource_preflight.ok_caption",
        "resource_preflight.details_expander",
        "visuals_top.roadmap_caption",
        "visuals_top.essential_mode_caption",
        "report_top.reports_essential_caption",
        "report_top.checks_essential_caption",
        "report_top.exports_essential_caption",
    }
    missing = checked_keys - set(ja)
    assert not missing, f"Missing Japanese copy keys: {sorted(missing)}"

    joined = "\n".join(ja[key] for key in sorted(checked_keys))
    awkward_fragments = [
        "初めての方へ",
        "\u521d\u5fc3\u8005",
        "簡易",
        "完全/ private",
        "hosted Streamlit UI",
        "View density",
        "再実行すること",
        "指定すること",
        "判断すること",
        "貼り付けてみる",
        "利用できる。",
        "使用すること",
    ]
    for fragment in awkward_fragments:
        assert fragment not in joined, f"Awkward Japanese UI fragment remains: {fragment}"

    assert "公開・共有用エクスポートモード" in ja["downloads.public_export_mode_label"]
    assert "正式な非識別化保証ではありません" in ja["downloads.public_export_mode_help"]
    assert "解釈ガイド" in ja["fit_details.guide_expander"]


def test_help_copy_frames_thresholds_as_diagnostic_guidance(locales):
    en = _flatten_leaves(locales["en"])
    ja = _flatten_leaves(locales["ja"])

    assert "not as an automatic decision rule" in en["help.interpretation_caution"]
    assert "screening heuristics" in en["help.interpretation_caution"]
    assert "自動判定ルールではありません" in ja["help.interpretation_caution"]
    assert "スクリーニング用の目安" in ja["help.interpretation_caution"]

    old_overclaims = [
        "The default and\n  recommended method",
        "Missing data is below **20%**",
        "Low is *desirable*",
        "No significant bias | None needed",
        "Report and interpret substantively",
        "デフォルトかつ推奨される推定法",
        "欠測が **20%** 未満",
        "低いことが *望ましい*",
        "有意な bias なし | 不要",
        "レポートに記載し、実質的に解釈する",
    ]
    help_text = "\n".join(
        value for key, value in {**en, **ja}.items() if key.startswith("help.")
    )
    for phrase in old_overclaims:
        assert phrase not in help_text

    assert "MFRM_Visual_Evidence_Binder.zip" in en["help.quick_start_body"]
    assert "visual_qa_preflight.csv" in ja["help.quick_start_body"]
    assert "Holm/BH" in en["help.analysis_workflow_body"]
    assert "Holm / BH" in ja["help.analysis_workflow_body"]
