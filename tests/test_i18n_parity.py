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
