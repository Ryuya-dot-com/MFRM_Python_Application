"""Integrity tests for the APA 7 reference library + citation scanner.

The Publication Document, sample-scenario About box, and help tabs
all lean on `_APA_REFERENCE_LIBRARY` + `_CITATION_TO_KEY`. If a
citation string in the narrative does not resolve to a library entry
it silently vanishes from the bibliography — which is exactly the
kind of silent drift pytest should catch.
"""

from __future__ import annotations

import re

import streamlit_app as app


def test_reference_library_is_non_empty():
    assert len(app._APA_REFERENCE_LIBRARY) >= 30


def test_citation_map_is_non_empty():
    assert len(app._CITATION_TO_KEY) >= 30


def test_every_citation_resolves_to_a_library_entry():
    missing: list[str] = []
    for citation, key in app._CITATION_TO_KEY.items():
        if key not in app._APA_REFERENCE_LIBRARY:
            missing.append(f"{citation!r} → {key!r}")
    assert not missing, f"unresolved citation mappings: {missing}"


def test_citation_keys_follow_lastname_year_convention():
    pattern = re.compile(r"^[A-Za-z][A-Za-z_-]*_\d{4}[a-z]?$")
    for key in app._APA_REFERENCE_LIBRARY:
        assert pattern.match(key), (
            f"reference key {key!r} does not follow `Lastname_YYYY` convention"
        )


def test_citation_parens_use_correct_form():
    """Every citation string should end with `, YYYY)` and start with `(Surname`.

    Accepts 1-, 2-, and 3-author forms:
      (Linacre, 2024)
      (Wright & Masters, 1982)
      (Bradlow, Wainer & Wang, 1999)
    """
    pattern = re.compile(r"^\([A-Z][^)]+, \d{4}\)$")
    bad = [c for c in app._CITATION_TO_KEY if not pattern.match(c)]
    assert not bad, f"malformed citation strings: {bad}"


def test_reference_entries_end_without_trailing_whitespace():
    for key, entry in app._APA_REFERENCE_LIBRARY.items():
        assert entry == entry.strip(), f"{key!r} has leading/trailing whitespace"


def test_reference_entries_include_year():
    """Every APA entry must contain a 4-digit year somewhere in the body."""
    year_pattern = re.compile(r"\b(19|20)\d{2}\b")
    missing_year = [k for k, v in app._APA_REFERENCE_LIBRARY.items() if not year_pattern.search(v)]
    assert not missing_year, f"entries without a 4-digit year: {missing_year}"


def test_collect_cited_references_finds_known_citation():
    text = "As shown by (Linacre, 2002), the fit statistics are robust."
    refs = app.collect_cited_references(text)
    assert any("Linacre" in r and "2002" in r for r in refs)


def test_collect_cited_references_deduplicates():
    text = "See (Linacre, 2002). Also (Linacre, 2002). And yet again (Linacre, 2002)."
    refs = app.collect_cited_references(text)
    linacre_hits = [r for r in refs if "Linacre" in r and "2002" in r]
    assert len(linacre_hits) == 1


def test_collect_cited_references_returns_sorted():
    text = (
        "See (Wright & Masters, 1982) and (Andrich, 1978). "
        "Also (Linacre, 2002)."
    )
    refs = app.collect_cited_references(text)
    assert refs == sorted(refs)


def test_collect_cited_references_ignores_unknown_citations():
    text = "Fake reference (Nobody, 9999) should be dropped."
    refs = app.collect_cited_references(text)
    assert not any("Nobody" in r for r in refs)


def test_sample_scenario_citations_all_resolve():
    """Every citation in every built-in scenario must map to a library entry."""
    for key, scenario in app.SAMPLE_DATA_SCENARIOS.items():
        for cite in scenario["citations"]:
            mapped = app._CITATION_TO_KEY.get(cite)
            assert mapped is not None, f"{key!r}: {cite!r} not in map"
            assert mapped in app._APA_REFERENCE_LIBRARY, (
                f"{key!r}: {cite!r} → {mapped!r} not in library"
            )
