"""Tests for the BibTeX / RIS export of the cited APA reference library.

The Streamlit app stores its bibliography as APA prose strings in
``_APA_REFERENCE_LIBRARY``. ``apa_entry_to_bibtex`` / ``apa_entry_to_ris``
parse each entry on the fly into a structured record, then serialise
it to the requested format. ``build_bibtex_from_cited`` /
``build_ris_from_cited`` walk the same ``_CITATION_TO_KEY`` map that
``build_apa_reference_list`` uses, so the three export formats stay in
sync for a given run.

The contract pinned here covers:

* Every entry in the live ``_APA_REFERENCE_LIBRARY`` round-trips through
  the parser without falling back to ``@misc`` — that's the regression
  net that catches a future APA prose entry whose format the parser
  cannot classify.
* The output for canonical entries matches expected BibTeX / RIS fields
  (author / year / title / journal / volume / number / pages / doi),
  including:
    - Single-author journal article (Andrich 1978)
    - Two-author journal article with hyphenated initials (Eckes & Jin 2021)
    - Three-author journal article (Lai, Wolfe & Vickers 2015)
    - Book with edition and DOI (Engelhard 2013)
    - Book chapter / proceedings (Hagberg et al. 2008)
    - Journal article with question-mark title (Linacre 2002)
    - Book with Winsteps.com publisher (Linacre 2007)
* Author conversion: ``Lastname, Initials`` pairs are recombined and
  joined with ``" and "`` (BibTeX) / one ``AU  -`` line each (RIS).
* Trailing periods are stripped from titles by BibTeX convention but
  ``?`` and ``!`` are preserved.
* Page ranges use the en-dash ``--`` for BibTeX (LaTeX convention).
* ``build_bibtex_from_cited`` / ``build_ris_from_cited`` emit the same
  set of references as ``build_apa_reference_list`` for the same
  ``text`` + ``always_include``.
* Empty input returns the empty string (no header, no entries).
"""

from __future__ import annotations

import pytest

import streamlit_app as app


# -----------------------------------------------------------------------------
# Coverage: every entry in the live library parses to a non-misc record
# -----------------------------------------------------------------------------


def test_every_apa_library_entry_parses_to_non_misc():
    """The parser must classify every shipped APA entry. A future entry
    whose format breaks the parser will trip this and prompt either a
    library-side rewording or a parser extension."""
    failures = []
    for key, apa in app._APA_REFERENCE_LIBRARY.items():
        bib = app.apa_entry_to_bibtex(key, apa)
        if bib.startswith("@misc"):
            failures.append(key)
    assert not failures, (
        f"{len(failures)} APA entries fell back to @misc: {failures}. "
        f"Either rephrase the APA prose or extend _parse_apa_entry."
    )


# -----------------------------------------------------------------------------
# Author conversion
# -----------------------------------------------------------------------------


def test_authors_single_author_preserves_initials_period():
    assert app._apa_authors_to_bibtex("Andrich, D.") == "Andrich, D."


def test_authors_two_authors_join_with_and():
    out = app._apa_authors_to_bibtex("Eckes, T., & Jin, K.-Y.")
    assert out == "Eckes, T. and Jin, K.-Y."


def test_authors_three_authors_join_with_and():
    out = app._apa_authors_to_bibtex("Lai, E. R., Wolfe, E. W., & Vickers, D.")
    assert out == "Lai, E. R. and Wolfe, E. W. and Vickers, D."


def test_authors_two_initial_pair():
    """Authors with two-initial given names (A. A., D. A.) must stay
    glued to the lastname, not split into a separate pseudo-author."""
    out = app._apa_authors_to_bibtex("Hagberg, A. A., Schult, D. A., & Swart, P. J.")
    assert out == "Hagberg, A. A. and Schult, D. A. and Swart, P. J."


def test_authors_jr_suffix_glued_to_previous():
    out = app._apa_authors_to_bibtex("Engelhard, G., Jr.")
    assert out == "Engelhard, G., Jr."


def test_authors_four_authors_join_with_and():
    out = app._apa_authors_to_bibtex(
        "Cronbach, L. J., Gleser, G. C., Nanda, H., & Rajaratnam, N."
    )
    assert out == (
        "Cronbach, L. J. and Gleser, G. C. and Nanda, H. and Rajaratnam, N."
    )


def test_authors_empty_input_returns_empty():
    assert app._apa_authors_to_bibtex("") == ""
    assert app._apa_authors_to_bibtex("   ") == ""


# -----------------------------------------------------------------------------
# Canonical BibTeX outputs
# -----------------------------------------------------------------------------


def test_bibtex_journal_article_with_doi():
    bib = app.apa_entry_to_bibtex(
        "Andrich_1978", app._APA_REFERENCE_LIBRARY["Andrich_1978"]
    )
    assert bib.startswith("@article{Andrich_1978,")
    assert "author = {Andrich, D.}" in bib
    assert "title = {A rating formulation for ordered response categories}" in bib
    assert "journal = {Psychometrika}" in bib
    assert "year = {1978}" in bib
    assert "volume = {43}" in bib
    assert "number = {4}" in bib
    assert "pages = {561--573}" in bib  # en-dash → BibTeX `--`
    assert "doi = {10.1007/BF02293814}" in bib


def test_bibtex_book_with_doi_and_subtitle():
    bib = app.apa_entry_to_bibtex(
        "Engelhard_2013", app._APA_REFERENCE_LIBRARY["Engelhard_2013"]
    )
    assert bib.startswith("@book{Engelhard_2013,")
    assert "author = {Engelhard, G., Jr.}" in bib
    assert "publisher = {Routledge}" in bib
    assert "doi = {10.4324/9780203073636}" in bib


def test_bibtex_book_winsteps_com_publisher():
    """The book publisher regex must accept ``Winsteps.com`` — periods
    inside the publisher value were a v0.2.14-beta regression."""
    bib = app.apa_entry_to_bibtex(
        "Linacre_2007", app._APA_REFERENCE_LIBRARY["Linacre_2007"]
    )
    assert "publisher = {Winsteps.com}" in bib


def test_bibtex_chapter_with_no_publisher():
    """Proceedings entries (Hagberg 2008) have no publisher after the
    page range — chapter regex must allow that."""
    bib = app.apa_entry_to_bibtex(
        "Hagberg_Schult_Swart_2008",
        app._APA_REFERENCE_LIBRARY["Hagberg_Schult_Swart_2008"],
    )
    assert bib.startswith("@incollection{Hagberg_Schult_Swart_2008,")
    assert "editor = {G. Varoquaux and T. Vaught and J. Millman}" in bib
    assert "booktitle = {Proceedings of the 7th Python in Science Conference}" in bib
    assert "pages = {11--15}" in bib


def test_bibtex_journal_article_question_mark_title():
    """Title ending in ``?`` (Linacre 2002) must keep the ``?`` and
    still parse as @article — the article terminator allows ``[.?!]``."""
    bib = app.apa_entry_to_bibtex(
        "Linacre_2002", app._APA_REFERENCE_LIBRARY["Linacre_2002"]
    )
    assert bib.startswith("@article{Linacre_2002,")
    assert "title = {What do infit and outfit, mean-square and standardized mean?}" in bib


# -----------------------------------------------------------------------------
# RIS output
# -----------------------------------------------------------------------------


def test_ris_journal_article_emits_one_author_per_line():
    ris = app.apa_entry_to_ris(
        "Lai_Wolfe_Vickers_2015",
        app._APA_REFERENCE_LIBRARY["Lai_Wolfe_Vickers_2015"],
    )
    assert ris.startswith("TY  - JOUR")
    assert ris.endswith("ER  - ")
    assert "AU  - Lai, E. R." in ris
    assert "AU  - Wolfe, E. W." in ris
    assert "AU  - Vickers, D." in ris
    assert "JO  - Educational and Psychological Measurement" in ris
    assert "SP  - 102" in ris
    assert "EP  - 125" in ris
    assert "DO  - 10.1177/0013164414530143" in ris


def test_ris_book_emits_pb_tag():
    ris = app.apa_entry_to_ris(
        "Linacre_1989", app._APA_REFERENCE_LIBRARY["Linacre_1989"]
    )
    assert "TY  - BOOK" in ris
    assert "PB  - MESA Press" in ris


def test_ris_chapter_emits_t2_and_a2_tags():
    ris = app.apa_entry_to_ris(
        "Hagberg_Schult_Swart_2008",
        app._APA_REFERENCE_LIBRARY["Hagberg_Schult_Swart_2008"],
    )
    assert "TY  - CHAP" in ris
    assert "T2  - Proceedings of the 7th Python in Science Conference" in ris
    assert "A2  - G. Varoquaux" in ris


# -----------------------------------------------------------------------------
# Bundle export: build_bibtex_from_cited / build_ris_from_cited
# -----------------------------------------------------------------------------


def test_build_bibtex_from_cited_includes_only_cited_references():
    text = "Estimation followed Andrich (1978) (Andrich, 1978)."
    bib = app.build_bibtex_from_cited(text)
    assert "@article{Andrich_1978," in bib
    # An uncited reference must not appear.
    assert "Eckes_Jin_2021" not in bib


def test_build_bibtex_from_cited_always_include_added():
    text = "no citation here"
    bib = app.build_bibtex_from_cited(text, always_include=["Andrich_1978"])
    assert "@article{Andrich_1978," in bib


def test_build_bibtex_from_cited_empty_returns_empty_string():
    """When neither the text nor always_include surface any reference,
    the BibTeX bundle is empty so the download button can be hidden."""
    assert app.build_bibtex_from_cited("no references at all") == ""


def test_build_bibtex_carries_header_with_app_version():
    bib = app.build_bibtex_from_cited(
        "Used the bootstrap (Efron & Tibshirani, 1993)."
    )
    assert bib.startswith("% BibTeX export")
    assert app.APP_VERSION in bib


def test_build_ris_from_cited_includes_only_cited_references():
    text = "Reliability was high (Andrich, 1978)."
    ris = app.build_ris_from_cited(text)
    assert "TY  - JOUR" in ris
    assert "AU  - Andrich, D." in ris


def test_build_ris_from_cited_empty_returns_empty_string():
    assert app.build_ris_from_cited("no references at all") == ""


# -----------------------------------------------------------------------------
# Parity: BibTeX bundle covers same keys as the APA reference list
# -----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "text, always_include",
    [
        ("The bootstrap (Efron & Tibshirani, 1993) underpins the CI.", None),
        ("RSM (Andrich, 1978) and PCM (Masters, 1982).", ["Linacre_1989"]),
    ],
)
def test_bibtex_and_apa_bundles_cover_same_keys(text, always_include):
    """Both bundles must include the same references — otherwise the
    user's APA reference list and BibTeX bundle drift."""
    apa_list = app.build_apa_reference_list(text, always_include=always_include)
    bib = app.build_bibtex_from_cited(text, always_include=always_include)
    # Count of @-entries in the BibTeX bundle.
    bib_count = sum(1 for line in bib.splitlines() if line.startswith("@"))
    assert bib_count == len(apa_list)
