"""Reading guidance survives each export and remains localized beside downloads."""
from io import BytesIO
import re
import xml.etree.ElementTree as ET

from docx import Document
import pandas as pd
import pytest
from reportlab.platypus import KeepTogether, Paragraph, SimpleDocTemplate
from streamlit.testing.v1 import AppTest

import streamlit_app as app


@pytest.mark.parametrize("include_details", [False, True])
@pytest.mark.parametrize("format", ["word", "pdf", "html"])
def test_all_reports_end_with_english_faq_and_its_references(format, include_details, monkeypatch):
    result = {"config": {"model": "RSM", "method": "MML", "estimate_population_sd": True}}
    diagnostics = {"measures": pd.DataFrame()}
    monkeypatch.setitem(app.st.session_state, "lang", "ja")
    monkeypatch.setattr(app, "_publication_figure_payloads", lambda *args: [])
    monkeypatch.setattr(app, "generate_method_appendix_text", lambda *args: "Method example (Andrich, 1978).")
    monkeypatch.setattr(app, "generate_manuscript_reporting_template", lambda *args: "Writing example.")
    # Inspect the actual PDF story while still building the bytes; no PDF reader dependency.
    pdf_text = []
    original_build = SimpleDocTemplate.build
    def capture(self, story, *args, **kwargs):
        for part in story:
            if isinstance(part, KeepTogether):
                pdf_text.extend(child.getPlainText() for child in part._content if isinstance(child, Paragraph))
            elif isinstance(part, Paragraph):
                pdf_text.append(part.getPlainText())
        return original_build(self, story, *args, **kwargs)
    monkeypatch.setattr(SimpleDocTemplate, "build", capture)
    data = getattr(app, f"build_publication_{format}_bytes")(
        result, diagnostics, include_details=include_details)
    if format == "word":
        text = "\n".join(p.text for p in Document(BytesIO(data)).paragraphs)
    elif format == "pdf":
        assert data.startswith(b"%PDF")
        text = "\n".join(pdf_text)
    else:
        root = ET.fromstring(data.decode().split("<!DOCTYPE html>\n", 1)[1])
        text = "\n".join(root.itertext())
        assert len(root.findall(".//section/h3")) == 4
    locale = app._load_locale("en")["report_faq"]
    for topic in ("overfit", "pcm_rsm", "wright", "bias"):
        assert text.count(locale[topic + "_question"]) == 1
    assert text.index("Results tables") < text.index("Questions and Answers") < text.index("References")
    assert ("Writing example." in text) == include_details
    assert "Interpretation is withheld for this calibration" in text
    assert "Bias/no-bias conclusions remain withheld" in text
    assert "PCA" not in text.split("Questions and Answers", 1)[1].split("References", 1)[0]
    assert "英語スピーチ" not in text
    assert "**" not in text
    for reference in app.build_apa_reference_list(app._publication_report_faq()):
        assert text.count(reference) == 1
        assert text.index(reference) > text.index("References")
    assert result["config"]["estimate_population_sd"] is True


def test_each_faq_citation_resolves_including_combined_citations():
    text = app._publication_report_faq()
    citations = [citation.strip() for group in re.findall(r"\(([^()]+)\)", text)
                 for citation in group.split(";") if re.search(r"\d{4}|n\.d\.", citation)]
    keys = {app._CITATION_TO_KEY[f"({citation})"] for citation in citations}
    refs = app.build_apa_reference_list(text + "\n" + text)
    assert len(keys) == 7
    assert set(refs) == {app._APA_REFERENCE_LIBRARY[key] for key in keys}
    for key in keys:
        if "(n.d." in app._APA_REFERENCE_LIBRARY[key]:
            reference = app._APA_REFERENCE_LIBRARY[key]
            # Online manuals retain source fields without turning retrieval dates into publication years.
            record = app._parse_apa_entry(key, reference)
            assert record["type"] == "manual" and record["year"] == ""
            assert record["url"] in reference and record["note"] == "Retrieved October 16, 2026"
            bib = app.apa_entry_to_bibtex(key, reference)
            ris = app.apa_entry_to_ris(key, reference)
            assert bib.startswith("@manual") and "year =" not in bib
            assert ris.startswith("TY  - ELEC") and "PY  -" not in ris
            assert record["url"] in bib and record["url"] in ris
            assert record["title"].replace("&", r"\&") in bib and record["title"] in ris


@pytest.mark.parametrize("lang", ["en", "ja"])
def test_faq_is_collapsed_localized_and_keeps_current_fit(lang):
    def view(lang):
        from unittest.mock import patch
        import streamlit as st
        import streamlit_app as app
        st.session_state["lang"] = lang
        st.session_state.setdefault("fit", {"analysis_id": "unchanged",
            "config": {"method": "MML", "estimate_population_sd": True}})
        st.session_state.setdefault("builds", 0)
        def build(*args, **kwargs):
            st.session_state["builds"] += 1
            return b"PDF"
        with patch.object(app, "mfrm_estimate", side_effect=AssertionError("Unexpected refit")), \
             patch.object(app, "build_publication_pdf_bytes", build):
            app._render_publication_document_section(st.session_state["fit"], {})
    at = AppTest.from_function(view, args=(lang,)).run(timeout=45)
    assert not at.exception
    locale = app._load_locale(lang)["report_faq"]
    panel = next(e for e in at.expander if e.label == locale["title"])
    assert not panel.proto.expanded
    assert panel.markdown[0].value == app._publication_report_faq(lang)
    assert panel.caption[0].value == locale["references_caption"]
    assert len(panel.caption) == 8  # Intro plus seven distinct sources, including combined citations.
    assert len(at.get("download_button")) == 1
    at.run()
    assert not at.exception and at.session_state["builds"] == 1
    assert at.session_state["fit"] == {"analysis_id": "unchanged",
        "config": {"method": "MML", "estimate_population_sd": True}}
