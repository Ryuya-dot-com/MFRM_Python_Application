"""Exercise output choices, inference holds, and current-result document caching."""

import pytest
from streamlit.testing.v1 import AppTest


@pytest.mark.parametrize("lang", ["en", "ja"])
@pytest.mark.parametrize("held", [False, True])
@pytest.mark.parametrize("gate_status", ["Ready", "Not ready", None])
def test_report_export_renders_only_the_selected_task(lang, held, gate_status):
    def view(lang, held, gate_status):
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app

        st.session_state["lang"] = lang
        st.session_state.setdefault("fit", {"analysis_id": "unchanged",
            "config": {"method": "MML" if held else "JMLE", "estimate_population_sd": held}})
        st.session_state["rendered"] = []
        gate = pd.DataFrame([{
            "GateArea": "Overall manuscript gate", "GateStatus": gate_status,
            "Evidence": "", "ManuscriptAction": "",
        }, {
            "GateArea": "Estimation", "GateStatus": "Ready",
            "Evidence": "Evidence retained", "ManuscriptAction": "Check the study context",
        }])
        if gate_status is None:
            gate = pd.DataFrame()
        def record(name):
            return lambda *args, **kwargs: st.session_state["rendered"].append(name)
        with patch.object(app, "build_publication_gate_summary", return_value=gate), \
             patch.object(app, "_render_publication_document_section", record("document")), \
             patch.object(app, "_render_manuscript_template_section", record("manuscript")), \
             patch.object(app, "_render_downloads", record("files")), \
             patch.object(app, "_render_guided_report_review_details", record("notes")), \
             patch.object(app, "show_report_section", record("advanced")):
            app._render_guided_report_export_section(
                st.session_state["fit"], {}, {}, pd.DataFrame(), pd.DataFrame(),
                None, None, generate_figures=False,
            )

    at = AppTest.from_function(view, args=(lang, held, gate_status)).run(timeout=45)
    assert not at.exception
    assert at.session_state["rendered"] == ["document"]
    assert len(at.warning) == int(held or gate_status != "Ready")
    assert not at.dataframe and not at.selectbox
    assert not at.subheader
    assert all(not item.value.startswith("guided.") for item in at.caption)
    assert all(not item.proto.expanded for item in at.expander)
    assert not any("{summary}" in item.value for item in at.markdown)
    for target in ("files", "review", "manuscript"):
        assert at.button(key=f"guided_report_next_{target}").label in at.button_group(key="guided_export_task").options
    at.button(key="guided_report_next_manuscript").click().run()
    assert not at.exception
    assert at.button_group(key="guided_export_task").value == "manuscript"
    assert at.session_state["rendered"] == ["manuscript"]
    assert len(at.warning) == int(held or gate_status != "Ready")
    at.button_group(key="guided_export_task").set_value("document").run()
    at.button(key="guided_report_next_files").click().run()
    assert not at.exception
    assert at.button_group(key="guided_export_task").value == "files"
    assert at.session_state["rendered"] == ["files"]
    at.button_group(key="guided_export_task").set_value("document").run()
    at.button(key="guided_report_next_review").click().run()
    assert not at.exception
    assert at.button_group(key="guided_export_task").value == "review"
    assert at.session_state["rendered"] == []
    if gate_status is not None:
        assert not at.expander[0].proto.expanded
        assert "Evidence retained" in [item.value for item in at.markdown]
        assert "Check the study context" in [item.value for item in at.markdown]
    else:
        assert not at.dataframe and at.info
    assert len(at.warning) == int(held or gate_status != "Ready")
    at.selectbox(key="guided_export_resource").select("work_notes").run()
    assert not at.exception
    assert at.session_state["rendered"] == ["notes"]
    at.selectbox(key="guided_export_resource").select("all_panels").run()
    assert not at.exception
    assert at.session_state["rendered"] == ["advanced"]
    assert at.session_state["fit"] == {"analysis_id": "unchanged",
        "config": {"method": "MML" if held else "JMLE", "estimate_population_sd": held}}


@pytest.mark.parametrize("previous_version", ["metric-definitions-v1", "report-faq-v1", "anchor-coverage-review-v1", "facet-equivalence-hold-v1", "uncertainty-output-hold-v1", "multiview-uncertainty-v1", "observed-writing-report-v1", "calibration-termination-hold-v1", "report-review-notes-v1", "report-first-read-qualification-v1", "report-connected-screening-v1"])
def test_document_formats_are_lazy_and_same_shape_changes_invalidate_cache(previous_version):
    def view(previous_version):
        import hashlib
        import pickle
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app

        st.session_state.setdefault("lang", "en")
        st.session_state.setdefault("builds", [])
        values = list(range(200))
        values[100] = st.session_state.get("middle_value", 100)
        diagnostics = {"measures": pd.DataFrame({"Estimate": values})}
        if not st.session_state.get("legacy_cache_seeded"):
            # A running session may still hold a report generated before the Q&A update.
            st.session_state["_publication_doc_pdf_fp"] = hashlib.sha256(pickle.dumps(
                (previous_version, {"config": {}}, diagnostics, None, "en", False), protocol=5,
            )).hexdigest()
            st.session_state["_publication_doc_pdf_bytes"] = b"old report without Q&A"
            st.session_state["legacy_cache_seeded"] = True
        def build(kind):
            def generate(*args, **kwargs):
                st.session_state["builds"].append(kind)
                if st.session_state.get("fail"):
                    raise ValueError("Export failed")
                return f"{kind}:{values[100]}:{kwargs['include_details']}".encode()
            return generate
        with patch.object(app, "build_publication_pdf_bytes", build("pdf")), \
             patch.object(app, "build_publication_word_bytes", build("word")), \
             patch.object(app, "build_publication_html_bytes", build("html")):
            app._render_publication_document_section({"config": {}}, diagnostics)

    at = AppTest.from_function(view, args=(previous_version,)).run(timeout=45)
    assert not at.exception
    assert at.session_state["builds"] == ["pdf"]
    assert len(at.get("download_button")) == 1
    assert not at.subheader
    assert all(not item.proto.expanded for item in at.expander)
    assert at.radio(key="publication_doc_format").label
    assert at.checkbox(key="publication_doc_details").label
    assert any("English" in item.value and "sharing" in item.value for item in at.caption)

    at.radio(key="publication_doc_format").set_value("word").run()
    assert at.session_state["builds"] == ["pdf", "word"]
    at.radio(key="publication_doc_format").set_value("pdf").run()
    assert at.session_state["builds"] == ["pdf", "word"]
    at.session_state["middle_value"] = 999
    at.run()
    assert at.session_state["builds"] == ["pdf", "word", "pdf"]
    assert at.session_state["_publication_doc_pdf_bytes"] == b"pdf:999:False"
    at.session_state["lang"] = "ja"
    at.run()
    assert at.session_state["builds"] == ["pdf", "word", "pdf", "pdf"]
    assert at.get("download_button")[0].proto.label == "PDFをダウンロード"
    assert not at.checkbox(key="publication_doc_details").value
    at.checkbox(key="publication_doc_details").check().run()
    assert at.session_state["_publication_doc_pdf_bytes"] == b"pdf:999:True"
    assert at.session_state["builds"] == ["pdf", "word", "pdf", "pdf", "pdf"]
    at.session_state["fail"] = True
    at.radio(key="publication_doc_format").set_value("html").run()
    assert not at.exception
    assert at.error and not at.get("download_button")
    at.session_state["fail"] = False
    at.run()
    assert not at.exception and not at.error
    assert len(at.get("download_button")) == 1


@pytest.mark.parametrize("lang", ["en", "ja"])
def test_single_csv_starts_with_summary_and_keeps_privacy_filter(lang):
    def view(lang):
        from collections import defaultdict
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app

        st.session_state["lang"] = lang
        frames = {
            "visualization_settings": pd.DataFrame({"Setting": ["default"]}),
            "summary": pd.DataFrame({"Converged": [True]}),
            "category_diagnostics": pd.DataFrame({"Category": [0, 1, 2], "Count": [2, 0, 3]}),
            "scorefile": pd.DataFrame({"Person": ["P1"], "Score": [2]}),
        }
        with patch.object(app, "collect_download_frames", return_value=(frames, defaultdict(pd.DataFrame))), \
             patch.object(app, "build_evidence_contract_text_assets", return_value={}):
            app._render_downloads({}, {}, {}, pd.DataFrame(), pd.DataFrame(), generate_figures=False)

    at = AppTest.from_function(view, args=(lang,)).run(timeout=45)
    assert not at.exception
    assert at.checkbox(key="downloads_public_export_mode").value
    assert at.selectbox(key="download_single_table").value == "summary"
    assert at.selectbox(key="download_single_table").options == (
        ["Analysis overview", "Category counts and fit"] if lang == "en" else
        ["分析の概要", "カテゴリーの使用度数と適合度"])
    assert "scorefile" not in at.selectbox(key="download_single_table").options
    csv_buttons = [b for b in at.get("download_button") if b.proto.id.endswith("-dl_selected_csv")]
    assert len(csv_buttons) == 1
    at.selectbox(key="download_single_table").select("category_diagnostics").run()
    assert not at.exception
    assert any("Count" in item.value and "Percent" in item.value for item in at.caption)
    assert any("mfrm_category_diagnostics.csv" in item.value and "3" in item.value for item in at.caption)
    at.checkbox(key="downloads_public_export_mode").uncheck().run()
    assert not at.exception
    # Changing privacy does not silently expose every technical table in the short list.
    assert "scorefile" not in at.selectbox(key="download_single_table").options
    at.checkbox(key="download_all_tables").check().run()
    assert "scorefile" in at.selectbox(key="download_single_table").options
    at.selectbox(key="download_single_table").select("scorefile").run()
    assert not at.exception
    at.checkbox(key="downloads_public_export_mode").check().run()
    assert not at.exception
    assert "scorefile" not in at.selectbox(key="download_single_table").options
    assert at.selectbox(key="download_single_table").value == "summary"
    at.selectbox(key="download_single_table").select("visualization_settings").run()
    at.checkbox(key="download_all_tables").uncheck().run()
    assert not at.exception
    assert at.selectbox(key="download_single_table").value == "summary"


@pytest.mark.parametrize("lang", ["en", "ja"])
def test_core_download_table_guidance_has_no_missing_translations(lang):
    import json
    from pathlib import Path
    import streamlit_app as app

    text = json.loads((Path(app.__file__).parent / "locales" / f"{lang}.json").read_text())
    for name in app._REPORT_TABLE_GUIDE_IDS:
        assert text["downloads"]["table_labels"][name].strip()
        assert text["downloads"]["table_reading"][name].strip()


@pytest.mark.parametrize("density,force_full", [("Essential", False), ("All panels", False), ("Essential", True)])
def test_advanced_report_renders_one_choice_and_restores_reading_position(density, force_full):
    def view(density, force_full):
        from contextlib import ExitStack
        from unittest.mock import patch
        import streamlit as st
        import streamlit_app as app
        st.radio("Language", ["en", "ja"], key="lang")
        st.session_state["app_view_density"] = density
        show = st.checkbox("Show report", value=True, key="show_report")
        st.session_state.setdefault("fit", {"analysis_id": "same-result", "config": {"method": "MML"}})
        st.session_state["rendered"] = []
        renderers = {
            "APA Report": "_render_apa_report", "Manuscript Template": "_render_manuscript_template_section",
            "Method Appendix": "_render_method_appendix_section", "Claim Guide": "_render_manuscript_claim_guide_section",
            "Tables": "_render_report_tables", "Reporting Checklist": "_render_reporting_checklist",
            "Additional analyses": "_render_additional_report_analyses",
            "Facet Equivalence": "_render_facet_equivalence", "Readiness": "_render_final_readiness_section",
            "Publication Document": "_render_publication_document_section",
        }
        def record(panel):
            return lambda *a, **kw: st.session_state["rendered"].append(panel)
        with ExitStack() as stack:
            stack.enter_context(patch.object(app, "mfrm_estimate", side_effect=AssertionError("Display must not refit")))
            for panel, renderer in renderers.items():
                stack.enter_context(patch.object(app, renderer, record(panel)))
            if show:
                app.show_report_section(st.session_state["fit"], {}, force_full=force_full)

    at = AppTest.from_function(view, args=(density, force_full)).run(timeout=30)
    assert not at.exception and not at.tabs
    assert at.session_state["rendered"] == ["APA Report"]
    panels = ["APA Report", "Claim Guide", "Tables", "Reporting Checklist", "Readiness", "Publication Document", "Additional analyses"]
    if force_full or density != "Essential":
        panels += ["Manuscript Template", "Method Appendix", "Facet Equivalence"]
    assert len(at.selectbox(key="report_detail_panel").options) == len(panels)
    for panel in panels:
        at.selectbox(key="report_detail_panel").set_value(panel).run()
        assert not at.exception
        assert at.session_state["rendered"] == [panel]
    at.selectbox(key="report_detail_panel").set_value("Tables").run()
    at.radio(key="lang").set_value("ja").run()
    assert not at.exception and at.session_state["rendered"] == ["Tables"]
    assert at.selectbox(key="report_detail_panel").value == "Tables"
    assert at.selectbox(key="report_detail_panel").label == "確認する内容"
    at.checkbox(key="show_report").uncheck().run()
    assert not at.selectbox and at.session_state["rendered"] == []
    at.checkbox(key="show_report").check().run()
    assert not at.exception and at.session_state["rendered"] == ["Tables"]
    assert at.session_state["fit"] == {"analysis_id": "same-result", "config": {"method": "MML"}}


def test_advanced_report_unavailable_selection_falls_back_in_essential_mode():
    def view():
        from contextlib import ExitStack
        from unittest.mock import patch
        import streamlit as st
        import streamlit_app as app
        st.radio("Density", ["All panels", "Essential"], key="app_view_density")
        with ExitStack() as stack:
            for name in ("_render_apa_report", "_render_manuscript_template_section", "_render_method_appendix_section",
                         "_render_manuscript_claim_guide_section", "_render_report_tables", "_render_reporting_checklist",
                         "_render_facet_equivalence", "_render_final_readiness_section", "_render_publication_document_section"):
                stack.enter_context(patch.object(app, name))
            app.show_report_section({}, {})
    at = AppTest.from_function(view).run(timeout=30)
    at.selectbox(key="report_detail_panel").set_value("Method Appendix").run()
    at.radio(key="app_view_density").set_value("Essential").run()
    assert not at.exception
    assert at.selectbox(key="report_detail_panel").value == "APA Report"
    assert at.session_state["report_detail_selection"] == "APA Report"
