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
    assert all(not item.value.startswith("guided.") for item in at.caption)
    at.button_group(key="guided_export_task").set_value("files").run()
    assert not at.exception
    assert at.session_state["rendered"] == ["files"]
    at.button_group(key="guided_export_task").set_value("review").run()
    assert not at.exception
    assert at.session_state["rendered"] == []
    if gate_status is not None:
        assert at.expander[0].proto.expanded
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


def test_document_formats_are_lazy_and_same_shape_changes_invalidate_cache():
    def view():
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app

        st.session_state.setdefault("lang", "en")
        st.session_state.setdefault("builds", [])
        values = list(range(200))
        values[100] = st.session_state.get("middle_value", 100)
        diagnostics = {"measures": pd.DataFrame({"Estimate": values})}
        def build(kind):
            def generate(*args):
                st.session_state["builds"].append(kind)
                if st.session_state.get("fail"):
                    raise ValueError("Export failed")
                return f"{kind}:{values[100]}".encode()
            return generate
        with patch.object(app, "build_publication_pdf_bytes", build("pdf")), \
             patch.object(app, "build_publication_word_bytes", build("word")), \
             patch.object(app, "build_publication_html_bytes", build("html")):
            app._render_publication_document_section({"config": {}}, diagnostics)

    at = AppTest.from_function(view).run(timeout=45)
    assert not at.exception
    assert at.session_state["builds"] == ["pdf"]
    assert len(at.get("download_button")) == 1

    at.radio(key="publication_doc_format").set_value("word").run()
    assert at.session_state["builds"] == ["pdf", "word"]
    at.radio(key="publication_doc_format").set_value("pdf").run()
    assert at.session_state["builds"] == ["pdf", "word"]
    at.session_state["middle_value"] = 999
    at.run()
    assert at.session_state["builds"] == ["pdf", "word", "pdf"]
    assert at.session_state["_publication_doc_pdf_bytes"] == b"pdf:999"
    at.session_state["lang"] = "ja"
    at.run()
    assert at.session_state["builds"] == ["pdf", "word", "pdf", "pdf"]
    assert at.get("download_button")[0].proto.label == "PDFをダウンロード"
    at.session_state["fail"] = True
    at.radio(key="publication_doc_format").set_value("html").run()
    assert not at.exception
    assert at.error and not at.get("download_button")
    at.session_state["fail"] = False
    at.run()
    assert not at.exception and not at.error
    assert len(at.get("download_button")) == 1


def test_single_csv_starts_with_summary_and_keeps_privacy_filter():
    def view():
        from collections import defaultdict
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app

        st.session_state["lang"] = "en"
        frames = {
            "visualization_settings": pd.DataFrame({"Setting": ["default"]}),
            "summary": pd.DataFrame({"Converged": [True]}),
            "scorefile": pd.DataFrame({"Person": ["P1"], "Score": [2]}),
        }
        with patch.object(app, "collect_download_frames", return_value=(frames, defaultdict(pd.DataFrame))), \
             patch.object(app, "build_evidence_contract_text_assets", return_value={}):
            app._render_downloads({}, {}, {}, pd.DataFrame(), pd.DataFrame(), generate_figures=False)

    at = AppTest.from_function(view).run(timeout=45)
    assert not at.exception
    assert at.checkbox(key="downloads_public_export_mode").value
    assert at.selectbox(key="download_single_table").value == "summary"
    assert "scorefile" not in at.selectbox(key="download_single_table").options
    csv_buttons = [b for b in at.get("download_button") if b.proto.label == "Download selected table (CSV)"]
    assert len(csv_buttons) == 1
    at.checkbox(key="downloads_public_export_mode").uncheck().run()
    assert not at.exception
    assert "scorefile" in at.selectbox(key="download_single_table").options
    at.selectbox(key="download_single_table").select("scorefile").run()
    assert not at.exception
    at.checkbox(key="downloads_public_export_mode").check().run()
    assert not at.exception
    assert "scorefile" not in at.selectbox(key="download_single_table").options
