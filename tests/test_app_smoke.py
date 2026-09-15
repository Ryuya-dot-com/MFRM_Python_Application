from pathlib import Path
import pytest
from streamlit.testing.v1 import AppTest


def test_app_initial_render():
    at = AppTest.from_file(Path(__file__).resolve().parents[1] / "streamlit_app.py").run(timeout=30)
    assert not at.exception
    assert any("MFRM | Rating analysis" in title.value for title in at.title)
    about = next(e for e in at.expander if e.label == "About this beta and its limitations")
    assert not about.proto.expanded
    assert any("source commit:" in caption.value for caption in about.caption)
    assert not any("Data privacy" in warning.value for warning in at.warning)
    assert not any("Keyboard shortcuts" in str(expander.label) for expander in at.expander)
    from mfrm_app import citation
    software = next(e for e in at.expander if e.label == "License & citation")
    assert not software.proto.expanded
    assert citation.APA in [item.value for item in software.code]
    assert {b.proto.label for b in software.get("download_button")} == {"BibTeX (.bib)", "CITATION.cff"}
    assert any("MIT License" in item.value for item in software.markdown)


def test_software_citation_matches_release_and_license():
    import streamlit_app as app
    from mfrm_app import citation

    assert citation.METADATA["version"] == app.APP_VERSION
    assert citation.METADATA["license"] == "MIT"
    assert (Path(app.__file__).parent / "LICENSE").read_text().startswith("MIT License\n")
    assert citation.METADATA["type"] == "software"
    assert citation.METADATA["version"] in citation.APA
    assert citation.METADATA["version"] in citation.BIBTEX
    assert "doi" not in citation.METADATA


def test_user_data_source_promotes_privacy_warning():
    at = AppTest.from_file(Path(__file__).resolve().parents[1] / "streamlit_app.py").run(timeout=30)
    at.button(key="onboarding_dismiss").click()
    at.run(timeout=30)

    assert not at.exception
    assert at.session_state["data_source_flat"] == "paste"
    assert any("pasted and uploaded rating files" in warning.value for warning in at.warning)


def test_two_level_source_picker_keeps_sample_detail_visible_and_stable() -> None:
    at = AppTest.from_file(Path(__file__).resolve().parents[1] / "streamlit_app.py").run(timeout=30)
    at.button(key="onboarding_skip_guide").click()
    at.run(timeout=30)

    assert at.selectbox(key="data_source_class").value == "sample"
    assert at.selectbox(key="data_source_scenario").value == "writing_essay"
    at.selectbox(key="data_source_scenario").set_value("clinical_osce")
    at.run(timeout=30)
    assert at.session_state["data_source_flat"] == "scenario:clinical_osce"

    at.selectbox(key="data_source_class").set_value("paste")
    at.run(timeout=30)
    assert at.session_state["data_source_flat"] == "paste"
    assert not any(item.key == "data_source_scenario" for item in at.selectbox)
    assert any("pasted and uploaded rating files" in warning.value for warning in at.warning)

    at.selectbox(key="data_source_class").set_value("sample")
    at.run(timeout=30)
    assert at.selectbox(key="data_source_scenario").value == "clinical_osce"
    assert at.session_state["data_source_flat"] == "scenario:clinical_osce"


def test_legacy_flat_source_state_migrates_into_two_level_picker() -> None:
    at = AppTest.from_file(Path(__file__).resolve().parents[1] / "streamlit_app.py").run(timeout=30)
    at.button(key="onboarding_skip_guide").click()
    at.run(timeout=30)

    at.session_state["data_source_flat"] = "upload"
    at.run(timeout=30)
    assert at.selectbox(key="data_source_class").value == "upload"
    assert at.session_state["data_source_flat"] == "upload"


def test_guided_setup_only_renders_first_run_decisions() -> None:
    at = AppTest.from_file(Path(__file__).resolve().parents[1] / "streamlit_app.py").run(timeout=30)
    at.button(key="onboarding_skip_guide").click()
    at.run(timeout=30)

    assert not at.exception
    assert at.radio(key="facets_mode_workflow_mode").value == "Guided defaults"
    assert at.selectbox(key="facets_mode_analysis_depth").options == [
        "Fast preview",
        "Standard (recommended)",
        "Full publication",
    ]
    selectbox_keys = {item.key for item in at.selectbox}
    checkbox_keys = {item.key for item in at.checkbox}
    number_input_keys = {item.key for item in at.number_input}
    radio_keys = {item.key for item in at.radio}
    assert "viz_ci_level" not in selectbox_keys
    assert "facet_regularization_mode" not in selectbox_keys
    assert not any(
        key and str(key).startswith("facets_mode_weight_col_")
        for key in selectbox_keys
    )
    assert "facets_mode_population_enabled" not in checkbox_keys
    assert "facets_mode_totalscore" not in checkbox_keys
    assert "facets_mode_maxit" not in number_input_keys
    assert "facets_mode_anchor_policy" not in radio_keys
    assert "Anchor constraints" not in {item.value for item in at.subheader}
    assert {"Analysis choices", "Analysis coverage"}.issubset(
        {item.value for item in at.subheader}
    )


def test_advanced_setup_restores_full_technical_controls() -> None:
    at = AppTest.from_file(Path(__file__).resolve().parents[1] / "streamlit_app.py").run(timeout=30)
    at.button(key="onboarding_skip_guide").click()
    at.run(timeout=30)
    at.radio(key="facets_mode_workflow_mode").set_value("Advanced controls")
    at.run(timeout=30)

    assert not at.exception
    assert "Custom" in at.selectbox(key="facets_mode_analysis_depth").options
    assert at.selectbox(key="viz_ci_level").value == 0.95
    assert at.selectbox(key="facet_regularization_mode").value == (
        "Off: unpenalized JMLE/MML"
    )
    assert at.checkbox(key="facets_mode_population_enabled").value is False
    assert at.number_input(key="facets_mode_maxit").value == 400
    assert at.radio(key="facets_mode_anchor_policy").value == "warn"
    assert at.checkbox(key="facets_mode_totalscore").value is True
    assert any(
        item.key and str(item.key).startswith("facets_mode_weight_col_")
        for item in at.selectbox
    )


def test_returning_to_guided_resets_hidden_visual_and_compute_customization() -> None:
    at = AppTest.from_file(Path(__file__).resolve().parents[1] / "streamlit_app.py").run(timeout=30)
    at.button(key="onboarding_skip_guide").click()
    at.run(timeout=30)
    at.radio(key="facets_mode_workflow_mode").set_value("Advanced controls")
    at.run(timeout=30)
    at.selectbox(key="viz_ci_level").set_value(0.66)
    at.selectbox(key="facets_mode_analysis_depth").set_value("Custom")
    at.run(timeout=30)

    at.radio(key="facets_mode_workflow_mode").set_value("Guided defaults")
    at.run(timeout=30)
    assert not at.exception
    assert at.session_state["viz_ci_level"] == 0.95
    assert at.selectbox(key="facets_mode_analysis_depth").value == (
        "Standard (recommended)"
    )
    assert not any(item.key == "viz_ci_level" for item in at.selectbox)
    assert "Custom" not in at.selectbox(key="facets_mode_analysis_depth").options


def test_compact_guided_setup_renders_in_japanese() -> None:
    at = AppTest.from_file(Path(__file__).resolve().parents[1] / "streamlit_app.py").run(timeout=30)
    at.button(key="onboarding_skip_guide").click()
    at.run(timeout=30)
    at.selectbox(key="facets_mode_analysis_depth").set_value("Full publication")
    at.run(timeout=30)
    at.radio(key="lang").set_value("ja")
    at.run(timeout=30)

    assert not at.exception
    subheaders = {item.value for item in at.subheader}
    assert {"モデルと推定方法", "分析内容の設定"}.issubset(subheaders)
    assert "Anchor 制約" not in subheaders
    assert any(
        "標準設定を使用しています" in str(item.value)
        for item in at.caption
    )

    source = at.selectbox(key="data_source_class")
    depth = at.selectbox(key="facets_mode_analysis_depth")
    assert source.value == "sample"
    assert depth.value == "Full publication"
    assert source.proto.set_value and source.proto.raw_value == "サンプルデータ"
    assert depth.proto.set_value and depth.proto.raw_value == "図・出力を含む（Full publication）"
    assert not any("Session State API" in warning.value for warning in at.warning)



def test_detected_weight_is_explained_and_not_offered_as_a_guided_facet() -> None:
    at = AppTest.from_string(
        """
import pandas as pd
import streamlit_app as app

data = pd.DataFrame({
    "Person": ["P1", "P1", "P2", "P2"],
    "Score": [1, 2, 2, 1],
    "Rater": ["R1", "R2", "R1", "R2"],
    "Task": ["T1", "T2", "T1", "T2"],
    "Weight": [1.0, 2.0, 1.0, 2.0],
})
app.run_facets_mode(app.load_core_namespace(), data)
""",
        default_timeout=30,
    ).run()

    assert not at.exception
    facet_picker = next(
        item for item in at.multiselect
        if item.key and str(item.key).startswith("facets_mode_facet_cols_")
    )
    assert "Weight" not in facet_picker.options
    assert any(
        "Guided defaults use equal weights" in str(item.value)
        for item in at.caption
    )
    assert not any(
        item.key and str(item.key).startswith("facets_mode_weight_col_")
        for item in at.selectbox
    )

    at.radio(key="facets_mode_workflow_mode").set_value("Advanced controls")
    at.run()
    assert not at.exception
    weight_picker = next(
        item for item in at.selectbox
        if item.key and str(item.key).startswith("facets_mode_weight_col_")
    )
    assert weight_picker.value == "Weight"


@pytest.mark.parametrize("lang", ["en", "ja"])
@pytest.mark.parametrize("builtin", [True, False])
def test_mapping_disclosure_keeps_actual_selections_and_warnings(lang, builtin):
    def view(lang, builtin):
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        st.session_state["lang"] = lang
        if builtin:
            st.session_state["_loaded_sample_scenario_key"] = "writing_essay"
        st.session_state.setdefault("facets_mode_output", {"analysis_id": "existing-fit"})
        data = pd.DataFrame({
            "Person": ["P1", "P1", "P2", "P2"], "Score": [1, 2, 2, 1],
            "Rater": ["R1", "R2", "R1", "R2"], "Task": ["T1", "T2", "T1", "T2"],
            "Alternate": [2, 1, 1, 2],
        })
        st.session_state["mapping_label"] = app.t("sidebar_estimation.column_mapping_subheader")
        st.session_state["facet_warning"] = app.t("sidebar_run_setup.needs_two_facets_warning")
        # Render the production setup, without a new estimation or the result body.
        app.run_facets_mode(app.load_core_namespace(), data, help_surface_active=True)

    at = AppTest.from_function(view, args=(lang, builtin)).run(timeout=45)
    assert not at.exception
    editor = next(e for e in at.sidebar.expander if e.label == at.session_state["mapping_label"])
    assert editor.proto.expanded is (not builtin)
    score = next(item for item in editor.selectbox if item.key.startswith("facets_mode_score_col_"))
    score.select("Alternate").run(timeout=45)
    assert not at.exception
    direct_captions = "\n".join(item.value for item in at.sidebar.caption)
    assert "Score" in direct_captions and "Alternate" in direct_captions
    facets = next(item for item in at.multiselect if item.key.startswith("facets_mode_facet_cols_"))
    facets.set_value(["Rater"]).run(timeout=45)
    assert not at.exception
    assert at.session_state["facet_warning"] in [w.value for w in at.sidebar.warning]
    assert at.radio(key="facets_mode_model_type").value == "RSM"
    assert at.radio(key="facets_mode_estimation_method").value == "JMLE"
    assert at.session_state["facets_mode_output"] == {"analysis_id": "existing-fit"}
