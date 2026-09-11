from streamlit.testing.v1 import AppTest


def test_app_initial_render():
    at = AppTest.from_file("streamlit_app.py").run(timeout=30)
    assert not at.exception
    assert any("MFRM | Rating analysis" in title.value for title in at.title)
    about = next(e for e in at.expander if e.label == "About this beta and its limitations")
    assert not about.proto.expanded
    assert any("source commit:" in caption.value for caption in about.caption)
    assert not any("Data privacy" in warning.value for warning in at.warning)
    assert not any("Keyboard shortcuts" in str(expander.label) for expander in at.expander)


def test_user_data_source_promotes_privacy_warning():
    at = AppTest.from_file("streamlit_app.py").run(timeout=30)
    at.button(key="onboarding_dismiss").click()
    at.run(timeout=30)

    assert not at.exception
    assert at.session_state["data_source_flat"] == "paste"
    assert any("pasted and uploaded rating files" in warning.value for warning in at.warning)


def test_two_level_source_picker_keeps_sample_detail_visible_and_stable() -> None:
    at = AppTest.from_file("streamlit_app.py").run(timeout=30)
    at.button(key="onboarding_skip_guide").click()
    at.run(timeout=30)

    assert at.radio(key="data_source_class").value == "sample"
    assert at.selectbox(key="data_source_scenario").value == "writing_essay"
    at.selectbox(key="data_source_scenario").set_value("clinical_osce")
    at.run(timeout=30)
    assert at.session_state["data_source_flat"] == "scenario:clinical_osce"

    at.radio(key="data_source_class").set_value("paste")
    at.run(timeout=30)
    assert at.session_state["data_source_flat"] == "paste"
    assert not any(item.key == "data_source_scenario" for item in at.selectbox)
    assert any("pasted and uploaded rating files" in warning.value for warning in at.warning)

    at.radio(key="data_source_class").set_value("sample")
    at.run(timeout=30)
    assert at.selectbox(key="data_source_scenario").value == "clinical_osce"
    assert at.session_state["data_source_flat"] == "scenario:clinical_osce"


def test_legacy_flat_source_state_migrates_into_two_level_picker() -> None:
    at = AppTest.from_file("streamlit_app.py").run(timeout=30)
    at.button(key="onboarding_skip_guide").click()
    at.run(timeout=30)

    at.session_state["data_source_flat"] = "upload"
    at.run(timeout=30)
    assert at.radio(key="data_source_class").value == "upload"
    assert at.session_state["data_source_flat"] == "upload"


def test_guided_setup_only_renders_first_run_decisions() -> None:
    at = AppTest.from_file("streamlit_app.py").run(timeout=30)
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
    at = AppTest.from_file("streamlit_app.py").run(timeout=30)
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
    at = AppTest.from_file("streamlit_app.py").run(timeout=30)
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
    at = AppTest.from_file("streamlit_app.py").run(timeout=30)
    at.button(key="onboarding_skip_guide").click()
    at.run(timeout=30)
    at.radio(key="lang").set_value("ja")
    at.run(timeout=30)

    assert not at.exception
    subheaders = {item.value for item in at.subheader}
    assert {"解析方法の選択", "解析内容"}.issubset(subheaders)
    assert "Anchor 制約" not in subheaders
    assert any(
        "技術設定には文書化された標準値" in str(item.value)
        for item in at.caption
    )


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
