"""Clean input audits collapse; excluded and invalid rows stay visible."""
import pytest
from streamlit.testing.v1 import AppTest


@pytest.mark.parametrize("case", ["clean", "missing", "fractional"])
def test_setup_audit_disclosure_keeps_problem_rows_and_downloads(case):
    def workspace(case):
        import numpy as np
        import streamlit as st
        import streamlit_app as app

        st.session_state["lang"] = "en"
        data = app.cached_sample_mfrm_data_by_key("writing_essay", seed=20240101).copy()
        if case == "missing":
            data.loc[data.index[0], "Score"] = np.nan
        elif case == "fractional":
            data["Score"] = data["Score"].astype(float)
            data.loc[data.index[0], "Score"] = .5
        app.render_analysis_setup_workspace(
            data, person_col="Person", score_col="Score", facet_cols=["Rater", "Task", "Criterion"],
            weight_col=None, workflow_mode="Guided defaults", model_type="RSM", est_method="JMLE",
            analysis_depth="Standard (recommended)",
        )

    at = AppTest.from_function(workspace, args=(case,)).run(timeout=45)
    assert not at.exception
    assert at.button(key="facets_mode_run_primary")
    assert at.button(key="facets_mode_run_primary").disabled == (case == "fractional")
    downloads = at.get("download_button")
    assert any(item.label == "Download full row audit CSV" for item in downloads)
    audit_folds = [item for item in at.expander if item.label == "Response-data audit"]
    mapping = next(item for item in at.expander if item.label == "Inspect likelihood-row mapping")
    assert mapping.proto.expanded
    if case == "clean":
        assert len(audit_folds) == 1 and not audit_folds[0].proto.expanded
        assert any("Every input row is usable" in item.value for item in audit_folds[0].success)
        assert not at.warning and not at.error
        quality = next(item for item in at.expander if item.label == "Show data quality detail")
        assert not quality.proto.expanded
        assert len(quality.dataframe) == 2
        assert not any("next action" in item.value.lower() for item in at.caption)
    else:
        assert not audit_folds
        banners = at.warning if case == "missing" else at.error
        expected = "excluded before likelihood calculation" if case == "missing" else "fractional score values"
        assert any(expected in item.value for item in banners)
        assert any(item.label == "Download excluded rows CSV" for item in downloads)
        assert any("next action" in item.value.lower() for item in at.caption)


@pytest.mark.parametrize("lang", ["en", "ja"])
def test_small_input_run_control_and_resource_limit(lang):
    def workspace(lang):
        import streamlit as st
        import streamlit_app as app
        from tests.test_readiness_report import _small_crossed_data

        st.session_state["lang"] = lang
        resource_block = st.checkbox("Resource block")
        app.render_analysis_setup_workspace(
            _small_crossed_data(), person_col="Person", score_col="Score",
            facet_cols=["Rater", "Task"], weight_col=None,
            workflow_mode="Guided defaults", model_type="RSM", est_method="JMLE",
            analysis_depth="Standard", run_disabled=resource_block,
        )
        st.session_state["count_note"] = app.t("input_readiness.counts_detail")
        st.session_state["count_action"] = app.t("input_readiness.action_review_sample_size")

    at = AppTest.from_function(workspace, args=(lang,)).run(timeout=45)
    assert not at.exception and not at.error
    assert at.warning and not at.button(key="facets_mode_run_primary").disabled
    text = "\n".join(frame.value.to_string() for frame in at.dataframe)
    assert at.session_state["count_note"] in text
    assert any(at.session_state["count_action"] in item.value for item in at.markdown)
    at.checkbox[0].check().run()
    assert not at.exception and at.button(key="facets_mode_run_primary").disabled


def test_full_app_blocks_invalid_input_in_both_modes_and_forced_triggers(monkeypatch):
    import streamlit_app as app
    from tests.test_readiness_report import _small_crossed_data

    calls = []

    def unexpected_estimate(*args, **kwargs):
        calls.append(True)
        raise AssertionError("Invalid input reached the estimator")

    core = app.load_core_namespace().copy()
    core["mfrm_estimate"] = unexpected_estimate
    monkeypatch.setattr(app, "load_core_namespace", lambda: core)
    at = AppTest.from_string("import streamlit_app as app\napp.main()", default_timeout=45).run()
    at.button(key="onboarding_dismiss").click().run()
    data = _small_crossed_data()
    at.text_area(key="paste_data_text").set_value(data.to_csv(index=False)).run()
    assert not at.exception and not at.button(key="facets_mode_run_primary").disabled
    data["Score"] = data.Score.astype(float)
    data.loc[0, "Score"] = .5
    at.text_area(key="paste_data_text").set_value(data.to_csv(index=False)).run()
    assert not at.exception and at.button(key="facets_mode_run_primary").disabled
    at.radio(key="facets_mode_workflow_mode").set_value("Advanced controls").run()
    assert not at.exception
    run_button = next(button for button in at.sidebar.button if button.label == app.t("sidebar_perf.run_button"))
    assert run_button.disabled
    for trigger in ("_facets_mode_force_rerun", "_onboarding_quickstart_fired"):
        at.session_state[trigger] = True
        at.run()
        assert not at.exception
        assert not calls
    at.text_area(key="paste_data_text").set_value(_small_crossed_data().to_csv(index=False)).run()
    assert not at.exception
    run_button = next(button for button in at.sidebar.button if button.label == app.t("sidebar_perf.run_button"))
    assert not run_button.disabled
    assert not calls
