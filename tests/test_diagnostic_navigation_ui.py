"""Recommended actions must open the named diagnostic and the correct PCA scope."""

import pytest
from streamlit.testing.v1 import AppTest


@pytest.mark.parametrize("lang", ["en", "ja"])
@pytest.mark.parametrize("target,panel", [
    ("Dimensionality", "dimensionality"),
    ("Bias / Interaction", "bias_interaction"),
    ("Categories / Steps", "categories_steps"),
])
def test_recommended_action_opens_its_diagnostic_without_changing_fit(lang, target, panel):
    def view(lang, target):
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app

        st.session_state["lang"] = lang
        st.session_state.setdefault("facets_mode_output", {"analysis_id": "same-fit"})
        st.session_state.setdefault("dimensionality_panel", "Rater")
        plan = pd.DataFrame([{
            "Priority": "P1", "SourceOrder": 1, "Area": "Assumptions",
            "Check": target, "Status": "Caution", "WhatItSays": "Inspect this evidence.",
            "NextAction": "Open the named diagnostic.", "DetailLocation": target,
        }])
        def render(name):
            return lambda *args, **kwargs: st.write(name)
        with patch.object(app, "_render_guided_start_section", render("start")), \
             patch.object(app, "_show_pca_panel", render("dimensionality")), \
             patch.object(app, "render_dimtest_panel"), \
             patch.object(app, "show_bias_section", render("bias_interaction")), \
             patch.object(app, "show_categories_section", render("categories_steps")):
            app._render_guided_essential_tabs(
                {}, pd.DataFrame(), st.session_state["facets_mode_output"], {"pca": {"eigenvalues": [3.74, 1.2]}}, {},
                pd.DataFrame(), pd.DataFrame(), ["Rater", "Task"], "Person", "Score",
                None, None, result_compute_pca=True, result_render_plots=False,
                result_generate_figures=False, action_plan=plan,
            )

    at = AppTest.from_function(view, args=(lang, target)).run(timeout=45)
    assert not at.exception
    assert at.button_group(key="guided_essential_section").value == "first_read"
    at.button(key="guided_action_hub_overview_primary").click().run(timeout=45)
    assert not at.exception
    assert at.session_state["guided_essential_section"] == "diagnostics"
    assert at.session_state["guided_diagnostics_panel"] == panel
    assert panel in [item.value for item in at.markdown]
    assert not any("Session State API" in item.value for item in at.warning)
    if panel == "dimensionality":
        assert at.session_state["dimensionality_panel"] == "overall"
    assert at.session_state["facets_mode_output"] == {"analysis_id": "same-fit"}
    at.button_group(key="guided_essential_section").set_value("first_read").run()
    assert not at.exception
    assert at.session_state["guided_essential_section"] == "first_read"
    assert at.session_state["facets_mode_output"] == {"analysis_id": "same-fit"}


@pytest.mark.parametrize("lang", ["en", "ja"])
def test_pca_opens_overall_before_rater_scope_and_preserves_cached_payload(lang):
    def view(lang):
        from unittest.mock import patch
        import streamlit as st
        import streamlit_app as app

        st.session_state["lang"] = lang
        pca = {"eigenvalues": [3.74, 1.2]}
        def render(payload, *, mode, **kwargs):
            assert payload == pca
            st.write("scope:" + mode + ":" + kwargs.get("facet_name", "all"))
        with patch.object(app, "_show_pca_panel", render), \
             patch.object(app, "render_dimtest_panel"):
            app.show_dimensionality_section({"pca": pca}, ["Rater", "Task"])
        assert pca == {"eigenvalues": [3.74, 1.2]}

    at = AppTest.from_function(view, args=(lang,)).run(timeout=45)
    assert not at.exception
    assert at.selectbox(key="dimensionality_panel").value == "overall"
    assert "scope:overall:all" in [item.value for item in at.markdown]
    assert all(not item.proto.expanded for item in at.expander)
    at.selectbox(key="dimensionality_panel").select("Rater").run()
    assert not at.exception
    assert "scope:facet:Rater" in [item.value for item in at.markdown]
