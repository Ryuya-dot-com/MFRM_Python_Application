"""Accepted MML facet constraints and rejection of ignored Person constraints."""
import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_analysis_journeys import _run
from tests.test_help_contextual_app_integration import _locale_value


def _args(engine="direct", model="RSM"):
    return dict(
        data=app._make_self_test_gradient_data(), person_col="Person",
        score_col="Score", facet_cols=["Rater", "Task"], model=model,
        method="MML", mml_engine=engine, step_facet="Task" if model != "RSM" else None,
        maxit=100, reltol=1e-6, quad_points=15,
    )


@pytest.mark.parametrize("engine", ["EM", "direct", "hybrid", "auto"])
@pytest.mark.parametrize("constraint", ["fixed", "group", "dummy"])
def test_person_constraints_stop_before_mml_optimization(engine, constraint, monkeypatch):
    args = _args(engine)
    if constraint == "fixed":
        args["anchor_df"] = pd.DataFrame({"Facet": ["Person"], "Level": ["P1"], "Anchor": [2.]})
    elif constraint == "group":
        args["group_anchor_df"] = pd.DataFrame({
            "Facet": ["Person"], "Level": ["P1"], "Group": ["Common"], "GroupValue": [2.],
        })
    else:
        args["dummy_facets"] = ["Person"]

    def forbidden(*args, **kwargs):
        pytest.fail("Unsupported Person constraints reached optimization")

    monkeypatch.setattr(app, "minimize", forbidden)
    with pytest.raises(ValueError, match="MML does not support Person anchors"):
        app.mfrm_estimate(**args)


@pytest.mark.parametrize("engine", ["EM", "direct", "hybrid", "auto"])
@pytest.mark.parametrize("model", ["RSM", "PCM", "GPCM"])
@pytest.mark.parametrize("constraint", ["fixed", "group"])
def test_non_person_constraints_reach_returned_mml_estimates(engine, model, constraint):
    args = _args(engine, model)
    if constraint == "fixed":
        args["anchor_df"] = pd.DataFrame({"Facet": ["Rater"], "Level": ["R1"], "Anchor": [-.4]})
    else:
        args["group_anchor_df"] = pd.DataFrame({
            "Facet": ["Rater", "Rater"], "Level": ["R1", "R2"],
            "Group": ["Reference", "Reference"], "GroupValue": [.2, .2],
        })
    result = app.mfrm_estimate(**args)
    estimates = result["facets"]["others"].query("Facet == 'Rater'").set_index("Level").Estimate
    assert np.isfinite(result["summary"].iloc[0].LogLik)
    if constraint == "fixed":
        assert estimates["R1"] == pytest.approx(-.4, abs=1e-12)
    else:
        assert estimates.loc[["R1", "R2"]].mean() == pytest.approx(.2, abs=1e-12)
    assert result["anchor_drift"]["summary"].iloc[0].ReviewCount == 0


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_mml_anchor_boundary_is_explained_and_invalid_run_has_no_result(lang):
    at = AppTest.from_string("import streamlit_app as app\napp.main()", default_timeout=120)
    _run(at)
    at.button(key="onboarding_dismiss").click()
    _run(at)
    at.text_area(key="paste_data_text").set_value(app._make_self_test_gradient_data().to_csv(index=False))
    _run(at)
    at.radio(key="facets_mode_workflow_mode").set_value("Advanced controls")
    _run(at)
    at.radio(key="facets_mode_estimation_method").set_value("MML")
    at.radio(key="lang").set_value(lang)
    _run(at)
    caption = _locale_value(lang, "sidebar_advanced.mml_anchor_scope_caption")
    assert caption in [c.value for c in at.caption]
    at.text_area(key="facets_mode_anchor_text").set_value("Facet,Level,Anchor\nPerson,P1,2.0")
    next(b for b in at.sidebar.button if b.label == _locale_value(lang, "sidebar_perf.run_button")).click()
    _run(at)
    assert _locale_value(lang, "sidebar_advanced.mml_person_anchor_error") in [e.value for e in at.error]
    assert "facets_mode_output" not in at.session_state or at.session_state["facets_mode_output"] is None
