"""End-to-end smoke tests across every built-in sample scenario.

The scenario registry has two distinct contracts:

1. Every built-in scenario must fit through the standalone Python engine.
2. Representative scenarios must still render through Streamlit AppTest after
   switching the scenario picker and clicking the run button.

Full AppTest renders every Streamlit tab, including publication/export helpers.
Running that UI path for every large scenario became a CI-timeout risk even
though the underlying estimator completes quickly. The engine test keeps all
scenario data shapes covered; the representative AppTest path keeps widget
constraint regressions visible.
"""

from __future__ import annotations

import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app


APPTEST_TIMEOUT = 180
APPTEST_SCENARIOS = ("writing_essay", "writing_with_missing")
WORKFLOW_CONTROL_KINDS = (
    "button",
    "checkbox",
    "radio",
    "selectbox",
    "slider",
    "text_area",
    "text_input",
    "number_input",
    "multiselect",
)
INITIAL_CONTROL_BUDGET = {
    "button": 4,
    "checkbox": 7,
    "radio": 7,
    "selectbox": 13,
    "slider": 4,
    "text_area": 2,
    "text_input": 1,
    "number_input": 9,
    "multiselect": 4,
}
RESULT_CONTROL_BUDGET = {
    "button": 9,
    "checkbox": 8,
    "radio": 7,
    "selectbox": 15,
    "slider": 4,
    "text_area": 2,
    "text_input": 1,
    "number_input": 9,
    "multiselect": 4,
}


def assert_control_topology_within_budget(
    app_test: AppTest,
    budget: dict[str, int],
) -> None:
    """Require an explicit IA budget change before adding rendered controls."""

    counts = {
        kind: len(getattr(app_test, kind))
        for kind in WORKFLOW_CONTROL_KINDS
    }
    excess = {
        kind: (counts[kind], maximum)
        for kind, maximum in budget.items()
        if counts[kind] > maximum
    }
    assert not excess, f"rendered workflow control budget exceeded: {excess}"


@pytest.mark.parametrize("scenario_key", list(app.SAMPLE_DATA_SCENARIOS.keys()))
def test_scenario_estimation_pipeline_completes(scenario_key: str):
    """Every built-in scenario should fit through the standalone engine."""
    df = app.sample_mfrm_data_by_key(scenario_key, seed=3)
    facet_cols = [c for c in df.columns if c not in {"Person", "Score"}]
    result = app.mfrm_estimate(
        data=df,
        person_col="Person",
        score_col="Score",
        facet_cols=facet_cols,
        model="RSM",
        method="JMLE",
        rating_min=int(df["Score"].min()),
        rating_max=int(df["Score"].max()),
        maxit=80,
        reltol=1e-4,
    )
    summary = result.get("summary")
    assert summary is not None and not summary.empty
    assert result.get("facets") is not None
    assert result.get("steps") is not None
    assert result.get("convergence") is not None


@pytest.mark.parametrize("scenario_key", APPTEST_SCENARIOS)
def test_representative_scenarios_render_without_streamlit_exception(scenario_key: str):
    """Representative UI path: load scenario, click Run, re-render cleanly."""
    at = AppTest.from_file("streamlit_app.py").run(timeout=APPTEST_TIMEOUT)
    assert not at.exception, (
        f"initial render already raised for {scenario_key!r}: "
        f"{[e.value for e in at.exception]}"
    )
    at.button(key="onboarding_skip_guide").click()
    at.run(timeout=APPTEST_TIMEOUT)

    assert at.selectbox(key="data_source_class").value == "sample"
    at.selectbox(key="data_source_scenario").set_value(scenario_key)

    at.run(timeout=APPTEST_TIMEOUT)
    assert not at.exception, (
        f"after scenario switch for {scenario_key!r}: "
        f"{[e.value for e in at.exception]}"
    )
    assert at.session_state["data_source_flat"] == f"scenario:{scenario_key}"
    assert_control_topology_within_budget(at, INITIAL_CONTROL_BUDGET)

    run_button = at.button(key="facets_mode_run_primary")
    assert run_button.label == "Run this analysis"
    assert not any(
        button.label and "Run FACETS-mode" in button.label
        for button in at.sidebar.button
    )
    run_button.click()
    at.run(timeout=APPTEST_TIMEOUT)
    assert not at.exception, (
        f"{scenario_key!r} crashed during or after Run: "
        f"{[e.value for e in at.exception]}"
    )
    subheaders = [str(item.value) for item in at.subheader]
    assert subheaders.count("Choose what you want to do now") == 1
    assert "First-read overview" not in subheaders
    assert any(
        str(item.value) == "**Interpretation readiness**"
        for item in at.markdown
    )
    visible_status_metrics = {
        (str(item.label), str(item.value))
        for item in at.metric
        if str(item.label) in {"Pause", "Caution", "Review", "OK"}
    }
    assert {label for label, _ in visible_status_metrics} == {
        "Pause", "Caution", "Review", "OK"
    }
    assert at.button(key="guided_action_hub_overview_primary").label == "Open recommended section"
    download_labels = {str(item.label) for item in at.get("download_button")}
    assert {"Download all results (ZIP)", "Download all results (Excel)"}.issubset(download_labels)
    assert any(
        "Stays visible while you scroll" in str(item.value)
        for item in at.caption
    )
    assert_control_topology_within_budget(at, RESULT_CONTROL_BUDGET)

    at.button(key="guided_action_hub_overview_primary").click()
    at.run(timeout=APPTEST_TIMEOUT)
    assert not at.exception, (
        f"{scenario_key!r} crashed after the recommended one-click route: "
        f"{[e.value for e in at.exception]}"
    )
    assert at.session_state["guided_essential_section"] != "start"


def test_scenario_names_match_registry():
    """Sanity: the scenario smoke tests still cover the intended registry."""
    keys = set(app.SAMPLE_DATA_SCENARIOS.keys())
    assert keys == {
        "writing_essay",
        "large_writing_pca",
        "speaking_test",
        "clinical_osce",
        "writing_with_missing",
        "music_peer_rating",
        "reading_testlet_binary",
    }, f"Scenario registry drifted: {keys}"


def test_guided_and_advanced_defaults_produce_the_same_estimates() -> None:
    """Progressive disclosure must not create a second statistical default."""

    at = AppTest.from_file("streamlit_app.py").run(timeout=APPTEST_TIMEOUT)
    at.button(key="onboarding_skip_guide").click()
    at.run(timeout=APPTEST_TIMEOUT)
    at.button(key="facets_mode_run_primary").click()
    at.run(timeout=APPTEST_TIMEOUT)
    assert not at.exception

    guided = at.session_state["facets_mode_output"]["result"]
    guided_params = guided["params"]
    guided_theta = np.asarray(guided_params["theta"], dtype=float).copy()
    guided_steps = np.asarray(guided_params["steps"], dtype=float).copy()
    guided_facets = {
        name: np.asarray(values, dtype=float).copy()
        for name, values in guided_params["facets"].items()
    }

    at.radio(key="facets_mode_workflow_mode").set_value("Advanced controls")
    at.run(timeout=APPTEST_TIMEOUT)
    sidebar_runs = [
        button for button in at.sidebar.button
        if button.label and "Run FACETS-mode" in button.label
    ]
    assert len(sidebar_runs) == 1
    sidebar_runs[0].click()
    at.run(timeout=APPTEST_TIMEOUT)
    assert not at.exception

    advanced = at.session_state["facets_mode_output"]["result"]
    advanced_params = advanced["params"]
    np.testing.assert_allclose(
        advanced_params["theta"], guided_theta, rtol=0, atol=1e-12
    )
    np.testing.assert_allclose(
        advanced_params["steps"], guided_steps, rtol=0, atol=1e-12
    )
    assert set(advanced_params["facets"]) == set(guided_facets)
    for name, expected in guided_facets.items():
        np.testing.assert_allclose(
            advanced_params["facets"][name], expected, rtol=0, atol=1e-12
        )
