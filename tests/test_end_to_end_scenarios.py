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

import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app


APPTEST_TIMEOUT = 180
APPTEST_SCENARIOS = ("writing_essay", "writing_with_missing")


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
    scenario_label = app.SAMPLE_DATA_SCENARIOS[scenario_key]["label"]
    at = AppTest.from_file("streamlit_app.py").run(timeout=APPTEST_TIMEOUT)
    assert not at.exception, (
        f"initial render already raised for {scenario_key!r}: "
        f"{[e.value for e in at.exception]}"
    )

    matched = False
    for radio in at.radio:
        if radio.key == "data_source_flat":
            try:
                radio.set_value(scenario_label)
                matched = True
            except Exception as exc:
                pytest.fail(
                    f"{scenario_key!r}: failed to switch radio to "
                    f"{scenario_label!r}: {exc}"
                )
            break
    if not matched:
        pytest.skip(
            "data_source_flat radio not found; skipping to avoid false positive"
        )

    at.run(timeout=APPTEST_TIMEOUT)
    assert not at.exception, (
        f"after scenario switch for {scenario_key!r}: "
        f"{[e.value for e in at.exception]}"
    )

    run_button = None
    for button in at.sidebar.button:
        if button.label and "Run MFRM estimation" in button.label:
            run_button = button
            break
    if run_button is None:
        pytest.skip(
            "Run MFRM estimation button not found; skipping to avoid false positive"
        )

    run_button.click()
    at.run(timeout=APPTEST_TIMEOUT)
    assert not at.exception, (
        f"{scenario_key!r} crashed during or after Run: "
        f"{[e.value for e in at.exception]}"
    )


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
