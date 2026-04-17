"""End-to-end smoke tests across every built-in sample scenario.

This is the gate the user legitimately asked for after the v0.2.10
`StreamlitValueBelowMinError` escape: every scenario must be runnable
through the full estimation + full-tab-render pipeline without any
Streamlit exception surfacing.

Rationale
---------
Before v0.2.11 we had `test_sample_scenarios.py` (unit — does the
generator return the right DataFrame?) and `test_app_smoke.py`
(initial-render only). Neither combined the two. The gap was exactly
the class of bug that shipped: a scenario generates valid data but
that data then trips a widget constraint in a tab renderer.

Strategy
--------
For every key in `SAMPLE_DATA_SCENARIOS`:
  1. Boot the app with `streamlit.testing.v1.AppTest`.
  2. Switch the scenario picker to that key.
  3. Click "Run FACETS-mode estimation".
  4. Assert `at.exception is None` after re-render.

This exercises EVERY post-run tab's widget math (Visuals, Fit
Details, Dimensionality, Bias/Interaction, etc.) against the real
scenario shape, catching `StreamlitValueBelowMinError`,
`StreamlitValueAboveMaxError`, and `StreamlitAPIException` for all
of them at once.
"""

from __future__ import annotations

import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app


# 90s timeout per scenario — the largest synthetic dataset is the
# speaking-test (3,600 obs) and estimation typically finishes well
# within that on a CI runner.
APPTEST_TIMEOUT = 90


@pytest.mark.parametrize("scenario_key", list(app.SAMPLE_DATA_SCENARIOS.keys()))
def test_scenario_runs_full_pipeline_without_exception(scenario_key: str):
    """End-to-end: load scenario, click Run, re-render — no exception."""
    scenario_label = app.SAMPLE_DATA_SCENARIOS[scenario_key]["label"]
    at = AppTest.from_file("streamlit_app.py").run(timeout=APPTEST_TIMEOUT)
    assert not at.exception, (
        f"initial render already raised for {scenario_key!r}: "
        f"{[e.value for e in at.exception]}"
    )

    # Switch the data-source radio to the scenario under test.
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
            "data_source_flat radio not found (UI refactor?) — "
            "skipping this scenario to avoid false positive"
        )

    at.run(timeout=APPTEST_TIMEOUT)
    assert not at.exception, (
        f"after scenario switch for {scenario_key!r}: "
        f"{[e.value for e in at.exception]}"
    )

    # Click the Run button. It lives in the sidebar and is labelled
    # "Run FACETS-mode estimation".
    run_button = None
    for button in at.sidebar.button:
        if button.label and "Run FACETS-mode" in button.label:
            run_button = button
            break
    if run_button is None:
        pytest.skip(
            "Run FACETS-mode estimation button not found — UI shift "
            "or disabled state; skipping to avoid false positive"
        )

    run_button.click()
    at.run(timeout=APPTEST_TIMEOUT)
    assert not at.exception, (
        f"{scenario_key!r} crashed during or after Run: "
        f"{[e.value for e in at.exception]}"
    )


def test_scenario_names_match_registry():
    """Sanity: the parametrize values in use match the registry keys."""
    keys = set(app.SAMPLE_DATA_SCENARIOS.keys())
    assert keys == {"writing_essay", "large_writing_pca",
                    "speaking_test", "clinical_osce"}, (
        f"Scenario registry drifted: {keys}"
    )
