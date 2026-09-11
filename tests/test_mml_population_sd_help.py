"""Tests for the Help / UX surfaces of the MML person population-SD feature.

The estimator was shipped first; these tests pin the user-facing explanation:
the Interpretation-Guide help table, the help-popover topic, the glossary
entry, and the convergence-panel summary that surfaces the estimated sigma
(with unqualified SE/CI withheld) for a free-SD run.
"""

from __future__ import annotations

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app


class _UIRecorder:
    """Stand-in for the module-level ``st`` that records rendered text."""

    def __init__(self, real):
        self._real = real
        self.messages: list[str] = []

    def __getattr__(self, name):
        return getattr(self._real, name)

    def metric(self, label="", value="", *a, **k):
        self.messages.append(f"{label}={value}")

    def caption(self, msg="", *a, **k):
        self.messages.append(str(msg))

    def info(self, msg="", *a, **k):
        self.messages.append(str(msg))

    def markdown(self, msg="", *a, **k):
        self.messages.append(str(msg))

    def dataframe(self, *a, **k):
        return None

    def columns(self, spec, *a, **k):
        n = spec if isinstance(spec, int) else len(spec)
        return [self for _ in range(n)]

    def _ctx(self):
        rec = self

        class _Ctx:
            def __enter__(self_inner):
                return rec

            def __exit__(self_inner, *exc):
                return False

        return _Ctx()

    def popover(self, *a, **k):
        return self._ctx()

    def expander(self, *a, **k):
        return self._ctx()


def _render(config: dict) -> list[str]:
    rec = _UIRecorder(app.st)
    saved = app.st
    app.st = rec
    try:
        app._render_population_sd_summary({"config": config})
    finally:
        app.st = saved
    return rec.messages


def test_help_table_has_fixed_and_free_rows():
    tbl = app.mml_population_sd_help_table()
    assert isinstance(tbl, pd.DataFrame)
    assert len(tbl) == 2
    assert len(tbl.columns) == 4
    joined = " ".join(tbl.astype(str).to_numpy().ravel().tolist())
    assert "EM" in joined  # free-SD row mentions the EM engine


def test_help_popover_topic_registered():
    topic = app._HELP_POPOVER_LIBRARY.get("mml_person_sd")
    assert topic is not None
    assert {"title", "what", "how", "watch"}.issubset(topic)
    assert "EM" in topic["how"] or "EM" in topic["watch"]


def test_glossary_entry_mentions_estimation():
    entry = app._MFRM_GLOSSARY.get("population prior sd", "")
    assert "estimat" in entry.lower()


def test_summary_renders_sigma_but_withholds_legacy_se_ci():
    msgs = _render({
        "method": "MML", "estimate_population_sd": True,
        "estimated_population_sd": 1.49, "population_sd_se": 0.09,
        "population_sd_ci": [1.32, 1.66],
        "population_sd_engine_notice": "overrode 'auto' to 'em'.",
    })
    assert any("1.49" in m for m in msgs)
    assert not any("0.09" in m or "1.32" in m or "1.66" in m for m in msgs)
    assert any("withheld" in m for m in msgs)
    assert any("overrode" in m for m in msgs)  # engine-override notice surfaced


def test_summary_handles_missing_se_gracefully():
    msgs = _render({
        "method": "MML", "estimate_population_sd": True,
        "estimated_population_sd": 0.05, "population_sd_se": float("nan"),
        "population_sd_ci": [float("nan"), float("nan")],
    })
    assert any("0.05" in m for m in msgs)  # sigma still shown even when SE is unavailable


@pytest.mark.parametrize("lang,hold", [("en", "withheld"), ("ja", "保留")])
def test_free_sd_summary_withholds_legacy_intervals_in_both_languages(lang, hold):
    def view(lang):
        import streamlit as st
        import streamlit_app as app

        st.session_state["lang"] = lang
        app._render_population_sd_summary({"config": {
            "method": "MML", "estimate_population_sd": True,
            "estimated_population_sd": 1.49, "population_sd_se": 0.09,
            "population_sd_ci": [1.32, 1.66],
        }})

    at = AppTest.from_function(view, args=(lang,)).run(timeout=30)
    assert not at.exception
    assert [metric.value for metric in at.metric] == ["1.49"]
    assert any(hold in caption.value for caption in at.caption)
    text = "\n".join(element.value for element in [*at.markdown, *at.caption])
    assert not any(value in text for value in ("0.09", "1.32", "1.66"))


def test_summary_quiet_caption_for_fixed_run():
    msgs = _render({"method": "MML", "estimate_population_sd": False, "population_prior_sd": 1.0})
    assert any("1.00" in m for m in msgs)


def test_summary_silent_for_non_mml_run():
    assert _render({"method": "JMLE"}) == []


def test_registered_self_test_passes():
    app._self_test_mml_population_sd_help_surfaces()
