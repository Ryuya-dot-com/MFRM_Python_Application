"""Rendered bilingual help stays readable without triggering estimation."""

import re

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_help_contextual_app_integration import _block_analysis_calls, _locale_value


PANELS = (
    "Quick Start", "Analysis Workflow", "Interpretation Guide", "Rater Effects",
    "Rating Scale Guide", "Glossary", "Reporting Guide", "Troubleshooting",
)
INTERNAL_COPY = re.compile(
    r"\b(?:G[0-6]|AnalysisID|EvidenceRecords?|ReasonCode|ComputationState|"
    r"StabilityState|backend|runner|manifest|pipeline)\b|実装上|適格性ゲート",
    re.IGNORECASE,
)


@pytest.mark.parametrize("language", ["ja", "en"])
def test_all_help_panels_keep_one_tutorial_and_hide_internal_copy(monkeypatch, language):
    calls = _block_analysis_calls(monkeypatch)
    for panel in PANELS:
        at = AppTest.from_string(
            "import streamlit_app as app\n"
            "app.show_tutorial()\n"
            "app.show_help_section(force_full=True)",
            default_timeout=120,
        )
        at.session_state["lang"] = language
        at.session_state["help_panel"] = panel
        at.run()
        assert not at.exception, (language, panel, at.exception)
        markdown = [element.value for element in at.markdown]
        assert markdown.count(_locale_value(language, "tutorial.section_estimation_method")) == 1
        prose = "\n".join(
            str(element.value)
            for kind in ("markdown", "caption", "info", "warning", "header", "subheader")
            for element in getattr(at, kind)
        )
        tables = "\n".join(element.value.to_string(index=False) for element in at.dataframe)
        assert not INTERNAL_COPY.search(prose + "\n" + tables), (language, panel)
        if panel == "Analysis Workflow":
            assert _locale_value(language, "help.analysis_workflow_body") in markdown
            assert len(at.expander) > 1
            assert all(not element.proto.expanded for element in at.expander)
        assert calls == {"estimate": 0, "refit": 0}


@pytest.mark.parametrize("language", ["ja", "en"])
def test_model_comparison_explains_holds_and_keeps_codes_in_export(monkeypatch, language):
    _block_analysis_calls(monkeypatch)
    monkeypatch.setattr(app, "_model_choice_can_compare", lambda result: (True, ""))
    monkeypatch.setattr(app, "_model_choice_cache_key", lambda result: "help-copy-test")
    exported = []
    original_to_csv = app.to_csv_bytes

    def record_export(frame):
        exported.append(frame.copy(deep=True))
        return original_to_csv(frame)

    monkeypatch.setattr(app, "to_csv_bytes", record_export)
    at = AppTest.from_string(
        'import streamlit_app as app\n'
        'app.render_model_choice_guidance({"config": {"method": "JMLE"}})',
        default_timeout=120,
    )
    at.session_state["lang"] = language
    at.session_state["model_choice_comparison::help-copy-test"] = {
        "available": True,
        "comparison": pd.DataFrame({"Model": ["RSM", "PCM"], "AIC": [12.0, 14.0]}),
        "lr_tests": pd.DataFrame({"Null": ["RSM"], "Alternative": ["PCM"], "p_value": [0.01]}),
        "recommendation": {"tier": "strong", "model": "RSM"},
        "caveat": "G2 gate: ReasonCode from an older cached comparison",
    }
    at.run()
    assert not at.exception
    prose = "\n".join(str(e.value) for kind in ("markdown", "caption", "warning") for e in getattr(at, kind))
    assert not INTERNAL_COPY.search(prose)
    assert _locale_value(language, "report_tables.model_choice_recommendation_withheld_template") in prose
    assert not at.success
    for table in at.dataframe:
        assert not {"ReasonCode", "RoadmapIssueID", "QualificationStatus", "PublicConclusionAllowed", "p_value"} & set(table.value.columns)
    assert len(exported) == 1
    assert {"ReasonCode", "RoadmapIssueID", "PublicConclusionAllowed", "p_value"} <= set(exported[0].columns)
    assert not exported[0]["PublicConclusionAllowed"].fillna(False).any()
