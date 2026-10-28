"""Person-level holds survive reporting and restoration without changing the fit."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app
from mfrm_app import basic_analysis as basic, engine, project, reliability
from mfrm_app import visual_maps


@pytest.fixture(scope="module")
def analysis():
    data = pd.DataFrame([
        dict(Person=f"P{i}", Rater=f"R{r}", Task=f"T{t}",
             Score=0 if i == 0 else 2 if i == 7 else (i + r + t) % 3)
        for i in range(8) for r in range(2) for t in range(2)
    ])
    result = basic.run_basic_analysis(data, person_col="Person", facet_cols=["Rater", "Task"],
                                     score_col="Score", rating_min=0, rating_max=2)
    assert result["result"]["opt"].success
    assert not basic.qualification_snapshot(result["result"])["calibration_hold_reason"]
    return result


@pytest.mark.parametrize("legacy", [False, True])
def test_person_hold_is_shared_and_keeps_raw_values(analysis, legacy):
    result = analysis["result"]
    raw = analysis["diagnostics"]["measures"].copy(deep=True)
    if legacy:
        raw = raw.drop(columns=[c for c in ("ExtremeScoreDirection", "FiniteJMLEEstimate",
            "PersonInferenceReady", "ReportableEstimate") if c in raw])
        raw["SE_Status"] = raw["CI_Status"] = "conditional_approximation"
    before = raw.copy(deep=True)
    qualified = engine.qualify_measure_uncertainty(raw, result)
    boundary = qualified.Facet.eq("Person") & qualified.Level.isin(["P0", "P7"])
    assert qualified.loc[boundary, "SE_Status"].eq("person_inference_withheld").all()
    assert qualified.loc[boundary, "CI_Status"].eq("person_inference_withheld").all()
    assert qualified.loc[~boundary, "SE_Status"].eq("conditional_approximation").all()
    pd.testing.assert_frame_equal(qualified[["Estimate", "SE", "CI_Lower", "CI_Upper"]],
                                  before[["Estimate", "SE", "CI_Lower", "CI_Upper"]])
    pd.testing.assert_frame_equal(raw, before)
    pd.testing.assert_frame_equal(engine.qualify_measure_uncertainty(qualified, result), qualified)


def assert_public_table(table):
    indexed = table.set_index(["Facet", "Level"])
    for person in ("P0", "P7"):
        row = indexed.loc[("Person", person)]
        assert row[["Estimate", "SE / SD", "Lower", "Upper", "Infit", "Outfit"]].eq("").all()
        assert row.Basis == "No finite JMLE"
    assert indexed.loc[("Person", "P1"), "Estimate"] != ""
    assert indexed.loc[("Rater", "R0"), "SE / SD"] != ""


def test_publication_and_desktop_tables_match_person_qualification(analysis):
    table, note = app._publication_measure_table(analysis["result"], analysis["diagnostics"])
    assert_public_table(table)
    assert "No finite JMLE" in note
    desktop = next(frame for title, frame, note in basic.current_tables(analysis) if title == "Measures and Fit")
    pd.testing.assert_frame_equal(table, desktop)
    _, compact, caption = next(item for item in app._publication_report_tables(
        analysis["result"], analysis["diagnostics"], {}, include_details=True)
        if item[0] == "Element Estimates and Fit")
    assert "Interval level" not in compact
    assert "Interval level = 0.95 for displayed intervals" in caption
    assert compact.set_index(["Facet", "Level"]).loc[("Person", "P0"), "Interval"] == ""


def test_person_reliability_is_held_without_dropping_people_or_valid_facets(analysis):
    measures = analysis["diagnostics"]["measures"]
    calculated = reliability.calculate_reliability(measures).set_index("Facet")
    person = calculated.loc["Person"]
    assert person.Levels == 8 and person.ValidLevels == 6
    assert person.ReliabilityStatus == "withheld"
    assert person.ReliabilityReasonCode == "stat.reliability.person_measure_withheld"
    assert calculated.loc[["Rater", "Task"], "ReliabilityStatus"].eq("diagnostic").all()
    # Reopen a legacy calculation: retain its raw numbers, but reapply today's hold.
    legacy_measures = measures.drop(columns=[c for c in ("FiniteJMLEEstimate", "PersonInferenceReady") if c in measures]).copy()
    legacy_measures["SE_Status"] = legacy_measures["CI_Status"] = "conditional_approximation"
    saved = reliability.calculate_reliability(legacy_measures)
    saved["ReliabilityCalculationVersion"] = "paired-reliability-v1"
    original = saved.copy(deep=True)
    diagnostics = dict(analysis["diagnostics"], reliability=saved)
    qualified = app.qualified_reliability_table(analysis["result"], diagnostics)
    assert qualified.set_index("Facet").loc["Person", "ReliabilityStatus"] == "withheld"
    pd.testing.assert_frame_equal(qualified[["RMSE", "Reliability", "Separation", "Strata"]],
                                  original[["RMSE", "Reliability", "Separation", "Strata"]])
    public, _ = app._publication_reliability_table(analysis["result"], diagnostics)
    public = public.set_index("Facet")
    assert public.loc["Person", ["RMSE", "Reliability", "Separation", "Strata"]].eq("").all()
    assert public.loc[["Rater", "Task"], "Use"].eq("Diagnostic").all()
    pd.testing.assert_frame_equal(saved, original)
    without_measures = reliability.qualify_saved_reliability(saved)
    assert without_measures.ReliabilityStatus.eq("unverified").all()


def test_missing_interior_estimate_is_not_mislabeled_as_an_endpoint(analysis):
    result = deepcopy(analysis["result"])
    persons = result["person_boundary"]["persons"]
    missing = persons.Person.eq("P1")
    persons.loc[missing, ["FiniteJMLEEstimate", "PersonInferenceReady"]] = False
    persons.loc[missing, "ReportableEstimate"] = np.nan
    table, _ = app._publication_measure_table(result, analysis["diagnostics"])
    row = table.set_index(["Facet", "Level"]).loc[("Person", "P1")]
    assert row.Basis == "Withheld"
    assert row[["Estimate", "SE / SD", "Lower", "Upper"]].eq("").all()


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_first_read_explains_person_hold_without_hiding_valid_facets(analysis, monkeypatch, lang):
    monkeypatch.setitem(app.st.session_state, "lang", lang)
    rows = app.build_first_read_guide_rows(analysis["result"], analysis["diagnostics"], localize=True)
    row = next(row for row in rows if row["Check"] == "3. Reliability / separation")
    assert row["Status"] != "OK"
    assert app.t("first_read.reliability_reasons.person_measure_withheld") in row["What it means"]
    assert "Rater" in row["What it means"] and "Task" in row["What it means"]


def test_saved_project_requalifies_without_refitting(analysis, tmp_path, monkeypatch):
    original = deepcopy(analysis)
    path = tmp_path / "境界得点.mfrm"
    project.save_project(path, analysis)
    def forbidden(*args, **kwargs):
        pytest.fail("Reading or reporting must not refit or recompute diagnostics")
    monkeypatch.setattr(engine, "mfrm_estimate", forbidden)
    monkeypatch.setattr(engine, "mfrm_diagnostics", forbidden)
    restored = project.load_project(path)
    assert restored["analysis_id"] == original["analysis_id"]
    assert restored["content_fingerprint"] == original["content_fingerprint"]
    for current in (analysis, restored):
        assert_public_table(basic.current_tables(current)[1][1])
        assert "No finite JMLE" in basic.basic_report_html(current)
    np.testing.assert_array_equal(restored["result"]["opt"].x, original["result"]["opt"].x)


def test_technical_figure_positions_carry_person_caution(analysis, monkeypatch):
    result = analysis["result"]
    assert "not finite JMLEs" in visual_maps.map_caution(result)
    assert "有限のJMLEではなく" in visual_maps.map_caution(result, lang="ja")
    normal = dict(result, person_boundary={"status": "no_extreme_persons"})
    assert visual_maps.map_caution(normal) == ""
    monkeypatch.setattr(app, "_plotly_or_matplotlib_png", lambda *args: b"image")
    figures = app._publication_figure_payloads(result, analysis["diagnostics"], include_details=True)
    captions = {name: caption for name, caption, _ in figures}
    for name in ("wright_map", "fit_scatter", "facet_distribution"):
        assert "not finite JMLEs" in captions[name]
        assert "technical review only" in captions[name]
