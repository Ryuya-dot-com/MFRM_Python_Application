"""Observed assignment counts and explicit Person-group joins, not DIF inference."""
from io import BytesIO

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from mfrm_app.design_assignment import build_group_assignment_coverage, validate_person_groups
from tests.test_analysis_journeys import _run
from tests.test_help_contextual_app_integration import _block_analysis_calls, _locale_value


def _school_data():
    data = app.sample_mfrm_data_by_key("two_school_linking")
    groups = app._two_school_sample_groups(data).set_index("Person").School
    return data, groups


def test_school_coverage_excludes_reference_recordings_and_counts_people_once():
    data, groups = _school_data()
    coverage = build_group_assignment_coverage(data, groups, facet="Rater")
    assert coverage["counts"].to_numpy().tolist() == [[20, 0], [20, 0], [0, 20], [0, 20]]
    assert list(coverage["counts"].index) == ["A_R1", "A_R2", "B_R1", "B_R2"]
    assert coverage["group_sizes"].to_dict() == {"A": 20, "B": 20}
    assert coverage["unlabeled_persons"] == 10
    assert coverage["levels_in_all_groups"] == 0
    changed = pd.concat([data, data]).assign(Score=999, Weight=-10)
    other = build_group_assignment_coverage(changed, groups, facet="Rater")
    pd.testing.assert_frame_equal(coverage["counts"], other["counts"])


def test_observed_cross_group_assignments_change_counts_without_implying_identification():
    data, groups = _school_data()
    added = pd.DataFrame({"Person": ["A01", "A01", "B01", "B01"],
                          "Rater": ["B_R1", "B_R2", "A_R1", "A_R2"]})
    coverage = build_group_assignment_coverage(pd.concat([data, added]), groups, facet="Rater")
    assert coverage["counts"].to_numpy().tolist() == [[20, 1], [20, 1], [1, 20], [1, 20]]
    assert coverage["levels_in_all_groups"] == 4
    assert not {"identified", "ready", "PublicConclusionAllowed"} & coverage.keys()


@pytest.mark.parametrize("index", [["P1", "P1"], [1, "1"], ["P1", None], ["P1", " "]])
def test_group_ids_must_be_unique_and_nonblank(index):
    with pytest.raises(ValueError, match="Person"):
        validate_person_groups(pd.Series(["A", "B"], index=index))


def test_missing_groups_stay_missing_and_ids_are_not_guessed():
    labels = validate_person_groups(pd.Series([" A ", " \t ", None, np.nan, "B"], index=["001", "2", "3", "4", " 5 "]))
    assert labels["001"] == "A"
    assert labels.loc[["2", "3", "4"]].isna().all()
    assert " 5 " in labels.index and "5" not in labels.index
    data = pd.DataFrame({"Person": ["001", "2", "3", "4", "5"], "Rater": ["R1"] * 5})
    coverage = build_group_assignment_coverage(data, labels, facet="Rater")
    assert coverage["group_sizes"].to_dict() == {"A": 1}
    assert coverage["unlabeled_persons"] == 4 and coverage["extra_persons"] == 1


def test_all_missing_groups_and_unusual_column_names_are_supported():
    data = pd.DataFrame({"Person": ["1", "2"], "_Group": ["r1", "r2"]})
    empty = build_group_assignment_coverage(data, {"1": None, "2": ""}, facet="_Group")
    assert empty["counts"].shape == (2, 0)
    assert empty["unlabeled_persons"] == 2
    grouped = build_group_assignment_coverage(data, {"1": "FacetLevel", "2": "B"}, facet="_Group")
    exported = grouped["counts"].stack().rename("Persons").reset_index()
    assert list(exported.columns) == ["FacetLevel", "Group", "Persons"]
    assert len(exported) == 4


def test_large_group_matrix_is_not_allocated(monkeypatch):
    size = 600
    data = pd.DataFrame({"Person": [f"P{i}" for i in range(size)], "Rater": [f"R{i}" for i in range(size)]})
    labels = pd.Series([f"G{i}" for i in range(size)], index=data.Person)
    def forbidden(*args, **kwargs):
        pytest.fail("A large dense group matrix was allocated")
    monkeypatch.setattr(pd.Series, "unstack", forbidden)
    result = build_group_assignment_coverage(data, labels, facet="Rater")
    assert result["counts"] is None
    assert len(result["long_counts"]) == size
    assert result["zero_cells"] == size * size - size
    assert result["levels_in_all_groups"] == 0


def test_classical_dif_api_drops_blank_groups_and_rejects_duplicate_ids():
    data, groups = _school_data()
    result = app.analyze_classical_dif(data, "Criterion", groups.fillna(""),
                                      min_group_n=5, fit_logistic=False, purify=False)
    assert set(result["settings"]["groups"]) == {"A", "B"}
    assert len(result["analysis_frame"]["frame"]) == 40
    assert not result["analysis_frame"]["frame"].Person.str.startswith("C").any()
    duplicate = pd.concat([groups, groups.iloc[:1]])
    assert "duplicate" in app.analyze_classical_dif(data, "Criterion", duplicate)["_skip_reason"]


@pytest.mark.parametrize("language", ["ja", "en"])
@pytest.mark.parametrize("file_format", ["csv", "xlsx"])
def test_group_upload_keeps_leading_zero_ids_and_waits_for_explicit_run(language, file_format, monkeypatch):
    data = pd.DataFrame({"Person": ["001", "002", "003", "004"],
                         "Rater": ["R1", "R1", "R2", "R2"],
                         "Criterion": ["C1", "C2", "C1", "C2"], "Score": [0, 1, 2, 3]})
    groups = pd.DataFrame({"Person": ["001", "002", "003", "004"], "School": ["A", "", "B", "B"]})
    upload = BytesIO()
    if file_format == "csv":
        upload.write(groups.to_csv(index=False).encode())
    else:
        groups.to_excel(upload, index=False)
    upload.name = "groups." + file_format
    def get_upload(*args, **kwargs):
        upload.seek(0)
        return upload
    monkeypatch.setattr(app.st, "file_uploader", get_upload)
    calls = _block_analysis_calls(monkeypatch)
    dif_calls = []
    def screen(prep, facet, labels, **kwargs):
        dif_calls.append(labels.copy())
        return {"_skip_reason": "Test stops before statistical calculation."}
    monkeypatch.setattr(app, "analyze_classical_dif", screen)
    at = AppTest.from_string('import streamlit as st\nimport streamlit_app as app\napp.show_classical_dif_section(st.session_state["result"], {})')
    at.session_state["lang"] = language
    at.session_state["result"] = {"prep": {"data": data}, "config": {"facet_names": ["Rater", "Criterion"]}}
    _run(at)
    assert not dif_calls
    assert at.dataframe[0].value.to_numpy().tolist() == [[1, 0], [0, 2]]
    assert at.dataframe[0].value.index.name == "Rater"
    assert at.dataframe[0].value.columns.name == "School"
    assert _locale_value(language, "dif_screening.assignment_missing_warning") in [w.value for w in at.warning]
    assert _locale_value(language, "dif_screening.assignment_boundary") in [w.value for w in at.info]
    at.button(key="dif_run").click()
    _run(at)
    assert len(dif_calls) == 1
    assert dif_calls[0]["001"] == "A" and pd.isna(dif_calls[0]["002"])
    assert calls == {"estimate": 0, "refit": 0}


@pytest.mark.parametrize("language", ["ja", "en"])
@pytest.mark.parametrize("ids", [["001", "001"], ["001", ""]])
def test_invalid_group_file_stops_before_screening(language, ids, monkeypatch):
    upload = BytesIO(pd.DataFrame({"Person": ids, "School": ["A", "B"]}).to_csv(index=False).encode())
    upload.name = "groups.csv"
    monkeypatch.setattr(app.st, "file_uploader", lambda *args, **kwargs: upload)
    def forbidden(*args, **kwargs):
        pytest.fail("Invalid group data reached screening")
    monkeypatch.setattr(app, "analyze_classical_dif", forbidden)
    at = AppTest.from_string('import streamlit as st\nimport streamlit_app as app\napp.show_classical_dif_section(st.session_state["result"], {})')
    at.session_state["lang"] = language
    at.session_state["result"] = {"prep": {"data": pd.DataFrame({"Person": ["001"], "Rater": ["R1"]})},
                                   "config": {"facet_names": ["Rater"]}}
    _run(at)
    assert _locale_value(language, "dif_screening.group_invalid_error") in [e.value for e in at.error]
    assert not at.button


def test_large_coverage_view_keeps_observed_counts_downloadable(monkeypatch):
    size = 300
    data = pd.DataFrame({"Person": [f"P{i}" for i in range(size)],
                         "Rater": [f"R{i}" for i in range(size)], "Score": 1})
    groups = pd.DataFrame({"Person": data.Person, "Group": [f"G{i}" for i in range(size)]})
    downloads = {}
    original = app.st.download_button
    def capture(label, data, *args, **kwargs):
        downloads[kwargs.get("file_name")] = data
        return original(label, data, *args, **kwargs)
    monkeypatch.setattr(app.st, "download_button", capture)
    def forbidden(*args, **kwargs):
        pytest.fail("Coverage review ran a DIF screen")
    monkeypatch.setattr(app, "analyze_classical_dif", forbidden)
    at = AppTest.from_string('import streamlit as st\nimport streamlit_app as app\napp.show_classical_dif_section(st.session_state["result"], {})')
    at.session_state["lang"] = "ja"
    at.session_state["result"] = {"prep": {"data": data}, "config": {"facet_names": ["Rater"]},
                                   "population": {"person_data": groups}}
    _run(at)
    assert _locale_value("ja", "dif_screening.assignment_large_info") in [i.value for i in at.info]
    assert not at.dataframe
    exported = pd.read_csv(BytesIO(downloads["dif_group_assignment_counts.csv"]))
    assert len(exported) == size and exported.Persons.eq(1).all()
    assert exported.Facet.eq("Rater").all() and exported.GroupVariable.eq("Group").all()


@pytest.mark.parametrize("facet", ["Missing", " Rater "])
def test_stale_or_space_padded_facet_names_do_not_crash_coverage(facet):
    data = pd.DataFrame({"Person": ["001", "002"], " Rater ": ["R1", "R2"], "Score": [0, 1]})
    groups = pd.DataFrame({"Person": ["001", "002"], "School": ["A", "B"]})
    at = AppTest.from_string('import streamlit as st\nimport streamlit_app as app\napp.show_classical_dif_section(st.session_state["result"], {})')
    at.session_state["lang"] = "en"
    at.session_state["result"] = {"prep": {"data": data}, "config": {"facet_names": [facet]},
                                   "population": {"person_data": groups}}
    _run(at)
    if facet == "Missing":
        assert _locale_value("en", "dif_screening.no_facets_info") in [i.value for i in at.info]
    else:
        assert at.selectbox(key="dif_assignment_facet").value == facet
        assert at.dataframe[0].value.to_numpy().tolist() == [[1, 0], [0, 1]]
