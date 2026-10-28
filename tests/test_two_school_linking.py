"""Teaching design and linking boundaries; not a statistical recovery study."""
from io import BytesIO

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_analysis_journeys import _run
from tests.test_help_contextual_app_integration import _block_analysis_calls, _locale_value


def _data():
    return app.sample_mfrm_data_by_key("two_school_linking")


def _prep(data):
    return app.prepare_mfrm_data(
        data, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Criterion"], rating_min=0, rating_max=4,
        keep_original=True,
    )


def _rank(data):
    prep = _prep(data)
    specs = app.prepare_constraint_specs(prep, noncenter_facet="Person")
    config = {
        "model": "RSM", "method": "JMLE", "n_cat": 5,
        "n_person": len(prep["levels"]["Person"]),
        "facet_names": prep["facet_names"],
        "facet_levels": {f: prep["levels"][f] for f in prep["facet_names"]},
        "facet_signs": {"Rater": -1, "Criterion": -1},
        "theta_spec": specs["theta_spec"], "facet_specs": specs["facet_specs"],
        "population_model": {"enabled": False, "n_params": 0},
    }
    return app.build_jmle_eta_identifiability_audit(prep, config)


def test_assignment_uses_same_reference_ids_and_local_school_raters():
    data = _data()
    assert len(data) == 360
    assert not data.duplicated(["Person", "Rater", "Task", "Criterion"]).any()
    assert data.groupby(data.Person.str[0]).size().to_dict() == {"A": 120, "B": 120, "C": 120}
    for person, rows in data.groupby("Person"):
        expected = ({"A_R1", "A_R2", "B_R1", "B_R2"} if person.startswith("C")
                    else {f"{person[0]}_R1", f"{person[0]}_R2"})
        assert set(rows.Rater) == expected
        assert rows.groupby("Rater").size().eq(3).all()


@pytest.mark.parametrize("case,components,nullity", [
    ("shared", 1, 0), ("removed", 2, 1), ("relabeled", 2, 1),
])
def test_only_genuinely_shared_recording_ids_bridge_the_schools(case, components, nullity):
    data = _data()
    common = data.Person.str.startswith("C")
    if case == "removed":
        data = data.loc[~common].copy()
    elif case == "relabeled":
        # Different IDs for the same recording destroy the modeled overlap.
        data.loc[common, "Person"] += "_" + data.loc[common, "Rater"].str[0]
    audit = _rank(data)
    assert bool(audit["eta_structurally_identified"]) == (case == "shared")
    assert audit["summary"].iloc[0].EtaStructuralNullity == nullity
    assert audit["connectivity"].set_index("Facet").loc["Rater", "Components"] == components


@pytest.mark.parametrize("constraint", ["none", "fixed", "group", "dummy"])
def test_anchor_counts_never_certify_valid_linking(constraint):
    kwargs = {}
    if constraint == "fixed":
        kwargs["anchor_df"] = pd.DataFrame({
            "Facet": ["Person"] * 10, "Level": [f"C{i:02d}" for i in range(1, 11)],
            "Anchor": [0.0] * 10,
        })
    elif constraint == "group":
        kwargs["group_anchor_df"] = pd.DataFrame({
            "Facet": ["Person"] * 10, "Level": [f"C{i:02d}" for i in range(1, 11)],
            "Group": ["Reference"] * 10, "GroupValue": [0.0] * 10,
        })
    elif constraint == "dummy":
        kwargs["dummy_facets"] = ["Rater"]
    audit = app.audit_mfrm_anchors(_prep(_data()), **kwargs)
    result = {"config": {"anchor_audit": audit}}
    result["equating_chain"] = app.build_equating_chain_summary(result)
    workflow = app.build_anchor_equating_workflow_plan(result)
    assert "Not linked" not in workflow.Status.tolist()
    assert workflow.set_index("Step").loc[2, "Status"] != "Ready"
    assert workflow.set_index("Step").loc[4, "Status"] != "Ready"
    assert "Linked" not in audit["summary"].Status.tolist()
    if constraint == "none":
        assert not result["equating_chain"]["available"]
        assert workflow.set_index("Step").loc[2, "Status"] == "Not applicable"
    else:
        public = app.build_equating_chain_readiness_summary(result).iloc[0]
        assert public.Status == "Review"
        assert "not validated links" in public.Boundary
        if constraint == "dummy":
            assert public.LinkedEdges == 0


def test_legacy_linked_label_still_requires_reference_review():
    result = {"equating_chain": {"available": True, "edges": pd.DataFrame({
        "Strength": ["Linked"], "CommonLinkedLevels": [10],
    })}}
    row = app.build_equating_chain_readiness_summary(result).iloc[0]
    assert row.LinkedEdges == 1
    assert row.Status == "Review"


def test_sample_selection_fit_lesson_and_language_preserve_results(monkeypatch):
    downloads = {}
    original_download = app.st.download_button

    def capture(label, data, *args, **kwargs):
        downloads[kwargs.get("file_name")] = data
        return original_download(label, data, *args, **kwargs)

    monkeypatch.setattr(app.st, "download_button", capture)
    at = AppTest.from_string("import streamlit_app as app\napp.main()", default_timeout=120)
    _run(at)
    at.button(key="onboarding_dismiss").click()
    _run(at)
    at.selectbox(key="data_source_class").set_value("sample")
    _run(at)
    at.selectbox(key="data_source_scenario").set_value("two_school_linking")
    at.radio(key="lang").set_value("ja")
    _run(at)
    pd.testing.assert_frame_equal(
        pd.read_csv(BytesIO(downloads["mfrm_sample_two_school_linking.csv"])), _data(),
    )
    group_file = pd.read_csv(BytesIO(downloads["two_school_person_groups.csv"]), dtype={"Person": str})
    assert len(group_file) == 50 and group_file.Person.is_unique
    assert group_file.School.value_counts().to_dict() == {"A": 20, "B": 20}
    assert group_file.loc[group_file.Person.str.startswith("C"), "School"].isna().all()
    assert at.radio(key="facets_mode_workflow_mode").value == "Guided defaults"
    pair_key = next(w.key for w in at.selectbox if w.key and w.key.startswith("facets_mode_bias_pair_"))
    assert at.selectbox(key=pair_key).value == "Rater x Criterion"
    at.button(key="facets_mode_run_primary").click()
    _run(at)
    output = at.session_state["facets_mode_output"]
    result = output["result"]
    assert result["config"]["method"] == "JMLE"
    assert result["config"]["model"] == "RSM"
    assert set(result["config"]["facet_names"]) == {"Rater", "Criterion"}
    assert bool(result["summary"].iloc[0].Converged)
    assert result["identifiability"]["status"] == "eta_identified"
    assert result["config"]["anchor_audit"]["summary"].AnchoredLevelsTotal.eq(0).all()
    original_measures = output["diagnostics"]["measures"].copy(deep=True)
    analysis_id = app.build_result_analysis_identity(result).analysis_id
    calls = _block_analysis_calls(monkeypatch)
    for lang in ("ja", "en"):
        at.radio(key="lang").set_value(lang)
        _run(at)
        title = _locale_value(lang, "sample_scenarios.two_school_linking.guide_title")
        guide = _locale_value(lang, "sample_scenarios.two_school_linking.guide")
        assert title in [e.label for e in at.expander]
        assert guide in [m.value for m in at.markdown]
        assert _locale_value(lang, "sample_scenarios.two_school_linking.assignment_next_step") in [e.value for e in at.info]
        at.button(key="two_school_open_dff").click()
        _run(at)
        assert at.session_state["guided_essential_section"] == "diagnostics"
        assert at.session_state["guided_diagnostics_panel"] == "bias_interaction"
        assert at.selectbox(key="bias_pair_selector").value == "Rater x Criterion"
        assert _locale_value(lang, "bias_interaction.reading_guide") in [m.value for m in at.markdown]
        details = next(e for e in at.expander if e.label == _locale_value(lang, "bias_interaction.technical_details_label"))
        assert not details.proto.expanded
        at.radio(key="bias_heatmap_metric").set_value("Standardized Bias")
        _run(at)
        at.button(key="mfrm_help_open_global").click()
        _run(at)
        at.selectbox(key=app._HELP_TOPIC_PICKER_KEY).set_value("help.results.differential_interaction")
        _run(at)
        help_details = next(e for e in at.expander if e.label == _locale_value(lang, "help_nav.details_heading"))
        assert not help_details.proto.expanded
        body = _locale_value(lang, "help_topics.results.differential_interaction.sections.external_group_design.body")
        assert body in [e.value for e in help_details.markdown]
        assert not any(code in body for code in ("model_identified", "partially_identified", "target_unidentified", "G3"))
        at.button(key="mfrm_help_return_to_source").click()
        _run(at)
        assert any(w.key == "bias_heatmap_metric" for w in at.radio), "\n".join(
            [w.value for w in at.warning] + [
                f"{w.key}={w.value}" for w in [*at.selectbox, *at.radio]
                if any(word in str(w.key) for word in ("guided", "panel", "workflow", "bias", "density"))
            ]
        )
        assert at.radio(key="bias_heatmap_metric").value == "Standardized Bias"
        assert at.selectbox(key="bias_pair_selector").value == "Rater x Criterion"
        assert at.session_state["facets_mode_output"] is output
        assert app.build_result_analysis_identity(result).analysis_id == analysis_id
        pd.testing.assert_frame_equal(output["diagnostics"]["measures"], original_measures, check_exact=True)
    assert calls == {"estimate": 0, "refit": 0}
    # A teaching default must not overwrite a deliberate choice on rerender.
    at.selectbox(key=pair_key).set_value("Person x Rater")
    _run(at)
    at.radio(key="lang").set_value("ja")
    _run(at)
    assert at.selectbox(key=pair_key).value == "Person x Rater"
    assert calls == {"estimate": 0, "refit": 0}
    # A different input must not inherit this sample's interpretation guide.
    at.selectbox(key="data_source_class").set_value("paste")
    _run(at)
    assert not any(e.label == _locale_value("ja", "sample_scenarios.two_school_linking.guide_title") for e in at.expander)
    assert not any(b.key == "two_school_open_dff" for b in at.button)
