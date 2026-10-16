"""Real fits retain uncertainty meaning across tables and publication formats."""
from io import BytesIO, StringIO

import numpy as np
import pandas as pd
import pytest
from docx import Document

import streamlit_app as app
from tests.test_jmle_facet_boundary import separated_data


@pytest.fixture(scope="module", params=["regular", "facet_boundary", "fixed_sd", "free_sd"])
def fitted_case(request):
    case = request.param
    data = separated_data() if case == "facet_boundary" else app._make_self_test_gradient_data()
    options = ({"method": "MML", "mml_engine": "EM", "quad_points": 15,
                "estimate_population_sd": case == "free_sd"}
               if case in {"fixed_sd", "free_sd"} else {"method": "JMLE"})
    result = app.mfrm_estimate(data, "Person", ["Rater", "Task"], "Score", maxit=100, **options)
    diagnostics = app.mfrm_diagnostics(result, compute_pca=False, interaction_pairs=[])
    return case, result, diagnostics


def test_real_measures_preserve_values_and_current_holds_in_csv_and_scaled_tables(fitted_case):
    case, result, diagnostics = fitted_case
    held = case in {"facet_boundary", "free_sd"}
    original = diagnostics["measures"].copy(deep=True)
    frames = app.build_result_bundle_frames(result, diagnostics)
    downloads, _ = app.collect_download_frames(result, diagnostics, {}, pd.DataFrame(), pd.DataFrame(),
                                              public_export_mode=False)
    apa = app._collect_apa_exportable_tables(result, diagnostics)["Combined measures and fit"]
    for table in (frames["measures"], downloads["measures"], apa):
        csv = pd.read_csv(StringIO(table.to_csv(index=False)))
        np.testing.assert_allclose(csv.Estimate, original.Estimate, equal_nan=True)
        np.testing.assert_allclose(csv.SE, original.SE, equal_nan=True)
        assert csv.SE_Status.eq("calibration_inference_withheld").all() == held
        if held:
            assert csv.CI_Status.eq("calibration_inference_withheld").all()
            assert csv.UncertaintyCaution.str.contains("inference is withheld").all()
    scaled = app.calc_facets_report_tbls(result, diagnostics, uscale=2.)
    for facet, table in scaled.items():
        source = original[original.Facet.eq(facet)].set_index("Level").reindex(table.Level)
        np.testing.assert_allclose(table.ModelSE, 2*source.SE, atol=.005, equal_nan=True)
        assert table.SE_Status.tolist() == source.SE_Status.tolist()
        assert table.SE_Method.tolist() == source.SE_Method.tolist()
    pd.testing.assert_frame_equal(diagnostics["measures"], original)


def test_real_publication_formats_keep_basis_and_suppress_held_intervals(fitted_case, monkeypatch):
    case, result, diagnostics = fitted_case
    held = case in {"facet_boundary", "free_sd"}
    # Figure rendering has its own tests; exercise the actual document builders and tables.
    monkeypatch.setattr(app, "_publication_figure_payloads", lambda *args: [])
    table, caption = app._publication_measure_table(result, diagnostics)
    assert len(table.columns) <= 10
    assert table.Basis.eq("Withheld").all() == held
    if held:
        assert table[["SE / SD", "Lower", "Upper"]].eq("").all().all()
        assert "technical terminal-point" in caption
    elif case == "fixed_sd":
        assert table.loc[table.Facet.eq("Person"), "Basis"].eq("Posterior SD").all()
        assert "not the sampling SE of EAP" in caption
    else:
        assert table.Basis.eq("Conditional").all()
    reliability, reliability_caption = app._publication_reliability_table(result, diagnostics)
    assert reliability.Use.eq("Withheld").all() == held
    if held:
        assert reliability[["RMSE", "Separation", "Strata", "Reliability"]].eq("").all().all()
        assert "Interpretation is withheld" in reliability_caption
    word = Document(BytesIO(app.build_publication_word_bytes(result, diagnostics)))
    assert [cell.text for cell in word.tables[0].rows[0].cells] == table.columns.tolist()
    basis_col = table.columns.get_loc("Basis")
    assert [row.cells[basis_col].text for row in word.tables[0].rows[1:]] == table.Basis.tolist()
    assert any(caption in paragraph.text for paragraph in word.paragraphs)
    assert any(reliability_caption in paragraph.text for paragraph in word.paragraphs)
    assert [cell.text for cell in word.tables[1].rows[0].cells] == reliability.columns.tolist()
    html = app.build_publication_html_bytes(result, diagnostics).decode()
    assert app._html_escape(caption) in html
    assert table.to_html(index=False, border=0) in html
    assert reliability.to_html(index=False, border=0) in html
    pdf_tables = []
    original_pdf_table = app._pdf_dataframe_flowable
    def record_table(frame, caption=None):
        pdf_tables.append((frame, caption))
        return original_pdf_table(frame, caption)
    monkeypatch.setattr(app, "_pdf_dataframe_flowable", record_table)
    pdf = app.build_publication_pdf_bytes(result, diagnostics)
    assert pdf.startswith(b"%PDF")
    pd.testing.assert_frame_equal(pdf_tables[0][0], table)
    pd.testing.assert_frame_equal(pdf_tables[1][0], reliability)
    assert caption in pdf_tables[0][1]


@pytest.mark.parametrize("config", [
    {"method": "MML", "estimate_population_sd": True},
    {"method": "JMLE", "jmle_facet_boundary_review": True},
])
def test_saved_positive_flags_do_not_override_current_hold_or_modify_source(config):
    source = pd.DataFrame({"Facet": ["Person"], "Level": ["P1"], "Estimate": [0.], "SE": [.2],
                           "CI_Lower": [-.392], "CI_Upper": [.392], "SE_Status": ["ok"], "CI_Status": ["ok"]})
    before = source.copy(deep=True)
    qualified = app.qualify_measure_uncertainty(source, {"config": config})
    assert qualified.SE_Status.eq("calibration_inference_withheld").all()
    assert qualified.CI_Status.eq("calibration_inference_withheld").all()
    assert qualified.SE.iloc[0] == .2
    table, _ = app._publication_measure_table({"config": config}, {"measures": source})
    assert table.Basis.iloc[0] == "Withheld"
    assert table["SE / SD"].iloc[0] == ""
    pd.testing.assert_frame_equal(source, before)


@pytest.mark.parametrize("lang", ["en", "ja"])
@pytest.mark.parametrize("case", ["facet_boundary", "free_sd"])
def test_actual_report_table_displays_uncertainty_statuses_and_hold(lang, case):
    from streamlit.testing.v1 import AppTest

    def view(lang, case):
        from contextlib import ExitStack
        from unittest.mock import patch
        import streamlit as st
        import streamlit_app as app
        from tests.test_jmle_facet_boundary import separated_data

        st.session_state["lang"] = lang
        data = separated_data() if case == "facet_boundary" else app._make_self_test_gradient_data()
        options = ({"method": "MML", "mml_engine": "EM", "quad_points": 15,
                    "estimate_population_sd": True} if case == "free_sd" else {"method": "JMLE"})
        result = app.mfrm_estimate(data, "Person", ["Rater", "Task"], "Score", maxit=100, **options)
        diagnostics = app.mfrm_diagnostics(result, compute_pca=False, interaction_pairs=[])
        # Optional planning and export widgets are independent of the report table.
        with ExitStack() as stack:
            for name in ("render_estimand_contract_panel", "render_model_choice_guidance",
                         "render_g_d_study_section", "render_design_network_section",
                         "render_rater_severity_section", "render_rater_halo_section",
                         "render_parameter_recovery_simulation", "render_apa_presets_section"):
                stack.enter_context(patch.object(app, name, return_value=None))
            app._render_report_tables(result, diagnostics)
        st.session_state["hold_message"] = app.t("guided.export_state_hold")

    at = AppTest.from_function(view, args=(lang, case)).run(timeout=60)
    assert not at.exception
    assert at.session_state["hold_message"] in [warning.value for warning in at.warning]
    table = next(frame.value for frame in at.dataframe if "CI_Status" in frame.value.columns)
    assert {"SE_Method", "SE_Status", "CI_Method", "UncertaintyCaution"}.issubset(table.columns)
    assert table.SE_Status.eq("calibration_inference_withheld").all()
    assert table.CI_Status.eq("calibration_inference_withheld").all()
