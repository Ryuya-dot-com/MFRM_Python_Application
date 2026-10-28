"""Definitions remain near decisions without changing values or hold states."""
import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app


@pytest.mark.parametrize('eap, variance, expected_g', [(False, 5., 2.), (False, .5, 0.), (True, 4., 2.)])
def test_strata_examples_use_actual_se_and_eap_formulas(eap, variance, expected_g):
    a = np.sqrt(variance)
    measures = pd.DataFrame({'Facet': ['Person']*3, 'Estimate': [-a, 0., a],
        'SE': [1.]*3, 'Infit': [1.]*3, 'Outfit': [1.]*3,
        'SE_Status': ['posterior_eap_sd' if eap else 'conditional_approximation']*3})
    before = measures.copy(deep=True)
    row = app.calc_reliability(measures).iloc[0]
    assert row.Separation == pytest.approx(expected_g)
    assert row.Strata == pytest.approx((4*expected_g+1)/3)
    assert row.Reliability == pytest.approx(expected_g**2/(1+expected_g**2))
    pd.testing.assert_frame_equal(measures, before)


@pytest.mark.parametrize('lang', ['en', 'ja'])
def test_visible_definitions_and_collapsed_details_match_unrounded_flags(lang):
    def view(lang):
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        st.session_state['lang'] = lang
        app.render_metric_reading_notes('fit')
        app.render_metric_reading_notes('reliability')
        # Definition >1 does not imply a >1.5 review flag; round only after decisions.
        app._show_misfit_flags(pd.DataFrame({'Infit': [.5, 1.2, 1.5, 1.50001, 2., 2.00001, .49999, float('nan')],
                                           'Outfit': [1.]*7 + [float('nan')]}))
    at = AppTest.from_function(view, args=(lang,)).run()
    assert not at.exception
    assert len(at.expander) == 2 and all(not e.proto.expanded for e in at.expander)
    text = '\n'.join(x.value for x in at.caption)
    assert 'Strata' in text and ('Underfit' in text or 'underfit' in text)
    assert '(4G+1)/3' in text
    details = '\n'.join(x.value for x in at.markdown)
    assert '1.2' in details and 'EAP' in details and '1/3' in details
    assert 'https://www.winsteps.com/facetman64/reliability.htm' in details
    assert [m.value for m in at.metric] == ['3/7', '1', '2', '1']
    if lang == 'ja':
        assert '能力が低いという意味ではありません' in details
        assert '実際に3群が存在する証拠でもありません' in details


@pytest.mark.parametrize('free_sd', [False, True])
def test_report_view_reuses_qualified_eap_values_and_definitions(monkeypatch, free_sd):
    # Isolate the existing Tables surface from its independent optional analyses.
    for name in ('render_estimand_contract_panel', 'render_model_choice_guidance', 'render_g_d_study_section',
                 'render_design_network_section', 'render_rater_severity_section', 'render_rater_halo_section',
                 'render_parameter_recovery_simulation', 'render_apa_presets_section', 'show_convergence_section'):
        monkeypatch.setattr(app, name, lambda *a, **kw: None)
    monkeypatch.setattr(app, 'qualified_estimation_summary', lambda *a: pd.DataFrame())
    monkeypatch.setattr(app, 'build_likelihood_information_criteria', lambda *a: pd.DataFrame())
    monkeypatch.setattr(app, 'build_output_qualification_table', lambda *a: pd.DataFrame())
    monkeypatch.setattr(app, 'build_weighting_policy_audit', lambda *a: pd.DataFrame())
    frames = []
    monkeypatch.setattr(app.st, 'dataframe', lambda frame, **kw: frames.append(frame.copy()))
    measures = pd.DataFrame({'Facet': ['Person']*3, 'Level': ['P1', 'P2', 'P3'],
        'Estimate': [-2., 0., 2.], 'SE': [1.]*3, 'Infit': [1.]*3, 'Outfit': [1.]*3,
        'SE_Status': ['posterior_eap_sd']*3})
    reliability = app.calc_reliability(measures)
    before = reliability.copy(deep=True)
    result = {'config': {'method': 'MML', 'estimate_population_sd': free_sd}}
    diagnostics = {'measures': measures, 'reliability': reliability}
    expected, note = app._publication_reliability_table(result, diagnostics)
    app._render_report_tables(result, diagnostics)
    rel_frames = [frame for frame in frames if 'Strata' in frame]
    assert len(rel_frames) == 1  # No second summary using a conflicting SE formula.
    pd.testing.assert_frame_equal(rel_frames[0], expected)
    assert 'For EAP person estimates' in note
    assert 'not an observed group count' in note
    assert expected[['RMSE', 'Separation', 'Strata', 'Reliability']].eq('').all().all()
    assert reliability.Strata.iloc[0] == 3.
    assert reliability.Reliability.iloc[0] == .8
    pd.testing.assert_frame_equal(reliability, before)
