"""Raw-score entry shares the qualified engine and never requires an MFRM fit."""
import io

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_g_d_study_boundaries import ratings, result


@pytest.mark.parametrize('two_facets', [False, True])
def test_raw_and_fit_entries_share_math_but_identify_score_basis(two_facets):
    data = ratings(two_facets)
    data['Person'] = data.Person.map({'-1': '001', '1': '002'})
    data['Score'] = data.Score * 0.25 + 0.125
    before = data.copy(deep=True)
    raw = app.compute_generalizability_from_data(
        data, object_facet='Person', random_facets=['Rater', 'Item'] if two_facets else ['Rater'],
        target_population='Students, comparable tasks, trained raters')
    fit = app.compute_generalizability_study(result(data))
    assert raw['available'] and fit['available']
    pd.testing.assert_frame_equal(raw['variance_components'], fit['variance_components'])
    pd.testing.assert_frame_equal(raw['d_study'], fit['d_study'])
    assert raw['observed_coefficients'] == fit['observed_coefficients']
    assert raw['design']['analysis_id'] != fit['design']['analysis_id']
    assert raw['design']['score_basis'] == 'original_observed_scores'
    assert fit['design']['score_basis'] == 'retained_preprocessed_scores'
    assert raw['design']['grand_mean'] == pytest.approx(data.Score.mean())
    assert not raw['bootstrap_ci']['inference_allowed']
    csv = app.g_d_study_export_frame(raw)
    assert csv.loc[csv.Field.eq('target_population'), 'Value'].item().startswith('Students')
    row = csv[csv._section.eq('observed_coefficients')].iloc[0]
    # Orthogonal reference already checked in test_g_d_study_boundaries, scaled by .25².
    assert row.RelativeErrorVariance == pytest.approx((6.375 if two_facets else 2) * .25**2)
    assert row.AbsoluteErrorVariance == pytest.approx((15.625 if two_facets else 5) * .25**2)
    pd.testing.assert_frame_equal(data, before)


def _run(at):
    at.run(timeout=45)
    assert not at.exception, [e.value for e in at.exception]


@pytest.mark.parametrize('lang', ['en', 'ja'])
def test_standalone_main_is_lazy_keeps_language_help_and_invalidates_changed_input(monkeypatch, lang):
    def forbidden(*a, **kw):
        raise AssertionError('Standalone entry must not require MFRM or bootstrap')
    for name in ('mfrm_estimate', 'load_core_namespace', '_bootstrap_g_phi_ci'):
        monkeypatch.setattr(app, name, forbidden)
    calls, exports = [], []
    original_compute = app.compute_generalizability_from_data
    original_download = app.st.download_button
    def compute(*a, **kw):
        calls.append(kw)
        return original_compute(*a, **kw)
    def download(*a, **kw):
        if kw.get('file_name') == 'mfrm_g_d_study.csv':
            exports.append(pd.read_csv(io.BytesIO(kw['data'])))
        return original_download(*a, **kw)
    monkeypatch.setattr(app, 'compute_generalizability_from_data', compute)
    monkeypatch.setattr(app.st, 'download_button', download)
    at = AppTest.from_string('import streamlit_app as app\napp.main()')
    _run(at)
    at.radio(key='lang').set_value(lang); _run(at)
    at.button(key='onboarding_gd').click(); _run(at)
    assert not calls and at.button(key='gd_run').disabled
    at.text_input(key='gd_population').set_value('Students, comparable tasks, trained raters')
    at.checkbox(key='gd_confirm').check(); _run(at)
    at.button(key='gd_run').click(); _run(at)
    assert len(calls) == 1 and len(at.metric) == 2
    analysis_id = at.session_state['_gd_result']['bundle']['design']['analysis_id']
    plan = next(w for w in at.selectbox if w.key.startswith('standalone_gd_plan_'))
    plan.set_value(0); _run(at)
    assert len(calls) == 1
    selected = exports[-1][exports[-1]._section.eq('selected_d_study')].iloc[0]
    assert selected.Rater == selected.Task == 1
    assert selected.RelativeErrorVariance > 0
    assert not any(c in exports[-1] for c in ('G_lower', 'Phi_lower'))
    at.button(key='mfrm_help_open_global').click(); _run(at)
    at.radio(key='lang').set_value('ja' if lang == 'en' else 'en'); _run(at)
    at.button(key='mfrm_help_return_to_source').click(); _run(at)
    assert len(calls) == 1 and len(at.metric) == 2
    assert at.session_state['_gd_result']['bundle']['design']['analysis_id'] == analysis_id
    assert not any('gd_entry.' in c.value for c in at.caption)
    at.text_input(key='gd_population').set_value('Changed population'); _run(at)
    assert not at.metric and not any(w.key == 'standalone_gd_download' for w in at.get('download_button'))
    at.text_input(key='gd_population').set_value('Students, comparable tasks, trained raters'); _run(at)
    assert len(at.metric) == 2 and len(calls) == 1
    at.button(key='gd_back').click(); _run(at)
    assert at.button(key='onboarding_gd')
    at.button(key='onboarding_gd').click(); _run(at)
    assert len(calls) == 1 and len(at.metric) == 2
    assert next(w for w in at.selectbox if w.key.startswith('standalone_gd_plan_')).value == 0
    assert 'facets_mode_output' not in at.session_state


@pytest.mark.parametrize('case,reason', [('missing', 'missing_cells'), ('duplicate', 'duplicate_cells'),
                                         ('nan', 'invalid_scores'), ('weight', 'unsupported_weights')])
def test_raw_input_errors_reject_instead_of_repair(case, reason):
    data = ratings()
    if case == 'missing':
        data = data.iloc[:-1]
    elif case == 'duplicate':
        data = pd.concat([data, data.iloc[[0]]], ignore_index=True)
    elif case == 'nan':
        data.loc[0, 'Score'] = float('nan')
    else:
        data['Weight'] = 2
    bundle = app.compute_generalizability_from_data(data, object_facet='Person', random_facets=['Rater', 'Item'])
    assert not bundle['available'] and bundle['reason_code'] == reason
    assert bundle['design']['score_basis'] == 'original_observed_scores'
    assert bundle['design']['n_input_rows'] == len(data)


def test_return_to_real_mfrm_workspace_preserves_input_result_and_never_refits(monkeypatch):
    from tests.test_analysis_journeys import _prepare, _fit
    from tests.test_help_contextual_app_integration import _block_analysis_calls
    at = _prepare('regular')
    fitted = _fit(at)
    original_text = at.text_area(key='paste_data_text').value
    original_scores = fitted['result']['prep']['data'].copy(deep=True)
    analysis_id = app.build_result_analysis_identity(fitted['result']).analysis_id
    calls = _block_analysis_calls(monkeypatch)
    at.button(key='sidebar_gd').click(); _run(at)
    assert at.sidebar.text_area(key='paste_data_text').value == original_text
    at.text_input(key='gd_population').set_value('Example population')
    at.checkbox(key='gd_confirm').check(); _run(at)
    at.button(key='gd_run').click(); _run(at)
    assert at.session_state['_gd_result']['bundle']['available']
    at.button(key='mfrm_help_open_global').click(); _run(at)
    at.button(key='mfrm_help_return_to_source').click(); _run(at)
    at.button(key='gd_back').click(); _run(at)
    assert at.main.text_area(key='paste_data_text').value == original_text
    assert at.number_input(key='facets_mode_maxit').value == 100
    after = at.session_state['facets_mode_output']['result']
    assert app.build_result_analysis_identity(after).analysis_id == analysis_id
    pd.testing.assert_frame_equal(after['prep']['data'], original_scores)
    assert calls == {"estimate": 0, "refit": 0}


def test_csv_upload_cache_survives_hidden_widget_and_can_be_cleared_and_reloaded(monkeypatch):
    original_upload = app.st.file_uploader
    def upload(*args, **kwargs):
        original_upload(*args, **kwargs)
        if app.st.session_state.pop('_test_upload_pending', False):
            file = io.BytesIO(b'Person,Rater,Score\n001,A,15.25\n001,B,17.25\n002,A,21.25\n002,B,27.25\n')
            file.name = 'scores.csv'
            file.size = len(file.getvalue())
            return file
        return None
    monkeypatch.setattr(app.st, 'file_uploader', upload)
    at = AppTest.from_string('import streamlit_app as app\napp.main()')
    _run(at)
    at.button(key='onboarding_gd').click(); _run(at)
    at.radio(key='gd_source').set_value('upload')
    at.session_state['_test_upload_pending'] = True
    _run(at)
    data = at.session_state['_gd_upload_table']
    assert list(data.Person.unique()) == ['001', '002']
    assert float(data.Score.iloc[0]) == 15.25
    at.text_input(key='gd_population').set_value('Students and trained raters')
    at.checkbox(key='gd_confirm').check(); _run(at)
    at.button(key='gd_run').click(); _run(at)
    assert at.session_state['_gd_result']['bundle']['observed_coefficients']['G'] == pytest.approx(30/32)
    at.button(key='gd_back').click(); _run(at)
    at.button(key='onboarding_gd').click(); _run(at)
    assert at.radio(key='gd_source').value == 'upload' and len(at.metric) == 2
    pd.testing.assert_frame_equal(at.session_state['_gd_upload_table'], data)
    at.button(key='gd_clear_upload').click(); _run(at)
    assert '_gd_upload_table' not in at.session_state and not at.metric
    assert at.get('file_uploader')[0].proto.id.endswith('gd_upload_1')
    at.session_state['_test_upload_pending'] = True
    _run(at)
    pd.testing.assert_frame_equal(at.session_state['_gd_upload_table'], data)


def test_paste_missing_cell_replaces_valid_result_with_explicit_refusal():
    at = AppTest.from_string('import streamlit_app as app\napp.main()')
    _run(at)
    at.button(key='onboarding_gd').click(); _run(at)
    at.radio(key='gd_source').set_value('paste'); _run(at)
    data = ratings(False).to_csv(index=False)
    at.text_area(key='gd_text').set_value(data); _run(at)
    at.text_input(key='gd_population').set_value('Students and trained raters')
    at.checkbox(key='gd_confirm').check(); _run(at)
    at.button(key='gd_run').click(); _run(at)
    assert len(at.metric) == 2
    at.text_area(key='gd_text').set_value(ratings(False).iloc[:-1].to_csv(index=False)); _run(at)
    assert not at.metric
    at.button(key='gd_run').click(); _run(at)
    assert not at.metric
    assert at.session_state['_gd_result']['bundle']['reason_code'] == 'missing_cells'
    assert any('combinations' in item.value for item in at.info)
