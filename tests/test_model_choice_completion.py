"""Unfinished fits cannot win an information-criterion comparison."""
from copy import deepcopy
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
from scipy.stats import chi2
from streamlit.testing.v1 import AppTest
import streamlit_app as app
from tests.test_information_criteria_units import legacy_result


def fits():
    out = {}
    for model, k, ll in [('RSM', 5, -123.), ('PCM', 7, -90.), ('GPCM', 9, -120.)]:
        result = legacy_result(k=k, ll=ll)
        result['config']['model'] = model
        result['config']['facet_names'] = ['Criterion']
        result['prep']['data']['Criterion'] = ['C1', 'C2']
        if model != 'RSM':
            result['config']['step_facet'] = 'Criterion'
        if model == 'GPCM':
            result['config']['slope_facet'] = 'Criterion'
        result['summary']['Model'] = model
        result['opt'] = SimpleNamespace(success=True, message='completed')
        out[model] = result
    return out


def test_completed_fits_have_correct_differences_weights_and_held_lr(monkeypatch):
    results = fits()
    monkeypatch.setattr(app, '_refit_candidate_model', lambda result, model, **kwargs: results[model])
    bundle = app.compute_model_choice_comparison(results['RSM'])
    table = bundle['comparison'].set_index('Model')
    assert table.ComparisonIncluded.all()
    for criterion in ['AIC', 'BIC']:
        expected_delta = table[criterion] - table[criterion].min()
        np.testing.assert_allclose(table['Delta' + criterion], expected_delta)
        relative = np.exp(-expected_delta / 2.)
        np.testing.assert_allclose(table[criterion + 'Weight'], relative / relative.sum())
        np.testing.assert_allclose(table[criterion + 'EvidenceRatio'], np.exp(expected_delta / 2.))
    lr = bundle['lr_tests']
    assert not lr.empty
    for _, row in lr.iterrows():
        expected_lr = 2. * (table.loc[row.Alternative, 'LogLik'] - table.loc[row.Null, 'LogLik'])
        assert row.ChiSq == pytest.approx(expected_lr)
        if np.isfinite(row.p):
            assert row.p == pytest.approx(chi2.sf(expected_lr, row.df))
    assert lr.QualificationStatus.eq('WITHHELD').all()
    assert not lr.PublicConclusionAllowed.any()
    assert set(lr.loc[lr.Alternative.eq('GPCM'), 'ReasonCode']) == {'model.model_002.gpcm_lrt_uncalibrated'}
    assert bundle['recommendation']['tier'] == 'withheld'


@pytest.mark.parametrize('termination', [False, None, 'true'])
def test_unfinished_lowest_ic_keeps_values_but_cannot_enter_comparison(monkeypatch, termination):
    results = fits()
    results['PCM']['opt'] = SimpleNamespace(success=termination, message='completion unconfirmed')
    results['PCM']['summary']['Converged'] = True  # An old label must not override the optimizer.
    original = deepcopy(results['PCM']['summary'])
    monkeypatch.setattr(app, '_refit_candidate_model', lambda result, model, **kwargs: results[model])
    bundle = app.compute_model_choice_comparison(results['RSM'])
    table = bundle['comparison'].set_index('Model')
    row = table.loc['PCM']
    assert row.FitStatus == 'refit_termination_unconfirmed'
    assert row.AIC == 194. and np.isnan(row.DeltaAIC) and np.isnan(row.AICWeight)
    assert not row.ComparisonIncluded and not row.AlgorithmTerminated
    assert table.loc['RSM', 'DeltaAIC'] == 0. and table.loc['GPCM', 'DeltaAIC'] == 2.
    assert table.AICWeight.sum() == pytest.approx(1.)
    assert bundle['technical_recommendation']['model'] == 'RSM'
    assert 'PCM' not in bundle['lr_tests'][['Null', 'Alternative']].to_numpy()
    assert bundle['recommendation']['tier'] == 'withheld'
    pd.testing.assert_frame_equal(results['PCM']['summary'], original)


def test_a_single_completed_fit_has_no_difference_weight_or_ranking(monkeypatch):
    results = fits()
    results['RSM']['opt'].success = False
    results['PCM']['summary']['LogLik'] = np.nan
    monkeypatch.setattr(app, '_refit_candidate_model', lambda result, model, **kwargs: results[model])
    bundle = app.compute_model_choice_comparison(results['RSM'])
    table = bundle['comparison'].set_index('Model')
    assert table.loc['RSM', 'FitStatus'] == 'current_termination_unconfirmed'
    assert table.loc['PCM', 'FitStatus'] == 'refit_nonfinite'
    assert table[['DeltaAIC', 'DeltaBIC', 'AICWeight', 'BICWeight']].isna().all().all()
    assert bundle['technical_recommendation']['model'] == ''
    assert bundle['lr_tests'].empty


def _view(result, bundle, lang):
    import streamlit as st
    import streamlit_app as app
    from unittest.mock import patch
    st.session_state['lang'] = lang
    if 'seeded' not in st.session_state:
        key = app._model_comparison_cache_key(result, 'Criterion', 'Criterion')
        st.session_state[f'model_choice_comparison::{key}'] = bundle
        st.session_state['seeded'] = True
    with patch.object(app, 'mfrm_estimate', side_effect=AssertionError('No automatic refit')):
        app.render_model_choice_guidance(result)


@pytest.mark.parametrize('lang', ['ja', 'en'])
def test_comparison_ui_has_readable_status_and_hides_internal_columns(monkeypatch, lang):
    results = fits(); results['PCM']['opt'].success = False
    monkeypatch.setattr(app, '_refit_candidate_model', lambda result, model, **kwargs: results[model])
    bundle = app.compute_model_choice_comparison(results['RSM'])
    at = AppTest.from_function(_view, args=(results['RSM'], bundle, lang)).run()
    assert not at.exception
    table = at.table[0].value
    assert not {'FitStatus', 'AICcStatus', 'AICWeight', 'ComparisonIncluded', 'AlgorithmTerminated'}.intersection(table.columns)
    assert ('終了未確認' if lang == 'ja' else 'completion unconfirmed') in table.to_string()
    assert not at.success


def test_old_comparison_is_not_reused_or_automatically_refitted():
    at = AppTest.from_function(_view, args=(fits()['RSM'], {'available': True}, 'en')).run()
    assert not at.exception and not at.table
    assert len(at.button) == 1
