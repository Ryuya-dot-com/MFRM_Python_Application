"""An optimizer's successful stop cannot rescue a confounded IC parameter count."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_model_choice_exports import view
from tests.test_model_choice_completion import fits


def fixture(confounded):
    result = fits()['RSM']
    pairs = [('A','A'), ('B','B')] if confounded else [('A','A'),('A','B'),('B','A'),('B','B')]
    result['prep']['data'] = pd.DataFrame([
        dict(Person=f'P{i}',Rater=r,Task=t,Score=i%3) for i,(r,t) in enumerate(pairs)])
    result['config'].update(facet_names=['Rater','Task'], method='MML',
        facet_specs={f:app.build_facet_constraint(['A','B'],centered=True) for f in ['Rater','Task']})
    return result


def test_mml_confounded_locations_stop_before_refitting(monkeypatch):
    result = fixture(True); original = deepcopy(result)
    monkeypatch.setattr(app,'_refit_candidate_model',lambda *a,**k:pytest.fail('No refits for an unidentified design'))
    bundle = app.compute_model_choice_comparison(result,step_facet='Task')
    assert not bundle['available'] and 'confounded' in bundle['reason']
    assert bundle['comparison'].empty and bundle['lr_tests'].empty
    pd.testing.assert_frame_equal(result['summary'],original['summary'])


def test_observed_crossing_restores_necessary_rank_condition():
    # Matching labels alone are harmless: it is the observed combinations that matter.
    assert app._model_choice_can_compare(fixture(False))[0]
    result = fixture(True)
    result['config']['facet_names'] = ['Rater']
    assert app._model_choice_can_compare(result)[0]


@pytest.mark.parametrize('identified',[False,np.bool_(False)])
def test_existing_jmle_rank_failure_is_respected(identified):
    result = fixture(False); result['config']['method']='JMLE'
    result['identifiability']={'eta_structurally_identified':identified}
    assert not app._model_choice_can_compare(result)[0]


@pytest.mark.parametrize('lang',['ja','en'])
def test_rank_failure_blocks_old_cached_comparison_in_both_languages(lang):
    at=AppTest.from_function(view)
    result=fixture(True)
    at.session_state['result']=result;at.session_state['lang']=lang
    cache=app._model_comparison_cache_key(result,'Task','Task')
    at.session_state[f'model_choice_comparison::{cache}']={'available':True,'old_marker':True}
    at.run()
    assert not at.exception and len(at.info)==1
    assert not at.table and not at.button and not at.selectbox
    assert 'report_tables.' not in at.info[0].value


def test_unknown_levels_and_excessive_rank_check_are_withheld():
    result=fixture(False)
    result['prep']['data'].loc[0,'Rater']='unknown'
    assert not app._model_choice_can_compare(result)[0]
    result=fixture(False)
    result['config']['facet_specs']['Rater']['n_params']=2_000_000
    assert not app._model_choice_can_compare(result)[0]


def test_missing_rows_or_facet_columns_do_not_crash_the_rank_guard():
    result=fixture(False); result['prep']['data']=result['prep']['data'].iloc[:0]
    assert not app._model_choice_can_compare(result)[0]
    result=fixture(False); result['prep']['data']=result['prep']['data'].drop(columns='Task')
    assert not app._model_choice_can_compare(result)[0]
