"""Audit decisions and displayed source cells follow row positions, not labels."""
import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


def audit(data, weighted=True):
    return app.build_response_data_audit(data, person_col='Person', facet_cols=['Rater','Task'],
                                        score_col='Score', weight_col='Weight' if weighted else None)


@pytest.mark.parametrize('index', [[8,3,3,1,7,5,2,0], list('aacdbefa'), list(reversed(range(8)))])
def test_mixed_rows_preserve_source_values_reasons_and_decisions(index):
    data = pd.DataFrame(dict(Person=['p0','', 'p2','p3','p4','p5','p6','p7'],
        Rater=['r0','r1','r2','r0','r1','r2','r0',''], Task=['t0']*8,
        Score=[0,1,None,'bad',1.5,2,3,''], Weight=[1,1,1,1,1,0,'bad',None]), index=index)
    original = data.copy(deep=True)
    result = audit(data)
    rows = result['rows']
    assert rows.OriginalIndex.tolist() == list(map(str,index))
    assert rows.Person.tolist() == data.Person.tolist()
    assert rows.Person.dtype == object
    assert rows.Score.tolist() == [0,1,'','bad',1.5,2,3,'']
    assert rows.IncludedInLikelihood.tolist() == [True]+[False]*7
    assert rows.LikelihoodAction.tolist() == ['included','excluded','excluded','excluded',
                                              'error_before_fit','excluded','excluded','excluded']
    assert rows.ExclusionReason.tolist() == ['', 'missing person','missing score','non-numeric score',
        'fractional score','non-positive weight','missing or non-numeric weight',
        'missing score; missing Rater; missing or non-numeric weight']
    assert rows.MissingColumns.tolist() == ['', 'Person','Score','Score','','','Weight','Score, Rater, Weight']
    assert result['excluded_rows'].InputRow.tolist() == list(range(2,9))
    expected = audit(data.reset_index(drop=True))
    pd.testing.assert_frame_equal(result['summary'], expected['summary'], check_exact=True)
    pd.testing.assert_frame_equal(data, original, check_exact=True)


def test_empty_missing_column_and_unweighted_audits_use_positions():
    data = pd.DataFrame(dict(Person=['p1','p0'], Rater=['r1','r0'], Task=['t0','t1'], Score=[1,0]), index=[12,12])
    result = audit(data, weighted=False)
    assert result['rows'].IncludedInLikelihood.all()
    assert result['rows'].Person.tolist() == ['p1','p0']
    assert result['excluded_rows'].empty
    missing = audit(data.drop(columns='Task'), weighted=False)['rows']
    assert missing.Person.tolist() == ['p1','p0']
    assert missing.LikelihoodAction.tolist() == ['error_before_fit']*2
    assert audit(data.iloc[:0], weighted=False)['rows'].empty


def test_duplicate_labels_do_not_change_prepared_likelihood_rows():
    data = pd.DataFrame(dict(Person=['p0','p0','p1','p1','p2'], Rater=['r0','r1','r0','r1','r0'],
                             Score=[0,1,2,1,np.nan]), index=['a','a','b','a','c'])
    kwargs = dict(person_col='Person',facet_cols=['Rater'],score_col='Score',rating_min=0,rating_max=2)
    repeated = app.prepare_mfrm_data(data, **kwargs)
    positional = app.prepare_mfrm_data(data.reset_index(drop=True), **kwargs)
    pd.testing.assert_frame_equal(repeated['data'].reset_index(drop=True),
                                  positional['data'].reset_index(drop=True), check_exact=True)
    assert repeated['data'].index.tolist() == ['a','a','b','a']
    assert repeated['levels'] == positional['levels']
    assert len(repeated['data']) == 4
