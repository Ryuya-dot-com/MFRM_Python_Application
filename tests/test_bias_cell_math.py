"""Closed-form scalar likelihood controls, failure denominators, and G3 holds."""
import numpy as np
import pandas as pd
import pytest
from scipy.special import expit
from scipy.stats import chi2

import streamlit_app as app


def binary_cell(scores, slope=1., weight=1., eta=0., bound=10.):
    y = np.asarray(scores, dtype=float)
    a = np.full(y.size, slope)
    w = np.full(y.size, weight)
    def probability(b):
        p = expit(a*(eta+b))
        return np.column_stack([1-p, p])
    def nll(b):
        z = a*(eta+b)
        return float(np.sum(w*(np.logaddexp(0., z)-y*z)))
    return app._fit_conditional_bias_cell(nll, probability, y, a, w, bound), nll


@pytest.mark.parametrize('slope', [.5, 1., 2.])
def test_binary_bias_agrees_with_closed_form_log_odds_and_information(slope):
    y = [0]*3+[1]*7
    fit, nll = binary_cell(y, slope=slope, eta=.2)
    expected = np.log(7/3)/slope-.2
    assert fit['status'] == 'interior'
    assert fit['estimate'] == pytest.approx(expected, abs=1e-11)
    assert fit['se'] == pytest.approx(1/np.sqrt(10*slope**2*.7*.3), rel=1e-12)
    assert abs(fit['score']) < 1e-10
    lr = 2*(nll(0)-nll(expected))
    assert fit['lr'] == pytest.approx(lr, abs=1e-12)
    assert fit['p'] == pytest.approx(chi2.sf(lr, 1), abs=1e-12)
    assert fit['ci']['status'] == 'ok'
    for end in ['lower','upper']:
        assert 2*(nll(fit['ci'][end])-nll(expected)) == pytest.approx(chi2.ppf(.95, 1), abs=1e-7)


@pytest.mark.parametrize('eta',[0.,1e-8,-1e-8,1e-6])
def test_near_null_lr_retains_small_signal_without_subtracting_nearly_equal_nlls(eta):
    fit,_=binary_cell([0]*5+[1]*5,eta=eta)
    # Exact LR = 2*n*log(cosh(eta/2)); this form avoids rounding cosh to 1.
    expected=20*np.log1p(2*np.sinh(eta/4)**2)
    assert fit['status']=='interior'
    assert fit['lr']==pytest.approx(expected,rel=2e-7,abs=5e-24)
    assert fit['p']==pytest.approx(chi2.sf(expected,1),abs=2e-15)


def test_conditional_probability_and_nll_disagreement_is_not_given_a_p_value():
    y=np.array([0.,1.]*5);a=np.ones(10);w=np.ones(10)
    def probability(b):
        p=expit(b+.1);return np.column_stack([np.full(10,1-p),np.full(10,p)])
    # This likelihood has the right maximizer but the wrong curvature/scale.
    def wrong_nll(b):return 2*float(np.sum(np.logaddexp(0.,b+.1)-y*(b+.1)))
    fit=app._fit_conditional_bias_cell(wrong_nll,probability,y,a,w,10.)
    assert fit['status']=='likelihood_failed' and np.isnan(fit['p'])


@pytest.mark.parametrize('score', [0, 1])
def test_extreme_binary_cells_have_no_finite_mle(score):
    fit, _ = binary_cell([score]*10)
    assert fit['status'] == 'no_finite_mle'
    assert fit['candidate'] == (-10 if score == 0 else 10)
    for key in ['estimate', 'se', 'p', 'lr']:
        assert np.isnan(fit[key])
    assert fit['probabilities'] is None


def test_interior_solution_outside_search_range_is_not_reported_as_mle():
    fit, _ = binary_cell([0]*3+[1]*7, bound=.1)
    assert fit['status'] == 'search_range_limited'
    assert fit['candidate'] == .1
    assert np.isnan(fit['estimate']) and np.isnan(fit['p'])


def test_arbitrary_weights_do_not_acquire_a_chi_square_reference():
    fit, _ = binary_cell([0]*3+[1]*7, weight=.5)
    assert fit['status'] == 'interior'
    assert fit['estimate'] == pytest.approx(np.log(7/3))
    assert fit['se'] == pytest.approx(1/np.sqrt(5*.7*.3))
    assert fit['reference_status'] == 'weighted_reference_unqualified'
    assert np.isnan(fit['p'])
    assert fit['ci']['status'] == 'not available'


@pytest.mark.parametrize('weight', [-1., float('nan'), 0.])
def test_invalid_or_zero_support_is_not_silently_dropped(weight):
    fit, _ = binary_cell([0, 1], weight=weight)
    assert fit['status'] == 'invalid_input'
    assert np.isnan(fit['estimate'])


def test_scalar_exception_is_retained_instead_of_using_zero_shift(monkeypatch):
    def broken(*args, **kwargs):
        raise RuntimeError('retained failure control')
    monkeypatch.setattr(app, 'root_scalar', broken)
    fit, _ = binary_cell([0, 1])
    assert fit['status'] == 'calculation_failed'
    assert 'retained failure control' in fit['error']
    assert fit['probabilities'] is None
    assert np.isnan(fit['p'])


def policy_table():
    return pd.DataFrame({
        'FacetA': ['Rater']*3+['Person'], 'FacetB': ['Task']*4,
        'FacetA_Level': ['R1', 'R2', 'R3', 'P1'], 'FacetB_Level': ['T1']*4,
        'Bias Size': [.6, .2, np.nan, .1], 'S.E.': [.2]*4,
        'ObsN': [10]*4, 'PositiveWeightObsN': [10]*4,
        'Prob.': [1.e-15]*4,  # Legacy t values must never determine the screen.
        'LR Prob.': [.02, .03, np.nan, .04],
        'Primary P Method': ['unit_weight_conditional_lr_v1']*4,
        'Reference Status': ['unit_weight_conditional_lr_candidate']*4,
        'Cell Status': ['interior', 'interior', 'solver_failed', 'interior'],
    })


def test_primary_lr_holm_uses_each_family_and_retains_failed_tests():
    table = app.build_dff_bias_screening_table(policy_table()).set_index('FacetA_Level')
    assert table.loc['R1', 'p'] == .02
    assert table.loc['R1', 'p_holm'] == pytest.approx(.06)
    assert table.loc['R2', 'p_holm'] == pytest.approx(.06)
    # BH .045 cannot create a statistical flag; Holm is the sole candidate rule.
    assert table.loc['R2', 'p_bh'] == pytest.approx(.045)
    assert table.loc['R2', 'EvidenceLevel'] == 'No flag'
    assert table.loc['R1', 'EvidenceLevel'] == 'Practical review'
    assert table.loc['P1', 'p_holm'] == .04
    assert table.loc['R3', 'EvidenceLevel'] == 'Inference unavailable'
    assert table.loc['R1', 'FamilySize'] == 3
    assert table.loc['P1', 'FamilySize'] == 1
    assert not table.PublicConclusionAllowed.any()


def test_legacy_and_all_missing_cells_never_become_no_flag_or_ready():
    table = policy_table().drop(columns=['Cell Status', 'Reference Status', 'Primary P Method'])
    dff = app.build_dff_bias_screening_table(table)
    assert dff.EvidenceLevel.eq('Inference unavailable').all()
    assert dff.p.isna().all()
    audit = app.build_bias_inference_audit(table)
    assert audit.Status.ne('Ready').all()
    assert audit.UnavailableCells.iloc[0] == 4
    assert not audit.PublicConclusionAllowed.any()
    publication = app._bias_table_for_publication({'cells': {'table': table}}, 'GPCM')
    assert not {'Prob.', 'LR Prob.', 'Profile CI Lower'} & set(publication.columns)
    assert publication['Cell Status'].eq('legacy_unverified').all()
    assert not publication.PublicConclusionAllowed.any()


@pytest.mark.parametrize('method,expected', [('holm',[.04,np.nan,.09,.09]), ('bh',[.04,np.nan,4*.04/3,4*.04/3])])
def test_adjustment_counts_missing_family_members(method, expected):
    np.testing.assert_allclose(app.adjust_p_values([.01,np.nan,.03,.04], method), expected, equal_nan=True)


@pytest.mark.parametrize('values,method', [([-.01], 'holm'), ([np.inf], 'bh'), ([[.1]], 'holm'), ([.1], 'typo')])
def test_invalid_multiplicity_input_is_explicit(values, method):
    with pytest.raises(ValueError):
        app.adjust_p_values(values, method)


@pytest.mark.parametrize('model', ['RSM','PCM','GPCM'])
def test_fitted_cells_share_scalar_reference_policy_and_retained_family(model):
    result = app.mfrm_estimate(app._make_self_test_gradient_data(), 'Person', ['Rater','Task'],
                               'Score', model=model, method='JMLE', maxit=2)
    diag = app.mfrm_diagnostics(result, compute_pca=False, interaction_pairs=[])
    out = app.estimate_bias_interaction(result, diag, 'Rater', 'Task', omit_extreme=False)
    table = out['table']
    interior = table['Cell Status'].eq('interior')
    assert interior.any()
    assert table.loc[interior, 'LR Prob.'].notna().all()
    np.testing.assert_allclose(table['Prob.'], table['LR Prob.'], equal_nan=True)
    assert table['t'].isna().all() and table['d.f.'].isna().all()
    assert table['Conditional Information'].loc[interior].gt(0).all()
    assert table.PublicConclusionAllowed.eq(False).all()
    assert out['chi_sq'].FixedProb.isna().all()
    # Explicitly exclude an extreme level but keep all originally observed cells.
    diag['obs'].loc[diag['obs'].Rater == 'R1', 'Observed'] = result['prep']['rating_min']
    excluded = app.estimate_bias_interaction(result, diag, 'Rater','Task', omit_extreme=True)['table']
    assert len(excluded) == len(table)
    assert excluded['Cell Status'].eq('excluded_extreme_level').any()
    assert app.build_dff_bias_screening_table(excluded).FamilySize.eq(len(table)).all()


def test_family_labels_cannot_merge_distinct_facet_pairs():
    table = policy_table().iloc[[0,3]].copy()
    table['FacetA'] = ['a', 'a x b']
    table['FacetB'] = ['b x c', 'c']
    dff = app.build_dff_bias_screening_table(table)
    assert dff.FamilySize.eq(1).all()
    assert dff.MultiplicityFamilyID.nunique() == 2


def test_duplicate_cells_do_not_silently_inflate_the_family():
    table = policy_table()
    with pytest.raises(ValueError, match='Duplicate'):
        app.build_dff_bias_screening_table(pd.concat([table,table.iloc[[0]]]))


@pytest.mark.parametrize('metric',['Bias Size (logits)','Standardized Bias'])
def test_heatmap_shows_descriptive_values_without_legacy_significance_markers(monkeypatch, metric):
    table=pd.DataFrame({'FacetA_Level':['R1','R1','R2'], 'FacetB_Level':['T1','T2','T1'],
                        'Bias Size':[.2,.3,np.nan], 'Standardized Bias':[1.,1.5,np.nan],
                        't':[100.,100.,100.], 'Prob.':[0.,0.,0.], 'Observd Count':[10.,20.,30.]})
    plots=[]
    monkeypatch.setattr(app.st,'radio',lambda *a,**k:metric)
    monkeypatch.setattr(app.st,'plotly_chart',lambda fig,**kwargs:plots.append(fig))
    app._draw_bias_heatmap(table,'Rater','Task')
    assert len(plots)==1
    heat=plots[0].data[0]
    assert heat.text is None
    assert all('Prob.' not in str(value) and 't: 100' not in str(value) for value in heat.customdata)
    np.testing.assert_allclose(heat.z, [[.2,.3] if metric=='Bias Size (logits)' else [1.,1.5]])


def test_calibration_on_same_cell_can_absorb_the_entire_bias_score():
    # Binary intercept calibrated on this very cell: eta_hat=logit(S/n).
    # Its conditional bias MLE and LR statistic are identically zero for
    # every non-extreme S, so a chi-square reference is not established by
    # positive scalar curvature. This is an estimand/calibration control.
    n = 10
    for successes in range(1,n):
        eta = np.log(successes/(n-successes))
        fit, _ = binary_cell([0]*(n-successes)+[1]*successes, eta=eta)
        assert fit['status'] == 'interior'
        assert fit['information'] > 0
        assert abs(fit['estimate']) < 1e-10
        assert abs(fit['lr']) < 1e-12
