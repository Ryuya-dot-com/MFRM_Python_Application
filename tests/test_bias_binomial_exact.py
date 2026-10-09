"""Closed-form probabilities check the deterministic screening counterexample."""
import pytest

from validation.bias_binomial_exact import enumerate_cell


def test_eight_balanced_binary_ratings_have_exact_non_nominal_chisquare_size():
    result=enumerate_cell(8,.5,0.)
    known=result['summaries']['known_baseline'];fitted=result['summaries']['same_data_intercept']
    assert result['probability_sum']==pytest.approx(1.)
    assert known['test_availability']==pytest.approx(254/256)
    assert known['rejection_per_generated_sample']==pytest.approx(16/256)
    assert known['rejection_given_available']==pytest.approx(16/254)
    assert known['coverage_per_generated_sample']==pytest.approx(238/256)
    assert known['coverage_given_interval']==pytest.approx(238/254)
    assert fitted['rejection_per_generated_sample']==0.
    assert fitted['coverage_per_generated_sample'] is None


def test_same_data_intercept_absorbs_nonzero_dgp_effect_and_has_zero_power():
    result=enumerate_cell(16,.3,.8)
    assert result['summaries']['same_data_intercept']['rejection_per_generated_sample']==0.
    assert result['summaries']['known_baseline']['rejection_per_generated_sample']>0.
    for row in result['rows'][1:-1]:
        fit=row['results']['same_data_intercept']
        assert fit['estimate']==pytest.approx(0.,abs=1e-10)
        assert fit['coverage_of_DGP_shift'] is None
