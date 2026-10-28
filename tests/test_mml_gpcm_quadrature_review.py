import copy
import pytest

from validation.mml_gpcm_quadrature_review import review


def evidence():
    rows = [dict(id=f'GPCM_{design}_Q{q}_{start}', all_checks_pass=True, expanded_se_relative_difference=1e-5)
            for design in ['regular', 'sparse', 'unbalanced']
            for q, start in [(121, 'from61'), (241, 'from61'), (321, 'from61'), (321, 'from31')]]
    fit = dict(comparisons=copy.deepcopy(rows), completed=12, planned=12, failures=[], R_returncode=0,
               cross_start_coordinate_difference=dict(regular=1e-8, sparse=1e-8, unbalanced=1e-8))
    moments = dict(comparisons=copy.deepcopy(rows), complete=True, closed_form=dict(passed=True))
    info = dict(comparisons=copy.deepcopy(rows), completed=12, planned=12, failures=[])
    return fit, moments, info


def test_coarse_failure_is_retained_while_highest_order_scope_stays_technical():
    fit, moments, info = evidence()
    info['comparisons'][0]['all_checks_pass'] = False
    result = review(fit, moments, info)
    assert result['complete'] and result['highest_order_all_development_checks_pass']
    assert not result['all_orders_all_checks_pass']
    assert not result['comparisons'][0]['all_checks_pass']
    assert not result['scientific_inference_ready']


@pytest.mark.parametrize('fault', ['missing_coarse', 'missing_highest', 'duplicate', 'fit_failure', 'failed_control', 'start_disagreement'])
def test_missing_or_failed_evidence_cannot_be_hidden_by_highest_order_summary(fault):
    fit, moments, info = evidence()
    if fault == 'missing_coarse': info['comparisons'].pop(0)
    elif fault == 'missing_highest': info['comparisons'].pop(2)
    elif fault == 'duplicate': info['comparisons'].append(copy.deepcopy(info['comparisons'][0]))
    elif fault == 'fit_failure': fit['failures'].append(dict(id='GPCM_regular_Q121_from61', error='failed'))
    elif fault == 'failed_control': moments['closed_form']['passed'] = False
    elif fault == 'start_disagreement': fit['cross_start_coordinate_difference']['regular'] = 1e-6
    result = review(fit, moments, info)
    assert not result['highest_order_all_development_checks_pass']
    assert not result['scientific_inference_ready']
