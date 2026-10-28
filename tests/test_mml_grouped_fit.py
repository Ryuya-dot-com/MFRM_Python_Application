"""No-fit failures, dependent slope bounds, and false solver-success controls."""
from types import SimpleNamespace

import numpy as np
import pytest

from tests.test_mml_grouped_likelihood import problem as fixture
from validation.mml_gpcm_joint_check import prepare_problem
from validation.mml_grouped_fit import fit_grouped_joint, normal_rule


def native_problem(missing=None):
    data, *_ = fixture('GPCM', 'regular')
    if missing is not None:
        data.loc[data.Score == missing, 'Score'] = (missing+1) % 4
    _, problem = prepare_problem(data, 31)
    return data, problem


def test_no_finite_mle_certificate_prevents_the_optimizer_from_running(monkeypatch):
    import validation.mml_grouped_fit as module
    data, problem = native_problem(missing=1)
    def forbidden(*args, **kwargs): raise AssertionError('Optimizer must not run on a certified ray')
    monkeypatch.setattr(module, 'minimize', forbidden)
    result = fit_grouped_joint(data, problem)
    assert result['status'] == 'prefit_no_finite_mle'
    assert result['point'] is result['covariance'] is None
    assert not result['runs']


def test_dependent_last_slope_is_constrained_and_solver_success_does_not_override_failures(monkeypatch):
    import validation.mml_grouped_fit as module
    data, problem = native_problem(); calls = []
    def false_success(fun, start, **kwargs):
        point = np.array(start); point[-3:-1] = 2.  # Last expanded log-slope = -4, outside [-3,3].
        constraint = kwargs['constraints'][0]
        assert np.any(constraint.A@point < constraint.lb) or np.any(constraint.A@point > constraint.ub)
        value, gradient = fun(point)
        calls.append(kwargs['options'])
        return SimpleNamespace(x=point, fun=value, jac=gradient, success=True, status=0, message='false success', nit=1)
    monkeypatch.setattr(module, 'minimize', false_success)
    result = fit_grouped_joint(data, problem, q=31)
    assert len(calls) == len(result['runs']) == 2
    assert not result['numerical_pass'] and result['covariance'] is None
    assert not result['checks']['interior'] and not result['checks']['constraints']
    assert all(run['success'] for run in result['runs'])


def test_large_order_rule_has_finite_standard_normal_moments_and_zero_tail_mass():
    z, w = normal_rule(641)
    assert np.isfinite(z).all() and np.isfinite(w).all() and np.any(w == 0)
    assert w.sum() == pytest.approx(1., abs=1e-15)
    assert w@z**2 == pytest.approx(1., abs=1e-13)
    assert w@z**4 == pytest.approx(3., abs=1e-12)
    for invalid in [True, 31.5, 0]:
        with pytest.raises(ValueError): normal_rule(invalid)
