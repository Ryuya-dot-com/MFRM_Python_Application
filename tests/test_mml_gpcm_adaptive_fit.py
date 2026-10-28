"""Acceptance controls that must survive apparent solver success."""
import copy

import numpy as np
import pytest
from scipy.optimize import LinearConstraint

from tests.test_mml_grouped_fit import native_problem
from validation.mml_gpcm_adaptive import AdaptiveGpcmLikelihood
from validation.mml_gpcm_adaptive_fit import (assess_candidates, boundary_diagnostics,
    fit_adaptive_gpcm, information_comparison, planned_starts, reference_diagnostics)
from validation.mml_gpcm_joint_check import joint_constraints


def test_dependent_slope_boundary_has_constrained_but_not_interior_stationarity():
    bounds=[(-3.,3.),(-3.,3.),(np.log(.05),np.log(10.))]
    constraints=[LinearConstraint([[-1.,-1.,0.]],-3.,3.)]
    at_boundary=boundary_diagnostics([1.5,1.5,0.],[-1.,-1.,0.],bounds,constraints)
    assert at_boundary['feasible'] and at_boundary['constrained_stationary']
    assert not at_boundary['interior'] and at_boundary['raw_score']==1.
    assert at_boundary['active_labels']==['linear_0_0:lower']
    opposite=boundary_diagnostics([1.5,1.5,0.],[1.,1.,0.],bounds,constraints)
    assert not opposite['constrained_stationary']
    outside=boundary_diagnostics([2.,2.,0.],[0.,0.,0.],bounds,constraints)
    assert not outside['feasible'] and outside['kkt_residual'] is None


def test_small_absolute_information_error_can_destroy_a_weak_contrast():
    result=information_comparison(np.diag([1.,1e-6]),np.diag([1.,2e-6]))
    assert result['checks']['absolute']
    assert result['relative_error']==pytest.approx(1.)
    assert not result['passed']


def test_near_boundary_stationarity_also_requires_complementarity():
    result=boundary_diagnostics([3.-5e-7],[-1e6],[(-3.,3.)],[])
    assert result['feasible'] and not result['interior']
    assert result['stationarity_residual']==0.
    assert result['complementarity_residual']==pytest.approx(.5)
    assert not result['constrained_stationary']


def test_relative_information_review_is_invariant_to_invertible_linear_coordinates():
    ref=np.array([[3.,.4],[.4,1.]])
    approximate=ref+np.array([[.0001,.00003],[.00003,-.0002]])
    change=np.array([[2.,.3],[0.,.5]])
    a=information_comparison(ref,approximate)
    b=information_comparison(change.T@ref@change,change.T@approximate@change)
    assert a['relative_error']==pytest.approx(b['relative_error'],abs=1e-14)
    for vector in [np.array([1.,0.]),np.array([0.,1.]),np.array([1.,-3.])]:
        ratio=np.sqrt((vector@np.linalg.solve(approximate,vector))/(vector@np.linalg.solve(ref,vector)))
        assert a['se_ratio_bounds'][0]<=ratio<=a['se_ratio_bounds'][1]


@pytest.mark.parametrize('matrix',[np.diag([1.,0.]),np.diag([1.,-1.]),np.diag([1.,1e-12]),
                                  np.array([[1.,.01],[0.,1.]])])
def test_nonpositive_ill_conditioned_or_asymmetric_information_is_not_repaired(matrix):
    result=information_comparison(np.eye(2),matrix)
    assert not result['passed'] and result['se_ratio_bounds'] is None


def test_reference_error_allowance_also_respects_the_weak_information_direction():
    continuous=dict(gradient=np.zeros(2),information=np.diag([1.,1e-5]),
        numeric_information_error_allowance=1e-6,tail_information_abs_bound=0.,
        numeric_gradient_error_allowance=1e-10,tail_gradient_abs_bound=0.,
        numeric_relative_mass_error_sum=1e-12,tail_relative_mass_bound_sum=0.)
    actual=reference_diagnostics(np.zeros(2),continuous,[(None,None)]*2,[])
    assert actual['boundary']['interior'] and actual['matrix']['passed']
    assert not actual['checks']['relative_information_allowance'] and not actual['passed']


def test_failed_missing_or_disagreeing_start_is_not_hidden_by_the_best_candidate():
    good=dict(point=np.zeros(2),runs=[dict(success=True),dict(success=True)],
              review=dict(passed=True,continuous=dict(nll=1.)))
    records=[copy.deepcopy(good) for _ in range(3)]
    assert assess_candidates(records,3)['passed']
    assert not assess_candidates(records,3)['global_optimum_established']
    records[-1]['runs'][0]['success']=False
    assert not assess_candidates(records,3)['passed']
    assert not assess_candidates(records[:2],3)['passed']
    records[-1]=copy.deepcopy(good);records[-1]['point'][0]=.01
    assert not assess_candidates(records,3)['checks']['point_agreement']
    records[-1]=dict(point=None,runs=[],failure=dict(message='numerical exception'))
    assert not assess_candidates(records,3)['passed']


def test_all_cold_starts_are_distinct_interior_and_keep_dependent_constraints():
    data,problem=native_problem(); evaluator=AdaptiveGpcmLikelihood(data,problem.config)
    starts=planned_starts(evaluator,problem)
    assert len({tuple(p) for p in starts})==3
    for point in starts:
        result=boundary_diagnostics(point,np.zeros_like(point),
            list(problem.structural_bounds)+[tuple(np.log(problem.sigma_bounds))],joint_constraints(problem))
        assert result['interior']


def test_structural_recession_prevents_fit_and_covariance(monkeypatch):
    import validation.mml_gpcm_adaptive_fit as module
    data,problem=native_problem(missing=1)
    def forbidden(*a,**kw):raise AssertionError('A certified ray must prevent estimation')
    monkeypatch.setattr(module,'minimize',forbidden)
    result=fit_adaptive_gpcm(data,problem,q=31)
    assert result['status']=='prefit_no_finite_mle'
    assert result['point'] is result['covariance'] is None
    assert not result['candidates'] and not result['scientific_inference_ready']


def test_failed_objective_point_and_all_planned_starts_are_retained(monkeypatch):
    data,problem=native_problem()
    def failed_evaluation(*args,**kwargs):
        raise FloatingPointError('Unrepresentable trial-point calculation')
    monkeypatch.setattr(AdaptiveGpcmLikelihood,'evaluate_log_weights',failed_evaluation)
    result=fit_adaptive_gpcm(data,problem,q=31)
    assert len(result['candidates'])==len(result['planned_starts'])==3
    for candidate,initial in zip(result['candidates'],result['planned_starts'],strict=True):
        np.testing.assert_array_equal(candidate['failed_evaluation_point'],initial)
        assert candidate['failure']['type']=='FloatingPointError'
    assert result['point'] is result['covariance'] is None
    assert not result['numerical_pass'] and not result['scientific_inference_ready']
