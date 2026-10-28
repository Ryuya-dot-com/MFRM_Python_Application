"""Development joint MML fits under the existing native parameterization.

Finite-GH stationarity, not continuous-integral or sampling qualification.
The fixed two-stage optimizer is retained in full; no result-selected retry.
"""
from functools import lru_cache

import numpy as np
from scipy.optimize import minimize
from scipy.special import roots_hermitenorm

from validation.mml_gpcm_joint_check import joint_constraints
from validation.mml_grouped_likelihood import GroupedMmlLikelihood


OPTIONS = (dict(maxiter=500, ftol=1e-12), dict(maxiter=500, ftol=1e-14))
LIMITS = dict(score=1e-4, newton_correction=1e-5, constraints=1e-10,
    interior=1e-6, symmetry=1e-10, condition=1e10, restart=1e-5, objective_worsening=1e-8)


def normal_rule(q):
    if isinstance(q, (bool, np.bool_)) or not isinstance(q, (int, np.integer)) or q < 3:
        raise ValueError('Quadrature order must be an integer of at least three')
    return _normal_rule_cached(int(q))


@lru_cache(maxsize=8)
def _normal_rule_cached(q):
    # SciPy remains finite at high orders where numpy.hermgauss can overflow.
    nodes, weights = roots_hermitenorm(q)
    weights /= np.sqrt(2*np.pi)
    nodes.setflags(write=False); weights.setflags(write=False)
    return nodes, weights


def fit_grouped_joint(data, problem, *, q=321):
    """Cold start and one mandatory polish, including all dependent constraints.

    RSM/PCM structural coordinates are unbounded; GPCM additionally bounds
    every expanded log-slope, including the dependent last coordinate. Sigma
    is bounded by the native problem. Only interior PD points return a joint
    covariance, still marked development-only. A certified structural ray is
    a no-fit failure, not an extreme finite estimate.
    """
    grouped = GroupedMmlLikelihood(data, problem.config)
    nodes, weights = normal_rule(q)
    audit = grouped.structural_recession_audit()
    result = dict(classification='DEVELOPMENT_ONLY', scientific_inference_ready=False,
        numerical_pass=False, q=int(q), recession=audit, runs=[], point=None, covariance=None)
    if audit['status'] == 'unavailable' or audit['no_finite_mle_certified']:
        return dict(result, status='prefit_'+audit['status'])
    if audit['structural_design_rank'] != audit['linear_dimension']:
        return dict(result, status='prefit_structural_rank_deficient')
    start = np.r_[problem.structural_start, np.log(problem.sigma_start)]
    bounds = list(problem.structural_bounds)+[tuple(np.log(problem.sigma_bounds))]
    constraints = joint_constraints(problem)

    def value_gradient(point):
        value = grouped.evaluate(point, nodes, weights, moments=False)
        return value['nll'], value['gradient']

    for options in OPTIONS:
        initial = start.copy()
        fit = minimize(value_gradient, initial, jac=True, method='SLSQP', bounds=bounds,
            constraints=constraints, options=options)
        result['runs'].append(dict(initial=initial, point=fit.x, nll=float(fit.fun), gradient=fit.jac,
            success=bool(fit.success), status=int(fit.status), message=str(fit.message), iterations=int(fit.nit)))
        start = fit.x
        if not np.isfinite(start).all():
            return dict(result, status='nonfinite_optimizer_point')
    evaluated = grouped.evaluate(start, nodes, weights, information=True)
    information = evaluated['information']
    symmetry = float(np.max(abs(information-information.T)))
    symmetric = (information+information.T)/2
    eigen = np.linalg.eigvalsh(symmetric)
    covariance = np.linalg.inv(symmetric) if eigen.min() > 0 else None
    distances = [float(v-lower) for v, (lower, _) in zip(start, bounds) if lower is not None]
    distances += [float(upper-v) for v, (_, upper) in zip(start, bounds) if upper is not None]
    for constraint in constraints:
        value = constraint.A@start
        distances.extend((value-constraint.lb).ravel()); distances.extend((constraint.ub-value).ravel())
    boundary_distance = float(min(distances))
    constraint_residual = float(problem.constraint_residual(start[:-1]))
    correction = float(np.max(abs(covariance@evaluated['gradient']))) if covariance is not None else None
    condition = float(eigen.max()/eigen.min()) if covariance is not None else None
    restart = float(np.max(abs(result['runs'][1]['point']-result['runs'][0]['point'])))
    checks = dict(optimizer_terminated=all(run['success'] for run in result['runs']),
        score=float(np.max(abs(evaluated['gradient']))) < LIMITS['score'],
        newton_correction=correction is not None and correction < LIMITS['newton_correction'],
        constraints=constraint_residual <= LIMITS['constraints'], interior=boundary_distance > LIMITS['interior'],
        positive_definite=covariance is not None, symmetric=symmetry < LIMITS['symmetry'],
        condition=condition is not None and condition < LIMITS['condition'], restart=restart < LIMITS['restart'],
        objective_nonworsening=result['runs'][1]['nll']-result['runs'][0]['nll'] <= LIMITS['objective_worsening'])
    numerical_pass = all(checks.values())
    result.update(status='finite_rule_stationary' if numerical_pass else 'numerical_checks_failed',
        point=start, evaluated=evaluated, information_eigenvalues=eigen, symmetry_residual=symmetry,
        condition_number=condition, newton_correction=correction, constraint_residual=constraint_residual,
        boundary_distance=boundary_distance, restart_difference=restart, checks=checks, numerical_pass=numerical_pass,
        covariance=covariance if numerical_pass else None)
    return result
