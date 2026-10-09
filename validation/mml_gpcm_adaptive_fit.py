"""Development bounded GPCM calibration with independent numerical gates.

Every planned start uses the same adaptive rule, including its polish.
Continuous-normal integration reviews the final points. No global optimum,
unrestricted GPCM, sampling performance, or public inference is certified.
"""
import time

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import minimize, nnls

from validation.mml_gpcm_adaptive import AdaptiveGpcmLikelihood, normal_log_rule
from validation.mml_gpcm_integral import GpcmIntegral
from validation.mml_gpcm_joint_check import joint_constraints
from validation.mml_grouped_fit import OPTIONS


LIMITS = dict(score=1e-4, correction=1e-5, feasibility=1e-10, interior=1e-6, complementarity=1e-8,
    symmetry=1e-6, condition=1e10, information=1e-3, relative_information=.005,
    information_step=1e-3, nll_per_person=1e-6, gradient=1e-4, moments=1e-4,
    reference_mass=1e-7, reference_tail_mass=1e-10,
    reference_information=1e-6, reference_gradient=1e-7, reference_nll=1e-7,
    reference_moments=1e-7, multistart_point=1e-4, multistart_nll=1e-6)
REFERENCE_SETTINGS = (dict(bound=16., rel_tol=1e-11, abs_tol=1e-13),
                      dict(bound=20., rel_tol=1e-12, abs_tol=1e-14))
DIFFERENCE_STEPS = (1e-4, 5e-5)


def boundary_diagnostics(point, gradient, bounds, constraints):
    """KKT residual for ALL linear inequalities, including the last slope.

    Margins have the form A*x-b >= 0, so stationarity is g-A.T*mu=0,
    mu>=0. Boundary stationarity does not authorize an interior Wald SE.
    """
    point, gradient = np.asarray(point, float), np.asarray(gradient, float)
    if point.ndim != 1 or gradient.shape != point.shape or not np.isfinite([point, gradient]).all():
        raise ValueError('Invalid point or gradient for boundary review')
    if len(bounds) != len(point):
        raise ValueError('Bounds must match the joint coordinates')
    rows, margins, labels = [], [], []
    def add(row, value, lower, upper, label):
        if lower is not None and np.isfinite(lower):
            rows.append(row); margins.append(value-lower); labels.append(label+':lower')
        if upper is not None and np.isfinite(upper):
            rows.append(-row); margins.append(upper-value); labels.append(label+':upper')
    for j, (lower, upper) in enumerate(bounds):
        add(np.eye(len(point))[j], point[j], lower, upper, 'coordinate_'+str(j))
    for i, constraint in enumerate(constraints):
        a = np.atleast_2d(constraint.A)
        for j, (row, value, lower, upper) in enumerate(zip(a, a@point,
                np.broadcast_to(constraint.lb, len(a)), np.broadcast_to(constraint.ub, len(a)), strict=True)):
            add(row, value, lower, upper, f'linear_{i}_{j}')
    margins = np.asarray(margins)
    feasible = bool(np.all(margins >= -LIMITS['feasibility']))
    active = margins <= LIMITS['interior']
    multipliers = np.zeros(int(active.sum()))
    residual = None
    if feasible:
        if active.any():
            matrix = np.asarray(rows)[active].T
            multipliers, _ = nnls(matrix, gradient)
            residual = gradient-matrix@multipliers
        else:
            residual = gradient
    score = None if residual is None else float(np.max(abs(residual)))
    complementarity = (float(np.max(abs(multipliers*margins[active]))) if active.any() else 0.) if feasible else None
    return dict(feasible=feasible, interior=feasible and not bool(active.any()),
        minimum_margin=float(margins.min()) if len(margins) else None,
        margins=margins, active_labels=np.asarray(labels)[active].tolist(), multipliers=multipliers,
        raw_score=float(np.max(abs(gradient))), stationarity_residual=score,
        complementarity_residual=complementarity,
        kkt_residual=max(score, complementarity) if feasible else None,
        constrained_stationary=feasible and score < LIMITS['score'] and complementarity < LIMITS['complementarity'])


def information_comparison(reference, approximate):
    """Relative perturbation in every contrast direction, without flooring.

    rho = ||L^-1 (H_approx-H_ref) L^-T||_2, H_ref=L L.T.
    For rho<1, SE_approx/SE_ref lies between (1+rho)^(-1/2) and
    (1-rho)^(-1/2) for every fixed linear contrast of the joint coordinates.
    """
    matrices = [np.asarray(x, float) for x in [reference, approximate]]
    if (matrices[0].ndim != 2 or matrices[0].shape[0] != matrices[0].shape[1]
        or matrices[1].shape != matrices[0].shape or not matrices[0].size
        or not all(np.isfinite(x).all() for x in matrices)):
        raise ValueError('Invalid joint information matrices')
    symmetry = [float(np.max(abs(x-x.T))) for x in matrices]
    symmetric = [(x+x.T)/2 for x in matrices]
    eigenvalues = [np.linalg.eigvalsh(x) for x in symmetric]
    positive = all(v[0] > 0 for v in eigenvalues)
    conditions = [float(v[-1]/v[0]) if v[0] > 0 else None for v in eigenvalues]
    absolute = float(np.max(abs(matrices[1]-matrices[0])))
    relative, se_range = None, None
    if positive and max(symmetry) < LIMITS['symmetry'] and max(conditions) < LIMITS['condition']:
        factor = np.linalg.cholesky(symmetric[0])
        left = solve_triangular(factor, symmetric[1]-symmetric[0], lower=True)
        whitened = solve_triangular(factor, left.T, lower=True).T
        relative = float(np.max(abs(np.linalg.eigvalsh((whitened+whitened.T)/2))))
        if relative < 1:
            se_range = [float(1/np.sqrt(1+relative)), float(1/np.sqrt(1-relative))]
    checks = dict(positive_definite=positive, symmetry=max(symmetry) < LIMITS['symmetry'],
        condition=all(c is not None and c < LIMITS['condition'] for c in conditions),
        absolute=absolute < LIMITS['information'],
        relative=relative is not None and relative < LIMITS['relative_information'])
    return dict(eigenvalues=eigenvalues, symmetry=symmetry, conditions=conditions,
        absolute_error=absolute, relative_error=relative, se_ratio_bounds=se_range,
        checks=checks, passed=all(checks.values()))


def planned_starts(evaluator, problem):
    """Three deterministic independent starts; no result-selected retries."""
    base = np.r_[problem.structural_start, np.log(problem.sigma_start)]
    offset = np.zeros_like(base)
    offset[:evaluator.linear_d] = .25*np.sin(np.arange(evaluator.linear_d)+1)
    expanded = np.linspace(-.4, .4, evaluator.slope_d+1)
    offset[evaluator.linear_d:-1] = expanded[:-1]
    offset[-1] = np.log(1.5)
    return [base, base+offset, base-offset]


def reference_diagnostics(point, continuous, bounds, constraints):
    boundary = boundary_diagnostics(point, continuous['gradient'], bounds, constraints)
    h = np.asarray(continuous['information'])
    matrix = information_comparison(h, h)
    correction = relative_allowance = correction_allowance = None
    if matrix['passed']:
        inverse = np.linalg.inv((h+h.T)/2)
        correction = float(np.max(abs(inverse@continuous['gradient'])))
        # A uniform entrywise error e has spectral norm at most dimension*e.
        # quad_vec supplies an estimate, not a rigorous interval enclosure.
        combined_information_error = (continuous['numeric_information_error_allowance']+
            continuous['tail_information_abs_bound']+2*continuous['numeric_gradient_error_allowance']*
            continuous['tail_gradient_abs_bound'])
        relative_allowance = len(point)*combined_information_error/matrix['eigenvalues'][0][0]
        correction_allowance = float(np.max(np.sum(abs(inverse), axis=1))*(
            continuous['numeric_gradient_error_allowance']+continuous['tail_gradient_abs_bound']))
    checks = dict(interior=boundary['interior'], stationary=boundary['constrained_stationary'],
        matrix=matrix['passed'], newton_correction=correction is not None and correction < LIMITS['correction'],
        mass=continuous['numeric_relative_mass_error_sum'] < LIMITS['reference_mass'],
        tail=continuous['tail_relative_mass_bound_sum'] < LIMITS['reference_tail_mass'],
        relative_information_allowance=relative_allowance is not None and relative_allowance < LIMITS['relative_information'],
        correction_allowance=correction_allowance is not None and correction_allowance < LIMITS['correction'])
    return dict(boundary=boundary, matrix=matrix, newton_correction=correction,
        relative_information_allowance=relative_allowance, correction_allowance=correction_allowance,
        checks=checks, passed=all(checks.values()))


def review_candidate(evaluator, reference, point, rule, bounds, constraints):
    actual = evaluator.evaluate_log_weights(point, *rule)
    finite_boundary = boundary_diagnostics(point, actual['gradient'], bounds, constraints)
    if not finite_boundary['feasible']:
        return dict(passed=False, status='infeasible', finite_boundary=finite_boundary)
    columns = []
    for h in DIFFERENCE_STEPS:
        columns.append(np.asarray([(evaluator.evaluate_log_weights(point+d, *rule, moments=False)['gradient']-
                                    evaluator.evaluate_log_weights(point-d, *rule, moments=False)['gradient'])/(2*h)
                                   for d in np.eye(len(point))*h]).T)
    information = (4*columns[1]-columns[0])/3
    continuous = reference.continuous(point, information=True, **REFERENCE_SETTINGS[0])
    reference_review = reference_diagnostics(point, continuous, bounds, constraints)
    comparison = information_comparison(continuous['information'], information)
    errors = {key:float(np.max(abs(actual[key]-continuous[key]))) for key in ['gradient','eap','sd']}
    errors['nll_per_person'] = abs(actual['nll']-continuous['nll'])/reference.n
    step_change = float(np.max(abs(columns[0]-columns[1])))
    checks = dict(interior=finite_boundary['interior'],
        finite_stationary=finite_boundary['constrained_stationary'],
        continuous_review=reference_review['passed'],
        information=comparison['passed'], information_step=step_change < LIMITS['information_step'],
        nll=errors['nll_per_person'] < LIMITS['nll_per_person'], gradient=errors['gradient'] < LIMITS['gradient'],
        moments=max(errors['eap'],errors['sd']) < LIMITS['moments'])
    return dict(passed=all(checks.values()), status='reviewed', finite=actual, continuous=continuous,
        finite_information=information, information_step_change=step_change, information_comparison=comparison,
        finite_boundary=finite_boundary, reference_review=reference_review, errors=errors, checks=checks)


def assess_candidates(candidates, expected):
    """One good candidate never erases a failed or missing planned start."""
    usable = [i for i,c in enumerate(candidates) if 'continuous' in c.get('review', {})]
    selected = min(usable, key=lambda i:candidates[i]['review']['continuous']['nll']) if usable else None
    point_spread = nll_spread = None
    complete = len(candidates) == expected and len(usable) == expected
    if complete:
        point_spread = float(np.max(np.ptp([c['point'] for c in candidates], axis=0)))
        nll_spread = float(np.ptp([c['review']['continuous']['nll'] for c in candidates]))
    checks = dict(all_starts_present=complete,
        all_candidates_pass=complete and all(c['review']['passed'] and
            len(c['runs']) == len(OPTIONS) and all(r['success'] for r in c['runs']) for c in candidates),
        point_agreement=point_spread is not None and point_spread < LIMITS['multistart_point'],
        nll_agreement=nll_spread is not None and nll_spread < LIMITS['multistart_nll'])
    return dict(selected=selected, point_spread=point_spread, nll_spread=nll_spread,
                checks=checks, passed=all(checks.values()), global_optimum_established=False)


def fit_adaptive_gpcm(data, problem, *, q=121):
    """Same-rule cold starts and mandatory polishes, then independent review.

    The native finite slope/SD bounds define the target; no clipping, flooring,
    missing-run deletion, or automatic increase in quadrature order is used.
    """
    evaluator = AdaptiveGpcmLikelihood(data, problem.config)
    reference = GpcmIntegral(data, problem.config)
    rule = normal_log_rule(q)
    bounds = list(problem.structural_bounds)+[tuple(np.log(problem.sigma_bounds))]
    constraints = joint_constraints(problem)
    starts = planned_starts(evaluator, problem)
    result = dict(classification='DEVELOPMENT_ONLY', scientific_inference_ready=False,
        numerical_pass=False, global_optimum_established=False, finite_mle_existence_established=False,
        q=q, planned_starts=starts, candidates=[], point=None, covariance=None,
        target='Native bounded GPCM: all expanded log-slopes [-3,3], population SD [.05,10]',
        optimizer_objective='Mean NLL for initial search; summed NLL for mandatory polish, same quadrature',
        optimizer_divisors=[evaluator.observations, 1.],
        options=OPTIONS, limits=LIMITS, reference_settings=REFERENCE_SETTINGS)
    # Do not silently relabel a different target as the native bounded study.
    if tuple(problem.config.get('gpcm_log_slope_bounds', ())) != (-3.,3.) or tuple(problem.sigma_bounds) != (.05,10.):
        raise ValueError('This development fitter retains the native slope and SD bounds')
    for point in starts:
        if not boundary_diagnostics(point, np.zeros(len(point)), bounds, constraints)['interior']:
            raise ValueError('A planned cold start is not feasible and interior')
    recession = evaluator.structural_recession_audit()
    result['recession'] = recession
    if recession['status'] == 'unavailable' or recession['no_finite_mle_certified']:
        return dict(result, status='prefit_'+recession['status'])
    if recession['structural_design_rank'] != recession['linear_dimension']:
        return dict(result, status='prefit_structural_rank_deficient')
    for index, initial in enumerate(starts):
        candidate = dict(index=index, initial=initial, runs=[], point=None)
        def value_gradient(point):
            try:
                value = evaluator.evaluate_log_weights(point, *rule, moments=False)
            except (ValueError, FloatingPointError, np.linalg.LinAlgError, RuntimeError):
                candidate['failed_evaluation_point'] = np.asarray(point).copy()
                raise
            # Positive constant scaling preserves the target. It prevents an
            # identity-Hessian first step from growing with response count.
            return value['nll']/divisor, value['gradient']/divisor
        try:
            point = initial.copy()
            for options, divisor in zip(OPTIONS, result['optimizer_divisors'], strict=True):
                started = time.perf_counter()
                fit = minimize(value_gradient, point, jac=True, method='SLSQP', bounds=bounds,
                               constraints=constraints, options=options)
                candidate['runs'].append(dict(initial=point, point=fit.x,
                    nll=float(fit.fun)*divisor, objective_value=float(fit.fun), objective_divisor=divisor,
                    success=bool(fit.success), status=int(fit.status), message=str(fit.message),
                    iterations=int(fit.nit), evaluations=int(fit.nfev), seconds=time.perf_counter()-started))
                point = np.asarray(fit.x, float)
                candidate['point'] = point
                if not boundary_diagnostics(point, np.zeros(len(point)), bounds, constraints)['feasible']:
                    raise ValueError('Optimizer returned an infeasible point')
            candidate['review'] = review_candidate(evaluator, reference, point, rule, bounds, constraints)
        except (ValueError, FloatingPointError, np.linalg.LinAlgError, RuntimeError) as exc:
            candidate['failure'] = dict(type=type(exc).__name__, message=str(exc))
        result['candidates'].append(candidate)
    assessment = assess_candidates(result['candidates'], len(starts))
    result['multistart'] = assessment
    refinement_pass = False
    if assessment['selected'] is not None:
        chosen = result['candidates'][assessment['selected']]
        result['point'] = chosen['point']
        try:
            tight = reference.continuous(chosen['point'], information=True, **REFERENCE_SETTINGS[1])
            previous = chosen['review']['continuous']
            changes = {k:float(np.max(abs(tight[k]-previous[k]))) for k in ['nll','gradient','information','eap','sd']}
            relative = information_comparison(tight['information'], previous['information'])
            tight_review = reference_diagnostics(chosen['point'], tight, bounds, constraints)
            checks = dict(nll=changes['nll'] < LIMITS['reference_nll'],
                gradient=changes['gradient'] < LIMITS['reference_gradient'],
                information=changes['information'] < LIMITS['reference_information'],
                moments=max(changes['eap'],changes['sd']) < LIMITS['reference_moments'],
                relative_information=relative['passed'], tighter_reference=tight_review['passed'])
            result['reference_refinement'] = dict(continuous=tight, changes=changes,
                reference_review=tight_review, checks=checks)
            refinement_pass = all(checks.values())
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            result['reference_refinement_failure'] = dict(type=type(exc).__name__, message=str(exc))
    result['numerical_pass'] = assessment['passed'] and refinement_pass
    result['status'] = 'numerically_reviewed_interior_candidate' if result['numerical_pass'] else 'numerical_review_failed'
    if result['numerical_pass']:
        # Use the reviewed matrix; no diagonal-only uncertainty or eigenvalue repair.
        h = chosen['review']['continuous']['information']
        result['covariance'] = np.linalg.inv((h+h.T)/2)
    return result
