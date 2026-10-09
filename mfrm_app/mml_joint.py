"""App-side transport of the audited finite-GH joint MML kernel.

The numerical kernel and two fixed SLSQP stages come from the implementation
sealed by amendment 7ab8190c672b1ada. Person-level dependence, the exact linear
recession certificate, all Louis cross terms and zero quadrature weights are
preserved. This module imports no validation drivers and grants no inference
qualification. Validation retains its separate sealed implementation.
"""
from fractions import Fraction
from functools import lru_cache
from math import lcm

import numpy as np
from scipy.optimize import linprog, minimize
from scipy.special import roots_hermitenorm


class _NativeRows:
    def __init__(self, data, config):
        if config["model"] not in {"RSM", "PCM"} or config["noncenter_facet"] not in config["facet_names"]:
            raise ValueError("Reference requires RSM/PCM with an explicit noncenter facet")
        if (config.get("population_model", {}).get("enabled") or config.get("dummy_facets")
            or config.get("facet_regularization_enabled") or config.get("facet_regularization", {}).get("enabled")):
            raise ValueError("Reference excludes regression, dummy facets and penalties")
        for facet, spec in config.get("facet_specs", {}).items():
            if (np.isfinite(spec["anchors"]).any() or any(g not in (None, "") for g in spec["groups"])
                or bool(spec["centered"]) != (facet != config["noncenter_facet"])):
                raise ValueError("Reference requires simple sum-to-zero/unconstrained facets without anchors")
        self.person_labels = sorted(data["Person"].astype(str).unique())
        self.person = np.array([self.person_labels.index(str(p)) for p in data["Person"]])
        self.n = len(self.person_labels)
        self.k = np.arange(config["n_cat"], dtype=float)
        raw_score = np.asarray(data["Score"] - config["rating_min"], dtype=float)
        if self.n == 0 or len(self.k) < 2 or not np.isin(raw_score, self.k).all() or ("Weight" in data and not (data["Weight"] == 1).all()):
            raise ValueError("Reference requires nonempty, valid, unit-weight observed rows")
        self.y = raw_score.astype(int)
        counts = {f: len(config["facet_levels"][f]) - int(f != config["noncenter_facet"])
                  for f in config["facet_names"]}
        step_count = 1 if config["model"] == "RSM" else len(config["facet_levels"][config["step_facet"]])
        self.d = sum(counts.values()) + step_count * (len(self.k) - 2)
        self.design = np.zeros((len(data), len(self.k), self.d))
        cursor = 0
        for facet in config["facet_names"]:
            levels = list(config["facet_levels"][facet])
            index = np.array([levels.index(str(v)) for v in data[facet]])
            for j in range(counts[facet]):
                effect = (index == j).astype(float)
                if facet != config["noncenter_facet"]:
                    effect -= (index == len(levels)-1)
                self.design[:, :, cursor+j] = config["facet_signs"][facet] * effect[:, None] * self.k
            cursor += counts[facet]
        step_index = np.zeros(len(data), dtype=int) if config["model"] == "RSM" else np.array([
            list(config["facet_levels"][config["step_facet"]]).index(str(v)) for v in data[config["step_facet"]]])
        for group in range(step_count):
            for j in range(len(self.k)-2):
                # Last step is the negative sum: the highest category has zero cumulative step.
                self.design[step_index == group, :, cursor+j] = -((self.k > j) & (self.k < len(self.k)-1)).astype(float)
            cursor += len(self.k)-2


class GroupedMmlLikelihood:
    def __init__(self, data, config):
        if config['model'] == 'GPCM':
            if config['slope_facet'] != config['step_facet']:
                raise ValueError('Native bounded GPCM requires a common step/slope facet')
            reference = _NativeRows(data, dict(config, model='PCM'))
            linear_d = reference.d
            levels = list(config['facet_levels'][config['slope_facet']])
            index = np.array([levels.index(str(v)) for v in data[config['slope_facet']]])
            basis = np.r_[np.eye(len(levels)-1), -np.ones((1, len(levels)-1))]
            slope_design = basis[index]
        else:
            reference = _NativeRows(data, config)
            linear_d = reference.d
            slope_design = np.empty((len(data), 0))
        features = np.c_[reference.design.reshape(len(data), -1), slope_design]
        _, first, group = np.unique(features, axis=0, return_index=True, return_inverse=True)
        self.design = reference.design[first].copy()
        self.slope_design = slope_design[first].copy()
        self.counts = np.zeros((reference.n, len(first), len(reference.k)))
        np.add.at(self.counts, (reference.person, group, reference.y), 1.)
        # Retain the original row likelihood: add no multinomial coefficients.
        self.cell_totals = self.counts.sum(axis=-1)
        self.k, self.person_labels = reference.k, reference.person_labels
        self.linear_d, self.slope_d = linear_d, self.slope_design.shape[1]
        self.dimension = linear_d+self.slope_d+1
        self.observations = len(data)

    def structural_recession_audit(self):
        """Sufficient certificate of nonexistence, never an existence proof.

        Holding positive slopes and sigma fixed, a direction v with
        (D_iy-D_ik)v >= 0 for every observed row/category makes every observed
        conditional probability nondecreasing, at every theta. One strict
        inequality makes its marginal likelihood strictly increase. Thus no
        finite unpenalized maximum exists. This covers only structural rays;
        absence of a ray says nothing about variance/slope boundaries or other
        paths. It is not the conditional-score-stratum CMLE existence test.

        LP floating-point output is only a candidate. The native design has
        integer entries; rational reconstruction supplies an EXACT integer
        witness, checked against every distinct inequality before certifying.
        """
        present = self.counts.sum(axis=0) > 0
        cell, observed = np.nonzero(present)
        contrast = self.design[cell, observed, None, :]-self.design[cell]
        contrast = np.unique(contrast.reshape(-1, self.linear_d), axis=0)
        result = dict(status='unavailable', no_finite_mle_certified=False,
            finite_mle_existence_established=False, scientific_inference_ready=False,
            scope='unbounded linear structural ray with slopes and population SD held fixed',
            inequalities=len(contrast), linear_dimension=self.linear_d,
            structural_design_rank=int(np.linalg.matrix_rank(contrast)), direction=None)
        if not np.array_equal(contrast, np.rint(contrast)):
            return dict(result, reason='noninteger_design_outside_exact_certificate_scope')
        fit = linprog(-contrast.sum(axis=0), A_ub=-contrast, b_ub=np.zeros(len(contrast)),
            bounds=[(-1., 1.)]*self.linear_d, method='highs',
            options=dict(primal_feasibility_tolerance=1e-9, dual_feasibility_tolerance=1e-9))
        result.update(lp_success=bool(fit.success), lp_status=int(fit.status), lp_message=str(fit.message))
        if not fit.success or fit.x is None or not np.isfinite(fit.x).all():
            return dict(result, reason='linear_program_did_not_return_a_finite_optimum')
        result['lp_candidate'] = fit.x.tolist()
        rational = [Fraction(float(x)).limit_denominator(10**6) for x in fit.x]
        denominator = lcm(*(x.denominator for x in rational))
        direction = [x.numerator*(denominator//x.denominator) for x in rational]
        margins = [sum(int(a)*b for a, b in zip(row, direction, strict=True)) for row in contrast]
        if min(margins) >= 0 and max(margins) > 0:
            result.update(status='no_finite_mle', no_finite_mle_certified=True,
                direction=direction, strict_inequalities=sum(m > 0 for m in margins),
                minimum_exact_margin=min(margins), maximum_exact_margin=max(margins))
        elif min(margins) < 0:
            result.update(reason='candidate_failed_exact_inequality_check')
        else:
            result.update(status='no_structural_ray_detected', reason='not_a_finite_existence_certificate')
        return result

    def evaluate(self, point, nodes, weights, *, information=False, moments=True):
        """Negative log likelihood and analytical joint derivatives at a point.

        Nodes are STANDARD normal coordinates. Their actual theta values move
        as exp(log_sigma)*nodes. Weights, including exactly-zero underflowed
        tail weights, are fixed; a positive normalized mass is required.
        Optimizer calls may omit EAP/SD with moments=False; the full posterior
        is still used for the gradient and Louis information without alteration.
        """
        point, nodes, weights = (np.asarray(x, dtype=float) for x in [point, nodes, weights])
        if point.shape != (self.dimension,) or not np.isfinite(point).all():
            raise ValueError('Invalid joint coordinates')
        if (nodes.ndim != 1 or not len(nodes) or weights.shape != nodes.shape or not np.isfinite(nodes).all()
            or not np.isfinite(weights).all() or np.any(weights < 0) or abs(weights.sum()-1) > 1e-12):
            raise ValueError('Invalid standardized quadrature or prior weights')
        a = np.exp(self.slope_design@point[self.linear_d:-1]); sigma = np.exp(point[-1])
        if not np.isfinite(a).all() or not np.isfinite(sigma) or sigma <= 0 or np.any(a <= 0):
            raise ValueError('Invalid slope or population scale')
        theta = sigma*nodes
        logits = a[None, :, None]*(self.design@point[:self.linear_d]+theta[:, None, None]*self.k)
        if not np.isfinite(logits).all():
            raise FloatingPointError('Nonfinite category logits')
        shifted = logits-logits.max(axis=-1, keepdims=True)
        # With finite logits, exp(shifted) is in [0,1] and each cell has a
        # unit term. Reuse its positive normalizer for log-p and p.
        exponential = np.exp(shifted)
        normalizer = exponential.sum(axis=-1, keepdims=True)
        log_probability = shifted-np.log(normalizer)
        probability = exponential/normalizer
        counts_flat = self.counts.reshape(len(self.counts), -1)
        log_likelihood = counts_flat @ log_probability.reshape(len(nodes), -1).T
        log_weight = np.full_like(weights, -np.inf)
        np.log(weights, out=log_weight, where=weights > 0)
        log_joint = log_likelihood+log_weight
        maximum = log_joint.max(axis=1, keepdims=True)
        relative_mass = np.exp(log_joint-maximum)
        mass = relative_mass.sum(axis=1, keepdims=True)
        marginal = (maximum+np.log(mass)).ravel()
        posterior = relative_mass/mass

        # Exchange finite sums for the gradient: posterior-weighted category
        # counts suffice. Person/node score arrays are needed only for Louis'
        # information below, not for every optimizer value/gradient call.
        weighted_count = (posterior.T @ counts_flat).reshape(len(nodes), *self.counts.shape[1:])
        total = weighted_count.sum(axis=-1)
        residual = total[:, :, None]*probability-weighted_count
        gradient = np.empty(self.dimension)
        gradient[:self.linear_d] = (residual.sum(axis=0)*a[:, None]).ravel() @ self.design.reshape(-1, self.linear_d)
        gradient[self.linear_d:-1] = (residual*logits).sum(axis=(0, 2)) @ self.slope_design
        gradient[-1] = np.sum((residual @ self.k)*theta[:, None]*a[None, :])
        result = dict(nll=float(-marginal.sum()), gradient=gradient,
            posterior=posterior, cells=len(self.design), observations=self.observations)
        if moments:
            mean_theta = posterior@theta
            variance = np.sum(posterior*(theta[None, :]-mean_theta[:, None])**2, axis=1)
            result.update(eap=mean_theta, sd=np.sqrt(variance))
        if information:
            derivative = np.zeros((*logits.shape, self.dimension))
            derivative[..., :self.linear_d] = a[None, :, None, None]*self.design
            derivative[..., self.linear_d:-1] = logits[..., None]*self.slope_design[None, :, None, :]
            derivative[..., -1] = theta[:, None, None]*a[None, :, None]*self.k
            mean_derivative = np.einsum('qik,qikd->qid', probability, derivative, optimize=True)
            score = np.einsum('ni,qid->nqd', self.cell_totals, mean_derivative, optimize=True)
            score -= np.einsum('nik,qikd->nqd', self.counts, derivative, optimize=True)
            mean_score = np.einsum('nq,nqd->nd', posterior, score, optimize=True)
            centered = derivative-mean_derivative[:, :, None, :]
            hessian = np.einsum('qi,qik,qikd,qike->de', total, probability, centered, centered, optimize=True)
            if self.slope_d:
                cross = np.einsum('qik,qikd,ie->de', residual, derivative[..., :self.linear_d], self.slope_design, optimize=True)
                hessian[:self.linear_d, self.linear_d:-1] += cross
                hessian[self.linear_d:-1, :self.linear_d] += cross.T
                hessian[self.linear_d:-1, self.linear_d:-1] += np.einsum(
                    'qik,qik,id,ie->de', residual, logits, self.slope_design, self.slope_design, optimize=True)
                sigma_cross = np.einsum('qik,qik,id->d', residual, derivative[..., -1], self.slope_design, optimize=True)
                hessian[self.linear_d:-1, -1] += sigma_cross
                hessian[-1, self.linear_d:-1] += sigma_cross
            hessian[-1, -1] += np.sum(residual*derivative[..., -1])
            centered_score = score-mean_score[:, None, :]
            hessian -= np.einsum('nq,nqd,nqe->de', posterior, centered_score, centered_score, optimize=True)
            result['information'] = hessian
        return result


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
    constraints = problem.joint_constraints

    def value_gradient(point):
        value = grouped.evaluate(point, nodes, weights, moments=False)
        return value['nll'], value['gradient']

    for options in OPTIONS:
        initial = start.copy()
        fit = minimize(value_gradient, initial, jac=True, method='SLSQP', bounds=bounds,
            constraints=constraints, options=options)
        result['runs'].append(dict(initial=initial, point=fit.x, nll=float(fit.fun), gradient=fit.jac,
            success=bool(fit.success), status=int(fit.status), message=str(fit.message), iterations=int(fit.nit), evaluations=int(fit.nfev)))
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
