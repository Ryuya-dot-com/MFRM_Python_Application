"""Exact finite-GH development likelihood using Person/cell/category counts.

Compressing repeated covariate patterns changes neither the likelihood nor its
derivatives. Complete-data scores remain grouped by Person before applying
Louis' missing-information identity. Supports simple identified, unit-weight
RSM/PCM and bounded GPCM, with normal mean zero and free log(population SD).
This reference does not change the app or assert continuous-integral accuracy.
"""
from fractions import Fraction
from math import lcm

import numpy as np
from scipy.optimize import linprog

from validation.mml_gpcm_integral import GpcmIntegral
from validation.mml_rasch_integral import RaschIntegral


class GroupedMmlLikelihood:
    def __init__(self, data, config):
        if config['model'] == 'GPCM':
            if config['slope_facet'] != config['step_facet']:
                raise ValueError('Native bounded GPCM requires a common step/slope facet')
            reference = GpcmIntegral(data, config)
            linear_d = reference.linear_d
            slope_design = reference.slope_design
        else:
            reference = RaschIntegral(data, config)
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
