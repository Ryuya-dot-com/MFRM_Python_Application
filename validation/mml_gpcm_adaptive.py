"""Development GPCM adaptive likelihood with the exact finite-sum gradient.

The original normal prior is retained. Posterior mode, local scale and
Jacobian depend on all free coordinates, including log(population SD).
Uses the existing simple, unit-weight, common step/slope-owner design.
No optimizer, covariance, continuous-integral accuracy or sampling claim.
"""

from functools import lru_cache

import numpy as np
from scipy.special import logsumexp, roots_hermitenorm

from validation.mml_grouped_likelihood import GroupedMmlLikelihood


def normal_log_rule(order):
    """Standard-normal GH nodes and log weights, including underflowed tails.

    For orthonormal probabilists' Hermite polynomials p_j, at a root x of
    He_n the normalized weight is 1/(n*p_(n-1)(x)**2). Scaled recurrence
    retains log weights that cannot first be represented as doubles.
    """
    if isinstance(order, (bool, np.bool_)) or not isinstance(order, (int, np.integer)) or order < 3:
        raise ValueError('Quadrature order must be an integer of at least three')
    return _normal_log_rule_cached(int(order))


@lru_cache(maxsize=8)
def _normal_log_rule_cached(order):
    nodes, _ = roots_hermitenorm(order)
    previous, current = np.ones(order), nodes.copy()
    log_scale = np.zeros(order)
    for degree in range(2, order):
        following = (nodes*current-np.sqrt(degree-1)*previous)/np.sqrt(degree)
        large = np.maximum(abs(following), abs(current)) > 1e100
        if large.any():
            following[large] *= 1e-100
            current[large] *= 1e-100
            log_scale[large] += 100*np.log(10.)
        previous, current = current, following
    log_weights = -np.log(order)-2*(np.log(abs(current))+log_scale)
    if not np.isfinite(log_weights).all() or abs(logsumexp(log_weights)) > 1e-12:
        raise FloatingPointError('Invalid logarithmic Hermite weights')
    nodes.setflags(write=False); log_weights.setflags(write=False)
    return nodes, log_weights


class AdaptiveGpcmLikelihood(GroupedMmlLikelihood):
    def __init__(self, data, config):
        if config['model'] != 'GPCM':
            raise ValueError('Adaptive development evaluator requires GPCM')
        super().__init__(data, config)
        self.observed_score = self.counts @ self.k

    def evaluate(self, point, nodes, weights, *, moments=True):
        """Differentiate the supplied floating-point rule, without clipping.

        Zero input weights contribute exactly zero. For high-order adaptive
        GH use evaluate_log_weights with normal_log_rule: a tiny base weight
        need not remain negligible after the adaptive density correction.
        """
        nodes, weights = (np.asarray(x, dtype=float) for x in [nodes, weights])
        if (weights.shape != nodes.shape or not np.isfinite(weights).all()
            or np.any(weights < 0) or abs(weights.sum()-1) > 1e-12):
            raise ValueError('Invalid standardized quadrature or prior weights')
        log_weights = np.full_like(weights, -np.inf)
        np.log(weights, out=log_weights, where=weights > 0)
        return self.evaluate_log_weights(point, nodes, log_weights, moments=moments)

    def evaluate_log_weights(self, point, nodes, log_weights, *, moments=True):
        """Exact gradient of the finite adaptive sum in standard-normal z.

        z = mode + scale*x. Count compression preserves individual Person
        likelihoods without multinomial factors. Log weights are normalized
        for the original rule; relocated weights are never renormalized.
        """
        point, nodes, log_weights = (np.asarray(x, dtype=float) for x in [point, nodes, log_weights])
        if point.shape != (self.dimension,) or not np.isfinite(point).all():
            raise ValueError('Invalid joint coordinates')
        if (nodes.ndim != 1 or not len(nodes) or log_weights.shape != nodes.shape
            or not np.isfinite(nodes).all() or np.isnan(log_weights).any()
            or np.isposinf(log_weights).any() or abs(logsumexp(log_weights)) > 1e-12):
            raise ValueError('Invalid standardized quadrature or log prior weights')
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            a = np.exp(self.slope_design @ point[self.linear_d:-1])
            sigma = np.exp(point[-1])
            if np.any(a <= 0) or sigma <= 0:
                raise FloatingPointError('Nonpositive representable slope or population scale')
            return self._evaluate(point, nodes, log_weights, a, sigma, moments)

    def _evaluate(self, point, nodes, log_weight, a, sigma, moments):
        base = self.design @ point[:self.linear_d]
        totals = self.cell_totals

        def kernel(z):
            logits = a[None, :, None]*(base + sigma*z[:, None, None]*self.k)
            probability = np.exp(logits-logsumexp(logits, axis=-1, keepdims=True))
            centered = self.k - (probability @ self.k)[..., None]
            variance = np.sum(probability*centered**2, axis=-1)
            residual = self.observed_score-totals*(probability @ self.k)
            score = sigma*(residual @ a)-z
            curvature = 1+sigma**2*((totals*variance) @ (a*a))
            if not np.isfinite(score).all() or not np.isfinite(curvature).all():
                raise FloatingPointError('Nonfinite posterior score or curvature')
            return score, curvature, probability, centered, variance, logits, residual

        # Since 0 <= E(K) <= Kmax, these bounds bracket every posterior mode.
        # h''(z) <= -1: the score residual also bounds the root error in z.
        lower = sigma*((self.observed_score-totals*self.k[-1]) @ a)
        upper = sigma*(self.observed_score @ a)
        mode = np.zeros(len(totals))
        for iteration in range(100):
            score, curvature, probability, centered, variance, logits, residual = kernel(mode)
            active = abs(score) > 1e-10
            if not active.any():
                break
            lower = np.where(active & (score > 0), mode, lower)
            upper = np.where(active & (score < 0), mode, upper)
            proposed = mode+score/curvature
            inside = (proposed > lower) & (proposed < upper)
            mode = np.where(active, np.where(inside, proposed, (lower+upper)/2), mode)
        else:
            raise FloatingPointError('Posterior mode did not meet the score tolerance')
        scale = 1/np.sqrt(curvature)

        def project(coefficient):
            # Logit derivatives at FIXED standardized z, not fixed theta.
            out = np.empty((len(totals), self.dimension))
            out[:, :self.linear_d] = np.einsum('nck,c,ckd->nd', coefficient, a,
                                              self.design, optimize=True)
            out[:, self.linear_d:-1] = np.sum(coefficient*logits, axis=-1) @ self.slope_design
            out[:, -1] = sigma*mode*((coefficient @ self.k) @ a)
            return out

        # Implicit differentiation of h_z(mode, psi)=0, H=-h_zz>0.
        mixed = -sigma*project(totals[..., None]*a[None, :, None]*probability*centered)
        mixed[:, self.linear_d:-1] += sigma*(residual*a) @ self.slope_design
        mixed[:, -1] += sigma*(residual @ a)
        mode_derivative = mixed/curvature[:, None]
        curvature_derivative = sigma**2*project(
            totals[..., None]*a[None, :, None]**2*probability*(centered**2-variance[..., None]))
        curvature_derivative[:, self.linear_d:-1] += 2*sigma**2*(totals*variance*a*a) @ self.slope_design
        curvature_derivative[:, -1] += 2*(curvature-1)
        third = np.sum(probability*centered**3, axis=-1)
        curvature_z = sigma**3*((totals*third) @ (a**3))
        log_scale_derivative = -(curvature_derivative+curvature_z[:, None]*mode_derivative)/(2*curvature[:, None])

        gradient = np.zeros(self.dimension)
        direct_gradient = np.zeros(self.dimension)
        nll = 0.
        eap, sd = [], []
        # Batch Persons to bound transient node/cell/category arrays; no
        # Person/node/parameter tensor is needed for the first derivative.
        for start in range(0, len(totals), 32):
            sl = slice(start, start+32)
            z = mode[sl, None]+scale[sl, None]*nodes
            theta = sigma*z
            logits = a[None, None, :, None]*(base+theta[..., None, None]*self.k)
            logp = logits-logsumexp(logits, axis=-1, keepdims=True)
            probability = np.exp(logp)
            ll = np.einsum('nck,nqck->nq', self.counts[sl], logp, optimize=True)
            joint = ll+log_weight-.5*z*z+.5*nodes*nodes+np.log(scale[sl, None])
            marginal = logsumexp(joint, axis=1)
            posterior = np.exp(joint-marginal[:, None])
            nll -= float(marginal.sum())
            residual = self.counts[sl, None]-totals[sl, None, :, None]*probability
            weighted = residual*posterior[..., None, None]
            fixed = np.empty(self.dimension)
            fixed[:self.linear_d] = (weighted.sum(axis=(0, 1))*a[:, None]).ravel() @ self.design.reshape(-1, self.linear_d)
            fixed[self.linear_d:-1] = np.sum(weighted*logits, axis=(0, 1, 3)) @ self.slope_design
            eta_score = (residual @ self.k) @ a
            fixed[-1] = np.sum(posterior*theta*eta_score)
            node_score = sigma*eta_score-z
            mode_factor = np.sum(posterior*node_score, axis=1)
            scale_factor = np.sum(posterior*(node_score*scale[sl, None]*nodes+1), axis=1)
            direct_gradient -= fixed
            gradient -= fixed+mode_factor @ mode_derivative[sl]+scale_factor @ log_scale_derivative[sl]
            if moments:
                mean = np.sum(posterior*theta, axis=1)
                eap.extend(mean)
                sd.extend(np.sqrt(np.sum(posterior*(theta-mean[:, None])**2, axis=1)))
        if not np.isfinite(nll) or not np.isfinite(gradient).all():
            raise FloatingPointError('Nonfinite adaptive objective or gradient')
        result = dict(nll=nll, gradient=gradient, fixed_z_gradient=direct_gradient,
                      mode=mode, scale=scale, mode_score=score,
                      mode_iterations=iteration+1, cells=len(self.design), observations=self.observations)
        if moments:
            if not np.isfinite(eap).all() or not np.isfinite(sd).all() or np.any(np.asarray(sd) <= 0):
                raise FloatingPointError('Invalid adaptive posterior moments')
            result.update(eap=np.asarray(eap), sd=np.asarray(sd))
        return result
