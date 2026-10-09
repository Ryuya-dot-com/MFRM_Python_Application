"""Independent continuous-normal GPCM likelihood and observed information.

Reuses only the reference PCM category-design construction. Derivatives keep
theta fixed and differentiate the normal density for log(sigma), unlike the
app's moving finite-GH nodes. Unit weights, simple native constraints only.
No app expansion, probability, score, information or quadrature functions.
"""

import numpy as np
from scipy.integrate import quad_vec
from scipy.optimize import brentq
from scipy.special import logsumexp, ndtr

from validation.mml_rasch_integral import RaschIntegral


class GpcmIntegral:
    def __init__(self, data, config):
        if config["model"] != "GPCM" or config.get("slope_facet") not in config["facet_names"]:
            raise ValueError("Reference requires GPCM with an explicit slope facet")
        base = RaschIntegral(data, dict(config, model="PCM"))
        self.design, self.y, self.k = base.design, base.y, base.k
        self.person, self.rows, self.n = base.person, base.rows, base.n
        self.person_labels = base.person_labels
        self.linear_d = base.d
        levels = list(config["facet_levels"][config["slope_facet"]])
        index = np.array([levels.index(str(v)) for v in data[config["slope_facet"]]])
        basis = np.r_[np.eye(len(levels)-1), -np.ones((1, len(levels)-1))]
        self.slope_design = basis[index]
        self.d = self.linear_d + len(levels)-1

    def _point(self, point):
        point = np.asarray(point, dtype=float)
        if point.shape != (self.d+1,) or not np.isfinite(point).all():
            raise ValueError("Invalid joint reference coordinates")
        slope = np.exp(self.slope_design @ point[self.linear_d:-1])
        sigma = np.exp(point[-1])
        if not np.isfinite(slope).all() or not np.isfinite(sigma) or sigma <= 0 or np.any(slope <= 0):
            raise ValueError("Nonfinite or nonpositive reference scale")
        return point, slope, sigma

    def conditional(self, point, z, *, information=False):
        point, slope, sigma = self._point(point)
        logits = slope[:, None] * (self.design @ point[:self.linear_d] + sigma*z*self.k)
        logp = logits - logsumexp(logits, axis=1, keepdims=True)
        p = np.exp(logp)
        d = np.concatenate((slope[:, None, None]*self.design,
                            logits[:, :, None]*self.slope_design[:, None, :]), axis=2)
        mean = np.einsum("ik,ikd->id", p, d)
        row_score = mean - d[np.arange(len(self.y)), self.y]
        score = np.zeros((self.n, self.d+1))
        np.add.at(score[:, :-1], self.person, row_score)
        score[:, -1] = 1-z*z  # Normal-density NLL derivative at fixed theta.
        ll = np.bincount(self.person, weights=logp[np.arange(len(self.y)), self.y], minlength=self.n)
        eta_score = np.bincount(self.person, weights=slope*(self.y-p@self.k), minlength=self.n)
        hessian = None
        if information:
            centered = d - mean[:, None, :]
            row_h = np.einsum("ik,ikd,ike->ide", p, centered, centered, optimize=True)
            # E[d² logit]-d² logit(y) is essential for nonlinear slopes.
            cross = row_score[:, :self.linear_d, None]*self.slope_design[:, None, :]
            row_h[:, :self.linear_d, self.linear_d:] += cross
            row_h[:, self.linear_d:, :self.linear_d] += cross.transpose(0, 2, 1)
            logit_residual = np.sum(p*logits, axis=1)-logits[np.arange(len(self.y)), self.y]
            row_h[:, self.linear_d:, self.linear_d:] += np.einsum(
                "i,ij,ik->ijk", logit_residual, self.slope_design, self.slope_design)
            hessian = np.zeros((self.n, self.d+1, self.d+1))
            np.add.at(hessian[:, :-1, :-1], self.person, row_h)
            hessian[:, -1, -1] = 2*z*z
        return ll, score, eta_score, hessian

    def _tail_numerators(self, point, bound):
        """Absolute polynomial bounds for score and H-score*score outside B.

        For any category feature, |E(feature)-feature(y)| <= its range;
        |Cov(f,g)| <= range(f)*range(g)/4. Slopes introduce degree-one
        ranges in |z|; the normal-density score has degree two. L(theta)<=1.
        """
        point, slope, sigma = self._point(point)
        width = self.d+1
        ranges = np.zeros((len(self.y), self.d, 2))
        ranges[:, :self.linear_d, 0] = slope[:, None]*np.ptp(self.design, axis=1)
        logit_range = slope*np.ptp(self.design @ point[:self.linear_d], axis=1)
        ranges[:, self.linear_d:, 0] = logit_range[:, None]*abs(self.slope_design)
        ranges[:, self.linear_d:, 1] = (slope*sigma*self.k[-1])[:, None]*abs(self.slope_design)
        score = np.zeros((self.n, width, 3))
        np.add.at(score[:, :-1, :2], self.person, ranges)
        score[:, -1, 0] = 1.; score[:, -1, 2] = 1.
        row_h = np.zeros((len(self.y), self.d, self.d, 3))
        for j in range(2):
            for k in range(2):
                row_h[:, :, :, j+k] += np.einsum("id,ie->ide", ranges[:, :, j], ranges[:, :, k])/4
        cross = ranges[:, :self.linear_d, 0, None]*abs(self.slope_design[:, None, :])
        row_h[:, :self.linear_d, self.linear_d:, 0] += cross
        row_h[:, self.linear_d:, :self.linear_d, 0] += cross.transpose(0, 2, 1)
        for degree, value in enumerate([logit_range, slope*sigma*self.k[-1]]):
            row_h[:, self.linear_d:, self.linear_d:, degree] += np.einsum(
                "i,ij,ik->ijk", value, abs(self.slope_design), abs(self.slope_design))
        second = np.zeros((self.n, width, width, 5))
        np.add.at(second[:, :-1, :-1, :3], self.person, row_h)
        second[:, -1, -1, 2] = 2.
        for j in range(3):
            for k in range(3):
                second[:, :, :, j+k] += np.einsum("pd,pe->pde", score[:, :, j], score[:, :, k])
        phi = np.exp(-bound*bound/2)/np.sqrt(2*np.pi)
        moments = [2*ndtr(-bound), 2*phi]
        for power in range(2, 5):
            moments.append(2*bound**(power-1)*phi+(power-1)*moments[power-2])
        return score @ moments[:3], second @ moments, moments[0]

    def continuous(self, point, *, bound=12., rel_tol=1e-10, abs_tol=1e-12, information=False):
        point, _, sigma = self._point(point)
        if not np.isfinite(bound) or bound <= 0 or min(rel_tol, abs_tol) <= 0:
            raise ValueError("Invalid integration settings")
        width = self.d+1
        modes = np.array([brentq(lambda z: sigma*self.conditional(point, z)[2][p]-z,
                                -bound, bound, xtol=1e-13) for p in range(self.n)])
        if np.max(abs(modes)) >= bound-.1:
            raise ValueError("Posterior mode reaches reference integration boundary")
        centers = np.array([self.conditional(point, z)[0][p]-z*z/2-np.log(2*np.pi)/2
                            for p, z in enumerate(modes)])

        def integrand(z):
            ll, score, _, complete_h = self.conditional(point, z, information=information)
            mass = np.exp(ll-z*z/2-np.log(2*np.pi)/2-centers)
            blocks = [np.ones((self.n, 1)), score]
            if information:
                blocks.append((complete_h-np.einsum("pi,pj->pij", score, score)).reshape(self.n, -1))
            blocks.extend([np.full((self.n, 1), z), np.full((self.n, 1), z*z)])
            return mass[:, None]*np.concatenate(blocks, axis=1)

        values, error, status = quad_vec(integrand, -bound, bound, epsrel=rel_tol, epsabs=abs_tol,
            norm="max", points=np.unique(np.r_[0., modes]), quadrature="gk21", workers=1, full_output=True)
        if not status.success or not np.isfinite(values).all() or np.any(values[:, 0] <= error):
            raise FloatingPointError("Adaptive integration error reaches a reference mass or integration failed")
        normalized = values / values[:, :1]
        log_mass = centers + np.log(values[:, 0])
        mean = normalized[:, 1:1+width]
        ez, ez2 = normalized[:, -2], normalized[:, -1]
        if np.any(ez2 <= ez*ez):
            raise FloatingPointError("Invalid posterior variance in reference integral")
        tail_score, tail_second, tail0 = self._tail_numerators(point, bound)
        tail_score_error = (tail_score+abs(mean)*tail0)*np.exp(-log_mass[:, None])
        numeric_score_error = error*(1+abs(mean))/(values[:, :1]-error)
        result = dict(nll=float(-log_mass.sum()), gradient=mean.sum(axis=0),
            eap=sigma*ez, sd=sigma*np.sqrt(ez2-ez*ez), sigma=float(sigma),
            density_log_sigma_score=float(np.sum(1-ez2)),
            numeric_relative_mass_error_sum=float(np.sum(error/values[:, 0])),
            tail_relative_mass_bound_sum=float(np.sum(tail0*np.exp(-log_mass))),
            tail_gradient_abs_bound=float(np.max(tail_score_error.sum(axis=0))),
            numeric_gradient_error_allowance=float(np.max(numeric_score_error.sum(axis=0))),
            quadrature_neval=int(status.neval), quadrature_status=int(status.status))
        if information:
            second = normalized[:, 1+width:1+width+width*width].reshape(self.n, width, width)
            result["information"] = np.sum(second+np.einsum("pi,pj->pij", mean, mean), axis=0)
            tail_error = (tail_second+abs(second)*tail0)*np.exp(-log_mass[:, None, None])
            numeric_error = error*(1+abs(second))/(values[:, :1, None]-error)
            def allowance(second_delta, score_delta):
                total = second_delta+np.einsum("pi,pj->pij", abs(mean), score_delta)
                total += np.einsum("pi,pj->pij", score_delta, abs(mean))+np.einsum("pi,pj->pij", score_delta, score_delta)
                return float(np.max(total.sum(axis=0)))
            result["tail_information_abs_bound"] = allowance(tail_error, tail_score_error)
            # quad_vec's error is estimated, not a rigorous enclosure.
            result["numeric_information_error_allowance"] = allowance(numeric_error, numeric_score_error)
        return result
