"""Independent unit-weight RSM/PCM reference for pre-fit development problems.

Normal-population integrals use explicit affine category logits and Louis'
identity. No app probability, expansion, derivative or quadrature code is used.
Anchors, latent regression, penalties and GPCM are outside this reference.
"""

import numpy as np
from scipy.integrate import quad_vec
from scipy.optimize import brentq
from scipy.special import logsumexp, ndtr


class RaschIntegral:
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
        self.observed = self.design[np.arange(len(data)), self.y]
        self.rows = [np.flatnonzero(self.person == p) for p in range(self.n)]

    def conditional(self, point, z, *, information=False):
        theta = np.exp(point[-1]) * z
        logits = self.design @ point[:-1] + theta * self.k
        logp = logits - logsumexp(logits, axis=1, keepdims=True)
        probability = np.exp(logp)
        residual = self.y - probability @ self.k
        structural = self.observed - np.einsum("ik,ikd->id", probability, self.design)
        row_score = np.column_stack((structural, theta * residual))
        score = np.zeros((self.n, self.d+1)); np.add.at(score, self.person, row_score)
        loglik = np.bincount(self.person, weights=logp[np.arange(len(self.y)), self.y], minlength=self.n)
        eta_score = np.bincount(self.person, weights=residual, minlength=self.n)
        if not information:
            return loglik, score, eta_score, None
        derivative = np.concatenate((self.design, np.broadcast_to(theta*self.k[None, :, None],
                                    (len(self.y), len(self.k), 1))), axis=2)
        mean = np.einsum("ik,ikd->id", probability, derivative)
        centered = derivative - mean[:, None, :]
        row_h = np.einsum("ik,ikd,ike->ide", probability, centered, centered)
        row_h[:, -1, -1] -= theta * residual  # d² theta / d log(sigma)² = theta.
        hessian = np.zeros((self.n, self.d+1, self.d+1)); np.add.at(hessian, self.person, row_h)
        return loglik, score, eta_score, hessian

    def finite(self, point, nodes, weights):
        evaluated = [self.conditional(point, z, information=True) for z in nodes]
        lp = np.array([v[0] for v in evaluated]).T + np.log(weights)
        marginal = logsumexp(lp, axis=1)
        posterior = np.exp(lp-marginal[:, None])
        score = np.array([v[1] for v in evaluated]).transpose(1, 0, 2)
        mean = np.einsum("pq,pqd->pd", posterior, score)
        centered = score - mean[:, None, :]
        info = np.einsum("pq,qpde->de", posterior, np.array([v[3] for v in evaluated]))
        info -= np.einsum("pq,pqd,pqe->de", posterior, centered, centered)
        ez = posterior @ nodes
        variance = np.sum(posterior * (nodes[None, :]-ez[:, None])**2, axis=1)
        sigma = np.exp(point[-1])
        return dict(nll=float(-marginal.sum()), gradient=-mean.sum(axis=0), information=info,
                    eap=sigma*ez, sd=sigma*np.sqrt(variance))

    def continuous(self, point, *, bound=12., rel_tol=1e-10, abs_tol=1e-12, information=False):
        point = np.asarray(point, dtype=float)
        if point.shape != (self.d+1,) or not np.isfinite(point).all():
            raise ValueError("Invalid joint reference coordinates")
        sigma = np.exp(point[-1]); width = self.d+1
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
            return mass[:, None] * np.concatenate(blocks, axis=1)

        values, error, status = quad_vec(integrand, -bound, bound, epsrel=rel_tol, epsabs=abs_tol,
            norm="max", points=np.unique(np.r_[0., modes]), quadrature="gk21", workers=1, full_output=True)
        if not status.success or not np.isfinite(values).all() or np.any(values[:, 0] <= 0):
            raise FloatingPointError("Adaptive reference integration did not return valid positive masses")
        normalized = values / values[:, :1]
        log_mass = centers + np.log(values[:, 0])
        mean_score = normalized[:, 1:1+width]
        ez, ez2 = normalized[:, -2], normalized[:, -1]
        if np.any(ez2 <= ez*ez):
            raise FloatingPointError("Invalid posterior variance in reference integral")
        result = dict(nll=float(-log_mass.sum()), gradient=-mean_score.sum(axis=0),
            eap=sigma*ez, sd=sigma*np.sqrt(ez2-ez*ez), sigma=float(sigma),
            density_log_sigma_score=float(-np.sum(ez2-1)),
            numeric_relative_mass_error_sum=float(np.sum(error/values[:, 0])),
            tail_relative_mass_bound_sum=float(np.sum(np.exp(np.log(2*ndtr(-bound))-log_mass))),
            quadrature_neval=int(status.neval), quadrature_status=int(status.status))
        # Tail enclosures use likelihood <= 1 and bounded category-feature
        # ranges; Gaussian moments control the log-SD derivative's z factors.
        phi = np.exp(-bound*bound/2)/np.sqrt(2*np.pi)
        tail0, tail1, tail2 = 2*ndtr(-bound), 2*phi, 2*(bound*phi+ndtr(-bound))
        ranges = np.ptp(self.design, axis=1)
        structural_range = np.array([ranges[rows].sum(axis=0) for rows in self.rows])
        count = np.array([len(rows) for rows in self.rows]); max_score = len(self.k)-1
        eta_range = count*max_score
        score_tail = np.column_stack((structural_range*tail0, sigma*eta_range*tail1))
        tail_score_error = (score_tail + np.abs(mean_score)*tail0)*np.exp(-log_mass[:, None])
        if np.any(values[:, 0] <= error):
            raise FloatingPointError("Estimated integration error reaches an integrated mass")
        score_error = error*(1+np.abs(mean_score))/(values[:, :1]-error)
        result["tail_gradient_abs_bound"] = float(np.max(tail_score_error.sum(axis=0)))
        result["numeric_gradient_error_allowance"] = float(np.max(score_error.sum(axis=0)))
        if information:
            second = normalized[:, 1+width:1+width+width*width].reshape(self.n, width, width)
            result["information"] = np.sum(second + np.einsum("pi,pj->pij", mean_score, mean_score), axis=0)
            second_tail = np.zeros_like(second)
            second_tail[:, :-1, :-1] = np.array([r.T@r/4 for r in (ranges[rows] for rows in self.rows)])*tail0
            second_tail[:, :-1, :-1] += np.einsum("pi,pj->pij", structural_range, structural_range)*tail0
            second_tail[:, :-1, -1] = sigma*structural_range*(max_score/4+eta_range[:, None])*tail1
            second_tail[:, -1, :-1] = second_tail[:, :-1, -1]
            second_tail[:, -1, -1] = sigma**2*(count*max_score**2/4+eta_range**2)*tail2 + sigma*eta_range*tail1
            tail_second_error = (second_tail+np.abs(second)*tail0)*np.exp(-log_mass[:, None, None])
            second_error = error*(1+np.abs(second))/(values[:, :1, None]-error)
            def information_allowance(second_delta, score_delta):
                d = second_delta + np.einsum("pi,pj->pij", np.abs(mean_score), score_delta)
                d += np.einsum("pi,pj->pij", score_delta, np.abs(mean_score))
                d += np.einsum("pi,pj->pij", score_delta, score_delta)
                return float(np.max(d.sum(axis=0)))
            result["tail_information_abs_bound"] = information_allowance(tail_second_error, tail_score_error)
            # QUADPACK/quad_vec error is an estimate, not a rigorous enclosure.
            result["numeric_information_error_allowance"] = information_allowance(second_error, score_error)
        return result
