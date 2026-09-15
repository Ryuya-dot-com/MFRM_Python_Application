"""Development adapter from a native app MML fit to the v2 numerical kernel.

This adapter is intentionally not wired into the default Streamlit estimator.
It lets engineering fixtures exercise the same likelihood used by the app
without changing legacy results or assigning prospective scientific gates.
Free-SD v2 accepts unregularized RSM/PCM. A separate opt-in fixed-SD PCM
development entry point returns candidates without replacing a fitted result.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
import json
import math
import operator
from typing import Any, Callable

import numpy as np
from scipy.special import logsumexp

import streamlit_app as app
from mfrm_app.mml_engine_v2 import (
    FreeSdStationarityRun, TwoStageFreeSdRun,
    run_free_sd_stationarity_v2, run_free_sd_two_stage,
)
from mfrm_app.mml_stationarity import (
    JointPolishOptions, audit_joint_gradient, information_diagnostics, polish_fixed_sd,
)


def _log_average_exp(log_weights: np.ndarray, change: np.ndarray) -> np.ndarray:
    # Same arithmetic as the immutable mml_pcm_likelihood_difference study;
    # importing that executable also imports its historical study dependencies.
    if not np.isfinite(change).all():
        raise FloatingPointError("Nonfinite log ratio")
    log_weights = log_weights - logsumexp(log_weights, axis=-1, keepdims=True)
    if np.max(np.abs(change)) <= 0.5:
        return np.log1p(np.sum(np.exp(log_weights) * np.expm1(change), axis=-1))
    return logsumexp(log_weights + change, axis=-1) - logsumexp(log_weights, axis=-1)


@dataclass(frozen=True)
class AppFreeSdProblem:
    config: dict[str, Any]
    sizes: Any
    idx: dict[str, Any]
    structural_start: tuple[float, ...]
    sigma_start: float
    quadrature_points: int
    observations: float
    structural_bounds: tuple[tuple[float | None, float | None], ...]
    sigma_bounds: tuple[float, float]

    def value(self, parameters: np.ndarray, sigma: float) -> float:
        quad = app.gauss_hermite_normal(self.quadrature_points, sd=sigma)
        return float(
            app.mfrm_loglik_mml(parameters, self.idx, self.config, self.sizes, quad)
        )

    def value_gradient(
        self,
        parameters: np.ndarray,
        sigma: float,
    ) -> tuple[float, np.ndarray]:
        quad = app.gauss_hermite_normal(self.quadrature_points, sd=sigma)
        value, gradient = app.mfrm_loglik_mml_value_grad(
            parameters,
            self.idx,
            self.config,
            self.sizes,
            quad,
        )
        return float(value), np.asarray(gradient, dtype=float)

    def joint_value_gradient(self, joint: np.ndarray) -> tuple[float, np.ndarray]:
        """Analytic finite-GH gradient in [structural, log(sigma)] coordinates."""
        quad = app.gauss_hermite_normal(self.quadrature_points, sd=float(np.exp(joint[-1])))
        value, gradient = app.mfrm_loglik_mml_value_grad(
            joint[:-1], self.idx, self.config, self.sizes, quad, include_log_sigma=True,
        )
        return float(value), np.asarray(gradient, dtype=float)

    def likelihood_difference(self, anchor: np.ndarray) -> Callable[[np.ndarray], float]:
        """NLL(point)-NLL(anchor) for the app's weighted, observed-row GH target.

        Response weights multiply conditional log probabilities before person
        integration. Fractional weights therefore define a powered likelihood.
        Missing rows are absent; an empty or zero-weight person contributes zero.
        """
        config, idx = copy.deepcopy(self.config), copy.deepcopy(self.idx)
        if config["model"] not in {"RSM", "PCM"} or config["method"] != "MML":
            raise ValueError("Likelihood differences require RSM/PCM MML")
        if config.get("facet_regularization_enabled") or config.get("facet_regularization", {}).get("enabled"):
            raise ValueError("Likelihood differences exclude penalized fits")
        shape = (sum(self.sizes.values()) + 1,)
        log_bounds = np.log(self.sigma_bounds)
        fixed_sd = not config.get('estimate_population_sd', False)
        fixed_log_sd = float(np.log(self.sigma_start)) if fixed_sd else None

        def checked(point):
            point = np.asarray(point, dtype=float)
            if point.shape != shape or not np.isfinite(point).all():
                raise ValueError("Expected finite joint coordinates with the identified dimension")
            if not log_bounds[0] <= point[-1] <= log_bounds[1]:
                raise ValueError("Population SD is outside its bounds")
            if fixed_sd and point[-1] != fixed_log_sd:
                raise ValueError("Fixed-SD likelihood differences cannot change population SD")
            return point

        anchor = checked(anchor).copy()
        person, score = idx["person"], idx["score_k"]
        weight = np.asarray(idx.get("weight", np.ones(len(score))), dtype=float)
        n_person, n_cat = config["n_person"], config["n_cat"]
        if (weight.shape != score.shape or not np.isfinite(weight).all()
            or np.any(weight < 0) or not np.any(weight > 0)):
            raise ValueError("Response weights must be finite, nonnegative and not all zero")
        if (person.shape != score.shape or np.any(person < 0) or np.any(person >= n_person)
            or np.any(score < 0) or np.any(score >= n_cat)):
            raise ValueError("Invalid observed person/category indices")
        rows = [np.flatnonzero(person == p) for p in range(n_person)]
        saved_rows = idx.get("rows_by_person", rows)
        if len(saved_rows) != n_person or any(not np.array_equal(a, b) for a, b in zip(rows, saved_rows)):
            raise ValueError("Person row groups disagree with observation indices")
        quad = app.gauss_hermite_normal(self.quadrature_points)
        z, w = quad["nodes"], quad["weights"]
        if (not np.isfinite(z).all() or not np.isfinite(w).all() or not np.all(w > 0)
            or abs(w.sum() - 1) > 1e-12 or abs(w @ z) > 1e-12
            or abs(w @ (z*z) - 1) > 1e-12):
            raise ValueError("GH rule requires finite nodes and positive standard-normal weights")
        theta = np.exp(anchor[-1]) * z
        if fixed_sd:
            # Preserve the exact fixed rule, including SD roundoff. The joint
            # helper is used only on its structural slice; SD is never fitted.
            theta = app.gauss_hermite_normal(self.quadrature_points, sd=self.sigma_start)['nodes']
        categories = np.arange(n_cat)

        # Expand changes with homogeneous constraints. Subtracting two expanded
        # vectors would lose tiny changes against nonzero anchors/group means.
        change_config = copy.deepcopy(config)
        for spec in change_config["facet_specs"].values():
            spec["anchors"][np.isfinite(spec["anchors"])] = 0.0
            spec["group_values"] = dict.fromkeys(spec["group_values"], 0.0)

        def offsets(parameters, cfg):
            expanded = app.expand_params(parameters, self.sizes, cfg)
            eta = app.compute_base_eta(idx, expanded, cfg) + app.compute_population_mu(expanded, cfg)[person]
            steps = expanded["steps"][None, :] if cfg["model"] == "RSM" else expanded["steps_mat"]
            cumulative = np.column_stack((np.zeros(len(steps)), np.cumsum(steps, axis=1)))
            cumulative = cumulative[0] if cfg["model"] == "RSM" else cumulative[idx["step_idx"]]
            return eta[:, None] * categories - cumulative

        # ponytail: cache O(rows * Q * categories); chunk by person if large-data
        # profiling shows this development adapter exceeds the memory budget.
        logits = offsets(anchor[:-1], config)[:, None, :] + theta[None, :, None] * categories
        log_category = logits - logsumexp(logits, axis=-1, keepdims=True)
        observed = (np.arange(len(score)), slice(None), score)
        log_posterior = np.tile(np.log(w), (n_person, 1))
        np.add.at(log_posterior, person, weight[:, None] * log_category[observed])

        def difference(point):
            change = checked(point) - anchor
            theta_change = theta * np.expm1(change[-1])
            logit_change = offsets(change[:-1], change_config)[:, None, :] + theta_change[None, :, None] * categories
            normalizer_change = _log_average_exp(log_category, logit_change)
            person_change = np.zeros((n_person, len(w)))
            np.add.at(person_change, person, weight[:, None] * (logit_change[observed] - normalizer_change))
            return -math.fsum(_log_average_exp(log_posterior, person_change))

        return difference

    def constraint_residual(self, parameters: np.ndarray) -> float:
        """Reconstruct exact expanded-coordinate identification constraints."""

        expanded = app.expand_params(parameters, self.sizes, self.config)
        residuals: list[float] = []
        for facet in self.config.get("facet_names", []):
            values = np.asarray(expanded["facets"][facet], dtype=float)
            spec = self.config["facet_specs"][facet]
            anchors = np.asarray(spec["anchors"], dtype=float)
            anchored = np.isfinite(anchors)
            if anchored.any():
                residuals.extend(np.abs(values[anchored] - anchors[anchored]).tolist())
            groups = list(spec.get("groups", []))
            for group in sorted({item for item in groups if item not in (None, "")}):
                members = np.array([item == group for item in groups], dtype=bool)
                target = float(spec.get("group_values", {}).get(group, 0.0))
                residuals.append(abs(float(np.mean(values[members])) - target))
            if bool(spec.get("centered", False)):
                residuals.append(abs(float(np.sum(values))))
        audit = app.build_parameterization_audit(self.config, self.sizes, expanded)
        if not audit.empty:
            audit_residuals = np.asarray(audit["ConstraintMaxAbsResidual"], dtype=float)
            residuals.extend(np.abs(audit_residuals).tolist())
            if not bool(audit["ExpandedBoundsSatisfied"].all()):
                return float("inf")
        if not residuals:
            return 0.0
        maximum = float(np.max(residuals))
        return maximum if np.isfinite(maximum) else float("inf")


def prepare_app_free_sd_problem(
    result: dict[str, Any],
    *,
    quadrature_points: int | None = None,
) -> AppFreeSdProblem:
    """Validate and adapt one already-fitted native free-SD MML result."""

    if not isinstance(result, dict):
        raise TypeError("result must be a native mfrm_estimate result dictionary")
    config_source = result.get("config")
    prep = result.get("prep")
    opt = result.get("opt")
    if not isinstance(config_source, dict) or not isinstance(prep, dict) or opt is None:
        raise ValueError("result lacks config, prep, or optimizer evidence")
    config = copy.deepcopy(config_source)
    if str(config.get("method", "")).upper() != "MML":
        raise ValueError("stationarity v2 requires an MML fit")
    if not bool(config.get("estimate_population_sd", False)):
        raise ValueError("stationarity v2 phase 1 requires free population SD")
    if str(config.get("model", "")).upper() not in {"RSM", "PCM"}:
        raise ValueError("stationarity v2 phase 1 supports RSM and PCM only")
    if bool(config.get("facet_regularization_enabled", False)) or bool(
        config.get("facet_regularization", {}).get("enabled", False)
    ):
        raise ValueError("stationarity v2 phase 1 excludes penalized fits")

    start = np.asarray(getattr(opt, "x", None), dtype=float)
    if start.ndim != 1 or start.size == 0 or not np.isfinite(start).all():
        raise ValueError("optimizer structural start is missing or nonfinite")
    summary = result.get("summary")
    if summary is None or not hasattr(summary, "iloc") or len(summary) != 1:
        raise ValueError("result lacks its one-row summary evidence")
    sigma_sources = {
        "config": config.get("estimated_population_sd"),
        "optimizer": getattr(opt, "estimated_population_sd", None),
        "summary": summary.iloc[0].get("EstimatedPopulationSD"),
    }
    try:
        sigma_values = {name: float(value) for name, value in sigma_sources.items()}
    except (TypeError, ValueError) as exc:
        raise ValueError("estimated population SD evidence is incomplete") from exc
    if any(not np.isfinite(value) or value <= 0 for value in sigma_values.values()):
        raise ValueError("estimated population SD evidence is nonpositive or nonfinite")
    sigma = sigma_values["config"]
    if any(
        abs(value - sigma) > 1e-12 * max(1.0, abs(value), abs(sigma))
        for value in sigma_values.values()
    ):
        raise ValueError("config, optimizer, and summary population SD evidence differs")
    raw_q = quadrature_points if quadrature_points is not None else config.get("quad_points")
    if isinstance(raw_q, (bool, np.bool_)):
        raise ValueError("quadrature_points must be an integer of at least 3")
    try:
        q = operator.index(raw_q)
    except TypeError as exc:
        raise ValueError("quadrature_points must be an integer of at least 3") from exc
    if q < 3:
        raise ValueError("quadrature_points must be an integer of at least 3")

    sizes = app.build_param_sizes(config)
    if int(sum(sizes.values())) != int(start.size):
        raise ValueError("optimizer coordinates do not match the identified parameterization")
    idx = app.build_indices(
        prep,
        step_facet=config.get("step_facet"),
        slope_facet=config.get("slope_facet"),
    )
    for value in idx.values():
        if isinstance(value, np.ndarray):
            value.setflags(write=False)
        elif isinstance(value, list):
            for item in value:
                if isinstance(item, np.ndarray):
                    item.setflags(write=False)
    weight = idx.get("weight")
    observations = float(np.sum(weight)) if weight is not None else float(len(idx["score_k"]))
    if not np.isfinite(observations) or observations <= 0:
        raise ValueError("fit has no positive finite observation count")
    bounds = app.build_optimizer_bounds(sizes, config)
    if bounds is None:
        bounds = [(None, None)] * start.size
    sigma_bounds = tuple(map(float, config.get("population_sd_bounds", (0.05, 10.0))))
    return AppFreeSdProblem(
        config=config,
        sizes=sizes,
        idx=idx,
        structural_start=tuple(map(float, start)),
        sigma_start=sigma,
        quadrature_points=q,
        observations=observations,
        structural_bounds=tuple(bounds),
        sigma_bounds=sigma_bounds,
    )


def run_app_free_sd_stationarity_v2(
    result: dict[str, Any],
    *,
    quadrature_points: int | None = None,
    options: JointPolishOptions | None = None,
) -> tuple[AppFreeSdProblem, FreeSdStationarityRun]:
    """Run the threshold-free v2 kernel on the app's exact finite-Q objective."""

    problem = prepare_app_free_sd_problem(result, quadrature_points=quadrature_points)
    run = run_free_sd_stationarity_v2(
        problem.structural_start,
        problem.sigma_start,
        problem.value,
        problem.value_gradient,
        observations=problem.observations,
        structural_bounds=problem.structural_bounds,
        sigma_bounds=problem.sigma_bounds,
        options=options,
        constraint_residual_function=problem.constraint_residual,
        joint_value_gradient=problem.joint_value_gradient,
    )
    return problem, run


def run_app_free_sd_two_stage(
    result: dict[str, Any],
    *,
    anchor_gradient_limit: float,
    refinement_options: JointPolishOptions,
    preliminary_options: JointPolishOptions | None = None,
    quadrature_points: int | None = None,
) -> tuple[AppFreeSdProblem, TwoStageFreeSdRun]:
    """Opt-in observed-row adapter; no default estimator or inference gate changes."""
    problem = prepare_app_free_sd_problem(result, quadrature_points=quadrature_points)
    run = run_free_sd_two_stage(
        problem.structural_start, problem.sigma_start, problem.value, problem.value_gradient,
        observations=problem.observations, difference_factory=problem.likelihood_difference,
        joint_value_gradient=problem.joint_value_gradient, anchor_gradient_limit=anchor_gradient_limit,
        refinement_options=refinement_options, preliminary_options=preliminary_options,
        structural_bounds=problem.structural_bounds, sigma_bounds=problem.sigma_bounds,
        constraint_residual_function=problem.constraint_residual,
    )
    return problem, run


def run_app_fixed_sd_polish(result, *, enabled=False, options=None):
    """Opt-in candidate/diagnostics for the restricted PCM; no in-place update.

    Disabled means no validation, likelihood evaluation, export or optimization.
    Enabled reuses the native-comparison model preflight, but runs Python only.
    Candidate scores contain EAP and posterior SD, not recalculated SEs or CIs.
    """
    if type(enabled) is not bool:
        raise ValueError('enabled must be a boolean')
    if not enabled:
        return None
    prepared = app.prepare_fixed_sd_pcm_check(result)
    cfg, prep, sizes, idx = (prepared[k] for k in ('config', 'prep', 'sizes', 'idx'))
    start, quad = prepared['coordinates'], prepared['quad']
    sigma, q = float(quad['sd']), len(quad['nodes'])
    assets = app.build_cross_engine_validation_bundle(result)
    problem = AppFreeSdProblem(cfg, sizes, idx, tuple(start), sigma, q,
        float(len(idx['score_k'])), ((None, None),)*len(start), (.05, 10.))
    value_gradient = lambda p: problem.value_gradient(p, sigma)
    initial_check = app.evaluate_fixed_sd_pcm_check(prepared)
    value, gradient = initial_check['raw_nll'], np.asarray(initial_check['gradient'])
    reported = prepared['reported_nll']
    log_sd = float(np.log(sigma))
    def factory(anchor):
        delta = problem.likelihood_difference(np.r_[anchor, log_sd])
        return lambda point: delta(np.r_[point, log_sd])
    run = polish_fixed_sd(start, value_gradient, factory, options=options)
    run.update(schema='mfrm_fixed_sd_candidate_v1', source='development_api',
        input_bundle_sha256=json.loads(assets['comparison_report.json'])['details']['input_bundle_sha256'],
        coordinate_blocks=dict(sizes), facet_levels=copy.deepcopy(cfg['facet_levels']),
        fixed_quadrature=dict(method='Gauss-Hermite', points=q, sigma=sigma,
            nodes=quad['nodes'].tolist(), weights=quad['weights'].tolist()),
        original_optimizer=dict(message=str(result['opt'].message), success=bool(result['opt'].success),
            reported_nll=reported, reevaluated_nll=float(value), gradient=gradient.tolist()),
        candidate_is_fitted_result=False)
    if run['completed']:
        point = np.asarray(run['candidate_coordinates'])
        run['information'] = [information_diagnostics(value_gradient, point, relative_step=h).to_dict() for h in (1e-4, 3e-5)]
        run['gradient_audit'] = audit_joint_gradient(lambda p:value_gradient(p)[0], value_gradient, point).to_dict()
        expanded = app.expand_params(point, sizes, cfg)
        score = app.compute_person_eap(idx, cfg, expanded, quad)
        score.insert(0, 'Person', prep['levels']['Person'])
        run['candidate_parameters'] = expanded
        run['candidate_person_scores'] = score
        run['constraint_residual'] = problem.constraint_residual(point)
    return run


__all__ = [
    "AppFreeSdProblem",
    "prepare_app_free_sd_problem",
    "run_app_free_sd_stationarity_v2",
    "run_app_free_sd_two_stage",
    "run_app_fixed_sd_polish",
]
