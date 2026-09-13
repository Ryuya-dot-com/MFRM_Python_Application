"""Development adapter from a native app MML fit to the v2 numerical kernel.

This adapter is intentionally not wired into the default Streamlit estimator.
It lets engineering fixtures exercise the same likelihood used by the app
without changing legacy results or assigning prospective scientific gates.
Only unregularized RSM/PCM free-population-SD fits are accepted in v2 phase 1.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
import operator
from typing import Any

import numpy as np

import streamlit_app as app
from mfrm_app.mml_engine_v2 import FreeSdStationarityRun, run_free_sd_stationarity_v2
from mfrm_app.mml_stationarity import JointPolishOptions


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


__all__ = [
    "AppFreeSdProblem",
    "prepare_app_free_sd_problem",
    "run_app_free_sd_stationarity_v2",
]
