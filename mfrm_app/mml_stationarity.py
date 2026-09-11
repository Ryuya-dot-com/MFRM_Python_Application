"""Streamlit-free numerical tools for joint free-SD MML stationarity.

This module does not define an MFRM likelihood and does not choose scientific
acceptance thresholds.  It supplies a deterministic numerical kernel that an
MML likelihood adapter can use to polish structural coordinates jointly with
``log(sigma)`` and to retain scale-explicit diagnostics.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable, Sequence

import numpy as np
from scipy.optimize import minimize


Array = np.ndarray
ValueFunction = Callable[[Array, float], float]
StructuralValueGradient = Callable[[Array, float], tuple[float, Array]]
JointValueGradient = Callable[[Array], tuple[float, Array]]


@dataclass(frozen=True)
class JointPolishOptions:
    """Deterministic optimizer settings; thresholds are intentionally absent."""

    maxiter: int = 500
    gtol: float = 1e-8
    ftol: float = 1e-15
    maxls: int = 50
    log_sigma_relative_step: float = 1e-5

    def validate(self) -> None:
        if isinstance(self.maxiter, bool) or int(self.maxiter) != self.maxiter or self.maxiter < 1:
            raise ValueError("maxiter must be a positive integer")
        if isinstance(self.maxls, bool) or int(self.maxls) != self.maxls or self.maxls < 1:
            raise ValueError("maxls must be a positive integer")
        for name in ("gtol", "ftol", "log_sigma_relative_step"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be positive and finite")


@dataclass(frozen=True)
class GradientAudit:
    coordinates: int
    analytical_supnorm: float
    finite_difference_supnorm: float
    maximum_absolute_difference: float
    analytical_gradient: tuple[float, ...]
    finite_difference_gradient: tuple[float, ...]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class InformationDiagnostics:
    coordinates: int
    finite: bool
    positive_definite: bool
    minimum_eigenvalue: float
    maximum_eigenvalue: float
    condition_number: float
    symmetry_residual_supnorm: float
    maximum_absolute_jacobian: float
    relative_symmetry_residual: float
    standardized_score_supnorm: float
    newton_correction_supnorm: float
    relative_step: float
    gradient_evaluations: int

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class JointPolishResult:
    structural_parameters: tuple[float, ...]
    sigma: float
    initial_joint_coordinates: tuple[float, ...]
    joint_coordinates: tuple[float, ...]
    initial_objective: float
    final_objective: float
    initial_scalar_objective: float
    final_scalar_objective: float
    objective_value_consistency_max_abs_difference: float
    objective_improvement: float
    maximum_coordinate_displacement: float
    final_gradient_supnorm: float
    final_projected_gradient_supnorm: float
    final_gradient_l2norm: float
    optimizer_success: bool
    optimizer_status: int
    optimizer_message: str
    optimizer_iterations: int
    function_evaluations: int
    gradient_evaluations: int
    objective_scale: str
    options: dict[str, object]
    gradient_method: str = "structural_analytic_log_sigma_central_difference"

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _vector(value: Sequence[float] | Array, label: str) -> Array:
    result = np.asarray(value, dtype=float)
    if result.ndim != 1 or result.size == 0 or not np.isfinite(result).all():
        raise ValueError(f"{label} must be a nonempty finite one-dimensional vector")
    return result


def _positive_sigma(value: float) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0:
        raise ValueError("sigma must be positive and finite")
    return result


def _relative_step(value: float, relative_step: float) -> float:
    relative_step = float(relative_step)
    if not np.isfinite(relative_step) or relative_step <= 0:
        raise ValueError("relative_step must be positive and finite")
    return relative_step * max(1.0, abs(float(value)))


def central_difference_gradient(
    objective: Callable[[Array], float],
    coordinates: Sequence[float] | Array,
    *,
    relative_step: float = 1e-5,
) -> Array:
    """Central-difference gradient with a deterministic coordinate step."""

    point = _vector(coordinates, "coordinates")
    gradient = np.empty_like(point)
    for index, value in enumerate(point):
        step = _relative_step(float(value), relative_step)
        plus = point.copy()
        minus = point.copy()
        plus[index] += step
        minus[index] -= step
        upper = float(objective(plus))
        lower = float(objective(minus))
        if not np.isfinite(upper) or not np.isfinite(lower):
            raise FloatingPointError(f"Nonfinite objective in coordinate {index}")
        gradient[index] = (upper - lower) / (2.0 * step)
    return gradient


def projected_gradient(
    coordinates: Sequence[float] | Array,
    gradient: Sequence[float] | Array,
    bounds: Sequence[tuple[float | None, float | None]],
    *,
    active_tolerance: float = 1e-10,
) -> Array:
    """Projected first-order residual for a box-constrained minimization."""

    point = _vector(coordinates, "coordinates")
    result = _vector(gradient, "gradient").copy()
    if result.size != point.size or len(bounds) != point.size:
        raise ValueError("coordinates, gradient, and bounds must have equal length")
    active_tolerance = float(active_tolerance)
    if not np.isfinite(active_tolerance) or active_tolerance < 0:
        raise ValueError("active_tolerance must be nonnegative and finite")
    for index, (lower, upper) in enumerate(bounds):
        if lower is not None and point[index] <= float(lower) + active_tolerance and result[index] > 0:
            result[index] = 0.0
        if upper is not None and point[index] >= float(upper) - active_tolerance and result[index] < 0:
            result[index] = 0.0
    return result


def make_joint_free_sd_functions(
    value_function: ValueFunction,
    structural_value_gradient: StructuralValueGradient,
    *,
    log_sigma_relative_step: float = 1e-5,
    joint_value_gradient: JointValueGradient | None = None,
) -> tuple[Callable[[Array], float], JointValueGradient]:
    """Adapt to ``[parameters, log(sigma)]``, optionally using a joint gradient.

    The supplied callback must differentiate the same finite-Q objective,
    including movement of the nodes with sigma. Otherwise log-SD is differenced.
    """

    _relative_step(0.0, log_sigma_relative_step)

    def objective(joint: Array) -> float:
        point = _vector(joint, "joint coordinates")
        sigma = _positive_sigma(float(np.exp(point[-1])))
        value = float(value_function(point[:-1], sigma))
        if not np.isfinite(value):
            raise FloatingPointError("Joint MML objective is nonfinite")
        return value

    def value_gradient(joint: Array) -> tuple[float, Array]:
        point = _vector(joint, "joint coordinates")
        sigma = _positive_sigma(float(np.exp(point[-1])))
        if joint_value_gradient is not None:
            value, gradient = joint_value_gradient(point)
            gradient = _vector(gradient, "joint gradient")
            if gradient.size != point.size:
                raise ValueError("Joint gradient dimension differs from joint coordinates")
            if not np.isfinite(float(value)):
                raise FloatingPointError("Joint MML objective is nonfinite")
            return float(value), gradient
        value, structural_gradient = structural_value_gradient(point[:-1], sigma)
        structural = _vector(structural_gradient, "structural gradient")
        if structural.size != point.size - 1:
            raise ValueError("Structural gradient dimension differs from structural coordinates")
        if not np.isfinite(float(value)):
            raise FloatingPointError("Joint MML objective is nonfinite")
        step = _relative_step(float(point[-1]), log_sigma_relative_step)
        plus = point.copy()
        minus = point.copy()
        plus[-1] += step
        minus[-1] -= step
        log_sigma_gradient = (objective(plus) - objective(minus)) / (2.0 * step)
        gradient = np.concatenate([structural, [log_sigma_gradient]])
        return float(value), gradient

    return objective, value_gradient


def audit_joint_gradient(
    objective: Callable[[Array], float],
    value_gradient: JointValueGradient,
    coordinates: Sequence[float] | Array,
    *,
    relative_step: float = 1e-5,
) -> GradientAudit:
    point = _vector(coordinates, "coordinates")
    _, analytical = value_gradient(point)
    analytical = _vector(analytical, "analytical gradient")
    numerical = central_difference_gradient(
        objective,
        point,
        relative_step=relative_step,
    )
    if analytical.size != numerical.size:
        raise ValueError("Analytical and finite-difference gradients differ in length")
    return GradientAudit(
        coordinates=int(point.size),
        analytical_supnorm=float(np.max(np.abs(analytical))),
        finite_difference_supnorm=float(np.max(np.abs(numerical))),
        maximum_absolute_difference=float(np.max(np.abs(analytical - numerical))),
        analytical_gradient=tuple(map(float, analytical)),
        finite_difference_gradient=tuple(map(float, numerical)),
    )


def information_diagnostics(
    value_gradient: JointValueGradient,
    coordinates: Sequence[float] | Array,
    *,
    relative_step: float = 1e-4,
) -> InformationDiagnostics:
    """Numerically evaluate local curvature without assigning a PASS threshold."""

    point = _vector(coordinates, "coordinates")
    relative_step = float(relative_step)
    _relative_step(0.0, relative_step)
    _, score = value_gradient(point)
    score = _vector(score, "gradient")
    if score.size != point.size:
        raise ValueError("Gradient dimension differs from coordinates")
    jacobian = np.empty((point.size, point.size), dtype=float)
    for index, value in enumerate(point):
        step = _relative_step(float(value), relative_step)
        plus = point.copy()
        minus = point.copy()
        plus[index] += step
        minus[index] -= step
        _, upper = value_gradient(plus)
        _, lower = value_gradient(minus)
        jacobian[:, index] = (np.asarray(upper) - np.asarray(lower)) / (2.0 * step)
    finite = bool(np.isfinite(jacobian).all())
    if not finite:
        return InformationDiagnostics(
            coordinates=int(point.size), finite=False, positive_definite=False,
            minimum_eigenvalue=np.nan, maximum_eigenvalue=np.nan,
            condition_number=np.inf, symmetry_residual_supnorm=np.inf,
            maximum_absolute_jacobian=np.inf, relative_symmetry_residual=np.inf,
            standardized_score_supnorm=np.inf, newton_correction_supnorm=np.inf,
            relative_step=relative_step,
            gradient_evaluations=int(1 + 2 * point.size),
        )
    information = 0.5 * (jacobian + jacobian.T)
    symmetry_residual = float(np.max(np.abs(jacobian - jacobian.T)))
    maximum_absolute_jacobian = float(np.max(np.abs(jacobian)))
    relative_symmetry_residual = float(
        symmetry_residual / max(1.0, maximum_absolute_jacobian)
    )
    eigenvalues = np.linalg.eigvalsh(information)
    minimum = float(np.min(eigenvalues))
    maximum = float(np.max(eigenvalues))
    positive_definite = bool(minimum > 0)
    condition = float(maximum / minimum) if positive_definite else np.inf
    diagonal = np.diag(information)
    standardized = (
        float(np.max(np.abs(score) / np.sqrt(diagonal)))
        if np.all(diagonal > 0)
        else np.inf
    )
    if positive_definite:
        correction = np.linalg.solve(information, score)
        newton_supnorm = float(np.max(np.abs(correction)))
    else:
        newton_supnorm = np.inf
    return InformationDiagnostics(
        coordinates=int(point.size),
        finite=True,
        positive_definite=positive_definite,
        minimum_eigenvalue=minimum,
        maximum_eigenvalue=maximum,
        condition_number=condition,
        symmetry_residual_supnorm=symmetry_residual,
        maximum_absolute_jacobian=maximum_absolute_jacobian,
        relative_symmetry_residual=relative_symmetry_residual,
        standardized_score_supnorm=standardized,
        newton_correction_supnorm=newton_supnorm,
        relative_step=relative_step,
        gradient_evaluations=int(1 + 2 * point.size),
    )


def polish_joint_free_sd(
    structural_start: Sequence[float] | Array,
    sigma_start: float,
    value_function: ValueFunction,
    structural_value_gradient: StructuralValueGradient,
    *,
    structural_bounds: Sequence[tuple[float | None, float | None]] | None = None,
    sigma_bounds: tuple[float, float] = (0.05, 10.0),
    objective_scale: str = "sum_negative_log_likelihood",
    options: JointPolishOptions | None = None,
    joint_value_gradient: JointValueGradient | None = None,
) -> JointPolishResult:
    """Jointly polish structural coordinates and ``log(sigma)`` once."""

    settings = options or JointPolishOptions()
    settings.validate()
    structural = _vector(structural_start, "structural_start")
    sigma = _positive_sigma(sigma_start)
    sigma_lower, sigma_upper = map(float, sigma_bounds)
    if not (
        np.isfinite(sigma_lower)
        and np.isfinite(sigma_upper)
        and 0 < sigma_lower < sigma_upper
        and sigma_lower <= sigma <= sigma_upper
    ):
        raise ValueError("sigma_bounds must contain sigma_start and be positive")
    if not objective_scale.strip():
        raise ValueError("objective_scale must be explicit")
    if structural_bounds is None:
        structural_bounds = [(None, None)] * structural.size
    if len(structural_bounds) != structural.size:
        raise ValueError("structural_bounds length differs from structural_start")
    bounds = [*structural_bounds, (np.log(sigma_lower), np.log(sigma_upper))]
    start = np.concatenate([structural, [np.log(sigma)]])
    objective, value_gradient = make_joint_free_sd_functions(
        value_function,
        structural_value_gradient,
        log_sigma_relative_step=settings.log_sigma_relative_step,
        joint_value_gradient=joint_value_gradient,
    )
    initial_value, _ = value_gradient(start)
    initial_scalar_value = objective(start)
    fitted = minimize(
        value_gradient,
        start,
        jac=True,
        method="L-BFGS-B",
        bounds=bounds,
        options={
            "maxiter": int(settings.maxiter),
            "gtol": float(settings.gtol),
            "ftol": float(settings.ftol),
            "maxls": int(settings.maxls),
        },
    )
    final_value, final_gradient = value_gradient(np.asarray(fitted.x, dtype=float))
    final_scalar_value = objective(np.asarray(fitted.x, dtype=float))
    projected = projected_gradient(fitted.x, final_gradient, bounds)
    return JointPolishResult(
        structural_parameters=tuple(map(float, fitted.x[:-1])),
        sigma=float(np.exp(fitted.x[-1])),
        initial_joint_coordinates=tuple(map(float, start)),
        joint_coordinates=tuple(map(float, fitted.x)),
        initial_objective=float(initial_value),
        final_objective=float(final_value),
        initial_scalar_objective=float(initial_scalar_value),
        final_scalar_objective=float(final_scalar_value),
        objective_value_consistency_max_abs_difference=float(
            max(
                abs(initial_value - initial_scalar_value),
                abs(final_value - final_scalar_value),
            )
        ),
        objective_improvement=float(initial_value - final_value),
        maximum_coordinate_displacement=float(np.max(np.abs(fitted.x - start))),
        final_gradient_supnorm=float(np.max(np.abs(final_gradient))),
        final_projected_gradient_supnorm=float(np.max(np.abs(projected))),
        final_gradient_l2norm=float(np.linalg.norm(final_gradient)),
        optimizer_success=bool(fitted.success),
        optimizer_status=int(fitted.status),
        optimizer_message=str(fitted.message),
        optimizer_iterations=int(fitted.nit),
        function_evaluations=int(fitted.nfev),
        gradient_evaluations=int(getattr(fitted, "njev", 0) or 0),
        objective_scale=objective_scale,
        options=asdict(settings),
        gradient_method=(
            "provided_joint_gradient" if joint_value_gradient is not None
            else "structural_analytic_log_sigma_central_difference"
        ),
    )


__all__ = [
    "GradientAudit",
    "InformationDiagnostics",
    "JointPolishOptions",
    "JointPolishResult",
    "audit_joint_gradient",
    "central_difference_gradient",
    "information_diagnostics",
    "make_joint_free_sd_functions",
    "polish_joint_free_sd",
    "projected_gradient",
]
