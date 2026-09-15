"""Threshold-explicit orchestration for stationarity-audited free-SD MML.

The legacy EM estimator is intentionally outside this module.  Callers may use
it as a warm start, then pass the resulting structural coordinates and latent
population SD here.  This module performs two joint optimizations of the
finite-quadrature objective in ``[structural parameters, log(sigma)]`` space.

Scientific qualification is deliberately separate from optimization.  The
engine records diagnostics without silently choosing tolerances.  This module
cannot emit ``InferenceReady=True``; only a future artifact-bound study-level
registrar may combine a frozen contract with independent Q31/Q61 evidence.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable, Sequence

import numpy as np

from mfrm_app.mml_stationarity import (
    GradientAudit,
    InformationDiagnostics,
    JointPolishOptions,
    JointPolishResult,
    JointValueGradient,
    StructuralValueGradient,
    ValueFunction,
    audit_joint_gradient,
    information_diagnostics,
    make_joint_free_sd_functions,
    polish_joint_free_sd,
    projected_gradient,
)


@dataclass(frozen=True)
class StationarityContract:
    """Predeclared numerical gates for a finite-quadrature MML solution."""

    max_projected_gradient_supnorm: float
    max_standardized_score_supnorm: float
    max_newton_correction_supnorm: float
    max_restart_improvement_per_observation: float
    max_restart_displacement: float
    max_gradient_fd_disagreement: float
    max_objective_value_disagreement: float
    max_information_relative_symmetry_residual: float
    max_information_condition_number: float
    max_constraint_residual: float
    objective_worsening_tolerance_total: float
    objective_worsening_tolerance_per_observation: float
    sigma_boundary_tolerance: float

    def validate(self) -> None:
        positive = (
            "max_projected_gradient_supnorm",
            "max_standardized_score_supnorm",
            "max_newton_correction_supnorm",
            "max_restart_improvement_per_observation",
            "max_restart_displacement",
            "max_gradient_fd_disagreement",
            "max_objective_value_disagreement",
            "max_information_relative_symmetry_residual",
            "max_information_condition_number",
            "max_constraint_residual",
        )
        for name in positive:
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be positive and finite")
        for name in (
            "objective_worsening_tolerance_total",
            "objective_worsening_tolerance_per_observation",
            "sigma_boundary_tolerance",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be nonnegative and finite")


@dataclass(frozen=True)
class FreeSdStationarityRun:
    """A threshold-free two-polish numerical record."""

    primary_polish: JointPolishResult
    restart_polish: JointPolishResult
    gradient_audit_h: GradientAudit
    gradient_audit_h_over_2: GradientAudit
    information: InformationDiagnostics
    observations: float
    sigma_lower_bound: float
    sigma_upper_bound: float
    final_sigma: float
    final_objective: float
    final_projected_gradient_supnorm: float
    gradient_step_agreement_supnorm: float
    restart_improvement_per_observation: float
    restart_displacement: float
    constraint_residual: float | None
    algorithm_terminated: bool
    finite_solution: bool
    objective_nonworsening: bool

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class TwoStageFreeSdRun:
    """Keep preliminary failures separate from any later refinement evidence."""

    preliminary: FreeSdStationarityRun
    refinement: FreeSdStationarityRun | None
    anchor_admitted: bool
    anchor_gradient_limit: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class StationarityAssessment:
    """Fail-closed interpretation of a run under an explicit contract."""

    algorithm_terminated: bool
    finite_solution: bool
    objective_nonworsening: bool
    sigma_interior: bool
    information_finite: bool
    information_positive_definite: bool
    information_symmetry_pass: bool
    information_condition_pass: bool
    constraint_pass: bool
    objective_value_consistency_pass: bool
    projected_gradient_pass: bool
    standardized_score_pass: bool
    newton_correction_pass: bool
    gradient_fd_agreement_pass: bool
    restart_improvement_pass: bool
    restart_displacement_pass: bool
    stationarity_pass: bool
    quadrature_sensitivity_pass: None
    inference_ready: bool
    status: str
    contract: dict[str, float]

    def to_dict(self) -> dict[str, object]:
        result = asdict(self)
        result.update(
            {
                "AlgorithmTerminated": self.algorithm_terminated,
                "StationarityPass": self.stationarity_pass,
                "QuadratureSensitivityPass": self.quadrature_sensitivity_pass,
                "InferenceReady": self.inference_ready,
            }
        )
        return result


def _validated_observations(value: float) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0:
        raise ValueError("observations must be positive and finite")
    return result


def _validated_sigma_bounds(bounds: tuple[float, float]) -> tuple[float, float]:
    lower, upper = map(float, bounds)
    if not (np.isfinite(lower) and np.isfinite(upper) and 0 < lower < upper):
        raise ValueError("sigma_bounds must be finite, positive, and ordered")
    return lower, upper


def _joint_bounds(
    structural_count: int,
    structural_bounds: Sequence[tuple[float | None, float | None]] | None,
    sigma_bounds: tuple[float, float],
) -> list[tuple[float | None, float | None]]:
    if structural_bounds is None:
        structural = [(None, None)] * structural_count
    else:
        structural = list(structural_bounds)
        if len(structural) != structural_count:
            raise ValueError("structural_bounds length differs from structural coordinates")
    lower, upper = _validated_sigma_bounds(sigma_bounds)
    return [*structural, (float(np.log(lower)), float(np.log(upper)))]


def run_free_sd_stationarity_v2(
    structural_start: Sequence[float] | np.ndarray,
    sigma_start: float,
    value_function: ValueFunction,
    structural_value_gradient: StructuralValueGradient,
    *,
    observations: float,
    structural_bounds: Sequence[tuple[float | None, float | None]] | None = None,
    sigma_bounds: tuple[float, float] = (0.05, 10.0),
    options: JointPolishOptions | None = None,
    constraint_residual_function: Callable[[np.ndarray], float] | None = None,
    joint_value_gradient: JointValueGradient | None = None,
    objective_difference: Callable[[np.ndarray], float] | None = None,
    objective_anchor: Sequence[float] | np.ndarray | None = None,
) -> FreeSdStationarityRun:
    """Polish once, restart once, and retain scale-explicit diagnostics."""

    n_observations = _validated_observations(observations)
    lower, upper = _validated_sigma_bounds(sigma_bounds)
    settings = options or JointPolishOptions()
    settings.validate()
    normalized_bounds = _joint_bounds(
        len(np.asarray(structural_start, dtype=float).reshape(-1)),
        structural_bounds,
        (lower, upper),
    )
    if any(bound != (None, None) for bound in normalized_bounds[:-1]):
        raise ValueError(
            "stationarity v2 phase 1 excludes bounded structural coordinates"
        )

    primary = polish_joint_free_sd(
        structural_start,
        sigma_start,
        value_function,
        structural_value_gradient,
        structural_bounds=structural_bounds,
        sigma_bounds=(lower, upper),
        options=settings,
        joint_value_gradient=joint_value_gradient,
        objective_difference=objective_difference,
        objective_anchor=objective_anchor,
    )
    restart = polish_joint_free_sd(
        primary.structural_parameters,
        primary.sigma,
        value_function,
        structural_value_gradient,
        structural_bounds=structural_bounds,
        sigma_bounds=(lower, upper),
        options=settings,
        joint_value_gradient=joint_value_gradient,
        objective_difference=objective_difference,
        objective_anchor=objective_anchor,
    )

    final = np.asarray(restart.joint_coordinates, dtype=float)
    objective_h, value_gradient_h = make_joint_free_sd_functions(
        value_function,
        structural_value_gradient,
        log_sigma_relative_step=settings.log_sigma_relative_step,
        joint_value_gradient=joint_value_gradient,
    )
    half_step = settings.log_sigma_relative_step / 2.0
    _, value_gradient_half = make_joint_free_sd_functions(
        value_function,
        structural_value_gradient,
        log_sigma_relative_step=half_step,
        joint_value_gradient=joint_value_gradient,
    )
    # Only the scalar audit changes origin; its derivative is still the raw gradient.
    if objective_difference is not None:
        objective_h = objective_difference
    audit_h = audit_joint_gradient(
        objective_h,
        value_gradient_h,
        final,
        relative_step=settings.log_sigma_relative_step,
    )
    audit_half = audit_joint_gradient(
        objective_h,
        value_gradient_half,
        final,
        relative_step=half_step,
    )
    _, gradient_h = value_gradient_h(final)
    _, gradient_half = value_gradient_half(final)
    joint_bounds = _joint_bounds(len(primary.structural_parameters), structural_bounds, (lower, upper))
    projected = projected_gradient(final, gradient_half, joint_bounds)
    information = information_diagnostics(value_gradient_half, final)
    constraint_residual = None
    if constraint_residual_function is not None:
        raw_residual = float(constraint_residual_function(final[:-1]))
        constraint_residual = raw_residual if np.isfinite(raw_residual) else None

    finite_solution = bool(
        np.isfinite(final).all()
        and np.isfinite(restart.final_objective)
        and np.isfinite(gradient_h).all()
        and np.isfinite(gradient_half).all()
    )
    algorithm_terminated = bool(
        primary.optimizer_success and restart.optimizer_success and finite_solution
    )
    objective_nonworsening = bool(
        primary.optimization_improvement >= 0
        and restart.optimization_improvement >= 0
    )
    return FreeSdStationarityRun(
        primary_polish=primary,
        restart_polish=restart,
        gradient_audit_h=audit_h,
        gradient_audit_h_over_2=audit_half,
        information=information,
        observations=n_observations,
        sigma_lower_bound=lower,
        sigma_upper_bound=upper,
        final_sigma=float(restart.sigma),
        final_objective=float(restart.final_objective),
        final_projected_gradient_supnorm=float(np.max(np.abs(projected))),
        gradient_step_agreement_supnorm=float(np.max(np.abs(gradient_h - gradient_half))),
        restart_improvement_per_observation=float(
            restart.optimization_improvement / n_observations
        ),
        restart_displacement=float(restart.maximum_coordinate_displacement),
        constraint_residual=constraint_residual,
        algorithm_terminated=algorithm_terminated,
        finite_solution=finite_solution,
        objective_nonworsening=objective_nonworsening,
    )


def run_free_sd_two_stage(
    structural_start: Sequence[float] | np.ndarray,
    sigma_start: float,
    value_function: ValueFunction,
    structural_value_gradient: StructuralValueGradient,
    *,
    observations: float,
    difference_factory: Callable[[np.ndarray], Callable[[np.ndarray], float]],
    joint_value_gradient: JointValueGradient,
    anchor_gradient_limit: float,
    refinement_options: JointPolishOptions,
    preliminary_options: JointPolishOptions | None = None,
    structural_bounds: Sequence[tuple[float | None, float | None]] | None = None,
    sigma_bounds: tuple[float, float] = (0.05, 10.0),
    constraint_residual_function: Callable[[np.ndarray], float] | None = None,
) -> TwoStageFreeSdRun:
    """Two raw polishes, then two fixed-anchor polishes when admission is met.

    A failed preliminary algorithm can leave a usable endpoint. Keep that failure
    in ``preliminary`` even if ``refinement`` later passes stationarity assessment.
    This function assigns neither stationarity nor scientific inference gates.
    """
    limit = float(anchor_gradient_limit)
    if not np.isfinite(limit) or limit <= 0:
        raise ValueError("anchor_gradient_limit must be positive and finite")
    refinement_options.validate()
    if not callable(difference_factory) or not callable(joint_value_gradient):
        raise ValueError("two-stage refinement requires difference and joint-gradient callbacks")
    common = dict(observations=observations, structural_bounds=structural_bounds,
                  sigma_bounds=sigma_bounds, constraint_residual_function=constraint_residual_function,
                  joint_value_gradient=joint_value_gradient)
    preliminary = run_free_sd_stationarity_v2(
        structural_start, sigma_start, value_function, structural_value_gradient,
        options=preliminary_options, **common,
    )
    endpoint = preliminary.restart_polish
    admitted = bool(preliminary.finite_solution and sigma_bounds[0] < endpoint.sigma < sigma_bounds[1]
                    and endpoint.final_gradient_supnorm <= limit)
    refinement = None
    if admitted:
        anchor = np.array(endpoint.joint_coordinates)
        difference = difference_factory(anchor.copy())
        refinement = run_free_sd_stationarity_v2(
            endpoint.structural_parameters, endpoint.sigma, value_function, structural_value_gradient,
            options=refinement_options, objective_difference=difference, objective_anchor=anchor, **common,
        )
    return TwoStageFreeSdRun(preliminary, refinement, admitted, limit)


def assess_free_sd_stationarity(
    run: FreeSdStationarityRun,
    contract: StationarityContract,
) -> StationarityAssessment:
    """Apply stationarity gates without claiming full inferential readiness.

    Quadrature sensitivity is a separate, artifact-bound study-level contract.
    A bare boolean is deliberately not accepted here; until that layer exists,
    this function always returns ``InferenceReady=False``.
    """

    contract.validate()
    def consistent(left: float, right: float) -> bool:
        left = float(left)
        right = float(right)
        return bool(
            np.isfinite(left)
            and np.isfinite(right)
            and abs(left - right) <= 1e-12 * max(1.0, abs(left), abs(right))
        )

    def validate_polish(result: JointPolishResult, label: str) -> tuple[np.ndarray, np.ndarray]:
        initial = np.asarray(result.initial_joint_coordinates, dtype=float)
        final = np.asarray(result.joint_coordinates, dtype=float)
        structural = np.asarray(result.structural_parameters, dtype=float)
        if (
            initial.ndim != 1
            or final.ndim != 1
            or structural.ndim != 1
            or initial.size == 0
            or initial.shape != final.shape
            or structural.size != final.size - 1
            or not np.isfinite(initial).all()
            or not np.isfinite(final).all()
            or not np.isfinite(structural).all()
        ):
            raise ValueError(f"{label} polish coordinate evidence is invalid")
        checks = (
            np.allclose(structural, final[:-1], rtol=0.0, atol=1e-12),
            consistent(result.sigma, float(np.exp(final[-1]))),
            consistent(
                result.objective_improvement,
                result.initial_objective - result.final_objective,
            ),
            consistent(
                result.maximum_coordinate_displacement,
                float(np.max(np.abs(final - initial))),
            ),
            consistent(
                result.objective_value_consistency_max_abs_difference,
                max(
                    abs(result.initial_objective - result.initial_scalar_objective),
                    abs(result.final_objective - result.final_scalar_objective),
                ),
            ),
        )
        if not all(checks):
            raise ValueError(f"{label} polish derived fields do not reconstruct")
        shift = result.objective_shift
        if shift is not None:
            anchor = np.asarray(shift.anchor_coordinates, dtype=float)
            values = np.array([shift.anchor_objective, shift.initial_difference,
                               shift.final_difference, shift.reconstruction_max_abs_difference])
            if (anchor.shape != final.shape or not np.isfinite(anchor).all()
                or not np.isfinite(values).all() or shift.reconstruction_max_abs_difference < 0):
                raise ValueError(f"{label} objective origin evidence is invalid")
            endpoint_error = max(
                abs(shift.anchor_objective + shift.initial_difference - result.initial_objective),
                abs(shift.anchor_objective + shift.final_difference - result.final_objective),
            )
            if endpoint_error > shift.reconstruction_max_abs_difference:
                raise ValueError(f"{label} objective origin does not reconstruct")
        return initial, final

    def validate_gradient_audit(audit: GradientAudit, label: str) -> tuple[np.ndarray, np.ndarray]:
        analytical = np.asarray(audit.analytical_gradient, dtype=float)
        numerical = np.asarray(audit.finite_difference_gradient, dtype=float)
        if (
            analytical.ndim != 1
            or analytical.size == 0
            or analytical.shape != numerical.shape
            or not np.isfinite(analytical).all()
            or not np.isfinite(numerical).all()
        ):
            raise ValueError(f"{label} gradient audit evidence is invalid")
        checks = (
            audit.coordinates == analytical.size,
            consistent(audit.analytical_supnorm, float(np.max(np.abs(analytical)))),
            consistent(
                audit.finite_difference_supnorm,
                float(np.max(np.abs(numerical))),
            ),
            consistent(
                audit.maximum_absolute_difference,
                float(np.max(np.abs(analytical - numerical))),
            ),
        )
        if not all(checks):
            raise ValueError(f"{label} gradient audit derived fields do not reconstruct")
        return analytical, numerical

    observations = _validated_observations(run.observations)
    sigma_lower, sigma_upper = _validated_sigma_bounds(
        (run.sigma_lower_bound, run.sigma_upper_bound)
    )
    primary = run.primary_polish
    restart = run.restart_polish
    primary_initial, primary_final = validate_polish(primary, "primary")
    restart_initial, final_coordinates = validate_polish(restart, "restart")
    left, right = primary.objective_shift, restart.objective_shift
    if ((left is None) != (right is None) or (left is not None and (
        left.anchor_coordinates != right.anchor_coordinates
        or not consistent(left.anchor_objective, right.anchor_objective)
        or not consistent(left.final_difference, right.initial_difference)
    ))):
        raise ValueError("stationarity objective origin lineage is inconsistent")
    audit_h, _ = validate_gradient_audit(run.gradient_audit_h, "h")
    audit_half, _ = validate_gradient_audit(run.gradient_audit_h_over_2, "h_over_2")
    if (
        audit_h.shape != final_coordinates.shape
        or audit_half.shape != final_coordinates.shape
        or not np.allclose(restart_initial, primary_final, rtol=0.0, atol=1e-12)
    ):
        raise ValueError("stationarity run coordinate lineage is inconsistent")

    derived_finite = bool(
        np.isfinite(final_coordinates).all()
        and np.isfinite(restart.final_objective)
        and np.isfinite(audit_h).all()
        and np.isfinite(audit_half).all()
    )
    derived_algorithm_terminated = bool(
        primary.optimizer_success and restart.optimizer_success and derived_finite
    )
    derived_objective_nonworsening = bool(
        primary.optimization_improvement >= 0
        and restart.optimization_improvement >= 0
    )
    derived_restart_improvement = float(
        primary.final_objective - restart.final_objective if right is None else
        left.final_difference - right.final_difference
    )
    derived_restart_per_observation = float(derived_restart_improvement / observations)
    derived_restart_displacement = float(
        np.max(np.abs(final_coordinates - primary_final))
    )
    derived_gradient_step_agreement = float(np.max(np.abs(audit_h - audit_half)))
    phase_one_bounds = [
        *[(None, None)] * (len(final_coordinates) - 1),
        (float(np.log(sigma_lower)), float(np.log(sigma_upper))),
    ]
    derived_projected_gradient_supnorm = float(
        np.max(
            np.abs(
                projected_gradient(
                    final_coordinates,
                    audit_half,
                    phase_one_bounds,
                )
            )
        )
    )
    consistency_checks = (
        run.algorithm_terminated == derived_algorithm_terminated,
        run.finite_solution == derived_finite,
        run.objective_nonworsening == derived_objective_nonworsening,
        consistent(run.final_sigma, restart.sigma),
        consistent(run.final_objective, restart.final_objective),
        consistent(run.restart_displacement, derived_restart_displacement),
        consistent(run.restart_improvement_per_observation, derived_restart_per_observation),
        consistent(run.gradient_step_agreement_supnorm, derived_gradient_step_agreement),
        consistent(
            run.final_projected_gradient_supnorm,
            derived_projected_gradient_supnorm,
        ),
    )
    if not all(consistency_checks):
        raise ValueError("stationarity run derived fields do not reconstruct")
    sigma_margin = float(contract.sigma_boundary_tolerance)
    sigma_interior = bool(
        run.final_sigma > sigma_lower + sigma_margin
        and run.final_sigma < sigma_upper - sigma_margin
    )
    information_condition_pass = bool(
        run.information.finite
        and np.isfinite(run.information.condition_number)
        and run.information.condition_number <= contract.max_information_condition_number
    )
    information_symmetry_pass = bool(
        np.isfinite(run.information.relative_symmetry_residual)
        and run.information.relative_symmetry_residual
        <= contract.max_information_relative_symmetry_residual
    )
    constraint_pass = bool(
        run.constraint_residual is not None
        and np.isfinite(run.constraint_residual)
        and run.constraint_residual <= contract.max_constraint_residual
    )
    projected_gradient_pass = bool(
        np.isfinite(derived_projected_gradient_supnorm)
        and derived_projected_gradient_supnorm
        <= contract.max_projected_gradient_supnorm
    )
    objective_value_consistency_pass = bool(
        primary.objective_value_consistency_max_abs_difference
        <= contract.max_objective_value_disagreement
        and all(p.objective_shift is None or p.objective_shift.reconstruction_max_abs_difference
                <= contract.max_objective_value_disagreement for p in (primary, restart))
        and restart.objective_value_consistency_max_abs_difference
        <= contract.max_objective_value_disagreement
    )
    standardized_score_pass = bool(
        np.isfinite(run.information.standardized_score_supnorm)
        and run.information.standardized_score_supnorm <= contract.max_standardized_score_supnorm
    )
    newton_correction_pass = bool(
        np.isfinite(run.information.newton_correction_supnorm)
        and run.information.newton_correction_supnorm <= contract.max_newton_correction_supnorm
    )
    gradient_fd_agreement_pass = bool(
        np.isfinite(run.gradient_step_agreement_supnorm)
        and run.gradient_step_agreement_supnorm <= contract.max_gradient_fd_disagreement
        and run.gradient_audit_h.maximum_absolute_difference
        <= contract.max_gradient_fd_disagreement
        and run.gradient_audit_h_over_2.maximum_absolute_difference
        <= contract.max_gradient_fd_disagreement
    )
    restart_improvement_pass = bool(
        np.isfinite(derived_restart_per_observation)
        and derived_restart_per_observation
        <= contract.max_restart_improvement_per_observation
        and derived_restart_per_observation
        >= -contract.objective_worsening_tolerance_per_observation
    )
    restart_displacement_pass = bool(
        np.isfinite(run.restart_displacement)
        and run.restart_displacement <= contract.max_restart_displacement
    )
    objective_nonworsening = bool(
        derived_objective_nonworsening
        or (
            primary.optimization_improvement
            >= -contract.objective_worsening_tolerance_total
            and restart.optimization_improvement
            >= -contract.objective_worsening_tolerance_total
        )
    )
    stationarity_pass = all(
        (
            derived_algorithm_terminated,
            derived_finite,
            objective_nonworsening,
            sigma_interior,
            run.information.finite,
            run.information.positive_definite,
            information_symmetry_pass,
            information_condition_pass,
            constraint_pass,
            objective_value_consistency_pass,
            projected_gradient_pass,
            standardized_score_pass,
            newton_correction_pass,
            gradient_fd_agreement_pass,
            restart_improvement_pass,
            restart_displacement_pass,
        )
    )
    inference_ready = False
    status = (
        "STATIONARITY_PASS_QUADRATURE_EVIDENCE_REQUIRED"
        if stationarity_pass
        else "STATIONARITY_NOT_QUALIFIED"
    )
    return StationarityAssessment(
        algorithm_terminated=derived_algorithm_terminated,
        finite_solution=derived_finite,
        objective_nonworsening=objective_nonworsening,
        sigma_interior=sigma_interior,
        information_finite=run.information.finite,
        information_positive_definite=run.information.positive_definite,
        information_symmetry_pass=information_symmetry_pass,
        information_condition_pass=information_condition_pass,
        constraint_pass=constraint_pass,
        objective_value_consistency_pass=objective_value_consistency_pass,
        projected_gradient_pass=projected_gradient_pass,
        standardized_score_pass=standardized_score_pass,
        newton_correction_pass=newton_correction_pass,
        gradient_fd_agreement_pass=gradient_fd_agreement_pass,
        restart_improvement_pass=restart_improvement_pass,
        restart_displacement_pass=restart_displacement_pass,
        stationarity_pass=stationarity_pass,
        quadrature_sensitivity_pass=None,
        inference_ready=inference_ready,
        status=status,
        contract={key: float(value) for key, value in asdict(contract).items()},
    )


__all__ = [
    "FreeSdStationarityRun",
    "StationarityAssessment",
    "StationarityContract",
    "TwoStageFreeSdRun",
    "assess_free_sd_stationarity",
    "run_free_sd_stationarity_v2",
    "run_free_sd_two_stage",
]
