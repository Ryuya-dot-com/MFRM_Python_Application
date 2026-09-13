"""Exact conditional maximum-likelihood calibration for Rasch-family MFRMs.

This module is deliberately Streamlit-free.  It implements the Phase 0
conditional likelihood core for additive RSM and rectangular PCM models.  The
person parameter is removed by conditioning on each person's raw score.  The
remaining structural effects use explicit sum-to-zero coordinates and an
explicit, caller-supplied rating support; observed maxima are never used to
truncate a virtual response unit's category support.

The public surface is intentionally small:

``prepare_cmle_design``
    Validate long-format data, construct the conditional design, group equal
    response-unit patterns, and perform a conditional-information rank audit.

``cmle_objective_value_grad``
    Evaluate the exact negative conditional log-likelihood and its analytical
    gradient by scaled dynamic programming.

``cmle_objective_value_grad_hessian``
    Additionally return the exact observed conditional information from the
    conditional sufficient-statistic covariance, without numerical
    differentiation.

``fit_cmle``
    Optimize an eligible design and return app-friendly structural tables,
    covariance evidence, convergence evidence, and category surfaces.

Phase 0 exclusions are fail-closed: estimated discriminations/GPCM, non-unit
row weights, unlabelled duplicate response units, latent regression,
regularization, and person estimation are not silently approximated.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import minimize


SUPPORTED_CMLE_MODELS = frozenset({"RSM", "PCM"})
CMLE_SCHEMA_VERSION = "native_exact_cmle_phase0_v1"
DEFAULT_RANK_AUDIT_MAX_BYTES = 512 * 1024 * 1024


class CMLEEligibilityError(ValueError):
    """Raised when a requested fit does not pass the CMLE pre-fit audit."""


@dataclass
class CMLEPattern:
    """Aggregated persons sharing the same multiset of response-unit designs."""

    signature: tuple[tuple[int, ...], ...]
    design: np.ndarray
    offset: np.ndarray
    score_frequencies: np.ndarray
    observed_statistics: np.ndarray
    observed_offset: float
    persons: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class CMLEStructuralParameter:
    """One expanded facet or threshold coordinate and its free-coordinate row."""

    parameter_type: str
    facet: str
    level: str
    step: int | None
    jacobian: np.ndarray
    offset: float = 0.0
    anchored: bool = False
    constraint: str = "sum_to_zero"


@dataclass
class CMLEDesign:
    """Prepared exact conditional-likelihood design."""

    data: pd.DataFrame
    person_col: str
    facet_cols: tuple[str, ...]
    score_col: str
    model: str
    step_facet: str | None
    rating_min: int
    rating_max: int
    n_categories: int
    facet_levels: dict[str, list[str]]
    facet_signs: dict[str, int]
    parameter_names: list[str]
    structural_parameters: list[CMLEStructuralParameter]
    row_design: np.ndarray
    row_offset: np.ndarray
    unit_codes: list[tuple[int, ...]]
    patterns: list[CMLEPattern]
    person_status: pd.DataFrame
    audit: dict[str, object]
    hard_anchors: pd.DataFrame
    response_unit_col: str | None = None

    @property
    def n_parameters(self) -> int:
        return len(self.parameter_names)


@dataclass
class _CMLEOptimizerEvaluationState:
    """Finite-domain evidence retained around the primary scipy optimizer."""

    invalid_evaluations: int = 0
    finite_evaluations: int = 0
    best_value: float = np.inf
    best_parameters: np.ndarray | None = None
    last_invalid_reason: str = ""


def _issue(code: str, severity: str, message: str, action: str) -> dict[str, str]:
    return {
        "Code": str(code),
        "Severity": str(severity),
        "Message": str(message),
        "Action": str(action),
    }


def _sum_zero_contrast(n_levels: int) -> np.ndarray:
    """Map ``n-1`` free coordinates to ``n`` effects summing exactly to zero."""
    n_levels = int(n_levels)
    if n_levels <= 1:
        return np.zeros((max(n_levels, 0), 0), dtype=float)
    out = np.zeros((n_levels, n_levels - 1), dtype=float)
    out[:-1, :] = np.eye(n_levels - 1, dtype=float)
    out[-1, :] = -1.0
    return out


def _finite_integer(value: object, name: str) -> int:
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite integer.") from exc
    if not np.isfinite(numeric) or abs(numeric - round(numeric)) > 1e-10:
        raise ValueError(f"{name} must be a finite integer.")
    return int(round(numeric))


def _normalize_hard_anchors(
    hard_anchors: pd.DataFrame | Sequence[Mapping[str, object]] | None,
    *,
    facet_cols: Sequence[str],
    facet_levels: dict[str, list[str]],
) -> pd.DataFrame:
    """Validate and canonicalize Phase-B fixed facet anchors."""
    columns = ["ParameterType", "Facet", "Level", "Step", "Value"]
    if hard_anchors is None:
        return pd.DataFrame(columns=columns)
    if isinstance(hard_anchors, pd.DataFrame):
        anchors = hard_anchors.copy()
    elif isinstance(hard_anchors, Sequence) and not isinstance(
        hard_anchors, (str, bytes)
    ):
        rows = list(hard_anchors)
        if not all(isinstance(row, Mapping) for row in rows):
            raise ValueError("hard_anchors sequence entries must be mappings.")
        anchors = pd.DataFrame(rows)
    else:
        raise ValueError("hard_anchors must be a pandas DataFrame or sequence of mappings.")
    if anchors.empty and len(anchors.columns) == 0:
        return pd.DataFrame(columns=columns)

    required = {"ParameterType", "Facet", "Level", "Value"}
    allowed = required | {"Step"}
    missing = sorted(required - set(anchors.columns))
    extra = sorted(set(anchors.columns) - allowed)
    if missing:
        raise ValueError("hard_anchors is missing columns: " + ", ".join(missing))
    if extra:
        raise ValueError("hard_anchors has unexpected columns: " + ", ".join(extra))
    if anchors[["ParameterType", "Facet", "Level", "Value"]].isna().any().any():
        raise ValueError("hard_anchors required fields cannot be missing.")
    parameter_types = anchors["ParameterType"].astype(str).str.strip().str.casefold()
    if not parameter_types.eq("facet").all():
        raise ValueError("Native hard-anchor v1 supports ParameterType='Facet' only.")
    if "Step" in anchors.columns and anchors["Step"].notna().any():
        raise ValueError("Facet hard anchors must not specify Step values.")

    anchors = anchors.assign(
        ParameterType="Facet",
        Facet=anchors["Facet"].astype(str),
        Level=anchors["Level"].astype(str),
        Value=pd.to_numeric(anchors["Value"], errors="coerce"),
    )
    if not np.isfinite(anchors["Value"].to_numpy(dtype=float)).all():
        raise ValueError("hard_anchors Value must be finite numeric.")
    if anchors.duplicated(["Facet", "Level"], keep=False).any():
        raise ValueError("hard_anchors contains duplicate Facet-Level keys.")
    unknown_facets = sorted(set(anchors["Facet"]) - set(facet_cols))
    if unknown_facets:
        raise ValueError("hard_anchors contains unknown facets: " + ", ".join(unknown_facets))
    unknown_levels = [
        f"{row.Facet}={row.Level}"
        for row in anchors.itertuples(index=False)
        if row.Level not in facet_levels[row.Facet]
    ]
    if unknown_levels:
        raise ValueError("hard_anchors contains unknown levels: " + ", ".join(unknown_levels))

    facet_order = {facet: index for index, facet in enumerate(facet_cols)}
    level_order = {
        (facet, level): index
        for facet in facet_cols
        for index, level in enumerate(facet_levels[facet])
    }
    anchors["__facet_order__"] = anchors["Facet"].map(facet_order)
    anchors["__level_order__"] = [
        level_order[(facet, level)]
        for facet, level in zip(anchors["Facet"], anchors["Level"])
    ]
    anchors = anchors.sort_values(
        ["__facet_order__", "__level_order__"], kind="stable"
    ).reset_index(drop=True)
    anchors["Step"] = pd.NA
    return anchors[columns]


def _parameter_layout(
    *,
    facet_levels: dict[str, list[str]],
    facet_cols: Sequence[str],
    facet_signs: dict[str, int],
    model: str,
    step_facet: str | None,
    n_categories: int,
    hard_anchor_values: Mapping[tuple[str, str], float],
) -> tuple[
    list[str],
    dict[str, tuple[slice, np.ndarray, np.ndarray]],
    dict[str, tuple[slice, np.ndarray]],
    list[CMLEStructuralParameter],
]:
    """Create free-coordinate slices and expanded-parameter Jacobian rows."""
    parameter_names: list[str] = []
    facet_blocks: dict[str, tuple[slice, np.ndarray, np.ndarray]] = {}
    step_blocks: dict[str, tuple[slice, np.ndarray]] = {}
    structural: list[CMLEStructuralParameter] = []
    cursor = 0

    for facet in facet_cols:
        levels = facet_levels[facet]
        anchored_levels = {
            level: float(hard_anchor_values[(facet, level)])
            for level in levels
            if (facet, level) in hard_anchor_values
        }
        if anchored_levels:
            unanchored_levels = [level for level in levels if level not in anchored_levels]
            contrast = np.zeros((len(levels), len(unanchored_levels)), dtype=float)
            unanchored_index = {level: index for index, level in enumerate(unanchored_levels)}
            for level_index, level in enumerate(levels):
                if level in unanchored_index:
                    contrast[level_index, unanchored_index[level]] = 1.0
            offsets = np.array(
                [anchored_levels.get(level, 0.0) for level in levels], dtype=float
            )
            free_labels = unanchored_levels
        else:
            contrast = _sum_zero_contrast(len(levels))
            offsets = np.zeros(len(levels), dtype=float)
            free_labels = levels[: contrast.shape[1]]
        width = contrast.shape[1]
        block_slice = slice(cursor, cursor + width)
        facet_blocks[facet] = (block_slice, contrast, offsets)
        parameter_names.extend(
            f"facet:{facet}:free:{level}" for level in free_labels
        )
        cursor += width

    n_steps = n_categories - 1
    step_contrast = _sum_zero_contrast(n_steps)
    if model == "RSM":
        step_levels = ["__shared__"]
    else:
        if step_facet is None:
            raise ValueError("PCM requires step_facet.")
        step_levels = facet_levels[step_facet]
    for level in step_levels:
        width = step_contrast.shape[1]
        block_slice = slice(cursor, cursor + width)
        step_blocks[level] = (block_slice, step_contrast)
        parameter_names.extend(
            f"step:{level}:free:{j + 1}" for j in range(width)
        )
        cursor += width

    n_parameters = cursor
    for facet in facet_cols:
        block_slice, contrast, offsets = facet_blocks[facet]
        facet_has_anchor = any(
            (facet, level) in hard_anchor_values for level in facet_levels[facet]
        )
        for level_index, level in enumerate(facet_levels[facet]):
            jac = np.zeros(n_parameters, dtype=float)
            jac[block_slice] = contrast[level_index]
            anchored = (facet, level) in hard_anchor_values
            structural.append(
                CMLEStructuralParameter(
                    parameter_type="Facet",
                    facet=facet,
                    level=level,
                    step=None,
                    jacobian=jac,
                    offset=float(offsets[level_index]),
                    anchored=anchored,
                    constraint=(
                        "hard_anchor"
                        if anchored
                        else "free_relative_to_hard_anchor"
                        if facet_has_anchor
                        else "sum_to_zero"
                    ),
                )
            )
    step_label = "Step" if model == "RSM" else str(step_facet)
    for level in step_levels:
        block_slice, contrast = step_blocks[level]
        for step_index in range(n_steps):
            jac = np.zeros(n_parameters, dtype=float)
            jac[block_slice] = contrast[step_index]
            structural.append(
                CMLEStructuralParameter(
                    parameter_type="Step",
                    facet=step_label,
                    level=level,
                    step=step_index + 1,
                    jacobian=jac,
                    offset=0.0,
                    anchored=False,
                    constraint="sum_to_zero",
                )
            )
    return parameter_names, facet_blocks, step_blocks, structural


def _build_row_design(
    *,
    data: pd.DataFrame,
    facet_cols: Sequence[str],
    facet_levels: dict[str, list[str]],
    facet_signs: dict[str, int],
    facet_blocks: dict[str, tuple[slice, np.ndarray, np.ndarray]],
    step_blocks: dict[str, tuple[slice, np.ndarray]],
    model: str,
    step_facet: str | None,
    n_categories: int,
    n_parameters: int,
) -> tuple[np.ndarray, np.ndarray, list[tuple[int, ...]]]:
    """Return per-row/per-category free-coordinate kernel design matrices."""
    n_rows = len(data)
    design = np.zeros((n_rows, n_categories, n_parameters), dtype=float)
    offset = np.zeros((n_rows, n_categories), dtype=float)
    level_maps = {
        facet: {level: index for index, level in enumerate(facet_levels[facet])}
        for facet in facet_cols
    }
    unit_codes: list[tuple[int, ...]] = []
    for row_number in range(n_rows):
        row = data.iloc[row_number]
        codes = tuple(level_maps[facet][str(row[facet])] for facet in facet_cols)
        unit_codes.append(codes)
        for category in range(n_categories):
            for facet, code in zip(facet_cols, codes):
                block_slice, contrast, facet_offsets = facet_blocks[facet]
                design[row_number, category, block_slice] += (
                    float(category) * float(facet_signs[facet]) * contrast[code]
                )
                offset[row_number, category] += (
                    float(category) * float(facet_signs[facet]) * facet_offsets[code]
                )
            step_level = (
                "__shared__" if model == "RSM" else str(row[str(step_facet)])
            )
            block_slice, contrast = step_blocks[step_level]
            if category > 0 and contrast.shape[1] > 0:
                design[row_number, category, block_slice] -= np.sum(
                    contrast[:category, :], axis=0
                )
    return design, offset, unit_codes


def _build_patterns(
    *,
    data: pd.DataFrame,
    person_col: str,
    score_internal_col: str,
    row_design: np.ndarray,
    row_offset: np.ndarray,
    unit_codes: Sequence[tuple[int, ...]],
    n_categories: int,
    n_parameters: int,
    response_unit_col: str | None,
) -> tuple[list[CMLEPattern], pd.DataFrame]:
    """Aggregate informative persons by observed response-unit design pattern."""
    pattern_map: dict[tuple[tuple[int, ...], ...], CMLEPattern] = {}
    status_rows: list[dict[str, object]] = []
    max_category = n_categories - 1

    for person, person_frame in data.groupby(person_col, sort=True, observed=True):
        row_indices = person_frame.index.to_numpy(dtype=int)

        def order_key(row_index: int) -> tuple[object, ...]:
            key: tuple[object, ...] = tuple(unit_codes[row_index])
            if response_unit_col is not None:
                key += (str(data.at[row_index, response_unit_col]),)
            return key

        ordered = np.array(sorted(row_indices.tolist(), key=order_key), dtype=int)
        scores = data.loc[ordered, score_internal_col].to_numpy(dtype=int)
        n_units = len(ordered)
        total_score = int(np.sum(scores))
        maximum_score = int(n_units * max_category)
        if n_units < 2:
            status = "insufficient_units"
        elif total_score == 0:
            status = "extreme_low"
        elif total_score == maximum_score:
            status = "extreme_high"
        else:
            status = "informative"
        status_rows.append(
            {
                "Person": str(person),
                "ResponseUnits": n_units,
                "RawScore": total_score,
                "MaximumScore": maximum_score,
                "ConditionalStatus": status,
            }
        )
        if status != "informative":
            continue

        signature = tuple(unit_codes[index] for index in ordered)
        person_design = row_design[ordered, :, :]
        person_offset = row_offset[ordered, :]
        observed = np.sum(
            person_design[np.arange(n_units), scores, :], axis=0
        )
        observed_offset = float(
            np.sum(person_offset[np.arange(n_units), scores])
        )
        if signature not in pattern_map:
            pattern_map[signature] = CMLEPattern(
                signature=signature,
                design=person_design.copy(),
                offset=person_offset.copy(),
                score_frequencies=np.zeros(maximum_score + 1, dtype=float),
                observed_statistics=np.zeros(n_parameters, dtype=float),
                observed_offset=0.0,
                persons=[],
            )
        pattern = pattern_map[signature]
        pattern.score_frequencies[total_score] += 1.0
        pattern.observed_statistics += observed
        pattern.observed_offset += observed_offset
        pattern.persons.append(str(person))

    return list(pattern_map.values()), pd.DataFrame(status_rows)


def _conditional_score_moments(
    log_kernel: np.ndarray,
    design: np.ndarray,
    *,
    compute_second: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Return score-wise log normalizers and sufficient-statistic moments.

    For each attainable raw score, the recurrence carries the polynomial
    coefficient ``Z``, its unnormalized first moment ``M = sum(w T)``, and,
    when requested, its second moment ``Q = sum(w T T')``.  The conditional
    mean is ``M / Z`` and the conditional covariance is
    ``Q / Z - (M / Z)(M / Z)'``.

    Per-unit log-kernel shifts and per-stage coefficient rescaling keep all
    recurrences in a finite numerical range.  A single scale applies to every
    score state at a stage, so neither moments nor conditional probabilities
    are altered.
    """
    log_kernel = np.asarray(log_kernel, dtype=float)
    design = np.asarray(design, dtype=float)
    if log_kernel.ndim != 2 or design.ndim != 3:
        raise ValueError("CMLE moment inputs must be a 2D kernel and 3D design.")
    if design.shape[:2] != log_kernel.shape:
        raise ValueError("CMLE kernel and design unit/category shapes differ.")
    if not np.all(np.isfinite(log_kernel)) or not np.all(np.isfinite(design)):
        raise FloatingPointError("CMLE moment inputs must be finite.")

    n_units, n_categories = log_kernel.shape
    n_parameters = design.shape[2]
    max_category = n_categories - 1
    max_score = n_units * max_category

    unit_shifts = np.max(log_kernel, axis=1)
    weights = np.exp(log_kernel - unit_shifts[:, None])
    coefficients = np.zeros(max_score + 1, dtype=float)
    first = np.zeros((max_score + 1, n_parameters), dtype=float)
    second = (
        np.zeros((max_score + 1, n_parameters, n_parameters), dtype=float)
        if compute_second
        else None
    )
    coefficients[0] = 1.0
    accumulated_log_scale = float(np.sum(unit_shifts))

    for unit in range(n_units):
        previous_length = unit * max_category + 1
        source_coefficients = coefficients[:previous_length]
        source_first = first[:previous_length, :]
        source_second = second[:previous_length, :, :] if second is not None else None
        next_length = (unit + 1) * max_category + 1
        next_coefficients = np.zeros(max_score + 1, dtype=float)
        next_first = np.zeros_like(first)
        next_second = np.zeros_like(second) if second is not None else None

        for category in range(n_categories):
            target = slice(category, category + previous_length)
            weight = float(weights[unit, category])
            statistic = design[unit, category, :]
            next_coefficients[target] += weight * source_coefficients
            next_first[target, :] += weight * (
                source_first + source_coefficients[:, None] * statistic[None, :]
            )
            if next_second is not None and source_second is not None:
                statistic_outer = np.outer(statistic, statistic)
                next_second[target, :, :] += weight * (
                    source_second
                    + source_first[:, :, None] * statistic[None, None, :]
                    + statistic[None, :, None] * source_first[:, None, :]
                    + source_coefficients[:, None, None]
                    * statistic_outer[None, :, :]
                )

        stage_scale = float(np.max(next_coefficients[:next_length]))
        if not np.isfinite(stage_scale) or stage_scale <= 0.0:
            raise FloatingPointError("CMLE score-polynomial scaling failed.")
        accumulated_log_scale += float(np.log(stage_scale))
        coefficients = next_coefficients / stage_scale
        first = next_first / stage_scale
        second = next_second / stage_scale if next_second is not None else None

    positive = coefficients > 0.0
    log_normalizers = np.full(max_score + 1, -np.inf, dtype=float)
    log_normalizers[positive] = (
        np.log(coefficients[positive]) + accumulated_log_scale
    )
    means = np.zeros_like(first)
    means[positive, :] = first[positive, :] / coefficients[positive, None]

    covariance: np.ndarray | None = None
    if second is not None:
        covariance = np.zeros_like(second)
        covariance[positive, :, :] = (
            second[positive, :, :] / coefficients[positive, None, None]
            - means[positive, :, None] * means[positive, None, :]
        )
        covariance = 0.5 * (covariance + np.swapaxes(covariance, 1, 2))
    return log_normalizers, means, covariance


def _cmle_objective(
    parameters: np.ndarray,
    design: CMLEDesign,
    *,
    compute_hessian: bool,
) -> tuple[float, np.ndarray, np.ndarray | None]:
    """Shared exact conditional objective, gradient, and information path."""
    parameters = np.asarray(parameters, dtype=float)
    if parameters.shape != (design.n_parameters,):
        raise ValueError(
            f"Expected {design.n_parameters} CMLE parameters, got {parameters.shape}."
        )
    log_likelihood = 0.0
    score_gradient = np.zeros(design.n_parameters, dtype=float)
    information = (
        np.zeros((design.n_parameters, design.n_parameters), dtype=float)
        if compute_hessian
        else None
    )
    for pattern in design.patterns:
        log_kernel = pattern.offset + np.einsum(
            "ukp,p->uk", pattern.design, parameters, optimize=True
        )
        log_normalizers, expected, covariance = _conditional_score_moments(
            log_kernel,
            pattern.design,
            compute_second=compute_hessian,
        )
        log_likelihood += float(
            pattern.observed_offset + pattern.observed_statistics @ parameters
        )
        score_gradient += pattern.observed_statistics
        used_scores = np.flatnonzero(pattern.score_frequencies > 0)
        frequencies = pattern.score_frequencies[used_scores]
        used_normalizers = log_normalizers[used_scores]
        if not np.all(np.isfinite(used_normalizers)):
            bad_score = int(used_scores[np.flatnonzero(~np.isfinite(used_normalizers))[0]])
            raise FloatingPointError(
                f"Conditional normalizer is non-finite for score {bad_score}."
            )
        log_likelihood -= float(frequencies @ used_normalizers)
        score_gradient -= np.einsum(
            "s,sp->p", frequencies, expected[used_scores, :], optimize=True
        )
        if information is not None and covariance is not None:
            information += np.einsum(
                "s,spq->pq", frequencies, covariance[used_scores, :, :], optimize=True
            )
    if information is not None:
        information = 0.5 * (information + information.T)
    return -float(log_likelihood), -score_gradient, information


def cmle_objective_value_grad(
    parameters: np.ndarray,
    design: CMLEDesign,
) -> tuple[float, np.ndarray]:
    """Exact negative conditional log-likelihood and analytical gradient."""
    value, gradient, _ = _cmle_objective(
        parameters, design, compute_hessian=False
    )
    return value, gradient


def cmle_objective_value_grad_hessian(
    parameters: np.ndarray,
    design: CMLEDesign,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Exact objective, gradient, and observed conditional information.

    The returned Hessian is the sum, over observed score strata, of the
    conditional covariance of the structural sufficient statistic.  It is
    therefore analytical for this linear Rasch-family kernel and does not use
    finite differences.
    """
    value, gradient, hessian = _cmle_objective(
        parameters, design, compute_hessian=True
    )
    if hessian is None:  # pragma: no cover - guaranteed by compute_hessian
        raise RuntimeError("CMLE analytical information was not computed.")
    return value, gradient, hessian


def _optimizer_value_grad(
    parameters: np.ndarray,
    design: CMLEDesign,
    state: _CMLEOptimizerEvaluationState,
) -> tuple[float, np.ndarray]:
    """Evaluate the exact objective while rejecting non-finite trial points.

    This adapter is used only by scipy's primary BFGS line search.  The public
    exact-objective functions retain their fail-closed behavior and continue
    to raise for a non-finite mathematical evaluation.  Returning positive
    infinity here lets the line search backtrack; it does not clip parameters,
    relax a readiness gate, or turn an invalid trial into likelihood evidence.
    """

    candidate = np.asarray(parameters, dtype=float)
    try:
        value, gradient = cmle_objective_value_grad(candidate, design)
        if not np.isfinite(value) or not np.all(np.isfinite(gradient)):
            raise FloatingPointError("CMLE objective or gradient is non-finite.")
    except FloatingPointError as exc:
        state.invalid_evaluations += 1
        state.last_invalid_reason = f"{type(exc).__name__}: {exc}"
        return np.inf, np.zeros_like(candidate)
    state.finite_evaluations += 1
    if value < state.best_value:
        state.best_value = float(value)
        state.best_parameters = candidate.copy()
    return float(value), np.asarray(gradient, dtype=float)


def _newton_polish(
    parameters: np.ndarray,
    design: CMLEDesign,
    *,
    stationarity_tolerance: float,
    max_steps: int,
) -> tuple[np.ndarray, float, np.ndarray, np.ndarray, list[dict[str, object]], str]:
    """Safeguarded exact-Newton refinement after the primary BFGS run.

    Polishing is attempted only while the analytical gradient exceeds the
    declared stationarity tolerance.  Each step requires full-rank positive-
    definite conditional information and an Armijo-accepted objective
    decrease.  Failure to meet any guard returns the best accepted point and a
    typed reason; it never converts a numerical failure into success by
    relaxing the tolerance.
    """
    parameters = np.asarray(parameters, dtype=float).copy()
    max_steps = max(0, int(max_steps))
    history: list[dict[str, object]] = []
    value, gradient, hessian = cmle_objective_value_grad_hessian(parameters, design)
    reason = (
        "polishing_disabled"
        if max_steps == 0
        else "maximum_polish_steps_reached"
    )

    for step_index in range(max_steps):
        gradient_sup = float(np.max(np.abs(gradient)))
        if gradient_sup <= float(stationarity_tolerance):
            reason = "stationarity_reached"
            break
        rank, nullity, _, rank_tolerance = _matrix_rank_evidence(hessian)
        minimum_eigenvalue = float(np.min(np.linalg.eigvalsh(hessian)))
        if nullity or rank != design.n_parameters:
            reason = "information_rank_deficient"
            break
        if minimum_eigenvalue <= rank_tolerance:
            reason = "information_not_positive_definite"
            break
        try:
            direction = np.linalg.solve(hessian, -gradient)
        except np.linalg.LinAlgError:
            reason = "newton_solve_failed"
            break
        directional_derivative = float(gradient @ direction)
        if not np.isfinite(directional_derivative) or directional_derivative >= 0.0:
            reason = "newton_direction_not_descent"
            break

        accepted = False
        acceptance = "none"
        step_size = 1.0
        candidate_value = np.nan
        for _ in range(24):
            candidate = parameters + step_size * direction
            candidate_value, candidate_gradient = cmle_objective_value_grad(
                candidate, design
            )
            if np.isfinite(candidate_value) and candidate_value <= (
                value + 1e-4 * step_size * directional_derivative
            ):
                accepted = True
                acceptance = "armijo"
                break
            objective_roundoff = (
                128.0
                * np.finfo(float).eps
                * max(1.0, abs(float(value)), abs(float(candidate_value)))
            )
            candidate_gradient_sup = float(np.max(np.abs(candidate_gradient)))
            if (
                np.isfinite(candidate_value)
                and candidate_value <= value + objective_roundoff
                and candidate_gradient_sup
                <= max(float(stationarity_tolerance), 0.5 * gradient_sup)
            ):
                accepted = True
                acceptance = "roundoff_gradient_reduction"
                break
            step_size *= 0.5
        history.append(
            {
                "Step": step_index + 1,
                "Accepted": accepted,
                "Acceptance": acceptance,
                "StepSize": float(step_size),
                "ObjectiveBefore": float(value),
                "ObjectiveAfter": float(candidate_value),
                "GradientSupNormBefore": gradient_sup,
            }
        )
        if not accepted:
            reason = "line_search_failed"
            break
        parameters = candidate
        value, gradient, hessian = cmle_objective_value_grad_hessian(
            parameters, design
        )
        reason = "maximum_polish_steps_reached"
    else:
        if float(np.max(np.abs(gradient))) <= float(stationarity_tolerance):
            reason = "stationarity_reached"

    return parameters, value, gradient, hessian, history, reason


def _matrix_rank_evidence(matrix: np.ndarray) -> tuple[int, int, float, float]:
    """Return rank, nullity, condition number, and rank tolerance."""
    matrix = np.asarray(matrix, dtype=float)
    if matrix.size == 0:
        return 0, 0, np.nan, np.nan
    singular = np.linalg.svd(matrix, compute_uv=False)
    largest = float(singular[0]) if singular.size else 0.0
    tolerance = max(1e-9, largest * max(matrix.shape) * 1e-8)
    rank = int(np.sum(singular > tolerance))
    nullity = int(matrix.shape[1] - rank)
    smallest = float(singular[-1]) if singular.size else 0.0
    condition = float(largest / smallest) if smallest > tolerance else np.inf
    return rank, nullity, condition, tolerance


def prepare_cmle_design(
    data: pd.DataFrame,
    *,
    person_col: str,
    facet_cols: Sequence[str],
    score_col: str,
    rating_min: int,
    rating_max: int,
    model: str = "RSM",
    step_facet: str | None = None,
    weight_col: str | None = None,
    response_unit_col: str | None = None,
    positive_facets: Iterable[str] | None = None,
    hard_anchors: pd.DataFrame | Sequence[Mapping[str, object]] | None = None,
    allow_duplicate_units: bool = False,
    rank_audit_max_parameters: int = 200,
    rank_audit_max_work: int = 50_000_000,
    rank_audit_max_bytes: int = DEFAULT_RANK_AUDIT_MAX_BYTES,
) -> CMLEDesign:
    """Prepare and audit an exact CMLE design from long-format rating data.

    ``rating_min`` and ``rating_max`` are mandatory so zero-count top or
    intermediate categories remain in the model support.  Row weights are
    audited but Phase 0 accepts only unit weights.
    """
    if not isinstance(data, pd.DataFrame) or data.empty:
        raise ValueError("CMLE requires a non-empty pandas DataFrame.")
    model = str(model or "RSM").strip().upper()
    if model not in SUPPORTED_CMLE_MODELS:
        raise ValueError(
            "Exact CMLE Phase 0 supports RSM and PCM only; estimated-slope "
            "GPCM does not retain the required raw-score sufficiency."
        )
    facet_cols = tuple(str(value) for value in facet_cols)
    if not facet_cols:
        raise ValueError("CMLE requires at least one non-person facet.")
    if len(set(facet_cols)) != len(facet_cols):
        raise ValueError("facet_cols contains duplicate names.")
    if model == "PCM":
        step_facet = str(step_facet) if step_facet is not None else facet_cols[0]
        if step_facet not in facet_cols:
            raise ValueError("PCM step_facet must be one of facet_cols.")
    else:
        step_facet = None

    rank_audit_max_parameters = _finite_integer(
        rank_audit_max_parameters, "rank_audit_max_parameters"
    )
    rank_audit_max_work = _finite_integer(
        rank_audit_max_work, "rank_audit_max_work"
    )
    rank_audit_max_bytes = _finite_integer(
        rank_audit_max_bytes, "rank_audit_max_bytes"
    )
    if min(
        rank_audit_max_parameters,
        rank_audit_max_work,
        rank_audit_max_bytes,
    ) <= 0:
        raise ValueError("CMLE rank-audit caps must be positive integers.")

    rating_min = _finite_integer(rating_min, "rating_min")
    rating_max = _finite_integer(rating_max, "rating_max")
    if rating_max <= rating_min:
        raise ValueError("rating_max must be larger than rating_min.")
    n_categories = rating_max - rating_min + 1

    required = [person_col, *facet_cols, score_col]
    if weight_col is not None:
        required.append(str(weight_col))
    if response_unit_col is not None:
        required.append(str(response_unit_col))
    if len(set(required)) != len(required):
        raise ValueError("Person, score, facet, weight, and response-unit roles must be distinct.")
    missing_columns = [column for column in required if column not in data.columns]
    if missing_columns:
        raise ValueError("CMLE input is missing columns: " + ", ".join(missing_columns))

    frame = data[required].copy()
    issues: list[dict[str, str]] = []
    score_numeric = pd.to_numeric(frame[score_col], errors="coerce")
    valid_required = frame[[person_col, *facet_cols]].notna().all(axis=1)
    if response_unit_col is not None:
        valid_required &= frame[response_unit_col].notna()
    valid_required &= score_numeric.notna() & np.isfinite(score_numeric)
    excluded_rows = int((~valid_required).sum())
    if excluded_rows:
        issues.append(
            _issue(
                "rows_excluded_missing_roles",
                "Caution",
                f"Excluded {excluded_rows} row(s) with missing/non-numeric required values.",
                "Represent planned missing responses as absent long-format rows and audit exclusions.",
            )
        )
    frame = frame.loc[valid_required].copy()
    score_numeric = score_numeric.loc[valid_required]
    if frame.empty:
        raise ValueError("No valid CMLE rows remain after required-value screening.")
    fractional = np.abs(score_numeric - np.round(score_numeric)) > 1e-10
    if fractional.any():
        raise ValueError("CMLE score values must be ordered integer category codes.")
    frame[score_col] = np.round(score_numeric).astype(int)
    outside = ~frame[score_col].between(rating_min, rating_max)
    if outside.any():
        values = sorted(frame.loc[outside, score_col].unique().tolist())
        raise ValueError(
            "Observed CMLE scores fall outside the declared rating support: "
            + ", ".join(map(str, values))
        )

    if weight_col is None:
        frame["__cmle_weight__"] = 1.0
    else:
        weights = pd.to_numeric(frame[weight_col], errors="coerce")
        frame["__cmle_weight__"] = weights
        if weights.isna().any() or (~np.isfinite(weights)).any() or (weights <= 0).any():
            issues.append(
                _issue(
                    "invalid_row_weights",
                    "Block",
                    "CMLE row weights contain missing, non-finite, or non-positive values.",
                    "Remove the weight role or supply valid unit weights for Phase 0.",
                )
            )
        elif not np.allclose(weights.to_numpy(dtype=float), 1.0, atol=1e-12, rtol=0):
            issues.append(
                _issue(
                    "nonunit_row_weights_unsupported",
                    "Block",
                    "Exact CMLE Phase 0 does not reinterpret observation weights as person-frequency weights.",
                    "Run an unweighted CMLE or retain the weighted analysis under its existing estimator.",
                )
            )

    frame[person_col] = frame[person_col].astype(str)
    for facet in facet_cols:
        frame[facet] = frame[facet].astype(str)
    if response_unit_col is not None:
        frame[response_unit_col] = frame[response_unit_col].astype(str)
    frame["__cmle_score__"] = frame[score_col] - rating_min
    frame.reset_index(drop=True, inplace=True)

    duplicate_key = [person_col, *facet_cols]
    if response_unit_col is not None:
        duplicate_key.append(response_unit_col)
    duplicate_rows = int(frame.duplicated(duplicate_key, keep=False).sum())
    if duplicate_rows and not allow_duplicate_units:
        issues.append(
            _issue(
                "duplicate_response_units",
                "Block",
                f"Found {duplicate_rows} row(s) in duplicated Person × response-unit cells.",
                "Add an explicit Occasion/Event response_unit_col or remove unintended duplicates.",
            )
        )
    elif duplicate_rows:
        issues.append(
            _issue(
                "duplicate_response_units_allowed",
                "Caution",
                f"Retained {duplicate_rows} duplicated response-unit row(s) by explicit override.",
                "Justify conditional independence of the repeated observations.",
            )
        )

    positive = {str(value) for value in (positive_facets or [])}
    unknown_positive = sorted(positive - set(facet_cols))
    if unknown_positive:
        raise ValueError("positive_facets contains unknown facets: " + ", ".join(unknown_positive))
    facet_signs = {facet: (1 if facet in positive else -1) for facet in facet_cols}
    facet_levels = {
        facet: sorted(frame[facet].astype(str).unique().tolist()) for facet in facet_cols
    }
    for facet, levels in facet_levels.items():
        if len(levels) < 2:
            issues.append(
                _issue(
                    "singleton_facet_fixed_zero",
                    "Caution",
                    f"Facet {facet!r} has one level and contributes no free contrast.",
                    "Remove the facet unless retaining its zero reference is substantively useful.",
                )
            )

    anchors = _normalize_hard_anchors(
        hard_anchors,
        facet_cols=facet_cols,
        facet_levels=facet_levels,
    )
    hard_anchor_values = {
        (str(row.Facet), str(row.Level)): float(row.Value)
        for row in anchors.itertuples(index=False)
    }

    observed_categories = set(frame[score_col].unique().tolist())
    unused_categories = sorted(set(range(rating_min, rating_max + 1)) - observed_categories)
    if unused_categories:
        issues.append(
            _issue(
                "declared_zero_count_categories_retained",
                "Caution",
                "Declared zero-count categories remain in the CMLE support: "
                + ", ".join(map(str, unused_categories)),
                "Review threshold information; do not truncate support to observed maxima.",
            )
        )

    parameter_names, facet_blocks, step_blocks, structural = _parameter_layout(
        facet_levels=facet_levels,
        facet_cols=facet_cols,
        facet_signs=facet_signs,
        model=model,
        step_facet=step_facet,
        n_categories=n_categories,
        hard_anchor_values=hard_anchor_values,
    )
    row_design, row_offset, unit_codes = _build_row_design(
        data=frame,
        facet_cols=facet_cols,
        facet_levels=facet_levels,
        facet_signs=facet_signs,
        facet_blocks=facet_blocks,
        step_blocks=step_blocks,
        model=model,
        step_facet=step_facet,
        n_categories=n_categories,
        n_parameters=len(parameter_names),
    )
    patterns, person_status = _build_patterns(
        data=frame,
        person_col=person_col,
        score_internal_col="__cmle_score__",
        row_design=row_design,
        row_offset=row_offset,
        unit_codes=unit_codes,
        n_categories=n_categories,
        n_parameters=len(parameter_names),
        response_unit_col=response_unit_col,
    )
    status_counts = person_status["ConditionalStatus"].value_counts().to_dict()
    informative_persons = int(status_counts.get("informative", 0))
    extreme_persons = int(status_counts.get("extreme_low", 0) + status_counts.get("extreme_high", 0))
    insufficient_persons = int(status_counts.get("insufficient_units", 0))
    if extreme_persons:
        issues.append(
            _issue(
                "extreme_persons_zero_structural_information",
                "Info",
                f"{extreme_persons} extreme-score Person(s) contribute zero conditional structural information.",
                "Retain their typed status; handle any later Person scoring separately.",
            )
        )
    if insufficient_persons:
        issues.append(
            _issue(
                "persons_with_fewer_than_two_units",
                "Caution",
                f"{insufficient_persons} Person(s) have fewer than two observed response units.",
                "They cannot contribute to exact conditional structural calibration.",
            )
        )
    if not parameter_names:
        issues.append(
            _issue(
                "no_free_structural_parameters",
                "Block",
                "The selected CMLE design has no free structural parameters.",
                "Include at least one varying facet or estimable threshold contrast.",
            )
        )
    if informative_persons == 0:
        issues.append(
            _issue(
                "no_informative_persons",
                "Block",
                "No Person has at least two units and a non-extreme raw score.",
                "Collect within-Person response variation before fitting CMLE.",
            )
        )

    audit: dict[str, object] = {
        "schema_version": CMLE_SCHEMA_VERSION,
        "eligible": False,
        "issues": issues,
        "persons_total": int(person_status.shape[0]),
        "persons_informative": informative_persons,
        "persons_extreme": extreme_persons,
        "persons_insufficient_units": insufficient_persons,
        "unique_missingness_design_patterns": int(len(patterns)),
        "parameters": int(len(parameter_names)),
        "conditional_rank": 0,
        "conditional_nullity": int(len(parameter_names)),
        "conditional_condition_number": np.nan,
        "conditional_min_eigenvalue": np.nan,
        "rank_tolerance": np.nan,
        "declared_rating_min": rating_min,
        "declared_rating_max": rating_max,
        "unused_declared_categories": unused_categories,
        "duplicate_rows": duplicate_rows,
        "excluded_rows": excluded_rows,
        "likelihood_scope": "structural_parameters_after_conditioning_out_person",
        "person_estimation_scope": "not_part_of_cmle_phase0",
        "hard_anchor_count": int(len(anchors)),
        "hard_anchored_facets": anchors["Facet"].drop_duplicates().tolist(),
        "hard_anchor_scale": "expanded_structural_estimate_before_facet_sign",
    }
    state_cells = int(
        sum(
            (pattern.design.shape[0] + 1)
            * (pattern.design.shape[0] * (n_categories - 1) + 1)
            for pattern in patterns
        )
    )
    rank_work = int(
        sum(
            pattern.design.shape[0]
            * n_categories
            * (pattern.design.shape[0] * (n_categories - 1) + 1)
            * max(1, len(parameter_names) ** 2)
            for pattern in patterns
        )
    )
    frame_bytes = int(frame.memory_usage(index=True, deep=True).sum())
    resident_design_bytes = int(
        row_design.nbytes
        + row_offset.nbytes
        + sum(
            pattern.design.nbytes
            + pattern.offset.nbytes
            + pattern.score_frequencies.nbytes
            + pattern.observed_statistics.nbytes
            for pattern in patterns
        )
    )
    workspace_bytes = int(
        max(
            (
                8
                * (pattern.design.shape[0] * (n_categories - 1) + 1)
                * (
                    4
                    + 4 * len(parameter_names)
                    + 6 * len(parameter_names) ** 2
                )
                + 8
                * pattern.design.shape[0]
                * n_categories
                * (len(parameter_names) + 1)
            )
            for pattern in patterns
        )
        if patterns
        else 0
    )
    peak_bytes = int(frame_bytes + resident_design_bytes + workspace_bytes)
    audit["conditional_dp_state_cells"] = state_cells
    audit["rank_audit_work_proxy"] = rank_work
    audit["rank_audit_method"] = "exact_conditional_second_moment_dp"
    audit["input_frame_bytes"] = frame_bytes
    audit["conditional_design_resident_bytes"] = resident_design_bytes
    audit["rank_audit_workspace_bytes_proxy"] = workspace_bytes
    audit["rank_audit_peak_bytes_proxy"] = peak_bytes
    audit["rank_audit_max_parameters"] = int(rank_audit_max_parameters)
    audit["rank_audit_max_work"] = int(rank_audit_max_work)
    audit["rank_audit_max_bytes"] = int(rank_audit_max_bytes)
    if len(parameter_names) > int(rank_audit_max_parameters):
        issues.append(
            _issue(
                "phase0_rank_audit_parameter_cap",
                "Block",
                f"CMLE Phase 0 rank audit has {len(parameter_names)} parameters; cap is {int(rank_audit_max_parameters)}.",
                "Reduce facet cardinality or validate a reviewed higher-cap local research run.",
            )
        )
    if rank_work > int(rank_audit_max_work):
        issues.append(
            _issue(
                "phase0_rank_audit_work_cap",
                "Block",
                f"CMLE Phase 0 rank-work proxy is {rank_work:,}; cap is {int(rank_audit_max_work):,}.",
                "Reduce response units/patterns or use a reviewed optimized rank-audit implementation.",
            )
        )
    if peak_bytes > int(rank_audit_max_bytes):
        issues.append(
            _issue(
                "phase0_rank_audit_memory_cap",
                "Block",
                "CMLE Phase 0 analytical-information peak-memory proxy is "
                f"{peak_bytes / (1024 ** 2):,.1f} MiB; cap is "
                f"{int(rank_audit_max_bytes) / (1024 ** 2):,.1f} MiB.",
                "Reduce response units, categories, facet cardinality, or "
                "validate a reviewed higher-memory local research run.",
            )
        )
    prepared = CMLEDesign(
        data=frame,
        person_col=person_col,
        facet_cols=facet_cols,
        score_col=score_col,
        model=model,
        step_facet=step_facet,
        rating_min=rating_min,
        rating_max=rating_max,
        n_categories=n_categories,
        facet_levels=facet_levels,
        facet_signs=facet_signs,
        parameter_names=parameter_names,
        structural_parameters=structural,
        row_design=row_design,
        row_offset=row_offset,
        unit_codes=unit_codes,
        patterns=patterns,
        person_status=person_status,
        audit=audit,
        hard_anchors=anchors,
        response_unit_col=response_unit_col,
    )

    has_block_before_rank = any(row["Severity"] == "Block" for row in issues)
    if parameter_names and informative_persons and not has_block_before_rank:
        try:
            _, _, hessian_zero = cmle_objective_value_grad_hessian(
                np.zeros(len(parameter_names)), prepared
            )
            rank, nullity, condition, tolerance = _matrix_rank_evidence(hessian_zero)
            audit["conditional_rank"] = rank
            audit["conditional_nullity"] = nullity
            audit["conditional_condition_number"] = condition
            audit["rank_tolerance"] = tolerance
            minimum_eigenvalue = float(np.min(np.linalg.eigvalsh(hessian_zero)))
            audit["conditional_min_eigenvalue"] = minimum_eigenvalue
            if nullity:
                issues.append(
                    _issue(
                        "conditional_information_rank_deficient",
                        "Block",
                        f"Conditional information rank is {rank}/{len(parameter_names)} (nullity {nullity}).",
                        "Add within-Person facet contrasts or remove conditionally unidentified effects.",
                    )
                )
            elif minimum_eigenvalue < -tolerance:
                issues.append(
                    _issue(
                        "conditional_information_not_positive_semidefinite",
                        "Block",
                        f"Conditional information has minimum eigenvalue {minimum_eigenvalue:.3g}.",
                        "Resolve the numerical/design inconsistency before optimization.",
                    )
                )
            elif np.isfinite(condition) and condition > 1e10:
                issues.append(
                    _issue(
                        "conditional_information_ill_conditioned",
                        "Caution",
                        f"Conditional information condition number is {condition:.3g}.",
                        "Inspect weak facet bridges and threshold/category exposure.",
                    )
                )
        except Exception as exc:  # pragma: no cover - defensive fail-closed guard
            issues.append(
                _issue(
                    "conditional_rank_audit_failed",
                    "Block",
                    f"Conditional rank audit failed: {exc}",
                    "Resolve the design or numerical failure before optimization.",
                )
            )
    audit["eligible"] = not any(row["Severity"] == "Block" for row in issues)
    audit["issues_table"] = pd.DataFrame(
        issues, columns=["Code", "Severity", "Message", "Action"]
    )
    return prepared


def audit_cmle_eligibility(*args, **kwargs) -> dict[str, object]:
    """Return a fail-closed CMLE audit without starting optimization."""
    try:
        design = prepare_cmle_design(*args, **kwargs)
    except Exception as exc:
        issues = [
            _issue(
                "cmle_input_invalid",
                "Block",
                str(exc),
                "Correct the input/model specification before fitting CMLE.",
            )
        ]
        return {
            "schema_version": CMLE_SCHEMA_VERSION,
            "eligible": False,
            "issues": issues,
            "issues_table": pd.DataFrame(issues),
        }
    return design.audit


def _expanded_parameter_table(
    design: CMLEDesign,
    parameters: np.ndarray,
    covariance: np.ndarray | None,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for structural in design.structural_parameters:
        estimate = float(structural.offset + structural.jacobian @ parameters)
        if structural.anchored:
            standard_error = 0.0
        elif covariance is None:
            standard_error = np.nan
        else:
            variance = float(structural.jacobian @ covariance @ structural.jacobian)
            standard_error = np.sqrt(max(variance, 0.0)) if np.isfinite(variance) else np.nan
        rows.append(
            {
                "ParameterType": structural.parameter_type,
                "Facet": structural.facet,
                "Level": structural.level,
                "Step": structural.step,
                "Estimate": estimate,
                "SE": standard_error,
                "Anchored": bool(structural.anchored),
                "Constraint": structural.constraint,
            }
        )
    return pd.DataFrame(rows)


def _category_surfaces(design: CMLEDesign, parameters: np.ndarray) -> pd.DataFrame:
    seen: set[tuple[int, ...]] = set()
    rows: list[dict[str, object]] = []
    for row_index, codes in enumerate(design.unit_codes):
        if codes in seen:
            continue
        seen.add(codes)
        labels = {
            facet: design.facet_levels[facet][code]
            for facet, code in zip(design.facet_cols, codes)
        }
        unit_id = "|".join(f"{facet}={labels[facet]}" for facet in design.facet_cols)
        for category in range(design.n_categories):
            kernel = float(
                design.row_offset[row_index, category]
                + design.row_design[row_index, category, :] @ parameters
            )
            rows.append(
                {
                    "VirtualUnit": unit_id,
                    **labels,
                    "Category": category + design.rating_min,
                    "InternalCategory": category,
                    "LogKernelWithoutPerson": kernel,
                    "CumulativeDifficulty": -kernel,
                }
            )
    return pd.DataFrame(rows)


def fit_cmle(
    data: pd.DataFrame,
    *,
    person_col: str,
    facet_cols: Sequence[str],
    score_col: str,
    rating_min: int,
    rating_max: int,
    model: str = "RSM",
    step_facet: str | None = None,
    weight_col: str | None = None,
    response_unit_col: str | None = None,
    positive_facets: Iterable[str] | None = None,
    hard_anchors: pd.DataFrame | Sequence[Mapping[str, object]] | None = None,
    allow_duplicate_units: bool = False,
    rank_audit_max_parameters: int = 200,
    rank_audit_max_work: int = 50_000_000,
    rank_audit_max_bytes: int = DEFAULT_RANK_AUDIT_MAX_BYTES,
    maxiter: int = 500,
    gtol: float = 1e-7,
    newton_polish_maxiter: int = 8,
    finite_mle_gate: bool = True,
) -> dict[str, object]:
    """Fit exact CMLE after a strict eligibility and conditional-rank audit."""
    if not isinstance(finite_mle_gate, (bool, np.bool_)):
        raise ValueError("finite_mle_gate must be boolean.")
    finite_mle_gate = bool(finite_mle_gate)
    newton_polish_maxiter = _finite_integer(
        newton_polish_maxiter, "newton_polish_maxiter"
    )
    if newton_polish_maxiter < 0:
        raise ValueError("newton_polish_maxiter must be non-negative.")
    design = prepare_cmle_design(
        data,
        person_col=person_col,
        facet_cols=facet_cols,
        score_col=score_col,
        rating_min=rating_min,
        rating_max=rating_max,
        model=model,
        step_facet=step_facet,
        weight_col=weight_col,
        response_unit_col=response_unit_col,
        positive_facets=positive_facets,
        hard_anchors=hard_anchors,
        allow_duplicate_units=allow_duplicate_units,
        rank_audit_max_parameters=rank_audit_max_parameters,
        rank_audit_max_work=rank_audit_max_work,
        rank_audit_max_bytes=rank_audit_max_bytes,
    )
    if not bool(design.audit.get("eligible", False)):
        blocking = [
            row["Message"]
            for row in design.audit.get("issues", [])
            if row.get("Severity") == "Block"
        ]
        raise CMLEEligibilityError("CMLE eligibility failed: " + " ".join(blocking))

    finite_mle_audit: dict[str, object] | None = None
    finite_mle_status = "not_evaluated"
    finite_mle_reason = "gate_explicitly_disabled"
    finite_mle_boundary: object = pd.NA
    finite_mle_qualified = False
    finite_mle_tolerance_grid = ""
    finite_mle_theoretical_configurations: object = np.nan
    finite_mle_oracle_state_cells: object = np.nan
    finite_mle_generated_constraints: object = np.nan
    finite_mle_oracle_calls: object = np.nan
    if finite_mle_gate:
        # Local import avoids a module cycle: the support-oracle module uses
        # CMLEDesign and the exact objective defined in this module.
        from mfrm_app.cmle_existence import audit_cmle_finite_mle_oracle

        finite_mle_audit = audit_cmle_finite_mle_oracle(design)
        finite_mle_summary = finite_mle_audit["summary"].iloc[0]
        finite_mle_status = str(finite_mle_summary["Status"])
        finite_mle_reason = str(finite_mle_summary["Reason"])
        finite_mle_boundary = finite_mle_summary.get(
            "BoundaryDetected", pd.NA
        )
        finite_mle_qualified = bool(
            finite_mle_summary.get("ExistenceQualified", False)
        )
        finite_mle_tolerance_grid = str(
            finite_mle_summary.get("ToleranceGrid", "")
        )
        finite_mle_theoretical_configurations = finite_mle_summary.get(
            "TheoreticalConfigurations", np.nan
        )
        finite_mle_oracle_state_cells = finite_mle_summary.get(
            "OracleStateCells", np.nan
        )
        finite_mle_generated_constraints = finite_mle_summary.get(
            "GeneratedConstraintsMax", np.nan
        )
        finite_mle_oracle_calls = finite_mle_summary.get(
            "OracleCallsTotal", np.nan
        )

    start = np.zeros(design.n_parameters, dtype=float)
    objective_trace: list[float] = []
    optimizer_state = _CMLEOptimizerEvaluationState()

    def callback(parameters: np.ndarray) -> None:
        try:
            value, _ = cmle_objective_value_grad(parameters, design)
        except FloatingPointError:  # pragma: no cover - scipy accepts finite points
            return
        if np.isfinite(value):
            objective_trace.append(float(value))

    result = minimize(
        _optimizer_value_grad,
        start,
        args=(design, optimizer_state),
        method="BFGS",
        jac=True,
        callback=callback,
        options={"maxiter": int(maxiter), "gtol": float(gtol)},
    )
    optimizer_parameters = np.asarray(result.x, dtype=float)
    optimizer_best_finite_fallback_used = False
    try:
        optimizer_nll, optimizer_gradient = cmle_objective_value_grad(
            optimizer_parameters, design
        )
        if not np.isfinite(optimizer_nll) or not np.all(
            np.isfinite(optimizer_gradient)
        ):
            raise FloatingPointError("Primary optimizer ended at a non-finite point.")
    except FloatingPointError:
        if optimizer_state.best_parameters is None:
            raise FloatingPointError(
                "CMLE primary optimizer did not retain a finite objective point."
            )
        optimizer_parameters = optimizer_state.best_parameters.copy()
        optimizer_nll, optimizer_gradient = cmle_objective_value_grad(
            optimizer_parameters, design
        )
        optimizer_best_finite_fallback_used = True
    optimizer_gradient_sup_norm = float(np.max(np.abs(optimizer_gradient)))
    stationarity_tolerance = max(1e-5, 10.0 * float(gtol))
    (
        parameters,
        final_nll,
        final_gradient,
        hessian,
        newton_polish_trace,
        newton_polish_reason,
    ) = _newton_polish(
        optimizer_parameters,
        design,
        stationarity_tolerance=stationarity_tolerance,
        max_steps=int(newton_polish_maxiter),
    )
    objective_trace.extend(
        float(row["ObjectiveAfter"])
        for row in newton_polish_trace
        if bool(row["Accepted"])
    )
    gradient_norm = float(np.linalg.norm(final_gradient))
    gradient_sup_norm = float(np.max(np.abs(final_gradient)))
    rank, nullity, condition, rank_tolerance = _matrix_rank_evidence(hessian)
    information_eigenvalues = np.linalg.eigvalsh(hessian)
    information_min_eigenvalue = float(np.min(information_eigenvalues))
    information_positive_definite = bool(
        information_min_eigenvalue > rank_tolerance
    )
    covariance: np.ndarray | None
    try:
        if nullity == 0 and information_positive_definite:
            covariance_raw = np.linalg.inv(hessian)
            covariance = 0.5 * (covariance_raw + covariance_raw.T)
        else:
            covariance = None
    except np.linalg.LinAlgError:
        covariance = None
    optimizer_success = bool(
        result.success
        and np.isfinite(optimizer_nll)
        and not optimizer_best_finite_fallback_used
    )
    converged = bool(
        np.isfinite(final_nll)
        and (optimizer_success or gradient_sup_norm <= stationarity_tolerance)
    )
    finite_mle_readiness_pass = bool(
        finite_mle_qualified if finite_mle_gate else True
    )
    inference_ready = bool(
        finite_mle_readiness_pass
        and converged
        and nullity == 0
        and information_positive_definite
        and covariance is not None
        and np.all(np.isfinite(covariance))
        and gradient_sup_norm <= stationarity_tolerance
    )
    readiness_reasons: list[str] = []
    if finite_mle_gate and not finite_mle_qualified:
        if finite_mle_status == "boundary_no_finite_cmle":
            readiness_reasons.append("finite_mle_boundary_no_finite_cmle")
        else:
            readiness_reasons.append(f"finite_mle_{finite_mle_status}")
    if not converged:
        readiness_reasons.append("stationarity_not_reached")
    if nullity:
        readiness_reasons.append("fitted_information_rank_deficient")
    if not information_positive_definite:
        readiness_reasons.append("fitted_information_not_positive_definite")
    if covariance is None or not np.all(np.isfinite(covariance)):
        readiness_reasons.append("covariance_unavailable")
    if gradient_sup_norm > stationarity_tolerance:
        readiness_reasons.append("gradient_above_readiness_tolerance")

    expanded = _expanded_parameter_table(design, parameters, covariance)
    facets = expanded.loc[expanded["ParameterType"] == "Facet"].copy()
    steps = expanded.loc[expanded["ParameterType"] == "Step"].copy()
    free_parameters = pd.DataFrame(
        {
            "Parameter": design.parameter_names,
            "Estimate": parameters,
            "Gradient": final_gradient,
            "SE": (
                np.sqrt(np.clip(np.diag(covariance), 0.0, np.inf))
                if covariance is not None
                else np.full(design.n_parameters, np.nan)
            ),
        }
    )
    conditional_loglik = -float(final_nll)
    summary = pd.DataFrame(
        [
            {
                "Model": design.model,
                "Method": "CMLE",
                "EstimatorLabel": "Native exact conditional maximum likelihood",
                "PersonsTotal": int(design.audit["persons_total"]),
                "PersonsInformative": int(design.audit["persons_informative"]),
                "PersonsExtreme": int(design.audit["persons_extreme"]),
                "Patterns": len(design.patterns),
                "KParams": design.n_parameters,
                "HardAnchors": int(len(design.hard_anchors)),
                "ConditionalLogLik": conditional_loglik,
                "ConditionalDeviance": -2.0 * conditional_loglik,
                "ConditionalAIC": 2.0 * design.n_parameters - 2.0 * conditional_loglik,
                "ConditionalAICScope": "compare only matched conditional likelihoods",
                "Converged": converged,
                "OptimizerSuccess": optimizer_success,
                "OptimizerGradientSupNorm": optimizer_gradient_sup_norm,
                "OptimizerFiniteEvaluations": int(
                    optimizer_state.finite_evaluations
                ),
                "OptimizerInvalidEvaluations": int(
                    optimizer_state.invalid_evaluations
                ),
                "OptimizerBestFiniteFallbackUsed": bool(
                    optimizer_best_finite_fallback_used
                ),
                "OptimizerLastInvalidReason": optimizer_state.last_invalid_reason,
                "NewtonPolishSteps": len(newton_polish_trace),
                "NewtonPolishAcceptedSteps": int(
                    sum(bool(row["Accepted"]) for row in newton_polish_trace)
                ),
                "NewtonPolishReason": newton_polish_reason,
                "FiniteMLEGateEnabled": finite_mle_gate,
                "FiniteMLEStatus": finite_mle_status,
                "FiniteMLEReason": finite_mle_reason,
                "FiniteMLEBoundaryDetected": finite_mle_boundary,
                "FiniteMLEExistenceQualified": finite_mle_qualified,
                "FiniteMLEToleranceGrid": finite_mle_tolerance_grid,
                "FiniteMLETheoreticalConfigurations": (
                    finite_mle_theoretical_configurations
                ),
                "FiniteMLEOracleStateCells": finite_mle_oracle_state_cells,
                "FiniteMLEGeneratedConstraintsMax": (
                    finite_mle_generated_constraints
                ),
                "FiniteMLEOracleCallsTotal": finite_mle_oracle_calls,
                "InferenceReady": inference_ready,
                "Iterations": int(getattr(result, "nit", 0)),
                "FunctionEvaluations": int(getattr(result, "nfev", 0)),
                "GradientNorm": gradient_norm,
                "GradientSupNorm": gradient_sup_norm,
                "StationarityTolerance": stationarity_tolerance,
                "InformationRank": rank,
                "InformationNullity": nullity,
                "InformationMethod": "exact conditional sufficient-statistic covariance",
                "InformationConditionNumber": condition,
                "InformationMinEigenvalue": information_min_eigenvalue,
                "InformationPositiveDefinite": information_positive_definite,
                "OptimizerMessage": str(result.message),
                "ReadinessReasons": ";".join(readiness_reasons),
            }
        ]
    )
    return {
        "summary": summary,
        "coefficients": free_parameters,
        "facets": {
            "person": pd.DataFrame(
                columns=["Person", "Estimate", "SE", "ScoringMethod"]
            ),
            "others": facets,
        },
        "steps": steps,
        "surfaces": _category_surfaces(design, parameters),
        "person_status": design.person_status.copy(),
        "covariance": covariance,
        "information": hessian,
        "audit": design.audit,
        "config": {
            "schema_version": CMLE_SCHEMA_VERSION,
            "model": design.model,
            "method": "CMLE",
            "person_col": person_col,
            "facet_names": list(design.facet_cols),
            "score_col": score_col,
            "step_facet": design.step_facet,
            "rating_min": design.rating_min,
            "rating_max": design.rating_max,
            "facet_signs": dict(design.facet_signs),
            "hard_anchors": design.hard_anchors.to_dict(orient="records"),
            "hard_anchor_scale": "expanded_structural_estimate_before_facet_sign",
            "response_unit_col": response_unit_col,
            "weight_col": weight_col,
            "person_scoring": "not_part_of_cmle_phase0",
            "information_method": "exact_conditional_second_moment_dp",
            "likelihood_comparability": (
                "conditional likelihood; do not compare generic AIC/BIC with JMLE or MML"
            ),
            "rank_tolerance": rank_tolerance,
            "newton_polish_maxiter": int(newton_polish_maxiter),
            "finite_mle_gate": finite_mle_gate,
            "finite_mle_status": finite_mle_status,
            "finite_mle_reason": finite_mle_reason,
            "rank_audit_max_parameters": int(rank_audit_max_parameters),
            "rank_audit_max_work": int(rank_audit_max_work),
            "rank_audit_max_bytes": int(rank_audit_max_bytes),
            "rank_audit_work_proxy": int(
                design.audit["rank_audit_work_proxy"]
            ),
            "rank_audit_peak_bytes_proxy": int(
                design.audit["rank_audit_peak_bytes_proxy"]
            ),
        },
        "optimizer": result,
        "objective_trace": objective_trace,
        "newton_polish_trace": pd.DataFrame(newton_polish_trace),
        "finite_mle_audit": finite_mle_audit,
        "design": design,
    }


__all__ = [
    "CMLEDesign",
    "CMLEEligibilityError",
    "CMLEPattern",
    "CMLE_SCHEMA_VERSION",
    "DEFAULT_RANK_AUDIT_MAX_BYTES",
    "SUPPORTED_CMLE_MODELS",
    "audit_cmle_eligibility",
    "cmle_objective_value_grad",
    "cmle_objective_value_grad_hessian",
    "fit_cmle",
    "prepare_cmle_design",
]
