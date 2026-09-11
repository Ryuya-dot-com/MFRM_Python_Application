"""Finite-CMLE existence diagnostics for the repository-only exact core.

For a full-rank conditional exponential family, a finite maximum exists only
when the observed sufficient statistic is in the relative interior of its
conditional convex support.  This module enumerates registered small/medium
score-stratum supports and searches their aggregate supporting cone by linear
programming.  It is a research audit, not yet a public estimator route.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import linprog

from mfrm_app.cmle import (
    CMLEDesign,
    cmle_objective_value_grad,
)


DEFAULT_EXISTENCE_TOLERANCES = (1e-10, 1e-9, 1e-8)
DEFAULT_DIRECTIONAL_RADII = (0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
DEFAULT_MAX_CONFIGURATIONS = 2_000_000
DEFAULT_MAX_CONSTRAINTS = 500_000
DEFAULT_MAX_BYTES = 512 * 1024 * 1024
DEFAULT_MAX_CUTTING_PLANE_ROUNDS = 10_000
DEFAULT_MAX_ORACLE_CUTS = 500_000
DEFAULT_MAX_ORACLE_STATE_CELLS = 50_000_000


class CMLEExistenceAuditError(RuntimeError):
    """Raised only for malformed internal existence-audit state."""


@dataclass(frozen=True)
class _StratumSupport:
    pattern_index: int
    raw_score: int
    persons: int
    configurations: int
    support: np.ndarray
    observed_sum: np.ndarray


@dataclass(frozen=True)
class _OracleStratum:
    pattern_index: int
    raw_score: int
    persons: int
    configurations: int
    state_cells: int
    design: np.ndarray
    observed_sum: np.ndarray


def _validate_tolerances(values: Iterable[float]) -> tuple[float, ...]:
    tolerances = tuple(float(value) for value in values)
    if (
        not tolerances
        or any(not np.isfinite(value) or value <= 0.0 for value in tolerances)
        or len(set(tolerances)) != len(tolerances)
    ):
        raise ValueError("existence tolerances must be unique positive finite values.")
    return tolerances


def _positive_integer(value: object, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a positive integer.")
    numeric = int(value)
    if float(value) != numeric or numeric < 1:
        raise ValueError(f"{name} must be a positive integer.")
    return numeric


@lru_cache(maxsize=None)
def _score_configuration_count(
    n_units: int, n_categories: int, raw_score: int
) -> int:
    """Count bounded category vectors without materializing them."""

    n_units = int(n_units)
    n_categories = int(n_categories)
    raw_score = int(raw_score)
    maximum = n_categories - 1
    if n_units < 1 or n_categories < 2 or not 0 <= raw_score <= n_units * maximum:
        return 0
    counts = np.zeros(raw_score + 1, dtype=object)
    counts[0] = 1
    for _ in range(n_units):
        updated = np.zeros(raw_score + 1, dtype=object)
        for score in range(raw_score + 1):
            if counts[score] == 0:
                continue
            for category in range(min(maximum, raw_score - score) + 1):
                updated[score + category] += counts[score]
        counts = updated
    return int(counts[raw_score])


@lru_cache(maxsize=None)
def _score_configurations(
    n_units: int, n_categories: int, raw_score: int
) -> np.ndarray:
    """Return every bounded category vector with the requested raw score."""

    n_units = int(n_units)
    n_categories = int(n_categories)
    raw_score = int(raw_score)
    maximum = n_categories - 1
    if n_units < 1 or n_categories < 2 or not 0 <= raw_score <= n_units * maximum:
        return np.empty((0, max(n_units, 0)), dtype=np.int16)
    rows: list[tuple[int, ...]] = []

    def visit(prefix: tuple[int, ...], units_left: int, score_left: int) -> None:
        if units_left == 0:
            if score_left == 0:
                rows.append(prefix)
            return
        lower = max(0, score_left - (units_left - 1) * maximum)
        upper = min(maximum, score_left)
        for category in range(lower, upper + 1):
            visit(prefix + (category,), units_left - 1, score_left - category)

    visit((), n_units, raw_score)
    return np.asarray(rows, dtype=np.int16).reshape(-1, n_units)


def _ordered_person_rows(design: CMLEDesign, person_frame: pd.DataFrame) -> np.ndarray:
    indices = person_frame.index.to_numpy(dtype=int)

    def order_key(row_index: int) -> tuple[object, ...]:
        key: tuple[object, ...] = tuple(design.unit_codes[row_index])
        if design.response_unit_col is not None:
            key += (str(design.data.at[row_index, design.response_unit_col]),)
        return key

    return np.asarray(sorted(indices.tolist(), key=order_key), dtype=int)


def _observed_by_pattern_score(
    design: CMLEDesign,
) -> dict[tuple[int, int], np.ndarray]:
    """Reconstruct score-specific observed sufficient-statistic sums."""

    groups = {
        str(person): frame
        for person, frame in design.data.groupby(
            design.person_col, sort=True, observed=True
        )
    }
    observed: dict[tuple[int, int], np.ndarray] = {}
    for pattern_index, pattern in enumerate(design.patterns):
        frequency_check = np.zeros_like(pattern.score_frequencies, dtype=int)
        for person in pattern.persons:
            if str(person) not in groups:
                raise CMLEExistenceAuditError(
                    f"Pattern Person {person!r} is absent from the prepared design."
                )
            ordered = _ordered_person_rows(design, groups[str(person)])
            scores = design.data.loc[ordered, "__cmle_score__"].to_numpy(dtype=int)
            signature = tuple(design.unit_codes[index] for index in ordered)
            if signature != pattern.signature:
                raise CMLEExistenceAuditError("Pattern signature reconstruction failed.")
            raw_score = int(np.sum(scores))
            statistic = np.sum(
                design.row_design[ordered, scores, :], axis=0
            )
            key = (pattern_index, raw_score)
            observed.setdefault(
                key, np.zeros(design.n_parameters, dtype=float)
            )
            observed[key] += statistic
            frequency_check[raw_score] += 1
        if not np.array_equal(
            frequency_check.astype(float), pattern.score_frequencies
        ):
            raise CMLEExistenceAuditError(
                "Score-specific observed frequency reconstruction failed."
            )
    return observed


def _enumerate_strata(
    design: CMLEDesign,
    *,
    max_configurations: int,
    max_constraints: int,
    max_bytes: int,
) -> tuple[list[_StratumSupport], pd.DataFrame, str | None]:
    observed = _observed_by_pattern_score(design)
    strata: list[_StratumSupport] = []
    rows: list[dict[str, object]] = []
    total_configurations = 0
    total_support = 0
    bytes_proxy = 0
    planned: list[tuple[int, object, int, int, int]] = []
    for pattern_index, pattern in enumerate(design.patterns):
        n_units = pattern.design.shape[0]
        used_scores = np.flatnonzero(pattern.score_frequencies > 0)
        for raw_score in used_scores:
            count = _score_configuration_count(
                n_units, design.n_categories, int(raw_score)
            )
            total_configurations += count
            planned.append(
                (pattern_index, pattern, n_units, int(raw_score), count)
            )
            rows.append(
                {
                    "PatternIndex": pattern_index,
                    "RawScore": int(raw_score),
                    "Persons": int(pattern.score_frequencies[int(raw_score)]),
                    "ResponseUnits": n_units,
                    "Configurations": count,
                    "UniqueSupportPoints": np.nan,
                    "ObservedStatistic": "",
                    "EnumerationBytesProxy": 0,
                }
            )
    if total_configurations > max_configurations:
        return [], pd.DataFrame(rows), "configuration_cap_exceeded"

    rows = []
    for pattern_index, pattern, n_units, raw_score, planned_count in planned:
        configurations = _score_configurations(
            n_units, design.n_categories, int(raw_score)
        )
        count = int(len(configurations))
        if count != planned_count:
            raise CMLEExistenceAuditError(
                "Score-configuration count and enumeration disagree."
            )
        statistic = np.zeros(
            (count, design.n_parameters), dtype=float
        )
        for unit in range(n_units):
            statistic += pattern.design[
                unit, configurations[:, unit], :
            ]
        support = np.unique(statistic, axis=0)
        total_support += len(support)
        if total_support > max_constraints:
            return strata, pd.DataFrame(rows), "constraint_cap_exceeded"
        bytes_proxy += int(
            configurations.nbytes + statistic.nbytes + support.nbytes
        )
        if bytes_proxy > max_bytes:
            return strata, pd.DataFrame(rows), "memory_cap_exceeded"
        key = (pattern_index, int(raw_score))
        observed_sum = observed.get(key)
        persons = int(pattern.score_frequencies[int(raw_score)])
        if observed_sum is None:
            raise CMLEExistenceAuditError(
                "Observed score stratum is missing after reconstruction."
            )
        strata.append(
            _StratumSupport(
                pattern_index=pattern_index,
                raw_score=int(raw_score),
                persons=persons,
                configurations=count,
                support=support,
                observed_sum=observed_sum,
            )
        )
        rows.append(
            {
                "PatternIndex": pattern_index,
                "RawScore": int(raw_score),
                "Persons": persons,
                "ResponseUnits": n_units,
                "Configurations": count,
                "UniqueSupportPoints": int(len(support)),
                "ObservedStatistic": json_vector(observed_sum),
                "EnumerationBytesProxy": int(bytes_proxy),
            }
        )
    return strata, pd.DataFrame(rows), None


def _oracle_strata(
    design: CMLEDesign,
    *,
    max_state_cells: int,
    max_bytes: int,
) -> tuple[list[_OracleStratum], pd.DataFrame, str | None]:
    """Build fixed-score oracle strata without materializing response vectors."""

    observed = _observed_by_pattern_score(design)
    strata: list[_OracleStratum] = []
    rows: list[dict[str, object]] = []
    total_state_cells = 0
    bytes_proxy = 0
    for pattern_index, pattern in enumerate(design.patterns):
        n_units = int(pattern.design.shape[0])
        used_scores = np.flatnonzero(pattern.score_frequencies > 0)
        for raw_score_value in used_scores:
            raw_score = int(raw_score_value)
            configurations = _score_configuration_count(
                n_units, design.n_categories, raw_score
            )
            state_cells = int((n_units + 1) * (raw_score + 1))
            total_state_cells += state_cells
            # Value, predecessor-score, and predecessor-category workspaces.
            bytes_proxy += int(state_cells * (8 + 4 + 2))
            key = (pattern_index, raw_score)
            observed_sum = observed.get(key)
            if observed_sum is None:
                raise CMLEExistenceAuditError(
                    "Observed score stratum is missing after reconstruction."
                )
            persons = int(pattern.score_frequencies[raw_score])
            rows.append(
                {
                    "PatternIndex": pattern_index,
                    "RawScore": raw_score,
                    "Persons": persons,
                    "ResponseUnits": n_units,
                    "TheoreticalConfigurations": configurations,
                    "OracleStateCells": state_cells,
                    "ObservedStatistic": json_vector(observed_sum),
                    "OracleBytesProxy": bytes_proxy,
                }
            )
            if total_state_cells > max_state_cells:
                return [], pd.DataFrame(rows), "oracle_state_cell_cap_exceeded"
            if bytes_proxy > max_bytes:
                return [], pd.DataFrame(rows), "oracle_memory_cap_exceeded"
            strata.append(
                _OracleStratum(
                    pattern_index=pattern_index,
                    raw_score=raw_score,
                    persons=persons,
                    configurations=configurations,
                    state_cells=state_cells,
                    design=pattern.design,
                    observed_sum=observed_sum,
                )
            )
    return strata, pd.DataFrame(rows), None


def _fixed_score_support_maximum(
    design: np.ndarray,
    raw_score: int,
    direction: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Maximize a linear sufficient statistic over one fixed-score support."""

    n_units, n_categories, n_parameters = design.shape
    raw_score = int(raw_score)
    maximum = n_categories - 1
    if not 0 <= raw_score <= n_units * maximum:
        raise ValueError("raw_score is outside the response-pattern support.")
    direction = np.asarray(direction, dtype=float)
    if direction.shape != (n_parameters,) or not np.all(np.isfinite(direction)):
        raise ValueError("direction must be finite and match the free-parameter width.")
    utility = np.einsum("ukp,p->uk", design, direction, optimize=True)
    values = np.full((n_units + 1, raw_score + 1), -np.inf, dtype=float)
    predecessor_score = np.full(
        (n_units + 1, raw_score + 1), -1, dtype=np.int32
    )
    predecessor_category = np.full(
        (n_units + 1, raw_score + 1), -1, dtype=np.int16
    )
    values[0, 0] = 0.0
    for unit in range(n_units):
        for previous_score in range(raw_score + 1):
            previous = values[unit, previous_score]
            if not np.isfinite(previous):
                continue
            upper = min(maximum, raw_score - previous_score)
            for category in range(upper + 1):
                next_score = previous_score + category
                candidate = float(previous + utility[unit, category])
                # Categories are visited in ascending order, giving deterministic
                # tie handling without using a tolerance to alter the optimum.
                if candidate > values[unit + 1, next_score]:
                    values[unit + 1, next_score] = candidate
                    predecessor_score[unit + 1, next_score] = previous_score
                    predecessor_category[unit + 1, next_score] = category
    maximum_value = float(values[n_units, raw_score])
    if not np.isfinite(maximum_value):
        raise CMLEExistenceAuditError(
            "Fixed-score support oracle did not reach the requested score."
        )
    categories = np.empty(n_units, dtype=np.int16)
    score = raw_score
    for unit in range(n_units, 0, -1):
        category = int(predecessor_category[unit, score])
        previous = int(predecessor_score[unit, score])
        if category < 0 or previous < 0:
            raise CMLEExistenceAuditError(
                "Fixed-score support oracle predecessor reconstruction failed."
            )
        categories[unit - 1] = category
        score = previous
    statistic = np.sum(
        design[np.arange(n_units), categories, :], axis=0
    )
    reconstructed_value = float(statistic @ direction)
    roundoff = 64.0 * np.finfo(float).eps * max(
        1.0, abs(maximum_value), abs(reconstructed_value)
    )
    if abs(reconstructed_value - maximum_value) > roundoff:
        raise CMLEExistenceAuditError(
            "Fixed-score support oracle value reconstruction failed."
        )
    return maximum_value, statistic, categories


def _fixed_score_support_value(
    design: np.ndarray,
    raw_score: int,
    direction: np.ndarray,
) -> float:
    """Return only the support value using a rolling vectorized DP workspace."""

    n_units, n_categories, n_parameters = design.shape
    raw_score = int(raw_score)
    maximum = n_categories - 1
    if not 0 <= raw_score <= n_units * maximum:
        raise ValueError("raw_score is outside the response-pattern support.")
    direction = np.asarray(direction, dtype=float)
    if direction.shape != (n_parameters,) or not np.all(np.isfinite(direction)):
        raise ValueError("direction must be finite and match the free-parameter width.")
    utility = np.einsum("ukp,p->uk", design, direction, optimize=True)
    values = np.full(raw_score + 1, -np.inf, dtype=float)
    values[0] = 0.0
    for unit in range(n_units):
        updated = np.full(raw_score + 1, -np.inf, dtype=float)
        for category in range(min(maximum, raw_score) + 1):
            candidates = values[: raw_score - category + 1] + utility[
                unit, category
            ]
            np.maximum(
                updated[category:], candidates, out=updated[category:]
            )
        values = updated
    result = float(values[raw_score])
    if not np.isfinite(result):
        raise CMLEExistenceAuditError(
            "Fixed-score support-value oracle did not reach the requested score."
        )
    return result


def cmle_score_stratum_support_maximum(
    design: CMLEDesign,
    *,
    pattern_index: int,
    raw_score: int,
    direction: Sequence[float],
) -> dict[str, object]:
    """Return an exact fixed-score support-function value for validation."""

    if not isinstance(design, CMLEDesign):
        raise TypeError("design must be a prepared CMLEDesign.")
    pattern_index = int(pattern_index)
    if not 0 <= pattern_index < len(design.patterns):
        raise ValueError("pattern_index is outside the prepared CMLE design.")
    vector = np.asarray(direction, dtype=float)
    value, statistic, categories = _fixed_score_support_maximum(
        design.patterns[pattern_index].design,
        int(raw_score),
        vector,
    )
    return {
        "PatternIndex": pattern_index,
        "RawScore": int(raw_score),
        "Maximum": value,
        "Statistic": statistic,
        "Categories": categories,
        "Configurations": _score_configuration_count(
            design.patterns[pattern_index].design.shape[0],
            design.n_categories,
            int(raw_score),
        ),
    }


def json_vector(values: Sequence[float]) -> str:
    return "[" + ",".join(f"{float(value):.17g}" for value in values) + "]"


def _support_inequalities(
    strata: Sequence[_StratumSupport], n_parameters: int
) -> tuple[np.ndarray, int, int]:
    raw_rows = []
    for stratum in strata:
        raw_rows.append(
            stratum.persons * stratum.support - stratum.observed_sum[None, :]
        )
    raw = (
        np.vstack(raw_rows)
        if raw_rows
        else np.empty((0, n_parameters), dtype=float)
    )
    unique = np.unique(raw, axis=0)
    scale = np.max(np.abs(unique), axis=1) if len(unique) else np.empty(0)
    nonzero = scale > 0.0
    normalized = unique[nonzero] / scale[nonzero, None]
    normalized = np.unique(normalized, axis=0)
    return normalized, int(len(raw)), int((~nonzero).sum())


def _boundary_search(
    inequalities: np.ndarray,
    *,
    n_parameters: int,
    tolerance: float,
) -> tuple[dict[str, object], np.ndarray | None]:
    best_objective = -np.inf
    best_direction: np.ndarray | None = None
    statuses: list[int] = []
    messages: list[str] = []
    for coordinate, sign in product(range(n_parameters), (-1.0, 1.0)):
        objective = np.zeros(n_parameters, dtype=float)
        objective[coordinate] = -sign
        result = linprog(
            objective,
            A_ub=inequalities if len(inequalities) else None,
            b_ub=np.zeros(len(inequalities), dtype=float) if len(inequalities) else None,
            bounds=[(-1.0, 1.0)] * n_parameters,
            method="highs",
            options={
                "primal_feasibility_tolerance": float(tolerance),
                "dual_feasibility_tolerance": float(tolerance),
                "ipm_optimality_tolerance": float(tolerance),
            },
        )
        statuses.append(int(result.status))
        messages.append(str(result.message))
        if not result.success or result.x is None:
            continue
        achieved = float(sign * result.x[coordinate])
        if achieved > best_objective:
            best_objective = achieved
            best_direction = np.asarray(result.x, dtype=float)
    solver_passed = bool(statuses and all(status == 0 for status in statuses))
    max_violation = np.nan
    direction_norm = 0.0
    boundary = False
    active_constraints = 0
    if best_direction is not None:
        direction_norm = float(np.max(np.abs(best_direction)))
        if direction_norm > 0.0:
            best_direction = best_direction / direction_norm
        residual = inequalities @ best_direction if len(inequalities) else np.empty(0)
        max_violation = float(np.max(residual)) if len(residual) else 0.0
        residual_tolerance = max(10.0 * float(tolerance), 1e-9)
        boundary = bool(
            solver_passed
            and best_objective >= 0.5
            and max_violation <= residual_tolerance
        )
        active_constraints = int(
            np.sum(np.abs(residual) <= residual_tolerance)
        )
    row = {
        "Tolerance": float(tolerance),
        "SolverPassed": solver_passed,
        "SolverStatuses": ";".join(str(value) for value in statuses),
        "SolverMessage": messages[0] if messages else "no_solver_result",
        "BestCoordinateObjective": float(best_objective),
        "DirectionMaxAbs": direction_norm,
        "MaxConstraintViolation": max_violation,
        "ActiveConstraints": active_constraints,
        "BoundaryDetected": boundary,
        "Direction": json_vector(best_direction) if best_direction is not None else "",
    }
    return row, best_direction if boundary else None


def _normalized_support_row(row: np.ndarray) -> np.ndarray | None:
    row = np.asarray(row, dtype=float)
    scale = float(np.max(np.abs(row))) if len(row) else 0.0
    if not np.isfinite(scale):
        raise CMLEExistenceAuditError("Support oracle returned a non-finite row.")
    if scale == 0.0:
        return None
    normalized = row / scale
    normalized[normalized == 0.0] = 0.0
    return normalized


def _oracle_pass(
    strata: Sequence[_OracleStratum],
    direction: np.ndarray,
) -> tuple[list[dict[str, object]], float, float]:
    """Return one most-violating valid support row per observed stratum."""

    violations: list[dict[str, object]] = []
    maximum_raw = -np.inf
    maximum_normalized = -np.inf
    for stratum in strata:
        value = _fixed_score_support_value(
            stratum.design, stratum.raw_score, direction
        )
        raw_violation = float(
            stratum.persons * value - stratum.observed_sum @ direction
        )
        roundoff = 128.0 * np.finfo(float).eps * max(1.0, abs(raw_violation))
        normalized = None
        categories = np.empty(0, dtype=np.int16)
        normalized_violation = 0.0
        # A maximizing response vector is needed only when a positive support
        # violation can add a cut. Certified/nonpositive oracle passes avoid
        # predecessor allocation and reconstruction entirely.
        if raw_violation > roundoff:
            reconstructed_value, statistic, categories = (
                _fixed_score_support_maximum(
                    stratum.design, stratum.raw_score, direction
                )
            )
            expected_violation = float(
                stratum.persons * reconstructed_value
                - stratum.observed_sum @ direction
            )
            if abs(raw_violation - expected_violation) > roundoff:
                raise CMLEExistenceAuditError(
                    "Support-oracle violation reconstruction failed."
                )
            row = stratum.persons * statistic - stratum.observed_sum
            normalized = _normalized_support_row(row)
            normalized_violation = (
                0.0 if normalized is None else float(normalized @ direction)
            )
        maximum_raw = max(maximum_raw, raw_violation)
        maximum_normalized = max(maximum_normalized, normalized_violation)
        violations.append(
            {
                "PatternIndex": stratum.pattern_index,
                "RawScore": stratum.raw_score,
                "RawViolation": raw_violation,
                "NormalizedViolation": normalized_violation,
                "Row": normalized,
                "Categories": categories,
            }
        )
    if not strata:
        maximum_raw = 0.0
        maximum_normalized = 0.0
    return violations, float(maximum_raw), float(maximum_normalized)


def _cut_key(row: np.ndarray) -> bytes:
    contiguous = np.ascontiguousarray(row, dtype=np.float64)
    return contiguous.tobytes()


def _boundary_search_oracle(
    strata: Sequence[_OracleStratum],
    *,
    n_parameters: int,
    tolerance: float,
    max_rounds: int,
    max_cuts: int,
) -> tuple[dict[str, object], np.ndarray | None, np.ndarray, pd.DataFrame]:
    """Certify the coordinate boundary LPs by exact constraint generation."""

    cuts: list[np.ndarray] = []
    cut_keys: set[bytes] = set()
    history: list[dict[str, object]] = []
    statuses: list[int] = []
    messages: list[str] = []
    best_objective = -np.inf
    best_direction: np.ndarray | None = None
    best_raw_violation = np.nan
    best_normalized_violation = np.nan
    oracle_calls = 0
    failure_reason = ""
    residual_tolerance = max(10.0 * float(tolerance), 1e-9)

    for coordinate, sign in product(range(n_parameters), (-1.0, 1.0)):
        certified = False
        for round_index in range(1, max_rounds + 1):
            inequality_matrix = (
                np.vstack(cuts)
                if cuts
                else np.empty((0, n_parameters), dtype=float)
            )
            objective = np.zeros(n_parameters, dtype=float)
            objective[coordinate] = -sign
            result = linprog(
                objective,
                A_ub=inequality_matrix if len(inequality_matrix) else None,
                b_ub=(
                    np.zeros(len(inequality_matrix), dtype=float)
                    if len(inequality_matrix)
                    else None
                ),
                bounds=[(-1.0, 1.0)] * n_parameters,
                method="highs",
                options={
                    "primal_feasibility_tolerance": float(tolerance),
                    "dual_feasibility_tolerance": float(tolerance),
                    "ipm_optimality_tolerance": float(tolerance),
                },
            )
            statuses.append(int(result.status))
            messages.append(str(result.message))
            if not result.success or result.x is None:
                failure_reason = "coordinate_lp_failed"
                history.append(
                    {
                        "Tolerance": float(tolerance),
                        "Coordinate": coordinate,
                        "Sign": sign,
                        "Round": round_index,
                        "SolverStatus": int(result.status),
                        "CutsBefore": len(cuts),
                        "CutsAdded": 0,
                        "MaxRawViolation": np.nan,
                        "MaxNormalizedViolation": np.nan,
                        "Certified": False,
                        "FailureReason": failure_reason,
                    }
                )
                break
            candidate = np.asarray(result.x, dtype=float)
            oracle_rows, raw_violation, normalized_violation = _oracle_pass(
                strata, candidate
            )
            oracle_calls += len(strata)
            added = 0
            if normalized_violation > residual_tolerance:
                for oracle_row in oracle_rows:
                    row = oracle_row["Row"]
                    if (
                        row is None
                        or float(oracle_row["NormalizedViolation"])
                        <= residual_tolerance
                    ):
                        continue
                    key = _cut_key(row)
                    if key in cut_keys:
                        continue
                    cuts.append(row)
                    cut_keys.add(key)
                    added += 1
                    if len(cuts) > max_cuts:
                        failure_reason = "oracle_cut_cap_exceeded"
                        break
                if not failure_reason and added == 0:
                    failure_reason = "oracle_positive_violation_without_new_cut"
            else:
                certified = True
            history.append(
                {
                    "Tolerance": float(tolerance),
                    "Coordinate": coordinate,
                    "Sign": sign,
                    "Round": round_index,
                    "SolverStatus": int(result.status),
                    "CutsBefore": len(cuts) - added,
                    "CutsAdded": added,
                    "MaxRawViolation": raw_violation,
                    "MaxNormalizedViolation": normalized_violation,
                    "Certified": certified,
                    "FailureReason": failure_reason,
                }
            )
            if failure_reason or certified:
                break
        if failure_reason:
            break
        if not certified:
            failure_reason = "oracle_round_cap_exceeded"
            break
        achieved = float(sign * candidate[coordinate])
        if achieved > best_objective:
            best_objective = achieved
            best_direction = candidate.copy()
            best_raw_violation = raw_violation
            best_normalized_violation = normalized_violation

    solver_passed = bool(
        not failure_reason and statuses and all(status == 0 for status in statuses)
    )
    direction_norm = 0.0
    boundary = False
    active_constraints = 0
    if best_direction is not None:
        direction_norm = float(np.max(np.abs(best_direction)))
        if direction_norm > 0.0:
            best_direction = best_direction / direction_norm
        _, best_raw_violation, best_normalized_violation = _oracle_pass(
            strata, best_direction
        )
        oracle_calls += len(strata)
        inequality_matrix = (
            np.vstack(cuts)
            if cuts
            else np.empty((0, n_parameters), dtype=float)
        )
        residual = (
            inequality_matrix @ best_direction
            if len(inequality_matrix)
            else np.empty(0)
        )
        active_constraints = int(
            np.sum(np.abs(residual) <= residual_tolerance)
        )
        boundary = bool(
            solver_passed
            and best_objective >= 0.5
            and best_normalized_violation <= residual_tolerance
        )
    inequality_matrix = (
        np.vstack(cuts)
        if cuts
        else np.empty((0, n_parameters), dtype=float)
    )
    row = {
        "Tolerance": float(tolerance),
        "SolverPassed": solver_passed,
        "SolverStatuses": ";".join(str(value) for value in statuses),
        "SolverMessage": messages[0] if messages else "no_solver_result",
        "BestCoordinateObjective": float(best_objective),
        "DirectionMaxAbs": direction_norm,
        "MaxRawOracleViolation": best_raw_violation,
        "MaxNormalizedOracleViolation": best_normalized_violation,
        "ActiveGeneratedConstraints": active_constraints,
        "GeneratedConstraints": len(cuts),
        "CuttingPlaneRounds": len(history),
        "OracleCalls": oracle_calls,
        "BoundaryDetected": boundary,
        "FailureReason": failure_reason,
        "Direction": (
            json_vector(best_direction) if best_direction is not None else ""
        ),
    }
    return (
        row,
        best_direction if boundary else None,
        inequality_matrix,
        pd.DataFrame(history),
    )


def _directional_trace(
    design: CMLEDesign,
    direction: np.ndarray,
    radii: Sequence[float],
) -> pd.DataFrame:
    rows = []
    for radius in radii:
        parameters = float(radius) * direction
        try:
            objective, gradient = cmle_objective_value_grad(parameters, design)
            rows.append(
                {
                    "Radius": float(radius),
                    "Objective": float(objective),
                    "ConditionalLogLik": -float(objective),
                    "ObjectiveDirectionalDerivative": float(gradient @ direction),
                    "Finite": bool(
                        np.isfinite(objective) and np.all(np.isfinite(gradient))
                    ),
                    "FailureReason": "",
                }
            )
        except FloatingPointError as exc:
            rows.append(
                {
                    "Radius": float(radius),
                    "Objective": np.nan,
                    "ConditionalLogLik": np.nan,
                    "ObjectiveDirectionalDerivative": np.nan,
                    "Finite": False,
                    "FailureReason": f"{type(exc).__name__}: {exc}",
                }
            )
    return pd.DataFrame(rows)


def audit_cmle_finite_mle(
    design: CMLEDesign,
    *,
    tolerances: Iterable[float] = DEFAULT_EXISTENCE_TOLERANCES,
    directional_radii: Sequence[float] = DEFAULT_DIRECTIONAL_RADII,
    max_configurations: int = DEFAULT_MAX_CONFIGURATIONS,
    max_constraints: int = DEFAULT_MAX_CONSTRAINTS,
    max_bytes: int = DEFAULT_MAX_BYTES,
) -> dict[str, object]:
    """Audit full-rank exact-CMLE finite-maximizer existence by convex support."""

    if not isinstance(design, CMLEDesign):
        raise TypeError("design must be a prepared CMLEDesign.")
    tolerance_grid = _validate_tolerances(tolerances)
    radii = tuple(float(value) for value in directional_radii)
    if (
        not radii
        or any(not np.isfinite(value) or value < 0.0 for value in radii)
        or any(right <= left for left, right in zip(radii, radii[1:]))
    ):
        raise ValueError("directional_radii must be finite, non-negative, and increasing.")
    max_configurations = _positive_integer(
        max_configurations, "max_configurations"
    )
    max_constraints = _positive_integer(max_constraints, "max_constraints")
    max_bytes = _positive_integer(max_bytes, "max_bytes")
    base = {
        "Method": "exact_conditional_convex_support_lp",
        "KParams": int(design.n_parameters),
        "PrefitRank": int(design.audit.get("conditional_rank", 0)),
        "PrefitNullity": int(design.audit.get("conditional_nullity", 0)),
        "ToleranceGrid": ";".join(f"{value:.17g}" for value in tolerance_grid),
        "EnumerationComplete": False,
        "LPComplete": False,
        "BoundaryDetected": pd.NA,
        "FiniteMLESupported": False,
        "ExistenceQualified": False,
        "Status": "unavailable",
        "Reason": "",
    }
    empty = pd.DataFrame()
    if not bool(design.audit.get("eligible", False)) or int(
        design.audit.get("conditional_nullity", 0)
    ):
        base.update(
            {
                "Status": "structural_nonidentification",
                "Reason": "exact_conditional_information_rank_deficient",
                "BoundaryDetected": pd.NA,
            }
        )
        return {
            "summary": pd.DataFrame([base]),
            "strata": empty,
            "lp_tolerances": empty,
            "directional_trace": empty,
            "direction": None,
            "inequalities": None,
        }
    try:
        strata, stratum_table, failure = _enumerate_strata(
            design,
            max_configurations=max_configurations,
            max_constraints=max_constraints,
            max_bytes=max_bytes,
        )
    except (CMLEExistenceAuditError, MemoryError) as exc:
        base["Reason"] = f"{type(exc).__name__}: {exc}"
        return {
            "summary": pd.DataFrame([base]),
            "strata": empty,
            "lp_tolerances": empty,
            "directional_trace": empty,
            "direction": None,
            "inequalities": None,
        }
    if failure is not None:
        base.update(
            {
                "Reason": failure,
                "Strata": int(len(stratum_table)),
                "Configurations": int(
                    pd.to_numeric(
                        stratum_table.get("Configurations", pd.Series(dtype=float)),
                        errors="coerce",
                    ).sum()
                ),
            }
        )
        return {
            "summary": pd.DataFrame([base]),
            "strata": stratum_table,
            "lp_tolerances": empty,
            "directional_trace": empty,
            "direction": None,
            "inequalities": None,
        }
    inequalities, raw_constraints, zero_constraints = _support_inequalities(
        strata, design.n_parameters
    )
    base.update(
        {
            "EnumerationComplete": True,
            "Strata": len(strata),
            "Configurations": int(sum(row.configurations for row in strata)),
            "UniqueSupportPoints": int(sum(len(row.support) for row in strata)),
            "RawConstraints": raw_constraints,
            "ZeroConstraintsRemoved": zero_constraints,
            "UniqueNormalizedConstraints": len(inequalities),
        }
    )
    lp_rows = []
    directions: list[np.ndarray | None] = []
    for tolerance in tolerance_grid:
        row, direction = _boundary_search(
            inequalities,
            n_parameters=design.n_parameters,
            tolerance=tolerance,
        )
        lp_rows.append(row)
        directions.append(direction)
    lp_table = pd.DataFrame(lp_rows)
    solver_complete = bool(lp_table["SolverPassed"].all())
    boundary_values = lp_table["BoundaryDetected"].astype(bool).tolist()
    base["LPComplete"] = solver_complete
    direction: np.ndarray | None = None
    if not solver_complete:
        base.update(
            {
                "Status": "lp_unavailable",
                "Reason": "one_or_more_lp_searches_failed",
                "BoundaryDetected": pd.NA,
            }
        )
    elif all(boundary_values):
        preferred = min(
            range(len(tolerance_grid)),
            key=lambda index: abs(np.log10(tolerance_grid[index]) + 9.0),
        )
        direction = directions[preferred]
        base.update(
            {
                "Status": "boundary_no_finite_cmle",
                "Reason": "observed_sufficient_statistic_on_convex_support_boundary",
                "BoundaryDetected": True,
            }
        )
    elif not any(boundary_values):
        base.update(
            {
                "Status": "interior_finite_cmle_supported",
                "Reason": "no_nonzero_supporting_direction_found",
                "BoundaryDetected": False,
                "FiniteMLESupported": True,
                "ExistenceQualified": True,
            }
        )
    else:
        base.update(
            {
                "Status": "tolerance_unstable",
                "Reason": "boundary_decision_changed_across_tolerance_grid",
                "BoundaryDetected": pd.NA,
            }
        )
    trace = (
        _directional_trace(design, direction, radii)
        if direction is not None
        else empty
    )
    if direction is not None and not trace.empty:
        finite = trace.loc[trace["Finite"]]
        differences = np.diff(finite["Objective"].to_numpy(dtype=float))
        monotone_tolerance = 1e-10 * max(
            1.0, float(np.max(np.abs(finite["Objective"])))
        )
        base["DirectionalTraceFinite"] = bool(trace["Finite"].all())
        base["DirectionalObjectiveNonincreasing"] = bool(
            np.all(differences <= monotone_tolerance)
        )
        base["DirectionalDerivativeAtLargestRadius"] = float(
            finite.iloc[-1]["ObjectiveDirectionalDerivative"]
        )
    else:
        base["DirectionalTraceFinite"] = pd.NA
        base["DirectionalObjectiveNonincreasing"] = pd.NA
        base["DirectionalDerivativeAtLargestRadius"] = np.nan
    return {
        "summary": pd.DataFrame([base]),
        "strata": stratum_table,
        "lp_tolerances": lp_table,
        "directional_trace": trace,
        "direction": direction,
        "inequalities": inequalities,
    }


def audit_cmle_finite_mle_oracle(
    design: CMLEDesign,
    *,
    tolerances: Iterable[float] = DEFAULT_EXISTENCE_TOLERANCES,
    directional_radii: Sequence[float] = DEFAULT_DIRECTIONAL_RADII,
    max_cutting_plane_rounds: int = DEFAULT_MAX_CUTTING_PLANE_ROUNDS,
    max_cuts: int = DEFAULT_MAX_ORACLE_CUTS,
    max_oracle_state_cells: int = DEFAULT_MAX_ORACLE_STATE_CELLS,
    max_bytes: int = DEFAULT_MAX_BYTES,
) -> dict[str, object]:
    """Audit finite-CMLE existence with an exact fixed-score support oracle."""

    if not isinstance(design, CMLEDesign):
        raise TypeError("design must be a prepared CMLEDesign.")
    tolerance_grid = _validate_tolerances(tolerances)
    radii = tuple(float(value) for value in directional_radii)
    if (
        not radii
        or any(not np.isfinite(value) or value < 0.0 for value in radii)
        or any(right <= left for left, right in zip(radii, radii[1:]))
    ):
        raise ValueError("directional_radii must be finite, non-negative, and increasing.")
    max_rounds = _positive_integer(
        max_cutting_plane_rounds, "max_cutting_plane_rounds"
    )
    max_cuts = _positive_integer(max_cuts, "max_cuts")
    max_state_cells = _positive_integer(
        max_oracle_state_cells, "max_oracle_state_cells"
    )
    max_bytes = _positive_integer(max_bytes, "max_bytes")
    base = {
        "Method": "exact_conditional_support_oracle_cutting_plane_lp",
        "KParams": int(design.n_parameters),
        "PrefitRank": int(design.audit.get("conditional_rank", 0)),
        "PrefitNullity": int(design.audit.get("conditional_nullity", 0)),
        "ToleranceGrid": ";".join(f"{value:.17g}" for value in tolerance_grid),
        "OraclePrepared": False,
        "OracleComplete": False,
        "LPComplete": False,
        "BoundaryDetected": pd.NA,
        "FiniteMLESupported": False,
        "ExistenceQualified": False,
        "Status": "unavailable",
        "Reason": "",
    }
    empty = pd.DataFrame()
    if not bool(design.audit.get("eligible", False)) or int(
        design.audit.get("conditional_nullity", 0)
    ):
        base.update(
            {
                "Status": "structural_nonidentification",
                "Reason": "exact_conditional_information_rank_deficient",
                "BoundaryDetected": pd.NA,
            }
        )
        return {
            "summary": pd.DataFrame([base]),
            "strata": empty,
            "lp_tolerances": empty,
            "cut_history": empty,
            "directional_trace": empty,
            "direction": None,
            "inequalities": None,
            "inequalities_by_tolerance": [],
        }
    try:
        strata, stratum_table, failure = _oracle_strata(
            design,
            max_state_cells=max_state_cells,
            max_bytes=max_bytes,
        )
    except (CMLEExistenceAuditError, MemoryError) as exc:
        base["Reason"] = f"{type(exc).__name__}: {exc}"
        return {
            "summary": pd.DataFrame([base]),
            "strata": empty,
            "lp_tolerances": empty,
            "cut_history": empty,
            "directional_trace": empty,
            "direction": None,
            "inequalities": None,
            "inequalities_by_tolerance": [],
        }
    theoretical_configurations = int(
        pd.to_numeric(
            stratum_table.get(
                "TheoreticalConfigurations", pd.Series(dtype=float)
            ),
            errors="coerce",
        ).sum()
    )
    state_cells = int(
        pd.to_numeric(
            stratum_table.get("OracleStateCells", pd.Series(dtype=float)),
            errors="coerce",
        ).sum()
    )
    base.update(
        {
            "Strata": len(stratum_table),
            "TheoreticalConfigurations": theoretical_configurations,
            "OracleStateCells": state_cells,
        }
    )
    if failure is not None:
        base["Reason"] = failure
        return {
            "summary": pd.DataFrame([base]),
            "strata": stratum_table,
            "lp_tolerances": empty,
            "cut_history": empty,
            "directional_trace": empty,
            "direction": None,
            "inequalities": None,
            "inequalities_by_tolerance": [],
        }
    base["OraclePrepared"] = True
    lp_rows: list[dict[str, object]] = []
    directions: list[np.ndarray | None] = []
    inequality_matrices: list[np.ndarray] = []
    histories: list[pd.DataFrame] = []
    try:
        for tolerance in tolerance_grid:
            row, direction, inequalities, history = _boundary_search_oracle(
                strata,
                n_parameters=design.n_parameters,
                tolerance=tolerance,
                max_rounds=max_rounds,
                max_cuts=max_cuts,
            )
            lp_rows.append(row)
            directions.append(direction)
            inequality_matrices.append(inequalities)
            histories.append(history)
    except (CMLEExistenceAuditError, MemoryError) as exc:
        base["Reason"] = f"{type(exc).__name__}: {exc}"
        return {
            "summary": pd.DataFrame([base]),
            "strata": stratum_table,
            "lp_tolerances": pd.DataFrame(lp_rows),
            "cut_history": (
                pd.concat(histories, ignore_index=True) if histories else empty
            ),
            "directional_trace": empty,
            "direction": None,
            "inequalities": None,
            "inequalities_by_tolerance": inequality_matrices,
        }
    lp_table = pd.DataFrame(lp_rows)
    cut_history = (
        pd.concat(histories, ignore_index=True) if histories else empty
    )
    solver_complete = bool(lp_table["SolverPassed"].all())
    boundary_values = lp_table["BoundaryDetected"].astype(bool).tolist()
    base.update(
        {
            "OracleComplete": solver_complete,
            "LPComplete": solver_complete,
            "GeneratedConstraintsMax": int(
                lp_table["GeneratedConstraints"].max()
            ),
            "GeneratedConstraintsTotal": int(
                lp_table["GeneratedConstraints"].sum()
            ),
            "CuttingPlaneRoundsTotal": int(
                lp_table["CuttingPlaneRounds"].sum()
            ),
            "OracleCallsTotal": int(lp_table["OracleCalls"].sum()),
        }
    )
    direction: np.ndarray | None = None
    preferred = min(
        range(len(tolerance_grid)),
        key=lambda index: abs(np.log10(tolerance_grid[index]) + 9.0),
    )
    if not solver_complete:
        reasons = sorted(
            {
                str(value)
                for value in lp_table["FailureReason"]
                if str(value)
            }
        )
        base.update(
            {
                "Status": "oracle_unavailable",
                "Reason": ";".join(reasons) or "one_or_more_oracle_lp_searches_failed",
                "BoundaryDetected": pd.NA,
            }
        )
    elif all(boundary_values):
        direction = directions[preferred]
        base.update(
            {
                "Status": "boundary_no_finite_cmle",
                "Reason": "observed_sufficient_statistic_on_convex_support_boundary",
                "BoundaryDetected": True,
            }
        )
    elif not any(boundary_values):
        base.update(
            {
                "Status": "interior_finite_cmle_supported",
                "Reason": "no_nonzero_supporting_direction_found",
                "BoundaryDetected": False,
                "FiniteMLESupported": True,
                "ExistenceQualified": True,
            }
        )
    else:
        base.update(
            {
                "Status": "tolerance_unstable",
                "Reason": "boundary_decision_changed_across_tolerance_grid",
                "BoundaryDetected": pd.NA,
            }
        )
    trace = (
        _directional_trace(design, direction, radii)
        if direction is not None
        else empty
    )
    if direction is not None and not trace.empty:
        finite = trace.loc[trace["Finite"]]
        differences = np.diff(finite["Objective"].to_numpy(dtype=float))
        monotone_tolerance = 1e-10 * max(
            1.0, float(np.max(np.abs(finite["Objective"])))
        )
        trace_finite = bool(trace["Finite"].all())
        trace_nonincreasing = bool(
            len(finite) == len(trace)
            and np.all(differences <= monotone_tolerance)
        )
        base["DirectionalTraceFinite"] = trace_finite
        base["DirectionalObjectiveNonincreasing"] = trace_nonincreasing
        base["DirectionalDerivativeAtLargestRadius"] = float(
            finite.iloc[-1]["ObjectiveDirectionalDerivative"]
        )
        if not trace_finite or not trace_nonincreasing:
            base.update(
                {
                    "Status": "directional_verification_failed",
                    "Reason": "reported_boundary_direction_failed_direct_objective_trace",
                    "BoundaryDetected": pd.NA,
                    "FiniteMLESupported": False,
                    "ExistenceQualified": False,
                }
            )
    else:
        base["DirectionalTraceFinite"] = pd.NA
        base["DirectionalObjectiveNonincreasing"] = pd.NA
        base["DirectionalDerivativeAtLargestRadius"] = np.nan
    return {
        "summary": pd.DataFrame([base]),
        "strata": stratum_table,
        "lp_tolerances": lp_table,
        "cut_history": cut_history,
        "directional_trace": trace,
        "direction": direction,
        "inequalities": inequality_matrices[preferred],
        "inequalities_by_tolerance": inequality_matrices,
    }


__all__ = [
    "CMLEExistenceAuditError",
    "DEFAULT_DIRECTIONAL_RADII",
    "DEFAULT_EXISTENCE_TOLERANCES",
    "DEFAULT_MAX_BYTES",
    "DEFAULT_MAX_CUTTING_PLANE_ROUNDS",
    "DEFAULT_MAX_CONFIGURATIONS",
    "DEFAULT_MAX_CONSTRAINTS",
    "DEFAULT_MAX_ORACLE_CUTS",
    "DEFAULT_MAX_ORACLE_STATE_CELLS",
    "audit_cmle_finite_mle",
    "audit_cmle_finite_mle_oracle",
    "cmle_score_stratum_support_maximum",
]
