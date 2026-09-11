"""App-specific Q-fit/Q-check adapter layered above frozen stationarity code."""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
from pathlib import Path
import operator
import platform
import sys
from typing import Any

import numpy as np
import scipy

import streamlit_app as app
import mfrm_app.mml_engine_v2 as mml_engine_v2_module
import mfrm_app.mml_qualification_batch as qualification_batch_module
import mfrm_app.mml_quadrature_sensitivity as quadrature_module
import mfrm_app.mml_stationarity as stationarity_module
from mfrm_app.mml_engine_v2 import (
    FreeSdStationarityRun,
    StationarityContract,
    run_free_sd_stationarity_v2,
)
from mfrm_app.mml_quadrature_sensitivity import (
    QuadratureSensitivityAssessment,
    QuadratureSensitivityContract,
    assess_quadrature_sensitivity,
)
from mfrm_app.mml_stationarity import JointPolishOptions
from mfrm_app.mml_qualification_batch import QuadratureObjectiveEvaluator
from validation.mml_free_sd_stationarity_adapter import (
    AppFreeSdProblem,
    prepare_app_free_sd_problem,
)


_IMPLEMENTATION_SOURCE_PATHS = {
    "streamlit_app_sha256": Path(app.__file__),
    "mml_engine_v2_sha256": Path(mml_engine_v2_module.__file__),
    "mml_quadrature_sensitivity_sha256": Path(quadrature_module.__file__),
    "mml_qualification_batch_sha256": Path(qualification_batch_module.__file__),
    "mml_stationarity_sha256": Path(stationarity_module.__file__),
    "stationarity_adapter_sha256": Path(__file__).with_name(
        "mml_free_sd_stationarity_adapter.py"
    ),
    "quadrature_adapter_sha256": Path(__file__),
}
_LOADED_IMPLEMENTATION_SOURCE_SHA256 = {
    name: hashlib.sha256(path.read_bytes()).hexdigest()
    for name, path in _IMPLEMENTATION_SOURCE_PATHS.items()
}


@dataclass(frozen=True)
class AppQuadratureSensitivityBundle:
    problem_digest: str
    primary_problem: AppFreeSdProblem
    sensitivity_problem: AppFreeSdProblem
    primary_run: FreeSdStationarityRun
    sensitivity_run: FreeSdStationarityRun
    assessment: QuadratureSensitivityAssessment
    objective_evaluator_implementation_sha256: str
    objective_evaluator: QuadratureObjectiveEvaluator


def _jsonable(value: Any) -> Any:
    """Canonical local-v2 encoding; reject unmodelled nonfinite values."""

    if isinstance(value, dict):
        return {
            str(key): _jsonable(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        array = np.asarray(value)
        if array.dtype.kind in "fc" and not np.isfinite(array).all():
            raise ValueError("nonfinite array cannot enter problem identity")
        if array.dtype.kind == "f":
            canonical = np.ascontiguousarray(array, dtype="<f8")
            dtype = "<f8"
        elif array.dtype.kind in "iu":
            canonical = np.ascontiguousarray(array, dtype="<i8")
            dtype = "<i8"
        elif array.dtype.kind == "b":
            canonical = np.ascontiguousarray(array, dtype=np.uint8)
            dtype = "bool-u1"
        elif array.dtype.kind in "US":
            return {
                "encoding": "unicode-values-v1",
                "shape": list(array.shape),
                "values": _jsonable(array.tolist()),
            }
        else:
            raise ValueError("unsupported array dtype in problem identity")
        return {
            "encoding": "canonical-numeric-array-v1",
            "dtype": dtype,
            "shape": list(canonical.shape),
            "sha256": hashlib.sha256(canonical.tobytes(order="C")).hexdigest(),
        }
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        result = float(value)
        if not np.isfinite(result):
            raise ValueError("nonfinite scalar cannot enter problem identity")
        return result
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if value is None or isinstance(value, (str, int)):
        return value
    raise TypeError(f"unsupported value in problem identity: {type(value).__name__}")


def _facet_specs_identity(config: dict[str, Any]) -> dict[str, Any]:
    specs = config.get("facet_specs")
    if not isinstance(specs, dict):
        raise ValueError("facet_specs are unavailable for problem identity")
    result: dict[str, Any] = {}
    for facet in config.get("facet_names", []):
        spec = specs.get(facet)
        if not isinstance(spec, dict):
            raise ValueError(f"facet spec is missing for {facet}")
        item = dict(spec)
        anchors = np.asarray(item.pop("anchors", []), dtype=float)
        if anchors.ndim != 1 or np.isinf(anchors).any():
            raise ValueError(f"facet anchors are malformed for {facet}")
        item["anchors"] = [
            {"state": "free"}
            if np.isnan(value)
            else {"state": "fixed", "value": float(value)}
            for value in anchors
        ]
        result[str(facet)] = _jsonable(item)
    return result


def _effective_facet_signs(config: dict[str, Any]) -> dict[str, int]:
    facets = tuple(map(str, config.get("facet_names", [])))
    raw = config.get("facet_signs")
    if not isinstance(raw, dict) or set(map(str, raw)) != set(facets):
        raise ValueError("facet_signs must exactly cover facet_names")
    signs: dict[str, int] = {}
    for facet in facets:
        value = raw.get(facet)
        if isinstance(value, (bool, np.bool_)):
            raise ValueError("facet signs must be numeric -1 or +1")
        try:
            numeric = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError("facet signs must be numeric -1 or +1") from exc
        if not np.isfinite(numeric) or numeric not in {-1.0, 1.0}:
            raise ValueError("facet signs must be numeric -1 or +1")
        signs[facet] = int(numeric)
    positive = {str(value) for value in config.get("positive_facets", [])}
    if positive != {facet for facet, sign in signs.items() if sign == 1}:
        raise ValueError("positive_facets and facet_signs are inconsistent")
    return signs


def _source_sha256(path: str | Path) -> str:
    source = Path(path)
    if not source.is_file():
        raise ValueError(f"objective source is unavailable: {source}")
    return hashlib.sha256(source.read_bytes()).hexdigest()


def app_free_sd_evaluator_implementation_sha256() -> str:
    """Bind every source and runtime component used by objective replay."""

    current_sources = {
        name: _source_sha256(path)
        for name, path in _IMPLEMENTATION_SOURCE_PATHS.items()
    }
    if current_sources != _LOADED_IMPLEMENTATION_SOURCE_SHA256:
        raise ValueError("objective evaluator source changed after module import")
    payload = {
        "schema_version": "app-free-sd-objective-evaluator-implementation-v1",
        "domain": "MFRM_APP_FREE_SD_OBJECTIVE_EVALUATOR_IMPLEMENTATION",
        "sources": current_sources,
        "runtime": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "byteorder": sys.byteorder,
        },
    }
    encoded = json.dumps(
        _jsonable(payload),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def app_free_sd_problem_digest(
    result: dict[str, Any],
    problem: AppFreeSdProblem,
) -> str:
    """Hash the exact Q-independent likelihood problem under local schema v2."""

    prep = result.get("prep")
    if not isinstance(prep, dict):
        raise ValueError("native result lacks prepared data for problem identity")
    data = prep.get("data")
    if data is None or not hasattr(data, "to_csv"):
        raise ValueError("prepared observation data are unavailable for problem identity")
    config = problem.config
    facet_signs = _effective_facet_signs(config)
    identity_columns = [
        name
        for name in (
            "Person",
            *config.get("facet_names", []),
            "Score",
            "Weight",
        )
        if name in data.columns
    ]
    observation_bytes = data[identity_columns].to_csv(
        index=False,
        lineterminator="\n",
    ).encode("utf-8")
    population = config.get("population_model", {})
    design = population.get("X")
    design_identity = None
    if design is not None:
        design_array = np.asarray(design, dtype=float)
        if design_array.ndim != 2 or not np.isfinite(design_array).all():
            raise ValueError("population design must be a finite two-dimensional array")
        design_identity = _jsonable(design_array)
    idx_identity = _jsonable(problem.idx)
    bounds_identity = _jsonable(
        {
            "structural": problem.structural_bounds,
            "sigma": problem.sigma_bounds,
        }
    )
    payload = {
        "schema_version": "app-free-sd-mml-problem-canonical-local-v2",
        "domain": "MFRM_APP_FREE_SD_FINITE_Q_LIKELIHOOD_PROBLEM",
        "objective": "unregularized_normal_population_marginal_negative_log_likelihood",
        "model": config.get("model"),
        "method": config.get("method"),
        "facet_names": config.get("facet_names"),
        "facet_levels": config.get("facet_levels"),
        "step_facet": config.get("step_facet"),
        "n_cat": config.get("n_cat"),
        "rating_min": config.get("rating_min"),
        "rating_max": config.get("rating_max"),
        "noncenter_facet": config.get("noncenter_facet"),
        "dummy_facets": sorted(map(str, config.get("dummy_facets", []) or [])),
        "positive_facets": sorted(
            map(str, config.get("positive_facets", []) or [])
        ),
        "facet_signs": facet_signs,
        "facet_specs": _facet_specs_identity(config),
        "parameter_sizes": list(problem.sizes.items()),
        "parameter_bounds": bounds_identity,
        "indices": idx_identity,
        "observations": problem.observations,
        "population": {
            "enabled": population.get("enabled"),
            "formula": population.get("formula"),
            "columns": population.get("columns"),
            "term_types": population.get("term_types"),
            "design": design_identity,
        },
        "observations_sha256": hashlib.sha256(observation_bytes).hexdigest(),
        "objective_evaluator_implementation_sha256": (
            app_free_sd_evaluator_implementation_sha256()
        ),
    }
    encoded = json.dumps(
        _jsonable(payload),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def run_app_free_sd_quadrature_sensitivity(
    result: dict[str, Any],
    *,
    primary_quadrature_points: int,
    sensitivity_quadrature_points: int,
    stationarity_contract: StationarityContract,
    sensitivity_contract: QuadratureSensitivityContract,
    options: JointPolishOptions | None = None,
) -> AppQuadratureSensitivityBundle:
    """Run and cross-evaluate two native app quadrature objectives."""

    stationarity_contract.validate()
    sensitivity_contract.validate()
    settings = options or JointPolishOptions()
    settings.validate()
    try:
        primary_points = operator.index(primary_quadrature_points)
        sensitivity_points = operator.index(sensitivity_quadrature_points)
    except TypeError as exc:
        raise ValueError("quadrature points must be exact integers") from exc
    if isinstance(primary_quadrature_points, (bool, np.bool_)) or isinstance(
        sensitivity_quadrature_points,
        (bool, np.bool_),
    ):
        raise ValueError("quadrature points must be exact integers")
    if (
        primary_points != sensitivity_contract.primary_quadrature_points
        or sensitivity_points != sensitivity_contract.sensitivity_quadrature_points
    ):
        raise ValueError("quadrature arguments differ from the sensitivity contract")
    implementation_start = app_free_sd_evaluator_implementation_sha256()
    primary_problem = prepare_app_free_sd_problem(
        result,
        quadrature_points=primary_points,
    )
    digest = app_free_sd_problem_digest(result, primary_problem)
    primary_run = run_free_sd_stationarity_v2(
        primary_problem.structural_start,
        primary_problem.sigma_start,
        primary_problem.value,
        primary_problem.value_gradient,
        observations=primary_problem.observations,
        structural_bounds=primary_problem.structural_bounds,
        sigma_bounds=primary_problem.sigma_bounds,
        options=settings,
        constraint_residual_function=primary_problem.constraint_residual,
    )
    sensitivity_base = prepare_app_free_sd_problem(
        result,
        quadrature_points=sensitivity_points,
    )
    if app_free_sd_problem_digest(result, sensitivity_base) != digest:
        raise ValueError("Q-fit and Q-check problem identities differ")
    sensitivity_problem = replace(
        sensitivity_base,
        structural_start=primary_run.restart_polish.structural_parameters,
        sigma_start=primary_run.restart_polish.sigma,
    )
    sensitivity_run = run_free_sd_stationarity_v2(
        sensitivity_problem.structural_start,
        sensitivity_problem.sigma_start,
        sensitivity_problem.value,
        sensitivity_problem.value_gradient,
        observations=sensitivity_problem.observations,
        structural_bounds=sensitivity_problem.structural_bounds,
        sigma_bounds=sensitivity_problem.sigma_bounds,
        options=settings,
        constraint_residual_function=sensitivity_problem.constraint_residual,
    )
    if (
        app_free_sd_problem_digest(result, primary_problem) != digest
        or app_free_sd_problem_digest(result, sensitivity_problem) != digest
    ):
        raise ValueError("Q-fit/Q-check problem identity changed during evaluation")
    assessment = assess_quadrature_sensitivity(
        problem_digest=digest,
        primary_quadrature_points=primary_problem.quadrature_points,
        sensitivity_quadrature_points=sensitivity_problem.quadrature_points,
        primary_run=primary_run,
        sensitivity_run=sensitivity_run,
        primary_value_function=primary_problem.value,
        sensitivity_value_function=sensitivity_problem.value,
        stationarity_contract=stationarity_contract,
        sensitivity_contract=sensitivity_contract,
    )
    if (
        app_free_sd_problem_digest(result, primary_problem) != digest
        or app_free_sd_problem_digest(result, sensitivity_problem) != digest
    ):
        raise ValueError("Q-fit/Q-check problem identity changed during assessment")
    if app_free_sd_evaluator_implementation_sha256() != implementation_start:
        raise ValueError("objective evaluator implementation changed during evaluation")
    objective_evaluator = QuadratureObjectiveEvaluator(
        problem_digest=digest,
        primary_quadrature_points=primary_problem.quadrature_points,
        sensitivity_quadrature_points=sensitivity_problem.quadrature_points,
        primary_value_function=primary_problem.value,
        sensitivity_value_function=sensitivity_problem.value,
    )
    objective_evaluator.validate()
    return AppQuadratureSensitivityBundle(
        problem_digest=digest,
        primary_problem=primary_problem,
        sensitivity_problem=sensitivity_problem,
        primary_run=primary_run,
        sensitivity_run=sensitivity_run,
        assessment=assessment,
        objective_evaluator_implementation_sha256=(
            implementation_start
        ),
        objective_evaluator=objective_evaluator,
    )


def resolve_app_free_sd_objective_evaluator(
    result: dict[str, Any],
    *,
    primary_quadrature_points: int,
    sensitivity_quadrature_points: int,
    expected_problem_digest: str,
    expected_implementation_sha256: str,
) -> QuadratureObjectiveEvaluator:
    """Auditor-side factory from raw app evidence, never from a batch record."""

    actual_implementation = app_free_sd_evaluator_implementation_sha256()
    if actual_implementation != expected_implementation_sha256:
        raise ValueError("objective evaluator implementation identity differs")
    primary = prepare_app_free_sd_problem(
        result,
        quadrature_points=primary_quadrature_points,
    )
    sensitivity = prepare_app_free_sd_problem(
        result,
        quadrature_points=sensitivity_quadrature_points,
    )
    if (
        app_free_sd_problem_digest(result, primary) != expected_problem_digest
        or app_free_sd_problem_digest(result, sensitivity) != expected_problem_digest
    ):
        raise ValueError("objective evaluator problem identity differs")
    evaluator = QuadratureObjectiveEvaluator(
        problem_digest=expected_problem_digest,
        primary_quadrature_points=primary.quadrature_points,
        sensitivity_quadrature_points=sensitivity.quadrature_points,
        primary_value_function=primary.value,
        sensitivity_value_function=sensitivity.value,
    )
    evaluator.validate()
    if app_free_sd_evaluator_implementation_sha256() != actual_implementation:
        raise ValueError("objective evaluator implementation changed during resolution")
    return evaluator


__all__ = [
    "AppQuadratureSensitivityBundle",
    "app_free_sd_evaluator_implementation_sha256",
    "app_free_sd_problem_digest",
    "resolve_app_free_sd_objective_evaluator",
    "run_app_free_sd_quadrature_sensitivity",
]
