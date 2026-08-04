"""Regression tests for the exact identified optimizer coordinates (S1/G1)."""

from __future__ import annotations

from collections import OrderedDict

import numpy as np
import pytest

import streamlit_app as app


def _config(model: str, method: str = "JMLE", n_cat: int = 4) -> dict:
    levels = ["C1", "C2", "C3"]
    facet_spec = app.build_facet_constraint(levels, centered=True)
    theta_spec = app.build_facet_constraint(["P1", "P2", "P3"], centered=False)
    return {
        "model": model,
        "method": method,
        "n_cat": n_cat,
        "n_person": 3,
        "facet_names": ["Criterion"],
        "facet_levels": {"Criterion": levels},
        "facet_specs": {"Criterion": facet_spec},
        "theta_spec": theta_spec,
        "step_facet": "Criterion" if model != "RSM" else None,
        "slope_facet": "Criterion" if model == "GPCM" else None,
        "gpcm_log_slope_bounds": (-3.0, 3.0),
    }


@pytest.mark.parametrize("full_size", [0, 1, 2, 5])
def test_sum_zero_expansion_has_exact_dimension_and_constraint(full_size):
    free = np.linspace(-0.6, 0.8, max(full_size - 1, 0))
    expanded = app.expand_sum_zero_free(free, full_size)
    jacobian = app.sum_zero_expansion_matrix(full_size)

    assert expanded.shape == (full_size,)
    assert jacobian.shape == (full_size, max(full_size - 1, 0))
    assert np.isclose(expanded.sum(), 0.0, atol=1e-14)
    assert np.allclose(expanded, jacobian @ free)


def test_sum_zero_free_roundtrip_matches_legacy_centering():
    legacy_redundant = np.array([1.7, -0.4, 0.2, 2.1])
    legacy_expanded = app.center_sum_zero(legacy_redundant)
    free = app.sum_zero_free_from_expanded(legacy_redundant)
    exact_expanded = app.expand_sum_zero_free(free, legacy_redundant.size)

    assert free.size == legacy_redundant.size - 1
    assert np.allclose(exact_expanded, legacy_expanded, atol=1e-14)


def test_sum_zero_gradient_is_exact_chain_rule():
    free = np.array([0.2, -0.4, 0.7])
    weights = np.array([1.3, -0.6, 0.2, 2.0])
    analytic = app.collapse_sum_zero_gradient(weights)
    eps = 1e-7
    numeric = np.empty_like(free)
    for j in range(free.size):
        plus = free.copy()
        minus = free.copy()
        plus[j] += eps
        minus[j] -= eps
        f_plus = weights @ app.expand_sum_zero_free(plus, 4)
        f_minus = weights @ app.expand_sum_zero_free(minus, 4)
        numeric[j] = (f_plus - f_minus) / (2.0 * eps)

    assert np.allclose(analytic, numeric, atol=1e-9)


@pytest.mark.parametrize(
    ("model", "expected_steps", "expected_slopes"),
    [
        ("RSM", 2, 0),
        ("PCM", 6, 0),
        ("GPCM", 6, 2),
    ],
)
def test_build_param_sizes_counts_only_identified_step_and_slope_coordinates(
    model, expected_steps, expected_slopes
):
    sizes = app.build_param_sizes(_config(model))

    assert sizes["steps"] == expected_steps
    assert sizes.get("log_slopes", 0) == expected_slopes


@pytest.mark.parametrize("model", ["RSM", "PCM", "GPCM"])
def test_expand_params_preserves_public_expanded_shapes(model):
    config = _config(model)
    sizes = app.build_param_sizes(config)
    par = np.linspace(-0.3, 0.4, sum(sizes.values()))
    params = app.expand_params(par, sizes, config)
    n_steps = config["n_cat"] - 1

    if model == "RSM":
        assert params["steps"].shape == (n_steps,)
        assert np.isclose(params["steps"].sum(), 0.0, atol=1e-14)
    else:
        assert params["steps_mat"].shape == (3, n_steps)
        assert np.allclose(params["steps_mat"].sum(axis=1), 0.0, atol=1e-14)
    if model == "GPCM":
        assert params["log_slopes"].shape == (3,)
        assert np.isclose(params["log_slopes"].sum(), 0.0, atol=1e-14)
        assert np.isclose(np.prod(params["slopes"]) ** (1.0 / 3.0), 1.0)


def test_gpcm_constraint_bounds_the_derived_log_slope_without_shrinking_free_box():
    config = _config("GPCM")
    sizes = app.build_param_sizes(config)
    slices = app._build_param_slices(sizes)
    constraints = app.build_optimizer_constraints(sizes, config)
    bounds = app.build_optimizer_bounds(sizes, config)

    assert len(constraints) == 1
    assert bounds[slices["log_slopes"].start] == (-3.0, 3.0)
    feasible = np.zeros(sum(sizes.values()))
    feasible[slices["log_slopes"]] = np.array([2.5, -1.0])
    derived = app.expand_sum_zero_free(feasible[slices["log_slopes"]], 3)
    assert np.all(derived >= -3.0) and np.all(derived <= 3.0)
    assert np.all(constraints[0].lb <= constraints[0].A @ feasible)
    assert np.all(constraints[0].A @ feasible <= constraints[0].ub)

    violates_derived_bound = feasible.copy()
    violates_derived_bound[slices["log_slopes"]] = np.array([3.0, 3.0])
    assert np.any(constraints[0].A @ violates_derived_bound > constraints[0].ub)


def test_split_params_covers_zero_and_nonzero_identified_blocks():
    sizes = OrderedDict([("theta", 0), ("facet", 2), ("steps", 2)])
    par = np.array([0.1, -0.2, 0.3, -0.4])
    parts = app.split_params(par, sizes)
    assert parts["theta"].size == 0
    assert np.allclose(parts["facet"], par[:2])
    assert np.allclose(parts["steps"], par[2:])


@pytest.mark.parametrize("method", ["JMLE", "MML"])
@pytest.mark.parametrize("model", ["RSM", "PCM", "GPCM"])
def test_supported_fixture_has_no_optimizer_coordinate_null_direction(model, method):
    data = app._make_self_test_gradient_data()
    result = app.mfrm_estimate(
        data,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        model=model,
        method=method,
        mml_engine="Direct",
        quad_points=5,
        maxit=40,
        reltol=1e-5,
    )
    config = result["config"]
    sizes = app.build_param_sizes(config)
    par = np.asarray(result["opt"].x, dtype=float)
    idx = app.build_indices(
        result["prep"],
        step_facet=config.get("step_facet"),
        slope_facet=config.get("slope_facet"),
    )
    if method == "JMLE":
        value_grad = lambda x: app.mfrm_loglik_jmle_value_grad(x, idx, config, sizes)
    else:
        quad = app.make_mml_quadrature(config, 5)
        value_grad = lambda x: app.mfrm_loglik_mml_value_grad(x, idx, config, sizes, quad)

    hessian = np.zeros((par.size, par.size), dtype=float)
    for j in range(par.size):
        step = 1e-5 * max(1.0, abs(float(par[j])))
        plus = par.copy()
        minus = par.copy()
        plus[j] += step
        minus[j] -= step
        _, grad_plus = value_grad(plus)
        _, grad_minus = value_grad(minus)
        hessian[:, j] = (grad_plus - grad_minus) / (2.0 * step)
    singular_values = np.linalg.svd((hessian + hessian.T) / 2.0, compute_uv=False)
    tolerance = max(hessian.shape) * np.finfo(float).eps * singular_values[0]
    rank = int(np.sum(singular_values > tolerance))
    audit = result["parameterization_audit"]

    assert par.size == sum(sizes.values())
    assert int(result["summary"]["KParams"].iloc[0]) == par.size
    assert config["parameterization"]["optimizer_dimension"] == par.size
    assert set(audit["Status"]) == {"PASS"}
    assert set(audit["ArtificialCoordinateNullDirections"]) == {0}
    if model == "GPCM" and method == "MML":
        quick_frames = app.build_result_bundle_frames(result, {})
        apa_tables = app._collect_apa_exportable_tables(result, {})
        assert "parameterization_audit" in quick_frames
        assert "Parameterization audit" in apa_tables
    assert rank == par.size, {
        "model": model,
        "method": method,
        "dimension": par.size,
        "rank": rank,
        "smallest_singular_value": singular_values[-1],
    }
