"""Development fixtures: an input-only problem must reproduce native MML."""

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app
from validation.mml_free_sd_quadrature_adapter import app_free_sd_problem_digest
from validation.mml_free_sd_stationarity_adapter import (
    prepare_app_free_sd_problem,
    prepare_app_free_sd_problem_from_data,
)


def observed_rows():
    # An engineering fixture, excluded from any future qualification cohort.
    rng = np.random.default_rng(2026100401)
    return pd.DataFrame([
        (f"P{p:02d}", f"R{r}", f"C{c}", int(rng.integers(3)))
        for p in range(16) for r in range(2) for c in range(3)
        if not (p % 3 == 0 and r == 1 and c == 2)
    ], columns=["Person", "Rater", "Criterion", "Score"])


def arguments(model="RSM"):
    return dict(person_col="Person", facet_cols=["Rater", "Criterion"], score_col="Score",
                rating_min=0, rating_max=2, model=model, noncenter_facet="Criterion",
                step_facet="Criterion" if model == "PCM" else None,
                positive_facets=("Rater",), quadrature_points=9)


@pytest.mark.parametrize("model", ["RSM", "PCM"])
def test_prefit_matches_native_problem_without_invoking_fit(model, monkeypatch):
    data, args = observed_rows(), arguments(model)
    original = data.copy(deep=True)
    native_args = {k: v for k, v in args.items() if k != "quadrature_points"}
    native = app.mfrm_estimate(data, **native_args, method="MML", keep_original=True,
                               quad_points=9, estimate_population_sd=True, maxit=300)
    fitted = prepare_app_free_sd_problem(native)

    def forbidden(*a, **kw):
        pytest.fail("Pre-fit problem construction attempted optimization")

    monkeypatch.setattr(app, "mfrm_estimate", forbidden)
    monkeypatch.setattr(app, "minimize", forbidden)
    prep, problem = prepare_app_free_sd_problem_from_data(data, **args)
    assert "estimated_population_sd" not in problem.config
    assert problem.sigma_start == 1.0
    assert problem.sizes == fitted.sizes
    assert problem.structural_bounds == fitted.structural_bounds
    assert app_free_sd_problem_digest({"prep": prep}, problem) == app_free_sd_problem_digest(native, fitted)
    assert app_free_sd_problem_digest({"prep": prep}, replace(problem, quadrature_points=31)) == app_free_sd_problem_digest(native, fitted)
    point = np.r_[np.linspace(-0.4, 0.7, sum(problem.sizes.values())), np.log(1.3)]
    nll, gradient = problem.joint_value_gradient(point)
    native_nll, native_gradient = fitted.joint_value_gradient(point)
    assert nll == native_nll
    np.testing.assert_array_equal(gradient, native_gradient)
    assert problem.constraint_residual(point[:-1]) < 1e-12
    np.testing.assert_allclose(problem.likelihood_difference(point)(point + 0.01),
                               problem.joint_value_gradient(point + 0.01)[0] - nll, atol=1e-12)
    pd.testing.assert_frame_equal(original, data)
    # The evaluator owns its inputs even if the caller edits the source/table.
    data.loc[:, "Score"] = 0
    prep["data"].loc[:, "score_k"] = 0
    assert problem.joint_value_gradient(point)[0] == nll


def test_prefit_preserves_declared_zero_count_category_and_input_identity():
    data = observed_rows()
    data.loc[data["Score"] == 1, "Score"] = 2
    prep, problem = prepare_app_free_sd_problem_from_data(data, **arguments())
    assert problem.config["n_cat"] == 3
    assert prep["unused_score_categories"] == [1]
    assert set(problem.idx["score_k"]) == {0, 2}
    first = app_free_sd_problem_digest({"prep": prep}, problem)
    args = arguments()
    args["positive_facets"] = ()
    prep2, other = prepare_app_free_sd_problem_from_data(data, **args)
    assert first != app_free_sd_problem_digest({"prep": prep2}, other)


@pytest.mark.parametrize("changes", [
    {"model": "GPCM"}, {"noncenter_facet": "Person"}, {"step_facet": "Criterion"},
    {"model": "PCM", "step_facet": None}, {"positive_facets": ("Unknown",)},
    {"quadrature_points": 9.5}, {"quadrature_points": True}, {"rating_min": None},
    {"rating_max": 2.0}, {"sigma_start": 0.05}, {"sigma_start": np.nan},
])
def test_prefit_rejects_ambiguous_or_unsupported_configuration(changes):
    with pytest.raises(ValueError):
        prepare_app_free_sd_problem_from_data(observed_rows(), **(arguments() | changes))


@pytest.mark.parametrize("bad_weight", [0.0, 0.5, -1.0, np.nan])
def test_prefit_rejects_nonunit_weights_before_row_exclusion(bad_weight):
    data = observed_rows().assign(Weight=1.0)
    data.loc[0, "Weight"] = bad_weight
    with pytest.raises(ValueError, match="unit weights"):
        prepare_app_free_sd_problem_from_data(data, **arguments(), weight_col="Weight")


def test_prefit_rejects_silent_exclusions_and_unselected_weight_column():
    data = observed_rows()
    data.loc[0, "Score"] = np.nan
    with pytest.raises(ValueError, match="valid observed rows"):
        prepare_app_free_sd_problem_from_data(data, **arguments())
    with pytest.raises(ValueError, match="explicitly"):
        prepare_app_free_sd_problem_from_data(observed_rows().assign(Weight=1.0), **arguments())
    prep, problem = prepare_app_free_sd_problem_from_data(observed_rows().assign(Weight=1.0),
                                                        **arguments(), weight_col="Weight")
    assert problem.observations == len(prep["data"])
