"""Prospectively registered native exact-CMLE hard-anchor contracts."""

from __future__ import annotations

from itertools import product

import numpy as np
import pandas as pd
import pytest
from scipy.special import logsumexp

from mfrm_app.cmle import (
    CMLEEligibilityError,
    cmle_objective_value_grad,
    cmle_objective_value_grad_hessian,
    fit_cmle,
    prepare_cmle_design,
)
from mfrm_app.cmle_person_scoring import score_cmle_persons_wle
from mfrm_app.cmle_wle_bootstrap import _refit_arguments, _row_category_kernels
from mfrm_app.cmle_wle_fit import compute_cmle_wle_person_fit
from mfrm_app.cmle_wle_uncertainty import cmle_wle_calibration_sensitivity
from mfrm_app.person_scoring import score_fixed_calibration_persons


def _frame() -> pd.DataFrame:
    scores = {
        "P1": [0, 1, 1, 2],
        "P2": [1, 0, 2, 1],
        "P3": [2, 1, 1, 0],
        "P4": [1, 2, 0, 1],
        "P5": [0, 2, 1, 1],
        "P6": [2, 0, 1, 1],
    }
    units = [("R1", "C1"), ("R1", "C2"), ("R2", "C1"), ("R2", "C2")]
    return pd.DataFrame(
        [
            (person, rater, criterion, score)
            for person, values in scores.items()
            for (rater, criterion), score in zip(units, values)
        ],
        columns=["Person", "Rater", "Criterion", "Score"],
    )


def _arguments(**overrides) -> dict[str, object]:
    arguments: dict[str, object] = {
        "person_col": "Person",
        "facet_cols": ["Rater", "Criterion"],
        "score_col": "Score",
        "rating_min": 0,
        "rating_max": 2,
        "model": "RSM",
    }
    arguments.update(overrides)
    return arguments


def _rater_anchor(value: float = -0.35) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "ParameterType": "Facet",
                "Facet": "Rater",
                "Level": "R1",
                "Value": value,
            }
        ]
    )


def _manual_anchored_nll(parameters: np.ndarray, design) -> float:
    log_likelihood = 0.0
    for pattern in design.patterns:
        log_kernel = pattern.offset + np.einsum(
            "ukp,p->uk", pattern.design, parameters
        )
        log_likelihood += float(
            pattern.observed_offset + pattern.observed_statistics @ parameters
        )
        for score in np.flatnonzero(pattern.score_frequencies > 0):
            values = []
            for response in product(
                range(design.n_categories), repeat=len(pattern.signature)
            ):
                if sum(response) == int(score):
                    values.append(
                        sum(
                            log_kernel[unit, category]
                            for unit, category in enumerate(response)
                        )
                    )
            log_likelihood -= float(pattern.score_frequencies[score]) * float(
                logsumexp(values)
            )
    return -log_likelihood


def _finite_difference(function, point: np.ndarray, step: float = 1e-6) -> np.ndarray:
    gradient = np.zeros_like(point, dtype=float)
    for index in range(len(point)):
        plus = point.copy()
        minus = point.copy()
        plus[index] += step
        minus[index] -= step
        gradient[index] = (function(plus) - function(minus)) / (2.0 * step)
    return gradient


def test_no_anchor_and_empty_anchor_inputs_preserve_legacy_coordinates_exactly():
    none = prepare_cmle_design(_frame(), **_arguments())
    empty = prepare_cmle_design(_frame(), **_arguments(hard_anchors=[]))
    parameters = np.array([0.31, -0.22, 0.47])

    assert none.parameter_names == [
        "facet:Rater:free:R1",
        "facet:Criterion:free:C1",
        "step:__shared__:free:1",
    ]
    assert none.parameter_names == empty.parameter_names
    assert np.array_equal(none.row_design, empty.row_design)
    assert np.array_equal(none.row_offset, np.zeros_like(none.row_offset))
    assert np.array_equal(none.row_offset, empty.row_offset)
    for left, right in zip(none.patterns, empty.patterns):
        assert np.array_equal(left.design, right.design)
        assert np.array_equal(left.offset, right.offset)
        assert left.observed_offset == right.observed_offset == 0.0
    for function in (cmle_objective_value_grad, cmle_objective_value_grad_hessian):
        left = function(parameters, none)
        right = function(parameters, empty)
        assert all(np.array_equal(np.asarray(a), np.asarray(b)) for a, b in zip(left, right))


def test_anchored_objective_matches_explicit_affine_enumeration_and_gradient():
    design = prepare_cmle_design(
        _frame(), **_arguments(hard_anchors=_rater_anchor(-0.35))
    )
    parameters = np.array([0.18, -0.22, 0.47])
    value, gradient = cmle_objective_value_grad(parameters, design)
    manual = _manual_anchored_nll(parameters, design)
    numeric = _finite_difference(
        lambda point: cmle_objective_value_grad(point, design)[0], parameters
    )

    assert value == pytest.approx(manual, abs=1e-10)
    assert np.max(np.abs(gradient - numeric)) < 1e-5
    assert np.any(design.row_offset != 0.0)
    assert bool(design.audit["eligible"])
    assert int(design.audit["conditional_nullity"]) == 0


def test_anchor_is_exact_in_outputs_and_anchor_row_order_is_invariant():
    rows = [
        {"ParameterType": "facet", "Facet": "Criterion", "Level": "C2", "Value": 0.15},
        {"ParameterType": "Facet", "Facet": "Rater", "Level": "R1", "Value": -0.35},
    ]
    first = fit_cmle(
        _frame(), **_arguments(hard_anchors=pd.DataFrame(rows)), gtol=1e-8
    )
    second = fit_cmle(
        _frame(), **_arguments(hard_anchors=pd.DataFrame(rows[::-1])), gtol=1e-8
    )
    first_facets = first["facets"]["others"].set_index(["Facet", "Level"])
    second_facets = second["facets"]["others"].set_index(["Facet", "Level"])

    assert first["design"].hard_anchors[["Facet", "Level"]].values.tolist() == [
        ["Rater", "R1"],
        ["Criterion", "C2"],
    ]
    assert first["design"].parameter_names == second["design"].parameter_names
    assert np.max(
        np.abs(first["coefficients"]["Estimate"] - second["coefficients"]["Estimate"])
    ) < 1e-12
    for key, expected in ((('Rater', 'R1'), -0.35), (('Criterion', 'C2'), 0.15)):
        row = first_facets.loc[key]
        assert float(row["Estimate"]) == expected
        assert float(row["SE"]) == 0.0
        assert bool(row["Anchored"])
        assert row["Constraint"] == "hard_anchor"
        assert float(second_facets.loc[key, "Estimate"]) == expected
    assert int(first["summary"].iloc[0]["HardAnchors"]) == 2
    assert int(first["summary"].iloc[0]["KParams"]) == 3


@pytest.mark.parametrize(
    ("anchors", "message"),
    [
        ([{"ParameterType": "Step", "Facet": "Rater", "Level": "R1", "Value": 0}], "Facet.*only"),
        ([{"ParameterType": "Facet", "Facet": "Rater", "Level": "R1", "Step": 1, "Value": 0}], "must not specify Step"),
        ([{"ParameterType": "Facet", "Facet": "Ghost", "Level": "R1", "Value": 0}], "unknown facets"),
        ([{"ParameterType": "Facet", "Facet": "Rater", "Level": "Ghost", "Value": 0}], "unknown levels"),
        ([{"ParameterType": "Facet", "Facet": "Rater", "Level": "R1", "Value": np.inf}], "finite numeric"),
        ([
            {"ParameterType": "Facet", "Facet": "Rater", "Level": "R1", "Value": 0},
            {"ParameterType": "Facet", "Facet": "Rater", "Level": "R1", "Value": 0},
        ], "duplicate"),
        ([{"ParameterType": "Facet", "Facet": "Rater", "Level": "R1", "Value": 0, "Typo": 1}], "unexpected columns"),
    ],
)
def test_invalid_anchor_contracts_fail_closed(anchors, message):
    with pytest.raises(ValueError, match=message):
        prepare_cmle_design(_frame(), **_arguments(hard_anchors=anchors))


def test_all_fixed_binary_design_fails_closed_when_no_free_parameter_remains():
    frame = _frame().assign(Score=lambda value: value["Score"].clip(upper=1))
    anchors = [
        {"ParameterType": "Facet", "Facet": facet, "Level": level, "Value": 0.0}
        for facet, levels in (("Rater", ["R1", "R2"]), ("Criterion", ["C1", "C2"]))
        for level in levels
    ]
    design = prepare_cmle_design(
        frame,
        **_arguments(rating_max=1, hard_anchors=anchors),
    )
    assert design.n_parameters == 0
    assert not bool(design.audit["eligible"])
    assert "no_free_structural_parameters" in set(
        design.audit["issues_table"]["Code"]
    )
    with pytest.raises(CMLEEligibilityError, match="no free structural parameters"):
        fit_cmle(frame, **_arguments(rating_max=1, hard_anchors=anchors))


def test_fixed_offsets_propagate_to_surfaces_wle_fit_bootstrap_and_draws():
    fit = fit_cmle(
        _frame(), **_arguments(hard_anchors=_rater_anchor(-0.35)), gtol=1e-8
    )
    design = fit["design"]
    parameters = fit["coefficients"].set_index("Parameter")["Estimate"].reindex(
        design.parameter_names
    ).to_numpy(dtype=float)
    intercepts = design.row_offset + np.einsum(
        "rkp,p->rk", design.row_design, parameters
    )
    _, bootstrap_kernels = _row_category_kernels(fit)
    assert np.max(np.abs(bootstrap_kernels - intercepts)) < 1e-12
    assert _refit_arguments(fit)["hard_anchors"].equals(design.hard_anchors)

    surfaces = fit["surfaces"]
    first_rows = design.data.reset_index().drop_duplicates(list(design.facet_cols))
    reconstructed = []
    for row in first_rows.itertuples(index=False):
        for category in range(design.n_categories):
            reconstructed.append(
                {
                    "Rater": row.Rater,
                    "Criterion": row.Criterion,
                    "InternalCategory": category,
                    "ExpectedKernel": intercepts[int(row.index), category],
                }
            )
    surface_check = surfaces.merge(
        pd.DataFrame(reconstructed),
        on=["Rater", "Criterion", "InternalCategory"],
        validate="one_to_one",
    )
    assert np.max(
        np.abs(surface_check["LogKernelWithoutPerson"] - surface_check["ExpectedKernel"])
    ) < 1e-12

    bridge = score_cmle_persons_wle(fit)
    direct = score_fixed_calibration_persons(
        design.data[design.person_col].astype(str),
        design.data[design.score_col].to_numpy(dtype=int) - design.rating_min,
        intercepts,
        np.tile(np.arange(design.n_categories), (len(design.data), 1)),
        person_levels=fit["person_status"]["Person"].astype(str).tolist(),
    )
    assert np.max(np.abs(bridge["Estimate"] - direct["Estimate"])) < 1e-12

    fit_result = compute_cmle_wle_person_fit(fit)
    expected_by_row = fit_result["observations"]["ExpectedInternalCategory"].to_numpy()
    assert np.isfinite(expected_by_row).all()
    sensitivity = cmle_wle_calibration_sensitivity(
        fit, n_draws=8, seed=20260810, covariance_scale=0.0
    )
    assert np.max(sensitivity["persons"]["CalibrationDrawMeanMinusBaseline"].abs()) < 1e-12

