"""Identification contracts for hard and group anchors in JMLE/MML facets."""

import numpy as np

import streamlit_app as app


def test_one_hard_anchor_identifies_origin_without_forcing_only_free_level_to_zero():
    spec = app.build_facet_constraint(
        ["R01", "R02"],
        anchors={"R01": 0.25},
        centered=True,
    )

    assert spec["requested_centered"] is True
    assert spec["absolute_origin_identified"] is True
    assert spec["origin_constraint"] == "absolute_anchor"
    assert spec["centered"] is False
    assert spec["n_params"] == 1
    assert np.allclose(
        app.expand_facet_with_constraints(np.array([0.80]), spec),
        [0.25, 0.80],
    )


def test_two_hard_anchors_leave_both_remaining_levels_estimable():
    spec = app.build_facet_constraint(
        ["R01", "R02", "R03", "R04"],
        anchors={"R01": -0.20, "R02": 0.10},
        centered=True,
    )
    free = np.array([0.35, -0.05])

    assert spec["n_params"] == 2
    assert np.allclose(
        app.expand_facet_with_constraints(free, spec),
        [-0.20, 0.10, 0.35, -0.05],
    )
    assert np.allclose(
        app.collapse_facet_gradient(np.array([9.0, 8.0, 0.7, -0.2]), spec),
        [0.7, -0.2],
    )


def test_all_hard_anchored_levels_have_no_free_coordinates():
    spec = app.build_facet_constraint(
        ["R01", "R02"],
        anchors={"R01": -0.20, "R02": 0.10},
        centered=True,
    )

    assert spec["n_params"] == 0
    assert np.allclose(
        app.expand_facet_with_constraints(np.array([]), spec),
        [-0.20, 0.10],
    )


def test_no_anchor_retains_existing_sum_to_zero_identification():
    spec = app.build_facet_constraint(
        ["R01", "R02", "R03"],
        centered=True,
    )
    expanded = app.expand_facet_with_constraints(np.array([0.20, -0.50]), spec)

    assert spec["absolute_origin_identified"] is False
    assert spec["origin_constraint"] == "sum_to_zero"
    assert spec["centered"] is True
    assert spec["n_params"] == 2
    assert np.allclose(expanded, [0.20, -0.50, 0.30])
    assert np.isclose(expanded.sum(), 0.0)


def test_absolute_group_anchor_frees_ungrouped_levels_from_extra_centering():
    spec = app.build_facet_constraint(
        ["R01", "R02", "R03", "R04"],
        groups={"R01": "A", "R02": "A"},
        group_values={"A": 0.20},
        centered=True,
    )
    # First coordinate identifies R01 within group A; R02 is derived so their
    # mean is 0.20. R03/R04 are both direct coordinates because group A fixes
    # the facet origin.
    expanded = app.expand_facet_with_constraints(
        np.array([0.10, 0.60, -0.40]),
        spec,
    )

    assert spec["absolute_origin_identified"] is True
    assert spec["origin_constraint"] == "absolute_group_anchor"
    assert spec["centered"] is False
    assert spec["n_params"] == 3
    assert np.allclose(expanded, [0.10, 0.30, 0.60, -0.40])
    assert np.isclose(expanded[:2].mean(), 0.20)
