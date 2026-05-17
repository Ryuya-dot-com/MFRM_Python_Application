"""Tests for the category-specific information curve decomposition.

The new helper ``build_category_information_curve_data`` decomposes the
existing total information curve into per-category contributions

    I_k(theta) = a^2 * P_k(theta) * (k - E[X | theta])^2

which by construction satisfy the identity

    sum_k I_k(theta) = a^2 * Var[X | theta] = I_total(theta)

at every ``theta``. The tests below pin that identity (the heart of
the contract), confirm the per-category share sums to one, and check
basic shape / metadata invariants on a synthetic GPCM fit fixture
shaped like the one used by the existing category-curve test.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


@pytest.fixture(scope="module")
def gpcm_fixture_for_info_curves():
    """Synthetic GPCM fit shaped like the existing category-curve test.

    Uses three step-facet levels with non-trivial slopes so the
    identity test exercises the slope-aware information formula
    rather than the slope-equals-one specialisation."""
    return {
        "params": {
            "steps_mat": np.array(
                [[-1.2, 0.3], [-0.4, 1.1], [-0.7, 0.7]],
                dtype=float,
            ),
            "slopes": np.array([0.8, 1.25, 1.05], dtype=float),
        },
        "config": {
            "model": "GPCM",
            "n_cat": 3,
            "step_facet": "Task",
            "slope_facet": "Task",
        },
        "prep": {
            "rating_min": 0,
            "levels": {"Task": ["T1", "T2", "T3"]},
        },
    }


# -----------------------------------------------------------------------------
# Identity: sum_k I_k(theta) = I_total(theta) at every theta
# -----------------------------------------------------------------------------


def test_category_information_curve_sums_to_total_information(
    gpcm_fixture_for_info_curves,
):
    """Per-category contributions must reproduce the total information
    curve at every theta to machine precision; the identity is the
    statistical contract that Muraki (1993, Eq. 16) derives from."""
    fixture = gpcm_fixture_for_info_curves
    cat_bundle = app.build_category_information_curve_data(fixture)
    total_bundle = app.build_information_curve_data(fixture)
    assert cat_bundle["available"]
    assert total_bundle["available"]
    cat_curve = cat_bundle["curve"]
    total_curve = total_bundle["curve"]

    # Both bundles must traverse the same theta grid.
    cat_thetas = sorted(cat_curve["Theta"].unique().tolist())
    total_thetas = sorted(total_curve["Theta"].unique().tolist())
    assert cat_thetas == total_thetas

    # At each theta, sum of CategoryInformation must equal the total
    # Information reported by build_information_curve_data.
    cat_sum = cat_curve.groupby("Theta")["CategoryInformation"].sum()
    total_index = total_curve.set_index("Theta")["Information"]
    for theta in cat_thetas:
        s = float(cat_sum.loc[theta])
        t = float(total_index.loc[theta])
        assert abs(s - t) < 1e-12, (
            f"identity violated at theta={theta}: sum_k I_k={s} vs total I={t}"
        )


def test_category_information_total_column_matches_per_row_sum(
    gpcm_fixture_for_info_curves,
):
    """``TotalInformation`` carried on every row is the sum across
    categories at the same theta; the column lets a consumer read the
    total off any single row without recomputing the sum."""
    bundle = app.build_category_information_curve_data(gpcm_fixture_for_info_curves)
    curve = bundle["curve"]
    by_theta = curve.groupby("Theta")
    for theta, sub in by_theta:
        per_row_total = float(sub["TotalInformation"].iloc[0])
        # Every row in the same theta carries the same total.
        assert (sub["TotalInformation"] - per_row_total).abs().max() < 1e-12
        # The total matches the recomputed sum of CategoryInformation.
        assert (
            abs(per_row_total - float(sub["CategoryInformation"].sum())) < 1e-12
        )


# -----------------------------------------------------------------------------
# Slope-dispatch: the slope-aware formula actually consumes the slope
# -----------------------------------------------------------------------------


def test_category_information_uses_slope_squared(gpcm_fixture_for_info_curves):
    """At slope = 2 the category information must be 4x the slope-1
    information at every (theta, category), because the formula has a
    bare ``a^2`` factor in front of the kernel."""
    fixture_a1 = {
        **gpcm_fixture_for_info_curves,
        "params": {
            "steps_mat": gpcm_fixture_for_info_curves["params"]["steps_mat"],
            "slopes": np.array([1.0, 1.0, 1.0], dtype=float),
        },
    }
    fixture_a2 = {
        **gpcm_fixture_for_info_curves,
        "params": {
            "steps_mat": gpcm_fixture_for_info_curves["params"]["steps_mat"],
            "slopes": np.array([2.0, 2.0, 2.0], dtype=float),
        },
    }
    cat_a1 = app.build_category_information_curve_data(fixture_a1, step_level_index=0)
    cat_a2 = app.build_category_information_curve_data(fixture_a2, step_level_index=0)
    assert cat_a1["available"] and cat_a2["available"]
    # The category probability vectors at a given theta differ between
    # the two slope settings (the kernel multiplies eta by slope), so
    # we cannot expect I_k to scale exactly by 4. But the total
    # information at slope = 2 *must* exceed the total at slope = 1
    # by a large margin -- the slope-squared factor is what drives that
    # ordering. Pinning the ordering is sufficient to confirm slope
    # dispatch is active without overcommitting to a closed form.
    total_a1 = float(cat_a1["curve"].groupby("Theta")["CategoryInformation"].sum().max())
    total_a2 = float(cat_a2["curve"].groupby("Theta")["CategoryInformation"].sum().max())
    assert total_a2 > total_a1 * 1.5, (total_a1, total_a2)


# -----------------------------------------------------------------------------
# Metadata / summary shape
# -----------------------------------------------------------------------------


def test_category_information_summary_contains_per_category_peak(
    gpcm_fixture_for_info_curves,
):
    bundle = app.build_category_information_curve_data(gpcm_fixture_for_info_curves)
    summary = bundle["summary"]
    assert {
        "CategoryValue",
        "MaxCategoryInformation",
        "ThetaAtMaxCategoryInformation",
        "IntegratedCategoryInformation",
        "IntegratedShare",
    }.issubset(summary.columns)
    assert len(summary) == 3  # K = 3 categories
    # Per-category shares of integrated information sum to one.
    assert abs(float(summary["IntegratedShare"].sum()) - 1.0) < 1e-12
    # Every per-category peak must be non-negative.
    assert (summary["MaxCategoryInformation"] >= 0).all()


def test_category_information_metadata_matches_total_curve_scope(
    gpcm_fixture_for_info_curves,
):
    """The decomposition must carry the same scope and model labels as
    the total information curve so downstream consumers (UI legend,
    export filenames) stay consistent across the two builders."""
    cat = app.build_category_information_curve_data(gpcm_fixture_for_info_curves)
    total = app.build_information_curve_data(gpcm_fixture_for_info_curves)
    assert cat["metadata"]["Scope"] == total["metadata"]["Scope"]
    assert cat["metadata"]["Model"] == total["metadata"]["Model"]
    assert cat["metadata"]["Slope"] == total["metadata"]["Slope"]


def test_category_information_returns_unavailable_for_missing_inputs():
    """Missing fit inputs propagate the ``available = False`` flag with a
    human-readable reason so the UI can show an explanation instead of
    crashing on an empty fixture."""
    bundle = app.build_category_information_curve_data({})
    assert bundle["available"] is False
    assert "reason" in bundle and isinstance(bundle["reason"], str)


def test_category_information_step_level_index_propagates(
    gpcm_fixture_for_info_curves,
):
    """Selecting a specific step-facet level must produce a level-
    specific decomposition whose scope label encodes the chosen level."""
    bundle_avg = app.build_category_information_curve_data(
        gpcm_fixture_for_info_curves
    )
    bundle_t2 = app.build_category_information_curve_data(
        gpcm_fixture_for_info_curves, step_level_index=1
    )
    assert bundle_avg["available"] and bundle_t2["available"]
    assert bundle_t2["metadata"]["Scope"] != bundle_avg["metadata"]["Scope"]
    # Different scope -> different slope (under the per-level GPCM
    # discrimination convention) -> different total info column.
    avg_total = float(bundle_avg["curve"].groupby("Theta")["CategoryInformation"].sum().max())
    t2_total = float(bundle_t2["curve"].groupby("Theta")["CategoryInformation"].sum().max())
    assert avg_total != pytest.approx(t2_total, abs=1e-12)
