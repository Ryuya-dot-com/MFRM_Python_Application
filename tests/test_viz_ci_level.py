"""Tests for the user-configurable visualization CI level.

The forest plot and EB-shrinkage error bars size their CI bands from
``get_viz_ci_level()`` (default 0.95). ``_ci_z(level)`` is the
two-sided normal critical value used to multiply SEs.

The contract pinned here covers:

* ``_ci_z`` matches the standard-normal quantile for the canonical
  levels (50, 66, 80, 89, 90, 95, 99) within 1e-3.
* ``_ci_z`` returns the 95 % value on NaN / out-of-range / negative
  input — must never raise.
* ``_ci_level_pct_label`` strips trailing zeros (0.95 → "95",
  0.665 → "66.5", 0.5 → "50").
* ``get_viz_ci_level`` defaults to 0.95 when the session has no value
  and clamps to the default on out-of-range values.
* ``VIZ_CI_LEVEL_OPTIONS`` covers the canonical levels (50 / 66 / 80 /
  89 / 90 / 95 / 99) and is sorted ascending — keeps the sidebar
  selectbox stable and predictable.
"""

from __future__ import annotations

import math

import pytest

import streamlit_app as app


# -----------------------------------------------------------------------------
# _ci_z — matches the standard-normal quantile within tolerance
# -----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "level, expected",
    [
        (0.50, 0.6744897501960817),
        (0.66, 0.9541652531461947),
        (0.80, 1.2815515655446004),
        (0.89, 1.5981931399228184),
        (0.90, 1.6448536269514722),
        (0.95, 1.959963984540054),
        (0.99, 2.5758293035489004),
    ],
)
def test_ci_z_matches_standard_normal_quantile(level: float, expected: float):
    z = app._ci_z(level)
    assert math.isclose(z, expected, rel_tol=0, abs_tol=1e-9), (
        f"_ci_z({level}) = {z}, expected {expected}"
    )


def test_ci_z_falls_back_on_nan_input():
    z = app._ci_z(float("nan"))
    assert math.isclose(z, 1.959963984540054, abs_tol=1e-9)


@pytest.mark.parametrize("bad_level", [0.0, 1.0, -0.1, 1.5, 99.0])
def test_ci_z_falls_back_on_out_of_range(bad_level: float):
    """Levels outside (0, 1) must fall back to the 95 % critical value
    rather than raise — a corrupted session value must not crash the
    plot pipeline."""
    z = app._ci_z(bad_level)
    assert math.isclose(z, 1.959963984540054, abs_tol=1e-9)


def test_ci_z_does_not_raise_on_string_input():
    """The fallback runs through np.isfinite, which would raise on
    strings. Confirm we coerce or fail soft — never propagate the
    TypeError to the plot."""
    try:
        z = app._ci_z("not a number")  # type: ignore[arg-type]
    except TypeError:
        # If the helper does not coerce, the calling site (the plot
        # function) sees the same TypeError and would crash. That is a
        # contract violation; pin the soft-fallback behaviour.
        pytest.fail(
            "_ci_z must fall back on non-numeric input rather than raise"
        )
    assert math.isclose(z, 1.959963984540054, abs_tol=1e-9)


# -----------------------------------------------------------------------------
# _ci_level_pct_label — strips trailing zeros
# -----------------------------------------------------------------------------


def test_ci_level_pct_label_strips_trailing_zero_on_round_values():
    assert app._ci_level_pct_label(0.95) == "95"
    assert app._ci_level_pct_label(0.50) == "50"
    assert app._ci_level_pct_label(0.66) == "66"
    assert app._ci_level_pct_label(0.89) == "89"
    assert app._ci_level_pct_label(0.99) == "99"


def test_ci_level_pct_label_preserves_fractional_part():
    """Non-integer percentages keep one decimal so labels stay
    readable (e.g. McElreath uses 0.89, but a custom 0.945 would render
    as ``94.5``, not ``94``)."""
    assert app._ci_level_pct_label(0.945) == "94.5"
    assert app._ci_level_pct_label(0.665) == "66.5"


# -----------------------------------------------------------------------------
# VIZ_CI_LEVEL_OPTIONS — stable, sorted, covers the canonical values
# -----------------------------------------------------------------------------


def test_viz_ci_level_options_cover_canonical_values():
    for lvl in (0.50, 0.66, 0.80, 0.89, 0.90, 0.95, 0.99):
        assert lvl in app.VIZ_CI_LEVEL_OPTIONS, (
            f"canonical CI level {lvl} missing from VIZ_CI_LEVEL_OPTIONS"
        )


def test_viz_ci_level_options_sorted_ascending():
    opts = list(app.VIZ_CI_LEVEL_OPTIONS)
    assert opts == sorted(opts), "VIZ_CI_LEVEL_OPTIONS must be ascending"


def test_viz_ci_level_default_is_in_options():
    assert app.VIZ_CI_LEVEL_DEFAULT in app.VIZ_CI_LEVEL_OPTIONS
    # And matches the classical 95 % convention.
    assert app.VIZ_CI_LEVEL_DEFAULT == 0.95


# -----------------------------------------------------------------------------
# get_viz_ci_level — defaults and clamping
# -----------------------------------------------------------------------------


def test_get_viz_ci_level_defaults_to_95_when_session_empty(monkeypatch):
    """Without a value in ``st.session_state["viz_ci_level"]``, the
    helper must return 0.95 so unchanged callsites keep their
    historical behaviour."""
    # st.session_state is a dict-like; we mimic an empty session.
    import streamlit as st
    monkeypatch.setitem(st.session_state, "viz_ci_level", None)
    # ``None`` round-trips to a TypeError on float(), which the helper
    # catches and falls back to the default.
    assert app.get_viz_ci_level() == app.VIZ_CI_LEVEL_DEFAULT


def test_get_viz_ci_level_returns_session_value():
    import streamlit as st
    st.session_state["viz_ci_level"] = 0.66
    try:
        assert app.get_viz_ci_level() == 0.66
    finally:
        # Restore for downstream tests.
        del st.session_state["viz_ci_level"]


def test_get_viz_ci_level_clamps_out_of_range():
    import streamlit as st
    st.session_state["viz_ci_level"] = 1.5
    try:
        assert app.get_viz_ci_level() == app.VIZ_CI_LEVEL_DEFAULT
    finally:
        del st.session_state["viz_ci_level"]
