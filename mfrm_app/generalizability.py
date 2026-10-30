"""UI-independent balanced crossed G/D calculation and qualified public output.

The observed-score entry does not fit an MFRM. Existing design restrictions,
variance handling, provenance and the interval hold apply in every frontend.
"""
from __future__ import annotations

from copy import deepcopy
import itertools
import json
import math

import numpy as np
import pandas as pd

from mfrm_app.engine import dataframe_fingerprint, stable_json_fingerprint


_GD_STUDY_DEFAULT_GRID = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
_GD_CALCULATION_VERSION = "balanced_anova_v2"
_GD_CI_REASON = "Person-only resampling has not qualified intervals for the declared random-facet G/D targets."
_RESERVED_GD_NAMES = {"Residual", "G", "Phi", "IsObservedDesign", "_section"}

GD_SECTION_ORDER = (
    "observed_coefficients",
    "selected_d_study",
    "variance_components",
    "d_study",
    "metadata",
)
GD_SECTION_NOTES = {
    "observed_coefficients": (
        "Observed-design point estimates on original scores: G for relative decisions, "
        "Phi for absolute decisions, and error variances / SEMs. Intervals for G and Phi are withheld."
    ),
    "selected_d_study": (
        "Projected G, Phi, error variances and SEMs at the selected D-Study counts "
        "assuming the same universe and fully crossed observations."
    ),
    "variance_components": (
        "Balanced random-effects ANOVA components. RawVariance retains signed method-of-moments estimates; "
        "Variance clips negative values at zero for G/Phi. Highest-order interaction is confounded with Residual."
    ),
    "d_study": (
        "Closed-form D-Study point projections across planned facet counts "
        "(IsObservedDesign=True marks the observed counts)."
    ),
    "metadata": (
        "Audit record of design checks, score basis, declared population and calculation settings "
        "(interval_inference_allowed=False)."
    ),
}

GD_VARIANCE_NEGATIVE_NOTE = (
    "Balanced random-effects ANOVA components. One or more RawVariance estimates are negative "
    "(possible sampling variability near zero or design misspecification) and were clipped at zero for G/Phi. "
    "Highest-order interaction is confounded with Residual."
)

_GD_INPUT_REASONS = {
    "invalid_result": "Result is not a fit dictionary.",
    "missing_data": "Fit does not carry a raw data frame for G-study.",
    "unsupported_facets": "One or two distinct non-object random facets are required for a G-study.",
    "invalid_columns": "Required columns are missing, duplicated, or use reserved component names.",
    "invalid_scores": "Every retained score must be finite; no rows were silently removed.",
    "invalid_labels": "Facet labels must be nonmissing, nonblank and unambiguous after text conversion.",
    "insufficient_levels": "Each facet, including the object of measurement, needs at least two levels.",
    "duplicate_cells": "Repeated cells are unsupported; define the observation unit before analysis.",
    "missing_cells": "The full crossed design is incomplete; missing ratings are not zero scores.",
    "unsupported_weights": "This G-study supports unit observation weights only.",
    "invalid_grid": "The maximum planned facet count must be a positive integer.",
    "invalid_components": "The variance decomposition did not produce finite components.",
    "undefined_coefficients": "G and Phi are undefined because usable object variance is absent.",
}


def _planned_facet_grid(n_max_per_facet: int) -> tuple[int, ...]:
    """Build a contiguous positive-integer D-Study grid up to n_max_per_facet.

    For standard planning ranges (<= 25), every integer 1..n_max_per_facet is
    included so custom counts (e.g., 5, 7, 9, 12) never fail grid lookup.
    For very large upper bounds, keep 1..15 exact and step cleanly up to n_max.
    """
    n_max = int(n_max_per_facet)
    if n_max <= 25:
        return tuple(range(1, n_max + 1))
    dense = list(range(1, 16))
    step = max(2, math.ceil((n_max - 15) / 10))
    tail = list(range(15 + step, n_max, step))
    if not tail or tail[-1] != n_max:
        tail.append(n_max)
    return tuple(dense + tail)


def validate_fixed_condition(fixed_condition, object_facet, score_col, random_facets):
    """Validate explicit fixed-facet conditioning settings for 3+ facet designs."""
    if fixed_condition is None:
        return None
    if not isinstance(fixed_condition, dict) or set(fixed_condition) != {"facet", "mode", "level"}:
        raise ValueError("Fixed-facet conditioning requires facet, mode and level fields")
    facet = fixed_condition.get("facet")
    mode = fixed_condition.get("mode")
    level = fixed_condition.get("level")
    if (not isinstance(facet, str) or not facet.strip() or ":" in facet or facet in _RESERVED_GD_NAMES
            or facet in {object_facet, score_col, *random_facets}):
        raise ValueError("Select a distinct fixed-facet column")
    if mode not in ("level", "sum", "mean"):
        raise ValueError("Fixed-facet conditioning mode must be level, sum or mean")
    if mode == "level":
        if not isinstance(level, str) or not level.strip():
            raise ValueError("Choose a level for fixed-facet slice conditioning")
    elif level is not None:
        raise ValueError("Composite fixed-facet modes aggregate across all levels without a single level")
    return dict(facet=facet, mode=mode, level=level)


def prepare_gd_input_data(data, settings):
    """Apply explicit fixed-facet conditioning only when elementary design rules hold."""
    fc = settings.get("fixed_condition") if isinstance(settings, dict) else None
    if fc is None:
        return data, "", None
    if not isinstance(data, pd.DataFrame) or data.empty:
        return data, "missing_data", None
    object_facet = settings["object_facet"]
    random_facets = list(settings["random_facets"])
    score_col = settings["score_col"]
    fixed_facet = fc["facet"]
    mode = fc["mode"]
    level = fc["level"]
    if len(random_facets) not in (1, 2) or len({object_facet, *random_facets}) != 1 + len(random_facets):
        return data, "unsupported_facets", None
    elementary_facets = [object_facet, *random_facets, fixed_facet]
    if (not data.columns.is_unique or score_col in elementary_facets
            or any(not f or ":" in f or f in _RESERVED_GD_NAMES for f in elementary_facets)
            or any(c not in data.columns for c in [*elementary_facets, score_col])):
        return data, "invalid_columns", None
    if "Weight" in data.columns:
        weights = pd.to_numeric(data["Weight"], errors="coerce")
        if not (np.isfinite(weights.to_numpy(dtype=float)).all() and weights.eq(1).all()):
            return data, "unsupported_weights", None
    f_labels = data[fixed_facet]
    if f_labels.isna().any() or f_labels.astype(str).str.strip().eq("").any():
        return data, "invalid_labels", None
    if f_labels.nunique() != f_labels.astype(str).nunique():
        return data, "invalid_labels", None
    f_str = f_labels.astype(str)
    if mode == "level":
        sliced = data.loc[f_str.eq(level)].copy()
        if sliced.empty:
            return data, "insufficient_levels", None
        return sliced, "", dict(facet=fixed_facet, mode="level", level=level, n_fixed_levels=1)

    scores = pd.to_numeric(data[score_col], errors="coerce")
    if not np.isfinite(scores.to_numpy(dtype=float)).all():
        return data, "invalid_scores", None
    work = data[elementary_facets].copy()
    for facet in elementary_facets:
        labels = work[facet]
        if labels.isna().any() or labels.astype(str).str.strip().eq("").any():
            return data, "invalid_labels", None
        if labels.nunique() != labels.astype(str).nunique():
            return data, "invalid_labels", None
        work[facet] = labels.astype(str)
    work[score_col] = scores.to_numpy(dtype=float)
    levels = {f: int(work[f].nunique()) for f in elementary_facets}
    if min(levels.values()) < 2:
        return data, "insufficient_levels", None
    if work.duplicated(elementary_facets).any():
        return data, "duplicate_cells", None
    if not _is_balanced_one_obs_per_cell(work, elementary_facets):
        return data, "missing_cells", None
    grouped = work.groupby([object_facet, *random_facets], sort=True, as_index=False)[score_col]
    agg = grouped.sum() if mode == "sum" else grouped.mean()
    if "Weight" in data.columns:
        agg["Weight"] = 1.0
    return agg, "", dict(facet=fixed_facet, mode=mode, level=None, n_fixed_levels=levels[fixed_facet])


def _prepare_g_d_data(data, object_facet, random_facets, score_col):
    """Validate the same observation table used by decomposition and reporting."""
    design = {"n_input_rows": len(data) if isinstance(data, pd.DataFrame) else 0,
              "n_observations": 0, "weight_basis": "unit_observation_weights"}
    if not isinstance(data, pd.DataFrame) or data.empty:
        return None, design, "missing_data"
    facets = [object_facet] + list(random_facets)
    if (len(random_facets) not in (1, 2)
            or any(not isinstance(f, str) for f in facets)
            or len(set(facets)) != len(facets)):
        return None, design, "unsupported_facets"
    if (not data.columns.is_unique or score_col in facets
            or any(not f or ":" in f or f in _RESERVED_GD_NAMES
                   for f in facets)
            or any(c not in data.columns for c in facets + [score_col])):
        return None, design, "invalid_columns"
    design["data_fingerprint"] = dataframe_fingerprint(data, length=64)
    work = data[facets + [score_col]].copy()
    work[score_col] = pd.to_numeric(work[score_col], errors="coerce")
    if not np.isfinite(work[score_col].to_numpy(dtype=float)).all():
        return None, design, "invalid_scores"
    if "Weight" in data.columns:
        weights = pd.to_numeric(data["Weight"], errors="coerce")
        if not (np.isfinite(weights.to_numpy(dtype=float)).all() and weights.eq(1).all()):
            return None, design, "unsupported_weights"
    for facet in facets:
        labels = work[facet]
        if labels.isna().any() or labels.astype(str).str.strip().eq("").any():
            return None, design, "invalid_labels"
        if labels.nunique() != labels.astype(str).nunique():
            return None, design, "invalid_labels"
        work[facet] = labels.astype(str)
    levels = {f: int(work[f].nunique()) for f in facets}
    design.update(object_facet=object_facet, random_facets=list(random_facets),
                  observed_n={f: levels[f] for f in random_facets}, levels=levels,
                  score_column=score_col, score_basis="retained_preprocessed_scores")
    design["n_expected_cells"] = math.prod(levels.values())
    if min(levels.values()) < 2:
        return None, design, "insufficient_levels"
    if work.duplicated(facets).any():
        return None, design, "duplicate_cells"
    if not _is_balanced_one_obs_per_cell(work, facets):
        return None, design, "missing_cells"
    design.update(n_observations=len(work), grand_mean=float(work[score_col].mean()))
    return work, design, ""


def _is_balanced_one_obs_per_cell(
    work: pd.DataFrame,
    facets: list[str],
) -> bool:
    """Detect whether the design has exactly one observation per cell."""
    if work.empty or not facets:
        return False
    if any(f not in work.columns for f in facets) or work[facets].isna().any().any():
        return False
    expected = 1
    for facet in facets:
        expected *= int(work[facet].nunique())
    cells = work.groupby(facets, observed=True).size()
    return bool(len(work) == expected and len(cells) == expected and (cells == 1).all())


def _crossed_anova_two_way_variance_components(
    work: pd.DataFrame,
    *,
    object_facet: str,
    random_facets: list[str],
    score_col: str,
    grand_mean: float,
    ss_total: float,
) -> pd.DataFrame | None:
    """Full random-effects 3-way ANOVA decomposition under one-obs-per-cell.

    Uses sequential orthogonal tensor projection (Brennan, 2001, Table 3.5)
    rather than subtracting large sums of squares, eliminating IEEE-754
    catastrophic cancellation when main-effect variance dominates residual error.
    """
    del grand_mean, ss_total  # Retained in signature for internal compatibility; centered tensor is used.
    if len(random_facets) != 2:
        return None
    p_facet = object_facet
    i_facet, j_facet = random_facets[0], random_facets[1]
    facets_3way = [p_facet, i_facet, j_facet]
    if not _is_balanced_one_obs_per_cell(work, facets_3way):
        return None

    p_cat = pd.Categorical(work[p_facet])
    i_cat = pd.Categorical(work[i_facet])
    j_cat = pd.Categorical(work[j_facet])
    n_p, n_i, n_j = len(p_cat.categories), len(i_cat.categories), len(j_cat.categories)
    if min(n_p, n_i, n_j) < 2:
        return None

    y_tensor = np.empty((n_p, n_i, n_j), dtype=np.float64)
    y_tensor[p_cat.codes, i_cat.codes, j_cat.codes] = work[score_col].to_numpy(dtype=np.float64)

    # Two-pass grand-mean centering for exact shift invariance.
    y0 = y_tensor - float(y_tensor.mean())
    y0 = y0 - float(y0.mean())

    # Sweep mutually orthogonal 1-way main effects in order of decreasing magnitude
    # so large main effects never pollute smaller main effects, interactions, or residuals.
    main_axes = {
        0: (1, 2),
        1: (0, 2),
        2: (0, 1),
    }
    init_ranges = {
        ax: float(np.ptp(y0.mean(axis=reduce_axes)))
        for ax, reduce_axes in main_axes.items()
    }
    y_swept = y0.copy()
    main_effects: dict[int, np.ndarray] = {
        0: np.zeros((n_p, 1, 1), dtype=np.float64),
        1: np.zeros((1, n_i, 1), dtype=np.float64),
        2: np.zeros((1, 1, n_j), dtype=np.float64),
    }
    sorted_main_axes = sorted(main_axes, key=lambda ax: -init_ranges[ax])
    for _ in range(2):
        for ax in sorted_main_axes:
            delta = y_swept.mean(axis=main_axes[ax], keepdims=True)
            main_effects[ax] = main_effects[ax] + delta
            y_swept = y_swept - delta

    eff_p, eff_i, eff_j = main_effects[0], main_effects[1], main_effects[2]
    ss_p = float(n_i * n_j * np.sum(eff_p ** 2))
    ss_i = float(n_p * n_j * np.sum(eff_i ** 2))
    ss_j = float(n_p * n_i * np.sum(eff_j ** 2))

    # Sweep mutually orthogonal 2-way interaction subspaces in order of decreasing magnitude.
    two_way_axes = {
        "pi": 2,
        "pj": 1,
        "ij": 0,
    }
    init_tw_ranges = {
        key: float(np.ptp(y_swept.mean(axis=reduce_ax)))
        for key, reduce_ax in two_way_axes.items()
    }
    tw_effects = {
        "pi": np.zeros((n_p, n_i, 1), dtype=np.float64),
        "pj": np.zeros((n_p, 1, n_j), dtype=np.float64),
        "ij": np.zeros((1, n_i, n_j), dtype=np.float64),
    }
    sorted_tw = sorted(two_way_axes, key=lambda k: -init_tw_ranges[k])
    for _ in range(2):
        for key in sorted_tw:
            delta = y_swept.mean(axis=two_way_axes[key], keepdims=True)
            tw_effects[key] = tw_effects[key] + delta
            y_swept = y_swept - delta

    ss_pi = float(n_j * np.sum(tw_effects["pi"] ** 2))
    ss_pj = float(n_i * np.sum(tw_effects["pj"] ** 2))
    ss_ij = float(n_p * np.sum(tw_effects["ij"] ** 2))
    ss_residual = float(np.sum(y_swept ** 2))

    # Degrees of freedom for each source (Brennan 2001 Table 3.5).
    df_p = n_p - 1
    df_i = n_i - 1
    df_j = n_j - 1
    df_pi = df_p * df_i
    df_pj = df_p * df_j
    df_ij = df_i * df_j
    df_residual = df_p * df_i * df_j

    def safe_div(num, den):
        return float(num / den) if (np.isfinite(num) and den > 0) else np.nan

    ms_p = safe_div(ss_p, df_p)
    ms_i = safe_div(ss_i, df_i)
    ms_j = safe_div(ss_j, df_j)
    ms_pi = safe_div(ss_pi, df_pi)
    ms_pj = safe_div(ss_pj, df_pj)
    ms_ij = safe_div(ss_ij, df_ij)
    ms_residual = safe_div(ss_residual, df_residual)

    sigma2_residual = max(ms_residual, 0.0)
    sigma2_pi = max((ms_pi - ms_residual) / max(n_j, 1), 0.0)
    sigma2_pj = max((ms_pj - ms_residual) / max(n_i, 1), 0.0)
    sigma2_ij = max((ms_ij - ms_residual) / max(n_p, 1), 0.0)
    sigma2_p = max(
        (ms_p - ms_pi - ms_pj + ms_residual) / max(n_i * n_j, 1), 0.0
    )
    sigma2_i = max(
        (ms_i - ms_pi - ms_ij + ms_residual) / max(n_p * n_j, 1), 0.0
    )
    sigma2_j = max(
        (ms_j - ms_pj - ms_ij + ms_residual) / max(n_p * n_i, 1), 0.0
    )

    rows = [
        {"Source": p_facet, "SumSquares": ss_p, "df": df_p,
         "MeanSquares": ms_p, "Variance": sigma2_p,
         "MeanObsPerLevel": float(n_i * n_j), "Levels": n_p,
         "Term": "main"},
        {"Source": i_facet, "SumSquares": ss_i, "df": df_i,
         "MeanSquares": ms_i, "Variance": sigma2_i,
         "MeanObsPerLevel": float(n_p * n_j), "Levels": n_i,
         "Term": "main"},
        {"Source": j_facet, "SumSquares": ss_j, "df": df_j,
         "MeanSquares": ms_j, "Variance": sigma2_j,
         "MeanObsPerLevel": float(n_p * n_i), "Levels": n_j,
         "Term": "main"},
        {"Source": f"{p_facet}:{i_facet}", "SumSquares": ss_pi, "df": df_pi,
         "MeanSquares": ms_pi, "Variance": sigma2_pi,
         "MeanObsPerLevel": float(n_j), "Levels": n_p * n_i,
         "Term": "two-way"},
        {"Source": f"{p_facet}:{j_facet}", "SumSquares": ss_pj, "df": df_pj,
         "MeanSquares": ms_pj, "Variance": sigma2_pj,
         "MeanObsPerLevel": float(n_i), "Levels": n_p * n_j,
         "Term": "two-way"},
        {"Source": f"{i_facet}:{j_facet}", "SumSquares": ss_ij, "df": df_ij,
         "MeanSquares": ms_ij, "Variance": sigma2_ij,
         "MeanObsPerLevel": float(n_p), "Levels": n_i * n_j,
         "Term": "two-way"},
        {"Source": "Residual", "SumSquares": ss_residual, "df": df_residual,
         "MeanSquares": ms_residual, "Variance": sigma2_residual,
         "MeanObsPerLevel": float("nan"),
         "Levels": None, "Term": "residual"},
    ]
    out = pd.DataFrame(rows)
    out["RawVariance"] = [
        (ms_p - ms_pi - ms_pj + ms_residual) / (n_i * n_j),
        (ms_i - ms_pi - ms_ij + ms_residual) / (n_p * n_j),
        (ms_j - ms_pj - ms_ij + ms_residual) / (n_p * n_i),
        (ms_pi - ms_residual) / n_j,
        (ms_pj - ms_residual) / n_i,
        (ms_ij - ms_residual) / n_p,
        ms_residual,
    ]
    total_var = float(np.nansum(out["Variance"]))
    out["ProportionVariance"] = (
        out["Variance"] / total_var if total_var > 0 else np.nan
    )
    return out[[
        "Source", "Term", "RawVariance", "Variance", "ProportionVariance",
        "SumSquares", "df", "MeanSquares", "MeanObsPerLevel", "Levels",
    ]]


def _crossed_anova_variance_components(
    obs_df: pd.DataFrame,
    *,
    object_facet: str,
    random_facets: list[str],
    score_col: str = "Score",
) -> pd.DataFrame:
    """Balanced crossed ANOVA for one or two non-object random facets.

    Exactly one observation per cell and at least two levels per factor
    are required. Unsupported designs return an empty table. RawVariance
    retains signed method-of-moments estimates; Variance clips at zero.
    Two random facets use the full three-way interaction decomposition;
    one random facet uses orthogonal two-way projection below.
    """
    work, _, reason = _prepare_g_d_data(obs_df, object_facet, random_facets, score_col)
    if reason:
        return pd.DataFrame()

    grand_mean = float(work[score_col].mean())
    ss_total = float(((work[score_col] - grand_mean) ** 2).sum())

    two_way = _crossed_anova_two_way_variance_components(
        work,
        object_facet=object_facet,
        random_facets=list(random_facets),
        score_col=score_col,
        grand_mean=grand_mean,
        ss_total=ss_total,
    )
    if two_way is not None and not two_way.empty:
        return two_way

    # Balanced one-random-facet (2-way) orthogonal tensor sweep.
    i_facet = random_facets[0]
    p_cat = pd.Categorical(work[object_facet])
    i_cat = pd.Categorical(work[i_facet])
    n_p, n_i = len(p_cat.categories), len(i_cat.categories)
    if min(n_p, n_i) < 2:
        return pd.DataFrame()

    y_mat = np.empty((n_p, n_i), dtype=np.float64)
    y_mat[p_cat.codes, i_cat.codes] = work[score_col].to_numpy(dtype=np.float64)
    y0 = y_mat - float(y_mat.mean())
    y0 = y0 - float(y0.mean())

    # Sweep larger main effect first to prevent cancellation.
    eff_p = np.zeros((n_p, 1), dtype=np.float64)
    eff_i = np.zeros((1, n_i), dtype=np.float64)
    y_swept = y0.copy()
    p_first = float(np.ptp(y0.mean(axis=1))) >= float(np.ptp(y0.mean(axis=0)))
    for _ in range(2):
        if p_first:
            dp = y_swept.mean(axis=1, keepdims=True)
            eff_p = eff_p + dp
            y_swept = y_swept - dp
            di = y_swept.mean(axis=0, keepdims=True)
            eff_i = eff_i + di
            y_swept = y_swept - di
        else:
            di = y_swept.mean(axis=0, keepdims=True)
            eff_i = eff_i + di
            y_swept = y_swept - di
            dp = y_swept.mean(axis=1, keepdims=True)
            eff_p = eff_p + dp
            y_swept = y_swept - dp

    ss_p = float(n_i * np.sum(eff_p ** 2))
    ss_i = float(n_p * np.sum(eff_i ** 2))
    ss_residual = float(np.sum(y_swept ** 2))

    df_p = n_p - 1
    df_i = n_i - 1
    df_residual = df_p * df_i

    ms_p = ss_p / df_p if df_p > 0 else np.nan
    ms_i = ss_i / df_i if df_i > 0 else np.nan
    ms_residual = ss_residual / df_residual if df_residual > 0 else np.nan

    rows = [
        {
            "Source": object_facet,
            "Term": "main",
            "SumSquares": ss_p,
            "df": df_p,
            "MeanSquares": ms_p,
            "MeanObsPerLevel": float(n_i),
            "Levels": n_p,
        },
        {
            "Source": i_facet,
            "Term": "main",
            "SumSquares": ss_i,
            "df": df_i,
            "MeanSquares": ms_i,
            "MeanObsPerLevel": float(n_p),
            "Levels": n_i,
        },
        {
            "Source": "Residual",
            "Term": "residual",
            "SumSquares": ss_residual,
            "df": df_residual,
            "MeanSquares": ms_residual,
            "MeanObsPerLevel": float("nan"),
            "Levels": None,
        },
    ]
    out = pd.DataFrame(rows)
    sigma2_e = ms_residual if np.isfinite(ms_residual) else np.nan
    out["RawVariance"] = [
        (ms_p - sigma2_e) / n_i if np.isfinite(ms_p) and np.isfinite(sigma2_e) else np.nan,
        (ms_i - sigma2_e) / n_p if np.isfinite(ms_i) and np.isfinite(sigma2_e) else np.nan,
        float(sigma2_e) if np.isfinite(sigma2_e) else np.nan,
    ]
    out["Variance"] = out["RawVariance"].clip(lower=0)
    total_var = float(np.nansum(out["Variance"]))
    out["ProportionVariance"] = (
        out["Variance"] / total_var if total_var > 0 else np.nan
    )
    return out[
        ["Source", "Term", "RawVariance", "Variance", "ProportionVariance", "SumSquares", "df",
         "MeanSquares", "MeanObsPerLevel", "Levels"]
    ]


def _g_phi_from_variance_components(
    var_components: pd.DataFrame,
    *,
    object_facet: str,
    facet_observations: dict[str, int] | None = None,
) -> dict:
    """Compute G / Phi from a variance-component table (Brennan, 2001, Eq. 3.18 / 3.19)."""
    if not isinstance(var_components, pd.DataFrame) or var_components.empty:
        return {"G": np.nan, "Phi": np.nan, "details": {}}

    sources = [str(s) for s in var_components["Source"].tolist()]
    variances = var_components["Variance"].to_numpy(dtype=float)
    lookup = {src: float(var) for src, var in zip(sources, variances)}
    sigma2_p = lookup.get(object_facet, np.nan)
    sigma2_e = lookup.get("Residual", np.nan)
    if not (np.isfinite(sigma2_p) and np.isfinite(sigma2_e) and sigma2_p > 0):
        return {"G": np.nan, "Phi": np.nan, "details": {"sigma2_p": sigma2_p, "sigma2_e": sigma2_e}}

    facet_obs = facet_observations or {}
    two_way_sources = (
        [src for src, term in zip(sources, var_components["Term"].tolist()) if term == "two-way"]
        if "Term" in var_components.columns else []
    )

    if two_way_sources:
        random_facets: list[str] = []
        for src in two_way_sources:
            for part in src.split(":"):
                if part != object_facet and part not in random_facets:
                    random_facets.append(part)
        n_a = max(int(facet_obs.get(random_facets[0], 1)), 1) if random_facets else 1
        n_b = max(int(facet_obs.get(random_facets[1], 1)), 1) if len(random_facets) > 1 else 1

        def lookup_two_way(pair: tuple[str, str]) -> float:
            target = set(pair)
            for src in two_way_sources:
                if set(src.split(":")) == target:
                    return float(lookup[src])
            return 0.0

        i_facet = random_facets[0] if random_facets else None
        j_facet = random_facets[1] if len(random_facets) > 1 else None
        sigma2_pi = lookup_two_way((object_facet, i_facet)) if i_facet else 0.0
        sigma2_pj = lookup_two_way((object_facet, j_facet)) if j_facet else 0.0
        sigma2_ij = lookup_two_way((i_facet, j_facet)) if i_facet and j_facet else 0.0
        sigma2_i = lookup.get(i_facet, 0.0) if i_facet else 0.0
        sigma2_j = lookup.get(j_facet, 0.0) if j_facet else 0.0

        sigma2_delta = (
            sigma2_pi / n_a + sigma2_pj / n_b + sigma2_e / (n_a * n_b)
        )
        sigma2_big_delta = (
            sigma2_delta
            + sigma2_i / n_a + sigma2_j / n_b + sigma2_ij / (n_a * n_b)
        )
        g_coef = sigma2_p / (sigma2_p + sigma2_delta)
        phi_coef = sigma2_p / (sigma2_p + sigma2_big_delta)
        return {
            "G": float(g_coef),
            "Phi": float(phi_coef),
            "details": {
                "sigma2_object": sigma2_p,
                "sigma2_residual": sigma2_e,
                "sigma2_relative_error": sigma2_delta,
                "sigma2_absolute_error": sigma2_big_delta,
                "facet_observations": dict(facet_obs),
                "decomposition": "full_3way_brennan_eq_3_18",
                "random_facets": list(random_facets),
            },
        }

    sigma2_main_per_obs = 0.0
    for source, variance in lookup.items():
        if source in {object_facet, "Residual"} or not np.isfinite(variance):
            continue
        n_per = facet_obs.get(source, 1)
        sigma2_main_per_obs += variance / max(int(n_per), 1)

    n_total = max(int(np.prod([max(facet_obs.get(s, 1), 1)
                               for s in lookup
                               if s not in {object_facet, "Residual"}])), 1)
    relative_err = sigma2_e / n_total
    g_coef = sigma2_p / (sigma2_p + relative_err)
    phi_coef = sigma2_p / (sigma2_p + sigma2_main_per_obs + relative_err)
    return {
        "G": float(g_coef),
        "Phi": float(phi_coef),
        "details": {
            "sigma2_object": sigma2_p,
            "sigma2_residual": sigma2_e,
            "sigma2_relative_error": relative_err,
            "sigma2_absolute_error": sigma2_main_per_obs + relative_err,
            "sigma2_main_per_obs": sigma2_main_per_obs,
            "n_total": n_total,
            "facet_observations": dict(facet_obs),
            "decomposition": "main_effects_only_approximation",
        },
    }


def _build_d_study_grid(
    var_components: pd.DataFrame,
    *,
    object_facet: str,
    random_facets: list[str],
    facet_grids: dict[str, tuple[int, ...]],
) -> pd.DataFrame:
    """Forecast G / Phi over the cartesian product of facet sample sizes."""
    if not isinstance(var_components, pd.DataFrame) or var_components.empty:
        return pd.DataFrame()
    if not random_facets:
        return pd.DataFrame()
    grids = [facet_grids.get(f, _GD_STUDY_DEFAULT_GRID) for f in random_facets]
    rows = []
    for combo in itertools.product(*grids):
        facet_obs = {f: int(n) for f, n in zip(random_facets, combo)}
        gp = _g_phi_from_variance_components(
            var_components, object_facet=object_facet, facet_observations=facet_obs,
        )
        row = {f: int(n) for f, n in zip(random_facets, combo)}
        row["G"] = float(gp["G"]) if np.isfinite(gp.get("G", np.nan)) else np.nan
        row["Phi"] = float(gp["Phi"]) if np.isfinite(gp.get("Phi", np.nan)) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _bootstrap_g_phi_ci(
    data: pd.DataFrame,
    *,
    object_facet: str,
    random_facets: list[str],
    score_col: str,
    facet_observations: dict[str, int],
    n_bootstrap: int = 200,
    confidence: float = 0.95,
    seed: int | None = 20260601,
    extra_designs: list[dict[str, int]] | None = None,
    _variance_fn=None,
    _coef_fn=None,
) -> dict:
    """Technical person-resampling calculation; not qualified for inference.

    Percentiles summarize finite replicates only. Each coefficient/design
    retains its success count and failure reasons against all attempts.
    This calculation does not establish coverage for random-facet G/D targets.
    """
    variance_fn = _variance_fn or _crossed_anova_variance_components
    coef_fn = _coef_fn or _g_phi_from_variance_components
    n_bootstrap = max(int(n_bootstrap), 1)
    confidence = float(min(max(confidence, 0.0), 1.0))
    base = {
        "available": False, "inference_allowed": False,
        "requested": True,
        "inference_reason": _GD_CI_REASON,
        "G_lower": np.nan, "G_upper": np.nan, "Phi_lower": np.nan, "Phi_upper": np.nan,
        "G_replicates": np.array([]), "Phi_replicates": np.array([]),
        "confidence": confidence, "n_bootstrap": n_bootstrap, "seed": seed,
        "n_attempted": 0, "n_success": 0, "n_success_phi": 0,
        "failure_counts": {"G": {}, "Phi": {}},
        "method": "cluster_bootstrap_on_persons", "extra_design_replicates": [],
    }
    data, _, reason = _prepare_g_d_data(data, object_facet, random_facets, score_col)
    if reason:
        return {**base, "reason": _GD_INPUT_REASONS[reason]}
    rng = np.random.default_rng(int(seed) if seed is not None else None)
    persons = data[object_facet].unique()
    rows_by_person = {p: data[data[object_facet] == p] for p in persons}
    designs = [facet_observations] + list(extra_designs or [])
    values = [{"G": [], "Phi": []} for _ in designs]
    failures = [{"G": {}, "Phi": {}} for _ in designs]

    def fail(counts, coefficient, why):
        counts[coefficient][why] = counts[coefficient].get(why, 0) + 1

    for _ in range(n_bootstrap):
        parts = []
        for k, person in enumerate(rng.choice(persons, size=len(persons), replace=True)):
            sub = rows_by_person[person].copy()
            sub[object_facet] = f"BOOT_{k:04d}_{person}"
            parts.append(sub)
        boot_df = pd.concat(parts, ignore_index=True)
        why = ""
        try:
            components = variance_fn(
                boot_df, object_facet=object_facet, random_facets=random_facets,
                score_col=score_col)
            if components.empty or not np.isfinite(components["Variance"]).all():
                why = "invalid_components"
        except Exception as exc:
            why = "decomposition_exception:" + type(exc).__name__
        if why:
            for counts in failures:
                for coefficient in ("G", "Phi"):
                    fail(counts, coefficient, why)
            continue
        for design, draws, counts in zip(designs, values, failures):
            try:
                coef = coef_fn(
                    components, object_facet=object_facet, facet_observations=design)
            except Exception as exc:
                for coefficient in ("G", "Phi"):
                    fail(counts, coefficient, "coefficient_exception:" + type(exc).__name__)
                continue
            for coefficient in ("G", "Phi"):
                value = float(coef.get(coefficient, np.nan))
                if np.isfinite(value):
                    draws[coefficient].append(value)
                else:
                    fail(counts, coefficient, "nonfinite_coefficient")

    def summary(draws, counts):
        out = {"n_attempted": n_bootstrap, "failure_counts": counts}
        for coefficient in ("G", "Phi"):
            arr = np.asarray(draws[coefficient], dtype=float)
            bounds = (np.quantile(arr, [(1-confidence)/2, (1+confidence)/2])
                      if len(arr) >= 2 else [np.nan, np.nan])
            out[coefficient + "_lower"], out[coefficient + "_upper"] = map(float, bounds)
            out[coefficient + "_replicates"] = arr
        out["n_success"] = len(draws["G"])
        out["n_success_phi"] = len(draws["Phi"])
        return out

    base.update(summary(values[0], failures[0]))
    base["available"] = base["n_success"] >= 2 or base["n_success_phi"] >= 2
    base["reason"] = ("" if base["available"] else
                      "Insufficient bootstrap replicates produced finite coefficients.")
    base["extra_design_replicates"] = [
        {"design": dict(design), **summary(draws, counts)}
        for design, draws, counts in zip(designs[1:], values[1:], failures[1:])
    ]
    return base


def compute_generalizability_study(
    res: dict,
    *,
    object_facet: str = "Person",
    random_facets: list[str] | None = None,
    n_max_per_facet: int = 10,
    score_col: str = "Score",
    fixed_condition: dict | None = None,
    bootstrap_ci: bool = False,
    n_bootstrap: int = 200,
    bootstrap_confidence: float = 0.95,
    bootstrap_seed: int | None = 20260601,
    _bootstrap_fn=None,
) -> dict:
    """Use retained MFRM scores; share the standalone G/D calculation and limits."""
    valid = isinstance(res, dict)
    config = res.get("config") if valid and isinstance(res.get("config"), dict) else {}
    prep = res.get("prep", {}) if valid else {}
    if random_facets is None:
        fixed_name = fixed_condition.get("facet") if isinstance(fixed_condition, dict) else None
        random_facets = [
            f for f in (config or {}).get("facet_names", [])
            if f != object_facet and f != fixed_name
        ]
    bundle = _compute_generalizability_study(
        prep.get("data") if isinstance(prep, dict) else None,
        object_facet=object_facet, random_facets=random_facets, score_col=score_col,
        n_max_per_facet=n_max_per_facet, score_basis="retained_preprocessed_scores",
        fixed_condition=fixed_condition,
        bootstrap_ci=bootstrap_ci, n_bootstrap=n_bootstrap,
        bootstrap_confidence=bootstrap_confidence, bootstrap_seed=bootstrap_seed,
        _bootstrap_fn=_bootstrap_fn)
    if not valid:
        bundle.update(reason_code="invalid_result", reason=_GD_INPUT_REASONS["invalid_result"])
    return bundle


def compute_generalizability_from_data(
    data: pd.DataFrame, *, object_facet: str, random_facets: list[str],
    score_col: str = "Score", n_max_per_facet: int = 10,
    target_population: str = "", fixed_condition: dict | None = None,
) -> dict:
    """Analyze observed scores directly, without MFRM preparation or estimation.

    Supports 1-2 random facets, with optional verified fixed-facet conditioning
    (single-level slice or complete-balanced sum/mean composite) for 3+ facet tables.
    """
    return _compute_generalizability_study(
        data, object_facet=object_facet, random_facets=random_facets, score_col=score_col,
        n_max_per_facet=n_max_per_facet, score_basis="original_observed_scores",
        target_population=target_population, fixed_condition=fixed_condition)


def _compute_generalizability_study(
    data: pd.DataFrame,
    *,
    object_facet: str = "Person",
    random_facets: list[str] | None = None,
    n_max_per_facet: int = 10,
    score_col: str = "Score",
    score_basis: str,
    target_population: str = "",
    fixed_condition: dict | None = None,
    bootstrap_ci: bool = False,
    n_bootstrap: int = 200,
    bootstrap_confidence: float = 0.95,
    bootstrap_seed: int | None = 20260601,
    _bootstrap_fn=None,
) -> dict:
    """Shared balanced ANOVA, design validation and qualification for both entries."""
    design: dict = {}

    def unavailable(code, components=None):
        return {
            "available": False, "reason_code": code, "reason": _GD_INPUT_REASONS[code],
            "calculation_version": _GD_CALCULATION_VERSION,
            "variance_components": components if components is not None else pd.DataFrame(),
            "observed_coefficients": {"G": np.nan, "Phi": np.nan},
            "d_study": pd.DataFrame(), "design": design,
            "caveat": "", "reference": "", "bootstrap_ci": {
                "available": False, "inference_allowed": False,
                "requested": bool(bootstrap_ci),
                "inference_reason": _GD_CI_REASON, "n_attempted": 0},
        }

    if not isinstance(random_facets, (list, tuple)):
        return unavailable("unsupported_facets")

    validated_fc = None
    fc_meta = None
    raw_fp = ""
    working_input = data
    if fixed_condition is not None:
        validated_fc = validate_fixed_condition(fixed_condition, object_facet, score_col, list(random_facets))
        raw_fp = dataframe_fingerprint(data, length=64) if isinstance(data, pd.DataFrame) else ""
        working_input, prep_reason, fc_meta = prepare_gd_input_data(
            data,
            dict(
                object_facet=object_facet,
                random_facets=list(random_facets),
                score_col=score_col,
                fixed_condition=validated_fc,
            ),
        )
        if prep_reason:
            design = dict(
                n_input_rows=len(data) if isinstance(data, pd.DataFrame) else 0,
                n_observations=0, weight_basis="unit_observation_weights",
                score_basis=score_basis, target_population=target_population,
                fixed_facets=[validated_fc["facet"]],
                fixed_condition=deepcopy(validated_fc),
                raw_data_fingerprint=raw_fp, data_fingerprint=raw_fp,
                generalization_basis="fixed_facet_conditioned_observed_scores",
            )
            return unavailable(prep_reason)

    work_df, design, reason = _prepare_g_d_data(working_input, object_facet, random_facets, score_col)
    design.update(score_basis=score_basis, target_population=target_population,
                  fixed_facets=[], generalization_basis="all_non_object_facets_random")
    if fc_meta is not None and validated_fc is not None:
        design.update(
            fixed_facets=[fc_meta["facet"]],
            fixed_condition=deepcopy(validated_fc),
            fixed_levels_aggregated=int(fc_meta["n_fixed_levels"]),
            raw_data_fingerprint=raw_fp,
            generalization_basis=(
                "single_fixed_level_slice" if fc_meta["mode"] == "level"
                else f"balanced_fixed_facet_{fc_meta['mode']}"
            ),
        )
    if reason:
        return unavailable(reason)
    if (isinstance(n_max_per_facet, (bool, np.bool_))
            or not isinstance(n_max_per_facet, (int, np.integer)) or n_max_per_facet < 1):
        return unavailable("invalid_grid")
    design["max_planned_count"] = int(n_max_per_facet)
    analysis_id = "gstudy_" + stable_json_fingerprint(
        {"version": _GD_CALCULATION_VERSION, "design": design}, length=24)
    design["analysis_id"] = analysis_id
    decomposition = ("full_3way_brennan_eq_3_18" if len(random_facets) == 2
                     else "balanced_2way_anova")
    design["decomposition"] = decomposition
    with np.errstate(over="ignore", invalid="ignore"):
        var_components = _crossed_anova_variance_components(
            work_df, object_facet=object_facet, random_facets=random_facets, score_col=score_col)
    if var_components.empty or not np.isfinite(
        var_components[["RawVariance", "Variance", "SumSquares", "MeanSquares"]].to_numpy(dtype=float)
    ).all():
        return unavailable("invalid_components")
    observed_n = design["observed_n"]
    observed_coefficients = _g_phi_from_variance_components(
        var_components, object_facet=object_facet, facet_observations=observed_n)
    if not all(np.isfinite(observed_coefficients[k]) for k in ("G", "Phi")):
        return unavailable("undefined_coefficients", var_components)
    planned_points = _planned_facet_grid(int(n_max_per_facet))
    facet_grids = {f: planned_points for f in random_facets}
    d_study = _build_d_study_grid(
        var_components, object_facet=object_facet,
        random_facets=random_facets, facet_grids=facet_grids)
    observed_mask = np.logical_and.reduce([d_study[f].eq(observed_n[f]) for f in random_facets])
    d_study["IsObservedDesign"] = observed_mask
    if not observed_mask.any():
        d_study = pd.concat([d_study, pd.DataFrame([{
            **observed_n, "G": observed_coefficients["G"], "Phi": observed_coefficients["Phi"],
            "IsObservedDesign": True,
        }])], ignore_index=True)
    observed_coefficients["details"]["decomposition"] = decomposition
    caveat = (
        "Balanced fully crossed random-effects ANOVA with one observation per cell; "
        "the highest-order interaction is confounded with residual error. "
        "RawVariance retains signed method-of-moments estimates; Variance clips them at zero "
        "for the reported G/Phi calculation. "
        + ("Scores are original observed scores, " if score_basis == "original_observed_scores"
           else "Scores are the retained preprocessed scores, ")
        + "with unit observation weights. G concerns relative decisions; Phi concerns absolute "
        "decisions. Neither coefficient alone establishes suitability for a use or a cut score. "
        "Intervals remain unqualified for the random-facet targets."
    )
    if fc_meta is not None:
        caveat += (
            f" Conditioned on fixed facet {fc_meta['facet']} = {fc_meta['level']}."
            if fc_meta["mode"] == "level"
            else (
                f" Aggregated ({fc_meta['mode']}) across {fc_meta['n_fixed_levels']} fixed "
                f"{fc_meta['facet']} levels after verifying complete balanced elementary cells."
            )
        )
    reference = "Cronbach, Gleser, Nanda, & Rajaratnam (1972); Brennan (2001)."

    bootstrap = {
        "available": False, "reason": "Bootstrap CI was not requested.",
        "requested": False,
        "G_lower": np.nan, "G_upper": np.nan,
        "Phi_lower": np.nan, "Phi_upper": np.nan,
        "confidence": float(bootstrap_confidence),
        "n_bootstrap": int(n_bootstrap),
        "n_success": 0, "n_success_phi": 0, "n_attempted": 0,
        "extra_design_replicates": [],
    }
    if bootstrap_ci:
        extra_designs: list[dict[str, int]] = []
        if isinstance(d_study, pd.DataFrame) and not d_study.empty:
            for _, row in d_study.iterrows():
                design_dict: dict[str, int] = {}
                for facet in random_facets:
                    if facet in row.index and pd.notna(row[facet]):
                        design_dict[facet] = int(row[facet])
                if design_dict:
                    extra_designs.append(design_dict)
        bootstrap_runner = _bootstrap_fn or _bootstrap_g_phi_ci
        bootstrap = bootstrap_runner(
            work_df,
            object_facet=object_facet,
            random_facets=list(random_facets),
            score_col=score_col,
            facet_observations=observed_n,
            n_bootstrap=int(n_bootstrap),
            confidence=float(bootstrap_confidence),
            seed=bootstrap_seed,
            extra_designs=extra_designs,
        )
        if (
            isinstance(d_study, pd.DataFrame)
            and not d_study.empty
            and bootstrap.get("available")
        ):
            extra_replicates = bootstrap.get("extra_design_replicates", [])
            ci_rows: dict[tuple, dict] = {}
            for entry in extra_replicates:
                key = tuple(int(entry["design"].get(f, 0)) for f in random_facets)
                ci_rows[key] = entry
            g_lows, g_highs, phi_lows, phi_highs = [], [], [], []
            for _, row in d_study.iterrows():
                key = tuple(int(row[f]) if f in row.index else 0 for f in random_facets)
                entry = ci_rows.get(key)
                if entry is None:
                    g_lows.append(np.nan)
                    g_highs.append(np.nan)
                    phi_lows.append(np.nan)
                    phi_highs.append(np.nan)
                else:
                    g_lows.append(float(entry["G_lower"]))
                    g_highs.append(float(entry["G_upper"]))
                    phi_lows.append(float(entry["Phi_lower"]))
                    phi_highs.append(float(entry["Phi_upper"]))
            d_study = d_study.copy()
            d_study["G_lower"] = g_lows
            d_study["G_upper"] = g_highs
            d_study["Phi_lower"] = phi_lows
            d_study["Phi_upper"] = phi_highs

    bootstrap.update(inference_allowed=False, inference_reason=_GD_CI_REASON)
    return {
        "available": True,
        "reason": "", "reason_code": "",
        "calculation_version": _GD_CALCULATION_VERSION,
        "variance_components": var_components,
        "observed_coefficients": observed_coefficients,
        "d_study": d_study,
        "design": design,
        "caveat": caveat,
        "reference": reference,
        "bootstrap_ci": bootstrap,
    }


def _g_d_error_row(coefficients: dict) -> dict:
    details = coefficients["details"]
    rel_var = float(details["sigma2_relative_error"])
    abs_var = float(details["sigma2_absolute_error"])
    return {
        "G": coefficients["G"],
        "Phi": coefficients["Phi"],
        "RelativeErrorVariance": rel_var,
        "AbsoluteErrorVariance": abs_var,
        "RelativeSEM": float(np.sqrt(max(rel_var, 0.0))) if np.isfinite(rel_var) else np.nan,
        "AbsoluteSEM": float(np.sqrt(max(abs_var, 0.0))) if np.isfinite(abs_var) else np.nan,
    }


def g_d_study_export_frame(bundle: dict) -> pd.DataFrame:
    """Self-describing public CSV; technical bootstrap bounds are held."""
    design = bundle.get("design", {})
    bootstrap = bundle.get("bootstrap_ci", {})
    metadata = {
        "available": bool(bundle.get("available")),
        "reason_code": bundle.get("reason_code", ""), "reason": bundle.get("reason", ""),
        "calculation_version": bundle.get("calculation_version", "unrecorded"),
        **design,
        "supported_design": "Fully crossed balanced; one or two non-object random facets; at least two levels per factor; one observation per cell; unit weights",
        "residual_basis": "Highest-order interaction confounded with residual error",
        "variance_policy": "RawVariance signed; Variance clipped at zero",
        "coefficient_targets": "G: relative decisions; Phi: absolute decisions",
        "error_variance_units": "score_squared",
        "sem_units": "score_units",
        "population_status": "declaration_not_verified",
        "interval_inference_allowed": False, "interval_hold_reason": _GD_CI_REASON,
        "bootstrap_requested": bool(bootstrap.get("requested", False)),
        "bootstrap_method": bootstrap.get("method", "not_run"),
        "bootstrap_n_requested": bootstrap.get("n_bootstrap", 0) if bootstrap.get("requested") else 0,
        "bootstrap_n_attempted": bootstrap.get("n_attempted", 0),
        "bootstrap_n_success_G": bootstrap.get("n_success", 0),
        "bootstrap_n_success_Phi": bootstrap.get("n_success_phi", 0),
        "bootstrap_failure_counts": bootstrap.get("failure_counts", {}),
        "bootstrap_seed": bootstrap.get("seed"),
        "bootstrap_confidence": bootstrap.get("confidence") if bootstrap.get("requested") else None,
        "scope": bundle.get("caveat", ""), "reference": bundle.get("reference", ""),
    }
    rows = [{"_section": "metadata", "Field": key,
             "Value": json.dumps(value, ensure_ascii=False, sort_keys=True)
             if isinstance(value, (dict, list)) else str(value)}
            for key, value in metadata.items()]
    frames = [pd.DataFrame(rows)]
    components = bundle.get("variance_components")
    if isinstance(components, pd.DataFrame) and not components.empty:
        frames.append(components.assign(_section="variance_components"))
    if bundle.get("available"):
        coef = bundle["observed_coefficients"]
        frames.append(pd.DataFrame([{
            "_section": "observed_coefficients", **design.get("observed_n", {}),
            **_g_d_error_row(coef), "IsObservedDesign": True}]))
        d_study = bundle.get("d_study")
        if isinstance(d_study, pd.DataFrame) and not d_study.empty:
            columns = list(design.get("random_facets", [])) + ["G", "Phi", "IsObservedDesign"]
            frames.append(d_study.reindex(columns=columns).assign(_section="d_study"))
    if bundle.get("available") and bundle.get("selected_plan"):
        frames.append(pd.DataFrame([{**bundle["selected_plan"], "_section": "selected_d_study"}]))
    return pd.concat(frames, ignore_index=True)
