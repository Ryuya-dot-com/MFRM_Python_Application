"""UI-independent balanced crossed G/D calculation and qualified public output.

The observed-score entry does not fit an MFRM. Existing design restrictions,
variance handling, provenance and the interval hold apply in every frontend.
"""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd

from mfrm_app.engine import dataframe_fingerprint, stable_json_fingerprint


_GD_STUDY_DEFAULT_GRID = (1, 2, 3, 4, 5, 6, 8, 10)
_GD_CALCULATION_VERSION = "balanced_anova_v2"
_GD_CI_REASON = "Person-only resampling has not qualified intervals for the declared random-facet G/D targets."
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
            or any(not f or ":" in f or f in {"Residual", "G", "Phi", "IsObservedDesign", "_section"}
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

    For a balanced p x i x j design with one observation per cell, all
    factors random, and Person as the object of measurement, the
    method-of-moments equations from Brennan (2001, Table 3.5) yield:

      sigma2_p   = (MS_p   - MS_pi - MS_pj + MS_pij) / (n_i n_j)
      sigma2_i   = (MS_i   - MS_pi - MS_ij + MS_pij) / (n_p n_j)
      sigma2_j   = (MS_j   - MS_pj - MS_ij + MS_pij) / (n_p n_i)
      sigma2_pi  = (MS_pi  - MS_pij) / n_j
      sigma2_pj  = (MS_pj  - MS_pij) / n_i
      sigma2_ij  = (MS_ij  - MS_pij) / n_p
      sigma2_pij = MS_pij              (residual; confounded with error)

    Returns the variance-component table with seven rows (three main
    effects + three two-way interactions + residual) when applicable,
    or ``None`` when the design does not satisfy the one-observation-
    per-cell balanced assumption. RawVariance retains signed estimates;
    Variance is clipped at zero for the coefficient calculation.
    """
    if len(random_facets) != 2:
        return None
    p_facet = object_facet
    i_facet, j_facet = random_facets[0], random_facets[1]
    facets_3way = [p_facet, i_facet, j_facet]
    if not _is_balanced_one_obs_per_cell(work, facets_3way):
        return None

    n_p = int(work[p_facet].nunique())
    n_i = int(work[i_facet].nunique())
    n_j = int(work[j_facet].nunique())
    if min(n_p, n_i, n_j) < 2:
        return None

    def cell_mean(group_cols):
        return work.groupby(group_cols, observed=False)[score_col].mean()

    # Main-effect SS (each weighted by the cell count of the
    # marginalised facets so the formula matches the balanced ANOVA).
    means_p = cell_mean(p_facet)
    ss_p = float(n_i * n_j * np.sum((means_p.to_numpy(dtype=float) - grand_mean) ** 2))
    means_i = cell_mean(i_facet)
    ss_i = float(n_p * n_j * np.sum((means_i.to_numpy(dtype=float) - grand_mean) ** 2))
    means_j = cell_mean(j_facet)
    ss_j = float(n_p * n_i * np.sum((means_j.to_numpy(dtype=float) - grand_mean) ** 2))

    # Two-way SS: SS_pi = n_j * sum (bar_Y_pi. - bar_Y_p.. - bar_Y_.i. + bar_Y_...)^2
    means_pi = cell_mean([p_facet, i_facet]).unstack().to_numpy()
    means_pj = cell_mean([p_facet, j_facet]).unstack().to_numpy()
    means_ij = cell_mean([i_facet, j_facet]).unstack().to_numpy()
    means_p_arr = means_p.to_numpy(dtype=float).reshape(-1, 1)
    means_i_arr = means_i.to_numpy(dtype=float).reshape(-1, 1)
    means_j_arr = means_j.to_numpy(dtype=float).reshape(-1, 1)

    diff_pi = (
        means_pi - means_p_arr - means_i.to_numpy(dtype=float).reshape(1, -1) + grand_mean
    )
    ss_pi = float(n_j * np.sum(diff_pi ** 2))
    diff_pj = (
        means_pj - means_p_arr - means_j.to_numpy(dtype=float).reshape(1, -1) + grand_mean
    )
    ss_pj = float(n_i * np.sum(diff_pj ** 2))
    diff_ij = (
        means_ij - means_i_arr - means_j.to_numpy(dtype=float).reshape(1, -1) + grand_mean
    )
    ss_ij = float(n_p * np.sum(diff_ij ** 2))

    # Residual (= three-way interaction confounded with error under
    # one obs per cell) by subtraction.
    ss_residual = max(ss_total - ss_p - ss_i - ss_j - ss_pi - ss_pj - ss_ij, 0.0)

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

    # Keep the existing zero-clipped coefficient calculation. Signed
    # method-of-moments estimates are retained separately as RawVariance.
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
    one random facet uses the two-way ANOVA below.
    """
    work, _, reason = _prepare_g_d_data(obs_df, object_facet, random_facets, score_col)
    if reason:
        return pd.DataFrame()
    n_obs = int(len(work))

    grand_mean = float(work[score_col].mean())
    ss_total = float(((work[score_col] - grand_mean) ** 2).sum())

    # Two random facets use the full three-way decomposition. Only the
    # already validated one-random-facet case proceeds to the two-way path.
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

    rows = []
    explained_ss = 0.0
    facet_terms = [object_facet] + list(random_facets)
    for facet in facet_terms:
        grouped = work.groupby(facet, observed=False)
        means = grouped[score_col].mean()
        ns = grouped[score_col].size().astype(float)
        k = int(len(means))
        if k < 2:
            continue
        ss_facet = float(np.sum(ns.to_numpy() * (means.to_numpy() - grand_mean) ** 2))
        df_facet = k - 1
        ms_facet = ss_facet / df_facet if df_facet > 0 else np.nan
        # Validation guarantees equal counts; the harmonic mean equals n_obs / k.
        ns_arr = ns.to_numpy(dtype=float)
        ns_valid = ns_arr[ns_arr > 0]
        if ns_valid.size > 0 and np.all(ns_valid > 0):
            harmonic_mean_n = float(ns_valid.size / np.sum(1.0 / ns_valid))
        else:
            harmonic_mean_n = float(n_obs / max(k, 1))
        mean_obs_per_level = harmonic_mean_n
        explained_ss += ss_facet
        rows.append({
            "Source": facet,
            "Term": "main",
            "SumSquares": ss_facet,
            "df": df_facet,
            "MeanSquares": ms_facet,
            "MeanObsPerLevel": mean_obs_per_level,
            "Levels": k,
        })
    ss_residual = max(float(ss_total - explained_ss), 0.0)
    df_residual = max(int(n_obs - sum(int(r["Levels"]) for r in rows) + len(rows) - 1), 1)
    ms_residual = ss_residual / df_residual if df_residual > 0 else np.nan
    rows.append({
        "Source": "Residual",
        "Term": "residual",
        "SumSquares": ss_residual,
        "df": df_residual,
        "MeanSquares": ms_residual,
        "MeanObsPerLevel": float("nan"),
        "Levels": None,
    })

    out = pd.DataFrame(rows)
    sigma2_e = ms_residual if np.isfinite(ms_residual) else np.nan
    variances = []
    for _, row in out.iterrows():
        source = row["Source"]
        if source == "Residual":
            variances.append(float(sigma2_e) if np.isfinite(sigma2_e) else np.nan)
            continue
        ms = float(row["MeanSquares"])
        m_per = float(row["MeanObsPerLevel"]) if pd.notna(row["MeanObsPerLevel"]) else np.nan
        if not (np.isfinite(ms) and np.isfinite(sigma2_e) and np.isfinite(m_per) and m_per > 0):
            variances.append(np.nan)
        else:
            variances.append((ms - sigma2_e) / m_per)
    out["RawVariance"] = variances
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
    """Compute G / Phi from a variance-component table.

    Two formulas are used depending on what the variance-component
    table carries:

    * **Full 3-way decomposition** (the table has explicit two-way
      interaction terms ``object:facetA``, ``object:facetB``,
      ``facetA:facetB``): the canonical Brennan (2001, Eq. 3.18)
      formulas with explicit interaction variances:

          sigma2(delta) = sigma2_pi / n_i + sigma2_pj / n_j + sigma2_pij / (n_i n_j)
          sigma2(Delta) = sigma2(delta) + sigma2_i / n_i + sigma2_j / n_j + sigma2_ij / (n_i n_j)
          G   = sigma2_p / (sigma2_p + sigma2(delta))
          Phi = sigma2_p / (sigma2_p + sigma2(Delta))

      Here ``sigma2_pij`` is the residual term (three-way interaction
      confounded with error under one observation per cell).

    * **One-random-facet design**: the non-object facet
      contributes ``sigma2_facet / n_facet`` to absolute error
      and the residual contributes ``sigma2_residual / n_total`` to
      both relative and absolute error. The person-by-facet interaction
      is confounded with residual error. The G-study entry point refuses
      unbalanced data and designs with more than two random facets.
    """
    if not isinstance(var_components, pd.DataFrame) or var_components.empty:
        return {"G": np.nan, "Phi": np.nan, "details": {}}

    # Lookup by Source name; the table may also carry a Term column
    # ("main" / "two-way" / "residual") to flag interaction rows.
    lookup = {
        str(row["Source"]): float(row["Variance"])
        for _, row in var_components.iterrows()
    }
    sigma2_p = lookup.get(object_facet, np.nan)
    sigma2_e = lookup.get("Residual", np.nan)
    if not (np.isfinite(sigma2_p) and np.isfinite(sigma2_e) and sigma2_p > 0):
        return {"G": np.nan, "Phi": np.nan, "details": {"sigma2_p": sigma2_p, "sigma2_e": sigma2_e}}

    facet_obs = facet_observations or {}
    has_term = "Term" in var_components.columns
    two_way_rows = (
        var_components[var_components["Term"] == "two-way"]
        if has_term else pd.DataFrame()
    )
    using_interactions = not two_way_rows.empty

    if using_interactions:
        # Full 3-way decomposition. Identify the two random facets
        # from the two-way Source labels that involve the object_facet.
        random_facets: list[str] = []
        for src in two_way_rows["Source"]:
            parts = [p for p in str(src).split(":") if p != object_facet]
            for p in parts:
                if p not in random_facets:
                    random_facets.append(p)
        n_a = max(int(facet_obs.get(random_facets[0], 1)), 1) if random_facets else 1
        n_b = max(int(facet_obs.get(random_facets[1], 1)), 1) if len(random_facets) > 1 else 1

        def lookup_two_way(facets: tuple[str, str]) -> float:
            for src, var in lookup.items():
                parts = set(str(src).split(":"))
                if parts == set(facets):
                    return float(var)
            return 0.0

        i_facet = random_facets[0] if random_facets else None
        j_facet = random_facets[1] if len(random_facets) > 1 else None
        sigma2_pi = lookup_two_way((object_facet, i_facet)) if i_facet else 0.0
        sigma2_pj = lookup_two_way((object_facet, j_facet)) if j_facet else 0.0
        sigma2_ij = lookup_two_way((i_facet, j_facet)) if i_facet and j_facet else 0.0
        sigma2_i = lookup.get(i_facet, 0.0) if i_facet else 0.0
        sigma2_j = lookup.get(j_facet, 0.0) if j_facet else 0.0

        # Brennan (2001) Eq. 3.18 / 3.19.
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

    # Main-effects-only fallback.
    sigma2_main_per_obs = 0.0
    for source, variance in lookup.items():
        if source in {object_facet, "Residual"} or not np.isfinite(variance):
            continue
        n_per = facet_obs.get(source, 1)
        sigma2_main_per_obs += variance / max(int(n_per), 1)

    # Total number of conditions-of-measurement observations per object.
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
    import itertools
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
) -> dict:
    """Technical person-resampling calculation; not qualified for inference.

    Percentiles summarize finite replicates only. Each coefficient/design
    retains its success count and failure reasons against all attempts.
    This calculation does not establish coverage for random-facet G/D targets.
    """
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
            # Each sampled cluster is a distinct measurement object.
            sub[object_facet] = f"BOOT_{k:04d}_{person}"
            parts.append(sub)
        boot_df = pd.concat(parts, ignore_index=True)
        why = ""
        try:
            components = _crossed_anova_variance_components(
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
                coef = _g_phi_from_variance_components(
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
    bootstrap_ci: bool = False,
    n_bootstrap: int = 200,
    bootstrap_confidence: float = 0.95,
    bootstrap_seed: int | None = 20260601,
) -> dict:
    """Use retained MFRM scores; share the standalone G/D calculation and limits."""
    valid = isinstance(res, dict)
    config = res.get("config") if valid and isinstance(res.get("config"), dict) else {}
    prep = res.get("prep", {}) if valid else {}
    if random_facets is None:
        random_facets = [f for f in (config or {}).get("facet_names", []) if f != object_facet]
    bundle = _compute_generalizability_study(
        prep.get("data") if isinstance(prep, dict) else None,
        object_facet=object_facet, random_facets=random_facets, score_col=score_col,
        n_max_per_facet=n_max_per_facet, score_basis="retained_preprocessed_scores",
        bootstrap_ci=bootstrap_ci, n_bootstrap=n_bootstrap,
        bootstrap_confidence=bootstrap_confidence, bootstrap_seed=bootstrap_seed)
    if not valid:
        bundle.update(reason_code="invalid_result", reason=_GD_INPUT_REASONS["invalid_result"])
    return bundle


def compute_generalizability_from_data(
    data: pd.DataFrame, *, object_facet: str, random_facets: list[str],
    score_col: str = "Score", n_max_per_facet: int = 10,
    target_population: str = "",
) -> dict:
    """Analyze observed scores directly, without MFRM preparation or estimation.

    All non-object facets are random; no fixed-facet or weighted model is fitted.
    The declared generalization population is provenance, not a verified assumption.
    """
    return _compute_generalizability_study(
        data, object_facet=object_facet, random_facets=random_facets, score_col=score_col,
        n_max_per_facet=n_max_per_facet, score_basis="original_observed_scores",
        target_population=target_population)


def _compute_generalizability_study(
    data: pd.DataFrame,
    *,
    object_facet: str = "Person",
    random_facets: list[str] | None = None,
    n_max_per_facet: int = 10,
    score_col: str = "Score",
    score_basis: str,
    target_population: str = "",
    bootstrap_ci: bool = False,
    n_bootstrap: int = 200,
    bootstrap_confidence: float = 0.95,
    bootstrap_seed: int | None = 20260601,
) -> dict:
    """Shared balanced ANOVA, design validation and qualification for both entries."""
    design = {}

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
    data, design, reason = _prepare_g_d_data(data, object_facet, random_facets, score_col)
    design.update(score_basis=score_basis, target_population=target_population,
                  fixed_facets=[], generalization_basis="all_non_object_facets_random")
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
            data, object_facet=object_facet, random_facets=random_facets, score_col=score_col)
    if var_components.empty or not np.isfinite(
        var_components[["RawVariance", "Variance", "SumSquares", "MeanSquares"]].to_numpy(dtype=float)
    ).all():
        return unavailable("invalid_components")
    observed_n = design["observed_n"]
    observed_coefficients = _g_phi_from_variance_components(
        var_components, object_facet=object_facet, facet_observations=observed_n)
    if not all(np.isfinite(observed_coefficients[k]) for k in ("G", "Phi")):
        return unavailable("undefined_coefficients", var_components)
    facet_grids = {
        f: tuple(g for g in _GD_STUDY_DEFAULT_GRID
                 if g <= min(max(observed_n[f] * 2, 4), n_max_per_facet))
        for f in random_facets
    }
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
        # Build the list of extra designs to evaluate per replicate.
        # Each d_study row contributes a facet_observations dict; the
        # bootstrap loop computes the variance components ONCE per
        # replicate and reuses them across all grid points, so the
        # marginal cost of the D-study CI bands is small (one G/Phi
        # closed-form evaluation per replicate per grid point).
        extra_designs: list[dict[str, int]] = []
        if isinstance(d_study, pd.DataFrame) and not d_study.empty:
            for _, row in d_study.iterrows():
                design_dict: dict[str, int] = {}
                for facet in random_facets:
                    if facet in row.index and pd.notna(row[facet]):
                        design_dict[facet] = int(row[facet])
                if design_dict:
                    extra_designs.append(design_dict)
        bootstrap = _bootstrap_g_phi_ci(
            data,
            object_facet=object_facet,
            random_facets=list(random_facets),
            score_col=score_col,
            facet_observations=observed_n,
            n_bootstrap=int(n_bootstrap),
            confidence=float(bootstrap_confidence),
            seed=bootstrap_seed,
            extra_designs=extra_designs,
        )
        # Merge per-design CIs into d_study DataFrame.
        if (
            isinstance(d_study, pd.DataFrame)
            and not d_study.empty
            and bootstrap.get("available")
        ):
            extra_replicates = bootstrap.get("extra_design_replicates", [])
            # Match by (n_facet1, n_facet2, ...) tuple.
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
    return {"G": coefficients["G"], "Phi": coefficients["Phi"],
            "RelativeErrorVariance": details["sigma2_relative_error"],
            "AbsoluteErrorVariance": details["sigma2_absolute_error"]}


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
