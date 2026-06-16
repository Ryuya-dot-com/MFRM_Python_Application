"""Resource preflight checks for hosted Streamlit runs."""

from __future__ import annotations

import os

import pandas as pd


ESTIMATION_HOSTED_MAX_OBS = 100_000
ESTIMATION_HOSTED_MAX_PARAMETERS = 30_000
ESTIMATION_HOSTED_MAX_MML_EVALS = 5_000_000
ESTIMATION_HOSTED_MAX_BIAS_CELLS = 250_000


def allow_large_hosted_run_override() -> bool:
    """Return whether hosted-resource hard blocks are explicitly overridden."""
    return str(os.environ.get("MFRM_ALLOW_LARGE_HOSTED_RUNS", "")).strip().lower() in {
        "1",
        "true",
        "yes",
        "on",
    }


def build_estimation_resource_preflight(
    *,
    data: pd.DataFrame,
    person_col: str,
    score_col: str,
    facet_cols: list[str],
    model_type: str,
    est_method: str,
    maxit: int,
    quad_points: int,
    bias_mode: str,
    bias_pairs_available: list[tuple[str, str]],
    selected_bias_pair: tuple[str, str] | None,
    compute_strict_marginal: bool,
    strict_marginal_pairwise: bool,
    strict_marginal_max_pair_cells: int,
    generate_figure_exports: bool,
    hosted_max_obs: int = ESTIMATION_HOSTED_MAX_OBS,
    hosted_max_parameters: int = ESTIMATION_HOSTED_MAX_PARAMETERS,
    hosted_max_mml_evals: int = ESTIMATION_HOSTED_MAX_MML_EVALS,
    hosted_max_bias_cells: int = ESTIMATION_HOSTED_MAX_BIAS_CELLS,
) -> dict:
    """Estimate hosted-runtime cost before the user starts a fit."""
    if not isinstance(data, pd.DataFrame) or data.empty:
        return {"block": False, "rows": []}

    n_obs = int(len(data))
    n_persons = int(data[person_col].nunique(dropna=True)) if person_col in data.columns else 0
    facet_level_counts = {
        str(facet): int(data[facet].nunique(dropna=True))
        for facet in facet_cols
        if facet in data.columns
    }
    score_unique = (
        int(pd.to_numeric(data[score_col], errors="coerce").dropna().nunique())
        if score_col in data.columns else 0
    )
    n_steps = 0
    if score_unique >= 2:
        if model_type in {"PCM", "GPCM"}:
            step_levels = max(facet_level_counts.values(), default=1)
            n_steps = max(0, score_unique - 1) * step_levels
            if model_type == "GPCM":
                n_steps += step_levels
        else:
            n_steps = max(0, score_unique - 1)
    approx_params = (
        (n_persons if est_method == "JMLE" else 0)
        + sum(facet_level_counts.values())
        + int(n_steps)
    )
    mml_evals = int(n_obs) * int(quad_points or 0) if est_method == "MML" else 0

    pair_counts = {"Person": n_persons, **facet_level_counts}
    if bias_mode == "All facet pairs":
        bias_pairs = list(bias_pairs_available or [])
    elif bias_mode == "Selected pair" and selected_bias_pair is not None:
        bias_pairs = [selected_bias_pair]
    else:
        bias_pairs = []
    bias_cells = int(sum(
        max(0, pair_counts.get(str(a), 0)) * max(0, pair_counts.get(str(b), 0))
        for a, b in bias_pairs
    ))

    rows: list[dict] = []

    def add(metric: str, value: object, limit: object, status: str, action: str) -> None:
        rows.append({
            "Metric": metric,
            "Value": value,
            "HostedLimit": limit,
            "Status": status,
            "Action": action,
        })

    add(
        "Observations",
        f"{n_obs:,}",
        f"{hosted_max_obs:,}",
        "Block" if n_obs > hosted_max_obs else ("Review" if n_obs > 50_000 else "OK"),
        "Sample rows, run Fast preview first, or run locally for large datasets.",
    )
    add(
        "Approx. free parameters",
        f"{approx_params:,}",
        f"{hosted_max_parameters:,}",
        "Block" if approx_params > hosted_max_parameters else ("Review" if approx_params > 12_000 else "OK"),
        "Drop high-cardinality facets or avoid treating row IDs as facets.",
    )
    if est_method == "MML":
        add(
            "MML quadrature evaluations",
            f"{mml_evals:,}",
            f"{hosted_max_mml_evals:,}",
            "Block" if mml_evals > hosted_max_mml_evals else ("Review" if mml_evals > 2_000_000 else "OK"),
            "Reduce quadrature points, use JMLE preview, or run MML locally.",
        )
    if bias_pairs:
        add(
            "Bias/interaction candidate cells",
            f"{bias_cells:,}",
            f"{hosted_max_bias_cells:,}",
            "Block" if bias_cells > hosted_max_bias_cells else ("Review" if bias_cells > 75_000 else "OK"),
            "Use Selected pair or Skip for the first run; avoid all-pair Person interactions on large data.",
        )
    if compute_strict_marginal and strict_marginal_pairwise:
        add(
            "Strict marginal pairwise cap",
            f"{strict_marginal_max_pair_cells:,}",
            "800 recommended",
            "Review" if strict_marginal_max_pair_cells > 2_000 else "OK",
            "Keep pairwise marginal diagnostics capped for hosted use.",
        )
    if maxit > 3_000 and n_obs > 25_000:
        add(
            "Max iterations",
            f"{maxit:,}",
            "3,000 with large data",
            "Review",
            "Use a smaller preview first; high maxit on large designs can pin the hosted worker.",
        )
    if generate_figure_exports and n_obs > 50_000:
        add(
            "Figure export bundle",
            "enabled",
            "large-data caution",
            "Review",
            "Disable figure export for preview runs; generate publication figures after a stable fit.",
        )

    hard_block = any(row["Status"] == "Block" for row in rows)
    override = allow_large_hosted_run_override()
    return {
        "block": bool(hard_block and not override),
        "override": override,
        "rows": rows,
        "approx_params": approx_params,
        "mml_evals": mml_evals,
        "bias_cells": bias_cells,
    }


def estimation_resource_preflight_status(preflight: dict) -> str:
    """Return the display severity for a preflight payload."""
    rows = preflight.get("rows", []) if isinstance(preflight, dict) else []
    if not rows:
        return "OK"
    statuses = {str(row.get("Status", "OK")) for row in rows if isinstance(row, dict)}
    if "Block" in statuses and not preflight.get("override"):
        return "Block"
    if "Review" in statuses or preflight.get("override"):
        return "Review"
    return "OK"


def should_render_estimation_resource_preflight(preflight: dict) -> bool:
    """Avoid alert fatigue: only surface resource budget details when action is needed."""
    return estimation_resource_preflight_status(preflight) != "OK"

