"""UI-independent conditional bias analysis and qualified public output.

Calibration is held fixed. Moving these calculations does not qualify tests,
intervals, pairwise claims or conclusions that bias is present or absent.
"""
from __future__ import annotations

from itertools import combinations
import json

import numpy as np
import pandas as pd
from scipy.optimize import root_scalar
from scipy.special import logsumexp
from scipy.stats import chi2

from mfrm_app import decision_stability as _decision_stability
from mfrm_app import jmle_boundary as _jmle_boundary
from mfrm_app import output_qualification as _output_qualification
from mfrm_app.engine import build_indices, category_prob_rsm, compute_eta, loglik_rsm


def get_extreme_levels(obs_df, facet_names, rating_min, rating_max):
    extreme_levels = {}
    for facet in facet_names:
        if facet not in obs_df.columns:
            extreme_levels[facet] = set()
            continue
        stat = (
            obs_df.groupby(facet, observed=False)["Observed"]
            .agg(MinScore="min", MaxScore="max")
            .reset_index()
        )
        extreme = stat[
            ((stat["MinScore"] == rating_min) & (stat["MaxScore"] == rating_min))
            | ((stat["MinScore"] == rating_max) & (stat["MaxScore"] == rating_max))
        ]
        extreme_levels[facet] = set(extreme[facet].astype(str).tolist())
    return extreme_levels


# =============================================================================
# Conditional likelihood slices for an additive bias shift
# Other calibration parameters remain fixed: this is not a nuisance profile.
# A chi-square(1) cutoff is a candidate asymptotic reference only for regular
# independent unit-weight observations with known calibration. Its validity
# after calibration on the same data requires G3/G4 operating-characteristic
# validation. A slope-squared information identity alone does not establish it.
# =============================================================================


def _safe_nll_value(nll, b):
    """Robust wrapper for evaluating a negative log-likelihood closure.

    Catches arbitrary exceptions (degenerate observation slices, log of
    zero from probabilities that collapse, etc.) and returns ``nan`` so
    the conditional likelihood search and the LR statistic can degrade
    gracefully instead of propagating a hard error from a single cell.
    """
    try:
        val = float(nll(b))
    except Exception:
        return float("nan")
    return val if np.isfinite(val) else float("nan")


def _profile_bias_ci(nll, estimate, nll_min, max_abs, level=0.95):
    """Conditional LR-cutoff interval (legacy helper name; no nuisance profiling).

    Parameters
    ----------
    nll : callable
        Negative log-likelihood as a function of the bias scalar.
    estimate : float
        The MLE ``b_hat`` for the bias parameter.
    nll_min : float
        ``nll(estimate)`` at the MLE. Pre-evaluated so the helper does
        not duplicate the cost of the optimisation.
    max_abs : float
        Half-width of the search bracket. The CI is computed inside
        ``[-max_abs, max_abs]``.
    level : float, default 0.95
        Nominal CI coverage. Maps to the chi-square cutoff
        ``chi2.ppf(level, df=1)``.

    Returns
    -------
    dict with keys ``lower``, ``upper``, ``level``, ``status``. ``status``
    is one of ``"ok"`` (both endpoints found via interior roots),
    ``"limited by search range"`` (one or both endpoints hit the
    ``max_abs`` boundary because the likelihood never fell far enough
    inside the bracket), or ``"not available"`` (the search did not
    return a finite endpoint).
    """
    if (not np.isfinite(estimate) or not np.isfinite(nll_min) or not np.isfinite(max_abs)
            or max_abs <= 0 or abs(estimate) > max_abs):
        return {"lower": float("nan"), "upper": float("nan"), "level": float(level), "status": "not available"}
    if not (0.0 < level < 1.0):
        return {"lower": float("nan"), "upper": float("nan"), "level": float(level), "status": "not available"}
    cutoff = float(chi2.ppf(level, df=1))

    def target(b):
        val = _safe_nll_value(nll, b)
        if not np.isfinite(val):
            return float("nan")
        return 2.0 * (val - nll_min) - cutoff

    def root_side(bound):
        if abs(bound - estimate) <= np.sqrt(np.finfo(float).eps):
            return {"value": float(bound), "status": "limited by search range"}
        t_bound = target(bound)
        if not np.isfinite(t_bound):
            return {"value": float("nan"), "status": "not available"}
        if t_bound <= 0:
            # The likelihood never drops by the cutoff inside the bracket,
            # so the CI extends to or beyond ``bound``. Report the bound
            # so the caller can flag the truncation in the status field.
            return {"value": float(bound), "status": "limited by search range"}
        lo = min(bound, estimate)
        hi = max(bound, estimate)
        try:
            sol = root_scalar(target, bracket=(lo, hi), method="brentq", xtol=1e-8)
        except Exception:
            return {"value": float("nan"), "status": "not available"}
        if not getattr(sol, "converged", False) or not np.isfinite(sol.root):
            return {"value": float("nan"), "status": "not available"}
        return {"value": float(sol.root), "status": "ok"}

    lower = root_side(-abs(float(max_abs)))
    upper = root_side(abs(float(max_abs)))
    statuses: list[str] = []
    for s in (lower["status"], upper["status"]):
        if s not in statuses:
            statuses.append(s)
    if statuses == ["ok"]:
        status = "ok"
    else:
        status = "; ".join(statuses)
    return {
        "lower": lower["value"],
        "upper": upper["value"],
        "level": float(level),
        "status": status,
    }


def _fit_conditional_bias_cell(nll, probability, scores, slopes, weights, max_abs):
    """Solve the concave, fixed-calibration scalar problem without a false MLE.

    d NLL/db = sum w*a*(E[X]-X), d² NLL/db² = sum w*a²*Var(X).
    Unit-weight chi-square calibration is an unqualified development candidate;
    fitted calibration uncertainty is not included. Other weights retain only
    point/curvature diagnostics. A search bound is never a finite interior MLE.
    """
    out = {"estimate": np.nan, "candidate": np.nan, "status": "invalid_input",
           "score": np.nan, "information": np.nan, "se": np.nan,
           "lr": np.nan, "p": np.nan, "reference_status": "not_available",
           "ci": {"lower": np.nan, "upper": np.nan, "level": .95,
                  "status": "not available"}, "probabilities": None, "error": ""}
    scores, slopes, weights = [np.asarray(x, dtype=float) for x in (scores, slopes, weights)]
    if (scores.ndim != 1 or not scores.size or slopes.shape != scores.shape
            or weights.shape != scores.shape or not np.isfinite(max_abs) or max_abs <= 0
            or not all(np.isfinite(x).all() for x in (scores, slopes, weights))
            or np.any(slopes <= 0) or np.any(weights < 0) or not np.any(weights > 0)):
        return out
    active = weights > 0

    def moments(b):
        probs = np.asarray(probability(b), dtype=float)
        if (probs.ndim != 2 or probs.shape[0] != scores.size or probs.shape[1] < 2
                or not np.isfinite(probs).all() or np.any(probs < 0)
                or not np.allclose(probs.sum(axis=1), 1., atol=1e-12, rtol=1e-12)):
            raise ValueError("Invalid category probabilities in bias cell")
        k = np.arange(probs.shape[1], dtype=float)
        mean = probs @ k
        # Central moment avoids cancellation; no small-variance row deletion.
        variance = np.sum(probs * (k[None, :] - mean[:, None])**2, axis=1)
        return probs, mean, variance

    def gradient(b):
        _, mean, _ = moments(b)
        return float(np.sum(weights * slopes * (mean-scores)))

    try:
        probs_null, mean_null, _ = moments(0.)
        max_score = probs_null.shape[1]-1
        if np.any(scores != np.floor(scores)) or np.any(scores < 0) or np.any(scores > max_score):
            return out
        # Positive slopes and finite thresholds give a strictly monotone score.
        if np.all(scores[active] == 0) or np.all(scores[active] == max_score):
            out.update(status="no_finite_mle", candidate=-max_abs if np.all(scores[active] == 0) else max_abs)
            return out
        lo, hi = gradient(-max_abs), gradient(max_abs)
        if not (np.isfinite(lo) and np.isfinite(hi)):
            raise ValueError("Nonfinite bias score at search bounds")
        if lo >= 0 or hi <= 0:
            out.update(status="search_range_limited", candidate=-max_abs if lo >= 0 else max_abs)
            return out
        solution = root_scalar(gradient, bracket=(-max_abs, max_abs), method="brentq", xtol=1e-12)
        if not solution.converged or not np.isfinite(solution.root):
            out["status"] = "solver_failed"
            return out
        estimate = float(solution.root)
        probs, mean, variance = moments(estimate)
        score = float(np.sum(weights*slopes*(mean-scores)))
        info = float(np.sum(weights*slopes**2*variance))
        out.update(candidate=estimate, score=score, information=info)
        # Numerical stationarity uses per weighted slope/category score scale.
        score_scale = max(1., float(np.sum(weights*slopes)*max_score))
        if not np.isfinite(info) or info <= 0 or abs(score) > 1e-9*score_scale:
            out["status"] = "stationarity_failed"
            return out
        nll_hat, nll_null = _safe_nll_value(nll, estimate), _safe_nll_value(nll, 0.)
        if not np.isfinite(nll_hat) or not np.isfinite(nll_null) or nll_hat > nll_null + 1e-10*max(1., abs(nll_null)):
            out["status"] = "likelihood_failed"
            return out
        lr_raw_nll = 2*(nll_null-nll_hat)
        lr = lr_raw_nll
        shift = slopes*estimate
        if np.max(np.abs(shift))*max_score <= .5:
            # Exponential tilting gives exactly the same likelihood ratio.
            # Centering and expm1/log1p avoid subtracting nearly equal NLLs
            # near the null; chi-square survival is sensitive to that loss.
            centered_shift = shift[:, None]*(np.arange(max_score+1)-mean_null[:, None])
            cumulant = np.log1p(np.sum(probs_null*np.expm1(centered_shift), axis=1))
            lr = float(2*np.sum(weights*(shift*(scores-mean_null)-cumulant)))
            allowance = 64*np.finfo(float).eps*max(1., abs(nll_null), abs(nll_hat)) + 1e-10*max(1., abs(lr_raw_nll))
            if not np.isfinite(lr) or abs(lr-lr_raw_nll) > allowance:
                out["status"] = "likelihood_failed"
                out["error"] = "Conditional probabilities and NLL likelihood ratios disagree"
                return out
        lr = max(0., lr)
        out.update(estimate=estimate, status="interior", se=1/np.sqrt(info),
                   lr=lr, lr_raw_nll=lr_raw_nll, probabilities=probs)
        if np.all(weights[active] == 1.):
            out.update(p=float(chi2.sf(lr, 1)), reference_status="unit_weight_conditional_lr_candidate",
                       ci=_profile_bias_ci(nll, estimate, nll_hat, max_abs, level=.95))
        else:
            out["reference_status"] = "weighted_reference_unqualified"
    except Exception as exc:
        out.update(status="calculation_failed", error=f"{type(exc).__name__}: {exc}",
                   estimate=np.nan, se=np.nan, lr=np.nan, p=np.nan, probabilities=None,
                   reference_status="not_available",
                   ci={"lower": np.nan, "upper": np.nan, "level": .95, "status": "not available"})
    return out


def _prepare_bias_cell_functions(model, eta, scores, weights, step_cum,
                                 step_cum_mat, step_idx, slopes, slope_idx):
    """Prepare fixed cell slices once; preserve the model's summation order."""
    if model == "RSM":
        return (lambda b: -loglik_rsm(eta+b, scores, step_cum, weight=weights),
                lambda b: category_prob_rsm(eta+b, step_cum))

    k_cat = step_cum_mat.shape[1]
    k = np.arange(k_cat, dtype=float if model == "GPCM" else int)
    groups = [np.where(step_idx == c)[0] for c in range(step_cum_mat.shape[0])]
    rows = np.concatenate(groups)
    ends = np.cumsum([len(group) for group in groups])
    slices = [slice(end-len(group), end) for group, end in zip(groups, ends) if len(group)]
    eta_rows, steps = eta[rows], step_cum_mat[step_idx[rows]]
    row_numbers, score_rows = np.arange(rows.size), scores[rows]
    weight_rows = weights[rows] if weights is not None else None
    slope_rows = np.asarray(slopes, dtype=float)[slope_idx[rows]] if model == "GPCM" else None

    def log_numerator(b):
        linear = np.outer(eta_rows+b, k)-steps
        return slope_rows[:, None]*linear if slope_rows is not None else linear

    def nll(b):
        log_num = log_numerator(b)
        diff = log_num[row_numbers, score_rows]-logsumexp(log_num, axis=1)
        if weight_rows is not None:
            diff = diff*weight_rows
        total = 0.0
        for group in slices:
            total += float(np.sum(diff[group]))
        return -total

    def probability(b):
        probs = np.zeros((len(eta), k_cat))
        log_num = log_numerator(b)
        probs[rows, :] = np.exp(log_num-logsumexp(log_num, axis=1)[:, None])
        return probs

    return nll, probability


def estimate_bias_interaction(
    res,
    diagnostics,
    facet_a,
    facet_b,
    max_abs=10.0,
    omit_extreme=True,
    max_iter=4,
    tol=1e-3,
):
    """Fit each scalar cell shift once with all calibration values fixed.

    ``max_iter`` and ``tol`` are retained as ignored legacy arguments: there
    is no outer updating algorithm. Repeating an unchanged scalar problem
    cannot establish convergence of a joint bias/calibration estimate.
    """
    # Returns either:
    #   - a full results dict (table, summary, chi-sq, ...) on success
    #   - a dict with a single "_skip_reason" key explaining why
    #     the computation was not possible. This lets the UI surface
    #     the cause instead of showing "not computed" silently.
    if res is None or diagnostics is None:
        return {"_skip_reason": "estimation result or diagnostics missing"}
    identifiability = res.get("identifiability", {}) if isinstance(res, dict) else {}
    identifiability_status = (
        str(identifiability.get("status", ""))
        if isinstance(identifiability, dict) else ""
    )
    if identifiability_status == "eta_rank_deficient":
        return {
            "_skip_reason": (
                "JMLE person/facet eta coordinates are structurally rank deficient; "
                "bias interaction analysis is withheld until the design is connected "
                "or otherwise identified."
            )
        }
    if identifiability_status == "audit_size_limit":
        return {
            "_skip_reason": (
                "JMLE structural identifiability was not established because the exact "
                "rank audit exceeded its technical size limit; bias interaction analysis "
                "is withheld."
            )
        }
    if _jmle_facet_boundary_hold(res.get("config")):
        return {"_skip_reason": _jmle_boundary.HOLD_MESSAGE}
    obs_df = diagnostics.get("obs")
    if obs_df is None or obs_df.empty:
        return {"_skip_reason": "observation table is empty (no usable rows for bias analysis)"}
    if facet_a == facet_b:
        return {
            "_skip_reason": (
                f"bias analysis requires two distinct facets, not `{facet_a}` × `{facet_a}`"
            )
        }

    facet_names = ["Person"] + res["config"]["facet_names"]
    if facet_a not in facet_names or facet_b not in facet_names:
        missing = [f for f in (facet_a, facet_b) if f not in facet_names]
        return {
            "_skip_reason": (
                f"facet(s) {missing} not in estimated facets {facet_names} "
                "(did you rename a column after running the estimation?)"
            )
        }

    if not np.isfinite(max_abs) or max_abs <= 0:
        raise ValueError("max_abs must be finite and positive")
    prep = res["prep"]
    config = res["config"]
    params = res["params"]
    idx = build_indices(prep, step_facet=config["step_facet"], slope_facet=config.get("slope_facet"))
    theta_hat = params["theta"] if config["method"] == "JMLE" else res["facets"]["person"]["Estimate"].to_numpy()
    eta_base = compute_eta(idx, params, config, theta_override=theta_hat if theta_hat.size else None)
    score_k = idx["score_k"]
    weight = idx.get("weight")
    step_idx = idx.get("step_idx")
    slope_idx = idx.get("slope_idx")

    if config["model"] == "RSM":
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
        step_cum_mat = None
    else:
        step_cum = None
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])

    if omit_extreme:
        extreme_levels = get_extreme_levels(
            obs_df,
            [facet_a, facet_b],
            prep["rating_min"],
            prep["rating_max"],
        )
    else:
        extreme_levels = {facet_a: set(), facet_b: set()}

    measures = diagnostics.get("measures", pd.DataFrame())
    meas_map = {}
    se_map = {}
    if not measures.empty:
        for row in measures.itertuples():
            meas_map[(row.Facet, str(row.Level))] = float(row.Estimate)
            se_map[(row.Facet, str(row.Level))] = float(row.SE) if pd.notna(row.SE) else np.nan

    level_map = {facet: list(map(str, prep["levels"][facet])) for facet in facet_names if facet in prep["levels"]}
    group_indices = obs_df.groupby([facet_a, facet_b], observed=False).indices
    groups = []
    for (lvl_a, lvl_b), idx_rows in group_indices.items():
        lvl_a_str = str(lvl_a)
        lvl_b_str = str(lvl_b)
        idx_rows = np.asarray(idx_rows, dtype=int)
        if idx_rows.size == 0:
            continue
        groups.append({
            "key": (lvl_a_str, lvl_b_str),
            "idx": idx_rows,
            "excluded_extreme": bool(omit_extreme and (
                lvl_a_str in extreme_levels.get(facet_a, set())
                or lvl_b_str in extreme_levels.get(facet_b, set()))),
        })

    def _slice_cell(idx_rows):
        """Slice the captured per-observation arrays to a single bias cell."""
        eta_sub = eta_base[idx_rows]
        score_k_sub = score_k[idx_rows]
        weight_sub = weight[idx_rows] if weight is not None else None
        step_idx_sub = step_idx[idx_rows] if step_idx is not None else None
        slope_idx_sub = slope_idx[idx_rows] if slope_idx is not None else step_idx_sub
        return eta_sub, score_k_sub, weight_sub, step_idx_sub, slope_idx_sub

    qualification = _output_qualification.bias_cell_qualification().to_record()
    rows = []
    seq = 1
    for g in groups:
        lvl_a_str, lvl_b_str = g["key"]
        idx_rows = g["idx"]
        eta_sub, score_k_sub, weight_sub, step_idx_sub, slope_idx_sub = _slice_cell(idx_rows)
        w = np.asarray(weight_sub, dtype=float) if weight_sub is not None else np.ones(len(idx_rows))
        slopes = (np.asarray(params["slopes"], dtype=float)[slope_idx_sub]
                  if config["model"] == "GPCM" else np.ones(len(idx_rows)))
        if g["excluded_extreme"]:
            fit = {"estimate": np.nan, "candidate": np.nan, "status": "excluded_extreme_level",
                   "se": np.nan, "score": np.nan, "information": np.nan,
                   "lr": np.nan, "p": np.nan, "reference_status": "not_available",
                   "probabilities": None, "error": "",
                   "ci": {"lower": np.nan, "upper": np.nan, "level": .95, "status": "not available"}}
        else:
            nll_sub, probability_sub = _prepare_bias_cell_functions(
                config["model"], eta_sub, score_k_sub, weight_sub, step_cum,
                step_cum_mat, step_idx_sub, params.get("slopes"), slope_idx_sub)
            fit = _fit_conditional_bias_cell(
                nll_sub, probability_sub, score_k_sub, slopes, w, max_abs)
        bias_hat, se = fit["estimate"], fit["se"]
        probs = fit["probabilities"]
        infit = outfit = np.nan
        if probs is not None:
            k = np.arange(probs.shape[1])
            mean = probs @ k
            variance = np.sum(probs*(k[None, :]-mean[:, None])**2, axis=1)
            residual_sq = (score_k_sub-mean)**2
            if np.sum(w*variance) > 0:
                infit = float(np.sum(w*residual_sq)/np.sum(w*variance))
            active = w > 0
            if np.all(variance[active] > 0):
                outfit = float(np.sum(w[active]*residual_sq[active]/variance[active])/np.sum(w[active]))
        profile_ci = fit["ci"]
        likelihood_basis = (
            "conditional likelihood slice for one additive bias shift; theta, "
            "facets, steps, and slopes fixed; no nuisance reoptimization"
        )
        obs_slice = obs_df.iloc[idx_rows]
        w_obs = obs_slice["Weight"].to_numpy(dtype=float) if "Weight" in obs_slice.columns else np.ones(len(obs_slice))
        obs_score = float(np.nansum(obs_slice["Observed"].to_numpy(dtype=float) * w_obs))
        exp_score = float(np.nansum(obs_slice["Expected"].to_numpy(dtype=float) * w_obs))
        obs_count = float(np.nansum(w_obs))
        obs_exp_avg = (obs_score - exp_score) / obs_count if obs_count > 0 else np.nan

        n_obs = int(len(obs_slice))
        standardized = bias_hat/se if np.isfinite(bias_hat) and np.isfinite(se) and se > 0 else np.nan

        rows.append({
            "Sq": seq,
            "Observd Score": obs_score,
            "Expctd Score": exp_score,
            "Observd Count": obs_count,
            "Obs-Exp Average": obs_exp_avg,
            "Bias Size": bias_hat,
            "S.E.": se,
            "Standardized Bias": standardized,
            "t": np.nan,
            "d.f.": np.nan,
            "Prob.": fit["p"],
            "LR ChiSq": fit["lr"],
            "LR d.f.": 1 if np.isfinite(fit["p"]) else np.nan,
            "LR Prob.": fit["p"],
            "Primary P Method": "unit_weight_conditional_lr_v1",
            "Reference Status": fit["reference_status"],
            "Cell Status": fit["status"],
            "Cell Error": fit["error"],
            "Bounded Candidate": fit["candidate"],
            "Search Bound": float(max_abs),
            "Score Residual": fit["score"],
            "Conditional Information": fit["information"],
            "PositiveWeightObsN": int(np.count_nonzero(w > 0)),
            "PersonN": int(obs_slice.loc[w > 0, "Person"].nunique()) if "Person" in obs_slice else np.nan,
            "Weight Sum": float(np.sum(w)),
            "Weight Sum Squares": float(np.sum(w**2)),
            "Profile CI Lower": profile_ci["lower"],
            "Profile CI Upper": profile_ci["upper"],
            "Profile CI Level": profile_ci["level"],
            "Profile CI Status": profile_ci["status"],
            "Likelihood Basis": likelihood_basis,
            "SE Basis": (
                "conditional information for one additive bias shift at fitted "
                "theta, facet, step, and slope values; main-effect uncertainty "
                "is not propagated"
            ),
            "Multiplicity Basis": (
                "Holm candidate within each observed facet-pair cell family, "
                "including failed/excluded cells in its size; BH descriptive only"
            ),
            "InferenceTier": "technical_only",
            **qualification,
            "Infit": infit,
            "Outfit": outfit,
            "ObsN": n_obs,
            "FacetA": facet_a,
            "FacetA_Level": lvl_a_str,
            "FacetA_Index": level_map.get(facet_a, []).index(lvl_a_str) + 1 if lvl_a_str in level_map.get(facet_a, []) else np.nan,
            "FacetA_Measure": meas_map.get((facet_a, lvl_a_str), np.nan),
            "FacetA_SE": se_map.get((facet_a, lvl_a_str), np.nan),
            "FacetA_Sign": 1 if facet_a == "Person" else config.get("facet_signs", {}).get(facet_a, -1),
            "FacetB": facet_b,
            "FacetB_Level": lvl_b_str,
            "FacetB_Index": level_map.get(facet_b, []).index(lvl_b_str) + 1 if lvl_b_str in level_map.get(facet_b, []) else np.nan,
            "FacetB_Measure": meas_map.get((facet_b, lvl_b_str), np.nan),
            "FacetB_SE": se_map.get((facet_b, lvl_b_str), np.nan),
            "FacetB_Sign": 1 if facet_b == "Person" else config.get("facet_signs", {}).get(facet_b, -1),
        })
        seq += 1

    bias_tbl = pd.DataFrame(rows)
    if bias_tbl.empty:
        return {
            "_skip_reason": (
                f"no valid `{facet_a}` × `{facet_b}` cells after filtering "
                "(all-extreme levels were removed and no valid observations remain)"
            )
        }

    numeric_cols = ["Observd Score", "Expctd Score", "Observd Count", "Obs-Exp Average", "Bias Size", "S.E."]
    mean_row = bias_tbl[numeric_cols].mean(numeric_only=True)
    sd_pop_row = bias_tbl[numeric_cols].std(ddof=0, numeric_only=True)
    sd_sample_row = bias_tbl[numeric_cols].std(ddof=1, numeric_only=True)
    summary_tbl = pd.DataFrame(
        [mean_row, sd_pop_row, sd_sample_row],
        index=["Mean (Count: {})".format(len(bias_tbl)), "S.D. (Population)", "S.D. (Sample)"],
    ).reset_index().rename(columns={"index": "Statistic"})

    se_bias = bias_tbl["S.E."].to_numpy()
    bias_vals = bias_tbl["Bias Size"].to_numpy()
    w_chi = np.where(np.isfinite(se_bias) & (se_bias > 0), 1 / (se_bias ** 2), np.nan)
    ok = np.isfinite(w_chi) & np.isfinite(bias_vals)
    fixed_chi = (
        np.sum(w_chi[ok] * bias_vals[ok] ** 2) - (np.sum(w_chi[ok] * bias_vals[ok]) ** 2) / np.sum(w_chi[ok])
        if np.sum(ok) >= 2
        else np.nan
    )
    # A descriptive dispersion; no calibrated omnibus reference distribution.
    fixed_df = np.nan
    fixed_prob = np.nan
    chi_tbl = pd.DataFrame({
        "FixedChiSq": [fixed_chi],
        "FixedDF": [fixed_df],
        "FixedProb": [fixed_prob],
        "ContributingCells": [int(np.sum(ok))],
        "Status": ["descriptive_only; shared calibration uncertainty unavailable"],
    })

    return {
        "facet_a": facet_a,
        "facet_b": facet_b,
        "table": bias_tbl,
        "summary": summary_tbl,
        "chi_sq": chi_tbl,
        "iteration": pd.DataFrame([{
            "Pass": 1, "PlannedCells": len(bias_tbl),
            "AttemptedCells": int(bias_tbl["Cell Status"].ne("excluded_extreme_level").sum()),
            "FittedCells": int(bias_tbl["Cell Status"].eq("interior").sum()),
            "FailedCells": int((~bias_tbl["Cell Status"].isin(["interior", "excluded_extreme_level"])).sum()),
            "ExcludedCells": int(bias_tbl["Cell Status"].eq("excluded_extreme_level").sum()),
            "CalibrationUpdated": False,
            "MaxAbsScoreResidual": bias_tbl["Score Residual"].abs().max(),
            "MaxAbsBias": bias_tbl["Bias Size"].abs().max(),
        }]),
        "calculation_mode": "single_pass_fixed_calibration",
    }


def calc_bias_pairwise(bias_tbl, target_facet, context_facet, *, target_sign=None):
    """Technical local contrasts conditional on the fitted calibration.

    Bias Size is an additive eta shift b. For a target with eta coefficient s,
    the local measure is m + b/s and the contrast is (b1-b2)/s. The shared m
    cancels exactly. Cell SEs describe disjoint cells with calibration held
    fixed; joint calibration uncertainty is not available from this table.
    No t reference distribution or pairwise p-value is assigned.
    """
    if bias_tbl is None or bias_tbl.empty:
        return pd.DataFrame()
    if target_facet == context_facet:
        raise ValueError("Target and context must be distinct facets.")
    forward = bias_tbl["FacetA"].eq(target_facet) & bias_tbl["FacetB"].eq(context_facet)
    reverse = bias_tbl["FacetB"].eq(target_facet) & bias_tbl["FacetA"].eq(context_facet)
    if forward.any() and reverse.any():
        raise ValueError("The requested facet pair has mixed orientations; supply one orientation.")
    if not (forward.any() or reverse.any()):
        return pd.DataFrame()

    target_prefix, context_prefix = ("FacetA", "FacetB") if forward.any() else ("FacetB", "FacetA")
    sub = bias_tbl.loc[forward if forward.any() else reverse]
    sign_column = f"{target_prefix}_Sign"
    if sign_column in sub:
        signs = pd.to_numeric(sub[sign_column], errors="coerce")
        if not signs.isin([-1, 1]).all() or signs.nunique() != 1:
            raise ValueError("Recorded target signs must be a single value of +1 or -1.")
        if target_sign is None:
            target_sign = float(signs.iloc[0])
        elif not signs.eq(target_sign).all():
            raise ValueError("target_sign conflicts with the recorded target sign.")
    if target_sign is None:
        raise ValueError("target_sign is required for a bias table without recorded facet signs.")
    if target_sign not in (-1, 1):
        raise ValueError("target_sign must be +1 or -1.")
    if target_facet == "Person" and target_sign != 1:
        raise ValueError("The Person target sign must be +1 in this model.")
    if sub.duplicated([f"{target_prefix}_Level", f"{context_prefix}_Level"]).any():
        raise ValueError("Duplicate target/context cells make the comparison ambiguous.")

    rows = []
    qualification = _output_qualification.bias_pairwise_qualification().to_record()
    for tgt_level, group_df in sub.groupby(f"{target_prefix}_Level"):
        contexts = group_df[f"{context_prefix}_Level"].unique().tolist()
        if len(contexts) < 2:
            continue
        if group_df[f"{target_prefix}_Measure"].nunique(dropna=False) != 1:
            raise ValueError("The shared target measure differs across context cells.")
        ctx_pairs = list(combinations(contexts, 2))
        for c1, c2 in ctx_pairs:
            r1 = group_df[group_df[f"{context_prefix}_Level"] == c1].iloc[0]
            r2 = group_df[group_df[f"{context_prefix}_Level"] == c2].iloc[0]

            tgt_measure = float(r1[f"{target_prefix}_Measure"]) if np.isfinite(r1[f"{target_prefix}_Measure"]) else np.nan
            tgt_se = float(r1[f"{target_prefix}_SE"]) if np.isfinite(r1[f"{target_prefix}_SE"]) else np.nan
            tgt_index = r1.get(f"{target_prefix}_Index", np.nan)

            bias1 = float(r1["Bias Size"]) if np.isfinite(r1["Bias Size"]) else np.nan
            bias2 = float(r2["Bias Size"]) if np.isfinite(r2["Bias Size"]) else np.nan
            bias_se1 = float(r1["S.E."]) if np.isfinite(r1["S.E."]) else np.nan
            bias_se2 = float(r2["S.E."]) if np.isfinite(r2["S.E."]) else np.nan

            local1 = tgt_measure + bias1 / target_sign if np.isfinite(tgt_measure) and np.isfinite(bias1) else np.nan
            local2 = tgt_measure + bias2 / target_sign if np.isfinite(tgt_measure) and np.isfinite(bias2) else np.nan
            se1 = bias_se1 / abs(target_sign) if np.isfinite(bias1) and np.isfinite(bias_se1) and bias_se1 > 0 else np.nan
            se2 = bias_se2 / abs(target_sign) if np.isfinite(bias2) and np.isfinite(bias_se2) and bias_se2 > 0 else np.nan

            # Cancel the shared term algebraically before floating-point evaluation.
            contrast = (bias1 - bias2) / target_sign if np.isfinite(bias1) and np.isfinite(bias2) else np.nan
            se_contrast = float(np.hypot(se1, se2)) if np.isfinite(se1) and np.isfinite(se2) else np.nan
            standardized = contrast / se_contrast if np.isfinite(contrast) and np.isfinite(se_contrast) and se_contrast > 0 else np.nan

            rows.append({
                "Target": tgt_level,
                "Target N": tgt_index,
                "Target Measure": tgt_measure,
                "Target S.E.": tgt_se,
                "Target Sign": target_sign,
                "Context1": c1,
                "Context1 N": r1.get(f"{context_prefix}_Index", np.nan),
                "Local Measure1": local1,
                "SE1": se1,
                "Obs-Exp Avg1": float(r1["Obs-Exp Average"]),
                "Count1": float(r1["Observd Count"]),
                "Context2": c2,
                "Context2 N": r2.get(f"{context_prefix}_Index", np.nan),
                "Local Measure2": local2,
                "SE2": se2,
                "Obs-Exp Avg2": float(r2["Obs-Exp Average"]),
                "Count2": float(r2["Observd Count"]),
                "Contrast": contrast,
                "SE": se_contrast,
                "Standardized Contrast": standardized,
                "t": np.nan,
                "d.f.": np.nan,
                "Prob.": np.nan,
                "SE Method": "Conditional information for disjoint cells; fitted calibration held fixed",
                "SE Status": "conditional_approximation" if np.isfinite(se_contrast) else "not_available",
                "Inference Detail": (
                    "The shared target effect cancels. Zero covariance between cell shifts "
                    "is assumed only conditional on fixed calibration and independent "
                    "observations. Joint calibration uncertainty and reference-distribution "
                    "calibration are unavailable; the standardized contrast is descriptive."
                ),
                **qualification,
            })

    return pd.DataFrame(rows)


def adjust_p_values(p_values, method: str = "holm") -> np.ndarray:
    """Retain missing tests in family size, as R p.adjust with explicit n=family size."""
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or np.isinf(p).any() or np.any((p < 0) | (p > 1)):
        raise ValueError("p-values must be a one-dimensional vector in [0, 1] or NaN")
    method_l = str(method).lower()
    if method_l not in {"holm", "holm-bonferroni", "holm_bonferroni", "bh", "fdr", "benjamini-hochberg", "benjamini_hochberg", "none"}:
        raise ValueError(f"Unknown multiplicity method: {method}")
    out = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p)
    if not np.any(ok):
        return out
    vals = p[ok]
    order = np.argsort(vals)
    sorted_p = vals[order]
    m, n_observed = len(p), len(vals)
    if method_l in {"bh", "fdr", "benjamini-hochberg", "benjamini_hochberg"}:
        adjusted_sorted = sorted_p * m / np.arange(1, n_observed+1, dtype=float)
        adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
    elif method_l in {"holm", "holm-bonferroni", "holm_bonferroni"}:
        adjusted_sorted = np.maximum.accumulate(sorted_p * (m-np.arange(n_observed)))
    else:
        adjusted_sorted = sorted_p
    adjusted = np.empty(n_observed, dtype=float)
    adjusted[order] = np.clip(adjusted_sorted, 0., 1.)
    out[ok] = adjusted
    return out


def _bias_screening_settings(bundle, *, alpha=None, min_n=None, practical_logit=None) -> dict:
    """Resolve explicit overrides, retained pair settings, then legacy defaults."""
    saved = bundle.get("screening_settings", {}) if isinstance(bundle, dict) else {}
    settings = {
        "alpha": saved.get("alpha", 0.05) if alpha is None else alpha,
        "min_n": saved.get("min_n", 5) if min_n is None else min_n,
        "practical_logit": saved.get("practical_logit", 0.50) if practical_logit is None else practical_logit,
    }
    if (not 0 < settings["alpha"] < 1 or settings["min_n"] < 1
            or not np.isfinite(settings["practical_logit"]) or settings["practical_logit"] < 0):
        raise ValueError("Invalid bias screening thresholds")
    return settings


def build_dff_bias_screening_table(
    bias_results: dict | pd.DataFrame | None,
    *,
    alpha: float | None = None,
    min_n: int | None = None,
    practical_logit: float | None = None,
) -> pd.DataFrame:
    """Create a DFF/bias screening table with multiplicity and sparse-cell flags."""
    if isinstance(bias_results, pd.DataFrame):
        tbl = bias_results.copy()
    elif isinstance(bias_results, dict):
        tbl = bias_results.get("table", pd.DataFrame()).copy()
    else:
        tbl = pd.DataFrame()
    if tbl.empty or "Bias Size" not in tbl.columns:
        return pd.DataFrame(columns=[
            "FacetPair", "FacetA", "FacetA_Level", "FacetB", "FacetB_Level",
            "ObsN", "BiasSize", "SE", "t", "p", "p_holm", "p_bh",
            "AbsBias", "Direction", "SparseCell", "EvidenceLevel", "Flag",
            "ClaimStatus", "LRChiSq", "LRProb", "ProfileCILower",
            "ProfileCIUpper", "ProfileCIStatus", "LikelihoodBasis",
            "SEBasis", "InferenceTier", "MultiplicityBasis",
            "MultiplicityFamily", "ReportingCaution",
            "ScreenAlpha", "ScreenMinN", "ScreenPracticalLogit",
            "Interpretation", "NextAction",
        ])
    settings = _bias_screening_settings(bias_results, alpha=alpha, min_n=min_n, practical_logit=practical_logit)
    alpha, min_n, practical_logit = settings["alpha"], settings["min_n"], settings["practical_logit"]
    tbl = tbl.reset_index(drop=True)

    def col_series(name: str, default=np.nan) -> pd.Series:
        if name in tbl.columns:
            return tbl[name]
        return pd.Series([default] * len(tbl))

    # Old tables lack the scalar-fit/reference audit. Never silently promote them.
    reference_ok = (col_series("Primary P Method", "").eq("unit_weight_conditional_lr_v1")
                    & col_series("Reference Status", "").eq("unit_weight_conditional_lr_candidate")
                    & col_series("Cell Status", "").eq("interior"))
    p_vals = pd.to_numeric(col_series("LR Prob."), errors="coerce").where(reference_ok).to_numpy(dtype=float)
    reference_ok &= np.isfinite(p_vals)
    bias_vals = pd.to_numeric(col_series("Bias Size"), errors="coerce")
    se_vals = pd.to_numeric(col_series("S.E."), errors="coerce")
    t_vals = pd.to_numeric(col_series("t"), errors="coerce")
    obs_source = col_series("PositiveWeightObsN") if "PositiveWeightObsN" in tbl else col_series("ObsN")
    obs_n = pd.to_numeric(obs_source, errors="coerce")
    out = pd.DataFrame({
        "FacetPair": col_series("FacetA", "FacetA").astype(str)
        + " x "
        + col_series("FacetB", "FacetB").astype(str),
        "FacetA": col_series("FacetA", "").astype(str),
        "FacetA_Level": col_series("FacetA_Level", "").astype(str),
        "FacetB": col_series("FacetB", "").astype(str),
        "FacetB_Level": col_series("FacetB_Level", "").astype(str),
        "ObsN": obs_n,
        "BiasSize": bias_vals,
        "SE": se_vals,
        "t": t_vals,
        "p": p_vals,
        "CellStatus": col_series("Cell Status", "legacy_unverified"),
        "ReferenceStatus": col_series("Reference Status", "not_available"),
        "PrimaryPMethod": col_series("Primary P Method", "not_available"),
        "InferenceAvailable": reference_ok,
        "PersonN": col_series("PersonN"),
    })
    out["p_holm"], out["p_bh"], out["FamilySize"] = np.nan, np.nan, 0
    # One family per facet pair, fixed to all observed cells before exclusion or
    # failed fits. Reversing a pair's orientation must not split the family.
    family = out.apply(lambda row: tuple(sorted([row.FacetA, row.FacetB])), axis=1)
    out["MultiplicityFamily"] = family.map(lambda pair: " x ".join(pair))
    out["MultiplicityFamilyID"] = family.map(lambda pair: json.dumps(pair, ensure_ascii=False))
    cells = out.apply(lambda row: tuple(sorted([(row.FacetA, row.FacetA_Level),
                                                (row.FacetB, row.FacetB_Level)])), axis=1)
    if cells.duplicated().any():
        raise ValueError("Duplicate bias cells in the same facet-pair family")
    for _, idx_family in out.groupby(family, sort=False).groups.items():
        values = out.loc[idx_family, "p"].to_numpy(dtype=float)
        out.loc[idx_family, "p_holm"] = adjust_p_values(values, "holm")
        out.loc[idx_family, "p_bh"] = adjust_p_values(values, "bh")
        out.loc[idx_family, "FamilySize"] = len(idx_family)
    for name, value in _output_qualification.bias_cell_qualification().to_record().items():
        out[name] = value
    out["AbsBias"] = out["BiasSize"].abs()
    out["Direction"] = np.select(
        [out["BiasSize"] > 0, out["BiasSize"] < 0],
        ["Observed higher than expected", "Observed lower than expected"],
        default="No signed departure",
    )

    sparse = pd.to_numeric(out["ObsN"], errors="coerce").fillna(0) < int(min_n)
    stat_holm = pd.to_numeric(out["p_holm"], errors="coerce") < float(alpha)
    # Holm is the sole candidate statistical flag. BH remains descriptive:
    # dependence through estimated shared calibration has not been qualified.
    practical = pd.to_numeric(out["AbsBias"], errors="coerce") >= float(practical_logit)
    out["SparseCell"] = sparse.astype(bool)
    levels = np.select(
        [
            ~reference_ok,
            sparse,
            stat_holm & practical,
            stat_holm,
            practical,
        ],
        [
            "Inference unavailable",
            "Sparse cell",
            "Strong review",
            "Statistical review",
            "Practical review",
        ],
        default="No flag",
    )
    out["EvidenceLevel"] = levels
    out["Flag"] = out["EvidenceLevel"].ne("No flag")
    out["ClaimStatus"] = "Do not claim"
    out["LRChiSq"] = pd.to_numeric(col_series("LR ChiSq"), errors="coerce")
    out["LRProb"] = pd.to_numeric(col_series("LR Prob."), errors="coerce")
    out["ProfileCILower"] = pd.to_numeric(col_series("Profile CI Lower"), errors="coerce")
    out["ProfileCIUpper"] = pd.to_numeric(col_series("Profile CI Upper"), errors="coerce")
    out["ProfileCIStatus"] = col_series("Profile CI Status", "").astype(str)
    out["LikelihoodBasis"] = col_series("Likelihood Basis", "").astype(str)
    out["SEBasis"] = col_series("SE Basis", "").astype(str)
    out["ScreenAlpha"] = alpha
    out["ScreenMinN"] = min_n
    out["ScreenPracticalLogit"] = practical_logit
    out["InferenceTier"] = col_series("InferenceTier", "").astype(str)
    out["MultiplicityBasis"] = col_series("Multiplicity Basis", "").astype(str)
    out["MultiplicityFamily"] += (
        " observed-cell family; alpha=" + f"{float(alpha):.3g}, min_n={int(min_n)}, "
        + f"practical_logit={float(practical_logit):.3g}; failed/excluded cells retained"
    )
    out["ReportingCaution"] = np.select(
        [
            out["EvidenceLevel"].eq("Inference unavailable"),
            out["EvidenceLevel"].eq("Sparse cell"),
            out["Flag"].astype(bool),
        ],
        [
            "Reference unavailable: do not interpret missing evidence as no bias.",
            "Sparse cells are design evidence, not stable bias estimates.",
            "Treat as a conditional screening flag; corroborate with design, fit, and substantive evidence.",
        ],
        default=(
            "No flag under this conditional screen; this is not evidence that "
            "bias is absent outside the tested cell family."
        ),
    )
    out["Interpretation"] = np.select(
        [
            out["EvidenceLevel"].eq("Inference unavailable"),
            out["EvidenceLevel"].eq("Sparse cell"),
            out["EvidenceLevel"].eq("Strong review"),
            out["EvidenceLevel"].eq("Statistical review"),
            out["EvidenceLevel"].eq("Practical review"),
        ],
        [
            "Scalar fit or reference distribution unavailable; see CellStatus and ReferenceStatus.",
            "Too few observations for a stable DFF/bias claim.",
            "Cell is statistically flagged after Holm correction and practically large.",
            "Cell is statistically flagged after multiplicity screening.",
            "Cell is practically large but not multiplicity-adjusted significant.",
        ],
        default="No flag in this development calculation; error-rate and power validation pending.",
    )
    out["NextAction"] = np.select(
        [
            out["EvidenceLevel"].eq("Inference unavailable"),
            out["EvidenceLevel"].eq("Sparse cell"),
            out["EvidenceLevel"].eq("Strong review"),
            out["EvidenceLevel"].eq("Statistical review"),
            out["EvidenceLevel"].eq("Practical review"),
        ],
        [
            "Resolve the recorded fit/scope issue; retain this cell in the planned family.",
            "Add observations or combine substantively defensible levels before interpreting.",
            "Review scoring rubrics, rater notes, content match, and rerun sensitivity checks.",
            "Inspect heatmap pattern and corroborate with fit, residual, and substantive evidence.",
            "Decide whether this logit departure is substantively meaningful for the use case.",
        ],
        default="No special action beyond routine review.",
    )
    return out.sort_values(["Flag", "AbsBias"], ascending=[False, False]).reset_index(drop=True)


def summarize_dff_bias_screening(dff_tbl: pd.DataFrame) -> pd.DataFrame:
    """Compact counts for the DFF/bias screening panel and exports."""
    if dff_tbl is None or dff_tbl.empty:
        return pd.DataFrame(columns=["Metric", "Value", "Interpretation"])
    levels = dff_tbl["EvidenceLevel"].astype(str) if "EvidenceLevel" in dff_tbl.columns else pd.Series(dtype=str)
    abs_bias = pd.to_numeric(dff_tbl.get("AbsBias", pd.Series(dtype=float)), errors="coerce")
    p_holm = pd.to_numeric(dff_tbl.get("p_holm", pd.Series(dtype=float)), errors="coerce")
    p_bh = pd.to_numeric(dff_tbl.get("p_bh", pd.Series(dtype=float)), errors="coerce")
    profile_status = dff_tbl["ProfileCIStatus"].astype(str) if "ProfileCIStatus" in dff_tbl.columns else pd.Series(dtype=str)
    flagged = dff_tbl[dff_tbl.get("Flag", pd.Series([False] * len(dff_tbl))).astype(bool)]
    strongest = "none"
    if len(flagged) and "FacetA_Level" in flagged.columns and "FacetB_Level" in flagged.columns:
        row = flagged.iloc[0]
        strongest = (
            f"{row.get('FacetA', 'FacetA')}={row.get('FacetA_Level', '')} x "
            f"{row.get('FacetB', 'FacetB')}={row.get('FacetB_Level', '')}"
        )
    alpha = float(dff_tbl["ScreenAlpha"].iloc[0]) if "ScreenAlpha" in dff_tbl else 0.05
    alpha_label = f"{alpha:g}".lstrip("0")
    return pd.DataFrame([
        {"Metric": "Cells screened", "Value": int(len(dff_tbl)), "Interpretation": "Facet-level combinations in this pair."},
        {"Metric": "Flagged cells", "Value": int(len(flagged)), "Interpretation": "Sparse, statistical, or practical review prompts."},
        {"Metric": "Strong review", "Value": int((levels == "Strong review").sum()), "Interpretation": "Holm-adjusted and practically large cells."},
        {"Metric": f"Holm p < {alpha_label}", "Value": int((p_holm < alpha).sum()), "Interpretation": "Holm candidate flags; input reference distributions remain unqualified."},
        {"Metric": f"BH q < {alpha_label}", "Value": int((p_bh < alpha).sum()), "Interpretation": "Descriptive BH values; not used as the primary flag."},
        {"Metric": "Sparse cells", "Value": int((levels == "Sparse cell").sum()), "Interpretation": "Cells below the minimum observation threshold."},
        {"Metric": "Max |bias|", "Value": float(abs_bias.max()) if len(abs_bias.dropna()) else np.nan, "Interpretation": "Largest absolute logit departure."},
        {"Metric": "Profile CI ok", "Value": int((profile_status == "ok").sum()), "Interpretation": "Cells with both conditional LR-cutoff endpoints found; coverage unqualified."},
        {"Metric": "Inference tier", "Value": "conditional screening", "Interpretation": "Other fitted parameters are held fixed; global uncertainty is not propagated."},
        {"Metric": "Strongest flagged cell", "Value": strongest, "Interpretation": "First flagged row after sorting by |bias|."},
    ])


def _iter_bias_result_bundles(
    bias_results: dict | pd.DataFrame | None,
) -> list[tuple[str, dict]]:
    """Normalize selected-pair and all-pair bias containers."""
    if isinstance(bias_results, pd.DataFrame):
        return [("Bias table", {"table": bias_results})]
    if not isinstance(bias_results, dict) or not bias_results:
        return []
    if any(k in bias_results for k in ("table", "bias_tbl", "_skip_reason")):
        label = "Bias table"
        facet_a = bias_results.get("facet_a")
        facet_b = bias_results.get("facet_b")
        if facet_a or facet_b:
            label = f"{facet_a or 'FacetA'} x {facet_b or 'FacetB'}"
        return [(label, bias_results)]
    out: list[tuple[str, dict]] = []
    for pair_label, bundle in bias_results.items():
        if isinstance(bundle, pd.DataFrame):
            out.append((str(pair_label), {"table": bundle}))
        elif isinstance(bundle, dict):
            out.append((str(pair_label), bundle))
    return out


def _common_scale_bias_status(result: dict | None, diagnostics: dict | None) -> tuple[str, int | None]:
    """Evaluate whether the observed design is connected for bias claims."""
    if not isinstance(diagnostics, dict):
        return "Not evaluated: diagnostics are unavailable.", None
    obs = diagnostics.get("obs", pd.DataFrame())
    if not isinstance(obs, pd.DataFrame) or obs.empty:
        return "Not evaluated: observation-level diagnostics are unavailable.", None
    config = result.get("config", {}) if isinstance(result, dict) else {}
    configured = list(config.get("facet_names", []) or [])
    facet_cols = ["Person"] + configured if "Person" in obs.columns else configured
    facet_cols = [c for c in facet_cols if c in obs.columns]
    if len(facet_cols) < 2:
        return "Not evaluated: fewer than two observed facet columns are available.", None
    try:
        subsets = calc_subsets(obs, facet_cols)
        summary = subsets.get("summary", pd.DataFrame()) if isinstance(subsets, dict) else pd.DataFrame()
    except Exception as exc:
        return f"Not evaluated: connectedness check failed ({exc}).", None
    if not isinstance(summary, pd.DataFrame) or summary.empty:
        return "Not evaluated: connectedness summary is empty.", None
    n_subsets = int(summary["Subset"].nunique()) if "Subset" in summary.columns else int(len(summary))
    if n_subsets > 1:
        return (
            f"Review: observed design has {n_subsets} disconnected subset(s); "
            "bias cells across subsets are not on one common scale.",
            n_subsets,
        )
    return (
        "Ready: observed design is connected for the fitted Person + facet graph.",
        1,
    )


def build_bias_decision_stability_audit(
    bias_results: dict | pd.DataFrame | None,
    *,
    alpha: float | None = None,
    min_n: int | None = None,
    practical_logit: float | None = None,
    display_decimals: int = _decision_stability.BIAS_DISPLAY_DECIMALS,
) -> pd.DataFrame:
    """Return cell-level floating/display boundary evidence for bias screens."""

    audits: list[pd.DataFrame] = []
    for pair_label, bundle in _iter_bias_result_bundles(bias_results):
        if isinstance(bundle, dict) and "_skip_reason" in bundle and "table" not in bundle:
            continue
        settings = _bias_screening_settings(bundle, alpha=alpha, min_n=min_n, practical_logit=practical_logit)
        dff_tbl = build_dff_bias_screening_table(bundle, **settings)
        if dff_tbl.empty:
            continue
        audit = _decision_stability.audit_bias_decision_stability(
            dff_tbl,
            alpha=settings["alpha"],
            practical_logit=settings["practical_logit"],
            display_decimals=display_decimals,
        )
        if audit.empty:
            continue
        audit.insert(0, "FacetPairBundle", str(pair_label))
        audits.append(audit)
    return pd.concat(audits, ignore_index=True) if audits else pd.DataFrame()


def build_bias_inference_audit(
    bias_results: dict | pd.DataFrame | None,
    result: dict | None = None,
    diagnostics: dict | None = None,
    *,
    alpha: float | None = None,
    min_n: int | None = None,
    practical_logit: float | None = None,
) -> pd.DataFrame:
    """Pair-level guardrails for conditional bias / local-interaction claims."""
    columns = [
        "FacetPair", "Status", "ClaimStatus", "CellsScreened", "FlaggedCells",
        "StrongReviewCells", "SparseCells", "HolmSignificantCells",
        "BHSignificantCells", "PracticalCells", "MaxAbsBias",
        "BoundarySensitiveCells", "BoundarySensitiveDecisions",
        "BoundarySensitivitySummary",
        "CommonScaleSubsets", "CommonScaleStatus", "InferenceScope",
        "MultiplicityScope", "ProfileCIStatusSummary", "InferenceTierSummary",
        "UnavailableCells", "QualificationStatus", "PublicConclusionAllowed",
        "ScreenAlpha", "ScreenMinN", "ScreenPracticalLogit",
        "EvidenceSummary", "RecommendedAction",
    ]
    bundles = _iter_bias_result_bundles(bias_results)
    if not bundles:
        return pd.DataFrame(columns=columns)

    common_status, common_subsets = _common_scale_bias_status(result, diagnostics)
    disconnected = str(common_status).startswith("Review:")
    rows: list[dict[str, object]] = []
    inference_scope = (
        "Conditional screening: fitted theta, facet, step, slope, and other "
        "main-effect values are held fixed; global model uncertainty is not propagated."
    )
    multiplicity_scope = (
        "Holm is the sole candidate statistical flag; BH is descriptive. "
        "Each observed facet-pair cell family includes failed/excluded cells. "
        "No whole-analysis error control or confirmatory qualification is claimed."
    )
    for pair_label, bundle in bundles:
        if isinstance(bundle, dict) and "_skip_reason" in bundle and "table" not in bundle:
            reason = str(bundle.get("_skip_reason", "Bias screen was skipped."))
            rows.append({
                "FacetPair": pair_label,
                "Status": "Missing",
                "ClaimStatus": "Do not claim",
                "CellsScreened": 0,
                "FlaggedCells": 0,
                "StrongReviewCells": 0,
                "SparseCells": 0,
                "HolmSignificantCells": 0,
                "BHSignificantCells": 0,
                "PracticalCells": 0,
                "MaxAbsBias": np.nan,
                "BoundarySensitiveCells": 0,
                "BoundarySensitiveDecisions": 0,
                "BoundarySensitivitySummary": "No bias decision statistics were available for a boundary audit.",
                "CommonScaleSubsets": common_subsets,
                "CommonScaleStatus": common_status,
                "InferenceScope": inference_scope,
                "MultiplicityScope": multiplicity_scope,
                "ProfileCIStatusSummary": "",
                "InferenceTierSummary": "",
                "EvidenceSummary": reason,
                "RecommendedAction": "Run a selected-pair or all-pair bias screen before writing no-bias claims.",
            })
            continue

        settings = _bias_screening_settings(bundle, alpha=alpha, min_n=min_n, practical_logit=practical_logit)
        dff_tbl = build_dff_bias_screening_table(bundle, **settings)
        if dff_tbl.empty:
            rows.append({
                "FacetPair": pair_label,
                "Status": "Missing",
                "ClaimStatus": "Do not claim",
                "CellsScreened": 0,
                "FlaggedCells": 0,
                "StrongReviewCells": 0,
                "SparseCells": 0,
                "HolmSignificantCells": 0,
                "BHSignificantCells": 0,
                "PracticalCells": 0,
                "MaxAbsBias": np.nan,
                "BoundarySensitiveCells": 0,
                "BoundarySensitiveDecisions": 0,
                "BoundarySensitivitySummary": "No bias decision statistics were available for a boundary audit.",
                "CommonScaleSubsets": common_subsets,
                "CommonScaleStatus": common_status,
                "InferenceScope": inference_scope,
                "MultiplicityScope": multiplicity_scope,
                "ProfileCIStatusSummary": "",
                "InferenceTierSummary": "",
                "EvidenceSummary": "No interpretable bias cells were exported.",
                "RecommendedAction": "Check sparse/extreme data and rerun the bias screen if local interaction is a reporting target.",
            })
            continue

        levels = dff_tbl["EvidenceLevel"].astype(str) if "EvidenceLevel" in dff_tbl.columns else pd.Series(dtype=str)
        flags = dff_tbl["Flag"].astype(bool) if "Flag" in dff_tbl.columns else pd.Series([False] * len(dff_tbl))
        p_holm = pd.to_numeric(dff_tbl.get("p_holm", pd.Series(dtype=float)), errors="coerce")
        p_bh = pd.to_numeric(dff_tbl.get("p_bh", pd.Series(dtype=float)), errors="coerce")
        abs_bias = pd.to_numeric(dff_tbl.get("AbsBias", pd.Series(dtype=float)), errors="coerce")
        profile_status = (
            dff_tbl["ProfileCIStatus"].astype(str)
            if "ProfileCIStatus" in dff_tbl.columns else pd.Series(dtype=str)
        )
        tiers = (
            dff_tbl["InferenceTier"].astype(str)
            if "InferenceTier" in dff_tbl.columns else pd.Series(dtype=str)
        )
        flagged_n = int(flags.sum())
        strong_n = int((levels == "Strong review").sum())
        sparse_n = int(dff_tbl["SparseCell"].sum())
        unavailable_n = int((~dff_tbl["InferenceAvailable"]).sum())
        holm_n = int((p_holm < float(settings["alpha"])).sum())
        bh_n = int((p_bh < float(settings["alpha"])).sum())
        practical_n = int((abs_bias >= float(settings["practical_logit"])).sum())
        max_abs = float(abs_bias.max()) if len(abs_bias.dropna()) else np.nan
        boundary_audit = _decision_stability.audit_bias_decision_stability(
            dff_tbl,
            alpha=settings["alpha"],
            practical_logit=settings["practical_logit"],
            display_decimals=_decision_stability.BIAS_DISPLAY_DECIMALS,
        )
        boundary_counts = _decision_stability.summarize_boundary_audit(boundary_audit)
        if boundary_audit.empty:
            boundary_sensitive_cells = 0
        else:
            boundary_mask = (
                boundary_audit["BoundaryStatus"].astype(str).isin(
                    {"numerical_boundary", "display_rounding_boundary"}
                )
                | ~boundary_audit["DisplayDecisionConsistent"].fillna(True).astype(bool)
            )
            boundary_sensitive_cells = int(
                boundary_audit.loc[boundary_mask, "SourceRow"].nunique()
            )
        boundary_sensitive_decisions = int(
            boundary_counts["numerical_boundary"]
            + boundary_counts["display_rounding_boundary"]
        )
        profile_summary = (
            dict(profile_status[profile_status.ne("")].value_counts())
            if len(profile_status) else {}
        )
        tier_summary = dict(tiers[tiers.ne("")].value_counts()) if len(tiers) else {}

        status = "Review"
        claim_status = "Do not claim"

        if disconnected:
            action = "Resolve linking/connectivity or restrict bias wording to each connected subset."
        elif unavailable_n:
            action = "Resolve unavailable cells before interpreting the conditional screen; missing evidence is not no bias."
        elif sparse_n:
            action = "Treat sparse cells as design prompts; add observations or combine defensible levels before strong claims."
        elif boundary_sensitive_cells:
            action = (
                "Inspect raw p-values and |bias| values near the configured thresholds; "
                "report sensitivity and do not infer the decision from rounded output."
            )
        elif strong_n:
            action = "Review scoring rubrics, rater notes, content match, fit, and sensitivity before reporting the flagged cells."
        elif flagged_n:
            action = "Use flagged cells as conditional screening evidence and corroborate with substantive design review."
        else:
            action = "Report only that no cell was flagged in the computed conditional screen; avoid broad no-bias claims."

        rows.append({
            "FacetPair": pair_label,
            "Status": status,
            "ClaimStatus": claim_status,
            "CellsScreened": int(len(dff_tbl)),
            "UnavailableCells": unavailable_n,
            "ScreenAlpha": settings["alpha"], "ScreenMinN": settings["min_n"],
            "ScreenPracticalLogit": settings["practical_logit"],
            "QualificationStatus": "TECHNICAL_ONLY",
            "PublicConclusionAllowed": False,
            "FlaggedCells": flagged_n,
            "StrongReviewCells": strong_n,
            "SparseCells": sparse_n,
            "HolmSignificantCells": holm_n,
            "BHSignificantCells": bh_n,
            "PracticalCells": practical_n,
            "MaxAbsBias": max_abs,
            "BoundarySensitiveCells": boundary_sensitive_cells,
            "BoundarySensitiveDecisions": boundary_sensitive_decisions,
            "BoundarySensitivitySummary": (
                f"{boundary_counts['numerical_boundary']} numerical-boundary, "
                f"{boundary_counts['display_rounding_boundary']} display-rounding-boundary, "
                f"{boundary_counts['display_decision_mismatch']} raw/display decision mismatch(es)."
            ),
            "CommonScaleSubsets": common_subsets,
            "CommonScaleStatus": common_status,
            "InferenceScope": inference_scope,
            "MultiplicityScope": multiplicity_scope,
            "ProfileCIStatusSummary": profile_summary,
            "InferenceTierSummary": tier_summary,
            "EvidenceSummary": (
                f"{flagged_n} flagged of {len(dff_tbl)} screened cell(s); "
                f"{strong_n} strong-review, {sparse_n} sparse, {unavailable_n} unavailable, "
                f"{holm_n} Holm-significant, {bh_n} BH-significant; "
                f"{boundary_sensitive_cells} threshold-sensitive cell(s)."
            ),
            "RecommendedAction": action,
        })
    audit = pd.DataFrame(rows, columns=columns)
    if not audit.empty:
        audit["QualificationStatus"] = "TECHNICAL_ONLY"
        audit["PublicConclusionAllowed"] = False
    return audit



def _jmle_facet_boundary_hold(config) -> bool:
    return bool(
        isinstance(config, dict)
        and str(config.get("method", "")).upper() == "JMLE"
        and config.get("jmle_facet_boundary_review", False)
    )


def calc_subsets(obs_df, facet_cols):
    if obs_df is None or obs_df.empty or not facet_cols:
        return {"summary": pd.DataFrame(), "nodes": pd.DataFrame()}

    df = obs_df[facet_cols].dropna()
    if df.empty:
        return {"summary": pd.DataFrame(), "nodes": pd.DataFrame()}

    nodes = []
    for facet in facet_cols:
        nodes.extend([f"{facet}:{val}" for val in df[facet].astype(str).unique()])
    nodes = list(dict.fromkeys(nodes))
    if not nodes:
        return {"summary": pd.DataFrame(), "nodes": pd.DataFrame()}

    parent = {node: node for node in nodes}

    def find_root(x):
        px = parent.get(x)
        if px is None:
            return None
        if px != x:
            parent[x] = find_root(px)
        return parent[x]

    def union_nodes(a, b):
        ra = find_root(a)
        rb = find_root(b)
        if ra is None or rb is None or ra == rb:
            return
        parent[rb] = ra

    for _, row in df.iterrows():
        row_nodes = [f"{facet}:{row[facet]}" for facet in facet_cols if pd.notna(row[facet])]
        if len(row_nodes) < 2:
            continue
        base = row_nodes[0]
        for node in row_nodes[1:]:
            union_nodes(base, node)

    comp_ids = {node: find_root(node) for node in nodes}
    comp_levels = list(dict.fromkeys(comp_ids.values()))
    comp_index = {comp: idx + 1 for idx, comp in enumerate(comp_levels)}

    node_tbl = pd.DataFrame({
        "Node": nodes,
        "Component": [comp_ids[n] for n in nodes],
        "Subset": [comp_index[comp_ids[n]] for n in nodes],
        "Facet": [n.split(":", 1)[0] for n in nodes],
        "Level": [n.split(":", 1)[1] for n in nodes],
    })

    facet_counts = (
        node_tbl.groupby(["Subset", "Facet"])["Level"]
        .nunique()
        .reset_index()
        .pivot(index="Subset", columns="Facet", values="Level")
        .fillna(0)
        .reset_index()
    )

    first_facet = facet_cols[0]
    row_subset = df[first_facet].astype(str).apply(lambda v: comp_index[find_root(f"{first_facet}:{v}")])
    obs_counts = row_subset.value_counts().rename_axis("Subset").reset_index(name="Observations")

    summary_tbl = facet_counts.merge(obs_counts, on="Subset", how="left")
    summary_tbl["Observations"] = summary_tbl["Observations"].fillna(0).astype(int)
    summary_tbl = summary_tbl.sort_values("Observations", ascending=False)

    return {"summary": summary_tbl, "nodes": node_tbl}


def _bias_table_for_publication(
    all_bias_results: dict | None,
    model: str | None,
) -> "pd.DataFrame | None":
    """Export descriptive cells with the current G3 hold, including legacy tables."""
    parts = []
    for pair_key, bundle in (all_bias_results or {}).items():
        if not isinstance(bundle, dict):
            continue
        table = bundle.get("table")
        if not isinstance(table, pd.DataFrame) or table.empty:
            continue
        rows = table.copy()
        rows.insert(0, "Pair", str(pair_key))
        rows["Use"] = "Withheld"
        if "Cell Status" not in rows:
            rows["Cell Status"] = "legacy_unverified"
        parts.append(rows)
    if not parts:
        return None
    out = pd.concat(parts, ignore_index=True)
    out["Cell Status"] = out["Cell Status"].map({
        "interior": "Interior solution", "legacy_unverified": "Not verified",
        "excluded_extreme_level": "Extreme level excluded", "invalid_input": "Invalid input",
        "no_finite_mle": "No finite estimate", "search_range_limited": "Search limit reached",
        "solver_failed": "Solver failed", "stationarity_failed": "Stationarity not met",
        "likelihood_failed": "Likelihood check failed", "calculation_failed": "Calculation failed",
    }).fillna("Needs review")
    cols = ["Pair", "FacetA_Level", "FacetB_Level", "Bias Size", "ObsN", "Cell Status", "Use"]
    return out[[c for c in cols if c in out]].rename(columns={
        "FacetA_Level": "Level A", "FacetB_Level": "Level B", "ObsN": "Ratings",
        "Cell Status": "Calculation",
    }).copy()


def _bias_table_caption(model: str | None, table_idx: int | None = None) -> str:
    """Public exports preserve the G3 boundary regardless of saved-table version."""
    return (
        (f"Table {table_idx}. " if table_idx is not None else "")
        + f"Descriptive {str(model or '').upper()} interaction estimates. "
        "Bias Size is an additive eta shift with calibration held fixed. "
        "Validation of tests and intervals is pending: these values are retained only in "
        "technical exports. Use = Withheld: these cells do not establish bias or absence of bias."
    )
