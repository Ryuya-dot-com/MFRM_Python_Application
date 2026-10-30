"""Mathematical equivalence tests for P3-B numerical engine vectorization.

Every vectorized kernel in ``mfrm_app.engine``, ``mfrm_app.joint_interaction``,
and ``mfrm_app.mml_fixed_sd_review`` is checked against standalone scalar-loop
reference implementations and full statistical invariants across:
- RSM, PCM, and GPCM response models
- JMLE and MML (EM, direct L-BFGS-B, observed-information covariance, EAP/PV)
- Unweighted and non-unit observation weights
- Anchored/grouped constraints, population covariates, and ``include_log_sigma``
- Zero-observation (empty) person levels and missing cells
- Bias interactions, Krippendorff's alpha (ordinal/nominal/interval), and
  residual PCA LOO + row-bootstrap stability
- Structural intervals and joint fixed-effect interaction Wald/Holm inference
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.special import logsumexp
from scipy.stats import chi2, t as t_dist

from mfrm_app import engine, joint_interaction, mml_fixed_sd_review, structural_intervals


# ---------------------------------------------------------------------------
# Standalone scalar-loop reference oracles (pre-vectorization mathematical spec)
# ---------------------------------------------------------------------------

def _ref_loglik_rsm(eta, score_k, step_cum, weight=None):
    if eta.size == 0:
        return 0.0
    k_cat = len(step_cum)
    log_num = np.outer(eta, np.arange(k_cat)) - step_cum
    log_denom = logsumexp(log_num, axis=1)
    diff = log_num[np.arange(len(eta)), score_k] - log_denom
    return float(np.sum(diff if weight is None else diff * weight))


def _ref_loglik_pcm(eta, score_k, step_cum_mat, criterion_idx, weight=None):
    if eta.size == 0:
        return 0.0
    total = 0.0
    k_cat = step_cum_mat.shape[1]
    for c_idx in range(step_cum_mat.shape[0]):
        rows = np.where(criterion_idx == c_idx)[0]
        if rows.size == 0:
            continue
        log_num = np.outer(eta[rows], np.arange(k_cat)) - step_cum_mat[c_idx]
        log_denom = logsumexp(log_num, axis=1)
        diff = log_num[np.arange(len(rows)), score_k[rows]] - log_denom
        total += float(np.sum(diff if weight is None else diff * weight[rows]))
    return total


def _ref_loglik_gpcm(eta, score_k, step_cum_mat, step_idx, slopes, slope_idx, weight=None):
    if eta.size == 0:
        return 0.0
    if slope_idx is None:
        slope_idx = step_idx
    total = 0.0
    k_cat = step_cum_mat.shape[1]
    k_vals = np.arange(k_cat, dtype=float)
    slopes = np.asarray(slopes, dtype=float)
    for c_idx in range(step_cum_mat.shape[0]):
        rows = np.where(step_idx == c_idx)[0]
        if rows.size == 0:
            continue
        linear = np.outer(eta[rows], k_vals) - step_cum_mat[c_idx]
        log_num = slopes[slope_idx[rows], None] * linear
        log_denom = logsumexp(log_num, axis=1)
        diff = log_num[np.arange(len(rows)), score_k[rows]] - log_denom
        total += float(np.sum(diff if weight is None else diff * weight[rows]))
    return total


def _ref_e_step_and_mml_nll(par, idx, config, sizes, quad):
    params = engine.expand_params(par, sizes, config)
    base_eta = engine.compute_base_eta(idx, params, config)
    rows_by_person = idx["rows_by_person"]
    n_persons = len(rows_by_person)
    q_len = len(quad["nodes"])
    log_w = engine._log_mml_quadrature_weights(quad)
    log_lik_mat = np.zeros((n_persons, q_len), dtype=float)
    pop_mu = engine.compute_population_mu(params, config)
    nonempty_ll = []

    if config["model"] == "RSM":
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
    else:
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])

    for j, rows in enumerate(rows_by_person):
        if len(rows) == 0:
            continue
        base = base_eta[rows]
        score_k = idx["score_k"][rows]
        w = idx.get("weight")[rows] if idx.get("weight") is not None else None
        mu_j = pop_mu[j] if j < len(pop_mu) else 0.0
        for q, theta_q in enumerate(quad["nodes"]):
            if config["model"] == "RSM":
                ll_val = _ref_loglik_rsm(mu_j + theta_q + base, score_k, step_cum, weight=w)
            elif config["model"] == "PCM":
                ll_val = _ref_loglik_pcm(
                    mu_j + theta_q + base, score_k, step_cum_mat, idx["step_idx"][rows], weight=w
                )
            else:
                s_idx_full = idx.get("slope_idx")
                s_idx = s_idx_full[rows] if s_idx_full is not None else idx["step_idx"][rows]
                ll_val = _ref_loglik_gpcm(
                    mu_j + theta_q + base,
                    score_k,
                    step_cum_mat,
                    idx["step_idx"][rows],
                    params["slopes"],
                    s_idx,
                    weight=w,
                )
            log_lik_mat[j, q] = ll_val
        nonempty_ll.append(logsumexp(log_w + log_lik_mat[j]))

    log_joint = log_w[None, :] + log_lik_mat
    log_marginal = logsumexp(log_joint, axis=1)
    post_weights = np.exp(log_joint - log_marginal[:, None])
    marginal_ll = float(np.sum(log_marginal))
    mml_nll = -float(np.sum(nonempty_ll)) if len(idx["score_k"]) > 0 else 0.0
    return post_weights, marginal_ll, mml_nll


def _ref_m_step_value_grad(par, idx, config, sizes, quad, post_weights, *, include_log_sigma=False):
    params = engine.expand_params(par, sizes, config)
    base_eta = engine.compute_base_eta(idx, params, config)
    score_k = idx["score_k"]
    weight = np.ones_like(score_k, dtype=float) if idx.get("weight") is None else np.asarray(idx["weight"], dtype=float)
    k_cat = int(config["n_cat"])
    k_vals = np.arange(k_cat, dtype=float)
    n_steps = max(k_cat - 1, 0)
    thresholds = np.arange(n_steps)
    facet_signs = config.get("facet_signs", {})
    facet_full_grads = {
        facet: np.zeros(len(config["facet_specs"][facet]["levels"]), dtype=float)
        for facet in config["facet_names"]
    }
    pop = config.get("population_model", {})
    pop_x = np.asarray(pop.get("X", np.zeros((int(config.get("n_person", 0)), 0))), dtype=float)
    pop_mu = engine.compute_population_mu(params, config)
    pop_mu_by_obs = pop_mu[idx["person"]] if pop_mu.size else np.zeros(len(score_k), dtype=float)
    pop_grad_full = np.zeros(pop_x.shape[1], dtype=float) if pop.get("enabled") and pop_x.ndim == 2 else np.array([], dtype=float)

    if config["model"] == "RSM":
        step_grad_full = np.zeros(n_steps, dtype=float)
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
    else:
        step_facet_n = len(config["facet_levels"][config["step_facet"]])
        step_grad_full = np.zeros((step_facet_n, n_steps), dtype=float)
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])
    if config["model"] == "GPCM":
        slope_idx = idx.get("slope_idx") if idx.get("slope_idx") is not None else idx["step_idx"]
        slope_obs = params["slopes"][slope_idx]
        slope_levels = len(config["facet_levels"][config["slope_facet"]])
        log_slope_grad_full = np.zeros(slope_levels, dtype=float)
    else:
        slope_idx = None
        slope_obs = np.ones(len(score_k), dtype=float)
        log_slope_grad_full = np.array([], dtype=float)

    nll = 0.0
    log_sigma_grad = 0.0
    for q, theta_q in enumerate(quad["nodes"]):
        eff_weight = weight * post_weights[idx["person"], q]
        eta = theta_q + pop_mu_by_obs + base_eta
        if config["model"] == "RSM":
            log_num = np.outer(eta, k_vals) - step_cum
            log_denom = logsumexp(log_num, axis=1)
            ll_obs = log_num[np.arange(len(eta)), score_k] - log_denom
            probs = np.exp(log_num - log_denom[:, None])
            nll -= float(np.sum(eff_weight * ll_obs))
            if n_steps > 0:
                observed_gt = score_k[:, None] > thresholds[None, :]
                prob_gt = np.flip(np.cumsum(np.flip(probs[:, 1:], axis=1), axis=1), axis=1)
                step_grad_full += np.sum(eff_weight[:, None] * (observed_gt - prob_gt), axis=0)
        elif config["model"] == "PCM":
            probs = np.zeros((len(eta), k_cat), dtype=float)
            ll_obs = np.zeros(len(eta), dtype=float)
            for c_idx in range(step_cum_mat.shape[0]):
                rows = np.where(idx["step_idx"] == c_idx)[0]
                if rows.size == 0:
                    continue
                log_num = np.outer(eta[rows], k_vals) - step_cum_mat[c_idx]
                log_denom = logsumexp(log_num, axis=1)
                ll_obs[rows] = log_num[np.arange(len(rows)), score_k[rows]] - log_denom
                probs[rows, :] = np.exp(log_num - log_denom[:, None])
                if n_steps > 0:
                    observed_gt = score_k[rows, None] > thresholds[None, :]
                    prob_gt = np.flip(np.cumsum(np.flip(probs[rows, 1:], axis=1), axis=1), axis=1)
                    step_grad_full[c_idx, :] += np.sum(eff_weight[rows, None] * (observed_gt - prob_gt), axis=0)
            nll -= float(np.sum(eff_weight * ll_obs))
        else:
            probs = np.zeros((len(eta), k_cat), dtype=float)
            ll_obs = np.zeros(len(eta), dtype=float)
            linear_all = np.zeros((len(eta), k_cat), dtype=float)
            for c_idx in range(step_cum_mat.shape[0]):
                rows = np.where(idx["step_idx"] == c_idx)[0]
                if rows.size == 0:
                    continue
                linear = np.outer(eta[rows], k_vals) - step_cum_mat[c_idx]
                log_num = slope_obs[rows, None] * linear
                log_denom = logsumexp(log_num, axis=1)
                ll_obs[rows] = log_num[np.arange(len(rows)), score_k[rows]] - log_denom
                probs[rows, :] = np.exp(log_num - log_denom[:, None])
                linear_all[rows, :] = linear
                if n_steps > 0:
                    observed_gt = score_k[rows, None] > thresholds[None, :]
                    prob_gt = np.flip(np.cumsum(np.flip(probs[rows, 1:], axis=1), axis=1), axis=1)
                    step_grad_full[c_idx, :] += np.sum(
                        eff_weight[rows, None] * slope_obs[rows, None] * (observed_gt - prob_gt),
                        axis=0,
                    )
            obs_linear = linear_all[np.arange(len(eta)), score_k]
            expected_linear = np.sum(probs * linear_all, axis=1)
            log_slope_grad_full += np.bincount(
                slope_idx,
                weights=-eff_weight * slope_obs * (obs_linear - expected_linear),
                minlength=slope_levels,
            )
            nll -= float(np.sum(eff_weight * ll_obs))

        expected = probs @ k_vals
        eta_grad = -eff_weight * slope_obs * (score_k - expected)
        if include_log_sigma:
            log_sigma_grad += theta_q * float(np.sum(eta_grad))
        for facet in config["facet_names"]:
            sign = facet_signs.get(facet, -1)
            facet_full_grads[facet] += np.bincount(
                idx["facets"][facet],
                weights=eta_grad * sign,
                minlength=len(config["facet_specs"][facet]["levels"]),
            )
        if pop_grad_full.size:
            person_eta_grad = np.bincount(idx["person"], weights=eta_grad, minlength=pop_x.shape[0])
            pop_grad_full += pop_x.T @ person_eta_grad

    grad_parts = []
    if sizes.get("population", 0) > 0:
        grad_parts.append(pop_grad_full)
    for facet in config["facet_names"]:
        if sizes.get(facet, 0) > 0:
            grad_parts.append(engine.collapse_facet_gradient(facet_full_grads[facet], config["facet_specs"][facet]))
    if n_steps > 0:
        if config["model"] == "RSM":
            grad_parts.append(engine.collapse_sum_zero_gradient(step_grad_full))
        else:
            grad_parts.append(np.vstack([
                engine.collapse_sum_zero_gradient(row) for row in step_grad_full
            ]).reshape(-1))
    if config["model"] == "GPCM" and sizes.get("log_slopes", 0) > 0:
        grad_parts.append(engine.collapse_sum_zero_gradient(log_slope_grad_full))

    grad = np.concatenate([g for g in grad_parts if g.size]) if grad_parts else np.array([], dtype=float)
    penalty_value, penalty_grad = engine.facet_regularization_value_grad(par, sizes, config)
    if config.get("facet_regularization", {}).get("enabled"):
        nll += penalty_value
        grad = grad + penalty_grad
    if include_log_sigma:
        grad = np.append(grad, log_sigma_grad)
    return nll, grad


def _ref_jmle_value_grad(par, idx, config, sizes):
    params = engine.expand_params(par, sizes, config)
    eta = engine.compute_eta(idx, params, config)
    score_k = idx["score_k"]
    weight = np.ones_like(score_k, dtype=float) if idx.get("weight") is None else np.asarray(idx["weight"], dtype=float)
    k_cat = int(config["n_cat"])
    k_vals = np.arange(k_cat, dtype=float)
    slope_obs = np.ones(len(eta), dtype=float)
    gpcm_linear = None

    if config["model"] == "RSM":
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
        log_num = np.outer(eta, k_vals) - step_cum
        log_denom = logsumexp(log_num, axis=1)
        ll_obs = log_num[np.arange(len(eta)), score_k] - log_denom
        probs = np.exp(log_num - log_denom[:, None])
    elif config["model"] == "PCM":
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])
        probs = np.zeros((len(eta), k_cat), dtype=float)
        ll_obs = np.zeros(len(eta), dtype=float)
        for c_idx in range(step_cum_mat.shape[0]):
            rows = np.where(idx["step_idx"] == c_idx)[0]
            if rows.size == 0:
                continue
            log_num = np.outer(eta[rows], k_vals) - step_cum_mat[c_idx]
            log_denom = logsumexp(log_num, axis=1)
            ll_obs[rows] = log_num[np.arange(len(rows)), score_k[rows]] - log_denom
            probs[rows, :] = np.exp(log_num - log_denom[:, None])
    else:
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])
        slope_idx = idx.get("slope_idx") if idx.get("slope_idx") is not None else idx["step_idx"]
        slope_obs = params["slopes"][slope_idx]
        probs = np.zeros((len(eta), k_cat), dtype=float)
        ll_obs = np.zeros(len(eta), dtype=float)
        gpcm_linear = np.zeros((len(eta), k_cat), dtype=float)
        for c_idx in range(step_cum_mat.shape[0]):
            rows = np.where(idx["step_idx"] == c_idx)[0]
            if rows.size == 0:
                continue
            linear = np.outer(eta[rows], k_vals) - step_cum_mat[c_idx]
            log_num = slope_obs[rows, None] * linear
            log_denom = logsumexp(log_num, axis=1)
            ll_obs[rows] = log_num[np.arange(len(rows)), score_k[rows]] - log_denom
            probs[rows, :] = np.exp(log_num - log_denom[:, None])
            gpcm_linear[rows, :] = linear

    nll = -float(np.sum(weight * ll_obs))
    expected = probs @ k_vals
    eta_grad = -weight * slope_obs * (score_k - expected)

    grad_parts = []
    if sizes.get("theta", 0) > 0:
        theta_levels = len(config["theta_spec"]["levels"])
        theta_full_grad = np.bincount(idx["person"], weights=eta_grad, minlength=theta_levels)
        grad_parts.append(engine.collapse_facet_gradient(theta_full_grad, config["theta_spec"]))

    facet_signs = config.get("facet_signs", {})
    for facet in config["facet_names"]:
        if sizes.get(facet, 0) == 0:
            continue
        sign = facet_signs.get(facet, -1)
        facet_levels = len(config["facet_specs"][facet]["levels"])
        facet_full_grad = np.bincount(idx["facets"][facet], weights=eta_grad * sign, minlength=facet_levels)
        grad_parts.append(engine.collapse_facet_gradient(facet_full_grad, config["facet_specs"][facet]))

    n_steps = max(k_cat - 1, 0)
    if n_steps > 0:
        thresholds = np.arange(n_steps)
        if config["model"] == "RSM":
            observed_gt = score_k[:, None] > thresholds[None, :]
            prob_gt = np.flip(np.cumsum(np.flip(probs[:, 1:], axis=1), axis=1), axis=1)
            step_full_grad = np.sum(weight[:, None] * (observed_gt - prob_gt), axis=0)
            grad_parts.append(engine.collapse_sum_zero_gradient(step_full_grad))
        else:
            step_facet_n = len(config["facet_levels"][config["step_facet"]])
            step_mat_grad = np.zeros((step_facet_n, n_steps), dtype=float)
            row_factor = weight if config["model"] == "PCM" else (weight * slope_obs)
            for c_idx in range(step_facet_n):
                rows = np.where(idx["step_idx"] == c_idx)[0]
                if rows.size == 0:
                    continue
                observed_gt = score_k[rows, None] > thresholds[None, :]
                prob_gt = np.flip(np.cumsum(np.flip(probs[rows, 1:], axis=1), axis=1), axis=1)
                step_mat_grad[c_idx, :] = np.sum(row_factor[rows, None] * (observed_gt - prob_gt), axis=0)
            grad_parts.append(np.vstack([
                engine.collapse_sum_zero_gradient(row) for row in step_mat_grad
            ]).reshape(-1))

    if config["model"] == "GPCM" and sizes.get("log_slopes", 0) > 0:
        slope_idx = idx.get("slope_idx") if idx.get("slope_idx") is not None else idx["step_idx"]
        obs_linear = gpcm_linear[np.arange(len(eta)), score_k]
        expected_linear = np.sum(probs * gpcm_linear, axis=1)
        slope_levels = len(config["facet_levels"][config["slope_facet"]])
        slope_full_grad = np.bincount(
            slope_idx,
            weights=-weight * slope_obs * (obs_linear - expected_linear),
            minlength=slope_levels,
        )
        grad_parts.append(engine.collapse_sum_zero_gradient(slope_full_grad))

    grad = np.concatenate([g for g in grad_parts if g.size]) if grad_parts else np.array([], dtype=float)
    penalty_value, penalty_grad = engine.facet_regularization_value_grad(par, sizes, config)
    if config.get("facet_regularization", {}).get("enabled"):
        nll += penalty_value
        grad = grad + penalty_grad
    return nll, grad


def _ref_calc_bias_interactions(obs_df, facet_cols, pairs=None, top_n=20):
    if len(facet_cols) < 2:
        return pd.DataFrame()
    if pairs is None:
        combos = [
            (facet_cols[i], facet_cols[j])
            for i in range(len(facet_cols))
            for j in range(i + 1, len(facet_cols))
        ]
    elif len(pairs) == 0:
        return pd.DataFrame()
    else:
        combos = pairs
    rows = []
    for pair1, pair2 in combos:
        grp = obs_df.groupby([pair1, pair2], observed=True)
        for (lvl1, lvl2), df in grp:
            if df.empty:
                continue
            w = engine.get_weights(df)
            n = np.nansum(w)
            observed_avg = engine.weighted_mean(df["Observed"].to_numpy(), w)
            expected_avg = engine.weighted_mean(df["Expected"].to_numpy(), w)
            mean_resid = engine.weighted_mean(df["Residual"].to_numpy(), w)
            mean_std_resid = engine.weighted_mean(df["StdResidual"].to_numpy(), w)
            mean_abs_std = engine.weighted_mean(np.abs(df["StdResidual"]).to_numpy(), w)
            chi_sq = np.nansum((df["StdResidual"] ** 2) * w)
            se_resid = np.sqrt(np.nansum(df["Var"] * w)) / n if n > 0 else np.nan
            se_std = 1 / np.sqrt(n) if n > 0 else np.nan
            df_chi = n - 1 if n > 1 else np.nan
            t_resid = mean_resid / se_resid if np.isfinite(se_resid) and se_resid > 0 else np.nan
            t_std = mean_std_resid / se_std if np.isfinite(se_std) and se_std > 0 else np.nan
            p_resid = 2 * t_dist.cdf(-abs(t_resid), df=df_chi) if np.isfinite(df_chi) and np.isfinite(t_resid) else np.nan
            p_std = 2 * t_dist.cdf(-abs(t_std), df=df_chi) if np.isfinite(df_chi) and np.isfinite(t_std) else np.nan
            chi_p = 1 - chi2.cdf(chi_sq, df=df_chi) if np.isfinite(chi_sq) and np.isfinite(df_chi) and df_chi > 0 else np.nan
            rows.append({
                "Pair": f"{pair1} x {pair2}",
                "Level": f"{lvl1} | {lvl2}",
                "N": n,
                "ObservedAverage": observed_avg,
                "ExpectedAverage": expected_avg,
                "Bias": mean_resid,
                "MeanResidual": mean_resid,
                "MeanStdResidual": mean_std_resid,
                "MeanAbsStdResidual": mean_abs_std,
                "ChiSq": chi_sq,
                "ChiDf": df_chi,
                "ChiP": chi_p,
                "SE_Residual": se_resid,
                "t_Residual": t_resid,
                "p_Residual": p_resid,
                "SE_StdResidual": se_std,
                "t_StdResidual": t_std,
                "p_StdResidual": p_std,
                "DF": df_chi,
            })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["AbsStd"] = out["MeanStdResidual"].abs()
    return out.sort_values("AbsStd", ascending=False).drop(columns=["AbsStd"]).head(top_n)


def _ref_krippendorff_alpha_from_wide(wide: pd.DataFrame, metric: str):
    pairable_values: list[float] = []
    observed_num = 0.0
    pairable_units = 0
    for _, row in wide.iterrows():
        vals = pd.to_numeric(row, errors="coerce").dropna().to_numpy(dtype=float)
        m = len(vals)
        if m < 2:
            continue
        pairable_units += 1
        pairable_values.extend(float(v) for v in vals)
        if metric == "nominal":
            observed_num += float((vals[:, None] != vals[None, :]).sum()) / float(m - 1)
        elif metric == "interval":
            observed_num += float(((vals[:, None] - vals[None, :]) ** 2).sum()) / float(m - 1)
    values = np.asarray(pairable_values, dtype=float)
    ratings_used = int(len(values))
    if ratings_used < 2 or pairable_units < 1:
        return np.nan, np.nan, np.nan
    if metric == "ordinal":
        categories, counts = np.unique(values, return_counts=True)
        order = np.argsort(categories)
        categories = categories[order]
        counts = counts[order].astype(float)
        n_cat = len(categories)
        dist = np.zeros((n_cat, n_cat), dtype=float)
        for i in range(n_cat):
            for j in range(i + 1, n_cat):
                span = float(np.sum(counts[i:j + 1]) - 0.5 * (counts[i] + counts[j]))
                dist[i, j] = dist[j, i] = span ** 2
        cat_index = {float(cat): idx for idx, cat in enumerate(categories)}
        observed_num = 0.0
        for _, row in wide.iterrows():
            vals = pd.to_numeric(row, errors="coerce").dropna().to_numpy(dtype=float)
            if len(vals) < 2:
                continue
            c_idx = [cat_index[float(v)] for v in vals]
            observed_num += float(dist[np.ix_(c_idx, c_idx)].sum()) / float(len(vals) - 1)
        expected_num = float((counts[:, None] * counts[None, :] * dist).sum())
    elif metric == "nominal":
        _, counts = np.unique(values, return_counts=True)
        expected_num = float(values.size ** 2 - np.sum(counts.astype(float) ** 2))
    else:
        mean_val = float(np.mean(values))
        expected_num = float(2.0 * values.size * np.sum((values - mean_val) ** 2))
    do = float(observed_num / ratings_used)
    de = float(expected_num / (ratings_used * (ratings_used - 1)))
    alpha = float(1.0 - do / de) if de > 0 else np.nan
    return alpha, do, de


# ---------------------------------------------------------------------------
# Synthetic test fixtures covering edge cases
# ---------------------------------------------------------------------------

def _synthetic_rating_data(seed: int = 20261030, weighted: bool = False) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    persons = [f"P{i:02d}" for i in range(1, 23)]
    raters = ["R1", "R2", "R3", "R4"]
    tasks = ["T1", "T2", "T3"]
    rows = []
    theta = rng.normal(0.0, 1.0, size=len(persons))
    rater_sev = np.array([-0.55, -0.10, 0.20, 0.45])
    task_diff = np.array([-0.35, 0.05, 0.30])
    steps = np.array([-0.9, -0.2, 0.3, 0.8])
    step_cum = np.concatenate([[0.0], np.cumsum(steps)])
    k_vals = np.arange(len(step_cum), dtype=float)
    for p_i, p in enumerate(persons):
        for r_i, r in enumerate(raters):
            for t_i, t in enumerate(tasks):
                # Deterministic sparse missingness to exercise ragged person counts
                if (p_i + r_i + t_i) % 7 == 0:
                    continue
                eta = theta[p_i] - rater_sev[r_i] - task_diff[t_i]
                logits = eta * k_vals - step_cum
                probs = np.exp(logits - logsumexp(logits))
                score = int(rng.choice(len(step_cum), p=probs)) + 1
                row = {
                    "Person": p,
                    "Rater": r,
                    "Task": t,
                    "Score": score,
                    "Group": "G1" if p_i % 2 == 0 else "G2",
                }
                if weighted:
                    row["Weight"] = float(0.5 + 0.5 * ((p_i + r_i) % 4 + 1))
                rows.append(row)
    return pd.DataFrame(rows)


@pytest.mark.parametrize("model", ["RSM", "PCM", "GPCM"])
@pytest.mark.parametrize("weighted", [False, True])
def test_kernel_level_parity_jmle_and_mml_with_empty_person_and_log_sigma(model: str, weighted: bool):
    df = _synthetic_rating_data(seed=101 if not weighted else 202, weighted=weighted)
    step_facet = None if model == "RSM" else "Task"
    slope_facet = "Task" if model == "GPCM" else None

    # Fit quick baseline to obtain a realistic config/sizes/idx structure
    res = engine.mfrm_estimate(
        df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        weight_col="Weight" if weighted else None,
        model=model,
        method="MML",
        step_facet=step_facet,
        slope_facet=slope_facet,
        quad_points=15,
        maxit=12,
    )
    cfg = dict(res["config"])
    base_sizes = dict(engine.build_param_sizes(cfg))
    prep = engine.prepare_mfrm_data(
        df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        weight_col="Weight" if weighted else None,
    )
    # Inject an unobserved / empty person level at the end to verify empty-person semantics
    extended_person_levels = list(prep["levels"]["Person"]) + ["P_EMPTY"]
    prep["data"] = prep["data"].copy()
    prep["data"]["Person"] = pd.Categorical(prep["data"]["Person"], categories=extended_person_levels)
    prep["levels"] = dict(prep["levels"], Person=extended_person_levels)
    prep["n_person"] = len(extended_person_levels)

    idx = engine.build_indices(prep, step_facet=step_facet, slope_facet=slope_facet)
    cfg_mml = dict(cfg)
    cfg_mml["n_person"] = len(extended_person_levels)
    # Also enable population covariates on unweighted run to test pop_grad contraction
    if not weighted:
        rng_x = np.random.default_rng(77)
        x_mat = rng_x.normal(0.0, 1.0, size=(cfg_mml["n_person"], 2))
        x_mat -= x_mat.mean(axis=0, keepdims=True)
        cfg_mml["population_model"] = {"enabled": True, "X": x_mat, "columns": ["cov1", "cov2"]}
        sizes_mml = {"population": 2, **base_sizes}
    else:
        sizes_mml = dict(base_sizes)

    rng = np.random.default_rng(303)
    par_len_mml = sum(sizes_mml.values())
    par_mml = rng.normal(0.0, 0.35, size=par_len_mml)
    quad = engine.gauss_hermite_normal(17, sd=1.15)

    ref_post, ref_mll, ref_nll = _ref_e_step_and_mml_nll(par_mml, idx, cfg_mml, sizes_mml, quad)
    params_mml = engine.expand_params(par_mml, sizes_mml, cfg_mml)
    got_post, got_mll = engine._e_step_posteriors(idx, cfg_mml, params_mml, quad)
    got_nll = engine.mfrm_loglik_mml(par_mml, idx, cfg_mml, sizes_mml, quad)

    np.testing.assert_allclose(got_post, ref_post, rtol=1e-12, atol=1e-12)
    assert abs(got_mll - ref_mll) <= 1e-11
    assert abs(got_nll - ref_nll) <= 1e-11
    # Empty person posterior equals prior weights and sums to 1
    np.testing.assert_allclose(got_post[-1], quad["weights"] / np.sum(quad["weights"]), atol=1e-13)

    # Check M-step value/grad both with and without include_log_sigma, and at a DIFFERENT par from E-step
    par_mstep = par_mml + rng.normal(0.0, 0.12, size=par_len_mml)
    for inc_sigma in (False, True):
        ref_q_nll, ref_q_grad = _ref_m_step_value_grad(
            par_mstep, idx, cfg_mml, sizes_mml, quad, ref_post, include_log_sigma=inc_sigma
        )
        got_q_nll, got_q_grad = engine._m_step_expected_ll_value_grad(
            par_mstep, idx, cfg_mml, sizes_mml, quad, got_post, include_log_sigma=inc_sigma
        )
        assert abs(got_q_nll - ref_q_nll) <= 1e-11
        np.testing.assert_allclose(got_q_grad, ref_q_grad, rtol=1e-11, atol=1e-11)

        # Direct single-pass value + gradient at par_mml
        _, ref_dir_grad = _ref_m_step_value_grad(
            par_mml, idx, cfg_mml, sizes_mml, quad, ref_post, include_log_sigma=inc_sigma
        )
        pen_val, _ = engine.facet_regularization_value_grad(par_mml, sizes_mml, cfg_mml)
        ref_dir_val = -ref_mll + pen_val
        got_dir_val, got_dir_grad = engine.mfrm_loglik_mml_value_grad(
            par_mml, idx, cfg_mml, sizes_mml, quad, include_log_sigma=inc_sigma
        )
        assert abs(got_dir_val - ref_dir_val) <= 1e-11
        np.testing.assert_allclose(got_dir_grad, ref_dir_grad, rtol=1e-11, atol=1e-11)

    # Also verify large-tensor chunked fallback produces identical results when forced
    if hasattr(engine, "_MML_MAX_TENSOR_ELEMENTS"):
        old_limit = engine._MML_MAX_TENSOR_ELEMENTS
        try:
            engine._MML_MAX_TENSOR_ELEMENTS = 64  # force per-node 2D path
            chunk_post, chunk_mll = engine._e_step_posteriors(idx, cfg_mml, params_mml, quad)
            chunk_val, chunk_grad = engine.mfrm_loglik_mml_value_grad(
                par_mml, idx, cfg_mml, sizes_mml, quad, include_log_sigma=True
            )
            np.testing.assert_allclose(chunk_post, ref_post, rtol=1e-12, atol=1e-12)
            assert abs(chunk_mll - ref_mll) <= 1e-11
            assert abs(chunk_val - (-ref_mll)) <= 1e-11
            np.testing.assert_allclose(chunk_grad, got_dir_grad, rtol=1e-11, atol=1e-11)
        finally:
            engine._MML_MAX_TENSOR_ELEMENTS = old_limit

    # Verify JMLE likelihood, probabilities, and analytical gradient
    res_jmle = engine.mfrm_estimate(
        df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        weight_col="Weight" if weighted else None,
        model=model,
        method="JMLE",
        step_facet=step_facet,
        slope_facet=slope_facet,
        maxit=10,
    )
    prep_j = engine.prepare_mfrm_data(
        df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        weight_col="Weight" if weighted else None,
    )
    idx_j = engine.build_indices(prep_j, step_facet=step_facet, slope_facet=slope_facet)
    sizes_j = engine.build_param_sizes(res_jmle["config"])
    par_j = rng.normal(0.0, 0.3, size=sum(sizes_j.values()))
    ref_j_nll, ref_j_grad = _ref_jmle_value_grad(par_j, idx_j, res_jmle["config"], sizes_j)
    got_j_nll, got_j_grad = engine.mfrm_loglik_jmle_value_grad(par_j, idx_j, res_jmle["config"], sizes_j)
    got_j_only = engine.mfrm_loglik_jmle(par_j, idx_j, res_jmle["config"], sizes_j)
    assert abs(got_j_nll - ref_j_nll) <= 1e-11
    assert abs(got_j_only - ref_j_nll) <= 1e-11
    np.testing.assert_allclose(got_j_grad, ref_j_grad, rtol=1e-11, atol=1e-11)


def test_diagnostics_hotspots_parity_bias_krippendorff_pca_and_inference():
    df = _synthetic_rating_data(seed=404, weighted=False)
    res = engine.mfrm_estimate(
        df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        model="PCM",
        method="MML",
        step_facet="Task",
        quad_points=15,
    )
    diag = engine.mfrm_diagnostics(res, whexact=False, compute_pca=True)
    obs_df = diag["obs"]

    # 1. Bias interactions parity (including weighted & unobserved category levels)
    obs_cat = obs_df.copy()
    obs_cat["Rater"] = pd.Categorical(obs_cat["Rater"], categories=["R1", "R2", "R3", "R4", "R_UNOBS"])
    obs_cat["Task"] = pd.Categorical(obs_cat["Task"], categories=["T1", "T2", "T3"])
    ref_bias = _ref_calc_bias_interactions(obs_cat, ["Rater", "Task"], top_n=50)
    got_bias = engine.calc_bias_interactions(obs_cat, ["Rater", "Task"], top_n=50)
    pd.testing.assert_frame_equal(got_bias.reset_index(drop=True), ref_bias.reset_index(drop=True), rtol=1e-12, atol=1e-12)

    obs_weighted = obs_df.copy()
    obs_weighted["Weight"] = np.linspace(0.5, 2.0, len(obs_weighted))
    ref_bias_w = _ref_calc_bias_interactions(obs_weighted, ["Rater", "Task"], top_n=20)
    got_bias_w = engine.calc_bias_interactions(obs_weighted, ["Rater", "Task"], top_n=20)
    pd.testing.assert_frame_equal(got_bias_w.reset_index(drop=True), ref_bias_w.reset_index(drop=True), rtol=1e-12, atol=1e-12)

    # 2. Krippendorff's alpha parity across ordinal, nominal, interval (unweighted & weighted)
    for candidate_obs in (obs_df, obs_weighted):
        by_facet = engine.calc_krippendorff_alpha_by_facet(candidate_obs, ["Person", "Rater", "Task"])
        assert not by_facet.empty
        for facet in ("Rater", "Task"):
            wide, _ = engine._agreement_score_wide(candidate_obs, ["Person", "Rater", "Task"], facet)
            for metric in ("ordinal", "nominal", "interval"):
                ref_a, ref_do, ref_de = _ref_krippendorff_alpha_from_wide(wide, metric)
                single = engine.calc_krippendorff_alpha(candidate_obs, ["Person", "Rater", "Task"], facet, metric=metric)
                row = single["summary"].iloc[0]
                row_bf = by_facet[(by_facet["Facet"] == facet) & (by_facet["DistanceMetric"] == metric)].iloc[0]
                np.testing.assert_allclose(row["KrippendorffAlpha"], ref_a, rtol=1e-12, atol=1e-12)
                np.testing.assert_allclose(row["ObservedDisagreement"], ref_do, rtol=1e-12, atol=1e-12)
                np.testing.assert_allclose(row["ExpectedDisagreement"], ref_de, rtol=1e-12, atol=1e-12)
                np.testing.assert_allclose(row_bf["KrippendorffAlpha"], ref_a, rtol=1e-12, atol=1e-12)

    # 3. Residual PCA LOO & Bootstrap stability parity (with missing cells)
    pca_bundle = diag["pca"]
    assert pca_bundle is not None
    resid_clean = pca_bundle["residual_matrix"].dropna(axis=1, how="all")
    full_eig = pca_bundle["eigenvalues"]
    full_load = pca_bundle["loadings"]
    loo = engine._pca_leave_one_column_stability(resid_clean, full_eig, full_load)
    boot = engine._pca_row_bootstrap_stability(resid_clean, full_eig, full_load)

    # Compare against reference DataFrame-based LOO and bootstrap loop
    ref_ev1, ref_shares, ref_lcorrs = [], [], []
    cols = list(resid_clean.columns)[: engine.PCA_STABILITY_LOO_MAX_COLUMNS]
    for col in cols:
        reduced = resid_clean.drop(columns=[col], errors="ignore")
        core = engine._pca_core_from_residual_matrix(reduced)
        eig = np.asarray(core["eigenvalues"], dtype=float)
        ref_ev1.append(float(eig[0]))
        ref_shares.append(float(eig[0] / reduced.shape[1]))
        corr = engine._abs_loading_correlation(full_load["PC1"].drop(index=[col], errors="ignore"), core["loadings"]["PC1"])
        if np.isfinite(corr):
            ref_lcorrs.append(float(corr))
    np.testing.assert_allclose(loo["LOOEV1Min"], np.min(ref_ev1), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(loo["LOOEV1Max"], np.max(ref_ev1), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(loo["LOOPC1MedianAbsLoadingCorrelation"], np.median(ref_lcorrs), rtol=1e-12, atol=1e-12)

    work = resid_clean.iloc[:, : engine.PCA_STABILITY_BOOTSTRAP_MAX_COLUMNS]
    full_pc1 = full_load["PC1"].reindex(work.columns)
    rng_b = np.random.default_rng(engine.PCA_STABILITY_BOOTSTRAP_SEED)
    b_ev1, b_corrs = [], []
    for _ in range(engine.PCA_STABILITY_BOOTSTRAP_REPS):
        b_idx = rng_b.integers(0, work.shape[0], size=work.shape[0])
        core_b = engine._pca_core_from_residual_matrix(work.iloc[b_idx].reset_index(drop=True))
        if not core_b:
            continue
        b_ev1.append(float(core_b["eigenvalues"][0]))
        c_b = engine._abs_loading_correlation(full_pc1, core_b["loadings"]["PC1"])
        if np.isfinite(c_b):
            b_corrs.append(float(c_b))
    np.testing.assert_allclose(boot["BootstrapEV1Mean"], np.mean(b_ev1), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(boot["BootstrapEV1SD"], np.std(b_ev1, ddof=1), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(boot["BootstrapEV1P05"], np.percentile(b_ev1, 5), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(boot["BootstrapEV1P95"], np.percentile(b_ev1, 95), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(boot["BootstrapPC1MedianAbsLoadingCorrelation"], np.median(b_corrs), rtol=1e-11, atol=1e-11)
    assert boot["BootstrapReplicates"] == len(b_ev1)

    # 4. Fixed-SD structural intervals & joint interaction end-to-end checks
    prepared = mml_fixed_sd_review.prepare_fixed_sd_check(res)
    q_check = mml_fixed_sd_review.evaluate_fixed_sd_quadrature_check(prepared)
    assert isinstance(q_check, dict) and q_check
    interval_rec = structural_intervals.compute_structural_intervals(
        res,
        facet="Rater",
        level=0.95,
    )
    interval_df = structural_intervals.read_structural_intervals(res, interval_rec)
    assert not interval_df.empty
    assert {"Estimate", "SE", "Lower", "Upper"}.issubset(interval_df.columns)

    joint_rec = joint_interaction.fit_joint_interaction(
        res,
        facet_a="Rater",
        facet_b="Task",
        level=0.95,
    )
    joint_df = joint_interaction.read_joint_interaction(res, joint_rec)
    assert not joint_df.empty
    assert {"Estimate", "SE", "Lower", "Upper", "WaldP", "HolmP"}.issubset(joint_df.columns)
