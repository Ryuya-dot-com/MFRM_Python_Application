import io
from itertools import combinations
from collections import OrderedDict

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize, root_scalar, minimize_scalar
from scipy.special import logsumexp
from scipy.stats import chi2, t as t_dist


# ---- math helpers ----
def gauss_hermite_normal(n: int):
    if n < 1:
        raise ValueError("n must be >= 1")
    nodes, weights = hermgauss(n)  # exp(-x^2) weights
    nodes = np.sqrt(2.0) * nodes
    weights = weights / np.sqrt(np.pi)
    return {"nodes": nodes, "weights": weights}


def center_sum_zero(x: np.ndarray):
    if x.size == 0:
        return x
    return x - np.mean(x)


def build_facet_constraint(levels, anchors=None, groups=None, group_values=None, centered=True):
    lvl = [str(x) for x in levels]
    anchors_vec = np.full(len(lvl), np.nan, dtype=float)
    if anchors:
        for k, v in anchors.items():
            if k in lvl and pd.notna(v):
                anchors_vec[lvl.index(k)] = float(v)

    groups_vec = [None] * len(lvl)
    if groups:
        for k, v in groups.items():
            if k in lvl and pd.notna(v):
                groups_vec[lvl.index(k)] = str(v)

    group_values_map = {}
    if group_values:
        for k, v in group_values.items():
            if pd.notna(k):
                group_values_map[str(k)] = float(v)

    spec = {
        "levels": lvl,
        "anchors": anchors_vec,
        "groups": groups_vec,
        "group_values": group_values_map,
        "centered": bool(centered),
    }
    spec["n_params"] = count_facet_params(spec)
    return spec


def count_facet_params(spec):
    anchors = spec["anchors"]
    groups = spec["groups"]
    free_idx = np.where(np.isnan(anchors))[0]
    if free_idx.size == 0:
        return 0

    n_params = 0
    group_ids = {
        groups[i]
        for i in free_idx
        if groups[i] not in (None, "", np.nan)
    }
    for gid in group_ids:
        group_levels = [i for i, g in enumerate(groups) if g == gid]
        free_in_group = [i for i in group_levels if np.isnan(anchors[i])]
        k = len(free_in_group)
        if k > 1:
            n_params += k - 1

    ungrouped_idx = [
        i for i in free_idx if groups[i] in (None, "", np.nan)
    ]
    m = len(ungrouped_idx)
    if m == 0:
        return n_params
    if spec["centered"]:
        n_params += max(m - 1, 0)
    else:
        n_params += m
    return n_params


def expand_facet_with_constraints(free, spec):
    out = np.array(spec["anchors"], dtype=float, copy=True)
    groups = spec["groups"]
    group_values = spec["group_values"]
    centered = bool(spec["centered"])
    free_idx = np.where(np.isnan(out))[0]
    if free_idx.size == 0:
        return out

    used = 0
    group_ids = {
        groups[i]
        for i in free_idx
        if groups[i] not in (None, "", np.nan)
    }
    for gid in group_ids:
        group_levels = [i for i, g in enumerate(groups) if g == gid]
        free_in_group = [i for i in group_levels if np.isnan(out[i])]
        if not free_in_group:
            continue
        group_value = group_values.get(gid, 0.0)
        anchor_sum = np.nansum(out[group_levels])
        target_sum = group_value * len(group_levels)
        k = len(free_in_group)
        if k == 1:
            out[free_in_group[0]] = target_sum - anchor_sum
        else:
            seg = np.array(free[used : used + k - 1], dtype=float)
            used += k - 1
            last_val = target_sum - anchor_sum - np.sum(seg)
            out[free_in_group] = np.concatenate([seg, [last_val]])

    ungrouped_idx = [
        i for i in free_idx if groups[i] in (None, "", np.nan)
    ]
    m = len(ungrouped_idx)
    if m == 0:
        return out
    if centered:
        if m == 1:
            out[ungrouped_idx[0]] = 0.0
        else:
            seg = np.array(free[used : used + m - 1], dtype=float)
            used += m - 1
            out[ungrouped_idx] = np.concatenate([seg, [-np.sum(seg)]])
    else:
        seg = np.array(free[used : used + m], dtype=float)
        used += m
        out[ungrouped_idx] = seg
    return out


def build_param_sizes(config):
    n_steps = max(config["n_cat"] - 1, 0)
    sizes = OrderedDict()
    sizes["theta"] = config["theta_spec"]["n_params"] if config["method"] == "JMLE" else 0
    for facet in config["facet_names"]:
        sizes[facet] = config["facet_specs"][facet]["n_params"]
    if config["model"] == "RSM":
        sizes["steps"] = n_steps
    else:
        if not config.get("step_facet") or config["step_facet"] not in config["facet_names"]:
            raise ValueError("PCM requires a valid step facet.")
        sizes["steps"] = len(config["facet_levels"][config["step_facet"]]) * n_steps
    return sizes


def split_params(par, sizes):
    out = {}
    idx = 0
    for name, k in sizes.items():
        if k == 0:
            out[name] = np.array([], dtype=float)
        else:
            out[name] = np.array(par[idx : idx + k], dtype=float)
            idx += k
    return out


def expand_params(par, sizes, config):
    parts = split_params(par, sizes)
    theta = (
        expand_facet_with_constraints(parts["theta"], config["theta_spec"])
        if config["method"] == "JMLE"
        else np.array([], dtype=float)
    )

    facets = {}
    for facet in config["facet_names"]:
        facets[facet] = expand_facet_with_constraints(parts[facet], config["facet_specs"][facet])

    if config["model"] == "RSM":
        steps = center_sum_zero(parts["steps"])
        steps_mat = None
    else:
        n_levels = len(config["facet_levels"][config["step_facet"]])
        n_steps = max(config["n_cat"] - 1, 0)
        if n_levels == 0 or n_steps == 0:
            steps_mat = np.zeros((n_levels, n_steps))
        else:
            steps_mat = np.array(parts["steps"], dtype=float).reshape((n_levels, n_steps))
            steps_mat = np.vstack([center_sum_zero(row) for row in steps_mat])
        steps = None

    return {"theta": theta, "facets": facets, "steps": steps, "steps_mat": steps_mat}


# ---- data preparation ----
def prepare_mfrm_data(
    data,
    person_col,
    facet_cols,
    score_col,
    rating_min=None,
    rating_max=None,
    weight_col=None,
    keep_original=False,
):
    required = [person_col, score_col] + list(facet_cols)
    if weight_col:
        required.append(weight_col)
    if len(set(required)) != len(required):
        raise ValueError("Person/score/facet columns must be distinct (no duplicates).")
    if not all(col in data.columns for col in required):
        raise ValueError("Some selected columns are not in the data.")
    if data.columns.duplicated().any():
        dupes = data.columns[data.columns.duplicated()].tolist()
        if any(col in dupes for col in required):
            raise ValueError("Selected columns include duplicate names in the data. Please rename columns to be unique.")
    if len(facet_cols) == 0:
        raise ValueError("Select at least one facet column.")

    cols = [person_col] + list(facet_cols) + [score_col]
    if weight_col:
        cols.append(weight_col)
    df = data[cols].copy()
    rename_map = {person_col: "Person", score_col: "Score"}
    if weight_col:
        rename_map[weight_col] = "Weight"
    df.rename(columns=rename_map, inplace=True)
    df["Person"] = df["Person"].astype(str)
    for facet in facet_cols:
        df[facet] = df[facet].astype(str)
    df["Score"] = pd.to_numeric(df["Score"], errors="coerce")
    if "Weight" in df.columns:
        df["Weight"] = pd.to_numeric(df["Weight"], errors="coerce")
    else:
        df["Weight"] = 1.0

    df = df.dropna()
    df = df[df["Weight"] > 0]
    df["Score"] = df["Score"].astype(int)

    if rating_min is None:
        rating_min = int(df["Score"].min())
    if rating_max is None:
        rating_max = int(df["Score"].max())

    if not keep_original:
        score_vals = np.sort(df["Score"].unique())
        expected_vals = np.arange(rating_min, rating_max + 1)
        if not np.array_equal(score_vals, expected_vals):
            mapping = {val: rating_min + i for i, val in enumerate(score_vals)}
            df["Score"] = df["Score"].map(mapping).astype(int)
            rating_max = rating_min + len(score_vals) - 1

    df["score_k"] = df["Score"] - rating_min

    df["Person"] = pd.Categorical(df["Person"])
    for facet in facet_cols:
        df[facet] = pd.Categorical(df[facet])

    facet_levels = {facet: list(df[facet].cat.categories) for facet in facet_cols}

    return {
        "data": df,
        "rating_min": rating_min,
        "rating_max": rating_max,
        "facet_names": list(facet_cols),
        "levels": {"Person": list(df["Person"].cat.categories), **facet_levels},
        "weight_col": "Weight",
    }


def build_indices(prep, step_facet=None):
    df = prep["data"]
    facets_idx = {facet: df[facet].cat.codes.to_numpy(dtype=int) for facet in prep["facet_names"]}
    step_idx = df[step_facet].cat.codes.to_numpy(dtype=int) if step_facet else None
    return {
        "person": df["Person"].cat.codes.to_numpy(dtype=int),
        "facets": facets_idx,
        "step_idx": step_idx,
        "score_k": df["score_k"].to_numpy(dtype=int),
        "weight": df["Weight"].to_numpy(dtype=float) if "Weight" in df.columns else np.ones(len(df)),
    }


def sample_mfrm_data(seed=20240131):
    rng = np.random.default_rng(seed)
    persons = [f"P{idx:02d}" for idx in range(1, 37)]
    raters = [f"R{idx}" for idx in range(1, 4)]
    tasks = [f"T{idx}" for idx in range(1, 5)]
    criteria = [f"C{idx}" for idx in range(1, 4)]

    rows = []
    for p in persons:
        for r in raters:
            for t in tasks:
                for c in criteria:
                    rows.append((p, r, t, c))
    df = pd.DataFrame(rows, columns=["Person", "Rater", "Task", "Criterion"])

    ability = rng.normal(0, 1, len(persons))
    rater_eff = np.array([-0.4, 0.0, 0.4])
    task_eff = np.linspace(-0.5, 0.5, len(tasks))
    crit_eff = np.array([-0.3, 0.0, 0.3])

    eta = (
        ability[[persons.index(p) for p in df["Person"]]]
        - rater_eff[[raters.index(r) for r in df["Rater"]]]
        - task_eff[[tasks.index(t) for t in df["Task"]]]
        - crit_eff[[criteria.index(c) for c in df["Criterion"]]]
    )
    raw = eta + rng.normal(0, 0.6, len(df))
    bins = [-np.inf, -1.0, -0.3, 0.3, 1.0, np.inf]
    score = pd.cut(raw, bins=bins, labels=[1, 2, 3, 4, 5]).astype(int)
    df["Score"] = score
    return df


def format_tab_template(df: pd.DataFrame):
    char_df = df.fillna("").astype(str)
    widths = [max([len(str(c)) for c in [col] + char_df[col].tolist()]) for col in char_df.columns]

    def format_row(row_vec):
        return "\t".join(str(val).ljust(width) for val, width in zip(row_vec, widths))

    header = format_row(list(char_df.columns))
    rows = [format_row(row) for row in char_df.to_numpy()] if len(char_df) else []
    return "\n".join([header] + rows)


def guess_col(cols, patterns, fallback=0):
    if not cols:
        return None
    lowered = [c.lower() for c in cols]
    for pattern in patterns:
        for idx, name in enumerate(lowered):
            if pattern in name:
                return cols[idx]
    return cols[min(fallback, len(cols) - 1)]


def truncate_label(value, width=28):
    text = str(value)
    if len(text) <= width:
        return text
    return text[: max(0, width - 1)] + "…"


# ---- likelihoods ----
def loglik_rsm(eta, score_k, step_cum, weight=None):
    if eta.size == 0:
        return 0.0
    k_cat = len(step_cum)
    eta_mat = np.outer(eta, np.arange(k_cat))
    log_num = eta_mat - step_cum
    log_denom = logsumexp(log_num, axis=1)
    log_num_obs = log_num[np.arange(len(eta)), score_k]
    diff = log_num_obs - log_denom
    if weight is None:
        return float(np.sum(diff))
    return float(np.sum(diff * weight))


def loglik_pcm(eta, score_k, step_cum_mat, criterion_idx, weight=None):
    if eta.size == 0:
        return 0.0
    total = 0.0
    k_cat = step_cum_mat.shape[1]
    for c_idx in range(step_cum_mat.shape[0]):
        rows = np.where(criterion_idx == c_idx)[0]
        if rows.size == 0:
            continue
        eta_c = eta[rows]
        step_cum = step_cum_mat[c_idx]
        eta_mat = np.outer(eta_c, np.arange(k_cat))
        log_num = eta_mat - step_cum
        log_denom = logsumexp(log_num, axis=1)
        log_num_obs = log_num[np.arange(len(rows)), score_k[rows]]
        diff = log_num_obs - log_denom
        if weight is None:
            total += float(np.sum(diff))
        else:
            total += float(np.sum(diff * weight[rows]))
    return total


def category_prob_rsm(eta, step_cum):
    if eta.size == 0:
        return np.zeros((0, len(step_cum)))
    k_cat = len(step_cum)
    eta_mat = np.outer(eta, np.arange(k_cat))
    log_num = eta_mat - step_cum
    log_denom = logsumexp(log_num, axis=1)
    return np.exp(log_num - log_denom[:, None])


def category_prob_pcm(eta, step_cum_mat, criterion_idx):
    if eta.size == 0:
        return np.zeros((0, step_cum_mat.shape[1]))
    k_cat = step_cum_mat.shape[1]
    probs = np.zeros((len(eta), k_cat))
    for c_idx in range(step_cum_mat.shape[0]):
        rows = np.where(criterion_idx == c_idx)[0]
        if rows.size == 0:
            continue
        step_cum = step_cum_mat[c_idx]
        eta_c = eta[rows]
        eta_mat = np.outer(eta_c, np.arange(k_cat))
        log_num = eta_mat - step_cum
        log_denom = logsumexp(log_num, axis=1)
        probs[rows, :] = np.exp(log_num - log_denom[:, None])
    return probs


def zstd_from_mnsq(mnsq, df, whexact=False):
    if not np.isfinite(mnsq) or not np.isfinite(df) or df <= 0:
        return np.nan
    if whexact:
        return (mnsq - 1) * np.sqrt(df / 2)
    return (mnsq ** (1 / 3) - (1 - 2 / (9 * df))) / np.sqrt(2 / (9 * df))


def compute_base_eta(idx, params, config):
    eta = np.zeros_like(idx["score_k"], dtype=float)
    facet_signs = config.get("facet_signs", {})
    for facet in config["facet_names"]:
        sign = facet_signs.get(facet, -1)
        eta += sign * params["facets"][facet][idx["facets"][facet]]
    return eta


def compute_eta(idx, params, config, theta_override=None):
    theta = theta_override if theta_override is not None else params["theta"]
    eta = theta[idx["person"]] if theta.size else np.zeros_like(idx["score_k"], dtype=float)
    return eta + compute_base_eta(idx, params, config)


def mfrm_loglik_jmle(par, idx, config, sizes):
    params = expand_params(par, sizes, config)
    eta = compute_eta(idx, params, config)
    if config["model"] == "RSM":
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
        ll = loglik_rsm(eta, idx["score_k"], step_cum, weight=idx.get("weight"))
    else:
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])
        ll = loglik_pcm(eta, idx["score_k"], step_cum_mat, idx["step_idx"], weight=idx.get("weight"))
    return -ll


def mfrm_loglik_mml(par, idx, config, sizes, quad):
    params = expand_params(par, sizes, config)
    n = len(idx["score_k"])
    if n == 0:
        return 0.0

    base_eta = compute_base_eta(idx, params, config)
    rows_by_person = {}
    for i, p_idx in enumerate(idx["person"]):
        rows_by_person.setdefault(p_idx, []).append(i)

    log_w = np.log(quad["weights"])

    if config["model"] == "RSM":
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
        ll_person = []
        for rows in rows_by_person.values():
            base = base_eta[rows]
            score_k = idx["score_k"][rows]
            w = idx.get("weight")[rows] if idx.get("weight") is not None else None
            ll_nodes = [loglik_rsm(theta + base, score_k, step_cum, weight=w) for theta in quad["nodes"]]
            ll_person.append(logsumexp(log_w + ll_nodes))
    else:
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])
        ll_person = []
        for rows in rows_by_person.values():
            base = base_eta[rows]
            score_k = idx["score_k"][rows]
            crit = idx["step_idx"][rows]
            w = idx.get("weight")[rows] if idx.get("weight") is not None else None
            ll_nodes = [loglik_pcm(theta + base, score_k, step_cum_mat, crit, weight=w) for theta in quad["nodes"]]
            ll_person.append(logsumexp(log_w + ll_nodes))

    return -float(np.sum(ll_person))


def compute_person_eap(idx, config, params, quad):
    n = len(idx["score_k"])
    if n == 0:
        return pd.DataFrame(columns=["Person", "Estimate", "SD"])

    base_eta = compute_base_eta(idx, params, config)
    rows_by_person = {}
    for i, p_idx in enumerate(idx["person"]):
        rows_by_person.setdefault(p_idx, []).append(i)

    log_w = np.log(quad["weights"])
    estimates = []
    if config["model"] == "RSM":
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
        for rows in rows_by_person.values():
            base = base_eta[rows]
            score_k = idx["score_k"][rows]
            w = idx.get("weight")[rows] if idx.get("weight") is not None else None
            ll_nodes = np.array([loglik_rsm(theta + base, score_k, step_cum, weight=w) for theta in quad["nodes"]])
            log_post = log_w + ll_nodes
            log_post = log_post - logsumexp(log_post)
            post_w = np.exp(log_post)
            eap = float(np.sum(quad["nodes"] * post_w))
            sd = float(np.sqrt(np.sum((quad["nodes"] - eap) ** 2 * post_w)))
            estimates.append((eap, sd))
    else:
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])
        for rows in rows_by_person.values():
            base = base_eta[rows]
            score_k = idx["score_k"][rows]
            crit = idx["step_idx"][rows]
            w = idx.get("weight")[rows] if idx.get("weight") is not None else None
            ll_nodes = np.array([loglik_pcm(theta + base, score_k, step_cum_mat, crit, weight=w) for theta in quad["nodes"]])
            log_post = log_w + ll_nodes
            log_post = log_post - logsumexp(log_post)
            post_w = np.exp(log_post)
            eap = float(np.sum(quad["nodes"] * post_w))
            sd = float(np.sqrt(np.sum((quad["nodes"] - eap) ** 2 * post_w)))
            estimates.append((eap, sd))

    est_mat = np.array(estimates)
    return pd.DataFrame({"Estimate": est_mat[:, 0], "SD": est_mat[:, 1]})


def normalize_anchor_df(df):
    if df is None or df.empty:
        return pd.DataFrame()
    nm = [c.lower() for c in df.columns]
    facet_col = next((i for i, c in enumerate(nm) if c in ("facet", "facets")), None)
    level_col = next((i for i, c in enumerate(nm) if c in ("level", "element", "label")), None)
    anchor_col = next((i for i, c in enumerate(nm) if c in ("anchor", "value", "measure")), None)
    if facet_col is None or level_col is None or anchor_col is None:
        return pd.DataFrame()
    out = pd.DataFrame({
        "Facet": df.iloc[:, facet_col].astype(str),
        "Level": df.iloc[:, level_col].astype(str),
        "Anchor": pd.to_numeric(df.iloc[:, anchor_col], errors="coerce"),
    })
    return out.dropna(subset=["Facet", "Level"])


def normalize_group_anchor_df(df):
    if df is None or df.empty:
        return pd.DataFrame()
    nm = [c.lower() for c in df.columns]
    facet_col = next((i for i, c in enumerate(nm) if c in ("facet", "facets")), None)
    level_col = next((i for i, c in enumerate(nm) if c in ("level", "element", "label")), None)
    group_col = next((i for i, c in enumerate(nm) if c in ("group", "subset")), None)
    value_col = next((i for i, c in enumerate(nm) if c in ("groupvalue", "value", "anchor")), None)
    if facet_col is None or level_col is None or group_col is None or value_col is None:
        return pd.DataFrame()
    out = pd.DataFrame({
        "Facet": df.iloc[:, facet_col].astype(str),
        "Level": df.iloc[:, level_col].astype(str),
        "Group": df.iloc[:, group_col].astype(str),
        "GroupValue": pd.to_numeric(df.iloc[:, value_col], errors="coerce"),
    })
    return out.dropna(subset=["Facet", "Level", "Group"])


def prepare_constraint_specs(prep, anchor_df=None, group_anchor_df=None, noncenter_facet="Person", dummy_facets=None):
    facet_names = prep["facet_names"]
    all_facets = ["Person"] + facet_names

    anchor_df = normalize_anchor_df(anchor_df) if anchor_df is not None else pd.DataFrame()
    anchor_df = anchor_df[anchor_df["Facet"].isin(all_facets)] if not anchor_df.empty else anchor_df

    group_anchor_df = normalize_group_anchor_df(group_anchor_df) if group_anchor_df is not None else pd.DataFrame()
    group_anchor_df = group_anchor_df[group_anchor_df["Facet"].isin(all_facets)] if not group_anchor_df.empty else group_anchor_df

    dummy_facets = set(dummy_facets or [])

    anchor_map = {facet: {} for facet in all_facets}
    group_map = {facet: {} for facet in all_facets}
    group_values = {facet: {} for facet in all_facets}

    for facet in all_facets:
        levels = prep["levels"].get(facet, [])
        if facet in dummy_facets:
            anchor_map[facet] = {level: 0.0 for level in levels}
            group_map[facet] = {}
            group_values[facet] = {}
            continue
        if not anchor_df.empty:
            df = anchor_df[anchor_df["Facet"] == facet]
            if not df.empty:
                anchors = {row.Level: row.Anchor for row in df.itertuples() if row.Level in levels}
                anchor_map[facet] = anchors
        if not group_anchor_df.empty:
            df = group_anchor_df[group_anchor_df["Facet"] == facet]
            if not df.empty:
                groups = {row.Level: row.Group for row in df.itertuples() if row.Level in levels}
                group_map[facet] = groups
                group_vals = (
                    df[["Group", "GroupValue"]]
                    .dropna(subset=["Group"])
                    .drop_duplicates(subset=["Group"])
                    .assign(GroupValue=lambda x: x["GroupValue"].fillna(0))
                )
                group_values[facet] = {row.Group: row.GroupValue for row in group_vals.itertuples()}

    theta_spec = build_facet_constraint(
        levels=prep["levels"]["Person"],
        anchors=anchor_map["Person"],
        groups=group_map["Person"],
        group_values=group_values["Person"],
        centered=noncenter_facet != "Person",
    )

    facet_specs = {}
    for facet in facet_names:
        facet_specs[facet] = build_facet_constraint(
            levels=prep["levels"][facet],
            anchors=anchor_map[facet],
            groups=group_map[facet],
            group_values=group_values[facet],
            centered=noncenter_facet != facet,
        )

    anchor_summary = pd.DataFrame({
        "Facet": all_facets,
        "AnchoredLevels": [len(anchor_map[f]) for f in all_facets],
        "GroupAnchors": [len(set(group_map[f].values())) if group_map[f] else 0 for f in all_facets],
        "DummyFacet": [f in dummy_facets for f in all_facets],
    })

    return {
        "theta_spec": theta_spec,
        "facet_specs": facet_specs,
        "anchor_summary": anchor_summary,
    }


def mfrm_estimate(
    data,
    person_col,
    facet_cols,
    score_col,
    rating_min=None,
    rating_max=None,
    weight_col=None,
    keep_original=False,
    model="RSM",
    method="JMLE",
    step_facet=None,
    anchor_df=None,
    group_anchor_df=None,
    noncenter_facet="Person",
    dummy_facets=None,
    positive_facets=None,
    quad_points=15,
    maxit=400,
    reltol=1e-6,
):
    prep = prepare_mfrm_data(
        data,
        person_col=person_col,
        facet_cols=facet_cols,
        score_col=score_col,
        rating_min=rating_min,
        rating_max=rating_max,
        weight_col=weight_col,
        keep_original=keep_original,
    )

    if model == "PCM":
        if step_facet is None:
            step_facet = prep["facet_names"][0]
        if step_facet not in prep["facet_names"]:
            raise ValueError("Selected step facet is not in the facet list.")
    else:
        step_facet = None

    idx = build_indices(prep, step_facet=step_facet)

    n_person = len(prep["levels"]["Person"])
    facet_levels = {f: prep["levels"][f] for f in prep["facet_names"]}
    n_cat = prep["rating_max"] - prep["rating_min"] + 1

    if noncenter_facet not in ["Person"] + prep["facet_names"]:
        noncenter_facet = "Person"

    dummy_facets = [f for f in (dummy_facets or []) if f in ["Person"] + prep["facet_names"]]

    positive_facets = set(positive_facets or [])
    facet_signs = {facet: (1 if facet in positive_facets else -1) for facet in prep["facet_names"]}

    config = {
        "model": model,
        "method": method,
        "n_person": n_person,
        "n_cat": n_cat,
        "facet_names": prep["facet_names"],
        "facet_levels": facet_levels,
        "step_facet": step_facet,
    }
    config["weight_col"] = weight_col if weight_col else None
    config["positive_facets"] = list(positive_facets)
    config["facet_signs"] = facet_signs

    constraint_specs = prepare_constraint_specs(
        prep=prep,
        anchor_df=anchor_df,
        group_anchor_df=group_anchor_df,
        noncenter_facet=noncenter_facet,
        dummy_facets=dummy_facets,
    )
    config["theta_spec"] = constraint_specs["theta_spec"]
    config["facet_specs"] = constraint_specs["facet_specs"]
    config["noncenter_facet"] = noncenter_facet
    config["dummy_facets"] = dummy_facets
    config["anchor_summary"] = constraint_specs["anchor_summary"]

    sizes = build_param_sizes(config)

    step_init = np.linspace(-1, 1, max(n_cat - 1, 0)) if n_cat > 1 else np.array([], dtype=float)
    facet_starts = np.concatenate([np.zeros(sizes[f]) for f in config["facet_names"]]) if config["facet_names"] else np.array([], dtype=float)
    start = np.concatenate([
        np.zeros(sizes["theta"]),
        facet_starts,
        step_init if model == "RSM" else np.tile(step_init, len(facet_levels[step_facet])),
    ])

    if method == "JMLE":
        opt = minimize(
            mfrm_loglik_jmle,
            start,
            args=(idx, config, sizes),
            method="BFGS",
            options={"maxiter": maxit, "gtol": reltol},
        )
    else:
        quad = gauss_hermite_normal(quad_points)
        opt = minimize(
            mfrm_loglik_mml,
            start,
            args=(idx, config, sizes, quad),
            method="BFGS",
            options={"maxiter": maxit, "gtol": reltol},
        )

    params = expand_params(opt.x, sizes, config)

    if method == "MML":
        quad = gauss_hermite_normal(quad_points)
        person_tbl = compute_person_eap(idx, config, params, quad)
        person_tbl.insert(0, "Person", prep["levels"]["Person"])
    else:
        person_tbl = pd.DataFrame({
            "Person": prep["levels"]["Person"],
            "Estimate": params["theta"],
        })

    facet_tbls = []
    for facet in config["facet_names"]:
        facet_tbls.append(pd.DataFrame({
            "Facet": facet,
            "Level": prep["levels"][facet],
            "Estimate": params["facets"][facet],
        }))
    facet_tbl = pd.concat(facet_tbls, ignore_index=True) if facet_tbls else pd.DataFrame(columns=["Facet", "Level", "Estimate"])

    if model == "RSM":
        step_tbl = pd.DataFrame({
            "Step": [f"Step_{i+1}" for i in range(n_cat - 1)],
            "Estimate": params["steps"],
        })
    else:
        rows = []
        for sf in prep["levels"][step_facet]:
            for i in range(n_cat - 1):
                rows.append((sf, f"Step_{i+1}", params["steps_mat"][prep["levels"][step_facet].index(sf), i]))
        step_tbl = pd.DataFrame(rows, columns=["StepFacet", "Step", "Estimate"])

    k_params = sum(sizes.values())
    loglik = -opt.fun
    if "Weight" in prep["data"].columns:
        n_obs = float(prep["data"]["Weight"].sum())
    else:
        n_obs = len(prep["data"])
    aic = 2 * k_params - 2 * loglik
    bic = np.log(n_obs) * k_params - 2 * loglik if n_obs > 0 else np.nan

    summary_tbl = pd.DataFrame({
        "Model": [model],
        "Method": [method],
        "N": [n_obs],
        "Persons": [n_person],
        "Facets": [len(config["facet_names"])],
        "Categories": [n_cat],
        "LogLik": [loglik],
        "AIC": [aic],
        "BIC": [bic],
        "Converged": [bool(opt.success)],
        "Iterations": [opt.nfev],
    })

    return {
        "summary": summary_tbl,
        "facets": {"person": person_tbl, "others": facet_tbl},
        "steps": step_tbl,
        "config": config,
        "prep": prep,
        "opt": opt,
        "params": params,
    }


def compute_obs_table(res):
    prep = res["prep"]
    config = res["config"]
    idx = build_indices(prep, step_facet=config["step_facet"])
    params = res["params"]
    theta_hat = params["theta"] if config["method"] == "JMLE" else res["facets"]["person"]["Estimate"].to_numpy()

    person_measure = theta_hat
    person_measure_by_row = person_measure[idx["person"]] if person_measure.size else np.zeros_like(idx["score_k"], dtype=float)

    eta = compute_eta(idx, params, config, theta_override=theta_hat if theta_hat.size else None)
    if config["model"] == "RSM":
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
        probs = category_prob_rsm(eta, step_cum)
    else:
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])
        probs = category_prob_pcm(eta, step_cum_mat, idx["step_idx"])

    k_vals = np.arange(probs.shape[1])
    expected_k = probs.dot(k_vals)
    var_k = probs.dot(k_vals ** 2) - expected_k ** 2
    var_k = np.where(var_k <= 1e-10, np.nan, var_k)
    resid_k = idx["score_k"] - expected_k
    std_sq = resid_k ** 2 / var_k

    df = prep["data"].copy()
    df["PersonMeasure"] = person_measure_by_row
    df["Observed"] = prep["rating_min"] + idx["score_k"]
    df["Expected"] = prep["rating_min"] + expected_k
    df["Var"] = var_k
    df["Residual"] = df["Observed"] - df["Expected"]
    df["StdResidual"] = df["Residual"] / np.sqrt(df["Var"])
    df["StdSq"] = std_sq
    return df


def compute_prob_matrix(res):
    prep = res["prep"]
    config = res["config"]
    idx = build_indices(prep, step_facet=config["step_facet"])
    params = res["params"]
    theta_hat = params["theta"] if config["method"] == "JMLE" else res["facets"]["person"]["Estimate"].to_numpy()
    eta = compute_eta(idx, params, config, theta_override=theta_hat if theta_hat.size else None)
    if config["model"] == "RSM":
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
        return category_prob_rsm(eta, step_cum)
    step_cum_mat = np.vstack([
        np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
    ])
    return category_prob_pcm(eta, step_cum_mat, idx["step_idx"])


def compute_scorefile(res):
    obs = compute_obs_table(res)
    if obs.empty:
        return pd.DataFrame()
    probs = compute_prob_matrix(res)
    if probs is None or probs.shape[0] != len(obs):
        return obs
    cat_vals = np.arange(res["prep"]["rating_min"], res["prep"]["rating_max"] + 1)
    prob_df = pd.DataFrame(probs, columns=[f"P_{v}" for v in cat_vals])
    max_idx = np.argmax(probs, axis=1)
    max_prob = probs[np.arange(len(probs)), max_idx]
    most_likely = cat_vals[max_idx]
    base_cols = ["Person"] + res["config"]["facet_names"]
    if "Weight" in obs.columns:
        base_cols.append("Weight")
    base_cols += ["Score", "Observed", "Expected", "Residual", "StdResidual", "Var", "PersonMeasure"]
    out = pd.concat([
        obs[base_cols],
        prob_df,
    ], axis=1)
    out["MostLikely"] = most_likely
    out["MaxProb"] = max_prob
    return out


def compute_residual_file(res):
    obs = compute_obs_table(res)
    if obs.empty:
        return pd.DataFrame()
    cols = ["Person", *res["config"]["facet_names"]]
    if "Weight" in obs.columns:
        cols.append("Weight")
    cols += ["Score", "Observed", "Expected", "Residual", "StdResidual", "Var", "StdSq", "PersonMeasure"]
    return obs[cols]


def get_weights(df):
    if "Weight" in df.columns:
        w = df["Weight"].to_numpy(dtype=float)
        w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
        return w
    return np.ones(len(df), dtype=float)


def calc_overall_fit(obs_df, whexact=False):
    w = get_weights(obs_df)
    infit = np.nansum(obs_df["StdSq"] * obs_df["Var"] * w) / np.nansum(obs_df["Var"] * w)
    outfit = np.nansum(obs_df["StdSq"] * w) / np.nansum(w)
    df_infit = np.nansum(obs_df["Var"] * w)
    df_outfit = np.nansum(w)
    return pd.DataFrame({
        "Infit": [infit],
        "Outfit": [outfit],
        "InfitZSTD": [zstd_from_mnsq(infit, df_infit, whexact=whexact)],
        "OutfitZSTD": [zstd_from_mnsq(outfit, df_outfit, whexact=whexact)],
        "DF_Infit": [df_infit],
        "DF_Outfit": [df_outfit],
    })


def calc_facet_fit(obs_df, facet_cols, whexact=False):
    rows = []
    for facet in facet_cols:
        grp = obs_df.groupby(facet, observed=False)
        for level, df in grp:
            w = get_weights(df)
            infit = np.nansum(df["StdSq"] * df["Var"] * w) / np.nansum(df["Var"] * w)
            outfit = np.nansum(df["StdSq"] * w) / np.nansum(w)
            df_infit = np.nansum(df["Var"] * w)
            df_outfit = np.nansum(w)
            rows.append({
                "Facet": facet,
                "Level": level,
                "N": np.nansum(w),
                "Infit": infit,
                "Outfit": outfit,
                "InfitZSTD": zstd_from_mnsq(infit, df_infit, whexact=whexact),
                "OutfitZSTD": zstd_from_mnsq(outfit, df_outfit, whexact=whexact),
                "DF_Infit": df_infit,
                "DF_Outfit": df_outfit,
            })
    return pd.DataFrame(rows)


def calc_facet_se(obs_df, facet_cols):
    rows = []
    for facet in facet_cols:
        grp = obs_df.groupby(facet, observed=False)
        for level, df in grp:
            w = get_weights(df)
            info = np.nansum(df["Var"] * w)
            se = 1 / np.sqrt(info) if info > 0 else np.nan
            rows.append({
                "Facet": facet,
                "Level": level,
                "N": np.nansum(w),
                "SE": se,
            })
    return pd.DataFrame(rows)


def calc_bias_facet(obs_df, facet_cols):
    rows = []
    for facet in facet_cols:
        grp = obs_df.groupby(facet, observed=False)
        for level, df in grp:
            w = get_weights(df)
            n = np.nansum(w)
            observed_avg = weighted_mean(df["Observed"].to_numpy(), w)
            expected_avg = weighted_mean(df["Expected"].to_numpy(), w)
            mean_resid = weighted_mean(df["Residual"].to_numpy(), w)
            mean_std_resid = weighted_mean(df["StdResidual"].to_numpy(), w)
            mean_abs_std = weighted_mean(np.abs(df["StdResidual"]).to_numpy(), w)
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
                "Facet": facet,
                "Level": level,
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
    return pd.DataFrame(rows)


def calc_bias_interactions(obs_df, facet_cols, pairs=None, top_n=20):
    if len(facet_cols) < 2:
        return pd.DataFrame()
    if pairs is None:
        combos = []
        for i in range(len(facet_cols)):
            for j in range(i + 1, len(facet_cols)):
                combos.append((facet_cols[i], facet_cols[j]))
    elif len(pairs) == 0:
        return pd.DataFrame()
    else:
        combos = pairs

    rows = []
    for pair1, pair2 in combos:
        grp = obs_df.groupby([pair1, pair2], observed=False)
        for (lvl1, lvl2), df in grp:
            w = get_weights(df)
            n = np.nansum(w)
            observed_avg = weighted_mean(df["Observed"].to_numpy(), w)
            expected_avg = weighted_mean(df["Expected"].to_numpy(), w)
            mean_resid = weighted_mean(df["Residual"].to_numpy(), w)
            mean_std_resid = weighted_mean(df["StdResidual"].to_numpy(), w)
            mean_abs_std = weighted_mean(np.abs(df["StdResidual"]).to_numpy(), w)
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
    out = out.sort_values("AbsStd", ascending=False).drop(columns=["AbsStd"]).head(top_n)
    return out


def weighted_mean(x, w):
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not np.any(ok):
        return np.nan
    w_sum = np.sum(w[ok])
    if w_sum <= 0:
        return np.nan
    return float(np.sum(x[ok] * w[ok]) / w_sum)


def build_histogram(values, bins=30):
    vals = np.array(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return pd.DataFrame(), {}
    hist, edges = np.histogram(vals, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    stats = {
        "mean": float(np.mean(vals)),
        "sd": float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
        "n": int(vals.size),
    }
    return pd.DataFrame({"Value": centers, "Count": hist}), stats


def format_fixed_width_table(
    df,
    columns,
    formats=None,
    right_align=None,
    max_col_width=16,
    min_col_width=6,
):
    if df is None or df.empty:
        return "No data"
    formats = formats or {}
    if right_align is None:
        right_align = {
            col for col in columns
            if col in df.columns and pd.api.types.is_numeric_dtype(df[col])
        }

    def fmt_val(col, val):
        if pd.isna(val):
            return ""
        fmt = formats.get(col)
        if callable(fmt):
            return fmt(val)
        if isinstance(fmt, str):
            try:
                return fmt.format(val)
            except Exception:
                return str(val)
        return str(val)

    str_cols = {}
    widths = {}
    for col in columns:
        if col not in df.columns:
            str_cols[col] = [""] * len(df)
            widths[col] = max(len(col), min_col_width)
            continue
        vals = [fmt_val(col, v) for v in df[col].tolist()]
        str_cols[col] = vals
        max_len = max([len(col)] + [len(v) for v in vals])
        widths[col] = max(min_col_width, min(max_len, max_col_width))

    def pad(col, text):
        text = text[: widths[col]]
        if col in right_align:
            return text.rjust(widths[col])
        return text.ljust(widths[col])

    header = " ".join([pad(col, col) for col in columns])
    lines = [header]
    for i in range(len(df)):
        row = " ".join([pad(col, str_cols[col][i]) for col in columns])
        lines.append(row)
    return "\n".join(lines)


def build_bias_fixed_text(
    table_df,
    summary_df,
    chi_df,
    facet_a,
    facet_b,
    columns,
    formats,
):
    if table_df is None or table_df.empty:
        return "No bias data"
    combined = table_df.copy()
    if summary_df is not None and not summary_df.empty:
        summary_dicts = [
            {k.replace("_", " "): v for k, v in row._asdict().items()}
            for row in summary_df.itertuples(index=False)
        ]
        sum_rows = []
        for row in summary_dicts:
            row_dict = {col: np.nan for col in columns}
            row_dict["Sq"] = row.get("Statistic", "")
            col_map = {
                "Observd Score": "Observd Score",
                "Expctd Score": "Expctd Score",
                "Observd Count": "Observd Count",
                "Obs-Exp Average": "Obs-Exp Average",
                "Bias Size": "Bias Size",
                "Model S.E.": "S.E.",
            }
            for out_col, src_col in col_map.items():
                if src_col in row:
                    row_dict[out_col] = row[src_col]
            sum_rows.append(row_dict)
        blank_row = {col: np.nan for col in columns}
        combined = pd.concat([combined, pd.DataFrame([blank_row] + sum_rows)], ignore_index=True)

    fixed_table = format_fixed_width_table(combined, columns, formats=formats, max_col_width=18)
    lines = [f"Bias/Interaction: {facet_a} x {facet_b}", "", fixed_table]
    if chi_df is not None and not chi_df.empty:
        chi = chi_df.iloc[0]
        chi_line = (
            f"Fixed (all = 0) chi-squared: {chi['FixedChiSq']:.2f}  "
            f"d.f.: {int(chi['FixedDF']) if pd.notna(chi['FixedDF']) else ''}  "
            f"significance (probability): {chi['FixedProb']:.4f}"
            if pd.notna(chi["FixedChiSq"])
            else "Fixed (all = 0) chi-squared: N/A"
        )
        lines.extend(["", chi_line])
    return "\n".join(lines)


def build_pairwise_fixed_text(pair_df, target_facet, context_facet, columns, formats):
    if pair_df is None or pair_df.empty:
        return "No pairwise data"
    fixed_table = format_fixed_width_table(pair_df, columns, formats=formats, max_col_width=18)
    lines = [
        f"Bias/Interaction Pairwise Report: Target={target_facet}  Context={context_facet}",
        "",
        fixed_table,
    ]
    return "\n".join(lines)


def get_extreme_levels(obs_df, facet_names, rating_min, rating_max):
    extreme_levels = {}
    for facet in facet_names:
        if facet not in obs_df.columns:
            extreme_levels[facet] = set()
            continue
        stat = (
            obs_df.groupby(facet)["Observed"]
            .agg(MinScore="min", MaxScore="max")
            .reset_index()
        )
        extreme = stat[
            ((stat["MinScore"] == rating_min) & (stat["MaxScore"] == rating_min))
            | ((stat["MinScore"] == rating_max) & (stat["MaxScore"] == rating_max))
        ]
        extreme_levels[facet] = set(extreme[facet].astype(str).tolist())
    return extreme_levels


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
    if res is None or diagnostics is None:
        return {}
    obs_df = diagnostics.get("obs")
    if obs_df is None or obs_df.empty:
        return {}
    if facet_a == facet_b:
        return {}

    facet_names = ["Person"] + res["config"]["facet_names"]
    if facet_a not in facet_names or facet_b not in facet_names:
        return {}

    prep = res["prep"]
    config = res["config"]
    params = res["params"]
    idx = build_indices(prep, step_facet=config["step_facet"])
    theta_hat = params["theta"] if config["method"] == "JMLE" else res["facets"]["person"]["Estimate"].to_numpy()
    eta_base = compute_eta(idx, params, config, theta_override=theta_hat if theta_hat.size else None)
    score_k = idx["score_k"]
    weight = idx.get("weight")
    step_idx = idx.get("step_idx")

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
    group_indices = obs_df.groupby([facet_a, facet_b]).indices
    groups = []
    for (lvl_a, lvl_b), idx_rows in group_indices.items():
        lvl_a_str = str(lvl_a)
        lvl_b_str = str(lvl_b)
        if omit_extreme:
            if lvl_a_str in extreme_levels.get(facet_a, set()) or lvl_b_str in extreme_levels.get(facet_b, set()):
                continue
        idx_rows = np.asarray(idx_rows, dtype=int)
        if idx_rows.size == 0:
            continue
        groups.append({
            "key": (lvl_a_str, lvl_b_str),
            "idx": idx_rows,
        })

    def estimate_bias_for_group(idx_rows):
        eta_sub = eta_base[idx_rows]
        score_k_sub = score_k[idx_rows]
        weight_sub = weight[idx_rows] if weight is not None else None
        step_idx_sub = step_idx[idx_rows] if step_idx is not None else None

        if config["model"] == "RSM":
            def nll(b):
                return -loglik_rsm(eta_sub + b, score_k_sub, step_cum, weight=weight_sub)
        else:
            def nll(b):
                return -loglik_pcm(eta_sub + b, score_k_sub, step_cum_mat, step_idx_sub, weight=weight_sub)

        try:
            opt = minimize_scalar(nll, bounds=(-max_abs, max_abs), method="bounded")
            if opt.success:
                return float(opt.x)
        except Exception:
            return np.nan
        return np.nan

    def iteration_metrics(bias_map):
        max_resid = 0.0
        max_resid_pct = np.nan
        max_categories = np.nan
        for g in groups:
            idx_rows = g["idx"]
            bias = bias_map.get(g["key"], 0.0)
            eta_sub = eta_base[idx_rows] + (bias if np.isfinite(bias) else 0.0)
            score_k_sub = score_k[idx_rows]
            step_idx_sub = step_idx[idx_rows] if step_idx is not None else None
            if config["model"] == "RSM":
                probs = category_prob_rsm(eta_sub, step_cum)
            else:
                probs = category_prob_pcm(eta_sub, step_cum_mat, step_idx_sub)
            k_vals = np.arange(probs.shape[1])
            expected_k = probs.dot(k_vals)
            expected_score = prep["rating_min"] + expected_k
            obs_score = obs_df.iloc[idx_rows]["Observed"].to_numpy(dtype=float)
            if "Weight" in obs_df.columns:
                w = obs_df.iloc[idx_rows]["Weight"].to_numpy(dtype=float)
                obs_score = obs_score * w
                expected_score = expected_score * w
            resid = obs_score - expected_score
            if resid.size:
                obs_sum = float(np.nansum(obs_score))
                exp_sum = float(np.nansum(expected_score))
                resid_sum = obs_sum - exp_sum
                if abs(resid_sum) >= abs(max_resid):
                    max_resid = resid_sum
                    max_resid_pct = (resid_sum / exp_sum * 100) if exp_sum != 0 else np.nan
                    max_categories = np.nan

        return {
            "max_resid": max_resid,
            "max_resid_pct": max_resid_pct,
            "max_resid_categories": max_categories,
        }

    bias_map = {g["key"]: 0.0 for g in groups}
    iter_rows = []
    for it in range(1, max_iter + 1):
        max_change_abs = 0.0
        max_change_signed = 0.0
        changes = []
        for g in groups:
            key = g["key"]
            bias_hat = estimate_bias_for_group(g["idx"])
            prev = bias_map.get(key, 0.0)
            if np.isfinite(bias_hat) and np.isfinite(prev):
                delta = bias_hat - prev
                if abs(delta) >= max_change_abs:
                    max_change_abs = abs(delta)
                    max_change_signed = delta
                changes.append(abs(delta))
            else:
                changes.append(np.nan)
            bias_map[key] = bias_hat
        resid_info = iteration_metrics(bias_map)
        iter_rows.append({
            "Iteration": it,
            "MaxScoreResidual": resid_info["max_resid"],
            "MaxScoreResidualPct": resid_info["max_resid_pct"],
            "MaxScoreResidualCategories": resid_info["max_resid_categories"],
            "MaxLogitChange": max_change_signed,
            "BiasCells": int(np.sum(np.array(changes, dtype=float) > tol)),
        })
        if max_change_abs < tol:
            break

    rows = []
    seq = 1
    for g in groups:
        lvl_a_str, lvl_b_str = g["key"]
        idx_rows = g["idx"]
        bias_hat = bias_map.get(g["key"], np.nan)
        eta_sub = eta_base[idx_rows]
        score_k_sub = score_k[idx_rows]
        weight_sub = weight[idx_rows] if weight is not None else None
        step_idx_sub = step_idx[idx_rows] if step_idx is not None else None

        if config["model"] == "RSM":
            probs = category_prob_rsm(eta_sub + (bias_hat if np.isfinite(bias_hat) else 0.0), step_cum)
        else:
            probs = category_prob_pcm(eta_sub + (bias_hat if np.isfinite(bias_hat) else 0.0), step_cum_mat, step_idx_sub)

        k_vals = np.arange(probs.shape[1])
        expected_k = probs.dot(k_vals)
        var_k = probs.dot(k_vals ** 2) - expected_k ** 2
        var_k = np.where(var_k <= 1e-10, np.nan, var_k)
        resid_k = score_k_sub - expected_k
        std_sq = resid_k ** 2 / var_k

        w = weight_sub if weight_sub is not None else np.ones(len(idx_rows))
        info = np.nansum(var_k * w)
        se = 1 / np.sqrt(info) if np.isfinite(info) and info > 0 else np.nan
        infit = np.nansum(std_sq * var_k * w) / np.nansum(var_k * w) if np.nansum(var_k * w) > 0 else np.nan
        outfit = np.nansum(std_sq * w) / np.nansum(w) if np.nansum(w) > 0 else np.nan

        obs_slice = obs_df.iloc[idx_rows]
        w_obs = obs_slice["Weight"].to_numpy(dtype=float) if "Weight" in obs_slice.columns else np.ones(len(obs_slice))
        obs_score = float(np.nansum(obs_slice["Observed"].to_numpy(dtype=float) * w_obs))
        exp_score = float(np.nansum(obs_slice["Expected"].to_numpy(dtype=float) * w_obs))
        obs_count = float(np.nansum(w_obs))
        obs_exp_avg = (obs_score - exp_score) / obs_count if obs_count > 0 else np.nan

        n_obs = int(len(obs_slice))
        df_t = max(n_obs - 1, 0)
        t_val = bias_hat / se if np.isfinite(bias_hat) and np.isfinite(se) and se > 0 else np.nan
        p_val = 2 * t_dist.cdf(-abs(t_val), df=df_t) if np.isfinite(t_val) and df_t > 0 else np.nan

        rows.append({
            "Sq": seq,
            "Observd Score": obs_score,
            "Expctd Score": exp_score,
            "Observd Count": obs_count,
            "Obs-Exp Average": obs_exp_avg,
            "Bias Size": bias_hat,
            "S.E.": se,
            "t": t_val,
            "d.f.": df_t,
            "Prob.": p_val,
            "Infit": infit,
            "Outfit": outfit,
            "ObsN": n_obs,
            "FacetA": facet_a,
            "FacetA_Level": lvl_a_str,
            "FacetA_Index": level_map.get(facet_a, []).index(lvl_a_str) + 1 if lvl_a_str in level_map.get(facet_a, []) else np.nan,
            "FacetA_Measure": meas_map.get((facet_a, lvl_a_str), np.nan),
            "FacetA_SE": se_map.get((facet_a, lvl_a_str), np.nan),
            "FacetB": facet_b,
            "FacetB_Level": lvl_b_str,
            "FacetB_Index": level_map.get(facet_b, []).index(lvl_b_str) + 1 if lvl_b_str in level_map.get(facet_b, []) else np.nan,
            "FacetB_Measure": meas_map.get((facet_b, lvl_b_str), np.nan),
            "FacetB_SE": se_map.get((facet_b, lvl_b_str), np.nan),
        })
        seq += 1

    bias_tbl = pd.DataFrame(rows)
    if bias_tbl.empty:
        return {}

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
    fixed_df = max(len(bias_tbl) - 1, 0)
    fixed_prob = 1 - chi2.cdf(fixed_chi, df=fixed_df) if np.isfinite(fixed_chi) and fixed_df > 0 else np.nan
    chi_tbl = pd.DataFrame({
        "FixedChiSq": [fixed_chi],
        "FixedDF": [fixed_df],
        "FixedProb": [fixed_prob],
    })

    return {
        "facet_a": facet_a,
        "facet_b": facet_b,
        "table": bias_tbl,
        "summary": summary_tbl,
        "chi_sq": chi_tbl,
        "iteration": pd.DataFrame(iter_rows),
    }


def calc_bias_pairwise(bias_tbl, target_facet, context_facet):
    if bias_tbl is None or bias_tbl.empty:
        return pd.DataFrame()

    use_a = bias_tbl["FacetA"].iloc[0] == target_facet
    use_b = bias_tbl["FacetB"].iloc[0] == target_facet
    if not (use_a or use_b):
        return pd.DataFrame()

    if use_a:
        target_prefix = "FacetA"
        context_prefix = "FacetB"
    else:
        target_prefix = "FacetB"
        context_prefix = "FacetA"

    sub = bias_tbl.copy()
    sub = sub[sub[f"{context_prefix}"].isin([context_facet])]
    if sub.empty:
        sub = bias_tbl.copy()

    rows = []
    for tgt_level, group_df in sub.groupby(f"{target_prefix}_Level"):
        contexts = group_df[f"{context_prefix}_Level"].unique().tolist()
        if len(contexts) < 2:
            continue
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

            local1 = tgt_measure + bias1 if np.isfinite(tgt_measure) and np.isfinite(bias1) else np.nan
            local2 = tgt_measure + bias2 if np.isfinite(tgt_measure) and np.isfinite(bias2) else np.nan
            se1 = np.sqrt(tgt_se ** 2 + bias_se1 ** 2) if np.isfinite(tgt_se) and np.isfinite(bias_se1) else np.nan
            se2 = np.sqrt(tgt_se ** 2 + bias_se2 ** 2) if np.isfinite(tgt_se) and np.isfinite(bias_se2) else np.nan

            contrast = local1 - local2 if np.isfinite(local1) and np.isfinite(local2) else np.nan
            se_contrast = np.sqrt(se1 ** 2 + se2 ** 2) if np.isfinite(se1) and np.isfinite(se2) else np.nan
            t_val = contrast / se_contrast if np.isfinite(contrast) and np.isfinite(se_contrast) and se_contrast > 0 else np.nan

            n1 = int(r1["ObsN"]) if np.isfinite(r1["ObsN"]) else 0
            n2 = int(r2["ObsN"]) if np.isfinite(r2["ObsN"]) else 0
            df_num = (se1 ** 2 + se2 ** 2) ** 2
            df_den = 0.0
            if n1 > 1 and np.isfinite(se1):
                df_den += (se1 ** 4) / (n1 - 1)
            if n2 > 1 and np.isfinite(se2):
                df_den += (se2 ** 4) / (n2 - 1)
            df_t_val = df_num / df_den if df_den > 0 else np.nan
            p_val = 2 * t_dist.cdf(-abs(t_val), df=df_t_val) if np.isfinite(t_val) and np.isfinite(df_t_val) and df_t_val > 0 else np.nan

            rows.append({
                "Target": tgt_level,
                "Target N": tgt_index,
                "Target Measure": tgt_measure,
                "Target S.E.": tgt_se,
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
                "t": t_val,
                "d.f.": df_t_val,
                "Prob.": p_val,
            })

    return pd.DataFrame(rows)


def safe_cor(x, y, w=None):
    ok = np.isfinite(x) & np.isfinite(y)
    if w is None:
        if not np.any(ok):
            return np.nan
        x = x[ok]
        y = y[ok]
        if len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
            return np.nan
        return float(np.corrcoef(x, y)[0, 1])

    ok = ok & np.isfinite(w) & (w > 0)
    if not np.any(ok):
        return np.nan
    x = x[ok]
    y = y[ok]
    w = w[ok]
    w_sum = np.sum(w)
    if w_sum <= 0:
        return np.nan
    mx = np.sum(w * x) / w_sum
    my = np.sum(w * y) / w_sum
    vx = np.sum(w * (x - mx) ** 2) / w_sum
    vy = np.sum(w * (y - my) ** 2) / w_sum
    if vx <= 0 or vy <= 0:
        return np.nan
    cov = np.sum(w * (x - mx) * (y - my)) / w_sum
    return float(cov / np.sqrt(vx * vy))


def weighted_mean_safe(x, w):
    ok = np.isfinite(x) & np.isfinite(w)
    if not np.any(ok):
        return np.nan
    return float(np.sum(x[ok] * w[ok]) / np.sum(w[ok]))


def calc_interrater_agreement(obs_df, facet_cols, rater_facet, res=None):
    if obs_df is None or obs_df.empty:
        return {"summary": pd.DataFrame(), "pairs": pd.DataFrame()}
    if rater_facet is None or rater_facet not in facet_cols:
        return {"summary": pd.DataFrame(), "pairs": pd.DataFrame()}
    context_cols = [c for c in facet_cols if c != rater_facet]
    if len(context_cols) == 0:
        return {"summary": pd.DataFrame(), "pairs": pd.DataFrame()}

    df = obs_df.copy()
    df["_context"] = df[context_cols].astype(str).agg("|".join, axis=1)
    base_cols = ["_context", rater_facet, "Observed"]
    if "Weight" in df.columns:
        base_cols.append("Weight")
    df = df[base_cols].copy()
    if "Weight" in df.columns:
        df_obs = (
            df.groupby(["_context", rater_facet])
            .apply(lambda g: weighted_mean(g["Observed"].to_numpy(), g["Weight"].to_numpy()))
            .reset_index(name="Score")
        )
    else:
        df_obs = df.groupby(["_context", rater_facet], as_index=False).agg(Score=("Observed", "mean"))
    if df_obs.empty:
        return {"summary": pd.DataFrame(), "pairs": pd.DataFrame()}

    prob_map = {}
    if res is not None:
        probs = compute_prob_matrix(res)
        if probs is not None and probs.shape[0] == len(obs_df):
            prob_cols = [f"_p{i}" for i in range(probs.shape[1])]
            df_probs = df[["_context", rater_facet]].copy()
            for i, col in enumerate(prob_cols):
                df_probs[col] = probs[:, i]
            if "Weight" in df.columns:
                df_probs["Weight"] = df["Weight"].to_numpy()
            for (ctx, r), g in df_probs.groupby(["_context", rater_facet]):
                if "Weight" in g.columns:
                    w = g["Weight"].to_numpy(dtype=float)
                    ok = np.isfinite(w) & (w > 0)
                    if not np.any(ok):
                        avg = np.full(len(prob_cols), np.nan)
                    else:
                        avg = np.average(g.loc[ok, prob_cols].to_numpy(), axis=0, weights=w[ok])
                else:
                    avg = np.nanmean(g[prob_cols].to_numpy(), axis=0)
                prob_map[(ctx, r)] = avg

    wide = df_obs.pivot(index="_context", columns=rater_facet, values="Score")
    if wide is None or wide.empty:
        return {"summary": pd.DataFrame(), "pairs": pd.DataFrame()}

    rater_cols = [c for c in wide.columns if c != "_context"]
    if len(rater_cols) < 2:
        return {"summary": pd.DataFrame(), "pairs": pd.DataFrame()}

    pairs = []
    for i in range(len(rater_cols)):
        for j in range(i + 1, len(rater_cols)):
            pairs.append((rater_cols[i], rater_cols[j]))

    pair_rows = []
    total_pairs = 0
    total_exact = 0
    total_expected = 0.0
    expected_available = False
    for r1, r2 in pairs:
        sub = wide[[r1, r2]].dropna()
        n_ok = int(len(sub))
        if n_ok == 0:
            pair_rows.append({
                "Rater1": r1,
                "Rater2": r2,
                "N": 0,
                "Exact": np.nan,
                "ExpectedExact": np.nan,
                "Adjacent": np.nan,
                "MeanDiff": np.nan,
                "MAD": np.nan,
                "Corr": np.nan,
            })
            continue

        v1 = sub[r1].to_numpy()
        v2 = sub[r2].to_numpy()
        diff = v1 - v2
        exact_count = int(np.sum(np.isclose(diff, 0)))

        exp_vals = []
        if prob_map:
            for ctx in sub.index:
                p1 = prob_map.get((ctx, r1))
                p2 = prob_map.get((ctx, r2))
                if p1 is None or p2 is None:
                    continue
                if np.any(np.isnan(p1)) or np.any(np.isnan(p2)):
                    continue
                exp_vals.append(float(np.sum(p1 * p2)))
        exp_mean = float(np.mean(exp_vals)) if exp_vals else np.nan
        exp_sum = float(np.sum(exp_vals)) if exp_vals else 0.0
        if exp_vals:
            expected_available = True

        pair_rows.append({
            "Rater1": r1,
            "Rater2": r2,
            "N": n_ok,
            "Exact": float(exact_count / n_ok),
            "ExpectedExact": exp_mean,
            "Adjacent": float(np.mean(np.abs(diff) <= 1)),
            "MeanDiff": float(np.mean(diff)),
            "MAD": float(np.mean(np.abs(diff))),
            "Corr": safe_cor(v1, v2),
        })

        total_pairs += n_ok
        total_exact += exact_count
        total_expected += exp_sum

    pair_tbl = pd.DataFrame(pair_rows)
    contexts_with_pairs = int(np.sum(np.sum(np.isfinite(wide.to_numpy()), axis=1) >= 2))

    if not expected_available:
        total_expected = np.nan
    summary_tbl = pd.DataFrame({
        "RaterFacet": [rater_facet],
        "Raters": [len(rater_cols)],
        "Pairs": [len(pair_tbl)],
        "Contexts": [contexts_with_pairs],
        "TotalPairs": [total_pairs],
        "ExactAgreements": [total_exact],
        "ExpectedAgreements": [total_expected if total_pairs > 0 else np.nan],
        "ExactAgreement": [total_exact / total_pairs if total_pairs > 0 else np.nan],
        "ExpectedExactAgreement": [total_expected / total_pairs if total_pairs > 0 else np.nan],
        "AdjacentAgreement": [weighted_mean_safe(pair_tbl["Adjacent"].to_numpy(), pair_tbl["N"].to_numpy()) if not pair_tbl.empty else np.nan],
        "MeanAbsDiff": [weighted_mean_safe(pair_tbl["MAD"].to_numpy(), pair_tbl["N"].to_numpy()) if not pair_tbl.empty else np.nan],
        "MeanCorr": [weighted_mean_safe(pair_tbl["Corr"].to_numpy(), pair_tbl["N"].to_numpy()) if not pair_tbl.empty else np.nan],
    })

    return {"summary": summary_tbl, "pairs": pair_tbl}


def calc_ptmea(obs_df, facet_cols):
    facet_cols = [f for f in facet_cols if f != "Person"]
    rows = []
    for facet in facet_cols:
        grp = obs_df.groupby(facet, observed=False)
        for level, df in grp:
            w = get_weights(df)
            rows.append({
                "Facet": facet,
                "Level": level,
                "PTMEA": safe_cor(df["Observed"].to_numpy(), df["PersonMeasure"].to_numpy(), w=w),
                "N": np.nansum(w),
            })
    return pd.DataFrame(rows)


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


def expected_score_from_eta(eta, step_cum, rating_min):
    if not np.isfinite(eta) or step_cum is None or len(step_cum) == 0:
        return np.nan
    probs = category_prob_rsm(np.array([eta], dtype=float), np.array(step_cum, dtype=float))
    k_vals = np.arange(probs.shape[1])
    return float(rating_min + np.sum(probs[0] * k_vals))


def estimate_eta_from_target(target, step_cum, rating_min, rating_max):
    if not np.isfinite(target) or step_cum is None or len(step_cum) == 0:
        return np.nan
    if target <= rating_min:
        return -np.inf
    if target >= rating_max:
        return np.inf

    def f(eta):
        return expected_score_from_eta(eta, step_cum, rating_min) - target

    low, high = -10.0, 10.0
    f_low = f(low)
    f_high = f(high)
    if not np.isfinite(f_low) or not np.isfinite(f_high) or f_low * f_high > 0:
        low, high = -20.0, 20.0
        f_low = f(low)
        f_high = f(high)
        if not np.isfinite(f_low) or not np.isfinite(f_high) or f_low * f_high > 0:
            return np.nan

    try:
        root = root_scalar(f, bracket=(low, high), method="brentq")
        return float(root.root) if root.converged else np.nan
    except Exception:
        return np.nan


def facet_anchor_status(facet, levels, config, extreme_levels=None):
    spec = config["theta_spec"] if facet == "Person" else config["facet_specs"].get(facet)
    if spec is None:
        return [""] * len(levels)
    anchors = spec.get("anchors")
    groups = spec.get("groups")
    status = []
    for lvl in levels:
        idx = spec["levels"].index(str(lvl)) if str(lvl) in spec["levels"] else None
        if idx is not None and anchors is not None and np.isfinite(anchors[idx]):
            status.append("A")
        elif idx is not None and groups is not None:
            grp_val = groups[idx]
            if grp_val not in (None, "", np.nan):
                if extreme_levels is not None and str(lvl) in extreme_levels:
                    status.append("X")
                else:
                    status.append("G")
            else:
                status.append("")
        else:
            status.append("")
    return status


def calc_facets_report_tbls(
    res,
    diagnostics,
    totalscore=True,
    umean=0,
    uscale=1,
    udecimals=2,
    omit_unobserved=False,
    xtreme=0,
):
    if res is None or diagnostics is None:
        return {}
    obs_df = diagnostics.get("obs")
    measures = diagnostics.get("measures")
    if obs_df is None or obs_df.empty or measures is None or measures.empty:
        return {}

    prep = res["prep"]
    config = res["config"]
    rating_min = prep["rating_min"]
    rating_max = prep["rating_max"]
    params = res["params"]

    theta_hat = params["theta"] if config["method"] == "JMLE" else res["facets"]["person"]["Estimate"].to_numpy()
    theta_mean = float(np.mean(theta_hat)) if len(theta_hat) > 0 else 0.0
    facet_means = {f: float(np.nanmean(params["facets"][f])) for f in config["facet_names"]}
    facet_signs = config.get("facet_signs", {f: -1 for f in config["facet_names"]})

    if config["model"] == "RSM":
        step_cum_common = np.concatenate([[0.0], np.cumsum(params["steps"])])
        step_cum_mean = step_cum_common
    else:
        step_mat = params["steps_mat"]
        if step_mat is None or len(step_mat) == 0:
            step_cum_common = np.array([])
            step_cum_mean = np.array([])
        else:
            step_mean = np.nanmean(step_mat, axis=0)
            step_cum_common = np.vstack([
                np.concatenate([[0.0], np.cumsum(row)]) for row in step_mat
            ])
            step_cum_mean = np.concatenate([[0.0], np.cumsum(step_mean)])

    facet_names = ["Person"] + config["facet_names"]
    facet_levels_all = {
        facet: (prep["levels"]["Person"] if facet == "Person" else prep["levels"][facet])
        for facet in facet_names
    }

    extreme_levels = {}
    for facet in facet_names:
        if facet not in obs_df.columns:
            extreme_levels[facet] = []
            continue
        stat = (
            obs_df.groupby(facet)["Observed"]
            .agg(MinScore="min", MaxScore="max")
            .reset_index()
        )
        extreme = stat[
            ((stat["MinScore"] == rating_min) & (stat["MaxScore"] == rating_min))
            | ((stat["MinScore"] == rating_max) & (stat["MaxScore"] == rating_max))
        ]
        extreme_levels[facet] = extreme[facet].astype(str).tolist()

    extreme_flags = {}
    for facet in facet_names:
        if facet in obs_df.columns:
            extreme_flags[facet] = obs_df[facet].astype(str).isin(extreme_levels[facet])
    if extreme_flags:
        extreme_flag_df = pd.DataFrame(extreme_flags)
        extreme_count = extreme_flag_df.sum(axis=1)
    else:
        extreme_count = pd.Series(0, index=obs_df.index, dtype=int)

    out = {}
    for facet in facet_names:
        if facet not in obs_df.columns:
            continue
        status_tbl = (
            obs_df.groupby(facet)["Observed"]
            .agg(MinScore="min", MaxScore="max", TotalCountAll="size")
            .reset_index()
        )
        if totalscore:
            score_source = obs_df
        else:
            flag = extreme_flags.get(facet)
            if flag is None:
                flag = pd.Series(False, index=obs_df.index)
            active_mask = (extreme_count == 0) | ((extreme_count == 1) & flag)
            score_source = obs_df.loc[active_mask]
        score_rows = []
        for level_val, df_lvl in score_source.groupby(facet):
            total_score = float(df_lvl["Observed"].sum())
            total_count = float(len(df_lvl))
            if "Weight" in df_lvl.columns:
                weightd_count = float(df_lvl["Weight"].sum())
                weightd_score = float((df_lvl["Observed"] * df_lvl["Weight"]).sum())
            else:
                weightd_count = total_count
                weightd_score = total_score
            observed_avg = weightd_score / weightd_count if weightd_count > 0 else np.nan
            score_rows.append({
                facet: level_val,
                "TotalScore": total_score,
                "TotalCount": total_count,
                "WeightdScore": weightd_score,
                "WeightdCount": weightd_count,
                "ObservedAverage": observed_avg,
            })
        score_tbl = pd.DataFrame(
            score_rows,
            columns=[facet, "TotalScore", "TotalCount", "WeightdScore", "WeightdCount", "ObservedAverage"],
        )
        level_tbl = pd.DataFrame({"Level": list(map(str, facet_levels_all[facet]))})

        score_tbl[facet] = score_tbl[facet].astype(str)
        status_tbl[facet] = status_tbl[facet].astype(str)

        tbl = level_tbl.merge(score_tbl, left_on="Level", right_on=facet, how="left")
        tbl = tbl.merge(status_tbl, left_on="Level", right_on=facet, how="left", suffixes=("", "_all"))
        tbl = tbl.drop(columns=[c for c in [facet, f"{facet}_all"] if c in tbl.columns])

        tbl["TotalScore"] = tbl["TotalScore"].fillna(0)
        tbl["TotalCount"] = tbl["TotalCount"].fillna(0)
        if "WeightdScore" in tbl.columns:
            tbl["WeightdScore"] = tbl["WeightdScore"].fillna(0)
        if "WeightdCount" in tbl.columns:
            tbl["WeightdCount"] = tbl["WeightdCount"].fillna(0)
        if "WeightdCount" in tbl.columns:
            tbl["ObservedAverage"] = tbl.apply(
                lambda r: np.nan if r["WeightdCount"] == 0 else r["ObservedAverage"], axis=1
            )
        else:
            tbl["ObservedAverage"] = tbl.apply(
                lambda r: np.nan if r["TotalCount"] == 0 else r["ObservedAverage"], axis=1
            )

        meas_tbl = measures[measures["Facet"] == facet].copy()
        if not meas_tbl.empty:
            meas_tbl["Level"] = meas_tbl["Level"].astype(str)
            tbl = tbl.merge(
                meas_tbl[
                    ["Level", "Estimate", "SE", "Infit", "Outfit", "InfitZSTD", "OutfitZSTD", "PTMEA"]
                ],
                on="Level",
                how="left",
            )

        tbl["Anchor"] = facet_anchor_status(
            facet,
            tbl["Level"].tolist(),
            config,
            extreme_levels=extreme_levels.get(facet),
        )
        tbl["Status"] = ""
        tbl.loc[tbl["TotalCountAll"].fillna(0) == 0, "Status"] = "No data"
        tbl.loc[
            (tbl["Status"] == "")
            & (tbl["MinScore"] == rating_min)
            & (tbl["MaxScore"] == rating_min),
            "Status",
        ] = "Minimum"
        tbl.loc[
            (tbl["Status"] == "")
            & (tbl["MinScore"] == rating_max)
            & (tbl["MaxScore"] == rating_max),
            "Status",
        ] = "Maximum"
        tbl.loc[(tbl["Status"] == "") & (tbl["TotalCountAll"] == 1), "Status"] = "One datum"

        sign = facet_signs.get(facet, -1)
        if facet == "Person":
            other_sum = float(np.sum([facet_signs[k] * v for k, v in facet_means.items()]))
            eta_m = tbl["Estimate"] + other_sum
            eta_z = tbl["Estimate"]
        else:
            other_sum = float(np.sum([
                facet_signs[k] * v for k, v in facet_means.items() if k != facet
            ]))
            eta_m = theta_mean + other_sum + sign * tbl["Estimate"]
            eta_z = sign * tbl["Estimate"]

        if config["model"] == "PCM" and config.get("step_facet"):
            step_levels = prep["levels"][config["step_facet"]]
            if facet == config["step_facet"] and len(step_levels) > 0 and len(step_cum_common) > 0:
                step_cum_list = []
                for lvl in tbl["Level"].tolist():
                    idx = step_levels.index(lvl) if lvl in step_levels else None
                    if idx is not None and idx < len(step_cum_common):
                        step_cum_list.append(step_cum_common[idx])
                    else:
                        step_cum_list.append(step_cum_mean)
            else:
                step_cum_list = [step_cum_mean for _ in range(len(tbl))]
        else:
            step_cum_list = [step_cum_common for _ in range(len(tbl))]

        tbl["FairM"] = [
            expected_score_from_eta(e, step, rating_min) for e, step in zip(eta_m, step_cum_list)
        ]
        tbl["FairZ"] = [
            expected_score_from_eta(e, step, rating_min) for e, step in zip(eta_z, step_cum_list)
        ]

        xtreme_target = np.where(
            tbl["Status"] == "Minimum",
            rating_min + xtreme,
            np.where(tbl["Status"] == "Maximum", rating_max - xtreme, np.nan),
        )
        xtreme_eta = [
            estimate_eta_from_target(t, step, rating_min, rating_max) if xtreme > 0 else np.nan
            for t, step in zip(xtreme_target, step_cum_list)
        ]

        measure_logit = tbl["Estimate"].copy()
        xtreme_mask = np.isfinite(xtreme_eta)
        if np.any(xtreme_mask):
            if facet == "Person":
                measure_logit = np.where(xtreme_mask, np.array(xtreme_eta) - other_sum, measure_logit)
            else:
                measure_logit = np.where(
                    xtreme_mask,
                    (np.array(xtreme_eta) - theta_mean - other_sum) / sign,
                    measure_logit,
                )

        scale_factor = uscale if np.isfinite(uscale) else 1
        scale_origin = umean if np.isfinite(umean) else 0

        tbl["Measure"] = np.where(np.isfinite(measure_logit), measure_logit * scale_factor + scale_origin, np.nan)
        tbl["ModelSE"] = np.where(np.isfinite(tbl["SE"]), np.abs(scale_factor) * tbl["SE"], np.nan)
        tbl["RealSE"] = np.where(
            np.isfinite(tbl["SE"]) & np.isfinite(tbl["Infit"]),
            np.abs(scale_factor) * tbl["SE"] * np.sqrt(np.maximum(tbl["Infit"], 1)),
            np.nan,
        )

        extreme_mask = tbl["Status"].isin(["Minimum", "Maximum"])
        tbl.loc[extreme_mask, ["Infit", "Outfit", "InfitZSTD", "OutfitZSTD"]] = np.nan
        if facet == "Person":
            tbl["PTMEA"] = np.nan

        tbl = tbl.rename(columns={
            "Infit": "InfitMnSq",
            "InfitZSTD": "InfitZStd",
            "Outfit": "OutfitMnSq",
            "OutfitZSTD": "OutfitZStd",
            "PTMEA": "PtMeaCorr",
        })

        tbl = tbl[[
            "TotalScore",
            "TotalCount",
            "WeightdScore",
            "WeightdCount",
            "ObservedAverage",
            "FairM",
            "FairZ",
            "Measure",
            "ModelSE",
            "RealSE",
            "InfitMnSq",
            "InfitZStd",
            "OutfitMnSq",
            "OutfitZStd",
            "PtMeaCorr",
            "Anchor",
            "Status",
            "Level",
            "TotalCountAll",
        ]]

        if omit_unobserved:
            tbl = tbl[tbl["TotalCountAll"] > 0]

        tbl = tbl.drop(columns=["TotalCountAll"]).sort_values(["Measure", "TotalCount"], ascending=[False, False])
        out[facet] = tbl

    return out


def calc_reliability(measure_df):
    rows = []
    if measure_df.empty:
        return pd.DataFrame()
    for facet, df in measure_df.groupby("Facet"):
        mv = np.nanvar(df["Estimate"], ddof=1)
        ev = np.nanmean(df["SE"] ** 2)
        rmse = np.sqrt(ev) if np.isfinite(ev) else np.nan
        if np.isfinite(mv) and np.isfinite(ev) and mv > 0 and ev > 0:
            tv = max(mv - ev, 0.0)
            separation = np.sqrt(tv / ev) if tv > 0 else 0.0
            reliability = tv / mv if mv > 0 else np.nan
            strata = (4 * separation + 1) / 3
        else:
            separation = np.nan
            reliability = np.nan
            strata = np.nan
        rows.append({
            "Facet": facet,
            "Levels": len(df),
            "SD": np.sqrt(mv) if np.isfinite(mv) else np.nan,
            "RMSE": rmse,
            "Separation": separation,
            "Strata": strata,
            "Reliability": reliability,
            "MeanInfit": np.nanmean(df["Infit"]),
            "MeanOutfit": np.nanmean(df["Outfit"]),
        })
    return pd.DataFrame(rows)


def calc_facets_chisq(measure_df):
    if measure_df.empty:
        return pd.DataFrame()
    rows = []
    for facet, df in measure_df.groupby("Facet"):
        levels = len(df)
        estimate = df["Estimate"].to_numpy()
        se = df["SE"].to_numpy()
        w = np.where(np.isfinite(se) & (se > 0), 1 / (se ** 2), np.nan)
        ok = np.isfinite(w) & np.isfinite(estimate)
        fixed_chi = (
            np.sum(w[ok] * estimate[ok] ** 2) - (np.sum(w[ok] * estimate[ok]) ** 2) / np.sum(w[ok])
            if np.sum(ok) >= 2
            else np.nan
        )
        fixed_df = max(levels - 1, 0)
        random_var = (
            (np.sum((estimate[ok] - np.nanmean(estimate[ok])) ** 2) / (np.sum(ok) - 1))
            - (np.sum(se[ok] ** 2) / np.sum(ok))
            if np.sum(ok) >= 2
            else np.nan
        )
        if np.isfinite(random_var) and random_var <= 0:
            random_var = np.nan
        if np.isfinite(random_var) and random_var > 0:
            w_r = np.where(np.isfinite(se) & (se > 0), 1 / (random_var + se ** 2), np.nan)
            ok_r = np.isfinite(w_r) & np.isfinite(estimate)
            random_chi = (
                np.sum(w_r[ok_r] * estimate[ok_r] ** 2) - (np.sum(w_r[ok_r] * estimate[ok_r]) ** 2) / np.sum(w_r[ok_r])
                if np.sum(ok_r) >= 2
                else np.nan
            )
        else:
            random_chi = np.nan

        random_df = max(levels - 2, 0)
        rows.append({
            "Facet": facet,
            "Levels": levels,
            "MeanMeasure": np.nanmean(estimate),
            "SD": np.nanstd(estimate, ddof=1),
            "FixedChiSq": fixed_chi,
            "FixedDF": fixed_df,
            "FixedProb": 1 - chi2.cdf(fixed_chi, df=fixed_df) if np.isfinite(fixed_chi) and fixed_df > 0 else np.nan,
            "RandomVar": random_var,
            "RandomChiSq": random_chi,
            "RandomDF": random_df,
            "RandomProb": 1 - chi2.cdf(random_chi, df=random_df) if np.isfinite(random_chi) and random_df > 0 else np.nan,
        })
    return pd.DataFrame(rows)


def ensure_positive_definite(mat, eps=1e-6):
    mat = np.array(mat, dtype=float, copy=True)
    mat = (mat + mat.T) / 2
    eigvals = np.linalg.eigvalsh(mat)
    min_eig = np.min(eigvals)
    if min_eig < eps:
        mat += np.eye(mat.shape[0]) * (eps - min_eig)
    return mat


def compute_pca_bundle(residual_matrix_wide):
    if residual_matrix_wide is None:
        return None
    if residual_matrix_wide.shape[0] < 2 or residual_matrix_wide.shape[1] < 2:
        return None

    residual_matrix_clean = residual_matrix_wide.dropna(axis=1, how="all")
    if residual_matrix_clean.shape[1] < 2:
        return None

    cor_df = residual_matrix_clean.corr(min_periods=2)
    cor_df = cor_df.fillna(0)
    np.fill_diagonal(cor_df.values, 1)

    cor_pd = ensure_positive_definite(cor_df.values)
    cor_pd_df = pd.DataFrame(cor_pd, index=cor_df.index, columns=cor_df.columns)

    eigvals, eigvecs = np.linalg.eigh(cor_pd)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    total = np.sum(eigvals)
    var_pct = (eigvals / total * 100) if total > 0 else np.zeros_like(eigvals)
    loadings = eigvecs * np.sqrt(np.maximum(eigvals, 0))
    loadings_df = pd.DataFrame(
        loadings,
        index=cor_pd_df.index,
        columns=[f"PC{i+1}" for i in range(loadings.shape[1])],
    )

    return {
        "eigenvalues": eigvals,
        "variance_pct": var_pct,
        "loadings": loadings_df,
        "cor_matrix": cor_pd_df,
        "residual_matrix": residual_matrix_wide,
    }


def compute_pca_overall(obs_df, facet_names):
    if obs_df is None or obs_df.empty or not facet_names:
        return None
    df_aug = obs_df.copy()
    df_aug["Person"] = df_aug["Person"].astype(str)
    df_aug["item_combination"] = df_aug[facet_names].astype(str).agg("_".join, axis=1)
    residual_matrix_prep = (
        df_aug.groupby(["Person", "item_combination"])["StdResidual"]
        .mean()
        .reset_index()
    )
    residual_matrix_wide = residual_matrix_prep.pivot(
        index="Person", columns="item_combination", values="StdResidual"
    )
    return compute_pca_bundle(residual_matrix_wide)


def compute_pca_by_facet(obs_df, facet_names):
    out = {}
    if obs_df is None or obs_df.empty:
        return out
    for facet in facet_names:
        df = obs_df.copy()
        df["Person"] = df["Person"].astype(str)
        df["Level"] = df[facet].astype(str)
        residual_matrix_prep = (
            df.groupby(["Person", "Level"])["StdResidual"]
            .mean()
            .reset_index()
        )
        residual_matrix_wide = residual_matrix_prep.pivot(
            index="Person", columns="Level", values="StdResidual"
        )
        out[facet] = compute_pca_bundle(residual_matrix_wide)
    return out


def calc_expected_category_counts(res):
    if res is None:
        return pd.DataFrame()
    prep = res["prep"]
    config = res["config"]
    idx = build_indices(prep, step_facet=config["step_facet"])
    params = res["params"]
    theta_hat = params["theta"] if config["method"] == "JMLE" else res["facets"]["person"]["Estimate"].to_numpy()
    eta = compute_eta(idx, params, config, theta_override=theta_hat if theta_hat.size else None)
    if config["model"] == "RSM":
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
        probs = category_prob_rsm(eta, step_cum)
    else:
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])
        probs = category_prob_pcm(eta, step_cum_mat, idx["step_idx"])
    if probs.size == 0:
        return pd.DataFrame()
    w = idx.get("weight")
    if w is None:
        exp_counts = np.nansum(probs, axis=0)
    else:
        exp_counts = np.nansum(probs * w[:, None], axis=0)
    total_exp = np.nansum(exp_counts)
    cat_vals = np.arange(prep["rating_min"], prep["rating_max"] + 1)
    return pd.DataFrame({
        "Category": cat_vals,
        "ExpectedCount": exp_counts,
        "ExpectedPercent": (100 * exp_counts / total_exp) if total_exp > 0 else np.nan,
    })


def calc_category_stats(obs_df, res=None, whexact=False):
    if obs_df.empty:
        return pd.DataFrame()
    total_n = np.nansum(get_weights(obs_df))
    rows = []
    for category, df_cat in obs_df.groupby("Observed"):
        w = get_weights(df_cat)
        count = np.nansum(w)
        var_w = np.nansum(df_cat["Var"] * w)
        infit = np.nansum(df_cat["StdSq"] * df_cat["Var"] * w) / var_w if var_w > 0 else np.nan
        outfit = np.nansum(df_cat["StdSq"] * w) / count if count > 0 else np.nan
        rows.append({
            "Category": category,
            "Count": count,
            "AvgPersonMeasure": weighted_mean(df_cat["PersonMeasure"].to_numpy(), w),
            "ExpectedAverage": weighted_mean(df_cat["Expected"].to_numpy(), w),
            "Infit": infit,
            "Outfit": outfit,
            "MeanResidual": weighted_mean(df_cat["Residual"].to_numpy(), w),
            "DF_Infit": var_w,
            "DF_Outfit": count,
        })
    summary = pd.DataFrame(rows)

    all_categories = (
        np.arange(res["prep"]["rating_min"], res["prep"]["rating_max"] + 1)
        if res is not None
        else np.sort(obs_df["Observed"].unique())
    )
    cat_tbl = pd.DataFrame({"Category": all_categories}).merge(summary, on="Category", how="left")
    cat_tbl["Count"] = cat_tbl["Count"].fillna(0)
    cat_tbl["Percent"] = np.where(total_n > 0, 100 * cat_tbl["Count"] / total_n, np.nan)
    cat_tbl["InfitZSTD"] = [
        zstd_from_mnsq(m, d, whexact=whexact) for m, d in zip(cat_tbl["Infit"], cat_tbl["DF_Infit"])
    ]
    cat_tbl["OutfitZSTD"] = [
        zstd_from_mnsq(m, d, whexact=whexact) for m, d in zip(cat_tbl["Outfit"], cat_tbl["DF_Outfit"])
    ]

    exp_tbl = calc_expected_category_counts(res) if res is not None else pd.DataFrame()
    if not exp_tbl.empty:
        cat_tbl = cat_tbl.merge(exp_tbl, on="Category", how="left")
        cat_tbl["DiffCount"] = cat_tbl["Count"] - cat_tbl["ExpectedCount"]
        cat_tbl["DiffPercent"] = cat_tbl["Percent"] - cat_tbl["ExpectedPercent"]

    cat_tbl["LowCount"] = cat_tbl["Count"] < 10
    cat_tbl["InfitFlag"] = cat_tbl["Infit"].apply(
        lambda v: np.nan if pd.isna(v) else (v < 0.5 or v > 1.5)
    )
    cat_tbl["OutfitFlag"] = cat_tbl["Outfit"].apply(
        lambda v: np.nan if pd.isna(v) else (v < 0.5 or v > 1.5)
    )
    cat_tbl["ZSTDFlag"] = (
        (cat_tbl["InfitZSTD"].abs() >= 2) | (cat_tbl["OutfitZSTD"].abs() >= 2)
    )
    return cat_tbl.sort_values("Category")


def calc_step_order(step_tbl):
    if step_tbl is None or step_tbl.empty:
        return pd.DataFrame()
    step_tbl = step_tbl.copy()
    step_tbl["StepIndex"] = step_tbl["Step"].str.extract(r"(\d+)").astype(float)
    if "StepFacet" not in step_tbl.columns:
        step_tbl["StepFacet"] = "Common"
    step_tbl = step_tbl.sort_values(["StepFacet", "StepIndex"])
    step_tbl["Spacing"] = step_tbl.groupby("StepFacet")["Estimate"].diff()
    step_tbl["Ordered"] = step_tbl["Spacing"].apply(lambda x: np.nan if pd.isna(x) else x > 0)
    return step_tbl


def category_warnings_text(cat_tbl, step_tbl=None):
    if cat_tbl is None or cat_tbl.empty:
        return "No category diagnostics available."
    msgs = []
    unused = cat_tbl[cat_tbl["Count"] == 0]
    if not unused.empty:
        msgs.append("Unused categories: " + ", ".join(map(str, unused["Category"].tolist())))
    low_counts = cat_tbl[cat_tbl["Count"] < 10]
    if not low_counts.empty:
        msgs.append("Low category counts (<10): " + ", ".join(map(str, low_counts["Category"].tolist())))
    if "DiffPercent" in cat_tbl.columns:
        diff_bad = cat_tbl[cat_tbl["DiffPercent"].abs() >= 5]
        if not diff_bad.empty:
            msgs.append("Observed vs expected % differs by >= 5: " + ", ".join(map(str, diff_bad["Category"].tolist())))
    if "InfitZSTD" in cat_tbl.columns and "OutfitZSTD" in cat_tbl.columns:
        zstd_bad = cat_tbl[(cat_tbl["InfitZSTD"].abs() >= 2) | (cat_tbl["OutfitZSTD"].abs() >= 2)]
        if not zstd_bad.empty:
            msgs.append("Category |ZSTD| >= 2: " + ", ".join(map(str, zstd_bad["Category"].tolist())))
    avg_tbl = cat_tbl.dropna(subset=["AvgPersonMeasure"]).sort_values("Category")
    if len(avg_tbl) >= 3 and not avg_tbl["AvgPersonMeasure"].is_monotonic_increasing:
        msgs.append("Category averages are not monotonic (Avg Measure by category).")
    if step_tbl is not None and not step_tbl.empty:
        disordered = step_tbl[step_tbl["Ordered"] == False]
        if not disordered.empty:
            labels = disordered.apply(lambda r: f"{r['StepFacet']}:{r['Step']}", axis=1).tolist()
            msgs.append("Disordered thresholds detected: " + ", ".join(labels))
    return "No major category warnings detected." if not msgs else "\n".join(msgs)


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def fmt_count(value):
    val = to_float(value)
    if not np.isfinite(val):
        return "NA"
    if abs(val - round(val)) < 1e-6:
        return str(int(round(val)))
    return f"{val:.0f}"


def fmt_num(value, decimals=2):
    val = to_float(value)
    if not np.isfinite(val):
        return "NA"
    return f"{val:.{decimals}f}"


def fmt_pvalue(value):
    val = to_float(value)
    if not np.isfinite(val):
        return "NA"
    if val < 0.001:
        return "< .001"
    return f"= {val:.3f}"


def describe_series(series):
    if series is None:
        return None
    arr = pd.to_numeric(series, errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return None
    return {
        "min": float(np.nanmin(arr)),
        "max": float(np.nanmax(arr)),
        "mean": float(np.nanmean(arr)),
        "sd": float(np.nanstd(arr, ddof=1)) if arr.size > 1 else np.nan,
    }


def build_apa_report_text(res, diagnostics, bias_results=None, context=None, whexact=False):
    context = context or {}
    summary = res["summary"].iloc[0] if "summary" in res and not res["summary"].empty else {}
    prep = res["prep"]
    config = res["config"]

    n_obs = to_float(summary.get("N", np.nan))
    n_person = to_float(summary.get("Persons", len(res["facets"]["person"])))
    n_cat = to_float(summary.get("Categories", config.get("n_cat", np.nan)))
    rating_min = to_float(prep.get("rating_min", np.nan))
    rating_max = to_float(prep.get("rating_max", np.nan))

    facet_names = list(config.get("facet_names", []))
    facet_levels = config.get("facet_levels", {})
    facet_counts = {f: len(facet_levels.get(f, [])) for f in facet_names}
    facets_text = ", ".join([f"{f} (n = {fmt_count(n)})" for f, n in facet_counts.items()]) if facet_counts else "no additional facets"

    assessment = context.get("assessment", "").strip()
    setting = context.get("setting", "").strip()
    rater_training = context.get("rater_training", "").strip()
    raters_per_response = context.get("raters_per_response", "").strip()
    scale_desc = context.get("scale_desc", "").strip()

    method_sentences = []
    if assessment:
        if setting:
            method_sentences.append(f"The analysis focused on {assessment} in {setting}.")
        else:
            method_sentences.append(f"The analysis focused on {assessment}.")
    method_sentences.append(
        "A many-facet Rasch model (MFRM) was fit to "
        f"{fmt_count(n_obs)} observations from {fmt_count(n_person)} persons scored on a "
        f"{fmt_count(n_cat)}-category scale ({fmt_count(rating_min)}-{fmt_count(rating_max)})."
    )
    if facet_names:
        method_sentences.append(f"The design included facets for {facets_text}.")
    else:
        method_sentences.append("No additional facets beyond Person were modeled.")
    if scale_desc:
        method_sentences.append(f"The rating scale was described as {scale_desc}.")
    if rater_training:
        method_sentences.append(f"Raters received {rater_training}.")
    if raters_per_response:
        method_sentences.append(f"Each response was scored by {raters_per_response} raters on average.")

    model = config.get("model", "RSM")
    method = config.get("method", "JMLE")
    model_sentence = f"The {model} specification was estimated using {method} in the Streamlit MFRM app."
    if model == "PCM" and config.get("step_facet"):
        model_sentence += f" The step structure varied by {config['step_facet']}."
    method_sentences.append(model_sentence)
    if config.get("weight_col"):
        method_sentences.append("Observation weights were applied as frequency counts.")

    anchor_summary = config.get("anchor_summary", pd.DataFrame())
    if isinstance(anchor_summary, pd.DataFrame) and not anchor_summary.empty:
        anchored_levels = int(anchor_summary["AnchoredLevels"].sum())
        group_anchors = int(anchor_summary["GroupAnchors"].sum())
        dummy_facets = anchor_summary.loc[anchor_summary["DummyFacet"], "Facet"].tolist()
        if anchored_levels > 0 or group_anchors > 0:
            method_sentences.append(
                f"Anchoring was used (anchored levels = {fmt_count(anchored_levels)}, "
                f"group anchors = {fmt_count(group_anchors)})."
            )
        if dummy_facets:
            method_sentences.append(f"Dummy-facet constraints were applied to: {', '.join(dummy_facets)}.")
    noncenter = config.get("noncenter_facet")
    if noncenter:
        method_sentences.append(f"The {noncenter} facet was specified as non-centered.")

    method_text = "Method.\n" + " ".join(method_sentences)

    results_sentences = []
    cat_tbl = calc_category_stats(diagnostics["obs"], res=res, whexact=whexact) if diagnostics else pd.DataFrame()
    step_order = calc_step_order(res.get("steps")) if res else pd.DataFrame()
    unused = int((cat_tbl["Count"] == 0).sum()) if not cat_tbl.empty else 0
    low_count = int((cat_tbl["Count"] < 10).sum()) if not cat_tbl.empty else 0
    disordered = step_order[step_order["Ordered"] == False] if not step_order.empty else pd.DataFrame()
    usage_label = "adequate" if unused == 0 and low_count == 0 else "uneven"
    threshold_text = "thresholds were ordered" if disordered.empty else (
        f"thresholds were disordered for {fmt_count(len(disordered))} step(s)"
    )
    results_sentences.append(
        f"Category usage was {usage_label} (unused categories = {fmt_count(unused)}, "
        f"low-count categories = {fmt_count(low_count)}), and {threshold_text}."
    )

    person_stats = describe_series(res["facets"]["person"]["Estimate"]) if res else None
    if person_stats:
        results_sentences.append(
            "Person measures ranged from "
            f"{fmt_num(person_stats['min'])} to {fmt_num(person_stats['max'])} logits "
            f"(M = {fmt_num(person_stats['mean'])}, SD = {fmt_num(person_stats['sd'])})."
        )

    facet_texts = []
    if res and not res["facets"]["others"].empty:
        for facet in facet_names:
            df_f = res["facets"]["others"][res["facets"]["others"]["Facet"] == facet]
            stats = describe_series(df_f["Estimate"]) if not df_f.empty else None
            if stats:
                facet_texts.append(
                    f"{facet} measures ranged from {fmt_num(stats['min'])} to "
                    f"{fmt_num(stats['max'])} logits (M = {fmt_num(stats['mean'])}, SD = {fmt_num(stats['sd'])})."
                )
    results_sentences.extend(facet_texts)

    overall_fit = diagnostics["overall_fit"].iloc[0] if diagnostics and not diagnostics["overall_fit"].empty else None
    if overall_fit is not None:
        infit = float(overall_fit.get("Infit", np.nan))
        outfit = float(overall_fit.get("Outfit", np.nan))
        fit_label = "acceptable" if 0.5 <= infit <= 1.5 and 0.5 <= outfit <= 1.5 else "elevated"
        results_sentences.append(
            f"Overall fit was {fit_label} (infit MnSq = {fmt_num(infit)}, "
            f"outfit MnSq = {fmt_num(outfit)})."
        )
    fit_tbl = diagnostics["fit"] if diagnostics else pd.DataFrame()
    if not fit_tbl.empty:
        misfit = (
            (fit_tbl["Infit"] < 0.5) | (fit_tbl["Infit"] > 1.5) |
            (fit_tbl["Outfit"] < 0.5) | (fit_tbl["Outfit"] > 1.5)
        )
        results_sentences.append(
            f"{fmt_count(int(misfit.sum()))} of {fmt_count(len(fit_tbl))} elements exceeded the 0.5-1.5 fit range."
        )

    rel_tbl = diagnostics["reliability"] if diagnostics else pd.DataFrame()
    if not rel_tbl.empty:
        rel_texts = []
        for _, row in rel_tbl.iterrows():
            rel_texts.append(
                f"{row['Facet']} reliability = {fmt_num(row['Reliability'])} "
                f"(separation = {fmt_num(row['Separation'])})."
            )
        results_sentences.append(" ".join(rel_texts))

    if bias_results and bias_results.get("table") is not None and not bias_results["table"].empty:
        bias_tbl = bias_results["table"].copy()
        bias_tbl = bias_tbl[np.isfinite(bias_tbl["t"])]
        if not bias_tbl.empty:
            idx = bias_tbl["t"].abs().idxmax()
            row = bias_tbl.loc[idx]
            results_sentences.append(
                f"Bias analysis for {bias_results['facet_a']} x {bias_results['facet_b']} "
                f"showed a largest contrast of {fmt_num(row['Bias Size'])} logits "
                f"(t = {fmt_num(row['t'])}, p {fmt_pvalue(row['Prob.'])})."
            )
    else:
        results_sentences.append("Bias analysis was not estimated in this run.")

    results_text = "Results.\n" + " ".join(results_sentences)
    return method_text + "\n\n" + results_text


def build_apa_table_figure_note_map(res, diagnostics, bias_results=None, context=None, whexact=False):
    context = context or {}
    summary = res["summary"].iloc[0] if "summary" in res and not res["summary"].empty else {}
    prep = res["prep"]
    config = res["config"]

    n_obs = to_float(summary.get("N", np.nan))
    n_person = to_float(summary.get("Persons", len(res["facets"]["person"])))
    n_cat = to_float(summary.get("Categories", config.get("n_cat", np.nan)))
    rating_min = to_float(prep.get("rating_min", np.nan))
    rating_max = to_float(prep.get("rating_max", np.nan))
    model = config.get("model", "RSM")
    method = config.get("method", "JMLE")
    rater_facet = context.get("rater_facet", "").strip()

    cat_tbl = calc_category_stats(diagnostics["obs"], res=res, whexact=whexact) if diagnostics else pd.DataFrame()
    step_order = calc_step_order(res.get("steps")) if res else pd.DataFrame()
    unused = int((cat_tbl["Count"] == 0).sum()) if not cat_tbl.empty else 0
    low_count = int((cat_tbl["Count"] < 10).sum()) if not cat_tbl.empty else 0
    disordered = step_order[step_order["Ordered"] == False] if not step_order.empty else pd.DataFrame()
    threshold_text = "ordered" if disordered.empty else f"disordered in {fmt_count(len(disordered))} step(s)"

    overall_fit = diagnostics["overall_fit"].iloc[0] if diagnostics and not diagnostics["overall_fit"].empty else None
    infit = float(overall_fit.get("Infit", np.nan)) if overall_fit is not None else np.nan
    outfit = float(overall_fit.get("Outfit", np.nan)) if overall_fit is not None else np.nan

    rel_tbl = diagnostics["reliability"] if diagnostics else pd.DataFrame()
    rater_rel = None
    if rater_facet and not rel_tbl.empty:
        match = rel_tbl[rel_tbl["Facet"] == rater_facet]
        if not match.empty:
            rater_rel = match.iloc[0]

    note_map = {}
    note_map["table1"] = (
        "Table 1. Facet summary\n"
        "Note. Measures are reported in logits; higher values indicate more of the modeled trait for that facet. "
        "SE = standard error; MnSq = mean-square fit. "
        f"Model = {model}; estimation = {method}; N = {fmt_count(n_obs)} observations "
        f"from {fmt_count(n_person)} persons on a {fmt_count(n_cat)}-category scale "
        f"({fmt_count(rating_min)}-{fmt_count(rating_max)})."
    )
    note_map["table2"] = (
        "Table 2. Rating scale diagnostics\n"
        "Note. Category counts and thresholds summarize scale functioning. "
        f"Thresholds were {threshold_text}; unused categories = {fmt_count(unused)}; "
        f"low-count categories (< 10) = {fmt_count(low_count)}."
    )
    fit_sentence = (
        f"Overall fit: infit MnSq = {fmt_num(infit)}, outfit MnSq = {fmt_num(outfit)}."
        if np.isfinite(infit) and np.isfinite(outfit)
        else "Overall fit indices are reported as mean infit and outfit MnSq."
    )
    note_map["table3"] = (
        "Table 3. Fit and reliability summary\n"
        "Note. Separation and reliability are based on observed variance and measurement error. "
        + fit_sentence
    )
    if rater_rel is not None:
        note_map["table3"] += (
            f" Rater facet ({rater_facet}) reliability = {fmt_num(rater_rel['Reliability'])}, "
            f"separation = {fmt_num(rater_rel['Separation'])}."
        )

    if bias_results and bias_results.get("table") is not None and not bias_results["table"].empty:
        bias_tbl = bias_results["table"]
        sig = bias_tbl[pd.to_numeric(bias_tbl["Prob."], errors="coerce") < 0.05]
        note_map["table4"] = (
            "Table 4. Bias/interaction effects\n"
            "Note. Bias contrasts are in logits and represent observed minus expected scores with main effects held fixed. "
            f"Significant interactions (p < .05) = {fmt_count(len(sig))}."
        )
    else:
        note_map["table4"] = (
            "Table 4. Bias/interaction effects\n"
            "Note. Bias contrasts are in logits and represent observed minus expected scores with main effects held fixed."
        )

    note_map["figure1"] = (
        "Figure 1. Wright map\n"
        "Note. Persons and facet elements are located on a common logit scale; higher values indicate higher ability or "
        "greater severity/difficulty depending on facet orientation."
    )
    note_map["figure2"] = (
        "Figure 2. Pathway map (measure vs. infit t)\n"
        "Note. Points show element measures and their standardized infit values. Extreme |ZSTD| values flag potential misfit."
    )
    note_map["figure3"] = (
        "Figure 3. Facet estimate distribution\n"
        "Note. Distributions summarize severity/difficulty spread within each facet."
    )
    note_map["figure4"] = (
        "Figure 4. Step/threshold estimates\n"
        "Note. Step ordering should generally increase; disordered thresholds suggest category structure issues."
    )
    note_map["figure5"] = (
        "Figure 5. Category probability curves\n"
        "Note. Curves show the most probable category across the latent continuum; well-functioning categories show "
        "distinct peaks in order."
    )
    note_map["figure6"] = (
        "Figure 6. Observed vs expected scores\n"
        "Note. Points summarize mean observed and expected scores by bin; deviations from the diagonal suggest local misfit."
    )
    note_map["figure7"] = (
        "Figure 7. Fit diagnostics (Infit vs Outfit)\n"
        "Note. Each point represents an element within a facet. Values near 1.0 indicate expected fit; "
        "values substantially above 1.0 suggest misfit."
    )
    note_map["figure8"] = (
        "Figure 8. Fit ZSTD distribution\n"
        "Note. Distributions of standardized fit help identify unusually large residuals across facets."
    )
    note_map["figure9"] = (
        "Figure 9. Misfit levels\n"
        "Note. Levels are ranked by maximum |ZSTD| to highlight potentially problematic elements."
    )

    return note_map


def build_apa_table_figure_notes(res, diagnostics, bias_results=None, context=None, whexact=False):
    note_map = build_apa_table_figure_note_map(
        res,
        diagnostics,
        bias_results=bias_results,
        context=context,
        whexact=whexact,
    )
    ordered_keys = [
        "table1",
        "table2",
        "table3",
        "table4",
        "figure1",
        "figure2",
        "figure3",
        "figure4",
        "figure5",
        "figure6",
        "figure7",
        "figure8",
        "figure9",
    ]
    return "\n\n".join([note_map[k] for k in ordered_keys if k in note_map])


def build_visual_warning_map(res, diagnostics, whexact=False, thresholds=None):
    warnings = {f"figure{i}": [] for i in range(1, 10)}
    if res is None or diagnostics is None:
        return warnings
    thresholds = thresholds or {}
    n_obs_min = thresholds.get("n_obs_min", 100)
    n_person_min = thresholds.get("n_person_min", 30)
    low_cat_min = thresholds.get("low_cat_min", 10)
    min_facet_levels = thresholds.get("min_facet_levels", 3)
    misfit_ratio_warn = thresholds.get("misfit_ratio_warn", 0.1)
    missing_fit_ratio_warn = thresholds.get("missing_fit_ratio_warn", 0.2)
    zstd2_ratio_warn = thresholds.get("zstd2_ratio_warn", 0.1)
    zstd3_ratio_warn = thresholds.get("zstd3_ratio_warn", 0.05)
    expected_var_min = thresholds.get("expected_var_min", 0.2)

    summary = res["summary"].iloc[0] if "summary" in res and not res["summary"].empty else {}
    n_obs = to_float(summary.get("N", np.nan))
    n_person = len(res["facets"]["person"]) if "facets" in res else 0
    if np.isfinite(n_obs) and n_obs < n_obs_min:
        warnings["figure1"].append(
            f"Small number of observations (N = {fmt_count(n_obs)} < {fmt_count(n_obs_min)})."
        )
        warnings["figure2"].append(
            f"Small number of observations (N = {fmt_count(n_obs)} < {fmt_count(n_obs_min)}); ZSTD values may be volatile."
        )
        warnings["figure6"].append(
            f"Small number of observations (N = {fmt_count(n_obs)} < {fmt_count(n_obs_min)}); bin averages may be noisy."
        )
    if n_person < n_person_min:
        warnings["figure1"].append(
            f"Small person sample (n = {fmt_count(n_person)} < {fmt_count(n_person_min)}); interpret spread cautiously."
        )

    facet_levels = res.get("config", {}).get("facet_levels", {})
    small_facets = [f for f, levels in facet_levels.items() if len(levels) < min_facet_levels]
    if small_facets:
        warnings["figure1"].append(
            "Facets with very few levels: " + ", ".join(small_facets) + "."
        )
        warnings["figure3"].append(
            "Facet distributions are based on few levels: " + ", ".join(small_facets) + "."
        )

    cat_tbl = calc_category_stats(diagnostics["obs"], res=res, whexact=whexact) if diagnostics else pd.DataFrame()
    if not cat_tbl.empty:
        unused = int((cat_tbl["Count"] == 0).sum())
        low_count = int((cat_tbl["Count"] < low_cat_min).sum())
        if unused > 0:
            warnings["figure5"].append(f"Unused categories detected (n = {fmt_count(unused)}).")
        if low_count > 0:
            warnings["figure5"].append(
                f"Low-count categories (< {fmt_count(low_cat_min)}) detected (n = {fmt_count(low_count)})."
            )

    step_order = calc_step_order(res.get("steps")) if res else pd.DataFrame()
    if not step_order.empty:
        disordered = step_order[step_order["Ordered"] == False]
        if not disordered.empty:
            warnings["figure4"].append(f"Disordered thresholds detected (n = {fmt_count(len(disordered))}).")
            warnings["figure5"].append("Disordered thresholds can distort category curves.")

    measures = diagnostics.get("measures", pd.DataFrame())
    if measures.empty:
        warnings["figure2"].append("Fit statistics are not available for this run.")
        warnings["figure7"].append("Fit statistics are not available for this run.")
        warnings["figure8"].append("ZSTD distributions are not available for this run.")
        warnings["figure9"].append("Misfit ranking requires fit statistics.")
        return warnings

    infit = pd.to_numeric(measures.get("Infit"), errors="coerce")
    outfit = pd.to_numeric(measures.get("Outfit"), errors="coerce")
    infit_z = pd.to_numeric(measures.get("InfitZSTD"), errors="coerce")
    outfit_z = pd.to_numeric(measures.get("OutfitZSTD"), errors="coerce")

    valid_fit = np.isfinite(infit) & np.isfinite(outfit)
    if valid_fit.mean() < 0.8:
        missing_pct = 100 * (1 - valid_fit.mean())
        if missing_pct / 100 >= missing_fit_ratio_warn:
            warnings["figure7"].append(
                f"Fit statistics missing for {missing_pct:.0f}% of elements."
            )

    misfit = (infit < 0.5) | (infit > 1.5) | (outfit < 0.5) | (outfit > 1.5)
    if np.isfinite(misfit).any():
        ratio = np.nanmean(misfit)
        if ratio > misfit_ratio_warn:
            warnings["figure7"].append(
                f"High proportion of misfit elements ({ratio * 100:.0f}%)."
            )

    zstd = np.nanmax(np.vstack([np.abs(infit_z), np.abs(outfit_z)]), axis=0)
    zstd = zstd[np.isfinite(zstd)]
    if zstd.size > 0:
        prop2 = np.mean(zstd >= 2)
        prop3 = np.mean(zstd >= 3)
        if prop2 > zstd2_ratio_warn:
            warnings["figure8"].append(
                f"Large share of |ZSTD| >= 2 ({prop2 * 100:.0f}%)."
            )
        if prop3 > zstd3_ratio_warn:
            warnings["figure8"].append(
                f"Notable |ZSTD| >= 3 ({prop3 * 100:.0f}%)."
            )

    obs = diagnostics.get("obs", pd.DataFrame())
    if not obs.empty and "Expected" in obs.columns:
        exp_var = np.nanvar(pd.to_numeric(obs["Expected"], errors="coerce"))
        if np.isfinite(exp_var) and exp_var < expected_var_min:
            warnings["figure6"].append(
                "Expected scores have limited spread; trends may be muted."
            )

    return warnings


def build_visual_summary_map(res, diagnostics, whexact=False, options=None):
    options = options or {}
    detail = str(options.get("detail", "standard")).lower()
    max_facet_ranges = int(options.get("max_facet_ranges", 4))
    top_misfit_n = int(options.get("top_misfit_n", 3))
    include_top_misfit = top_misfit_n > 0
    summaries = {f"figure{i}": [] for i in range(1, 10)}
    if res is None or diagnostics is None:
        return summaries

    summary = res["summary"].iloc[0] if "summary" in res and not res["summary"].empty else {}
    n_obs = to_float(summary.get("N", np.nan))
    n_person = len(res["facets"]["person"]) if "facets" in res else 0
    if np.isfinite(n_obs):
        summaries["figure1"].append(f"Observations: N = {fmt_count(n_obs)}.")
    summaries["figure1"].append(f"Persons: n = {fmt_count(n_person)}.")

    person_stats = describe_series(res["facets"]["person"]["Estimate"])
    if person_stats:
        summaries["figure1"].append(
            "Person range "
            f"{fmt_num(person_stats['min'])} to {fmt_num(person_stats['max'])} "
            f"(M = {fmt_num(person_stats['mean'])}, SD = {fmt_num(person_stats['sd'])})."
        )
        if detail == "detailed" and np.isfinite(person_stats.get("sd", np.nan)):
            summaries["figure1"].append(
                f"Person spread (SD) = {fmt_num(person_stats['sd'])} logits."
            )

    if not res["facets"]["others"].empty:
        facet_stats = []
        for facet, df in res["facets"]["others"].groupby("Facet"):
            stats = describe_series(df["Estimate"])
            if stats:
                facet_stats.append(
                    f"{facet}: n = {fmt_count(len(df))}, range {fmt_num(stats['min'])} to {fmt_num(stats['max'])}"
                )
        if facet_stats:
            summaries["figure1"].append("Facet ranges: " + "; ".join(facet_stats[:max_facet_ranges]) + ".")
            if len(facet_stats) > max_facet_ranges:
                summaries["figure1"].append("Additional facets omitted for brevity.")

    measures = diagnostics.get("measures", pd.DataFrame())
    if not measures.empty:
        infit_z = pd.to_numeric(measures.get("InfitZSTD"), errors="coerce")
        valid = measures[np.isfinite(infit_z)]
        if not valid.empty:
            count2 = int(np.sum(np.abs(infit_z) >= 2))
            summaries["figure2"].append(
                f"Elements with |Infit ZSTD| >= 2: {fmt_count(count2)} of {fmt_count(len(valid))}."
            )
            if include_top_misfit:
                top = valid.assign(Abs=np.abs(infit_z)).sort_values("Abs", ascending=False).head(top_misfit_n)
                if not top.empty:
                    labels = [
                        f"{r['Facet']}: {truncate_label(r['Level'], 20)} (|Z|={fmt_num(r['Abs'])})"
                        for _, r in top.iterrows()
                    ]
                    summaries["figure2"].append("Largest |Z|: " + "; ".join(labels) + ".")

    if not res["facets"]["others"].empty:
        summaries["figure3"].append("Distributions show the spread of severity/difficulty within each facet.")
        if detail == "detailed":
            facet_ranges = []
            for facet, df in res["facets"]["others"].groupby("Facet"):
                stats = describe_series(df["Estimate"])
                if stats:
                    facet_ranges.append((facet, stats["max"] - stats["min"]))
            if facet_ranges:
                facet_ranges.sort(key=lambda x: x[1], reverse=True)
                top_facet, top_range = facet_ranges[0]
                summaries["figure3"].append(
                    f"Widest facet spread: {top_facet} (range = {fmt_num(top_range)})."
                )

    step_tbl = res.get("steps", pd.DataFrame())
    if step_tbl is not None and not step_tbl.empty:
        step_count = len(step_tbl)
        step_order = calc_step_order(step_tbl)
        disordered = step_order[step_order["Ordered"] == False] if not step_order.empty else pd.DataFrame()
        summaries["figure4"].append(f"Steps estimated: {fmt_count(step_count)}.")
        summaries["figure4"].append(f"Disordered steps: {fmt_count(len(disordered))}.")
        if detail == "detailed" and not step_order.empty:
            spacing = step_order["Spacing"].dropna()
            if not spacing.empty:
                summaries["figure4"].append(
                    f"Average step spacing: {fmt_num(float(np.nanmean(spacing)))}."
                )

    cat_tbl = calc_category_stats(diagnostics["obs"], res=res, whexact=whexact) if diagnostics else pd.DataFrame()
    if not cat_tbl.empty:
        used = int(np.sum(cat_tbl["Count"] > 0))
        total = len(cat_tbl)
        max_pct = np.nanmax(cat_tbl["Percent"]) if "Percent" in cat_tbl.columns else np.nan
        summaries["figure5"].append(
            f"Categories used: {fmt_count(used)} of {fmt_count(total)}."
        )
        if np.isfinite(max_pct):
            summaries["figure5"].append(f"Largest category share: {fmt_num(max_pct, 1)}%.")
        if detail == "detailed":
            unused = int(np.sum(cat_tbl["Count"] == 0))
            if unused > 0:
                summaries["figure5"].append(f"Unused categories: {fmt_count(unused)}.")

    obs = diagnostics.get("obs", pd.DataFrame())
    if not obs.empty and "Observed" in obs.columns and "Expected" in obs.columns:
        resid = pd.to_numeric(obs["Observed"], errors="coerce") - pd.to_numeric(obs["Expected"], errors="coerce")
        if "Weight" in obs.columns:
            w = pd.to_numeric(obs["Weight"], errors="coerce").fillna(0).to_numpy(dtype=float)
            w = np.where(w > 0, w, 0.0)
            mean_resid = np.nansum(resid * w) / np.nansum(w) if np.nansum(w) > 0 else np.nan
            mae = np.nansum(np.abs(resid) * w) / np.nansum(w) if np.nansum(w) > 0 else np.nan
        else:
            mean_resid = np.nanmean(resid)
            mae = np.nanmean(np.abs(resid))
        summaries["figure6"].append(f"Mean residual: {fmt_num(mean_resid)}.")
        summaries["figure6"].append(f"Mean absolute residual: {fmt_num(mae)}.")
        if detail == "detailed":
            obs_num = pd.to_numeric(obs["Observed"], errors="coerce")
            exp_num = pd.to_numeric(obs["Expected"], errors="coerce")
            ok = np.isfinite(obs_num) & np.isfinite(exp_num)
            if np.sum(ok) > 1:
                corr = np.corrcoef(obs_num[ok], exp_num[ok])[0, 1]
                summaries["figure6"].append(f"Observed-expected correlation: {fmt_num(corr)}.")

    if not measures.empty:
        infit = pd.to_numeric(measures.get("Infit"), errors="coerce")
        outfit = pd.to_numeric(measures.get("Outfit"), errors="coerce")
        ok = np.isfinite(infit) & np.isfinite(outfit)
        if ok.any():
            misfit = (infit < 0.5) | (infit > 1.5) | (outfit < 0.5) | (outfit > 1.5)
            summaries["figure7"].append(
                f"Misfit elements (0.5-1.5 rule): {fmt_count(np.nansum(misfit))} of {fmt_count(np.nansum(ok))}."
            )
            if detail == "detailed":
                summaries["figure7"].append(
                    f"Mean infit = {fmt_num(np.nanmean(infit))}, mean outfit = {fmt_num(np.nanmean(outfit))}."
                )

        zstd_full = np.nanmax(
            np.vstack([np.abs(pd.to_numeric(measures.get("InfitZSTD"), errors="coerce")),
                       np.abs(pd.to_numeric(measures.get("OutfitZSTD"), errors="coerce"))]),
            axis=0,
        )
        zstd_valid = zstd_full[np.isfinite(zstd_full)]
        if zstd_valid.size > 0:
            summaries["figure8"].append(f"|ZSTD| >= 2: {fmt_count(np.sum(zstd_valid >= 2))}.")
            summaries["figure8"].append(f"|ZSTD| >= 3: {fmt_count(np.sum(zstd_valid >= 3))}.")

        if zstd_valid.size > 0 and include_top_misfit:
            tmp = measures.copy()
            tmp["AbsZSTD"] = zstd_full
            top = tmp[np.isfinite(tmp["AbsZSTD"])].sort_values("AbsZSTD", ascending=False).head(top_misfit_n)
            if not top.empty:
                labels = [
                    f"{r['Facet']}: {truncate_label(r['Level'], 20)} (|Z|={fmt_num(r['AbsZSTD'])})"
                    for _, r in top.iterrows()
                ]
                summaries["figure9"].append("Top misfit: " + "; ".join(labels) + ".")

    return summaries


def render_note_text(note_text):
    if not note_text:
        return
    parts = note_text.split("\n", 1)
    title = parts[0].strip()
    body = parts[1].strip() if len(parts) > 1 else ""
    if title:
        st.markdown(f"**{title}**")
    if body:
        st.markdown(body)


CONTEXT_NOTES = {
    "data": [
        "One row per observation; required columns are Person, Score, and at least two facets.",
        "Scores must be ordered categories; use Keep original values if categories are non-contiguous.",
    ],
    "person": [
        "Person measures are logits; higher values indicate higher performance.",
        "Extreme-only scores can yield unstable measures; interpret with caution.",
    ],
    "facets": [
        "Higher rater measures indicate stricter severity; higher task measures indicate greater difficulty.",
        "Compare facet ranges to detect uneven severity/difficulty.",
    ],
    "report": [
        "Measures and fit are reported in logits; use Umean/Uscale to present in user units.",
        "Omit unobserved hides elements with zero observations.",
    ],
    "fit": [
        "Infit/Outfit near 1.0 indicate expected fit; > 1.0 suggests misfit; < 1.0 suggests overfit.",
        "ZSTD can inflate with large samples; use as a flag, not a verdict.",
    ],
    "reliability": [
        "Separation reliability reflects how well elements are distinguished within a facet.",
        "For raters, lower reliability can indicate more consistent severity across raters.",
        "Agreement compares observed vs expected matches under identical conditions.",
    ],
    "bias": [
        "Bias/interaction captures systematic deviations between facet pairs.",
        "Interpret Bias Size alongside t/p and Obs-Exp averages for practical impact.",
    ],
    "steps": [
        "Ordered thresholds are expected; disordered steps often indicate category issues.",
        "Sparse categories can destabilize step estimates.",
    ],
    "categories": [
        "Low or unused categories can distort thresholds and curves.",
        "Category probability curves should show distinct, ordered peaks.",
    ],
    "dimensionality": [
        "Large PC1 variance or low PC1/PC2 ratio can indicate multidimensionality.",
    ],
    "subsets": [
        "Disconnected subsets cannot be compared without anchors or linking data.",
    ],
}


def render_context_notes(key):
    bullets = CONTEXT_NOTES.get(key, [])
    if not bullets:
        return
    st.markdown("**Interpretation notes**")
    st.markdown("\n".join([f"- {b}" for b in bullets]))


def build_apa_table_figure_captions(res, diagnostics, bias_results=None, context=None):
    context = context or {}
    assessment = context.get("assessment", "").strip()
    facet_pair = ""
    if bias_results and bias_results.get("facet_a") and bias_results.get("facet_b"):
        facet_pair = f"{bias_results['facet_a']} x {bias_results['facet_b']}"

    assessment_phrase = f" for {assessment}" if assessment else ""

    blocks = []
    blocks.append(
        "Table 1\n"
        f"Facet Summary (Measures, SE, Fit, Reliability){assessment_phrase}"
    )
    blocks.append(
        "Table 2\n"
        "Rating Scale Diagnostics (Category Counts and Thresholds)"
    )
    blocks.append(
        "Table 3\n"
        "Fit and Reliability Summary"
    )
    if facet_pair:
        blocks.append(
            "Table 4\n"
            f"Bias/Interaction Effects for {facet_pair}"
        )
    else:
        blocks.append(
            "Table 4\n"
            "Bias/Interaction Effects"
        )

    blocks.append(
        "Figure 1\n"
        f"Wright Map of Person and Facet Measures{assessment_phrase}"
    )
    blocks.append(
        "Figure 2\n"
        "Pathway Map (Measure vs. Infit t)"
    )
    blocks.append(
        "Figure 3\n"
        "Facet Estimate Distribution"
    )
    blocks.append(
        "Figure 4\n"
        "Step/Threshold Estimates"
    )
    blocks.append(
        "Figure 5\n"
        "Category Probability Curves"
    )
    blocks.append(
        "Figure 6\n"
        "Observed vs. Expected Scores"
    )
    blocks.append(
        "Figure 7\n"
        "Fit Diagnostics (Infit vs Outfit)"
    )
    blocks.append(
        "Figure 8\n"
        "Fit ZSTD Distribution"
    )
    blocks.append(
        "Figure 9\n"
        "Misfit Levels (Max |ZSTD|)"
    )

    return "\n\n".join(blocks)


def mfrm_diagnostics(res, interaction_pairs=None, top_n_interactions=20, whexact=False):
    obs_df = compute_obs_table(res)
    facet_cols = ["Person"] + res["config"]["facet_names"]
    overall_fit = calc_overall_fit(obs_df, whexact=whexact)
    fit_tbl = calc_facet_fit(obs_df, facet_cols, whexact=whexact)
    se_tbl = calc_facet_se(obs_df, facet_cols)
    bias_tbl = calc_bias_facet(obs_df, facet_cols)
    interaction_tbl = calc_bias_interactions(obs_df, facet_cols, pairs=interaction_pairs, top_n=top_n_interactions)
    ptmea_tbl = calc_ptmea(obs_df, facet_cols)

    person_tbl = res["facets"]["person"].copy()
    person_tbl["Facet"] = "Person"
    person_tbl["Level"] = person_tbl["Person"].astype(str)
    person_tbl["SE"] = person_tbl["SD"] if "SD" in person_tbl.columns else np.nan

    facet_tbl = res["facets"]["others"].copy()
    facet_tbl["Level"] = facet_tbl["Level"].astype(str)
    facet_tbl["SE"] = np.nan

    measures = pd.concat([
        person_tbl[["Facet", "Level", "Estimate", "SE"]],
        facet_tbl[["Facet", "Level", "Estimate", "SE"]],
    ], ignore_index=True)

    measures = measures.merge(se_tbl, on=["Facet", "Level"], how="left", suffixes=("", "_calc"))
    measures["SE"] = measures["SE"].fillna(measures["SE_calc"])
    measures = measures.drop(columns=["SE_calc"])
    if "N" in measures.columns:
        measures = measures.rename(columns={"N": "N_SE"})
    fit_tbl = fit_tbl.rename(columns={"N": "N_Fit"}) if "N" in fit_tbl.columns else fit_tbl
    bias_tbl = bias_tbl.rename(columns={"N": "N_Bias"}) if "N" in bias_tbl.columns else bias_tbl
    ptmea_tbl = ptmea_tbl.rename(columns={"N": "N_PTMEA"}) if "N" in ptmea_tbl.columns else ptmea_tbl
    measures = measures.merge(fit_tbl, on=["Facet", "Level"], how="left")
    measures = measures.merge(bias_tbl, on=["Facet", "Level"], how="left")
    measures = measures.merge(ptmea_tbl, on=["Facet", "Level"], how="left")
    measures["CI_Lower"] = measures["Estimate"] - 1.96 * measures["SE"]
    measures["CI_Upper"] = measures["Estimate"] + 1.96 * measures["SE"]

    reliability_tbl = calc_reliability(measures)

    return {
        "obs": obs_df,
        "overall_fit": overall_fit,
        "measures": measures,
        "fit": fit_tbl,
        "reliability": reliability_tbl,
        "bias": bias_tbl,
        "interactions": interaction_tbl,
    }


def read_flexible_table(text_value, file_input, header=True):
    if file_input is not None:
        name = file_input.name.lower()
        sep = "\t" if name.endswith((".tsv", ".txt")) else ","
        return pd.read_csv(file_input, sep=sep, header=0 if header else None, dtype=str)
    if text_value is None or not str(text_value).strip():
        return pd.DataFrame()
    text_value = str(text_value).strip()
    if "\t" in text_value:
        sep = "\t"
    elif ";" in text_value:
        sep = ";"
    else:
        sep = ","
    return pd.read_csv(io.StringIO(text_value), sep=sep, header=0 if header else None, dtype=str)


def facet_summary_table(df, facet_cols):
    rows = []
    for col in ["Person"] + list(facet_cols):
        values = df[col].dropna().astype(str)
        examples = ", ".join(values.unique()[:5])
        rows.append({
            "Facet": col,
            "Levels": values.nunique(),
            "Missing": df[col].isna().sum(),
            "Examples": examples,
        })
    return pd.DataFrame(rows)


def response_distribution(df):
    if "Weight" in df.columns:
        dist = df.groupby("Score")["Weight"].sum().sort_index()
    else:
        dist = df["Score"].value_counts().sort_index()
    return pd.DataFrame({"Score": dist.index.astype(int), "Count": dist.values})


# ---- Streamlit UI ----
st.set_page_config(page_title="MFRM Estimation (Streamlit)", layout="wide")

st.title("MFRM Estimation (Streamlit)")

sample_df = sample_mfrm_data()

template_demo = format_tab_template(sample_df.head(24))
template_toy = format_tab_template(sample_df.head(8))

def ensure_session_state():
    if "pasted_data" not in st.session_state:
        st.session_state.pasted_data = template_demo

ensure_session_state()

st.sidebar.header("Data source")
data_source = st.sidebar.radio(
    "Source",
    ["Sample data (built-in)", "Paste table", "Upload file"],
    index=0,
)

anchor_text = None
anchor_file = None
group_anchor_text = None
group_anchor_file = None

if data_source == "Upload file":
    file = st.sidebar.file_uploader("Upload CSV/TSV", type=["csv", "tsv", "txt"])
    header = st.sidebar.checkbox("Header", value=True)
    sep = st.sidebar.selectbox("Separator", [",", "\t", ";"], index=0)
    df_raw = pd.read_csv(file, sep=sep, header=0 if header else None) if file else pd.DataFrame()
    if not header and not df_raw.empty:
        df_raw.columns = [f"Column_{i+1}" for i in range(df_raw.shape[1])]
elif data_source == "Paste table":
    st.sidebar.caption("Paste data into the template. Tab-delimited is recommended.")
    col1, col2, col3 = st.sidebar.columns(3)
    if col1.button("Clear"):
        st.session_state.pasted_data = ""
    if col2.button("Demo template"):
        st.session_state.pasted_data = template_demo
    if col3.button("Toy template"):
        st.session_state.pasted_data = template_toy
    pasted = st.sidebar.text_area("Paste table data", key="pasted_data", height=200)
    paste_header = st.sidebar.checkbox("Header row", value=True)
    paste_sep = st.sidebar.selectbox("Separator", ["\t", ",", ";"], index=0)
    df_raw = read_flexible_table(pasted, None, header=paste_header)
    if not paste_header and not df_raw.empty:
        df_raw.columns = [f"Column_{i+1}" for i in range(df_raw.shape[1])]
else:
    st.sidebar.info("Sample data is synthetic and included for quick demonstration.")
    df_raw = sample_df.copy()

if df_raw is None:
    df_raw = pd.DataFrame()

if not df_raw.empty:
    st.sidebar.subheader("Column selection")
    cols = list(df_raw.columns)
    person_guess = guess_col(cols, ["person", "examinee", "student", "id"], fallback=0)
    score_guess = guess_col(cols, ["score", "rating", "result"], fallback=min(1, len(cols) - 1))

    person_col = st.sidebar.selectbox(
        "Person column",
        cols,
        index=cols.index(person_guess) if person_guess in cols else 0,
        help="Identifier for the unit being measured (e.g., student, examinee).",
    )
    score_col = st.sidebar.selectbox(
        "Score column",
        cols,
        index=cols.index(score_guess) if score_guess in cols else min(1, len(cols) - 1),
        help="Observed rating category for each row.",
    )
    weight_choice = st.sidebar.selectbox(
        "Weight column (optional)",
        ["(None)"] + cols,
        index=0,
        help="Optional frequency/weight per observation. Must be positive.",
    )
    weight_col = None if weight_choice == "(None)" else weight_choice

    facet_choices = [c for c in cols if c != weight_col]
    default_facets = [c for c in facet_choices if c not in (person_col, score_col)]
    facet_cols = st.sidebar.multiselect(
        "Facet columns (2+)",
        facet_choices,
        default=default_facets[:3] if len(default_facets) >= 2 else default_facets,
        help="Facets beyond Person and Score (e.g., Rater, Task, Criterion).",
    )
    if person_col in facet_cols or score_col in facet_cols:
        st.sidebar.warning("Facet columns should not include Person or Score. Please deselect them.")

    st.sidebar.subheader("Model")
    model = st.sidebar.radio(
        "Model",
        ["RSM", "PCM"],
        index=0,
        help="RSM = shared step structure; PCM = step structure varies by a selected facet.",
    )
    method = st.sidebar.radio(
        "Estimation",
        ["JMLE", "MML"],
        index=0,
        help="JMLE is the traditional many-facet approach; MML provides EAP/SD (experimental).",
    )
    keep_original = st.sidebar.checkbox(
        "Keep original category values (K)",
        value=False,
        help="Use if your rating categories are non-contiguous (e.g., 1,3,5).",
    )

    step_facet = None
    if model == "PCM":
        if facet_cols:
            step_facet = st.sidebar.selectbox(
                "Step facet",
                facet_cols,
                index=0,
                help="Facet that defines distinct step structures under PCM.",
            )
        else:
            step_facet = None

    quad_points = st.sidebar.number_input(
        "Quadrature points (MML)",
        min_value=5,
        max_value=61,
        value=15,
        step=2,
        help="Higher values improve accuracy but increase computation.",
    ) if method == "MML" else 15
    maxit = st.sidebar.number_input(
        "Max iterations",
        min_value=50,
        max_value=2000,
        value=400,
        step=50,
        help="Maximum optimization iterations.",
    )

    st.sidebar.subheader("Constraints (optional)")
    noncenter_facet = st.sidebar.selectbox(
        "Non-centered facet",
        ["Person"] + list(facet_cols),
        index=0,
        help="Select exactly one facet to float (not centered).",
    )
    dummy_facets = st.sidebar.multiselect(
        "Dummy facets",
        ["Person"] + list(facet_cols),
        default=[],
        help="All elements anchored at 0 (use for adjustment facets).",
    )
    positive_facets = st.sidebar.multiselect(
        "Positive facets (higher values increase score)",
        list(facet_cols),
        default=[],
        help="By default facets are negatively oriented (higher = more severe).",
    )

    with st.sidebar.expander("Anchors (CSV/TSV)"):
        anchor_file = st.file_uploader("Anchor table file", type=["csv", "tsv", "txt"], key="anchor_file")
        anchor_text = st.text_area(
            "Anchor table text",
            height=120,
            key="anchor_text",
            help="Columns: Facet, Level, Anchor (logit).",
        )

    with st.sidebar.expander("Group anchors (CSV/TSV)"):
        group_anchor_file = st.file_uploader("Group anchor table file", type=["csv", "tsv", "txt"], key="group_anchor_file")
        group_anchor_text = st.text_area(
            "Group anchor table text",
            height=120,
            key="group_anchor_text",
            help="Columns: Facet, Level, Group, GroupValue.",
        )

    st.sidebar.subheader("Reporting (Facets-style)")
    report_totalscore_label = st.sidebar.radio(
        "Total scores reported",
        ["Yes (Total Score/Count)", "No (Obsvd Score/Count)"],
        index=0,
        help="Whether extreme-score observations are included in totals.",
    )
    report_totalscore = report_totalscore_label.startswith("Yes")
    report_omit_unobserved_label = st.sidebar.radio(
        "Omit unobserved elements",
        ["No", "Yes"],
        index=0,
        help="Hide elements with zero observations in the report.",
    )
    report_omit_unobserved = report_omit_unobserved_label == "Yes"
    report_xtreme = st.sidebar.number_input(
        "Xtreme correction (fraction)",
        value=0.0,
        min_value=0.0,
        step=0.1,
        help="Finite correction for minimum/maximum scores.",
    )
    report_umean = st.sidebar.number_input(
        "Umean (origin)",
        value=0.0,
        step=0.1,
        help="User-scale origin (reported measure = logit * uscale + umean).",
    )
    report_uscale = st.sidebar.number_input(
        "Uscale (units per logit)",
        value=1.0,
        step=0.1,
        help="User-scale slope (units per logit).",
    )
    report_udecimals = st.sidebar.number_input(
        "Udecimals (report precision)",
        value=2,
        min_value=0,
        step=1,
        help="Decimal places in report tables.",
    )

    st.sidebar.subheader("Fit statistics")
    whexact = st.sidebar.checkbox(
        "WHEXACT (Exact ZSTD)",
        value=False,
        help="Use exact standardization for ZSTD (slower).",
    )

    st.sidebar.subheader("Bias/Interaction (Facets-style)")
    bias_run = st.sidebar.checkbox(
        "Estimate bias/interaction",
        value=False,
        help="Re-estimate bias terms with main effects fixed (Facets-style).",
    )
    bias_pair = "(None)"
    bias_max_abs = 10.0
    bias_omit_extreme = True
    if bias_run:
        facets_all = ["Person"] + list(facet_cols)
        pair_options = [f"{a} x {b}" for a, b in combinations(facets_all, 2)]
        bias_pair = st.sidebar.selectbox(
            "Bias pair",
            ["(None)"] + pair_options,
            index=0,
            help="Choose two facets to estimate interaction/bias terms.",
        )
        bias_max_abs = st.sidebar.number_input(
            "Max |bias| (logit)",
            value=10.0,
            step=1.0,
            help="Bounds for bias estimation (logits).",
        )
        bias_omit_extreme = st.sidebar.checkbox(
            "Omit extreme elements",
            value=True,
            help="Exclude elements with only minimum or maximum scores.",
        )

    st.sidebar.subheader("Run")
    run_button = st.sidebar.button("Run estimation", type="primary")
else:
    st.warning("No data loaded. Please paste or upload data.")
    run_button = False

if run_button:
    if len(facet_cols) < 2:
        st.error("Select at least two facet columns (Person + 2+ facets).")
    else:
        with st.spinner("Running estimation..."):
            try:
                anchor_df = read_flexible_table(anchor_text, anchor_file)
                group_anchor_df = read_flexible_table(group_anchor_text, group_anchor_file)
                res = mfrm_estimate(
                    df_raw,
                    person_col=person_col,
                    facet_cols=facet_cols,
                    score_col=score_col,
                    weight_col=weight_col,
                    keep_original=keep_original,
                    model=model,
                    method=method,
                    step_facet=step_facet,
                    anchor_df=anchor_df,
                    group_anchor_df=group_anchor_df,
                    noncenter_facet=noncenter_facet,
                    dummy_facets=dummy_facets,
                    positive_facets=positive_facets,
                    quad_points=int(quad_points),
                    maxit=int(maxit),
                )
                diagnostics = mfrm_diagnostics(res, whexact=whexact)
                st.session_state["results"] = res
                st.session_state["diagnostics"] = diagnostics
                st.session_state["whexact"] = whexact
                bias_results = None
                if bias_run and bias_pair != "(None)":
                    try:
                        facet_a, facet_b = bias_pair.split(" x ", 1)
                        bias_results = estimate_bias_interaction(
                            res,
                            diagnostics,
                            facet_a,
                            facet_b,
                            max_abs=float(bias_max_abs),
                            omit_extreme=bool(bias_omit_extreme),
                        )
                    except Exception:
                        bias_results = None
                st.session_state["bias_results"] = bias_results
            except Exception as exc:
                st.error(f"Estimation failed: {exc}")

res = st.session_state.get("results")
diagnostics = st.session_state.get("diagnostics")

if res:
    tabs = st.tabs([
        "Data",
        "Summary",
        "Results",
        "Person",
        "Facets",
        "Report",
        "Fit",
        "Reliability",
        "Bias",
        "Interactions",
        "Steps",
        "Categories",
        "Dimensionality",
        "Subsets",
        "Visuals",
        "Downloads",
        "Help",
    ])

    with tabs[0]:
        st.caption("Preview your data and verify the facets and score distribution.")
        st.subheader("Data preview")
        st.dataframe(res["prep"]["data"].head(200))
        st.subheader("Facet summary")
        st.dataframe(facet_summary_table(res["prep"]["data"], res["config"]["facet_names"]))
        st.subheader("Response distribution")
        dist = response_distribution(res["prep"]["data"])
        st.plotly_chart(px.bar(dist, x="Score", y="Count", title="Observed category counts"), use_container_width=True)
        render_context_notes("data")

    with tabs[1]:
        st.caption("Model/estimation summary and anchor diagnostics.")
        st.subheader("Estimation summary")
        st.dataframe(res["summary"])
        st.subheader("Anchor summary")
        st.dataframe(res["config"]["anchor_summary"])

    with tabs[2]:
        st.caption("High-level checks: reliability and overall model summary.")
        st.subheader("Quick overview")
        st.dataframe(res["summary"])
        st.subheader("Reliability snapshot")
        st.dataframe(diagnostics["reliability"])
        if not diagnostics["reliability"].empty:
            fig = px.bar(diagnostics["reliability"], x="Facet", y="Reliability", title="Reliability by facet")
            st.plotly_chart(fig, use_container_width=True)

    with tabs[3]:
        st.caption("Person (ratee) measures and distribution.")
        st.subheader("Person measures")
        st.dataframe(res["facets"]["person"])
        if not res["facets"]["person"].empty:
            fig = px.histogram(res["facets"]["person"], x="Estimate", nbins=30, title="Person measure distribution")
            st.plotly_chart(fig, use_container_width=True)
        render_context_notes("person")

    with tabs[4]:
        st.caption("Facet measures for raters/tasks/criteria, etc.")
        st.subheader("Facet measures")
        st.dataframe(res["facets"]["others"])
        render_context_notes("facets")

    with tabs[5]:
        st.caption("Facets-style measurement report (app-friendly).")
        st.subheader("Measurement report (Facets-style)")
        report_tbls = calc_facets_report_tbls(
            res,
            diagnostics,
            totalscore=report_totalscore,
            umean=report_umean,
            uscale=report_uscale,
            udecimals=report_udecimals,
            omit_unobserved=report_omit_unobserved,
            xtreme=report_xtreme,
        )
        if not report_tbls:
            st.info("No report available.")
        else:
            facets = list(report_tbls.keys())
            selected_facet = st.selectbox("Facet", facets, index=0, key="facets_report_select")
            tbl = report_tbls[selected_facet].copy()
            decimals = int(report_udecimals) if report_udecimals is not None else 2
            weight_is_int = True
            if "Weight" in res["prep"]["data"].columns:
                w = res["prep"]["data"]["Weight"].to_numpy(dtype=float)
                weight_is_int = np.allclose(w, np.round(w))
            for col in tbl.select_dtypes(include=[np.number]).columns:
                if col in ("TotalScore", "TotalCount"):
                    tbl[col] = tbl[col].round(0).astype("Int64")
                elif col in ("WeightdScore", "WeightdCount") and weight_is_int:
                    tbl[col] = tbl[col].round(0).astype("Int64")
                else:
                    tbl[col] = tbl[col].round(decimals)
            st.dataframe(tbl)
            st.download_button(
                f"Download {selected_facet} report (CSV)",
                tbl.to_csv(index=False).encode("utf-8"),
                file_name=f"mfrm_facets_report_{selected_facet}.csv",
            )
            combined = []
            for facet, df in report_tbls.items():
                tmp = df.copy()
                tmp.insert(0, "Facet", facet)
                combined.append(tmp)
            combined_df = pd.concat(combined, ignore_index=True)
            st.download_button(
                "Download all facets report (CSV)",
                combined_df.to_csv(index=False).encode("utf-8"),
                file_name="mfrm_facets_report_all.csv",
            )
            note_map = build_apa_table_figure_note_map(
                res,
                diagnostics,
                bias_results=st.session_state.get("bias_results"),
                context={},
                whexact=st.session_state.get("whexact", False),
            )
            render_note_text(note_map.get("table1", ""))
            render_context_notes("report")

    with tabs[6]:
        st.caption("Overall and facet-level fit statistics.")
        st.subheader("Overall fit")
        st.dataframe(diagnostics["overall_fit"])
        st.subheader("Facet fit")
        st.dataframe(diagnostics["fit"])
        if not diagnostics["fit"].empty:
            fig = px.scatter(
                diagnostics["fit"],
                x="Infit",
                y="Outfit",
                color="Facet",
                hover_data=["Level", "InfitZSTD", "OutfitZSTD"],
                title="Infit vs Outfit",
            )
            st.plotly_chart(fig, use_container_width=True)
        render_context_notes("fit")

    with tabs[7]:
        st.caption("Reliability, chi-square tests, and inter-rater agreement.")
        st.subheader("Reliability")
        st.dataframe(diagnostics["reliability"])
        st.subheader("Reliability chi-square")
        st.dataframe(calc_facets_chisq(diagnostics["measures"]))
        st.subheader("Inter-rater agreement")
        facet_names = res["config"]["facet_names"]
        rater_note_facet = ""
        if not facet_names:
            st.info("Select facet columns to enable agreement checks.")
        else:
            agreement_facet = st.selectbox(
                "Rater facet",
                facet_names,
                index=0,
                key="agreement_facet",
            )
            rater_note_facet = agreement_facet
            agreement = calc_interrater_agreement(
                diagnostics["obs"],
                ["Person"] + facet_names,
                agreement_facet,
                res=res,
            )
            if agreement["summary"].empty:
                st.info("Agreement summary not available.")
            else:
                st.dataframe(agreement["summary"])
                st.dataframe(agreement["pairs"])
        note_map = build_apa_table_figure_note_map(
            res,
            diagnostics,
            bias_results=st.session_state.get("bias_results"),
            context={"rater_facet": rater_note_facet},
            whexact=st.session_state.get("whexact", False),
        )
        render_note_text(note_map.get("table3", ""))
        render_context_notes("reliability")

    with tabs[8]:
        st.caption("Bias/interaction analyses with fixed main effects.")
        st.subheader("Bias/Interaction (Facets-style)")
        bias_results = st.session_state.get("bias_results")
        if bias_results and bias_results.get("table") is not None and not bias_results["table"].empty:
            facet_a = bias_results["facet_a"]
            facet_b = bias_results["facet_b"]
            tbl = bias_results["table"].copy()
            tbl_display = tbl[[
                "Sq",
                "Observd Score",
                "Expctd Score",
                "Observd Count",
                "Obs-Exp Average",
                "Bias Size",
                "S.E.",
                "t",
                "d.f.",
                "Prob.",
                "Infit",
                "Outfit",
                "FacetA_Index",
                "FacetA_Level",
                "FacetA_Measure",
                "FacetB_Index",
                "FacetB_Level",
                "FacetB_Measure",
            ]].rename(columns={
                "S.E.": "Model S.E.",
                "Infit": "Infit MnSq",
                "Outfit": "Outfit MnSq",
                "FacetA_Index": f"{facet_a} N",
                "FacetA_Level": facet_a,
                "FacetA_Measure": f"{facet_a} measr",
                "FacetB_Index": f"{facet_b} N",
                "FacetB_Level": facet_b,
                "FacetB_Measure": f"{facet_b} measr",
            })
            st.dataframe(tbl_display)
            if "iteration" in bias_results and bias_results["iteration"] is not None:
                st.subheader("Bias iteration report (convergence)")
                iter_tbl = bias_results["iteration"].copy()
                iter_display = iter_tbl[[
                    "Iteration",
                    "MaxScoreResidual",
                    "MaxScoreResidualPct",
                    "MaxScoreResidualCategories",
                    "MaxLogitChange",
                    "BiasCells",
                ]].rename(columns={
                    "MaxScoreResidual": "Max. Score Residual",
                    "MaxScoreResidualPct": "%",
                    "MaxScoreResidualCategories": "Categories",
                    "MaxLogitChange": "Max. Logit Change",
                    "BiasCells": "Cells",
                })
                iter_display["Iteration"] = iter_display["Iteration"].apply(lambda v: f"BIAS {int(v)}" if pd.notna(v) else v)
                st.dataframe(iter_display)

            st.subheader("Bias summary (facet-pair overview)")
            st.dataframe(bias_results["summary"])
            size_hist, size_stats = build_histogram(tbl["Bias Size"].to_numpy(), bins=20)
            if not size_hist.empty:
                fig = px.bar(size_hist, x="Value", y="Count", title="Bias size distribution")
                for mult in [0, 1, 2, 3]:
                    if mult == 0:
                        pos = size_stats["mean"]
                        fig.add_vline(x=pos, line_dash="solid", line_color="black")
                    else:
                        fig.add_vline(x=size_stats["mean"] + mult * size_stats["sd"], line_dash="dash", line_color="gray")
                        fig.add_vline(x=size_stats["mean"] - mult * size_stats["sd"], line_dash="dash", line_color="gray")
                st.plotly_chart(fig, use_container_width=True)

            t_hist, t_stats = build_histogram(tbl["t"].to_numpy(), bins=20)
            if not t_hist.empty:
                fig = px.bar(t_hist, x="Value", y="Count", title="Bias significance (t) distribution")
                for mult in [0, 1, 2, 3]:
                    if mult == 0:
                        pos = t_stats["mean"]
                        fig.add_vline(x=pos, line_dash="solid", line_color="black")
                    else:
                        fig.add_vline(x=t_stats["mean"] + mult * t_stats["sd"], line_dash="dash", line_color="gray")
                        fig.add_vline(x=t_stats["mean"] - mult * t_stats["sd"], line_dash="dash", line_color="gray")
                st.plotly_chart(fig, use_container_width=True)

            st.subheader("Fixed (all = 0) chi-square")
            st.dataframe(bias_results["chi_sq"])

            st.subheader("Pairwise bias report (target vs context)")
            target_facet = st.selectbox("Target facet", [facet_a, facet_b], index=0, key="bias_target_facet")
            context_facet = facet_b if target_facet == facet_a else facet_a
            pairwise_tbl = calc_bias_pairwise(bias_results["table"], target_facet, context_facet)
            if pairwise_tbl.empty:
                st.info("Pairwise bias report not available for the current selection.")
            else:
                pair_display = pairwise_tbl[[
                    "Target N",
                    "Target",
                    "Target Measure",
                    "Target S.E.",
                    "Context1 N",
                    "Context1",
                    "Local Measure1",
                    "SE1",
                    "Obs-Exp Avg1",
                    "Count1",
                    "Context2 N",
                    "Context2",
                    "Local Measure2",
                    "SE2",
                    "Obs-Exp Avg2",
                    "Count2",
                    "Contrast",
                    "SE",
                    "t",
                    "d.f.",
                    "Prob.",
                ]].rename(columns={
                    "Target Measure": "Target Measr",
                    "Target S.E.": "Target S.E.",
                    "Local Measure1": "Context1 Measr",
                    "Local Measure2": "Context2 Measr",
                    "SE1": "Context1 S.E.",
                    "SE2": "Context2 S.E.",
                    "Obs-Exp Avg1": "Obs-Exp Avg1",
                    "Obs-Exp Avg2": "Obs-Exp Avg2",
                })
                st.dataframe(pair_display)

            with st.expander("Fixed-width output (Facets-style text)"):
                def fmt_int(v):
                    return "" if pd.isna(v) else f"{int(round(v))}"

                def fmt_2(v):
                    return "" if pd.isna(v) else f"{float(v):.2f}"

                def fmt_3(v):
                    return "" if pd.isna(v) else f"{float(v):.3f}"

                bias_cols = [
                    "Sq",
                    "Observd Score",
                    "Expctd Score",
                    "Observd Count",
                    "Obs-Exp Average",
                    "Bias Size",
                    "Model S.E.",
                    "t",
                    "d.f.",
                    "Prob.",
                    "Infit MnSq",
                    "Outfit MnSq",
                    f"{facet_a} N",
                    facet_a,
                    f"{facet_a} measr",
                    f"{facet_b} N",
                    facet_b,
                    f"{facet_b} measr",
                ]
                bias_formats = {
                    "Sq": "{}",
                    "Observd Score": "{:.2f}",
                    "Expctd Score": "{:.2f}",
                    "Observd Count": "{:.0f}",
                    "Obs-Exp Average": "{:.2f}",
                    "Bias Size": "{:.2f}",
                    "Model S.E.": "{:.2f}",
                    "t": "{:.2f}",
                    "d.f.": "{:.0f}",
                    "Prob.": "{:.4f}",
                    "Infit MnSq": "{:.2f}",
                    "Outfit MnSq": "{:.2f}",
                    f"{facet_a} N": "{:.0f}",
                    f"{facet_a} measr": "{:.2f}",
                    f"{facet_b} N": "{:.0f}",
                    f"{facet_b} measr": "{:.2f}",
                }
                bias_fixed = format_fixed_width_table(
                    tbl_display,
                    bias_cols,
                    formats=bias_formats,
                    max_col_width=18,
                )
                bias_fixed_text = build_bias_fixed_text(
                    tbl_display,
                    bias_results["summary"],
                    bias_results["chi_sq"],
                    facet_a,
                    facet_b,
                    bias_cols,
                    bias_formats,
                )
                st.code(bias_fixed_text, language="text")

                if not pairwise_tbl.empty:
                    pair_cols = [
                        "Target N",
                        "Target",
                        "Target Measr",
                        "Target S.E.",
                        "Context1 N",
                        "Context1",
                        "Context1 Measr",
                        "Context1 S.E.",
                        "Obs-Exp Avg1",
                        "Count1",
                        "Context2 N",
                        "Context2",
                        "Context2 Measr",
                        "Context2 S.E.",
                        "Obs-Exp Avg2",
                        "Count2",
                        "Contrast",
                        "SE",
                        "t",
                        "d.f.",
                        "Prob.",
                    ]
                    pair_formats = {
                        "Target N": "{:.0f}",
                        "Target Measr": "{:.2f}",
                        "Target S.E.": "{:.2f}",
                        "Context1 N": "{:.0f}",
                        "Context1 Measr": "{:.2f}",
                        "Context1 S.E.": "{:.2f}",
                        "Obs-Exp Avg1": "{:.2f}",
                        "Count1": "{:.0f}",
                        "Context2 N": "{:.0f}",
                        "Context2 Measr": "{:.2f}",
                        "Context2 S.E.": "{:.2f}",
                        "Obs-Exp Avg2": "{:.2f}",
                        "Count2": "{:.0f}",
                        "Contrast": "{:.2f}",
                        "SE": "{:.2f}",
                        "t": "{:.2f}",
                        "d.f.": "{:.0f}",
                        "Prob.": "{:.4f}",
                    }
                    pair_fixed_text = build_pairwise_fixed_text(
                        pair_display,
                        target_facet,
                        context_facet,
                        pair_cols,
                        pair_formats,
                    )
                    st.code(pair_fixed_text, language="text")
            note_map = build_apa_table_figure_note_map(
                res,
                diagnostics,
                bias_results=bias_results,
                context={},
                whexact=st.session_state.get("whexact", False),
            )
            render_note_text(note_map.get("table4", ""))
            render_context_notes("bias")
        else:
            st.info("Run bias estimation in the sidebar to see Facets-style bias tables.")
        with st.expander("Legacy residual bias (main model)"):
            st.dataframe(diagnostics["bias"])

    with tabs[9]:
        st.caption("Largest interaction residuals (legacy diagnostics).")
        st.subheader("Top interactions")
        st.dataframe(diagnostics["interactions"])

    with tabs[10]:
        st.caption("Step/threshold estimates and ordering.")
        st.subheader("Steps / thresholds")
        st.dataframe(res["steps"])
        step_order = calc_step_order(res["steps"])
        if not step_order.empty:
            st.subheader("Step ordering")
            st.dataframe(step_order)
        render_context_notes("steps")

    with tabs[11]:
        st.caption("Category counts, fit, and threshold warnings.")
        st.subheader("Category diagnostics")
        cat_tbl = calc_category_stats(
            diagnostics["obs"],
            res,
            whexact=st.session_state.get("whexact", False),
        )
        st.dataframe(cat_tbl)
        st.subheader("Category warnings")
        st.text(category_warnings_text(cat_tbl, calc_step_order(res["steps"])))
        if not cat_tbl.empty:
            st.plotly_chart(px.bar(cat_tbl, x="Category", y="Count", title="Category counts"), use_container_width=True)
            st.plotly_chart(px.bar(cat_tbl, x="Category", y="Percent", title="Category percent"), use_container_width=True)
        note_map = build_apa_table_figure_note_map(
            res,
            diagnostics,
            bias_results=st.session_state.get("bias_results"),
            context={},
            whexact=st.session_state.get("whexact", False),
        )
        render_note_text(note_map.get("table2", ""))
        render_context_notes("categories")

    with tabs[12]:
        st.caption("Residual PCA for dimensionality checks.")
        st.subheader("Residual PCA (Dimensionality)")
        facet_choices = ["Overall"] + res["config"]["facet_names"]
        pca_facet = st.selectbox("Facet for residual PCA", facet_choices, index=0, key="pca_facet")
        if pca_facet == "Overall":
            pca_bundle = compute_pca_overall(diagnostics["obs"], res["config"]["facet_names"])
        else:
            pca_bundle = compute_pca_by_facet(diagnostics["obs"], [pca_facet]).get(pca_facet)

        if pca_bundle is None:
            st.info("PCA results not available. Check data variability and sample size.")
        else:
            eigenvalues = pca_bundle["eigenvalues"]
            var_pct = pca_bundle["variance_pct"]
            loadings = pca_bundle["loadings"]

            max_components = min(20, len(eigenvalues))
            scree_df = pd.DataFrame({
                "Component": np.arange(1, max_components + 1),
                "Eigenvalue": eigenvalues[:max_components],
            })
            fig = px.line(
                scree_df,
                x="Component",
                y="Eigenvalue",
                markers=True,
                title=f"Scree plot ({pca_facet})",
            )
            fig.add_hline(y=1, line_dash="dash", line_color="gray")
            st.plotly_chart(fig, use_container_width=True)

            n_components = min(10, len(eigenvalues))
            eigen_df = pd.DataFrame({
                "Component": [f"PC{i}" for i in range(1, n_components + 1)],
                "Eigenvalue": eigenvalues[:n_components],
                "Variance_Pct": var_pct[:n_components],
                "Cumulative_Pct": np.cumsum(var_pct[:n_components]),
            })
            st.subheader("Eigenvalues and variance explained")
            st.dataframe(eigen_df.round(3))

            if "PC1" in loadings.columns:
                top_loadings = (
                    loadings["PC1"]
                    .abs()
                    .sort_values(ascending=False)
                    .head(20)
                    .index
                )
                loadings_plot = loadings.loc[top_loadings, "PC1"].sort_values()
                fig = px.bar(
                    loadings_plot,
                    x=loadings_plot.values,
                    y=loadings_plot.index,
                    orientation="h",
                    title="PC1 loadings (top 20)",
                )
                st.plotly_chart(fig, use_container_width=True)

            if "PC1" in loadings.columns and "PC2" in loadings.columns:
                loadings_df = loadings[["PC1", "PC2"]].copy()
                loadings_df["Distance"] = np.sqrt(loadings_df["PC1"] ** 2 + loadings_df["PC2"] ** 2)
                top_biplot = loadings_df.sort_values("Distance", ascending=False).head(30)
                fig = px.scatter(
                    top_biplot,
                    x="PC1",
                    y="PC2",
                    text=top_biplot.index,
                    title="Loadings biplot (top 30)",
                )
                fig.update_traces(textposition="top center")
                fig.add_hline(y=0, line_dash="dash", line_color="gray")
                fig.add_vline(x=0, line_dash="dash", line_color="gray")
                st.plotly_chart(fig, use_container_width=True)

            pc1_var = var_pct[0] if len(var_pct) > 0 else np.nan
            eig_ratio = eigenvalues[0] / eigenvalues[1] if len(eigenvalues) > 1 else np.nan
            eig_above_1 = int(np.sum(eigenvalues > 1))
            assessment = (
                "Potential multidimensionality"
                if np.isfinite(pc1_var) and (pc1_var > 20 or (np.isfinite(eig_ratio) and eig_ratio < 3))
                else "Acceptable unidimensionality"
            )
            dim_df = pd.DataFrame({
                "Metric": [
                    "PC1 variance",
                    "Eigenvalue ratio (PC1/PC2)",
                    "Eigenvalues > 1",
                    "Assessment",
                ],
                "Value": [
                    f"{pc1_var:.1f}%" if np.isfinite(pc1_var) else "N/A",
                    f"{eig_ratio:.2f}" if np.isfinite(eig_ratio) else "N/A",
                    eig_above_1,
                    assessment,
                ],
            })
            st.subheader("Unidimensionality assessment")
            st.dataframe(dim_df)
        render_context_notes("dimensionality")

    with tabs[13]:
        st.caption("Connectivity diagnostics for disjoint subsets.")
        st.subheader("Connectivity subsets")
        subsets = calc_subsets(diagnostics["obs"], ["Person"] + res["config"]["facet_names"])
        if subsets["summary"].empty:
            st.info("No subset diagnostics available.")
        else:
            st.dataframe(subsets["summary"])
            st.subheader("Subset nodes")
            st.dataframe(subsets["nodes"])
        render_context_notes("subsets")

    with tabs[14]:
        st.caption("Visual summaries for measures, fit, and categories.")
        st.markdown(
            "**Workflow** 1) Read the APA title/note. 2) Check warnings and summary cues. 3) Interpret the chart."
        )
        with st.expander("Notes and warning options", expanded=False):
            preset_defs = {
                "Balanced (default)": {
                    "n_obs_min": 100,
                    "n_person_min": 30,
                    "low_cat_min": 10,
                    "min_facet_levels": 3,
                    "misfit_ratio_warn": 0.10,
                    "missing_fit_ratio_warn": 0.20,
                    "zstd2_ratio_warn": 0.10,
                    "zstd3_ratio_warn": 0.05,
                    "expected_var_min": 0.2,
                },
                "Classroom / low stakes": {
                    "n_obs_min": 50,
                    "n_person_min": 20,
                    "low_cat_min": 5,
                    "min_facet_levels": 2,
                    "misfit_ratio_warn": 0.20,
                    "missing_fit_ratio_warn": 0.30,
                    "zstd2_ratio_warn": 0.20,
                    "zstd3_ratio_warn": 0.10,
                    "expected_var_min": 0.1,
                },
                "High-stakes": {
                    "n_obs_min": 200,
                    "n_person_min": 50,
                    "low_cat_min": 20,
                    "min_facet_levels": 5,
                    "misfit_ratio_warn": 0.05,
                    "missing_fit_ratio_warn": 0.10,
                    "zstd2_ratio_warn": 0.05,
                    "zstd3_ratio_warn": 0.02,
                    "expected_var_min": 0.3,
                },
            }
            preset_name = st.selectbox(
                "Warning preset",
                list(preset_defs.keys()),
                index=0,
                key="warn_preset",
            )
            if st.button("Apply preset", key="warn_apply_preset"):
                preset = preset_defs[preset_name]
                st.session_state["warn_n_obs"] = preset["n_obs_min"]
                st.session_state["warn_n_person"] = preset["n_person_min"]
                st.session_state["warn_low_cat"] = preset["low_cat_min"]
                st.session_state["warn_min_facet_levels"] = preset["min_facet_levels"]
                st.session_state["warn_misfit_ratio"] = preset["misfit_ratio_warn"]
                st.session_state["warn_missing_fit_ratio"] = preset["missing_fit_ratio_warn"]
                st.session_state["warn_zstd2_ratio"] = preset["zstd2_ratio_warn"]
                st.session_state["warn_zstd3_ratio"] = preset["zstd3_ratio_warn"]
                st.session_state["warn_expected_var"] = preset["expected_var_min"]
            col1, col2 = st.columns(2)
            with col1:
                warn_n_obs = st.number_input(
                    "Warn if observations (N) <",
                    min_value=0,
                    value=st.session_state.get("warn_n_obs", 100),
                    step=10,
                    key="warn_n_obs",
                )
                warn_n_person = st.number_input(
                    "Warn if persons <",
                    min_value=0,
                    value=st.session_state.get("warn_n_person", 30),
                    step=5,
                    key="warn_n_person",
                )
                warn_low_cat = st.number_input(
                    "Low category count threshold",
                    min_value=1,
                    value=st.session_state.get("warn_low_cat", 10),
                    step=1,
                    key="warn_low_cat",
                )
                warn_min_facet_levels = st.number_input(
                    "Min facet levels for stability",
                    min_value=1,
                    value=st.session_state.get("warn_min_facet_levels", 3),
                    step=1,
                    key="warn_min_facet_levels",
                )
            with col2:
                warn_misfit_ratio = st.slider(
                    "Warn if misfit proportion >",
                    min_value=0.0,
                    max_value=1.0,
                    value=st.session_state.get("warn_misfit_ratio", 0.10),
                    step=0.05,
                    key="warn_misfit_ratio",
                )
                warn_missing_fit_ratio = st.slider(
                    "Warn if missing fit stats >",
                    min_value=0.0,
                    max_value=1.0,
                    value=st.session_state.get("warn_missing_fit_ratio", 0.20),
                    step=0.05,
                    key="warn_missing_fit_ratio",
                )
                warn_zstd2_ratio = st.slider(
                    "Warn if |ZSTD|>=2 proportion >",
                    min_value=0.0,
                    max_value=1.0,
                    value=st.session_state.get("warn_zstd2_ratio", 0.10),
                    step=0.05,
                    key="warn_zstd2_ratio",
                )
                warn_zstd3_ratio = st.slider(
                    "Warn if |ZSTD|>=3 proportion >",
                    min_value=0.0,
                    max_value=1.0,
                    value=st.session_state.get("warn_zstd3_ratio", 0.05),
                    step=0.05,
                    key="warn_zstd3_ratio",
                )
                warn_expected_var = st.number_input(
                    "Warn if expected score variance <",
                    min_value=0.0,
                    value=st.session_state.get("warn_expected_var", 0.2),
                    step=0.1,
                    key="warn_expected_var",
                )
            st.markdown("---")
            summary_detail = st.selectbox(
                "Summary detail",
                ["Brief", "Standard", "Detailed"],
                index=1,
                key="summary_detail",
            )
            summary_max_facets = st.slider(
                "Max facet ranges to list",
                min_value=1,
                max_value=8,
                value=4,
                step=1,
                key="summary_max_facets",
            )
            summary_top_misfit = st.slider(
                "Top misfit levels to list",
                min_value=0,
                max_value=10,
                value=3,
                step=1,
                key="summary_top_misfit",
            )
            show_summaries = st.checkbox("Show auto summaries", value=True, key="show_summaries")

        warn_thresholds = {
            "n_obs_min": warn_n_obs,
            "n_person_min": warn_n_person,
            "low_cat_min": warn_low_cat,
            "min_facet_levels": warn_min_facet_levels,
            "misfit_ratio_warn": warn_misfit_ratio,
            "missing_fit_ratio_warn": warn_missing_fit_ratio,
            "zstd2_ratio_warn": warn_zstd2_ratio,
            "zstd3_ratio_warn": warn_zstd3_ratio,
            "expected_var_min": warn_expected_var,
        }
        note_map = build_apa_table_figure_note_map(
            res,
            diagnostics,
            bias_results=st.session_state.get("bias_results"),
            context={},
            whexact=st.session_state.get("whexact", False),
        )
        warning_map = build_visual_warning_map(
            res,
            diagnostics,
            whexact=st.session_state.get("whexact", False),
            thresholds=warn_thresholds,
        )
        summary_map = build_visual_summary_map(
            res,
            diagnostics,
            whexact=st.session_state.get("whexact", False),
            options={
                "detail": summary_detail,
                "max_facet_ranges": summary_max_facets,
                "top_misfit_n": summary_top_misfit,
            },
        )

        st.subheader("Warnings summary")
        warnings_rows = []
        for key, msgs in warning_map.items():
            if not msgs:
                continue
            note_text = note_map.get(key, "")
            fig_title = note_text.split("\n", 1)[0].strip() if note_text else key
            for msg in msgs:
                warnings_rows.append({"Figure": fig_title, "Warning": msg})
        if not warnings_rows:
            st.success("No warnings triggered with the current thresholds.")
        else:
            warn_df = pd.DataFrame(warnings_rows)
            st.warning(f"{len(warnings_rows)} warning(s) across {warn_df['Figure'].nunique()} figure(s).")
            st.dataframe(warn_df)

        def render_apa_figure_note(key):
            note_text = note_map.get(key, "")
            if note_text:
                parts = note_text.split("\n", 1)
                title = parts[0].strip()
                note = parts[1].strip() if len(parts) > 1 else ""
                st.markdown(f"**{title}**")
                if note:
                    st.markdown(note)
            warn_msgs = warning_map.get(key, [])
            if warn_msgs:
                st.warning("Warnings:\n" + "\n".join([f"- {msg}" for msg in warn_msgs]))
            if show_summaries:
                summary_lines = summary_map.get(key, [])
                if summary_lines:
                    st.info("Summary:\n" + "\n".join([f"- {line}" for line in summary_lines]))
        vtabs = st.tabs([
            "Wright Map",
            "Pathway Map",
            "Facet Distribution",
            "Steps/Thresholds",
            "Category Curves",
            "Observed vs Expected",
            "Fit Scatter",
            "Fit ZSTD",
            "Misfit Levels",
        ])

        with vtabs[0]:
            render_apa_figure_note("figure1")
            st.subheader("Wright map (side-by-side)")
            theta_hat = (
                res["params"]["theta"]
                if res["config"]["method"] == "JMLE"
                else res["facets"]["person"]["Estimate"].to_numpy()
            )
            person_df = pd.DataFrame({"Estimate": theta_hat})
            facet_rows = []
            for facet, values in res["params"]["facets"].items():
                facet_rows.append(pd.DataFrame({"Type": facet, "Estimate": values}))
            facet_df = pd.concat(facet_rows, ignore_index=True) if facet_rows else pd.DataFrame(columns=["Type", "Estimate"])

            fig = make_subplots(rows=1, cols=2, subplot_titles=("Persons", "Facets"))
            fig.add_trace(
                go.Histogram(
                    x=person_df["Estimate"],
                    histnorm="probability density",
                    marker_color="#6baed6",
                    opacity=0.6,
                    name="Persons",
                ),
                row=1,
                col=1,
            )

            if not facet_df.empty:
                facet_names = facet_df["Type"].unique().tolist()
                y_map = {name: idx for idx, name in enumerate(facet_names)}
                rng = np.random.default_rng(42)
                facet_df["y"] = facet_df["Type"].map(y_map) + rng.normal(scale=0.08, size=len(facet_df))
                for facet_name in facet_names:
                    sub = facet_df[facet_df["Type"] == facet_name]
                    fig.add_trace(
                        go.Scatter(
                            x=sub["Estimate"],
                            y=sub["y"],
                            mode="markers",
                            name=str(facet_name),
                            marker=dict(size=6, opacity=0.85),
                            showlegend=True,
                            hovertemplate="Facet: %{text}<br>Estimate: %{x:.3f}<extra></extra>",
                            text=[facet_name] * len(sub),
                        ),
                        row=1,
                        col=2,
                    )
                fig.update_yaxes(
                    tickvals=list(y_map.values()),
                    ticktext=list(y_map.keys()),
                    row=1,
                    col=2,
                    title_text="Facet",
                )

            fig.update_xaxes(title_text="Logit scale", row=1, col=1)
            fig.update_xaxes(title_text="Logit scale", row=1, col=2)
            fig.update_layout(height=360, legend_orientation="h", legend_y=-0.15)
            st.plotly_chart(fig, use_container_width=True)

        with vtabs[1]:
            render_apa_figure_note("figure2")
            st.subheader("Pathway map (measure vs. infit t)")
            plot_df = diagnostics["measures"].copy()
            plot_df = plot_df.dropna(subset=["InfitZSTD", "Estimate"])
            if plot_df.empty:
                st.info("No pathway plot available.")
            else:
                plot_df["Label"] = plot_df["Level"].map(lambda x: truncate_label(x, width=16))
                plot_df["Tooltip"] = (
                    "Facet: " + plot_df["Facet"].astype(str)
                    + "<br>Level: " + plot_df["Level"].astype(str)
                    + "<br>Measure: " + plot_df["Estimate"].round(3).astype(str)
                    + "<br>Infit ZSTD: " + plot_df["InfitZSTD"].round(2).astype(str)
                )
                scatter_kwargs = dict(
                    data_frame=plot_df,
                    x="InfitZSTD",
                    y="Estimate",
                    color="Facet",
                    custom_data=["Tooltip"],
                )
                if "SE" in plot_df.columns:
                    scatter_kwargs["error_y"] = "SE"
                fig = px.scatter(**scatter_kwargs)
                fig.update_traces(hovertemplate="%{customdata[0]}<extra></extra>")
                fig.add_vline(x=-2, line_dash="dash", line_color="gray")
                fig.add_vline(x=2, line_dash="dash", line_color="gray")

                label_df = plot_df[plot_df["InfitZSTD"].abs() >= 2].copy()
                label_df = label_df.sort_values("InfitZSTD", key=lambda s: s.abs(), ascending=False).head(20)
                if not label_df.empty:
                    fig.add_trace(
                        go.Scatter(
                            x=label_df["InfitZSTD"],
                            y=label_df["Estimate"],
                            mode="text",
                            text=label_df["Label"],
                            textposition="top center",
                            showlegend=False,
                        )
                    )
                fig.update_layout(xaxis_title="Infit ZSTD", yaxis_title="Measure", height=360)
                st.plotly_chart(fig, use_container_width=True)

        with vtabs[2]:
            render_apa_figure_note("figure3")
            st.subheader("Facet estimate distribution")
            facet_rows = []
            for facet, values in res["params"]["facets"].items():
                facet_rows.append(pd.DataFrame({"Facet": facet, "Estimate": values}))
            facet_df = pd.concat(facet_rows, ignore_index=True) if facet_rows else pd.DataFrame()
            if facet_df.empty:
                st.info("No facet estimates available.")
            else:
                fig = px.box(
                    facet_df,
                    x="Facet",
                    y="Estimate",
                    color="Facet",
                    points="all",
                )
                fig.update_layout(showlegend=False, height=360)
                st.plotly_chart(fig, use_container_width=True)

        with vtabs[3]:
            render_apa_figure_note("figure4")
            st.subheader("Step / threshold estimates")
            step_tbl = res["steps"].copy()
            if step_tbl.empty:
                st.info("No step estimates available.")
            else:
                step_tbl["StepIndex"] = pd.to_numeric(
                    step_tbl["Step"].str.extract(r"(\d+)")[0],
                    errors="coerce",
                )
                if "StepFacet" not in step_tbl.columns:
                    step_tbl["StepFacet"] = "Common"
                plot_tbl = step_tbl.dropna(subset=["StepIndex"]).copy()
                plot_tbl["StepIndex"] = plot_tbl["StepIndex"].astype(int)
                fig = px.line(
                    plot_tbl,
                    x="StepIndex",
                    y="Estimate",
                    color="StepFacet",
                    markers=True,
                )
                fig.update_layout(xaxis_title="Step", yaxis_title="Estimate", height=360)
                st.plotly_chart(fig, use_container_width=True)

        with vtabs[4]:
            render_apa_figure_note("figure5")
            st.subheader("Category probability curves")
            theta_grid = np.linspace(-4, 4, 81)
            params = res["params"]
            if res["config"]["model"] == "RSM":
                step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
                probs = category_prob_rsm(theta_grid, step_cum)
            else:
                step_mean = np.nanmean(params["steps_mat"], axis=0) if params["steps_mat"] is not None else np.array([])
                step_cum = np.concatenate([[0.0], np.cumsum(step_mean)]) if step_mean.size else np.array([0.0])
                probs = category_prob_rsm(theta_grid, step_cum)
            fig = go.Figure()
            for k in range(probs.shape[1]):
                fig.add_trace(
                    go.Scatter(
                        x=theta_grid,
                        y=probs[:, k],
                        mode="lines",
                        name=f"Cat {k + res['prep']['rating_min']}",
                    )
                )
            fig.update_layout(title="Category probability curves", xaxis_title="Theta", yaxis_title="Probability", height=360)
            st.plotly_chart(fig, use_container_width=True)

        with vtabs[5]:
            render_apa_figure_note("figure6")
            st.subheader("Observed vs expected score")
            cols = ["Observed", "Expected"]
            if "Weight" in diagnostics["obs"].columns:
                cols.append("Weight")
            obs_exp = diagnostics["obs"][cols].dropna()
            if obs_exp.empty:
                st.info("No observed/expected data available.")
            else:
                try:
                    obs_exp["bin"] = pd.qcut(obs_exp["Expected"], 10, duplicates="drop")
                except ValueError:
                    obs_exp["bin"] = pd.cut(obs_exp["Expected"], 10)
                if "Weight" in obs_exp.columns:
                    plot_df = (
                        obs_exp.groupby("bin")
                        .apply(lambda g: pd.Series({
                            "Expected": weighted_mean(g["Expected"].to_numpy(), g["Weight"].to_numpy()),
                            "Observed": weighted_mean(g["Observed"].to_numpy(), g["Weight"].to_numpy()),
                            "n": np.nansum(g["Weight"].to_numpy()),
                        }))
                        .reset_index(drop=True)
                    )
                else:
                    plot_df = (
                        obs_exp.groupby("bin")
                        .agg(Expected=("Expected", "mean"), Observed=("Observed", "mean"), n=("Observed", "size"))
                        .reset_index(drop=True)
                    )
                fig = px.scatter(
                    plot_df,
                    x="Expected",
                    y="Observed",
                    size="n",
                    size_max=10,
                )
                fig.add_trace(
                    go.Scatter(
                        x=plot_df["Expected"],
                        y=plot_df["Observed"],
                        mode="lines",
                        line=dict(color="#1f77b4"),
                        showlegend=False,
                    )
                )
                fig.add_trace(
                    go.Scatter(
                        x=[plot_df["Expected"].min(), plot_df["Expected"].max()],
                        y=[plot_df["Expected"].min(), plot_df["Expected"].max()],
                        mode="lines",
                        line=dict(color="gray", dash="dash"),
                        showlegend=False,
                    )
                )
                fig.update_layout(xaxis_title="Expected score", yaxis_title="Observed score", height=360)
                st.plotly_chart(fig, use_container_width=True)

        with vtabs[6]:
            render_apa_figure_note("figure7")
            st.subheader("Fit scatter (Infit vs Outfit)")
            fit_df = diagnostics["measures"].copy()
            fit_df = fit_df.dropna(subset=["Infit", "Outfit"])
            if fit_df.empty:
                st.info("No fit statistics available.")
            else:
                facets = sorted(fit_df["Facet"].unique())
                selected = st.multiselect("Facets to show", facets, default=facets, key="fit_scatter_facets")
                if selected:
                    fit_df = fit_df[fit_df["Facet"].isin(selected)]
                if "N" not in fit_df.columns:
                    fit_df["N"] = 1
                fig = px.scatter(
                    fit_df,
                    x="Infit",
                    y="Outfit",
                    color="Facet",
                    size="N",
                    hover_data=["Level", "Infit", "Outfit", "N"],
                )
                fig.add_vline(x=1, line_dash="dash", line_color="gray")
                fig.add_hline(y=1, line_dash="dash", line_color="gray")
                fig.update_layout(xaxis_title="Infit MNSQ", yaxis_title="Outfit MNSQ", height=360)
                st.plotly_chart(fig, use_container_width=True)

        with vtabs[7]:
            render_apa_figure_note("figure8")
            st.subheader("ZSTD distribution by facet")
            zstd_df = diagnostics["measures"].copy()
            zstd_df = zstd_df[["Facet", "InfitZSTD", "OutfitZSTD"]].melt(
                id_vars=["Facet"],
                value_vars=["InfitZSTD", "OutfitZSTD"],
                var_name="Statistic",
                value_name="ZSTD",
            )
            zstd_df = zstd_df[np.isfinite(zstd_df["ZSTD"])]
            if zstd_df.empty:
                st.info("No ZSTD data available.")
            else:
                fig = px.histogram(
                    zstd_df,
                    x="ZSTD",
                    color="Statistic",
                    facet_col="Facet",
                    facet_col_wrap=3,
                    histnorm="probability density",
                    opacity=0.5,
                )
                fig.update_layout(height=360)
                st.plotly_chart(fig, use_container_width=True)

        with vtabs[8]:
            render_apa_figure_note("figure9")
            st.subheader("Top misfit levels (max |ZSTD|)")
            misfit_df = diagnostics["measures"].copy()
            if misfit_df.empty:
                st.info("No misfit statistics available.")
            else:
                top_n = st.slider("Top N levels", min_value=5, max_value=100, value=20, step=5)
                threshold = st.slider("Misfit |ZSTD| threshold", min_value=0.0, max_value=5.0, value=2.0, step=0.5)
                compare = st.radio("Threshold mode", [">=", "<="], index=0, horizontal=True)
                misfit_df["AbsZSTD"] = np.nanmax(
                    np.vstack([misfit_df["InfitZSTD"].abs().to_numpy(), misfit_df["OutfitZSTD"].abs().to_numpy()]),
                    axis=0,
                )
                misfit_df = misfit_df[np.isfinite(misfit_df["AbsZSTD"])]
                if compare == ">=":
                    misfit_df = misfit_df[misfit_df["AbsZSTD"] >= threshold]
                else:
                    misfit_df = misfit_df[misfit_df["AbsZSTD"] <= threshold]
                misfit_df = misfit_df.sort_values("AbsZSTD", ascending=False).head(top_n)
                if misfit_df.empty:
                    st.info("No misfit levels for the current threshold.")
                else:
                    misfit_df["Label"] = misfit_df.apply(
                        lambda r: truncate_label(f"{r['Facet']}: {r['Level']}", 32),
                        axis=1,
                    )
                    fig = px.bar(
                        misfit_df[::-1],
                        x="AbsZSTD",
                        y="Label",
                        color="Facet",
                        orientation="h",
                    )
                    fig.add_vline(x=2, line_dash="dash", line_color="gray")
                    fig.update_layout(xaxis_title="|ZSTD|", yaxis_title="", height=420)
                    st.plotly_chart(fig, use_container_width=True)

    with tabs[15]:
        st.caption("Export tables and diagnostic files.")
        st.subheader("Downloads")
        summary_csv = res["summary"].to_csv(index=False).encode("utf-8")
        measures_csv = diagnostics["measures"].to_csv(index=False).encode("utf-8")
        fit_csv = diagnostics["fit"].to_csv(index=False).encode("utf-8")
        steps_csv = res["steps"].to_csv(index=False).encode("utf-8")
        scorefile_csv = compute_scorefile(res).to_csv(index=False).encode("utf-8")
        residuals_csv = compute_residual_file(res).to_csv(index=False).encode("utf-8")
        report_tbls = calc_facets_report_tbls(
            res,
            diagnostics,
            totalscore=report_totalscore,
            umean=report_umean,
            uscale=report_uscale,
            udecimals=report_udecimals,
            omit_unobserved=report_omit_unobserved,
            xtreme=report_xtreme,
        )
        agreement_facet = st.session_state.get("agreement_facet")
        agreement_tbls = None
        if agreement_facet:
            agreement_tbls = calc_interrater_agreement(
                diagnostics["obs"],
                ["Person"] + res["config"]["facet_names"],
                agreement_facet,
                res=res,
            )
        subset_tbls = calc_subsets(diagnostics["obs"], ["Person"] + res["config"]["facet_names"])

        st.download_button("Summary (CSV)", summary_csv, file_name="mfrm_summary.csv")
        st.download_button("Measures (CSV)", measures_csv, file_name="mfrm_measures.csv")
        st.download_button("Fit statistics (CSV)", fit_csv, file_name="mfrm_fit.csv")
        st.download_button("Steps (CSV)", steps_csv, file_name="mfrm_steps.csv")
        if report_tbls:
            combined = []
            for facet, df in report_tbls.items():
                tmp = df.copy()
                tmp.insert(0, "Facet", facet)
                combined.append(tmp)
            combined_df = pd.concat(combined, ignore_index=True)
            report_csv = combined_df.to_csv(index=False).encode("utf-8")
            st.download_button("Facets-style report (CSV)", report_csv, file_name="mfrm_facets_report_all.csv")
        if agreement_tbls and not agreement_tbls["summary"].empty:
            st.download_button(
                "Agreement summary (CSV)",
                agreement_tbls["summary"].to_csv(index=False).encode("utf-8"),
                file_name="mfrm_agreement_summary.csv",
            )
            st.download_button(
                "Agreement pairs (CSV)",
                agreement_tbls["pairs"].to_csv(index=False).encode("utf-8"),
                file_name="mfrm_agreement_pairs.csv",
            )
        if subset_tbls and not subset_tbls["summary"].empty:
            st.download_button(
                "Subsets summary (CSV)",
                subset_tbls["summary"].to_csv(index=False).encode("utf-8"),
                file_name="mfrm_subsets_summary.csv",
            )
            st.download_button(
                "Subsets nodes (CSV)",
                subset_tbls["nodes"].to_csv(index=False).encode("utf-8"),
                file_name="mfrm_subsets_nodes.csv",
            )
        bias_results = st.session_state.get("bias_results")
        if bias_results and bias_results.get("table") is not None and not bias_results["table"].empty:
            st.download_button(
                "Bias/Interaction report (CSV)",
                bias_results["table"].to_csv(index=False).encode("utf-8"),
                file_name="mfrm_bias_report.csv",
            )
            st.download_button(
                "Bias summary (CSV)",
                bias_results["summary"].to_csv(index=False).encode("utf-8"),
                file_name="mfrm_bias_summary.csv",
            )
            if bias_results.get("iteration") is not None and not bias_results["iteration"].empty:
                st.download_button(
                    "Bias iteration report (CSV)",
                    bias_results["iteration"].to_csv(index=False).encode("utf-8"),
                    file_name="mfrm_bias_iteration.csv",
                )
            pairwise_tbl = calc_bias_pairwise(
                bias_results["table"],
                bias_results["facet_a"],
                bias_results["facet_b"],
            )
            if not pairwise_tbl.empty:
                st.download_button(
                    "Bias pairwise report (CSV)",
                    pairwise_tbl.to_csv(index=False).encode("utf-8"),
                    file_name="mfrm_bias_pairwise.csv",
                )
            # Fixed-width text outputs
            facet_a = bias_results["facet_a"]
            facet_b = bias_results["facet_b"]
            tbl = bias_results["table"].copy()
            tbl_display = tbl[[
                "Sq",
                "Observd Score",
                "Expctd Score",
                "Observd Count",
                "Obs-Exp Average",
                "Bias Size",
                "S.E.",
                "t",
                "d.f.",
                "Prob.",
                "Infit",
                "Outfit",
                "FacetA_Index",
                "FacetA_Level",
                "FacetA_Measure",
                "FacetB_Index",
                "FacetB_Level",
                "FacetB_Measure",
            ]].rename(columns={
                "S.E.": "Model S.E.",
                "Infit": "Infit MnSq",
                "Outfit": "Outfit MnSq",
                "FacetA_Index": f"{facet_a} N",
                "FacetA_Level": facet_a,
                "FacetA_Measure": f"{facet_a} measr",
                "FacetB_Index": f"{facet_b} N",
                "FacetB_Level": facet_b,
                "FacetB_Measure": f"{facet_b} measr",
            })
            bias_cols = [
                "Sq",
                "Observd Score",
                "Expctd Score",
                "Observd Count",
                "Obs-Exp Average",
                "Bias Size",
                "Model S.E.",
                "t",
                "d.f.",
                "Prob.",
                "Infit MnSq",
                "Outfit MnSq",
                f"{facet_a} N",
                facet_a,
                f"{facet_a} measr",
                f"{facet_b} N",
                facet_b,
                f"{facet_b} measr",
            ]
            bias_formats = {
                "Sq": "{}",
                "Observd Score": "{:.2f}",
                "Expctd Score": "{:.2f}",
                "Observd Count": "{:.0f}",
                "Obs-Exp Average": "{:.2f}",
                "Bias Size": "{:.2f}",
                "Model S.E.": "{:.2f}",
                "t": "{:.2f}",
                "d.f.": "{:.0f}",
                "Prob.": "{:.4f}",
                "Infit MnSq": "{:.2f}",
                "Outfit MnSq": "{:.2f}",
                f"{facet_a} N": "{:.0f}",
                f"{facet_a} measr": "{:.2f}",
                f"{facet_b} N": "{:.0f}",
                f"{facet_b} measr": "{:.2f}",
            }
            bias_fixed = build_bias_fixed_text(
                tbl_display,
                bias_results["summary"],
                bias_results["chi_sq"],
                facet_a,
                facet_b,
                bias_cols,
                bias_formats,
            )
            st.download_button(
                "Bias/Interaction report (fixed-width TXT)",
                bias_fixed.encode("utf-8"),
                file_name="mfrm_bias_report_fixed.txt",
            )
            if not pairwise_tbl.empty:
                pair_display = pairwise_tbl[[
                    "Target N",
                    "Target",
                    "Target Measure",
                    "Target S.E.",
                    "Context1 N",
                    "Context1",
                    "Local Measure1",
                    "SE1",
                    "Obs-Exp Avg1",
                    "Count1",
                    "Context2 N",
                    "Context2",
                    "Local Measure2",
                    "SE2",
                    "Obs-Exp Avg2",
                    "Count2",
                    "Contrast",
                    "SE",
                    "t",
                    "d.f.",
                    "Prob.",
                ]].rename(columns={
                    "Target Measure": "Target Measr",
                    "Local Measure1": "Context1 Measr",
                    "SE1": "Context1 S.E.",
                    "Local Measure2": "Context2 Measr",
                    "SE2": "Context2 S.E.",
                })
                pair_cols = [
                    "Target N",
                    "Target",
                    "Target Measr",
                    "Target S.E.",
                    "Context1 N",
                    "Context1",
                    "Context1 Measr",
                    "Context1 S.E.",
                    "Obs-Exp Avg1",
                    "Count1",
                    "Context2 N",
                    "Context2",
                    "Context2 Measr",
                    "Context2 S.E.",
                    "Obs-Exp Avg2",
                    "Count2",
                    "Contrast",
                    "SE",
                    "t",
                    "d.f.",
                    "Prob.",
                ]
                pair_formats = {
                    "Target N": "{:.0f}",
                    "Target Measr": "{:.2f}",
                    "Target S.E.": "{:.2f}",
                    "Context1 N": "{:.0f}",
                    "Context1 Measr": "{:.2f}",
                    "Context1 S.E.": "{:.2f}",
                    "Obs-Exp Avg1": "{:.2f}",
                    "Count1": "{:.0f}",
                    "Context2 N": "{:.0f}",
                    "Context2 Measr": "{:.2f}",
                    "Context2 S.E.": "{:.2f}",
                    "Obs-Exp Avg2": "{:.2f}",
                    "Count2": "{:.0f}",
                    "Contrast": "{:.2f}",
                    "SE": "{:.2f}",
                    "t": "{:.2f}",
                    "d.f.": "{:.0f}",
                    "Prob.": "{:.4f}",
                }
                pair_fixed = build_pairwise_fixed_text(
                    pair_display,
                    bias_results["facet_a"],
                    bias_results["facet_b"],
                    pair_cols,
                    pair_formats,
                )
                st.download_button(
                    "Bias pairwise report (fixed-width TXT)",
                    pair_fixed.encode("utf-8"),
                    file_name="mfrm_bias_pairwise_fixed.txt",
                )
        st.download_button("Scorefile (CSV)", scorefile_csv, file_name="mfrm_scorefile.csv")
        st.download_button("Residuals (CSV)", residuals_csv, file_name="mfrm_residuals.csv")
    with tabs[16]:
        st.subheader("Help & guidance")
        help_tabs = st.tabs(["Basics", "Interpretation", "Troubleshooting", "Glossary", "Reporting (APA)", "Examples"])

        with help_tabs[0]:
            st.markdown(
                """
                **Quick start**
                1. Choose a data source (sample, paste, or upload).
                2. Select Person, Score, and facet columns (2+ facets).
                3. Pick a model (RSM or PCM) and estimation (JMLE or MML).
                4. Optional: anchors, group anchors, non-centering, dummy facets.
                5. Click **Run estimation**.

                **Data format**
                - One row per observation.
                - Required: Person, Score, and 2+ facet columns (e.g., Rater, Task, Criterion).
                - Optional: Weight (positive numeric).
                - Scores must be ordered categories (0 or positive integers). If you use labels, recode them to ordered numbers.
                - Use "Keep original category values (K)" if your categories are non-contiguous or if you want to retain unobserved intermediate categories.

                **Model choices**
                - **RSM**: shared step structure across facets.
                - **PCM**: step structure varies by a chosen facet.
                - **JMLE**: closest to traditional many-facet output.
                - **MML**: EAP/SD for person estimates (experimental).

                **Design and sample size (rules of thumb)**
                - Aim for **30+ observations per element** and **10+ per category** for stable estimates.
                - Ensure facets are sufficiently crossed to avoid disjoint subsets.
                - Use shared raters/tasks or anchors to link separate administrations or cohorts.
                - If you need subgroup summaries (e.g., classes, schools), add a dummy facet and report it as a group.
                """
            )

        with help_tabs[1]:
            st.markdown(
                """
                **What the measures mean**
                - Measures are on a logit scale. Higher **Person** measures indicate higher performance; higher **Rater** measures indicate stricter or more severe raters.
                - "Fair" scores adjust for the average severity/difficulty of other facets, which helps compare elements on a common scale.

                **Fit statistics (Infit/Outfit)**
                - Values close to 1.0 indicate expected behavior. Values **> 1** suggest unexpected ratings (misfit); values **< 1** suggest overly predictable ratings (overfit).
                - A commonly used working range is about **0.5-1.5**; high-stakes uses often prefer a tighter band.
                - Infit is more sensitive to in-range patterns; outfit is more sensitive to outliers.
                - Examine high mean-squares first; low values can reflect redundancy or overfit.
                - ZSTD tests exact fit and can inflate with large samples; treat it as a flag, not a verdict.

                **Rating scale diagnostics**
                - Check category counts. Very low counts or unused categories often produce unstable thresholds.
                - Category probability curves should show distinct peaks; disordered thresholds can indicate problematic categories or rater confusion.

                **Reliability vs. agreement**
                - **Separation/Reliability** describes how distinguishable elements are within a facet (ordering reproducibility).
                  - For **Persons**, higher reliability is desirable (better differentiation).
                  - For **Raters**, lower reliability can be desirable if the goal is consistent severity across raters.
                - **Agreement** indices show how often raters agree under identical conditions; very high agreement can indicate "rating machines" and inflate separation reliability.
                - Model vs. Real reliability provide upper/lower bounds when unexpectedness is small/large.

                **Bias / interaction**
                - Bias terms capture systematic differences between pairs of facets (e.g., a rater being lenient for a specific task).
                - Bias estimation fixes the main-effect measures, then analyzes residuals to detect interaction patterns.
                - Interpret **Bias Size** with **t / p** and the **Obs-Exp Average**; large, significant values may indicate rater effects such as differential leniency or task-specific bias.

                **Research-based rater effects (common in performance assessment)**
                - **Severity/leniency**: consistent harshness or generosity in scoring.
                - **Central tendency**: avoidance of extreme categories (often low fit or compressed category usage).
                - **Randomness**: inconsistent scoring (often high fit statistics).
                - **Halo**: ratings across criteria look too similar for a given person.
                - **Differential leniency**: rater is lenient only for specific tasks or criteria (bias interaction).
                """
            )

        with help_tabs[2]:
            st.markdown(
                """
                **Common errors and fixes**
                - **Duplicate column names** or choosing Person/Score as a facet -> rename columns and reselect facets.
                - **Non-numeric scores** -> ensure Score is numeric and ordered.
                - **Sparse categories** -> collapse categories or collect more data.
                - **Disordered thresholds** -> review category labels and rater training.
                - **Extreme scores only** -> consider Xtreme correction; interpret with caution.
                - **Disconnected subsets** -> add linking data or anchors.

                **If estimates do not converge**
                - Check that exactly one facet is non-centered.
                - Reduce model complexity (fewer facets, simpler rating scale) and re-run.
                - Consider anchoring or removing problematic elements.

                **If fit looks problematic**
                - High mean-squares often reflect miscodes, ambiguous rubrics, or outlier raters/tasks; investigate those first.
                - If only a few elements misfit, compare runs with/without them and consider anchoring them back in.
                - Central tendency can show up as low infit/outfit and heavy use of middle categories.
                - For bias tables, set a practical effect-size threshold (e.g., half a category) in addition to p-values.
                """
            )

        with help_tabs[3]:
            st.markdown("**Glossary**")
            term_query = st.text_input("Search terms", value="", help="Type a keyword, e.g., infit, anchor, subset.")
            glossary = [
                {
                    "category": "Core Concepts",
                    "term": "Facet",
                    "definition": "A set of elements such as persons, raters, tasks, or criteria that contribute variability.",
                    "interpretation": "Choose facets that represent the main sources of variation you want to adjust for.",
                },
                {
                    "category": "Core Concepts",
                    "term": "Element",
                    "definition": "A member within a facet, identified by a label and index.",
                    "interpretation": "Each element receives a separate measure and fit statistics.",
                },
                {
                    "category": "Core Concepts",
                    "term": "Observation",
                    "definition": "A recorded rating or count for a specific combination of elements.",
                    "interpretation": "One row in the data corresponds to one observation.",
                },
                {
                    "category": "Core Concepts",
                    "term": "Measure (logit)",
                    "definition": "A linear estimate for an element on the logit scale.",
                    "interpretation": "Higher person measures indicate higher performance; higher rater measures indicate greater severity.",
                },
                {
                    "category": "Core Concepts",
                    "term": "Standard Error (SE)",
                    "definition": "An estimate of precision for an element's measure.",
                    "interpretation": "Smaller SE means more precise measurement.",
                },
                {
                    "category": "Model and Estimation",
                    "term": "Many-Facet Rasch Model (MFRM)",
                    "definition": "A Rasch model that estimates measures for multiple facets simultaneously.",
                    "interpretation": "Use MFRM when ratings depend on more than persons and items.",
                },
                {
                    "category": "Model and Estimation",
                    "term": "Rating Scale Model (RSM)",
                    "definition": "A model with a common set of thresholds across elements.",
                    "interpretation": "Use when the rating scale functions similarly across the chosen step facet.",
                },
                {
                    "category": "Model and Estimation",
                    "term": "Partial Credit Model (PCM)",
                    "definition": "A model with thresholds that vary by a specified facet.",
                    "interpretation": "Use when different tasks or criteria have distinct step structures.",
                },
                {
                    "category": "Model and Estimation",
                    "term": "JMLE",
                    "definition": "Joint maximum likelihood estimation for facet measures and steps.",
                    "interpretation": "Standard for many-facet output and familiar FACETS-style tables.",
                },
                {
                    "category": "Model and Estimation",
                    "term": "MML / EAP",
                    "definition": "Marginal maximum likelihood with EAP estimates for persons.",
                    "interpretation": "Use when you want posterior means and SDs for persons.",
                },
                {
                    "category": "Rating Scale and Thresholds",
                    "term": "Rating Scale",
                    "definition": "An ordered set of categories used for scoring.",
                    "interpretation": "Categories should be used meaningfully and in order.",
                },
                {
                    "category": "Rating Scale and Thresholds",
                    "term": "Category",
                    "definition": "An ordered score level on the rating scale (e.g., 1-5).",
                    "interpretation": "Sparse or unused categories can destabilize thresholds.",
                },
                {
                    "category": "Rating Scale and Thresholds",
                    "term": "Rasch-Andrich Thresholds (Step Calibrations)",
                    "definition": "Threshold values between adjacent categories on the latent scale.",
                    "interpretation": "Ordered thresholds suggest categories function as intended.",
                },
                {
                    "category": "Rating Scale and Thresholds",
                    "term": "Step Calibration",
                    "definition": "Another name for Rasch-Andrich thresholds in rating/partial credit scales.",
                    "interpretation": "Use to check whether category transitions are distinct.",
                },
                {
                    "category": "Rating Scale and Thresholds",
                    "term": "General vs Specific Scale",
                    "definition": "General uses one set of thresholds; Specific estimates separate thresholds per model statement.",
                    "interpretation": "Specific allows item- or task-specific category structures.",
                },
                {
                    "category": "Rating Scale and Thresholds",
                    "term": "Category Recoding",
                    "definition": "Mapping multiple raw values to a single category during input.",
                    "interpretation": "Use to collapse sparse categories or harmonize labels.",
                },
                {
                    "category": "Rating Scale and Thresholds",
                    "term": "Anchored Category",
                    "definition": "A category whose threshold is fixed at a preset value.",
                    "interpretation": "Anchoring helps equate scales across analyses.",
                },
                {
                    "category": "Rating Scale and Thresholds",
                    "term": "Disordered Thresholds",
                    "definition": "Thresholds that do not increase in order across categories.",
                    "interpretation": "Consider collapsing or revising categories.",
                },
                {
                    "category": "Fit and Diagnostics",
                    "term": "Infit MnSq",
                    "definition": "Information-weighted fit mean-square statistic.",
                    "interpretation": "Values near 1.0 indicate expected fit; >1 suggests underfit; <1 suggests overfit.",
                },
                {
                    "category": "Fit and Diagnostics",
                    "term": "Outfit MnSq",
                    "definition": "Outlier-sensitive fit mean-square statistic.",
                    "interpretation": "High outfit can indicate unexpected outliers or unusual ratings.",
                },
                {
                    "category": "Fit and Diagnostics",
                    "term": "ZSTD",
                    "definition": "Standardized (z) fit statistics for infit and outfit mean-squares.",
                    "interpretation": "Large |ZSTD| flags unexpected responses; it can inflate with large samples.",
                },
                {
                    "category": "Fit and Diagnostics",
                    "term": "WHEXACT",
                    "definition": "Wilson-Hilferty exact standardization option for ZSTD.",
                    "interpretation": "Use when you want the exact transformation rather than the approximation.",
                },
                {
                    "category": "Fit and Diagnostics",
                    "term": "Residual (Yni)",
                    "definition": "Observed step value minus expected value for an observation.",
                    "interpretation": "Large residuals indicate unexpected ratings.",
                },
                {
                    "category": "Fit and Diagnostics",
                    "term": "Standardized Residual (Zni)",
                    "definition": "A residual scaled by its model-based standardization.",
                    "interpretation": "Used to flag unexpected observations in Table 4.",
                },
                {
                    "category": "Fit and Diagnostics",
                    "term": "Observed Score",
                    "definition": "The raw score or count recorded in the data.",
                    "interpretation": "Use with expected score to compute residuals and fit.",
                },
                {
                    "category": "Fit and Diagnostics",
                    "term": "Expected Score",
                    "definition": "The model-implied average score given the measures and thresholds.",
                    "interpretation": "Comparing observed vs expected highlights misfit.",
                },
                {
                    "category": "Diagnostics and Formulas",
                    "term": "RMSE",
                    "definition": "Root mean-square error, the average of standard errors.",
                    "interpretation": "Used with True SD to compute separation and reliability.",
                },
                {
                    "category": "Diagnostics and Formulas",
                    "term": "Adj (True) SD",
                    "definition": "Standard deviation adjusted for measurement error.",
                    "interpretation": "Represents the estimated true spread of measures.",
                },
                {
                    "category": "Diagnostics and Formulas",
                    "term": "Separation (formula)",
                    "definition": "Separation = True SD / RMSE.",
                    "interpretation": "Higher values mean better differentiation among elements.",
                },
                {
                    "category": "Diagnostics and Formulas",
                    "term": "Reliability (formula)",
                    "definition": "Reliability = True variance / Observed variance.",
                    "interpretation": "The Rasch equivalent of KR-20 or Cronbach alpha.",
                },
                {
                    "category": "Diagnostics and Formulas",
                    "term": "Strata (formula)",
                    "definition": "Strata = (4 * Separation + 1) / 3.",
                    "interpretation": "Approximate number of statistically distinct levels.",
                },
                {
                    "category": "Diagnostics and Formulas",
                    "term": "Real S.E. (formula)",
                    "definition": "Real S.E. = Model S.E. * sqrt(max(Infit MnSq, 1)).",
                    "interpretation": "Inflates error when misfit is present.",
                },
                {
                    "category": "Diagnostics and Formulas",
                    "term": "Unexpected (standardized residual threshold)",
                    "definition": "Cutoff for listing observations in Table 4.",
                    "interpretation": "About 5% of residuals exceed |2| and about 1% exceed |3| when data fit the model.",
                },
                {
                    "category": "Reliability and Agreement",
                    "term": "Separation",
                    "definition": "The spread of measures relative to measurement error.",
                    "interpretation": "Higher separation means elements are more distinguishable.",
                },
                {
                    "category": "Reliability and Agreement",
                    "term": "Reliability (Separation)",
                    "definition": "The proportion of observed variance not due to measurement error.",
                    "interpretation": "Higher values indicate more stable ordering of elements.",
                },
                {
                    "category": "Reliability and Agreement",
                    "term": "Strata",
                    "definition": "Approximate number of statistically distinct levels implied by separation.",
                    "interpretation": "Larger strata implies clearer distinctions among elements.",
                },
                {
                    "category": "Reliability and Agreement",
                    "term": "Model Reliability vs Real Reliability",
                    "definition": "Model treats all unexpectedness as random; Real treats it as misfit.",
                    "interpretation": "True reliability is expected to fall between the two.",
                },
                {
                    "category": "Reliability and Agreement",
                    "term": "Model S.E. vs Real S.E.",
                    "definition": "Model S.E. assumes Rasch-consistent randomness; Real S.E. allows extra randomness.",
                    "interpretation": "Real S.E. is typically larger in the presence of misfit.",
                },
                {
                    "category": "Reliability and Agreement",
                    "term": "Agreement (Exact/Expected)",
                    "definition": "Observed exact agreement compared to expected agreement under the model.",
                    "interpretation": "Much higher-than-expected agreement can indicate dependence.",
                },
                {
                    "category": "Fair Averages and Scores",
                    "term": "Fair Average",
                    "definition": "Expected average raw score if other facets are held at mean or zero levels.",
                    "interpretation": "Use to compare elements as if they faced the same conditions.",
                },
                {
                    "category": "Fair Averages and Scores",
                    "term": "Fair(M)",
                    "definition": "Fair Average computed using mean measures of other facets.",
                    "interpretation": "This is the default fair-average reference.",
                },
                {
                    "category": "Fair Averages and Scores",
                    "term": "Fair(Z)",
                    "definition": "Fair Average computed using the zero/origin of other facets.",
                    "interpretation": "Use when you want a criterion-referenced baseline.",
                },
                {
                    "category": "Fair Averages and Scores",
                    "term": "Fair Score",
                    "definition": "A raw-score metric version of the measure for communication.",
                    "interpretation": "Useful for audiences who prefer the original rating scale.",
                },
                {
                    "category": "Bias and Interaction",
                    "term": "Bias/Interaction",
                    "definition": "A facet-by-facet interaction beyond main effects.",
                    "interpretation": "Indicates systematic differences for a specific facet pair.",
                },
                {
                    "category": "Bias and Interaction",
                    "term": "Bias Size",
                    "definition": "Estimated magnitude of the interaction effect (logits).",
                    "interpretation": "Larger absolute values indicate stronger interaction.",
                },
                {
                    "category": "Bias and Interaction",
                    "term": "Obs-Exp Average",
                    "definition": "Observed average score minus expected average score.",
                    "interpretation": "Positive values indicate observed scores exceed expectations.",
                },
                {
                    "category": "Anchoring and Linking",
                    "term": "Anchor (Element)",
                    "definition": "A fixed element value used to link analyses.",
                    "interpretation": "Anchors stabilize measures across administrations.",
                },
                {
                    "category": "Anchoring and Linking",
                    "term": "Anchorfile",
                    "definition": "An output specification file containing anchored measures.",
                    "interpretation": "Use to reuse measures as anchors in later runs.",
                },
                {
                    "category": "Anchoring and Linking",
                    "term": "Anchoring Rating Scales",
                    "definition": "Fixing Rasch-Andrich thresholds to equate scales across analyses.",
                    "interpretation": "Ensures comparability of category structures between runs.",
                },
                {
                    "category": "Anchoring and Linking",
                    "term": "Group Anchor",
                    "definition": "A fixed group mean rather than individual elements.",
                    "interpretation": "Useful when individual anchors are unstable.",
                },
                {
                    "category": "Anchoring and Linking",
                    "term": "Dummy Facet",
                    "definition": "A facet with all elements anchored at 0 used for interactions.",
                    "interpretation": "Supports interaction analyses without shifting the scale.",
                },
                {
                    "category": "Anchoring and Linking",
                    "term": "Non-centered Facet",
                    "definition": "A facet allowed to float rather than being centered at zero.",
                    "interpretation": "Typically used for the person facet.",
                },
                {
                    "category": "Extreme Scores and Status Flags",
                    "term": "Extreme Score",
                    "definition": "Only minimum or maximum possible scores for an element.",
                    "interpretation": "Extreme-only elements yield unstable measures.",
                },
                {
                    "category": "Extreme Scores and Status Flags",
                    "term": "Xtreme (Extreme Score Adjustment)",
                    "definition": "A fractional score adjustment used to estimate finite measures for extremes.",
                    "interpretation": "Use when extreme-only elements would otherwise be infinite.",
                },
                {
                    "category": "Extreme Scores and Status Flags",
                    "term": "Minimum / Maximum",
                    "definition": "Status flags indicating minimum or maximum possible extreme scores.",
                    "interpretation": "Fit statistics are not computed for these elements.",
                },
                {
                    "category": "Extreme Scores and Status Flags",
                    "term": "Unmeasurable",
                    "definition": "Status flag for elements not connected to estimable elements.",
                    "interpretation": "Indicates a disconnected subset requiring anchors or linking data.",
                },
                {
                    "category": "Extreme Scores and Status Flags",
                    "term": "Anchored (Status)",
                    "definition": "Status flag for elements fixed at an anchored value.",
                    "interpretation": "Anchored elements do not update during estimation.",
                },
                {
                    "category": "Extreme Scores and Status Flags",
                    "term": "Facet Orientation (+/-)",
                    "definition": "Plus or minus signs indicate whether a facet is positively or negatively oriented.",
                    "interpretation": "Check orientation before comparing directions across facets.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Output Tables Menu",
                    "definition": "Menu that exports post-analysis tables to files or screen.",
                    "interpretation": "Use to generate Table 4-8 reports.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 1 Specifications Summary",
                    "definition": "Summary of the analysis specifications.",
                    "interpretation": "Check that the input setup matches your intended design.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 2 Data Summary",
                    "definition": "Summary of the data used in the analysis.",
                    "interpretation": "Use to verify counts and missingness.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 3 Iteration Report",
                    "definition": "Iteration history of the main analysis.",
                    "interpretation": "Use to confirm convergence and stability.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 4 Unexpected Responses",
                    "definition": "Outlying observations with large standardized residuals.",
                    "interpretation": "Look for repeated elements to diagnose data issues or rater effects.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 5 Measurable Data Summary",
                    "definition": "Summary statistics for measurable data.",
                    "interpretation": "Use to check overall data coverage.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 6.0 All Facet Summary (Wright Map)",
                    "definition": "Vertical rulers that position each facet's elements by measure.",
                    "interpretation": "Compare facet ranges and alignments on a common scale.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 6.0.0 Disjoint Element Listing",
                    "definition": "Elements in disconnected subsets.",
                    "interpretation": "Use to identify linking gaps across facets.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 6.2 Facet Statistics Graphic",
                    "definition": "Graphical description of facet statistics.",
                    "interpretation": "Quick visual check of facet distributions.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 7 Measurement Report",
                    "definition": "Facet element measures with SE, fit, and fair averages.",
                    "interpretation": "Primary output for reporting measures.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 7 Reliability and Chi-square",
                    "definition": "Reliability, separation, strata, RMSE, and chi-square statistics.",
                    "interpretation": "Use to judge reproducibility of ordering.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 7 Agreement Statistics",
                    "definition": "Observed vs expected agreement for rater pairs.",
                    "interpretation": "High exact agreement can indicate rater dependence.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 8 Rating Scale Structure",
                    "definition": "Category statistics and threshold diagnostics.",
                    "interpretation": "Evaluate category functioning and ordering.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 8.1 Polytomous / Dichotomous",
                    "definition": "Rating scale or dichotomy statistics for scale structure.",
                    "interpretation": "Use with Table 8 to judge scale quality.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Table 8 Bar Chart",
                    "definition": "Scale structure bar-chart output.",
                    "interpretation": "Quick visual check of category functioning.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Total Score (Table 7)",
                    "definition": "Observed raw score including extreme elements.",
                    "interpretation": "Shown when Totalscore=Yes.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Obsvd Score (Table 7)",
                    "definition": "Observed raw score after removing extreme elements.",
                    "interpretation": "Shown when Totalscore=No.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Weightd Score (Table 7)",
                    "definition": "Observed raw score multiplied by weights.",
                    "interpretation": "Shown when weighting is active.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Total Count (Table 7)",
                    "definition": "Observed response count including extremes.",
                    "interpretation": "Count of all responses for the element.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Obsvd Count (Table 7)",
                    "definition": "Observed response count after removing extremes.",
                    "interpretation": "Use to check effective data size.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Weightd Count (Table 7)",
                    "definition": "Observed response count multiplied by weights.",
                    "interpretation": "Shown when weighting is active.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Fair(M) Average (Table 7)",
                    "definition": "Fair average using mean measures of other facets.",
                    "interpretation": "Standardizes ratings to average conditions.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Fair(Z) Average (Table 7)",
                    "definition": "Fair average using facet zero points as baseline.",
                    "interpretation": "Use for criterion-referenced interpretation.",
                },
                {
                    "category": "Output Tables (FACETS)",
                    "term": "Anchoring Status (A/G/X)",
                    "definition": "A=anchored, G=group anchor including extremes, X=group anchor excluding extremes.",
                    "interpretation": "Indicates how element measures were fixed.",
                },
            ]
            q = term_query.strip().lower()
            if q:
                filtered = [
                    g for g in glossary
                    if q in g["term"].lower()
                    or q in g["definition"].lower()
                    or q in g["interpretation"].lower()
                    or q in g["category"].lower()
                ]
            else:
                filtered = glossary
            if not filtered:
                st.info("No glossary entries match your search.")
            else:
                ordered_categories = []
                for g in filtered:
                    if g["category"] not in ordered_categories:
                        ordered_categories.append(g["category"])
                for category in ordered_categories:
                    st.markdown(f"**{category}**")
                    for g in [x for x in filtered if x["category"] == category]:
                        st.markdown(
                            f"- **{g['term']}**: {g['definition']}  \n"
                            f"  Interpretation: {g['interpretation']}"
                        )

        with help_tabs[4]:
            st.markdown("**APA 7 report builder (auto-generated from current results)**")
            with st.expander("Report inputs (optional)", expanded=True):
                assessment_name = st.text_input("Assessment name", value="", key="report_assessment")
                assessment_setting = st.text_input("Setting / context", value="", key="report_setting")
                rater_training = st.text_input("Rater training (short phrase)", value="", key="report_training")
                raters_per_response = st.text_input("Raters per response (e.g., 2 raters)", value="", key="report_raters_per")
                scale_desc = st.text_input("Scale description (e.g., 0-4 analytic scale)", value="", key="report_scale_desc")
                rater_facet_name = st.text_input("Primary rater facet name (optional)", value="", key="report_rater_facet")

            bias_results = st.session_state.get("bias_results")
            report_text = build_apa_report_text(
                res,
                diagnostics,
                bias_results=bias_results,
                context={
                    "assessment": assessment_name,
                    "setting": assessment_setting,
                    "rater_training": rater_training,
                    "raters_per_response": raters_per_response,
                    "scale_desc": scale_desc,
                },
                whexact=st.session_state.get("whexact", False),
            )
            st.text_area(
                "APA-style Method and Results (auto-generated)",
                report_text,
                height=420,
            )
            st.download_button(
                "Download APA Method/Results (TXT)",
                report_text.encode("utf-8"),
                file_name="mfrm_apa_method_results.txt",
            )

            table_notes = build_apa_table_figure_notes(
                res,
                diagnostics,
                bias_results=bias_results,
                context={"rater_facet": rater_facet_name},
                whexact=st.session_state.get("whexact", False),
            )
            st.text_area(
                "APA table and figure notes (auto-generated)",
                table_notes,
                height=380,
            )
            st.download_button(
                "Download APA table/figure notes (TXT)",
                table_notes.encode("utf-8"),
                file_name="mfrm_apa_table_figure_notes.txt",
            )

            captions = build_apa_table_figure_captions(
                res,
                diagnostics,
                bias_results=bias_results,
                context={"assessment": assessment_name},
            )
            st.text_area(
                "APA table and figure titles (auto-generated)",
                captions,
                height=300,
            )
            st.download_button(
                "Download APA table/figure titles (TXT)",
                captions.encode("utf-8"),
                file_name="mfrm_apa_table_figure_titles.txt",
            )
            st.markdown("---")

            st.markdown(
                """
                **Reporting in APA 7 style (many-facet Rasch / MFRM)**
                Use the templates below to draft concise Method and Results sections. Replace brackets with your actual values from this app.

                **Method (APA 7 template)**
                - **Participants/targets**: "The analysis included [P] persons (e.g., students) who completed [assessment name] in [setting]."
                - **Raters and tasks**: "Ratings were provided by [R] raters across [T] tasks/[K] criteria using a [scale description] rubric."
                - **Procedure**: "Raters received [training/calibration], and each response was scored by [n] raters per task on a [min-max] category scale."
                - **Data structure**: "The data file contained one row per observation with columns for Person, Score, and facet variables (e.g., Rater, Task, Criterion)."
                - **Model and estimation**: "A many-facet Rasch model was fit using a [RSM/PCM] specification with [JMLE/MML] estimation in [software/app]."
                - **Anchoring/extremes**: "Anchors were [not used/used for facet X], and extreme scores were [handled with Xtreme correction/excluded/not adjusted]."

                **Results (APA 7 template)**
                - **Rating scale functioning**: "Category usage was [adequate/uneven], with thresholds [ordered/disordered] across [k] categories."
                - **Facet measures**: "Person measures ranged from [min] to [max] logits (M = [x.xx], SD = [x.xx]). Rater measures ranged from [min] to [max], indicating [severity spread]."
                - **Fit**: "Mean infit MnSq was [x.xx] (SD = [x.xx]) and mean outfit MnSq was [x.xx] (SD = [x.xx]). [n] elements showed notable misfit."
                - **Reliability/separation**: "Person separation reliability was [x.xx], indicating [interpretation]. Rater separation reliability was [x.xx], indicating [interpretation]."
                - **Bias/interaction**: "Bias analysis identified [facet A x facet B] interactions; the largest contrast was [x.xx] logits (t = [x.xx], p = [value])."
                - **Interpretation**: "Overall, results suggest [rating quality conclusion], with practical implications for [training, rubric revision, or task design]."

                **Step-by-step reporting flow**
                1. **Design and data**: Describe the assessment, facets, scale, and observation counts.
                2. **Model and estimation**: State MFRM model (RSM/PCM), estimation (JMLE/MML), and software settings.
                3. **Rating scale functioning**: Summarize category counts and threshold ordering.
                4. **Facet measures and reliability**: Report measures, SE, separation, and reliability per facet.
                5. **Fit**: Summarize mean infit/outfit and identify notable misfit elements.
                6. **Bias/interaction**: Report meaningful interactions with effect sizes and significance.
                7. **Interpretation**: Relate findings to scoring quality and practical decisions.

                **Table/figure suggestions (APA-friendly)**
                - Table: Facet summary (measure, SE, infit/outfit, separation, reliability).
                - Table: Rater fit summary with counts and severity range.
                - Table: Significant bias interactions (effect size, t, p, Obs-Exp).
                - Figure: Wright map (persons vs. facets).
                - Figure: Category probability curves for the rating scale.
                """
            )

        with help_tabs[5]:
            st.markdown(
                """
                **Interpretation examples**
                - **Leniency/Severity**: A rater has a much lower measure than peers -> more lenient.
                - **Central tendency**: Rater rarely uses extreme categories -> low counts at ends, disordered thresholds.
                - **Randomness**: Rater shows high infit/outfit -> inconsistent use of the scale.
                - **Halo**: Ratings are overly consistent across traits -> trait-level differences shrink.
                - **Differential leniency**: A rater is lenient only for one task -> bias/interaction flagged with large t.

                **Practical workflow**
                1. Start with RSM + JMLE.
                2. Check category counts and step ordering.
                3. Review fit (Infit/Outfit) at facet and element levels.
                4. Inspect rater effects using bias/interaction.
                5. If needed, refine scale or rater training and re-run.

                **Q&A (English education examples)**
                - **Q: Our speaking test has multiple raters. One rater looks much stricter. What should we check?**
                  **A:** Compare rater measures and fit. A higher rater measure indicates more severe scoring. If fit is also poor, the rater may be inconsistent. Review their score distribution and consider calibration training.
                - **Q: Writing Task A seems harder than Task B. How does that show up?**
                  **A:** Task A should have a higher measure (more difficult). If person measures drop sharply only on Task A, check for task-specific bias or poorly aligned prompts.
                - **Q: A rubric has five categories, but almost nobody uses the top band. Is that a problem?**
                  **A:** Sparse categories can lead to unstable thresholds. Consider collapsing adjacent bands or revising descriptors so raters use the full scale.
                - **Q: Two raters are consistent overall, but disagree on “Grammar” only. What does that imply?**
                  **A:** Look for a significant bias/interaction between Rater and Criterion=Grammar. That suggests a systematic difference in how grammar is interpreted.
                - **Q: Some students only receive minimum scores on a short diagnostic quiz. How should we interpret them?**
                  **A:** Extreme scores can yield infinite measures. Use Xtreme correction and interpret cautiously, or add more items to better differentiate low performers.
                - **Q: We want to compare two cohorts (spring vs fall) with different prompts.**
                  **A:** Use anchors (shared raters or common prompts) or group anchoring to link the scales before comparing person measures.
                - **Q: Speaking ratings are consistent, but raters avoid the top and bottom categories.**
                  **A:** That suggests central tendency. Check category counts and thresholds, and consider rubric revisions or targeted rater training.
                - **Q: In an analytic writing rubric, all criteria scores move together.**
                  **A:** That can signal halo. Look for high correlations across criteria and investigate whether descriptors are distinct enough.
                - **Q: A new rater is erratic only on Task B, not Task A.**
                  **A:** Check rater-by-task bias and the rater's fit for Task B; task-specific guidance may be needed.
                - **Q: Two teachers score classroom participation similarly, but one is harsher for high-level students.**
                  **A:** Inspect bias for rater-by-proficiency group interactions and consider whether expectations differ by level.
                - **Q: A reading rubric works well overall, but one criterion shows misfit.**
                  **A:** Review the wording of that criterion and check if it behaves differently across tasks; misfit can indicate construct drift.

                **Source notes (paraphrased)**
                Guidance above draws on the FACETS manual and the rater-effects overview by Myford & Wolfe (2003), which emphasize systematic rater effects, the interpretation of fit statistics, and careful use of separation/Agreement indices.
                """
            )
else:
    st.info("Run estimation to see results.")
