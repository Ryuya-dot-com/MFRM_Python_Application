"""Shared Python estimator, diagnostics and result identity; no UI runtime.

Moved from the Streamlit entrypoint without changing numerical algorithms.
The Web app imports these same functions; desktop services call them directly.
"""
from __future__ import annotations

import hashlib
import json
import re
import time
from collections import OrderedDict
from collections.abc import Mapping
from datetime import date, datetime, timedelta
from enum import Enum
from functools import lru_cache
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import LinearConstraint, minimize, root_scalar
from scipy.special import logsumexp, roots_hermitenorm
from scipy.stats import chi2, norm as _norm, t as t_dist

from mfrm_app import evidence as _evidence
from mfrm_app import jmle_boundary as _jmle_boundary
from mfrm_app import output_qualification as _output_qualification

APP_VERSION = "0.2.19-beta"


APP_RELEASE_LABEL = "standalone Python beta"


MML_COVARIANCE_AUTO_MAX_PARAMS = 180


PCA_STABILITY_MIN_PERSONS = 30


PCA_STABILITY_MIN_COLUMNS = 3


PCA_STABILITY_MAX_MISSING_SHARE = 0.50


PCA_STABILITY_MIN_PAIRWISE_N = 10


PCA_STABILITY_LOO_MAX_COLUMNS = 40


PCA_STABILITY_BOOTSTRAP_REPS = 40


PCA_STABILITY_BOOTSTRAP_MAX_COLUMNS = 40


PCA_STABILITY_BOOTSTRAP_SEED = 20260513


PCA_STABILITY_MAX_LOO_EV1_SHARE_DELTA = 0.10


PCA_STABILITY_MAX_BOOTSTRAP_EV1_CV = 0.20


PCA_STABILITY_MIN_LOADING_CORR = 0.70


def _make_reason(key: str, **kwargs) -> dict:
    """Return a serializable reason payload that defers locale resolution.

    Diagnostics that surface user-visible messages (notably ``pca_reason``)
    must keep working when the sidebar language is toggled *after* the
    estimation run. Storing the i18n key + format kwargs — rather than the
    rendered string — lets the Dimensionality tab re-translate on every
    rerun without re-running the estimator.
    """
    return {"key": key, "kwargs": dict(kwargs)}


class _LazyKeyLabel:
    """Semantic fragment; the UI resolves its locale when rendering a reason."""

    __slots__ = ("mode", "facet")

    def __init__(self, *, mode: str, facet: str | None = None) -> None:
        self.mode = mode
        self.facet = facet

    def __str__(self) -> str:
        return "item-combination" if self.mode == "overall" else f"{self.facet} levels"

    def __format__(self, spec: str) -> str:
        return str(self)


def _fingerprint_content(value: object):
    """Lossless typed content for supported settings; never use display reprs."""
    encode = _fingerprint_content
    if value is pd.NA:
        return ["pandas.NA"]
    if value is pd.NaT:
        return ["pandas.NaT"]
    if isinstance(value, np.ndarray):
        dtype = value.dtype.descr if value.dtype.fields else value.dtype.str
        content = (encode(value.ravel(order="C").tolist()) if value.dtype.hasobject
                   else value.tobytes(order="C").hex())
        return ["ndarray", encode(dtype), list(value.shape), content]
    if isinstance(value, np.generic):
        return ["numpy.scalar", encode(np.asarray(value))]
    if isinstance(value, pd.CategoricalDtype):
        return ["category", encode(value.categories), value.ordered]
    if isinstance(value, pd.api.extensions.ExtensionDtype):
        return [type(value).__name__, str(value),
                encode(getattr(value, "storage", None)), encode(value.na_value)]
    if isinstance(value, np.dtype):
        return ["dtype", encode(value.descr if value.fields else value.str)]
    if isinstance(value, pd.MultiIndex):
        return ["MultiIndex", encode(value.names), encode(list(value.levels)),
                encode(list(value.codes))]
    if isinstance(value, pd.Index):
        return ["Index", encode(value.dtype), encode(value.name), encode(value.tolist())]
    if isinstance(value, pd.Series):
        return ["Series", encode(value.name), encode(value.index), encode(value.dtype),
                encode(value.tolist())]
    if isinstance(value, pd.DataFrame):
        return ["DataFrame", encode(value.index), encode(value.columns),
                [encode(value.iloc[:, i]) for i in range(value.shape[1])]]
    if isinstance(value, Enum):
        return ["enum", type(value).__module__, type(value).__qualname__, encode(value.value)]
    if value is None or type(value) in (str, bool, int):
        return value
    if type(value) is float:
        return ["float", value.hex()]
    if isinstance(value, Mapping):
        pairs = [[encode(key), encode(item)] for key, item in value.items()]
        pairs.sort(key=lambda pair: json.dumps(pair[0], ensure_ascii=False, separators=(",", ":")))
        return ["mapping", pairs]
    if isinstance(value, (list, tuple, set, frozenset)):
        items = [encode(item) for item in value]
        if isinstance(value, (set, frozenset)):
            items.sort(key=lambda item: json.dumps(item, ensure_ascii=False, separators=(",", ":")))
        return [type(value).__name__, items]
    if isinstance(value, (pd.Timestamp, datetime, date)):
        return [type(value).__name__, value.isoformat(), str(getattr(value, "tzinfo", None))]
    if isinstance(value, pd.Timedelta):
        return ["Timedelta", value.isoformat()]
    if isinstance(value, timedelta):
        return ["timedelta", value.days, value.seconds, value.microseconds]
    if isinstance(value, Path):
        return ["path", str(value)]
    if isinstance(value, bytes):
        return ["bytes", value.hex()]
    raise TypeError(f"Unsupported fingerprint content type: {type(value).__name__}")


def stable_json_fingerprint(payload: object, length: int = 16) -> str:
    """Versioned SHA-256 of complete typed content, independent of display options."""
    text = json.dumps(
        ["mfrm_content_fingerprint_v2", _fingerprint_content(payload)],
        ensure_ascii=False, allow_nan=False, separators=(",", ":"),
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[: int(length)]


def dataframe_fingerprint(df: pd.DataFrame | None, length: int = 16) -> str | None:
    """Fingerprint DataFrame content without storing the DataFrame in Streamlit cache."""
    if df is None:
        return None
    if not isinstance(df, pd.DataFrame):
        raise TypeError("dataframe_fingerprint requires a DataFrame or None")

    return stable_json_fingerprint(df, length=length)


def config_fingerprint(config: dict | None, length: int = 16) -> str:
    """Fingerprint analysis settings while excluding derived fingerprint fields."""
    if not isinstance(config, dict):
        return stable_json_fingerprint({}, length=length)
    excluded = {
        "run_fingerprint",
        "analysis_config_fingerprint",
        "config_export_fingerprint",
    }
    payload = {k: v for k, v in config.items() if k not in excluded}
    return stable_json_fingerprint(payload, length=length)


FINAL_PCA_EIGENVALUE_READY = 2.0


FINAL_PCA_EIGENVALUE_REVIEW = 3.0


@lru_cache(maxsize=8)
def _standard_normal_quadrature(n: int):
    """Normalized standard-normal GH rule, including zero tail weights."""
    nodes, weights = roots_hermitenorm(n)
    weights /= np.sqrt(2.0 * np.pi)
    nodes.setflags(write=False)
    weights.setflags(write=False)
    return nodes, weights


def gauss_hermite_normal(n: int, sd: float = 1.0):
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError("n must be an integer >= 1")
    sd = float(sd if sd is not None else 1.0)
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("sd must be a positive finite value")
    nodes, weights = _standard_normal_quadrature(int(n))
    with np.errstate(over="ignore"):
        scaled_nodes = nodes * sd
    if not np.isfinite(scaled_nodes).all():
        raise ValueError("sd produces nonfinite quadrature nodes")
    # Existing callers receive independent arrays; mutation cannot poison the
    # cached rule or a later fit at another population SD.
    return {"nodes": scaled_nodes, "weights": weights.copy(), "sd": sd}


def _log_mml_quadrature_weights(quad):
    """Log of normalized prior masses; exact zeros stay at minus infinity."""
    nodes = np.asarray(quad["nodes"], dtype=float)
    weights = np.asarray(quad["weights"], dtype=float)
    if (nodes.ndim != 1 or nodes.size == 0 or weights.shape != nodes.shape
            or not np.isfinite(nodes).all() or not np.isfinite(weights).all()
            or np.any(weights < 0) or abs(weights.sum() - 1.0) > 1e-12):
        raise ValueError("Invalid normalized MML quadrature")
    result = np.full_like(weights, -np.inf)
    np.log(weights, out=result, where=weights > 0)
    return result


def get_population_prior_sd(config: dict | None) -> float:
    if not isinstance(config, dict):
        return 1.0
    raw = config.get("population_prior_sd", 1.0)
    sd = float(raw if raw is not None else 1.0)
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("population_prior_sd must be a positive finite value.")
    return sd


def make_mml_quadrature(config: dict, quad_points: int | None = None) -> dict:
    points = config.get("quad_points") if quad_points is None else quad_points
    if points is None:
        points = 15
    return gauss_hermite_normal(points, sd=get_population_prior_sd(config))


def center_sum_zero(x: np.ndarray):
    if x.size == 0:
        return x
    return x - np.mean(x)


def expand_sum_zero_free(free: np.ndarray, full_size: int) -> np.ndarray:
    """Expand ``full_size - 1`` coordinates to an exact sum-zero vector.

    The final coordinate is derived as the negative sum of the free block.  In
    contrast to centring a redundant ``full_size`` block, this parameterization
    has no artificial null direction in the optimizer or observed information.
    """
    full_size = int(full_size)
    if full_size < 0:
        raise ValueError("full_size must be non-negative.")
    free = np.asarray(free, dtype=float).reshape(-1)
    expected = max(full_size - 1, 0)
    if free.size != expected:
        raise ValueError(
            f"Expected {expected} free sum-zero coordinate(s) for full_size={full_size}, "
            f"got {free.size}."
        )
    if full_size == 0:
        return np.array([], dtype=float)
    if full_size == 1:
        return np.array([0.0], dtype=float)
    return np.concatenate([free, [-float(np.sum(free))]])


def sum_zero_free_from_expanded(expanded: np.ndarray) -> np.ndarray:
    """Return free coordinates for an arbitrary vector after sum-zero centring."""
    expanded = np.asarray(expanded, dtype=float).reshape(-1)
    if expanded.size <= 1:
        return np.array([], dtype=float)
    return center_sum_zero(expanded)[:-1]


def collapse_sum_zero_gradient(full_grad: np.ndarray) -> np.ndarray:
    """Apply the exact sum-zero expansion Jacobian transpose to a gradient."""
    full_grad = np.asarray(full_grad, dtype=float).reshape(-1)
    if full_grad.size <= 1:
        return np.array([], dtype=float)
    return full_grad[:-1] - full_grad[-1]


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

    requested_centered = bool(centered)
    has_hard_anchor = bool(np.isfinite(anchors_vec).any())
    has_group_anchor = any(
        value not in (None, "") and not bool(pd.isna(value))
        for value in groups_vec
    )
    # A hard element anchor or an absolute group-mean constraint already fixes
    # this facet's origin.  Applying sum-to-zero again to only the remaining
    # unanchored levels removes one estimable coordinate and produces a
    # different, over-constrained likelihood (not FACETS element-anchor
    # semantics).  Sum-to-zero remains the identification rule only when no
    # absolute origin constraint is present.
    absolute_origin_identified = has_hard_anchor or has_group_anchor
    effective_centered = requested_centered and not absolute_origin_identified

    spec = {
        "levels": lvl,
        "anchors": anchors_vec,
        "groups": groups_vec,
        "group_values": group_values_map,
        "centered": effective_centered,
        "requested_centered": requested_centered,
        "absolute_origin_identified": absolute_origin_identified,
        "origin_constraint": (
            "absolute_anchor"
            if has_hard_anchor else
            "absolute_group_anchor"
            if has_group_anchor else
            "sum_to_zero"
            if requested_centered else
            "free"
        ),
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
    group_ids = sorted({
        groups[i]
        for i in free_idx
        if groups[i] not in (None, "", np.nan)
    })
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


def collapse_facet_gradient(full_grad, spec):
    """Map full level gradients back to the constrained free-parameter vector."""
    full_grad = np.asarray(full_grad, dtype=float)
    anchors = spec["anchors"]
    groups = spec["groups"]
    free_idx = np.where(np.isnan(anchors))[0]
    if free_idx.size == 0:
        return np.array([], dtype=float)

    pieces = []
    group_ids = sorted({
        groups[i]
        for i in free_idx
        if groups[i] not in (None, "", np.nan)
    })
    grouped_free = set()
    for gid in group_ids:
        group_levels = [i for i, g in enumerate(groups) if g == gid]
        free_in_group = [i for i in group_levels if np.isnan(anchors[i])]
        grouped_free.update(free_in_group)
        if len(free_in_group) > 1:
            last = free_in_group[-1]
            pieces.append(full_grad[free_in_group[:-1]] - full_grad[last])

    ungrouped_idx = [
        i for i in free_idx
        if i not in grouped_free and groups[i] in (None, "", np.nan)
    ]
    m = len(ungrouped_idx)
    if m > 0:
        if spec["centered"]:
            if m > 1:
                last = ungrouped_idx[-1]
                pieces.append(full_grad[ungrouped_idx[:-1]] - full_grad[last])
        else:
            pieces.append(full_grad[ungrouped_idx])

    if not pieces:
        return np.array([], dtype=float)
    return np.concatenate([np.asarray(p, dtype=float) for p in pieces])


def parse_population_formula(population_formula: str | None) -> dict:
    """Parse a small, explicit population formula for latent regression."""
    text = str(population_formula or "").strip()
    if not text:
        return {"enabled": False, "intercept": True, "terms": [], "formula": ""}
    rhs = text.split("~", 1)[1].strip() if "~" in text else text
    if not rhs:
        raise ValueError("population_formula must have a right-hand side, e.g. '~ grade + ses'.")
    if any(op in rhs for op in ("*", ":", "^", "/", "(", ")")):
        raise ValueError(
            "population_formula currently supports main effects only, e.g. '~ grade + ses'. "
            "Interactions and transformations are not yet enabled."
        )
    intercept = True
    terms: list[str] = []
    for raw in rhs.split("+"):
        term = raw.strip()
        if not term or term == "1":
            continue
        if term in {"0", "-1"}:
            intercept = False
            continue
        if term.startswith("`") and term.endswith("`") and len(term) >= 2:
            term = term[1:-1].strip()
        if not term:
            continue
        terms.append(term)
    if not intercept and not terms:
        raise ValueError("population_formula cannot remove the intercept without adding any covariate terms.")
    return {"enabled": True, "intercept": intercept, "terms": terms, "formula": text}


def parse_population_covariate_terms(raw_terms) -> list[str]:
    """Parse comma/newline-separated population covariate override terms."""
    if raw_terms is None:
        return []
    if isinstance(raw_terms, (list, tuple, set)):
        items = raw_terms
    else:
        items = re.split(r"[,;\n]+", str(raw_terms))
    out: list[str] = []
    seen: set[str] = set()
    for item in items:
        term = str(item).strip()
        if term.startswith("`") and term.endswith("`") and len(term) >= 2:
            term = term[1:-1].strip()
        if term and term not in seen:
            seen.add(term)
            out.append(term)
    return out


POPULATION_TYPE_SUMMARY_COLUMNS = [
    "Term",
    "InferredType",
    "Override",
    "N",
    "UniqueValues",
    "NumericConvertible",
    "IntegerLike",
    "ExampleValues",
    "ReviewFlag",
    "Recommendation",
]


def infer_population_covariate_type(
    term: str,
    series: pd.Series,
    forced_categorical: bool = False,
    forced_numeric: bool = False,
) -> dict:
    """Infer how a latent-regression covariate will be encoded and flag risky numeric codes."""
    clean = pd.Series(series).dropna()
    text_values = clean.astype(str)
    numeric = pd.to_numeric(clean, errors="coerce")
    n = int(len(clean))
    unique_values = int(text_values.nunique(dropna=True)) if n else 0
    numeric_convertible = bool(n > 0 and numeric.notna().all())
    finite_numeric = numeric.dropna().to_numpy(dtype=float) if numeric_convertible else np.array([], dtype=float)
    integer_like = bool(
        finite_numeric.size > 0
        and np.all(np.isfinite(finite_numeric))
        and np.allclose(finite_numeric, np.round(finite_numeric), rtol=0.0, atol=1e-9)
    )
    if forced_numeric:
        inferred_type = "numeric" if numeric_convertible else "numeric_invalid"
        override = "numeric"
    elif forced_categorical:
        inferred_type = "categorical"
        override = "categorical"
    elif numeric_convertible:
        inferred_type = "numeric"
        override = "none"
    else:
        inferred_type = "categorical"
        override = "none"

    examples = text_values.drop_duplicates().head(8).tolist()
    review_flag = False
    recommendation = "No action needed."
    if inferred_type == "numeric_invalid":
        review_flag = True
        recommendation = (
            f"`{term}` was forced numeric but contains non-numeric values. "
            "Clean the person_data column or remove it from the numeric override."
        )
    elif (
        inferred_type == "numeric"
        and override == "none"
        and integer_like
        and 2 <= unique_values <= 10
    ):
        review_flag = True
        recommendation = (
            f"`{term}` is being treated as numeric, but it has {unique_values} integer-like levels. "
            "If these numbers are labels, IDs, grade groups, or ordered categories rather than a continuous scale, "
            "add the term to Force categorical covariates."
        )
    elif inferred_type == "categorical" and unique_values > max(20, min(50, n // 2 if n else 20)):
        review_flag = True
        recommendation = (
            f"`{term}` is categorical with {unique_values} levels. "
            "High-cardinality dummy coding can be unstable; consider simplifying the grouping or using a numeric scale if appropriate."
        )

    return {
        "Term": str(term),
        "InferredType": inferred_type,
        "Override": override,
        "N": n,
        "UniqueValues": unique_values,
        "NumericConvertible": numeric_convertible,
        "IntegerLike": integer_like,
        "ExampleValues": ", ".join(examples),
        "ReviewFlag": bool(review_flag),
        "Recommendation": recommendation,
    }


def build_population_model(
    prep,
    person_data=None,
    person_id_col=None,
    population_formula=None,
    standardize_numeric=False,
    categorical_terms=None,
    numeric_terms=None,
) -> dict:
    """Build the person-level latent regression design matrix for MML."""
    parsed = parse_population_formula(population_formula)
    standardize_numeric = bool(standardize_numeric)
    categorical_set = set(parse_population_covariate_terms(categorical_terms))
    numeric_set = set(parse_population_covariate_terms(numeric_terms))
    overlap = sorted(categorical_set & numeric_set)
    if overlap:
        raise ValueError(
            "Population covariate terms cannot be forced to both categorical and numeric: "
            + ", ".join(overlap)
        )
    if not parsed["enabled"]:
        return {
            "enabled": False,
            "formula": "",
            "intercept": True,
            "terms": [],
            "columns": [],
            "X": np.zeros((len(prep["levels"]["Person"]), 0), dtype=float),
            "n_params": 0,
            "person_data": pd.DataFrame({"Person": prep["levels"]["Person"]}),
            "standardize_numeric": False,
            "numeric_transforms": {},
            "categorical_terms": [],
            "numeric_terms": [],
            "term_types": {},
            "type_summary": pd.DataFrame(columns=POPULATION_TYPE_SUMMARY_COLUMNS),
            "warnings": [],
        }
    unknown_type_terms = sorted((categorical_set | numeric_set) - set(parsed["terms"]))
    if unknown_type_terms:
        raise ValueError(
            "Population covariate type override term(s) are not in population_formula: "
            + ", ".join(unknown_type_terms)
        )

    person_levels = [str(x) for x in prep["levels"]["Person"]]
    need_person_data = bool(parsed["terms"])
    if person_data is None or not isinstance(person_data, pd.DataFrame) or person_data.empty:
        if need_person_data:
            raise ValueError("person_data is required when population_formula includes covariates.")
        joined = pd.DataFrame({"Person": person_levels})
    else:
        pdf = person_data.copy()
        if person_id_col is None:
            person_id_col = "Person" if "Person" in pdf.columns else pdf.columns[0]
        if person_id_col not in pdf.columns:
            raise ValueError(f"person_id_col '{person_id_col}' was not found in person_data.")
        pdf = pdf.rename(columns={person_id_col: "Person"})
        pdf["Person"] = pdf["Person"].astype(str)
        if pdf["Person"].duplicated().any():
            examples = pdf.loc[pdf["Person"].duplicated(), "Person"].head(5).tolist()
            raise ValueError("person_data contains duplicate person IDs: " + ", ".join(map(str, examples)))
        joined = pd.DataFrame({"Person": person_levels}).merge(pdf, on="Person", how="left", indicator=True)
        missing = joined.loc[joined["_merge"] != "both", "Person"].astype(str).head(10).tolist()
        if missing and need_person_data:
            raise ValueError(
                "person_data is missing fitted persons required by population_formula: "
                + ", ".join(missing)
            )
        if "_merge" in joined.columns:
            joined = joined.drop(columns=["_merge"])

    columns: list[str] = []
    matrix_parts: list[np.ndarray] = []
    warnings: list[str] = []
    numeric_transforms: dict[str, dict[str, float | bool]] = {}
    term_types: dict[str, str] = {}
    type_summary_rows: list[dict] = []
    if parsed["intercept"]:
        columns.append("Intercept")
        matrix_parts.append(np.ones((len(person_levels), 1), dtype=float))

    for term in parsed["terms"]:
        if term not in joined.columns:
            raise ValueError(f"population_formula term '{term}' was not found in person_data.")
        series = joined[term]
        if series.isna().any():
            missing_people = joined.loc[series.isna(), "Person"].astype(str).head(10).tolist()
            raise ValueError(
                f"population covariate '{term}' has missing values for fitted persons: "
                + ", ".join(missing_people)
            )
        numeric = pd.to_numeric(series, errors="coerce")
        forced_numeric = term in numeric_set
        forced_categorical = term in categorical_set
        type_info = infer_population_covariate_type(
            term,
            series,
            forced_categorical=forced_categorical,
            forced_numeric=forced_numeric,
        )
        type_summary_rows.append(type_info)
        if type_info.get("ReviewFlag"):
            warnings.append(str(type_info.get("Recommendation", "")))
        if forced_numeric or (not forced_categorical and type_info.get("InferredType") == "numeric"):
            if numeric.isna().any():
                raise ValueError(
                    f"population covariate '{term}' was forced numeric but contains non-numeric value(s)."
                )
            values = numeric.to_numpy(dtype=float)
            center = 0.0
            scale = 1.0
            if standardize_numeric:
                center = float(np.mean(values)) if values.size else 0.0
                scale = float(np.std(values, ddof=0)) if values.size else 1.0
                if not np.isfinite(scale) or scale <= 1e-12:
                    warnings.append(
                        f"population covariate '{term}' has near-zero SD; it was centered without rescaling."
                    )
                    scale = 1.0
                values = (values - center) / scale
            numeric_transforms[term] = {
                "center": float(center),
                "scale": float(scale),
                "standardized": bool(standardize_numeric),
            }
            columns.append(term)
            matrix_parts.append(values.reshape(-1, 1))
            term_types[term] = "numeric"
        else:
            cat = pd.Categorical(series.astype(str), categories=sorted(series.astype(str).unique()))
            dummies = pd.get_dummies(cat, prefix=term, drop_first=True, dtype=float)
            if dummies.shape[1] == 0:
                warnings.append(f"population covariate '{term}' has one level and was omitted.")
                term_types[term] = "categorical_omitted_one_level"
                continue
            columns.extend(list(dummies.columns))
            matrix_parts.append(dummies.to_numpy(dtype=float))
            term_types[term] = "categorical"

    X = np.column_stack(matrix_parts) if matrix_parts else np.zeros((len(person_levels), 0), dtype=float)
    if X.size and not np.all(np.isfinite(X)):
        raise ValueError("population model matrix contains non-finite values.")
    type_summary = pd.DataFrame(type_summary_rows, columns=POPULATION_TYPE_SUMMARY_COLUMNS)
    return {
        "enabled": True,
        "formula": parsed["formula"],
        "intercept": bool(parsed["intercept"]),
        "terms": parsed["terms"],
        "columns": columns,
        "X": X,
        "n_params": int(X.shape[1]),
        "person_data": joined,
        "standardize_numeric": bool(standardize_numeric),
        "numeric_transforms": numeric_transforms,
        "categorical_terms": sorted(categorical_set),
        "numeric_terms": sorted(numeric_set),
        "term_types": term_types,
        "type_summary": type_summary,
        "warnings": warnings,
    }


def compute_population_mu(params, config):
    pop = config.get("population_model", {})
    if not pop or not pop.get("enabled"):
        return np.zeros(int(config.get("n_person", 0)), dtype=float)
    X = np.asarray(pop.get("X", np.zeros((int(config.get("n_person", 0)), 0))), dtype=float)
    beta = np.asarray(params.get("population", np.array([], dtype=float)), dtype=float)
    if X.size == 0 or beta.size == 0:
        return np.zeros(int(config.get("n_person", X.shape[0] if X.ndim == 2 else 0)), dtype=float)
    return X @ beta


def build_param_sizes(config):
    n_steps = max(config["n_cat"] - 1, 0)
    n_step_free = max(n_steps - 1, 0)
    sizes = OrderedDict()
    sizes["theta"] = config["theta_spec"]["n_params"] if config["method"] == "JMLE" else 0
    pop = config.get("population_model", {})
    sizes["population"] = int(pop.get("n_params", 0)) if config["method"] == "MML" and pop.get("enabled") else 0
    for facet in config["facet_names"]:
        sizes[facet] = config["facet_specs"][facet]["n_params"]
    if config["model"] == "RSM":
        sizes["steps"] = n_step_free
    elif config["model"] in {"PCM", "GPCM"}:
        if not config.get("step_facet") or config["step_facet"] not in config["facet_names"]:
            raise ValueError(f"{config['model']} requires a valid step facet.")
        sizes["steps"] = len(config["facet_levels"][config["step_facet"]]) * n_step_free
        if config["model"] == "GPCM":
            if not config.get("slope_facet") or config["slope_facet"] not in config["facet_names"]:
                raise ValueError("GPCM requires a valid slope facet.")
            if config["slope_facet"] != config["step_facet"]:
                raise ValueError("This bounded GPCM implementation requires slope_facet == step_facet.")
            sizes["log_slopes"] = max(
                len(config["facet_levels"][config["slope_facet"]]) - 1,
                0,
            )
    else:
        raise ValueError("Model must be one of: RSM, PCM, GPCM.")
    return sizes


def _constraint_expansion_matrix_for_identifiability(spec):
    """Exact affine-expansion Jacobian for one constrained facet block."""
    n_params = int(spec.get("n_params", 0))
    n_levels = len(spec.get("levels", []))
    if n_params == 0:
        return np.zeros((n_levels, 0), dtype=float)
    origin = expand_facet_with_constraints(np.zeros(n_params, dtype=float), spec)
    matrix = np.zeros((n_levels, n_params), dtype=float)
    for coordinate in range(n_params):
        basis = np.zeros(n_params, dtype=float)
        basis[coordinate] = 1.0
        matrix[:, coordinate] = expand_facet_with_constraints(basis, spec) - origin
    return matrix


def _constraint_free_labels_for_identifiability(spec):
    """Label exact free coordinates and their constraint-derived counterparts."""
    anchors = np.asarray(spec.get("anchors", []), dtype=float)
    groups = list(spec.get("groups", []))
    levels = [str(value) for value in spec.get("levels", [])]
    free_idx = np.where(np.isnan(anchors))[0].tolist()
    rows = []
    grouped_free = set()

    def has_group(value):
        if value is None or value == "":
            return False
        try:
            return not bool(pd.isna(value))
        except (TypeError, ValueError):
            return True

    group_ids = sorted({
        groups[index]
        for index in free_idx
        if has_group(groups[index])
    })
    for group_id in group_ids:
        group_levels = [index for index, value in enumerate(groups) if value == group_id]
        free_in_group = [index for index in group_levels if np.isnan(anchors[index])]
        grouped_free.update(free_in_group)
        if len(free_in_group) > 1:
            derived = levels[free_in_group[-1]]
            for index in free_in_group[:-1]:
                rows.append({
                    "Level": levels[index],
                    "DerivedCounterpart": derived,
                    "CoordinateMeaning": f"group {group_id}: final level derived",
                })
    ungrouped = [
        index for index in free_idx
        if index not in grouped_free and not has_group(groups[index])
    ]
    if bool(spec.get("centered", False)):
        if len(ungrouped) > 1:
            derived = levels[ungrouped[-1]]
            for index in ungrouped[:-1]:
                rows.append({
                    "Level": levels[index],
                    "DerivedCounterpart": derived,
                    "CoordinateMeaning": "final unanchored level derived by centering constraint",
                })
    else:
        for index in ungrouped:
            rows.append({
                "Level": levels[index],
                "DerivedCounterpart": "",
                "CoordinateMeaning": "direct unconstrained level",
            })
    return rows


def _person_facet_connectivity_summary(prep):
    """Connected-component summaries for each observed Person x facet graph."""
    data = prep.get("data", pd.DataFrame())
    rows = []
    if not isinstance(data, pd.DataFrame) or data.empty:
        return pd.DataFrame()
    for facet in prep.get("facet_names", []):
        adjacency = {}
        pairs = data[["Person", facet]].drop_duplicates()
        for pair in pairs.itertuples(index=False, name=None):
            person_node = f"Person::{pair[0]}"
            facet_node = f"{facet}::{pair[1]}"
            adjacency.setdefault(person_node, set()).add(facet_node)
            adjacency.setdefault(facet_node, set()).add(person_node)
        seen = set()
        component_sizes = []
        for node in sorted(adjacency):
            if node in seen:
                continue
            stack = [node]
            component = set()
            while stack:
                current = stack.pop()
                if current in component:
                    continue
                component.add(current)
                seen.add(current)
                stack.extend(adjacency[current] - component)
            component_sizes.append(len(component))
        per_person = data.groupby("Person", observed=True)[facet].nunique()
        per_level = data.groupby(facet, observed=True)["Person"].nunique()
        rows.append({
            "Facet": str(facet),
            "Components": int(len(component_sizes)),
            "MinimumRatersOrLevelsPerPerson": int(per_person.min()) if len(per_person) else 0,
            "MedianRatersOrLevelsPerPerson": float(per_person.median()) if len(per_person) else np.nan,
            "MaximumRatersOrLevelsPerPerson": int(per_person.max()) if len(per_person) else 0,
            "MinimumPersonsPerFacetLevel": int(per_level.min()) if len(per_level) else 0,
            "MaximumPersonsPerFacetLevel": int(per_level.max()) if len(per_level) else 0,
            "ComponentNodeSizes": "|".join(map(str, sorted(component_sizes, reverse=True))),
        })
    return pd.DataFrame(rows)


def build_jmle_eta_identifiability_audit(
    prep,
    config,
    sizes=None,
    idx=None,
    max_dense_cells=5_000_000,
):
    """Audit exact structural rank of JMLE person/facet eta coordinates.

    This audit covers the additive coordinates entering the linear predictor
    ``eta`` after anchors and centering constraints are expanded. RSM/PCM/GPCM
    step and GPCM slope identification are intentionally outside this scope.
    """
    scope = "eta_person_facet_blocks_only"
    empty_summary = pd.DataFrame([{
        "Available": False,
        "Scope": scope,
        "Status": "not_applicable_mml",
        "EtaDesignRows": 0,
        "EtaDesignColumns": 0,
        "EtaDesignRank": np.nan,
        "EtaStructuralNullity": np.nan,
        "EtaStructurallyIdentified": pd.NA,
        "InferenceReady": pd.NA,
        "Action": "Use the estimator-specific MML qualification path.",
    }])
    if str(config.get("method", "")).upper() != "JMLE":
        return {
            "available": False,
            "scope": scope,
            "status": "not_applicable_mml",
            "eta_structurally_identified": None,
            "inference_ready": None,
            "summary": empty_summary,
            "coordinate_null_weight": pd.DataFrame(),
            "null_space_block_energy": pd.DataFrame(),
            "connectivity": _person_facet_connectivity_summary(prep),
        }
    sizes = build_param_sizes(config) if sizes is None else sizes
    idx = build_indices(
        prep,
        step_facet=config.get("step_facet"),
        slope_facet=config.get("slope_facet"),
    ) if idx is None else idx
    pieces = []
    metadata = []
    global_coordinate = 0

    def add_block(block, spec, observation_index, sign=1.0):
        nonlocal global_coordinate
        matrix = _constraint_expansion_matrix_for_identifiability(spec)
        if matrix.shape[1] == 0:
            return
        pieces.append(float(sign) * matrix[np.asarray(observation_index, dtype=int)])
        labels = _constraint_free_labels_for_identifiability(spec)
        if len(labels) != matrix.shape[1]:
            raise RuntimeError(
                f"Identifiability labels for {block} have {len(labels)} rows; "
                f"expected {matrix.shape[1]}."
            )
        for local_coordinate, label in enumerate(labels):
            metadata.append({
                "GlobalCoordinate": global_coordinate,
                "Block": str(block),
                "BlockCoordinate": int(local_coordinate),
                **label,
            })
            global_coordinate += 1

    n_rows = len(idx.get("person", []))
    specs = [config["theta_spec"], *(config["facet_specs"][facet] for facet in config["facet_names"])]
    n_columns = sum(spec["n_params"] for spec in specs)
    largest_expansion = max((len(spec["levels"]) * spec["n_params"] for spec in specs), default=0)
    connectivity = _person_facet_connectivity_summary(prep)
    if max(int(n_rows) * int(n_columns), largest_expansion) > int(max_dense_cells):
        summary = pd.DataFrame([{
            "Available": False,
            "Scope": scope,
            "Status": "audit_size_limit",
            "EtaDesignRows": int(n_rows),
            "EtaDesignColumns": int(n_columns),
            "EtaDesignRank": np.nan,
            "EtaStructuralNullity": np.nan,
            "EtaStructurallyIdentified": False,
            "InferenceReady": False,
            "Action": (
                "The exact dense rank audit exceeded its technical size limit. "
                "Do not interpret optimizer convergence as identification; use a "
                "validated large-design audit before inference."
            ),
        }])
        return {
            "available": False,
            "scope": scope,
            "status": "audit_size_limit",
            "eta_structurally_identified": False,
            "inference_ready": False,
            "summary": summary,
            "coordinate_null_weight": pd.DataFrame(metadata),
            "null_space_block_energy": pd.DataFrame(),
            "connectivity": connectivity,
        }
    add_block("theta", config["theta_spec"], idx["person"], 1.0)
    facet_signs = config.get("facet_signs", {})
    for facet in config.get("facet_names", []):
        add_block(
            facet,
            config["facet_specs"][facet],
            idx["facets"][facet],
            facet_signs.get(facet, -1),
        )
    design = np.column_stack(pieces) if pieces else np.zeros((n_rows, 0), dtype=float)
    # Economy SVD avoids allocating an n_columns x n_columns Vh for a wide
    # design that still falls below the dense-cell guard. Rank remains exact
    # under the registered tolerance. A complete right-null basis is available
    # from the economy form only when rows >= columns.
    _, singular_values, vh = np.linalg.svd(design, full_matrices=False)
    largest = float(singular_values[0]) if len(singular_values) else 0.0
    tolerance = max(design.shape, default=0) * np.finfo(float).eps * largest
    rank = int(np.sum(singular_values > tolerance))
    nullity = int(design.shape[1] - rank)
    identified = nullity == 0
    null_basis_complete = bool(design.shape[0] >= design.shape[1])
    null_vectors = (
        vh[rank:, :]
        if nullity > 0 and null_basis_complete
        else np.zeros((0, design.shape[1]))
    )
    null_weight = (
        np.sum(null_vectors ** 2, axis=0)
        if nullity > 0 and null_basis_complete
        else (
            np.zeros(design.shape[1])
            if nullity == 0 else np.full(design.shape[1], np.nan)
        )
    )
    coordinate_null = pd.DataFrame(metadata)
    coordinate_null["NullProjectionWeight"] = null_weight
    if nullity > 0 and not null_basis_complete:
        coordinate_null["NullVulnerable"] = pd.Series(
            pd.array([pd.NA] * len(coordinate_null), dtype="boolean"),
            index=coordinate_null.index,
        )
    else:
        coordinate_null["NullVulnerable"] = coordinate_null["NullProjectionWeight"].gt(1e-10)
    energy_rows = []
    for vector_index, vector in enumerate(null_vectors):
        denominator = float(np.sum(vector ** 2))
        for block, block_frame in coordinate_null.groupby("Block", sort=False):
            coordinates = block_frame["GlobalCoordinate"].to_numpy(dtype=int)
            energy_rows.append({
                "NullVector": int(vector_index + 1),
                "Block": str(block),
                "SquaredEnergyShare": (
                    float(np.sum(vector[coordinates] ** 2) / denominator)
                    if denominator > 0 else np.nan
                ),
                "MaximumAbsLoading": float(np.max(np.abs(vector[coordinates]))),
            })
    status = "eta_identified" if identified else "eta_rank_deficient"
    action = (
        "No exact person/facet eta null direction was detected; retain the separate numerical and model-fit checks."
        if identified else
        "Do not use this JMLE fit for inference or bias decisions. Add cross-facet connectivity or anchors, redesign the model, or use an explicitly qualified alternative estimator."
    )
    summary = pd.DataFrame([{
        "Available": True,
        "Scope": scope,
        "Status": status,
        "EtaDesignRows": int(design.shape[0]),
        "EtaDesignColumns": int(design.shape[1]),
        "EtaDesignRank": rank,
        "EtaStructuralNullity": nullity,
        "EtaStructurallyIdentified": bool(identified),
        "InferenceReady": bool(identified),
        "NullSpaceBasisComplete": null_basis_complete,
        "RankTolerance": tolerance,
        "LargestSingularValue": largest,
        "SmallestSingularValue": (
            float(singular_values[-1]) if len(singular_values) else np.nan
        ),
        "Action": action,
    }])
    return {
        "available": True,
        "scope": scope,
        "status": status,
        "eta_structurally_identified": bool(identified),
        "inference_ready": bool(identified),
        "summary": summary,
        "coordinate_null_weight": coordinate_null,
        "null_space_block_energy": pd.DataFrame(energy_rows),
        "connectivity": connectivity,
    }


def build_jmle_person_boundary_audit(
    prep,
    config,
    estimates=None,
    *,
    fit_inference_ready=False,
):
    """Classify JMLE Person score boundaries without estimate cutoffs.

    An all-minimum or all-maximum response pattern has no finite unbounded
    JMLE Person MLE. The fitted optimizer/constraint value is retained for
    reproducibility, while ``ReportableEstimate`` is withheld. Classification
    uses retained integer scores only, never the magnitude or display rounding
    of the terminal estimate.
    """
    scope = "jmle_person_score_boundary_only"
    method = str(config.get("method", "")).upper()
    if method != "JMLE":
        summary = pd.DataFrame([{
            "Available": False,
            "Scope": scope,
            "Status": "not_applicable_mml",
            "Persons": int(config.get("n_person", 0)),
            "ExtremePersons": 0,
            "AllMinimumPersons": 0,
            "AllMaximumPersons": 0,
            "FiniteJMLEPersons": pd.NA,
            "PersonMeasureInferenceReady": pd.NA,
            "Action": "Use the estimator-specific MML/EAP person-score qualification path.",
        }])
        return {
            "available": False,
            "scope": scope,
            "status": "not_applicable_mml",
            "person_measure_inference_ready": None,
            "summary": summary,
            "persons": pd.DataFrame(),
        }

    data = prep.get("data", pd.DataFrame())
    levels = [str(value) for value in prep.get("levels", {}).get("Person", [])]
    rating_min = int(prep.get("rating_min", 0))
    rating_max = int(prep.get("rating_max", max(int(config.get("n_cat", 1)) - 1, 0)))
    estimate_values = np.asarray(
        estimates if estimates is not None else np.full(len(levels), np.nan),
        dtype=float,
    ).reshape(-1)
    if len(estimate_values) != len(levels):
        raise ValueError(
            "JMLE Person boundary audit estimate count does not match Person levels."
        )
    estimate_map = dict(zip(levels, estimate_values, strict=False))
    score_source = pd.DataFrame({
        "Person": data.get("Person", pd.Series(dtype=str)).astype(str),
        "Score": pd.to_numeric(data.get("Score"), errors="coerce"),
        "Weight": pd.to_numeric(
            data.get("Weight", pd.Series(1.0, index=data.index)),
            errors="coerce",
        ),
    })
    score_source = score_source.loc[
        score_source["Score"].notna()
        & score_source["Weight"].notna()
        & score_source["Weight"].gt(0)
    ].copy()
    score_source["WeightedScore"] = score_source["Score"] * score_source["Weight"]
    score_stats = (
        score_source.groupby("Person", observed=True)
        .agg(
            ObservedCount=("Score", "size"),
            EffectiveWeight=("Weight", "sum"),
            ObservedScoreMin=("Score", "min"),
            ObservedScoreMax=("Score", "max"),
            RawScoreTotal=("Score", "sum"),
            WeightedScoreTotal=("WeightedScore", "sum"),
        )
        .reindex(levels)
    )
    rows = []
    for person in levels:
        score_row = score_stats.loc[person]
        observed_count = int(score_row["ObservedCount"]) if pd.notna(score_row["ObservedCount"]) else 0
        effective_weight = float(score_row["EffectiveWeight"]) if observed_count else 0.0
        observed_min = float(score_row["ObservedScoreMin"]) if observed_count else np.nan
        observed_max = float(score_row["ObservedScoreMax"]) if observed_count else np.nan
        all_minimum = bool(observed_count > 0 and observed_min == rating_min and observed_max == rating_min)
        all_maximum = bool(observed_count > 0 and observed_min == rating_max and observed_max == rating_max)
        extreme = bool(all_minimum or all_maximum)
        direction = "all_minimum" if all_minimum else "all_maximum" if all_maximum else "interior"
        estimate = float(estimate_map.get(person, np.nan))
        finite_jmle = bool(not extreme and np.isfinite(estimate))
        person_ready = bool(finite_jmle and fit_inference_ready)
        rows.append({
            "Person": person,
            "ObservedCount": observed_count,
            "EffectiveWeight": effective_weight,
            "ObservedScoreMin": observed_min,
            "ObservedScoreMax": observed_max,
            "RawScoreTotal": float(score_row["RawScoreTotal"]) if observed_count else np.nan,
            "WeightedScoreTotal": float(score_row["WeightedScoreTotal"]) if observed_count else np.nan,
            "MinimumPossibleWeightedTotal": rating_min * effective_weight,
            "MaximumPossibleWeightedTotal": rating_max * effective_weight,
            "ExtremeScorePattern": extreme,
            "ExtremeScoreDirection": direction,
            "FiniteJMLEEstimate": finite_jmle,
            "PersonInferenceReady": person_ready,
            "Estimate": estimate,
            "ReportableEstimate": estimate if person_ready else np.nan,
            "EstimateRole": (
                "optimizer_or_constraint_value_not_finite_person_mle"
                if extreme else
                "finite_jmle_estimate" if person_ready else
                "finite_boundary_status_but_fit_not_inference_ready"
            ),
        })
    persons = pd.DataFrame(rows)
    extreme_n = int(persons["ExtremeScorePattern"].sum()) if not persons.empty else 0
    all_min_n = int(persons["ExtremeScoreDirection"].eq("all_minimum").sum()) if not persons.empty else 0
    all_max_n = int(persons["ExtremeScoreDirection"].eq("all_maximum").sum()) if not persons.empty else 0
    finite_n = int(persons["FiniteJMLEEstimate"].sum()) if not persons.empty else 0
    person_ready = bool(not persons.empty and persons["PersonInferenceReady"].all())
    status = "extreme_persons_present" if extreme_n else "no_extreme_persons"
    summary = pd.DataFrame([{
        "Available": True,
        "Scope": scope,
        "Status": status,
        "Persons": int(len(persons)),
        "ExtremePersons": extreme_n,
        "AllMinimumPersons": all_min_n,
        "AllMaximumPersons": all_max_n,
        "FiniteJMLEPersons": finite_n,
        "PersonMeasureInferenceReady": person_ready,
        "Action": (
            "Do not report finite JMLE Person measures for boundary-score Persons; "
            "use ReportableEstimate and retain Estimate only as a technical terminal value."
            if extreme_n else
            "No score-boundary Person was detected; retain structural, numerical, and model-fit checks."
        ),
    }])
    return {
        "available": True,
        "scope": scope,
        "status": status,
        "person_measure_inference_ready": person_ready,
        "summary": summary,
        "persons": persons,
    }


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
        steps = expand_sum_zero_free(parts["steps"], max(config["n_cat"] - 1, 0))
        steps_mat = None
    else:
        n_levels = len(config["facet_levels"][config["step_facet"]])
        n_steps = max(config["n_cat"] - 1, 0)
        n_step_free = max(n_steps - 1, 0)
        if n_levels == 0 or n_steps == 0:
            steps_mat = np.zeros((n_levels, n_steps))
        else:
            free_mat = np.array(parts["steps"], dtype=float).reshape((n_levels, n_step_free))
            steps_mat = np.vstack([
                expand_sum_zero_free(row, n_steps)
                for row in free_mat
            ])
        steps = None

    log_slopes = np.array([], dtype=float)
    slopes = np.array([], dtype=float)
    if config["model"] == "GPCM":
        n_slope_levels = len(config["facet_levels"][config["slope_facet"]])
        log_slopes = expand_sum_zero_free(
            np.asarray(parts.get("log_slopes", []), dtype=float),
            n_slope_levels,
        )
        slopes = np.exp(log_slopes)

    return {
        "theta": theta,
        "population": np.asarray(parts.get("population", []), dtype=float),
        "facets": facets,
        "steps": steps,
        "steps_mat": steps_mat,
        "log_slopes": log_slopes,
        "slopes": slopes,
    }


def build_optimizer_bounds(sizes, config):
    """Box bounds for free optimizer coordinates."""
    bounds = []
    log_slope_bounds = tuple(config.get("gpcm_log_slope_bounds", (-3.0, 3.0)))
    for name, k in sizes.items():
        if name == "log_slopes" and config.get("model") == "GPCM":
            bounds.extend([log_slope_bounds] * int(k))
        else:
            bounds.extend([(None, None)] * int(k))
    return bounds if bounds else None


def build_optimizer_constraints(sizes, config):
    """Linear constraints that keep every expanded GPCM log-slope bounded.

    The first ``L - 1`` log-slopes are free and box-bounded.  The last is
    ``-sum(free)``; this constraint applies the same configured bounds to that
    derived coordinate without shrinking the feasible identified parameter
    space.
    """
    if config.get("model") != "GPCM" or int(sizes.get("log_slopes", 0)) <= 0:
        return []
    lower, upper = tuple(config.get("gpcm_log_slope_bounds", (-3.0, 3.0)))
    if not (np.isfinite(lower) and np.isfinite(upper) and float(lower) < float(upper)):
        raise ValueError("gpcm_log_slope_bounds must be finite and strictly increasing.")
    slices = _build_param_slices(sizes)
    row = np.zeros(int(sum(sizes.values())), dtype=float)
    row[slices["log_slopes"]] = 1.0
    return [LinearConstraint(row, -float(upper), -float(lower))]


def optimizer_method_and_options(sizes, config, maxiter, gtol, ftol):
    """Select an optimizer compatible with the identified bounded space."""
    # The exact n-1 coordinates remove flat directions, so a loose relative
    # objective stop can otherwise fire while the transformed gradient still
    # has visible movement.  Keep the user's gradient tolerance but cap the
    # objective tolerance at a numerically modest 1e-9.
    effective_ftol = min(float(ftol), 1e-9)
    if build_optimizer_constraints(sizes, config):
        return "SLSQP", {"maxiter": int(maxiter), "ftol": effective_ftol}
    return "L-BFGS-B", {
        "maxiter": int(maxiter),
        "gtol": float(gtol),
        "ftol": effective_ftol,
    }


def build_parameterization_audit(config, sizes, params) -> pd.DataFrame:
    """Describe identified optimizer blocks and verify expanded constraints."""
    model = str(config.get("model", ""))
    n_steps = max(int(config.get("n_cat", 0)) - 1, 0)
    n_step_levels = (
        1
        if model == "RSM"
        else len(config.get("facet_levels", {}).get(config.get("step_facet"), []))
    )
    if model == "RSM":
        step_matrix = np.asarray(params.get("steps", []), dtype=float).reshape(1, -1)
    else:
        step_matrix = np.asarray(
            params.get("steps_mat", np.zeros((n_step_levels, n_steps))),
            dtype=float,
        ).reshape(n_step_levels, n_steps)
    step_constraint = (
        float(np.max(np.abs(np.sum(step_matrix, axis=1))))
        if step_matrix.size
        else 0.0
    )
    step_expanded = int(n_step_levels * n_steps)
    step_free = int(sizes.get("steps", 0))
    rows = [{
        "SchemaVersion": "mfrm_exact_identified_v1",
        "Block": "step_thresholds",
        "ExpandedCoordinates": step_expanded,
        "FreeCoordinates": step_free,
        "ExpansionRank": step_free,
        "ArtificialCoordinateNullDirections": 0,
        "Constraint": "sum(step) = 0 within each step-facet level",
        "ConstraintMaxAbsResidual": step_constraint,
        "ExpandedBoundsSatisfied": True,
        "Status": "PASS" if step_constraint <= 1e-10 else "FAIL",
        "RoadmapIssueID": "STAT-001",
    }]

    if model == "GPCM":
        log_slopes = np.asarray(params.get("log_slopes", []), dtype=float)
        lower, upper = tuple(config.get("gpcm_log_slope_bounds", (-3.0, 3.0)))
        slope_constraint = abs(float(np.sum(log_slopes))) if log_slopes.size else 0.0
        bounds_ok = bool(
            np.all(log_slopes >= float(lower) - 1e-10)
            and np.all(log_slopes <= float(upper) + 1e-10)
        )
        slope_free = int(sizes.get("log_slopes", 0))
        rows.append({
            "SchemaVersion": "mfrm_exact_identified_v1",
            "Block": "gpcm_log_slopes",
            "ExpandedCoordinates": int(log_slopes.size),
            "FreeCoordinates": slope_free,
            "ExpansionRank": slope_free,
            "ArtificialCoordinateNullDirections": 0,
            "Constraint": "sum(log_slope) = 0; every expanded log_slope within configured bounds",
            "ConstraintMaxAbsResidual": slope_constraint,
            "ExpandedBoundsSatisfied": bounds_ok,
            "Status": "PASS" if slope_constraint <= 1e-10 and bounds_ok else "FAIL",
            "RoadmapIssueID": "STAT-001; STAT-003",
        })
    return pd.DataFrame(rows)


def _coerce_bool(value, default: bool = True) -> bool:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return default
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"", "na", "nan", "none"}:
            return default
        if text in {"true", "t", "yes", "y", "1", "enabled", "on"}:
            return True
        if text in {"false", "f", "no", "n", "0", "disabled", "off"}:
            return False
    return bool(value)


def normalize_facet_regularization_rows(facet_regularization) -> pd.DataFrame:
    """Normalize optional facet-regularization input to a row table."""
    columns = ["Facet", "Level", "ParameterClass", "Mean", "SD", "Enabled", "Source"]
    if facet_regularization is None or facet_regularization is False:
        return pd.DataFrame(columns=columns)
    if isinstance(facet_regularization, pd.DataFrame):
        raw = facet_regularization.copy()
    elif isinstance(facet_regularization, dict):
        if not _coerce_bool(facet_regularization.get("enabled", True), default=True):
            return pd.DataFrame(columns=columns)
        if "specs" in facet_regularization:
            raw = pd.DataFrame(facet_regularization.get("specs") or [])
        elif "rows" in facet_regularization:
            raw = pd.DataFrame(facet_regularization.get("rows") or [])
        else:
            facets = facet_regularization.get("facets") or facet_regularization.get("facet") or []
            if isinstance(facets, str):
                facets = [facets]
            sd = facet_regularization.get("sd", facet_regularization.get("SD", facet_regularization.get("prior_sd", np.nan)))
            mean = facet_regularization.get("mean", facet_regularization.get("Mean", 0.0))
            raw = pd.DataFrame([
                {
                    "Facet": facet,
                    "Level": "*",
                    "ParameterClass": "facet",
                    "Mean": mean,
                    "SD": sd,
                    "Enabled": True,
                    "Source": facet_regularization.get("source", "dict"),
                }
                for facet in facets
            ])
    elif isinstance(facet_regularization, (list, tuple)):
        raw = pd.DataFrame(list(facet_regularization))
    else:
        raise ValueError("facet_regularization must be None, a dict, a DataFrame, or a list of row dicts.")
    if raw.empty:
        return pd.DataFrame(columns=columns)
    out = raw.copy()
    rename = {}
    for col in out.columns:
        lower = str(col).strip().lower()
        if lower in {"facet", "facetname"}:
            rename[col] = "Facet"
        elif lower in {"level", "element"}:
            rename[col] = "Level"
        elif lower in {"parameterclass", "parameter_class", "class"}:
            rename[col] = "ParameterClass"
        elif lower in {"mean", "prior_mean", "target"}:
            rename[col] = "Mean"
        elif lower in {"sd", "prior_sd", "sigma", "scale"}:
            rename[col] = "SD"
        elif lower in {"enabled", "enable", "active"}:
            rename[col] = "Enabled"
        elif lower in {"source"}:
            rename[col] = "Source"
    out = out.rename(columns=rename)
    for col in columns:
        if col not in out.columns:
            if col == "Level":
                out[col] = "*"
            elif col == "ParameterClass":
                out[col] = "facet"
            elif col == "Mean":
                out[col] = 0.0
            elif col == "Enabled":
                out[col] = True
            elif col == "Source":
                out[col] = "user"
            else:
                out[col] = np.nan
    return out[columns]


def prepare_facet_regularization(facet_regularization, config: dict, sizes: OrderedDict) -> dict:
    """Validate facet regularization and align it to constrained facet levels.

    This implements a Gaussian MAP-style penalty on non-person facet main
    effects only. It does not make facets random effects and it does not alter
    the MML person population distribution.
    """
    rows = normalize_facet_regularization_rows(facet_regularization)
    facet_names = list(config.get("facet_names", []))
    means: dict[str, list[float]] = {}
    sds: dict[str, list[float]] = {}
    active: dict[str, list[bool]] = {}
    for facet in facet_names:
        levels = list(config.get("facet_levels", {}).get(facet, []))
        means[facet] = [0.0] * len(levels)
        sds[facet] = [np.nan] * len(levels)
        active[facet] = [False] * len(levels)

    audit_rows: list[dict] = []
    if rows.empty:
        return {
            "enabled": False,
            "requested": False,
            "mode": "off",
            "scope": "none",
            "spec": pd.DataFrame(columns=["Facet", "Level", "ParameterClass", "Mean", "SD", "Enabled", "Source"]),
            "audit": pd.DataFrame(columns=["Facet", "Level", "Status", "Message", "Mean", "SD", "Source"]),
            "runtime": {
                "enabled": False,
                "means": means,
                "sds": sds,
                "active": active,
                "fingerprint": stable_json_fingerprint({"enabled": False}),
            },
        }

    for i, row in rows.iterrows():
        enabled = _coerce_bool(row.get("Enabled"), default=True)
        facet = str(row.get("Facet", "")).strip()
        level = str(row.get("Level", "*")).strip() or "*"
        parameter_class = str(row.get("ParameterClass", "facet")).strip().lower() or "facet"
        source = str(row.get("Source", "user")).strip() or "user"
        if not enabled or not facet:
            continue
        if parameter_class not in {"facet", "main", "severity", "difficulty"}:
            audit_rows.append({
                "Facet": facet, "Level": level, "Status": "ignored",
                "Message": f"Unsupported ParameterClass '{parameter_class}'. Only facet main effects are currently supported.",
                "Mean": row.get("Mean", np.nan), "SD": row.get("SD", np.nan), "Source": source,
            })
            continue
        if facet == "Person":
            audit_rows.append({
                "Facet": facet, "Level": level, "Status": "ignored",
                "Message": "Person ability distribution is controlled by Person population SD (MML), not facet regularization.",
                "Mean": row.get("Mean", np.nan), "SD": row.get("SD", np.nan), "Source": source,
            })
            continue
        if facet not in facet_names:
            raise ValueError(f"facet_regularization refers to unknown facet '{facet}'.")
        mean = float(row.get("Mean", 0.0) if pd.notna(row.get("Mean", 0.0)) else 0.0)
        sd = float(row.get("SD", np.nan))
        if not np.isfinite(sd) or sd <= 0:
            raise ValueError(f"facet_regularization SD must be positive and finite for facet '{facet}', level '{level}'.")
        levels = list(config.get("facet_levels", {}).get(facet, []))
        spec = config.get("facet_specs", {}).get(facet, {})
        anchors = np.asarray(spec.get("anchors", np.full(len(levels), np.nan)), dtype=float)
        if sizes.get(facet, 0) <= 0:
            audit_rows.append({
                "Facet": facet, "Level": level, "Status": "ignored",
                "Message": "Facet has no free parameters after anchors/dummy constraints.",
                "Mean": mean, "SD": sd, "Source": source,
            })
            continue
        target_levels = levels if level == "*" or level.lower() == "all" else [level]
        for target_level in target_levels:
            if target_level not in levels:
                raise ValueError(f"facet_regularization refers to unknown level '{target_level}' in facet '{facet}'.")
            j = levels.index(target_level)
            if j < anchors.size and np.isfinite(anchors[j]):
                audit_rows.append({
                    "Facet": facet, "Level": target_level, "Status": "ignored",
                    "Message": "Anchored or dummy-constrained level is fixed and is not regularized.",
                    "Mean": mean, "SD": sd, "Source": source,
                })
                continue
            status = "applied"
            message = "Gaussian MAP-style regularization applied to this free facet level."
            if active[facet][j]:
                status = "applied_override"
                message = "Duplicate regularization row; later row overrides earlier mean/SD."
            means[facet][j] = mean
            sds[facet][j] = sd
            active[facet][j] = True
            audit_rows.append({
                "Facet": facet, "Level": target_level, "Status": status,
                "Message": message,
                "Mean": mean, "SD": sd, "Source": source,
            })

    audit = pd.DataFrame(audit_rows, columns=["Facet", "Level", "Status", "Message", "Mean", "SD", "Source"])
    applied = audit[audit["Status"].astype(str).str.startswith("applied")].copy() if not audit.empty else pd.DataFrame()
    runtime_payload = {"enabled": not applied.empty, "means": means, "sds": sds, "active": active}
    runtime_payload["fingerprint"] = stable_json_fingerprint(runtime_payload)
    return {
        "enabled": bool(not applied.empty),
        "requested": bool(len(rows) > 0),
        "mode": "normal_map_penalty" if not applied.empty else "requested_no_free_parameters",
        "scope": "non-person facet main effects",
        "spec": rows.copy(),
        "audit": audit,
        "runtime": runtime_payload,
    }


def facet_regularization_value_grad(par, sizes, config):
    """Return Gaussian facet regularization penalty and gradient."""
    runtime = config.get("facet_regularization", {})
    if not isinstance(runtime, dict) or not runtime.get("enabled"):
        return 0.0, np.zeros(len(par), dtype=float)
    params = expand_params(par, sizes, config)
    means = runtime.get("means", {})
    sds = runtime.get("sds", {})
    active = runtime.get("active", {})
    value = 0.0
    grad_parts = []
    for name, k in sizes.items():
        k = int(k)
        if k <= 0:
            continue
        if name in config.get("facet_names", []):
            full = np.asarray(params["facets"].get(name, np.array([], dtype=float)), dtype=float)
            mean = np.asarray(means.get(name, [0.0] * len(full)), dtype=float)
            sd = np.asarray(sds.get(name, [np.nan] * len(full)), dtype=float)
            is_active = np.asarray(active.get(name, [False] * len(full)), dtype=bool)
            valid = is_active & np.isfinite(sd) & (sd > 0)
            full_grad = np.zeros(len(full), dtype=float)
            if valid.any():
                diff = full[valid] - mean[valid]
                value += 0.5 * float(np.sum((diff / sd[valid]) ** 2))
                full_grad[valid] = diff / (sd[valid] ** 2)
            collapsed = collapse_facet_gradient(full_grad, config["facet_specs"][name])
            if collapsed.size != k:
                raise RuntimeError(f"Facet regularization gradient size mismatch for {name}: expected {k}, got {collapsed.size}.")
            grad_parts.append(collapsed)
        else:
            grad_parts.append(np.zeros(k, dtype=float))
    grad = np.concatenate(grad_parts) if grad_parts else np.zeros(len(par), dtype=float)
    if grad.size != len(par):
        raise RuntimeError(f"Facet regularization gradient size mismatch: expected {len(par)}, got {grad.size}.")
    return value, grad


def facet_regularization_penalty_summary(par, sizes, config) -> pd.DataFrame:
    """Per-facet penalty contribution at the fitted estimates."""
    runtime = config.get("facet_regularization", {})
    if not isinstance(runtime, dict) or not runtime.get("enabled"):
        return pd.DataFrame(columns=["Facet", "PenalizedLevels", "PenaltyValue", "MeanSD", "Scope"])
    params = expand_params(par, sizes, config)
    rows = []
    for facet in config.get("facet_names", []):
        full = np.asarray(params["facets"].get(facet, np.array([], dtype=float)), dtype=float)
        sd = np.asarray(runtime.get("sds", {}).get(facet, [np.nan] * len(full)), dtype=float)
        mean = np.asarray(runtime.get("means", {}).get(facet, [0.0] * len(full)), dtype=float)
        is_active = np.asarray(runtime.get("active", {}).get(facet, [False] * len(full)), dtype=bool)
        valid = is_active & np.isfinite(sd) & (sd > 0)
        penalty = 0.0
        if valid.any():
            penalty = 0.5 * float(np.sum(((full[valid] - mean[valid]) / sd[valid]) ** 2))
        rows.append({
            "Facet": facet,
            "PenalizedLevels": int(valid.sum()),
            "PenaltyValue": penalty,
            "MeanSD": float(np.nanmean(sd[valid])) if valid.any() else np.nan,
            "Scope": "free facet main effects only",
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out.loc[len(out)] = {
            "Facet": "TOTAL",
            "PenalizedLevels": int(out["PenalizedLevels"].sum()),
            "PenaltyValue": float(out["PenaltyValue"].sum()),
            "MeanSD": np.nan,
            "Scope": "normal MAP-style penalty; not full Bayesian sampling",
        }
    return out


def _blank_or_na_mask(series: pd.Series) -> pd.Series:
    """Return True for true NA or whitespace-only cells."""
    s = pd.Series(series)
    text = s.astype(object).where(~s.isna(), "").map(lambda x: str(x).strip())
    return s.isna() | text.eq("")


def normalize_response_input_data(
    data: pd.DataFrame,
    *,
    person_col: str,
    facet_cols: list[str],
    score_col: str,
    weight_col: str | None = None,
) -> pd.DataFrame:
    """Normalize blank selected response cells to NA before model preparation."""
    out = data.copy()
    selected = [person_col, score_col, *list(facet_cols)]
    if weight_col:
        selected.append(weight_col)
    for col in dict.fromkeys(selected):
        if col in out.columns:
            mask = _blank_or_na_mask(out[col])
            if mask.any():
                out.loc[mask, col] = np.nan
    return out


def build_response_data_audit(
    data: pd.DataFrame,
    *,
    person_col: str,
    facet_cols: list[str],
    score_col: str,
    weight_col: str | None = None,
) -> dict:
    """Describe which input rows enter the likelihood and why rows are excluded."""
    facet_cols = list(facet_cols or [])
    selected = [person_col, score_col, *facet_cols]
    if weight_col:
        selected.append(weight_col)
    selected_unique = list(dict.fromkeys([c for c in selected if c is not None]))

    if data is None or not isinstance(data, pd.DataFrame):
        empty = pd.DataFrame()
        summary = pd.DataFrame([
            {
                "Item": "Input rows",
                "Count": 0,
                "Percent": 0.0,
                "Interpretation": "No input data were available.",
            }
        ])
        return {"summary": summary, "rows": empty, "excluded_rows": empty}

    n = len(data)
    audit = pd.DataFrame({
        "InputRow": np.arange(1, n + 1, dtype=int),
        "OriginalIndex": data.index.astype(str),
    })
    for col in selected_unique:
        if col in data.columns:
            # Audit rows follow input positions, including repeated index labels.
            audit[col] = data[col].astype(object).where(~data[col].isna(), "").reset_index(drop=True)

    missing_required_cols = [c for c in selected_unique if c not in data.columns]
    if missing_required_cols:
        audit["IncludedInLikelihood"] = False
        audit["LikelihoodAction"] = "error_before_fit"
        audit["ExclusionReason"] = "Selected column(s) missing: " + ", ".join(missing_required_cols)
        audit["MissingColumns"] = ", ".join(missing_required_cols)
        summary = pd.DataFrame([
            {
                "Item": "Input rows",
                "Count": int(n),
                "Percent": 100.0 if n else 0.0,
                "Interpretation": "Selected columns are missing; estimation cannot start.",
            }
        ])
        return {"summary": summary, "rows": audit, "excluded_rows": audit.copy()}

    person_missing = _blank_or_na_mask(data[person_col])
    score_blank = _blank_or_na_mask(data[score_col])
    score_numeric = pd.to_numeric(data[score_col], errors="coerce")
    score_non_numeric = (~score_blank) & score_numeric.isna()
    score_non_finite = score_numeric.notna() & ~np.isfinite(score_numeric)
    score_tol = np.sqrt(np.finfo(float).eps)
    score_fractional = (
        score_numeric.notna()
        & np.isfinite(score_numeric)
        & ((score_numeric - np.round(score_numeric)).abs() > score_tol)
    )

    facet_missing_by_col: dict[str, pd.Series] = {}
    facet_missing_any = pd.Series(False, index=data.index)
    for facet in facet_cols:
        mask = _blank_or_na_mask(data[facet])
        facet_missing_by_col[facet] = mask
        facet_missing_any = facet_missing_any | mask

    if weight_col:
        weight_blank = _blank_or_na_mask(data[weight_col])
        weight_numeric = pd.to_numeric(data[weight_col], errors="coerce")
        weight_non_numeric = (~weight_blank) & weight_numeric.isna()
        weight_missing_or_non_numeric = weight_blank | weight_non_numeric
        weight_non_positive = weight_numeric.notna() & (weight_numeric <= 0)
        weight_non_finite = weight_numeric.notna() & ~np.isfinite(weight_numeric)
    else:
        weight_numeric = pd.Series(np.ones(n, dtype=float), index=data.index)
        weight_missing_or_non_numeric = pd.Series(False, index=data.index)
        weight_non_positive = pd.Series(False, index=data.index)
        weight_non_finite = pd.Series(False, index=data.index)

    score_missing_or_non_numeric = score_blank | score_non_numeric
    excluded = (
        person_missing
        | score_missing_or_non_numeric
        | facet_missing_any
        | weight_missing_or_non_numeric
        | weight_non_positive
    )
    error_before_fit = score_fractional | score_non_finite | weight_non_finite
    included = (~excluded) & (~error_before_fit)

    missing_columns_per_row = [""] * n
    reasons_per_row = [""] * n
    included_array = included.to_numpy(dtype=bool)
    actions = np.where(error_before_fit.to_numpy(dtype=bool), "error_before_fit",
                       np.where(included_array, "included", "excluded"))
    # Read each mask once; ordinary included rows need no exclusion text.
    person_mask, score_blank_mask, score_non_numeric_mask, fractional_mask = (
        mask.to_numpy(dtype=bool) for mask in
        (person_missing, score_blank, score_non_numeric, score_fractional)
    )
    facet_masks = {facet: mask.to_numpy(dtype=bool) for facet, mask in facet_missing_by_col.items()}
    weight_missing_mask = weight_missing_or_non_numeric.to_numpy(dtype=bool)
    weight_non_positive_mask = weight_non_positive.to_numpy(dtype=bool)
    score_non_finite_mask = score_non_finite.to_numpy(dtype=bool)
    weight_non_finite_mask = weight_non_finite.to_numpy(dtype=bool)
    for row in np.flatnonzero(~included_array):
        missing_cols = []
        reasons = []
        if person_mask[row]:
            missing_cols.append(person_col)
            reasons.append("missing person")
        if score_blank_mask[row]:
            missing_cols.append(score_col)
            reasons.append("missing score")
        elif score_non_numeric_mask[row]:
            missing_cols.append(score_col)
            reasons.append("non-numeric score")
        if fractional_mask[row]:
            reasons.append("fractional score")
        if score_non_finite_mask[row]:
            reasons.append("non-finite score")
        for facet, mask in facet_masks.items():
            if mask[row]:
                missing_cols.append(facet)
                reasons.append(f"missing {facet}")
        if weight_col:
            if weight_missing_mask[row]:
                missing_cols.append(weight_col)
                reasons.append("missing or non-numeric weight")
            elif weight_non_finite_mask[row]:
                reasons.append("non-finite weight")
            elif weight_non_positive_mask[row]:
                reasons.append("non-positive weight")
        missing_columns_per_row[row] = ", ".join(dict.fromkeys(missing_cols))
        reasons_per_row[row] = "; ".join(dict.fromkeys(reasons))

    audit["ScoreNumeric"] = score_numeric.to_numpy()
    if weight_col:
        audit["WeightNumeric"] = weight_numeric.to_numpy()
    audit["PersonMissing"] = person_missing.to_numpy(dtype=bool)
    audit["ScoreMissingOrNonNumeric"] = score_missing_or_non_numeric.to_numpy(dtype=bool)
    audit["ScoreFractional"] = score_fractional.to_numpy(dtype=bool)
    audit["ScoreNonFinite"] = score_non_finite_mask
    audit["FacetMissing"] = facet_missing_any.to_numpy(dtype=bool)
    if weight_col:
        audit["WeightMissingOrNonNumeric"] = weight_missing_or_non_numeric.to_numpy(dtype=bool)
        audit["WeightNonPositive"] = weight_non_positive.to_numpy(dtype=bool)
        audit["WeightNonFinite"] = weight_non_finite_mask
    audit["IncludedInLikelihood"] = included.to_numpy(dtype=bool)
    audit["LikelihoodAction"] = actions
    audit["MissingColumns"] = missing_columns_per_row
    audit["ExclusionReason"] = reasons_per_row

    included_count = int(included.sum())
    excluded_count = int(excluded.sum())
    fractional_count = int(score_fractional.sum())
    not_used_count = int((~included).sum())
    included_weight = float(weight_numeric[included].sum()) if included_count else 0.0

    def pct(count: int) -> float:
        return round(100.0 * int(count) / n, 2) if n else 0.0

    summary_rows = [
        {
            "Item": "Input rows",
            "Count": int(n),
            "Percent": 100.0 if n else 0.0,
            "Interpretation": "Rows supplied by the user after optional missing-code recoding.",
        },
        {
            "Item": "Rows included in likelihood",
            "Count": included_count,
            "Percent": pct(included_count),
            "Interpretation": "Rows with non-missing person, facets, integer score, and positive weight.",
        },
        {
            "Item": "Rows not used in likelihood",
            "Count": not_used_count,
            "Percent": pct(not_used_count),
            "Interpretation": "Rows excluded as missing/zero-weight plus rows that trigger pre-fit errors.",
        },
        {
            "Item": "Rows excluded before likelihood",
            "Count": excluded_count,
            "Percent": pct(excluded_count),
            "Interpretation": "Rows treated as structurally absent, not as zero scores.",
        },
        {
            "Item": "Effective included weight",
            "Count": included_weight,
            "Percent": np.nan,
            "Interpretation": "Sum of positive weights for rows entering the likelihood.",
        },
        {
            "Item": "Missing/non-numeric score",
            "Count": int(score_missing_or_non_numeric.sum()),
            "Percent": pct(int(score_missing_or_non_numeric.sum())),
            "Interpretation": "Rows excluded because the score is absent or cannot be parsed.",
        },
        {
            "Item": "Missing person",
            "Count": int(person_missing.sum()),
            "Percent": pct(int(person_missing.sum())),
            "Interpretation": "Rows excluded because the measured object is absent.",
        },
        {
            "Item": "Missing selected facet cell",
            "Count": int(facet_missing_any.sum()),
            "Percent": pct(int(facet_missing_any.sum())),
            "Interpretation": "Rows excluded because at least one selected facet level is absent.",
        },
    ]
    if weight_col:
        summary_rows.extend([
            {
                "Item": "Missing/non-numeric weight",
                "Count": int(weight_missing_or_non_numeric.sum()),
                "Percent": pct(int(weight_missing_or_non_numeric.sum())),
                "Interpretation": "Rows excluded because the selected weight cannot be used.",
            },
            {
                "Item": "Non-positive weight",
                "Count": int(weight_non_positive.sum()),
                "Percent": pct(int(weight_non_positive.sum())),
                "Interpretation": "Rows excluded because weight <= 0 contributes no likelihood.",
            },
        ])
    if fractional_count:
        summary_rows.append({
            "Item": "Fractional score error",
            "Count": fractional_count,
            "Percent": pct(fractional_count),
            "Interpretation": "The run stops until these scores are recoded to integer categories.",
        })
    for label, mask in [("Non-finite score error", score_non_finite),
                        ("Non-finite weight error", weight_non_finite)]:
        count = int(mask.sum())
        if count:
            summary_rows.append({
                "Item": label, "Count": count, "Percent": pct(count),
                "Interpretation": "The run stops until infinite input values are corrected.",
            })

    summary = pd.DataFrame(summary_rows)
    excluded_rows = audit.loc[~audit["IncludedInLikelihood"]].copy()
    return {"summary": summary, "rows": audit, "excluded_rows": excluded_rows}


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

    audit_bundle = build_response_data_audit(
        data,
        person_col=person_col,
        facet_cols=list(facet_cols),
        score_col=score_col,
        weight_col=weight_col,
    )
    data = normalize_response_input_data(
        data,
        person_col=person_col,
        facet_cols=list(facet_cols),
        score_col=score_col,
        weight_col=weight_col,
    )

    cols = [person_col] + list(facet_cols) + [score_col]
    if weight_col:
        cols.append(weight_col)
    df = data[cols].copy()
    rename_map = {person_col: "Person", score_col: "Score"}
    if weight_col:
        rename_map[weight_col] = "Weight"
    df = df.rename(columns=rename_map)

    # v0.2.13-beta guard: if the rename produces duplicate column names
    # (e.g. a facet named "Score" alongside a score_col renamed to
    # "Score"), `df["Score"]` would return a DataFrame and
    # `pd.to_numeric(…)` raises a cryptic TypeError. Detect the
    # collision up front and explain the fix.
    if df.columns.duplicated().any():
        dups = df.columns[df.columns.duplicated()].tolist()
        raise ValueError(
            f"Column rename produced duplicates: {dups!r}. This typically "
            "means one of your facet columns has the same name as the "
            "Person or Score role after renaming. Rename or drop the "
            "conflicting facet column before fitting."
        )

    score_num = pd.to_numeric(df["Score"], errors="coerce")
    if (score_num.notna() & ~np.isfinite(score_num)).any():
        raise ValueError("Score must contain finite integer category codes; infinite values were found.")
    score_tol = np.sqrt(np.finfo(float).eps)
    fractional_score = score_num.notna() & np.isfinite(score_num) & ((score_num - np.round(score_num)).abs() > score_tol)
    if fractional_score.any():
        examples = df.loc[fractional_score, "Score"].astype(str).drop_duplicates().head(5).tolist()
        raise ValueError(
            "Score must contain ordered integer category codes. "
            f"Fractional value(s) were found: {', '.join(examples)}. "
            "Recode the score column before fitting."
        )
    df["Score"] = score_num

    if "Weight" in df.columns:
        df["Weight"] = pd.to_numeric(df["Weight"], errors="coerce")
        if (df["Weight"].notna() & ~np.isfinite(df["Weight"])).any():
            raise ValueError("Weights must be finite; infinite values were found.")
    else:
        df["Weight"] = 1.0

    df = df.dropna()
    df = df[df["Weight"] > 0]
    if df.empty:
        raise ValueError("No valid observations remain after removing missing values and non-positive weights.")

    df["Person"] = df["Person"].astype(str)
    for facet in facet_cols:
        df[facet] = df[facet].astype(str)
    df["Score"] = df["Score"].astype(int)
    observed_score_values = np.sort(df["Score"].unique())
    if len(observed_score_values) < 2:
        raise ValueError("MFRM requires at least two distinct score categories.")

    rating_min_supplied = rating_min is not None
    rating_max_supplied = rating_max is not None
    if rating_min is None:
        rating_min = int(df["Score"].min())
    else:
        rating_min_num = float(rating_min)
        if not np.isfinite(rating_min_num) or abs(rating_min_num - round(rating_min_num)) > score_tol:
            raise ValueError("rating_min must be a single finite integer category value.")
        rating_min = int(round(rating_min_num))
    if rating_max is None:
        rating_max = int(df["Score"].max())
    else:
        rating_max_num = float(rating_max)
        if not np.isfinite(rating_max_num) or abs(rating_max_num - round(rating_max_num)) > score_tol:
            raise ValueError("rating_max must be a single finite integer category value.")
        rating_max = int(round(rating_max_num))
    if rating_max <= rating_min:
        raise ValueError("rating_max must be larger than rating_min.")

    expected_vals = np.arange(rating_min, rating_max + 1)
    out_of_range = observed_score_values[(observed_score_values < rating_min) | (observed_score_values > rating_max)]
    if out_of_range.size:
        raise ValueError(
            "Observed score categories fall outside the supplied rating range: "
            + ", ".join(map(str, out_of_range.tolist()))
            + ". Adjust rating_min/rating_max or recode the score column."
        )

    score_messages = []
    preserve_score_support = bool(keep_original)
    if not keep_original:
        score_vals = np.sort(df["Score"].unique())
        observed_contiguous = np.array_equal(score_vals, np.arange(score_vals[0], score_vals[-1] + 1))
        explicit_rating_range = rating_min_supplied or rating_max_supplied
        boundary_only_gap = bool(explicit_rating_range and observed_contiguous and np.all(np.isin(score_vals, expected_vals)))
        if not np.array_equal(score_vals, expected_vals) and not boundary_only_gap:
            mapping = {val: rating_min + i for i, val in enumerate(score_vals)}
            df["Score"] = df["Score"].map(mapping).astype(int)
            recoded_vals = np.arange(rating_min, rating_min + len(score_vals))
            score_messages.append(
                "Observed score categories were non-consecutive "
                f"({', '.join(map(str, score_vals.tolist()))}) and were recoded internally to "
                f"{', '.join(map(str, recoded_vals.tolist()))}. "
                "Use Keep original category values to retain zero-count intermediate categories."
            )
            rating_max = rating_min + len(score_vals) - 1
            expected_vals = np.arange(rating_min, rating_max + 1)
        elif boundary_only_gap:
            preserve_score_support = True

    df["score_k"] = df["Score"] - rating_min
    unused_score_categories = sorted(set(range(rating_min, rating_max + 1)) - set(df["Score"].unique().tolist()))
    if unused_score_categories:
        score_messages.append(
            "Zero-count score categories retained in the fitted support: "
            + ", ".join(map(str, unused_score_categories))
            + ". Adjacent thresholds should be treated as weakly identified."
        )

    if preserve_score_support:
        score_map = pd.DataFrame({
            "OriginalScore": np.arange(rating_min, rating_max + 1, dtype=int),
            "InternalScore": np.arange(rating_min, rating_max + 1, dtype=int),
        })
    else:
        observed_after = np.sort(df["Score"].unique())
        score_map = pd.DataFrame({
            "OriginalScore": observed_score_values,
            "InternalScore": observed_after,
        })

    df["Person"] = pd.Categorical(df["Person"])
    for facet in facet_cols:
        df[facet] = pd.Categorical(df[facet])

    facet_levels = {facet: list(df[facet].cat.categories) for facet in facet_cols}

    n_obs = int(df["Weight"].sum()) if "Weight" in df.columns and df["Weight"].notna().any() else len(df)
    n_person = df["Person"].nunique()

    return {
        "data": df,
        "rating_min": rating_min,
        "rating_max": rating_max,
        "score_map": score_map,
        "unused_score_categories": unused_score_categories,
        "score_messages": score_messages,
        "response_data_audit": audit_bundle,
        "row_audit": audit_bundle.get("rows", pd.DataFrame()),
        "excluded_rows": audit_bundle.get("excluded_rows", pd.DataFrame()),
        "audit_summary": audit_bundle.get("summary", pd.DataFrame()),
        "keep_original": bool(keep_original),
        "score_range_explicit": bool(rating_min_supplied or rating_max_supplied),
        "facet_names": list(facet_cols),
        "levels": {"Person": list(df["Person"].cat.categories), **facet_levels},
        "weight_col": "Weight",
        "n_obs": n_obs,
        "n_person": n_person,
    }


_MML_MAX_TENSOR_ELEMENTS = 4_000_000
_GROUP_REDUCTION_CACHE: OrderedDict[tuple[int, int], tuple[object, dict[str, object]]] = OrderedDict()
_GROUP_REDUCTION_CACHE_MAX = 64


def _build_group_reduction(group_idx: np.ndarray, n_groups: int) -> dict[str, object]:
    group_arr = np.asarray(group_idx, dtype=int)
    n_groups = int(max(n_groups, 0))
    if group_arr.size == 0 or n_groups == 0:
        return {
            "n_groups": n_groups,
            "n_obs": int(group_arr.size),
            "order": np.empty(0, dtype=int),
            "nonempty": np.empty(0, dtype=int),
            "reduceat": np.empty(0, dtype=int),
            "counts": np.zeros(n_groups, dtype=int),
        }
    valid_mask = (group_arr >= 0) & (group_arr < n_groups)
    if np.all(valid_mask):
        valid_rows = None
        g_valid = group_arr
    else:
        valid_rows = np.flatnonzero(valid_mask)
        g_valid = group_arr[valid_rows]
    order_local = np.argsort(g_valid, kind="mergesort")
    order = order_local if valid_rows is None else valid_rows[order_local]
    counts = np.bincount(g_valid, minlength=n_groups)
    nonempty = np.flatnonzero(counts > 0)
    reduceat = np.concatenate(([0], np.cumsum(counts[nonempty])[:-1])) if nonempty.size else np.empty(0, dtype=int)
    return {
        "n_groups": n_groups,
        "n_obs": int(group_arr.size),
        "order": order,
        "nonempty": nonempty,
        "reduceat": reduceat,
        "counts": counts,
    }


def _cached_group_reduction(group_arr, n_groups: int) -> dict[str, object]:
    n_groups = int(max(n_groups, 0))
    key = (id(group_arr), n_groups)
    entry = _GROUP_REDUCTION_CACHE.get(key)
    if entry is not None:
        arr_ref, red = entry
        if arr_ref() is group_arr and red["n_obs"] == len(group_arr):
            _GROUP_REDUCTION_CACHE.move_to_end(key)
            return red
        _GROUP_REDUCTION_CACHE.pop(key, None)
    red = _build_group_reduction(group_arr, n_groups)
    try:
        import weakref
        arr_ref = weakref.ref(group_arr)
    except TypeError:
        return red
    _GROUP_REDUCTION_CACHE[key] = (arr_ref, red)
    while len(_GROUP_REDUCTION_CACHE) > _GROUP_REDUCTION_CACHE_MAX:
        _GROUP_REDUCTION_CACHE.popitem(last=False)
    return red


def _get_person_reduction(idx: dict, n_persons: int) -> dict[str, object]:
    return _cached_group_reduction(idx["person"], n_persons)


def _get_step_reduction(idx: dict, n_step_groups: int) -> dict[str, object]:
    return _cached_group_reduction(idx["step_idx"], n_step_groups)


def _sum_2d_by_group(values_2d: np.ndarray, reduction: dict[str, object]) -> np.ndarray:
    n_groups = int(reduction["n_groups"])
    m_cols = int(values_2d.shape[1])
    out = np.zeros((n_groups, m_cols), dtype=float)
    nonempty = reduction["nonempty"]
    if nonempty.size == 0 or values_2d.shape[0] == 0:
        return out
    ordered = values_2d[reduction["order"]]
    out[nonempty] = np.add.reduceat(ordered, reduction["reduceat"], axis=0)
    return out


def _build_step_cum_mat(steps_mat: np.ndarray) -> np.ndarray:
    arr = np.asarray(steps_mat, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[0] == 0:
        return np.zeros((0, arr.shape[1] + 1), dtype=float)
    if arr.shape[1] == 0:
        return np.zeros((arr.shape[0], 1), dtype=float)
    return np.hstack([np.zeros((arr.shape[0], 1), dtype=float), np.cumsum(arr, axis=1)])


def _log_denom_and_probs(log_num: np.ndarray, *, axis: int = -1, need_probs: bool = True):
    max_log = np.max(log_num, axis=axis, keepdims=True)
    if not np.all(np.isfinite(max_log)):
        log_denom = logsumexp(log_num, axis=axis)
        if not need_probs:
            return log_denom, None
        probs = np.exp(log_num - np.expand_dims(log_denom, axis=axis))
        return log_denom, probs
    shifted = np.exp(log_num - max_log)
    sum_exp = np.sum(shifted, axis=axis, keepdims=True)
    log_denom = np.squeeze(max_log + np.log(sum_exp), axis=axis)
    if not need_probs:
        return log_denom, None
    return log_denom, shifted / sum_exp


def build_indices(prep, step_facet=None, slope_facet=None):
    df = prep["data"]
    person_idx = df["Person"].cat.codes.to_numpy(dtype=int)
    facets_idx = {facet: df[facet].cat.codes.to_numpy(dtype=int) for facet in prep["facet_names"]}
    step_idx = df[step_facet].cat.codes.to_numpy(dtype=int) if step_facet else None
    slope_idx = df[slope_facet].cat.codes.to_numpy(dtype=int) if slope_facet else None
    n_persons = len(prep["levels"]["Person"])
    person_red = _cached_group_reduction(person_idx, n_persons)
    rows_by_person = [np.empty(0, dtype=int) for _ in range(n_persons)]
    if person_red["nonempty"].size:
        splits = np.cumsum(person_red["counts"][person_red["nonempty"]])[:-1]
        for p_i, chunk in zip(person_red["nonempty"], np.split(person_red["order"], splits)):
            rows_by_person[int(p_i)] = chunk
    if step_facet and step_idx is not None:
        _cached_group_reduction(step_idx, len(prep["levels"][step_facet]))
    return {
        "person": person_idx,
        "facets": facets_idx,
        "step_idx": step_idx,
        "slope_idx": slope_idx,
        "score_k": df["score_k"].to_numpy(dtype=int),
        "weight": df["Weight"].to_numpy(dtype=float) if "Weight" in df.columns else np.ones(len(df)),
        "rows_by_person": rows_by_person,
    }


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


def loglik_gpcm(eta, score_k, step_cum_mat, step_idx, slopes, slope_idx, weight=None):
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
        slope_obs = slopes[slope_idx[rows]]
        log_num = slope_obs[:, None] * linear
        log_denom = logsumexp(log_num, axis=1)
        diff = log_num[np.arange(len(rows)), score_k[rows]] - log_denom
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


def category_prob_gpcm(eta, step_cum_mat, step_idx, slopes, slope_idx):
    if eta.size == 0:
        return np.zeros((0, step_cum_mat.shape[1]))
    if slope_idx is None:
        slope_idx = step_idx
    k_cat = step_cum_mat.shape[1]
    k_vals = np.arange(k_cat, dtype=float)
    probs = np.zeros((len(eta), k_cat))
    slopes = np.asarray(slopes, dtype=float)
    for c_idx in range(step_cum_mat.shape[0]):
        rows = np.where(step_idx == c_idx)[0]
        if rows.size == 0:
            continue
        linear = np.outer(eta[rows], k_vals) - step_cum_mat[c_idx]
        slope_obs = slopes[slope_idx[rows]]
        log_num = slope_obs[:, None] * linear
        log_denom = logsumexp(log_num, axis=1)
        probs[rows, :] = np.exp(log_num - log_denom[:, None])
    return probs


def _resolve_zstd_transform(whexact=False, zstd_transform=None):
    """Retain the historical boolean without implying FACETS WHEXACT parity."""
    allowed = {"wilson_hilferty", "wilson_hilferty_no_offset", "linear_normal"}
    resolved = zstd_transform or ("linear_normal" if whexact else "wilson_hilferty")
    if resolved not in allowed:
        raise ValueError(f"Unknown ZSTD transform: {resolved}")
    if whexact and resolved != "linear_normal":
        raise ValueError("Use zstd_transform alone; whexact=True is the legacy linear approximation.")
    return resolved


def _zstd_value(mnsq, df, transform):
    if transform == "linear_normal":
        return (mnsq - 1) * np.sqrt(df / 2)
    if transform == "wilson_hilferty_no_offset":
        return (mnsq ** (1 / 3) - 1) / np.sqrt(2 / (9 * df))
    return (mnsq ** (1 / 3) - (1 - 2 / (9 * df))) / np.sqrt(2 / (9 * df))


def zstd_from_mnsq(mnsq, df, whexact=False, *, zstd_transform=None):
    """Wilson-Hilferty (1931) standardised MnSq, engine convention.

    The engine convention guards against ``df < 1``: in that regime the
    ``1 - 2/(9 df)`` centring term becomes dominant and the resulting
    ZSTD can flip sign relative to the underlying MnSq, producing
    spurious large values. This legacy convention returns NaN there. The
    FACETS-convention helper ``zstd_from_mnsq_facets`` relaxes this to
    ``df > 0`` (the Wright-Masters Welch-Satterthwaite d.f. can be
    fractional by design) and applies an explicit cap for the
    pathological cells.

    Prefer ``zstd_transform``: ``wilson_hilferty`` includes the offset,
    ``wilson_hilferty_no_offset`` omits it, and ``linear_normal`` uses
    the historical linear approximation. The legacy ``whexact`` boolean
    remains unchanged (False = full WH, True = linear); it is not the
    FACETS Yes/No setting.
    """
    transform = _resolve_zstd_transform(whexact, zstd_transform)
    if not np.isfinite(mnsq) or mnsq < 0 or not np.isfinite(df) or df < 1:
        return np.nan
    return _zstd_value(mnsq, df, transform)


_FIT_DF_METHOD_CHOICES = ("engine", "facets", "both")


DEFAULT_FIT_DF_METHOD = "facets"


DEFAULT_FACETS_ZSTD_CAP = 9.0


def _resolve_fit_df_method(fit_df_method) -> str:
    """Coerce ``fit_df_method`` to a known token, defaulting to FACETS."""
    token = str(fit_df_method or DEFAULT_FIT_DF_METHOD).strip().lower()
    if token in _FIT_DF_METHOD_CHOICES:
        return token
    return DEFAULT_FIT_DF_METHOD


def fit_df_model_metadata(model: str | None, fit_df_method: str | None = None) -> dict[str, str]:
    """Describe the model scope of the selected fit-d.f. convention.

    FACETS/Wright-Masters d.f. is the primary reporting convention for RSM
    and PCM. The same fourth-moment calculation is available for bounded
    GPCM fits, but is explicitly labelled as a FACETS-style approximation
    because free discrimination parameters fall outside the strict FACETS
    Rasch-family parameterisation.
    """
    model_token = str(model or "RSM").strip().upper()
    method = _resolve_fit_df_method(fit_df_method)
    if method == "engine":
        return {
            "applicability": "engine_homogeneous_variance_approximation",
            "warning": "",
        }
    if model_token == "GPCM":
        return {
            "applicability": "facets_style_approximation_for_gpcm",
            "warning": (
                "GPCM uses free discrimination parameters, so its FACETS/Wright-Masters "
                "d.f. and ZSTD are a FACETS-style approximation. Interpret the FACETS-primary "
                "columns with the engine sidecars and do not claim exact FACETS equivalence."
            ),
        }
    return {
        "applicability": "facets_wright_masters_rasch_family",
        "warning": "",
    }


def zstd_from_mnsq_facets(
    mnsq,
    df,
    whexact: bool = False,
    cap: float = DEFAULT_FACETS_ZSTD_CAP,
    *,
    zstd_transform=None,
):
    """FACETS-style ZSTD with explicit cap on extreme values.

    Implements the same Wilson-Hilferty cube-root transformation as
    ``zstd_from_mnsq`` but with the FACETS convention:

    * ``df > 0`` (rather than ``>= 1``) is accepted because the
      Wright-Masters Welch-Satterthwaite d.f. can be fractional;
    * the result is clamped to ``[-cap, +cap]`` (FACETS default 9.0).

    For FACETS WHEXACT=Yes use ``zstd_transform="wilson_hilferty"``;
    for No use ``"wilson_hilferty_no_offset"``. Omit ``whexact`` when
    using these names. Its historical True value means linear normal
    approximation and is preserved for compatibility.
    """
    # WHEXACT in FACETS has different boolean semantics. Use the named
    # transform: full WH corresponds to Yes, no-offset WH corresponds to No.
    transform = _resolve_zstd_transform(whexact, zstd_transform)
    if not np.isfinite(mnsq) or not np.isfinite(df) or mnsq <= 0 or df <= 0:
        return np.nan
    z = _zstd_value(mnsq, df, transform)
    if cap is None or not np.isfinite(cap) or cap <= 0:
        return float(z)
    return float(np.clip(z, -float(cap), float(cap)))


def _facets_fit_df_terms(var, fourth, weight):
    """Return the Wright-Masters denominators for the FACETS d.f.

    The denominators are ``sum w * (M4 - sigma^4)`` for Infit and
    ``sum w * (M4 / sigma^4 - 1)`` for Outfit, evaluated only over
    observations with strictly positive variance and finite fourth
    moment. ``weight`` defaults to all-ones if ``None`` is passed.

    The linear-``w`` form matches the mfrmr 0.2.0 reference exactly
    and is the formula reported in the Wright & Masters (1982,
    Eq. 4.27) literature. The independence-based derivation under
    non-unit weights would substitute ``w^2`` for ``w``; for the
    typical MFRM case of unit observation weights the two forms are
    identical, and the linear-``w`` form is what FACETS / Wright-
    Masters report. The bundle's ``settings`` field names the
    formula in plain English so downstream consumers can name it
    correctly in a manuscript.
    """
    var = np.asarray(var, dtype=float)
    fourth = np.asarray(fourth, dtype=float)
    if weight is None:
        weight = np.ones_like(var)
    else:
        weight = np.asarray(weight, dtype=float)
    if weight.shape != var.shape:
        weight = np.broadcast_to(weight, var.shape).astype(float, copy=True)
    safe_var = np.where(np.isfinite(var) & (var > 1e-12), var, np.nan)
    valid = (
        np.isfinite(safe_var)
        & np.isfinite(fourth)
        & np.isfinite(weight)
        & (weight > 0)
    )
    denom_infit = float(np.nansum(weight[valid] * (fourth[valid] - safe_var[valid] ** 2)))
    denom_outfit = float(
        np.nansum(weight[valid] * (fourth[valid] / (safe_var[valid] ** 2) - 1))
    )
    return denom_infit, denom_outfit


def _facets_fit_df(sum_var_w: float, sum_w: float, denom_infit: float, denom_outfit: float):
    """Return ``(DF_Infit_FACETS, DF_Outfit_FACETS)`` from cell totals."""
    if (
        np.isfinite(denom_infit)
        and denom_infit > 1e-12
        and np.isfinite(sum_var_w)
        and sum_var_w > 0
    ):
        df_infit = 2.0 * (sum_var_w ** 2) / denom_infit
    else:
        df_infit = np.nan
    if (
        np.isfinite(denom_outfit)
        and denom_outfit > 1e-12
        and np.isfinite(sum_w)
        and sum_w > 0
    ):
        df_outfit = 2.0 * (sum_w ** 2) / denom_outfit
    else:
        df_outfit = np.nan
    return float(df_infit), float(df_outfit)


def fit_zstd_transform_label(whexact: bool = False, *, zstd_transform=None) -> str:
    """Human-readable label for the ZSTD transform used in a report."""
    return {"linear_normal": "linear normal approximation", "wilson_hilferty": "Wilson-Hilferty",
            "wilson_hilferty_no_offset": "Wilson-Hilferty without offset"}[_resolve_zstd_transform(whexact, zstd_transform)]


def _apply_fit_df_method(
    tbl: pd.DataFrame,
    fit_df_method: str = DEFAULT_FIT_DF_METHOD,
    whexact: bool = False,
    facets_zstd_cap: float = DEFAULT_FACETS_ZSTD_CAP,
    *,
    zstd_transform=None,
) -> pd.DataFrame:
    """Dispatch a fit table between engine, FACETS, or both d.f. conventions.

    Expects ``tbl`` to carry ``Infit``, ``Outfit``, ``InfitZSTD``,
    ``OutfitZSTD``, ``DF_Infit``, ``DF_Outfit``, ``DF_Infit_FACETS``,
    ``DF_Outfit_FACETS``. Returns a frame whose primary columns reflect
    the requested method:

    * ``"engine"`` (explicit legacy convention): drops every FACETS
      column and the status fields, leaving the engine values as-is.
    * ``"facets"`` (default): primary ``InfitZSTD`` / ``OutfitZSTD`` / ``DF_*``
      become the FACETS values; original engine values are preserved
      under ``*_ENGINE`` suffix columns; FACETS values also under
      ``*_FACETS`` for symmetry.
    * ``"both"``: keeps both ``*_ENGINE`` and ``*_FACETS`` suffix
      columns; primary columns continue to hold the engine values.

    Adds ``FitDfMethod``, ``FitZSTDTransform``, and ``FitZSTDCap``
    metadata columns in every mode so downstream consumers
    can read the convention off the row.
    """
    method = _resolve_fit_df_method(fit_df_method)
    required = (
        "Infit", "Outfit", "InfitZSTD", "OutfitZSTD",
        "DF_Infit", "DF_Outfit", "DF_Infit_FACETS", "DF_Outfit_FACETS",
    )
    if not all(c in tbl.columns for c in required):
        return tbl

    out = tbl.copy()
    out["DF_Infit_ENGINE"] = out["DF_Infit"]
    out["DF_Outfit_ENGINE"] = out["DF_Outfit"]
    out["InfitZSTD_ENGINE"] = out["InfitZSTD"]
    out["OutfitZSTD_ENGINE"] = out["OutfitZSTD"]
    out["InfitZSTD_FACETS"] = [
        zstd_from_mnsq_facets(m, d, whexact=whexact, cap=facets_zstd_cap, zstd_transform=zstd_transform)
        for m, d in zip(out["Infit"], out["DF_Infit_FACETS"])
    ]
    out["OutfitZSTD_FACETS"] = [
        zstd_from_mnsq_facets(m, d, whexact=whexact, cap=facets_zstd_cap, zstd_transform=zstd_transform)
        for m, d in zip(out["Outfit"], out["DF_Outfit_FACETS"])
    ]

    if method == "facets":
        out["DF_Infit"] = out["DF_Infit_FACETS"]
        out["DF_Outfit"] = out["DF_Outfit_FACETS"]
        out["InfitZSTD"] = out["InfitZSTD_FACETS"]
        out["OutfitZSTD"] = out["OutfitZSTD_FACETS"]
        out["FitDfMethod"] = "facets_wright_masters"
        out["FitDfFormula"] = "linear_w_wright_masters_eq_4_27"
    elif method == "both":
        out["FitDfMethod"] = "engine_primary_facets_available"
        out["FitDfFormula"] = "linear_w_wright_masters_eq_4_27"
    else:  # engine
        drop_cols = [
            "DF_Infit_ENGINE", "DF_Outfit_ENGINE",
            "InfitZSTD_ENGINE", "OutfitZSTD_ENGINE",
            "DF_Infit_FACETS", "DF_Outfit_FACETS",
            "InfitZSTD_FACETS", "OutfitZSTD_FACETS",
        ]
        out = out.drop(columns=[c for c in drop_cols if c in out.columns])
        out["FitDfMethod"] = "engine_homogeneous_variance_approximation"
        out["FitZSTDTransform"] = fit_zstd_transform_label(whexact, zstd_transform=zstd_transform)
        out["FitZSTDCap"] = np.nan
        return out

    out["FitZSTDTransform"] = fit_zstd_transform_label(whexact, zstd_transform=zstd_transform)
    out["FitZSTDCap"] = (
        float(facets_zstd_cap)
        if (facets_zstd_cap is not None and np.isfinite(facets_zstd_cap) and facets_zstd_cap > 0)
        else np.nan
    )
    return out


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
    elif config["model"] == "PCM":
        step_cum_mat = _build_step_cum_mat(params["steps_mat"])
        ll = loglik_pcm(eta, idx["score_k"], step_cum_mat, idx["step_idx"], weight=idx.get("weight"))
    else:
        step_cum_mat = _build_step_cum_mat(params["steps_mat"])
        ll = loglik_gpcm(
            eta,
            idx["score_k"],
            step_cum_mat,
            idx["step_idx"],
            params["slopes"],
            idx.get("slope_idx"),
            weight=idx.get("weight"),
        )
    return -ll


def mfrm_loglik_jmle_value_grad(par, idx, config, sizes):
    params = expand_params(par, sizes, config)
    eta = compute_eta(idx, params, config)
    score_k = idx["score_k"]
    weight = idx.get("weight")
    if weight is None:
        weight = np.ones_like(score_k, dtype=float)
    else:
        weight = np.asarray(weight, dtype=float)

    n_obs = len(eta)
    k_cat = int(config["n_cat"])
    k_vals = np.arange(k_cat, dtype=float)
    slope_obs = np.ones(n_obs, dtype=float)
    gpcm_linear = None

    if n_obs == 0:
        ll_obs = np.zeros(0, dtype=float)
        probs = np.zeros((0, k_cat), dtype=float)
        if config["model"] == "GPCM":
            gpcm_linear = np.zeros((0, k_cat), dtype=float)
    elif config["model"] == "RSM":
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
        log_num = np.outer(eta, k_vals) - step_cum
        log_denom, probs = _log_denom_and_probs(log_num, axis=1, need_probs=True)
        ll_obs = log_num[np.arange(n_obs), score_k] - log_denom
    elif config["model"] == "PCM":
        step_cum_mat = _build_step_cum_mat(params["steps_mat"])
        log_num = np.outer(eta, k_vals) - step_cum_mat[idx["step_idx"]]
        log_denom, probs = _log_denom_and_probs(log_num, axis=1, need_probs=True)
        ll_obs = log_num[np.arange(n_obs), score_k] - log_denom
    else:
        step_cum_mat = _build_step_cum_mat(params["steps_mat"])
        slope_idx = idx.get("slope_idx")
        if slope_idx is None:
            slope_idx = idx["step_idx"]
        slope_obs = params["slopes"][slope_idx]
        gpcm_linear = np.outer(eta, k_vals) - step_cum_mat[idx["step_idx"]]
        log_num = slope_obs[:, None] * gpcm_linear
        log_denom, probs = _log_denom_and_probs(log_num, axis=1, need_probs=True)
        ll_obs = log_num[np.arange(n_obs), score_k] - log_denom

    nll = -float(np.sum(weight * ll_obs))
    expected = probs @ k_vals
    eta_grad = -weight * slope_obs * (score_k - expected)

    grad_parts = []
    if sizes.get("theta", 0) > 0:
        theta_levels = len(config["theta_spec"]["levels"])
        theta_full_grad = np.bincount(idx["person"], weights=eta_grad, minlength=theta_levels)
        grad_parts.append(collapse_facet_gradient(theta_full_grad, config["theta_spec"]))

    facet_signs = config.get("facet_signs", {})
    for facet in config["facet_names"]:
        if sizes.get(facet, 0) == 0:
            continue
        sign = facet_signs.get(facet, -1)
        facet_levels = len(config["facet_specs"][facet]["levels"])
        facet_full_grad = np.bincount(
            idx["facets"][facet],
            weights=eta_grad * sign,
            minlength=facet_levels,
        )
        grad_parts.append(collapse_facet_gradient(facet_full_grad, config["facet_specs"][facet]))

    n_steps = max(k_cat - 1, 0)
    if n_steps > 0:
        thresholds = np.arange(n_steps)
        observed_gt = score_k[:, None] > thresholds[None, :]
        prob_gt = np.cumsum(probs[:, :0:-1], axis=1)[:, ::-1]
        if config["model"] == "RSM":
            step_full_grad = np.sum(weight[:, None] * (observed_gt - prob_gt), axis=0)
            grad_parts.append(collapse_sum_zero_gradient(step_full_grad))
        else:
            step_facet_n = len(config["facet_levels"][config["step_facet"]])
            step_red = _get_step_reduction(idx, step_facet_n)
            row_factor = weight if config["model"] == "PCM" else (weight * slope_obs)
            step_mat_grad = _sum_2d_by_group(row_factor[:, None] * (observed_gt - prob_gt), step_red)
            grad_parts.append(np.vstack([
                collapse_sum_zero_gradient(row)
                for row in step_mat_grad
            ]).reshape(-1))

    if config["model"] == "GPCM" and sizes.get("log_slopes", 0) > 0:
        slope_idx = idx.get("slope_idx")
        if slope_idx is None:
            slope_idx = idx["step_idx"]
        obs_linear = gpcm_linear[np.arange(n_obs), score_k] if n_obs else np.zeros(0, dtype=float)
        expected_linear = np.sum(probs * gpcm_linear, axis=1) if n_obs else np.zeros(0, dtype=float)
        slope_levels = len(config["facet_levels"][config["slope_facet"]])
        slope_full_grad = np.bincount(
            slope_idx,
            weights=-weight * slope_obs * (obs_linear - expected_linear),
            minlength=slope_levels,
        )
        grad_parts.append(collapse_sum_zero_gradient(slope_full_grad))

    grad = np.concatenate([g for g in grad_parts if g.size]) if grad_parts else np.array([], dtype=float)
    if grad.size != len(par):
        raise RuntimeError(f"Internal gradient size mismatch: expected {len(par)}, got {grad.size}.")
    penalty_value, penalty_grad = facet_regularization_value_grad(par, sizes, config)
    if config.get("facet_regularization", {}).get("enabled"):
        nll += penalty_value
        grad = grad + penalty_grad
    return nll, grad


def _mml_person_count(idx: dict, config: dict) -> int:
    rows_by_person = idx.get("rows_by_person")
    if rows_by_person is not None:
        return len(rows_by_person)
    return int(config.get("n_person", 0))


def _mml_pop_mu_per_person(params: dict, config: dict, n_persons: int) -> np.ndarray:
    pop_mu_raw = compute_population_mu(params, config)
    pop_mu_person = np.zeros(int(max(n_persons, 0)), dtype=float)
    if pop_mu_raw.size and n_persons > 0:
        n_copy = min(int(n_persons), int(pop_mu_raw.size))
        pop_mu_person[:n_copy] = pop_mu_raw[:n_copy]
    return pop_mu_person


def _compute_mml_person_node_loglik(
    idx: dict,
    config: dict,
    params: dict,
    quad: dict,
    *,
    n_persons: int | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    if n_persons is None:
        n_persons = _mml_person_count(idx, config)
    person_red = _get_person_reduction(idx, n_persons)
    nodes = np.asarray(quad["nodes"], dtype=float)
    q_len = len(nodes)
    score_k = idx["score_k"]
    n_obs = len(score_k)
    if n_obs == 0 or n_persons == 0 or q_len == 0:
        return np.zeros((n_persons, q_len), dtype=float), person_red

    base_eta = compute_base_eta(idx, params, config)
    pop_mu_person = _mml_pop_mu_per_person(params, config, n_persons)
    base_shift = pop_mu_person[idx["person"]] + base_eta
    weight = idx.get("weight")
    w_arr = None if weight is None else np.asarray(weight, dtype=float)
    k_cat = int(config["n_cat"])
    k_vals = np.arange(k_cat, dtype=float)

    if config["model"] == "RSM":
        step_cum_obs = np.concatenate([[0.0], np.cumsum(params["steps"])])[None, :]
        slope_obs = None
    elif config["model"] == "PCM":
        step_cum_obs = _build_step_cum_mat(params["steps_mat"])[idx["step_idx"]]
        slope_obs = None
    else:
        step_cum_obs = _build_step_cum_mat(params["steps_mat"])[idx["step_idx"]]
        slope_idx = idx.get("slope_idx") if idx.get("slope_idx") is not None else idx["step_idx"]
        slope_obs = np.asarray(params["slopes"], dtype=float)[slope_idx]

    if n_obs * q_len * k_cat <= _MML_MAX_TENSOR_ELEMENTS:
        eta_nq = base_shift[:, None] + nodes[None, :]
        linear = eta_nq[:, :, None] * k_vals[None, None, :] - step_cum_obs[:, None, :]
        log_num = linear if slope_obs is None else (slope_obs[:, None, None] * linear)
        log_denom, _ = _log_denom_and_probs(log_num, axis=2, need_probs=False)
        obs_log_num = np.take_along_axis(log_num, score_k[:, None, None], axis=2)[:, :, 0]
        ll_obs_nq = obs_log_num - log_denom
        w_ll_nq = ll_obs_nq if w_arr is None else (w_arr[:, None] * ll_obs_nq)
        return _sum_2d_by_group(w_ll_nq, person_red), person_red

    log_lik_mat = np.zeros((n_persons, q_len), dtype=float)
    row_arange = np.arange(n_obs)
    step_2d = step_cum_obs[0] if config["model"] == "RSM" else step_cum_obs
    for q_i, theta_q in enumerate(nodes):
        eta_q = theta_q + base_shift
        linear_q = np.outer(eta_q, k_vals) - step_2d
        log_num_q = linear_q if slope_obs is None else (slope_obs[:, None] * linear_q)
        log_denom_q, _ = _log_denom_and_probs(log_num_q, axis=1, need_probs=False)
        diff_q = log_num_q[row_arange, score_k] - log_denom_q
        if w_arr is not None:
            diff_q = diff_q * w_arr
        log_lik_mat[:, q_i] = np.bincount(idx["person"], weights=diff_q, minlength=n_persons)
    return log_lik_mat, person_red


def mfrm_loglik_mml(par, idx, config, sizes, quad):
    """Marginal log-likelihood (negative) for MML — used in M-step."""
    n = len(idx["score_k"])
    if n == 0:
        return 0.0
    params = expand_params(par, sizes, config)
    log_lik_mat, person_red = _compute_mml_person_node_loglik(idx, config, params, quad)
    nonempty = person_red["nonempty"]
    if nonempty.size == 0:
        return 0.0
    log_w = _log_mml_quadrature_weights(quad)
    ll_person = logsumexp(log_w[None, :] + log_lik_mat[nonempty], axis=1)
    return -float(np.sum(ll_person))


def mfrm_loglik_mml_value_grad(par, idx, config, sizes, quad, *, include_log_sigma=False):
    """Direct MML negative marginal log-likelihood with analytical gradient.

    The gradient uses Fisher's identity: at the current parameter vector, the
    gradient of the marginal log-likelihood equals the posterior expectation of
    the complete-data score. When the 3D observation × node × category tensor
    fits within ``_MML_MAX_TENSOR_ELEMENTS``, probabilities and log-likelihoods
    are evaluated once and reused for both the E-step posteriors and the M-step
    score expectation.

    Development RSM/PCM/GPCM callers may append d(NLL)/d(log sigma). This assumes
    quad["nodes"] = sigma * z with fixed standard-normal nodes and weights;
    population means are separate and do not move with sigma.
    """
    if include_log_sigma and config["model"] not in {"RSM", "PCM", "GPCM"}:
        raise ValueError("Analytic log-SD gradient requires a native response model")
    params = expand_params(par, sizes, config)
    n_obs = len(idx["score_k"])
    q_len = len(quad["nodes"])
    k_cat = int(config["n_cat"])
    if n_obs > 0 and q_len > 0 and (n_obs * q_len * k_cat <= _MML_MAX_TENSOR_ELEMENTS):
        marginal_ll, _, grad = _mml_joint_or_mstep_3d(
            par,
            params,
            idx,
            config,
            sizes,
            quad,
            post_weights=None,
            include_log_sigma=include_log_sigma,
        )
        penalty_value, _ = facet_regularization_value_grad(par, sizes, config)
        return -marginal_ll + penalty_value, grad

    post_weights, marginal_ll = _e_step_posteriors(idx, config, params, quad)
    _, grad = _m_step_expected_ll_value_grad(
        par, idx, config, sizes, quad, post_weights, include_log_sigma=include_log_sigma,
    )
    penalty_value, _ = facet_regularization_value_grad(par, sizes, config)
    return -marginal_ll + penalty_value, grad


def mfrm_direct_mml(start, idx, config, sizes, quad, maxit=400, reltol=1e-6):
    """Direct analytical-gradient MML optimization."""
    method, options = optimizer_method_and_options(
        sizes,
        config,
        maxiter=maxit,
        gtol=reltol,
        ftol=reltol,
    )
    opt = minimize(
        mfrm_loglik_mml_value_grad,
        np.array(start, dtype=float, copy=True),
        args=(idx, config, sizes, quad),
        jac=True,
        method=method,
        bounds=build_optimizer_bounds(sizes, config),
        constraints=build_optimizer_constraints(sizes, config),
        options=options,
    )
    opt.mml_engine = "direct"
    raw_nll = mfrm_loglik_mml(opt.x, idx, config, sizes, quad)
    penalty_value, _ = facet_regularization_value_grad(opt.x, sizes, config)
    opt.raw_fun = raw_nll
    opt.penalty_value = penalty_value
    opt.ll_trace = [-float(raw_nll)] if np.isfinite(raw_nll) else []
    opt.gradient_norm = float(np.linalg.norm(opt.jac)) if getattr(opt, "jac", None) is not None else np.nan
    return opt


def mfrm_hybrid_mml(start, idx, config, sizes, quad, maxit=400, reltol=1e-6):
    """Hybrid MML: short EM warm start followed by direct MML optimization."""
    warm_maxit = max(2, min(10, int(maxit // 4))) if maxit >= 4 else max(1, int(maxit))
    direct_maxit = max(1, int(maxit) - warm_maxit)
    warm = mfrm_em_mml(start, idx, config, sizes, quad, maxit=warm_maxit, reltol=reltol)
    opt = mfrm_direct_mml(warm.x, idx, config, sizes, quad, maxit=direct_maxit, reltol=reltol)
    opt.mml_engine = "hybrid"
    opt.em_warmup_failure = getattr(warm, "m_step_failure", None)
    opt.em_warmup_iterations = int(getattr(warm, "nit", warm_maxit))
    opt.direct_iterations = int(getattr(opt, "nit", 0))
    opt.nit = opt.em_warmup_iterations + opt.direct_iterations
    opt.nfev = int(getattr(warm, "nfev", 0)) + int(getattr(opt, "nfev", 0))
    opt.ll_trace = list(getattr(warm, "ll_trace", [])) + list(getattr(opt, "ll_trace", []))
    return opt


def mfrm_auto_mml(start, idx, config, sizes, quad, maxit=400, reltol=1e-6):
    """Auto MML: try hybrid first, then fall back to EM if needed."""
    hybrid = mfrm_hybrid_mml(start, idx, config, sizes, quad, maxit=maxit, reltol=reltol)
    hybrid.auto_attempted = ["hybrid"]
    hybrid.auto_selected = "hybrid"
    if bool(getattr(hybrid, "success", False)) and np.isfinite(getattr(hybrid, "fun", np.nan)):
        return hybrid

    fallback = mfrm_em_mml(start, idx, config, sizes, quad, maxit=maxit, reltol=reltol)
    fallback.mml_engine = "em"
    fallback.auto_attempted = ["hybrid", "em"]
    fallback.auto_selected = "em"
    fallback.auto_fallback_reason = str(getattr(hybrid, "message", "hybrid did not converge") or "hybrid did not converge")

    if bool(getattr(fallback, "success", False)):
        return fallback
    if np.isfinite(getattr(fallback, "fun", np.nan)) and (
        not np.isfinite(getattr(hybrid, "fun", np.nan)) or fallback.fun <= hybrid.fun
    ):
        return fallback

    hybrid.auto_attempted = ["hybrid", "em"]
    hybrid.auto_selected = "hybrid"
    hybrid.auto_fallback_reason = "EM fallback did not improve convergence or likelihood."
    return hybrid


def _mml_inference_hold_reason(config) -> str:
    return _output_qualification.mml_inference_hold_reason(config)


def _free_sd_mml_inference_withheld(config) -> bool:
    return bool(
        isinstance(config, dict)
        and str(config.get("method", "")).upper() == "MML"
        and config.get("estimate_population_sd", False)
    )


def qualify_measure_uncertainty(measures, result):
    """Carry current person and calibration holds without altering source values."""
    if not isinstance(measures, pd.DataFrame) or measures.empty:
        return measures
    out = measures.copy()
    person = out.get("Facet", pd.Series("", index=out.index)).eq("Person")
    # Legacy combined tables may omit flags retained in the person audit.
    boundary = (result.get("person_boundary") or {}) if isinstance(result, dict) else {}
    persons = boundary.get("persons") if isinstance(boundary, dict) else None
    if isinstance(persons, pd.DataFrame) and "Person" in persons and "Level" in out:
        lookup = persons.drop_duplicates("Person").set_index("Person")
        matched = person & out.Level.isin(lookup.index)
        for column in ("ExtremeScoreDirection", "FiniteJMLEEstimate", "PersonInferenceReady", "ReportableEstimate"):
            if column in lookup and matched.any():
                if column not in out:
                    out[column] = pd.Series(None, index=out.index, dtype=object)
                out.loc[matched, column] = out.loc[matched, "Level"].map(lookup[column])
    person_held = pd.Series(False, index=out.index)
    for column in ("FiniteJMLEEstimate", "PersonInferenceReady"):
        if column in out:
            person_held |= person & out[column].eq(False).fillna(False)
    if person_held.any():
        out.loc[person_held, "SE_Status"] = "person_inference_withheld"
        out.loc[person_held, "CI_Status"] = "person_inference_withheld"
        out.loc[person_held, "UncertaintyReasonCode"] = "stat.jmle.person_inference_withheld"
        out.loc[person_held, "UncertaintyCaution"] = (
            "Person inference is withheld; retained estimates, SEs and intervals are technical terminal-point values."
        )
    hold = _output_qualification.calibration_inference_hold_message(result)
    if hold:
        out["SE_Status"] = "calibration_inference_withheld"
        out["CI_Status"] = "calibration_inference_withheld"
        out["UncertaintyReasonCode"] = _output_qualification.calibration_inference_hold_reason(result)
        if any(column in out for column in ("PersonInferenceReady", "ReportableEstimate")):
            out = qualified_person_boundary_table(out, result)
        for prefix in ("FairM", "FairZ"):
            status_column = f"{prefix}_SE_Status"
            if status_column in out or f"{prefix}_CI_Status" in out:
                out[status_column] = "calibration_inference_withheld"
                out[f"{prefix}_CI_Status"] = "calibration_inference_withheld"
                qualification_column = f"{prefix}_SE_QualificationStatus"
                if qualification_column in out:
                    out[qualification_column] = "WITHHELD"
                public_column = f"{prefix}_SE_PublicConclusionAllowed"
                if public_column in out:
                    out[public_column] = False
        note = f"Technical terminal-point calculations only; inference is withheld. {hold}"
        out["UncertaintyCaution"] = note
    return out


def qualified_estimation_summary(result) -> pd.DataFrame:
    """Report explicit IC units and inference holds on a copy of saved evidence."""
    summary = result.get("summary") if isinstance(result, dict) else None
    if not isinstance(summary, pd.DataFrame):
        return pd.DataFrame()
    out = summary.copy()
    if not out.empty:
        if "Converged" in out:
            stopped = _output_qualification.optimizer_termination_established(result)
            if not out["Converged"].map(lambda value: isinstance(value, (bool, np.bool_))
                                       and bool(value) == stopped).all() and "SavedConverged" not in out:
                out["SavedConverged"] = out["Converged"]
            out["Converged"] = stopped
        config = result.get("config", {})
        config = config if isinstance(config, dict) else {}
        method = str(config.get("method", out.iloc[0].get("Method", ""))).upper()
        n_bic = _summary_metric(out, "Persons" if method == "MML" else "N")
        if method == "MML" and (not np.isfinite(n_bic) or n_bic <= 0 or not n_bic.is_integer()):
            n_bic = np.nan
        out["BICSampleSize"] = n_bic
        out["BICSampleSizeBasis"] = (
            "independent_person_marginals" if method == "MML" else "conditional_response_rows_technical"
        )
        if method == "MML" and not np.isfinite(n_bic):
            out["BICSampleSizeBasis"] = "unavailable_person_count"
        out["BICScope"] = (
            "fixed parameter dimension; bounded responses per person; unit weights; regular interior fit required"
            if method == "MML" else "incidental person parameters; ordinary BIC qualification unavailable"
        )
        if method == "MML":
            # Legacy snapshots used response exposure in the BIC penalty.
            # Keep that value for audit; never infer Person count from N.
            if "BICSampleSize" not in summary and "BIC" in summary and "SavedBIC" not in out:
                out["SavedBIC"] = summary["BIC"]
            k = _summary_metric(out, "KParams", _safe_float(config.get("parameter_count")))
            ll = _summary_metric(out, "LogLik")
            out["BIC"] = -2*ll + k*np.log(n_bic) if np.isfinite(n_bic) and n_bic > 0 else np.nan
            n_obs = _summary_metric(out, "N")
            out["BICPerObs"] = out["BIC"]/n_obs if np.isfinite(n_obs) and n_obs > 0 else np.nan
    if not out.empty and (reason := _output_qualification.calibration_inference_hold_reason(
            dict(result, config=dict(config, method=method)))):
        out["InferenceReady"] = False
        out["InferenceReadinessReason"] = reason
        if "PersonMeasureInferenceReady" in out:
            out["PersonMeasureInferenceReady"] = False
    return out


def qualified_person_boundary_table(table, result):
    """Retain score-boundary facts without trusting saved person-use permission."""
    if not isinstance(table, pd.DataFrame) or table.empty:
        return table
    out = table.copy()
    if reason := _output_qualification.calibration_inference_hold_reason(result):
        for column in ("PersonInferenceReady", "PersonMeasureInferenceReady"):
            if column in out:
                out[column] = False
        if "ReportableEstimate" in out:
            out["ReportableEstimate"] = np.nan
        if "EstimateRole" in out:
            out.loc[out.EstimateRole.eq("finite_jmle_estimate"), "EstimateRole"] = "finite_boundary_status_but_fit_not_inference_ready"
        out["InferenceReadinessReason"] = reason
    return out


def build_convergence_summary(opt, config: dict, loglik: float, elapsed_seconds: float | None = None) -> pd.DataFrame:
    """One-row optimizer convergence table for reporting and downloads."""
    if opt is None:
        return pd.DataFrame()
    ll_trace = list(getattr(opt, "ll_trace", []) or [])
    ll_start = float(ll_trace[0]) if ll_trace and np.isfinite(ll_trace[0]) else np.nan
    ll_final = float(loglik) if np.isfinite(loglik) else np.nan
    ll_change = ll_final - ll_start if np.isfinite(ll_final) and np.isfinite(ll_start) else np.nan
    grad_norm = getattr(opt, "gradient_norm", np.nan)
    if not np.isfinite(grad_norm) and getattr(opt, "jac", None) is not None:
        try:
            grad_norm = float(np.linalg.norm(opt.jac))
        except Exception:
            grad_norm = np.nan
    rows = [{
        "Method": config.get("method"),
        "RequestedMmlEngine": config.get("mml_engine_requested"),
        "ResolvedMmlEngine": config.get("mml_engine"),
        "Optimizer": config.get("optimizer"),
        "FacetRegularization": "On" if config.get("facet_regularization_enabled") else "Off",
        "Converged": _output_qualification.optimizer_termination_established({"opt": opt}),
        "Message": str(getattr(opt, "message", "")),
        "Iterations": int(getattr(opt, "nit", 0) or 0),
        "FunctionEvals": int(getattr(opt, "nfev", 0) or 0),
        "WarmupEMIterations": getattr(opt, "em_warmup_iterations", np.nan),
        "DirectIterations": getattr(opt, "direct_iterations", np.nan),
        "FinalLogLik": ll_final,
        "LogLikStart": ll_start,
        "LogLikChange": ll_change,
        "RegularizationPenalty": float(getattr(opt, "penalty_value", np.nan))
        if np.isfinite(getattr(opt, "penalty_value", np.nan)) else np.nan,
        "FinalPenalizedObjective": (
            -ll_final + float(getattr(opt, "penalty_value", 0.0))
            if np.isfinite(ll_final) and np.isfinite(getattr(opt, "penalty_value", 0.0)) else np.nan
        ),
        "GradientNorm": float(grad_norm) if np.isfinite(grad_norm) else np.nan,
        "GradientScope": (
            "structural coordinates and log(population SD)"
            if config.get("mml_engine") == "joint" else (
                "structural coordinates only; population SD excluded"
                if _free_sd_mml_inference_withheld(config) else "optimizer coordinates"
            )
        ),
        "ElapsedSeconds": float(elapsed_seconds) if elapsed_seconds is not None and np.isfinite(elapsed_seconds) else np.nan,
        "AutoAttempted": ", ".join(getattr(opt, "auto_attempted", []) or []),
        "AutoSelected": getattr(opt, "auto_selected", ""),
        "AutoFallbackReason": getattr(opt, "auto_fallback_reason", ""),
    }]
    if _mml_inference_hold_reason(config):
        rows[0].update(QuadratureAccuracyQualified=False, ScientificInferenceReady=False)
    if hasattr(opt, "em_relative_loglik_tolerance"):
        rows[0].update(
            EMRelativeLogLikTolerance=opt.em_relative_loglik_tolerance,
            EMStoppingRelativeLogLikChange=opt.em_stopping_relative_loglik_change,
        )
        if getattr(opt, "em_m_step_options", None):
            controls = opt.em_m_step_options
            rows[0].update(EMMStepOptimizer=controls["method"],
                          EMMStepFunctionTolerance=controls["ftol"],
                          EMMStepGradientTolerance=controls.get("gtol", np.nan),
                          EMMStepMaxIterations=controls["maxiter"])
    if getattr(opt, "m_step_failure", None):
        rows[0].update(MStepFailureReason=opt.m_step_failure["reason"],
                       MStepOptimizerStatus=opt.m_step_failure["status"],
                       MStepOptimizerMessage=opt.m_step_failure["message"])
    if getattr(opt, "em_warmup_failure", None):
        rows[0]["EMWarmupFailure"] = str(opt.em_warmup_failure)
    if config.get("mml_engine") == "joint":
        audit = getattr(opt, "joint_mml", {})
        rows[0].update(JointStationarityPass=bool(audit.get("numerical_pass", False)),
            JointNumericalStatus=audit.get("status", "unavailable"),
            AlgorithmTerminated=all(r.get("success", False) for r in audit.get("runs", []))
                if audit.get("runs") else False,
            QuadratureAccuracyQualified=False, ScientificInferenceReady=False)
    return pd.DataFrame(rows)


def _safe_float(value, default: float = np.nan) -> float:
    """Coerce a scalar to float; return default for missing/non-finite failures."""
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if np.isfinite(out) else default


def _summary_metric(summary: pd.DataFrame | None, name: str, default: float = np.nan) -> float:
    if not isinstance(summary, pd.DataFrame) or summary.empty or name not in summary.columns:
        return default
    return _safe_float(summary.iloc[0].get(name), default=default)


def build_likelihood_information_criteria(result: dict | None) -> pd.DataFrame:
    """Single-run likelihood and information-criteria table for reporting."""
    if not isinstance(result, dict):
        return pd.DataFrame()
    summary = qualified_estimation_summary(result)
    if not isinstance(summary, pd.DataFrame) or summary.empty:
        return pd.DataFrame()
    row = summary.iloc[0]
    config = result.get("config", {}) if isinstance(result.get("config", {}), dict) else {}
    loglik = _summary_metric(summary, "LogLik")
    n_eff = _summary_metric(summary, "N")
    n_bic = _summary_metric(summary, "BICSampleSize")
    k_params = _summary_metric(summary, "KParams")
    if not np.isfinite(k_params):
        k_params = _summary_metric(summary, "Parameters")
    if not np.isfinite(k_params):
        k_params = _safe_float(config.get("parameter_count"))
    deviance = _summary_metric(summary, "Deviance")
    if not np.isfinite(deviance) and np.isfinite(loglik):
        deviance = -2.0 * loglik
    aic = _summary_metric(summary, "AIC")
    if not np.isfinite(aic) and np.isfinite(deviance) and np.isfinite(k_params):
        aic = deviance + 2.0 * k_params
    bic = _summary_metric(summary, "BIC")
    if not np.isfinite(bic) and np.isfinite(deviance) and np.isfinite(k_params) and np.isfinite(n_bic) and n_bic > 0:
        bic = deviance + np.log(n_bic) * k_params
    penalty_value = _summary_metric(summary, "RegularizationPenalty", default=0.0)
    penalized_objective = _summary_metric(summary, "PenalizedObjective")
    if not np.isfinite(penalized_objective) and np.isfinite(deviance):
        penalized_objective = -loglik + penalty_value if np.isfinite(loglik) else np.nan

    def per_obs(value: float) -> float:
        return float(value / n_eff) if np.isfinite(value) and np.isfinite(n_eff) and n_eff > 0 else np.nan

    model = str(row.get("Model", config.get("model", "")))
    method = str(row.get("Method", config.get("method", "")))
    rows = [
        {
            "Metric": "N",
            "Value": n_eff,
            "PreferredDirection": "context",
            "Formula": "effective likelihood rows or total positive weight",
            "Interpretation": "Response exposure for descriptive normalization; not the MML independent sample size.",
        },
        {
            "Metric": "BICSampleSize",
            "Value": n_bic,
            "PreferredDirection": "context",
            "Formula": str(row.get("BICSampleSizeBasis", "unavailable")),
            "Interpretation": str(row.get("BICScope", "unavailable")),
        },
        {
            "Metric": "KParams",
            "Value": k_params,
            "PreferredDirection": "context",
            "Formula": "free estimated parameters counted for AIC/BIC",
            "Interpretation": "Penalty term size; larger models are penalized more strongly.",
        },
        {
            "Metric": "LogLik",
            "Value": loglik,
            "PreferredDirection": "higher",
            "Formula": "maximized log-likelihood",
            "Interpretation": "Higher is better only under the same likelihood definition and data.",
        },
        {
            "Metric": "Deviance",
            "Value": deviance,
            "PreferredDirection": "lower",
            "Formula": "-2 * LogLik",
            "Interpretation": "Lower is better under comparable likelihoods.",
        },
        {
            "Metric": "RegularizationPenalty",
            "Value": penalty_value,
            "PreferredDirection": "context",
            "Formula": "0.5 * sum(((facet effect - mean) / SD)^2)",
            "Interpretation": "Penalty added to the optimized objective; zero for unpenalized runs.",
        },
        {
            "Metric": "PenalizedObjective",
            "Value": penalized_objective,
            "PreferredDirection": "lower",
            "Formula": "-LogLik + RegularizationPenalty",
            "Interpretation": "Objective minimized when facet regularization is active; not an ordinary log-likelihood.",
        },
        {
            "Metric": "AIC",
            "Value": aic,
            "PreferredDirection": "lower",
            "Formula": "2 * KParams - 2 * unpenalized LogLik",
            "Interpretation": "Lower AIC favors the model after a lighter complexity penalty.",
        },
        {
            "Metric": "BIC",
            "Value": bic,
            "PreferredDirection": "lower",
            "Formula": "log(BICSampleSize) * KParams - 2 * unpenalized LogLik",
            "Interpretation": "Technical criterion under the stated sample-size scope; comparison qualification is separate.",
        },
        {
            "Metric": "LogLikPerObs",
            "Value": per_obs(loglik),
            "PreferredDirection": "higher",
            "Formula": "LogLik / N",
            "Interpretation": "Scale-normalized likelihood index for descriptive review.",
        },
        {
            "Metric": "AICPerObs",
            "Value": per_obs(aic),
            "PreferredDirection": "lower",
            "Formula": "AIC / N",
            "Interpretation": "Scale-normalized information criterion for descriptive review.",
        },
        {
            "Metric": "BICPerObs",
            "Value": per_obs(bic),
            "PreferredDirection": "lower",
            "Formula": "BIC / N",
            "Interpretation": "Scale-normalized information criterion for descriptive review.",
        },
    ]
    out = pd.DataFrame(rows)
    out.insert(0, "Method", method)
    out.insert(0, "Model", model)
    return out


def _e_step_posteriors(idx, config, params, quad):
    """E-step: compute posterior weights for each person × quadrature node.

    Returns
    -------
    post_weights : ndarray, shape (n_persons, Q)
        Posterior probabilities  P(θ_q | x_j, δ)  for person j at node q.
    marginal_ll : float
        Current marginal log-likelihood  Σ_j log Σ_q w_q L(x_j|θ_q,δ).

    Reference: Bock, R. D. & Aitkin, M. (1981). Marginal maximum
    likelihood estimation of item parameters: Application of an EM
    algorithm.  Psychometrika, 46(4), 443–459.
    """
    n_persons = _mml_person_count(idx, config)
    log_w = _log_mml_quadrature_weights(quad)           # (Q,)
    log_lik_mat, _ = _compute_mml_person_node_loglik(
        idx, config, params, quad, n_persons=n_persons,
    )

    # log posterior (unnormalized) = log w_q + log L(x_j | θ_q, δ)
    log_joint = log_w[None, :] + log_lik_mat           # (n_persons, Q)
    log_marginal = logsumexp(log_joint, axis=1)        # (n_persons,)
    log_post = log_joint - log_marginal[:, None]       # normalise per person
    post_weights = np.exp(log_post)                    # (n_persons, Q)

    marginal_ll = float(np.sum(log_marginal))
    return post_weights, marginal_ll


def _assemble_mml_gradient_from_contractions(
    par,
    idx,
    config,
    sizes,
    *,
    nll: float,
    eff_w_sum: np.ndarray,
    bar_probs: np.ndarray,
    slope_obs: np.ndarray,
    slope_row_weights: np.ndarray | None,
    log_sigma_grad: float,
    include_log_sigma: bool,
):
    score_k = idx["score_k"]
    k_cat = int(config["n_cat"])
    k_vals = np.arange(k_cat, dtype=float)
    n_steps = max(k_cat - 1, 0)
    facet_signs = config.get("facet_signs", {})

    eta_grad_sum = -slope_obs * (eff_w_sum * score_k - (bar_probs @ k_vals))

    pop = config.get("population_model", {})
    pop_X = np.asarray(pop.get("X", np.zeros((int(config.get("n_person", 0)), 0))), dtype=float)
    if pop.get("enabled") and pop_X.ndim == 2 and pop_X.shape[1] > 0:
        person_eta_grad = np.bincount(idx["person"], weights=eta_grad_sum, minlength=pop_X.shape[0])
        pop_grad_full = pop_X.T @ person_eta_grad
    else:
        pop_grad_full = np.array([], dtype=float)

    facet_full_grads = {}
    for facet in config["facet_names"]:
        sign = facet_signs.get(facet, -1)
        facet_full_grads[facet] = np.bincount(
            idx["facets"][facet],
            weights=eta_grad_sum * sign,
            minlength=len(config["facet_specs"][facet]["levels"]),
        )

    if n_steps > 0:
        thresholds = np.arange(n_steps)
        observed_gt = score_k[:, None] > thresholds[None, :]
        bar_prob_gt = np.cumsum(bar_probs[:, :0:-1], axis=1)[:, ::-1]
        if config["model"] == "RSM":
            step_grad_full = np.sum(eff_w_sum[:, None] * observed_gt - bar_prob_gt, axis=0)
        else:
            step_facet_n = len(config["facet_levels"][config["step_facet"]])
            step_red = _get_step_reduction(idx, step_facet_n)
            row_step = eff_w_sum[:, None] * observed_gt - bar_prob_gt
            if config["model"] == "GPCM":
                row_step = slope_obs[:, None] * row_step
            step_grad_full = _sum_2d_by_group(row_step, step_red)
    else:
        step_grad_full = np.zeros(0, dtype=float)

    if config["model"] == "GPCM" and slope_row_weights is not None:
        slope_idx = idx.get("slope_idx") if idx.get("slope_idx") is not None else idx["step_idx"]
        slope_levels = len(config["facet_levels"][config["slope_facet"]])
        log_slope_grad_full = np.bincount(slope_idx, weights=slope_row_weights, minlength=slope_levels)
    else:
        log_slope_grad_full = np.array([], dtype=float)

    grad_parts = []
    if sizes.get("population", 0) > 0:
        grad_parts.append(pop_grad_full)
    for facet in config["facet_names"]:
        if sizes.get(facet, 0) > 0:
            grad_parts.append(collapse_facet_gradient(facet_full_grads[facet], config["facet_specs"][facet]))
    if n_steps > 0:
        if config["model"] == "RSM":
            grad_parts.append(collapse_sum_zero_gradient(step_grad_full))
        else:
            grad_parts.append(np.vstack([
                collapse_sum_zero_gradient(row)
                for row in step_grad_full
            ]).reshape(-1))
    if config["model"] == "GPCM" and sizes.get("log_slopes", 0) > 0:
        grad_parts.append(collapse_sum_zero_gradient(log_slope_grad_full))

    grad = np.concatenate([g for g in grad_parts if g.size]) if grad_parts else np.array([], dtype=float)
    if grad.size != len(par):
        raise RuntimeError(f"Internal MML gradient size mismatch: expected {len(par)}, got {grad.size}.")
    penalty_value, penalty_grad = facet_regularization_value_grad(par, sizes, config)
    if config.get("facet_regularization", {}).get("enabled"):
        nll += penalty_value
        grad = grad + penalty_grad
    if include_log_sigma:
        grad = np.append(grad, log_sigma_grad)
    return nll, grad


def _mml_joint_or_mstep_3d(
    par,
    params,
    idx,
    config,
    sizes,
    quad,
    *,
    post_weights: np.ndarray | None,
    include_log_sigma: bool = False,
):
    score_k = idx["score_k"]
    n_obs = len(score_k)
    weight_raw = idx.get("weight")
    w_arr = np.ones(n_obs, dtype=float) if weight_raw is None else np.asarray(weight_raw, dtype=float)
    nodes = np.asarray(quad["nodes"], dtype=float)
    k_cat = int(config["n_cat"])
    k_vals = np.arange(k_cat, dtype=float)

    n_persons = (
        _mml_person_count(idx, config)
        if post_weights is None
        else max(_mml_person_count(idx, config), int(post_weights.shape[0]))
    )
    base_eta = compute_base_eta(idx, params, config)
    pop_mu_person = _mml_pop_mu_per_person(params, config, n_persons)
    base_shift = pop_mu_person[idx["person"]] + base_eta

    if config["model"] == "RSM":
        step_cum_obs = np.concatenate([[0.0], np.cumsum(params["steps"])])[None, :]
        slope_obs = np.ones(n_obs, dtype=float)
    elif config["model"] == "PCM":
        step_cum_obs = _build_step_cum_mat(params["steps_mat"])[idx["step_idx"]]
        slope_obs = np.ones(n_obs, dtype=float)
    else:
        step_cum_obs = _build_step_cum_mat(params["steps_mat"])[idx["step_idx"]]
        slope_idx = idx.get("slope_idx") if idx.get("slope_idx") is not None else idx["step_idx"]
        slope_obs = np.asarray(params["slopes"], dtype=float)[slope_idx]

    eta_nq = base_shift[:, None] + nodes[None, :]
    linear = eta_nq[:, :, None] * k_vals[None, None, :] - step_cum_obs[:, None, :]
    log_num = linear if config["model"] != "GPCM" else (slope_obs[:, None, None] * linear)
    log_denom, probs = _log_denom_and_probs(log_num, axis=2, need_probs=True)
    obs_log_num = np.take_along_axis(log_num, score_k[:, None, None], axis=2)[:, :, 0]
    ll_obs_nq = obs_log_num - log_denom

    marginal_ll = np.nan
    if post_weights is None:
        person_red = _get_person_reduction(idx, n_persons)
        w_ll_nq = ll_obs_nq if weight_raw is None else (w_arr[:, None] * ll_obs_nq)
        log_lik_mat = _sum_2d_by_group(w_ll_nq, person_red)
        log_w = _log_mml_quadrature_weights(quad)
        log_joint = log_w[None, :] + log_lik_mat
        log_marginal = logsumexp(log_joint, axis=1)
        post_weights = np.exp(log_joint - log_marginal[:, None])
        marginal_ll = float(np.sum(log_marginal))

    eff_w_nq = w_arr[:, None] * post_weights[idx["person"], :]
    expected_nll = -float(np.sum(eff_w_nq * ll_obs_nq))
    eff_w_sum = np.sum(eff_w_nq, axis=1)
    bar_probs = np.einsum("nq,nqk->nk", eff_w_nq, probs, optimize=True)

    log_sigma_grad = 0.0
    if include_log_sigma:
        expected_nq = probs @ k_vals
        eta_grad_nq = -eff_w_nq * slope_obs[:, None] * (score_k[:, None] - expected_nq)
        log_sigma_grad = float(np.sum(np.sum(eta_grad_nq, axis=0) * nodes))

    slope_row_weights = None
    if config["model"] == "GPCM":
        obs_lin_nq = np.take_along_axis(linear, score_k[:, None, None], axis=2)[:, :, 0]
        exp_lin_nq = np.sum(probs * linear, axis=2)
        slope_row_weights = -slope_obs * np.sum(eff_w_nq * (obs_lin_nq - exp_lin_nq), axis=1)

    expected_nll, grad = _assemble_mml_gradient_from_contractions(
        par,
        idx,
        config,
        sizes,
        nll=expected_nll,
        eff_w_sum=eff_w_sum,
        bar_probs=bar_probs,
        slope_obs=slope_obs,
        slope_row_weights=slope_row_weights,
        log_sigma_grad=log_sigma_grad,
        include_log_sigma=include_log_sigma,
    )
    return marginal_ll, expected_nll, grad


def _m_step_expected_ll_value_grad(
    par, idx, config, sizes, quad, post_weights, *, include_log_sigma=False,
):
    params = expand_params(par, sizes, config)
    score_k = idx["score_k"]
    n_obs = len(score_k)
    nodes = np.asarray(quad["nodes"], dtype=float)
    q_len = len(nodes)
    k_cat = int(config["n_cat"])

    if n_obs > 0 and q_len > 0 and (n_obs * q_len * k_cat <= _MML_MAX_TENSOR_ELEMENTS):
        _, expected_nll, grad = _mml_joint_or_mstep_3d(
            par,
            params,
            idx,
            config,
            sizes,
            quad,
            post_weights=post_weights,
            include_log_sigma=include_log_sigma,
        )
        return expected_nll, grad

    base_eta = compute_base_eta(idx, params, config)
    weight = idx.get("weight")
    w_arr = np.ones_like(score_k, dtype=float) if weight is None else np.asarray(weight, dtype=float)
    k_vals = np.arange(k_cat, dtype=float)
    n_persons = max(_mml_person_count(idx, config), int(post_weights.shape[0]))
    pop_mu_person = _mml_pop_mu_per_person(params, config, n_persons)
    base_shift = pop_mu_person[idx["person"]] + base_eta if n_obs else np.zeros(0, dtype=float)

    if config["model"] == "RSM":
        step_2d = np.concatenate([[0.0], np.cumsum(params["steps"])])
        slope_obs = np.ones(n_obs, dtype=float)
        slope_row_weights = None
    elif config["model"] == "PCM":
        step_2d = _build_step_cum_mat(params["steps_mat"])[idx["step_idx"]] if n_obs else np.zeros((0, k_cat))
        slope_obs = np.ones(n_obs, dtype=float)
        slope_row_weights = None
    else:
        step_2d = _build_step_cum_mat(params["steps_mat"])[idx["step_idx"]] if n_obs else np.zeros((0, k_cat))
        slope_idx = idx.get("slope_idx") if idx.get("slope_idx") is not None else idx["step_idx"]
        slope_obs = np.asarray(params["slopes"], dtype=float)[slope_idx] if n_obs else np.zeros(0, dtype=float)
        slope_row_weights = np.zeros(n_obs, dtype=float)

    nll = 0.0
    log_sigma_grad = 0.0
    eff_w_sum = np.zeros(n_obs, dtype=float)
    bar_probs = np.zeros((n_obs, k_cat), dtype=float)
    row_arange = np.arange(n_obs)

    for q_i, theta_q in enumerate(nodes):
        if n_obs == 0:
            break
        eff_weight = w_arr * post_weights[idx["person"], q_i]
        eff_w_sum += eff_weight
        eta = theta_q + base_shift
        linear = np.outer(eta, k_vals) - step_2d
        log_num = linear if config["model"] != "GPCM" else (slope_obs[:, None] * linear)
        log_denom, probs = _log_denom_and_probs(log_num, axis=1, need_probs=True)
        ll_obs = log_num[row_arange, score_k] - log_denom
        nll -= float(np.sum(eff_weight * ll_obs))
        bar_probs += eff_weight[:, None] * probs
        if include_log_sigma:
            expected = probs @ k_vals
            eta_grad = -eff_weight * slope_obs * (score_k - expected)
            log_sigma_grad += float(theta_q) * float(np.sum(eta_grad))
        if config["model"] == "GPCM":
            obs_linear = linear[row_arange, score_k]
            expected_linear = np.sum(probs * linear, axis=1)
            slope_row_weights -= eff_weight * slope_obs * (obs_linear - expected_linear)

    return _assemble_mml_gradient_from_contractions(
        par,
        idx,
        config,
        sizes,
        nll=nll,
        eff_w_sum=eff_w_sum,
        bar_probs=bar_probs,
        slope_obs=slope_obs,
        slope_row_weights=slope_row_weights,
        log_sigma_grad=log_sigma_grad,
        include_log_sigma=include_log_sigma,
    )


def mfrm_em_mml(start, idx, config, sizes, quad, maxit=200, reltol=1e-6):
    """EM algorithm for MML estimation of MFRM parameters.

    Implements the Bock & Aitkin (1981) EM algorithm:
      E-step — compute posterior  P(θ_q | x_j, δ^(t))  for each
               person j at each quadrature node q.
      M-step — maximise expected complete-data log-likelihood
               Q(δ | δ^(t)) = Σ_j Σ_q r_jq log L(x_j | θ_q, δ)
               w.r.t. item/facet parameters δ using analytical-gradient
               L-BFGS-B.

    Parameters
    ----------
    start : ndarray
        Initial parameter vector (facet + step parameters only; no θ).
    idx : dict
        Observation indices from build_indices().
    config : dict
        Model configuration.
    sizes : OrderedDict
        Parameter block sizes.
    quad : dict
        Gauss-Hermite quadrature nodes and weights.
    maxit : int
        Maximum number of EM iterations.
    reltol : float
        Convergence criterion on relative change in marginal LL.

    Returns
    -------
    result : dict with keys
        'x'        — final parameter vector
        'fun'      — final negative marginal log-likelihood
        'success'  — whether convergence was achieved
        'nit'      — number of EM iterations
        'nfev'     — total function evaluations across all M-steps
        'll_trace' — marginal LL at each iteration

    References
    ----------
    Bock, R. D. & Aitkin, M. (1981). Marginal maximum likelihood
        estimation of item parameters: Application of an EM algorithm.
        Psychometrika, 46(4), 443–459.
    Bock, R. D. & Mislevy, R. J. (1982). Adaptive EAP estimation of
        ability in a microcomputer environment. Applied Psychological
        Measurement, 6(4), 431–444.
    """
    par = np.array(start, dtype=float, copy=True)
    params = expand_params(par, sizes, config)

    # Optional free latent-population-SD estimation (Bock-Aitkin variance
    # M-step). Opt-in via config['estimate_population_sd']; when off, the quad
    # below stays the caller's fixed-SD object and this engine is unchanged.
    estimate_sd = bool(config.get("estimate_population_sd", False))
    if estimate_sd:
        quad = dict(quad)  # local mutable copy; never mutate the caller's quad
        n_quad_points = len(quad["nodes"])
        sd_lo, sd_hi = config.get("population_sd_bounds", (0.05, 10.0))
        sigma = float(quad.get("sd", get_population_prior_sd(config)))
        persons_with_obs = np.unique(np.asarray(idx["person"]))
        n_eff = max(int(persons_with_obs.size), 1)
        sigma_trace = [sigma]

    prev_ll = -np.inf
    converged = False
    ll_trace = []
    total_nfev = 0
    m_step_failure = None
    rel_change = np.nan
    m_step_options = None

    for it in range(maxit):
        # ── E-step ──────────────────────────────────────────────
        post_weights, marginal_ll = _e_step_posteriors(idx, config, params, quad)
        ll_trace.append(marginal_ll)

        # ── Convergence check (relative change in marginal LL) ─
        if it > 0:
            rel_change = abs(marginal_ll - prev_ll) / (abs(prev_ll) + 1e-10)
            if rel_change < reltol:
                converged = True
                break
        prev_ll = marginal_ll

        # ── M-step (analytical gradient; constrained for GPCM) ─
        method, options = optimizer_method_and_options(
            sizes,
            config,
            maxiter=50,
            gtol=1e-5,
            ftol=1e-7,
        )
        m_step_options = {"method": method, **options}
        opt = minimize(
            _m_step_expected_ll_value_grad,
            par,
            args=(idx, config, sizes, quad, post_weights),
            jac=True,
            method=method,
            bounds=build_optimizer_bounds(sizes, config),
            constraints=build_optimizer_constraints(sizes, config),
            options=options,
        )
        total_nfev += opt.nfev
        candidate = np.asarray(opt.x, dtype=float)
        failure_reason = ""
        if not opt.success:
            failure_reason = "stat.mml_em.m_step_optimizer_failed"
        elif candidate.shape != par.shape or not np.all(np.isfinite(candidate)):
            failure_reason = "stat.mml_em.m_step_nonfinite"
        else:
            bounds = build_optimizer_bounds(sizes, config) or []
            feasible = all((lo is None or value >= lo - 1e-8) and
                           (hi is None or value <= hi + 1e-8)
                           for value, (lo, hi) in zip(candidate, bounds))
            for constraint in build_optimizer_constraints(sizes, config):
                value = constraint.A @ candidate
                feasible = feasible and bool(np.all(value >= constraint.lb - 1e-8)
                                             and np.all(value <= constraint.ub + 1e-8))
            if not feasible:
                failure_reason = "stat.mml_em.m_step_infeasible"
            else:
                old_q, _ = _m_step_expected_ll_value_grad(par, idx, config, sizes, quad, post_weights)
                new_q, new_grad = _m_step_expected_ll_value_grad(candidate, idx, config, sizes, quad, post_weights)
                total_nfev += 2
                if not (np.isfinite(old_q) and np.isfinite(new_q) and np.all(np.isfinite(new_grad))):
                    failure_reason = "stat.mml_em.m_step_nonfinite"
                elif new_q > old_q + 1e-10 * (1.0 + abs(old_q)):
                    failure_reason = "stat.mml_em.m_step_objective_increased"
        if failure_reason:
            m_step_failure = dict(reason=failure_reason, iteration=it + 1,
                                  status=getattr(opt, "status", None), message=str(opt.message))
            break
        par = candidate
        params = expand_params(par, sizes, config)

        # ── Population-variance M-step (Bock-Aitkin) ───────────
        # sigma^2 = (1/n_eff) Σ_j Σ_q post[j,q] · node_q^2, where node_q is the
        # quadrature deviation (theta_jq − mu_j). Rebuild the grid at the new
        # sigma so the next E-step integrates over N(mu_j, sigma^2).
        # Regridding does not guarantee monotonicity of the finite-Q objective.
        # Empty-person rows are excluded.
        if estimate_sd:
            nodes_sq = (quad["nodes"] ** 2)[None, :]
            sigma2_new = float(np.sum(post_weights[persons_with_obs] * nodes_sq) / n_eff)
            sigma = float(np.clip(np.sqrt(max(sigma2_new, sd_lo ** 2)), sd_lo, sd_hi))
            quad = gauss_hermite_normal(n_quad_points, sd=sigma)
            sigma_trace.append(sigma)

    # Final E-step to get the converged marginal LL
    _, final_ll = _e_step_posteriors(idx, config, params, quad)
    ll_trace.append(final_ll)

    # Use a serializable result so reports and saved runs can reuse the fit.
    from scipy.optimize import OptimizeResult

    result = OptimizeResult()
    result.x = par
    penalty_value, _ = facet_regularization_value_grad(par, sizes, config)
    result.raw_fun = -final_ll
    result.penalty_value = penalty_value
    result.fun = -final_ll + penalty_value
    result.success = converged
    result.nit = len(ll_trace) - 1
    result.nfev = total_nfev
    result.ll_trace = ll_trace
    result.estimated_population_sd = float(sigma) if estimate_sd else None
    result.sigma_trace = sigma_trace if estimate_sd else None
    result.final_quad_sd = float(quad.get("sd")) if estimate_sd else None
    result.m_step_failure = m_step_failure
    result.em_relative_loglik_tolerance = float(reltol)
    # The final trace entry repeats the stopping point; its zero difference is not the tested change.
    result.em_stopping_relative_loglik_change = float(rel_change) if converged else np.nan
    result.em_m_step_options = m_step_options
    result.status = 2 if m_step_failure else (0 if converged else 1)
    result.message = "EM relative log-likelihood convergence reached." if converged else "EM reached max iterations before relative log-likelihood convergence."
    if m_step_failure:
        result.message = (f"EM stopped at M-step {m_step_failure['iteration']}: "
                          f"{m_step_failure['reason']}. {m_step_failure['message']}")
    return result


def _weighted_quantile_from_grid(nodes: np.ndarray, weights: np.ndarray, probs=(0.05, 0.5, 0.95)) -> list[float]:
    order = np.argsort(nodes)
    x = np.asarray(nodes, dtype=float)[order]
    w = np.asarray(weights, dtype=float)[order]
    total = float(np.sum(w))
    if total <= 0 or not np.isfinite(total):
        return [np.nan for _ in probs]
    cdf = np.cumsum(w / total)
    return [float(x[min(np.searchsorted(cdf, p, side="left"), len(x) - 1)]) for p in probs]


def compute_person_posterior_outputs(
    idx,
    config,
    params,
    quad,
    person_levels=None,
    n_plausible_values=0,
    seed=20260411,
):
    """Person posterior summaries conditional on fitted MML calibration.

    The quadrature represents the configured normal population distribution.
    MAPGrid maximizes posterior density over its nodes, not quadrature mass.
    """
    if config.get("method") != "MML":
        return {
            "available": False,
            "reason": "Posterior scoring and plausible values are defined for MML runs only.",
            "scores": pd.DataFrame(),
            "plausible_values": pd.DataFrame(),
        }

    n = len(idx["score_k"])
    n_person = int(config.get("n_person", 0))
    if n_person <= 0:
        return {
            "available": False,
            "reason": "No persons available for posterior scoring.",
            "scores": pd.DataFrame(),
            "plausible_values": pd.DataFrame(),
        }

    if person_levels is None:
        person_levels = [str(i) for i in range(n_person)]

    log_w = _log_mml_quadrature_weights(quad)
    rng = np.random.default_rng(int(seed))
    n_pv = max(0, int(n_plausible_values or 0))
    score_rows: list[dict] = []
    pv_rows: list[dict] = []
    nodes = np.asarray(quad["nodes"], dtype=float)
    pop_mu_person = _mml_pop_mu_per_person(params, config, n_person)
    ll_mat, person_red = _compute_mml_person_node_loglik(
        idx, config, params, quad, n_persons=n_person,
    )
    theta_grid = pop_mu_person[:, None] + nodes[None, :]
    log_post = log_w[None, :] + ll_mat
    log_post = log_post - logsumexp(log_post, axis=1)[:, None]
    post_w_mat = np.exp(log_post)
    eap_arr = np.sum(theta_grid * post_w_mat, axis=1)
    sd_arr = np.sqrt(np.sum((theta_grid - eap_arr[:, None]) ** 2 * post_w_mat, axis=1))
    prior_sd = get_population_prior_sd(config)
    log_density = ll_mat - 0.5 * (nodes[None, :] / prior_sd) ** 2
    map_idx = np.argmax(log_density, axis=1)
    obs_rows_arr = person_red["counts"]
    weight_col_arr = idx.get("weight")
    if weight_col_arr is not None:
        obs_weight_arr = np.bincount(
            idx["person"],
            weights=np.asarray(weight_col_arr, dtype=float),
            minlength=n_person,
        )
    else:
        obs_weight_arr = obs_rows_arr.astype(float)

    for p_idx in range(n_person):
        theta_nodes = theta_grid[p_idx]
        post_w = post_w_mat[p_idx]
        q05, q50, q95 = _weighted_quantile_from_grid(theta_nodes, post_w, probs=(0.05, 0.5, 0.95))
        person_id = person_levels[p_idx] if p_idx < len(person_levels) else str(p_idx)
        score_rows.append({
            "Person": person_id,
            "PosteriorMean": float(eap_arr[p_idx]),
            "PosteriorSD": float(sd_arr[p_idx]),
            "PosteriorMedian": q50,
            "PosteriorQ05": q05,
            "PosteriorQ95": q95,
            "PosteriorMAPGrid": float(theta_nodes[int(map_idx[p_idx])]),
            "PopulationMean": float(pop_mu_person[p_idx]),
            "ObservedRows": int(obs_rows_arr[p_idx]),
            "ObservedWeight": float(obs_weight_arr[p_idx]),
        })
        if n_pv > 0:
            draws = rng.choice(theta_nodes, size=n_pv, replace=True, p=post_w / np.sum(post_w))
            pv_rows.append({
                "Person": person_id,
                **{f"PV_{i + 1}": float(draws[i]) for i in range(n_pv)},
            })

    return {
        "available": True,
        "reason": "Posterior scores computed from the MML quadrature posterior.",
        "scores": pd.DataFrame(score_rows),
        "plausible_values": pd.DataFrame(pv_rows),
        "n_plausible_values": n_pv,
        "seed": int(seed),
        "interpretation": (
            "PosteriorMean is the EAP person score. All posterior summaries and "
            "plausible values condition on fitted calibration and population parameters; "
            "calibration uncertainty is not propagated. PosteriorMAPGrid maximizes "
            "posterior density over the quadrature nodes, not over the continuous "
            "ability scale. Plausible values are draws from the conditional quadrature "
            "posterior for group-level uncertainty propagation, not individual ranking."
        ),
    }


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


def audit_mfrm_anchors(
    prep,
    anchor_df=None,
    group_anchor_df=None,
    dummy_facets=None,
    min_common_anchors=2,
    min_obs_per_element=2,
    min_obs_per_category=1,
):
    """Audit anchor inputs before they are used for scale constraints."""
    all_facets = ["Person"] + prep["facet_names"]
    dummy_facets = set(dummy_facets or [])
    levels_map = {facet: set(map(str, prep["levels"].get(facet, []))) for facet in all_facets}

    data = prep.get("data", pd.DataFrame())
    if not data.empty:
        element_counts = {}
        for facet in all_facets:
            if facet in data.columns:
                element_counts[facet] = data.groupby(facet, observed=False)["Score"].size().to_dict()
            else:
                element_counts[facet] = {}
        cat_counts = data["Score"].value_counts().to_dict() if "Score" in data.columns else {}
    else:
        element_counts = {facet: {} for facet in all_facets}
        cat_counts = {}

    issues: list[dict] = []
    anchor_norm = normalize_anchor_df(anchor_df) if anchor_df is not None else pd.DataFrame()
    group_norm = normalize_group_anchor_df(group_anchor_df) if group_anchor_df is not None else pd.DataFrame()

    if anchor_df is not None and not getattr(anchor_df, "empty", True) and anchor_norm.empty:
        issues.append({
            "Severity": "error",
            "Type": "anchor_format",
            "Facet": "",
            "Level": "",
            "Message": "Anchor table could not be parsed. Required columns: Facet, Level, Anchor.",
            "Action": "Rename columns or paste a CSV with Facet,Level,Anchor.",
        })
    if group_anchor_df is not None and not getattr(group_anchor_df, "empty", True) and group_norm.empty:
        issues.append({
            "Severity": "error",
            "Type": "group_anchor_format",
            "Facet": "",
            "Level": "",
            "Message": "Group-anchor table could not be parsed. Required columns: Facet, Level, Group, GroupValue.",
            "Action": "Rename columns or paste a CSV with Facet,Level,Group,GroupValue.",
        })

    valid_anchor_rows = []
    if not anchor_norm.empty:
        dup_mask = anchor_norm.duplicated(subset=["Facet", "Level"], keep="last")
        for row in anchor_norm[dup_mask].itertuples(index=False):
            issues.append({
                "Severity": "warning",
                "Type": "duplicate_anchor",
                "Facet": row.Facet,
                "Level": row.Level,
                "Message": "Duplicate anchor row detected; the last value will be used.",
                "Action": "Keep only one anchor value per Facet x Level for reproducibility.",
            })
        anchor_work = anchor_norm.drop_duplicates(subset=["Facet", "Level"], keep="last")
        for row in anchor_work.itertuples(index=False):
            facet = str(row.Facet)
            level = str(row.Level)
            if facet not in all_facets:
                issues.append({
                    "Severity": "warning",
                    "Type": "unknown_facet",
                    "Facet": facet,
                    "Level": level,
                    "Message": "Anchor references a facet that is not in the current model; it will be ignored.",
                    "Action": "Check spelling or include that facet in the model.",
                })
                continue
            if level not in levels_map[facet]:
                issues.append({
                    "Severity": "warning",
                    "Type": "unknown_level",
                    "Facet": facet,
                    "Level": level,
                    "Message": "Anchor references a level not observed in the current data; it will be ignored.",
                    "Action": "Check spelling or use a dataset containing this level.",
                })
                continue
            if not np.isfinite(row.Anchor):
                issues.append({
                    "Severity": "warning",
                    "Type": "invalid_anchor_value",
                    "Facet": facet,
                    "Level": level,
                    "Message": "Anchor value is missing or non-finite; it will be ignored.",
                    "Action": "Provide a finite logit value.",
                })
                continue
            n_level = int(element_counts.get(facet, {}).get(level, 0))
            if n_level < int(min_obs_per_element):
                issues.append({
                    "Severity": "warning",
                    "Type": "low_anchor_observations",
                    "Facet": facet,
                    "Level": level,
                    "Message": f"Anchored level has only {n_level} observation(s).",
                    "Action": "Treat this anchor cautiously or add more linked observations.",
                })
            valid_anchor_rows.append({"Facet": facet, "Level": level, "Anchor": float(row.Anchor)})

    valid_group_rows = []
    if not group_norm.empty:
        dup_mask = group_norm.duplicated(subset=["Facet", "Level"], keep="last")
        for row in group_norm[dup_mask].itertuples(index=False):
            issues.append({
                "Severity": "warning",
                "Type": "duplicate_group_anchor",
                "Facet": row.Facet,
                "Level": row.Level,
                "Message": "Duplicate group-anchor row detected; the last value will be used.",
                "Action": "Keep only one group assignment per Facet x Level.",
            })
        group_work = group_norm.drop_duplicates(subset=["Facet", "Level"], keep="last")
        for row in group_work.itertuples(index=False):
            facet = str(row.Facet)
            level = str(row.Level)
            group = str(row.Group)
            if facet not in all_facets:
                issues.append({
                    "Severity": "warning",
                    "Type": "unknown_group_facet",
                    "Facet": facet,
                    "Level": level,
                    "Message": "Group anchor references a facet that is not in the current model; it will be ignored.",
                    "Action": "Check spelling or include that facet in the model.",
                })
                continue
            if level not in levels_map[facet]:
                issues.append({
                    "Severity": "warning",
                    "Type": "unknown_group_level",
                    "Facet": facet,
                    "Level": level,
                    "Message": "Group anchor references a level not observed in the current data; it will be ignored.",
                    "Action": "Check spelling or use a dataset containing this level.",
                })
                continue
            if not np.isfinite(row.GroupValue):
                issues.append({
                    "Severity": "warning",
                    "Type": "invalid_group_anchor_value",
                    "Facet": facet,
                    "Level": level,
                    "Message": "Group-anchor value is missing or non-finite; it will be ignored.",
                    "Action": "Provide a finite group-anchor logit value; do not rely on an implicit zero.",
                })
                continue
            valid_group_rows.append({
                "Facet": facet,
                "Level": level,
                "Group": group,
                "GroupValue": float(row.GroupValue),
            })

    valid_anchor_df = pd.DataFrame(valid_anchor_rows, columns=["Facet", "Level", "Anchor"])
    valid_group_df = pd.DataFrame(valid_group_rows, columns=["Facet", "Level", "Group", "GroupValue"])

    summary_rows = []
    for facet in all_facets:
        n_levels = len(levels_map[facet])
        fixed_levels = (
            set(valid_anchor_df.loc[valid_anchor_df["Facet"] == facet, "Level"].astype(str))
            if not valid_anchor_df.empty else set()
        )
        group_levels = (
            set(valid_group_df.loc[valid_group_df["Facet"] == facet, "Level"].astype(str))
            if not valid_group_df.empty else set()
        )
        n_anch = len(fixed_levels)
        n_group_levels = len(group_levels)
        n_groups = int(valid_group_df.loc[valid_group_df["Facet"] == facet, "Group"].nunique()) if not valid_group_df.empty else 0
        anchored_levels = set(levels_map[facet]) if facet in dummy_facets else fixed_levels | group_levels
        anchored_total = len(anchored_levels)
        anchor_share = float(anchored_total / n_levels) if n_levels else np.nan
        unanchored_levels = max(0, n_levels - anchored_total)
        if facet in dummy_facets:
            status = "Dummy"
            action = "Facet fixed at 0 by dummy-facet setting."
        elif anchored_total == 0:
            status = "Unanchored"
            action = "No fixed/group anchors supplied. Joint calibration can link groups through shared observations; inspect design identification separately."
        elif anchored_total < int(min_common_anchors):
            status = "Below count threshold"
            action = f"Review the configured count threshold ({min_common_anchors}) and reference quality; no count alone establishes valid linking."
            issues.append({
                "Severity": "warning",
                "Type": "few_common_anchors",
                "Facet": facet,
                "Level": "",
                "Message": f"Only {anchored_total} anchored/common level(s) available for this facet.",
                "Action": action,
            })
        else:
            status = "Count threshold met"
            action = "The configured count threshold is met. Content match, identification, fit and empirical reference stability still require review."
        summary_rows.append({
            "Facet": facet,
            "Levels": n_levels,
            "FixedAnchors": n_anch,
            "GroupAnchoredLevels": n_group_levels,
            "AnchoredLevelsTotal": anchored_total,
            "UnanchoredLevels": unanchored_levels,
            "AnchorShare": anchor_share,
            "Groups": n_groups,
            "DummyFacet": facet in dummy_facets,
            "Status": status,
            "Action": action,
            "CoverageInterpretation": (
                "Descriptive coverage only: no universal anchor percentage is sufficient. "
                "Review common-anchor count, content match, observations, drift, and connectedness together."
            ),
        })

    for cat, cnt in sorted(cat_counts.items()):
        if int(cnt) < int(min_obs_per_category):
            issues.append({
                "Severity": "warning",
                "Type": "low_category_count",
                "Facet": "Score",
                "Level": str(cat),
                "Message": f"Score category {cat} has only {int(cnt)} observation(s).",
                "Action": "Threshold and anchor-based linking diagnostics may be unstable for sparse categories.",
            })

    issues_df = pd.DataFrame(issues, columns=["Severity", "Type", "Facet", "Level", "Message", "Action"])
    summary_df = pd.DataFrame(summary_rows)
    n_error = int((issues_df["Severity"] == "error").sum()) if not issues_df.empty else 0
    n_warning = int((issues_df["Severity"] == "warning").sum()) if not issues_df.empty else 0
    overall_status = "error" if n_error else ("warning" if n_warning else "ok")
    linking_review = {
        "overall_status": overall_status,
        "message": (
            "Anchor audit found error(s). Fix anchor inputs before estimating."
            if n_error else
            "Anchor audit found warning(s). Review before making linking claims."
            if n_warning else
            "Anchor audit found no obvious anchor-input issues."
        ),
        "summary": summary_df,
        "issues": issues_df,
        "valid_anchors": valid_anchor_df,
        "valid_group_anchors": valid_group_df,
        "settings": {
            "min_common_anchors": int(min_common_anchors),
            "min_obs_per_element": int(min_obs_per_element),
            "min_obs_per_category": int(min_obs_per_category),
        },
    }
    return linking_review


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


def _combined_measure_table_for_linking(result: dict) -> pd.DataFrame:
    facets = result.get("facets", {})
    frames = []
    person_tbl = facets.get("person", pd.DataFrame()) if isinstance(facets, dict) else pd.DataFrame()
    if isinstance(person_tbl, pd.DataFrame) and {"Person", "Estimate"}.issubset(person_tbl.columns):
        frames.append(pd.DataFrame({
            "Facet": "Person",
            "Level": person_tbl["Person"].astype(str),
            "Estimate": pd.to_numeric(person_tbl["Estimate"], errors="coerce"),
        }))
    other_tbl = facets.get("others", pd.DataFrame()) if isinstance(facets, dict) else pd.DataFrame()
    if isinstance(other_tbl, pd.DataFrame) and {"Facet", "Level", "Estimate"}.issubset(other_tbl.columns):
        frames.append(pd.DataFrame({
            "Facet": other_tbl["Facet"].astype(str),
            "Level": other_tbl["Level"].astype(str),
            "Estimate": pd.to_numeric(other_tbl["Estimate"], errors="coerce"),
        }))
    if not frames:
        return pd.DataFrame(columns=["Facet", "Level", "Estimate"])
    return pd.concat(frames, ignore_index=True)


def build_anchor_drift_review(result: dict, tolerance: float | None = None, warning_abs: float = 0.05) -> dict:
    """Post-fit anchor consistency screen.

    Because this app uses anchor inputs as hard constraints, this is not an
    empirical unconstrained drift estimate. It verifies that the fitted scale
    respects the supplied baseline anchors and group-anchor target means.
    """
    config = result.get("config", {})
    audit = config.get("anchor_audit", {})
    if not isinstance(audit, dict):
        return {
            "available": False,
            "reason": "No anchor audit is available.",
            "summary": pd.DataFrame(),
            "anchor_drift": pd.DataFrame(),
            "group_drift": pd.DataFrame(),
        }
    valid_anchors = audit.get("valid_anchors", pd.DataFrame())
    valid_groups = audit.get("valid_group_anchors", pd.DataFrame())
    has_anchors = isinstance(valid_anchors, pd.DataFrame) and not valid_anchors.empty
    has_groups = isinstance(valid_groups, pd.DataFrame) and not valid_groups.empty
    if not has_anchors and not has_groups:
        return {
            "available": False,
            "reason": "No fixed anchors or group anchors were supplied for this run.",
            "summary": pd.DataFrame(),
            "anchor_drift": pd.DataFrame(),
            "group_drift": pd.DataFrame(),
            "scope": "post_fit_constraint_check_only",
        }

    tol = float(tolerance) if tolerance is not None else max(1e-6, 10.0 * float(config.get("reltol") or 1e-6))
    measures = _combined_measure_table_for_linking(result)
    anchor_tbl = pd.DataFrame()
    if has_anchors:
        anchor_tbl = valid_anchors.copy()
        anchor_tbl["Facet"] = anchor_tbl["Facet"].astype(str)
        anchor_tbl["Level"] = anchor_tbl["Level"].astype(str)
        anchor_tbl["Anchor"] = pd.to_numeric(anchor_tbl["Anchor"], errors="coerce")
        anchor_tbl = anchor_tbl.merge(measures, on=["Facet", "Level"], how="left")
        anchor_tbl["Drift"] = anchor_tbl["Estimate"] - anchor_tbl["Anchor"]
        anchor_tbl["AbsDrift"] = anchor_tbl["Drift"].abs()
        anchor_tbl["Tolerance"] = tol
        anchor_tbl["Status"] = "Constraint ok"
        anchor_tbl.loc[anchor_tbl["Estimate"].isna(), "Status"] = "Missing estimate"
        anchor_tbl.loc[
            anchor_tbl["Estimate"].notna() & (anchor_tbl["AbsDrift"] > tol) & (anchor_tbl["AbsDrift"] <= float(warning_abs)),
            "Status",
        ] = "Small numerical drift"
        anchor_tbl.loc[
            anchor_tbl["Estimate"].notna() & (anchor_tbl["AbsDrift"] > float(warning_abs)),
            "Status",
        ] = "Review"
        anchor_tbl["Scope"] = "Constrained anchor consistency check"
        anchor_tbl["Action"] = np.where(
            anchor_tbl["Status"].eq("Constraint ok"),
            "No action needed; this verifies the hard anchor constraint.",
            "Check anchor spelling, duplicate inputs, and optimizer convergence before making linking claims.",
        )

    group_tbl = pd.DataFrame()
    if has_groups:
        group_rows = valid_groups.copy()
        group_rows["Facet"] = group_rows["Facet"].astype(str)
        group_rows["Level"] = group_rows["Level"].astype(str)
        group_rows["Group"] = group_rows["Group"].astype(str)
        group_rows["GroupValue"] = pd.to_numeric(group_rows["GroupValue"], errors="coerce").fillna(0.0)
        group_rows = group_rows.merge(measures, on=["Facet", "Level"], how="left")
        group_tbl = (
            group_rows.groupby(["Facet", "Group"], dropna=False)
            .agg(
                GroupValue=("GroupValue", "first"),
                Levels=("Level", "nunique"),
                EstimatedMean=("Estimate", "mean"),
                EstimatedSD=("Estimate", "std"),
                MissingEstimates=("Estimate", lambda x: int(pd.to_numeric(x, errors="coerce").isna().sum())),
            )
            .reset_index()
        )
        group_tbl["MeanDrift"] = group_tbl["EstimatedMean"] - group_tbl["GroupValue"]
        group_tbl["AbsMeanDrift"] = group_tbl["MeanDrift"].abs()
        group_tbl["Tolerance"] = tol
        group_tbl["Status"] = "Group target ok"
        group_tbl.loc[group_tbl["MissingEstimates"] > 0, "Status"] = "Missing estimate"
        group_tbl.loc[
            group_tbl["EstimatedMean"].notna() & (group_tbl["AbsMeanDrift"] > tol) & (group_tbl["AbsMeanDrift"] <= float(warning_abs)),
            "Status",
        ] = "Small numerical drift"
        group_tbl.loc[
            group_tbl["EstimatedMean"].notna() & (group_tbl["AbsMeanDrift"] > float(warning_abs)),
            "Status",
        ] = "Review"
        group_tbl["Scope"] = "Group-anchor mean consistency check"
        group_tbl["Action"] = np.where(
            group_tbl["Status"].eq("Group target ok"),
            "No action needed; the group target mean is respected.",
            "Check group-anchor values and convergence before making group-linking claims.",
        )

    max_anchor = float(pd.to_numeric(anchor_tbl.get("AbsDrift", pd.Series(dtype=float)), errors="coerce").max()) if not anchor_tbl.empty else np.nan
    max_group = float(pd.to_numeric(group_tbl.get("AbsMeanDrift", pd.Series(dtype=float)), errors="coerce").max()) if not group_tbl.empty else np.nan
    review_count = 0
    if not anchor_tbl.empty:
        review_count += int(anchor_tbl["Status"].isin(["Review", "Missing estimate"]).sum())
    if not group_tbl.empty:
        review_count += int(group_tbl["Status"].isin(["Review", "Missing estimate"]).sum())
    summary = pd.DataFrame([{
        "Available": True,
        "ExplicitAnchors": int(len(valid_anchors)) if has_anchors else 0,
        "GroupAnchorLevels": int(len(valid_groups)) if has_groups else 0,
        "MaxAbsAnchorDrift": max_anchor,
        "MaxAbsGroupMeanDrift": max_group,
        "Tolerance": tol,
        "ReviewCount": review_count,
        "Status": "Ready" if review_count == 0 else "Review",
        "Scope": "post_fit_constraint_check_only",
        "Action": (
            "Anchor constraints are internally consistent; this does not prove empirical invariance."
            if review_count == 0 else
            "Resolve anchor consistency warnings before using this run for linking."
        ),
    }])
    return {
        "available": True,
        "reason": "Anchor drift review is a post-fit constraint consistency check for the supplied anchors.",
        "summary": summary,
        "anchor_drift": anchor_tbl,
        "group_drift": group_tbl,
        "scope": "post_fit_constraint_check_only",
        "empirical_unconstrained_drift_available": False,
    }


def build_equating_chain_summary(result: dict) -> dict:
    """Single-run linking-chain summary from anchor and group-anchor inputs."""
    config = result.get("config", {})
    audit = config.get("anchor_audit", {})
    if not isinstance(audit, dict):
        return {
            "available": False,
            "reason": "No anchor audit is available.",
            "summary": pd.DataFrame(),
            "edges": pd.DataFrame(),
        }
    audit_summary = audit.get("summary", pd.DataFrame())
    if not isinstance(audit_summary, pd.DataFrame) or audit_summary.empty:
        return {
            "available": False,
            "reason": "Anchor summary is unavailable.",
            "summary": pd.DataFrame(),
            "edges": pd.DataFrame(),
        }
    min_common = int(audit.get("settings", {}).get("min_common_anchors", 2))
    rows = []
    edges = []
    for row in audit_summary.itertuples(index=False):
        facet = str(getattr(row, "Facet", ""))
        levels = int(getattr(row, "Levels", 0) or 0)
        fixed = int(getattr(row, "FixedAnchors", 0) or 0)
        group_levels = int(getattr(row, "GroupAnchoredLevels", 0) or 0)
        groups = int(getattr(row, "Groups", 0) or 0)
        dummy = bool(getattr(row, "DummyFacet", False))
        common = int(
            getattr(
                row,
                "AnchoredLevelsTotal",
                fixed + group_levels + (levels if dummy else 0),
            ) or 0
        )
        coverage = (common / levels) if levels > 0 else np.nan
        if dummy:
            link_type = "Dummy facet"
        elif fixed and group_levels:
            link_type = "Fixed + group anchors"
        elif fixed:
            link_type = "Fixed anchors"
        elif group_levels:
            link_type = "Group anchors"
        else:
            link_type = "None"
        if dummy:
            strength = "Dummy constraint"
        elif common >= min_common and common > 0:
            strength = "Count threshold met"
        elif common > 0:
            strength = "Below count threshold"
        else:
            strength = "Unanchored"
        action = (
            "No fixed/group anchors supplied for this facet. Joint-calibration connectivity is assessed separately."
            if strength == "Unanchored" else
            "Fixing a facet to zero is a modeling constraint, not evidence of a stable reference baseline."
            if dummy else
            "Review reference quality, content match, identification and empirical stability; the count threshold does not establish valid linking."
        )
        rows.append({
            "Facet": facet,
            "Levels": levels,
            "LinkType": link_type,
            "CommonLinkedLevels": common,
            "FixedAnchors": fixed,
            "GroupAnchoredLevels": group_levels,
            "Groups": groups,
            "Coverage": coverage,
            "MinimumCommonAnchors": min_common,
            "Strength": strength,
            "Action": action,
        })
        if common > 0:
            edges.append({
                "From": "Current run",
                "To": f"Anchor baseline via {facet}",
                "Facet": facet,
                "LinkType": link_type,
                "CommonLinkedLevels": common,
                "Strength": strength,
                "Action": action,
            })

    summary = pd.DataFrame(rows)
    edges_df = pd.DataFrame(edges, columns=[
        "From", "To", "Facet", "LinkType", "CommonLinkedLevels", "Strength", "Action",
    ])
    available = not edges_df.empty
    return {
        "available": available,
        "reason": (
            "Single-run anchor chain summary built from fixed anchors and group anchors."
            if available else "No fixed/group anchors or dummy constraints were supplied. Shared observations for joint calibration are assessed separately."
        ),
        "summary": summary,
        "edges": edges_df,
        "scope": "single_run_anchor_chain_summary",
        "full_multi_run_chain_available": False,
        "interpretation": (
            "This describes supplied anchor constraints and their counts. It does not "
            "establish empirical reference stability, validate scale linking, assess "
            "joint-calibration connectivity or reconstruct a historical equating chain."
        ),
    }


def build_anchor_equating_workflow_plan(result: dict) -> pd.DataFrame:
    """Operational anchor/equating checklist for the current Streamlit run."""
    config = result.get("config", {}) if isinstance(result, dict) else {}
    audit = config.get("anchor_audit", {}) if isinstance(config, dict) else {}
    drift = result.get("anchor_drift", {}) if isinstance(result, dict) else {}
    chain = result.get("equating_chain", {}) if isinstance(result, dict) else {}
    rows: list[dict] = []

    if not isinstance(audit, dict):
        return pd.DataFrame([{
            "Step": 1,
            "Check": "Anchor inputs",
            "Status": "Missing",
            "Evidence": "Anchor audit was not recorded for this run.",
            "Action": "Run the model again with current anchor settings before making linking claims.",
        }])

    audit_status = str(audit.get("overall_status", "unknown"))
    audit_summary = audit.get("summary", pd.DataFrame())
    audit_issues = audit.get("issues", pd.DataFrame())
    n_issues = len(audit_issues) if isinstance(audit_issues, pd.DataFrame) else 0
    rows.append({
        "Step": 1,
        "Check": "Anchor input validation",
        "Status": "Ready" if audit_status == "ok" else "Review",
        "Evidence": f"{audit.get('message', 'No audit message.')} Issues: {n_issues}.",
        "Action": "Fix format, spelling, duplicate, sparse, or missing-anchor issues before claiming a linked scale.",
    })

    linked = weak = unanchored = 0
    if isinstance(audit_summary, pd.DataFrame) and not audit_summary.empty:
        chain_summary = chain.get("summary", pd.DataFrame()) if isinstance(chain, dict) else pd.DataFrame()
        if isinstance(chain_summary, pd.DataFrame) and not chain_summary.empty and "Strength" in chain_summary.columns:
            linked = int(chain_summary["Strength"].astype(str).isin(["Count threshold met", "Linked"]).sum())
            weak = int(chain_summary["Strength"].astype(str).isin(["Below count threshold", "Weak link"]).sum())
            unanchored = int(chain_summary["Strength"].astype(str).eq("Unanchored").sum())
    common_status = "Review" if linked or weak else "Not applicable"
    rows.append({
        "Step": 2,
        "Check": "Common-anchor coverage",
        "Status": common_status,
        "Evidence": f"Facets meeting the count threshold: {linked}; below threshold: {weak}; without anchors: {unanchored}.",
        "Action": "Counts do not validate linking. Review reference stability and content equivalence. Joint calibration can connect groups without fixed anchors.",
    })

    if isinstance(drift, dict) and drift.get("available"):
        drift_summary = drift.get("summary", pd.DataFrame())
        drift_status = "Review"
        evidence = drift.get("reason", "Anchor drift review available.")
        if isinstance(drift_summary, pd.DataFrame) and not drift_summary.empty:
            drift_status = str(drift_summary.get("Status", pd.Series(["Review"])).iloc[0])
            evidence = (
                f"Max fixed-anchor drift: {drift_summary.get('MaxAbsAnchorDrift', pd.Series([np.nan])).iloc[0]}; "
                f"max group mean drift: {drift_summary.get('MaxAbsGroupMeanDrift', pd.Series([np.nan])).iloc[0]}."
            )
        rows.append({
            "Step": 3,
            "Check": "Hard-anchor consistency",
            "Status": "Ready" if drift_status == "Ready" else "Review",
            "Evidence": evidence,
            "Action": "Resolve constraint or convergence warnings. Agreement with imposed values is not an empirical test of reference stability.",
        })
    else:
        rows.append({
            "Step": 3,
            "Check": "Hard-anchor consistency",
            "Status": "Not applicable",
            "Evidence": drift.get("reason", "No fixed anchors or group anchors were supplied.") if isinstance(drift, dict) else "No drift review.",
            "Action": "No hard-anchor check is needed for joint calibration without fixed anchors; inspect shared observations and identification instead.",
        })

    chain_available = bool(isinstance(chain, dict) and chain.get("available"))
    rows.append({
        "Step": 4,
        "Check": "Equating-chain scope",
        "Status": "Review" if chain_available else "Not applicable",
        "Evidence": (
            chain.get("interpretation", chain.get("reason", "No chain summary."))
            if isinstance(chain, dict) else "No chain summary."
        ),
        "Action": "This is a constraint inventory, not validation of a linked scale. Retain independent reference-stability and historical evidence.",
    })
    rows.append({
        "Step": 5,
        "Check": "Next-run anchor package",
        "Status": "Optional",
        "Evidence": "Downloads include current estimated anchors, anchor audit tables, drift tables, and equating summary when available.",
        "Action": "Export anchors after a vetted run, version the file, and reuse unchanged labels for the next administration.",
    })
    return pd.DataFrame(rows)


def _population_sd_conditional_curvature_scale(loglik, sigma, bounds) -> float:
    """Diagnostic 1/sqrt(-d²ell/dsigma²) with all nuisance parameters fixed.

    This is not a profile standard error. A finite value does not establish
    joint stationarity, nuisance-adjusted uncertainty, or interval coverage.
    """
    sigma = float(sigma)
    lower, upper = map(float, bounds)
    h = max(1e-3, 0.01 * sigma)
    if not (
        np.isfinite([sigma, lower, upper]).all()
        and 0 < lower < sigma - h < sigma + h < upper
    ):
        return float("nan")
    try:
        values = [float(loglik(sigma + offset)) for offset in (h, 0.0, -h)]
        information = -(values[0] - 2.0 * values[1] + values[2]) / (h * h)
        if np.isfinite(values).all() and np.isfinite(information) and information > 0:
            return float(1.0 / np.sqrt(information))
    except (ArithmeticError, ValueError):
        pass
    return float("nan")


def mfrm_joint_mml(start, prep, idx, config, sizes, quad):
    """Explicit free-SD joint fit; numerical evidence never releases inference."""
    from types import SimpleNamespace
    from scipy.optimize import OptimizeResult
    from mfrm_app.mml_joint import fit_grouped_joint, OPTIONS, LIMITS

    if config.get("method") != "MML" or not config.get("estimate_population_sd"):
        raise ValueError("Joint MML requires MML with estimate_population_sd=True.")
    if sizes.get("theta", 0) or sizes.get("population", 0):
        raise ValueError("Joint MML excludes fitted Person coordinates and latent regression.")
    theta_spec = config.get("theta_spec", {})
    if (np.isfinite(theta_spec.get("anchors", [])).any()
        or any(g not in (None, "") for g in theta_spec.get("groups", []))):
        raise ValueError("Joint MML excludes Person anchors and group anchors.")
    if tuple(config.get("population_sd_bounds", ())) != (0.05, 10.0):
        raise ValueError("Joint MML requires the audited population-SD bounds (0.05, 10.0).")
    sigma_start = float(config["population_prior_sd"])
    if not 0.05 <= sigma_start <= 10.0:
        raise ValueError("Joint MML initial population SD must lie within (0.05, 10.0).")
    if config["model"] == "GPCM" and tuple(config.get("gpcm_log_slope_bounds", ())) != (-3.0, 3.0):
        raise ValueError("Joint GPCM requires the audited expanded log-slope bounds (-3.0, 3.0).")

    def residual(parameters):
        expanded = expand_params(parameters, sizes, config)
        values = [abs(float(np.sum(expanded["facets"][f]))) for f in config["facet_names"]
                  if config["facet_specs"][f]["centered"]]
        audit = build_parameterization_audit(config, sizes, expanded)
        if not bool(audit["ExpandedBoundsSatisfied"].all()):
            return float("inf")
        values.extend(audit["ConstraintMaxAbsResidual"].to_numpy(dtype=float))
        return float(max(values, default=0.0))

    structural_bounds = build_optimizer_bounds(sizes, config)
    constraints = [LinearConstraint(np.column_stack([np.atleast_2d(c.A),
        np.zeros((np.atleast_2d(c.A).shape[0], 1))]), c.lb, c.ub)
        for c in build_optimizer_constraints(sizes, config)]
    problem = SimpleNamespace(config=config, structural_start=np.asarray(start, dtype=float),
        sigma_start=sigma_start, sigma_bounds=config["population_sd_bounds"],
        structural_bounds=structural_bounds if structural_bounds is not None else [(None,None)]*len(start),
        joint_constraints=constraints, constraint_residual=residual)
    raw = fit_grouped_joint(prep["data"], problem, q=len(quad["nodes"]))
    if raw["point"] is None:
        raise ValueError("Joint MML did not produce a finite candidate: "+raw["status"])
    point = np.asarray(raw["point"], dtype=float)
    evaluated = raw["evaluated"]
    opt = OptimizeResult(x=point[:-1].copy(), fun=evaluated["nll"], jac=evaluated["gradient"][:-1].copy(),
        success=bool(raw["numerical_pass"]), status=0 if raw["numerical_pass"] else 1,
        message=raw["status"], nit=sum(r["iterations"] for r in raw["runs"]),
        nfev=sum(r["evaluations"] for r in raw["runs"]))
    opt.mml_engine = "joint"
    opt.estimated_population_sd = float(np.exp(point[-1]))
    opt.gradient_norm = float(np.linalg.norm(evaluated["gradient"]))
    opt.joint_gradient = evaluated["gradient"].copy()
    opt.joint_mml = dict(raw, evaluated={k:v for k,v in evaluated.items() if k != "posterior"},
        coordinate_order="native structural coordinates, then log(population SD)",
        input_sha256=dataframe_fingerprint(prep["data"], length=64),
        optimizer_options=[dict(o) for o in OPTIONS], numerical_limits=dict(LIMITS),
        quadrature_accuracy_qualified=False, scientific_inference_ready=False)
    return opt


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
    slope_facet=None,
    person_data=None,
    person_id_col=None,
    population_formula=None,
    population_standardize_numeric=False,
    population_categorical_terms=None,
    population_numeric_terms=None,
    anchor_df=None,
    group_anchor_df=None,
    noncenter_facet="Person",
    dummy_facets=None,
    positive_facets=None,
    quad_points=15,
    population_prior_sd=1.0,
    estimate_population_sd=False,
    facet_regularization=None,
    maxit=400,
    reltol=1e-6,
    mml_engine="EM",
    anchor_policy="warn",
    min_common_anchors=2,
    min_obs_per_element=2,
    min_obs_per_category=1,
    compute_plausible_values=False,
    n_plausible_values=5,
    plausible_seed=20260411,
):
    response_data_audit = build_response_data_audit(
        data,
        person_col=person_col,
        facet_cols=list(facet_cols),
        score_col=score_col,
        weight_col=weight_col,
    )
    data_for_fit = normalize_response_input_data(
        data,
        person_col=person_col,
        facet_cols=list(facet_cols),
        score_col=score_col,
        weight_col=weight_col,
    )
    prep = prepare_mfrm_data(
        data_for_fit,
        person_col=person_col,
        facet_cols=facet_cols,
        score_col=score_col,
        rating_min=rating_min,
        rating_max=rating_max,
        weight_col=weight_col,
        keep_original=keep_original,
    )
    prep["response_data_audit"] = response_data_audit
    prep["row_audit"] = response_data_audit.get("rows", pd.DataFrame())
    prep["excluded_rows"] = response_data_audit.get("excluded_rows", pd.DataFrame())
    prep["audit_summary"] = response_data_audit.get("summary", pd.DataFrame())

    model = str(model or "RSM").upper()
    if model not in {"RSM", "PCM", "GPCM"}:
        raise ValueError("Model must be one of: RSM, PCM, GPCM.")
    population_prior_sd = float(population_prior_sd if population_prior_sd is not None else 1.0)
    if not np.isfinite(population_prior_sd) or population_prior_sd <= 0:
        raise ValueError("population_prior_sd must be a positive finite value.")

    if model in {"PCM", "GPCM"}:
        if step_facet is None:
            step_facet = prep["facet_names"][0]
        if step_facet not in prep["facet_names"]:
            raise ValueError("Selected step facet is not in the facet list.")
        if model == "GPCM":
            if slope_facet is None:
                slope_facet = step_facet
            if slope_facet not in prep["facet_names"]:
                raise ValueError("Selected slope facet is not in the facet list.")
            if slope_facet != step_facet:
                raise ValueError("Bounded GPCM currently requires slope_facet == step_facet.")
        else:
            slope_facet = None
    else:
        step_facet = None
        slope_facet = None

    idx = build_indices(prep, step_facet=step_facet, slope_facet=slope_facet)

    n_person = len(prep["levels"]["Person"])
    facet_levels = {f: prep["levels"][f] for f in prep["facet_names"]}
    n_cat = prep["rating_max"] - prep["rating_min"] + 1
    if str(population_formula or "").strip() and method != "MML":
        raise ValueError("population_formula / latent regression is currently available for MML only.")
    population_model = (
        build_population_model(
            prep,
            person_data=person_data,
            person_id_col=person_id_col,
            population_formula=population_formula,
            standardize_numeric=population_standardize_numeric,
            categorical_terms=population_categorical_terms,
            numeric_terms=population_numeric_terms,
        )
        if method == "MML" else build_population_model(prep)
    )

    if noncenter_facet not in ["Person"] + prep["facet_names"]:
        noncenter_facet = "Person"

    dummy_facets = [f for f in (dummy_facets or []) if f in ["Person"] + prep["facet_names"]]

    positive_facets = set(positive_facets or [])
    facet_signs = {facet: (1 if facet in positive_facets else -1) for facet in prep["facet_names"]}
    anchor_policy = str(anchor_policy or "warn").lower()
    if anchor_policy not in {"warn", "error", "silent"}:
        anchor_policy = "warn"
    anchor_audit = audit_mfrm_anchors(
        prep,
        anchor_df=anchor_df,
        group_anchor_df=group_anchor_df,
        dummy_facets=dummy_facets,
        min_common_anchors=min_common_anchors,
        min_obs_per_element=min_obs_per_element,
        min_obs_per_category=min_obs_per_category,
    )
    if anchor_policy == "error" and anchor_audit.get("overall_status") in {"error", "warning"}:
        issues = anchor_audit.get("issues", pd.DataFrame())
        msg = "Anchor audit failed under policy='error'."
        if isinstance(issues, pd.DataFrame) and not issues.empty:
            first = issues.iloc[0]
            msg += f" First issue: {first.get('Message', '')} Action: {first.get('Action', '')}"
        raise ValueError(msg)
    anchor_df_for_fit = anchor_audit.get("valid_anchors", pd.DataFrame())
    group_anchor_df_for_fit = anchor_audit.get("valid_group_anchors", pd.DataFrame())

    config = {
        "model": model,
        "method": method,
        "app_version": APP_VERSION,
        "release_label": APP_RELEASE_LABEL,
        "runtime_scope": "standalone Python; no external MFRM engine called",
        "fit_df_method": DEFAULT_FIT_DF_METHOD,
        "fit_zstd_transform": fit_zstd_transform_label(False),
        "facets_zstd_cap": DEFAULT_FACETS_ZSTD_CAP,
        "n_person": n_person,
        "n_cat": n_cat,
        "facet_names": prep["facet_names"],
        "facet_levels": facet_levels,
        "step_facet": step_facet,
        "slope_facet": slope_facet,
        "population_model": population_model,
    }
    config["weight_col"] = weight_col if weight_col else None
    config["keep_original"] = bool(keep_original)
    config["rating_min"] = prep["rating_min"]
    config["rating_max"] = prep["rating_max"]
    config["input_rating_min"] = rating_min
    config["input_rating_max"] = rating_max
    config["score_range_explicit"] = bool(prep.get("score_range_explicit", False))
    config["positive_facets"] = list(positive_facets)
    config["facet_signs"] = facet_signs
    config["maxit"] = maxit
    config["reltol"] = reltol
    config["quad_points"] = quad_points if method == "MML" else None
    config["population_prior_sd"] = float(population_prior_sd) if method == "MML" else None
    config["population_formula"] = population_model.get("formula", "") if method == "MML" else ""
    config["population_standardize_numeric"] = bool(population_model.get("standardize_numeric", False)) if method == "MML" else False
    config["population_categorical_terms"] = population_model.get("categorical_terms", []) if method == "MML" else []
    config["population_numeric_terms"] = population_model.get("numeric_terms", []) if method == "MML" else []
    config["population_term_types"] = population_model.get("term_types", {}) if method == "MML" else {}
    config["person_id_col"] = person_id_col if method == "MML" and population_model.get("enabled") else None
    config["anchor_policy"] = anchor_policy
    config["compute_plausible_values"] = bool(compute_plausible_values) if method == "MML" else False
    config["n_plausible_values"] = int(n_plausible_values) if method == "MML" and compute_plausible_values else 0
    config["plausible_seed"] = int(plausible_seed) if method == "MML" and compute_plausible_values else None
    if model == "GPCM":
        config["gpcm_spec"] = {
            "scope": "bounded",
            "step_facet": step_facet,
            "slope_facet": slope_facet,
            "slope_identification": (
                "positive slopes; L-1 free log slopes with the final coordinate "
                "derived by exact sum-to-zero; geometric mean slope fixed to 1"
            ),
            "log_slope_bounds": [-3.0, 3.0],
        }
        config["gpcm_log_slope_bounds"] = (-3.0, 3.0)
    mml_engine_raw = str(mml_engine or "EM").strip().lower()
    config["mml_engine_requested"] = str(mml_engine) if method == "MML" else None
    if mml_engine_raw == "joint":
        mml_engine_resolved = "joint"
        if method != "MML" or not estimate_population_sd:
            raise ValueError("Joint MML requires MML with estimate_population_sd=True.")
        if ((anchor_df is not None and len(anchor_df))
            or (group_anchor_df is not None and len(group_anchor_df))):
            raise ValueError("Joint MML excludes individual and group anchors.")
    elif "auto" in mml_engine_raw:
        mml_engine_resolved = "auto"
    elif "hybrid" in mml_engine_raw:
        mml_engine_resolved = "hybrid"
    elif "direct" in mml_engine_raw:
        mml_engine_resolved = "direct"
    else:
        mml_engine_resolved = "em"
    # Existing free-SD paths retain their EM override. The explicit joint path
    # uses the audited two-stage solver and preserves the same inference hold.
    config["estimate_population_sd"] = bool(estimate_population_sd) if method == "MML" else False
    config["population_sd_bounds"] = (0.05, 10.0)
    if config["estimate_population_sd"] and mml_engine_resolved != "joint":
        maxit = max(int(maxit), 300)
        if mml_engine_resolved != "em":
            config["mml_engine_requested_before_sd_override"] = mml_engine_resolved
            config["population_sd_engine_notice"] = (
                "Population-SD estimation requires the EM engine; overrode "
                f"requested '{mml_engine_resolved}' to 'em'."
            )
            mml_engine_resolved = "em"
    config["mml_engine"] = mml_engine_resolved if method == "MML" else None
    exact_optimizer = "SLSQP" if model == "GPCM" else "L-BFGS-B"
    if method == "JMLE":
        config["optimizer"] = (
            "SLSQP with exact identified coordinates and bounded log slopes"
            if model == "GPCM"
            else "L-BFGS-B with exact identified step coordinates"
        )
    elif mml_engine_resolved == "joint":
        from mfrm_app.mml_joint import OPTIONS, LIMITS
        config["optimizer"] = "Joint structural/log-population-SD MML; two fixed SLSQP stages"
        config["joint_optimizer_options"] = [dict(o) for o in OPTIONS]
        config["joint_numerical_limits"] = dict(LIMITS)
        config["maxit_requested"] = maxit
        config["reltol_requested"] = reltol
        config["maxit"] = OPTIONS[-1]["maxiter"]
        config["reltol"] = OPTIONS[-1]["ftol"]
    elif mml_engine_resolved == "direct":
        config["optimizer"] = f"Direct MML with analytical-gradient {exact_optimizer}"
    elif mml_engine_resolved == "hybrid":
        config["optimizer"] = f"Hybrid MML: EM warm start + analytical-gradient {exact_optimizer}"
    elif mml_engine_resolved == "auto":
        config["optimizer"] = "Auto MML: hybrid first, EM fallback if needed"
    else:
        config["optimizer"] = (
            f"EM with analytical-gradient {exact_optimizer} M-step and free population SD"
            if config.get("estimate_population_sd")
            else f"EM with analytical-gradient {exact_optimizer} M-step"
        )

    constraint_specs = prepare_constraint_specs(
        prep=prep,
        anchor_df=anchor_df_for_fit,
        group_anchor_df=group_anchor_df_for_fit,
        noncenter_facet=noncenter_facet,
        dummy_facets=dummy_facets,
    )
    config["theta_spec"] = constraint_specs["theta_spec"]
    config["facet_specs"] = constraint_specs["facet_specs"]
    config["noncenter_facet"] = noncenter_facet
    config["dummy_facets"] = dummy_facets
    config["anchor_summary"] = constraint_specs["anchor_summary"]
    config["anchor_audit"] = anchor_audit

    # MML integrates person locations; theta_spec does not constrain EAP scores.
    if method == "MML" and (
        np.isfinite(config["theta_spec"]["anchors"]).any()
        or any(group not in (None, "") for group in config["theta_spec"]["groups"])
    ):
        raise ValueError(
            "MML does not support Person anchors, Person group anchors, or a dummy "
            "Person facet. Person scores are posterior estimates, not fixed "
            "parameters. Remove these anchor constraints, or use JMLE when fixed Person "
            "measures match the intended model. Non-Person facet anchors remain "
            "available in ordinary MML engines."
        )

    sizes = build_param_sizes(config)
    n_steps = max(n_cat - 1, 0)
    n_step_levels = 1 if model == "RSM" else len(facet_levels[step_facet])
    n_slope_levels = len(facet_levels[slope_facet]) if model == "GPCM" else 0
    config["parameterization"] = {
        "schema_version": "mfrm_exact_identified_v1",
        "step_coordinate": "last_derived_negative_sum",
        "step_expanded_count": int(n_step_levels * n_steps),
        "step_free_count": int(sizes.get("steps", 0)),
        "log_slope_coordinate": (
            "last_derived_negative_sum_with_expanded_bounds"
            if model == "GPCM"
            else "not_applicable"
        ),
        "log_slope_expanded_count": int(n_slope_levels),
        "log_slope_free_count": int(sizes.get("log_slopes", 0)),
        "optimizer_dimension": int(sum(sizes.values())) + int(method == "MML" and mml_engine_resolved == "joint"),
        "roadmap_issue_ids": ["STAT-001", "STAT-003"],
    }
    regularization_bundle = prepare_facet_regularization(facet_regularization, config, sizes)
    config["facet_regularization"] = regularization_bundle["runtime"]
    config["facet_regularization_scope"] = regularization_bundle["scope"]
    config["facet_regularization_mode"] = regularization_bundle["mode"]
    config["facet_regularization_objective_mode"] = regularization_bundle["mode"]
    config["facet_regularization_enabled"] = bool(regularization_bundle["enabled"])
    config["facet_regularization_fingerprint"] = regularization_bundle["runtime"].get("fingerprint")

    identifiability_bundle = build_jmle_eta_identifiability_audit(
        prep,
        config,
        sizes=sizes,
        idx=idx,
    )
    config["identifiability_scope"] = identifiability_bundle.get("scope")
    config["identifiability_status"] = identifiability_bundle.get("status")
    config["eta_structurally_identified"] = identifiability_bundle.get(
        "eta_structurally_identified"
    )
    identifiability_summary = identifiability_bundle.get("summary", pd.DataFrame())
    if isinstance(identifiability_summary, pd.DataFrame) and not identifiability_summary.empty:
        identifiability_row = identifiability_summary.iloc[0]
        config["eta_design_rank"] = identifiability_row.get("EtaDesignRank")
        config["eta_design_columns"] = identifiability_row.get("EtaDesignColumns")
        config["eta_structural_nullity"] = identifiability_row.get("EtaStructuralNullity")

    step_init_expanded = (
        np.linspace(-1, 1, n_steps)
        if n_steps > 0
        else np.array([], dtype=float)
    )
    step_init = sum_zero_free_from_expanded(step_init_expanded)
    facet_starts = np.concatenate([np.zeros(sizes[f]) for f in config["facet_names"]]) if config["facet_names"] else np.array([], dtype=float)
    start = np.concatenate([
        np.zeros(sizes["theta"]),
        np.zeros(sizes.get("population", 0)),
        facet_starts,
        step_init if model == "RSM" else np.tile(step_init, len(facet_levels[step_facet])),
        np.zeros(sizes.get("log_slopes", 0)),
    ])

    fit_start_time = time.perf_counter()
    if method == "JMLE":
        optimizer_method, optimizer_options = optimizer_method_and_options(
            sizes,
            config,
            maxiter=maxit,
            gtol=reltol,
            ftol=reltol,
        )
        opt = minimize(
            mfrm_loglik_jmle_value_grad,
            start,
            args=(idx, config, sizes),
            jac=True,
            method=optimizer_method,
            bounds=build_optimizer_bounds(sizes, config),
            constraints=build_optimizer_constraints(sizes, config),
            options=optimizer_options,
        )
    else:
        quad = make_mml_quadrature(config, quad_points)
        if mml_engine_resolved == "joint":
            opt = mfrm_joint_mml(start, prep, idx, config, sizes, quad)
        elif mml_engine_resolved == "direct":
            opt = mfrm_direct_mml(start, idx, config, sizes, quad,
                                  maxit=maxit, reltol=reltol)
        elif mml_engine_resolved == "hybrid":
            opt = mfrm_hybrid_mml(start, idx, config, sizes, quad,
                                  maxit=maxit, reltol=reltol)
        elif mml_engine_resolved == "auto":
            opt = mfrm_auto_mml(start, idx, config, sizes, quad,
                                maxit=maxit, reltol=reltol)
        else:
            opt = mfrm_em_mml(start, idx, config, sizes, quad,
                              maxit=maxit, reltol=reltol)
            opt.mml_engine = "em"
        config["mml_engine"] = getattr(opt, "mml_engine", mml_engine_resolved)
        if mml_engine_resolved == "auto":
            config["mml_engine_auto_selected"] = getattr(opt, "auto_selected", config["mml_engine"])
            config["mml_engine_auto_attempted"] = getattr(opt, "auto_attempted", ["hybrid"])
            config["mml_engine_auto_fallback_reason"] = getattr(opt, "auto_fallback_reason", "")
            if config["mml_engine"] == "em":
                config["optimizer"] = "Auto MML selected EM after hybrid fallback"
            elif config["mml_engine"] == "hybrid":
                config["optimizer"] = "Auto MML selected hybrid"
        # Thread the estimated population SD back into config so the final
        # marginal-NLL, EAP, and plausible-value quadratures (all built via
        # make_mml_quadrature, which reads population_prior_sd) use the fitted
        # scale. The user's original input is preserved separately.
        if config.get("estimate_population_sd"):
            est_sd = getattr(opt, "estimated_population_sd", None)
            if est_sd is not None and np.isfinite(est_sd) and est_sd > 0:
                config["estimated_population_sd"] = float(est_sd)
                config["population_prior_sd_input"] = float(population_prior_sd)
                config["population_prior_sd"] = float(est_sd)
                config["sigma_trace"] = list(getattr(opt, "sigma_trace", []) or [])
    elapsed_seconds = time.perf_counter() - fit_start_time

    params = expand_params(opt.x, sizes, config)

    # Keep fixed-nuisance curvature as a technical diagnostic, never an SE/CI.
    if config.get("estimate_population_sd") and config.get("estimated_population_sd"):
        sigma_hat = float(config["estimated_population_sd"])
        n_quad = int(quad_points or config.get("quad_points") or 15)

        def _mll_at_sigma(s):
            _, ll = _e_step_posteriors(
                idx, config, params, gauss_hermite_normal(n_quad, sd=s)
            )
            return float(ll)

        config["population_sd_conditional_curvature_scale"] = (
            _population_sd_conditional_curvature_scale(
                _mll_at_sigma, sigma_hat,
                config.get("population_sd_bounds", (0.05, 10.0)),
            )
        )
        config["population_sd_se"] = None
        config["population_sd_ci"] = None
        config["population_sd_inference_ready"] = False
        config["population_sd_se_basis"] = (
            "SE/CI withheld: fixed-nuisance curvature is not profile information. "
            "Joint stationarity, nuisance-adjusted uncertainty, and coverage "
            "have not been qualified."
        )

    if method == "MML":
        quad = make_mml_quadrature(config, quad_points)
        posterior_outputs = compute_person_posterior_outputs(
            idx,
            config,
            params,
            quad,
            person_levels=prep["levels"]["Person"],
            n_plausible_values=(int(n_plausible_values) if compute_plausible_values else 0),
            seed=int(plausible_seed),
        )
        person_columns = ["Person", "PosteriorMean", "PosteriorSD"]
        if population_model.get("enabled"):
            person_columns.append("PopulationMean")
        person_tbl = posterior_outputs["scores"][person_columns].rename(
            columns={"PosteriorMean": "Estimate", "PosteriorSD": "SD"}
        ).copy()
    else:
        person_tbl = pd.DataFrame({
            "Person": prep["levels"]["Person"],
            "Estimate": params["theta"],
        })
        posterior_outputs = {
            "available": False,
            "reason": "Posterior scoring and plausible values are defined for MML runs only.",
            "scores": pd.DataFrame(),
            "plausible_values": pd.DataFrame(),
        }

    facet_boundary = _jmle_boundary.facet_score_boundary_audit(prep, config)
    facet_boundary_hold = bool(facet_boundary["RequiresReview"].any())
    config["jmle_facet_boundary_review"] = facet_boundary_hold

    fit_inference_ready_for_person = bool(getattr(opt, "success", False)) and not facet_boundary_hold and (
        bool(identifiability_bundle.get("inference_ready", False))
        if method == "JMLE" else True
    )
    person_boundary_bundle = build_jmle_person_boundary_audit(
        prep,
        config,
        estimates=params["theta"] if method == "JMLE" else None,
        fit_inference_ready=fit_inference_ready_for_person,
    )
    config["person_boundary_scope"] = person_boundary_bundle.get("scope")
    config["person_boundary_status"] = person_boundary_bundle.get("status")
    if method == "JMLE":
        person_boundary_rows = person_boundary_bundle.get("persons", pd.DataFrame())
        if isinstance(person_boundary_rows, pd.DataFrame) and not person_boundary_rows.empty:
            annotation_columns = [
                "Person",
                "ExtremeScorePattern",
                "ExtremeScoreDirection",
                "FiniteJMLEEstimate",
                "PersonInferenceReady",
                "ReportableEstimate",
                "EstimateRole",
            ]
            person_tbl = person_tbl.merge(
                person_boundary_rows[annotation_columns],
                on="Person",
                how="left",
                validate="one_to_one",
            )

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

    if model == "GPCM":
        slope_tbl = pd.DataFrame({
            "SlopeFacet": slope_facet,
            "Level": prep["levels"][slope_facet],
            "LogEstimate": params["log_slopes"],
            "Estimate": params["slopes"],
        })
    else:
        slope_tbl = pd.DataFrame(columns=["SlopeFacet", "Level", "LogEstimate", "Estimate"])

    if method == "MML" and population_model.get("enabled") and sizes.get("population", 0) > 0:
        population_tbl = pd.DataFrame({
            "Term": population_model.get("columns", []),
            "Estimate": params["population"],
        })
        person_population_tbl = population_model.get("person_data", pd.DataFrame()).copy()
        pop_mu = compute_population_mu(params, config)
        if isinstance(person_population_tbl, pd.DataFrame) and not person_population_tbl.empty:
            person_population_tbl["PopulationMean"] = pop_mu
    else:
        population_tbl = pd.DataFrame(columns=["Term", "Estimate"])
        person_population_tbl = pd.DataFrame()
    transform_rows = []
    if method == "MML" and population_model.get("enabled"):
        for term, transform in population_model.get("numeric_transforms", {}).items():
            if isinstance(transform, dict):
                transform_rows.append({
                    "Term": term,
                    "Center": float(transform.get("center", 0.0)),
                    "Scale": float(transform.get("scale", 1.0)),
                    "Standardized": bool(transform.get("standardized", False)),
                })
    population_transform_tbl = pd.DataFrame(
        transform_rows,
        columns=["Term", "Center", "Scale", "Standardized"],
    )

    k_params = int(sum(sizes.values())) + (1 if config.get("estimate_population_sd") else 0)
    config["parameter_count"] = k_params
    if method == "MML":
        final_quad_for_ll = make_mml_quadrature(config, quad_points)
        final_raw_nll = mfrm_loglik_mml(opt.x, idx, config, sizes, final_quad_for_ll)
    else:
        final_raw_nll = mfrm_loglik_jmle(opt.x, idx, config, sizes)
    penalty_value, _ = facet_regularization_value_grad(opt.x, sizes, config)
    penalized_objective = final_raw_nll + penalty_value
    loglik = -final_raw_nll
    opt.raw_fun = final_raw_nll
    opt.penalty_value = penalty_value
    opt.penalized_objective = penalized_objective
    if not np.isfinite(getattr(opt, "gradient_norm", np.nan)):
        try:
            if method == "MML":
                final_quad = make_mml_quadrature(config, quad_points)
                _, final_grad = mfrm_loglik_mml_value_grad(opt.x, idx, config, sizes, final_quad)
            else:
                _, final_grad = mfrm_loglik_jmle_value_grad(opt.x, idx, config, sizes)
            opt.gradient_norm = float(np.linalg.norm(final_grad))
        except Exception:
            opt.gradient_norm = np.nan
    opt.elapsed_seconds = float(elapsed_seconds)
    convergence_tbl = build_convergence_summary(opt, config, loglik, elapsed_seconds)
    if "Weight" in prep["data"].columns:
        n_obs = float(prep["data"]["Weight"].sum())
    else:
        n_obs = len(prep["data"])
    aic = 2 * k_params - 2 * loglik
    bic_sample_size = n_person if method == "MML" else n_obs
    bic = np.log(bic_sample_size) * k_params - 2 * loglik if bic_sample_size > 0 else np.nan
    deviance = -2 * loglik
    loglik_per_obs = loglik / n_obs if n_obs > 0 else np.nan
    aic_per_obs = aic / n_obs if n_obs > 0 else np.nan
    bic_per_obs = bic / n_obs if n_obs > 0 else np.nan

    n_iter = getattr(opt, "nit", opt.nfev)
    eta_ready = (
        bool(identifiability_bundle.get("inference_ready", False))
        if method == "JMLE" else True
    )
    inference_ready = bool(opt.success) and eta_ready and not facet_boundary_hold
    summary_tbl = pd.DataFrame({
        "Model": [model],
        "Method": [method],
        "EstimatorLabel": [f"Penalized {method}" if regularization_bundle["enabled"] else method],
        "MML Engine": [config.get("mml_engine") if method == "MML" else None],
        "PopulationPriorSD": [config.get("population_prior_sd") if method == "MML" else None],
        "FacetRegularization": ["On" if regularization_bundle["enabled"] else "Off"],
        "Optimizer": [config.get("optimizer")],
        "N": [n_obs],
        "Persons": [n_person],
        "Facets": [len(config["facet_names"])],
        "Categories": [n_cat],
        "KParams": [k_params],
        "PopulationSDMode": [
            "estimated" if config.get("estimate_population_sd")
            else ("fixed" if config.get("method") == "MML" else "n/a")
        ],
        "EstimatedPopulationSD": [config.get("estimated_population_sd")],
        "LogLik": [loglik],
        "Deviance": [deviance],
        "RegularizationPenalty": [penalty_value],
        "PenalizedObjective": [penalized_objective],
        "ObjectiveBasis": [
            "negative log-likelihood + Gaussian facet penalty"
            if regularization_bundle["enabled"] else "negative log-likelihood"
        ],
        "LogLikBasis": ["unpenalized likelihood evaluated at fitted estimates"],
        "AIC": [aic],
        "BIC": [bic],
        "BICSampleSize": [bic_sample_size],
        "LogLikPerObs": [loglik_per_obs],
        "AICPerObs": [aic_per_obs],
        "BICPerObs": [bic_per_obs],
        "Converged": [bool(opt.success)],
        "InferenceReady": [inference_ready],
        "EtaStructurallyIdentified": [
            identifiability_bundle.get("eta_structurally_identified")
        ],
        "EtaDesignRank": [config.get("eta_design_rank")],
        "EtaDesignColumns": [config.get("eta_design_columns")],
        "EtaStructuralNullity": [config.get("eta_structural_nullity")],
        "IdentifiabilityStatus": [identifiability_bundle.get("status")],
        "ExtremeJMLEPersons": [
            int(person_boundary_bundle.get("summary", pd.DataFrame()).iloc[0]["ExtremePersons"])
            if method == "JMLE" else pd.NA
        ],
        "FiniteJMLEPersonEstimates": [
            int(person_boundary_bundle.get("summary", pd.DataFrame()).iloc[0]["FiniteJMLEPersons"])
            if method == "JMLE" else pd.NA
        ],
        "PersonMeasureInferenceReady": [
            person_boundary_bundle.get("person_measure_inference_ready")
        ],
        "PersonBoundaryStatus": [person_boundary_bundle.get("status")],
        "Iterations": [n_iter],
        "GradientNorm": [getattr(opt, "gradient_norm", np.nan)],
        "ElapsedSeconds": [elapsed_seconds],
    })
    parameterization_audit = build_parameterization_audit(config, sizes, params)

    result = {
        "summary": summary_tbl,
        "facets": {"person": person_tbl, "others": facet_tbl},
        "steps": step_tbl,
        "slopes": slope_tbl,
        "population": {
            "enabled": bool(population_model.get("enabled", False)),
            "formula": population_model.get("formula", ""),
            "coefficients": population_tbl,
            "person_data": person_population_tbl,
            "numeric_transforms": population_transform_tbl,
            "term_types": dict(population_model.get("term_types", {})),
            "type_summary": population_model.get(
                "type_summary",
                pd.DataFrame(columns=POPULATION_TYPE_SUMMARY_COLUMNS),
            ),
            "categorical_terms": list(population_model.get("categorical_terms", [])),
            "numeric_terms": list(population_model.get("numeric_terms", [])),
            "warnings": list(population_model.get("warnings", [])),
        },
        "config": config,
        "prep": prep,
        "opt": opt,
        "params": params,
        "parameterization_audit": parameterization_audit,
        "identifiability": identifiability_bundle,
        "person_boundary": person_boundary_bundle,
        "facet_boundary": facet_boundary,
        "convergence": convergence_tbl,
        "posterior": posterior_outputs,
        "regularization": {
            "settings": regularization_bundle["spec"],
            "audit": regularization_bundle["audit"],
            "penalty_summary": facet_regularization_penalty_summary(opt.x, sizes, config),
        },
    }
    # Qualification needs the optimizer record and identification audit, not just the table.
    result["summary"] = qualified_estimation_summary(result)
    if method == "MML" and config.get("mml_engine") == "joint":
        result["joint_mml"] = opt.joint_mml
    result["likelihood_information"] = build_likelihood_information_criteria(result)
    result["anchor_drift"] = build_anchor_drift_review(result)
    result["equating_chain"] = build_equating_chain_summary(result)
    result["anchor_equating_workflow"] = build_anchor_equating_workflow_plan(result)
    config["anchor_drift_scope"] = result["anchor_drift"].get("scope")
    config["equating_chain_scope"] = result["equating_chain"].get("scope")
    return result


def compute_obs_table(res):
    prep = res["prep"]
    config = res["config"]
    idx = build_indices(prep, step_facet=config["step_facet"], slope_facet=config.get("slope_facet"))
    params = res["params"]
    theta_hat = params["theta"] if config["method"] == "JMLE" else res["facets"]["person"]["Estimate"].to_numpy()

    person_measure = theta_hat
    person_measure_by_row = person_measure[idx["person"]] if person_measure.size else np.zeros_like(idx["score_k"], dtype=float)

    eta = compute_eta(idx, params, config, theta_override=theta_hat if theta_hat.size else None)
    if config["model"] == "RSM":
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
        probs = category_prob_rsm(eta, step_cum)
    elif config["model"] == "PCM":
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])
        probs = category_prob_pcm(eta, step_cum_mat, idx["step_idx"])
    else:
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])
        probs = category_prob_gpcm(
            eta,
            step_cum_mat,
            idx["step_idx"],
            params["slopes"],
            idx.get("slope_idx"),
        )

    k_vals = np.arange(probs.shape[1])
    expected_k = probs.dot(k_vals)
    var_k = probs.dot(k_vals ** 2) - expected_k ** 2
    var_k = np.where(var_k <= 1e-10, np.nan, var_k)
    resid_k = idx["score_k"] - expected_k
    std_sq = resid_k ** 2 / var_k
    # E[(X - E[X])^4 | theta_i] for the polytomous model, in closed form:
    # M4_i = sum_k P_k * (k - E[X])^4. Used by the FACETS / Wright-Masters
    # Welch-Satterthwaite d.f. for ZSTD reporting (Wright & Masters, 1982,
    # Eqs. 4.20 and 4.27).
    diff_k = k_vals[np.newaxis, :] - expected_k[:, np.newaxis]
    fourth_central_moment = np.sum(probs * (diff_k ** 4), axis=1)

    df = prep["data"].copy()
    df["PersonMeasure"] = person_measure_by_row
    df["Observed"] = prep["rating_min"] + idx["score_k"]
    df["Expected"] = prep["rating_min"] + expected_k
    df["Var"] = var_k
    df["FourthCentralMoment"] = fourth_central_moment
    df["Residual"] = df["Observed"] - df["Expected"]
    df["StdResidual"] = df["Residual"] / np.sqrt(df["Var"])
    df["StdSq"] = std_sq

    # ---- Per-observation diagnostics for polytomous person-fit ----
    # These six columns are the inputs to the Drasgow (1985) lz statistic
    # and the Snijders (2001) lz* correction (Snijders, 2001, Eq. 16). The
    # math is summarised in the lz / lz* helpers below; the columns are
    # added here so a single pass over the per-observation
    # probability matrix produces every term they need.
    n_rows = probs.shape[0]
    score_k_int = idx["score_k"].astype(int)
    pr_observed = np.full(n_rows, np.nan)
    if n_rows > 0:
        valid_score = (score_k_int >= 0) & (score_k_int < probs.shape[1])
        rows_idx = np.where(valid_score)[0]
        pr_observed[rows_idx] = probs[rows_idx, score_k_int[rows_idx]]
    # Use clipping to avoid log(0); the clip floor mirrors the R reference
    # (.Machine$double.eps in R, np.finfo(float).eps here).
    eps = float(np.finfo(float).eps)
    log_p_mat = np.log(np.clip(probs, eps, None))
    # E[log P]_i = sum_k P_k * log P_k   (per-item negentropy; the R code
    # calls this ItemEntropy for parity with mfrmr 0.2.0, even though
    # information-theoretically it is the negative of Shannon entropy).
    item_entropy = np.sum(probs * log_p_mat, axis=1)
    # Var[log P]_i = sum_k P_k * (log P_k)^2  -  (E[log P])^2
    item_var_logp = np.sum(probs * log_p_mat ** 2, axis=1) - item_entropy ** 2
    item_var_logp = np.maximum(item_var_logp, 0.0)
    # Slope per observation (1.0 for RSM / PCM; estimated for GPCM under
    # the geometric-mean-one identification).
    if config["model"] == "GPCM":
        slope_idx_arr = idx.get("slope_idx")
        if slope_idx_arr is None:
            slope_idx_arr = idx["step_idx"]
        slope_obs = np.asarray(params["slopes"], dtype=float)[slope_idx_arr]
    else:
        slope_obs = np.ones(n_rows, dtype=float)
    # ScoreInformation_i = a^2 * Var[X|theta]   (Muraki, 1993, Eq. 16).
    score_information = (slope_obs ** 2) * var_k
    # ObservedScoreDerivative_i = d log P_i / d theta evaluated at the
    # observed category: a_i * (X_i - E[X|theta]). At the JML estimate
    # of theta_n the sum across a person's observations is zero
    # (definition of the MLE); for any individual observation it can
    # take either sign.
    observed_score_derivative = slope_obs * (idx["score_k"].astype(float) - expected_k)
    # ItemLogPScoreCov_i = Cov[log P, score]_i
    #   = a_i * Cov[log P, X | theta]_i
    #   = a_i * (E[X * log P]_i - E[X]_i * E[log P]_i)
    # where E[X * log P]_i = sum_k k * P_k * log P_k.
    item_e_k_logp = np.sum(k_vals * probs * log_p_mat, axis=1)
    item_logp_score_cov = slope_obs * (item_e_k_logp - expected_k * item_entropy)

    df["PrObserved"] = pr_observed
    df["ItemEntropy"] = item_entropy
    df["ItemVarLogP"] = item_var_logp
    df["ItemLogPScoreCov"] = item_logp_score_cov
    df["ScoreInformation"] = score_information
    df["ObservedScoreDerivative"] = observed_score_derivative
    return df


def compute_prob_matrix(res):
    prep = res["prep"]
    config = res["config"]
    idx = build_indices(prep, step_facet=config["step_facet"], slope_facet=config.get("slope_facet"))
    params = res["params"]
    theta_hat = params["theta"] if config["method"] == "JMLE" else res["facets"]["person"]["Estimate"].to_numpy()
    eta = compute_eta(idx, params, config, theta_override=theta_hat if theta_hat.size else None)
    if config["model"] == "RSM":
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
        return category_prob_rsm(eta, step_cum)
    step_cum_mat = np.vstack([
        np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
    ])
    if config["model"] == "PCM":
        return category_prob_pcm(eta, step_cum_mat, idx["step_idx"])
    return category_prob_gpcm(
        eta,
        step_cum_mat,
        idx["step_idx"],
        params["slopes"],
        idx.get("slope_idx"),
    )


def compute_scorefile(res):
    obs = compute_obs_table(res)
    if obs.empty:
        return pd.DataFrame()
    probs = compute_prob_matrix(res)
    if probs is None or probs.shape[0] != len(obs):
        return obs
    cat_vals = np.arange(res["prep"]["rating_min"], res["prep"]["rating_max"] + 1)
    # Probabilities follow observation order; retained row labels can have
    # gaps or duplicates after filtering. Preserve them during concatenation.
    prob_df = pd.DataFrame(probs, index=obs.index, columns=[f"P_{v}" for v in cat_vals])
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


def get_weights(df):
    if "Weight" in df.columns:
        w = df["Weight"].to_numpy(dtype=float)
        w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
        return w
    return np.ones(len(df), dtype=float)


def calc_overall_fit(
    obs_df,
    whexact: bool = False,
    fit_df_method: str = DEFAULT_FIT_DF_METHOD,
    facets_zstd_cap: float = DEFAULT_FACETS_ZSTD_CAP,
    *,
    zstd_transform=None,
):
    """Overall Infit / Outfit mean-square with optional FACETS d.f.

    The primary ``InfitZSTD`` / ``OutfitZSTD`` columns use the FACETS /
    Wright-Masters Welch-Satterthwaite d.f. by default and retain the
    historical engine values under ``*_ENGINE`` sidecar columns. Pass
    ``fit_df_method="engine"`` for the explicit legacy convention, or
    ``"both"`` for an engine-primary comparison view. The
    ``facets_zstd_cap`` argument (default 9.0) bounds the FACETS ZSTD
    on extreme cells to match the FACETS software convention.
    """
    w = get_weights(obs_df)
    var = pd.to_numeric(obs_df["Var"], errors="coerce").to_numpy(dtype=float)
    fourth = (
        pd.to_numeric(obs_df.get("FourthCentralMoment"), errors="coerce").to_numpy(dtype=float)
        if "FourthCentralMoment" in obs_df.columns
        else np.full(len(obs_df), np.nan)
    )
    infit = np.nansum(obs_df["StdSq"] * obs_df["Var"] * w) / np.nansum(obs_df["Var"] * w)
    outfit = np.nansum(obs_df["StdSq"] * w) / np.nansum(w)
    df_infit = float(np.nansum(obs_df["Var"] * w))
    df_outfit = float(np.nansum(w))
    denom_infit, denom_outfit = _facets_fit_df_terms(var, fourth, w)
    df_infit_facets, df_outfit_facets = _facets_fit_df(df_infit, df_outfit, denom_infit, denom_outfit)
    tbl = pd.DataFrame({
        "Infit": [infit],
        "Outfit": [outfit],
        "InfitZSTD": [zstd_from_mnsq(infit, df_infit, whexact=whexact, zstd_transform=zstd_transform)],
        "OutfitZSTD": [zstd_from_mnsq(outfit, df_outfit, whexact=whexact, zstd_transform=zstd_transform)],
        "DF_Infit": [df_infit],
        "DF_Outfit": [df_outfit],
        "DF_Infit_FACETS": [df_infit_facets],
        "DF_Outfit_FACETS": [df_outfit_facets],
    })
    return _apply_fit_df_method(
        tbl,
        fit_df_method=fit_df_method,
        whexact=whexact,
        zstd_transform=zstd_transform,
        facets_zstd_cap=facets_zstd_cap,
    )


def calc_facet_fit(
    obs_df,
    facet_cols,
    whexact: bool = False,
    fit_df_method: str = DEFAULT_FIT_DF_METHOD,
    facets_zstd_cap: float = DEFAULT_FACETS_ZSTD_CAP,
    *,
    zstd_transform=None,
):
    """Per-facet/-level Infit / Outfit mean-square with optional FACETS d.f."""
    rows = []
    has_fourth = "FourthCentralMoment" in obs_df.columns
    for facet in facet_cols:
        grp = obs_df.groupby(facet, observed=False)
        for level, df in grp:
            w = get_weights(df)
            infit = np.nansum(df["StdSq"] * df["Var"] * w) / np.nansum(df["Var"] * w)
            outfit = np.nansum(df["StdSq"] * w) / np.nansum(w)
            df_infit = float(np.nansum(df["Var"] * w))
            df_outfit = float(np.nansum(w))
            var = pd.to_numeric(df["Var"], errors="coerce").to_numpy(dtype=float)
            fourth = (
                pd.to_numeric(df["FourthCentralMoment"], errors="coerce").to_numpy(dtype=float)
                if has_fourth
                else np.full(len(df), np.nan)
            )
            denom_infit, denom_outfit = _facets_fit_df_terms(var, fourth, w)
            df_infit_facets, df_outfit_facets = _facets_fit_df(
                df_infit, df_outfit, denom_infit, denom_outfit
            )
            rows.append({
                "Facet": facet,
                "Level": level,
                "N": np.nansum(w),
                "Infit": infit,
                "Outfit": outfit,
                "InfitZSTD": zstd_from_mnsq(infit, df_infit, whexact=whexact, zstd_transform=zstd_transform),
                "OutfitZSTD": zstd_from_mnsq(outfit, df_outfit, whexact=whexact, zstd_transform=zstd_transform),
                "DF_Infit": df_infit,
                "DF_Outfit": df_outfit,
                "DF_Infit_FACETS": df_infit_facets,
                "DF_Outfit_FACETS": df_outfit_facets,
            })
    return _apply_fit_df_method(
        pd.DataFrame(rows),
        fit_df_method=fit_df_method,
        whexact=whexact,
        zstd_transform=zstd_transform,
        facets_zstd_cap=facets_zstd_cap,
    )


def calc_facet_se(obs_df, facet_cols):
    rows = []
    for facet in facet_cols:
        grp = obs_df.groupby(facet, observed=False)
        for level, df in grp:
            w = get_weights(df)
            # GPCM location information is a^2 Var(K); Rasch has a = 1.
            info = np.nansum(df.get("ScoreInformation", df["Var"]) * w)
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
    if len(facet_cols) < 2 or obs_df is None:
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

    w_all = np.asarray(get_weights(obs_df), dtype=float)
    w_finite_mask = np.isfinite(w_all)
    w_pos_mask = w_finite_mask & (w_all > 0)
    w_nansafe = np.where(w_finite_mask, w_all, 0.0)

    obs_arr = obs_df["Observed"].to_numpy(dtype=float)
    exp_arr = obs_df["Expected"].to_numpy(dtype=float)
    res_arr = obs_df["Residual"].to_numpy(dtype=float)
    std_arr = obs_df["StdResidual"].to_numpy(dtype=float)
    abs_std_arr = np.abs(std_arr)
    var_arr = obs_df["Var"].to_numpy(dtype=float)

    def _cell_weighted_mean(values: np.ndarray, cell_codes: np.ndarray, valid_rows_mask: np.ndarray, n_cells: int) -> np.ndarray:
        ok = valid_rows_mask & w_pos_mask & np.isfinite(values)
        out = np.full(n_cells, np.nan, dtype=float)
        if not np.any(ok):
            return out
        c_ok = cell_codes[ok]
        w_ok = w_all[ok]
        den = np.bincount(c_ok, weights=w_ok, minlength=n_cells)
        num = np.bincount(c_ok, weights=values[ok] * w_ok, minlength=n_cells)
        pos = den > 0
        out[pos] = num[pos] / den[pos]
        return out

    def _cell_nansum_weighted(values: np.ndarray, cell_codes: np.ndarray, valid_rows_mask: np.ndarray, n_cells: int) -> np.ndarray:
        prod = values * w_all
        ok = valid_rows_mask & np.isfinite(prod)
        if not np.any(ok):
            return np.zeros(n_cells, dtype=float)
        return np.bincount(cell_codes[ok], weights=prod[ok], minlength=n_cells)

    def _encode_axis(series: pd.Series) -> tuple[np.ndarray, list[object]]:
        if isinstance(series.dtype, pd.CategoricalDtype):
            return series.cat.codes.to_numpy(dtype=np.int64), list(series.cat.categories)
        codes, uniques = pd.factorize(series, sort=True)
        return codes.astype(np.int64), list(uniques)

    encoded_axes: dict[str, tuple[np.ndarray, list[object]]] = {}
    part_frames: list[pd.DataFrame] = []
    for pair1, pair2 in combos:
        if pair1 not in encoded_axes:
            encoded_axes[pair1] = _encode_axis(obs_df[pair1])
        if pair2 not in encoded_axes:
            encoded_axes[pair2] = _encode_axis(obs_df[pair2])
        c1, lvl1_list = encoded_axes[pair1]
        c2, lvl2_list = encoded_axes[pair2]
        n2 = max(len(lvl2_list), 1)
        valid_mask = (c1 >= 0) & (c2 >= 0)
        if not np.any(valid_mask):
            continue
        raw_keys = c1[valid_mask] * np.int64(n2) + c2[valid_mask]
        unique_keys, codes_valid = np.unique(raw_keys, return_inverse=True)
        n_cells = len(unique_keys)
        if n_cells == 0:
            continue
        cell_codes = np.zeros(len(obs_df), dtype=int)
        cell_codes[valid_mask] = codes_valid
        u1_idx = (unique_keys // np.int64(n2)).astype(int)
        u2_idx = (unique_keys % np.int64(n2)).astype(int)
        level_labels = [f"{lvl1_list[i1]} | {lvl2_list[i2]}" for i1, i2 in zip(u1_idx, u2_idx)]

        n_arr = (
            np.bincount(cell_codes[valid_mask], weights=w_nansafe[valid_mask], minlength=n_cells)
            if np.any(valid_mask)
            else np.zeros(n_cells, dtype=float)
        )
        observed_avg = _cell_weighted_mean(obs_arr, cell_codes, valid_mask, n_cells)
        expected_avg = _cell_weighted_mean(exp_arr, cell_codes, valid_mask, n_cells)
        mean_resid = _cell_weighted_mean(res_arr, cell_codes, valid_mask, n_cells)
        mean_std_resid = _cell_weighted_mean(std_arr, cell_codes, valid_mask, n_cells)
        mean_abs_std = _cell_weighted_mean(abs_std_arr, cell_codes, valid_mask, n_cells)
        chi_sq = _cell_nansum_weighted(std_arr ** 2, cell_codes, valid_mask, n_cells)
        var_sum = _cell_nansum_weighted(var_arr, cell_codes, valid_mask, n_cells)

        pos_n = n_arr > 0
        se_resid = np.full(n_cells, np.nan, dtype=float)
        se_std = np.full(n_cells, np.nan, dtype=float)
        se_resid[pos_n] = np.sqrt(var_sum[pos_n]) / n_arr[pos_n]
        se_std[pos_n] = 1.0 / np.sqrt(n_arr[pos_n])
        df_chi = np.where(n_arr > 1.0, n_arr - 1.0, np.nan)

        t_resid = np.where(np.isfinite(se_resid) & (se_resid > 0), mean_resid / se_resid, np.nan)
        t_std = np.where(np.isfinite(se_std) & (se_std > 0), mean_std_resid / se_std, np.nan)

        p_resid = np.full(n_cells, np.nan, dtype=float)
        ok_pr = np.isfinite(df_chi) & np.isfinite(t_resid)
        if np.any(ok_pr):
            p_resid[ok_pr] = 2.0 * t_dist.cdf(-np.abs(t_resid[ok_pr]), df=df_chi[ok_pr])

        p_std = np.full(n_cells, np.nan, dtype=float)
        ok_ps = np.isfinite(df_chi) & np.isfinite(t_std)
        if np.any(ok_ps):
            p_std[ok_ps] = 2.0 * t_dist.cdf(-np.abs(t_std[ok_ps]), df=df_chi[ok_ps])

        chi_p = np.full(n_cells, np.nan, dtype=float)
        ok_cp = np.isfinite(chi_sq) & np.isfinite(df_chi) & (df_chi > 0)
        if np.any(ok_cp):
            chi_p[ok_cp] = 1.0 - chi2.cdf(chi_sq[ok_cp], df=df_chi[ok_cp])

        part_frames.append(pd.DataFrame({
            "Pair": f"{pair1} x {pair2}",
            "Level": level_labels,
            "N": n_arr,
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
        }))

    if not part_frames:
        return pd.DataFrame()
    out = pd.concat(part_frames, ignore_index=True)
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


def compute_empirical_bayes_shrinkage(
    measures: pd.DataFrame,
    *,
    min_levels: int = 3,
    shrink_person: bool = False,
    facet_prior_sd: float | None = None,
) -> dict:
    """Post-hoc normal-normal plug-in shrinkage for facet measures."""
    empty = {
        "available": False,
        "reason": "Measure table is unavailable.",
        "measures": pd.DataFrame(),
        "report": pd.DataFrame(),
    }
    required = {"Facet", "Level", "Estimate", "SE"}
    if not isinstance(measures, pd.DataFrame) or measures.empty or not required.issubset(measures.columns):
        return empty
    rows: list[pd.DataFrame] = []
    report_rows: list[dict[str, object]] = []
    for facet, df in measures.groupby("Facet", observed=False):
        if str(facet) == "Person" and not shrink_person:
            continue
        sub = df.copy()
        sub["Estimate"] = pd.to_numeric(sub["Estimate"], errors="coerce")
        sub["SE"] = pd.to_numeric(sub["SE"], errors="coerce")
        valid = sub["Estimate"].notna() & sub["SE"].notna() & np.isfinite(sub["SE"]) & (sub["SE"] > 0)
        k_valid = int(valid.sum())
        sub["ShrinkageTau2"] = np.nan
        sub["ShrinkageFactor"] = 0.0
        sub["ShrunkEstimate"] = sub["Estimate"]
        sub["ShrunkSE"] = sub["SE"]
        sub["ShrinkageStatus"] = "not enough valid levels" if k_valid < int(min_levels) else "not shrunk"
        if k_valid >= int(min_levels):
            est = sub.loc[valid, "Estimate"].to_numpy(dtype=float)
            se = sub.loc[valid, "SE"].to_numpy(dtype=float)
            if facet_prior_sd is not None and np.isfinite(float(facet_prior_sd)) and float(facet_prior_sd) > 0:
                tau2 = float(facet_prior_sd) ** 2
                tau_source = "user-specified prior SD"
            else:
                tau2 = max(0.0, float(np.mean(est ** 2) - np.mean(se ** 2)))
                tau_source = "method-of-moments EB"
            denom = tau2 + se ** 2
            b = np.where(denom > 0, se ** 2 / denom, 1.0)
            b = np.clip(b, 0.0, 1.0)
            shrunk = (1.0 - b) * est
            shrunk_se = np.sqrt(np.maximum((1.0 - b) * se ** 2, 0.0))
            sub.loc[valid, "ShrinkageTau2"] = tau2
            sub.loc[valid, "ShrinkageFactor"] = b
            sub.loc[valid, "ShrunkEstimate"] = shrunk
            sub.loc[valid, "ShrunkSE"] = shrunk_se
            sub.loc[valid, "ShrinkageStatus"] = "shrunk"
            report_rows.append({
                "Facet": str(facet),
                "ValidLevels": k_valid,
                "Tau2": tau2,
                "TauSource": tau_source,
                "MeanShrinkageFactor": float(np.mean(b)),
                "MaxShrinkageFactor": float(np.max(b)),
                "EffectiveDF": float(np.sum(1.0 - b)),
                "Caution": (
                    "Posterior SE treats tau2 as known; raw estimates remain the primary fitted estimates. "
                    "This is post-hoc shrinkage, not a random-effects model."
                ),
            })
        else:
            report_rows.append({
                "Facet": str(facet),
                "ValidLevels": k_valid,
                "Tau2": np.nan,
                "TauSource": "not estimated",
                "MeanShrinkageFactor": np.nan,
                "MaxShrinkageFactor": np.nan,
                "EffectiveDF": np.nan,
                "Caution": f"At least {int(min_levels)} valid levels with positive SE are required.",
            })
        rows.append(sub)
    if not rows:
        return {
            "available": False,
            "reason": "No non-person facet measures with usable SE were available.",
            "measures": pd.DataFrame(),
            "report": pd.DataFrame(report_rows),
        }
    out = pd.concat(rows, ignore_index=True)
    out["ShrinkageDelta"] = pd.to_numeric(out["ShrunkEstimate"], errors="coerce") - pd.to_numeric(out["Estimate"], errors="coerce")
    return qualify_eb_shrinkage({
        "available": True,
        "reason": "Post-hoc empirical-Bayes shrinkage computed for non-person facet measures.",
        "measures": out,
        "report": pd.DataFrame(report_rows),
    })


def qualify_eb_shrinkage(shrinkage: dict, result: dict | None = None) -> dict:
    """Qualify fresh and saved plug-in calculations without changing their values."""
    out = dict(shrinkage) if isinstance(shrinkage, dict) else {}
    caution = (
        "ShrunkSE is a plug-in normal posterior SD treating the shrinkage variance as known; "
        "its uncertainty and calibration covariance are not propagated. At an estimated zero "
        "variance, zero ShrunkSE is not evidence of exact equality or zero uncertainty. "
        "The prior mean is fixed at zero; changing the logit origin while retaining this prior "
        "changes the shrinkage assumption and can change contrasts. "
        "Point estimates are technical sensitivity calculations; interval inference is unqualified."
    )
    for key in ("measures", "report"):
        frame = out.get(key)
        if not isinstance(frame, pd.DataFrame) or frame.empty:
            continue
        frame = (qualify_measure_uncertainty(frame, result) if key == "measures" and result is not None
                 else frame.copy())
        frame["ShrinkageCI_Status"] = "unqualified_plugin_uncertainty"
        frame["ShrinkagePriorMean"] = 0.0
        frame["ShrinkageReasonCode"] = "stat.eb_shrinkage.uncertainty_unqualified"
        frame["ShrinkageCaution"] = caution
        out[key] = frame
    return out


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


KRIPPENDORFF_ALPHA_METRICS: tuple[str, ...] = ("ordinal", "nominal", "interval")


def _agreement_score_wide(
    obs_df: pd.DataFrame,
    facet_cols: list[str],
    target_facet: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Return unit-by-element score matrix for agreement-style diagnostics."""
    meta: dict[str, object] = {
        "TargetFacet": target_facet,
        "ContextColumns": "",
        "Aggregation": "mean",
        "RowsUsed": 0,
    }
    if obs_df is None or obs_df.empty or target_facet not in obs_df.columns:
        return pd.DataFrame(), meta
    available_facets = [c for c in facet_cols if c in obs_df.columns]
    context_cols = [c for c in available_facets if c != target_facet]
    if not context_cols:
        return pd.DataFrame(), meta

    if "Observed" not in obs_df.columns:
        return pd.DataFrame(), meta
    ctx_series = obs_df[context_cols[0]].astype(str)
    for col in context_cols[1:]:
        ctx_series = ctx_series.str.cat(obs_df[col].astype(str), sep="|")
    df = pd.DataFrame({
        "_context": ctx_series,
        target_facet: obs_df[target_facet],
        "_score_for_agreement": pd.to_numeric(obs_df["Observed"], errors="coerce"),
    })
    if "Weight" in obs_df.columns:
        df["Weight"] = obs_df["Weight"]
    df = df.dropna(subset=["_score_for_agreement", target_facet]).copy()
    if df.empty:
        return pd.DataFrame(), meta

    if "Weight" in df.columns:
        df["Weight"] = pd.to_numeric(df["Weight"], errors="coerce")
        df = df.dropna(subset=["Weight"])
        df = df[df["Weight"] > 0].copy()
        if df.empty:
            return pd.DataFrame(), meta
        df["_weighted_score_for_agreement"] = df["_score_for_agreement"] * df["Weight"]
        df_obs = (
            df.groupby(["_context", target_facet], observed=True, as_index=False)
            .agg(WeightedObserved=("_weighted_score_for_agreement", "sum"), Weight=("Weight", "sum"))
        )
        df_obs["Score"] = np.where(
            df_obs["Weight"].to_numpy(dtype=float) > 0,
            df_obs["WeightedObserved"].to_numpy(dtype=float) / df_obs["Weight"].to_numpy(dtype=float),
            np.nan,
        )
        meta["Aggregation"] = "weighted_mean"
    else:
        df_obs = (
            df.groupby(["_context", target_facet], observed=True, as_index=False)
            .agg(Score=("_score_for_agreement", "mean"))
        )
        meta["Aggregation"] = "mean"

    wide = df_obs.pivot(index="_context", columns=target_facet, values="Score")
    wide = wide.apply(pd.to_numeric, errors="coerce")
    meta["ContextColumns"] = "; ".join(context_cols)
    meta["RowsUsed"] = int(len(df))
    return wide, meta


def _krippendorff_alpha_status(alpha_value: float) -> tuple[str, str]:
    if not np.isfinite(alpha_value):
        return (
            "Unavailable",
            "Alpha could not be estimated from the shared contexts in this run.",
        )
    if alpha_value >= 0.80:
        return (
            "Ready",
            "High agreement for this selected facet under the chosen distance metric.",
        )
    if alpha_value >= 0.67:
        return (
            "Caution",
            "Moderate agreement; usable for some summaries, but inspect pairwise patterns and MFRM fit.",
        )
    return (
        "Review",
        "Low agreement; inspect scoring design, facet fit, category functioning, and possible rater or task drift.",
    )


def _krippendorff_ordinal_distance_matrix(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    categories, counts = np.unique(values, return_counts=True)
    order = np.argsort(categories)
    categories = categories[order]
    counts = counts[order].astype(float)
    n_cat = len(categories)
    if n_cat <= 1:
        return categories, np.zeros((n_cat, n_cat), dtype=float)
    cum = np.concatenate(([0.0], np.cumsum(counts[:-1])))
    midpoints = cum + 0.5 * counts
    dist = (midpoints[:, None] - midpoints[None, :]) ** 2
    return categories, dist


def _krippendorff_alpha_from_wide(
    wide: pd.DataFrame,
    meta: dict[str, object],
    target_facet: str,
    *,
    metric: str = "ordinal",
    include_coverage: bool = True,
) -> dict[str, pd.DataFrame]:
    metric = str(metric or "ordinal").lower()
    if metric not in KRIPPENDORFF_ALPHA_METRICS:
        metric = "ordinal"
    empty_summary = pd.DataFrame(columns=[
        "Facet", "DistanceMetric", "KrippendorffAlpha", "Status", "PlainLanguage",
        "Units", "Elements", "PairableUnits", "RatingsUsed", "MissingCells",
        "UnitsWithSingleRating", "ObservedDisagreement", "ExpectedDisagreement", "ContextColumns",
        "Aggregation", "Caveat",
    ])
    if wide.empty or wide.shape[1] < 2:
        return {"summary": empty_summary, "coverage": pd.DataFrame()}

    mat = wide.to_numpy(dtype=float)
    finite_mask = np.isfinite(mat)
    m_per_unit = np.sum(finite_mask, axis=1)
    units_with_single_rating = int(np.sum(m_per_unit == 1))
    pairable_mask = m_per_unit >= 2
    pairable_units = int(np.sum(pairable_mask))
    mat_pair = mat[pairable_mask]
    mask_pair = finite_mask[pairable_mask]
    m_pair_int = m_per_unit[pairable_mask]
    m_pair = m_pair_int.astype(float)
    values = mat_pair[mask_pair]
    ratings_used = int(values.size)
    units = int(mat.shape[0])
    elements = int(mat.shape[1])
    missing_cells = int(mat.size - int(np.sum(finite_mask)))
    caveat = (
        "Experimental agreement screen. It complements MFRM fit/reliability; "
        "it does not replace model-based severity, bias, or category diagnostics."
    )
    if ratings_used < 2 or pairable_units < 1:
        status, plain = _krippendorff_alpha_status(np.nan)
        summary = pd.DataFrame([{
            "Facet": target_facet,
            "DistanceMetric": metric,
            "KrippendorffAlpha": np.nan,
            "Status": status,
            "PlainLanguage": plain,
            "Units": units,
            "Elements": elements,
            "PairableUnits": pairable_units,
            "UnitsWithSingleRating": units_with_single_rating,
            "RatingsUsed": ratings_used,
            "MissingCells": missing_cells,
            "ObservedDisagreement": np.nan,
            "ExpectedDisagreement": np.nan,
            "ContextColumns": meta.get("ContextColumns", ""),
            "Aggregation": meta.get("Aggregation", "mean"),
            "Caveat": caveat,
        }])
        return {"summary": summary, "coverage": pd.DataFrame()}

    if metric == "interval":
        safe = np.where(mask_pair, mat_pair, 0.0)
        row_mean = np.sum(safe, axis=1) / m_pair
        sq_dev = np.sum(np.where(mask_pair, (mat_pair - row_mean[:, None]) ** 2, 0.0), axis=1)
        observed_num = float(np.sum(2.0 * m_pair * sq_dev / (m_pair - 1.0)))
        mean_value = float(np.mean(values))
        expected_num = float(2.0 * ratings_used * np.sum((values - mean_value) ** 2))
    elif metric == "nominal":
        categories, inv, counts_int = np.unique(values, return_inverse=True, return_counts=True)
        counts = counts_int.astype(float)
        n_cat = len(categories)
        unit_ids = np.repeat(np.arange(pairable_units, dtype=np.int64), m_pair_int)
        pair_key = unit_ids * np.int64(max(n_cat, 1)) + inv.astype(np.int64)
        _, key_idx, cell_counts = np.unique(pair_key, return_index=True, return_counts=True)
        sum_n2 = np.bincount(
            unit_ids[key_idx],
            weights=cell_counts.astype(float) ** 2,
            minlength=pairable_units,
        )
        observed_num = float(np.sum((m_pair ** 2 - sum_n2) / (m_pair - 1.0)))
        expected_num = float(ratings_used ** 2 - np.sum(counts ** 2))
    else:
        # Ordinal squared rank-span distance d^2(c, k) = (mid_c - mid_k)^2 is exact
        # Euclidean squared distance on cumulative marginal midpoints mid_c.
        categories, inv, counts_int = np.unique(values, return_inverse=True, return_counts=True)
        counts = counts_int.astype(float)
        cum = np.concatenate(([0.0], np.cumsum(counts[:-1])))
        midpoints = cum + 0.5 * counts
        z_vals = midpoints[inv]
        unit_ids = np.repeat(np.arange(pairable_units, dtype=int), m_pair_int)
        unit_z_sum = np.bincount(unit_ids, weights=z_vals, minlength=pairable_units)
        unit_z_mean = unit_z_sum / m_pair
        sq_dev_z = np.bincount(
            unit_ids,
            weights=(z_vals - unit_z_mean[unit_ids]) ** 2,
            minlength=pairable_units,
        )
        observed_num = float(np.sum(2.0 * m_pair * sq_dev_z / (m_pair - 1.0)))
        mean_z = float(np.mean(z_vals))
        expected_num = float(2.0 * ratings_used * np.sum((z_vals - mean_z) ** 2))

    observed_disagreement = float(observed_num / ratings_used) if ratings_used else np.nan
    expected_disagreement = (
        float(expected_num / (ratings_used * (ratings_used - 1)))
        if ratings_used > 1 else np.nan
    )
    alpha_value = (
        float(1.0 - observed_disagreement / expected_disagreement)
        if np.isfinite(observed_disagreement)
        and np.isfinite(expected_disagreement)
        and expected_disagreement > 0
        else np.nan
    )
    status, plain = _krippendorff_alpha_status(alpha_value)

    if include_coverage:
        col_counts = np.sum(finite_mask, axis=0)
        col_sums = np.sum(np.where(finite_mask, mat, 0.0), axis=0)
        col_means = np.where(col_counts > 0, col_sums / np.maximum(col_counts, 1), np.nan)
        col_missing = units - col_counts
        coverage_df = pd.DataFrame({
            "Facet": target_facet,
            "Element": list(wide.columns),
            "Ratings": col_counts.astype(int),
            "MeanScore": col_means.astype(float),
            "MissingUnits": col_missing.astype(int),
        })
    else:
        coverage_df = pd.DataFrame()

    summary = pd.DataFrame([{
        "Facet": target_facet,
        "DistanceMetric": metric,
        "KrippendorffAlpha": alpha_value,
        "Status": status,
        "PlainLanguage": plain,
        "Units": units,
        "Elements": elements,
        "PairableUnits": pairable_units,
        "UnitsWithSingleRating": units_with_single_rating,
        "RatingsUsed": ratings_used,
        "MissingCells": missing_cells,
        "ObservedDisagreement": observed_disagreement,
        "ExpectedDisagreement": expected_disagreement,
        "ContextColumns": meta.get("ContextColumns", ""),
        "Aggregation": meta.get("Aggregation", "mean"),
        "Caveat": caveat,
    }])
    return {"summary": summary, "coverage": coverage_df}


def calc_krippendorff_alpha(
    obs_df: pd.DataFrame,
    facet_cols: list[str],
    target_facet: str,
    *,
    metric: str = "ordinal",
) -> dict[str, pd.DataFrame]:
    """Compute an experimental Krippendorff alpha table for one selected facet."""
    wide, meta = _agreement_score_wide(obs_df, facet_cols, target_facet)
    return _krippendorff_alpha_from_wide(
        wide,
        meta,
        target_facet,
        metric=metric,
        include_coverage=True,
    )


def calc_krippendorff_alpha_by_facet(
    obs_df: pd.DataFrame,
    facet_cols: list[str],
    *,
    metrics: tuple[str, ...] = KRIPPENDORFF_ALPHA_METRICS,
) -> pd.DataFrame:
    """Compute experimental Krippendorff alpha summaries for each non-person facet."""
    rows: list[pd.DataFrame] = []
    target_facets = [facet for facet in facet_cols if facet != "Person"]
    for facet in target_facets:
        wide, meta = _agreement_score_wide(obs_df, facet_cols, facet)
        for metric in metrics:
            result = _krippendorff_alpha_from_wide(
                wide,
                meta,
                facet,
                metric=metric,
                include_coverage=False,
            )
            summary = result.get("summary", pd.DataFrame())
            if isinstance(summary, pd.DataFrame) and not summary.empty:
                rows.append(summary)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


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


def expected_score_from_eta(eta, step_cum, rating_min, slope=1.0):
    """Expected score under the slope-aware GPCM / PCM / RSM kernel.

    At ``slope == 1`` the body reduces byte-for-byte to the PCM/RSM Linacre
    expected score; for ``slope > 0`` it implements the GPCM kernel of
    Muraki (1992, Eqs. 2-3, 10) in log-space via ``logsumexp`` for numerical
    stability. A non-finite or non-positive slope returns ``np.nan`` rather
    than silently falling back to slope = 1, so a degenerate fit cannot
    masquerade as a successful one.
    """

    if not np.isfinite(eta) or step_cum is None or len(step_cum) == 0:
        return np.nan
    step_cum = np.asarray(step_cum, dtype=float)
    if step_cum.ndim != 1:
        return np.nan
    if not np.isfinite(slope) or slope <= 0:
        return np.nan
    slope = float(slope)
    k_vals = np.arange(len(step_cum), dtype=float)
    log_num = slope * (eta * k_vals - step_cum)
    if not np.all(np.isfinite(log_num)):
        return np.nan
    log_denom = logsumexp(log_num)
    probs = np.exp(log_num - log_denom)
    return float(rating_min + np.sum(probs * k_vals))


def estimate_eta_from_target(target, step_cum, rating_min, rating_max, slope=1.0):
    if not np.isfinite(target) or step_cum is None or len(step_cum) == 0:
        return np.nan
    if target <= rating_min:
        return -np.inf
    if target >= rating_max:
        return np.inf

    def f(eta):
        return expected_score_from_eta(eta, step_cum, rating_min, slope=slope) - target

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
    measures = qualify_measure_uncertainty(diagnostics.get("measures"), res)
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
        slope_by_step_level = np.array([], dtype=float)
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
        if config["model"] == "GPCM":
            slope_by_step_level = np.asarray(params.get("slopes", []), dtype=float)
        else:
            slope_by_step_level = np.array([], dtype=float)

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
            obs_df.groupby(facet, observed=False)["Observed"]
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
            obs_df.groupby(facet, observed=False)["Observed"]
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
        for level_val, df_lvl in score_source.groupby(facet, observed=False):
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
            fit_sidecars = [
                "DF_Infit", "DF_Outfit",
                "DF_Infit_ENGINE", "DF_Outfit_ENGINE",
                "DF_Infit_FACETS", "DF_Outfit_FACETS",
                "InfitZSTD_ENGINE", "OutfitZSTD_ENGINE",
                "InfitZSTD_FACETS", "OutfitZSTD_FACETS",
                "FitDfMethod", "FitDfFormula", "FitZSTDTransform",
                "FitZSTDCap", "FitDfApplicability", "FitDfWarning",
            ]
            measure_cols = [
                "Level", "Estimate", "SE", "Infit", "Outfit",
                "InfitZSTD", "OutfitZSTD", "PTMEA",
                *[col for col in fit_sidecars if col in meas_tbl.columns],
                *[col for col in ("SE_Method", "SE_Status", "CI_Method", "CI_Status", "UncertaintyCaution")
                  if col in meas_tbl.columns],
            ]
            tbl = tbl.merge(
                meas_tbl[measure_cols],
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

        # Build per-row step-cumulative thresholds and per-row slopes for the
        # slope-aware GPCM Linacre fair-average. The contract here matches the
        # mfrmr R reference implementation:
        #   * Step facet's own rows use that level's own cumulative thresholds.
        #   * Slope facet's own rows use that level's own discrimination.
        #   * Any other facet row uses (a) the column-mean cumulative
        #     thresholds and (b) slope = 1, which is the discrimination
        #     geometric mean under the sum-to-zero log-slope identification.
        # The "slope = 1" choice (rather than the arithmetic mean of the
        # estimated slopes) makes the GPCM fair-average reduce exactly to the
        # PCM Linacre construction on non-slope-facet rows, restoring continuity
        # with the simpler models.
        slope_list = [1.0 for _ in range(len(tbl))]
        if config["model"] in {"PCM", "GPCM"} and config.get("step_facet"):
            step_levels = prep["levels"][config["step_facet"]]
            if facet == config["step_facet"] and len(step_levels) > 0 and len(step_cum_common) > 0:
                step_cum_list = []
                slope_list = []
                for lvl in tbl["Level"].tolist():
                    idx = step_levels.index(lvl) if lvl in step_levels else None
                    if idx is not None and idx < len(step_cum_common):
                        step_cum_list.append(step_cum_common[idx])
                        if (
                            config["model"] == "GPCM"
                            and idx < len(slope_by_step_level)
                            and np.isfinite(slope_by_step_level[idx])
                            and slope_by_step_level[idx] > 0
                        ):
                            slope_list.append(float(slope_by_step_level[idx]))
                        else:
                            slope_list.append(1.0)
                    else:
                        step_cum_list.append(step_cum_mean)
                        slope_list.append(1.0)
            else:
                step_cum_list = [step_cum_mean for _ in range(len(tbl))]
                slope_list = [1.0 for _ in range(len(tbl))]
        else:
            step_cum_list = [step_cum_common for _ in range(len(tbl))]

        tbl["FairM"] = [
            expected_score_from_eta(e, step, rating_min, slope=slope)
            for e, step, slope in zip(eta_m, step_cum_list, slope_list)
        ]
        tbl["FairZ"] = [
            expected_score_from_eta(e, step, rating_min, slope=slope)
            for e, step, slope in zip(eta_z, step_cum_list, slope_list)
        ]

        xtreme_target = np.where(
            tbl["Status"] == "Minimum",
            rating_min + xtreme,
            np.where(tbl["Status"] == "Maximum", rating_max - xtreme, np.nan),
        )
        xtreme_eta = [
            estimate_eta_from_target(t, step, rating_min, rating_max, slope=slope) if xtreme > 0 else np.nan
            for t, step, slope in zip(xtreme_target, step_cum_list, slope_list)
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
        extreme_fit_cols = [
            col for col in (
                "Infit", "Outfit", "InfitZSTD", "OutfitZSTD",
                "InfitZSTD_ENGINE", "OutfitZSTD_ENGINE",
                "InfitZSTD_FACETS", "OutfitZSTD_FACETS",
            )
            if col in tbl.columns
        ]
        tbl.loc[extreme_mask, extreme_fit_cols] = np.nan
        if facet == "Person":
            tbl["PTMEA"] = np.nan

        tbl = tbl.rename(columns={
            "Infit": "InfitMnSq",
            "InfitZSTD": "InfitZStd",
            "Outfit": "OutfitMnSq",
            "OutfitZSTD": "OutfitZStd",
            "PTMEA": "PtMeaCorr",
        })

        display_cols = [
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
            *[col for col in ("SE_Method", "SE_Status", "CI_Method", "CI_Status", "UncertaintyCaution")
              if col in tbl.columns],
            "InfitMnSq",
            "InfitZStd",
            "OutfitMnSq",
            "OutfitZStd",
            "PtMeaCorr",
            *[
                col for col in (
                    "DF_Infit", "DF_Outfit",
                    "DF_Infit_ENGINE", "DF_Outfit_ENGINE",
                    "DF_Infit_FACETS", "DF_Outfit_FACETS",
                    "InfitZSTD_ENGINE", "OutfitZSTD_ENGINE",
                    "InfitZSTD_FACETS", "OutfitZSTD_FACETS",
                    "FitDfMethod", "FitDfFormula", "FitZSTDTransform",
                    "FitZSTDCap", "FitDfApplicability", "FitDfWarning",
                )
                if col in tbl.columns
            ],
            "Anchor",
            "Status",
            "Level",
            "TotalCountAll",
        ]
        tbl = tbl[display_cols]

        if omit_unobserved:
            tbl = tbl[tbl["TotalCountAll"] > 0]

        tbl = tbl.drop(columns=["TotalCountAll"]).sort_values(["Measure", "TotalCount"], ascending=[False, False])
        out[facet] = tbl

    return out


def _invert_information_matrix(info):
    """Invert a symmetric observed-information matrix with eigenvalue regularization.

    Returns a tuple ``(cov, regularized, rank)`` where ``cov`` is the
    inverse (or ``None`` if inversion failed), ``regularized`` is ``True``
    when any eigenvalue was at or below the tolerance and had to be floored,
    and ``rank`` is the count of eigenvalues strictly above the tolerance.
    The eigenvalue floor is ``max|lambda| * sqrt(eps)`` (``sqrt(eps)``
    for a zero matrix); this is not NumPy's default matrix-rank tolerance.
    """
    if info is None:
        return None, False, 0
    info = np.asarray(info, dtype=float)
    if info.size == 0 or info.ndim != 2 or info.shape[0] != info.shape[1]:
        return None, False, 0
    if not np.all(np.isfinite(info)):
        return None, False, 0
    info_sym = (info + info.T) / 2.0
    try:
        eigvals, eigvecs = np.linalg.eigh(info_sym)
    except np.linalg.LinAlgError:
        return None, False, 0
    if eigvals.size == 0:
        return None, False, 0
    max_val = float(np.max(np.abs(eigvals)))
    eps = float(np.finfo(float).eps)
    tol = max_val * np.sqrt(eps) if max_val > 0 else np.sqrt(eps)
    rank = int(np.sum(eigvals > tol))
    regularized = bool(np.any(eigvals <= tol))
    safe_eigvals = np.maximum(eigvals, tol)
    cov = eigvecs @ np.diag(1.0 / safe_eigvals) @ eigvecs.T
    cov = (cov + cov.T) / 2.0
    return cov, regularized, rank


def _information_matrix_spectrum(info) -> dict[str, object]:
    """Numerical spectrum diagnostics for an observed-information matrix."""
    out = {
        "HessianFinite": False,
        "ParamCount": 0,
        "Rank": 0,
        "RankDeficiency": np.nan,
        "MinEigenvalue": np.nan,
        "MaxEigenvalue": np.nan,
        "MinPositiveEigenvalue": np.nan,
        "EigFloorTolerance": np.nan,
        "RegularizedEigenvalues": np.nan,
        "ConditionNumber": np.nan,
    }
    if info is None:
        return out
    arr = np.asarray(info, dtype=float)
    if arr.size == 0 or arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        return out
    out["ParamCount"] = int(arr.shape[0])
    if not np.all(np.isfinite(arr)):
        return out
    info_sym = (arr + arr.T) / 2.0
    try:
        eigvals = np.linalg.eigvalsh(info_sym)
    except np.linalg.LinAlgError:
        return out
    if eigvals.size == 0:
        return out
    max_abs = float(np.max(np.abs(eigvals)))
    eps = float(np.finfo(float).eps)
    tol = max_abs * np.sqrt(eps) if max_abs > 0 else np.sqrt(eps)
    positive = eigvals[eigvals > tol]
    rank = int(positive.size)
    min_positive = float(np.min(positive)) if positive.size else np.nan
    max_positive = float(np.max(positive)) if positive.size else np.nan
    condition = (
        float(max_positive / min_positive)
        if np.isfinite(max_positive) and np.isfinite(min_positive) and min_positive > 0
        else np.inf
    )
    out.update({
        "HessianFinite": True,
        "Rank": rank,
        "RankDeficiency": int(arr.shape[0] - rank),
        "MinEigenvalue": float(np.min(eigvals)),
        "MaxEigenvalue": float(np.max(eigvals)),
        "MinPositiveEigenvalue": min_positive,
        "EigFloorTolerance": float(tol),
        "RegularizedEigenvalues": int(np.sum(eigvals <= tol)),
        "ConditionNumber": condition,
    })
    return out


def _mml_covariance_scope_issue(config):
    """Reasons the current structural Hessian is not a sampling covariance."""
    if str(config.get("method", "")).upper() != "MML":
        return ""
    reasons = []
    if config.get("estimate_population_sd"):
        reasons.append(
            "The joint finite-quadrature information is retained for technical audit; "
            "quadrature accuracy and sampling uncertainty are not yet qualified."
            if config.get("mml_engine") == "joint" else
            "The structural Hessian holds the estimated population SD fixed; "
            "joint nuisance-adjusted covariance is not implemented."
        )
    if config.get("facet_regularization_enabled") or (config.get("facet_regularization") or {}).get("enabled"):
        reasons.append(
            "The Hessian includes a facet regularization penalty; inverse "
            "penalized curvature is not an unpenalized sampling covariance. "
            "Penalized sampling uncertainty is not implemented."
        )
    return " ".join(reasons)


def build_mml_covariance_audit(covariance: dict | None, result: dict | None = None) -> pd.DataFrame:
    """Machine-readable audit for MML observed-information covariance quality."""
    cov_info = covariance if isinstance(covariance, dict) else {}
    config = (result or {}).get("config", {}) if isinstance(result, dict) else {}
    scope_issue = _mml_covariance_scope_issue(config)
    if scope_issue:
        cov_info = {"status": "unsupported_covariance_scope", "detail": scope_issue}
    par = _get_opt_par(result or {}) if isinstance(result, dict) else None
    hessian = cov_info.get("hessian")
    spectrum = _information_matrix_spectrum(hessian)
    param_count = int(spectrum.get("ParamCount") or 0)
    if param_count == 0 and par is not None:
        param_count = int(np.asarray(par).size)
    status = str(cov_info.get("status", "not_available"))
    detail = str(cov_info.get("detail", "MML covariance audit was not available."))
    cov = cov_info.get("cov")
    cov_arr = np.asarray(cov, dtype=float) if cov is not None else np.asarray([])
    cov_finite = bool(cov_arr.size and cov_arr.ndim == 2 and np.all(np.isfinite(cov_arr)))
    cov_diag = np.diag(cov_arr) if cov_finite and cov_arr.shape[0] == cov_arr.shape[1] else np.asarray([])
    condition = float(spectrum.get("ConditionNumber", np.nan))
    rank_def = spectrum.get("RankDeficiency", np.nan)
    regularized_eigs = spectrum.get("RegularizedEigenvalues", np.nan)
    resolved_rank = int(cov_info.get("rank", spectrum.get("Rank", 0)) or 0)
    output_qualification = _output_qualification.covariance_qualification(
        status=status,
        rank=resolved_rank,
        param_count=param_count,
    )

    if status == "ok":
        claim_status = "Ready"
        interpretation = "Observed-information covariance was inverted without eigenvalue regularization."
        action = "Report SE/CI status and retain the covariance audit with the analysis archive."
    elif status == "regularized":
        claim_status = "Report with caveat"
        interpretation = "Observed-information covariance required eigenvalue regularization; SEs are regularization-aware large-sample approximations."
        action = "Report the regularized status, condition number, and sensitivity caveat before treating SE/CI as calibration evidence."
    elif status == "skipped_large_parameter_count":
        claim_status = "Report with caveat"
        interpretation = "Routine observed-information covariance was skipped by the parameter-count gate."
        action = "Use conditional SEs for screening, or rerun a smaller/targeted validation if structural SEs are required."
    elif status == "not_applicable":
        claim_status = "Not applicable"
        interpretation = "Observed-information covariance is only defined for the current MML path."
        action = "Use the method-specific SE/CI status columns instead of claiming structural MML covariance."
    else:
        claim_status = "Do not claim"
        interpretation = "Observed-information covariance was unavailable or failed diagnostics."
        action = "Do not claim structural covariance-based SE/CI; report fallback SE_Status and investigate the covariance failure."

    if np.isfinite(condition) and condition > 1e8 and claim_status == "Ready":
        claim_status = "Report with caveat"
        interpretation += " The Hessian condition number is high."
        action = "Report the high condition number and consider a sensitivity check or simpler parameterization."
    if np.isfinite(rank_def) and float(rank_def) > 0 and claim_status == "Ready":
        claim_status = "Report with caveat"

    # S0 fail-closed boundary: numerical covariance values remain available
    # for technical audit, but neither eigenvalue flooring nor full numerical
    # rank is sufficient evidence for confirmatory SE/CI claims. Coverage is a
    # separate G2 requirement.
    if output_qualification.status == _output_qualification.OutputUse.WITHHELD:
        claim_status = "Do not claim"
        interpretation += " Public inferential use is withheld by the S0 output-qualification gate."
        action = output_qualification.user_action
    elif output_qualification.status == _output_qualification.OutputUse.TECHNICAL_ONLY:
        claim_status = "Technical only"
        interpretation += " Numerical full rank does not yet establish interval coverage."
        action = output_qualification.user_action
    elif output_qualification.status == _output_qualification.OutputUse.NOT_APPLICABLE:
        claim_status = "Not applicable"
        action = output_qualification.user_action

    row = {
        "Area": "MML observed-information covariance",
        "Status": status,
        "ClaimStatus": claim_status,
        "Method": "numerical Jacobian of analytical MML gradient; eigendecomposition inverse",
        "Model": str(config.get("model", "")),
        "FitMethod": str(config.get("method", "")),
        "ParamCount": param_count,
        "AutoParamLimit": int(MML_COVARIANCE_AUTO_MAX_PARAMS),
        "Rank": resolved_rank,
        "RankDeficiency": rank_def,
        "HessianFinite": bool(spectrum.get("HessianFinite", False)),
        "MinEigenvalue": spectrum.get("MinEigenvalue", np.nan),
        "MaxEigenvalue": spectrum.get("MaxEigenvalue", np.nan),
        "MinPositiveEigenvalue": spectrum.get("MinPositiveEigenvalue", np.nan),
        "EigFloorTolerance": spectrum.get("EigFloorTolerance", np.nan),
        "RegularizedEigenvalues": regularized_eigs,
        "ConditionNumber": condition,
        "Regularized": bool(cov_info.get("regularized", False)),
        "CovarianceFinite": cov_finite,
        "CovarianceMinDiagonal": float(np.nanmin(cov_diag)) if cov_diag.size else np.nan,
        "CovarianceMaxDiagonal": float(np.nanmax(cov_diag)) if cov_diag.size else np.nan,
        "Detail": detail,
        "Interpretation": interpretation,
        "RecommendedAction": action,
        "QualificationSchemaVersion": output_qualification.schema_version,
        "QualificationStatus": output_qualification.status.value,
        "ReasonCode": output_qualification.reason_code.value,
        "RoadmapIssueID": output_qualification.roadmap_issue_id,
        "NextGate": output_qualification.next_gate,
        "RawTechnicalExportAllowed": output_qualification.raw_technical_export_allowed,
        "PublicConclusionAllowed": output_qualification.public_conclusion_allowed,
    }
    return pd.DataFrame([row])


def _compute_mml_observed_information(
    par, idx, config, sizes, quad, rel_step=1e-5, *, include_log_sigma=False,
):
    """Curvature of the MML objective in the supplied free coordinates.

    Computed as the numerical Jacobian of the analytical gradient
    returned by ``mfrm_loglik_mml_value_grad`` via central differences at
    relative step size ``rel_step``. Symmetrised before return; downstream
    inversion is responsible for regularising near-singular spectra. An
    enabled facet penalty is included, so the public covariance helper
    restricts use to unpenalized fits with a prespecified population SD.

    Development RSM/PCM/GPCM callers can include log(SD) as the last coordinate.
    In that case quad must be the rule at the supplied SD. Its nodes move
    with SD, while standardized nodes and weights remain fixed. The full
    matrix includes structural/log-SD cross derivatives; it does not itself
    certify stationarity, integration accuracy, or sampling covariance.
    """
    par = np.asarray(par, dtype=float)
    if not np.isfinite(rel_step) or rel_step <= 0:
        raise ValueError("rel_step must be finite and positive")
    if include_log_sigma:
        if config["model"] not in {"RSM", "PCM", "GPCM"}:
            raise ValueError("Joint log-SD information requires a native response model")
        if par.ndim != 1 or par.size != sum(sizes.values()) + 1 or not np.isfinite(par).all():
            raise ValueError("Expected finite structural coordinates followed by log(SD)")
        sigma = float(np.exp(par[-1]))
        if not np.isfinite(sigma) or sigma <= 0:
            raise ValueError("Population SD must be finite and positive")
        standard_nodes = np.asarray(quad["nodes"], dtype=float) / sigma

    def gradient(point):
        if include_log_sigma:
            moving_quad = {**quad, "nodes": standard_nodes * np.exp(point[-1])}
            return mfrm_loglik_mml_value_grad(
                point[:-1], idx, config, sizes, moving_quad, include_log_sigma=True,
            )[1]
        return mfrm_loglik_mml_value_grad(point, idx, config, sizes, quad)[1]

    P = par.size
    if P == 0:
        return np.zeros((0, 0), dtype=float)
    H = np.zeros((P, P), dtype=float)
    for j in range(P):
        h = rel_step * max(1.0, abs(float(par[j])))
        par_p = par.copy()
        par_m = par.copy()
        par_p[j] += h
        par_m[j] -= h
        grad_p = gradient(par_p)
        grad_m = gradient(par_m)
        H[:, j] = (np.asarray(grad_p, dtype=float) - np.asarray(grad_m, dtype=float)) / (2.0 * h)
    return (H + H.T) / 2.0


def _build_param_slices(sizes):
    """Return a dict mapping each parameter section name to its slice into ``par``."""
    out = {}
    cursor = 0
    for name, count in sizes.items():
        out[name] = slice(cursor, cursor + int(count))
        cursor += int(count)
    return out


def _get_opt_par(res):
    """Robustly extract the optimisation parameter vector from a fit result."""
    opt = res.get("opt") if isinstance(res, dict) else None
    if opt is None:
        return None
    if hasattr(opt, "x"):
        return np.asarray(getattr(opt, "x"), dtype=float)
    if isinstance(opt, dict):
        for key in ("x", "par"):
            if key in opt:
                return np.asarray(opt[key], dtype=float)
    return None


def compute_mml_parameter_covariance(res, rel_step=1e-5):
    """Compute observed-information covariance for unpenalized fixed-SD MML.

    Returns a dict with the keys ``cov`` (P x P ndarray or ``None``),
    ``hessian`` (P x P ndarray or ``None``), ``sizes`` (the OrderedDict from
    ``build_param_sizes``), ``param_slices`` (section-name -> slice),
    ``status`` (``"ok" | "regularized" | "not_applicable" | "fallback"``
    or ``"unsupported_covariance_scope"``),
    ``detail`` (human-readable description), ``regularized`` (bool),
    ``rank`` (int), and a one-row ``audit`` table with Hessian spectrum
    and covariance quality diagnostics.
    """
    config = (res or {}).get("config", {}) or {}
    try:
        sizes = build_param_sizes(config) if config.get("facet_names") else OrderedDict()
    except Exception:
        sizes = OrderedDict()
    param_slices = _build_param_slices(sizes)

    def envelope(status, detail, cov=None, hessian=None, regularized=False, rank=0):
        out = {
            "cov": cov,
            "hessian": hessian,
            "sizes": sizes,
            "param_slices": param_slices,
            "status": status,
            "detail": detail,
            "regularized": regularized,
            "rank": rank,
        }
        out["audit"] = build_mml_covariance_audit(out, res)
        return out

    if config.get("method") != "MML":
        return envelope(
            "not_applicable",
            "MML observed-information covariance is only computed for MML fits.",
        )

    scope_issue = _mml_covariance_scope_issue(config)
    if scope_issue:
        return envelope("unsupported_covariance_scope", scope_issue)

    par = _get_opt_par(res)
    if par is None or par.size == 0:
        return envelope(
            "fallback",
            "Optimised parameter vector was not available; observed-information covariance was not computed.",
        )

    try:
        idx = build_indices(
            res["prep"],
            step_facet=config.get("step_facet"),
            slope_facet=config.get("slope_facet"),
        )
    except Exception as exc:
        return envelope(
            "fallback",
            f"Could not rebuild observation indices for covariance: {exc!s}",
        )

    quad_points = int(config.get("quad_points") or 15)
    try:
        quad = gauss_hermite_normal(max(1, quad_points), sd=get_population_prior_sd(config))
    except Exception as exc:
        return envelope(
            "fallback",
            f"Could not build Gauss-Hermite quadrature for covariance: {exc!s}",
        )

    try:
        hess = _compute_mml_observed_information(par, idx, config, sizes, quad, rel_step=rel_step)
    except Exception as exc:
        return envelope(
            "fallback",
            f"Observed-information Hessian computation raised: {exc!s}",
        )

    if not np.all(np.isfinite(hess)):
        return envelope(
            "fallback",
            "Observed-information matrix contained non-finite entries; covariance was not computed.",
            hessian=hess,
        )

    cov, regularized, rank = _invert_information_matrix(hess)
    if cov is None:
        return envelope(
            "fallback",
            "Observed-information matrix could not be inverted; covariance was not computed.",
            hessian=hess,
        )

    if regularized:
        detail = (
            "MML parameter covariance from observed information; a near-singular "
            "Hessian was regularized during inversion."
        )
        status = "regularized"
    else:
        detail = "MML parameter covariance from observed information."
        status = "ok"

    return envelope(status, detail, cov=cov, hessian=hess, regularized=regularized, rank=rank)


def _finite_difference_gradient(fn, par, rel_step=1e-5):
    """Central-difference gradient of a scalar ``fn(par)``.

    Falls back to one-sided differences if the symmetric perturbation
    returns a non-finite value at either end (e.g. the per-row fair-average
    crosses a rating bound). Returns NaN for entries that cannot be
    estimated at all.
    """
    par = np.asarray(par, dtype=float)
    P = par.size
    if P == 0:
        return np.zeros(0, dtype=float)
    grad = np.full(P, np.nan, dtype=float)
    try:
        base_val = float(fn(par))
    except Exception:
        base_val = float("nan")
    for j in range(P):
        h = rel_step * max(1.0, abs(float(par[j])))
        par_p = par.copy()
        par_m = par.copy()
        par_p[j] += h
        par_m[j] -= h
        try:
            f_p = float(fn(par_p))
        except Exception:
            f_p = float("nan")
        try:
            f_m = float(fn(par_m))
        except Exception:
            f_m = float("nan")
        if np.isfinite(f_p) and np.isfinite(f_m):
            grad[j] = (f_p - f_m) / (2.0 * h)
        elif np.isfinite(f_p) and np.isfinite(base_val):
            grad[j] = (f_p - base_val) / h
        elif np.isfinite(f_m) and np.isfinite(base_val):
            grad[j] = (base_val - f_m) / h
    return grad


def _mml_covariance_for_diagnostics(res, max_params=MML_COVARIANCE_AUTO_MAX_PARAMS):
    """Observed-information covariance with an automatic size gate for diagnostics.

    ``compute_mml_parameter_covariance`` remains the public helper. This
    wrapper is used during normal diagnostics so large MML runs do not silently
    become much slower; when skipped, downstream tables keep conditional SEs and
    carry an explicit ``skipped_large_parameter_count`` status.
    """
    config = (res or {}).get("config", {}) or {}
    try:
        sizes = build_param_sizes(config) if config.get("facet_names") else OrderedDict()
    except Exception:
        sizes = OrderedDict()
    param_slices = _build_param_slices(sizes)

    def envelope(status, detail):
        out = {
            "cov": None,
            "hessian": None,
            "sizes": sizes,
            "param_slices": param_slices,
            "status": status,
            "detail": detail,
            "regularized": False,
            "rank": 0,
        }
        out["audit"] = build_mml_covariance_audit(out, res)
        return out

    if config.get("method") != "MML":
        return envelope(
            "not_applicable",
            "MML observed-information covariance is only computed for MML fits.",
        )
    scope_issue = _mml_covariance_scope_issue(config)
    if scope_issue:
        return envelope("unsupported_covariance_scope", scope_issue)
    par = _get_opt_par(res)
    if par is not None and max_params is not None and int(par.size) > int(max_params):
        return envelope(
            "skipped_large_parameter_count",
            (
                f"MML observed-information covariance was skipped for routine "
                f"diagnostics because the parameter vector has {int(par.size)} "
                f"parameters (auto limit {int(max_params)}). The measure table "
                "keeps conditional information SEs. Reduce the supported design "
                "or rerun within the covariance limit when structural SEs are "
                "required; otherwise report the conditional-SE scope and mark "
                "structural covariance as not assessable."
            ),
        )
    return compute_mml_parameter_covariance(res)


def _compute_measure_value_from_par(res, par, facet, level):
    """Return a non-person facet measure as a scalar function of ``par``."""
    config = (res or {}).get("config", {}) or {}
    if config.get("method") != "MML" or str(facet) == "Person":
        return float("nan")
    prep = (res or {}).get("prep", {}) or {}
    facet = str(facet)
    level = str(level)
    if facet not in list(config.get("facet_names", [])):
        return float("nan")
    try:
        sizes = build_param_sizes(config)
        params = expand_params(np.asarray(par, dtype=float), sizes, config)
    except Exception:
        return float("nan")
    levels = [str(x) for x in (prep.get("levels", {}).get(facet) or [])]
    if level not in levels:
        return float("nan")
    idx = levels.index(level)
    values = np.asarray(params.get("facets", {}).get(facet, np.array([])), dtype=float)
    if idx >= values.size:
        return float("nan")
    value = float(values[idx])
    return value if np.isfinite(value) else float("nan")


def compute_mml_structural_measure_se_table(res, covariance=None, ci_level=0.95):
    """Delta-method SEs for MML non-person facet estimates.

    Person EAP estimates are posterior summaries and are intentionally excluded:
    they are not coordinates of the marginal-likelihood parameter vector. The
    returned table is keyed by ``Facet`` / ``Level`` and is merged into the main
    measure table by ``annotate_measure_uncertainty`` when finite SEs are
    available.
    """
    config = (res or {}).get("config", {}) or {}
    if config.get("method") != "MML":
        return pd.DataFrame()
    scope_issue = _mml_covariance_scope_issue(config)
    if scope_issue:
        covariance = {"status": "unsupported_covariance_scope", "detail": scope_issue}
    elif covariance is None:
        covariance = _mml_covariance_for_diagnostics(res)
    prep = (res or {}).get("prep", {}) or {}
    par = _get_opt_par(res)
    status = str((covariance or {}).get("status", "fallback"))
    detail = str((covariance or {}).get("detail", "MML covariance unavailable."))
    cov = (covariance or {}).get("cov")
    rows: list[dict[str, object]] = []
    z = float(_norm.ppf((1.0 + float(ci_level)) / 2.0)) if 0.0 < float(ci_level) < 1.0 else 1.959963984540054

    for facet in list(config.get("facet_names", [])):
        levels = [str(x) for x in (prep.get("levels", {}).get(facet) or [])]
        for level in levels:
            row = {
                "Facet": str(facet),
                "Level": str(level),
                "StructuralSE": np.nan,
                "StructuralCI_Lower": np.nan,
                "StructuralCI_Upper": np.nan,
                "StructuralCI_Level": float(ci_level),
                "StructuralSE_Status": status,
                "StructuralSE_Method": "MML observed-information delta method",
                "StructuralSE_Detail": detail,
            }
            if status in {"ok", "regularized"} and cov is not None and par is not None:
                cov_arr = np.asarray(cov, dtype=float)
                if cov_arr.ndim == 2 and cov_arr.shape[0] == cov_arr.shape[1] == par.size:
                    fn = lambda p, f=facet, lv=level: _compute_measure_value_from_par(res, p, f, lv)
                    grad = _finite_difference_gradient(fn, par)
                    if grad.size == cov_arr.shape[0] and np.all(np.isfinite(grad)):
                        var = float(grad @ cov_arr @ grad)
                        if -1e-10 < var < 0:
                            var = 0.0
                        estimate = _compute_measure_value_from_par(res, par, facet, level)
                        if np.isfinite(var) and var >= 0 and np.isfinite(estimate):
                            se = float(np.sqrt(var))
                            row.update({
                                "StructuralSE": se,
                                "StructuralCI_Lower": estimate - z * se,
                                "StructuralCI_Upper": estimate + z * se,
                                "StructuralSE_Status": status,
                                "StructuralSE_Detail": detail,
                            })
                        else:
                            row["StructuralSE_Status"] = "not_available"
                            row["StructuralSE_Detail"] = "Delta-method variance was not finite and non-negative."
                    else:
                        row["StructuralSE_Status"] = "not_available"
                        row["StructuralSE_Detail"] = "Finite-difference gradient was not available."
            rows.append(row)
    return pd.DataFrame(rows)


def _compute_step_value_from_par(res, par, step_facet_level, step_index):
    """Return a reported step / threshold estimate as a scalar function of ``par``.

    Mirrors :func:`_compute_measure_value_from_par` for the rating-scale step
    thresholds. The reported step is the centred value (``center_sum_zero``)
    that the fit surfaces, so differentiating it through ``par`` carries the
    sum-to-zero identification automatically -- no explicit Jacobian is needed.
    """
    config = (res or {}).get("config", {}) or {}
    if config.get("method") != "MML":
        return float("nan")
    n_cat = int(config.get("n_cat", 0) or 0)
    n_steps = max(n_cat - 1, 0)
    step_index = int(step_index)
    if step_index < 0 or step_index >= n_steps:
        return float("nan")
    try:
        sizes = build_param_sizes(config)
        params = expand_params(np.asarray(par, dtype=float), sizes, config)
    except Exception:
        return float("nan")
    if config.get("model") == "RSM":
        steps = np.asarray(params.get("steps"), dtype=float)
        if steps.ndim != 1 or step_index >= steps.size:
            return float("nan")
        value = float(steps[step_index])
    else:
        prep = (res or {}).get("prep", {}) or {}
        step_facet = config.get("step_facet")
        levels = [str(x) for x in (prep.get("levels", {}).get(step_facet) or [])]
        sf = str(step_facet_level)
        steps_mat = params.get("steps_mat")
        if sf not in levels or steps_mat is None:
            return float("nan")
        steps_mat = np.asarray(steps_mat, dtype=float)
        level_idx = levels.index(sf)
        if steps_mat.ndim != 2 or level_idx >= steps_mat.shape[0] or step_index >= steps_mat.shape[1]:
            return float("nan")
        value = float(steps_mat[level_idx, step_index])
    return value if np.isfinite(value) else float("nan")


def compute_mml_step_se_table(res, covariance=None, ci_level=0.95):
    """Delta-method SEs for MML rating-scale step / threshold estimates.

    Uses the same MML observed-information delta method as
    :func:`compute_mml_structural_measure_se_table`, applied to the step
    thresholds. Returns a table keyed by ``Step`` (RSM) or ``StepFacet`` /
    ``Step`` (PCM / GPCM) so it can be merged onto the fitted steps table.
    Empty for non-MML fits; SE columns are ``NaN`` with a status when the
    marginal-likelihood covariance is unavailable.
    """
    config = (res or {}).get("config", {}) or {}
    if config.get("method") != "MML":
        return pd.DataFrame()
    n_cat = int(config.get("n_cat", 0) or 0)
    n_steps = max(n_cat - 1, 0)
    if n_steps <= 0:
        return pd.DataFrame()
    model = config.get("model")
    prep = (res or {}).get("prep", {}) or {}
    scope_issue = _mml_covariance_scope_issue(config)
    if scope_issue:
        covariance = {"status": "unsupported_covariance_scope", "detail": scope_issue}
    elif covariance is None:
        covariance = _mml_covariance_for_diagnostics(res)
    par = _get_opt_par(res)
    status = str((covariance or {}).get("status", "fallback"))
    detail = str((covariance or {}).get("detail", "MML covariance unavailable."))
    cov = (covariance or {}).get("cov")
    z = float(_norm.ppf((1.0 + float(ci_level)) / 2.0)) if 0.0 < float(ci_level) < 1.0 else 1.959963984540054

    targets: list[tuple[object, int, dict[str, object]]] = []
    if model == "RSM":
        for i in range(n_steps):
            targets.append((None, i, {"Step": f"Step_{i + 1}"}))
    else:
        step_facet = config.get("step_facet")
        for sf in [str(x) for x in (prep.get("levels", {}).get(step_facet) or [])]:
            for i in range(n_steps):
                targets.append((sf, i, {"StepFacet": sf, "Step": f"Step_{i + 1}"}))

    cov_arr = np.asarray(cov, dtype=float) if cov is not None else None
    par_arr = np.asarray(par, dtype=float) if par is not None else None
    usable = (
        status in {"ok", "regularized"}
        and cov_arr is not None
        and par_arr is not None
        and cov_arr.ndim == 2
        and cov_arr.shape[0] == cov_arr.shape[1] == par_arr.size
    )

    rows: list[dict[str, object]] = []
    for sf, i, key in targets:
        row = dict(key)
        row.update({
            "SE": np.nan,
            "CI_Lower": np.nan,
            "CI_Upper": np.nan,
            "CI_Level": float(ci_level),
            "SE_Method": "MML observed-information delta method",
            "CI_Method": "Normal approximation from observed-information delta-method SE; conditional on model and calibration scope",
            "SE_Status": status,
            "SE_CovarianceStatus": status,
            "SE_Detail": detail,
        })
        if usable:
            fn = lambda p, _sf=sf, _i=i: _compute_step_value_from_par(res, p, _sf, _i)
            grad = _finite_difference_gradient(fn, par_arr)
            if grad.size == cov_arr.shape[0] and np.all(np.isfinite(grad)):
                var = float(grad @ cov_arr @ grad)
                if -1e-10 < var < 0:
                    var = 0.0
                estimate = _compute_step_value_from_par(res, par_arr, sf, i)
                if np.isfinite(var) and var >= 0 and np.isfinite(estimate):
                    se = float(np.sqrt(var))
                    row.update({
                        "SE": se,
                        "CI_Lower": estimate - z * se,
                        "CI_Upper": estimate + z * se,
                    })
                else:
                    row["SE_Status"] = "not_available"
                    row["SE_Detail"] = "Delta-method variance was not finite and non-negative."
            else:
                row["SE_Status"] = "not_available"
                row["SE_Detail"] = "Finite-difference gradient was not available."
        rows.append(row)
    return qualify_measure_uncertainty(pd.DataFrame(rows), res)


def annotate_measure_uncertainty(measures, res, obs_df=None, *, ci_level=0.95):
    """Attach SE/CI provenance columns to the measure table."""
    if not isinstance(measures, pd.DataFrame) or measures.empty:
        return measures, {
            "covariance": {"status": "not_available", "detail": "Measure table is empty."},
            "structural_se": pd.DataFrame(),
            "summary": pd.DataFrame(),
        }
    out = measures.copy()
    config = (res or {}).get("config", {}) or {}
    method = str(config.get("method", ""))
    z = float(_norm.ppf((1.0 + float(ci_level)) / 2.0)) if 0.0 < float(ci_level) < 1.0 else 1.959963984540054

    out["SE_Method"] = "conditional information approximation: 1/sqrt(sum w*a^2*Var), other fitted parameters held fixed"
    out["SE_Status"] = "conditional_approximation"
    out["SE_Detail"] = "Does not propagate uncertainty from other estimated person, facet, step, slope, or anchor parameters."
    out["CI_Level"] = float(ci_level)
    out["CI_Method"] = "Wald normal CI from the reported SE"
    out["CI_Status"] = "conditional_approximation"
    out["UncertaintyCaution"] = "Use for screening unless the method/status column says structural covariance was used."

    if method == "MML" and "Facet" in out.columns:
        person_mask = out["Facet"].astype(str) == "Person"
        out.loc[person_mask, "SE_Method"] = "MML latent-ability posterior SD from quadrature"
        out.loc[person_mask, "SE_Status"] = "posterior_eap_sd"
        out.loc[person_mask, "SE_Detail"] = (
            "Person SE is the latent-ability posterior SD conditional on fitted "
            "facet, step, slope, and population parameters. It is not the sampling "
            "SE of the EAP estimator."
        )
        out.loc[person_mask, "CI_Method"] = "Normal approximation: EAP +/- z * conditional posterior SD"
        out.loc[person_mask, "CI_Status"] = "posterior_interval_approximation"
        out.loc[person_mask, "UncertaintyCaution"] = (
            "The interval uses a normal approximation to the conditional posterior, "
            "not posterior quantiles; calibration uncertainty is not propagated."
        )

        covariance = _mml_covariance_for_diagnostics(res)
        structural = compute_mml_structural_measure_se_table(
            res,
            covariance=covariance,
            ci_level=ci_level,
        )
        if isinstance(structural, pd.DataFrame) and not structural.empty:
            out = out.merge(structural, on=["Facet", "Level"], how="left")
            structural_ok = (
                out.get("StructuralSE").notna()
                & out.get("StructuralSE_Status").astype(str).isin({"ok", "regularized"})
            )
            out.loc[structural_ok, "SE"] = out.loc[structural_ok, "StructuralSE"]
            out.loc[structural_ok, "CI_Lower"] = out.loc[structural_ok, "StructuralCI_Lower"]
            out.loc[structural_ok, "CI_Upper"] = out.loc[structural_ok, "StructuralCI_Upper"]
            out.loc[structural_ok, "SE_Method"] = out.loc[structural_ok, "StructuralSE_Method"]
            out.loc[structural_ok, "SE_Status"] = out.loc[structural_ok, "StructuralSE_Status"]
            out.loc[structural_ok, "SE_Detail"] = out.loc[structural_ok, "StructuralSE_Detail"]
            out.loc[structural_ok, "CI_Method"] = "Wald normal CI using MML observed-information delta-method SE"
            out.loc[structural_ok, "CI_Status"] = out.loc[structural_ok, "StructuralSE_Status"]
            out.loc[structural_ok, "UncertaintyCaution"] = (
                "Non-person MML facet interval propagates structural covariance "
                "for fitted model parameters; it remains a large-sample Wald CI."
            )
            helper_cols = [
                "StructuralSE", "StructuralCI_Lower", "StructuralCI_Upper",
                "StructuralCI_Level", "StructuralSE_Status",
                "StructuralSE_Method", "StructuralSE_Detail",
            ]
            out = out.drop(columns=[c for c in helper_cols if c in out.columns])
    else:
        covariance = {"status": "not_applicable", "detail": "MML covariance is not applicable to this run."}
        structural = pd.DataFrame()

    covariance_audit = (
        covariance.get("audit")
        if isinstance(covariance, dict) and isinstance(covariance.get("audit"), pd.DataFrame)
        else build_mml_covariance_audit(covariance, res)
    )
    if isinstance(covariance_audit, pd.DataFrame) and not covariance_audit.empty:
        ca = covariance_audit.iloc[0]
        out["SE_CovarianceStatus"] = str(ca.get("Status", "not_available"))
        out["SE_CovarianceClaimStatus"] = str(ca.get("ClaimStatus", "Do not claim"))
        out["SE_CovarianceConditionNumber"] = ca.get("ConditionNumber", np.nan)
        out["SE_CovarianceRank"] = ca.get("Rank", np.nan)
        out["SE_CovarianceRankDeficiency"] = ca.get("RankDeficiency", np.nan)
        out["SE_CovarianceRegularizedEigenvalues"] = ca.get("RegularizedEigenvalues", np.nan)
        out["SE_CovarianceDetail"] = str(ca.get("Detail", ""))
    else:
        out["SE_CovarianceStatus"] = "not_available"
        out["SE_CovarianceClaimStatus"] = "Do not claim"
        out["SE_CovarianceConditionNumber"] = np.nan
        out["SE_CovarianceRank"] = np.nan
        out["SE_CovarianceRankDeficiency"] = np.nan
        out["SE_CovarianceRegularizedEigenvalues"] = np.nan
        out["SE_CovarianceDetail"] = ""

    est = pd.to_numeric(out.get("Estimate", pd.Series(dtype=float)), errors="coerce")
    se = pd.to_numeric(out.get("SE", pd.Series(dtype=float)), errors="coerce")
    needs_ci = out["CI_Lower"].isna() | out["CI_Upper"].isna() if {"CI_Lower", "CI_Upper"}.issubset(out.columns) else pd.Series(True, index=out.index)
    out.loc[needs_ci, "CI_Lower"] = est.loc[needs_ci] - z * se.loc[needs_ci]
    out.loc[needs_ci, "CI_Upper"] = est.loc[needs_ci] + z * se.loc[needs_ci]

    out = qualify_measure_uncertainty(out, res)

    summary_rows = []
    if "SE_Status" in out.columns:
        for status, sub in out.groupby("SE_Status", dropna=False, observed=False):
            summary_rows.append({
                "Area": "Measure SE/CI",
                "Status": str(status),
                "Rows": int(len(sub)),
                "Method": "; ".join(sorted(set(sub["SE_Method"].astype(str))))[:500],
                "Caution": "; ".join(sorted(set(sub["UncertaintyCaution"].astype(str))))[:500],
            })
    if method == "MML":
        summary_rows.append({
            "Area": "MML observed-information covariance",
            "Status": str(covariance.get("status", "unknown")),
            "Rows": int(len(out)),
            "Method": "numerical Jacobian of analytical MML gradient; eigendecomposition inverse",
            "Caution": str(covariance.get("detail", "")),
            "ClaimStatus": (
                str(covariance_audit.iloc[0].get("ClaimStatus", "Do not claim"))
                if isinstance(covariance_audit, pd.DataFrame) and not covariance_audit.empty
                else "Do not claim"
            ),
            "ConditionNumber": (
                covariance_audit.iloc[0].get("ConditionNumber", np.nan)
                if isinstance(covariance_audit, pd.DataFrame) and not covariance_audit.empty
                else np.nan
            ),
            "QualificationStatus": (
                str(covariance_audit.iloc[0].get("QualificationStatus", "WITHHELD"))
                if isinstance(covariance_audit, pd.DataFrame) and not covariance_audit.empty
                else "WITHHELD"
            ),
            "ReasonCode": (
                str(covariance_audit.iloc[0].get("ReasonCode", "stat.stat_003.covariance_unavailable"))
                if isinstance(covariance_audit, pd.DataFrame) and not covariance_audit.empty
                else "stat.stat_003.covariance_unavailable"
            ),
            "PublicConclusionAllowed": (
                bool(covariance_audit.iloc[0].get("PublicConclusionAllowed", False))
                if isinstance(covariance_audit, pd.DataFrame) and not covariance_audit.empty
                else False
            ),
        })
    return out, {
        "covariance": covariance,
        "covariance_audit": covariance_audit,
        "structural_se": structural,
        "summary": pd.DataFrame(summary_rows),
    }


def calc_reliability(measure_df):
    """Compute descriptive indices only from complete paired measure inputs."""
    from mfrm_app.reliability import calculate_reliability
    return calculate_reliability(measure_df)


def compute_marginal_fit_diagnostics(
    res,
    include_pairwise=False,
    max_pair_cells=400,
    min_expected=1e-8,
):
    """MML-only marginal observed-vs-expected category diagnostics.

    This integrates category probabilities over the population quadrature
    distribution, conditional on the observed design rows. It is a marginal
    distribution screen, not a posterior/person-residual diagnostic.
    """
    config = res.get("config", {})
    if config.get("method") != "MML":
        return {
            "available": False,
            "reason": "Strict marginal diagnostics are currently defined for MML runs only.",
            "summary": pd.DataFrame(),
            "counts": pd.DataFrame(),
            "pairwise": {"available": False, "reason": "MML run required.", "summary": pd.DataFrame(), "counts": pd.DataFrame()},
        }

    prep = res.get("prep", {})
    params = res.get("params", {})
    if not prep or not params:
        return {
            "available": False,
            "reason": "Fitted data or parameters are unavailable.",
            "summary": pd.DataFrame(),
            "counts": pd.DataFrame(),
            "pairwise": {"available": False, "reason": "Fitted data or parameters are unavailable.", "summary": pd.DataFrame(), "counts": pd.DataFrame()},
        }

    idx = build_indices(prep, step_facet=config.get("step_facet"), slope_facet=config.get("slope_facet"))
    n = len(idx["score_k"])
    n_cat = int(config.get("n_cat", prep["rating_max"] - prep["rating_min"] + 1))
    if n == 0 or n_cat < 2:
        return {
            "available": False,
            "reason": "Not enough observations or score categories for marginal diagnostics.",
            "summary": pd.DataFrame(),
            "counts": pd.DataFrame(),
            "pairwise": {"available": False, "reason": "Not enough data.", "summary": pd.DataFrame(), "counts": pd.DataFrame()},
        }

    quad_points = int(config.get("quad_points") or 15)
    quad = make_mml_quadrature(config, quad_points)
    base_eta = compute_base_eta(idx, params, config)
    pop_mu = compute_population_mu(params, config)
    pop_mu_by_obs = pop_mu[idx["person"]] if pop_mu.size else np.zeros(n, dtype=float)
    k_vals = np.arange(n_cat)
    marginal_probs = np.zeros((n, n_cat), dtype=float)

    if config.get("model") == "RSM":
        step_cum = np.concatenate([[0.0], np.cumsum(params["steps"])])
        for theta_q, w_q in zip(quad["nodes"], quad["weights"]):
            eta = theta_q + pop_mu_by_obs + base_eta
            log_num = np.outer(eta, k_vals) - step_cum
            log_denom = logsumexp(log_num, axis=1, keepdims=True)
            marginal_probs += float(w_q) * np.exp(log_num - log_denom)
    elif config.get("model") == "PCM":
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])
        step_idx = idx["step_idx"]
        for theta_q, w_q in zip(quad["nodes"], quad["weights"]):
            eta = theta_q + pop_mu_by_obs + base_eta
            for c_idx in range(step_cum_mat.shape[0]):
                rows = np.where(step_idx == c_idx)[0]
                if rows.size == 0:
                    continue
                log_num = np.outer(eta[rows], k_vals) - step_cum_mat[c_idx]
                log_denom = logsumexp(log_num, axis=1, keepdims=True)
                marginal_probs[rows, :] += float(w_q) * np.exp(log_num - log_denom)
    else:
        step_cum_mat = np.vstack([
            np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]
        ])
        for theta_q, w_q in zip(quad["nodes"], quad["weights"]):
            eta = theta_q + pop_mu_by_obs + base_eta
            marginal_probs += float(w_q) * category_prob_gpcm(
                eta,
                step_cum_mat,
                idx["step_idx"],
                params["slopes"],
                idx.get("slope_idx"),
            )

    row_weight = idx.get("weight")
    if row_weight is None:
        row_weight = np.ones(n, dtype=float)
    else:
        row_weight = np.asarray(row_weight, dtype=float)

    obs_k = np.asarray(idx["score_k"], dtype=int)
    rating_min = int(prep["rating_min"])
    facet_names = list(config.get("facet_names", []))
    data = prep.get("data", pd.DataFrame())

    def aggregate(assignments: pd.DataFrame, scope: str, label_cols: list[str]):
        count_rows = []
        summary_rows = []
        if assignments.empty:
            return pd.DataFrame(), pd.DataFrame()
        grouped = assignments.groupby(label_cols, dropna=False, observed=False, sort=False) if label_cols else [((), assignments)]
        for key, sub in grouped:
            row_idx = sub["_row"].to_numpy(dtype=int)
            if row_idx.size == 0:
                continue
            key_tuple = key if isinstance(key, tuple) else (key,)
            label = {col: key_tuple[i] if i < len(key_tuple) else "" for i, col in enumerate(label_cols)}
            chi_sq = 0.0
            df_eff = 0
            total_obs = 0.0
            total_exp = 0.0
            for k in range(n_cat):
                obs = float(np.sum(row_weight[row_idx] * (obs_k[row_idx] == k)))
                exp = float(np.sum(row_weight[row_idx] * marginal_probs[row_idx, k]))
                resid = obs - exp
                std = resid / np.sqrt(exp) if exp > min_expected else np.nan
                contrib = (resid * resid / exp) if exp > min_expected else np.nan
                if np.isfinite(contrib):
                    chi_sq += contrib
                    df_eff += 1
                total_obs += obs
                total_exp += exp
                count_rows.append({
                    "Scope": scope,
                    **label,
                    "Category": rating_min + k,
                    "Observed": obs,
                    "Expected": exp,
                    "Residual": resid,
                    "StdResidual": std,
                    "ChiSquareContribution": contrib,
                })
            df = max(df_eff - 1, 0)
            p_val = float(chi2.sf(chi_sq, df)) if df > 0 and np.isfinite(chi_sq) else np.nan
            summary_rows.append({
                "Scope": scope,
                **label,
                "Rows": int(row_idx.size),
                "ObservedTotal": total_obs,
                "ExpectedTotal": total_exp,
                "ChiSquare": chi_sq,
                "df": df,
                "p": p_val,
                "MaxAbsStdResidual": float(np.nanmax(np.abs([r["StdResidual"] for r in count_rows[-n_cat:]]))),
                "Status": "Review" if (np.isfinite(p_val) and p_val < 0.01) or (np.nanmax(np.abs([r["StdResidual"] for r in count_rows[-n_cat:]])) >= 3) else "OK",
            })
        return pd.DataFrame(summary_rows), pd.DataFrame(count_rows)

    row_base = pd.DataFrame({"_row": np.arange(n, dtype=int)})
    overall_summary, overall_counts = aggregate(row_base, "Overall", [])

    summary_parts = [overall_summary]
    count_parts = [overall_counts]
    for facet in facet_names:
        assignments = pd.DataFrame({
            "_row": np.arange(n, dtype=int),
            "Facet": facet,
            "Level": data[facet].astype(str).to_numpy() if facet in data.columns else idx["facets"][facet].astype(str),
        })
        s_tbl, c_tbl = aggregate(assignments, f"Facet:{facet}", ["Facet", "Level"])
        summary_parts.append(s_tbl)
        count_parts.append(c_tbl)

    pair_summary = pd.DataFrame()
    pair_counts = pd.DataFrame()
    pair_reason = "Pairwise marginal diagnostics were not requested."
    pair_available = False
    if include_pairwise and len(facet_names) >= 2:
        pair_summary_parts = []
        pair_count_parts = []
        n_cells = 0
        for fa, fb in combinations(facet_names, 2):
            if fa not in data.columns or fb not in data.columns:
                continue
            n_cells += int(data[[fa, fb]].drop_duplicates().shape[0]) * n_cat
            if n_cells > int(max_pair_cells):
                pair_reason = (
                    f"Pairwise marginal diagnostics skipped after exceeding {max_pair_cells} expected cells. "
                    "Increase the Custom limit only when you need a final, slower screening run."
                )
                break
            assignments = pd.DataFrame({
                "_row": np.arange(n, dtype=int),
                "FacetA": fa,
                "LevelA": data[fa].astype(str).to_numpy(),
                "FacetB": fb,
                "LevelB": data[fb].astype(str).to_numpy(),
            })
            s_tbl, c_tbl = aggregate(assignments, f"Pair:{fa} x {fb}", ["FacetA", "LevelA", "FacetB", "LevelB"])
            pair_summary_parts.append(s_tbl)
            pair_count_parts.append(c_tbl)
        if pair_summary_parts:
            pair_summary = pd.concat(pair_summary_parts, ignore_index=True)
            pair_counts = pd.concat(pair_count_parts, ignore_index=True)
            pair_available = True
            pair_reason = "Pairwise marginal diagnostics computed."

    summary = pd.concat([p for p in summary_parts if isinstance(p, pd.DataFrame) and not p.empty], ignore_index=True)
    counts = pd.concat([p for p in count_parts if isinstance(p, pd.DataFrame) and not p.empty], ignore_index=True)
    n_review = int((summary.get("Status", pd.Series(dtype=str)) == "Review").sum()) if not summary.empty else 0

    return {
        "available": True,
        "reason": "Marginal diagnostics computed for an MML run.",
        "summary": summary,
        "counts": counts,
        "n_review": n_review,
        "quad_points": quad_points,
        "interpretation": (
            "Review rows with small p-values or |StdResidual| >= 3. "
            "This is a marginal distribution screen; interpret together with residual fit, category functioning, and design review."
        ),
        "pairwise": {
            "available": pair_available,
            "reason": pair_reason,
            "summary": pair_summary,
            "counts": pair_counts,
            "max_pair_cells": int(max_pair_cells),
        },
    }


def ensure_positive_definite(mat, eps=1e-6):
    mat = np.array(mat, dtype=float, copy=True)
    mat = (mat + mat.T) / 2
    eigvals = np.linalg.eigvalsh(mat)
    min_eig = np.min(eigvals)
    if min_eig < eps:
        mat += np.eye(mat.shape[0]) * (eps - min_eig)
    return mat


def _pca_core_from_residual_matrix(residual_matrix: pd.DataFrame | None) -> dict | None:
    """Compute the correlation PCA core without recursive stability checks."""
    if residual_matrix is None or residual_matrix.shape[0] < 2 or residual_matrix.shape[1] < 2:
        return None
    residual_matrix_clean = residual_matrix.dropna(axis=1, how="all")
    if residual_matrix_clean.shape[1] < 2:
        return None
    cor_df = residual_matrix_clean.corr(min_periods=2)
    cor_df = cor_df.fillna(0)
    cor_arr = cor_df.to_numpy(copy=True)
    np.fill_diagonal(cor_arr, 1)

    cor_pd = ensure_positive_definite(cor_arr)
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
        "residual_matrix_clean": residual_matrix_clean,
    }


def _abs_loading_correlation(reference: pd.Series, candidate: pd.Series) -> float:
    joined = pd.concat([reference, candidate], axis=1, join="inner").dropna()
    if joined.shape[0] < 2:
        return np.nan
    x = pd.to_numeric(joined.iloc[:, 0], errors="coerce")
    y = pd.to_numeric(joined.iloc[:, 1], errors="coerce")
    ok = x.notna() & y.notna()
    if int(ok.sum()) < 2:
        return np.nan
    x_arr = x[ok].to_numpy(dtype=float)
    y_arr = y[ok].to_numpy(dtype=float)
    if np.nanstd(x_arr) <= 0 or np.nanstd(y_arr) <= 0:
        return np.nan
    return float(abs(np.corrcoef(x_arr, y_arr)[0, 1]))


def _abs_loading_corr_arrays(x_arr: np.ndarray, y_arr: np.ndarray) -> float:
    ok = np.isfinite(x_arr) & np.isfinite(y_arr)
    n_ok = int(np.sum(ok))
    if n_ok < 2:
        return np.nan
    x = x_arr[ok]
    y = y_arr[ok]
    dx = x - np.mean(x)
    dy = y - np.mean(y)
    vx = float(np.dot(dx, dx))
    vy = float(np.dot(dy, dy))
    if vx <= 0.0 or vy <= 0.0:
        return np.nan
    return float(min(abs(float(np.dot(dx, dy)) / np.sqrt(vx * vy)), 1.0))


def _pairwise_corr_matrix_numpy(mat: np.ndarray) -> np.ndarray:
    """Pairwise Pearson correlation matching DataFrame.corr(min_periods=2).fillna(0) with unit diagonal."""
    n_rows, n_cols = mat.shape
    mask = np.isfinite(mat)
    if np.all(mask):
        if n_rows < 2:
            cor = np.zeros((n_cols, n_cols), dtype=float)
            np.fill_diagonal(cor, 1.0)
            return cor
        centered = mat - np.mean(mat, axis=0, keepdims=True)
        ss = np.sum(centered ** 2, axis=0)
        denom = np.sqrt(np.outer(ss, ss))
        cov_num = centered.T @ centered
        cor = np.divide(cov_num, denom, out=np.zeros_like(cov_num), where=denom > 0)
        cor = np.clip(cor, -1.0, 1.0)
        np.fill_diagonal(cor, 1.0)
        return cor

    mf = mask.astype(float)
    col_counts = np.sum(mf, axis=0)
    col_means = np.divide(
        np.sum(np.where(mask, mat, 0.0), axis=0),
        col_counts,
        out=np.zeros(n_cols, dtype=float),
        where=col_counts > 0,
    )
    xc = np.where(mask, mat - col_means[None, :], 0.0)
    n_ab = mf.T @ mf
    sx = xc.T @ mf
    sxx = (xc ** 2).T @ mf
    sxy = xc.T @ xc
    cov_num = sxy - np.divide(sx * sx.T, n_ab, out=np.zeros_like(sxy), where=n_ab > 0)
    var_a = np.maximum(sxx - np.divide(sx ** 2, n_ab, out=np.zeros_like(sxx), where=n_ab > 0), 0.0)
    denom = np.sqrt(var_a * var_a.T)
    cor = np.divide(cov_num, denom, out=np.zeros_like(cov_num), where=(n_ab >= 2.0) & (denom > 0.0))
    cor = np.clip(cor, -1.0, 1.0)
    np.fill_diagonal(cor, 1.0)
    return cor


def _pc1_ev1_and_loadings_from_cor(cor_arr: np.ndarray) -> tuple[float, np.ndarray]:
    cor_pd = ensure_positive_definite(cor_arr)
    eigvals, eigvecs = np.linalg.eigh(cor_pd)
    top_idx = int(np.argsort(eigvals)[::-1][0])
    ev1 = float(eigvals[top_idx])
    pc1 = eigvecs[:, top_idx] * np.sqrt(max(ev1, 0.0))
    return ev1, pc1


def _pca_leave_one_column_stability(
    residual_matrix_clean: pd.DataFrame,
    full_eigenvalues: np.ndarray,
    full_loadings: pd.DataFrame,
) -> dict[str, object]:
    if (
        residual_matrix_clean is None
        or residual_matrix_clean.shape[1] < 3
        or full_eigenvalues is None
        or len(full_eigenvalues) == 0
        or not isinstance(full_loadings, pd.DataFrame)
        or "PC1" not in full_loadings.columns
    ):
        return {
            "LOOColumnsEvaluated": 0,
            "LOOEV1Min": np.nan,
            "LOOEV1Max": np.nan,
            "LOOEV1MaxAbsDelta": np.nan,
            "LOOEV1ShareMin": np.nan,
            "LOOEV1ShareMax": np.nan,
            "LOOEV1ShareMaxAbsDelta": np.nan,
            "LOOPC1MedianAbsLoadingCorrelation": np.nan,
            "LOOPC1MinAbsLoadingCorrelation": np.nan,
        }
    clean_df = residual_matrix_clean.dropna(axis=1, how="all")
    n_total_cols = int(clean_df.shape[1])
    if n_total_cols < 3 or clean_df.shape[0] < 2:
        return {
            "LOOColumnsEvaluated": 0,
            "LOOEV1Min": np.nan,
            "LOOEV1Max": np.nan,
            "LOOEV1MaxAbsDelta": np.nan,
            "LOOEV1ShareMin": np.nan,
            "LOOEV1ShareMax": np.nan,
            "LOOEV1ShareMaxAbsDelta": np.nan,
            "LOOPC1MedianAbsLoadingCorrelation": np.nan,
            "LOOPC1MinAbsLoadingCorrelation": np.nan,
        }
    mat_full = clean_df.to_numpy(dtype=float)
    raw_cor_full = _pairwise_corr_matrix_numpy(mat_full)
    full_pc1_aligned = pd.to_numeric(full_loadings["PC1"].reindex(clean_df.columns), errors="coerce").to_numpy(dtype=float)

    n_eval = min(n_total_cols, int(PCA_STABILITY_LOO_MAX_COLUMNS))
    ev1_values: list[float] = []
    ev1_shares: list[float] = []
    loading_corrs: list[float] = []
    full_ev1 = float(full_eigenvalues[0])
    full_ev1_share = float(full_ev1 / residual_matrix_clean.shape[1]) if residual_matrix_clean.shape[1] else np.nan
    all_idx = np.arange(n_total_cols)

    for drop_j in range(n_eval):
        keep_idx = all_idx[all_idx != drop_j]
        if keep_idx.size < 2:
            continue
        sub_cor = raw_cor_full[np.ix_(keep_idx, keep_idx)]
        ev1, pc1 = _pc1_ev1_and_loadings_from_cor(sub_cor)
        ev1_values.append(ev1)
        ev1_shares.append(float(ev1 / keep_idx.size))
        corr = _abs_loading_corr_arrays(full_pc1_aligned[keep_idx], pc1)
        if np.isfinite(corr):
            loading_corrs.append(float(corr))

    ev_arr = np.asarray(ev1_values, dtype=float)
    share_arr = np.asarray(ev1_shares, dtype=float)
    corr_arr = np.asarray(loading_corrs, dtype=float)
    return {
        "LOOColumnsEvaluated": int(len(ev_arr)),
        "LOOEV1Min": float(np.nanmin(ev_arr)) if ev_arr.size else np.nan,
        "LOOEV1Max": float(np.nanmax(ev_arr)) if ev_arr.size else np.nan,
        "LOOEV1MaxAbsDelta": float(np.nanmax(np.abs(ev_arr - full_ev1))) if ev_arr.size else np.nan,
        "LOOEV1ShareMin": float(np.nanmin(share_arr)) if share_arr.size else np.nan,
        "LOOEV1ShareMax": float(np.nanmax(share_arr)) if share_arr.size else np.nan,
        "LOOEV1ShareMaxAbsDelta": float(np.nanmax(np.abs(share_arr - full_ev1_share))) if share_arr.size and np.isfinite(full_ev1_share) else np.nan,
        "LOOPC1MedianAbsLoadingCorrelation": float(np.nanmedian(corr_arr)) if corr_arr.size else np.nan,
        "LOOPC1MinAbsLoadingCorrelation": float(np.nanmin(corr_arr)) if corr_arr.size else np.nan,
    }


def _pca_row_bootstrap_stability(
    residual_matrix_clean: pd.DataFrame,
    full_eigenvalues: np.ndarray,
    full_loadings: pd.DataFrame,
) -> dict[str, object]:
    if (
        residual_matrix_clean is None
        or residual_matrix_clean.shape[0] < 4
        or residual_matrix_clean.shape[1] < 2
        or full_eigenvalues is None
        or len(full_eigenvalues) == 0
        or not isinstance(full_loadings, pd.DataFrame)
        or "PC1" not in full_loadings.columns
    ):
        return {
            "BootstrapReplicates": 0,
            "BootstrapEV1Mean": np.nan,
            "BootstrapEV1SD": np.nan,
            "BootstrapEV1CV": np.nan,
            "BootstrapEV1P05": np.nan,
            "BootstrapEV1P95": np.nan,
            "BootstrapPC1MedianAbsLoadingCorrelation": np.nan,
            "BootstrapPC1MinAbsLoadingCorrelation": np.nan,
            "BootstrapEV1ThresholdCrossing": False,
        }
    work = residual_matrix_clean
    if work.shape[1] > PCA_STABILITY_BOOTSTRAP_MAX_COLUMNS:
        work = work.iloc[:, :PCA_STABILITY_BOOTSTRAP_MAX_COLUMNS]
    work_arr = work.to_numpy(dtype=float)
    full_pc1_arr = pd.to_numeric(full_loadings["PC1"].reindex(work.columns), errors="coerce").to_numpy(dtype=float)
    rng = np.random.default_rng(PCA_STABILITY_BOOTSTRAP_SEED)
    ev1_values: list[float] = []
    loading_corrs: list[float] = []
    n_rows = int(work_arr.shape[0])
    for _ in range(PCA_STABILITY_BOOTSTRAP_REPS):
        b_idx = rng.integers(0, n_rows, size=n_rows)
        sample = work_arr[b_idx]
        nonempty_cols = np.any(np.isfinite(sample), axis=0)
        if int(np.sum(nonempty_cols)) < 2:
            continue
        sample_sub = sample[:, nonempty_cols]
        cor_boot = _pairwise_corr_matrix_numpy(sample_sub)
        ev1, pc1 = _pc1_ev1_and_loadings_from_cor(cor_boot)
        ev1_values.append(ev1)
        corr = _abs_loading_corr_arrays(full_pc1_arr[nonempty_cols], pc1)
        if np.isfinite(corr):
            loading_corrs.append(float(corr))
    ev_arr = np.asarray(ev1_values, dtype=float)
    corr_arr = np.asarray(loading_corrs, dtype=float)
    ev_mean = float(np.nanmean(ev_arr)) if ev_arr.size else np.nan
    ev_sd = float(np.nanstd(ev_arr, ddof=1)) if ev_arr.size > 1 else np.nan
    ev_cv = float(ev_sd / ev_mean) if np.isfinite(ev_sd) and np.isfinite(ev_mean) and ev_mean > 0 else np.nan
    ev_p05 = float(np.nanpercentile(ev_arr, 5)) if ev_arr.size else np.nan
    ev_p95 = float(np.nanpercentile(ev_arr, 95)) if ev_arr.size else np.nan
    full_ev1 = float(full_eigenvalues[0])
    threshold_cross = False
    for threshold in (FINAL_PCA_EIGENVALUE_READY, FINAL_PCA_EIGENVALUE_REVIEW):
        if np.isfinite(ev_p05) and np.isfinite(ev_p95) and ev_p05 < threshold <= ev_p95:
            threshold_cross = True
        if np.isfinite(full_ev1) and full_ev1 < threshold <= ev_p95:
            threshold_cross = True
        if np.isfinite(full_ev1) and full_ev1 >= threshold and ev_p05 < threshold:
            threshold_cross = True
    return {
        "BootstrapReplicates": int(len(ev_arr)),
        "BootstrapEV1Mean": ev_mean,
        "BootstrapEV1SD": ev_sd,
        "BootstrapEV1CV": ev_cv,
        "BootstrapEV1P05": ev_p05,
        "BootstrapEV1P95": ev_p95,
        "BootstrapPC1MedianAbsLoadingCorrelation": float(np.nanmedian(corr_arr)) if corr_arr.size else np.nan,
        "BootstrapPC1MinAbsLoadingCorrelation": float(np.nanmin(corr_arr)) if corr_arr.size else np.nan,
        "BootstrapEV1ThresholdCrossing": bool(threshold_cross),
    }


def compute_pca_bundle(residual_matrix_wide):
    if residual_matrix_wide is None:
        return None
    if residual_matrix_wide.shape[0] < 2 or residual_matrix_wide.shape[1] < 2:
        return None

    residual_matrix_clean = residual_matrix_wide.dropna(axis=1, how="all")
    if residual_matrix_clean.shape[1] < 2:
        return None

    valid_mask = residual_matrix_clean.notna().astype(int)
    pairwise_n = valid_mask.T.dot(valid_mask)
    if pairwise_n.shape[0] > 1:
        off_diag = pairwise_n.to_numpy(dtype=float)
        off_diag = off_diag[~np.eye(off_diag.shape[0], dtype=bool)]
        min_pairwise_n = float(np.nanmin(off_diag)) if off_diag.size else np.nan
        median_pairwise_n = float(np.nanmedian(off_diag)) if off_diag.size else np.nan
    else:
        min_pairwise_n = np.nan
        median_pairwise_n = np.nan
    missing_share = float(residual_matrix_clean.isna().to_numpy().mean())
    n_persons, n_columns = residual_matrix_clean.shape
    stability_flags = []
    if n_persons < PCA_STABILITY_MIN_PERSONS:
        stability_flags.append(f"fewer than {PCA_STABILITY_MIN_PERSONS} persons")
    if n_columns < PCA_STABILITY_MIN_COLUMNS:
        stability_flags.append(f"fewer than {PCA_STABILITY_MIN_COLUMNS} residual columns")
    if missing_share > PCA_STABILITY_MAX_MISSING_SHARE:
        stability_flags.append(f"{missing_share:.0%} matrix missingness")
    if np.isfinite(min_pairwise_n) and min_pairwise_n < PCA_STABILITY_MIN_PAIRWISE_N:
        stability_flags.append(f"minimum pairwise overlap {min_pairwise_n:.0f}")

    core = _pca_core_from_residual_matrix(residual_matrix_clean)
    if core is None:
        return None
    eigvals = np.asarray(core["eigenvalues"], dtype=float)
    var_pct = np.asarray(core["variance_pct"], dtype=float)
    loadings_df = core["loadings"]
    cor_pd_df = core["cor_matrix"]

    loo = _pca_leave_one_column_stability(residual_matrix_clean, eigvals, loadings_df)
    boot = _pca_row_bootstrap_stability(residual_matrix_clean, eigvals, loadings_df)
    sensitivity_flags: list[str] = []
    loo_delta = float(loo.get("LOOEV1MaxAbsDelta", np.nan))
    loo_share_delta = float(loo.get("LOOEV1ShareMaxAbsDelta", np.nan))
    loo_corr = float(loo.get("LOOPC1MinAbsLoadingCorrelation", np.nan))
    boot_cv = float(boot.get("BootstrapEV1CV", np.nan))
    boot_corr = float(boot.get("BootstrapPC1MinAbsLoadingCorrelation", np.nan))
    if np.isfinite(loo_share_delta) and loo_share_delta > PCA_STABILITY_MAX_LOO_EV1_SHARE_DELTA:
        sensitivity_flags.append(f"leave-one-column EV1 share shift {loo_share_delta:.2f}")
    if np.isfinite(loo_corr) and loo_corr < PCA_STABILITY_MIN_LOADING_CORR:
        sensitivity_flags.append(f"leave-one-column loading correlation {loo_corr:.2f}")
    if np.isfinite(boot_cv) and boot_cv > PCA_STABILITY_MAX_BOOTSTRAP_EV1_CV:
        sensitivity_flags.append(f"bootstrap EV1 CV {boot_cv:.2f}")
    if np.isfinite(boot_corr) and boot_corr < PCA_STABILITY_MIN_LOADING_CORR:
        sensitivity_flags.append(f"bootstrap loading correlation {boot_corr:.2f}")
    if bool(boot.get("BootstrapEV1ThresholdCrossing", False)):
        sensitivity_flags.append("bootstrap EV1 interval crosses an interpretation threshold")
    stability_flags.extend(sensitivity_flags)
    stability_status = "Review" if stability_flags else "Stable screen"
    stability_caution = (
        "Residual PCA may be unstable because " + "; ".join(stability_flags) + ". "
        "Use loadings as exploratory prompts and corroborate with fit, content, and marginal diagnostics."
        if stability_flags else
        "Residual matrix has adequate basic overlap for a screening PCA; still interpret as a diagnostic, not proof."
    )
    stability_row = {
        "PCAStabilityStatus": stability_status,
        "Persons": int(n_persons),
        "ResidualColumns": int(n_columns),
        "MissingShare": missing_share,
        "MinPairwiseOverlap": min_pairwise_n,
        "MedianPairwiseOverlap": median_pairwise_n,
        "FullEV1": float(eigvals[0]) if eigvals.size else np.nan,
        "FullEV2": float(eigvals[1]) if eigvals.size > 1 else np.nan,
        "FullPC1VariancePct": float(var_pct[0]) if var_pct.size else np.nan,
        "SensitivityFlagSummary": "; ".join(sensitivity_flags) if sensitivity_flags else "none",
        "Rule": (
            f"Review if persons < {PCA_STABILITY_MIN_PERSONS}, columns < "
            f"{PCA_STABILITY_MIN_COLUMNS}, missing share > "
            f"{PCA_STABILITY_MAX_MISSING_SHARE:.0%}, or min pairwise overlap < "
            f"{PCA_STABILITY_MIN_PAIRWISE_N}; also review if leave-one-column "
            f"EV1 share shift > {PCA_STABILITY_MAX_LOO_EV1_SHARE_DELTA:.0%}, bootstrap EV1 CV > "
            f"{PCA_STABILITY_MAX_BOOTSTRAP_EV1_CV:.0%}, loading correlation < "
            f"{PCA_STABILITY_MIN_LOADING_CORR}, or bootstrap EV1 crosses EV=2/3 thresholds."
        ),
        "Caution": stability_caution,
    }
    stability_row.update(loo)
    stability_row.update(boot)
    stability_table = pd.DataFrame([stability_row])

    return {
        "eigenvalues": eigvals,
        "variance_pct": var_pct,
        "loadings": loadings_df,
        "cor_matrix": cor_pd_df,
        "residual_matrix": residual_matrix_wide,
        "stability": stability_table.iloc[0].to_dict(),
        "stability_table": stability_table,
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


def diagnose_pca_skip_reason(obs_df, facet_names, *, mode="overall", facet=None):
    """Return a human-readable reason why residual PCA could not be computed.

    Called by ``mfrm_diagnostics`` whenever ``compute_pca_overall`` /
    ``compute_pca_by_facet`` returns None. The Dimensionality tab now
    surfaces this reason instead of the silent "PCA of residuals is
    not available for this estimation" message.

    Parameters
    ----------
    obs_df : pd.DataFrame | None
        Observation table from ``compute_obs_table(res)``.
    facet_names : list[str] | None
        Facet names from ``res["config"]["facet_names"]``.
    mode : {"overall", "facet"}
        "overall" concatenates all facets into ``item_combination``;
        "facet" keys the residual matrix by a single facet level.
    facet : str | None
        Required when ``mode="facet"``.

    Returns
    -------
    dict
        A reason payload of the form ``{"key": str, "kwargs": dict}``
        suitable for ``_render_reason``. The payload is locale-agnostic
        and re-translates on every Streamlit rerun, so the Dimensionality
        tab matches the active sidebar language even when the user
        toggles language *after* the run. Callers must treat the value
        as opaque and never match it against fixed English literals.
    """
    if obs_df is None:
        return _make_reason("dimensionality.skip_reason_obs_none")
    if not isinstance(obs_df, pd.DataFrame):
        return _make_reason(
            "dimensionality.skip_reason_obs_type_template",
            type_name=type(obs_df).__name__,
        )
    if obs_df.empty:
        return _make_reason("dimensionality.skip_reason_obs_empty")
    if mode == "overall" and not facet_names:
        return _make_reason("dimensionality.skip_reason_no_facets")
    if mode == "facet":
        if not facet:
            return _make_reason("dimensionality.skip_reason_facet_not_provided")
        if facet not in obs_df.columns:
            return _make_reason(
                "dimensionality.skip_reason_facet_not_in_obs_template",
                facet=facet,
            )
    if "Person" not in obs_df.columns:
        return _make_reason("dimensionality.skip_reason_person_missing")
    if "StdResidual" not in obs_df.columns:
        return _make_reason("dimensionality.skip_reason_stdresidual_missing")

    try:
        df_aug = obs_df.copy()
        df_aug["Person"] = df_aug["Person"].astype(str)
        if mode == "overall":
            df_aug["_col_key"] = df_aug[facet_names].astype(str).agg("_".join, axis=1)
        else:
            df_aug["_col_key"] = df_aug[facet].astype(str)
        residual_matrix_prep = (
            df_aug.groupby(["Person", "_col_key"])["StdResidual"].mean().reset_index()
        )
        residual_matrix_wide = residual_matrix_prep.pivot(
            index="Person", columns="_col_key", values="StdResidual"
        )
    except Exception as exc:  # pragma: no cover - defensive
        return _make_reason(
            "dimensionality.skip_reason_construction_error_template",
            error_type=type(exc).__name__,
            error=str(exc),
        )

    n_persons, n_cols = residual_matrix_wide.shape
    if n_persons < 2:
        return _make_reason(
            "dimensionality.skip_reason_too_few_persons_template",
            n_persons=n_persons,
        )
    if n_cols < 2:
        # ``key_label`` itself is a translatable fragment, so we record the
        # mode discriminator and the facet kwarg and let _render_reason
        # resolve the final concatenation at display time.
        if mode == "overall":
            return _make_reason(
                "dimensionality.skip_reason_too_few_cols_template",
                n_cols=n_cols,
                key_label=_LazyKeyLabel(mode="overall"),
            )
        return _make_reason(
            "dimensionality.skip_reason_too_few_cols_template",
            n_cols=n_cols,
            key_label=_LazyKeyLabel(mode="facet", facet=facet),
        )
    clean = residual_matrix_wide.dropna(axis=1, how="all")
    if clean.shape[1] < 2:
        return _make_reason(
            "dimensionality.skip_reason_after_dropna_template",
            n_remaining=clean.shape[1],
        )

    nan_share = float(residual_matrix_wide.isna().to_numpy().mean())
    if nan_share > 0.95:
        return _make_reason(
            "dimensionality.skip_reason_too_sparse_template",
            nan_pct=f"{nan_share:.0%}",
        )

    return _make_reason("dimensionality.skip_reason_unknown")


def mfrm_diagnostics(
    res,
    interaction_pairs=None,
    top_n_interactions=20,
    whexact=False,
    compute_pca=True,
    compute_marginal=False,
    marginal_pairwise=False,
    marginal_max_pair_cells=400,
    compute_eb_shrinkage=False,
    fit_df_method: str = DEFAULT_FIT_DF_METHOD,
    facets_zstd_cap: float = DEFAULT_FACETS_ZSTD_CAP,
    *,
    zstd_transform=None,
):
    zstd_transform = _resolve_zstd_transform(whexact, zstd_transform)
    fit_df_method = _resolve_fit_df_method(fit_df_method)
    model_family = str(res.get("config", {}).get("model", "RSM")).upper()
    fit_df_metadata = fit_df_model_metadata(model_family, fit_df_method)
    obs_df = compute_obs_table(res)
    facet_cols = ["Person"] + res["config"]["facet_names"]
    overall_fit = calc_overall_fit(
        obs_df,
        whexact=whexact,
        zstd_transform=zstd_transform,
        fit_df_method=fit_df_method,
        facets_zstd_cap=facets_zstd_cap,
    )
    fit_tbl = calc_facet_fit(
        obs_df,
        facet_cols,
        whexact=whexact,
        zstd_transform=zstd_transform,
        fit_df_method=fit_df_method,
        facets_zstd_cap=facets_zstd_cap,
    )
    if fit_df_method != "engine":
        overall_fit["FitDfApplicability"] = fit_df_metadata["applicability"]
        fit_tbl["FitDfApplicability"] = fit_df_metadata["applicability"]
        if fit_df_metadata["warning"]:
            overall_fit["FitDfWarning"] = fit_df_metadata["warning"]
            fit_tbl["FitDfWarning"] = fit_df_metadata["warning"]
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
    person_boundary_columns = [
        column for column in (
            "ExtremeScorePattern",
            "ExtremeScoreDirection",
            "FiniteJMLEEstimate",
            "PersonInferenceReady",
            "ReportableEstimate",
            "EstimateRole",
        )
        if column in person_tbl.columns
    ]
    if person_boundary_columns:
        measures = measures.merge(
            person_tbl[["Facet", "Level", *person_boundary_columns]],
            on=["Facet", "Level"],
            how="left",
            validate="one_to_one",
        )

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
    ci_z = float(_norm.ppf(0.975))
    measures["CI_Lower"] = measures["Estimate"] - ci_z * measures["SE"]
    measures["CI_Upper"] = measures["Estimate"] + ci_z * measures["SE"]
    measures, uncertainty_bundle = annotate_measure_uncertainty(
        measures,
        res,
        obs_df=obs_df,
        ci_level=0.95,
    )

    # Attach MML observed-information delta-method SEs / CIs to the fitted step
    # thresholds, reusing the covariance already computed for the facet
    # measures. MML only; a no-op for JMLE fits and when the steps table
    # already carries an SE column.
    try:
        steps_tbl = res.get("steps") if isinstance(res, dict) else None
        if (
            isinstance(steps_tbl, pd.DataFrame)
            and not steps_tbl.empty
            and "Estimate" in steps_tbl.columns
            and "SE" not in steps_tbl.columns
        ):
            step_se_tbl = compute_mml_step_se_table(
                res,
                covariance=uncertainty_bundle.get("covariance") if isinstance(uncertainty_bundle, dict) else None,
                ci_level=0.95,
            )
            if isinstance(step_se_tbl, pd.DataFrame) and not step_se_tbl.empty:
                key_cols = [c for c in ("StepFacet", "Step") if c in steps_tbl.columns and c in step_se_tbl.columns]
                if key_cols:
                    res["steps"] = steps_tbl.merge(step_se_tbl, on=key_cols, how="left")
    except Exception:
        pass

    reliability_tbl = calc_reliability(measures)
    krippendorff_alpha_tbl = calc_krippendorff_alpha_by_facet(
        obs_df,
        ["Person"] + list(res["config"]["facet_names"]),
    )

    # PCA of standardised residuals
    facet_names = res["config"]["facet_names"]
    pca_overall = None
    pca_by_facet = {}
    pca_reason: dict | None = None
    # Reason payloads are stored as dicts ({"key": ..., "kwargs": ...}) so
    # the Dimensionality tab can re-translate them on every rerun even
    # when the user toggles the sidebar language after the run completes.
    pca_by_facet_reasons: dict[str, dict] = {}
    if compute_pca:
        try:
            pca_overall = compute_pca_overall(obs_df, facet_names)
            if pca_overall is None:
                pca_reason = diagnose_pca_skip_reason(obs_df, facet_names, mode="overall")
        except Exception as exc:
            pca_overall = None
            pca_reason = _make_reason(
                "dimensionality.exception_overall_template",
                error_type=type(exc).__name__,
                error=str(exc),
            )
        try:
            pca_by_facet = compute_pca_by_facet(obs_df, facet_names)
            # Record a per-facet skip reason so each Dimensionality sub-tab
            # can show *why* its panel is empty even when overall PCA succeeded.
            for facet in facet_names or []:
                bundle = pca_by_facet.get(facet) if isinstance(pca_by_facet, dict) else None
                if bundle is None:
                    pca_by_facet_reasons[facet] = diagnose_pca_skip_reason(
                        obs_df, facet_names, mode="facet", facet=facet,
                    )
        except Exception as exc:
            pca_by_facet = {}
            pca_by_facet_reasons = {
                facet: _make_reason(
                    "dimensionality.exception_facet_template",
                    error_type=type(exc).__name__,
                    error=str(exc),
                )
                for facet in (facet_names or [])
            }
    else:
        pca_reason = _make_reason("dimensionality.disabled_by_depth_setting")

    marginal_fit = {
        "available": False,
        "reason": "Strict marginal diagnostics were skipped for this run.",
        "summary": pd.DataFrame(),
        "counts": pd.DataFrame(),
        "pairwise": {"available": False, "reason": "Strict marginal diagnostics were skipped.", "summary": pd.DataFrame(), "counts": pd.DataFrame()},
    }
    if compute_marginal:
        try:
            marginal_fit = compute_marginal_fit_diagnostics(
                res,
                include_pairwise=bool(marginal_pairwise),
                max_pair_cells=int(marginal_max_pair_cells),
            )
        except Exception as marg_exc:
            marginal_fit = {
                "available": False,
                "reason": f"Strict marginal diagnostics failed: {marg_exc}",
                "summary": pd.DataFrame(),
                "counts": pd.DataFrame(),
                "pairwise": {"available": False, "reason": f"Strict marginal diagnostics failed: {marg_exc}", "summary": pd.DataFrame(), "counts": pd.DataFrame()},
            }

    if compute_eb_shrinkage:
        shrinkage = compute_empirical_bayes_shrinkage(measures)
    else:
        shrinkage = {
            "available": False,
            "reason": "Empirical-Bayes shrinkage was not requested for this run.",
            "measures": pd.DataFrame(),
            "report": pd.DataFrame(),
        }

    return {
        "obs": obs_df,
        "overall_fit": overall_fit,
        "measures": measures,
        "fit": fit_tbl,
        "reliability": reliability_tbl,
        "krippendorff_alpha": krippendorff_alpha_tbl,
        "bias": bias_tbl,
        "interactions": interaction_tbl,
        "pca": pca_overall,
        "pca_by_facet": pca_by_facet,
        "pca_reason": pca_reason,
        "pca_by_facet_reasons": pca_by_facet_reasons,
        "pca_enabled": bool(compute_pca),
        "marginal_fit": marginal_fit,
        "marginal_fit_enabled": bool(compute_marginal),
        "eb_shrinkage": shrinkage,
        "eb_shrinkage_enabled": bool(compute_eb_shrinkage),
        "uncertainty": uncertainty_bundle,
        "fit_df_method": fit_df_method,
        "fit_df_applicability": fit_df_metadata["applicability"],
        "fit_df_warning": fit_df_metadata["warning"],
        "fit_df_model": model_family,
        "fit_zstd_transform_id": zstd_transform,
        "fit_zstd_transform": fit_zstd_transform_label(whexact, zstd_transform=zstd_transform),
        "facets_zstd_cap": (
            float(facets_zstd_cap)
            if (facets_zstd_cap is not None and np.isfinite(facets_zstd_cap) and facets_zstd_cap > 0)
            else None
        ),
    }


def _result_identity_settings(result: dict) -> dict[str, object]:
    config = result.get("config", {}) if isinstance(result, dict) else {}
    config = config if isinstance(config, dict) else {}
    config_fp = config_fingerprint(config, length=64)
    opt = result.get("opt")
    fitted_state = {
        "coordinates": getattr(opt, "x", None),
        "success": getattr(opt, "success", None),
        "params": result.get("params"),
    }
    summary = result.get("summary")
    if isinstance(summary, pd.DataFrame):
        fitted_state["summary"] = summary.drop(columns=["ElapsedSeconds"], errors="ignore")
    # Bind the estimates consumed by reports, excluding subsequently added SE/CI.
    facets = result.get("facets", {})
    population = result.get("population", {})
    for name, frame in (
        ("person", facets.get("person")), ("others", facets.get("others")),
        ("steps", result.get("steps")), ("slopes", result.get("slopes")),
        ("population", population.get("coefficients")),
    ):
        if isinstance(frame, pd.DataFrame):
            columns = [c for c in ("Person", "Facet", "Level", "StepFacet", "Step",
                                   "SlopeFacet", "Term", "Estimate", "LogEstimate", "SD")
                       if c in frame.columns]
            fitted_state[name] = frame[columns]

    def finite_number(value):
        try:
            number = float(value)
        except (TypeError, ValueError):
            return None
        return number if np.isfinite(number) else None

    return {
        "schema_version": "mfrm_result_identity_settings_v2",
        "base_config_fingerprint": str(config_fp),
        "normalized_response_fingerprint": _result_normalized_response_fingerprint(result),
        # Planned sensitivity variants bind this supplied baseline fit; their
        # future refit is checked separately by _mml_prior_refit_attestation.
        "source_fitted_state_fingerprint": stable_json_fingerprint(fitted_state, length=64),
        "model": str(config.get("model", "unknown")),
        "method": str(config.get("method", "unknown")),
        "facet_names": [str(value) for value in config.get("facet_names", [])],
        "rating_min": finite_number(config.get("rating_min")),
        "rating_max": finite_number(config.get("rating_max")),
        "step_facet": None if config.get("step_facet") is None else str(config.get("step_facet")),
        "slope_facet": None if config.get("slope_facet") is None else str(config.get("slope_facet")),
        "mml_engine": None if config.get("mml_engine") is None else str(config.get("mml_engine")),
        "quad_points": finite_number(config.get("quad_points")),
        "population_prior_sd": finite_number(config.get("population_prior_sd")),
        "estimate_population_sd": bool(config.get("estimate_population_sd", False)),
        "facet_regularization_fingerprint": (
            None
            if config.get("facet_regularization_fingerprint") is None
            else str(config.get("facet_regularization_fingerprint"))
        ),
    }


def _result_input_data_fingerprint(result: dict) -> str:
    config = result.get("config", {}) if isinstance(result, dict) else {}
    config = config if isinstance(config, dict) else {}
    recorded = str(config.get("input_data_fingerprint") or "").strip()
    if recorded:
        return recorded
    prep = result.get("prep", {}) if isinstance(result, dict) else {}
    prep = prep if isinstance(prep, dict) else {}
    data = prep.get("data")
    if isinstance(data, pd.DataFrame):
        return f"dataframe_sha256:{dataframe_fingerprint(data, length=64)}"
    metadata = {
        "n_obs": prep.get("n_obs"),
        "n_person": prep.get("n_person"),
        "facet_names": prep.get("facet_names", config.get("facet_names", [])),
    }
    return f"metadata_only_sha256:{stable_json_fingerprint(metadata, length=64)}"


def _result_normalized_response_fingerprint(result: dict) -> str | None:
    """Fingerprint the normalized likelihood input, independent of upload columns."""
    prep = result.get("prep", {}) if isinstance(result, dict) else {}
    prep = prep if isinstance(prep, dict) else {}
    data = prep.get("data")
    if not isinstance(data, pd.DataFrame) or data.empty:
        return None
    fingerprint = dataframe_fingerprint(data, length=64)
    return f"normalized_dataframe_sha256:{fingerprint}" if fingerprint else None


def build_result_analysis_identity(
    result: dict,
    *,
    analysis_type: str = "mfrm.fit",
    settings_overrides: dict[str, object] | None = None,
    variant_id: str = "baseline",
    baseline: _evidence.AnalysisIdentity | None = None,
) -> _evidence.AnalysisIdentity:
    """Build the deterministic identity used by readiness and sensitivity ledgers."""
    settings = _result_identity_settings(result)
    if settings_overrides:
        settings.update(settings_overrides)
    config = result.get("config", {}) if isinstance(result, dict) else {}
    config = config if isinstance(config, dict) else {}
    return _evidence.build_analysis_identity(
        input_data_fingerprint=_result_input_data_fingerprint(result),
        resolved_settings=settings,
        analysis_type=analysis_type,
        engine_version=str(config.get("app_version") or APP_VERSION),
        variant_id=variant_id,
        baseline=baseline,
        seed=None,
    )
