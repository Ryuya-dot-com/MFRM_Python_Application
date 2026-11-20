"""Fail-closed qualification contracts for statistically sensitive outputs.

The numerical engine can retain intermediate values for debugging and migration
while the public application withholds conclusions that have not cleared their
validation gate.  This module keeps that distinction machine-readable and free
of Streamlit dependencies.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from typing import Iterable

import numpy as np
import pandas as pd


OUTPUT_QUALIFICATION_SCHEMA_VERSION = "mfrm_output_qualification_v1"
FREE_SD_MML_INFERENCE_HOLD_REASON = "stat.mml_free_sd.numerical_qualification_pending"
FIXED_SD_MML_INFERENCE_HOLD_REASON = "stat.mml_fixed_sd.numerical_qualification_pending"
OPTIMIZER_TERMINATION_HOLD_REASON = "stat.estimation.optimizer_termination_unverified"
JMLE_RANK_DEFICIENT_HOLD_REASON = "stat.jmle.eta_rank_deficient"
JMLE_AUDIT_SIZE_LIMIT_HOLD_REASON = "stat.jmle.identifiability_audit_size_limit"
CATEGORY_UNSUPPORTED_HOLD_REASON = "stat.category.unsupported_internal_step_contrast"
UNOBSERVED_BOUNDARY_CATEGORY_HOLD_REASON = "stat.category.unobserved_boundary_category"
UNOBSERVED_OWNER_ENDPOINT_HOLD_REASON = "stat.category.unobserved_owner_endpoint"
MML_FACET_BOUNDARY_HOLD_REASON = "stat.mml.extreme_free_facet_location"
MML_DISCONNECTED_HOLD_REASON = "stat.mml.population_assumption_linked"


def _empty_boundary_audit() -> dict[str, object]:
    return {
        "active": False,
        "unsupported_internal": False,
        "unobserved_global_boundary": False,
        "unobserved_owner_endpoints": set(),
        "unbounded_step_owners": set(),
        "unbounded_facet_levels": set(),
        "contaminated_facets": set(),
        "disconnected_facets": set(),
        "hypergraph_components": 1,
        "observed_categories": {},
    }


def _connected_component_count(adjacency: dict) -> int:
    seen = set()
    comps = 0
    for node in adjacency:
        if node in seen:
            continue
        comps += 1
        stack = [node]
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            stack.extend(adjacency.get(cur, ()) - seen)
    return comps


def _unlinked_anchor_levels(data, config):
    """Test level estimability after removing an unrestricted intercept for each person.

    A fixed population distribution can link otherwise disconnected cohorts.
    Within-person differences remove that assumption. A target is design-linked
    exactly when its contrast is in the row space of those differences.
    """
    from .engine import _constraint_expansion_matrix_for_identifiability

    facets = config["facet_names"]
    specs = config["facet_specs"]
    dimension = sum(specs[f]["n_params"] for f in facets)
    observations = data[["Person", *facets]].drop_duplicates()
    n_levels = sum(len(specs[f]["levels"]) for f in facets)
    if max(len(observations), n_levels, dimension) * dimension > 5_000_000:
        return None
    design = np.zeros((len(observations), dimension))
    targets = np.zeros((n_levels, dimension))
    labels = []
    column = row = 0
    for facet in facets:
        spec = specs[facet]
        block = _constraint_expansion_matrix_for_identifiability(spec)
        levels = pd.Index(spec["levels"]).astype(str)
        indices = levels.get_indexer(observations[facet].astype(str))
        if (indices < 0).any():
            return None
        cols = slice(column, column + block.shape[1])
        design[:, cols] = config.get("facet_signs", {}).get(facet, -1) * block[indices]
        targets[row:row + len(levels), cols] = block
        labels.extend((facet, level) for level in levels)
        row += len(levels)
        column = cols.stop
    if not dimension:
        return set()
    people, names = pd.factorize(observations.Person.astype(str))
    _, first = np.unique(people, return_index=True)
    differences = design - design[first][people]
    _, singular, vh = np.linalg.svd(differences, full_matrices=False)
    tol = max(differences.shape) * np.finfo(float).eps
    basis = vh[singular > tol * singular[0]]

    def unlinked(matrix):
        residual = matrix - (matrix @ basis.T) @ basis
        return np.linalg.norm(residual, axis=1) > tol * np.maximum(1., np.linalg.norm(matrix, axis=1))

    unresolved = {label for label, held in zip(labels, unlinked(targets)) if held}
    unresolved.update(("Person", name) for name, held in zip(names, unlinked(design[first])) if held)
    return unresolved


def audit_boundary_and_support(result) -> dict[str, object]:
    """Audit category support, free endpoint facets, and cohort connectivity from retained data."""
    if not isinstance(result, dict) or not isinstance(result.get("params"), dict):
        return _empty_boundary_audit()
    prep = result.get("prep") if isinstance(result.get("prep"), dict) else {}
    config = result.get("config") if isinstance(result.get("config"), dict) else {}
    data = prep.get("data")
    if not isinstance(data, pd.DataFrame) or data.empty or "Score" not in data.columns:
        return _empty_boundary_audit()

    if "Weight" in data.columns:
        w = pd.to_numeric(data["Weight"], errors="coerce")
        df = data.loc[np.isfinite(w) & (w > 0)]
    else:
        df = data
    if df.empty:
        return _empty_boundary_audit()

    bounds = pd.to_numeric(pd.Series([prep.get("rating_min"), prep.get("rating_max")]), errors="coerce").to_numpy(dtype=float)
    if not (
        np.isfinite(bounds).all()
        and np.equal(bounds, np.floor(bounds)).all()
        and 0 < bounds[1] - bounds[0] < 1000
    ):
        return _empty_boundary_audit()
    rmin, rmax = int(bounds[0]), int(bounds[1])
    scores_num = pd.to_numeric(df["Score"], errors="coerce")
    valid_s = np.isfinite(scores_num)
    df = df.loc[valid_s]
    if df.empty:
        return _empty_boundary_audit()
    scores = scores_num.loc[valid_s].astype(int)

    facet_names = [str(f) for f in (config.get("facet_names") or prep.get("facet_names") or ())]
    levels_cfg = prep.get("levels") if isinstance(prep.get("levels"), dict) else {}
    observed_categories: dict[tuple[str, str], int] = {}
    for f_col in ["Person", *facet_names]:
        for lv in (levels_cfg.get(f_col) or ()):
            observed_categories[(f_col, str(lv))] = 0
        if f_col in df.columns:
            nunq = scores.groupby(df[f_col].astype(str), observed=True).nunique()
            for lv, cnt in nunq.items():
                observed_categories[(f_col, str(lv))] = int(cnt)

    model = str(config.get("model", "RSM")).upper()
    method = str(config.get("method", "JMLE")).upper()
    step_facet = str(config.get("step_facet") or "") if model in {"PCM", "GPCM"} else ""

    global_vc = scores.value_counts()
    unsupported_internal = any(int(global_vc.get(c, 0)) == 0 for c in range(rmin + 1, rmax))
    unobserved_global_boundary = any(int(global_vc.get(c, 0)) == 0 for c in (rmin, rmax))
    unobserved_owner_endpoints: set[tuple[str, str]] = set()
    unbounded_step_owners: set[tuple[str, str]] = set()

    if model in {"PCM", "GPCM"} and step_facet and step_facet in df.columns:
        unsupported_internal = False
        step_levels = [str(x) for x in (levels_cfg.get(step_facet) or df[step_facet].astype(str).unique())]
        sf_col = df[step_facet].astype(str)
        for lv in step_levels:
            vc = scores.loc[sf_col.eq(lv)].value_counts()
            n_dist = int((vc > 0).sum())
            has_internal_gap = any(int(vc.get(c, 0)) == 0 for c in range(rmin + 1, rmax))
            has_endpoint_gap = any(int(vc.get(c, 0)) == 0 for c in (rmin, rmax))
            if n_dist < 2 or has_internal_gap:
                unsupported_internal = True
                unbounded_step_owners.add((step_facet, lv))
            elif has_endpoint_gap:
                unobserved_owner_endpoints.add((step_facet, lv))
                unbounded_step_owners.add((step_facet, lv))

    facet_specs = config.get("facet_specs") if isinstance(config.get("facet_specs"), dict) else {}
    theta_spec = config.get("theta_spec") if isinstance(config.get("theta_spec"), dict) else {}
    reg_enabled = bool((config.get("facet_regularization") or {}).get("enabled", False))

    person_free_shift = (
        method == "JMLE"
        and bool(theta_spec)
        and not bool(theta_spec.get("centered", False))
        and not _spec_has_explicit_anchor(theta_spec)
    )
    pop_free_shift = method == "MML" and bool((config.get("population_model") or {}).get("enabled", False))
    uncentered_free_facets = {
        fname
        for fname in facet_names
        if isinstance(facet_specs.get(fname), dict)
        and not bool(facet_specs[fname].get("centered", False))
        and not _spec_has_explicit_anchor(facet_specs[fname])
    }
    common_shift = bool(person_free_shift or pop_free_shift or uncentered_free_facets)

    unbounded_facet_levels: set[tuple[str, str]] = set()
    contaminated_facets: set[str] = set()
    if not reg_enabled:
        for fname in facet_names:
            spec = facet_specs.get(fname)
            if not isinstance(spec, dict) or fname not in df.columns:
                continue
            levels = [str(x) for x in (spec.get("levels") or ())]
            raw_anchors = spec.get("anchors")
            anchors_arr = (
                np.asarray(raw_anchors, dtype=float)
                if raw_anchors is not None and len(raw_anchors) == len(levels)
                else np.full(len(levels), np.nan)
            )
            free = ~np.isfinite(anchors_arr)
            raw_groups = spec.get("groups") or [None] * len(levels)
            groups = [
                None if (g is None or g == "" or (isinstance(g, float) and np.isnan(g)) or pd.isna(g)) else str(g)
                for g in raw_groups
            ]
            group_counts: dict[str | None, int] = {}
            for is_free, grp in zip(free, groups):
                if is_free:
                    group_counts[grp] = group_counts.get(grp, 0) + 1
            centered = bool(spec.get("centered", False))
            agg = scores.groupby(df[fname].astype(str), observed=True).agg(["size", "min", "max"])

            bucket_min: dict[str, list[str]] = {}
            bucket_max: dict[str, list[str]] = {}
            for lv, is_free, grp in zip(levels, free, groups):
                varies = bool(is_free) and (
                    group_counts.get(grp, 0) > 1 if (grp is not None or centered) else True
                )
                if not varies or lv not in agg.index:
                    continue
                row_s = agg.loc[lv]
                if int(row_s["size"]) <= 0 or row_s["min"] != row_s["max"]:
                    continue
                bucket = f"grp:{grp}" if grp is not None else ("centered" if centered else "uncentered")
                if int(row_s["min"]) == rmin:
                    bucket_min.setdefault(bucket, []).append(lv)
                elif int(row_s["max"]) == rmax:
                    bucket_max.setdefault(bucket, []).append(lv)

            for bucket in set(bucket_min) | set(bucket_max):
                mins = bucket_min.get(bucket, [])
                maxs = bucket_max.get(bucket, [])
                is_uncentered = bucket == "uncentered"
                opposite_pair = bool(mins and maxs)
                if method == "JMLE" or is_uncentered or common_shift or opposite_pair:
                    for lv in [*mins, *maxs]:
                        unbounded_facet_levels.add((fname, lv))
                    if not is_uncentered:
                        contaminated_facets.add(fname)
                    if method == "MML" and common_shift and not is_uncentered:
                        contaminated_facets.update(uncentered_free_facets)

    disconnected_facets: set[str] = set()
    hypergraph_components = 1
    if "Person" in df.columns:
        persons = df["Person"].astype(str)
        all_edges = {}
        for fname in facet_names:
            if fname not in df.columns:
                continue
            edges = {}
            for person, level in zip(persons, df[fname].astype(str)):
                pn, fn = ("Person", person), (fname, level)
                edges.setdefault(pn, set()).add(fn)
                edges.setdefault(fn, set()).add(pn)
            if _connected_component_count(edges) > 1:
                disconnected_facets.add(fname)
            for node, neighbors in edges.items():
                all_edges.setdefault(node, set()).update(neighbors)
        hypergraph_components = max(1, _connected_component_count(all_edges))

    disconnected_levels = None
    if (method == "MML" and (disconnected_facets or hypergraph_components > 1)
            and any(_spec_has_explicit_anchor(spec) for spec in facet_specs.values())
            and not _spec_has_explicit_anchor(theta_spec)):
        disconnected_levels = _unlinked_anchor_levels(df, config)
        if disconnected_levels is not None:
            disconnected_facets = {f for f, _ in disconnected_levels if f != "Person"}
            if not disconnected_levels:
                hypergraph_components = 1

    return {
        "active": True,
        "unsupported_internal": unsupported_internal,
        "unobserved_global_boundary": unobserved_global_boundary,
        "unobserved_owner_endpoints": unobserved_owner_endpoints,
        "unbounded_step_owners": unbounded_step_owners,
        "unbounded_facet_levels": unbounded_facet_levels,
        "contaminated_facets": contaminated_facets,
        "disconnected_facets": disconnected_facets,
        "disconnected_levels": disconnected_levels,
        "hypergraph_components": hypergraph_components,
        "observed_categories": observed_categories,
    }


def boundary_and_support_hold_reason(result, *, for_screening_block: bool = False) -> str:
    """Return a data-support, endpoint-boundary, or cohort-connectivity hold code when active."""
    b = audit_boundary_and_support(result)
    if not b["active"]:
        return ""
    config = result.get("config") if isinstance(result, dict) and isinstance(result.get("config"), dict) else {}
    method = str(config.get("method", "JMLE")).upper()
    if b["unsupported_internal"]:
        return CATEGORY_UNSUPPORTED_HOLD_REASON
    if b["unbounded_facet_levels"]:
        if method == "MML":
            return MML_FACET_BOUNDARY_HOLD_REASON
        from mfrm_app.jmle_boundary import HOLD_REASON
        return HOLD_REASON
    if b["unobserved_owner_endpoints"]:
        return UNOBSERVED_OWNER_ENDPOINT_HOLD_REASON
    if b["unobserved_global_boundary"]:
        return UNOBSERVED_BOUNDARY_CATEGORY_HOLD_REASON
    if method == "MML":
        if int(b["hypergraph_components"]) > 1 and (not for_screening_block or b.get("disconnected_levels") is None):
            return MML_DISCONNECTED_HOLD_REASON
        if not for_screening_block and b["disconnected_facets"]:
            return MML_DISCONNECTED_HOLD_REASON
    return ""


def boundary_and_support_hold_message(reason: str) -> str:
    if reason == CATEGORY_UNSUPPORTED_HOLD_REASON:
        return (
            "One or more internal score categories (or step-facet levels with fewer than two observed categories) "
            "have zero observations; adjacent step thresholds are not empirically identified without recoding or collapsing."
        )
    if reason == UNOBSERVED_OWNER_ENDPOINT_HOLD_REASON:
        return (
            "A PCM/GPCM step-facet level has zero observations in an endpoint rubric category; "
            "its endpoint step parameter has no finite unpenalized estimate."
        )
    if reason == UNOBSERVED_BOUNDARY_CATEGORY_HOLD_REASON:
        return (
            "An endpoint rubric category (minimum or maximum) has zero observations across the retained data; "
            "the corresponding boundary step threshold is unbounded without recoding the active score range."
        )
    if reason == MML_FACET_BOUNDARY_HOLD_REASON:
        return (
            "A free MML facet level has only minimum or only maximum scores along an unpenalized recession direction; "
            "optimizer termination does not establish a finite structural estimate."
        )
    if reason == MML_DISCONNECTED_HOLD_REASON:
        return (
            "Rating cohorts are disconnected across unanchored facets; cross-cohort MML scale alignment "
            "relies solely on the assumed common person distribution rather than shared observations or anchors."
        )
    return ""


def _spec_has_explicit_anchor(spec) -> bool:
    if not isinstance(spec, dict):
        return False
    anchors = spec.get("anchors")
    if anchors is not None and np.isfinite(np.asarray(anchors, dtype=float)).any():
        return True
    groups = spec.get("groups") or ()
    return any(g not in (None, "") and not (isinstance(g, float) and np.isnan(g)) for g in groups)


def identifiability_audit_state(result) -> str:
    """Distinguish proven rank deficiency from an unverified size-limited dense audit."""
    result = result or {}
    audit = result.get("identifiability")
    if not isinstance(audit, dict):
        return ""
    status = str(audit.get("status") or "")
    if status == "audit_size_limit":
        conn = audit.get("connectivity")
        config = result.get("config") if isinstance(result.get("config"), dict) else {}
        person_anchored = _spec_has_explicit_anchor(config.get("theta_spec"))
        facet_specs = config.get("facet_specs") if isinstance(config.get("facet_specs"), dict) else {}
        if conn is not None and hasattr(conn, "iterrows") and {"Facet", "Components"}.issubset(getattr(conn, "columns", ())):
            for _, row in conn.iterrows():
                try:
                    comp = int(row["Components"])
                except (TypeError, ValueError):
                    comp = 1
                if comp > 1 and not person_anchored and not _spec_has_explicit_anchor(facet_specs.get(str(row["Facet"]))):
                    return "eta_rank_deficient"
    return status


def mml_inference_hold_reason(config) -> str:
    if not isinstance(config, dict) or str(config.get("method", "")).upper() != "MML":
        return ""
    return (FREE_SD_MML_INFERENCE_HOLD_REASON if config.get("estimate_population_sd", False)
            else FIXED_SD_MML_INFERENCE_HOLD_REASON)


def mml_inference_hold_message(config) -> str:
    if not mml_inference_hold_reason(config):
        return ""
    if config.get("estimate_population_sd", False):
        return "Free-population-SD MML numerical accuracy and sampling uncertainty remain unqualified."
    return (
        "Fixed-population-SD MML general measure and reliability tables do not establish formal inference "
        "(SupportsFormalInference=FALSE; general SEs and intervals are screening_only)."
    )


def jmle_inference_hold_reason(result) -> str:
    result = result or {}
    config = result.get("config", {})
    if (isinstance(config, dict) and str(config.get("method", "")).upper() == "JMLE"
            and config.get("jmle_facet_boundary_review", False)):
        from mfrm_app.jmle_boundary import HOLD_REASON
        return HOLD_REASON
    state = identifiability_audit_state(result)
    if state == "eta_rank_deficient":
        return JMLE_RANK_DEFICIENT_HOLD_REASON
    if state == "audit_size_limit":
        return JMLE_AUDIT_SIZE_LIMIT_HOLD_REASON
    return ""


def jmle_inference_hold_message(result) -> str:
    reason = jmle_inference_hold_reason(result)
    if not reason:
        return ""
    if reason == JMLE_RANK_DEFICIENT_HOLD_REASON:
        return (
            "JMLE linear predictor is structurally rank-deficient (unidentified); "
            "add cross-facet connectivity or anchors before interpreting estimates or uncertainty."
        )
    if reason == JMLE_AUDIT_SIZE_LIMIT_HOLD_REASON:
        return (
            "Exact dense JMLE rank audit exceeded its matrix size limit (identification unverified by dense SVD); "
            "formal inference is not established, and conditional SEs/intervals are screening_only (SupportsFormalInference=FALSE)."
        )
    from mfrm_app.jmle_boundary import HOLD_MESSAGE
    return HOLD_MESSAGE


def optimizer_termination_established(result) -> bool:
    """Require an explicitly recorded Boolean; truthy legacy values are not evidence."""
    stopped = getattr((result or {}).get("opt"), "success", None)
    return isinstance(stopped, (bool, np.bool_)) and bool(stopped)


def calibration_inference_hold_reason(result) -> str:
    """Requalify current and saved results, separating non-convergence, rank deficiency and audit limits."""
    result = result or {}
    config = result.get("config", {})
    if isinstance(config, dict) and str(config.get("method", "")).upper() == "MML" and config.get("estimate_population_sd", False):
        return FREE_SD_MML_INFERENCE_HOLD_REASON
    jmle_reason = jmle_inference_hold_reason(result)
    if jmle_reason and jmle_reason != JMLE_AUDIT_SIZE_LIMIT_HOLD_REASON:
        return jmle_reason
    if "opt" in result and not optimizer_termination_established(result):
        return OPTIMIZER_TERMINATION_HOLD_REASON
    if b_reason := boundary_and_support_hold_reason(result, for_screening_block=False):
        return b_reason
    if reason := mml_inference_hold_reason(config):
        return reason
    if jmle_reason:
        return jmle_reason
    return "" if optimizer_termination_established(result) else OPTIMIZER_TERMINATION_HOLD_REASON


def calibration_inference_hold_message(result) -> str:
    result = result or {}
    reason = calibration_inference_hold_reason(result)
    if not reason:
        return ""
    if reason == OPTIMIZER_TERMINATION_HOLD_REASON:
        return "Successful optimizer termination has not been established; review the numerical diagnostics."
    if msg := boundary_and_support_hold_message(reason):
        return msg
    if reason in {FREE_SD_MML_INFERENCE_HOLD_REASON, FIXED_SD_MML_INFERENCE_HOLD_REASON}:
        return mml_inference_hold_message(result.get("config"))
    return jmle_inference_hold_message(result)


def calibration_screening_block_reason(result) -> str:
    """Return a non-empty reason only when even exploratory screening_only SEs/intervals must be withheld."""
    result = result or {}
    if "opt" in result and not optimizer_termination_established(result):
        return OPTIMIZER_TERMINATION_HOLD_REASON
    config = result.get("config", {})
    if isinstance(config, dict) and str(config.get("method", "")).upper() == "MML" and config.get("estimate_population_sd", False):
        return FREE_SD_MML_INFERENCE_HOLD_REASON
    if (isinstance(config, dict) and str(config.get("method", "")).upper() == "JMLE"
            and config.get("jmle_facet_boundary_review", False)):
        from mfrm_app.jmle_boundary import HOLD_REASON
        return HOLD_REASON
    if identifiability_audit_state(result) == "eta_rank_deficient":
        return JMLE_RANK_DEFICIENT_HOLD_REASON
    if b_reason := boundary_and_support_hold_reason(result, for_screening_block=True):
        return b_reason
    return "" if optimizer_termination_established(result) else OPTIMIZER_TERMINATION_HOLD_REASON


def calibration_screening_block_message(result) -> str:
    reason = calibration_screening_block_reason(result)
    if not reason:
        return ""
    if reason == OPTIMIZER_TERMINATION_HOLD_REASON:
        return "Successful optimizer termination has not been established; review the numerical diagnostics."
    if reason == FREE_SD_MML_INFERENCE_HOLD_REASON:
        return "Free-population-SD MML numerical accuracy and sampling uncertainty remain unqualified."
    if reason == JMLE_RANK_DEFICIENT_HOLD_REASON:
        return jmle_inference_hold_message(result)
    if msg := boundary_and_support_hold_message(reason):
        return msg
    from mfrm_app.jmle_boundary import HOLD_MESSAGE
    return HOLD_MESSAGE


def run_qualification_state(result) -> dict[str, object]:
    """Shared run-level diagnostic state across tables, figures, HTML, CSV, and reloaded archives."""
    result = result or {}
    converged = optimizer_termination_established(result)
    ident_state = identifiability_audit_state(result)
    formal_reason = calibration_inference_hold_reason(result)
    formal_message = calibration_inference_hold_message(result)
    block_reason = calibration_screening_block_reason(result)
    block_message = calibration_screening_block_message(result)
    config = result.get("config", {}) if isinstance(result.get("config"), dict) else {}
    method = str(config.get("method", "")).upper()
    if block_reason == OPTIMIZER_TERMINATION_HOLD_REASON:
        state = "non_converged"
    elif block_reason == JMLE_RANK_DEFICIENT_HOLD_REASON:
        state = "unidentifiable"
    elif block_reason == FREE_SD_MML_INFERENCE_HOLD_REASON:
        state = "free_sd_unqualified"
    elif block_reason in {
        CATEGORY_UNSUPPORTED_HOLD_REASON,
        UNOBSERVED_BOUNDARY_CATEGORY_HOLD_REASON,
        UNOBSERVED_OWNER_ENDPOINT_HOLD_REASON,
    }:
        state = "unsupported_category_structure"
    elif block_reason == MML_DISCONNECTED_HOLD_REASON or formal_reason == MML_DISCONNECTED_HOLD_REASON:
        state = "disconnected_mml_cohorts"
    elif block_reason:
        state = "facet_score_boundary"
    elif ident_state == "audit_size_limit":
        state = "audit_size_limit_unverified"
    elif method == "MML":
        state = "fixed_sd_mml_exploratory"
    else:
        state = "jmle_conditional_screening"
    return {
        "state": state,
        "converged": converged,
        "identifiability_state": ident_state,
        "formal_hold_reason": formal_reason,
        "formal_hold_message": formal_message,
        "screening_block_reason": block_reason,
        "screening_block_message": block_message,
        "supports_formal_inference": False,
    }


class OutputUse(str, Enum):
    """Allowed use of an output under the active remediation baseline."""

    READY = "READY"
    CONDITIONAL = "CONDITIONAL"
    SCREENING_ONLY = "SCREENING_ONLY"
    TECHNICAL_ONLY = "TECHNICAL_ONLY"
    WITHHELD = "WITHHELD"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class QualificationReason(str, Enum):
    """Stable reason codes linked to the statistical remediation ledger."""

    PARAMETER_COUNT_UNIDENTIFIED = "stat.stat_002.parameter_count_unidentified"
    MODEL_SELECTION_SCOPE_UNVALIDATED = (
        "model.model_001.selection_scope_unvalidated"
    )
    JMLE_AUTO_RECOMMENDATION_WITHHELD = (
        "model.model_001.jmle_auto_recommendation_withheld"
    )
    RSM_PCM_LRT_REGULARITY_UNVALIDATED = (
        "model.model_001.rsm_pcm_lrt_regularity_unvalidated"
    )
    GPCM_LRT_UNCALIBRATED = "model.model_002.gpcm_lrt_uncalibrated"
    COVARIANCE_RANK_DEFICIENT = "stat.stat_003.covariance_rank_deficient"
    COVARIANCE_COVERAGE_UNVALIDATED = (
        "stat.stat_003.covariance_coverage_unvalidated"
    )
    COVARIANCE_UNAVAILABLE = "stat.stat_003.covariance_unavailable"
    COVARIANCE_NOT_APPLICABLE = "stat.stat_003.covariance_not_applicable"
    POPULATION_SD_UNCERTAINTY_UNQUALIFIED = (
        "stat.stat_003.population_sd_uncertainty_unqualified"
    )
    BIAS_CELL_UNCALIBRATED = "bias.bias_003.cell_reference_uncalibrated"
    BIAS_PAIRWISE_UNVERIFIED = "bias.bias_001_002.pairwise_unverified"
    FACET_EQUIVALENCE_UNVALIDATED = "stat.stat_003.facet_equivalence_unvalidated"


@dataclass(frozen=True)
class OutputQualification:
    """One exportable qualification decision for a named output."""

    output_id: str
    status: OutputUse
    reason_code: QualificationReason
    roadmap_issue_id: str
    next_gate: str
    raw_technical_export_allowed: bool
    public_conclusion_allowed: bool
    user_action: str
    schema_version: str = OUTPUT_QUALIFICATION_SCHEMA_VERSION

    def to_record(self) -> dict[str, object]:
        record = asdict(self)
        record["status"] = self.status.value
        record["reason_code"] = self.reason_code.value
        return {
            "SchemaVersion": record["schema_version"],
            "OutputID": record["output_id"],
            "QualificationStatus": record["status"],
            "ReasonCode": record["reason_code"],
            "RoadmapIssueID": record["roadmap_issue_id"],
            "NextGate": record["next_gate"],
            "RawTechnicalExportAllowed": record["raw_technical_export_allowed"],
            "PublicConclusionAllowed": record["public_conclusion_allowed"],
            "UserAction": record["user_action"],
        }


def _qualification(
    output_id: str,
    status: OutputUse,
    reason: QualificationReason,
    issue: str,
    gate: str,
    action: str,
    *,
    raw_export: bool = True,
    public_conclusion: bool = False,
) -> OutputQualification:
    return OutputQualification(
        output_id=output_id,
        status=status,
        reason_code=reason,
        roadmap_issue_id=issue,
        next_gate=gate,
        raw_technical_export_allowed=raw_export,
        public_conclusion_allowed=public_conclusion,
        user_action=action,
    )


def model_choice_qualifications(method: object) -> tuple[OutputQualification, ...]:
    """Qualification matrix for IC, recommendations, and nested LR outputs."""
    method_label = str(method or "").upper()
    recommendation_reason = (
        QualificationReason.JMLE_AUTO_RECOMMENDATION_WITHHELD
        if method_label == "JMLE"
        else QualificationReason.MODEL_SELECTION_SCOPE_UNVALIDATED
    )
    recommendation_issue = "MODEL-001"
    return (
        _qualification(
            "model_choice.information_criteria",
            OutputUse.TECHNICAL_ONLY,
            QualificationReason.MODEL_SELECTION_SCOPE_UNVALIDATED,
            "MODEL-001",
            "G2",
            "Exact constraint coordinates are counted; observed-design identification and estimator/comparison scope still require validation before ranking models.",
        ),
        _qualification(
            "model_choice.automatic_recommendation",
            OutputUse.WITHHELD,
            recommendation_reason,
            recommendation_issue,
            "G2",
            "Do not generate or report an automatic preferred-model conclusion.",
        ),
        _qualification(
            "model_choice.lrt.rsm_pcm",
            OutputUse.WITHHELD,
            QualificationReason.RSM_PCM_LRT_REGULARITY_UNVALIDATED,
            "MODEL-001",
            "G2",
            "Do not report a chi-square p-value or retain/reject decision.",
        ),
        _qualification(
            "model_choice.lrt.pcm_gpcm",
            OutputUse.WITHHELD,
            QualificationReason.GPCM_LRT_UNCALIBRATED,
            "MODEL-002",
            "G2",
            "Do not use an ordinary chi-square reference distribution.",
        ),
        _qualification(
            "model_choice.lrt.rsm_gpcm",
            OutputUse.WITHHELD,
            QualificationReason.GPCM_LRT_UNCALIBRATED,
            "MODEL-002",
            "G2",
            "Do not use an ordinary chi-square reference distribution.",
        ),
    )


def covariance_qualification(
    *,
    status: object,
    rank: object = None,
    param_count: object = None,
) -> OutputQualification:
    """Qualify MML structural covariance without changing retained numerics."""
    status_label = str(status or "not_available").lower()
    try:
        rank_value = int(rank)
        param_value = int(param_count)
    except (TypeError, ValueError, OverflowError):
        rank_value = 0
        param_value = 0

    if status_label == "not_applicable":
        return _qualification(
            "mml.structural_covariance",
            OutputUse.NOT_APPLICABLE,
            QualificationReason.COVARIANCE_NOT_APPLICABLE,
            "STAT-003",
            "G2",
            "Use the estimation-method-specific uncertainty status instead.",
            raw_export=False,
        )
    # A missing/skipped Hessian does not establish a deficient rank.
    if (
        status_label == "regularized"
        or (status_label == "ok" and param_value > 0 and rank_value < param_value)
    ):
        return _qualification(
            "mml.structural_covariance",
            OutputUse.WITHHELD,
            QualificationReason.COVARIANCE_RANK_DEFICIENT,
            "STAT-003",
            "G1",
            "Do not use the regularized covariance for confirmatory SE or CI claims.",
        )
    if status_label == "ok" and param_value > 0 and rank_value == param_value:
        return _qualification(
            "mml.structural_covariance",
            OutputUse.TECHNICAL_ONLY,
            QualificationReason.COVARIANCE_COVERAGE_UNVALIDATED,
            "STAT-003",
            "G2",
            "Retain for audit; wait for coverage validation before confirmatory use.",
        )
    return _qualification(
        "mml.structural_covariance",
        OutputUse.WITHHELD,
        QualificationReason.COVARIANCE_UNAVAILABLE,
        "STAT-003",
        "G2",
        "Do not claim structural covariance-based SE or CI.",
        raw_export=False,
    )


def population_sd_uncertainty_qualification() -> OutputQualification:
    """Withhold free-SD inference, including legacy fixed-nuisance SE/CI."""
    return _qualification(
        "mml.population_sd_uncertainty",
        OutputUse.WITHHELD,
        QualificationReason.POPULATION_SD_UNCERTAINTY_UNQUALIFIED,
        "STAT-003",
        "G2",
        "Do not report population-SD SE/CI until joint stationarity, "
        "nuisance-adjusted uncertainty, and coverage are qualified. "
        "Fixed-nuisance curvature is a technical diagnostic, not profile information.",
    )


def bias_cell_qualification() -> OutputQualification:
    """Retain scalar likelihood diagnostics without promoting uncalibrated claims."""
    return _qualification(
        "bias.cell_inference", OutputUse.TECHNICAL_ONLY,
        QualificationReason.BIAS_CELL_UNCALIBRATED, "BIAS-003", "G3",
        "Conditional LR and cutoff intervals hold calibration fixed. "
        "Wait for error-rate, power, coverage, and calibration-uncertainty validation "
        "before reporting bias or no-bias conclusions.",
    )


def bias_pairwise_qualification() -> OutputQualification:
    """Keep conditional local contrasts technical until G3 inference is qualified."""
    return _qualification(
        "bias.pairwise_local_measure",
        OutputUse.WITHHELD,
        QualificationReason.BIAS_PAIRWISE_UNVERIFIED,
        "BIAS-001;BIAS-002",
        "G3",
        "Use conditional pairwise contrasts for technical audit only. "
        "Boundary behavior, joint calibration uncertainty, and reference-distribution "
        "and error-rate validation remain unqualified.",
    )


def facet_equivalence_qualification() -> OutputQualification:
    """Withhold equivalence claims until their uncertainty and decision rules are validated."""
    return _qualification(
        "facet.equivalence", OutputUse.WITHHELD,
        QualificationReason.FACET_EQUIVALENCE_UNVALIDATED, "STAT-003", "G2",
        "Specify the contrast, substantive equivalence margin and claim family. "
        "Validate joint uncertainty and the decision rule before reporting equivalence; "
        "nonsignificance and agreement among unvalidated methods do not establish it.",
        raw_export=False,
    )


def records(
    qualifications: Iterable[OutputQualification],
) -> list[dict[str, object]]:
    """Convert qualifications to a stable table-ready record sequence."""
    return [qualification.to_record() for qualification in qualifications]
