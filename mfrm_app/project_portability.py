"""UI-independent bidirectional `.mfrm` portability between Web and Desktop.

Exports eligible Web MFRM and observed-score G/D-Study runs to validated
schema-2 `.mfrm` archives without refitting, and loads/validates `.mfrm`
bytes in Web without running any estimator.
"""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from mfrm_app import basic_analysis as basic, desktop_outputs, engine, project


_MML_ENGINE_DISPLAY = {"em": "EM", "hybrid": "Hybrid", "direct": "Direct"}


def _extract_anchors_from_config(cfg: dict, *, person_col: str, facet_cols: list[str]) -> tuple[list[dict], list[dict]]:
    fixed_raw, group_raw = [], []
    for facet in ["Person", *facet_cols]:
        spec = cfg.get("theta_spec") if facet == "Person" else cfg.get("facet_specs", {}).get(facet)
        if not isinstance(spec, dict):
            continue
        levels = list(map(str, spec.get("levels", [])))
        anchors = spec.get("anchors", [])
        for idx, val in enumerate(anchors):
            if idx < len(levels) and val is not None and np.isfinite(val):
                fixed_raw.append({"Facet": facet, "Level": levels[idx], "Anchor": float(val)})
        groups = spec.get("groups", [])
        gvals = spec.get("group_values", {}) or {}
        for idx, gname in enumerate(groups):
            if idx < len(levels) and gname not in (None, "") and not pd.isna(gname):
                g_str = str(gname)
                if g_str in gvals and np.isfinite(gvals[g_str]):
                    group_raw.append({
                        "Facet": facet,
                        "Level": levels[idx],
                        "Group": g_str,
                        "GroupValue": float(gvals[g_str]),
                    })
    norm_a = basic.normalize_anchor_records(fixed_raw, person_col=person_col, facet_cols=facet_cols)
    norm_g = basic.normalize_group_anchor_records(
        group_raw,
        person_col=person_col,
        facet_cols=facet_cols,
        fixed_pairs={(r["Facet"], r["Level"]) for r in norm_a},
    )
    return norm_a, norm_g


def check_web_mfrm_portability(
    result: dict,
    *,
    input_data: pd.DataFrame | None = None,
    person_col: str | None = None,
    score_col: str | None = None,
) -> dict:
    """Return `{'eligible': bool, 'reasons': list[str]}` for `.mfrm` export."""
    reasons: list[str] = []
    if not isinstance(result, dict):
        return {"eligible": False, "reasons": ["No fitted MFRM result is available."]}
    cfg = result.get("config")
    prep = result.get("prep")
    if not isinstance(cfg, dict) or not isinstance(prep, dict):
        return {"eligible": False, "reasons": ["Result configuration or prepared data record is incomplete."]}

    model = str(cfg.get("model", "")).upper()
    method = str(cfg.get("method", "")).upper()
    if model not in ("RSM", "PCM"):
        reasons.append(f"Model '{model}' is a Web research option; `.mfrm` projects support RSM and PCM.")
    if method not in ("JMLE", "MML"):
        reasons.append(f"Estimation method '{method}' is outside the JMLE / fixed-SD MML `.mfrm` archive contract.")
    if bool(cfg.get("estimate_population_sd")):
        reasons.append("Free population SD MML is on research hold and outside the `.mfrm` archive contract.")
    if method == "MML":
        c_eng = str(cfg.get("mml_engine", "")).lower()
        if c_eng not in _MML_ENGINE_DISPLAY:
            reasons.append(f"MML engine '{c_eng}' is outside the fixed-SD MML archive contract.")
        qp = cfg.get("quad_points")
        prior_sd = cfg.get("population_prior_sd")
        if not isinstance(qp, int) or qp < 3 or prior_sd is None or not np.isfinite(prior_sd) or float(prior_sd) <= 0:
            reasons.append("MML quadrature points (>= 3) or positive finite fixed population SD is missing.")

    if cfg.get("population_formula") or cfg.get("population_data_fingerprint"):
        reasons.append("Latent-regression covariates are outside the Desktop `.mfrm` archive contract.")
    if cfg.get("weight_col"):
        reasons.append("Differential observation weights are outside the unit-weight `.mfrm` archive contract.")

    prep_data = prep.get("data")
    if not isinstance(prep_data, pd.DataFrame) or prep_data.empty:
        reasons.append("Retained observation table (`prep['data']`) is missing.")
    elif "Weight" in prep_data.columns:
        w = pd.to_numeric(prep_data["Weight"], errors="coerce")
        if w.isna().any() or not np.allclose(w.to_numpy(dtype=float), 1.0):
            reasons.append("Observation weights differ from 1.0; `.mfrm` projects require unit weights.")

    specs = [cfg.get("theta_spec"), *(cfg.get("facet_specs") or {}).values()]
    for spec in specs:
        if isinstance(spec, dict):
            if bool(spec.get("regularized")) or float(spec.get("penalty_precision", 0.0) or 0.0) != 0.0:
                reasons.append("Facet shrinkage/regularization is a Web research option outside `.mfrm` archives.")
                break
            if bool(spec.get("dummy")):
                reasons.append("Dummy-facet suppression (`Status D`) is outside the Desktop `.mfrm` archive contract.")
                break

    if cfg.get("noncenter_facet", "Person") != "Person":
        reasons.append("Desktop `.mfrm` archives require `Person` as the non-centered reference facet.")
    signs = cfg.get("facet_signs") or {}
    if any(int(v) != -1 for v in signs.values()):
        reasons.append("Desktop `.mfrm` archives use standard `-1` subtractive facet orientation for non-Person facets.")

    rmin, rmax = cfg.get("rating_min"), cfg.get("rating_max")
    if type(rmin) is not int or type(rmax) is not int or rmin >= rmax:
        reasons.append("Integer minimum and maximum score bounds (`rating_min < rating_max`) are required.")

    score_map = prep.get("score_map")
    if isinstance(score_map, pd.DataFrame) and {"OriginalScore", "InternalScore"}.issubset(score_map.columns):
        orig = pd.to_numeric(score_map["OriginalScore"], errors="coerce")
        internal = pd.to_numeric(score_map["InternalScore"], errors="coerce")
        if orig.isna().any() or internal.isna().any() or not np.array_equal(orig.to_numpy(), internal.to_numpy()):
            reasons.append(
                "Unobserved score gaps were collapsed during recoding; enable 'Keep original category values' "
                "before exporting a `.mfrm` archive."
            )

    facet_names = list(cfg.get("facet_names") or [])
    if not facet_names:
        reasons.append("At least one non-Person facet is required.")
    if model == "PCM" and cfg.get("step_facet") not in facet_names:
        reasons.append("PCM `.mfrm` archives require a valid `step_facet`.")

    if person_col and score_col:
        cols = [str(person_col), *facet_names, str(score_col)]
        if len(set(cols)) != len(cols):
            reasons.append("Selected person, facet, and score column names must be distinct.")
        elif input_data is not None and any(c not in input_data.columns for c in cols):
            reasons.append("Selected person, facet, and score columns must exist distinctly in the input table.")

    return {"eligible": len(reasons) == 0, "reasons": reasons}


def build_mfrm_project_from_web(
    result: dict,
    diagnostics: dict,
    report_tables: dict | None = None,
    scorefile: pd.DataFrame | None = None,
    *,
    input_data: pd.DataFrame | None = None,
    person_col: str | None = None,
    score_col: str | None = None,
    all_bias_results: dict | None = None,
    umean: float = 0.0,
    uscale: float = 1.0,
) -> bytes:
    """Serialize an eligible Web MFRM run into validated schema-2 `.mfrm` bytes without refitting."""
    chk = check_web_mfrm_portability(
        result, input_data=input_data, person_col=person_col, score_col=score_col
    )
    if not chk["eligible"]:
        raise project.ProjectError("; ".join(chk["reasons"]))

    cfg = result["config"]
    prep = result["prep"]
    facet_cols = list(cfg["facet_names"])
    if isinstance(input_data, pd.DataFrame) and person_col and score_col and all(
        c in input_data.columns for c in [person_col, *facet_cols, score_col]
    ):
        archive_input = input_data.copy(deep=True)
        p_col, s_col = str(person_col), str(score_col)
    elif person_col and score_col and len({str(person_col), *facet_cols, str(score_col)}) == 2 + len(facet_cols):
        p_col, s_col = str(person_col), str(score_col)
        archive_input = (
            prep["data"][["Person", *facet_cols, "Score"]]
            .rename(columns={"Person": p_col, "Score": s_col})
            .copy(deep=True)
        )
    else:
        p_col, s_col = "Person", "Score"
        archive_input = prep["data"][["Person", *facet_cols, "Score"]].copy(deep=True)

    f_umean = float(umean)
    f_uscale = float(uscale)
    if not np.isfinite(f_umean) or not np.isfinite(f_uscale) or f_uscale <= 0:
        raise project.ProjectError("umean must be finite and uscale must be positive finite")
    if f_umean == 0.0:
        f_umean = 0.0

    compute_pca = bool(diagnostics.get("pca_enabled", False)) if isinstance(diagnostics, dict) else False
    norm_anchors, norm_groups = _extract_anchors_from_config(cfg, person_col=p_col, facet_cols=facet_cols)

    fit_settings: dict = {
        "person_col": p_col,
        "facet_cols": facet_cols,
        "score_col": s_col,
        "rating_min": int(cfg["rating_min"]),
        "rating_max": int(cfg["rating_max"]),
        "maxit": int(cfg.get("maxit", 400)),
        "reltol": float(cfg.get("reltol", 1e-6)),
        "model": str(cfg["model"]),
        "method": str(cfg["method"]),
        "keep_original": True,
    }
    if cfg["model"] == "PCM":
        fit_settings["step_facet"] = str(cfg["step_facet"])
    if cfg["method"] == "MML":
        fit_settings.update(
            quad_points=int(cfg["quad_points"]),
            population_prior_sd=float(cfg["population_prior_sd"]),
            estimate_population_sd=False,
            mml_engine=_MML_ENGINE_DISPLAY[str(cfg["mml_engine"]).lower()],
        )
    if norm_anchors:
        fit_settings["anchors"] = norm_anchors
    if norm_groups:
        fit_settings["group_anchors"] = norm_groups
    if compute_pca:
        fit_settings["compute_pca"] = True
    if (f_umean, f_uscale) != (0.0, 1.0):
        fit_settings["umean"] = f_umean
        fit_settings["uscale"] = f_uscale

    clean_result = {k: v for k, v in result.items() if not str(k).startswith("_")}
    analysis_id = engine.build_result_analysis_identity(clean_result).analysis_id

    if not isinstance(report_tables, dict) or not report_tables:
        tables = engine.calc_facets_report_tbls(clean_result, diagnostics, umean=f_umean, uscale=f_uscale)
    else:
        # Ensure stored report_tables match (umean, uscale) with zero refitting.
        tables = engine.calc_facets_report_tbls(clean_result, diagnostics, umean=f_umean, uscale=f_uscale)

    sfile = scorefile if isinstance(scorefile, pd.DataFrame) else engine.compute_scorefile(clean_result)

    bias_entries = []
    if isinstance(all_bias_results, dict):
        allowed_facets = {"Person", *facet_cols}
        for bres in all_bias_results.values():
            if not isinstance(bres, dict):
                continue
            fa, fb = bres.get("facet_a"), bres.get("facet_b")
            if not fa or not fb:
                btbl = bres.get("table")
                if isinstance(btbl, pd.DataFrame) and not btbl.empty and {"FacetA", "FacetB"}.issubset(btbl.columns):
                    fa = str(btbl["FacetA"].iloc[0])
                    fb = str(btbl["FacetB"].iloc[0])
            if fa in allowed_facets and fb in allowed_facets and fa != fb:
                b_settings = {
                    "facet_a": str(fa),
                    "facet_b": str(fb),
                    "max_abs": float(bres.get("max_abs", 10.0)),
                    "omit_extreme": bool(bres.get("omit_extreme", True)),
                }
                btbl = bres.get("table")
                if (isinstance(btbl, pd.DataFrame) and not btbl.empty) or bres.get("_skip_reason"):
                    bias_entries.append({
                        "analysis_id": analysis_id,
                        "settings": b_settings,
                        "result": deepcopy(bres),
                    })

    interval_entries = []
    si_rec = result.get("_structural_intervals")
    if isinstance(si_rec, dict) and cfg["method"] == "MML":
        targets = si_rec.get("targets") or []
        if targets:
            t_facet = str(targets[0]["Facet"])
            diffs = [r for r in targets if r.get("Kind") == "difference"]
            diff_spec = [str(diffs[0]["Level"]), str(diffs[0]["Reference"])] if diffs else None
            if t_facet in facet_cols:
                interval_entries.append({
                    "analysis_id": analysis_id,
                    "settings": {"facet": t_facet, "difference": diff_spec, "level": float(si_rec["level"])},
                    "record": deepcopy(si_rec),
                    "status": "computed",
                    "reason": "",
                })

    joint_entry = None
    ji_rec = result.get("_joint_interaction")
    if isinstance(ji_rec, dict) and cfg["method"] == "MML" and isinstance(ji_rec.get("settings"), dict):
        joint_entry = {
            "analysis_id": analysis_id,
            "settings": deepcopy(ji_rec["settings"]),
            "record": deepcopy(ji_rec),
            "status": "computed",
            "reason": "",
        }

    has_extras = bool(
        norm_anchors or norm_groups or compute_pca or (f_umean, f_uscale) != (0.0, 1.0)
        or bias_entries or interval_entries or joint_entry is not None
        or cfg["method"] != "JMLE" or cfg["model"] != "RSM"
    )
    analysis = {
        "kind": "mfrm-analysis-v2" if has_extras else "basic-rsm-jmle-v1",
        "analysis_id": analysis_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source": basic._source_versions(),
        "input_data": archive_input,
        "fit_settings": fit_settings,
        "result": clean_result,
        "diagnostics": diagnostics,
        "report_tables": tables,
        "scorefile": sfile,
        "qualification_at_save": basic.qualification_snapshot(clean_result),
        "scope": {
            "includes_input": True,
            "optimizer_checkpoint": False,
            "pca": compute_pca,
            "marginal_fit": False,
            "bias_refits": bool(bias_entries),
        },
    }
    if bias_entries or interval_entries or has_extras:
        analysis["bias_analyses"] = bias_entries
        analysis["interval_analyses"] = interval_entries
    if joint_entry is not None:
        analysis["joint_analysis"] = joint_entry

    analysis["content_fingerprint"] = project.content_fingerprint(analysis)
    with tempfile.TemporaryDirectory(prefix="mfrm-web-export-") as td:
        out_path = Path(td) / "analysis.mfrm"
        project.save_project(out_path, analysis)
        return out_path.read_bytes()


def build_gd_project_from_web(
    bundle: dict,
    *,
    input_data: pd.DataFrame,
    n_max_per_facet: int = 10,
) -> bytes:
    """Serialize a Web observed-score G/D-Study bundle into validated `.mfrm` bytes without recomputation."""
    if not isinstance(bundle, dict) or not isinstance(input_data, pd.DataFrame) or input_data.empty:
        raise project.ProjectError("Valid G/D bundle and input table are required")
    design = bundle.get("design") or {}
    obj = design.get("object_facet")
    rfacets = list(design.get("random_facets") or [])
    scol = design.get("score_column") or design.get("score_col")
    if not obj or not rfacets or not scol:
        raise project.ProjectError("G/D design metadata is incomplete")

    entry_bundle = deepcopy(bundle)
    selected_counts = None
    if entry_bundle.get("available"):
        if isinstance(entry_bundle.get("selected_plan"), dict) and all(f in entry_bundle["selected_plan"] for f in rfacets):
            selected_counts = {f: int(entry_bundle["selected_plan"][f]) for f in rfacets}
            if selected_counts == dict(design.get("observed_n") or {}):
                selected_counts = None
        else:
            from mfrm_app import generalizability as gd_mod
            obs_counts = {f: int(design["observed_n"][f]) for f in rfacets}
            coefs = gd_mod._g_phi_from_variance_components(
                entry_bundle["variance_components"], object_facet=str(obj), facet_observations=obs_counts
            )
            entry_bundle["selected_plan"] = {**obs_counts, **gd_mod._g_d_error_row(coefs)}

    d_tbl = bundle.get("d_study")
    inferred_n_max = int(design.get("max_planned_count") or n_max_per_facet)
    if isinstance(d_tbl, pd.DataFrame) and not d_tbl.empty and rfacets[0] in d_tbl.columns:
        inferred_n_max = max(inferred_n_max, int(pd.to_numeric(d_tbl[rfacets[0]], errors="coerce").max()))

    gd_settings: dict = {
        "object_facet": str(obj),
        "random_facets": [str(f) for f in rfacets],
        "score_col": str(scol),
        "n_max_per_facet": inferred_n_max,
        "target_population": str(design.get("target_population", "")),
        "selected_plan": selected_counts,
    }
    if design.get("fixed_condition") is not None:
        gd_settings["fixed_condition"] = deepcopy(design["fixed_condition"])

    calc_ver = bundle["calculation_version"]
    gd_id = project.gd_analysis_identity(input_data, gd_settings, calc_ver)
    entry = {
        "settings": gd_settings,
        "result": entry_bundle,
        "analysis_id": gd_id,
    }
    analysis = {
        "kind": "gd-analysis-v2",
        "analysis_id": gd_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source": basic._source_versions(),
        "input_data": input_data.copy(deep=True),
        "gd_analysis": entry,
        "scope": {"includes_input": True, "optimizer_checkpoint": False, "mfrm_fitted": False},
    }
    analysis["content_fingerprint"] = project.content_fingerprint(analysis)
    with tempfile.TemporaryDirectory(prefix="mfrm-gd-export-") as td:
        out_path = Path(td) / "gd_analysis.mfrm"
        project.save_project(out_path, analysis)
        return out_path.read_bytes()


def load_project_from_bytes(payload: bytes) -> dict:
    """Validate and decode a `.mfrm` archive from raw bytes with zero estimation."""
    if not isinstance(payload, (bytes, bytearray)):
        raise project.ProjectError("Project payload must be raw bytes")
    if len(payload) > project.MAX_PROJECT_BYTES:
        raise project.ProjectError("Project exceeds the project size limit (128 MiB)")
    with tempfile.TemporaryDirectory(prefix="mfrm-web-import-") as td:
        in_path = Path(td) / "uploaded.mfrm"
        in_path.write_bytes(bytes(payload))
        analysis = project.load_project(in_path)
    return {
        "analysis": analysis,
        "tables": desktop_outputs.current_tables(analysis),
        "figure_choices": desktop_outputs.figure_choices(analysis),
        "request_settings": desktop_outputs.request_settings(analysis),
    }
