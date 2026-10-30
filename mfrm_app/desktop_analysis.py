"""Desktop calculation requests and retained results, independent of either UI."""
from copy import deepcopy
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mfrm_app import basic_analysis as basic, bias_analysis as bias, engine, reporting_facts
from mfrm_app import generalizability as gd, project, structural_intervals as intervals, joint_interaction as joint


def _bias_requests(requests, facets):
    out, seen = [], set()
    for request in requests or []:
        if not isinstance(request, dict) or set(request) - {"facet_a", "facet_b", "max_abs", "omit_extreme"}:
            raise ValueError("Unknown bias request fields")
        a, b = request.get("facet_a"), request.get("facet_b")
        if a not in facets or b not in facets or a == b or frozenset((a, b)) in seen:
            raise ValueError("Select each distinct fitted bias facet pair once")
        bound, omit = request.get("max_abs", 10.), request.get("omit_extreme", True)
        if isinstance(bound, bool) or not np.isfinite(bound) or bound <= 0 or type(omit) is not bool:
            raise ValueError("Bias search bound must be positive and extreme-level handling explicit")
        seen.add(frozenset((a, b)))
        out.append(dict(facet_a=a, facet_b=b, max_abs=float(bound), omit_extreme=omit))
    return out


def _interval_requests(requests, facets, method):
    out = []
    for request in requests or []:
        if not isinstance(request, dict) or set(request) - {"facet", "difference", "level"}:
            raise ValueError("Unknown interval request fields")
        facet, difference, level = request.get("facet"), request.get("difference"), request.get("level", .95)
        if method != "MML" or facet not in facets:
            raise ValueError("Specified intervals require a fitted structural facet and fixed-SD MML")
        if isinstance(level, bool) or not np.isfinite(level) or not 0 < level < 1:
            raise ValueError("Interval level must be between zero and one")
        if difference is not None and (not isinstance(difference, (list, tuple)) or len(difference) != 2
                                       or any(not isinstance(x, str) for x in difference)
                                       or difference[0] == difference[1]):
            raise ValueError("Select two distinct levels for an A-minus-B difference")
        out.append(dict(facet=facet, difference=list(difference) if difference is not None else None,
                        level=float(level)))
    return out


_RESERVED_GD_NAMES = gd._RESERVED_GD_NAMES
_validate_fixed_condition = gd.validate_fixed_condition
prepare_gd_input_data = gd.prepare_gd_input_data


def _gd_settings(settings):
    if not isinstance(settings, dict) or set(settings) - {
        "object_facet", "random_facets", "score_col", "n_max_per_facet",
        "target_population", "selected_plan", "fixed_condition",
    }:
        raise ValueError("Unknown G/D settings")
    raw = deepcopy(settings)
    fixed_raw = raw.pop("fixed_condition", None)
    out = dict(score_col="Score", n_max_per_facet=10, target_population="", selected_plan=None)
    out.update(raw)
    if (not isinstance(out.get("object_facet"), str)
            or not isinstance(out.get("random_facets"), (list, tuple))
            or any(not isinstance(f, str) for f in out["random_facets"])
            or not isinstance(out["score_col"], str) or not isinstance(out["target_population"], str)
            or type(out["n_max_per_facet"]) is not int or out["n_max_per_facet"] < 1):
        raise ValueError("G/D requires named columns, random facets and a positive integer plan limit")
    out["random_facets"] = list(out["random_facets"])
    selected = out["selected_plan"]
    if selected is not None and (not isinstance(selected, dict) or set(selected) != set(out["random_facets"])
                                 or any(type(n) is not int or n < 1 for n in selected.values())):
        raise ValueError("A selected D plan must specify positive integer counts for every random facet")
    fixed_validated = _validate_fixed_condition(
        fixed_raw, out["object_facet"], out["score_col"], out["random_facets"]
    )
    if fixed_validated is not None:
        out["fixed_condition"] = fixed_validated
    return out


def preview_mfrm_connectivity(
    data,
    *,
    person_col: str,
    facet_cols: list[str],
    score_col: str,
    anchors=None,
    group_anchors=None,
) -> dict:
    """Fast, optimization-free pre-flight check of connectivity, sparsity, and anchor linking."""
    facet_cols = list(facet_cols or [])
    if not isinstance(data, pd.DataFrame) or data.empty:
        return dict(
            ok=False, status="missing_data", nullity=None, density_pct=0.0, disconnected_facets=[],
            summary_ja="データを読み込んでください。",
            summary_en="Load a rating dataset first.",
        )
    cols = [person_col, *facet_cols, score_col]
    if not person_col or not score_col or not facet_cols or len(set(cols)) != len(cols) or any(c not in data.columns for c in cols):
        return dict(
            ok=False, status="incomplete_selection", nullity=None, density_pct=0.0, disconnected_facets=[],
            summary_ja="対象者・得点・1つ以上の相を重複なく選択すると、連結性（未接続サブセット）と充填率を診断します。",
            summary_en="Select distinct person, score and facet columns to inspect connectivity and cell density.",
        )
    audit = engine.build_response_data_audit(data, person_col=person_col, facet_cols=facet_cols, score_col=score_col)
    retained = data.loc[audit["rows"]["IncludedInLikelihood"].to_numpy(dtype=bool)]
    excluded = len(audit["excluded_rows"])
    if retained.empty or pd.to_numeric(retained[score_col], errors="coerce").nunique() < 2:
        return dict(
            ok=False, status="insufficient_data", nullity=None, density_pct=0.0, disconnected_facets=[],
            summary_ja="有効な評定行または得点段階（2種類以上）が不足しています。",
            summary_en="At least two distinct valid numeric score values are required.",
        )
    observed_levels = {
        "Person": set(retained[person_col].dropna().astype(str).str.strip()),
        **{f: set(retained[f].dropna().astype(str).str.strip()) for f in facet_cols},
    }
    try:
        norm_a = basic.normalize_anchor_records(
            anchors, person_col=person_col, facet_cols=facet_cols, observed_levels=observed_levels
        )
        norm_g = basic.normalize_group_anchor_records(
            group_anchors,
            person_col=person_col,
            facet_cols=facet_cols,
            observed_levels=observed_levels,
            fixed_pairs={(r["Facet"], r["Level"]) for r in norm_a},
        )
    except ValueError as exc:
        return dict(
            ok=False, status="anchor_error", nullity=None, density_pct=0.0, disconnected_facets=[],
            summary_ja=f"アンカー設定エラー：{exc}",
            summary_en=f"Anchor configuration error: {exc}",
        )

    prep = engine.prepare_mfrm_data(
        retained, person_col=person_col, facet_cols=facet_cols, score_col=score_col, keep_original=True
    )
    conn = engine._person_facet_connectivity_summary(prep)
    disconnected = (
        [f"{r.Facet} ({int(r.Components)})" for r in conn.itertuples(index=False) if int(r.Components) > 1]
        if not conn.empty else []
    )
    disconnected_names = (
        [str(r.Facet) for r in conn.itertuples(index=False) if int(r.Components) > 1]
        if not conn.empty else []
    )

    n_persons = max(1, len(prep["levels"]["Person"]))
    total_grid = n_persons
    for f in facet_cols:
        total_grid *= max(1, len(prep["levels"][f]))
    density = 100.0 * len(prep["data"]) / float(total_grid)

    anchor_df = pd.DataFrame(norm_a) if norm_a else None
    group_df = pd.DataFrame(norm_g) if norm_g else None
    specs = engine.prepare_constraint_specs(prep, anchor_df=anchor_df, group_anchor_df=group_df, noncenter_facet="Person")
    cfg_check = {
        "method": "JMLE",
        "model": "RSM",
        "n_cat": max(2, int(prep["rating_max"] - prep["rating_min"] + 1)),
        "facet_names": list(prep["facet_names"]),
        "facet_levels": {f: list(prep["levels"][f]) for f in prep["facet_names"]},
        "theta_spec": specs["theta_spec"],
        "facet_specs": specs["facet_specs"],
        "facet_signs": {f: -1 for f in prep["facet_names"]},
        "step_facet": None,
        "slope_facet": None,
    }
    ident = engine.build_jmle_eta_identifiability_audit(prep, cfg_check)
    nullity = (
        int(ident["summary"]["EtaStructuralNullity"].iloc[0])
        if ident.get("available") and not ident["summary"].empty
        else (len(disconnected_names))
    )

    excl_ja = f"・除外 {excluded:,}行" if excluded else ""
    excl_en = f", {excluded:,} excluded" if excluded else ""
    anch_count = len(norm_a) + len(norm_g)
    comp_counts = {str(r.Facet): int(r.Components) for r in conn.itertuples(index=False)} if not conn.empty else {}
    n_inc = len(prep["data"])

    if nullity > 0:
        disc_str = ", ".join(disconnected) if disconnected else "Person/Facet"
        return dict(
            ok=False,
            status="disconnected",
            nullity=nullity,
            included=n_inc,
            excluded=excluded,
            density_pct=density,
            component_counts=comp_counts,
            disconnected_facets=disconnected_names,
            summary_ja=(
                f"⚠ 未接続サブセット検出：Person × [{disc_str}] が分断されており、"
                f"線形予測子に {nullity} 自由度の未識別（階数欠損）があります（有効 {n_inc:,}行{excl_ja}・充填率 {density:.1f}%）。"
                "このままでは分断グループ間の厳しさ・能力を同一尺度で直接比較できません。共通評価者の追加、または「アンカー・等化」タブでの固定値指定を検討してください。"
            ),
            summary_en=(
                f"⚠ Disconnected subsets detected: Person × [{disc_str}] splits into separate components with "
                f"structural nullity = {nullity} ({n_inc:,} valid rows{excl_en}, {density:.1f}% cell density). "
                "Subsets cannot be directly compared on one scale without shared ratings or anchors on the Anchors / Equating tab."
            ),
        )

    if disconnected_names or anch_count > 0:
        obs_note_ja = (
            f"観測上の Person × ({', '.join(disconnected_names)}) は複数サブセットですが、指定アンカー（{anch_count}件）により単一ロジット尺度へ識別（階数欠損 0）されています"
            if disconnected_names else
            f"全相が単一ネットワークで連結され、指定アンカー（{anch_count}件）も整合しています（階数欠損 0）"
        )
        obs_note_en = (
            f"Observational Person × ({', '.join(disconnected_names)}) subsets are linked onto one identified logit scale via {anch_count} anchor entry/entries (nullity = 0)"
            if disconnected_names else
            f"All facets form a single connected network with {anch_count} active anchor entry/entries (nullity = 0)"
        )
        return dict(
            ok=True,
            status="anchored_ok",
            nullity=0,
            included=n_inc,
            excluded=excluded,
            density_pct=density,
            component_counts=comp_counts,
            disconnected_facets=disconnected_names,
            summary_ja=f"⚓ アンカー・連結確認OK：{obs_note_ja}（有効 {n_inc:,}行{excl_ja}・充填率 {density:.1f}%）。",
            summary_en=f"⚓ Anchor & connectivity OK: {obs_note_en} ({n_inc:,} valid rows{excl_en}, {density:.1f}% density).",
        )

    sparse_ja = "・スパース連結設計" if density < 60.0 else ""
    sparse_en = ", sparse linked design" if density < 60.0 else ""
    return dict(
        ok=True,
        status="ok",
        nullity=0,
        included=n_inc,
        excluded=excluded,
        density_pct=density,
        component_counts=comp_counts,
        disconnected_facets=[],
        summary_ja=(
            f"✓ 連結・識別確認OK：すべての相が単一ネットワークで連結されています"
            f"（有効 {n_inc:,}行{excl_ja}・充填率 {density:.1f}%{sparse_ja}・未接続サブセット 0・階数欠損 0）。"
        ),
        summary_en=(
            f"✓ Connectivity & identification OK: All facets form one connected network "
            f"({n_inc:,} valid rows{excl_en}, {density:.1f}% cell density{sparse_en}, 0 disconnected subsets, nullity = 0)."
        ),
    )


def _connectivity_and_anchor_table(analysis) -> tuple[str, pd.DataFrame, str] | None:
    """Build retained connectivity & anchor linking summary without recomputation."""
    if "result" not in analysis:
        return None
    settings = analysis.get("fit_settings", {})
    anchors = settings.get("anchors") or []
    group_anchors = settings.get("group_anchors") or []
    ident = analysis["result"].get("identifiability", {})
    conn = ident.get("connectivity", pd.DataFrame())
    has_disconnected = isinstance(conn, pd.DataFrame) and not conn.empty and (conn["Components"] > 1).any()
    if not anchors and not group_anchors and not has_disconnected:
        return None

    rows = []
    if isinstance(conn, pd.DataFrame) and not conn.empty:
        for r in conn.itertuples(index=False):
            rows.append({
                "Category": "Connectivity",
                "Facet": str(r.Facet),
                "Target": "Person × " + str(r.Facet),
                "Constraint": "Connected" if int(r.Components) == 1 else f"Disconnected ({int(r.Components)} subsets)",
                "Value": float(r.Components),
                "Detail": f"Min levels/person={int(r.MinimumRatersOrLevelsPerPerson)}, Min persons/level={int(r.MinimumPersonsPerFacetLevel)}",
            })
    for a in anchors:
        rows.append({
            "Category": "FixedAnchor",
            "Facet": a["Facet"],
            "Target": a["Level"],
            "Constraint": "Fixed logit (Status A)",
            "Value": float(a["Anchor"]),
            "Detail": "Element fixed at pre-calibrated logit value",
        })
    for g in group_anchors:
        rows.append({
            "Category": "GroupAnchor",
            "Facet": g["Facet"],
            "Target": g["Level"],
            "Constraint": f"Group {g['Group']} mean (Status G)",
            "Value": float(g["GroupValue"]),
            "Detail": f"Mean of group {g['Group']} fixed at {g['GroupValue']:g} logits",
        })
    note = (
        "Retained bipartite Person × Facet connectivity counts and explicit fixed/group anchor constraints "
        "used to identify the logit scale. Anchor equating assumes invariance across administrations; review fit and drift."
    )
    return ("Connectivity and Scale Linking", pd.DataFrame(rows), note)


def _blank_na_columns(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for col in out.columns[out.isna().any()]:
        out[col] = out[col].astype(object).where(out[col].notna(), "")
    return out


def _reliability_summary_table(analysis) -> tuple[str, pd.DataFrame, str] | None:
    if "result" not in analysis:
        return None
    diagnostics = analysis.get("diagnostics")
    if not isinstance(diagnostics, dict) or not isinstance(diagnostics.get("reliability"), pd.DataFrame):
        return None
    locale_path = Path(__file__).resolve().parents[1] / "locales" / "en.json"
    rel_note = json.loads(locale_path.read_text(encoding="utf-8"))["metric_notes"]["reliability"]
    table, note = reporting_facts.publication_reliability_table(
        analysis["result"], diagnostics, rel_note=rel_note
    )
    if table.empty:
        return None
    note = note.replace(
        "Full notes and original technical values are retained in the reliability CSV.",
        "Full notes and original technical values are retained in this project.",
    )
    return ("Separation and Reliability by Facet", table, note)


def _agreement_and_alpha_tables(analysis) -> list[tuple[str, pd.DataFrame, str]]:
    if "result" not in analysis:
        return []
    result = analysis["result"]
    diagnostics = analysis.get("diagnostics")
    if not isinstance(diagnostics, dict):
        return []
    ia_bundle = diagnostics.get("interrater_agreement")
    ia_summary = ia_bundle.get("summary") if isinstance(ia_bundle, dict) else None
    ia_pairs = ia_bundle.get("pairs_by_facet", {}) if isinstance(ia_bundle, dict) else {}
    alpha_tbl = diagnostics.get("krippendorff_alpha")
    has_ia = isinstance(ia_summary, pd.DataFrame) and not ia_summary.empty
    has_alpha = isinstance(alpha_tbl, pd.DataFrame) and not alpha_tbl.empty
    if not has_ia and not has_alpha:
        return []

    cal_hold = bool(reporting_facts.current_manuscript_hold(result))
    facet_names = list(result.get("config", {}).get("facet_names", []) or [])
    rows = []
    for facet in facet_names:
        ia_row = (
            ia_summary.loc[ia_summary["RaterFacet"].astype(str) == str(facet)].head(1)
            if has_ia and "RaterFacet" in ia_summary.columns
            else pd.DataFrame()
        )
        alpha_sub = (
            alpha_tbl.loc[alpha_tbl["Facet"].astype(str) == str(facet)]
            if has_alpha and "Facet" in alpha_tbl.columns
            else pd.DataFrame()
        )
        if ia_row.empty and alpha_sub.empty:
            continue

        def _alpha_val(metric_name: str) -> float:
            if alpha_sub.empty or "DistanceMetric" not in alpha_sub.columns:
                return np.nan
            m = alpha_sub.loc[alpha_sub["DistanceMetric"].astype(str).str.lower() == metric_name]
            if m.empty or "KrippendorffAlpha" not in m.columns:
                return np.nan
            v = pd.to_numeric(m["KrippendorffAlpha"].iloc[0], errors="coerce")
            return round(float(v), 3) if pd.notna(v) and np.isfinite(v) else np.nan

        pairable_units = np.nan
        if not alpha_sub.empty and "PairableUnits" in alpha_sub.columns:
            pu = pd.to_numeric(alpha_sub["PairableUnits"].iloc[0], errors="coerce")
            if pd.notna(pu) and np.isfinite(pu):
                pairable_units = int(pu)

        if not ia_row.empty:
            r0 = ia_row.iloc[0]
            exact_pct = round(100.0 * float(r0["ExactAgreement"]), 2) if pd.notna(r0.get("ExactAgreement")) else np.nan
            exp_pct = (
                np.nan
                if cal_hold or pd.isna(r0.get("ExpectedExactAgreement"))
                else round(100.0 * float(r0["ExpectedExactAgreement"]), 2)
            )
            adj_pct = round(100.0 * float(r0["AdjacentAgreement"]), 2) if pd.notna(r0.get("AdjacentAgreement")) else np.nan
            mad = round(float(r0["MeanAbsDiff"]), 3) if pd.notna(r0.get("MeanAbsDiff")) else np.nan
            corr = round(float(r0["MeanCorr"]), 3) if pd.notna(r0.get("MeanCorr")) else np.nan
            elements = int(r0["Elements"]) if pd.notna(r0.get("Elements")) else np.nan
            contexts = int(r0["Contexts"]) if pd.notna(r0.get("Contexts")) else np.nan
            total_pairs = int(r0["TotalPairs"]) if pd.notna(r0.get("TotalPairs")) else np.nan
        else:
            r_al = alpha_sub.iloc[0]
            exact_pct = exp_pct = adj_pct = mad = corr = contexts = total_pairs = np.nan
            elements = int(r_al["Elements"]) if pd.notna(r_al.get("Elements")) else np.nan

        rows.append({
            "Facet": str(facet),
            "Elements": elements,
            "SharedContexts": contexts,
            "ComparedPairs": total_pairs,
            "ExactPct": exact_pct,
            "ExpectedExactPct": exp_pct,
            "AdjacentPct": adj_pct,
            "MeanAbsDiff": mad,
            "MeanCorr": corr,
            "AlphaOrdinal": _alpha_val("ordinal"),
            "AlphaInterval": _alpha_val("interval"),
            "AlphaNominal": _alpha_val("nominal"),
            "PairableUnits": pairable_units,
        })

    out_tables: list[tuple[str, pd.DataFrame, str]] = []
    if rows:
        summary_df = _blank_na_columns(pd.DataFrame(rows))
        note = (
            "Observed exact and adjacent (+/-1 category) agreement percentages across shared rating contexts, "
            "model-expected exact agreement percentage under independence given fitted category probabilities, "
            "and Krippendorff's alpha (ordinal, interval, nominal). Agreement describes raw-score consistency "
            "across overlapping assignments; it does not adjust for severity differences or replace MFRM separation."
        )
        if cal_hold:
            note += " Model-expected exact agreement is withheld while calibration inference is on hold."
        out_tables.append(("Inter-rater Agreement and Krippendorff's Alpha", summary_df, note))

    if isinstance(ia_pairs, dict) and ia_pairs:
        from mfrm_app.design_assignment import infer_rater_facet

        cfg = result.get("config", {})
        roles = cfg.get("facet_roles", {})
        mapping = infer_rater_facet(
            facet_names,
            configured_rater_facet=(cfg.get("rater_facet") or (roles.get("rater") if isinstance(roles, dict) else None)),
        )
        primary = mapping.get("rater_facet")
        if not primary or primary not in ia_pairs or not isinstance(ia_pairs[primary], pd.DataFrame) or ia_pairs[primary].empty:
            primary = next(
                (f for f in facet_names if isinstance(ia_pairs.get(f), pd.DataFrame) and not ia_pairs[f].empty),
                None,
            )
        if primary is not None:
            raw_pairs = ia_pairs[primary].copy()
            pair_df = pd.DataFrame({
                "Element1": raw_pairs["Element1"].astype(str),
                "Element2": raw_pairs["Element2"].astype(str),
                "SharedContexts": pd.to_numeric(raw_pairs["N"], errors="coerce").astype("Int64"),
                "ExactPct": (100.0 * pd.to_numeric(raw_pairs["Exact"], errors="coerce")).round(2),
                "ExpectedExactPct": (
                    pd.Series(np.nan, index=raw_pairs.index, dtype=float)
                    if cal_hold
                    else (100.0 * pd.to_numeric(raw_pairs["ExpectedExact"], errors="coerce")).round(2)
                ),
                "AdjacentPct": (100.0 * pd.to_numeric(raw_pairs["Adjacent"], errors="coerce")).round(2),
                "MeanDiff": pd.to_numeric(raw_pairs["MeanDiff"], errors="coerce").round(3),
                "MAD": pd.to_numeric(raw_pairs["MAD"], errors="coerce").round(3),
                "Corr": pd.to_numeric(raw_pairs["Corr"], errors="coerce").round(3),
            })
            pair_note = (
                f"Pairwise observed and model-expected agreement across shared contexts for {primary}. "
                "Pairs without shared contexts are omitted."
            )
            if cal_hold:
                pair_note += " Model-expected exact agreement is withheld while calibration inference is on hold."
            out_tables.append((f"Pairwise Agreement: {primary}", _blank_na_columns(pair_df), pair_note))

    return out_tables


def _rescaled_measures_table(analysis) -> tuple[str, pd.DataFrame, str] | None:
    if "result" not in analysis:
        return None
    settings = analysis.get("fit_settings", {})
    umean = float(settings.get("umean", 0.0))
    uscale = float(settings.get("uscale", 1.0))
    if (umean, uscale) == (0.0, 1.0):
        return None
    result = analysis["result"]
    diagnostics = analysis.get("diagnostics", {})
    pub, _ = reporting_facts.publication_measure_table(result, diagnostics, fit_note="")
    if pub.empty:
        return None
    raw_measures = engine.qualify_measure_uncertainty(diagnostics.get("measures"), result)
    est_visible = pub["Estimate"].astype(str).ne("")
    se_visible = pub["SE / SD"].astype(str).ne("") if "SE / SD" in pub.columns else pd.Series(False, index=pub.index)
    lo_visible = pub["Lower"].astype(str).ne("") if "Lower" in pub.columns else pd.Series(False, index=pub.index)
    hi_visible = pub["Upper"].astype(str).ne("") if "Upper" in pub.columns else pd.Series(False, index=pub.index)

    raw_est = pd.to_numeric(raw_measures["Estimate"], errors="coerce")
    raw_se = pd.to_numeric(raw_measures.get("SE", pd.Series(np.nan, index=raw_measures.index)), errors="coerce")
    raw_lo = pd.to_numeric(raw_measures.get("CI_Lower", pd.Series(np.nan, index=raw_measures.index)), errors="coerce")
    raw_hi = pd.to_numeric(raw_measures.get("CI_Upper", pd.Series(np.nan, index=raw_measures.index)), errors="coerce")

    out = pd.DataFrame({
        "Facet": pub["Facet"].astype(str),
        "Level": pub["Level"].astype(str),
        "LogitEstimate": raw_est.where(est_visible, np.nan).round(3),
        "ScaledMeasure": (raw_est * uscale + umean).where(est_visible, np.nan).round(3),
        "ScaledSE": (raw_se * uscale).where(se_visible, np.nan).round(3),
        "ScaledLower": (raw_lo * uscale + umean).where(lo_visible, np.nan).round(3),
        "ScaledUpper": (raw_hi * uscale + umean).where(hi_visible, np.nan).round(3),
    })

    report_tables = analysis.get("report_tables", {})
    if isinstance(report_tables, dict) and report_tables:
        aux_frames = []
        for f_name, f_tbl in report_tables.items():
            if isinstance(f_tbl, pd.DataFrame) and "Level" in f_tbl.columns:
                sub = pd.DataFrame({
                    "Facet": str(f_name),
                    "Level": f_tbl["Level"].astype(str),
                    "FairM": pd.to_numeric(f_tbl.get("FairM"), errors="coerce").round(3),
                    "ObservedAverage": pd.to_numeric(f_tbl.get("ObservedAverage"), errors="coerce").round(3),
                    "Anchor": f_tbl.get("Anchor", pd.Series("", index=f_tbl.index)).fillna("").astype(str),
                })
                aux_frames.append(sub)
        if aux_frames:
            aux_all = pd.concat(aux_frames, ignore_index=True)
            out = out.merge(aux_all, on=["Facet", "Level"], how="left")

    out["Basis"] = pub["Basis"].astype(str).to_numpy()
    note = (
        f"Linear score rescaling (Umean = {umean:g}, Uscale = {uscale:g}): applies M' = Umean + Uscale * M and "
        "SE' = Uscale * SE to canonical logit estimates (separation and reliability ratios are affine-invariant). "
        "Baseline Measures and Fit, Wright maps, specified intervals and joint interactions remain in canonical logits. "
        "Qualification Basis and withheld cells match Measures and Fit."
    )
    return ("Rescaled Measures (Umean / Uscale)", _blank_na_columns(out), note)


def _pca_tables(analysis) -> list[tuple[str, pd.DataFrame, str]]:
    if "result" not in analysis:
        return []
    diagnostics = analysis.get("diagnostics", {})
    if not isinstance(diagnostics, dict) or not diagnostics.get("pca_enabled"):
        return []
    pca = diagnostics.get("pca")
    if not isinstance(pca, dict) or len(pca.get("eigenvalues", [])) == 0:
        unavail = pd.DataFrame([{
            "Status": "Unavailable",
            "Reason": "Fewer than 2 estimable facet-combination residual columns or insufficient residual variance.",
        }])
        return [(
            "Residual PCA: Eigenvalues and Stability",
            unavail,
            "Optional standardized-residual PCA was requested, but the retained rating design did not produce a valid residual correlation matrix.",
        )]

    ev = np.asarray(pca["eigenvalues"], dtype=float)
    vp = np.asarray(pca["variance_pct"], dtype=float)
    cum = np.cumsum(vp)
    stab = pca.get("stability", {}) if isinstance(pca.get("stability"), dict) else {}
    k = min(10, len(ev))
    ev_rows = []
    for idx in range(k):
        row = {
            "Component": f"PC{idx + 1}",
            "Eigenvalue": round(float(ev[idx]), 3),
            "VariancePct": round(float(vp[idx]), 2),
            "CumulativePct": round(float(cum[idx]), 2),
            "StabilityStatus": str(stab.get("PCAStabilityStatus", "")) if idx == 0 else "",
            "ResidualColumns": int(stab["ResidualColumns"]) if idx == 0 and pd.notna(stab.get("ResidualColumns")) else np.nan,
            "MissingSharePct": round(100.0 * float(stab["MissingShare"]), 1) if idx == 0 and pd.notna(stab.get("MissingShare")) else np.nan,
            "LOOEV1Range": (
                f"{float(stab['LOOEV1Min']):.2f}–{float(stab['LOOEV1Max']):.2f}"
                if idx == 0 and pd.notna(stab.get("LOOEV1Min")) and pd.notna(stab.get("LOOEV1Max"))
                else ""
            ),
            "BootstrapEV1CI90": (
                f"{float(stab['BootstrapEV1P05']):.2f}–{float(stab['BootstrapEV1P95']):.2f}"
                if idx == 0 and pd.notna(stab.get("BootstrapEV1P05")) and pd.notna(stab.get("BootstrapEV1P95"))
                else ""
            ),
            "Flags": str(stab.get("SensitivityFlagSummary", "")) if idx == 0 else "",
        }
        ev_rows.append(row)

    ev_note = (
        "Principal components of standardized observation residuals across facet-combination columns "
        "(sum of eigenvalues equals residual column count). First-contrast eigenvalues near or above 2.0 "
        "serve as an inspection prompt for secondary dimensions or local dependence, not an automatic pass/fail test; "
        "sparse cell overlap and mean-zero residual imputation can inflate contrast eigenvalues."
    )
    out = [("Residual PCA: Eigenvalues and Stability", _blank_na_columns(pd.DataFrame(ev_rows)), ev_note)]

    loadings = pca.get("loadings")
    if isinstance(loadings, pd.DataFrame) and not loadings.empty and "PC1" in loadings.columns:
        pc1 = pd.to_numeric(loadings["PC1"], errors="coerce")
        pc2 = pd.to_numeric(loadings["PC2"], errors="coerce") if "PC2" in loadings.columns else pd.Series(np.nan, index=loadings.index)
        pole = np.where(pc1 >= 0.40, "Positive (>= +0.40)", np.where(pc1 <= -0.40, "Negative (<= -0.40)", ""))
        load_df = pd.DataFrame({
            "Combination": loadings.index.astype(str),
            "PC1": pc1.round(3).to_numpy(),
            "PC2": pc2.round(3).to_numpy(),
            "AbsPC1": pc1.abs().round(3).to_numpy(),
            "ContrastPole": pole,
        }).sort_values("PC1", ascending=False, kind="mergesort").reset_index(drop=True)
        load_note = (
            "Standardized-residual loadings on the first two principal contrasts, ordered from the positive pole "
            "to the negative pole of PC1 to identify contrasting tasks, criteria, or raters."
        )
        out.append(("Residual PCA: First-Contrast Loadings", _blank_na_columns(load_df), load_note))

    return out


def preview_gd_design(data, *, object_facet, score_col, random_facets, fixed_condition=None):
    """Non-destructive pre-flight design check for the Desktop G/D-Study setup."""
    random_facets = list(random_facets or [])
    if not isinstance(data, pd.DataFrame) or data.empty:
        return dict(
            ok=False, status="missing_data", reason_code="missing_data", observed_n={}, extra_facets=[],
            summary_ja="データを読み込んでください。",
            summary_en="Load a rating dataset first.",
        )
    if not object_facet or not score_col or object_facet not in data.columns or score_col not in data.columns or object_facet == score_col:
        return dict(
            ok=False, status="incomplete_selection", reason_code="invalid_columns", observed_n={}, extra_facets=[],
            summary_ja="測定対象と観測得点にそれぞれ別の列を選択してください。",
            summary_en="Select distinct columns for the object of measurement and observed score.",
        )
    if len(random_facets) == 0:
        return dict(
            ok=False, status="incomplete_selection", reason_code="unsupported_facets", observed_n={}, extra_facets=[],
            summary_ja="ランダム相を1〜2列選択してください。",
            summary_en="Select 1 or 2 random-facet columns.",
        )
    if len(random_facets) > 2:
        return dict(
            ok=False, status="incomplete_selection", reason_code="unsupported_facets", observed_n={}, extra_facets=[],
            summary_ja=f"単変量G/D-Studyはランダム相1〜2列に対応しています（現在 {len(random_facets)} 列選択中）。残りの相は下の「固定相の扱い」で指定できます。",
            summary_en=f"Univariate G/D-Study supports 1–2 random facets ({len(random_facets)} selected). Handle a third facet under Fixed-facet handling below.",
        )
    try:
        validated_fc = _validate_fixed_condition(fixed_condition, object_facet, score_col, random_facets)
    except ValueError as exc:
        return dict(
            ok=False, status="incomplete_selection", reason_code="invalid_fixed_condition", observed_n={}, extra_facets=[],
            summary_ja=f"固定相の指定を確認してください：{exc}",
            summary_en=f"Check fixed-facet conditioning: {exc}",
        )
    pseudo_settings = dict(
        object_facet=object_facet, score_col=score_col, random_facets=random_facets,
        n_max_per_facet=10, target_population="", selected_plan=None,
    )
    if validated_fc is not None:
        pseudo_settings["fixed_condition"] = validated_fc
    prepared_df, prep_reason, fc_meta = prepare_gd_input_data(data, pseudo_settings)
    if prep_reason:
        reason_en = gd._GD_INPUT_REASONS.get(prep_reason, prep_reason)
        reason_ja_map = {
            "missing_cells": f"固定相 {validated_fc['facet']} を含む基本セルに欠測があるため集約できません（完全均衡交差が必要）。",
            "duplicate_cells": f"固定相 {validated_fc['facet']} を含めても同一セルに重複行があります。採点回等の列を確認してください。",
            "insufficient_levels": f"固定相 {validated_fc['facet']} または選択した相の水準数が不足しています（2水準以上必要）。",
            "invalid_scores": "観測得点に空欄または数値でない値が含まれています。",
            "invalid_labels": "相のラベルに空欄または表記揺れがあります。",
            "unsupported_weights": "Weight列はすべて1である必要があります。",
        }
        return dict(
            ok=False, status="incomplete", reason_code=prep_reason, observed_n={}, extra_facets=[],
            summary_ja="設計確認：" + reason_ja_map.get(prep_reason, reason_en),
            summary_en="Design check: " + reason_en,
        )

    work, design, reason = gd._prepare_g_d_data(prepared_df, object_facet, random_facets, score_col)
    if not reason and work is not None:
        levels = design["levels"]
        dims = " × ".join(f"{f} {levels[f]}" for f in [object_facet, *random_facets])
        n_cells = int(design["n_observations"])
        if fc_meta is None:
            ja = f"設計確認OK：{dims}（{n_cells:,}セル・欠測0・重複0）"
            en = f"Balanced crossed design ready: {dims} ({n_cells:,} cells, 0 missing, 0 duplicates)"
        elif fc_meta["mode"] == "level":
            ja = f"設計確認OK（{fc_meta['facet']} = {fc_meta['level']} を抽出）：{dims}（{n_cells:,}セル・欠測0）"
            en = f"Balanced slice ready ({fc_meta['facet']} = {fc_meta['level']}): {dims} ({n_cells:,} cells, 0 missing)"
        else:
            mode_ja = "合計" if fc_meta["mode"] == "sum" else "平均"
            ja = f"設計確認OK（{fc_meta['facet']} {fc_meta['n_fixed_levels']}水準の完全均衡を確認し{mode_ja}得点へ集約）：{dims}（{n_cells:,}セル）"
            en = f"Balanced {fc_meta['mode']} composite ready across {fc_meta['n_fixed_levels']} {fc_meta['facet']} levels: {dims} ({n_cells:,} cells)"
        return dict(
            ok=True, status="ready", reason_code="", observed_n=dict(design["observed_n"]),
            extra_facets=[], summary_ja=ja, summary_en=en,
        )

    if reason == "duplicate_cells" and validated_fc is None:
        used = {object_facet, score_col, *random_facets, "Weight"}
        extra = []
        for col in data.columns:
            if col not in used:
                n_lv = int(data[col].dropna().astype(str).str.strip().replace("", np.nan).dropna().nunique())
                if n_lv >= 2:
                    extra.append((col, n_lv))
        if extra:
            desc_ja = "、".join(f"{c}（{n}水準）" for c, n in extra)
            desc_en = ", ".join(f"{c} ({n} levels)" for c, n in extra)
            return dict(
                ok=False, status="duplicate_extra_facets", reason_code="duplicate_cells", observed_n={}, extra_facets=extra,
                summary_ja=f"同一セルに複数の評定があります（未選択の相：{desc_ja}）。下の「固定相の扱い」で1水準の抽出または完全均衡の合計／平均集約を指定してください。",
                summary_en=f"Multiple ratings per cell (unselected facet: {desc_en}). Choose single-level slice or verified balanced sum/mean composite under Fixed-facet handling below.",
            )

    expected = design.get("n_expected_cells", 0)
    obs = len(prepared_df) if isinstance(prepared_df, pd.DataFrame) else 0
    reason_en = gd._GD_INPUT_REASONS.get(reason, reason)
    reason_ja_map = {
        "duplicate_cells": "同一セルに複数の観測があります。反復評定を区別するか、固定相の扱いを指定してください。",
        "missing_cells": f"欠測セルがあるため完全均衡交差を満たしません（観測 {obs:,} / 期待 {expected:,} セル；欠測は0点とみなしません）。",
        "insufficient_levels": "測定対象および各ランダム相にはそれぞれ2水準以上が必要です。",
        "invalid_scores": "観測得点に空欄または数値でない値が含まれています。",
        "invalid_labels": "相のラベルに空欄または表記揺れがあります。",
        "unsupported_weights": "Weight列はすべて1である必要があります。",
    }
    return dict(
        ok=False, status="incomplete", reason_code=reason, observed_n={}, extra_facets=[],
        summary_ja="設計確認：" + reason_ja_map.get(reason, reason_en),
        summary_en="Design check: " + reason_en,
    )


def _run_gd(data, settings):
    arguments = {k: v for k, v in settings.items() if k != "selected_plan"}
    result = gd.compute_generalizability_from_data(data, **arguments)
    if result["available"]:
        counts = settings["selected_plan"] or result["design"]["observed_n"]
        coefficients = gd._g_phi_from_variance_components(
            result["variance_components"],
            object_facet=settings["object_facet"],
            facet_observations=counts,
        )
        plans = result["d_study"]
        in_grid = np.logical_and.reduce([plans[f].eq(n) for f, n in counts.items()]).any()
        if not in_grid:
            extra_row = pd.DataFrame([{
                **counts,
                "G": float(coefficients["G"]),
                "Phi": float(coefficients["Phi"]),
                "IsObservedDesign": bool(counts == result["design"]["observed_n"]),
            }])
            result["d_study"] = (
                pd.concat([plans, extra_row], ignore_index=True)
                .sort_values(list(settings["random_facets"]), kind="mergesort")
                .reset_index(drop=True)
            )
        result["selected_plan"] = {**counts, **gd._g_d_error_row(coefficients)}
    return dict(settings=settings, result=result,
                analysis_id=project.gd_analysis_identity(data, settings, result["calculation_version"]))


def run_analysis(data, *, workflow="mfrm", bias_requests=None, interval_requests=None,
                 gd_settings=None, joint_settings=None, progress=None, **fit_settings):
    """Compute once in a worker and archive settings, failure states and provenance.

    An unsupported G/D design is a retained unavailable result. Interval review
    failures do not discard a fitted calibration or silently refit it.
    """
    progress = progress or (lambda phase: None)
    if not isinstance(data, pd.DataFrame):
        raise ValueError("An input data frame is required")
    if workflow not in {"mfrm", "gd"}:
        raise ValueError("Choose the MFRM or observed-score G/D workflow")
    gd_settings = _gd_settings(gd_settings) if gd_settings is not None else None
    if workflow == "gd":
        if fit_settings or bias_requests or interval_requests or joint_settings is not None or gd_settings is None:
            raise ValueError("Standalone G/D requires G/D settings, without MFRM requests")
        progress("generalizability")
        entry = _run_gd(data, gd_settings)
        analysis = dict(kind="gd-analysis-v2", analysis_id=entry["analysis_id"],
            created_at=datetime.now(timezone.utc).isoformat(), source=basic._source_versions(),
            input_data=data.copy(deep=True), gd_analysis=entry,
            scope=dict(includes_input=True, optimizer_checkpoint=False, mfrm_fitted=False))
    else:
        facets = list(fit_settings.get("facet_cols", []))
        biases = _bias_requests(bias_requests, ["Person", *facets])
        requested_intervals = _interval_requests(interval_requests, facets, fit_settings.get("method", "JMLE"))
        joint_settings = joint.settings(joint_settings, facets, fit_settings.get("method", "JMLE"))
        if joint_settings is not None and fit_settings.get('population_prior_sd', 1.) != 1.:
            raise ValueError('Joint interactions require a fixed N(0,1) population')
        analysis = basic.run_basic_analysis(data, **fit_settings, progress=progress)
        if analysis["kind"] != "mfrm-analysis-v2" and (
            biases or requested_intervals or gd_settings is not None or joint_settings is not None
            or fit_settings.get("anchors") or fit_settings.get("group_anchors")
        ):
            analysis["kind"] = "mfrm-analysis-v2"
        analysis["bias_analyses"], analysis["interval_analyses"] = [], []
        for settings in biases:
            progress("bias")
            output = bias.estimate_bias_interaction(analysis["result"], analysis["diagnostics"], **settings)
            analysis["bias_analyses"].append(dict(analysis_id=analysis["analysis_id"], settings=settings, result=output))
        analysis["scope"]["bias_refits"] = bool(biases)
        for settings in requested_intervals:
            progress("intervals")
            try:
                record = intervals.compute_structural_intervals(analysis["result"], **settings)
                entry = dict(record=record, status="computed", reason="")
            except ValueError as exc:
                entry = dict(record=None, status="unavailable", reason=str(exc))
            except Exception as exc:
                entry = dict(record=None, status="failed", reason=f"{type(exc).__name__}: {exc}")
            analysis["interval_analyses"].append(dict(analysis_id=analysis["analysis_id"], settings=settings, **entry))
        if joint_settings is not None:
            progress('joint_interaction')
            try:
                record = joint.fit_joint_interaction(analysis['result'], **joint_settings)
                entry = dict(record=record, status='computed', reason='')
            except ValueError as exc:
                entry = dict(record=None, status='unavailable', reason=str(exc))
            except Exception as exc:
                entry = dict(record=None, status='failed', reason=f'{type(exc).__name__}: {exc}')
            analysis['joint_analysis'] = dict(analysis_id=analysis['analysis_id'], settings=joint_settings, **entry)
        if gd_settings is not None:
            progress("generalizability")
            analysis["gd_analysis"] = _run_gd(data, gd_settings)
    analysis["content_fingerprint"] = project.content_fingerprint(analysis)
    project.validate_analysis(analysis)
    return analysis


def current_additional_tables(analysis):
    """Reapply output guards to retained records, without estimation or ANOVA."""
    project.validate_analysis(analysis)
    tables = []
    for entry in analysis.get("bias_analyses", []):
        settings, result = entry["settings"], entry["result"]
        pair = settings["facet_a"] + " × " + settings["facet_b"]
        table = bias._bias_table_for_publication({pair: result}, analysis["result"]["config"]["model"])
        if table is None:
            table = pd.DataFrame([dict(Use="Unavailable", Reason=result.get("_skip_reason", "No computed cells"))])
        tables.append(("Bias: " + pair, table, bias._bias_table_caption(analysis["result"]["config"]["model"])))
    for entry in analysis.get("interval_analyses", []):
        try:
            if entry["status"] != "computed":
                raise ValueError(entry["reason"])
            table, note = intervals.publication_table(analysis["result"], entry["record"])
        except ValueError as exc:
            table, note = pd.DataFrame([dict(Use="Withheld", Reason=str(exc))]), intervals.CAUTION
        tables.append(("Specified intervals: " + entry["settings"]["facet"], table, note))
    if entry := analysis.get("gd_analysis"):
        frame = gd.g_d_study_export_frame(entry["result"])
        by_section = {
            section: table.drop(columns="_section").dropna(axis=1, how="all")
            for section, table in frame.groupby("_section", sort=False)
        }
        ordered_sections = [s for s in gd.GD_SECTION_ORDER if s in by_section]
        ordered_sections += [s for s in by_section if s not in gd.GD_SECTION_ORDER]
        for section in ordered_sections:
            sec_table = by_section[section]
            note = gd.GD_SECTION_NOTES.get(section, gd._GD_CI_REASON)
            if (
                section == "variance_components"
                and "RawVariance" in sec_table.columns
                and (pd.to_numeric(sec_table["RawVariance"], errors="coerce") < -1e-12).any()
            ):
                note = gd.GD_VARIANCE_NEGATIVE_NOTE
            tables.append(("G/D: " + section, sec_table, note))
    if entry := analysis.get('joint_analysis'):
        try:
            if entry['status'] != 'computed':
                raise ValueError(entry['reason'])
            table, note = joint.publication_table(analysis['result'], entry['record'])
        except ValueError as exc:
            table, note = pd.DataFrame([dict(Use='Withheld', Reason=str(exc))]), joint.CAUTION
        tables.append(('Joint interaction: '+entry['settings']['facet_a']+' × '+entry['settings']['facet_b'], table, note))
        if entry['status'] == 'computed' and 'CIEligible' in table:
            tables.append(('Joint model: main effects and steps',
                pd.DataFrame(entry['record']['main_effects_and_steps']),
                'Point estimates from the separate joint model; baseline maps and tables retain the original main-effects fit.'))
    if linking := _connectivity_and_anchor_table(analysis):
        tables.append(linking)
    if rel_tbl := _reliability_summary_table(analysis):
        tables.append(rel_tbl)
    tables.extend(_agreement_and_alpha_tables(analysis))
    if rescaled := _rescaled_measures_table(analysis):
        tables.append(rescaled)
    tables.extend(_pca_tables(analysis))
    return tables
