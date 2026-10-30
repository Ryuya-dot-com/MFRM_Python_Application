"""Desktop calculation requests and retained results, independent of either UI."""
from copy import deepcopy
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from mfrm_app import basic_analysis as basic, bias_analysis as bias
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
        if biases or requested_intervals or gd_settings is not None or joint_settings is not None:
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
    return tables
