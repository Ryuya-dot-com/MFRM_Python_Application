"""UI-independent calibration and the original basic desktop reports.

RSM/PCM with JMLE or fixed-SD MML; the basic native figures still serve RSM.
Supplementary calculations and project requests live in desktop_analysis.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import html
import importlib.metadata
import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from mfrm_app import engine, output_qualification, project, reporting_facts
from mfrm_app._resources import metric_note


def qualification_snapshot(result):
    state = output_qualification.run_qualification_state(result)
    return {
        "schema": output_qualification.OUTPUT_QUALIFICATION_SCHEMA_VERSION,
        "optimizer_terminated": output_qualification.optimizer_termination_established(result),
        "calibration_hold_reason": output_qualification.calibration_inference_hold_reason(result),
        "calibration_hold_message": output_qualification.calibration_inference_hold_message(result),
        "screening_block_reason": state["screening_block_reason"],
        "identifiability_state": state["identifiability_state"],
        "run_state": state["state"],
        "supports_formal_inference": False,
        "facet_equivalence": output_qualification.facet_equivalence_qualification().to_record(),
        "bias_cell": output_qualification.bias_cell_qualification().to_record(),
        "bias_pairwise": output_qualification.bias_pairwise_qualification().to_record(),
    }


def _source_versions():
    build_record = Path(__file__).with_name("desktop_build.json")
    if build_record.exists():
        source_hashes = json.loads(build_record.read_text(encoding="utf-8"))["source_sha256"]
    else:
        modules = (engine, output_qualification, reporting_facts, project)
        paths = [Path(module.__file__) for module in modules] + [Path(__file__)]
        for name in ("desktop_analysis", "desktop_examples", "desktop_outputs", "bias_analysis", "generalizability",
                     "reliability", "structural_intervals", "joint_interaction", "project_portability",
                     "mml_fixed_sd_review", "mml_joint"):
            candidate = Path(__file__).with_name(name + ".py")
            if candidate.exists():
                paths.append(candidate)
        source_hashes = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    return {
        "app_version": engine.APP_VERSION,
        "source_sha256": source_hashes,
        "packages": {name: importlib.metadata.version(name) for name in ("numpy", "pandas", "scipy")},
    }


def normalize_anchor_records(records, *, person_col: str, facet_cols: list[str], observed_levels: dict[str, set[str]] | None = None) -> list[dict]:
    """Validate and canonicalize fixed-logit anchor dicts (`Facet`, `Level`, `Anchor`)."""
    if not records:
        return []
    if not isinstance(records, (list, tuple)):
        raise ValueError("Anchors must be provided as a list of {Facet, Level, Anchor} mappings")
    allowed = ["Person", *facet_cols]
    out, seen = [], set()
    for item in records:
        if not isinstance(item, dict) or set(item) - {"Facet", "Level", "Anchor"} or not {"Facet", "Level", "Anchor"}.issubset(item):
            raise ValueError("Each fixed anchor requires Facet, Level and Anchor")
        raw_facet, level, val = item["Facet"], item["Level"], item["Anchor"]
        if not isinstance(raw_facet, str) or not isinstance(level, str) or not raw_facet.strip() or not level.strip():
            raise ValueError("Anchor Facet and Level must be nonblank strings")
        facet = "Person" if raw_facet.strip() in {"Person", person_col} else raw_facet.strip()
        level = level.strip()
        if facet not in allowed:
            raise ValueError(f"Anchor facet {raw_facet!r} is not in the selected model facets")
        if observed_levels is not None and level not in observed_levels.get(facet, set()):
            raise ValueError(f"Anchor level {level!r} was not observed in facet {facet!r}")
        if isinstance(val, (bool, np.bool_)):
            raise ValueError("Anchor values must be finite numbers")
        try:
            fval = float(val)
        except (TypeError, ValueError) as exc:
            raise ValueError("Anchor values must be finite numbers") from exc
        if not np.isfinite(fval):
            raise ValueError("Anchor values must be finite numbers")
        if (facet, level) in seen:
            raise ValueError(f"Duplicate anchor specified for {facet}: {level}")
        seen.add((facet, level))
        out.append({"Facet": facet, "Level": level, "Anchor": fval})
    return sorted(out, key=lambda r: (allowed.index(r["Facet"]), r["Level"]))


def normalize_group_anchor_records(
    records,
    *,
    person_col: str,
    facet_cols: list[str],
    observed_levels: dict[str, set[str]] | None = None,
    fixed_pairs: set[tuple[str, str]] | None = None,
    fixed_records: list[dict] | None = None,
) -> list[dict]:
    """Validate and canonicalize group-mean anchor dicts (`Facet`, `Level`, `Group`, `GroupValue`)."""
    if not records:
        return []
    if not isinstance(records, (list, tuple)):
        raise ValueError("Group anchors must be provided as a list of {Facet, Level, Group, GroupValue} mappings")
    allowed = ["Person", *facet_cols]
    fixed_pairs = set(fixed_pairs or set())
    if fixed_records:
        fixed_pairs |= {(r["Facet"], r["Level"]) for r in fixed_records}
    out, seen, group_means = [], set(), {}
    for item in records:
        req = {"Facet", "Level", "Group", "GroupValue"}
        if not isinstance(item, dict) or set(item) - req or not req.issubset(item):
            raise ValueError("Each group anchor requires Facet, Level, Group and GroupValue")
        raw_facet, level, group, val = item["Facet"], item["Level"], item["Group"], item["GroupValue"]
        if (not isinstance(raw_facet, str) or not isinstance(level, str) or not isinstance(group, str)
                or not raw_facet.strip() or not level.strip() or not group.strip()):
            raise ValueError("Group-anchor Facet, Level and Group must be nonblank strings")
        facet = "Person" if raw_facet.strip() in {"Person", person_col} else raw_facet.strip()
        level, group = level.strip(), group.strip()
        if facet not in allowed:
            raise ValueError(f"Group-anchor facet {raw_facet!r} is not in the selected model facets")
        if observed_levels is not None and level not in observed_levels.get(facet, set()):
            raise ValueError(f"Group-anchor level {level!r} was not observed in facet {facet!r}")
        if (facet, level) in fixed_pairs:
            raise ValueError(f"{facet}: {level} cannot have both a fixed anchor and a group anchor")
        if (facet, level) in seen:
            raise ValueError(f"Duplicate group anchor specified for {facet}: {level}")
        if isinstance(val, (bool, np.bool_)):
            raise ValueError("Group-anchor values must be finite numbers")
        try:
            fval = float(val)
        except (TypeError, ValueError) as exc:
            raise ValueError("Group-anchor values must be finite numbers") from exc
        if not np.isfinite(fval):
            raise ValueError("Group-anchor values must be finite numbers")
        prev = group_means.get((facet, group))
        if prev is not None and abs(prev - fval) > 1e-12:
            raise ValueError(f"Group {group!r} in facet {facet!r} has conflicting GroupValue numbers ({prev:g} vs {fval:g})")
        group_means[(facet, group)] = fval
        seen.add((facet, level))
        out.append({"Facet": facet, "Level": level, "Group": group, "GroupValue": fval})
    counts_by_group = {}
    for row in out:
        key = (row["Facet"], row["Group"])
        counts_by_group[key] = counts_by_group.get(key, 0) + 1
    for (facet, group), cnt in counts_by_group.items():
        if cnt < 2:
            raise ValueError(f"Group anchor {group!r} in facet {facet!r} needs at least 2 levels (use a fixed anchor for 1 level)")
    return sorted(out, key=lambda r: (allowed.index(r["Facet"]), r["Group"], r["Level"]))


def normalize_scale_orientation(
    noncenter_facet="Person",
    positive_facets=None,
    *,
    person_col: str,
    facet_cols: list[str],
) -> tuple[str, list[str]]:
    """Validate and canonicalize `noncenter_facet` and `positive_facets` against model facets."""
    allowed = ["Person", *facet_cols]
    if noncenter_facet is None:
        norm_nc = "Person"
    else:
        if not isinstance(noncenter_facet, str) or not noncenter_facet.strip():
            raise ValueError("noncenter_facet must be 'Person' or a selected facet column name")
        raw_nc = noncenter_facet.strip()
        norm_nc = "Person" if raw_nc in {"Person", person_col} else raw_nc
        if norm_nc not in allowed:
            raise ValueError(f"noncenter_facet {noncenter_facet!r} is not in the selected model facets")

    if positive_facets is None:
        return norm_nc, []
    if isinstance(positive_facets, (str, bytes)) or not isinstance(positive_facets, (list, tuple, set, frozenset)):
        raise ValueError("positive_facets must be a list or tuple of selected facet names")
    pos_set: set[str] = set()
    for item in positive_facets:
        if not isinstance(item, str) or not item.strip():
            raise ValueError("Each positive_facets entry must be a nonblank facet name")
        raw_f = item.strip()
        mapped = "Person" if raw_f in {"Person", person_col} else raw_f
        if mapped not in allowed:
            raise ValueError(f"positive_facets entry {item!r} is not in the selected model facets")
        if mapped != "Person":
            pos_set.add(mapped)
    return norm_nc, [f for f in facet_cols if f in pos_set]


def _finish_analysis(result, data, settings, progress):
    compute_pca = bool(settings.get("compute_pca", False))
    umean = float(settings.get("umean", 0.0))
    uscale = float(settings.get("uscale", 1.0))
    has_rescaling = (umean, uscale) != (0.0, 1.0)
    has_scale_override = (
        settings.get("noncenter_facet", "Person") != "Person"
        or bool(settings.get("positive_facets"))
    )
    progress("diagnostics")
    diagnostics = engine.mfrm_diagnostics(result, compute_pca=compute_pca, compute_marginal=False,
                                         compute_eb_shrinkage=False)
    progress("tables")
    tables = engine.calc_facets_report_tbls(result, diagnostics, umean=umean, uscale=uscale)
    has_anchors = bool(settings.get("anchors") or settings.get("group_anchors"))
    analysis = {
        "kind": ("basic-rsm-jmle-v1"
                 if settings["model"] == "RSM" and settings["method"] == "JMLE"
                 and not has_anchors and not compute_pca and not has_rescaling and not has_scale_override
                 else "mfrm-analysis-v2"),
        "analysis_id": engine.build_result_analysis_identity(result).analysis_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source": _source_versions(),
        "input_data": data.copy(deep=True),
        "fit_settings": dict(settings),
        "result": result,
        "diagnostics": diagnostics,
        "report_tables": tables,
        "scorefile": engine.compute_scorefile(result),
        "qualification_at_save": qualification_snapshot(result),
        "scope": {"includes_input": True, "optimizer_checkpoint": False,
                  "pca": compute_pca, "marginal_fit": False, "bias_refits": False},
    }
    analysis["content_fingerprint"] = project.content_fingerprint(analysis)
    return analysis


def _is_int(value) -> bool:
    return isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_))


def run_basic_analysis(data, *, person_col, facet_cols, score_col, rating_min, rating_max,
                       maxit=400, reltol=1e-6, model="RSM", method="JMLE", step_facet=None,
                       quad_points=None, population_prior_sd=1.0, keep_original=True,
                       estimate_population_sd=False, mml_engine="EM",
                       noncenter_facet="Person", positive_facets=None,
                       anchors=None, group_anchors=None,
                       compute_pca=False, umean=0.0, uscale=1.0, progress=None):
    """Fit once and retain qualified tables. Failure/unknown termination stays so.

    Progress identifies phases, not a fabricated percentage. A desktop worker
    owns cancellation; the estimator is never run in the Tk event loop.
    """
    if not isinstance(data, pd.DataFrame) or data.empty:
        raise ValueError("A nonempty long-form ratings table is required")
    facet_cols = list(facet_cols)
    columns = [person_col, *facet_cols, score_col]
    if not facet_cols or len(set(columns)) != len(columns) or not data.columns.is_unique:
        raise ValueError("Select distinct person, facet and score columns")
    if any(column not in data for column in columns):
        raise ValueError("A selected column is missing from the input")
    if not _is_int(rating_min) or not _is_int(rating_max) or int(rating_min) >= int(rating_max):
        raise ValueError("Specify integer minimum/maximum scores with minimum < maximum")
    rating_min, rating_max = int(rating_min), int(rating_max)
    if not _is_int(maxit) or int(maxit) < 1 or isinstance(reltol, (bool, np.bool_)) or not np.isfinite(reltol) or float(reltol) <= 0:
        raise ValueError("Invalid iteration limit or stopping tolerance")
    maxit, reltol = int(maxit), float(reltol)
    if model not in {"RSM", "PCM"} or method not in {"JMLE", "MML"}:
        raise ValueError("Desktop analysis supports RSM/PCM with JMLE or fixed-SD MML")
    if not isinstance(keep_original, (bool, np.bool_)) or not bool(keep_original) or not isinstance(estimate_population_sd, (bool, np.bool_)) or bool(estimate_population_sd) or mml_engine != "EM":
        raise ValueError("Desktop requests retain original scores and use fixed-SD EM for MML")
    keep_original, estimate_population_sd = True, False
    if model == "PCM" and step_facet not in facet_cols:
        raise ValueError("PCM requires an explicitly selected step facet")
    if model == "RSM" and step_facet is not None:
        raise ValueError("RSM uses common steps; do not select a PCM step facet")
    if isinstance(population_prior_sd, (bool, np.bool_)) or not np.isfinite(population_prior_sd) or float(population_prior_sd) <= 0:
        raise ValueError("MML requires at least three quadrature points and a positive finite fixed SD")
    if method == "JMLE":
        if (quad_points is not None and (not _is_int(quad_points) or int(quad_points) != 15)) or float(population_prior_sd) != 1.0:
            raise ValueError("Population SD and quadrature settings apply only to MML")
        quad_points = 15
    else:
        resolved_qp = 31 if quad_points is None else quad_points
        if not _is_int(resolved_qp) or int(resolved_qp) < 3:
            raise ValueError("MML requires at least three quadrature points and a positive finite fixed SD")
        quad_points = int(resolved_qp)
    norm_noncenter, norm_positive = normalize_scale_orientation(
        noncenter_facet=noncenter_facet,
        positive_facets=positive_facets,
        person_col=person_col,
        facet_cols=facet_cols,
    )
    if not isinstance(compute_pca, (bool, np.bool_)):
        raise ValueError("Residual PCA option must be True or False")
    compute_pca = bool(compute_pca)
    if isinstance(umean, (bool, np.bool_)) or isinstance(uscale, (bool, np.bool_)):
        raise ValueError("Scale origin (umean) and multiplier (uscale) must be finite numbers with uscale > 0")
    try:
        f_umean = float(umean)
        f_uscale = float(uscale)
    except (TypeError, ValueError) as exc:
        raise ValueError("Scale origin (umean) and multiplier (uscale) must be finite numbers with uscale > 0") from exc
    if not np.isfinite(f_umean) or not np.isfinite(f_uscale) or f_uscale <= 0:
        raise ValueError("Scale origin (umean) and multiplier (uscale) must be finite numbers with uscale > 0")
    if f_umean == 0.0:
        f_umean = 0.0

    audit_rows = engine.build_response_data_audit(data, person_col=person_col, facet_cols=facet_cols, score_col=score_col)
    retained = data.loc[audit_rows["rows"]["IncludedInLikelihood"].to_numpy(dtype=bool)]
    observed_levels = {
        "Person": set(retained[person_col].dropna().astype(str).str.strip()) if not retained.empty else set(),
        **{
            f: set(retained[f].dropna().astype(str).str.strip()) if not retained.empty else set()
            for f in facet_cols
        },
    }
    norm_anchors = normalize_anchor_records(
        anchors, person_col=person_col, facet_cols=facet_cols, observed_levels=observed_levels
    )
    norm_groups = normalize_group_anchor_records(
        group_anchors,
        person_col=person_col,
        facet_cols=facet_cols,
        observed_levels=observed_levels,
        fixed_pairs={(r["Facet"], r["Level"]) for r in norm_anchors},
    )

    settings = dict(person_col=person_col, facet_cols=facet_cols, score_col=score_col,
                    rating_min=rating_min, rating_max=rating_max, maxit=maxit, reltol=reltol,
                    model=model, method=method, keep_original=True)
    if model == "PCM":
        settings["step_facet"] = step_facet
    if method == "MML":
        settings.update(quad_points=quad_points, population_prior_sd=float(population_prior_sd),
                        estimate_population_sd=False, mml_engine="EM")
    if norm_noncenter != "Person":
        settings["noncenter_facet"] = norm_noncenter
    if norm_positive:
        settings["positive_facets"] = norm_positive
    if norm_anchors:
        settings["anchors"] = norm_anchors
    if norm_groups:
        settings["group_anchors"] = norm_groups
    if compute_pca:
        settings["compute_pca"] = True
    if (f_umean, f_uscale) != (0.0, 1.0):
        settings["umean"] = f_umean
        settings["uscale"] = f_uscale

    estimate_kwargs = {
        k: v for k, v in settings.items()
        if k not in {"anchors", "group_anchors", "compute_pca", "umean", "uscale"}
    }
    if norm_anchors:
        estimate_kwargs["anchor_df"] = pd.DataFrame(norm_anchors)
    if norm_groups:
        estimate_kwargs["group_anchor_df"] = pd.DataFrame(norm_groups)

    progress = progress or (lambda phase: None)
    progress("estimating")
    result = engine.mfrm_estimate(data=data, **estimate_kwargs)
    return _finish_analysis(result, data, settings, progress)


def current_tables(analysis, *, extended: bool = False):
    """Reapply current display qualification without refitting or rewriting provenance."""
    project.validate_analysis(analysis)
    result = analysis["result"]
    facts = reporting_facts.observed_report_tables(result)
    tables = [(title, table, note) for title, table, note in facts]
    fit_note = metric_note("fit")
    measures, note = reporting_facts.publication_measure_table(
        result, analysis["diagnostics"], fit_note=fit_note, extended=extended
    )
    note = note.replace("Method details and source values: measures CSV.",
                        "Method details and source values are retained in this project.")
    tables.insert(1, ("Measures and Fit", measures, note))
    return tables


def category_figure(analysis):
    """Representative RSM probability curves at the net predictor, not raw scores."""
    from mfrm_app.visual_maps import category_probability_figure

    result = analysis["result"]
    if result["config"]["model"] != "RSM":
        raise ValueError("PCM category curves require a selected step-facet level")
    steps = np.asarray(result["params"]["steps"], dtype=float)
    if not np.isfinite(steps).all():
        raise ValueError("Category curves require finite step estimates")
    return category_probability_figure(result, steps, title="Rating scale · category probabilities", figsize=(7, 4))


def basic_figure(analysis, kind, *, facet=None, statistic="Infit"):
    from mfrm_app.visual_maps import wright_figure, pathway_figure
    if kind == "wright":
        return wright_figure(analysis["result"])
    if kind == "pathway":
        return pathway_figure(analysis["result"], analysis["diagnostics"], facet=facet, statistic=statistic)
    if kind == "category":
        return category_figure(analysis)
    raise ValueError(f"Unknown figure: {kind}")


def figure_bytes(analysis, kind, *, facet=None, statistic="Infit", format="svg", figure=None):
    """Export the displayed coordinates, without estimation or diagnostics."""
    from matplotlib import rc_context
    if format not in {"png", "svg"}:
        raise ValueError("Choose PNG or SVG")
    figure = figure if figure is not None else basic_figure(analysis, kind, facet=facet, statistic=statistic)
    output = io.BytesIO()
    with rc_context({"svg.hashsalt": analysis["analysis_id"] + kind}):
        figure.savefig(output, format=format, dpi=300,
                       metadata={"Date": None} if format == "svg" else None)
    content = output.getvalue()
    if format == "svg":
        # Fixed ids such as axes_1 also collide when three inline SVGs share a page.
        svg = content.decode("utf-8")
        ids = re.findall(r'\bid="([^"]+)"', svg)
        for old in sorted(set(ids), key=len, reverse=True):
            new = kind + "_" + old
            svg = svg.replace(f'id="{old}"', f'id="{new}"')
            svg = svg.replace(f'href="#{old}"', f'href="#{new}"')
            svg = svg.replace(f'url(#{old})', f'url(#{new})')
        content = svg.encode("utf-8")
    return content


def basic_report_html(analysis, *, pathway_statistic="Infit", pathway_facet=None):
    """Standalone HTML with qualified shared tables and three offline SVGs."""
    if analysis["kind"] != "basic-rsm-jmle-v1":
        raise ValueError("The basic HTML report supports the original RSM/JMLE workflow only")
    tables = current_tables(analysis)
    result = analysis["result"]
    current = qualification_snapshot(result)
    original = analysis["qualification_at_save"]
    escape = html.escape
    paragraphs = ["RSM / JMLE. Input data are retained in the project file.",
                  "Residual PCA, marginal diagnostics and bias refits were not run in this workflow.",
                  reporting_facts.optimizer_termination_text(result),
                  reporting_facts.recorded_method_details(result),
                  current["calibration_hold_message"],
                  "Automatic facet equivalence conclusions remain withheld."]
    if current != original:
        paragraphs.append("The conditions for interpreting this result have changed since it was saved. "
                          "The current report uses the updated conditions; the original record is retained in the project file.")
    parts = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>MFRM basic analysis</title>',
             '<style>body{font:16px sans-serif;margin:2em;max-width:1200px}table{border-collapse:collapse;font-size:12px}',
             'th,td{padding:.4em;border-bottom:1px solid #bbb}section{overflow-x:auto}svg{max-width:100%;height:auto}</style>',
             '<h1>MFRM basic analysis</h1>', '<p>Analysis ID: ' + escape(analysis["analysis_id"]) + '</p>']
    parts.extend('<p>' + escape(text) + '</p>' for text in paragraphs if text)
    parts.append('<h2>Estimation summary</h2>' + engine.qualified_estimation_summary(result).to_html(index=False, escape=True))
    for title, table, note in tables:
        parts.extend(['<h2>' + escape(title) + '</h2><section>', table.to_html(index=False, escape=True),
                      '</section><p>' + escape(note) + '</p>'])
    from mfrm_app.visual_maps import MAP_GUIDES
    from mfrm_app.desktop_guidance import WEB_APP_URL
    for kind, (title, guide) in MAP_GUIDES.items():
        options = {"statistic": pathway_statistic, "facet": pathway_facet} if kind == "pathway" else {}
        svg = figure_bytes(analysis, kind, **options).decode("utf-8")
        parts.extend(['<h2>' + escape(title) + '</h2>', '<p lang="ja">' + escape(guide) + '</p>', svg[svg.index('<svg'):]])
    parts.append('<p><a href="' + escape(WEB_APP_URL, quote=True) + '">Web版の理論・ヘルプ</a>：'
                 '入口からHelp／Tutorialを選択してください。リンクに入力データ・結果は含まれません。</p>')
    parts.append('<p>These curves describe the fitted rating scale, conditional on its assumptions; '
                 'they do not certify model fit or uncertainty coverage.</p></html>')
    return '\n'.join(parts)
