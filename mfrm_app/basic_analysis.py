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


def qualification_snapshot(result):
    return {
        "schema": output_qualification.OUTPUT_QUALIFICATION_SCHEMA_VERSION,
        "optimizer_terminated": output_qualification.optimizer_termination_established(result),
        "calibration_hold_reason": output_qualification.calibration_inference_hold_reason(result),
        "calibration_hold_message": output_qualification.calibration_inference_hold_message(result),
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
        paths += [Path(__file__).with_name(name + ".py") for name in
                  ("desktop_analysis", "desktop_outputs", "bias_analysis", "generalizability", "structural_intervals",
                   "mml_fixed_sd_review", "mml_joint")]
        source_hashes = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    return {
        "app_version": engine.APP_VERSION,
        "source_sha256": source_hashes,
        "packages": {name: importlib.metadata.version(name) for name in ("numpy", "pandas", "scipy")},
    }


def _finish_analysis(result, data, settings, progress):
    progress("diagnostics")
    diagnostics = engine.mfrm_diagnostics(result, compute_pca=False, compute_marginal=False,
                                         compute_eb_shrinkage=False)
    progress("tables")
    tables = engine.calc_facets_report_tbls(result, diagnostics)
    analysis = {
        "kind": ("basic-rsm-jmle-v1" if settings["model"] == "RSM" and settings["method"] == "JMLE"
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
                  "pca": False, "marginal_fit": False, "bias_refits": False},
    }
    analysis["content_fingerprint"] = project.content_fingerprint(analysis)
    return analysis


def run_basic_analysis(data, *, person_col, facet_cols, score_col, rating_min, rating_max,
                       maxit=400, reltol=1e-6, model="RSM", method="JMLE", step_facet=None,
                       quad_points=15, population_prior_sd=1.0, keep_original=True,
                       estimate_population_sd=False, mml_engine="EM", progress=None):
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
    if (type(rating_min) is not int or type(rating_max) is not int or rating_min >= rating_max):
        raise ValueError("Specify integer minimum/maximum scores with minimum < maximum")
    if type(maxit) is not int or maxit < 1 or not np.isfinite(reltol) or reltol <= 0:
        raise ValueError("Invalid iteration limit or stopping tolerance")
    if model not in {"RSM", "PCM"} or method not in {"JMLE", "MML"}:
        raise ValueError("Desktop analysis supports RSM/PCM with JMLE or fixed-SD MML")
    if keep_original is not True or estimate_population_sd is not False or mml_engine != "EM":
        raise ValueError("Desktop requests retain original scores and use fixed-SD EM for MML")
    if model == "PCM" and step_facet not in facet_cols:
        raise ValueError("PCM requires an explicitly selected step facet")
    if model == "RSM" and step_facet is not None:
        raise ValueError("RSM uses common steps; do not select a PCM step facet")
    if (type(quad_points) is not int or quad_points < 3
            or isinstance(population_prior_sd, (bool, np.bool_))
            or not np.isfinite(population_prior_sd) or population_prior_sd <= 0):
        raise ValueError("MML requires at least three quadrature points and a positive finite fixed SD")
    if method == "JMLE" and (quad_points != 15 or population_prior_sd != 1.0):
        raise ValueError("Population SD and quadrature settings apply only to MML")
    settings = dict(person_col=person_col, facet_cols=facet_cols, score_col=score_col,
                    rating_min=rating_min, rating_max=rating_max, maxit=maxit, reltol=reltol,
                    model=model, method=method, keep_original=True)
    if model == "PCM":
        settings["step_facet"] = step_facet
    if method == "MML":
        settings.update(quad_points=quad_points, population_prior_sd=float(population_prior_sd),
                        estimate_population_sd=False, mml_engine="EM")
    progress = progress or (lambda phase: None)
    progress("estimating")
    result = engine.mfrm_estimate(data=data, **settings)
    return _finish_analysis(result, data, settings, progress)


def current_tables(analysis):
    """Reapply current display qualification without refitting or rewriting provenance."""
    project.validate_analysis(analysis)
    result = analysis["result"]
    facts = reporting_facts.observed_report_tables(result)
    tables = [(title, table, note) for title, table, note in facts]
    locale_path = Path(__file__).resolve().parents[1] / "locales" / "en.json"
    fit_note = json.loads(locale_path.read_text(encoding="utf-8"))["metric_notes"]["fit"]
    measures, note = reporting_facts.publication_measure_table(result, analysis["diagnostics"], fit_note=fit_note)
    note = note.replace("Method details and source values: measures CSV.",
                        "Method details and source values are retained in this project.")
    tables.insert(1, ("Measures and Fit", measures, note))
    return tables


def category_figure(analysis):
    """Representative RSM probability curves at the net predictor, not raw scores."""
    from matplotlib.figure import Figure

    result = analysis["result"]
    if result["config"]["model"] != "RSM":
        raise ValueError("PCM category curves require a selected step-facet level")
    steps = np.asarray(result["params"]["steps"], dtype=float)
    if not np.isfinite(steps).all():
        raise ValueError("Category curves require finite step estimates")
    eta = np.linspace(min(-4., float(steps.min()) - 2.), max(4., float(steps.max()) + 2.), 201)
    probabilities = engine.category_prob_rsm(eta, np.r_[0., np.cumsum(steps)])
    figure = Figure(figsize=(7, 4), layout="constrained")
    axes = figure.subplots()
    for index, values in enumerate(probabilities.T):
        axes.plot(eta, values, label=str(index + result["prep"]["rating_min"]))
    axes.set(xlabel="Net predictor (logits)", ylabel="Category probability", ylim=(0, 1))
    axes.legend(title="Score", ncols=min(6, probabilities.shape[1]))
    from mfrm_app.visual_maps import style_figure, map_caution
    figure.suptitle("Rating scale · category probabilities", fontsize=15, fontweight="bold", color="#17344b")
    figure.supxlabel(map_caution(result), fontsize=8, color="#526274")
    return style_figure(figure)


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
