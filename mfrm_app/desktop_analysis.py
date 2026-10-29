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


def _gd_settings(settings):
    if not isinstance(settings, dict) or set(settings) - {
        "object_facet", "random_facets", "score_col", "n_max_per_facet", "target_population", "selected_plan"
    }:
        raise ValueError("Unknown G/D settings")
    out = dict(score_col="Score", n_max_per_facet=10, target_population="", selected_plan=None)
    out.update(deepcopy(settings))
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
    return out


def _run_gd(data, settings):
    arguments = {k: v for k, v in settings.items() if k != "selected_plan"}
    result = gd.compute_generalizability_from_data(data, **arguments)
    if result["available"]:
        counts = settings["selected_plan"] or result["design"]["observed_n"]
        plans = result["d_study"]
        if not np.logical_and.reduce([plans[f].eq(n) for f, n in counts.items()]).any():
            raise ValueError("The selected D plan is outside the computed plan grid")
        coefficients = gd._g_phi_from_variance_components(result["variance_components"],
            object_facet=settings["object_facet"], facet_observations=counts)
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
        for section, table in frame.groupby("_section", sort=False):
            tables.append(("G/D: " + section, table.drop(columns="_section").dropna(axis=1, how="all"),
                           gd._GD_CI_REASON))
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
