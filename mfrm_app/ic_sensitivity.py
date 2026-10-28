"""Optional fixed-SD MML sensitivity checks and refits; no inference qualification."""
from copy import deepcopy
from itertools import combinations
import json
import operator

import numpy as np
import pandas as pd

from . import engine
from .output_qualification import optimizer_termination_established


MODELS = ("RSM", "PCM", "GPCM")
SCHEMA = "mfrm_ic_same_point_v1"
REFINEMENT_SCHEMA = "mfrm_ic_refinement_v1"


def capture_inputs(fits):
    """Retain only the inputs needed to check the three candidate points later."""
    captured = {}
    for model in MODELS:
        fit = fits[model]
        cfg, prep = fit["config"], fit["prep"]
        if (cfg["method"] != "MML" or cfg["model"] != model or cfg.get("estimate_population_sd")
                or cfg.get("population_model", {}).get("enabled") or cfg.get("anchors_active")
                or fit.get("anchors", {}).get("active") or cfg.get("facet_regularization_enabled")):
            raise ValueError("unsupported_scope")
        q = cfg["quad_points"]
        if isinstance(q, (bool, np.bool_)) or not 1 <= operator.index(q) <= 121:
            raise ValueError("unsupported_scope")
        data = prep["data"]
        if not 0 < len(data) <= 1000 or not np.all(data["Weight"] == 1):
            raise ValueError("unsupported_scope")
        facets = cfg["facet_names"]
        if prep["facet_names"] != facets or cfg["n_person"] != len(prep["levels"]["Person"]):
            raise ValueError("inconsistent_inputs")
        for f in ["Person", *facets]:
            if list(data[f].cat.categories) != prep["levels"][f] or (data[f].cat.codes < 0).any():
                raise ValueError("inconsistent_inputs")
        if not np.array_equal(data["score_k"], data["Score"]-prep["rating_min"]):
            raise ValueError("inconsistent_inputs")
        point = np.asarray(fit["opt"].x, dtype=float)
        k = sum(engine.build_param_sizes(cfg).values())
        row = fit["summary"].iloc[0]
        n = cfg["n_person"]
        if (k < 1 or point.shape != (k,) or not np.isfinite(point).all() or row["KParams"] != k
                or row["Persons"] != n or row["N"] != len(data) or not np.isfinite(row["LogLik"])):
            raise ValueError("inconsistent_inputs")
        # Exclude bulky diagnostics and fitted result tables from the retained bundle.
        small_prep = {key:prep[key] for key in ["data", "levels", "facet_names", "rating_min", "rating_max"]}
        entry = deepcopy(dict(config=cfg, prep=small_prep, point=point,
            reported_nll=-float(row["LogLik"]), K=k, persons=n,
            analysis_id=engine.build_result_analysis_identity(fit).analysis_id,
            terminated=optimizer_termination_established(fit)))
        entry["point_id"] = engine.stable_json_fingerprint(entry, length=64)
        captured[model] = entry
    return captured


def input_id(inputs):
    """Bind cached measurements to all retained data, settings and coordinates."""
    return engine.stable_json_fingerprint(dict(schema=SCHEMA, inputs=inputs), length=64)


def _prepared_inputs(inputs):
    """Check every candidate before allowing evaluation or optimization."""
    if set(inputs) != set(MODELS):
        raise ValueError("incomplete_candidates")
    prepared, signatures = [], []
    for model in MODELS:
        entry = inputs[model]
        payload = {key:value for key,value in entry.items() if key != "point_id"}
        if engine.stable_json_fingerprint(payload, length=64) != entry["point_id"]:
            raise ValueError("changed_inputs")
        cfg, prep, point = entry["config"], entry["prep"], entry["point"]
        q = operator.index(cfg["quad_points"])
        signatures.append(engine.stable_json_fingerprint(dict(
            data=prep["data"], levels=prep["levels"], rating_min=prep["rating_min"], rating_max=prep["rating_max"],
            facets=cfg["facet_specs"], signs=cfg["facet_signs"], q=q, sd=cfg["population_prior_sd"])))
        if len(set(signatures)) != 1:
            raise ValueError("different_likelihoods")
        sizes = engine.build_param_sizes(cfg)
        idx = engine.build_indices(prep, step_facet=cfg.get("step_facet"), slope_facet=cfg.get("slope_facet"))
        prepared.append((model, entry, sizes, idx))
    return prepared


def evaluate(inputs):
    """Reevaluate each unchanged point at its original and two finer GH rules.

    Changes are numerical sensitivity, not error bounds. A small difference
    between two rules is not proof of integration accuracy or a maximum.
    """
    rows = []
    for model, entry, sizes, idx in _prepared_inputs(inputs):
        cfg, point = entry["config"], entry["point"]
        q = operator.index(cfg["quad_points"])
        for order in [q, max(61, 2*q+1), max(121, 4*q+1)]:
            rule = engine.make_mml_quadrature(cfg, order)
            nll, gradient = engine.mfrm_loglik_mml_value_grad(point, idx, cfg, sizes, rule)
            if not np.isfinite(nll) or not np.isfinite(gradient).all() or np.shape(gradient) != np.shape(point):
                raise ValueError("nonfinite_evaluation")
            if order == q and abs(nll-entry["reported_nll"]) > 1e-8:
                raise ValueError("reported_likelihood_mismatch")
            rows.append(dict(Model=model, Points=order, OriginalPoints=q, LogLik=-float(nll),
                AIC=2*nll+2*entry["K"], BIC=2*nll+np.log(entry["persons"])*entry["K"],
                ICChange=2*(nll-entry["reported_nll"]), GradientSupnorm=float(np.max(np.abs(gradient))),
                KParams=entry["K"], BICSampleSize=entry["persons"], FixedPopulationSD=cfg["population_prior_sd"],
                AlgorithmTerminated=entry["terminated"], CandidateAnalysisID=entry["analysis_id"], PointID=entry["point_id"],
                Refitted=False, ContinuousReferenceChecked=False, InferenceQualified=False))
    table = pd.DataFrame(rows)
    finest = table.groupby("Model", sort=False).tail(1).set_index("Model")
    pairs = []
    for first, second in combinations(MODELS, 2):
        a, b = finest.loc[first], finest.loc[second]
        for criterion in ["AIC", "BIC"]:
            fine = b[criterion]-a[criterion]
            shift = b.ICChange-a.ICChange
            pairs.append(dict(First=first, Second=second, Criterion=criterion,
                OriginalDifference=fine-shift, FinerRuleDifference=fine, DifferenceChange=shift,
                FirstPointID=a.PointID, SecondPointID=b.PointID,
                Refitted=False, ContinuousReferenceChecked=False, InferenceQualified=False))
    return dict(schema=SCHEMA, input_id=input_id(inputs), evaluations=table, pairs=pd.DataFrame(pairs),
        refitted=False, continuous_reference_checked=False, inference_qualified=False)


def refine(inputs):
    """Refit on a finer GH rule, retaining every outcome and an extra same-point check.

    Uses the existing constrained direct-MML optimizer. These are technical
    candidates, not replacement fits or precision/model-selection certificates.
    """
    prepared = _prepared_inputs(inputs)
    # Replay all original likelihoods before any model is optimized.
    for model, entry, sizes, idx in prepared:
        cfg = entry["config"]
        original, _ = engine.mfrm_loglik_mml_value_grad(entry["point"], idx, cfg, sizes, engine.make_mml_quadrature(cfg))
        if not np.isfinite(original) or abs(original-entry["reported_nll"]) > 1e-8:
            raise ValueError("reported_likelihood_mismatch")
    rows, points = [], {}
    for model, entry, sizes, idx in prepared:
        cfg = deepcopy(entry["config"])
        order = max(121, 4*operator.index(cfg["quad_points"])+1)
        check_order = 2*order+1
        cfg.update(quad_points=order, mml_engine="direct", maxit=400, reltol=1e-12)
        rule = engine.make_mml_quadrature(cfg)
        k, n = entry["K"], entry["persons"]
        row = dict(Model=model, FitPoints=order, CheckPoints=check_order, KParams=k, BICSampleSize=n,
            FixedPopulationSD=cfg["population_prior_sd"], OriginalPointID=entry["point_id"],
            OriginalCandidateAnalysisID=entry["analysis_id"], OriginalAlgorithmTerminated=entry["terminated"],
            OriginalAIC=2*entry["reported_nll"]+2*k, OriginalBIC=2*entry["reported_nll"]+np.log(n)*k,
            CalculationStatus="failed", CalculationReason="", AlgorithmTerminated=False, Feasible=False, Boundary=None,
            AIC=np.nan, BIC=np.nan, ICChangeFromOriginal=np.nan, CheckICChange=np.nan,
            GradientSupnorm=np.nan, CheckGradientSupnorm=np.nan,
            Refitted=True, OriginalResultReplaced=False, ContinuousReferenceChecked=False, InferenceQualified=False)
        try:
            method, options = engine.optimizer_method_and_options(sizes, cfg, 400, 1e-12, 1e-12)
            row.update(Optimizer=method, OptimizerOptions=json.dumps(options, sort_keys=True))
            initial, _ = engine.mfrm_loglik_mml_value_grad(entry["point"], idx, cfg, sizes, rule)
            fitted = engine.mfrm_direct_mml(entry["point"], idx, cfg, sizes, rule, maxit=400, reltol=1e-12)
            point = np.asarray(fitted.x, dtype=float)
            row.update(InitialFineNLL=float(initial), AlgorithmTerminated=optimizer_termination_established({"opt":fitted}),
                OptimizerMessage=str(fitted.message), Iterations=int(fitted.nit), FunctionEvals=int(fitted.nfev),
                MaxIterations=400, RequestedTolerance=1e-12)
            if point.shape != entry["point"].shape or not np.isfinite(point).all():
                raise ValueError("invalid_candidate")
            points[model] = point.tolist()
            row["CandidatePointID"] = engine.stable_json_fingerprint(dict(config=cfg, point=point, parent=entry["point_id"]), length=64)
            margins = []
            for value, (lower, upper) in zip(point, engine.build_optimizer_bounds(sizes,cfg) or []):
                if lower is not None: margins.append(value-lower)
                if upper is not None: margins.append(upper-value)
            for constraint in engine.build_optimizer_constraints(sizes,cfg):
                value = constraint.A @ point
                margins.extend(np.ravel(value-constraint.lb))
                margins.extend(np.ravel(constraint.ub-value))
            row["Feasible"] = bool(all(m >= -1e-8 for m in margins))
            row["Boundary"] = bool(any(m <= 1e-6 for m in margins))
            if not row["Feasible"]:
                raise ValueError("infeasible_candidate")
            nll, grad = engine.mfrm_loglik_mml_value_grad(point, idx, cfg, sizes, rule)
            checked, check_grad = engine.mfrm_loglik_mml_value_grad(point, idx, cfg, sizes, engine.make_mml_quadrature(cfg,check_order))
            if not (np.isfinite([initial,nll,checked]).all() and np.isfinite(grad).all() and np.isfinite(check_grad).all()
                    and np.shape(grad) == point.shape and np.shape(check_grad) == point.shape):
                raise ValueError("nonfinite_evaluation")
            row.update(LogLik=-float(nll), CheckLogLik=-float(checked), AIC=2*nll+2*k, BIC=2*nll+np.log(n)*k,
                ICChangeFromOriginal=2*(nll-entry["reported_nll"]), CheckICChange=2*(checked-nll),
                GradientSupnorm=float(np.max(abs(grad))), CheckGradientSupnorm=float(np.max(abs(check_grad))),
                ObjectiveImproved=bool(nll <= initial+1e-8))
            row["CalculationStatus"] = "completed" if row["AlgorithmTerminated"] and row["ObjectiveImproved"] else "unfinished"
            if not row["ObjectiveImproved"]:
                row["CalculationReason"] = "objective_worsened"
            elif not row["AlgorithmTerminated"]:
                row["CalculationReason"] = "termination_unconfirmed"
        except (ValueError, TypeError, KeyError, AttributeError, IndexError, ArithmeticError, RuntimeError) as exc:
            row["Failure"] = str(exc)
            row["CalculationReason"] = str(exc) if str(exc) in (
                "invalid_candidate", "infeasible_candidate", "nonfinite_evaluation"
            ) else "calculation_failed"
        rows.append(row)
    table = pd.DataFrame(rows).set_index("Model", drop=False)
    pairs = []
    for first, second in combinations(MODELS,2):
        a,b = table.loc[first], table.loc[second]
        included = a.CalculationStatus == b.CalculationStatus == "completed"
        for criterion in ["AIC","BIC"]:
            pairs.append(dict(First=first, Second=second, Criterion=criterion, Included=bool(included),
                Difference=b[criterion]-a[criterion] if included else np.nan,
                CheckedDifference=b[criterion]+b.CheckICChange-a[criterion]-a.CheckICChange if included else np.nan,
                Refitted=True, OriginalResultReplaced=False, ContinuousReferenceChecked=False, InferenceQualified=False))
    return dict(schema=REFINEMENT_SCHEMA, input_id=input_id(inputs), candidates=table.reset_index(drop=True),
        points=points, pairs=pd.DataFrame(pairs), refitted=True, original_result_replaced=False,
        continuous_reference_checked=False, inference_qualified=False)
