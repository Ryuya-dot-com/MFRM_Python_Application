"""Observed reporting facts and conservative manuscript boundaries; no estimation."""

import numpy as np
import pandas as pd
from mfrm_app import output_qualification


def current_manuscript_hold(result: dict) -> str:
    """Apply current policy to direct and reopened manuscript exports alike."""
    return output_qualification.calibration_inference_hold_message(result)


def observed_report_tables(result: dict) -> list[tuple[str, pd.DataFrame, str]]:
    """Use retained rows, not saved display summaries, without exposing identifiers."""
    prep = result.get("prep", {})
    config = result.get("config", {})
    data = prep.get("data")
    if not isinstance(data, pd.DataFrame) or data.empty or not {"Person", "Score"}.issubset(data):
        return []
    scores = pd.to_numeric(data["Score"], errors="coerce")
    if scores.isna().any() or not np.isfinite(scores).all() or data["Person"].isna().any():
        return []
    rows = [("Retained ratings", len(data)), ("Retained person identifiers", data["Person"].nunique())]
    for facet in ["Person", *list(config.get("facet_names", prep.get("facet_names", [])) or [])]:
        if facet not in data:
            continue
        counts = data.groupby(facet, observed=True).size()
        if facet != "Person":
            rows.append((f"{facet} levels", len(counts)))
        if not counts.empty:
            rows.append((f"Ratings per {facet} level", f"{int(counts.min())}–{int(counts.max())}"))
        if facet != "Person":
            incidence = data[["Person", facet]].drop_duplicates().groupby("Person", observed=True)[facet].nunique()
            rows.append((f"Unique {facet} levels per person", f"{int(incidence.min())}–{int(incidence.max())}"))
    audit = prep.get("response_data_audit", {})
    audit = audit.get("summary") if isinstance(audit, dict) else None
    if isinstance(audit, pd.DataFrame) and {"Item", "Count"}.issubset(audit):
        for name in ("Input rows", "Rows not used in likelihood"):
            value = pd.to_numeric(audit.loc[audit.Item.eq(name), "Count"], errors="coerce")
            if len(value) == 1 and np.isfinite(value.iloc[0]):
                rows.append((name, int(value.iloc[0])))
    design_note = (
        "Counts and workload ranges describe retained rating rows. Repeated ratings are not independent "
        "participants. Unobserved assignments cannot be classified as planned or accidental missingness "
        "from retained rows alone. Unique facet levels per person describe observed assignment overlap, "
        "not same-task agreement, connectivity, or measurement precision."
    )
    from mfrm_app.design_assignment import infer_rater_facet
    roles = config.get("facet_roles", {})
    mapping = infer_rater_facet(config.get("facet_names", []) or [], configured_rater_facet=(
        config.get("rater_facet") or (roles.get("rater") if isinstance(roles, dict) else None)))
    rater = mapping["rater_facet"]
    if mapping["mapping_confirmed"] and rater in data:
        incidence = data[["Person", rater]].dropna().drop_duplicates()
        if 0 < incidence[rater].nunique() <= 100:
            from scipy.sparse import csr_matrix
            persons, _ = pd.factorize(incidence["Person"])
            raters, levels = pd.factorize(incidence[rater])
            matrix = csr_matrix((np.ones(len(incidence), dtype=np.int64), (persons, raters)))
            shared = (matrix.T @ matrix).toarray()[np.triu_indices(len(levels), k=1)]
            positive = shared[shared > 0]
            rows.append((f"{rater} pairs sharing persons", len(positive)))
            rows.append((f"Shared persons per overlapping {rater} pair",
                         f"{positive.min()}–{positive.max()}" if len(positive) else "none"))
        elif incidence[rater].nunique() > 100:
            design_note += " Pair overlap counts are omitted above 100 rater levels; use the separate design-network diagnostics."
    design = pd.DataFrame(rows, columns=["Characteristic", "Value"])
    design["Value"] = design["Value"].astype(str)

    bounds = pd.to_numeric(pd.Series([prep.get("rating_min"), prep.get("rating_max")]), errors="coerce").to_numpy(dtype=float)
    # A restored range must not allocate an unbounded report table.
    support_recorded = (np.isfinite(bounds).all() and np.equal(bounds, np.floor(bounds)).all()
                        and 0 <= bounds[1] - bounds[0] < 1000
                        and np.abs(bounds).max() < np.iinfo(np.int64).max // 2)
    declared = np.arange(int(bounds[0]), int(bounds[1]) + 1) if support_recorded else np.array([])
    outside_support = bool(support_recorded and not scores.isin(declared).all())
    support = np.union1d(declared, scores.unique())
    counts = scores.value_counts().reindex(support, fill_value=0)
    categories = pd.DataFrame({"Category": support, "Count": counts.to_numpy(dtype=int)})
    categories["Percent"] = (100 * categories["Count"] / len(data)).round(2)
    category_note = (
        "Categories use the analyzed score coding. Count is the unweighted number of retained ratings; Percent uses all retained ratings as its "
        "denominator. Counts alone do not "
        "establish category functioning or justify merging categories."
    )
    categories.attrs["ReasonCode"] = ("report.score_support_unavailable" if not support_recorded else
                                       "report.observed_outside_recorded_support" if outside_support else "")
    if not support_recorded:
        category_note += " Recorded score support was unavailable or could not be expanded safely; only observed categories are listed."
    elif outside_support:
        category_note += " Observed scores outside the recorded support are retained here; review the saved score-range metadata."
    else:
        category_note += " All recorded score categories are listed, including zero-count categories when present."
    weights = pd.to_numeric(data.get("Weight", pd.Series(1., index=data.index)), errors="coerce")
    if not weights.isna().any() and np.isfinite(weights).all() and (weights > 0).all() and not weights.eq(1.).all():
        sums = weights.groupby(scores, observed=True).sum().reindex(support, fill_value=0.)
        categories["Weight sum"] = sums.to_numpy()
        categories["Weighted percent"] = (100 * sums.to_numpy() / weights.sum()).round(2)
        category_note += " Weight sum and Weighted percent report the separate likelihood weighting."
    return [("Analysis Design and Retained Ratings", design, design_note),
            ("Observed Category Use", categories, category_note)]


def observed_score_support_text(result: dict) -> str:
    """Recompute unused categories from retained ratings and bounded recorded support."""
    tables = observed_report_tables(result)
    if not tables:
        return ""
    categories = tables[1][1]
    if categories.attrs.get("ReasonCode"):
        return ("Recorded score support requires review; observed category counts are retained, "
                "but the full set of unused categories cannot be established reliably.")
    unused = categories.loc[categories.Count.eq(0), "Category"].astype(str).tolist()
    return (f"{categories.Category.min()} to {categories.Category.max()}; {len(categories)} categories; "
            f"unused categories retained: {', '.join(unused) if unused else 'none'}")


def recorded_method_details(result: dict) -> str:
    """Describe recorded settings, never filling unknown settings with defaults."""
    config = result.get("config", {})
    sentences = []
    settings = [("optimizer", "optimizer"), ("mml_engine", "MML engine"), ("maxit", "maximum iterations"),
                ("reltol", "stopping tolerance")]
    recorded = [f"{label} = {config[key]}" for key, label in settings if config.get(key) is not None]
    if recorded:
        sentences.append("Recorded estimation settings were " + "; ".join(recorded) + ".")
    signs = config.get("facet_signs", {})
    if isinstance(signs, dict) and signs:
        labels = [f"{name}: {int(sign):+d}" for name, sign in signs.items() if sign in (-1, 1)]
        if labels:
            sentences.append("Facet coefficients in the response predictor were " + ", ".join(labels)
                             + "; a positive coefficient increases the expected score as that measure increases.")
    specs = config.get("facet_specs", {})
    if isinstance(specs, dict) and specs:
        centered = [str(name) for name, spec in specs.items() if isinstance(spec, dict) and spec.get("centered")]
        if centered:
            sentences.append("Sum-to-zero constraints were applied to " + ", ".join(centered) + ".")
    noncenter = config.get("noncenter_facet")
    if noncenter:
        sentences.append(f"The selected non-centered facet was {noncenter}.")
    if config.get("step_facet"):
        sentences.append(f"Category-step parameters were indexed by {config['step_facet']}.")
    if config.get("weight_col"):
        sentences.append(f"Likelihood weights used the recorded column {config['weight_col']}; counts of ratings remain unweighted.")
    if str(config.get("method", "")).upper() == "MML":
        if config.get("quad_points") is not None:
            sentences.append(f"The recorded MML quadrature size was {config['quad_points']} nodes.")
        if config.get("estimate_population_sd"):
            sentences.append("The person population SD was estimated; its recorded terminal value was "
                             f"{config.get('estimated_population_sd', 'not recorded')}.")
        elif config.get("population_prior_sd") is not None:
            sentences.append(f"The person population SD was fixed at {config['population_prior_sd']}.")
        if config.get("population_formula"):
            sentences.append(f"The population regression formula was {config['population_formula']}.")
    return " ".join(sentences)


def optimizer_termination_text(result: dict) -> str:
    """A stopping flag is not a stationarity or uncertainty certificate."""
    opt = result.get("opt")
    stopped = getattr(opt, "success", None)
    if not isinstance(stopped, (bool, np.bool_)):
        text = "The optimizer termination status was not recorded."
    elif bool(stopped):
        text = "The optimizer reported that its stopping rule was reached."
    else:
        text = "The optimizer did not report successful termination."
    iterations = getattr(opt, "nit", None)
    if iterations is not None:
        text += f" Recorded iterations: {iterations}."
    text += " Optimizer termination alone does not establish stationarity, numerical accuracy, or valid inference."
    return text


def observed_results_text(result: dict) -> str:
    tables = observed_report_tables(result)
    if not tables:
        return ""
    design = tables[0][1].set_index("Characteristic")["Value"]
    categories = tables[1][1]
    unused = categories.loc[categories.Count.eq(0), "Category"].astype(str).tolist()
    category_text = (f"Observed scores occupied {int(categories.Count.gt(0).sum())} categories; "
                     "the recorded score support requires review. "
                     if categories.attrs.get("ReasonCode") else
                     f"{int(categories.Count.gt(0).sum())} of {len(categories)} recorded score categories were observed. "
                     f"Zero-count categories: {', '.join(unused) if unused else 'none'}. ")
    return (
        f"The retained data contained {design['Retained ratings']} ratings from "
        f"{design['Retained person identifiers']} person identifiers. "
        + category_text +
        "These are descriptive counts of the analyzed ratings, not evidence of model fit."
    )


def publication_measure_table(result: dict, diagnostics: dict, *, fit_note: str) -> tuple[pd.DataFrame, str]:
    """Keep uncertainty meaning visible within the Word/PDF column limits."""
    from mfrm_app.engine import qualify_measure_uncertainty

    measures = qualify_measure_uncertainty(diagnostics.get("measures"), result)
    if not isinstance(measures, pd.DataFrame) or measures.empty:
        return pd.DataFrame(), ""
    columns = [c for c in ("Facet", "Level", "Estimate", "Infit", "Outfit", "SE", "CI_Lower", "CI_Upper", "CI_Level")
               if c in measures]
    out = measures[columns].copy().round(3)
    status = measures.get("SE_Status", pd.Series("", index=measures.index)).astype(str)
    basis = status.map({"conditional_approximation": "Conditional", "posterior_eap_sd": "Posterior SD",
                        "ok": "Technical"}).fillna("Unverified")
    ci_status = measures.get("CI_Status", pd.Series("", index=measures.index)).astype(str)
    held = status.isin({"calibration_inference_withheld", "regularized", "rank_deficient",
                       "unsupported_covariance_scope"}) | ci_status.eq("calibration_inference_withheld")
    held |= bool(current_manuscript_hold(result))
    person_held = status.eq("person_inference_withheld") | ci_status.eq("person_inference_withheld")
    for column in ("FiniteJMLEEstimate", "PersonInferenceReady"):
        if column in measures:
            person_held |= measures.Facet.eq("Person") & measures[column].eq(False).fillna(False)
    held |= person_held
    basis.loc[held] = "Withheld"
    if "FiniteJMLEEstimate" in measures:
        endpoint = measures.get("ExtremeScoreDirection", pd.Series("", index=measures.index)).isin(
            {"all_minimum", "all_maximum"})
        basis.loc[person_held & endpoint & measures.FiniteJMLEEstimate.eq(False).fillna(False)] = "No finite JMLE"
    for column in ("SE", "CI_Lower", "CI_Upper", "CI_Level"):
        if column in out:
            out.loc[held, column] = np.nan
    for column in ("Estimate", "Infit", "Outfit"):
        if column in out:
            out.loc[person_held, column] = np.nan
    out = out.rename(columns={"SE": "SE / SD", "CI_Lower": "Lower", "CI_Upper": "Upper",
                              "CI_Level": "Interval level"})
    out["Basis"] = basis
    notes = {
        "Conditional": "Conditional: other fitted parameters are held fixed; intervals are screening approximations.",
        "Posterior SD": "Posterior SD: conditional posterior spread, not the sampling SE of EAP; intervals are normal approximations, not posterior quantiles.",
        "Technical": "Technical: structural covariance values are retained for audit; interval coverage has not been qualified.",
        "Unverified": "Unverified: the uncertainty method is not established; do not treat these numbers as validated SEs or intervals.",
        "Withheld": "Withheld: interpretation is on hold; estimates are technical terminal-point values. SE/SD and interval cells are blank.",
        "No finite JMLE": "No finite JMLE: all ratings for this person are at a scale endpoint. The person estimate, fit, SE and interval cells are blank; optimizer terminal values remain in the technical record.",
    }
    caption = (
        "Estimates and interval bounds are in logits. "
        "SE = standard error; SD = standard deviation. "
        + fit_note + " "
        "These are not automatic pass/fail judgments. "
    ) + " ".join(
        notes[label] for label in notes if basis.eq(label).any()
    ) + " Method details and source values: measures CSV."
    return out.fillna(""), caption
