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


def publication_measure_table(
    result: dict,
    diagnostics: dict,
    *,
    fit_note: str,
    extended: bool = False,
) -> tuple[pd.DataFrame, str]:
    """Shared qualified measure and fit table across tables, HTML, CSV, Word/PDF, and reloaded projects."""
    from mfrm_app.engine import qualify_measures_for_reporting

    diagnostics = diagnostics or {}
    measures = diagnostics.get("measures")
    if isinstance(measures, pd.DataFrame) and not measures.empty:
        for key in ("Facet", "Level"):
            if key not in measures or measures[key].astype("string").str.strip().fillna("").eq("").any():
                return pd.DataFrame(), f"Measure estimates are withheld because saved {key} identifiers are missing."
    qualified = qualify_measures_for_reporting(
        measures,
        result,
        obs_df=diagnostics.get("obs"),
    )
    if not isinstance(qualified, pd.DataFrame) or qualified.empty:
        return pd.DataFrame(), ""

    person = qualified.get("Facet", pd.Series("", index=qualified.index)).astype(str).eq("Person")
    person_held = pd.Series(False, index=qualified.index)
    for column in ("FiniteJMLEEstimate", "PersonInferenceReady"):
        if column in qualified.columns:
            person_held |= person & qualified[column].eq(False).fillna(False)
    unbounded_facet = qualified.get("UnboundedFacetEstimate", pd.Series(False, index=qualified.index)).eq(True).fillna(False)
    level_held = person_held | unbounded_facet
    # Only mask Estimate/Fit on level-specific boundary holds, not run-level screening-eligible holds.
    se_eligible = qualified["SEUse"].eq("screening_only")
    ci_eligible = qualified["CIUse"].eq("screening_only")

    est_col = pd.to_numeric(qualified.get("Estimate", pd.Series(np.nan, index=qualified.index)), errors="coerce").where(~level_held, np.nan).round(3)
    model_se_col = pd.to_numeric(qualified.get("ModelSE", pd.Series(np.nan, index=qualified.index)), errors="coerce").where(se_eligible, np.nan).round(3)
    real_se_col = pd.to_numeric(qualified.get("RealSE", pd.Series(np.nan, index=qualified.index)), errors="coerce").where(se_eligible, np.nan).round(3)
    lo_col = pd.to_numeric(qualified.get("CI_Lower", pd.Series(np.nan, index=qualified.index)), errors="coerce").where(ci_eligible, np.nan).round(3)
    hi_col = pd.to_numeric(qualified.get("CI_Upper", pd.Series(np.nan, index=qualified.index)), errors="coerce").where(ci_eligible, np.nan).round(3)
    lvl_col = pd.to_numeric(qualified.get("CI_Level", pd.Series(0.95, index=qualified.index)), errors="coerce").where(ci_eligible, np.nan).round(3)
    infit_col = pd.to_numeric(qualified.get("Infit", pd.Series(np.nan, index=qualified.index)), errors="coerce").where(~level_held, np.nan).round(3)
    outfit_col = pd.to_numeric(qualified.get("Outfit", pd.Series(np.nan, index=qualified.index)), errors="coerce").where(~level_held, np.nan).round(3)
    basis = qualified["Basis"].astype(str)

    if not extended:
        out = pd.DataFrame({
            "Facet": qualified["Facet"].astype(str),
            "Level": qualified["Level"].astype(str),
            "Estimate": est_col,
        })
        if "Infit" in qualified.columns:
            out["Infit"] = infit_col
        if "Outfit" in qualified.columns:
            out["Outfit"] = outfit_col
        if "SE" in qualified.columns or "ModelSE" in qualified.columns:
            out["SE / SD"] = model_se_col
        if "CI_Lower" in qualified.columns or extended:
            out["Lower"] = lo_col
        if "CI_Upper" in qualified.columns or extended:
            out["Upper"] = hi_col
        if "CI_Level" in qualified.columns:
            out["Interval level"] = lvl_col
        out["Basis"] = basis
    else:
        infit_z = pd.to_numeric(qualified.get("InfitZSTD", pd.Series(np.nan, index=qualified.index)), errors="coerce").where(~level_held, np.nan).round(3)
        outfit_z = pd.to_numeric(qualified.get("OutfitZSTD", pd.Series(np.nan, index=qualified.index)), errors="coerce").where(~level_held, np.nan).round(3)
        df_infit = pd.to_numeric(qualified.get("DF_Infit", pd.Series(np.nan, index=qualified.index)), errors="coerce").where(~level_held, np.nan).round(3)
        df_outfit = pd.to_numeric(qualified.get("DF_Outfit", pd.Series(np.nan, index=qualified.index)), errors="coerce").where(~level_held, np.nan).round(3)
        out = pd.DataFrame({
            "Facet": qualified["Facet"].astype(str),
            "Level": qualified["Level"].astype(str),
            "Estimate": est_col,
            "ModelSE": model_se_col,
            "RealSE": real_se_col,
            "SE / SD": model_se_col,
            "Lower": lo_col,
            "Upper": hi_col,
            "Interval level": lvl_col,
            "Infit": infit_col,
            "Outfit": outfit_col,
            "InfitZSTD": infit_z,
            "OutfitZSTD": outfit_z,
            "DF_Infit": df_infit,
            "DF_Outfit": df_outfit,
            "FitStatus": qualified["FitStatus"].astype(str),
            "FitReason": qualified["FitReason"].astype(str),
            "Basis": basis,
            "SEUse": qualified["SEUse"].astype(str),
            "CIUse": qualified["CIUse"].astype(str),
            "SupportsFormalInference": qualified["SupportsFormalInference"].astype(bool),
            "StatusReason": qualified["StatusReason"].astype(str),
        })
        if "ObservedCategories" in qualified.columns and qualified["ObservedCategories"].notna().any():
            out["ObservedCategories"] = qualified["ObservedCategories"].astype("Int64")
        if "SavedEstimate" in qualified.columns and unbounded_facet.any():
            sv_est = pd.to_numeric(qualified["SavedEstimate"], errors="coerce").round(3)
            if (unbounded_facet & sv_est.notna()).any():
                out["SavedEstimate"] = sv_est.where(unbounded_facet, np.nan)
        if "SavedUncertaintyReasonCode" in qualified.columns:
            saved_rc = qualified["SavedUncertaintyReasonCode"].fillna("").astype(str)
            if (saved_rc.ne("") & saved_rc.ne(out["StatusReason"])).any():
                out["SavedReasonCode"] = saved_rc
        if "SavedSE_Status" in qualified.columns:
            saved_st = qualified["SavedSE_Status"].fillna("").astype(str)
            if (saved_st.ne("") & saved_st.ne(qualified["EffectiveSE_Status"].astype(str))).any():
                out["SavedSE_Status"] = saved_st
        for raw_z, saved_z in (("InfitZSTD", "SavedInfitZSTD"), ("OutfitZSTD", "SavedOutfitZSTD")):
            if saved_z in qualified.columns:
                sv = pd.to_numeric(qualified[saved_z], errors="coerce").round(3)
                cur = pd.to_numeric(out[raw_z], errors="coerce")
                changed = sv.notna() & (cur.isna() | ~np.isclose(sv.fillna(0.0), cur.fillna(0.0), atol=1e-6))
                if changed.any():
                    out[saved_z] = sv

    notes = {
        "Conditional": (
            "Conditional: other fitted parameters are held fixed; "
            "exploratory normal bands are screening_only (SupportsFormalInference=FALSE)."
        ),
        "Conditional (audit size limit)": (
            "Conditional (audit size limit): exact dense JMLE rank audit exceeded its size limit "
            "(stat.jmle.identifiability_audit_size_limit); conditional SEs and intervals are screening_only (SupportsFormalInference=FALSE)."
        ),
        "Conditional (fallback)": (
            "Conditional (fallback): observation-table conditional SE fallback; "
            "intervals are screening_only (SupportsFormalInference=FALSE)."
        ),
        "Posterior SD": (
            "Posterior SD: conditional posterior spread given fitted parameters, not the sampling SE of EAP; "
            "normal bands are screening_only (SupportsFormalInference=FALSE)."
        ),
        "Observed info (fixed-SD)": (
            "Observed info (fixed-SD): MML observed-information delta-method SE under fixed population SD; "
            "general-table normal bands are screening_only (SupportsFormalInference=FALSE; use Specified intervals for target-checked pointwise intervals)."
        ),
        "Fixed constraint": "Fixed constraint: parameter fixed by anchor or single-level constraint; SE and interval are not applicable.",
        "Unverified": "Unverified: the uncertainty method is not established; SE/SD and interval cells are blank.",
        "Unavailable": "Unavailable: finite positive standard error was not available; SE/SD and interval cells are blank.",
        "Withheld (non-converged)": (
            "Withheld (non-converged): optimizer termination was not established "
            "(stat.estimation.optimizer_termination_unverified); estimates are technical terminal-point values and SE/SD and interval cells are blank."
        ),
        "Withheld (unidentifiable)": (
            "Withheld (unidentifiable): linear predictor is structurally rank-deficient "
            "(stat.jmle.eta_rank_deficient); estimates are technical terminal-point values and SE/SD and interval cells are blank."
        ),
        "Withheld (rank-deficient covariance)": (
            "Withheld (rank-deficient covariance): MML observed-information covariance is rank-deficient or regularized; "
            "SE/SD and interval cells are blank."
        ),
        "Withheld": (
            "Withheld: interpretation is on hold; estimates are technical terminal-point values. SE/SD and interval cells are blank."
        ),
        "No finite JMLE": (
            "No finite JMLE: all ratings for this person are at a scale endpoint. "
            "The person estimate, fit, SE and interval cells are blank; optimizer terminal values remain in the technical record."
        ),
    }
    fit_statuses = set(qualified.get("FitStatus", pd.Series(dtype=str)).astype(str))
    fit_guard_note = ""
    if "zstd_unavailable_boundary_proxy" in fit_statuses:
        fit_guard_note += " ZSTD is blank where scope includes JMLE extreme-score boundary proxies."
    if "zstd_unavailable_low_df" in fit_statuses:
        fit_guard_note += " ZSTD is blank where degrees of freedom are below 1."
    if "unavailable_zero_df" in fit_statuses:
        fit_guard_note += " Fit is blank where a free JMLE parameter has a single observation (0 residual degrees of freedom)."
    if "unavailable_boundary" in fit_statuses:
        fit_guard_note += " Estimate and fit are blank where a free facet or step-facet level has an unbounded endpoint or category-gap pattern."
    if "unavailable" in fit_statuses:
        fit_guard_note += " Fit is blank where response variance, residuals or weights are unavailable (no partial-row fit)."
    if extended and "ObservedCategories" in qualified.columns and (qualified["ObservedCategories"] == 1).fillna(False).any():
        fit_guard_note += " Levels with ObservedCategories = 1 used only a single score category (zero raw score variance)."

    reason_set = set(qualified.get("StatusReason", pd.Series(dtype=str)).astype(str))
    hold_detail = ""
    if output_qualification.OPTIMIZER_TERMINATION_HOLD_REASON in reason_set:
        hold_detail += " Hold cause: non-converged optimizer (stat.estimation.optimizer_termination_unverified)."
    if output_qualification.JMLE_RANK_DEFICIENT_HOLD_REASON in reason_set:
        hold_detail += " Hold cause: structurally rank-deficient / unidentifiable JMLE design (stat.jmle.eta_rank_deficient)."
    if "stat.stat_003.covariance_rank_deficient" in reason_set:
        hold_detail += " Hold cause: rank-deficient or regularized MML covariance (stat.stat_003.covariance_rank_deficient)."
    for b_code in (
        output_qualification.CATEGORY_UNSUPPORTED_HOLD_REASON,
        output_qualification.UNOBSERVED_OWNER_ENDPOINT_HOLD_REASON,
        output_qualification.UNOBSERVED_BOUNDARY_CATEGORY_HOLD_REASON,
        output_qualification.MML_FACET_BOUNDARY_HOLD_REASON,
        output_qualification.MML_DISCONNECTED_HOLD_REASON,
    ):
        if b_code in reason_set:
            hold_detail += f" Hold cause: {output_qualification.boundary_and_support_hold_message(b_code)} ({b_code})"

    caption = (
        "Estimates and interval bounds are in logits. "
        "SE = standard error; SD = standard deviation; RealSE = ModelSE * sqrt(max(Infit, 1)). "
        "Displayed general-table SEs and intervals are exploratory (Use = screening_only; SupportsFormalInference = FALSE). "
        + fit_note + " "
        "These are not automatic pass/fail judgments."
        + fit_guard_note + " "
        + " ".join(notes[label] for label in notes if basis.eq(label).any())
        + hold_detail
        + " Method details and source values: measures CSV."
    )
    for column in out.columns[out.isna().any()]:
        if out[column].dtype == bool:
            continue
        out[column] = out[column].astype(object).where(out[column].notna(), "")
    return out, caption


def qualified_reliability_table(result: dict, diagnostics: dict) -> pd.DataFrame:
    """Qualify current or saved 1-row reliability without replacing saved numbers."""
    from mfrm_app.engine import qualify_measure_uncertainty, qualify_measures_for_reporting
    from mfrm_app.reliability import qualify_saved_reliability

    result = result or {}
    diagnostics = diagnostics or {}
    reason = output_qualification.calibration_screening_block_reason(result)
    measures = diagnostics.get("measures")
    if isinstance(measures, pd.DataFrame) and not measures.empty:
        if reason or "SE_Method" not in measures.columns:
            measures = qualify_measure_uncertainty(measures, result)
        else:
            q = qualify_measures_for_reporting(measures, result, obs_df=diagnostics.get("obs"))
            measures = measures.copy()
            measures["SE_Status"] = q["EffectiveSE_Status"]
    return qualify_saved_reliability(
        diagnostics.get("reliability"),
        measures,
        calibration_hold=bool(reason),
        hold_reason=reason,
    )


def publication_reliability_table(
    result: dict,
    diagnostics: dict,
    *,
    rel_note: str,
    four_series: bool = False,
) -> tuple[pd.DataFrame, str]:
    """Qualified facet separation/reliability table (compact 1-row or full 4-series Model/Real × Sample/Population)."""
    from mfrm_app.reliability import calculate_reliability_series

    result = result or {}
    diagnostics = diagnostics or {}
    if not four_series:
        reliability = qualified_reliability_table(result, diagnostics)
        if not isinstance(reliability, pd.DataFrame) or reliability.empty:
            return pd.DataFrame(), ""
        columns = [c for c in ("Facet", "Levels", "SD", "RMSE", "Separation", "Strata", "Reliability")
                   if c in reliability]
        out = reliability[columns].copy().round(3)
        block_reason = output_qualification.calibration_screening_block_reason(result)
        held = reliability["ReliabilityStatus"].ne("diagnostic") | bool(block_reason)
        out["Use"] = np.where(held, "Withheld", "Diagnostic")
        caption = "Reliability and separation by facet (Use = screening_only / Diagnostic; SupportsFormalInference = FALSE). "
        if held.any():
            for column in ("RMSE", "Separation", "Strata", "Reliability"):
                if column in out:
                    out.loc[held, column] = np.nan
            caption += (
                "Interpretation is withheld for rows whose calibration, paired inputs or uncertainty basis are unqualified, "
                "or whose saved calculation cannot be verified. Precision indices are intentionally blank; "
                "SD is the spread of technical terminal-point estimates. "
            )
        else:
            caption += "These descriptive indices depend on the uncertainty method, observed spread and rating design; they do not establish accuracy for other designs. "
        caption += rel_note + " "
        caption += (
            "SD = sample standard deviation of estimates; RMSE = root mean square of SEs (posterior SDs for EAP). "
            "For SE-based measures, G = sqrt(max(SD^2 - RMSE^2, 0)) / RMSE. "
            "Reliability = G^2 / (1 + G^2). "
        )
        basis = reliability.get("ReliabilityBasis", pd.Series(dtype=str)).astype(str)
        if basis.str.contains("EAP", regex=False).any():
            caption += (
                "For EAP person estimates, G = SD(EAP) / RMSE and reliability = "
                "Var(EAP) / (Var(EAP) + mean(posterior SD^2)); this is not directly equivalent to FACETS separation reliability. "
            )
        caption += "Full notes and original technical values are retained in the reliability CSV."
        for column in out.columns[out.isna().any()]:
            out[column] = out[column].astype(object).where(out[column].notna(), "")
        return out, caption

    series = calculate_reliability_series(
        diagnostics.get("measures"),
        result,
        obs_df=diagnostics.get("obs"),
        saved_reliability=diagnostics.get("reliability"),
    )
    if series.empty:
        return pd.DataFrame(), ""

    held = series["Use"].ne("screening_only")
    out = pd.DataFrame({
        "Facet": series["Facet"].astype(str),
        "SEMode": series["SEMode"].astype(str),
        "SEColumn": series["SEColumn"].astype(str),
        "DistributionBasis": series["DistributionBasis"].astype(str),
        "Levels": pd.to_numeric(series["Levels"], errors="coerce").astype("Int64"),
        "EstimateAvailable": pd.to_numeric(series["EstimateAvailable"], errors="coerce").astype("Int64"),
        "SEAvailable": pd.to_numeric(series["SEAvailable"], errors="coerce").astype("Int64"),
        "ExcludedLevels": pd.to_numeric(series["ExcludedLevels"], errors="coerce").astype("Int64"),
        "AggregationTarget": series["AggregationTarget"].astype(str),
        "SD": pd.to_numeric(series["SD"], errors="coerce").round(3),
        "RMSE": pd.to_numeric(series["RMSE"], errors="coerce").where(~held, np.nan).round(3),
        "TrueSD": pd.to_numeric(series["TrueSD"], errors="coerce").where(~held, np.nan).round(3),
        "Separation": pd.to_numeric(series["Separation"], errors="coerce").where(~held, np.nan).round(3),
        "Strata": pd.to_numeric(series["Strata"], errors="coerce").where(~held, np.nan).round(3),
        "Reliability": pd.to_numeric(series["Reliability"], errors="coerce").where(~held, np.nan).round(3),
        "Basis": series["ReliabilityBasisCode"].map({
            "eap_posterior_sd": "EAP posterior SD",
            "se_variance_correction": "SE variance subtraction",
        }).fillna("Unverified"),
        "Use": series["Use"].astype(str),
        "SupportsFormalInference": series["SupportsFormalInference"].astype(bool),
        "StatusReason": series["StatusReason"].astype(str),
    })
    if "SavedReliability" in series.columns and series["SavedReliability"].notna().any():
        out["SavedReliability"] = pd.to_numeric(series["SavedReliability"], errors="coerce").round(3)
    if "SavedReliabilityReasonCode" in series.columns:
        saved_rc = series["SavedReliabilityReasonCode"].fillna("").astype(str)
        if (saved_rc.ne("") & saved_rc.ne(out["StatusReason"])).any():
            out["SavedReasonCode"] = saved_rc

    caption = (
        "Four descriptive precision series per facet: SEMode = model (ModelSE) vs fit_adjusted "
        "(RealSE = ModelSE * sqrt(max(Infit, 1))) crossed with DistributionBasis = sample (ddof=1) vs population (ddof=0). "
        "All derived separation/reliability indices are exploratory summaries (Use = screening_only; SupportsFormalInference = FALSE) "
        "and do not propagate inter-estimate covariance. "
        "EstimateAvailable, SEAvailable, ExcludedLevels and AggregationTarget state the exact calculation scope; "
        "JMLE extreme-score boundary persons without finite estimates are excluded explicitly rather than silently dropped. "
    )
    if held.any():
        caption += (
            "Precision cells (RMSE, TrueSD, Separation, Strata, Reliability) are blank where Use = withheld; "
            "inspect StatusReason for the exact cause. "
        )
    if rel_note:
        caption += rel_note + " "
    if series["ReliabilityBasisCode"].eq("eap_posterior_sd").any():
        caption += (
            "For MML EAP person estimates, Separation = SD(EAP) / RMSE and Reliability = Var(EAP) / (Var(EAP) + RMSE^2); "
            "conditional posterior SD is not the sampling SE of EAP. "
        )
    caption += "Full notes and original technical values are retained in this project."
    for column in out.columns[out.isna().any()]:
        if out[column].dtype == bool:
            continue
        out[column] = out[column].astype(object).where(out[column].notna(), "")
    return out, caption
