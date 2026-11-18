"""Descriptive reliability with complete paired inputs, 4-series Model/Real × Sample/Population summaries, and saved provenance."""
from __future__ import annotations

import numpy as np
import pandas as pd

from mfrm_app import output_qualification


CALCULATION_VERSION = "paired-reliability-v2"
SERIES_CALCULATION_VERSION = "paired-reliability-4series-v1"
PRECISION_COLUMNS = ("RMSE", "Separation", "Strata", "Reliability")


def calculate_reliability(measures):
    rows = []
    if not isinstance(measures, pd.DataFrame) or measures.empty:
        return pd.DataFrame()
    for facet, frame in measures.groupby("Facet"):
        estimate = pd.to_numeric(frame["Estimate"], errors="coerce")
        se = pd.to_numeric(frame["SE"], errors="coerce")
        valid = np.isfinite(estimate) & np.isfinite(se) & se.gt(0)
        with np.errstate(over="ignore", invalid="ignore"):
            mv = float(estimate.var(ddof=1)) if len(frame) >= 2 and np.isfinite(estimate).all() else np.nan
            ev = float(np.mean(se ** 2)) if valid.all() else np.nan
        methods = frame.get("SE_Method", pd.Series("", index=frame.index)).fillna("").astype(str).str.lower()
        statuses = frame.get("SE_Status", pd.Series("", index=frame.index)).fillna("").astype(str)
        person_held = pd.Series(False, index=frame.index)
        if facet == "Person":
            person_held |= statuses.eq("person_inference_withheld")
            for column in ("FiniteJMLEEstimate", "PersonInferenceReady"):
                if column in frame:
                    person_held |= frame[column].eq(False).fillna(False)
        posterior = methods.str.contains("latent-ability posterior sd", regex=False)
        se_based = methods.str.startswith(("conditional information approximation", "mml observed-information delta method"))
        legacy_posterior = statuses.eq("posterior_eap_sd")
        eap = bool(posterior.all() or (methods.eq("").all() and legacy_posterior.all()))
        known_basis = bool(posterior.all() or se_based.all())
        mixed = bool((posterior.any() and not posterior.all()) or
                     (methods.eq("").all() and legacy_posterior.any() and not legacy_posterior.all()))
        status, reason, note = "diagnostic", "stat.reliability.descriptive_only", ""
        separation = reliability = strata = np.nan
        if not valid.all():
            status, reason = "unassessable", "stat.reliability.invalid_measure_pairs"
            note = f"{int((~valid).sum())} level(s) have a non-finite estimate or a non-finite/non-positive SE; no levels were silently dropped."
        elif len(frame) < 2:
            status, reason = "unassessable", "stat.reliability.insufficient_levels"
            note = "At least two complete estimate/SE pairs are required."
        elif mixed:
            status, reason = "unassessable", "stat.reliability.mixed_uncertainty_basis"
            note = "Posterior SD and SE bases are mixed within this facet."
        elif not np.isfinite(mv) or not np.isfinite(ev) or ev <= 0:
            status, reason = "unassessable", "stat.reliability.invalid_variance"
            note = "Observed or error variance is not finite and usable."
        elif mv == 0 and not eap:
            status, reason = "unassessable", "stat.reliability.zero_observed_variance"
            note = "Identical estimates give zero observed variance; SE-based reliability has a zero denominator."
        else:
            signal = mv if eap else max(mv - ev, 0.)
            separation = float(np.sqrt(signal / ev))
            reliability = float(mv / (mv + ev) if eap else signal / mv)
            strata = float((4 * separation + 1) / 3)
            if signal == 0:
                note = "Zero observed EAP variance; empirical reliability = 0." if eap else "Measurement variance ≤ error variance; separation = 0."
            if statuses.eq("calibration_inference_withheld").any():
                status, reason = "withheld", "stat.reliability.calibration_inference_withheld"
            elif not known_basis or not statuses.isin({"posterior_eap_sd"} if eap else {"conditional_approximation", "ok"}).all():
                status, reason = "unverified", "stat.reliability.uncertainty_basis_unverified"
                note += " SE method/status provenance is incomplete or unqualified; values are technical only."
        if person_held.any() and not statuses.eq("calibration_inference_withheld").any():
            status, reason = "withheld", "stat.reliability.person_measure_withheld"
            note += (f" {int(person_held.sum())} person measure(s) have no finite JMLE or are withheld; "
                     "retained calculations are technical only. No persons were dropped from the calculation.")
        counts = statuses.value_counts().to_dict()
        status_summary = "; ".join(f"{key}={value}" for key, value in sorted(counts.items()))
        basis = ("Empirical reliability for EAP measures: Var(EAP) / (Var(EAP) + mean posterior SD^2). "
                 "Conditional posterior SD is not the sampling SE of EAP; this is not directly comparable to JMLE-style separation reliability."
                 if eap else "Reliability uses observed estimate variance minus mean squared SE; interpretation follows the stated uncertainty basis.")
        fit_means = {}
        for source, target in (("Infit", "MeanInfit"), ("Outfit", "MeanOutfit")):
            values = pd.to_numeric(frame.get(source, pd.Series(dtype=float)), errors="coerce")
            fit_means[target] = values[np.isfinite(values)].mean()
        rows.append(dict(Facet=facet, Levels=len(frame), ValidLevels=int((valid & ~person_held).sum()),
            SD=np.sqrt(mv) if np.isfinite(mv) else np.nan,
            RMSE=np.sqrt(ev) if np.isfinite(ev) else np.nan,
            Separation=separation, Strata=strata, Reliability=reliability, **fit_means,
            ReliabilityNote=f"{note} SE status: {status_summary}".strip(),
            SE_StatusSummary=status_summary, ReliabilityBasis=basis,
            ReliabilityBasisCode=("eap_posterior_sd" if eap else "se_variance_correction") if known_basis else "unverified",
            ReliabilityStatus=status, ReliabilityReasonCode=reason,
            ReliabilityCalculationVersion=CALCULATION_VERSION))
    return pd.DataFrame(rows)


def qualify_saved_reliability(saved, measures=None, *, calibration_hold=False, hold_reason=""):
    """Revalidate saved values without replacing them with new calculations."""
    if not isinstance(saved, pd.DataFrame) or saved.empty:
        return pd.DataFrame()
    out = saved.copy()
    if "ReliabilityStatus" in out.columns and "SavedReliabilityStatus" not in out.columns:
        out["SavedReliabilityStatus"] = out["ReliabilityStatus"]
    if "ReliabilityReasonCode" in out.columns and "SavedReliabilityReasonCode" not in out.columns:
        out["SavedReliabilityReasonCode"] = out["ReliabilityReasonCode"]
    current = calculate_reliability(measures) if isinstance(measures, pd.DataFrame) and {"Facet", "Estimate", "SE"}.issubset(measures) else pd.DataFrame()
    numeric = ("Levels", "SD", *PRECISION_COLUMNS)
    for index, row in out.iterrows():
        status, reason = "unverified", "stat.reliability.saved_calculation_unverified"
        match = current.loc[current.Facet.eq(row.get("Facet"))] if not current.empty else pd.DataFrame()
        if not match.empty:
            checked = match.iloc[0]
            status, reason = checked.ReliabilityStatus, checked.ReliabilityReasonCode
            if status == "diagnostic":
                actual = pd.to_numeric(row.reindex(numeric), errors="coerce").to_numpy(dtype=float)
                expected = pd.to_numeric(checked.reindex(numeric), errors="coerce").to_numpy(dtype=float)
                if not np.isclose(actual, expected, rtol=1e-10, atol=1e-12, equal_nan=True).all():
                    status, reason = "unverified", "stat.reliability.saved_values_mismatch"
        elif not isinstance(measures, pd.DataFrame) and row.get("ReliabilityCalculationVersion") == CALCULATION_VERSION:
            status = row.get("ReliabilityStatus", status)
            reason = row.get("ReliabilityReasonCode", reason)
            if status == "diagnostic":
                values = pd.to_numeric(row.reindex(numeric), errors="coerce").to_numpy(dtype=float)
                if (not np.isfinite(values).all() or values[0] < 2 or values[1] < 0 or values[2] <= 0
                    or values[3] < 0 or values[4] < 1 / 3
                    or row.get("ReliabilityBasisCode") not in {"eap_posterior_sd", "se_variance_correction"}
                    or not 0 <= row.get("Reliability", np.nan) <= 1):
                    status, reason = "unverified", "stat.reliability.saved_calculation_unverified"
        if calibration_hold:
            status, reason = "withheld", hold_reason or "stat.reliability.calibration_inference_withheld"
        out.loc[index, "ReliabilityStatus"] = status
        out.loc[index, "ReliabilityReasonCode"] = reason
    return out


def calculate_reliability_series(
    measures: pd.DataFrame | None,
    result: dict | None = None,
    *,
    obs_df: pd.DataFrame | None = None,
    saved_reliability: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Compute the 4 reliability series (Model/Real SE × Sample/Population) with explicit counts and saved provenance."""
    from mfrm_app import engine

    result = result or {}
    run_state = output_qualification.run_qualification_state(result)
    block_reason = str(run_state["screening_block_reason"])
    ident_state = str(run_state["identifiability_state"])
    method = str((result.get("config") or {}).get("method", "")).upper() if isinstance(result.get("config"), dict) else ""

    saved_lookup: dict[str, pd.Series] = {}
    if isinstance(saved_reliability, pd.DataFrame) and not saved_reliability.empty and "Facet" in saved_reliability.columns:
        for _, srow in saved_reliability.iterrows():
            saved_lookup[str(srow["Facet"])] = srow

    if not isinstance(measures, pd.DataFrame) or measures.empty:
        if not saved_lookup:
            return pd.DataFrame()
        qualified_saved = qualify_saved_reliability(
            saved_reliability,
            None,
            calibration_hold=bool(block_reason),
            hold_reason=block_reason,
        )
        fallback_rows: list[dict[str, object]] = []
        for _, srow in qualified_saved.iterrows():
            facet = str(srow["Facet"])
            s_status = str(srow.get("SavedReliabilityStatus", srow.get("ReliabilityStatus", "")))
            s_reason = str(srow.get("SavedReliabilityReasonCode", srow.get("ReliabilityReasonCode", "")))
            cur_status = str(srow.get("ReliabilityStatus", "unverified"))
            cur_reason = str(srow.get("ReliabilityReasonCode", "stat.reliability.saved_calculation_unverified"))
            levels = int(srow["Levels"]) if pd.notna(srow.get("Levels")) else 0
            valid_n = int(srow["ValidLevels"]) if pd.notna(srow.get("ValidLevels")) else levels
            for se_mode, se_col in (("model", "ModelSE"), ("fit_adjusted", "RealSE")):
                for dist_basis in ("sample", "population"):
                    is_primary = (se_mode == "model" and dist_basis == "sample")
                    usable = is_primary and cur_status == "diagnostic" and not block_reason
                    fallback_rows.append({
                        "Facet": facet,
                        "SEMode": se_mode,
                        "SEColumn": se_col,
                        "DistributionBasis": dist_basis,
                        "Levels": levels,
                        "EstimateAvailable": valid_n if is_primary else np.nan,
                        "SEAvailable": valid_n if usable else (0 if not is_primary else np.nan),
                        "ExcludedLevels": max(levels - valid_n, 0) if is_primary else np.nan,
                        "AggregationTarget": (
                            f"saved_paired_levels ({valid_n}/{levels})"
                            if is_primary
                            else "unavailable (source measures table not retained)"
                        ),
                        "SD": float(srow["SD"]) if is_primary and pd.notna(srow.get("SD")) else np.nan,
                        "RMSE": float(srow["RMSE"]) if is_primary and pd.notna(srow.get("RMSE")) else np.nan,
                        "TrueSD": np.nan,
                        "Separation": float(srow["Separation"]) if is_primary and pd.notna(srow.get("Separation")) else np.nan,
                        "Strata": float(srow["Strata"]) if is_primary and pd.notna(srow.get("Strata")) else np.nan,
                        "Reliability": float(srow["Reliability"]) if is_primary and pd.notna(srow.get("Reliability")) else np.nan,
                        "ReliabilityBasisCode": str(srow.get("ReliabilityBasisCode", "unverified")) if is_primary else "unavailable",
                        "Use": "screening_only" if usable else "withheld",
                        "SupportsFormalInference": False,
                        "ReliabilityStatus": "diagnostic" if usable else ("unassessable" if not is_primary else cur_status),
                        "StatusReason": (
                            cur_reason if is_primary else "stat.reliability.source_measures_unavailable"
                        ),
                        "SavedReliability": (
                            float(srow["Reliability"])
                            if is_primary and s_status in {"diagnostic", ""} and not usable and pd.notna(srow.get("Reliability"))
                            else np.nan
                        ),
                        "SavedReliabilityStatus": s_status,
                        "SavedReliabilityReasonCode": s_reason,
                        "CalculationVersion": SERIES_CALCULATION_VERSION,
                    })
        return pd.DataFrame(fallback_rows)

    q_measures = engine.qualify_measures_for_reporting(measures, result, obs_df=obs_df)
    config_facets = list((result.get("config") or {}).get("facet_names") or []) if isinstance(result.get("config"), dict) else []
    observed_facets = [str(f) for f in q_measures["Facet"].dropna().unique()]
    ordered_facets = list(dict.fromkeys([*(["Person"] if "Person" in observed_facets else []), *[f for f in config_facets if f in observed_facets], *observed_facets]))

    rows: list[dict[str, object]] = []
    for facet in ordered_facets:
        sub = q_measures.loc[q_measures["Facet"].astype(str).eq(facet)].copy()
        if sub.empty:
            continue
        n_total = int(len(sub))
        fixed_mask = sub.get("Fixed", pd.Series(False, index=sub.index)).eq(True).fillna(False)
        n_fixed = int(fixed_mask.sum())

        endpoint = sub.get("ExtremeScoreDirection", pd.Series("", index=sub.index)).astype(str).isin({"all_minimum", "all_maximum"})
        extreme_mask = (
            sub["Facet"].astype(str).eq("Person")
            & endpoint
            & sub.get("FiniteJMLEEstimate", pd.Series(True, index=sub.index)).eq(False).fillna(False)
        )
        n_extreme = int(extreme_mask.sum())

        other_person_held = (
            sub["Facet"].astype(str).eq("Person")
            & (
                sub.get("FiniteJMLEEstimate", pd.Series(True, index=sub.index)).eq(False).fillna(False)
                | sub.get("PersonInferenceReady", pd.Series(True, index=sub.index)).eq(False).fillna(False)
            )
            & ~extreme_mask
        )

        unbounded_facet = sub.get("UnboundedFacetEstimate", pd.Series(False, index=sub.index)).eq(True).fillna(False)
        n_unbounded_facet = int((unbounded_facet & ~fixed_mask & ~extreme_mask).sum())

        target_mask = ~fixed_mask & ~extreme_mask
        n_target = int(target_mask.sum())
        est = pd.to_numeric(sub["Estimate"], errors="coerce")
        est_ok = target_mask & np.isfinite(est) & ~other_person_held & ~unbounded_facet
        n_est_avail = int(est_ok.sum())
        n_excluded = int(n_total - n_est_avail)

        methods = sub.loc[target_mask, "SE_Method"].fillna("").astype(str).str.lower() if "SE_Method" in sub.columns else pd.Series("", index=sub.loc[target_mask].index)
        eff_statuses = sub.loc[target_mask, "EffectiveSE_Status"].fillna("").astype(str) if "EffectiveSE_Status" in sub.columns else pd.Series("", index=sub.loc[target_mask].index)
        posterior = methods.str.contains("latent-ability posterior sd", regex=False)
        se_based = methods.str.startswith(("conditional information approximation", "mml observed-information delta"))
        legacy_post = eff_statuses.eq("posterior_eap_sd")
        eap = bool(n_target > 0 and (posterior.all() or (methods.eq("").all() and legacy_post.all())))
        mixed = bool(n_target > 0 and ((posterior.any() and not posterior.all()) or (methods.eq("").all() and legacy_post.any() and not legacy_post.all())))
        known_basis = bool(eap or (n_target > 0 and se_based.all()))
        basis_code = ("eap_posterior_sd" if eap else "se_variance_correction") if known_basis else "unverified"

        srow = saved_lookup.get(facet)
        raw_saved_rel = float(srow["Reliability"]) if srow is not None and pd.notna(srow.get("Reliability")) else np.nan
        saved_status = str(srow.get("ReliabilityStatus", "")) if srow is not None else ""
        saved_reason = str(srow.get("ReliabilityReasonCode", "")) if srow is not None else ""
        # Check saved value integrity when a 1-row saved record was supplied without block/extreme changes.
        saved_mismatch = False
        if srow is not None and saved_status == "diagnostic" and n_extreme == 0 and n_fixed == 0:
            cur_legacy = calculate_reliability(measures.loc[measures["Facet"].astype(str).eq(facet)])
            if not cur_legacy.empty and cur_legacy.iloc[0]["ReliabilityStatus"] == "diagnostic":
                actual_arr = pd.to_numeric(srow.reindex(("Levels", "SD", *PRECISION_COLUMNS)), errors="coerce").to_numpy(dtype=float)
                expect_arr = pd.to_numeric(cur_legacy.iloc[0].reindex(("Levels", "SD", *PRECISION_COLUMNS)), errors="coerce").to_numpy(dtype=float)
                if not np.isclose(actual_arr, expect_arr, rtol=1e-10, atol=1e-12, equal_nan=True).all():
                    saved_mismatch = True

        if n_fixed > 0:
            agg_target = f"unavailable ({n_fixed}/{n_total} fixed by constraints)"
        elif n_extreme > 0:
            agg_target = f"finite_non_extreme_levels ({n_est_avail}/{n_total}; {n_extreme} extreme boundary excluded)"
        elif n_unbounded_facet > 0:
            agg_target = f"finite_supported_levels ({n_est_avail}/{n_total}; {n_unbounded_facet} unbounded boundary excluded)"
        else:
            agg_target = f"all_levels ({n_est_avail}/{n_total})"

        infit_target = pd.to_numeric(sub.loc[target_mask, "Infit"], errors="coerce") if "Infit" in sub.columns else pd.Series(np.nan, index=sub.loc[target_mask].index)

        for se_mode, se_col in (("model", "ModelSE"), ("fit_adjusted", "RealSE")):
            se_series = pd.to_numeric(sub.get(se_col, pd.Series(np.nan, index=sub.index)), errors="coerce")
            se_ok = est_ok & np.isfinite(se_series) & (se_series > 0) & sub["SEUse"].eq("screening_only")
            n_se_avail = int(se_ok.sum())

            for dist_basis, ddof_val in (("sample", 1), ("population", 0)):
                sd_val = rmse_val = true_sd_val = sep_val = strata_val = rel_val = np.nan
                if n_est_avail >= 2 and n_est_avail == n_target and n_fixed == 0:
                    with np.errstate(over="ignore", invalid="ignore"):
                        mv = float(est.loc[est_ok].var(ddof=ddof_val))
                        sd_val = float(np.sqrt(mv)) if np.isfinite(mv) and mv >= 0 else np.nan
                else:
                    mv = np.nan

                if block_reason:
                    status, reason, use_label = "withheld", block_reason, "withheld"
                elif saved_mismatch:
                    status, reason, use_label = "unverified", "stat.reliability.saved_values_mismatch", "withheld"
                elif n_fixed > 0:
                    status, reason, use_label = "unassessable", "stat.reliability.fixed_levels_present", "withheld"
                elif other_person_held.any() or n_est_avail < n_target or n_se_avail < n_est_avail:
                    if (
                        se_mode == "fit_adjusted"
                        and n_est_avail == n_target
                        and "Infit" in sub.columns
                        and (~np.isfinite(infit_target)).any()
                    ):
                        status, reason, use_label = "unassessable", "stat.reliability.fit_unavailable_for_real_se", "withheld"
                    elif (sub.loc[target_mask, "StatusReason"] == "stat.stat_003.covariance_rank_deficient").any():
                        status, reason, use_label = "withheld", "stat.stat_003.covariance_rank_deficient", "withheld"
                    else:
                        status, reason, use_label = "unassessable", "stat.reliability.invalid_measure_pairs", "withheld"
                elif n_est_avail < 2:
                    status, reason, use_label = "unassessable", "stat.reliability.insufficient_levels", "withheld"
                elif mixed:
                    status, reason, use_label = "unassessable", "stat.reliability.mixed_uncertainty_basis", "withheld"
                else:
                    ev = float(np.mean(se_series.loc[se_ok] ** 2))
                    rmse_val = float(np.sqrt(ev)) if np.isfinite(ev) and ev > 0 else np.nan
                    if not np.isfinite(mv) or not np.isfinite(ev) or ev <= 0:
                        status, reason, use_label = "unassessable", "stat.reliability.invalid_variance", "withheld"
                    elif mv == 0 and not eap:
                        status, reason, use_label = "unassessable", "stat.reliability.zero_observed_variance", "withheld"
                    elif not known_basis:
                        status, reason, use_label = "unverified", "stat.reliability.uncertainty_basis_unverified", "withheld"
                    else:
                        signal = mv if eap else max(mv - ev, 0.0)
                        true_sd_val = float(np.sqrt(signal))
                        sep_val = float(true_sd_val / rmse_val)
                        rel_val = float(mv / (mv + ev) if eap else signal / mv)
                        strata_val = float((4.0 * sep_val + 1.0) / 3.0)
                        status, use_label = "diagnostic", "screening_only"
                        if eap:
                            reason = "stat.mml_fixed_sd.eap_empirical_reliability_screening_only"
                        elif method == "MML":
                            reason = "stat.mml_fixed_sd.se_reliability_screening_only"
                        elif method == "JMLE" and ident_state == "audit_size_limit":
                            reason = output_qualification.JMLE_AUDIT_SIZE_LIMIT_HOLD_REASON
                        else:
                            reason = "stat.reliability.descriptive_screening_only"

                show_saved_rel = (
                    se_mode == "model"
                    and dist_basis == "sample"
                    and np.isfinite(raw_saved_rel)
                    and saved_status in {"diagnostic", ""}
                    and (status != "diagnostic" or not np.isclose(raw_saved_rel, rel_val, atol=1e-6))
                )
                rows.append({
                    "Facet": facet,
                    "SEMode": se_mode,
                    "SEColumn": se_col,
                    "DistributionBasis": dist_basis,
                    "Levels": n_total,
                    "EstimateAvailable": n_est_avail,
                    "SEAvailable": n_se_avail,
                    "ExcludedLevels": n_excluded,
                    "AggregationTarget": agg_target,
                    "SD": sd_val,
                    "RMSE": rmse_val,
                    "TrueSD": true_sd_val,
                    "Separation": sep_val,
                    "Strata": strata_val,
                    "Reliability": rel_val,
                    "ReliabilityBasisCode": basis_code,
                    "Use": use_label,
                    "SupportsFormalInference": False,
                    "ReliabilityStatus": status,
                    "StatusReason": reason,
                    "SavedReliability": raw_saved_rel if show_saved_rel else np.nan,
                    "SavedReliabilityStatus": saved_status,
                    "SavedReliabilityReasonCode": saved_reason,
                    "CalculationVersion": SERIES_CALCULATION_VERSION,
                })
    return pd.DataFrame(rows)
