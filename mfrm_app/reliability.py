"""Descriptive reliability with complete paired inputs and explicit provenance."""
import numpy as np
import pandas as pd


CALCULATION_VERSION = "paired-reliability-v2"
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
