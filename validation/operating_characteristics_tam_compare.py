#!/usr/bin/env python3
"""Normalize TAM MML evidence against Python, sirt, and structural gates."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


Q21 = "TAM_MML_Q21_SENSITIVITY"
Q61 = "TAM_MML_Q61_PRIMARY"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _pair_modes(
    frame: pd.DataFrame,
    *,
    keys: list[str],
    estimate: str,
    parameter_type: str | pd.Series,
) -> pd.DataFrame:
    q21 = frame.loc[
        frame["Mode"].eq(Q21), keys + [estimate, "IncludedInSummary"]
    ].rename(
        columns={estimate: "EstimateQ21", "IncludedInSummary": "IncludedQ21"}
    )
    q61 = frame.loc[
        frame["Mode"].eq(Q61), keys + [estimate, "IncludedInSummary"]
    ].rename(
        columns={estimate: "EstimateQ61", "IncludedInSummary": "IncludedQ61"}
    )
    out = q21.merge(q61, on=keys, how="outer", validate="one_to_one", indicator=True)
    out["ParameterType"] = parameter_type if isinstance(parameter_type, str) else ""
    if not isinstance(parameter_type, str):
        # Facet type is stable across modes and can be reconstructed from the key.
        out["ParameterType"] = "Facet: " + out["Facet"].astype(str)
    out["BothPresent"] = out["_merge"].eq("both")
    out["BothIncluded"] = (
        out["BothPresent"]
        & out["IncludedQ21"].fillna(False).astype(bool)
        & out["IncludedQ61"].fillna(False).astype(bool)
    )
    out["DifferenceQ61MinusQ21"] = out["EstimateQ61"] - out["EstimateQ21"]
    out["AbsDifference"] = out["DifferenceQ61MinusQ21"].abs()
    return out.drop(columns="_merge")


def build_quadrature_evidence(
    facets: pd.DataFrame,
    steps: pd.DataFrame,
    persons: pd.DataFrame,
    surfaces: pd.DataFrame,
    runs: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    blocks = [
        _pair_modes(
            facets,
            keys=["RunId", "Facet", "Level"],
            estimate="EstimateAligned",
            parameter_type=facets["Facet"],
        ),
        _pair_modes(
            steps,
            keys=["RunId", "Step"],
            estimate="Estimate",
            parameter_type="Common step",
        ),
        _pair_modes(
            persons,
            keys=["RunId", "Person"],
            estimate="Estimate",
            parameter_type="Person EAP",
        ),
        _pair_modes(
            surfaces,
            keys=["RunId", "GeneralizedItem", "Category"],
            estimate="CumulativeDifficulty",
            parameter_type="Cumulative response-surface difficulty",
        ),
    ]
    metadata = runs.loc[
        runs["Mode"].eq(Q61),
        ["RunId", "ConditionId", "Design", "TruthBias", "Replicate", "Seed"],
    ]
    pairs = pd.concat(blocks, ignore_index=True, sort=False).merge(
        metadata, on="RunId", how="left", validate="many_to_one"
    )
    rows: list[dict[str, object]] = []
    for (parameter_type, design), group in pairs.groupby(
        ["ParameterType", "Design"], sort=False
    ):
        for scope, mask in [
            ("all_jointly_returned", group["BothPresent"]),
            ("jointly_analysis_eligible", group["BothIncluded"]),
        ]:
            scoped = group.loc[mask]
            run_max = scoped.groupby("RunId")["AbsDifference"].max()
            rows.append(
                {
                    "ParameterType": parameter_type,
                    "Design": design,
                    "Scope": scope,
                    "RunIds": int(scoped["RunId"].nunique()),
                    "Pairs": int(len(scoped)),
                    "MeanAbsDifference": float(scoped["AbsDifference"].mean()) if len(scoped) else np.nan,
                    "MaxAbsDifference": float(scoped["AbsDifference"].max()) if len(scoped) else np.nan,
                    "RunsMaxDifferenceLE1e1": int((run_max <= 1e-1).sum()),
                    "RunsMaxDifferenceLE1e2": int((run_max <= 1e-2).sum()),
                    "RunsMaxDifferenceLE1e3": int((run_max <= 1e-3).sum()),
                    "ThresholdDenominator": int(len(run_max)),
                    "ThresholdBasis": "unrounded retained floating-point values",
                }
            )

    run_pairs = runs.loc[runs["Mode"].eq(Q21)].merge(
        runs.loc[runs["Mode"].eq(Q61)],
        on="RunId",
        suffixes=("Q21", "Q61"),
        validate="one_to_one",
    )
    for metric in ["PopulationSD", "LogLik", "EAPReliability"]:
        run_pairs[f"{metric}DifferenceQ61MinusQ21"] = (
            run_pairs[f"{metric}Q61"] - run_pairs[f"{metric}Q21"]
        )
        run_pairs[f"Abs{metric}Difference"] = run_pairs[
            f"{metric}DifferenceQ61MinusQ21"
        ].abs()
    run_rows = []
    for design, group in run_pairs.groupby("DesignQ61", sort=False):
        eligible = group["AnalysisEligibleQ21"].astype(bool) & group[
            "AnalysisEligibleQ61"
        ].astype(bool)
        for scope, mask in [
            ("all_jointly_returned", pd.Series(True, index=group.index)),
            ("jointly_analysis_eligible", eligible),
        ]:
            scoped = group.loc[mask]
            run_rows.append(
                {
                    "ParameterType": "Run-level population/loglik",
                    "Design": design,
                    "Scope": scope,
                    "RunIds": int(len(scoped)),
                    "Pairs": int(len(scoped)),
                    "MeanAbsDifference": float(scoped["AbsPopulationSDDifference"].mean()) if len(scoped) else np.nan,
                    "MaxAbsDifference": float(scoped["AbsPopulationSDDifference"].max()) if len(scoped) else np.nan,
                    "RunsMaxDifferenceLE1e1": int((scoped["AbsPopulationSDDifference"] <= 1e-1).sum()),
                    "RunsMaxDifferenceLE1e2": int((scoped["AbsPopulationSDDifference"] <= 1e-2).sum()),
                    "RunsMaxDifferenceLE1e3": int((scoped["AbsPopulationSDDifference"] <= 1e-3).sum()),
                    "ThresholdDenominator": int(len(scoped)),
                    "ThresholdBasis": "unrounded retained floating-point values",
                    "MaxAbsLogLikDifference": float(scoped["AbsLogLikDifference"].max()) if len(scoped) else np.nan,
                    "MaxAbsEAPReliabilityDifference": float(scoped["AbsEAPReliabilityDifference"].max()) if len(scoped) else np.nan,
                }
            )
    summary = pd.concat([pd.DataFrame(rows), pd.DataFrame(run_rows)], ignore_index=True)
    return pairs, summary, run_pairs


def build_python_tam_agreement(
    python_recovery: pd.DataFrame,
    tam_recovery: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tam = tam_recovery.loc[tam_recovery["Mode"].eq(Q61)].copy()
    keys = ["RunId", "Facet", "Level", "ComparisonScale"]
    left = python_recovery[keys + [
        "ConditionId", "Design", "TruthBias", "Replicate", "Seed",
        "TruthAligned", "EstimateAligned", "SE", "ErrorAligned",
        "IncludedInSummary",
    ]].rename(
        columns={
            "EstimateAligned": "EstimatePythonJMLE",
            "SE": "SEPythonJMLE",
            "ErrorAligned": "ErrorPythonJMLE",
            "IncludedInSummary": "IncludedPythonJMLE",
        }
    )
    right = tam[keys + [
        "EstimateAligned", "SE", "ErrorAligned", "IncludedInSummary",
        "Anchored", "DerivedConstraint", "CoverageEligible",
    ]].rename(
        columns={
            "EstimateAligned": "EstimateTAMMML",
            "SE": "SETAMMML",
            "ErrorAligned": "ErrorTAMMML",
            "IncludedInSummary": "IncludedTAMMML",
        }
    )
    pairs = left.merge(right, on=keys, how="outer", validate="one_to_one", indicator=True)
    pairs["BothPresent"] = pairs["_merge"].eq("both")
    pairs["BothIncluded"] = (
        pairs["BothPresent"]
        & pairs["IncludedPythonJMLE"].fillna(False).astype(bool)
        & pairs["IncludedTAMMML"].fillna(False).astype(bool)
    )
    pairs["DifferencePythonJMLEMinusTAMMML"] = (
        pairs["EstimatePythonJMLE"] - pairs["EstimateTAMMML"]
    )
    pairs["AbsDifference"] = pairs["DifferencePythonJMLEMinusTAMMML"].abs()
    pairs["InterpretationBoundary"] = (
        "same additive RSM facet surface but JMLE versus normal-population MML Person treatment; not parity"
    )
    pairs = pairs.drop(columns="_merge")

    rows = []
    for (design, facet), group in pairs.groupby(["Design", "Facet"], sort=False):
        included = group.loc[group["BothIncluded"]]
        rows.append(
            {
                "Design": design,
                "Facet": facet,
                "RunIds": int(group["RunId"].nunique()),
                "BothIncludedRunIds": int(included["RunId"].nunique()),
                "Pairs": int(len(group)),
                "BothIncludedPairs": int(len(included)),
                "MeanAbsDifference": float(included["AbsDifference"].mean()) if len(included) else np.nan,
                "MaxAbsDifference": float(included["AbsDifference"].max()) if len(included) else np.nan,
                "PythonJMLERMSE": float(np.sqrt(np.mean(np.square(included["ErrorPythonJMLE"])))) if len(included) else np.nan,
                "TAMMMLRMSE": float(np.sqrt(np.mean(np.square(included["ErrorTAMMML"])))) if len(included) else np.nan,
                "InterpretationBoundary": "different Person treatment; design-specific sensitivity, not parity",
            }
        )
    return pairs, pd.DataFrame(rows)


def build_tam_sirt_rater_agreement(
    tam_recovery: pd.DataFrame,
    sirt_recovery: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tam = tam_recovery.loc[
        tam_recovery["Mode"].eq(Q61) & tam_recovery["Facet"].eq("Rater")
    ]
    sirt = sirt_recovery.loc[sirt_recovery["Mode"].eq("SIRT_MML_Q61_PRIMARY")]
    keys = ["RunId", "Level", "ComparisonScale"]
    left = tam[keys + [
        "ConditionId", "Design", "TruthBias", "Replicate", "Seed",
        "EstimateAligned", "IncludedInSummary",
    ]].rename(
        columns={
            "EstimateAligned": "EstimateTAMMML",
            "IncludedInSummary": "IncludedTAMMML",
        }
    )
    right = sirt[keys + ["EstimateAligned", "IncludedInSummary"]].rename(
        columns={
            "EstimateAligned": "EstimateSirtMML",
            "IncludedInSummary": "IncludedSirtMML",
        }
    )
    pairs = left.merge(right, on=keys, how="outer", validate="one_to_one", indicator=True)
    pairs["BothPresent"] = pairs["_merge"].eq("both")
    pairs["BothIncluded"] = (
        pairs["BothPresent"]
        & pairs["IncludedTAMMML"].fillna(False).astype(bool)
        & pairs["IncludedSirtMML"].fillna(False).astype(bool)
    )
    pairs["DifferenceTAMMinusSirt"] = pairs["EstimateTAMMML"] - pairs["EstimateSirtMML"]
    pairs["AbsDifference"] = pairs["DifferenceTAMMinusSirt"].abs()
    pairs["InterpretationBoundary"] = (
        "both MML, but TAM uses additive RSM facets and sirt uses virtual-item-specific PCM thresholds"
    )
    pairs = pairs.drop(columns="_merge")
    rows = []
    for design, group in pairs.groupby("Design", sort=False):
        included = group.loc[group["BothIncluded"]]
        rows.append(
            {
                "Design": design,
                "RunIds": int(group["RunId"].nunique()),
                "BothIncludedRunIds": int(included["RunId"].nunique()),
                "Pairs": int(len(group)),
                "BothIncludedPairs": int(len(included)),
                "MeanAbsDifference": float(included["AbsDifference"].mean()) if len(included) else np.nan,
                "MaxAbsDifference": float(included["AbsDifference"].max()) if len(included) else np.nan,
                "InterpretationBoundary": "MML comparison with different item-threshold structure and constraint mechanics",
            }
        )
    return pairs, pd.DataFrame(rows)


def build_anchor_sensitivity(
    tam_recovery: pd.DataFrame,
    tam_runs: pd.DataFrame,
) -> pd.DataFrame:
    primary = tam_recovery.loc[tam_recovery["Mode"].eq(Q61)]
    clean = primary.loc[primary["Design"].eq("balanced_large_anchors")]
    drift = primary.loc[primary["Design"].eq("anchor_drift")]
    keys = ["TruthBias", "Replicate", "Seed", "Facet", "Level"]
    pairs = clean.merge(drift, on=keys, suffixes=("Clean", "Drift"), validate="one_to_one")
    pairs["EstimateShiftDriftMinusClean"] = pairs["EstimateDrift"] - pairs["EstimateClean"]
    pairs["ShiftErrorFromInjected0p25"] = pairs["EstimateShiftDriftMinusClean"] - 0.25
    pairs["ExpectedConstraintRole"] = np.select(
        [
            pairs["Facet"].eq("Rater") & pairs["AnchoredClean"].astype(bool),
            pairs["Facet"].eq("Rater") & ~pairs["AnchoredClean"].astype(bool),
        ],
        ["fixed anchor receives +0.25", "unanchored Rater compensates under sum-zero constraint"],
        default="downstream facet sensitivity",
    )
    run_clean = tam_runs.loc[
        tam_runs["Mode"].eq(Q61) & tam_runs["Design"].eq("balanced_large_anchors")
    ]
    run_drift = tam_runs.loc[
        tam_runs["Mode"].eq(Q61) & tam_runs["Design"].eq("anchor_drift")
    ]
    run_pairs = run_clean.merge(
        run_drift,
        on=["TruthBias", "Replicate", "Seed"],
        suffixes=("Clean", "Drift"),
        validate="one_to_one",
    )
    run_pairs["DevianceShiftDriftMinusClean"] = run_pairs["DevianceDrift"] - run_pairs["DevianceClean"]
    run_pairs["PopulationSDShiftDriftMinusClean"] = run_pairs["PopulationSDDrift"] - run_pairs["PopulationSDClean"]
    return pairs.merge(
        run_pairs[[
            "TruthBias", "Replicate", "Seed", "DevianceShiftDriftMinusClean",
            "PopulationSDShiftDriftMinusClean",
        ]],
        on=["TruthBias", "Replicate", "Seed"],
        validate="many_to_one",
    )


def build_omitted_bias_sensitivity(tam_recovery: pd.DataFrame) -> pd.DataFrame:
    primary = tam_recovery.loc[tam_recovery["Mode"].eq(Q61)]
    null = primary.loc[primary["TruthBias"].eq(0)].rename(
        columns={
            "RunId": "RunIdNull", "EstimateAligned": "EstimateNull",
            "IncludedInSummary": "IncludedNull",
        }
    )
    alternative = primary.loc[primary["TruthBias"].ne(0)].rename(
        columns={
            "RunId": "RunIdAlternative", "TruthBias": "TruthBiasAlternative",
            "EstimateAligned": "EstimateAlternative",
            "IncludedInSummary": "IncludedAlternative",
        }
    )
    keys = ["Design", "Replicate", "Seed", "Facet", "Level"]
    pairs = null[keys + ["RunIdNull", "EstimateNull", "IncludedNull"]].merge(
        alternative[keys + [
            "RunIdAlternative", "TruthBiasAlternative", "EstimateAlternative",
            "IncludedAlternative",
        ]],
        on=keys,
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    pairs["PairEligible"] = (
        pairs["_merge"].eq("both")
        & pairs["IncludedNull"].fillna(False).astype(bool)
        & pairs["IncludedAlternative"].fillna(False).astype(bool)
    )
    pairs["EstimateChangeAlternativeMinusNull"] = (
        pairs["EstimateAlternative"] - pairs["EstimateNull"]
    )
    pairs["AbsEstimateChange"] = pairs["EstimateChangeAlternativeMinusNull"].abs()
    pairs["InterpretationBoundary"] = (
        "omitted Rater x Task interaction leakage into additive TAM facets; not bias detection"
    )
    return pairs.drop(columns="_merge")


def build_variance_boundary_audit(tam_runs: pd.DataFrame) -> pd.DataFrame:
    q21 = tam_runs.loc[tam_runs["Mode"].eq(Q21)]
    q61 = tam_runs.loc[tam_runs["Mode"].eq(Q61)]
    out = q21.merge(q61, on="RunId", suffixes=("Q21", "Q61"), validate="one_to_one")
    for mode in ["Q21", "Q61"]:
        out[f"StoredOffsetFromLowerBound{mode}"] = (
            out[f"PopulationVariance{mode}"] - out[f"MinimumPopulationVariance{mode}"]
        )
        out[f"NaiveExactBoundaryDecision{mode}"] = (
            out[f"PopulationVariance{mode}"] <= out[f"MinimumPopulationVariance{mode}"]
        )
        out[f"TolerantBoundaryDecision{mode}"] = out[
            f"PopulationVarianceAtLowerBound{mode}"
        ].astype(bool)
        out[f"NaiveVsTolerantMismatch{mode}"] = (
            out[f"NaiveExactBoundaryDecision{mode}"]
            != out[f"TolerantBoundaryDecision{mode}"]
        )
    out["InterpretationBoundary"] = (
        "variance-boundary decisions use the retained 1e-8 numerical band; display rounding is not authoritative"
    )
    return out


def build_progress_threshold_summary(tam_runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for mode, group in tam_runs.groupby("Mode", sort=False):
        parsed = group.loc[group["TerminalProgressParsed"]]
        for threshold in [1e-4, 1e-5, 1e-6]:
            rows.append(
                {
                    "Mode": mode,
                    "Threshold": threshold,
                    "RunNumerator": int((parsed["TerminalMaxParameterChange"] <= threshold).sum()),
                    "RunDenominator": int(len(parsed)),
                    "EvidenceBasis": (
                        "parsed TAM progress text; primary 1e-5 item/deviance decision is also enforced by source-loop exit"
                    ),
                    "PrecisionBoundary": (
                        "progress values are formatted text and are sensitivity evidence, not higher-precision gradients"
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_sparse_identification(
    tam_runs: pd.DataFrame,
    sirt_runs: pd.DataFrame,
    python_runs: pd.DataFrame,
    cmle_runs: pd.DataFrame,
    mfrmr_runs: pd.DataFrame,
) -> pd.DataFrame:
    tam = tam_runs.loc[tam_runs["Mode"].eq(Q61) & tam_runs["Design"].eq("sparse_missing")]
    sirt = sirt_runs.loc[
        sirt_runs["Mode"].eq("SIRT_MML_Q61_PRIMARY")
        & sirt_runs["Design"].eq("sparse_missing"),
        ["RunId", "AnalysisEligible"],
    ].rename(columns={"AnalysisEligible": "SirtMMLAnalysisEligible"})
    python = python_runs.loc[python_runs["Design"].eq("sparse_missing"), [
        "RunId", "FitReturned", "Converged", "InferenceReady",
    ]].rename(columns=lambda name: name if name == "RunId" else f"PythonJMLE{name}")
    cmle = cmle_runs.loc[cmle_runs["Design"].eq("sparse_missing"), [
        "RunId", "ConditionalDesignEligible", "ConditionalRank", "ConditionalNullity",
    ]].rename(columns=lambda name: name if name == "RunId" else f"PythonCMLE{name}")
    mfrmr = mfrmr_runs.loc[
        mfrmr_runs["Design"].eq("sparse_missing") & mfrmr_runs["Mode"].eq("MFRMR_JML_STRICT"),
        ["RunId", "FitReturned", "Converged", "InferenceReady", "FailureReason"],
    ].rename(columns=lambda name: name if name == "RunId" else f"Mfrmr{name}")
    out = tam.merge(sirt, on="RunId", validate="one_to_one").merge(
        python, on="RunId", validate="one_to_one"
    ).merge(cmle, on="RunId", validate="one_to_one").merge(
        mfrmr, on="RunId", validate="one_to_one"
    )
    out["IdentificationBoundary"] = (
        "TAM and sirt MML can use common-population assumptions across sparse facet rows; this is not observed-design rank"
    )
    return out


def build_first_read_summary(
    runs: pd.DataFrame,
    quadrature_pairs: pd.DataFrame,
    python_tam_summary: pd.DataFrame,
    tam_sirt_summary: pd.DataFrame,
    anchors: pd.DataFrame,
    variance_audit: pd.DataFrame,
    sparse: pd.DataFrame,
) -> pd.DataFrame:
    eligible_q = quadrature_pairs.loc[quadrature_pairs["BothIncluded"]]
    qmax = eligible_q.groupby("ParameterType")["AbsDifference"].max()
    balanced_rater = python_tam_summary.loc[
        python_tam_summary["Design"].eq("balanced_small")
        & python_tam_summary["Facet"].eq("Rater")
    ].iloc[0]
    balanced_mml = tam_sirt_summary.loc[
        tam_sirt_summary["Design"].eq("balanced_small")
    ].iloc[0]
    anchored = anchors.loc[anchors["ExpectedConstraintRole"].eq("fixed anchor receives +0.25")]
    compensating = anchors.loc[
        anchors["ExpectedConstraintRole"].eq(
            "unanchored Rater compensates under sum-zero constraint"
        )
    ]
    mismatch_count = int(
        variance_audit[["NaiveVsTolerantMismatchQ21", "NaiveVsTolerantMismatchQ61"]]
        .to_numpy(dtype=bool)
        .sum()
    )
    primary = runs.loc[runs["Mode"].eq(Q61)]
    return pd.DataFrame(
        [
            {
                "Order": 1,
                "Section": "Evidence scope",
                "Status": "Pilot only",
                "Headline": "Two-replicate TAM MML smoke sensitivity",
                "Evidence": "Estimator boundaries and failures are retained; operating characteristics are not estimated.",
                "UserAction": "Do not use these rows as performance or sample-size evidence.",
            },
            {
                "Order": 2,
                "Section": "Run accounting",
                "Status": "Review",
                "Headline": f"{int(runs['FitReturned'].sum())}/{len(runs)} returned; {int(runs['AnalysisEligible'].sum())}/{len(runs)} eligible",
                "Evidence": "Two q21 sparse fits reached the configured Person-variance lower boundary.",
                "UserAction": "Separate returned fits from inference-ready fits.",
            },
            {
                "Order": 3,
                "Section": "Estimator boundary",
                "Status": "Sensitivity only",
                "Headline": "TAM MML is not Python JMLE parity",
                "Evidence": f"Balanced-small Rater MAE={balanced_rater['MeanAbsDifference']:.6g}; maximum={balanced_rater['MaxAbsDifference']:.6g} logits.",
                "UserAction": "Attribute remaining differences to Person treatment before alleging implementation error.",
            },
            {
                "Order": 4,
                "Section": "MML structure",
                "Status": "Review",
                "Headline": "TAM and sirt are close only in some designs",
                "Evidence": f"Balanced-small TAM/sirt Rater MAE={balanced_mml['MeanAbsDifference']:.6g}; anchor constraints diverge under contamination.",
                "UserAction": "Check item-threshold structure and identification constraints for every comparison.",
            },
            {
                "Order": 5,
                "Section": "Quadrature",
                "Status": "Caution",
                "Headline": "Output sensitivity differs by parameter type",
                "Evidence": f"Eligible maxima: Rater={qmax['Facet: Rater']:.6g}; Person EAP={qmax['Person EAP']:.6g}; surface={qmax['Cumulative response-surface difficulty']:.6g}.",
                "UserAction": "Inspect q21/q61 sensitivity before accepting Person or surface output.",
            },
            {
                "Order": 6,
                "Section": "Floating boundary",
                "Status": "Caution",
                "Headline": f"{mismatch_count} naive-versus-tolerant variance-boundary mismatches",
                "Evidence": "TAM stored the lower-bound cases as 0.0010000001; exact <=0.001 alone would miss them.",
                "UserAction": "Use the retained numerical band and never classify from displayed rounding.",
            },
            {
                "Order": 7,
                "Section": "Anchor contamination",
                "Status": "Caution",
                "Headline": "Sum-zero constraints redistribute anchor drift",
                "Evidence": f"Anchors shifted {anchored['EstimateShiftDriftMinusClean'].mean():.6g}; unanchored Raters shifted {compensating['EstimateShiftDriftMinusClean'].mean():.6g} logits on average.",
                "UserAction": "Audit the constraint, not only whether fixed anchors were reproduced.",
            },
            {
                "Order": 8,
                "Section": "Sparse identification",
                "Status": "Caution",
                "Headline": f"Primary TAM eligible for {int(sparse['AnalysisEligible'].sum())}/{len(sparse)} sparse RunIds",
                "Evidence": "MML population assumptions can return estimates where exact CMLE and strict mfrmr fail rank gates.",
                "UserAction": "Do not relabel assumption-based identification as observed connectivity.",
            },
            {
                "Order": 9,
                "Section": "Public application surface",
                "Status": "Withheld",
                "Headline": "No TAM or cross-engine result button is enabled",
                "Evidence": f"Primary TAM eligible runs={int(primary['AnalysisEligible'].sum())}/{len(primary)} under a smoke-only profile.",
                "UserAction": "Complete study-depth, platform, and independent UX review gates first.",
            },
        ]
    ).assign(PublicSurfaceEnabled=False)


def plot_python_tam_agreement(pairs: pd.DataFrame, path: Path) -> None:
    data = pairs.loc[pairs["BothPresent"]].copy()
    colors = {"Rater": "#dc2626", "Task": "#2563eb", "Criterion": "#059669"}
    fig, ax = plt.subplots(figsize=(8.3, 7.0))
    for facet, group in data.groupby("Facet", sort=False):
        included = group.loc[group["BothIncluded"]]
        excluded = group.loc[~group["BothIncluded"]]
        ax.scatter(
            included["EstimatePythonJMLE"], included["EstimateTAMMML"],
            s=42, alpha=0.68, color=colors.get(facet, "#6b7280"), label=facet,
        )
        if len(excluded):
            ax.scatter(
                excluded["EstimatePythonJMLE"], excluded["EstimateTAMMML"],
                marker="x", s=58, linewidth=1.6,
                color=colors.get(facet, "#6b7280"),
            )
    values = data[["EstimatePythonJMLE", "EstimateTAMMML"]].to_numpy()
    limits = [float(np.nanmin(values)) - 0.15, float(np.nanmax(values)) + 0.15]
    ax.plot(limits, limits, "--", color="#111827", linewidth=1.2, label="Identity reference")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_xlabel("Python additive RSM JMLE facet estimate (logits)")
    ax.set_ylabel("TAM additive RSM MML facet estimate (logits)")
    ax.set_title("Facet estimates: same response surface, different Person treatment")
    ax.grid(alpha=0.2)
    handles, labels = ax.get_legend_handles_labels()
    if (~data["BothIncluded"]).any():
        handles.append(Line2D([0], [0], marker="x", linestyle="none", color="#4b5563"))
        labels.append("returned but excluded")
    ax.legend(handles, labels, frameon=False, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_quadrature_sensitivity(pairs: pd.DataFrame, path: Path) -> None:
    order = [
        "Facet: Criterion", "Facet: Task", "Facet: Rater", "Common step",
        "Cumulative response-surface difficulty", "Person EAP",
    ]
    data = pairs.loc[pairs["BothPresent"] & pairs["ParameterType"].isin(order)]
    fig, ax = plt.subplots(figsize=(12.3, 6.5))
    rng = np.random.default_rng(20260809)
    for index, parameter_type in enumerate(order):
        block = data.loc[data["ParameterType"].eq(parameter_type)].dropna(subset=["AbsDifference"])
        included = block.loc[block["BothIncluded"]]
        excluded = block.loc[~block["BothIncluded"]]
        ax.scatter(
            np.full(len(included), index) + rng.uniform(-0.14, 0.14, len(included)),
            included["AbsDifference"], s=18, alpha=0.45,
        )
        if len(excluded):
            ax.scatter(
                np.full(len(excluded), index) + rng.uniform(-0.14, 0.14, len(excluded)),
                excluded["AbsDifference"], marker="x", s=28, color="#dc2626", alpha=0.8,
            )
        if len(included):
            ax.scatter(index, included["AbsDifference"].median(), marker="D", s=68, color="#111827", zorder=4)
    for value, style in [(1e-1, ":"), (1e-2, "--"), (1e-3, ":")]:
        ax.axhline(value, linestyle=style, linewidth=1.1, label=f"{value:g} logit reference")
    ax.set_yscale("log")
    ax.set_xticks(
        range(len(order)),
        ["Criterion", "Task", "Rater", "Common\nstep", "Cumulative\nsurface", "Person\nEAP"],
    )
    ax.set_ylabel("Absolute q61 - q21 estimate difference")
    ax.set_title("TAM quadrature sensitivity and retained exclusions")
    ax.grid(alpha=0.2, which="both", axis="y")
    handles, labels = ax.get_legend_handles_labels()
    handles.extend([
        Line2D([0], [0], marker="D", linestyle="none", color="#111827"),
        Line2D([0], [0], marker="x", linestyle="none", color="#dc2626"),
    ])
    labels.extend(["eligible median", "returned but excluded"])
    ax.legend(handles, labels, frameon=False, fontsize=8.3, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_anchor_transmission(anchor: pd.DataFrame, path: Path) -> None:
    data = anchor.loc[anchor["Facet"].eq("Rater")].copy()
    order = sorted(data["Level"].unique())
    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    positions = {level: index for index, level in enumerate(order)}
    for level, group in data.groupby("Level", sort=False):
        x = positions[level]
        ax.scatter(
            np.full(len(group), x) + np.linspace(-0.08, 0.08, len(group)),
            group["EstimateShiftDriftMinusClean"],
            s=48,
            color="#d97706" if group["AnchoredClean"].all() else "#2563eb",
            alpha=0.8,
        )
    ax.axhline(0.25, linestyle="--", color="#d97706", label="injected anchor shift")
    ax.axhline(0, linestyle=":", color="#111827", label="no shift")
    ax.set_xticks(range(len(order)), order)
    ax.set_ylabel("Drift - clean TAM Rater estimate (logits)")
    ax.set_title("TAM sum-zero constraint redistributes contaminated anchors")
    ax.grid(alpha=0.2, axis="y")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_results(
    output_dir: Path,
    runs: pd.DataFrame,
    q_pairs: pd.DataFrame,
    python_tam: pd.DataFrame,
    tam_sirt: pd.DataFrame,
    anchors: pd.DataFrame,
    omitted: pd.DataFrame,
    variance_audit: pd.DataFrame,
    sparse: pd.DataFrame,
) -> None:
    eligible_q = q_pairs.loc[q_pairs["BothIncluded"]]
    qmax = eligible_q.groupby("ParameterType")["AbsDifference"].max()
    balanced_rater = python_tam.loc[
        python_tam["Design"].eq("balanced_small") & python_tam["Facet"].eq("Rater")
    ].iloc[0]
    sparse_rater = python_tam.loc[
        python_tam["Design"].eq("sparse_missing") & python_tam["Facet"].eq("Rater")
    ].iloc[0]
    balanced_sirt = tam_sirt.loc[tam_sirt["Design"].eq("balanced_small")].iloc[0]
    drift_sirt = tam_sirt.loc[tam_sirt["Design"].eq("anchor_drift")].iloc[0]
    anchored = anchors.loc[anchors["ExpectedConstraintRole"].eq("fixed anchor receives +0.25")]
    compensating = anchors.loc[
        anchors["ExpectedConstraintRole"].eq(
            "unanchored Rater compensates under sum-zero constraint"
        )
    ]
    balanced_omitted = omitted.loc[
        omitted["PairEligible"] & omitted["Design"].eq("balanced_small")
        & omitted["Facet"].eq("Rater")
    ]
    mismatch = int(
        variance_audit[["NaiveVsTolerantMismatchQ21", "NaiveVsTolerantMismatchQ61"]]
        .to_numpy(dtype=bool)
        .sum()
    )
    text = f"""# TAM tam.mml.mfr operating-characteristics smoke sensitivity

## First read

- TAM uses the registered additive RSM surface directly: Criterion item + Rater + Task + common step. It is still normal-population MML, not Python JMLE or exact CMLE parity.
- All {int(runs['FitReturned'].sum())}/{len(runs)} q21/q61 fits returned and met the retained loop/progress convergence audit. {int(runs['AnalysisEligible'].sum())}/{len(runs)} were analysis eligible. Two q21 sparse fits were excluded because Person variance was on TAM's configured lower boundary; all {int(runs.loc[runs['Mode'].eq(Q61), 'AnalysisEligible'].sum())}/16 primary q61 fits were eligible.
- The variance boundary is a concrete floating-point case: TAM retained `0.0010000001` against a configured `0.001` lower bound. A naive exact comparison missed {mismatch} boundary decisions; the prespecified `1e-8` numerical band caught them. Display rounding is not used.
- In balanced-small data, Python JMLE versus primary TAM MML Rater estimates differed by `{balanced_rater['MeanAbsDifference']:.4g}` logits on average and `{balanced_rater['MaxAbsDifference']:.4g}` at most. The additive surface matches, but Person treatment does not.
- Sparse data again exposed the estimand difference. Primary TAM was eligible for {int(sparse['AnalysisEligible'].sum())}/4 RunIds while exact CMLE and strict mfrmr rejected all four structurally. Python/TAM Rater MAE was `{sparse_rater['MeanAbsDifference']:.4g}` and the maximum was `{sparse_rater['MaxAbsDifference']:.4g}` logits.
- TAM and sirt were close for balanced-small Raters (MAE `{balanced_sirt['MeanAbsDifference']:.4g}`, maximum `{balanced_sirt['MaxAbsDifference']:.4g}`), but anchor drift increased their MAE to `{drift_sirt['MeanAbsDifference']:.4g}` because their constraint and item-threshold mechanics differ.
- Contaminated anchors were reproduced exactly and shifted by `{anchored['EstimateShiftDriftMinusClean'].mean():.4g}` logits. Under TAM's Rater sum-zero constraint, unanchored Raters compensated by `{compensating['EstimateShiftDriftMinusClean'].mean():.4g}` logits on average. Exact anchor reproduction is therefore not robustness.
- Among jointly eligible q21/q61 runs, maximum differences were `{qmax['Facet: Rater']:.4g}` logits for Rater, `{qmax['Cumulative response-surface difficulty']:.4g}` for cumulative response-surface difficulty, and `{qmax['Person EAP']:.4g}` for Person EAP. Returned boundary fits remain visible under the all-returned scope.
- The omitted local Rater x Task interaction moved balanced-small TAM Rater main effects by `{balanced_omitted['AbsEstimateChange'].mean():.4g}` logits on average and `{balanced_omitted['AbsEstimateChange'].max():.4g}` at most. This is leakage, not bias detection.
- Derived last-level facet SE coverage is withheld. TAM 4.3-25 expands those SEs from a diagonal matrix of free-xi SEs, without a retained full covariance. Information criteria are retained only for within-TAM audit and are not compared across estimators.

## Interpretation boundary

This remains a two-replicate smoke study. It does not establish package superiority, stable recovery, coverage, a sample-size rule, or MML/JMLE equivalence. The public Streamlit surface remains withheld.

## Retained evidence

- `tam_runs.csv`, `tam_facet_recovery.csv`, `tam_step_estimates.csv`, `tam_person_estimates.csv`, `tam_surface_estimates.csv`
- `python_tam_facet_pairs.csv`, `python_tam_facet_agreement.csv`
- `tam_sirt_rater_pairs.csv`, `tam_sirt_rater_agreement.csv`
- `tam_quadrature_pairs.csv`, `tam_quadrature_summary.csv`, `tam_quadrature_run_pairs.csv`
- `tam_anchor_sensitivity.csv`, `tam_omitted_bias_sensitivity.csv`, `tam_sparse_identification.csv`
- `tam_variance_boundary_audit.csv`, `tam_progress_threshold_summary.csv`
- `tam_first_read_summary.csv` is the ordered future-UI projection; its public status remains `Withheld`.
- `python_tam_facet_agreement.png`, `tam_quadrature_sensitivity.png`, `tam_anchor_constraint_transmission.png`
- `tam_adapter_identity.csv` records bundle, TAM MFR/progress functions, and adapter-script hashes.
"""
    (output_dir / "TAM_RESULTS.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    input_dir = args.input.resolve()
    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    runs = pd.read_csv(input_dir / "tam_runs.csv")
    facets = pd.read_csv(input_dir / "tam_facet_recovery.csv")
    steps = pd.read_csv(input_dir / "tam_step_estimates.csv")
    persons = pd.read_csv(input_dir / "tam_person_estimates.csv")
    surfaces = pd.read_csv(input_dir / "tam_surface_estimates.csv")
    identity = pd.read_csv(input_dir / "tam_adapter_identity.csv").iloc[0]
    if sha256_file(input_dir / "generated_bundle_files.csv") != str(
        identity["BundleInventorySHA256"]
    ):
        raise RuntimeError("TAM adapter is stale for the current generated bundle.")
    adapter_path = Path(__file__).with_name("operating_characteristics_tam.R")
    if sha256_file(adapter_path) != str(identity["AdapterScriptSHA256"]):
        raise RuntimeError("TAM adapter script changed after its retained run.")

    python_recovery = pd.read_csv(input_dir / "parameter_recovery.csv")
    python_runs = pd.read_csv(input_dir / "runs.csv")
    cmle_runs = pd.read_csv(input_dir / "python_cmle_runs.csv")
    mfrmr_runs = pd.read_csv(input_dir / "mfrmr_runs.csv")
    sirt_runs = pd.read_csv(input_dir / "sirt_runs.csv")
    sirt_recovery = pd.read_csv(input_dir / "sirt_rater_recovery.csv")

    q_pairs, q_summary, q_run_pairs = build_quadrature_evidence(
        facets, steps, persons, surfaces, runs
    )
    python_tam_pairs, python_tam_summary = build_python_tam_agreement(
        python_recovery, facets
    )
    tam_sirt_pairs, tam_sirt_summary = build_tam_sirt_rater_agreement(
        facets, sirt_recovery
    )
    anchor = build_anchor_sensitivity(facets, runs)
    omitted = build_omitted_bias_sensitivity(facets)
    variance = build_variance_boundary_audit(runs)
    progress = build_progress_threshold_summary(runs)
    sparse = build_sparse_identification(
        runs, sirt_runs, python_runs, cmle_runs, mfrmr_runs
    )
    first_read = build_first_read_summary(
        runs, q_pairs, python_tam_summary, tam_sirt_summary,
        anchor, variance, sparse,
    )

    outputs = {
        "tam_quadrature_pairs.csv": q_pairs,
        "tam_quadrature_summary.csv": q_summary,
        "tam_quadrature_run_pairs.csv": q_run_pairs,
        "python_tam_facet_pairs.csv": python_tam_pairs,
        "python_tam_facet_agreement.csv": python_tam_summary,
        "tam_sirt_rater_pairs.csv": tam_sirt_pairs,
        "tam_sirt_rater_agreement.csv": tam_sirt_summary,
        "tam_anchor_sensitivity.csv": anchor,
        "tam_omitted_bias_sensitivity.csv": omitted,
        "tam_variance_boundary_audit.csv": variance,
        "tam_progress_threshold_summary.csv": progress,
        "tam_sparse_identification.csv": sparse,
        "tam_first_read_summary.csv": first_read,
    }
    for filename, frame in outputs.items():
        frame.to_csv(output_dir / filename, index=False)
    plot_python_tam_agreement(
        python_tam_pairs, output_dir / "python_tam_facet_agreement.png"
    )
    plot_quadrature_sensitivity(
        q_pairs, output_dir / "tam_quadrature_sensitivity.png"
    )
    plot_anchor_transmission(
        anchor, output_dir / "tam_anchor_constraint_transmission.png"
    )
    write_results(
        output_dir, runs, q_pairs, python_tam_summary, tam_sirt_summary,
        anchor, omitted, variance, sparse,
    )


if __name__ == "__main__":
    main()
