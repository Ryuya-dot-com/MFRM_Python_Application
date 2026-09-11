"""Post-result decomposition helpers for the known-truth MnSq pilot."""

from __future__ import annotations

import hashlib
import math

import numpy as np
import pandas as pd

from mfrm_app.fit_threshold_operating_characteristics import (
    PERSON_MEASURE_TRUE,
    PERSON_MEASURE_WLE,
    fit_rule_flags,
)


RULE_DEPENDENCIES = {
    "either_upper": ("AcceptableUpper",),
    "either_distorting": ("NoisyUpper",),
    "either_nonacceptable": ("OverfitUpper", "AcceptableUpper"),
    "either_overfit": ("OverfitUpper",),
}
THRESHOLD_COLUMNS = ("OverfitUpper", "AcceptableUpper", "NoisyUpper")


def wle_person_recovery(persons: pd.DataFrame) -> pd.DataFrame:
    """Return one WLE-minus-generating-theta recovery row per Person."""

    required = {
        "ConditionId",
        "Replicate",
        "Person",
        "PersonMeasureSource",
        "PersonMeasure",
        "TrueTheta",
        "TruthGroup",
        "WLEExtremeScorePattern",
        "Infit",
        "Outfit",
    }
    missing = sorted(required - set(persons.columns))
    if missing:
        raise ValueError("Known-truth Person rows are missing columns: " + ", ".join(missing))
    work = persons.loc[persons["PersonMeasureSource"].eq(PERSON_MEASURE_WLE)].copy()
    if work.duplicated(["ConditionId", "Replicate", "Person"]).any():
        raise ValueError("WLE recovery identity contains duplicate Persons.")
    work["WLEError"] = work["PersonMeasure"] - work["TrueTheta"]
    work["AbsoluteWLEError"] = work["WLEError"].abs()
    work["SquaredWLEError"] = work["WLEError"] ** 2
    work["WLEExactExtremePattern"] = np.fromiter(
        (
            False if pd.isna(value) else bool(value)
            for value in work["WLEExtremeScorePattern"]
        ),
        dtype=bool,
        count=len(work),
    )
    flags = fit_rule_flags(work["Infit"], work["Outfit"])
    work["CanonicalEitherUpperFlag"] = flags["either_upper"]
    return work


def replicate_wle_recovery(recovery: pd.DataFrame) -> pd.DataFrame:
    """Compute replicate-specific WLE recovery metrics."""

    rows: list[dict[str, object]] = []
    groups = ["ConditionId", "Replicate", "TruthGroup", "WLEExactExtremePattern"]
    for identity, frame in recovery.groupby(groups, sort=False, dropna=False):
        error = frame["WLEError"].to_numpy(dtype=float)
        rows.append(
            {
                **dict(zip(groups, identity)),
                "Persons": int(len(frame)),
                "MeanBias": float(np.mean(error)),
                "MeanAbsoluteError": float(np.mean(np.abs(error))),
                "RootMeanSquaredError": float(np.sqrt(np.mean(error**2))),
                "MedianError": float(np.median(error)),
                "IndependentMonteCarloUnit": "replicate",
            }
        )
    return pd.DataFrame(rows)


def summarize_wle_recovery(replicate: pd.DataFrame) -> pd.DataFrame:
    """Summarize recovery metrics over independent replicates."""

    metrics = ("MeanBias", "MeanAbsoluteError", "RootMeanSquaredError", "MedianError")
    groups = ["ConditionId", "TruthGroup", "WLEExactExtremePattern"]
    rows: list[dict[str, object]] = []
    for identity, frame in replicate.groupby(groups, sort=False, dropna=False):
        for metric in metrics:
            values = pd.to_numeric(frame[metric], errors="coerce")
            values = values[np.isfinite(values)]
            sd = float(values.std(ddof=1)) if len(values) > 1 else np.nan
            rows.append(
                {
                    **dict(zip(groups, identity)),
                    "Metric": metric,
                    "ReplicatesAvailable": int(len(values)),
                    "MeanReplicateMetric": float(values.mean()) if len(values) else np.nan,
                    "BetweenReplicateSD": sd,
                    "ReplicateMetricMCSE": float(sd / math.sqrt(len(values)))
                    if len(values) > 1
                    else np.nan,
                    "ConfirmatoryRecoveryClaim": False,
                }
            )
    return pd.DataFrame(rows)


def recovery_by_fit_flag(recovery: pd.DataFrame) -> pd.DataFrame:
    """Describe WLE error by the canonical either-upper flag, without causality."""

    rows: list[dict[str, object]] = []
    groups = ["ConditionId", "TruthGroup", "CanonicalEitherUpperFlag"]
    for identity, frame in recovery.groupby(groups, sort=False, dropna=False):
        rows.append(
            {
                **dict(zip(groups, identity)),
                "Persons": int(len(frame)),
                "MeanBias": float(frame["WLEError"].mean()),
                "MeanAbsoluteError": float(frame["AbsoluteWLEError"].mean()),
                "RootMeanSquaredError": float(
                    np.sqrt(frame["SquaredWLEError"].mean())
                ),
                "DescriptiveAssociationOnly": True,
            }
        )
    return pd.DataFrame(rows)


def canonical_measure_source_pairs(persons: pd.DataFrame) -> pd.DataFrame:
    """Pair WLE and generating-theta canonical either-upper flags by Person."""

    required = {
        "ConditionId",
        "Replicate",
        "Person",
        "TruthGroup",
        "PersonMeasureSource",
        "Infit",
        "Outfit",
    }
    missing = sorted(required - set(persons.columns))
    if missing:
        raise ValueError("Known-truth Person rows are missing columns: " + ", ".join(missing))
    parts = []
    for source in (PERSON_MEASURE_WLE, PERSON_MEASURE_TRUE):
        frame = persons.loc[persons["PersonMeasureSource"].eq(source)].copy()
        flags = fit_rule_flags(frame["Infit"], frame["Outfit"])
        frame["EitherUpperFlag"] = flags["either_upper"]
        parts.append(
            frame[
                [
                    "ConditionId",
                    "Replicate",
                    "Person",
                    "TruthGroup",
                    "EitherUpperFlag",
                ]
            ].rename(columns={"EitherUpperFlag": source})
        )
    keys = ["ConditionId", "Replicate", "Person", "TruthGroup"]
    paired = parts[0].merge(parts[1], on=keys, how="outer", validate="one_to_one", indicator=True)
    paired["BothPresent"] = paired["_merge"].eq("both")
    wle_flag = paired[PERSON_MEASURE_WLE].fillna(False).astype(bool)
    true_flag = paired[PERSON_MEASURE_TRUE].fillna(False).astype(bool)
    paired["FlagComparison"] = np.select(
        [
            ~wle_flag & ~true_flag,
            wle_flag & ~true_flag,
            ~wle_flag & true_flag,
            wle_flag & true_flag,
        ],
        ["both_unflagged", "wle_only", "generating_theta_only", "both_flagged"],
        default="unavailable",
    )
    paired["FlagDisagreement"] = wle_flag != true_flag
    return paired.drop(columns="_merge")


def summarize_measure_source_pairs(pairs: pd.DataFrame) -> pd.DataFrame:
    """Summarize the four WLE-versus-generating-theta flag cells."""

    rows: list[dict[str, object]] = []
    groups = ["ConditionId", "TruthGroup"]
    cells = ("both_unflagged", "wle_only", "generating_theta_only", "both_flagged")
    for identity, frame in pairs.groupby(groups, sort=False, dropna=False):
        counts = frame["FlagComparison"].value_counts().to_dict()
        rows.append(
            {
                **dict(zip(groups, identity)),
                "Persons": int(len(frame)),
                "BothUnflagged": int(counts.get(cells[0], 0)),
                "WLEOnly": int(counts.get(cells[1], 0)),
                "GeneratingThetaOnly": int(counts.get(cells[2], 0)),
                "BothFlagged": int(counts.get(cells[3], 0)),
                "FlagDisagreements": int(frame["FlagDisagreement"].sum()),
                "IdentityUnavailable": int((~frame["BothPresent"]).sum()),
                "MeasureSourceSensitivityOnly": True,
            }
        )
    return pd.DataFrame(rows)


def _vector_hash(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype=np.int64).tobytes()).hexdigest()


def threshold_dimension_audit(rates: pd.DataFrame) -> pd.DataFrame:
    """Audit nominal triplets against each rule's effective threshold inputs."""

    identity = [
        "ConditionId",
        "Replicate",
        "PersonMeasureSource",
        "TruthGroup",
    ]
    required = set(identity) | set(THRESHOLD_COLUMNS) | {"Rule", "Flagged", "PersonsEligible"}
    missing = sorted(required - set(rates.columns))
    if missing:
        raise ValueError("Threshold-rate rows are missing columns: " + ", ".join(missing))
    rows: list[dict[str, object]] = []
    for rule, relevant in RULE_DEPENDENCIES.items():
        frame = rates.loc[rates["Rule"].eq(rule)].copy()
        irrelevant = [column for column in THRESHOLD_COLUMNS if column not in relevant]
        group_columns = [*identity, *relevant]
        invariant = frame.groupby(group_columns, dropna=False).agg(
            FlagValues=("Flagged", "nunique"),
            EligibleValues=("PersonsEligible", "nunique"),
        )
        mismatches = int(
            (invariant["FlagValues"].ne(1) | invariant["EligibleValues"].ne(1)).sum()
        )
        deduplicated = frame.drop_duplicates(
            [*identity, *relevant, "Flagged", "PersonsEligible"]
        )
        vector_hashes = []
        for _, configuration in deduplicated.groupby(list(relevant), sort=True):
            ordered = configuration.sort_values(identity)
            vector_hashes.append(_vector_hash(ordered["Flagged"].to_numpy()))
        rows.append(
            {
                "Rule": rule,
                "RelevantThresholds": ";".join(relevant),
                "IrrelevantThresholds": ";".join(irrelevant),
                "NominalThresholdTriplets": int(
                    frame[list(THRESHOLD_COLUMNS)].drop_duplicates().shape[0]
                ),
                "EffectiveInputConfigurations": int(
                    frame[list(relevant)].drop_duplicates().shape[0]
                ),
                "ObservedDistinctFlagCountVectors": int(len(set(vector_hashes))),
                "IrrelevantDimensionMismatchGroups": mismatches,
                "Passed": mismatches == 0,
            }
        )
    return pd.DataFrame(rows)


__all__ = [
    "RULE_DEPENDENCIES",
    "canonical_measure_source_pairs",
    "recovery_by_fit_flag",
    "replicate_wle_recovery",
    "summarize_measure_source_pairs",
    "summarize_wle_recovery",
    "threshold_dimension_audit",
    "wle_person_recovery",
]
