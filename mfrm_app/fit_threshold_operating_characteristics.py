"""Known-truth pilot contracts for fixed-calibration Person MnSq rules.

The simulation in this module deliberately conditions on a known working
calibration.  It isolates Person scoring and fit-rule behavior; it does not
represent CMLE calibration uncertainty, anchors, or cross-engine equivalence.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
import math
from typing import Any

import numpy as np
import pandas as pd

from mfrm_app.cmle_wle_fit import compute_fixed_calibration_person_fit
from mfrm_app.fit_threshold_sensitivity import (
    CANONICAL_FIT_THRESHOLDS,
    classify_fit_mnsq_array,
    validate_fit_thresholds,
)
from mfrm_app.person_scoring import score_fixed_calibration_persons


FIT_RULES = (
    "either_upper",
    "either_distorting",
    "either_nonacceptable",
    "either_overfit",
)
PERSON_MEASURE_WLE = "fit_sample_wle_primary"
PERSON_MEASURE_TRUE = "generating_theta_sensitivity"


def rsm_category_intercepts(
    unit_difficulties: Sequence[float],
    category_thresholds: Sequence[float],
) -> np.ndarray:
    """Return row-by-category RSM log kernels without the Person term."""

    difficulties = np.asarray(unit_difficulties, dtype=float).reshape(-1)
    thresholds = np.asarray(category_thresholds, dtype=float).reshape(-1)
    if not len(difficulties) or not np.isfinite(difficulties).all():
        raise ValueError("unit_difficulties must be a non-empty finite sequence.")
    if not len(thresholds) or not np.isfinite(thresholds).all():
        raise ValueError("category_thresholds must be a non-empty finite sequence.")
    if not np.all(np.diff(thresholds) > 0.0):
        raise ValueError("category_thresholds must be strictly increasing.")
    categories = np.arange(len(thresholds) + 1, dtype=float)
    cumulative = np.r_[0.0, np.cumsum(thresholds)]
    return -difficulties[:, None] * categories[None, :] - cumulative[None, :]


def _sample_from_uniform(probabilities: np.ndarray, uniforms: np.ndarray) -> np.ndarray:
    probabilities = np.asarray(probabilities, dtype=float)
    uniforms = np.asarray(uniforms, dtype=float).reshape(-1)
    if probabilities.ndim != 2 or len(uniforms) != len(probabilities):
        raise ValueError("probabilities and uniforms must be row-aligned.")
    if not np.isfinite(probabilities).all() or np.any(probabilities < 0.0):
        raise ValueError("probabilities must be finite and non-negative.")
    row_sums = probabilities.sum(axis=1)
    if not np.allclose(row_sums, 1.0, rtol=0.0, atol=1e-12):
        raise ValueError("probability rows must sum to one.")
    if np.any((uniforms < 0.0) | (uniforms >= 1.0) | ~np.isfinite(uniforms)):
        raise ValueError("uniforms must lie in [0, 1).")
    cumulative = np.cumsum(probabilities, axis=1)
    cumulative[:, -1] = 1.0
    return np.sum(uniforms[:, None] >= cumulative, axis=1).astype(int)


def _probabilities(intercepts: np.ndarray, theta: np.ndarray) -> np.ndarray:
    categories = np.arange(intercepts.shape[1], dtype=float)
    logits = intercepts + theta[:, None] * categories[None, :]
    logits -= np.max(logits, axis=1, keepdims=True)
    exponentiated = np.exp(logits)
    return exponentiated / exponentiated.sum(axis=1, keepdims=True)


def _positive_integer(value: Any, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a positive integer.")
    numeric = int(value)
    if float(value) != numeric or numeric < 1:
        raise ValueError(f"{name} must be a positive integer.")
    return numeric


def simulate_known_truth_person_fit_sample(
    condition: Mapping[str, Any],
    *,
    seed: int,
    persons: int,
    rater_effects: Sequence[float],
    criterion_effects: Sequence[float],
    category_thresholds: Sequence[float],
    threshold_deviations: Mapping[str, Sequence[float]] | None = None,
) -> dict[str, object]:
    """Generate one paired known-truth response sample for a frozen condition."""

    n_persons = _positive_integer(persons, "persons")
    raters = np.asarray(rater_effects, dtype=float).reshape(-1)
    criteria = np.asarray(criterion_effects, dtype=float).reshape(-1)
    if not len(raters) or not len(criteria):
        raise ValueError("At least one Rater and Criterion are required.")
    if not np.isfinite(raters).all() or not np.isfinite(criteria).all():
        raise ValueError("Rater and Criterion effects must be finite.")
    thresholds = np.asarray(category_thresholds, dtype=float).reshape(-1)
    categories = int(condition["Categories"])
    if categories != len(thresholds) + 1:
        raise ValueError("Condition Categories differs from category_thresholds.")
    observations = _positive_integer(
        condition["ObservationsPerPerson"], "ObservationsPerPerson"
    )
    n_units = len(raters) * len(criteria)
    if observations > n_units:
        raise ValueError("ObservationsPerPerson exceeds the available unit count.")
    mechanism = str(condition["Mechanism"])
    supported = {
        "clean",
        "uniform_random_replacement",
        "previous_response_copy",
        "affected_person_pcm_threshold_heterogeneity",
    }
    if mechanism not in supported:
        raise ValueError(f"Unsupported known-truth mechanism: {mechanism}")
    affected_share = float(condition["AffectedShare"])
    mechanism_probability = float(condition["MechanismProbability"])
    if not 0.0 <= affected_share <= 1.0 or not 0.0 <= mechanism_probability <= 1.0:
        raise ValueError("AffectedShare and MechanismProbability must lie in [0, 1].")
    if mechanism == "clean" and (affected_share != 0.0 or mechanism_probability != 0.0):
        raise ValueError("Clean conditions must have zero affected share and probability.")

    unit_rows = [
        (rater_index, criterion_index)
        for rater_index in range(len(raters))
        for criterion_index in range(len(criteria))
    ]
    unit_difficulty = np.asarray(
        [raters[rater] + criteria[criterion] for rater, criterion in unit_rows]
    )
    unit_working_intercepts = rsm_category_intercepts(unit_difficulty, thresholds)
    sequence = np.random.SeedSequence(int(seed))
    theta_rng, administration_rng, response_rng, affected_rng, mechanism_rng, category_rng = [
        np.random.default_rng(child) for child in sequence.spawn(6)
    ]
    theta_by_person = theta_rng.normal(size=n_persons)
    affected_count = int(round(n_persons * affected_share))
    affected_indices = set(
        affected_rng.permutation(n_persons)[:affected_count].astype(int).tolist()
    )

    person_values: list[str] = []
    person_indices: list[int] = []
    selected_units: list[int] = []
    width = max(4, len(str(n_persons)))
    for person_index in range(n_persons):
        units = (
            np.arange(n_units, dtype=int)
            if observations == n_units
            else np.sort(
                administration_rng.choice(n_units, size=observations, replace=False)
            )
        )
        person_values.extend([f"P{person_index + 1:0{width}d}"] * observations)
        person_indices.extend([person_index] * observations)
        selected_units.extend(units.astype(int).tolist())
    person_index_array = np.asarray(person_indices, dtype=int)
    unit_index_array = np.asarray(selected_units, dtype=int)
    row_theta = theta_by_person[person_index_array]
    working_intercepts = unit_working_intercepts[unit_index_array].copy()
    true_intercepts = working_intercepts.copy()
    affected = np.fromiter(
        (index in affected_indices for index in person_index_array),
        dtype=bool,
        count=len(person_index_array),
    )

    if mechanism == "affected_person_pcm_threshold_heterogeneity":
        if threshold_deviations is None:
            raise ValueError("Threshold heterogeneity requires threshold_deviations.")
        for criterion_index in range(len(criteria)):
            criterion = f"C{criterion_index + 1:02d}"
            if criterion not in threshold_deviations:
                raise ValueError(f"Missing threshold deviations for {criterion}.")
            deviations = np.asarray(threshold_deviations[criterion], dtype=float)
            if deviations.shape != thresholds.shape or not np.isfinite(deviations).all():
                raise ValueError("Threshold deviations must match the threshold vector.")
            altered = thresholds + deviations
            if not np.all(np.diff(altered) > 0.0):
                raise ValueError("Altered thresholds must remain strictly increasing.")
            rows = affected & np.fromiter(
                (unit_rows[unit][1] == criterion_index for unit in unit_index_array),
                dtype=bool,
                count=len(unit_index_array),
            )
            if rows.any():
                true_intercepts[rows] = rsm_category_intercepts(
                    unit_difficulty[unit_index_array[rows]], altered
                )

    uniforms = response_rng.random(len(person_index_array))
    baseline = _sample_from_uniform(
        _probabilities(working_intercepts, row_theta), uniforms
    )
    observed = (
        _sample_from_uniform(_probabilities(true_intercepts, row_theta), uniforms)
        if mechanism == "affected_person_pcm_threshold_heterogeneity"
        else baseline.copy()
    )
    mechanism_applied = np.zeros(len(observed), dtype=bool)
    if mechanism == "uniform_random_replacement":
        apply = affected & (mechanism_rng.random(len(observed)) < mechanism_probability)
        replacement = category_rng.integers(0, categories, size=len(observed))
        observed[apply] = replacement[apply]
        mechanism_applied = apply
    elif mechanism == "previous_response_copy":
        for person_index in range(n_persons):
            rows = np.flatnonzero(person_index_array == person_index)
            if person_index not in affected_indices:
                continue
            for position in range(1, len(rows)):
                if mechanism_rng.random() < mechanism_probability:
                    observed[rows[position]] = observed[rows[position - 1]]
                    mechanism_applied[rows[position]] = True
    elif mechanism == "affected_person_pcm_threshold_heterogeneity":
        mechanism_applied = affected.copy()

    rater_index = np.asarray([unit_rows[unit][0] for unit in unit_index_array])
    criterion_index = np.asarray([unit_rows[unit][1] for unit in unit_index_array])
    responses = pd.DataFrame(
        {
            "Person": person_values,
            "PersonIndex": person_index_array,
            "Rater": [f"R{index + 1:02d}" for index in rater_index],
            "Criterion": [f"C{index + 1:02d}" for index in criterion_index],
            "UnitIndex": unit_index_array,
            "ObservedCategory": observed,
            "BaselineCategory": baseline,
            "TrueTheta": row_theta,
            "Affected": affected,
            "TruthGroup": np.where(affected, "affected", "clean"),
            "Mechanism": mechanism,
            "MechanismApplied": mechanism_applied,
            "ResponseModified": observed != baseline,
        }
    )
    slopes = np.tile(np.arange(categories, dtype=float), (len(responses), 1))
    return {
        "responses": responses,
        "working_intercepts": working_intercepts,
        "true_intercepts": true_intercepts,
        "category_slopes": slopes,
        "person_theta": pd.DataFrame(
            {
                "Person": [f"P{index + 1:0{width}d}" for index in range(n_persons)],
                "TrueTheta": theta_by_person,
                "Affected": [index in affected_indices for index in range(n_persons)],
            }
        ),
        "audit": {
            "Seed": int(seed),
            "Persons": n_persons,
            "Rows": int(len(responses)),
            "Categories": categories,
            "AffectedPersons": affected_count,
            "MechanismAppliedRows": int(mechanism_applied.sum()),
            "ModifiedRows": int((observed != baseline).sum()),
            "CategorySupportViolations": int(
                np.sum((observed < 0) | (observed >= categories))
            ),
        },
    }


def score_known_truth_person_fit(sample: Mapping[str, object]) -> pd.DataFrame:
    """Score WLE and generating-theta Person fit under the working calibration."""

    responses = sample["responses"]
    if not isinstance(responses, pd.DataFrame) or responses.empty:
        raise ValueError("sample responses must be a non-empty DataFrame.")
    intercepts = np.asarray(sample["working_intercepts"], dtype=float)
    slopes = np.asarray(sample["category_slopes"], dtype=float)
    persons = responses["Person"].astype(str).to_numpy(dtype=object)
    observed = responses["ObservedCategory"].to_numpy(dtype=int)
    levels = list(dict.fromkeys(persons.tolist()))
    wle = score_fixed_calibration_persons(
        persons,
        observed,
        intercepts,
        slopes,
        person_levels=levels,
        method="WLE",
    )
    if not wle["Status"].eq("ok").all() or not np.isfinite(wle["Estimate"]).all():
        raise RuntimeError("At least one fixed-calibration WLE was unavailable.")
    truth = (
        responses.groupby("Person", sort=False)
        .agg(
            TrueTheta=("TrueTheta", "first"),
            Affected=("Affected", "first"),
            TruthGroup=("TruthGroup", "first"),
            Mechanism=("Mechanism", "first"),
        )
        .reset_index()
    )
    sources = (
        (PERSON_MEASURE_WLE, wle.set_index("Person")["Estimate"]),
        (PERSON_MEASURE_TRUE, truth.set_index("Person")["TrueTheta"]),
    )
    parts: list[pd.DataFrame] = []
    for source, estimates in sources:
        row_theta = pd.Series(persons).map(estimates).to_numpy(dtype=float)
        result = compute_fixed_calibration_person_fit(
            persons,
            observed,
            intercepts,
            row_theta,
            person_levels=levels,
        )["persons"].copy()
        result["PersonMeasureSource"] = source
        result["PersonMeasure"] = result["Person"].map(estimates)
        result = result.merge(truth, on="Person", how="left", validate="one_to_one")
        if source == PERSON_MEASURE_WLE:
            result = result.merge(
                wle[
                    [
                        "Person",
                        "StandardError",
                        "ExtremeScorePattern",
                        "Status",
                        "AdjustedScoreResidual",
                    ]
                ].rename(
                    columns={
                        "StandardError": "WLEStandardError",
                        "ExtremeScorePattern": "WLEExtremeScorePattern",
                        "Status": "WLEStatus",
                        "AdjustedScoreResidual": "WLEAdjustedScoreResidual",
                    }
                ),
                on="Person",
                how="left",
                validate="one_to_one",
            )
        else:
            result["WLEStandardError"] = np.nan
            result["WLEExtremeScorePattern"] = pd.NA
            result["WLEStatus"] = "not_applicable_true_theta_sensitivity"
            result["WLEAdjustedScoreResidual"] = np.nan
        parts.append(result)
    return pd.concat(parts, ignore_index=True)


def fit_rule_flags(
    infit: Sequence[float],
    outfit: Sequence[float],
    thresholds: Sequence[float] = CANONICAL_FIT_THRESHOLDS,
) -> dict[str, np.ndarray]:
    """Evaluate the four frozen raw MnSq rules with exact endpoints."""

    lower, acceptable, noisy = validate_fit_thresholds(*thresholds)
    infit_values = np.asarray(infit, dtype=float).reshape(-1)
    outfit_values = np.asarray(outfit, dtype=float).reshape(-1)
    if infit_values.shape != outfit_values.shape:
        raise ValueError("Infit and Outfit must have equal lengths.")
    finite = np.isfinite(infit_values) & np.isfinite(outfit_values)
    return {
        "eligible": finite,
        "either_upper": finite
        & ((infit_values > acceptable) | (outfit_values > acceptable)),
        "either_distorting": finite
        & ((infit_values > noisy) | (outfit_values > noisy)),
        "either_nonacceptable": finite
        & (
            (infit_values < lower)
            | (infit_values > acceptable)
            | (outfit_values < lower)
            | (outfit_values > acceptable)
        ),
        "either_overfit": finite
        & ((infit_values < lower) | (outfit_values < lower)),
    }


def replicate_fit_rule_rates(
    persons: pd.DataFrame,
    *,
    threshold_triplets: Iterable[Sequence[float]],
    identity_columns: Sequence[str] = ("ConditionId", "Replicate", "PersonMeasureSource"),
) -> pd.DataFrame:
    """Aggregate rule rates within the Person clusters of one or more replicates."""

    required = set(identity_columns) | {"TruthGroup", "Infit", "Outfit"}
    missing = sorted(required - set(persons.columns))
    if missing:
        raise ValueError("Person-fit rows are missing columns: " + ", ".join(missing))
    triplets = [validate_fit_thresholds(*values) for values in threshold_triplets]
    if not triplets or len(triplets) != len(set(triplets)):
        raise ValueError("threshold_triplets must be non-empty and unique.")
    rows: list[dict[str, object]] = []
    group_columns = [*identity_columns, "TruthGroup"]
    group_key: str | list[str] = group_columns[0] if len(group_columns) == 1 else group_columns
    for identity, frame in persons.groupby(group_key, sort=False, dropna=False):
        values = identity if isinstance(identity, tuple) else (identity,)
        identity_row = dict(zip(group_columns, values))
        for thresholds in triplets:
            flags = fit_rule_flags(frame["Infit"], frame["Outfit"], thresholds)
            eligible = flags["eligible"]
            denominator = int(eligible.sum())
            for rule in FIT_RULES:
                numerator = int(flags[rule].sum())
                rows.append(
                    {
                        **identity_row,
                        "OverfitUpper": thresholds[0],
                        "AcceptableUpper": thresholds[1],
                        "NoisyUpper": thresholds[2],
                        "CanonicalThresholds": thresholds
                        == CANONICAL_FIT_THRESHOLDS,
                        "Rule": rule,
                        "PersonsAttempted": int(len(frame)),
                        "PersonsEligible": denominator,
                        "PersonsUnavailable": int(len(frame) - denominator),
                        "Flagged": numerator,
                        "Rate": float(numerator / denominator)
                        if denominator
                        else np.nan,
                        "IndependentMonteCarloUnit": "replicate",
                        "PersonRowsClustered": True,
                    }
                )
    return pd.DataFrame(rows)


def summarize_replicate_fit_rates(
    replicate_rates: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize replicate-specific rates without a naive Person-level Wilson CI."""

    identity = [
        "ConditionId",
        "PersonMeasureSource",
        "TruthGroup",
        "OverfitUpper",
        "AcceptableUpper",
        "NoisyUpper",
        "CanonicalThresholds",
        "Rule",
    ]
    missing = sorted(set(identity + ["Replicate", "PersonsAttempted", "PersonsEligible", "PersonsUnavailable", "Flagged", "Rate"]) - set(replicate_rates.columns))
    if missing:
        raise ValueError("Replicate-rate rows are missing columns: " + ", ".join(missing))
    rows: list[dict[str, object]] = []
    for values, frame in replicate_rates.groupby(identity, sort=False, dropna=False):
        rates = pd.to_numeric(frame["Rate"], errors="coerce")
        available = rates[np.isfinite(rates)]
        replicate_sd = float(available.std(ddof=1)) if len(available) > 1 else np.nan
        total_eligible = int(frame["PersonsEligible"].sum())
        total_flagged = int(frame["Flagged"].sum())
        truth_group = str(values[2])
        rows.append(
            {
                **dict(zip(identity, values)),
                "ReplicatesAttempted": int(frame["Replicate"].nunique()),
                "ReplicatesAvailable": int(len(available)),
                "PersonsAttempted": int(frame["PersonsAttempted"].sum()),
                "PersonsEligible": total_eligible,
                "PersonsUnavailable": int(frame["PersonsUnavailable"].sum()),
                "Flagged": total_flagged,
                "MeanReplicateRate": float(available.mean())
                if len(available)
                else np.nan,
                "BetweenReplicateSD": replicate_sd,
                "ReplicateRateMCSE": float(replicate_sd / math.sqrt(len(available)))
                if len(available) > 1
                else np.nan,
                "PooledPersonRateDescriptive": float(total_flagged / total_eligible)
                if total_eligible
                else np.nan,
                "RateLabel": (
                    "false_positive_rate_pilot"
                    if truth_group == "clean"
                    else "detection_rate_pilot_not_confirmatory_power"
                ),
                "IndependentBinomialWilsonIntervalAuthorized": False,
                "ConfirmatoryPerformanceClaim": False,
            }
        )
    return pd.DataFrame(rows)


def rounded_rule_disagreements(
    persons: pd.DataFrame,
    *,
    decimals: Iterable[int],
    thresholds: Sequence[float] = CANONICAL_FIT_THRESHOLDS,
) -> pd.DataFrame:
    """Return raw-versus-rounded rule and class disagreement counts."""

    required = {"ConditionId", "Replicate", "PersonMeasureSource", "TruthGroup", "Infit", "Outfit"}
    missing = sorted(required - set(persons.columns))
    if missing:
        raise ValueError("Person-fit rows are missing columns: " + ", ".join(missing))
    precisions = [int(value) for value in decimals]
    if not precisions or len(precisions) != len(set(precisions)) or any(value < 0 for value in precisions):
        raise ValueError("decimals must be a non-empty unique sequence of non-negative integers.")
    threshold_values = validate_fit_thresholds(*thresholds)
    rows: list[dict[str, object]] = []
    group_columns = ["ConditionId", "Replicate", "PersonMeasureSource", "TruthGroup"]
    for identity, frame in persons.groupby(group_columns, sort=False, dropna=False):
        raw_infit = frame["Infit"].to_numpy(dtype=float)
        raw_outfit = frame["Outfit"].to_numpy(dtype=float)
        raw_flags = fit_rule_flags(raw_infit, raw_outfit, threshold_values)
        raw_infit_class = classify_fit_mnsq_array(raw_infit, threshold_values)
        raw_outfit_class = classify_fit_mnsq_array(raw_outfit, threshold_values)
        for precision in precisions:
            rounded_infit = np.round(raw_infit, precision)
            rounded_outfit = np.round(raw_outfit, precision)
            rounded_flags = fit_rule_flags(
                rounded_infit, rounded_outfit, threshold_values
            )
            rounded_infit_class = classify_fit_mnsq_array(
                rounded_infit, threshold_values
            )
            rounded_outfit_class = classify_fit_mnsq_array(
                rounded_outfit, threshold_values
            )
            for rule in FIT_RULES:
                rows.append(
                    {
                        **dict(zip(group_columns, identity)),
                        "DisplayDecimals": precision,
                        "Rule": rule,
                        "Persons": int(len(frame)),
                        "RawFlagged": int(raw_flags[rule].sum()),
                        "RoundedFlagged": int(rounded_flags[rule].sum()),
                        "RuleDisagreements": int(
                            np.sum(raw_flags[rule] != rounded_flags[rule])
                        ),
                        "InfitClassMismatches": int(
                            np.sum(raw_infit_class != rounded_infit_class)
                        ),
                        "OutfitClassMismatches": int(
                            np.sum(raw_outfit_class != rounded_outfit_class)
                        ),
                        "RawDecisionRetained": True,
                    }
                )
    return pd.DataFrame(rows)


__all__ = [
    "FIT_RULES",
    "PERSON_MEASURE_TRUE",
    "PERSON_MEASURE_WLE",
    "fit_rule_flags",
    "replicate_fit_rule_rates",
    "rounded_rule_disagreements",
    "rsm_category_intercepts",
    "score_known_truth_person_fit",
    "simulate_known_truth_person_fit_sample",
    "summarize_replicate_fit_rates",
]
