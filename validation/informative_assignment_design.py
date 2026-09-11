#!/usr/bin/env python3
"""Locked observation masks for the informative-assignment PCM screen."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from validation.estimand_distribution_study import FIT_MODEL
from validation.facets_pcm_boundary_pilot import (
    apply_observation_design,
    constrained_adjacent_design_audit,
    person_rater_components,
)
from validation.facets_pcm_known_truth_smoke import RATER_TRUTH


DESIGNS = (
    "complete",
    "planned_connected",
    "ability_severity_aligned_connected",
)
ALIGNED_PAIRS = (
    ("R01", "R02"),
    ("R01", "R03"),
    ("R02", "R04"),
    ("R03", "R04"),
)
EXPECTED_PAIR_MEANS = (-0.30, -0.15, 0.15, 0.30)


def ability_quartiles(frame: pd.DataFrame) -> pd.DataFrame:
    """Return one stable, exactly balanced Theta quartile per Person."""

    required = {"Person", "Theta"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Assignment frame is missing columns: {sorted(required - set(frame))}")
    persons = frame[["Person", "Theta"]].copy()
    theta_counts = persons.groupby("Person", sort=False)["Theta"].nunique(dropna=False)
    if not theta_counts.eq(1).all():
        raise ValueError("Each Person must have one finite Theta")
    persons = persons.drop_duplicates("Person")
    persons["Theta"] = pd.to_numeric(persons["Theta"], errors="raise")
    if not np.isfinite(persons["Theta"]).all():
        raise ValueError("Theta values must be finite")
    if len(persons) % 4 != 0:
        raise ValueError("Person count must be divisible by four")
    persons = persons.sort_values(["Theta", "Person"], kind="mergesort").reset_index(drop=True)
    persons["AbilityRank"] = np.arange(1, len(persons) + 1, dtype=int)
    persons["AbilityQuartile"] = (
        np.arange(len(persons), dtype=int) // (len(persons) // 4) + 1
    )
    return persons


def apply_informative_assignment_design(frame: pd.DataFrame, design: str) -> pd.DataFrame:
    """Apply one of the three frozen masks without inspecting realized scores."""

    design = str(design)
    if design in {"complete", "planned_connected"}:
        return apply_observation_design(frame, design)
    if design != "ability_severity_aligned_connected":
        raise ValueError(f"Unknown informative-assignment design: {design}")
    quartiles = ability_quartiles(frame)
    pair_map = {
        person: set(ALIGNED_PAIRS[int(quartile) - 1])
        for person, quartile in quartiles[["Person", "AbilityQuartile"]].itertuples(index=False)
    }
    keep = [
        str(rater) in pair_map[str(person)]
        for person, rater in frame[["Person", "Rater"]].itertuples(index=False)
    ]
    return frame.loc[np.asarray(keep, dtype=bool)].copy().reset_index(drop=True)


def informative_assignment_audit(frame: pd.DataFrame, design: str) -> dict[str, Any]:
    """Audit density, exposure, graph/rank, and the registered assignment mechanism."""

    output = apply_informative_assignment_design(frame, design)
    quartiles = ability_quartiles(frame)
    person_raters = (
        output.groupby("Person", sort=True)["Rater"]
        .agg(lambda values: "|".join(sorted(set(map(str, values)))))
        .rename("RaterPair")
        .reset_index()
    )
    person_raters["AssignedMeanSeverity"] = person_raters["RaterPair"].map(
        lambda value: float(np.mean([RATER_TRUTH[level] for level in value.split("|")]))
    )
    person_assignment = quartiles.merge(
        person_raters, on="Person", how="left", validate="one_to_one"
    )
    if person_assignment["RaterPair"].isna().any():
        raise RuntimeError("At least one Person has no retained Rater pair")
    quartile_summary = (
        person_assignment.groupby("AbilityQuartile", as_index=False)
        .agg(
            Persons=("Person", "size"),
            MeanAssignedSeverity=("AssignedMeanSeverity", "mean"),
            MeanTheta=("Theta", "mean"),
        )
    )
    pair_counts = (
        person_assignment.groupby(["AbilityQuartile", "RaterPair"])
        .size()
        .rename("Persons")
        .reset_index()
    )
    rater_person_counts = (
        output[["Person", "Rater"]].drop_duplicates().groupby("Rater").size().to_dict()
    )
    rater_row_counts = output.groupby("Rater").size().to_dict()
    person_rater_counts = output.groupby("Person")["Rater"].nunique()
    rank = constrained_adjacent_design_audit(output, FIT_MODEL)
    components = person_rater_components(output)
    if person_assignment["AssignedMeanSeverity"].nunique() < 2:
        correlation = np.nan
    else:
        correlation = spearmanr(
            person_assignment["Theta"], person_assignment["AssignedMeanSeverity"]
        ).statistic
    sparse = design != "complete"
    expected_rows = 960 if sparse else 1920
    invariants = {
        "Rows": int(len(output)) == expected_rows,
        "Persons": int(output["Person"].nunique()) == 80,
        "Raters": int(output["Rater"].nunique()) == 4,
        "RatersPerPerson": person_rater_counts.eq(2 if sparse else 4).all(),
        "PersonsPerRater": (
            set(map(int, rater_person_counts.values())) == ({40} if sparse else {80})
        ),
        "RatingRowsPerRater": (
            set(map(int, rater_row_counts.values())) == ({240} if sparse else {480})
        ),
        "PersonRaterComponents": len(components) == 1,
        "ConstrainedPCMNullity": int(rank["Nullity"]) == 0,
        "QuartilesBalanced": quartile_summary["Persons"].eq(20).all(),
    }
    if design == "ability_severity_aligned_connected":
        observed_pairs = person_assignment.groupby("AbilityQuartile")["RaterPair"].unique()
        invariants["RegisteredPairByQuartile"] = all(
            list(values) == ["|".join(ALIGNED_PAIRS[index - 1])]
            for index, values in observed_pairs.items()
        )
        invariants["RegisteredMeanSeverityByQuartile"] = np.allclose(
            quartile_summary["MeanAssignedSeverity"], EXPECTED_PAIR_MEANS, atol=1e-12
        )
    return {
        "Design": design,
        "Rows": int(len(output)),
        "PersonRaterComponents": len(components),
        "ConstrainedPCMRank": int(rank["Rank"]),
        "ConstrainedPCMNullity": int(rank["Nullity"]),
        "RaterPersonCounts": {str(key): int(value) for key, value in rater_person_counts.items()},
        "RaterRowCounts": {str(key): int(value) for key, value in rater_row_counts.items()},
        "ThetaAssignedSeveritySpearman": float(correlation),
        "QuartileSummary": quartile_summary.to_dict(orient="records"),
        "QuartilePairCounts": pair_counts.to_dict(orient="records"),
        "Invariants": {key: bool(value) for key, value in invariants.items()},
        "AllInvariantsPass": all(bool(value) for value in invariants.values()),
    }


__all__ = [
    "ALIGNED_PAIRS",
    "DESIGNS",
    "EXPECTED_PAIR_MEANS",
    "ability_quartiles",
    "apply_informative_assignment_design",
    "informative_assignment_audit",
]
