"""Person-coordinate generator contracts for assignment sensitivity.

The assignment design is constructed from a source fitted Person ordering.
Two response-generation questions are deliberately separated:

``source_fitted``
    Hold the fitted Person measures fixed.  For MML these are posterior/EAP
    summaries and are not known truth.

``mml_population_rank_preserving``
    Draw a fresh iid normal Person population using the fitted MML population
    SD, then map its order statistics to the strictly ordered source EAP
    ranks.  This preserves the rank-to-assignment relation used to construct
    the counterfactual while varying population-scale realizations.  It is not
    an estimate of the assignment mechanism or an unconditional new-sample
    experiment.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np
import pandas as pd


GENERATOR_SCHEMA_VERSION = "assignment_person_generator_v1"
SOURCE_FITTED_MODE = "source_fitted"
MML_RANK_PRESERVING_MODE = "mml_population_rank_preserving"
SUPPORTED_MODES = frozenset({SOURCE_FITTED_MODE, MML_RANK_PRESERVING_MODE})


def build_assignment_person_generation_plan(
    *,
    method: str,
    person_scores: Mapping[str, float],
    mode: str = SOURCE_FITTED_MODE,
    population_sd: float | None = None,
    population_model_enabled: bool = False,
) -> dict[str, object]:
    """Validate and document one Person-coordinate generation mode."""

    method = str(method).upper()
    mode = str(mode).lower()
    if mode not in SUPPORTED_MODES:
        raise ValueError(
            "mode must be 'source_fitted' or 'mml_population_rank_preserving'."
        )
    clean_scores = {str(key): float(value) for key, value in person_scores.items()}
    people = sorted(clean_scores)
    values = np.asarray([clean_scores[person] for person in people], dtype=float)
    gate_rows: list[dict[str, object]] = []

    def gate(name: str, passed: bool, evidence: str, action: str) -> None:
        gate_rows.append(
            {
                "Gate": name,
                "Passed": bool(passed),
                "Evidence": evidence,
                "ActionIfFailed": action,
            }
        )

    finite = bool(len(values) >= 2 and np.isfinite(values).all())
    gate(
        "Finite source Person coordinates",
        finite,
        f"finite={int(np.isfinite(values).sum())}/{len(values)}",
        "Resolve missing/nonreportable Person coordinates before response generation.",
    )
    gate(
        "Supported source estimator",
        method in {"JMLE", "MML"},
        f"method={method or 'not available'}",
        "Use a qualified JMLE or MML source fit.",
    )

    is_population_mode = mode == MML_RANK_PRESERVING_MODE
    gate(
        "Generator-estimator compatibility",
        (method == "MML") if is_population_mode else method in {"JMLE", "MML"},
        f"mode={mode}; method={method}",
        "Rank-preserving population generation is available only for MML.",
    )
    gate(
        "No latent-regression population model",
        (not population_model_enabled) if is_population_mode else True,
        f"population_model_enabled={bool(population_model_enabled)}",
        "A covariate-stratified generator contract is required for latent regression.",
    )
    sd_value = pd.to_numeric(pd.Series([population_sd]), errors="coerce").iloc[0]
    sd_ready = bool(np.isfinite(sd_value) and float(sd_value) > 0)
    gate(
        "Finite positive MML population SD",
        sd_ready if is_population_mode else True,
        f"population_sd={sd_value if np.isfinite(sd_value) else 'not available'}",
        "Fit MML with a finite positive fixed or estimated population SD.",
    )
    distinct = int(pd.Series(values).nunique(dropna=True)) if finite else 0
    strict_order = bool(finite and distinct == len(values))
    sorted_values = np.sort(values) if finite else np.array([], dtype=float)
    minimum_gap = (
        float(np.min(np.diff(sorted_values))) if len(sorted_values) >= 2 else np.nan
    )
    gate(
        "Strict source Person ordering",
        strict_order if is_population_mode else True,
        f"distinct={distinct}/{len(values)}; minimum_adjacent_gap={minimum_gap}",
        "Resolve tied source EAP coordinates; their within-tie assignment rank is unidentified.",
    )

    gates = pd.DataFrame(gate_rows)
    available = bool(not gates.empty and gates["Passed"].all())
    failed = gates.loc[~gates["Passed"], "Gate"].astype(str).tolist()
    contract = pd.DataFrame(
        [
            {
                "SchemaVersion": GENERATOR_SCHEMA_VERSION,
                "Mode": mode,
                "Method": method,
                "AssignmentCoordinate": (
                    "source fitted Person rank; MML uses posterior/EAP rank"
                    if method == "MML"
                    else "source fitted JMLE Person rank"
                ),
                "ResponsePersonCoordinate": (
                    "normal order statistics mapped to source EAP ranks"
                    if is_population_mode
                    else "source fitted Person measure"
                ),
                "PopulationSD": float(sd_value) if is_population_mode and sd_ready else np.nan,
                "GeneratedPersonTruthKnownWithinReplicate": is_population_mode,
                "SourceFittedPersonIsKnownTruth": False,
                "RankRelationPreserved": is_population_mode,
                "CommonPersonCoordinateAcrossScenariosWithinReplicate": True,
                "UnconditionalNewSample": False,
                "AssignmentMechanismIdentified": False,
                "AllowedClaim": (
                    "MML normal-population rank-preserving fitted-assignment sensitivity"
                    if is_population_mode
                    else "local fitted-Person conditional sensitivity"
                ),
            }
        ]
    )
    ordered_people = [
        person
        for person, _ in sorted(clean_scores.items(), key=lambda item: (item[1], item[0]))
    ]
    return {
        "available": available,
        "reason": (
            "The requested Person generator is ready."
            if available
            else "Blocked by: " + "; ".join(failed)
        ),
        "schema_version": GENERATOR_SCHEMA_VERSION,
        "mode": mode,
        "method": method,
        "population_sd": float(sd_value) if sd_ready else np.nan,
        "gates": gates,
        "contract": contract,
        "ordered_people": ordered_people,
        "source_coordinates": clean_scores,
    }


def draw_assignment_person_coordinates(
    plan: Mapping[str, object],
    *,
    rng: np.random.Generator,
    replicate: int,
) -> dict[str, object]:
    """Draw one auditable coordinate vector for a validated plan."""

    if not isinstance(plan, Mapping) or not plan.get("available"):
        raise ValueError("An available Person-generation plan is required.")
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng must be a numpy.random.Generator.")
    mode = str(plan.get("mode"))
    ordered_people = [str(person) for person in plan.get("ordered_people", [])]
    source = {
        str(key): float(value)
        for key, value in dict(plan.get("source_coordinates", {})).items()
    }
    if mode == SOURCE_FITTED_MODE:
        coordinates = dict(source)
        random_population = False
    elif mode == MML_RANK_PRESERVING_MODE:
        population_sd = float(plan.get("population_sd", np.nan))
        draws = np.sort(rng.normal(0.0, population_sd, size=len(ordered_people)))
        coordinates = {
            person: float(draw) for person, draw in zip(ordered_people, draws)
        }
        random_population = True
    else:
        raise ValueError(f"Unsupported Person-generation mode: {mode}")

    frame = pd.DataFrame(
        {
            "Replicate": int(replicate),
            "Person": ordered_people,
            "SourceCoordinate": [source[person] for person in ordered_people],
            "GeneratedCoordinate": [coordinates[person] for person in ordered_people],
        }
    )
    frame["SourceRank"] = frame["SourceCoordinate"].rank(method="average")
    frame["GeneratedRank"] = frame["GeneratedCoordinate"].rank(method="average")
    rank_correlation = float(
        frame["SourceCoordinate"].corr(frame["GeneratedCoordinate"], method="spearman")
    )
    generated = frame["GeneratedCoordinate"].to_numpy(dtype=float)
    summary = pd.DataFrame(
        [
            {
                "Replicate": int(replicate),
                "Mode": mode,
                "Persons": len(frame),
                "RandomPopulationDraw": random_population,
                "GeneratorPopulationSD": (
                    float(plan.get("population_sd", np.nan))
                    if random_population
                    else np.nan
                ),
                "GeneratedMean": float(np.mean(generated)),
                "GeneratedSD": float(np.std(generated, ddof=1)) if len(generated) > 1 else np.nan,
                "GeneratedMinimum": float(np.min(generated)),
                "GeneratedMaximum": float(np.max(generated)),
                "SourceGeneratedRankSpearman": rank_correlation,
                "RankPreservedExactly": bool(np.isclose(rank_correlation, 1.0)),
                "CommonAcrossAssignmentScenarios": True,
            }
        ]
    )
    return {
        "coordinates": coordinates,
        "draws": frame,
        "summary": summary,
    }
