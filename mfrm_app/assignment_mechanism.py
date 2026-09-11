"""Known-mechanism assignment generators for repository validation.

This module defines an exponential-family distribution over bipartite
Person-Rater graphs with fixed Person and Rater degrees::

    Pr(G | degrees, connected) proportional to exp(gamma * T(G))

where ``T(G)`` is the sum of products of standardized Person ability and
Rater severity coordinates over observed edges.  ``gamma=0`` is uniform over
the conditioned state space; positive gamma favors ability-severity alignment.

Small state spaces can be enumerated exactly.  Larger spaces use a symmetric
Metropolis degree-preserving 2-switch kernel.  This is repository-only
validation infrastructure: it generates a known conditional assignment
mechanism and does not estimate a mechanism from observed ratings.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import combinations
import time
from typing import Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.special import logsumexp


MECHANISM_SCHEMA_VERSION = "degree_conditioned_assignment_mechanism_v1"
MAX_ENUMERATION_STATES = 200_000
MAX_ENUMERATION_NODES = 2_000_000
MAX_CHAIN_STEPS = 5_000_000
MAX_DP_STATES = 2_000_000
MAX_EXACT_INDEPENDENT_SAMPLES = 100_000
MAX_DENSE_DP_CELLS = 20_000_000


def _standardized_coordinates(
    values: Mapping[str, float], levels: list[str], label: str
) -> dict[str, float]:
    series = pd.Series({level: values.get(level, np.nan) for level in levels}, dtype=float)
    finite = np.isfinite(series.to_numpy(dtype=float))
    if not finite.all():
        missing = series.index[~finite].astype(str).tolist()
        raise ValueError(f"Non-finite {label} coordinate for: " + ", ".join(missing[:10]))
    centered = series - float(series.mean())
    sd = float(np.sqrt(np.mean(centered.to_numpy(dtype=float) ** 2)))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError(f"{label} coordinates must contain at least two distinct values.")
    standardized = centered / sd
    return {level: float(standardized.loc[level]) for level in levels}


def _direct_overlap_components(edges: set[tuple[str, str]], raters: list[str]) -> int:
    if not raters:
        return 0
    people_by_rater = {
        rater: {person for person, candidate in edges if candidate == rater}
        for rater in raters
    }
    adjacency = {rater: set() for rater in raters}
    for left, right in combinations(raters, 2):
        if people_by_rater[left] & people_by_rater[right]:
            adjacency[left].add(right)
            adjacency[right].add(left)
    remaining = set(raters)
    components = 0
    while remaining:
        components += 1
        stack = [remaining.pop()]
        while stack:
            current = stack.pop()
            unseen = adjacency[current] & remaining
            remaining.difference_update(unseen)
            stack.extend(unseen)
    return components


def _edge_signature(edges: set[tuple[str, str]]) -> str:
    return "|".join(f"{person}::{rater}" for person, rater in sorted(edges))


def _assignment_statistic(
    edges: set[tuple[str, str]],
    person_z: Mapping[str, float],
    rater_z: Mapping[str, float],
) -> float:
    return float(sum(person_z[person] * rater_z[rater] for person, rater in edges))


def _assignment_correlation(
    edges: set[tuple[str, str]],
    person_z: Mapping[str, float],
    rater_z: Mapping[str, float],
) -> float:
    ordered = sorted(edges)
    left = np.asarray([person_z[person] for person, _ in ordered], dtype=float)
    right = np.asarray([rater_z[rater] for _, rater in ordered], dtype=float)
    if len(ordered) < 2 or np.std(left) <= 0 or np.std(right) <= 0:
        return float("nan")
    return float(np.corrcoef(left, right)[0, 1])


def _prepare_degree_contract(
    *,
    person_degrees: Mapping[str, int],
    rater_degrees: Mapping[str, int],
    person_coordinates: Mapping[str, float],
    rater_coordinates: Mapping[str, float],
    gamma: float,
) -> dict[str, object]:
    people = sorted(str(level) for level in person_degrees)
    raters = sorted(str(level) for level in rater_degrees)
    p_degree = {str(key): int(value) for key, value in person_degrees.items()}
    r_degree = {str(key): int(value) for key, value in rater_degrees.items()}
    if len(people) < 2 or len(raters) < 2:
        raise ValueError("At least two Persons and two Raters are required.")
    if any(value < 1 or value > len(raters) for value in p_degree.values()):
        raise ValueError("Each Person degree must lie in [1, number of Raters].")
    if any(value < 1 or value > len(people) for value in r_degree.values()):
        raise ValueError("Each Rater degree must lie in [1, number of Persons].")
    edges_person = int(sum(p_degree.values()))
    edges_rater = int(sum(r_degree.values()))
    if edges_person != edges_rater:
        raise ValueError(
            f"Person and Rater degree totals differ: {edges_person} != {edges_rater}."
        )
    gamma = float(gamma)
    if not np.isfinite(gamma):
        raise ValueError("gamma must be finite.")
    person_z = _standardized_coordinates(person_coordinates, people, "Person")
    rater_z = _standardized_coordinates(rater_coordinates, raters, "Rater")
    return {
        "people": people,
        "raters": raters,
        "person_degrees": p_degree,
        "rater_degrees": r_degree,
        "person_z": person_z,
        "rater_z": rater_z,
        "gamma": gamma,
        "edge_count": edges_person,
    }


def enumerate_degree_conditioned_assignments(
    *,
    person_degrees: Mapping[str, int],
    rater_degrees: Mapping[str, int],
    person_coordinates: Mapping[str, float],
    rater_coordinates: Mapping[str, float],
    gamma: float,
    require_connected: bool = True,
    max_states: int = MAX_ENUMERATION_STATES,
    max_nodes: int = MAX_ENUMERATION_NODES,
) -> dict[str, object]:
    """Enumerate and normalize a small fixed-degree assignment state space."""

    contract = _prepare_degree_contract(
        person_degrees=person_degrees,
        rater_degrees=rater_degrees,
        person_coordinates=person_coordinates,
        rater_coordinates=rater_coordinates,
        gamma=gamma,
    )
    people = contract["people"]
    raters = contract["raters"]
    p_degree = contract["person_degrees"]
    initial_remaining = tuple(contract["rater_degrees"][rater] for rater in raters)
    states: list[set[tuple[str, str]]] = []
    visited_nodes = 0

    def visit(
        person_index: int,
        remaining: tuple[int, ...],
        edges: set[tuple[str, str]],
    ) -> None:
        nonlocal visited_nodes
        visited_nodes += 1
        if visited_nodes > int(max_nodes):
            raise RuntimeError(
                f"Enumeration exceeded max_nodes={int(max_nodes)} before completion."
            )
        if person_index == len(people):
            if any(remaining):
                return
            if require_connected and _direct_overlap_components(edges, raters) != 1:
                return
            states.append(set(edges))
            if len(states) > int(max_states):
                raise RuntimeError(
                    f"Enumeration exceeded max_states={int(max_states)}."
                )
            return

        person = people[person_index]
        degree = p_degree[person]
        remaining_people = len(people) - person_index - 1
        available = [index for index, value in enumerate(remaining) if value > 0]
        for selected in combinations(available, degree):
            next_remaining = list(remaining)
            for index in selected:
                next_remaining[index] -= 1
            if any(value < 0 or value > remaining_people for value in next_remaining):
                continue
            next_edges = set(edges)
            next_edges.update((person, raters[index]) for index in selected)
            visit(person_index + 1, tuple(next_remaining), next_edges)

    visit(0, initial_remaining, set())
    if not states:
        raise ValueError("No assignment graph satisfies the requested degree/connectivity contract.")

    rows: list[dict[str, object]] = []
    for state_id, edges in enumerate(states, start=1):
        statistic = _assignment_statistic(
            edges, contract["person_z"], contract["rater_z"]
        )
        rows.append(
            {
                "StateId": state_id,
                "EdgeSignature": _edge_signature(edges),
                "Statistic": statistic,
                "AssignmentCorrelation": _assignment_correlation(
                    edges, contract["person_z"], contract["rater_z"]
                ),
                "LogWeight": contract["gamma"] * statistic,
                "Connected": _direct_overlap_components(edges, raters) == 1,
            }
        )
    table = pd.DataFrame(rows).sort_values("EdgeSignature").reset_index(drop=True)
    table["StateId"] = np.arange(1, len(table) + 1, dtype=int)
    log_normalizer = float(logsumexp(table["LogWeight"].to_numpy(dtype=float)))
    table["Probability"] = np.exp(table["LogWeight"] - log_normalizer)
    expected_statistic = float(np.dot(table["Probability"], table["Statistic"]))
    expected_correlation = float(
        np.dot(table["Probability"], table["AssignmentCorrelation"])
    )
    variance = float(
        np.dot(
            table["Probability"],
            (table["Statistic"] - expected_statistic) ** 2,
        )
    )
    summary = pd.DataFrame(
        [
            {
                "SchemaVersion": MECHANISM_SCHEMA_VERSION,
                "Gamma": contract["gamma"],
                "Persons": len(people),
                "Raters": len(raters),
                "Edges": contract["edge_count"],
                "RequireConnected": bool(require_connected),
                "States": len(table),
                "EnumerationNodes": visited_nodes,
                "LogNormalizer": log_normalizer,
                "ExpectedStatistic": expected_statistic,
                "SDStatistic": float(np.sqrt(max(variance, 0.0))),
                "ExpectedAssignmentCorrelation": expected_correlation,
                "ProbabilitySum": float(table["Probability"].sum()),
            }
        ]
    )
    return {
        "available": True,
        "reason": "The conditioned assignment state space was enumerated exactly.",
        "schema_version": MECHANISM_SCHEMA_VERSION,
        "states": table,
        "summary": summary,
        "contract": contract,
        "claim_boundary": (
            "Probabilities are exact only for the declared fixed-degree state space and optional direct-overlap "
            "connectivity condition. Gamma is a known generator setting, not estimated from observed assignments."
        ),
    }


def _initial_positive_sequence_ess(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    n = len(values)
    if n < 3 or not np.isfinite(values).all():
        return np.nan
    centered = values - float(np.mean(values))
    variance = float(np.dot(centered, centered) / n)
    if variance <= 0:
        return float(n)
    rho_sum = 0.0
    for lag in range(1, min(n // 2, 2_000) + 1):
        covariance = float(np.dot(centered[:-lag], centered[lag:]) / (n - lag))
        rho = covariance / variance
        if not np.isfinite(rho) or rho <= 0:
            break
        rho_sum += rho
    tau = max(1.0, 1.0 + 2.0 * rho_sum)
    return float(n / tau)


def sample_degree_conditioned_assignment_chain(
    *,
    initial_edges: Iterable[tuple[str, str]],
    person_coordinates: Mapping[str, float],
    rater_coordinates: Mapping[str, float],
    gamma: float,
    n_samples: int,
    burnin: int = 2_000,
    thin: int = 10,
    seed: int = 20260811,
    require_connected: bool = True,
) -> dict[str, object]:
    """Sample fixed-degree graphs using a symmetric Metropolis 2-switch."""

    edges = {(str(person), str(rater)) for person, rater in initial_edges}
    if not edges:
        raise ValueError("initial_edges must not be empty.")
    person_degrees = pd.Series([person for person, _ in edges]).value_counts().to_dict()
    rater_degrees = pd.Series([rater for _, rater in edges]).value_counts().to_dict()
    contract = _prepare_degree_contract(
        person_degrees=person_degrees,
        rater_degrees=rater_degrees,
        person_coordinates=person_coordinates,
        rater_coordinates=rater_coordinates,
        gamma=gamma,
    )
    people = contract["people"]
    raters = contract["raters"]
    if {person for person, _ in edges} != set(people):
        raise ValueError("initial_edges do not cover every declared Person coordinate.")
    if {rater for _, rater in edges} != set(raters):
        raise ValueError("initial_edges do not cover every declared Rater coordinate.")
    if require_connected and _direct_overlap_components(edges, raters) != 1:
        raise ValueError("initial_edges are not connected under direct Rater overlap.")
    n_samples = int(n_samples)
    burnin = int(burnin)
    thin = int(thin)
    if n_samples < 1 or burnin < 0 or thin < 1:
        raise ValueError("n_samples>=1, burnin>=0, and thin>=1 are required.")
    total_steps = burnin + n_samples * thin
    if total_steps > MAX_CHAIN_STEPS:
        raise ValueError(
            f"Requested chain steps {total_steps} exceed {MAX_CHAIN_STEPS}."
        )

    rng = np.random.default_rng(int(seed))
    statistic = _assignment_statistic(
        edges, contract["person_z"], contract["rater_z"]
    )
    proposals = accepted = invalid = disconnected = 0
    sample_rows: list[dict[str, object]] = []
    for step in range(1, total_steps + 1):
        proposals += 1
        ordered_edges = sorted(edges)
        selected = rng.choice(len(ordered_edges), size=2, replace=False)
        (p1, r1), (p2, r2) = ordered_edges[int(selected[0])], ordered_edges[int(selected[1])]
        if p1 == p2 or r1 == r2 or (p1, r2) in edges or (p2, r1) in edges:
            invalid += 1
        else:
            candidate = set(edges)
            candidate.remove((p1, r1))
            candidate.remove((p2, r2))
            candidate.add((p1, r2))
            candidate.add((p2, r1))
            if require_connected and _direct_overlap_components(candidate, raters) != 1:
                disconnected += 1
            else:
                delta = (
                    contract["person_z"][p1] * contract["rater_z"][r2]
                    + contract["person_z"][p2] * contract["rater_z"][r1]
                    - contract["person_z"][p1] * contract["rater_z"][r1]
                    - contract["person_z"][p2] * contract["rater_z"][r2]
                )
                log_acceptance = min(0.0, contract["gamma"] * float(delta))
                if np.log(rng.random()) < log_acceptance:
                    edges = candidate
                    statistic += float(delta)
                    accepted += 1
        if step > burnin and (step - burnin) % thin == 0:
            # Recompute rather than trusting accumulated deltas; the residual
            # is retained as a numerical audit of long chains.
            recomputed = _assignment_statistic(
                edges, contract["person_z"], contract["rater_z"]
            )
            residual = float(statistic - recomputed)
            statistic = recomputed
            sample_rows.append(
                {
                    "Sample": len(sample_rows) + 1,
                    "Step": step,
                    "Statistic": statistic,
                    "AssignmentCorrelation": _assignment_correlation(
                        edges, contract["person_z"], contract["rater_z"]
                    ),
                    "EdgeSignature": _edge_signature(edges),
                    "StatisticUpdateResidual": residual,
                }
            )
    samples = pd.DataFrame(sample_rows)
    diagnostics = pd.DataFrame(
        [
            {
                "SchemaVersion": MECHANISM_SCHEMA_VERSION,
                "Gamma": contract["gamma"],
                "Seed": int(seed),
                "Burnin": burnin,
                "Thin": thin,
                "Samples": len(samples),
                "TotalSteps": total_steps,
                "Proposals": proposals,
                "Accepted": accepted,
                "InvalidSwitchProposals": invalid,
                "DisconnectedSwitchProposals": disconnected,
                "AcceptanceRateAllProposals": accepted / proposals,
                "AcceptanceRateAdmissibleProposals": (
                    accepted / (proposals - invalid - disconnected)
                    if proposals > invalid + disconnected else np.nan
                ),
                "UniqueSampledStates": int(samples["EdgeSignature"].nunique()),
                "StatisticESS": _initial_positive_sequence_ess(
                    samples["Statistic"].to_numpy(dtype=float)
                ),
                "MaximumStatisticUpdateResidual": float(
                    samples["StatisticUpdateResidual"].abs().max()
                ),
                "FinalOverlapComponents": _direct_overlap_components(edges, raters),
                "DegreeMarginsPreserved": True,
                "ProposalKernel": "uniform unordered edge-pair symmetric 2-switch",
            }
        ]
    )
    available = bool(
        len(samples) == n_samples
        and diagnostics.iloc[0]["MaximumStatisticUpdateResidual"] <= 1e-10
        and (not require_connected or diagnostics.iloc[0]["FinalOverlapComponents"] == 1)
    )
    return {
        "available": available,
        "reason": (
            "The requested conditioned assignment chain completed."
            if available else "The conditioned assignment chain failed its numerical/invariant audit."
        ),
        "schema_version": MECHANISM_SCHEMA_VERSION,
        "samples": samples,
        "diagnostics": diagnostics,
        "final_edges": set(edges),
        "contract": contract,
        "claim_boundary": (
            "Metropolis samples target the declared fixed-degree, optionally connected exponential-family "
            "assignment distribution. Finite-chain mixing must be checked against an exact oracle or multiple-chain "
            "diagnostics before performance claims. Gamma is known by construction, not estimated."
        ),
    }


def sample_degree_conditioned_assignment_heatbath_chain(
    *,
    initial_edges: Iterable[tuple[str, str]],
    person_coordinates: Mapping[str, float],
    rater_coordinates: Mapping[str, float],
    gamma: float,
    n_samples: int,
    burnin: int = 2_000,
    thin: int = 10,
    seed: int = 20260811,
    require_connected: bool = True,
) -> dict[str, object]:
    """Random-scan Person-pair heat-bath sampler for the same target law.

    For a uniformly selected pair of Persons, common Raters remain fixed and
    the union of their exclusive Raters is repartitioned while each Person's
    degree is held fixed.  Every connected repartition is then sampled from its
    exact conditional exponential-family probability.  This Curveball-style
    block update is a Gibbs kernel, not a different assignment estimand.
    """

    edges = {(str(person), str(rater)) for person, rater in initial_edges}
    if not edges:
        raise ValueError("initial_edges must not be empty.")
    person_degrees = pd.Series([person for person, _ in edges]).value_counts().to_dict()
    rater_degrees = pd.Series([rater for _, rater in edges]).value_counts().to_dict()
    contract = _prepare_degree_contract(
        person_degrees=person_degrees,
        rater_degrees=rater_degrees,
        person_coordinates=person_coordinates,
        rater_coordinates=rater_coordinates,
        gamma=gamma,
    )
    people = contract["people"]
    raters = contract["raters"]
    if {person for person, _ in edges} != set(people):
        raise ValueError("initial_edges do not cover every declared Person coordinate.")
    if {rater for _, rater in edges} != set(raters):
        raise ValueError("initial_edges do not cover every declared Rater coordinate.")
    if require_connected and _direct_overlap_components(edges, raters) != 1:
        raise ValueError("initial_edges are not connected under direct Rater overlap.")
    n_samples = int(n_samples)
    burnin = int(burnin)
    thin = int(thin)
    if n_samples < 1 or burnin < 0 or thin < 1:
        raise ValueError("n_samples>=1, burnin>=0, and thin>=1 are required.")
    total_steps = burnin + n_samples * thin
    if total_steps > MAX_CHAIN_STEPS:
        raise ValueError(
            f"Requested chain steps {total_steps} exceed {MAX_CHAIN_STEPS}."
        )

    rng = np.random.default_rng(int(seed))
    statistic = _assignment_statistic(edges, contract["person_z"], contract["rater_z"])
    updates = nontrivial_updates = moved_updates = 0
    disconnected_candidates = candidate_states_total = 0
    sample_rows: list[dict[str, object]] = []
    for step in range(1, total_steps + 1):
        updates += 1
        selected_people = rng.choice(len(people), size=2, replace=False)
        p1 = people[int(selected_people[0])]
        p2 = people[int(selected_people[1])]
        neighbors1 = {rater for person, rater in edges if person == p1}
        neighbors2 = {rater for person, rater in edges if person == p2}
        common = neighbors1 & neighbors2
        exclusive1 = neighbors1 - common
        exclusive2 = neighbors2 - common
        union = sorted(exclusive1 | exclusive2)
        old_pair_contribution = float(
            sum(contract["person_z"][p1] * contract["rater_z"][rater] for rater in neighbors1)
            + sum(contract["person_z"][p2] * contract["rater_z"][rater] for rater in neighbors2)
        )
        candidates: list[tuple[set[tuple[str, str]], float]] = []
        for selected_raters in combinations(union, len(exclusive1)):
            selected_set = set(selected_raters)
            next_neighbors1 = common | selected_set
            next_neighbors2 = common | (set(union) - selected_set)
            candidate = {
                edge for edge in edges if edge[0] not in {p1, p2}
            }
            candidate.update((p1, rater) for rater in next_neighbors1)
            candidate.update((p2, rater) for rater in next_neighbors2)
            if require_connected and _direct_overlap_components(candidate, raters) != 1:
                disconnected_candidates += 1
                continue
            next_pair_contribution = float(
                sum(
                    contract["person_z"][p1] * contract["rater_z"][rater]
                    for rater in next_neighbors1
                )
                + sum(
                    contract["person_z"][p2] * contract["rater_z"][rater]
                    for rater in next_neighbors2
                )
            )
            candidates.append(
                (candidate, statistic + next_pair_contribution - old_pair_contribution)
            )
        if not candidates:
            raise RuntimeError("Person-pair heat-bath update has no connected candidate.")
        candidate_states_total += len(candidates)
        if len(candidates) > 1:
            nontrivial_updates += 1
        log_weights = np.asarray(
            [contract["gamma"] * value for _, value in candidates], dtype=float
        )
        probabilities = np.exp(log_weights - logsumexp(log_weights))
        selected_index = int(rng.choice(len(candidates), p=probabilities))
        selected_edges, selected_statistic = candidates[selected_index]
        if selected_edges != edges:
            moved_updates += 1
        edges = selected_edges
        statistic = float(selected_statistic)
        if step > burnin and (step - burnin) % thin == 0:
            recomputed = _assignment_statistic(
                edges, contract["person_z"], contract["rater_z"]
            )
            residual = float(statistic - recomputed)
            statistic = recomputed
            sample_rows.append(
                {
                    "Sample": len(sample_rows) + 1,
                    "Step": step,
                    "Statistic": statistic,
                    "AssignmentCorrelation": _assignment_correlation(
                        edges, contract["person_z"], contract["rater_z"]
                    ),
                    "EdgeSignature": _edge_signature(edges),
                    "StatisticUpdateResidual": residual,
                }
            )
    samples = pd.DataFrame(sample_rows)
    diagnostics = pd.DataFrame(
        [
            {
                "SchemaVersion": MECHANISM_SCHEMA_VERSION,
                "Gamma": contract["gamma"],
                "Seed": int(seed),
                "Burnin": burnin,
                "Thin": thin,
                "Samples": len(samples),
                "TotalSteps": total_steps,
                "Updates": updates,
                "NontrivialUpdates": nontrivial_updates,
                "MovedUpdates": moved_updates,
                "MovementRateAllUpdates": moved_updates / updates,
                "MovementRateNontrivialUpdates": (
                    moved_updates / nontrivial_updates if nontrivial_updates else np.nan
                ),
                "MeanConnectedCandidateStates": candidate_states_total / updates,
                "DisconnectedCandidateStates": disconnected_candidates,
                "UniqueSampledStates": int(samples["EdgeSignature"].nunique()),
                "StatisticESS": _initial_positive_sequence_ess(
                    samples["Statistic"].to_numpy(dtype=float)
                ),
                "MaximumStatisticUpdateResidual": float(
                    samples["StatisticUpdateResidual"].abs().max()
                ),
                "FinalOverlapComponents": _direct_overlap_components(edges, raters),
                "DegreeMarginsPreserved": True,
                "ProposalKernel": "uniform Person-pair connected heat-bath Curveball repartition",
            }
        ]
    )
    available = bool(
        len(samples) == n_samples
        and diagnostics.iloc[0]["MaximumStatisticUpdateResidual"] <= 1e-10
        and (not require_connected or diagnostics.iloc[0]["FinalOverlapComponents"] == 1)
    )
    return {
        "available": available,
        "reason": (
            "The requested conditioned assignment heat-bath chain completed."
            if available
            else "The conditioned assignment heat-bath chain failed its numerical/invariant audit."
        ),
        "schema_version": MECHANISM_SCHEMA_VERSION,
        "samples": samples,
        "diagnostics": diagnostics,
        "final_edges": set(edges),
        "contract": contract,
        "claim_boundary": (
            "This random-scan Person-pair heat-bath kernel targets the same declared fixed-degree, optionally "
            "connected assignment distribution. Mixing still requires an exact oracle or multiple-chain "
            "qualification; gamma remains known by construction rather than estimated."
        ),
    }


def sample_degree_conditioned_assignment_dp(
    *,
    person_degrees: Mapping[str, int],
    rater_degrees: Mapping[str, int],
    person_coordinates: Mapping[str, float],
    rater_coordinates: Mapping[str, float],
    gamma: float,
    n_samples: int,
    seed: int = 20260811,
    require_connected: bool = True,
    max_dp_states: int = MAX_DP_STATES,
    max_total_rejections: int = 100_000,
    retain_edge_signatures: bool = True,
) -> dict[str, object]:
    """Draw independent fixed-margin graphs using exact dynamic programming.

    The backward recursion integrates all remaining Rater-degree allocations.
    Sequential draws therefore follow the exact unconditioned fixed-margin
    exponential-family law.  Rejecting disconnected draws yields independent
    samples from that law conditioned on direct Rater-overlap connectivity.
    This route is intended for a small number of Rater levels; the DP state cap
    fails closed when the remaining-margin state space becomes too large.
    """

    contract = _prepare_degree_contract(
        person_degrees=person_degrees,
        rater_degrees=rater_degrees,
        person_coordinates=person_coordinates,
        rater_coordinates=rater_coordinates,
        gamma=gamma,
    )
    people = contract["people"]
    raters = contract["raters"]
    p_degree = contract["person_degrees"]
    initial_remaining = tuple(contract["rater_degrees"][rater] for rater in raters)
    n_samples = int(n_samples)
    max_dp_states = int(max_dp_states)
    max_total_rejections = int(max_total_rejections)
    if n_samples < 1 or n_samples > MAX_EXACT_INDEPENDENT_SAMPLES:
        raise ValueError(
            f"n_samples must lie in [1, {MAX_EXACT_INDEPENDENT_SAMPLES}]."
        )
    if max_dp_states < 1 or max_total_rejections < 0:
        raise ValueError("Positive max_dp_states and nonnegative max_total_rejections are required.")

    suffix_edges = np.zeros(len(people) + 1, dtype=int)
    for index in range(len(people) - 1, -1, -1):
        suffix_edges[index] = suffix_edges[index + 1] + p_degree[people[index]]
    states_evaluated = 0

    @lru_cache(maxsize=None)
    def log_partition(person_index: int, remaining: tuple[int, ...]) -> float:
        nonlocal states_evaluated
        states_evaluated += 1
        if states_evaluated > max_dp_states:
            raise RuntimeError(
                f"Dynamic program exceeded max_dp_states={max_dp_states}."
            )
        if sum(remaining) != int(suffix_edges[person_index]):
            return float("-inf")
        remaining_people = len(people) - person_index
        if any(value < 0 or value > remaining_people for value in remaining):
            return float("-inf")
        if person_index == len(people):
            return 0.0 if not any(remaining) else float("-inf")
        person = people[person_index]
        degree = p_degree[person]
        available = [index for index, value in enumerate(remaining) if value > 0]
        values: list[float] = []
        for selected in combinations(available, degree):
            next_remaining = list(remaining)
            for rater_index in selected:
                next_remaining[rater_index] -= 1
            continuation = log_partition(person_index + 1, tuple(next_remaining))
            if not np.isfinite(continuation):
                continue
            local_statistic = contract["person_z"][person] * sum(
                contract["rater_z"][raters[rater_index]]
                for rater_index in selected
            )
            values.append(contract["gamma"] * local_statistic + continuation)
        return float(logsumexp(values)) if values else float("-inf")

    log_normalizer_unconditioned = log_partition(0, initial_remaining)
    if not np.isfinite(log_normalizer_unconditioned):
        raise ValueError("No graph satisfies the fixed Person/Rater degree contract.")

    @lru_cache(maxsize=None)
    def conditional_options(
        person_index: int, remaining: tuple[int, ...]
    ) -> tuple[tuple[tuple[int, ...], tuple[int, ...], float, float], ...]:
        """Cache exact full-conditionals for recurrent remaining-margin states."""

        person = people[person_index]
        degree = p_degree[person]
        available = [index for index, value in enumerate(remaining) if value > 0]
        candidates: list[tuple[tuple[int, ...], tuple[int, ...], float, float]] = []
        for selected in combinations(available, degree):
            next_remaining = list(remaining)
            for rater_index in selected:
                next_remaining[rater_index] -= 1
            continuation = log_partition(person_index + 1, tuple(next_remaining))
            if not np.isfinite(continuation):
                continue
            local_statistic = contract["person_z"][person] * sum(
                contract["rater_z"][raters[rater_index]]
                for rater_index in selected
            )
            candidates.append(
                (
                    selected,
                    tuple(next_remaining),
                    float(local_statistic),
                    float(contract["gamma"] * local_statistic + continuation),
                )
            )
        if not candidates:
            return ()
        log_weights = np.asarray([candidate[3] for candidate in candidates], dtype=float)
        probabilities = np.exp(log_weights - logsumexp(log_weights))
        return tuple(
            (selected, next_remaining, local_statistic, float(probability))
            for (selected, next_remaining, local_statistic, _), probability in zip(
                candidates, probabilities, strict=True
            )
        )

    edge_count = float(contract["edge_count"])
    person_edge_mean = float(
        sum(
            contract["person_z"][person] * p_degree[person]
            for person in people
        ) / edge_count
    )
    rater_edge_mean = float(
        sum(
            contract["rater_z"][rater] * contract["rater_degrees"][rater]
            for rater in raters
        ) / edge_count
    )
    person_edge_sd = float(
        np.sqrt(
            sum(
                p_degree[person]
                * (contract["person_z"][person] - person_edge_mean) ** 2
                for person in people
            ) / edge_count
        )
    )
    rater_edge_sd = float(
        np.sqrt(
            sum(
                contract["rater_degrees"][rater]
                * (contract["rater_z"][rater] - rater_edge_mean) ** 2
                for rater in raters
            ) / edge_count
        )
    )

    rng = np.random.default_rng(int(seed))
    sample_rows: list[dict[str, object]] = []
    connected_rejections = 0
    final_edges: set[tuple[str, str]] = set()
    maximum_residual = 0.0
    while len(sample_rows) < n_samples:
        remaining = initial_remaining
        edges: set[tuple[str, str]] = set()
        statistic = 0.0
        for person_index, person in enumerate(people):
            candidates = conditional_options(person_index, remaining)
            if not candidates:
                raise RuntimeError("Exact sequential sampler reached an empty conditional allocation.")
            draw = float(rng.random())
            cumulative = 0.0
            selected, remaining, local_statistic, _ = candidates[-1]
            for candidate in candidates:
                cumulative += candidate[3]
                if draw <= cumulative:
                    selected, remaining, local_statistic, _ = candidate
                    break
            edges.update((person, raters[rater_index]) for rater_index in selected)
            statistic += local_statistic
        if any(remaining):
            raise RuntimeError("Exact sequential sampler did not exhaust Rater margins.")
        if require_connected and _direct_overlap_components(edges, raters) != 1:
            connected_rejections += 1
            if connected_rejections > max_total_rejections:
                raise RuntimeError(
                    "Connected rejection sampling exceeded max_total_rejections."
                )
            continue
        recomputed = _assignment_statistic(edges, contract["person_z"], contract["rater_z"])
        residual = float(statistic - recomputed)
        maximum_residual = max(maximum_residual, abs(residual))
        row: dict[str, object] = {
            "Sample": len(sample_rows) + 1,
            "Statistic": recomputed,
            "AssignmentCorrelation": (
                (recomputed / edge_count - person_edge_mean * rater_edge_mean)
                / (person_edge_sd * rater_edge_sd)
            ),
            "StatisticUpdateResidual": residual,
        }
        if retain_edge_signatures:
            row["EdgeSignature"] = _edge_signature(edges)
        sample_rows.append(row)
        final_edges = edges

    samples = pd.DataFrame(sample_rows)
    total_attempts = n_samples + connected_rejections
    diagnostics = pd.DataFrame(
        [
            {
                "SchemaVersion": MECHANISM_SCHEMA_VERSION,
                "Gamma": contract["gamma"],
                "Seed": int(seed),
                "Samples": len(samples),
                "ExactIndependentSamples": True,
                "DPStatesEvaluated": states_evaluated,
                "CachedConditionalStates": conditional_options.cache_info().currsize,
                "DPStateCap": max_dp_states,
                "DPLogNormalizerUnconditioned": log_normalizer_unconditioned,
                "ConnectedRejections": connected_rejections,
                "TotalDrawAttempts": total_attempts,
                "ConnectedAcceptanceRate": n_samples / total_attempts,
                "UniqueSampledStates": (
                    int(samples["EdgeSignature"].nunique())
                    if "EdgeSignature" in samples else np.nan
                ),
                "StatisticESSDiagnostic": _initial_positive_sequence_ess(
                    samples["Statistic"].to_numpy(dtype=float)
                ),
                "MaximumStatisticUpdateResidual": maximum_residual,
                "FinalOverlapComponents": _direct_overlap_components(final_edges, raters),
                "DegreeMarginsPreserved": True,
                "SamplingKernel": "exact fixed-margin DP with connected rejection sampling",
            }
        ]
    )
    available = bool(
        len(samples) == n_samples
        and maximum_residual <= 1e-10
        and (not require_connected or diagnostics.iloc[0]["FinalOverlapComponents"] == 1)
    )
    return {
        "available": available,
        "reason": (
            "Exact independent fixed-margin assignment samples completed."
            if available
            else "Exact independent assignment samples failed their invariant audit."
        ),
        "schema_version": MECHANISM_SCHEMA_VERSION,
        "samples": samples,
        "diagnostics": diagnostics,
        "final_edges": final_edges,
        "contract": contract,
        "claim_boundary": (
            "Dynamic programming is exact for the declared fixed-margin exponential family before "
            "connectivity conditioning; rejection of disconnected draws yields exact independent conditional "
            "samples. The state-space cap limits this route to small Rater dimensions, and gamma remains known."
        ),
    }


def sample_degree_conditioned_assignment_dense_dp(
    *,
    person_degrees: Mapping[str, int],
    rater_degrees: Mapping[str, int],
    person_coordinates: Mapping[str, float],
    rater_coordinates: Mapping[str, float],
    gamma: float,
    n_samples: int,
    seed: int = 20260811,
    require_connected: bool = True,
    max_dense_dp_cells: int = MAX_DENSE_DP_CELLS,
    max_total_rejections: int = 100_000,
    retain_edge_signatures: bool = True,
) -> dict[str, object]:
    """Draw exact independent graphs using a vectorized dense-margin DP.

    The last Rater margin is implied by the remaining edge total, so each
    backward layer has one axis for every other Rater.  Vectorized shifts add
    every feasible Rater subset for the current Person.  The target law and
    connected rejection step are identical to
    :func:`sample_degree_conditioned_assignment_dp`; this implementation only
    changes how the partition function is evaluated.
    """

    contract = _prepare_degree_contract(
        person_degrees=person_degrees,
        rater_degrees=rater_degrees,
        person_coordinates=person_coordinates,
        rater_coordinates=rater_coordinates,
        gamma=gamma,
    )
    people = contract["people"]
    raters = contract["raters"]
    p_degree = contract["person_degrees"]
    r_degree = contract["rater_degrees"]
    n_samples = int(n_samples)
    max_dense_dp_cells = int(max_dense_dp_cells)
    max_total_rejections = int(max_total_rejections)
    if n_samples < 1 or n_samples > MAX_EXACT_INDEPENDENT_SAMPLES:
        raise ValueError(
            f"n_samples must lie in [1, {MAX_EXACT_INDEPENDENT_SAMPLES}]."
        )
    if max_dense_dp_cells < 1 or max_total_rejections < 0:
        raise ValueError(
            "Positive max_dense_dp_cells and nonnegative max_total_rejections are required."
        )
    axis_raters = raters[:-1]
    shape = tuple(r_degree[rater] + 1 for rater in axis_raters)
    cells_per_layer = int(np.prod(shape, dtype=np.int64))
    total_cells = int((len(people) + 1) * cells_per_layer)
    if total_cells > max_dense_dp_cells:
        raise ValueError(
            f"Dense DP requires {total_cells} cells, exceeding max_dense_dp_cells="
            f"{max_dense_dp_cells}."
        )

    started = time.perf_counter()
    layers = np.full((len(people) + 1, *shape), -np.inf, dtype=np.float64)
    layers[(len(people), *([0] * len(axis_raters)))] = 0.0
    rater_indices = tuple(range(len(raters)))
    for person_index in range(len(people) - 1, -1, -1):
        person = people[person_index]
        current = layers[person_index]
        following = layers[person_index + 1]
        for selected in combinations(rater_indices, p_degree[person]):
            increments = tuple(
                1 if rater_index in selected else 0
                for rater_index in range(len(axis_raters))
            )
            source_slices = tuple(
                slice(0, dimension - increment)
                for dimension, increment in zip(shape, increments, strict=True)
            )
            target_slices = tuple(
                slice(increment, dimension)
                for dimension, increment in zip(shape, increments, strict=True)
            )
            local_statistic = contract["person_z"][person] * sum(
                contract["rater_z"][raters[rater_index]]
                for rater_index in selected
            )
            candidate = following[source_slices] + contract["gamma"] * local_statistic
            np.logaddexp(current[target_slices], candidate, out=current[target_slices])
    initial_remaining = [r_degree[rater] for rater in raters]
    initial_index = tuple(initial_remaining[:-1])
    log_normalizer_unconditioned = float(layers[(0, *initial_index)])
    if not np.isfinite(log_normalizer_unconditioned):
        raise ValueError("No graph satisfies the dense fixed-margin degree contract.")
    partition_seconds = time.perf_counter() - started

    edge_count = float(contract["edge_count"])
    person_edge_mean = float(
        sum(
            contract["person_z"][person] * p_degree[person]
            for person in people
        ) / edge_count
    )
    rater_edge_mean = float(
        sum(
            contract["rater_z"][rater] * r_degree[rater]
            for rater in raters
        ) / edge_count
    )
    person_edge_sd = float(
        np.sqrt(
            sum(
                p_degree[person]
                * (contract["person_z"][person] - person_edge_mean) ** 2
                for person in people
            ) / edge_count
        )
    )
    rater_edge_sd = float(
        np.sqrt(
            sum(
                r_degree[rater]
                * (contract["rater_z"][rater] - rater_edge_mean) ** 2
                for rater in raters
            ) / edge_count
        )
    )

    rng = np.random.default_rng(int(seed))
    sample_rows: list[dict[str, object]] = []
    connected_rejections = 0
    maximum_statistic_residual = 0.0
    maximum_conditional_log_residual = 0.0
    final_edges: set[tuple[str, str]] = set()
    while len(sample_rows) < n_samples:
        remaining = list(initial_remaining)
        edges: set[tuple[str, str]] = set()
        statistic = 0.0
        for person_index, person in enumerate(people):
            available = [index for index, value in enumerate(remaining) if value > 0]
            candidates: list[tuple[tuple[int, ...], list[int], float, float]] = []
            for selected in combinations(available, p_degree[person]):
                next_remaining = list(remaining)
                for rater_index in selected:
                    next_remaining[rater_index] -= 1
                continuation = float(
                    layers[(person_index + 1, *tuple(next_remaining[:-1]))]
                )
                if not np.isfinite(continuation):
                    continue
                local_statistic = contract["person_z"][person] * sum(
                    contract["rater_z"][raters[rater_index]]
                    for rater_index in selected
                )
                candidates.append(
                    (
                        selected,
                        next_remaining,
                        float(local_statistic),
                        float(contract["gamma"] * local_statistic + continuation),
                    )
                )
            if not candidates:
                raise RuntimeError(
                    "Dense exact sampler reached an empty conditional allocation."
                )
            log_weights = np.asarray([candidate[3] for candidate in candidates], dtype=float)
            log_total = float(logsumexp(log_weights))
            current_log_partition = float(
                layers[(person_index, *tuple(remaining[:-1]))]
            )
            maximum_conditional_log_residual = max(
                maximum_conditional_log_residual,
                abs(log_total - current_log_partition),
            )
            probabilities = np.exp(log_weights - log_total)
            selected_index = int(rng.choice(len(candidates), p=probabilities))
            selected, remaining, local_statistic, _ = candidates[selected_index]
            edges.update((person, raters[rater_index]) for rater_index in selected)
            statistic += local_statistic
        if any(remaining):
            raise RuntimeError("Dense exact sampler did not exhaust Rater margins.")
        if require_connected and _direct_overlap_components(edges, raters) != 1:
            connected_rejections += 1
            if connected_rejections > max_total_rejections:
                raise RuntimeError(
                    "Dense connected rejection sampling exceeded max_total_rejections."
                )
            continue
        recomputed = _assignment_statistic(
            edges, contract["person_z"], contract["rater_z"]
        )
        residual = float(statistic - recomputed)
        maximum_statistic_residual = max(
            maximum_statistic_residual, abs(residual)
        )
        row: dict[str, object] = {
            "Sample": len(sample_rows) + 1,
            "Statistic": recomputed,
            "AssignmentCorrelation": (
                (recomputed / edge_count - person_edge_mean * rater_edge_mean)
                / (person_edge_sd * rater_edge_sd)
            ),
            "StatisticUpdateResidual": residual,
        }
        if retain_edge_signatures:
            row["EdgeSignature"] = _edge_signature(edges)
        sample_rows.append(row)
        final_edges = edges
    samples = pd.DataFrame(sample_rows)
    total_attempts = n_samples + connected_rejections
    diagnostics = pd.DataFrame(
        [
            {
                "SchemaVersion": MECHANISM_SCHEMA_VERSION,
                "Gamma": contract["gamma"],
                "Seed": int(seed),
                "Samples": len(samples),
                "ExactIndependentSamples": True,
                "DenseDPShape": "x".join(map(str, shape)),
                "DenseDPCells": total_cells,
                "DenseDPFiniteStates": int(np.isfinite(layers).sum()),
                "DenseDPMemoryBytes": int(layers.nbytes),
                "DenseDPCellCap": max_dense_dp_cells,
                "DenseDPPartitionSeconds": partition_seconds,
                "DPLogNormalizerUnconditioned": log_normalizer_unconditioned,
                "MaximumConditionalLogResidual": maximum_conditional_log_residual,
                "ConnectedRejections": connected_rejections,
                "TotalDrawAttempts": total_attempts,
                "ConnectedAcceptanceRate": n_samples / total_attempts,
                "UniqueSampledStates": (
                    int(samples["EdgeSignature"].nunique())
                    if "EdgeSignature" in samples
                    else np.nan
                ),
                "StatisticESSDiagnostic": _initial_positive_sequence_ess(
                    samples["Statistic"].to_numpy(dtype=float)
                ),
                "MaximumStatisticUpdateResidual": maximum_statistic_residual,
                "FinalOverlapComponents": _direct_overlap_components(final_edges, raters),
                "DegreeMarginsPreserved": True,
                "SamplingKernel": (
                    "exact vectorized dense fixed-margin DP with connected rejection sampling"
                ),
            }
        ]
    )
    available = bool(
        len(samples) == n_samples
        and maximum_statistic_residual <= 1e-10
        and maximum_conditional_log_residual <= 1e-10
        and (
            not require_connected
            or diagnostics.iloc[0]["FinalOverlapComponents"] == 1
        )
    )
    return {
        "available": available,
        "reason": (
            "Exact independent dense-DP assignment samples completed."
            if available
            else "Dense-DP assignment samples failed their invariant audit."
        ),
        "schema_version": MECHANISM_SCHEMA_VERSION,
        "samples": samples,
        "diagnostics": diagnostics,
        "final_edges": final_edges,
        "contract": contract,
        "claim_boundary": (
            "The dense DP is numerically equivalent to the recursive fixed-margin partition "
            "function within its registered tolerance. Its cell cap makes the small-Rater "
            "applicability boundary explicit; connected rejection remains exact and gamma known."
        ),
    }


def compare_chain_to_exact_oracle(
    chain: Mapping[str, object], oracle: Mapping[str, object]
) -> dict[str, object]:
    """Compare sampled state frequencies and moments to an exact enumeration."""

    samples = chain.get("samples", pd.DataFrame()) if isinstance(chain, Mapping) else pd.DataFrame()
    states = oracle.get("states", pd.DataFrame()) if isinstance(oracle, Mapping) else pd.DataFrame()
    if not isinstance(samples, pd.DataFrame) or samples.empty:
        raise ValueError("Chain samples are unavailable.")
    if not isinstance(states, pd.DataFrame) or states.empty:
        raise ValueError("Exact oracle states are unavailable.")
    empirical = samples["EdgeSignature"].value_counts(normalize=True).rename("EmpiricalProbability")
    comparison = states[["EdgeSignature", "Probability", "Statistic"]].merge(
        empirical,
        left_on="EdgeSignature",
        right_index=True,
        how="left",
        validate="one_to_one",
    )
    comparison["EmpiricalProbability"] = comparison["EmpiricalProbability"].fillna(0.0)
    comparison["AbsoluteProbabilityError"] = (
        comparison["EmpiricalProbability"] - comparison["Probability"]
    ).abs()
    total_variation = float(0.5 * comparison["AbsoluteProbabilityError"].sum())
    exact_mean = float(np.dot(comparison["Probability"], comparison["Statistic"]))
    empirical_mean = float(samples["Statistic"].mean())
    summary = pd.DataFrame(
        [
            {
                "ExactStates": len(comparison),
                "SampledStates": int((comparison["EmpiricalProbability"] > 0).sum()),
                "Samples": len(samples),
                "TotalVariationDistance": total_variation,
                "MaximumAbsoluteProbabilityError": float(
                    comparison["AbsoluteProbabilityError"].max()
                ),
                "ExactExpectedStatistic": exact_mean,
                "EmpiricalMeanStatistic": empirical_mean,
                "MeanStatisticError": empirical_mean - exact_mean,
                "ChainProbabilityOutsideOracle": float(
                    (~samples["EdgeSignature"].isin(set(states["EdgeSignature"]))).mean()
                ),
            }
        ]
    )
    return {
        "comparison": comparison,
        "summary": summary,
    }


def materialize_exchangeable_assignment_design(
    data: pd.DataFrame,
    target_edges: Iterable[tuple[str, str]],
    *,
    person_col: str = "Person",
    rater_col: str = "Rater",
    facet_cols: Iterable[str] = (),
    context_cols: Iterable[str] = (),
    weight_col: str | None = "Weight",
    rater_mapping_confirmed: bool = False,
) -> dict[str, object]:
    """Apply a sampled graph to equal-context blocks without copying outcomes."""

    # Local import avoids adding mechanism code to the public runner's import
    # path while reusing its deliberately strict exchangeability gate.
    from mfrm_app.assignment_sensitivity import (
        evaluate_fixed_density_perturbation_feasibility,
    )

    facet_cols = [str(column) for column in facet_cols]
    context_cols = [str(column) for column in context_cols]
    preflight = evaluate_fixed_density_perturbation_feasibility(
        data,
        person_col=person_col,
        rater_col=rater_col,
        context_cols=context_cols,
        weight_col=weight_col,
        rater_mapping_confirmed=rater_mapping_confirmed,
    )
    empty = {
        "design": pd.DataFrame(),
        "assignment_map": pd.DataFrame(),
        "invariants": pd.DataFrame(),
    }
    if not preflight.get("available"):
        return {**preflight, **empty}

    source_edges = {
        (str(person), str(rater)) for person, rater in preflight["edges"]
    }
    target = {(str(person), str(rater)) for person, rater in target_edges}
    people = list(preflight["persons"])
    raters = list(preflight["raters"])
    target_people = {person for person, _ in target}
    target_raters = {rater for _, rater in target}

    def degree(edges: set[tuple[str, str]], side: int) -> pd.Series:
        levels = people if side == 0 else raters
        counts = pd.Series(
            [edge[side] for edge in edges], dtype="object"
        ).value_counts()
        return counts.reindex(levels, fill_value=0).astype(int).sort_index()

    structural_rows: list[tuple[str, bool, str]] = [
        (
            "Target Person levels",
            target_people == set(people),
            f"target={len(target_people)}; source={len(people)}",
        ),
        (
            "Target Rater levels",
            target_raters == set(raters),
            f"target={len(target_raters)}; source={len(raters)}",
        ),
        (
            "Target edge count",
            len(target) == len(source_edges),
            f"target={len(target)}; source={len(source_edges)}",
        ),
        (
            "Per-Person Rater degree",
            degree(target, 0).equals(degree(source_edges, 0)),
            "exact equality",
        ),
        (
            "Per-Rater Person exposure",
            degree(target, 1).equals(degree(source_edges, 1)),
            "exact equality",
        ),
        (
            "Connected direct Rater overlap",
            _direct_overlap_components(target, raters) == 1,
            f"components={_direct_overlap_components(target, raters)}",
        ),
    ]
    structural = pd.DataFrame(
        structural_rows, columns=["Invariant", "Passed", "Evidence"]
    )
    if not structural["Passed"].all():
        return {
            **preflight,
            **empty,
            "available": False,
            "reason": "Target graph violates the fixed-density source contract.",
            "invariants": structural,
        }

    mapping: dict[tuple[str, str], str] = {}
    for person in people:
        source_raters = {rater for candidate, rater in source_edges if candidate == person}
        target_raters_for_person = {rater for candidate, rater in target if candidate == person}
        common = sorted(source_raters & target_raters_for_person)
        removed = sorted(source_raters - target_raters_for_person)
        added = sorted(target_raters_for_person - source_raters)
        if len(removed) != len(added):
            raise RuntimeError(f"Target mapping is inconsistent for Person {person}.")
        for rater in common:
            mapping[(person, rater)] = rater
        for source_rater, target_rater in zip(removed, added):
            mapping[(person, source_rater)] = target_rater

    work_cols = list(dict.fromkeys([person_col, rater_col, *facet_cols, *context_cols]))
    if weight_col and weight_col in data.columns:
        work_cols.append(weight_col)
    source = data[work_cols].copy().reset_index(drop=True)
    for column in dict.fromkeys([person_col, rater_col, *facet_cols, *context_cols]):
        source[column] = source[column].astype(str)
    source["_SourceRow"] = np.arange(len(source), dtype=int)
    counterfactual = source.copy()
    counterfactual[rater_col] = [
        mapping[(str(person), str(rater))]
        for person, rater in zip(source[person_col], source[rater_col])
    ]
    materialized_edges = {
        tuple(row)
        for row in counterfactual[[person_col, rater_col]].drop_duplicates().itertuples(
            index=False, name=None
        )
    }

    before_rater_rows = source.groupby(rater_col, observed=True).size().sort_index()
    after_rater_rows = counterfactual.groupby(rater_col, observed=True).size().sort_index()
    context_group = [rater_col, *context_cols]
    before_context = source.groupby(context_group, observed=True).size().sort_index()
    after_context = counterfactual.groupby(context_group, observed=True).size().sort_index()
    changed_edges = len(source_edges.symmetric_difference(materialized_edges)) // 2
    materialized_rows = [
        ("Materialized target graph", materialized_edges == target, "exact edge equality"),
        ("Response-row count", len(source) == len(counterfactual), f"{len(source)} -> {len(counterfactual)}"),
        ("Per-Rater response-row exposure", before_rater_rows.equals(after_rater_rows), "exact equality"),
        ("Per-Rater context-cell rows", before_context.equals(after_context), "exact equality"),
        ("Score-free design", "Score" not in counterfactual.columns, "observed outcome omitted"),
        ("Source-row identity retained", counterfactual["_SourceRow"].equals(source["_SourceRow"]), "exact equality"),
    ]
    invariants = pd.concat(
        [
            structural,
            pd.DataFrame(materialized_rows, columns=["Invariant", "Passed", "Evidence"]),
        ],
        ignore_index=True,
    )
    assignment_map = pd.DataFrame(
        [
            {
                "Person": person,
                "SourceRater": source_rater,
                "CounterfactualRater": target_rater,
                "Changed": source_rater != target_rater,
            }
            for (person, source_rater), target_rater in sorted(mapping.items())
        ]
    )
    available = bool(invariants["Passed"].all())
    return {
        **preflight,
        "available": available,
        "reason": (
            "The sampled graph was materialized as a score-free exchangeable-block design."
            if available else "The materialized design failed at least one invariant."
        ),
        "schema_version": MECHANISM_SCHEMA_VERSION,
        "design": counterfactual if available else pd.DataFrame(),
        "assignment_map": assignment_map,
        "invariants": invariants,
        "changed_edges": changed_edges,
        "claim_boundary": (
            "Only exchangeable equal-context Person-Rater blocks are materialized. Observed outcomes are not "
            "inputs or outputs; scores must be generated from an explicit response model afterward."
        ),
    }
