"""Outcome-blind assignment diagnostics and estimator estimand contracts.

The functions in this module deliberately keep three questions separate:

1. Is the observed person-by-rater assignment structure describable?
2. Can the assignment *mechanism* be identified from that structure? (No.)
3. Which likelihood/estimand does each estimator target?

No score, fitted ability, fitted severity, residual, or fit statistic is used by
the design audit.  This makes it suitable for display before users interpret
estimated measures and prevents a descriptive topology check from being
misreported as an MNAR/informative-assignment test.
"""

from __future__ import annotations

from itertools import combinations
from typing import Iterable, Mapping

import numpy as np
import pandas as pd


ASSIGNMENT_EVIDENCE_VERSION = "informative_assignment_confirmatory100_20260811"

_RATER_KEYWORDS = (
    "rater",
    "scorer",
    "judge",
    "evaluator",
    "assessor",
    "reader",
    "marker",
    "observer",
    "評定者",
    "採点者",
    "評価者",
    "判定者",
    "審査員",
    "観察者",
)


def _clean_names(values: Iterable[object]) -> list[str]:
    out: list[str] = []
    for value in values:
        name = str(value).strip()
        if name and name not in out and name != "Person":
            out.append(name)
    return out


def infer_rater_facet(
    facet_names: Iterable[object],
    *,
    configured_rater_facet: object | None = None,
) -> dict[str, object]:
    """Return a transparent rater-role mapping without silently certifying it.

    A configured role is preferred.  A keyword match is considered a plausible
    automatic mapping.  If neither is available, the first facet is returned as
    an *unconfirmed display candidate* so the UI can still explain what must be
    reviewed; it is never labelled as a confirmed rater mapping.
    """

    facets = _clean_names(facet_names)
    configured = str(configured_rater_facet or "").strip()
    if configured and configured in facets:
        return {
            "rater_facet": configured,
            "mapping_basis": "configured_role",
            "mapping_confirmed": True,
            "candidate_only": False,
        }

    for facet in facets:
        lowered = facet.casefold()
        if any(keyword.casefold() in lowered for keyword in _RATER_KEYWORDS):
            return {
                "rater_facet": facet,
                "mapping_basis": "name_keyword",
                "mapping_confirmed": True,
                "candidate_only": False,
            }

    if facets:
        return {
            "rater_facet": facets[0],
            "mapping_basis": "first_facet_unconfirmed",
            "mapping_confirmed": False,
            "candidate_only": True,
        }
    return {
        "rater_facet": None,
        "mapping_basis": "not_available",
        "mapping_confirmed": False,
        "candidate_only": False,
    }


def _gini(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0 or np.all(values == 0):
        return float("nan")
    values = np.sort(np.clip(values, 0.0, None))
    n = values.size
    return float((2.0 * np.dot(np.arange(1, n + 1), values) / (n * values.sum())) - (n + 1) / n)


def _component_labels(nodes: list[str], edges: list[tuple[str, str]]) -> dict[str, int]:
    parent = {node: node for node in nodes}

    def find(node: str) -> str:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left: str, right: str) -> None:
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for left, right in edges:
        union(left, right)
    roots = {root: idx + 1 for idx, root in enumerate(sorted({find(node) for node in nodes}))}
    return {node: roots[find(node)] for node in nodes}


def _empty_assignment_bundle(reason: str, mapping: Mapping[str, object] | None = None) -> dict[str, object]:
    mapping = dict(mapping or {})
    rater_facet = mapping.get("rater_facet")
    summary = pd.DataFrame(
        [
            {
                "AuditStatus": "Not auditable",
                "Reason": reason,
                "RaterFacet": rater_facet or "not mapped",
                "RaterMappingBasis": mapping.get("mapping_basis", "not_available"),
                "RaterMappingConfirmed": bool(mapping.get("mapping_confirmed", False)),
                "AssignmentMechanismStatus": "Not identified from observed assignments",
                "FixedDensityStressTestStatus": "Not run for this dataset",
                "OutcomeBlindAudit": True,
                "LikelihoodRowsOnly": True,
                "ClaimBoundary": (
                    "No MCAR, MAR, MNAR, informative-assignment, bias, or fairness conclusion "
                    "can be drawn from this audit."
                ),
            }
        ]
    )
    return {
        "available": False,
        "reason": reason,
        "mapping": mapping,
        "summary": summary,
        "rater_exposure": pd.DataFrame(),
        "rater_overlap": pd.DataFrame(),
        "sensitivity_plan": build_fixed_density_sensitivity_plan(False, reason=reason),
    }


def build_assignment_design_audit(
    data: pd.DataFrame,
    *,
    person_col: str = "Person",
    facet_names: Iterable[object] = (),
    rater_facet: str | None = None,
    rater_mapping_basis: str | None = None,
    rater_mapping_confirmed: bool | None = None,
) -> dict[str, object]:
    """Audit person-by-rater assignment topology without using outcomes.

    Counts are based on unique person-rater assignments, rather than raw response
    rows, so repeated tasks/criteria do not masquerade as additional assignment
    overlap.  Frequency weights are intentionally not used because they cannot
    reconstruct the topology of compressed individual assignments.
    """

    if not isinstance(data, pd.DataFrame) or data.empty:
        return _empty_assignment_bundle("Fitted likelihood-row data are unavailable.")

    configured = rater_facet
    mapping = infer_rater_facet(facet_names, configured_rater_facet=configured)
    if rater_mapping_basis is not None:
        mapping["mapping_basis"] = str(rater_mapping_basis)
    if rater_mapping_confirmed is not None:
        mapping["mapping_confirmed"] = bool(rater_mapping_confirmed)
        mapping["candidate_only"] = not bool(rater_mapping_confirmed)
    rater_facet = mapping.get("rater_facet")
    if not isinstance(rater_facet, str) or not rater_facet:
        return _empty_assignment_bundle("No non-person facet is available for rater-role review.", mapping)
    if person_col not in data.columns or rater_facet not in data.columns:
        return _empty_assignment_bundle(
            f"Required assignment columns are missing: person={person_col!r}, rater={rater_facet!r}.",
            mapping,
        )

    assignments = data[[person_col, rater_facet]].copy()
    assignments = assignments.dropna(subset=[person_col, rater_facet])
    assignments[person_col] = assignments[person_col].astype(str)
    assignments[rater_facet] = assignments[rater_facet].astype(str)
    assignments = assignments.loc[
        assignments[person_col].str.strip().ne("") & assignments[rater_facet].str.strip().ne("")
    ].drop_duplicates()
    if assignments.empty:
        return _empty_assignment_bundle("No nonblank person-rater assignments remain in likelihood rows.", mapping)

    people = sorted(assignments[person_col].unique().tolist())
    raters = sorted(assignments[rater_facet].unique().tolist())
    if len(raters) < 2:
        return _empty_assignment_bundle(
            "At least two rater levels are required for an assignment-overlap audit.",
            mapping,
        )

    person_loads = assignments.groupby(person_col, observed=True)[rater_facet].nunique().astype(int)
    rater_people = {
        rater: set(assignments.loc[assignments[rater_facet].eq(rater), person_col].tolist())
        for rater in raters
    }
    rater_counts = np.asarray([len(rater_people[rater]) for rater in raters], dtype=float)

    overlap_rows: list[dict[str, object]] = []
    overlap_edges: list[tuple[str, str]] = []
    for left, right in combinations(raters, 2):
        shared = len(rater_people[left] & rater_people[right])
        union_n = len(rater_people[left] | rater_people[right])
        jaccard = float(shared / union_n) if union_n else float("nan")
        if shared > 0:
            overlap_edges.append((left, right))
        overlap_rows.append(
            {
                "RaterA": left,
                "RaterB": right,
                "SharedPersons": int(shared),
                "UnionPersons": int(union_n),
                "JaccardOverlap": jaccard,
                "DirectOverlap": bool(shared > 0),
            }
        )
    overlap = pd.DataFrame(overlap_rows)
    component_by_rater = _component_labels(raters, overlap_edges)
    component_count = len(set(component_by_rater.values()))

    exposure_rows: list[dict[str, object]] = []
    for rater in raters:
        direct_neighbors = {
            right for left, right in overlap_edges if left == rater
        } | {
            left for left, right in overlap_edges if right == rater
        }
        shared_counts = [
            int(row["SharedPersons"])
            for row in overlap_rows
            if rater in {row["RaterA"], row["RaterB"]} and int(row["SharedPersons"]) > 0
        ]
        exposure_rows.append(
            {
                "Rater": rater,
                "PersonsAssessed": int(len(rater_people[rater])),
                "DirectOverlapRaters": int(len(direct_neighbors)),
                "TotalSharedPersonLinks": int(sum(shared_counts)),
                "OverlapComponent": int(component_by_rater[rater]),
                "IsolatedInDirectOverlapGraph": bool(len(direct_neighbors) == 0),
            }
        )
    exposure = pd.DataFrame(exposure_rows)

    possible_assignment_cells = len(people) * len(raters)
    possible_rater_pairs = len(raters) * (len(raters) - 1) // 2
    overlapping_pair_count = len(overlap_edges)
    positive_overlap = pd.to_numeric(
        overlap.loc[overlap["DirectOverlap"], "SharedPersons"], errors="coerce"
    )
    mean_exposure = float(rater_counts.mean()) if rater_counts.size else float("nan")
    exposure_cv = (
        float(rater_counts.std(ddof=0) / mean_exposure)
        if rater_counts.size and mean_exposure > 0
        else float("nan")
    )
    mapping_confirmed = bool(mapping.get("mapping_confirmed", False))
    topology_connected = component_count == 1
    if not mapping_confirmed:
        audit_status = "Review rater mapping"
        reason = "The first facet is only an unconfirmed rater-role candidate."
    elif not topology_connected:
        audit_status = "Review direct rater overlap"
        reason = f"The direct rater-overlap graph has {component_count} components."
    else:
        audit_status = "Descriptive audit ready"
        reason = "Rater role is mapped and the direct rater-overlap graph is connected."
    stress_eligible = bool(mapping_confirmed and topology_connected)
    stress_status = (
        "Eligible but not run for this dataset"
        if stress_eligible
        else "Not ready; confirm the rater role and connected overlap first"
    )

    summary = pd.DataFrame(
        [
            {
                "AuditStatus": audit_status,
                "Reason": reason,
                "RaterFacet": rater_facet,
                "RaterMappingBasis": mapping.get("mapping_basis"),
                "RaterMappingConfirmed": mapping_confirmed,
                "LikelihoodRows": int(len(data)),
                "UniquePersonRaterAssignments": int(len(assignments)),
                "Persons": int(len(people)),
                "Raters": int(len(raters)),
                "AssignmentCellDensity": float(len(assignments) / possible_assignment_cells),
                "MedianRatersPerPerson": float(person_loads.median()),
                "MinRatersPerPerson": int(person_loads.min()),
                "MaxRatersPerPerson": int(person_loads.max()),
                "CoRatedPersonShare": float((person_loads >= 2).mean()),
                "MinPersonsPerRater": int(rater_counts.min()),
                "MedianPersonsPerRater": float(np.median(rater_counts)),
                "MaxPersonsPerRater": int(rater_counts.max()),
                "RaterExposureCV": exposure_cv,
                "RaterExposureGini": _gini(rater_counts),
                "RaterOverlapComponents": int(component_count),
                "OverlappingRaterPairs": int(overlapping_pair_count),
                "PossibleRaterPairs": int(possible_rater_pairs),
                "OverlappingRaterPairShare": float(overlapping_pair_count / possible_rater_pairs),
                "MedianSharedPersonsPositivePairs": (
                    float(positive_overlap.median()) if not positive_overlap.empty else 0.0
                ),
                "AssignmentMechanismStatus": "Not identified from observed assignments",
                "FixedDensityStressTestStatus": stress_status,
                "OutcomeBlindAudit": True,
                "LikelihoodRowsOnly": True,
                "ClaimBoundary": (
                    "Topology and exposure are descriptive. They do not establish MCAR, MAR, MNAR, "
                    "informative assignment, estimator bias, or fairness."
                ),
            }
        ]
    )
    return {
        "available": True,
        "reason": reason,
        "mapping": mapping,
        "summary": summary,
        "rater_exposure": exposure,
        "rater_overlap": overlap,
        "sensitivity_plan": build_fixed_density_sensitivity_plan(
            stress_eligible,
            reason=(
                "Rater-role confirmation and a connected direct-overlap graph are required."
                if not stress_eligible else ""
            ),
        ),
    }


def build_fixed_density_sensitivity_plan(eligible: bool, *, reason: str = "") -> pd.DataFrame:
    """Return a prospective, non-executing sensitivity-analysis contract."""

    status = "Eligible; not run" if eligible else "Not ready"
    prerequisite = (
        "Observed person-rater topology is available."
        if eligible
        else (reason or "A defensible person-rater mapping and topology are required.")
    )
    rows = [
        {
            "Step": 1,
            "Analysis": "Freeze estimands and reporting rules",
            "Status": status,
            "RequiredControl": "Name JMLE, MML, and exact CMLE likelihood bases before fitting.",
            "Purpose": "Prevent estimator ranking and estimand drift.",
            "ClaimBoundary": "Never compare cross-basis likelihood, deviance, AIC, or BIC.",
        },
        {
            "Step": 2,
            "Analysis": "Construct fixed-density counterfactual assignments",
            "Status": status,
            "RequiredControl": (
                "Preserve included-row count, per-person load, rater exposure counts, and connectedness; "
                "vary only the strength/direction of assignment alignment."
            ),
            "Purpose": "Separate assignment dependence from simple sparsity or exposure imbalance.",
            "ClaimBoundary": "A counterfactual is a sensitivity scenario, not an estimate of the actual mechanism.",
        },
        {
            "Step": 3,
            "Analysis": "Refit within each estimator basis",
            "Status": status,
            "RequiredControl": "Hold model, anchors, scale, optimizer settings, and random seeds fixed within method.",
            "Purpose": "Measure local parameter movement and severity compression under assignment perturbation.",
            "ClaimBoundary": "Compare shifts within JMLE, within MML, or within exact CMLE only.",
        },
        {
            "Step": 4,
            "Analysis": "Report local sensitivity with denominators",
            "Status": status,
            "RequiredControl": "Report convergence, failures, extremes, constraints, and uncertainty for every scenario.",
            "Purpose": "Make sensitivity evidence auditable and failure-aware.",
            "ClaimBoundary": (
                "Do not transplant the external simulation effect size to the current dataset or label observed "
                "missingness as MNAR."
            ),
        },
    ]
    out = pd.DataFrame(rows)
    out["Prerequisite"] = prerequisite
    out["EvidenceVersion"] = ASSIGNMENT_EVIDENCE_VERSION
    return out


def build_informative_assignment_evidence_register() -> pd.DataFrame:
    """Return product-facing validation evidence without transporting effects.

    These rows are a compact mirror of the preregistered fresh-100 study.  They
    justify offering a sensitivity workflow; they are not correction factors or
    priors for a user's dataset.
    """

    common = {
        "EvidenceVersion": ASSIGNMENT_EVIDENCE_VERSION,
        "EvidenceSource": "validation/informative_assignment_confirmatory100_assessment_20260811.json",
        "Replicates": 100,
        "DesignContrast": "ability-severity aligned minus planned; equal density",
        "UseInCurrentDataset": "Rationale for sensitivity analysis only; do not transport as a correction.",
    }
    rows = [
        {
            **common,
            "Endpoint": "MML rater RMSE; free population SD",
            "Estimate": 0.09863414655507183,
            "CI95Lower": 0.0843219376607637,
            "CI95Upper": 0.11294635544937996,
            "MultiplicityRole": "Preregistered primary",
        },
        {
            **common,
            "Endpoint": "MML rater RMSE; fixed population SD",
            "Estimate": 0.06850620166504019,
            "CI95Lower": 0.05656143016767359,
            "CI95Upper": 0.08045097316240679,
            "MultiplicityRole": "Holm-confirmed secondary",
        },
        {
            **common,
            "Endpoint": "MML estimated population SD; free-SD fit",
            "Estimate": -0.11839387581077312,
            "CI95Lower": -0.1339436174586849,
            "CI95Upper": -0.10284413416286134,
            "MultiplicityRole": "Holm-confirmed secondary",
        },
        {
            **common,
            "Endpoint": "MML rater severity compression slope; free-SD fit",
            "Estimate": -0.48932000810870263,
            "CI95Lower": -0.5303870617183142,
            "CI95Upper": -0.4482529544990911,
            "MultiplicityRole": "Holm-confirmed secondary",
        },
        {
            **common,
            "Endpoint": "MML rater severity compression slope; fixed-SD fit",
            "Estimate": -0.38383509644697145,
            "CI95Lower": -0.41824896631225433,
            "CI95Upper": -0.34942122658168856,
            "MultiplicityRole": "Holm-confirmed secondary",
        },
    ]
    return pd.DataFrame(rows)


def build_estimand_contract(result_or_config: Mapping[str, object] | None) -> pd.DataFrame:
    """Describe JMLE/MML/exact-CMLE side by side without making a ranking.

    The selected run is marked, but likelihood values remain comparable only
    within the row's stated basis and like-for-like parameterization.
    """

    payload = dict(result_or_config or {})
    config_obj = payload.get("config", payload)
    config = dict(config_obj) if isinstance(config_obj, Mapping) else {}
    current_method = str(config.get("method", "JMLE") or "JMLE").upper()
    if current_method in {"EXACT CMLE", "CMLE", "CML"}:
        current_method = "EXACT_CMLE"
    free_sd = bool(config.get("estimate_population_sd", False))
    mml_scale = "estimated Gaussian population SD" if free_sd else "fixed Gaussian population SD"

    rows = [
        {
            "Method": "JMLE",
            "CurrentRun": current_method == "JMLE",
            "Availability": "Public Streamlit estimator",
            "EstimandClass": "Fixed-person joint calibration",
            "LikelihoodBasis": "Joint response likelihood with person and facet parameters",
            "PersonTreatment": "Each person is an incidental fixed parameter; boundary totals are flagged",
            "PopulationScale": "No random-person population distribution is estimated",
            "PersonOutput": "Finite JMLE person MLEs; boundary-score persons require separate handling",
            "FACETSRelation": (
                "External FACETS 4.5.0 comparator in qualified RSM/PCM JMLE scopes; "
                "the Python run does not call FACETS"
            ),
            "LikelihoodComparisonRule": "Compare only like-for-like JMLE fits on the same rows and parameterization",
            "EvidenceTier": "FACETS-calibrated JMLE scope plus Python known-truth evidence",
            "FitPrecisionBoundary": (
                "Do not recompute or adjudicate fit from FACETS two-decimal display values; use native precision "
                "or an explicitly tolerance-qualified export"
            ),
            "ClaimBoundary": "No cross-basis likelihood/AIC/BIC comparison or estimator ranking versus MML/CMLE",
        },
        {
            "Method": "MML",
            "CurrentRun": current_method == "MML",
            "Availability": "Public Streamlit estimator",
            "EstimandClass": "Gaussian-population marginal calibration",
            "LikelihoodBasis": "Person-marginal response likelihood",
            "PersonTreatment": "Person ability is integrated over the population model; person reports use posterior/EAP scores",
            "PopulationScale": mml_scale,
            "PersonOutput": "Posterior/EAP summaries, not JMLE person MLEs",
            "FACETSRelation": "Complementary estimator; FACETS JMLE is not a gold standard for the MML estimand",
            "LikelihoodComparisonRule": "Compare only like-for-like MML marginal fits with the same rows, scale, and integration basis",
            "EvidenceTier": "Python-native MML validation and fresh known-truth sensitivity evidence",
            "FitPrecisionBoundary": "Use unrounded Python probabilities and residual moments for fit calculations",
            "ClaimBoundary": "No FACETS parity claim and no cross-basis likelihood/AIC/BIC ranking versus JMLE/CMLE",
        },
        {
            "Method": "EXACT_CMLE",
            "CurrentRun": current_method == "EXACT_CMLE",
            "Availability": "Repository-only guarded workflow; not in the public estimator selector",
            "EstimandClass": "Person-total conditional facet calibration",
            "LikelihoodBasis": "Exact conditional likelihood given person sufficient totals",
            "PersonTreatment": "Person nuisance parameters are eliminated from facet calibration",
            "PopulationScale": "No Gaussian person-population distribution is required",
            "PersonOutput": "No native unconditional person estimate; WLE/person scoring is a downstream sidecar",
            "FACETSRelation": "Complementary conditional estimator; not a FACETS parity target",
            "LikelihoodComparisonRule": "Compare only like-for-like exact-CMLE conditional fits on the same conditional sample",
            "EvidenceTier": "Exact enumeration, finite-MLE existence, and independent cross-engine evidence",
            "FitPrecisionBoundary": "Retain full native precision; downstream WLE/fit evidence has its own qualification",
            "ClaimBoundary": "No cross-basis likelihood/AIC/BIC comparison or estimator ranking versus JMLE/MML",
        },
    ]
    return pd.DataFrame(rows)
