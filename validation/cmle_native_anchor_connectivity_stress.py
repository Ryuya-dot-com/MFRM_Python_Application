#!/usr/bin/env python3
"""Run the prospectively registered native-CMLE anchor-connectivity stress."""

from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import combinations
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle import fit_cmle, prepare_cmle_design  # noqa: E402
from mfrm_app.cmle_person_scoring import score_cmle_persons_wle  # noqa: E402
from mfrm_app.cmle_wle_fit import compute_cmle_wle_person_fit  # noqa: E402
from mfrm_app.fit_threshold_operating_characteristics import fit_rule_flags  # noqa: E402
from mfrm_app.operating_characteristics import frame_fingerprint  # noqa: E402


DEFAULT_PLAN = ROOT / "validation/cmle_native_anchor_connectivity_plan_20260810.json"
DEFAULT_OUTPUT = ROOT / "validation/cmle_native_anchor_connectivity_20260810"
RATERS = tuple(f"R{i:02d}" for i in range(1, 7))
CRITERIA = tuple(f"C{i:02d}" for i in range(1, 5))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, lineterminator="\n", float_format="%.17g")


def validate_plan(path: Path) -> dict[str, object]:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("study_id") != "cmle_native_anchor_connectivity_v1":
        raise ValueError("Unexpected anchor-connectivity plan.")
    mismatches = [
        relative
        for relative, digest in plan["input_identity"].items()
        if not (ROOT / relative).exists() or sha256_file(ROOT / relative) != digest
    ]
    if mismatches:
        raise ValueError(f"Anchor-connectivity input identity failed: {mismatches}")
    if plan["integrity_gates"]["performance_value_success_gate"] is not None:
        raise ValueError("Connectivity debugging stress cannot have a performance gate.")
    expected = (
        len(plan["topologies"])
        * int(plan["simulation"]["replicates"])
        * len(plan["anchor_scenarios"])
    )
    if expected != int(plan["fit_contract"]["expected_audits"]):
        raise ValueError("Expected audit count differs from the frozen factorial.")
    return plan


def _probabilities(theta: np.ndarray, difficulty: np.ndarray, steps: np.ndarray) -> np.ndarray:
    categories = np.arange(len(steps) + 1, dtype=float)
    cumulative = np.r_[0.0, np.cumsum(steps)]
    logits = (theta - difficulty)[:, None] * categories[None, :] - cumulative[None, :]
    logits -= logits.max(axis=1, keepdims=True)
    exp_logits = np.exp(logits)
    return exp_logits / exp_logits.sum(axis=1, keepdims=True)


def _draw_categories(probabilities: np.ndarray, uniforms: np.ndarray) -> np.ndarray:
    cumulative = np.cumsum(probabilities, axis=1)
    cumulative[:, -1] = 1.0
    return np.sum(uniforms[:, None] >= cumulative, axis=1).astype(int)


def generate_one(
    topology: dict[str, object], replicate: int, plan: dict[str, object]
) -> pd.DataFrame:
    constants = plan["simulation"]
    persons = int(constants["persons_per_dataset"])
    root_seed = int(constants["root_seed"])
    topology_order = int(topology["order"])
    theta_rng = np.random.default_rng(np.random.SeedSequence([root_seed, replicate, 0]))
    administration_rng = np.random.default_rng(
        np.random.SeedSequence([root_seed, replicate, topology_order, 1])
    )
    response_rng = np.random.default_rng(
        np.random.SeedSequence([root_seed, replicate, topology_order, 2])
    )
    theta = theta_rng.normal(size=persons)
    bridge_count = int(topology["bridge_persons"])
    person_order = administration_rng.permutation(persons)
    bridge_people = person_order[:bridge_count]
    nonbridge_people = person_order[bridge_count:]

    edge_assignments: dict[int, tuple[str, str]] = {}
    edges = [tuple(edge) for edge in topology["edges"]]
    if bridge_count:
        edge_vector = [edge for edge in edges for _ in range(int(topology["persons_per_edge"]))]
        if len(edge_vector) != bridge_count:
            raise ValueError("Frozen bridge count does not match edge support.")
        administration_rng.shuffle(edge_vector)
        edge_assignments = {
            int(person): tuple(edge)
            for person, edge in zip(bridge_people, edge_vector, strict=True)
        }
    base_vector = np.resize(np.asarray(RATERS, dtype=object), len(nonbridge_people))
    administration_rng.shuffle(base_vector)
    base_assignments = {
        int(person): str(rater)
        for person, rater in zip(nonbridge_people, base_vector, strict=True)
    }

    rows: list[dict[str, object]] = []
    width = len(str(persons))
    for person_index in range(persons):
        assigned = (
            edge_assignments[person_index]
            if person_index in edge_assignments
            else (base_assignments[person_index],)
        )
        for rater in assigned:
            for criterion in CRITERIA:
                rows.append(
                    {
                        "TopologyOrder": topology_order,
                        "Topology": str(topology["topology"]),
                        "TopologyFamily": str(topology["family"]),
                        "Replicate": int(replicate),
                        "RunId": f"{topology['topology']}::rep-{replicate:05d}",
                        "Person": f"P{person_index + 1:0{width}d}",
                        "PersonIndex": int(person_index),
                        "Rater": str(rater),
                        "Criterion": str(criterion),
                        "BridgePerson": bool(person_index in edge_assignments),
                        "TrueTheta": float(theta[person_index]),
                    }
                )
    frame = pd.DataFrame(rows)
    expected_rows = int(topology["response_rows_per_replicate"])
    if len(frame) != expected_rows:
        raise ValueError(f"Generated {len(frame)} rows; plan requires {expected_rows}.")
    rater_truth = constants["rater_truth"]
    criterion_truth = constants["criterion_truth"]
    step_truth = np.asarray(list(constants["step_truth"].values()), dtype=float)
    difficulty = np.asarray(
        [float(rater_truth[r]) + float(criterion_truth[c]) for r, c in zip(frame["Rater"], frame["Criterion"], strict=True)]
    )
    probabilities = _probabilities(
        frame["TrueTheta"].to_numpy(dtype=float), difficulty, step_truth
    )
    frame["ObservedCategory"] = _draw_categories(
        probabilities, response_rng.random(len(frame))
    )
    return frame


def generate_all(plan: dict[str, object]) -> pd.DataFrame:
    frames = [
        generate_one(topology, replicate, plan)
        for topology in plan["topologies"]
        for replicate in range(1, int(plan["simulation"]["replicates"]) + 1)
    ]
    combined = pd.concat(frames, ignore_index=True)
    if len(combined) != int(plan["simulation"]["unique_response_rows"]):
        raise ValueError("Generated response-row count differs from plan.")
    if combined["RunId"].nunique() != int(plan["simulation"]["unique_datasets"]):
        raise ValueError("Generated dataset count differs from plan.")
    return combined


def anchor_levels(plan: dict[str, object], replicate: int, scenario: str) -> list[str]:
    if scenario == "a0_unanchored":
        return []
    if scenario == "a6_all_correct":
        return list(RATERS)
    frozen = plan["frozen_random_anchor_levels"][str(int(replicate))]
    levels = list(frozen[scenario])
    expected = next(
        int(row["anchor_count"])
        for row in plan["anchor_scenarios"]
        if row["scenario"] == scenario
    )
    if len(levels) != expected or len(set(levels)) != expected:
        raise ValueError("Frozen random anchor assignment is malformed.")
    return levels


def hard_anchors(plan: dict[str, object], levels: list[str]) -> pd.DataFrame | None:
    if not levels:
        return None
    truth = plan["simulation"]["rater_truth"]
    return pd.DataFrame(
        [
            {
                "ParameterType": "Facet",
                "Facet": "Rater",
                "Level": level,
                "Value": float(truth[level]),
            }
            for level in levels
        ]
    )


class _UnionFind:
    def __init__(self, values: tuple[str, ...]):
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: str, right: str) -> None:
        a, b = self.find(left), self.find(right)
        if a != b:
            self.parent[b] = a


def realized_graph(run_frame: pd.DataFrame) -> dict[str, object]:
    totals = (
        run_frame.groupby("Person", sort=False)
        .agg(Total=("ObservedCategory", "sum"), Rows=("ObservedCategory", "size"))
    )
    informative = set(
        totals.index[(totals["Total"] > 0) & (totals["Total"] < 3 * totals["Rows"])]
    )
    support: dict[tuple[str, str], int] = defaultdict(int)
    union = _UnionFind(RATERS)
    for person, group in run_frame.groupby("Person", sort=False):
        if person not in informative:
            continue
        levels = sorted(set(group["Rater"].astype(str)))
        for left, right in combinations(levels, 2):
            edge = (left, right)
            support[edge] += 1
            union.union(left, right)
    grouped: dict[str, list[str]] = defaultdict(list)
    for rater in RATERS:
        grouped[union.find(rater)].append(rater)
    components = sorted((sorted(values) for values in grouped.values()), key=lambda x: x[0])
    return {
        "informative_persons": len(informative),
        "extreme_persons": int(run_frame["Person"].nunique() - len(informative)),
        "edge_support": dict(sorted(support.items())),
        "components": components,
    }


def graph_prediction(components: list[list[str]], anchors: list[str]) -> bool:
    anchored = set(anchors)
    if not anchored:
        return len(components) == 1
    for component in components:
        free = set(component) - anchored
        if free and not (set(component) & anchored):
            return False
    return True


def truth_for(plan: dict[str, object], row: pd.Series) -> float:
    if row["ParameterType"] == "Facet":
        return float(plan["simulation"][f"{str(row['Facet']).lower()}_truth"][str(row["Level"])])
    return float(plan["simulation"]["step_truth"][str(int(row["Step"]))])


def analyze_one(
    run_frame: pd.DataFrame,
    topology: dict[str, object],
    scenario: dict[str, object],
    plan: dict[str, object],
    graph: dict[str, object],
) -> tuple[dict[str, object], dict[str, object] | None, pd.DataFrame, pd.DataFrame]:
    replicate = int(run_frame["Replicate"].iloc[0])
    scenario_name = str(scenario["scenario"])
    levels = anchor_levels(plan, replicate, scenario_name)
    anchors = hard_anchors(plan, levels)
    analysis = run_frame[["Person", "Rater", "Criterion", "ObservedCategory"]].copy()
    started = time.perf_counter()
    design = prepare_cmle_design(
        analysis,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="ObservedCategory",
        rating_min=0,
        rating_max=3,
        model="RSM",
        hard_anchors=anchors,
    )
    audit = design.audit
    issues = audit["issues_table"]
    predicted = graph_prediction(graph["components"], levels)
    audit_row = {
        "TopologyOrder": int(topology["order"]),
        "Topology": str(topology["topology"]),
        "TopologyFamily": str(topology["family"]),
        "Replicate": replicate,
        "RunId": str(run_frame["RunId"].iloc[0]),
        "ScenarioOrder": int(scenario["order"]),
        "Scenario": scenario_name,
        "AnchorCount": len(levels),
        "AnchorLevels": ";".join(levels),
        "ResponseRows": len(analysis),
        "InputFingerprint": frame_fingerprint(analysis),
        "InformativePersonsGraph": int(graph["informative_persons"]),
        "ExtremePersonsGraph": int(graph["extreme_persons"]),
        "RealizedEdges": len(graph["edge_support"]),
        "RealizedComponents": len(graph["components"]),
        "ComponentMembership": "|".join(",".join(x) for x in graph["components"]),
        "MinimumPositiveEdgeSupport": min(graph["edge_support"].values()) if graph["edge_support"] else 0,
        "GraphPredictedRaterIdentified": predicted,
        "ExactEligible": bool(audit["eligible"]),
        "GraphExactAgreement": predicted == bool(audit["eligible"]),
        "KParams": int(design.n_parameters),
        "ConditionalRank": int(audit["conditional_rank"]),
        "ConditionalNullity": int(audit["conditional_nullity"]),
        "ConditionalConditionNumber": float(audit["conditional_condition_number"]),
        "RankTolerance": float(audit["rank_tolerance"]),
        "RankWorkProxy": int(audit["rank_audit_work_proxy"]),
        "RankPeakBytesProxy": int(audit["rank_audit_peak_bytes_proxy"]),
        "IssueCodes": ";".join(issues["Code"].astype(str)),
        "FitAttempted": bool(audit["eligible"]),
        "AuditElapsedSeconds": float(time.perf_counter() - started),
    }
    if not bool(audit["eligible"]):
        return audit_row, None, pd.DataFrame(), pd.DataFrame()

    fit_started = time.perf_counter()
    try:
        fit = fit_cmle(
            analysis,
            person_col="Person",
            facet_cols=["Rater", "Criterion"],
            score_col="ObservedCategory",
            rating_min=0,
            rating_max=3,
            model="RSM",
            hard_anchors=anchors,
            gtol=float(plan["fit_contract"]["gtol"]),
            maxiter=int(plan["fit_contract"]["maxiter"]),
            newton_polish_maxiter=int(plan["fit_contract"]["newton_polish_maxiter"]),
        )
        summary = fit["summary"].iloc[0]
        fitted_anchors = fit["facets"]["others"].loc[fit["facets"]["others"]["Anchored"]]
        expected_values = {
            level: float(plan["simulation"]["rater_truth"][level]) for level in levels
        }
        anchor_exact = len(fitted_anchors) == len(levels) and all(
            float(row.Estimate) == expected_values[str(row.Level)] and float(row.SE) == 0.0
            for row in fitted_anchors.itertuples(index=False)
        )
        fit_row = {
            **{key: audit_row[key] for key in (
                "TopologyOrder", "Topology", "TopologyFamily", "Replicate", "RunId",
                "ScenarioOrder", "Scenario", "AnchorCount", "AnchorLevels",
                "MinimumPositiveEdgeSupport", "RealizedComponents", "InputFingerprint"
            )},
            "Attempted": True,
            "Returned": True,
            "FailureReason": "",
            "AnchorExact": anchor_exact,
            "Converged": bool(summary["Converged"]),
            "InferenceReady": bool(summary["InferenceReady"]),
            "KParams": int(summary["KParams"]),
            "ConditionalLogLik": float(summary["ConditionalLogLik"]),
            "InformationRank": int(summary["InformationRank"]),
            "InformationNullity": int(summary["InformationNullity"]),
            "InformationConditionNumber": float(summary["InformationConditionNumber"]),
            "GradientSupNorm": float(summary["GradientSupNorm"]),
            "OptimizerFiniteEvaluations": int(
                summary.get("OptimizerFiniteEvaluations", 0)
            ),
            "OptimizerInvalidEvaluations": int(
                summary.get("OptimizerInvalidEvaluations", 0)
            ),
            "OptimizerBestFiniteFallbackUsed": bool(
                summary.get("OptimizerBestFiniteFallbackUsed", False)
            ),
            "OptimizerLastInvalidReason": str(
                summary.get("OptimizerLastInvalidReason", "")
            ),
            "FiniteMLEGateEnabled": bool(
                summary.get("FiniteMLEGateEnabled", False)
            ),
            "FiniteMLEStatus": str(
                summary.get("FiniteMLEStatus", "not_evaluated")
            ),
            "FiniteMLEReason": str(summary.get("FiniteMLEReason", "")),
            "FiniteMLEBoundaryDetected": summary.get(
                "FiniteMLEBoundaryDetected", pd.NA
            ),
            "FiniteMLEExistenceQualified": bool(
                summary.get("FiniteMLEExistenceQualified", False)
            ),
            "FiniteMLETheoreticalConfigurations": summary.get(
                "FiniteMLETheoreticalConfigurations", np.nan
            ),
            "FiniteMLEGeneratedConstraintsMax": summary.get(
                "FiniteMLEGeneratedConstraintsMax", np.nan
            ),
            "FiniteMLEOracleCallsTotal": summary.get(
                "FiniteMLEOracleCallsTotal", np.nan
            ),
            "ReadinessReasons": str(summary["ReadinessReasons"]),
            "ElapsedSeconds": float(time.perf_counter() - fit_started),
        }
        structural = pd.concat([fit["facets"]["others"], fit["steps"]], ignore_index=True)
        structural.insert(0, "Scenario", scenario_name)
        structural.insert(0, "RunId", audit_row["RunId"])
        structural.insert(0, "Replicate", replicate)
        structural.insert(0, "Topology", audit_row["Topology"])
        structural["Truth"] = structural.apply(lambda row: truth_for(plan, row), axis=1)
        structural["Error"] = structural["Estimate"] - structural["Truth"]
        structural["AbsoluteError"] = structural["Error"].abs()
        structural["SquaredError"] = structural["Error"] ** 2
        persons = pd.DataFrame()
        if bool(summary["InferenceReady"]):
            wle = score_cmle_persons_wle(fit)[
                ["Person", "Estimate", "StandardError", "ExtremeScorePattern", "Status"]
            ].rename(columns={
                "Estimate": "WLEEstimate",
                "StandardError": "WLEStandardError",
                "ExtremeScorePattern": "WLEExactExtreme",
                "Status": "WLEStatus",
            })
            person_fit = compute_cmle_wle_person_fit(fit)["persons"]
            truth = run_frame.groupby("Person", sort=False)["TrueTheta"].first().reset_index()
            persons = (
                wle.merge(person_fit[["Person", "Infit", "Outfit", "PersonFitReady"]], on="Person", validate="one_to_one")
                .merge(truth, on="Person", validate="one_to_one")
            )
            persons.insert(0, "Scenario", scenario_name)
            persons.insert(0, "RunId", audit_row["RunId"])
            persons.insert(0, "Replicate", replicate)
            persons.insert(0, "Topology", audit_row["Topology"])
            persons["WLEError"] = persons["WLEEstimate"] - persons["TrueTheta"]
            raw = fit_rule_flags(persons["Infit"], persons["Outfit"])
            persons["FitEligibleRaw"] = raw["eligible"]
            persons["EitherUpperRaw"] = raw["either_upper"]
            for decimals in plan["raw_fit_contract"]["display_counterfactual_decimals"]:
                rounded = fit_rule_flags(
                    np.round(persons["Infit"].to_numpy(dtype=float), int(decimals)),
                    np.round(persons["Outfit"].to_numpy(dtype=float), int(decimals)),
                )["either_upper"]
                persons[f"EitherUpperRounded{int(decimals)}"] = rounded
                persons[f"RawRounded{int(decimals)}Disagreement"] = rounded != raw["either_upper"]
            persons["EitherDistanceTo1.5"] = np.minimum(
                np.abs(persons["Infit"] - 1.5), np.abs(persons["Outfit"] - 1.5)
            )
        return audit_row, fit_row, structural, persons
    except Exception as exc:
        fit_row = {
            **{key: audit_row[key] for key in (
                "TopologyOrder", "Topology", "TopologyFamily", "Replicate", "RunId",
                "ScenarioOrder", "Scenario", "AnchorCount", "AnchorLevels",
                "MinimumPositiveEdgeSupport", "RealizedComponents", "InputFingerprint"
            )},
            "Attempted": True,
            "Returned": False,
            "FailureReason": f"{type(exc).__name__}: {str(exc)[:500]}",
            "AnchorExact": False,
            "Converged": False,
            "InferenceReady": False,
            "KParams": np.nan,
            "ConditionalLogLik": np.nan,
            "InformationRank": np.nan,
            "InformationNullity": np.nan,
            "InformationConditionNumber": np.nan,
            "GradientSupNorm": np.nan,
            "OptimizerFiniteEvaluations": 0,
            "OptimizerInvalidEvaluations": 0,
            "OptimizerBestFiniteFallbackUsed": False,
            "OptimizerLastInvalidReason": "",
            "FiniteMLEGateEnabled": True,
            "FiniteMLEStatus": "fit_exception",
            "FiniteMLEReason": "fit_exception_before_status_return",
            "FiniteMLEBoundaryDetected": pd.NA,
            "FiniteMLEExistenceQualified": False,
            "FiniteMLETheoreticalConfigurations": np.nan,
            "FiniteMLEGeneratedConstraintsMax": np.nan,
            "FiniteMLEOracleCallsTotal": np.nan,
            "ReadinessReasons": "fit_exception",
            "ElapsedSeconds": float(time.perf_counter() - fit_started),
        }
        return audit_row, fit_row, pd.DataFrame(), pd.DataFrame()


def aggregate_audits(audits: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, group in audits.groupby(
        ["TopologyOrder", "Topology", "ScenarioOrder", "Scenario"], sort=True
    ):
        finite_condition = pd.to_numeric(group["ConditionalConditionNumber"], errors="coerce")
        finite_condition = finite_condition[np.isfinite(finite_condition)]
        rows.append({
            "TopologyOrder": keys[0],
            "Topology": keys[1],
            "ScenarioOrder": keys[2],
            "Scenario": keys[3],
            "Replicates": len(group),
            "GraphPredictedIdentified": int(group["GraphPredictedRaterIdentified"].sum()),
            "ExactEligible": int(group["ExactEligible"].sum()),
            "GraphExactAgreement": int(group["GraphExactAgreement"].sum()),
            "MeanNullity": float(group["ConditionalNullity"].mean()),
            "MeanConditionNumberFinite": float(finite_condition.mean()) if len(finite_condition) else np.nan,
            "MaxConditionNumberFinite": float(finite_condition.max()) if len(finite_condition) else np.nan,
            "MeanMinimumPositiveEdgeSupport": float(group["MinimumPositiveEdgeSupport"].mean()),
        })
    return pd.DataFrame(rows)


def aggregate_fits(fits: pd.DataFrame, persons: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if fits.empty:
        return pd.DataFrame()
    for keys, group in fits.groupby(
        ["TopologyOrder", "Topology", "ScenarioOrder", "Scenario"], sort=True
    ):
        subset = persons.loc[
            persons["Topology"].eq(keys[1]) & persons["Scenario"].eq(keys[3])
        ] if not persons.empty else pd.DataFrame()
        rows.append({
            "TopologyOrder": keys[0],
            "Topology": keys[1],
            "ScenarioOrder": keys[2],
            "Scenario": keys[3],
            "AttemptedFits": len(group),
            "ReturnedFits": int(group["Returned"].sum()),
            "InferenceReadyFits": int(group["InferenceReady"].sum()),
            "MeanFittedConditionNumber": float(pd.to_numeric(group["InformationConditionNumber"], errors="coerce").mean()),
            "PersonRows": len(subset),
            "WLERMSE": float(np.sqrt(np.mean(np.square(subset["WLEError"])))) if len(subset) else np.nan,
            "EitherUpperRawRate": float(subset["EitherUpperRaw"].mean()) if len(subset) else np.nan,
            "Rounded3Disagreements": int(subset["RawRounded3Disagreement"].sum()) if len(subset) else 0,
            "Rounded6Disagreements": int(subset["RawRounded6Disagreement"].sum()) if len(subset) else 0,
        })
    return pd.DataFrame(rows)


def render_figures(audits: pd.DataFrame, fits: pd.DataFrame, output: Path) -> None:
    topology_order = audits[["TopologyOrder", "Topology"]].drop_duplicates().sort_values("TopologyOrder")
    scenario_order = audits[["ScenarioOrder", "Scenario"]].drop_duplicates().sort_values("ScenarioOrder")
    matrix = audits.pivot_table(
        index="Topology", columns="Scenario", values="ExactEligible", aggfunc="mean"
    ).reindex(index=topology_order["Topology"], columns=scenario_order["Scenario"])
    fig, ax = plt.subplots(figsize=(10.5, 5.7))
    image = ax.imshow(matrix.to_numpy(dtype=float), vmin=0, vmax=1, cmap="RdYlGn", aspect="auto")
    ax.set_xticks(range(len(matrix.columns)), [x.replace("_", "\n") for x in matrix.columns], fontsize=8)
    ax.set_yticks(range(len(matrix.index)), [x.replace("_", " ") for x in matrix.index], fontsize=9)
    for i in range(len(matrix.index)):
        for j in range(len(matrix.columns)):
            ax.text(j, i, f"{matrix.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8)
    ax.set_title("Exact CMLE prefit eligibility rate (10 replicates)")
    ax.set_xlabel("Prospectively frozen correct-anchor scenario")
    ax.set_ylabel("Rater bridge topology")
    fig.colorbar(image, ax=ax, label="Eligibility rate")
    fig.tight_layout()
    fig.savefig(output / "connectivity_eligibility_heatmap.png", dpi=180)
    plt.close(fig)

    ready = fits.loc[fits["InferenceReady"] & np.isfinite(fits["InformationConditionNumber"])].copy()
    fig, ax = plt.subplots(figsize=(9.5, 5.7))
    for scenario, group in ready.groupby("Scenario", sort=False):
        ax.scatter(
            group["MinimumPositiveEdgeSupport"],
            group["InformationConditionNumber"],
            alpha=0.68,
            s=28,
            label=scenario,
        )
    ax.set_yscale("log")
    ax.set_xlabel("Minimum realized informative-Person support over positive Rater edges")
    ax.set_ylabel("Fitted information condition number (log scale)")
    ax.set_title("Connectivity support, anchor content, and numerical conditioning")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(output / "connectivity_condition_number.png", dpi=180)
    plt.close(fig)


def _markdown_table(frame: pd.DataFrame) -> str:
    """Render a small evidence table without pandas' optional tabulate dependency."""

    def display(value: object) -> str:
        if pd.isna(value):
            return ""
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.4g}"
        return str(value).replace("|", "\\|")

    header = "| " + " | ".join(str(column) for column in frame.columns) + " |"
    separator = "| " + " | ".join("---" for _ in frame.columns) + " |"
    rows = [
        "| " + " | ".join(display(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *rows])


def report_text(decision: dict[str, object], audit_summary: pd.DataFrame) -> str:
    table = _markdown_table(audit_summary[[
        "Topology", "Scenario", "Replicates", "GraphPredictedIdentified",
        "ExactEligible", "GraphExactAgreement", "MeanNullity",
        "MeanConditionNumberFinite",
    ]])
    return f"""# Native CMLE anchor-connectivity stress

> Prospectively registered ten-replicate debugging evidence. It does not establish a universal anchor share, performance probability, cross-engine equivalence, or public UI readiness.

## Contract

- Contract passed: `{decision['contract_passed']}`
- Audits: `{decision['audits']}/350`
- Prefit eligible / attempted fits: `{decision['prefit_eligible']}/{decision['fit_attempts']}`
- Returned / inference-ready fits: `{decision['returned_fits']}/{decision['inference_ready_fits']}`
- Graph prediction versus exact eligibility agreement: `{decision['graph_exact_agreements']}/{decision['audits']}`
- Raw versus 3-decimal disagreements: `{decision['rounded3_disagreements']}`; 6-decimal: `{decision['rounded6_disagreements']}`

## Rank and conditioning summary

{table}

## Boundary

The graph rule is a transparent Rater-block prediction, not a substitute for exact conditional-information rank. Nominal bridge edges are counted only when an observed non-extreme Person supplies the within-Person contrast. Correct anchors can fix component origins but cannot create empirical contrasts within a component. All fit decisions use raw values; rounded columns are counterfactual audits only. Public UI and estimator-performance claims remain withheld.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    plan_path = args.plan.resolve()
    output = args.output.resolve()
    plan = validate_plan(plan_path)
    if output.exists() and any(output.iterdir()):
        if not args.overwrite:
            raise FileExistsError(f"Output is non-empty: {output}; pass --overwrite explicitly.")
        shutil.rmtree(output)
    output.mkdir(parents=True, exist_ok=True)
    shutil.copy2(plan_path, output / plan_path.name)

    generated = generate_all(plan)
    response_path = output / "connectivity_responses.csv"
    write_csv(generated, response_path)
    reloaded = pd.read_csv(response_path, float_precision="round_trip")
    source_fingerprint = frame_fingerprint(generated)
    reloaded_fingerprint = frame_fingerprint(reloaded)
    source_roundtrip_exact = source_fingerprint == reloaded_fingerprint
    if not source_roundtrip_exact:
        raise ValueError("Generated source failed the exact round-trip fingerprint.")

    audit_rows: list[dict[str, object]] = []
    fit_rows: list[dict[str, object]] = []
    structural_frames: list[pd.DataFrame] = []
    person_frames: list[pd.DataFrame] = []
    graph_edge_rows: list[dict[str, object]] = []
    graph_component_rows: list[dict[str, object]] = []
    topologies = {str(row["topology"]): row for row in plan["topologies"]}
    for run_id, run_frame in reloaded.groupby("RunId", sort=False):
        topology = topologies[str(run_frame["Topology"].iloc[0])]
        graph = realized_graph(run_frame)
        for (left, right), support in graph["edge_support"].items():
            graph_edge_rows.append({
                "Topology": str(topology["topology"]),
                "Replicate": int(run_frame["Replicate"].iloc[0]),
                "RunId": str(run_id),
                "RaterLeft": left,
                "RaterRight": right,
                "InformativePersonSupport": int(support),
            })
        for component_index, component in enumerate(graph["components"], start=1):
            for rater in component:
                graph_component_rows.append({
                    "Topology": str(topology["topology"]),
                    "Replicate": int(run_frame["Replicate"].iloc[0]),
                    "RunId": str(run_id),
                    "Component": component_index,
                    "Rater": rater,
                    "InformativePersons": int(graph["informative_persons"]),
                })
        for scenario in plan["anchor_scenarios"]:
            audit, fit, structural, persons = analyze_one(
                run_frame, topology, scenario, plan, graph
            )
            audit_rows.append(audit)
            if fit is not None:
                fit_rows.append(fit)
            if not structural.empty:
                structural_frames.append(structural)
            if not persons.empty:
                person_frames.append(persons)

    audits = pd.DataFrame(audit_rows).sort_values(
        ["TopologyOrder", "Replicate", "ScenarioOrder"]
    ).reset_index(drop=True)
    fits = pd.DataFrame(fit_rows).sort_values(
        ["TopologyOrder", "Replicate", "ScenarioOrder"]
    ).reset_index(drop=True)
    structural = pd.concat(structural_frames, ignore_index=True) if structural_frames else pd.DataFrame()
    persons = pd.concat(person_frames, ignore_index=True) if person_frames else pd.DataFrame()
    graph_edges = pd.DataFrame(graph_edge_rows)
    graph_components = pd.DataFrame(graph_component_rows)
    audit_summary = aggregate_audits(audits)
    fit_summary = aggregate_fits(fits, persons)

    outputs = {
        "connectivity_audit_ledger.csv": audits,
        "connectivity_audit_summary.csv": audit_summary,
        "connectivity_fit_ledger.csv": fits,
        "connectivity_fit_summary.csv": fit_summary,
        "connectivity_graph_edges.csv": graph_edges,
        "connectivity_graph_components.csv": graph_components,
        "connectivity_structural_results.csv": structural,
        "connectivity_person_results.csv": persons,
    }
    for name, frame in outputs.items():
        write_csv(frame, output / name)
    render_figures(audits, fits, output)

    raw_recomputed = True
    if not persons.empty:
        recomputed = fit_rule_flags(persons["Infit"], persons["Outfit"])["either_upper"]
        raw_recomputed = bool(np.array_equal(recomputed, persons["EitherUpperRaw"].to_numpy(dtype=bool)))
    expected_audits = int(plan["fit_contract"]["expected_audits"])
    prefit_eligible = int(audits["ExactEligible"].sum())
    anchor_exact = bool(fits.loc[fits["Returned"], "AnchorExact"].all()) if len(fits) else False
    same_byte = bool(audits.groupby("RunId")["InputFingerprint"].nunique().eq(1).all())
    all_eligible_attempted = prefit_eligible == len(fits) and bool(audits.loc[audits["ExactEligible"], "FitAttempted"].all())
    contract_passed = bool(
        len(audits) == expected_audits
        and source_roundtrip_exact
        and same_byte
        and all_eligible_attempted
        and anchor_exact
        and raw_recomputed
    )
    decision = {
        "schema_version": "mfrm-cmle-native-anchor-connectivity-decision-v1",
        "study_id": plan["study_id"],
        "contract_passed": contract_passed,
        "audits": len(audits),
        "prefit_eligible": prefit_eligible,
        "fit_attempts": len(fits),
        "returned_fits": int(fits["Returned"].sum()),
        "inference_ready_fits": int(fits["InferenceReady"].sum()),
        "graph_exact_agreements": int(audits["GraphExactAgreement"].sum()),
        "graph_exact_disagreements": int((~audits["GraphExactAgreement"]).sum()),
        "same_byte_pairing_passed": same_byte,
        "source_roundtrip_exact": source_roundtrip_exact,
        "source_fingerprint": source_fingerprint,
        "reloaded_fingerprint": reloaded_fingerprint,
        "retained_response_sha256": sha256_file(response_path),
        "plan_sha256": sha256_file(plan_path),
        "anchor_exact_passed": anchor_exact,
        "raw_fit_recomputation_passed": raw_recomputed,
        "rounded3_disagreements": int(persons["RawRounded3Disagreement"].sum()) if len(persons) else 0,
        "rounded6_disagreements": int(persons["RawRounded6Disagreement"].sum()) if len(persons) else 0,
        "within_0.0005_of_1.5": int((persons["EitherDistanceTo1.5"] <= 0.0005).sum()) if len(persons) else 0,
        "within_0.001_of_1.5": int((persons["EitherDistanceTo1.5"] <= 0.001).sum()) if len(persons) else 0,
        "performance_value_success_gate": None,
        "universal_anchor_share_ready": False,
        "cross_engine_claims_ready": False,
        "public_ui_ready": False,
        "output_sha256": {
            name: sha256_file(output / name) for name in outputs
        },
    }
    decision_path = output / "connectivity_decision.json"
    decision_path.write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "CMLE_NATIVE_ANCHOR_CONNECTIVITY_STRESS.md").write_text(
        report_text(decision, audit_summary), encoding="utf-8"
    )
    print(json.dumps(decision, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
