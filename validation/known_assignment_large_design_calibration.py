"""Outcome-blind large-design mixing and assignment-dose calibration."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import scipy

from mfrm_app.assignment_mechanism import (
    MECHANISM_SCHEMA_VERSION,
    materialize_exchangeable_assignment_design,
    sample_degree_conditioned_assignment_chain,
)


PLAN_PATH = ROOT / "validation" / "known_assignment_large_design_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "known_assignment_large_design_20260811"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _planned_edges(people: list[str], raters: list[str]) -> set[tuple[str, str]]:
    edges: set[tuple[str, str]] = set()
    for index, person in enumerate(people):
        edges.add((person, raters[index % len(raters)]))
        edges.add((person, raters[(index + 1) % len(raters)]))
    return edges


def _exchangeable_rows(edges: set[tuple[str, str]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for person, rater in sorted(edges):
        for task in ("T01", "T02", "T03"):
            for criterion in ("C01", "C02"):
                rows.append(
                    {
                        "Person": person,
                        "Rater": rater,
                        "Task": task,
                        "Criterion": criterion,
                        "Weight": 1.0,
                    }
                )
    return pd.DataFrame(rows)


def _split_rhat(chains: list[np.ndarray]) -> float:
    if len(chains) < 2:
        return float("nan")
    length = min(len(chain) for chain in chains)
    half = length // 2
    if half < 2:
        return float("nan")
    split = np.asarray(
        [part for chain in chains for part in (chain[:half], chain[-half:])],
        dtype=float,
    )
    n = split.shape[1]
    within = float(np.mean(np.var(split, axis=1, ddof=1)))
    if within <= 0:
        return 1.0
    between = float(n * np.var(np.mean(split, axis=1), ddof=1))
    variance = ((n - 1) / n) * within + between / n
    return float(np.sqrt(variance / within))


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite frozen output: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    design = plan["design"]
    people = [f"P{index:03d}" for index in range(1, int(design["persons"]) + 1)]
    raters = [f"R{index:02d}" for index in range(1, int(design["raters"]) + 1)]
    rng = np.random.default_rng(int(design["person_seed"]))
    abilities = rng.normal(size=len(people))
    abilities = (abilities - abilities.mean()) / abilities.std(ddof=0) * 0.8
    person_coordinates = dict(zip(people, abilities, strict=True))
    rater_coordinates = dict(
        zip(raters, map(float, design["rater_coordinates"]), strict=True)
    )
    initial_edges = _planned_edges(people, raters)
    source_rows = _exchangeable_rows(initial_edges)

    start_edges: list[set[tuple[str, str]]] = []
    start_rows: list[dict[str, object]] = []
    for start_number, seed in enumerate(plan["chains"]["start_seeds"], start=1):
        start_chain = sample_degree_conditioned_assignment_chain(
            initial_edges=initial_edges,
            person_coordinates=person_coordinates,
            rater_coordinates=rater_coordinates,
            gamma=float(plan["chains"]["start_generation_gamma"]),
            n_samples=1,
            burnin=int(plan["chains"]["start_burnin"]),
            thin=1,
            seed=int(seed),
            require_connected=True,
        )
        start_edges.append(start_chain["final_edges"])
        start_rows.append(
            {
                "Start": start_number,
                "Seed": int(seed),
                "EdgeSignature": "|".join(
                    f"{person}::{rater}"
                    for person, rater in sorted(start_chain["final_edges"])
                ),
            }
        )
    if len({row["EdgeSignature"] for row in start_rows}) != len(start_rows):
        raise RuntimeError("Independent start generation returned duplicate graph states.")

    diagnostic_parts: list[pd.DataFrame] = []
    trace_parts: list[pd.DataFrame] = []
    materialization_parts: list[pd.DataFrame] = []
    traces_by_gamma: dict[float, list[np.ndarray]] = {}
    for gamma_index, gamma_value in enumerate(plan["candidate_gammas"]):
        gamma = float(gamma_value)
        traces_by_gamma[gamma] = []
        for chain_index, edges in enumerate(start_edges, start=1):
            seed = int(plan["chains"]["chain_seed_base"]) + gamma_index * 100 + chain_index
            chain_id = f"gamma_{gamma:+.1f}_chain_{chain_index}"
            chain = sample_degree_conditioned_assignment_chain(
                initial_edges=edges,
                person_coordinates=person_coordinates,
                rater_coordinates=rater_coordinates,
                gamma=gamma,
                n_samples=int(plan["chains"]["samples"]),
                burnin=int(plan["chains"]["burnin"]),
                thin=int(plan["chains"]["thin"]),
                seed=seed,
                require_connected=True,
            )
            diagnostics = chain["diagnostics"].copy()
            diagnostics.insert(0, "ChainId", chain_id)
            diagnostics["Available"] = bool(chain["available"])
            diagnostic_parts.append(diagnostics)
            trace = chain["samples"].drop(columns="EdgeSignature").copy()
            trace.insert(0, "ChainId", chain_id)
            trace.insert(1, "Gamma", gamma)
            trace_parts.append(trace)
            traces_by_gamma[gamma].append(
                chain["samples"]["Statistic"].to_numpy(dtype=float)
            )
            materialized = materialize_exchangeable_assignment_design(
                source_rows,
                chain["final_edges"],
                facet_cols=["Rater", "Task", "Criterion"],
                context_cols=["Task", "Criterion"],
                rater_mapping_confirmed=True,
            )
            audit = materialized["invariants"].copy()
            audit.insert(0, "ChainId", chain_id)
            audit["ScoreColumnAbsent"] = "Score" not in materialized["design"].columns
            materialization_parts.append(audit)

    diagnostics_table = pd.concat(diagnostic_parts, ignore_index=True)
    trace_table = pd.concat(trace_parts, ignore_index=True)
    materialization_table = pd.concat(materialization_parts, ignore_index=True)
    gamma_rows: list[dict[str, object]] = []
    for gamma in sorted(traces_by_gamma):
        subset = trace_table.loc[trace_table["Gamma"].eq(gamma)]
        gamma_rows.append(
            {
                "Gamma": gamma,
                "Chains": len(traces_by_gamma[gamma]),
                "SamplesPerChain": min(map(len, traces_by_gamma[gamma])),
                "SplitRhatStatistic": _split_rhat(traces_by_gamma[gamma]),
                "MeanStatistic": float(subset["Statistic"].mean()),
                "MeanAssignmentCorrelation": float(
                    subset["AssignmentCorrelation"].mean()
                ),
                "SDAssignmentCorrelationAcrossSamples": float(
                    subset["AssignmentCorrelation"].std(ddof=1)
                ),
            }
        )
    gamma_summary = pd.DataFrame(gamma_rows).sort_values("Gamma").reset_index(drop=True)

    gates = plan["qualification_gates"]
    base_checks = {
        "all_chains_available": bool(diagnostics_table["Available"].all()),
        "degree_margins": bool(diagnostics_table["DegreeMarginsPreserved"].all()),
        "connectivity": bool(
            diagnostics_table["FinalOverlapComponents"].eq(
                gates["final_overlap_components"]
            ).all()
        ),
        "numerical_residual": bool(
            diagnostics_table["MaximumStatisticUpdateResidual"].le(
                gates["maximum_statistic_update_residual"]
            ).all()
        ),
        "split_rhat": bool(
            gamma_summary["SplitRhatStatistic"].le(gates["split_rhat_max"]).all()
        ),
        "ess": bool(
            diagnostics_table["StatisticESS"].ge(
                gates["minimum_per_chain_statistic_ess"]
            ).all()
        ),
        "admissible_acceptance": bool(
            diagnostics_table["AcceptanceRateAdmissibleProposals"].ge(
                gates["minimum_admissible_acceptance_rate"]
            ).all()
        ),
        "neutral_correlation": bool(
            abs(
                float(
                    gamma_summary.loc[
                        gamma_summary["Gamma"].eq(0.0),
                        "MeanAssignmentCorrelation",
                    ].iloc[0]
                )
            ) <= gates["neutral_absolute_mean_assignment_correlation_max"]
        ),
        "monotone_correlation": bool(
            np.all(
                np.diff(
                    gamma_summary["MeanAssignmentCorrelation"].to_numpy(dtype=float)
                ) > 0
            )
        ),
        "materialization": bool(materialization_table["Passed"].all()),
        "score_absent": bool(materialization_table["ScoreColumnAbsent"].all()),
    }

    selection = plan["dose_selection_rule"]
    selected_gamma: float | None = None
    dose_rows: list[dict[str, object]] = []
    for magnitude in sorted(
        map(float, selection["eligible_absolute_gammas"]), reverse=True
    ):
        negative = gamma_summary.loc[gamma_summary["Gamma"].eq(-magnitude)].iloc[0]
        positive = gamma_summary.loc[gamma_summary["Gamma"].eq(magnitude)].iloc[0]
        pair_diagnostics = diagnostics_table.loc[
            diagnostics_table["Gamma"].abs().eq(magnitude)
        ]
        correlations = [
            abs(float(negative["MeanAssignmentCorrelation"])),
            abs(float(positive["MeanAssignmentCorrelation"])),
        ]
        checks = {
            "split_rhat": max(
                float(negative["SplitRhatStatistic"]),
                float(positive["SplitRhatStatistic"]),
            ) <= gates["split_rhat_max"],
            "ess": float(pair_diagnostics["StatisticESS"].min())
            >= gates["minimum_per_chain_statistic_ess"],
            "acceptance": float(
                pair_diagnostics["AcceptanceRateAdmissibleProposals"].min()
            ) >= gates["minimum_admissible_acceptance_rate"],
            "correlation_range": min(correlations)
            >= selection["minimum_absolute_mean_assignment_correlation"]
            and max(correlations)
            <= selection["maximum_absolute_mean_assignment_correlation"],
            "signs": float(negative["MeanAssignmentCorrelation"]) < 0
            < float(positive["MeanAssignmentCorrelation"]),
            "symmetry": abs(
                float(negative["MeanAssignmentCorrelation"])
                + float(positive["MeanAssignmentCorrelation"])
            ) <= selection["positive_and_negative_symmetry_error_max"],
        }
        eligible = all(checks.values())
        dose_rows.append(
            {
                "AbsoluteGamma": magnitude,
                "NegativeMeanCorrelation": float(
                    negative["MeanAssignmentCorrelation"]
                ),
                "PositiveMeanCorrelation": float(
                    positive["MeanAssignmentCorrelation"]
                ),
                "SymmetryError": abs(
                    float(negative["MeanAssignmentCorrelation"])
                    + float(positive["MeanAssignmentCorrelation"])
                ),
                **{f"Gate_{key}": bool(value) for key, value in checks.items()},
                "Eligible": eligible,
            }
        )
        if selected_gamma is None and eligible:
            selected_gamma = magnitude
    dose_table = pd.DataFrame(dose_rows)
    base_checks["dose_selected"] = selected_gamma is not None
    assessment = {
        "schema_version": "known_assignment_large_design_assessment_v1",
        "decision": "pass" if all(base_checks.values()) else "fail",
        "selected_absolute_gamma": selected_gamma,
        "checks": base_checks,
        "minimum_statistic_ess": float(diagnostics_table["StatisticESS"].min()),
        "maximum_split_rhat": float(gamma_summary["SplitRhatStatistic"].max()),
        "minimum_admissible_acceptance_rate": float(
            diagnostics_table["AcceptanceRateAdmissibleProposals"].min()
        ),
        "claim_boundary": plan["claim_boundary"],
    }

    OUTPUT_DIR.mkdir(parents=False, exist_ok=False)
    artifacts = {
        "person_coordinates.csv": pd.DataFrame(
            {"Person": people, "Theta": abilities}
        ),
        "start_states.csv": pd.DataFrame(start_rows),
        "chain_diagnostics.csv": diagnostics_table,
        "chain_trace.csv": trace_table,
        "gamma_summary.csv": gamma_summary,
        "dose_selection.csv": dose_table,
        "materialization_audit.csv": materialization_table,
    }
    for name, table in artifacts.items():
        table.to_csv(OUTPUT_DIR / name, index=False)
    assessment_path = OUTPUT_DIR / "assessment.json"
    assessment_path.write_text(
        json.dumps(assessment, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    script_path = Path(__file__).resolve()
    manifest = {
        "schema_version": "known_assignment_large_design_manifest_v1",
        "mechanism_schema_version": MECHANISM_SCHEMA_VERSION,
        "plan_sha256": _sha256(PLAN_PATH),
        "runner_sha256": _sha256(script_path),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "artifacts": {
            name: _sha256(OUTPUT_DIR / name)
            for name in [*artifacts, assessment_path.name]
        },
    }
    (OUTPUT_DIR / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(assessment, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
