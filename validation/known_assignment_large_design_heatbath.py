"""Frozen heat-bath remediation of failed large-design assignment mixing."""

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
    sample_degree_conditioned_assignment_heatbath_chain,
)
from validation.known_assignment_large_design_calibration import (
    _exchangeable_rows,
    _planned_edges,
    _split_rhat,
)


PLAN_PATH = ROOT / "validation" / "known_assignment_large_design_heatbath_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "known_assignment_large_design_heatbath_20260811"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _resolve(relative: str) -> Path:
    return ROOT / Path(relative)


def _parse_edges(signature: str) -> set[tuple[str, str]]:
    return {
        tuple(value.split("::", maxsplit=1))
        for value in str(signature).split("|")
    }


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite frozen output: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    locked = plan["locked_inputs"]
    for key, hash_key in (
        ("parent_plan", "parent_plan_sha256"),
        ("start_states", "start_states_sha256"),
    ):
        path = _resolve(locked[key])
        if _sha256(path).lower() != str(locked[hash_key]).lower():
            raise ValueError(f"Locked remediation input hash mismatch: {key}")
    failed_path = _resolve(plan["trigger"]["failed_assessment"])
    if _sha256(failed_path).lower() != plan["trigger"]["failed_assessment_sha256"]:
        raise ValueError("Failed parent assessment hash mismatch")

    parent = json.loads(_resolve(locked["parent_plan"]).read_text(encoding="utf-8"))
    coordinate_table = pd.read_csv(_resolve(locked["person_coordinates"]))
    person_coordinates = dict(
        coordinate_table[["Person", "Theta"]].itertuples(index=False, name=None)
    )
    raters = [f"R{index:02d}" for index in range(1, 5)]
    rater_coordinates = dict(
        zip(raters, map(float, parent["design"]["rater_coordinates"]), strict=True)
    )
    people = sorted(person_coordinates)
    initial_edges = _planned_edges(people, raters)
    source_rows = _exchangeable_rows(initial_edges)
    starts = pd.read_csv(_resolve(locked["start_states"]))
    start_edges = [_parse_edges(value) for value in starts["EdgeSignature"]]

    diagnostic_parts: list[pd.DataFrame] = []
    trace_parts: list[pd.DataFrame] = []
    materialization_parts: list[pd.DataFrame] = []
    traces_by_gamma: dict[float, list[np.ndarray]] = {}
    for gamma_index, gamma_value in enumerate(locked["candidate_gammas"]):
        gamma = float(gamma_value)
        traces_by_gamma[gamma] = []
        for chain_index, edges in enumerate(start_edges, start=1):
            seed = int(locked["chain_seed_base"]) + gamma_index * 100 + chain_index
            chain_id = f"gamma_{gamma:+.1f}_chain_{chain_index}"
            chain = sample_degree_conditioned_assignment_heatbath_chain(
                initial_edges=edges,
                person_coordinates=person_coordinates,
                rater_coordinates=rater_coordinates,
                gamma=gamma,
                n_samples=int(locked["samples"]),
                burnin=int(locked["burnin"]),
                thin=int(locked["thin"]),
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

    diagnostics = pd.concat(diagnostic_parts, ignore_index=True)
    traces = pd.concat(trace_parts, ignore_index=True)
    materialization = pd.concat(materialization_parts, ignore_index=True)
    gamma_rows: list[dict[str, object]] = []
    for gamma in sorted(traces_by_gamma):
        subset = traces.loc[traces["Gamma"].eq(gamma)]
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
    checks = {
        "all_chains_available": bool(diagnostics["Available"].all()),
        "degree_margins": bool(diagnostics["DegreeMarginsPreserved"].all()),
        "connectivity": bool(
            diagnostics["FinalOverlapComponents"].eq(
                gates["final_overlap_components"]
            ).all()
        ),
        "numerical_residual": bool(
            diagnostics["MaximumStatisticUpdateResidual"].le(
                gates["maximum_statistic_update_residual"]
            ).all()
        ),
        "split_rhat": bool(
            gamma_summary["SplitRhatStatistic"].le(gates["split_rhat_max"]).all()
        ),
        "ess": bool(
            diagnostics["StatisticESS"].ge(
                gates["minimum_per_chain_statistic_ess"]
            ).all()
        ),
        "movement": bool(
            diagnostics["MovementRateAllUpdates"].ge(
                gates["minimum_all_update_movement_rate"]
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
        "materialization": bool(materialization["Passed"].all()),
        "score_absent": bool(materialization["ScoreColumnAbsent"].all()),
    }

    selection = plan["dose_selection_rule"]
    selected_gamma: float | None = None
    dose_rows: list[dict[str, object]] = []
    for magnitude in sorted(
        map(float, selection["eligible_absolute_gammas"]), reverse=True
    ):
        negative = gamma_summary.loc[gamma_summary["Gamma"].eq(-magnitude)].iloc[0]
        positive = gamma_summary.loc[gamma_summary["Gamma"].eq(magnitude)].iloc[0]
        pair_diagnostics = diagnostics.loc[diagnostics["Gamma"].abs().eq(magnitude)]
        correlations = [
            abs(float(negative["MeanAssignmentCorrelation"])),
            abs(float(positive["MeanAssignmentCorrelation"])),
        ]
        pair_checks = {
            "split_rhat": max(
                float(negative["SplitRhatStatistic"]),
                float(positive["SplitRhatStatistic"]),
            ) <= gates["split_rhat_max"],
            "ess": float(pair_diagnostics["StatisticESS"].min())
            >= gates["minimum_per_chain_statistic_ess"],
            "movement": float(pair_diagnostics["MovementRateAllUpdates"].min())
            >= gates["minimum_all_update_movement_rate"],
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
        eligible = all(pair_checks.values())
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
                **{f"Gate_{key}": bool(value) for key, value in pair_checks.items()},
                "Eligible": eligible,
            }
        )
        if selected_gamma is None and eligible:
            selected_gamma = magnitude
    checks["dose_selected"] = selected_gamma is not None
    dose_selection = pd.DataFrame(dose_rows)
    assessment = {
        "schema_version": "known_assignment_large_design_heatbath_assessment_v1",
        "decision": "pass" if all(checks.values()) else "fail",
        "selected_absolute_gamma": selected_gamma,
        "checks": checks,
        "minimum_statistic_ess": float(diagnostics["StatisticESS"].min()),
        "maximum_split_rhat": float(gamma_summary["SplitRhatStatistic"].max()),
        "minimum_movement_rate": float(diagnostics["MovementRateAllUpdates"].min()),
        "claim_boundary": plan["claim_boundary"],
    }

    OUTPUT_DIR.mkdir(parents=False, exist_ok=False)
    artifacts = {
        "chain_diagnostics.csv": diagnostics,
        "chain_trace.csv": traces,
        "gamma_summary.csv": gamma_summary,
        "dose_selection.csv": dose_selection,
        "materialization_audit.csv": materialization,
    }
    for name, table in artifacts.items():
        table.to_csv(OUTPUT_DIR / name, index=False)
    assessment_path = OUTPUT_DIR / "assessment.json"
    assessment_path.write_text(
        json.dumps(assessment, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    dependencies = {
        "plan": PLAN_PATH,
        "runner": Path(__file__).resolve(),
        "mechanism": ROOT / "mfrm_app" / "assignment_mechanism.py",
        "parent_runner_helpers": ROOT
        / "validation"
        / "known_assignment_large_design_calibration.py",
    }
    manifest = {
        "schema_version": "known_assignment_large_design_heatbath_manifest_v1",
        "mechanism_schema_version": MECHANISM_SCHEMA_VERSION,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "dependency_sha256": {
            name: _sha256(path) for name, path in dependencies.items()
        },
        "artifact_sha256": {
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
