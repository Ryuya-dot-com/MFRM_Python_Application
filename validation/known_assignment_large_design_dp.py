"""Exact-DP remediation and outcome-blind assignment-dose selection."""

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
    _initial_positive_sequence_ess,
    materialize_exchangeable_assignment_design,
    sample_degree_conditioned_assignment_dp,
)
from validation.known_assignment_large_design_calibration import (
    _exchangeable_rows,
    _planned_edges,
    _split_rhat,
)


PLAN_PATH = ROOT / "validation" / "known_assignment_large_design_dp_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "known_assignment_large_design_dp_20260811"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _resolve(relative: str) -> Path:
    return ROOT / Path(relative)


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite frozen output: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    locked = plan["locked_inputs"]
    if _sha256(_resolve(locked["base_plan"])).lower() != locked["base_plan_sha256"]:
        raise ValueError("Base plan hash mismatch")
    failure_path = _resolve(plan["trigger"]["failed_heatbath_assessment"])
    if _sha256(failure_path).lower() != plan["trigger"]["failed_heatbath_assessment_sha256"]:
        raise ValueError("Failed heat-bath assessment hash mismatch")
    base = json.loads(_resolve(locked["base_plan"]).read_text(encoding="utf-8"))
    coordinates = pd.read_csv(_resolve(locked["person_coordinates"]))
    person_coordinates = dict(
        coordinates[["Person", "Theta"]].itertuples(index=False, name=None)
    )
    people = sorted(person_coordinates)
    raters = [f"R{index:02d}" for index in range(1, 5)]
    rater_coordinates = dict(
        zip(raters, map(float, base["design"]["rater_coordinates"]), strict=True)
    )
    person_degrees = {person: 2 for person in people}
    rater_degrees = {rater: 40 for rater in raters}
    source_edges = _planned_edges(people, raters)
    source_rows = _exchangeable_rows(source_edges)
    exact = plan["exact_sampling"]
    batches = int(exact["batches_per_gamma"])
    batch_size = int(exact["independent_draws_per_batch"])
    total_samples = batches * batch_size

    diagnostics_parts: list[pd.DataFrame] = []
    trace_parts: list[pd.DataFrame] = []
    materialization_parts: list[pd.DataFrame] = []
    batch_rows: list[dict[str, object]] = []
    for gamma_index, gamma_value in enumerate(locked["candidate_gammas"]):
        gamma = float(gamma_value)
        sampled = sample_degree_conditioned_assignment_dp(
            person_degrees=person_degrees,
            rater_degrees=rater_degrees,
            person_coordinates=person_coordinates,
            rater_coordinates=rater_coordinates,
            gamma=gamma,
            n_samples=total_samples,
            seed=int(exact["seed_base"]) + gamma_index,
            require_connected=True,
            max_dp_states=int(exact["max_dp_states"]),
            max_total_rejections=int(exact["max_total_connected_rejections"]),
            retain_edge_signatures=bool(exact["edge_signatures_retained"]),
        )
        diagnostic = sampled["diagnostics"].copy()
        diagnostic.insert(0, "RunId", f"gamma_{gamma:+.1f}")
        diagnostic["Available"] = bool(sampled["available"])
        diagnostics_parts.append(diagnostic)
        trace = sampled["samples"].copy()
        trace.insert(0, "Gamma", gamma)
        trace["Batch"] = (np.arange(len(trace), dtype=int) // batch_size) + 1
        trace["BatchSample"] = (np.arange(len(trace), dtype=int) % batch_size) + 1
        trace_parts.append(trace)
        for batch, batch_frame in trace.groupby("Batch", sort=True):
            values = batch_frame["Statistic"].to_numpy(dtype=float)
            batch_rows.append(
                {
                    "Gamma": gamma,
                    "Batch": int(batch),
                    "Samples": len(batch_frame),
                    "MeanStatistic": float(np.mean(values)),
                    "MeanAssignmentCorrelation": float(
                        batch_frame["AssignmentCorrelation"].mean()
                    ),
                    "StatisticESSDiagnostic": _initial_positive_sequence_ess(values),
                }
            )
        materialized = materialize_exchangeable_assignment_design(
            source_rows,
            sampled["final_edges"],
            facet_cols=["Rater", "Task", "Criterion"],
            context_cols=["Task", "Criterion"],
            rater_mapping_confirmed=True,
        )
        audit = materialized["invariants"].copy()
        audit.insert(0, "RunId", f"gamma_{gamma:+.1f}")
        audit["ScoreColumnAbsent"] = "Score" not in materialized["design"].columns
        materialization_parts.append(audit)

    diagnostics = pd.concat(diagnostics_parts, ignore_index=True)
    traces = pd.concat(trace_parts, ignore_index=True)
    batches_table = pd.DataFrame(batch_rows)
    materialization = pd.concat(materialization_parts, ignore_index=True)
    gamma_rows: list[dict[str, object]] = []
    for gamma, group in traces.groupby("Gamma", sort=True):
        split = [
            part["Statistic"].to_numpy(dtype=float)
            for _, part in group.groupby("Batch", sort=True)
        ]
        gamma_rows.append(
            {
                "Gamma": float(gamma),
                "IndependentSamples": len(group),
                "SplitRhatStatistic": _split_rhat(split),
                "MeanStatistic": float(group["Statistic"].mean()),
                "MeanAssignmentCorrelation": float(
                    group["AssignmentCorrelation"].mean()
                ),
                "SDAssignmentCorrelation": float(
                    group["AssignmentCorrelation"].std(ddof=1)
                ),
            }
        )
    gamma_summary = pd.DataFrame(gamma_rows).sort_values("Gamma").reset_index(drop=True)

    gates = plan["qualification_gates"]
    checks = {
        "all_runs_available": bool(diagnostics["Available"].all()),
        "exact_independent": bool(diagnostics["ExactIndependentSamples"].all()),
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
        "dp_cap": bool(
            diagnostics["DPStatesEvaluated"].lt(diagnostics["DPStateCap"]).all()
        ),
        "connected_acceptance": bool(
            diagnostics["ConnectedAcceptanceRate"].ge(
                gates["minimum_connected_acceptance_rate"]
            ).all()
        ),
        "split_rhat": bool(
            gamma_summary["SplitRhatStatistic"].le(gates["split_rhat_max"]).all()
        ),
        "batch_ess": bool(
            batches_table["StatisticESSDiagnostic"].ge(
                gates["minimum_batch_statistic_ess_diagnostic"]
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
        correlations = [
            abs(float(negative["MeanAssignmentCorrelation"])),
            abs(float(positive["MeanAssignmentCorrelation"])),
        ]
        pair_checks = {
            "split_rhat": max(
                float(negative["SplitRhatStatistic"]),
                float(positive["SplitRhatStatistic"]),
            ) <= gates["split_rhat_max"],
            "batch_ess": float(
                batches_table.loc[
                    batches_table["Gamma"].abs().eq(magnitude),
                    "StatisticESSDiagnostic",
                ].min()
            ) >= gates["minimum_batch_statistic_ess_diagnostic"],
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
        "schema_version": "known_assignment_large_design_dp_assessment_v1",
        "decision": "pass" if all(checks.values()) else "fail",
        "selected_absolute_gamma": selected_gamma,
        "checks": checks,
        "minimum_batch_ess_diagnostic": float(
            batches_table["StatisticESSDiagnostic"].min()
        ),
        "maximum_split_rhat": float(gamma_summary["SplitRhatStatistic"].max()),
        "maximum_dp_states": int(diagnostics["DPStatesEvaluated"].max()),
        "minimum_connected_acceptance_rate": float(
            diagnostics["ConnectedAcceptanceRate"].min()
        ),
        "claim_boundary": plan["claim_boundary"],
    }

    OUTPUT_DIR.mkdir(parents=False, exist_ok=False)
    artifacts = {
        "run_diagnostics.csv": diagnostics,
        "independent_trace.csv": traces,
        "batch_diagnostics.csv": batches_table,
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
        "schema_version": "known_assignment_large_design_dp_manifest_v1",
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
