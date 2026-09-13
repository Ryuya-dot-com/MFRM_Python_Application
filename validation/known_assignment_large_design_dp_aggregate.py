"""Aggregate immutable exact-DP gamma shards under the frozen gates."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

from validation.known_assignment_large_design_calibration import _split_rhat
from validation.known_assignment_large_design_dp_shard import _sha256, _slug


PLAN_PATH = ROOT / "validation" / "known_assignment_large_design_dp_plan_20260811.json"
AMENDMENT_PATH = ROOT / "validation" / "known_assignment_large_design_dp_execution_amendment_20260811.json"
SHARD_ROOT = ROOT / "validation" / "known_assignment_large_design_dp_shards_20260811"
OUTPUT_DIR = ROOT / "validation" / "known_assignment_large_design_dp_20260811"


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite frozen aggregate: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    gammas = [float(value) for value in plan["locked_inputs"]["candidate_gammas"]]
    diagnostics_parts: list[pd.DataFrame] = []
    trace_parts: list[pd.DataFrame] = []
    batch_parts: list[pd.DataFrame] = []
    materialization_parts: list[pd.DataFrame] = []
    shard_manifests: dict[str, str] = {}
    for gamma in gammas:
        root = SHARD_ROOT / _slug(gamma)
        manifest_path = root / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if float(manifest["gamma"]) != gamma:
            raise ValueError(f"Shard gamma mismatch: {root}")
        for name, expected in manifest["artifact_sha256"].items():
            if _sha256(root / name) != expected:
                raise ValueError(f"Shard artifact hash mismatch: {root / name}")
        diagnostics_parts.append(pd.read_csv(root / "run_diagnostics.csv"))
        trace_parts.append(pd.read_csv(root / "independent_trace.csv"))
        batch_parts.append(pd.read_csv(root / "batch_diagnostics.csv"))
        materialization_parts.append(pd.read_csv(root / "materialization_audit.csv"))
        shard_manifests[_slug(gamma)] = _sha256(manifest_path)
    diagnostics = pd.concat(diagnostics_parts, ignore_index=True)
    traces = pd.concat(trace_parts, ignore_index=True)
    batches = pd.concat(batch_parts, ignore_index=True)
    materialization = pd.concat(materialization_parts, ignore_index=True)

    gamma_rows: list[dict[str, object]] = []
    for gamma, group in traces.groupby("Gamma", sort=True):
        split = [
            frame["Statistic"].to_numpy(dtype=float)
            for _, frame in group.groupby("Batch", sort=True)
        ]
        gamma_rows.append(
            {
                "Gamma": float(gamma),
                "IndependentSamples": len(group),
                "SplitRhatStatistic": _split_rhat(split),
                "MeanStatistic": float(group["Statistic"].mean()),
                "MeanAssignmentCorrelation": float(group["AssignmentCorrelation"].mean()),
                "SDAssignmentCorrelation": float(group["AssignmentCorrelation"].std(ddof=1)),
            }
        )
    gamma_summary = pd.DataFrame(gamma_rows).sort_values("Gamma").reset_index(drop=True)
    gates = plan["qualification_gates"]
    checks = {
        "all_runs_available": bool(diagnostics["Available"].all()),
        "exact_independent": bool(diagnostics["ExactIndependentSamples"].all()),
        "degree_margins": bool(diagnostics["DegreeMarginsPreserved"].all()),
        "connectivity": bool(diagnostics["FinalOverlapComponents"].eq(gates["final_overlap_components"]).all()),
        "numerical_residual": bool(diagnostics["MaximumStatisticUpdateResidual"].le(gates["maximum_statistic_update_residual"]).all()),
        "dp_cap": bool(diagnostics["DPStatesEvaluated"].lt(diagnostics["DPStateCap"]).all()),
        "connected_acceptance": bool(diagnostics["ConnectedAcceptanceRate"].ge(gates["minimum_connected_acceptance_rate"]).all()),
        "split_rhat": bool(gamma_summary["SplitRhatStatistic"].le(gates["split_rhat_max"]).all()),
        "batch_ess": bool(batches["StatisticESSDiagnostic"].ge(gates["minimum_batch_statistic_ess_diagnostic"]).all()),
        "neutral_correlation": bool(abs(float(gamma_summary.loc[gamma_summary["Gamma"].eq(0.0), "MeanAssignmentCorrelation"].iloc[0])) <= gates["neutral_absolute_mean_assignment_correlation_max"]),
        "monotone_correlation": bool(np.all(np.diff(gamma_summary["MeanAssignmentCorrelation"].to_numpy(dtype=float)) > 0)),
        "materialization": bool(materialization["Passed"].astype(bool).all()),
        "score_absent": bool(materialization["ScoreColumnAbsent"].astype(bool).all()),
    }

    selection = plan["dose_selection_rule"]
    selected_gamma: float | None = None
    dose_rows: list[dict[str, object]] = []
    for magnitude in sorted(map(float, selection["eligible_absolute_gammas"]), reverse=True):
        negative = gamma_summary.loc[gamma_summary["Gamma"].eq(-magnitude)].iloc[0]
        positive = gamma_summary.loc[gamma_summary["Gamma"].eq(magnitude)].iloc[0]
        correlations = [abs(float(negative["MeanAssignmentCorrelation"])), abs(float(positive["MeanAssignmentCorrelation"]))]
        pair_checks = {
            "split_rhat": max(float(negative["SplitRhatStatistic"]), float(positive["SplitRhatStatistic"])) <= gates["split_rhat_max"],
            "batch_ess": float(batches.loc[batches["Gamma"].abs().eq(magnitude), "StatisticESSDiagnostic"].min()) >= gates["minimum_batch_statistic_ess_diagnostic"],
            "correlation_range": min(correlations) >= selection["minimum_absolute_mean_assignment_correlation"] and max(correlations) <= selection["maximum_absolute_mean_assignment_correlation"],
            "signs": float(negative["MeanAssignmentCorrelation"]) < 0 < float(positive["MeanAssignmentCorrelation"]),
            "symmetry": abs(float(negative["MeanAssignmentCorrelation"]) + float(positive["MeanAssignmentCorrelation"])) <= selection["positive_and_negative_symmetry_error_max"],
        }
        eligible = all(pair_checks.values())
        dose_rows.append(
            {
                "AbsoluteGamma": magnitude,
                "NegativeMeanCorrelation": float(negative["MeanAssignmentCorrelation"]),
                "PositiveMeanCorrelation": float(positive["MeanAssignmentCorrelation"]),
                "SymmetryError": abs(float(negative["MeanAssignmentCorrelation"]) + float(positive["MeanAssignmentCorrelation"])),
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
        "minimum_batch_ess_diagnostic": float(batches["StatisticESSDiagnostic"].min()),
        "maximum_split_rhat": float(gamma_summary["SplitRhatStatistic"].max()),
        "maximum_dp_states": int(diagnostics["DPStatesEvaluated"].max()),
        "minimum_connected_acceptance_rate": float(diagnostics["ConnectedAcceptanceRate"].min()),
        "claim_boundary": plan["claim_boundary"],
    }

    OUTPUT_DIR.mkdir(parents=False, exist_ok=False)
    artifacts = {
        "run_diagnostics.csv": diagnostics,
        "independent_trace.csv": traces,
        "batch_diagnostics.csv": batches,
        "gamma_summary.csv": gamma_summary,
        "dose_selection.csv": dose_selection,
        "materialization_audit.csv": materialization,
    }
    for name, table in artifacts.items():
        table.to_csv(OUTPUT_DIR / name, index=False)
    assessment_path = OUTPUT_DIR / "assessment.json"
    assessment_path.write_text(json.dumps(assessment, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    dependencies = {
        "plan": PLAN_PATH,
        "amendment": AMENDMENT_PATH,
        "aggregator": Path(__file__).resolve(),
        "shard_runner": ROOT / "validation" / "known_assignment_large_design_dp_shard.py",
        "mechanism": ROOT / "mfrm_app" / "assignment_mechanism.py",
    }
    manifest = {
        "schema_version": "known_assignment_large_design_dp_aggregate_manifest_v1",
        "dependency_sha256": {name: _sha256(path) for name, path in dependencies.items()},
        "shard_manifest_sha256": shard_manifests,
        "artifact_sha256": {name: _sha256(OUTPUT_DIR / name) for name in [*artifacts, assessment_path.name]},
    }
    (OUTPUT_DIR / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(assessment, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
