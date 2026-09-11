"""Frozen equivalence and performance qualification of the dense exact DP."""

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
    compare_chain_to_exact_oracle,
    enumerate_degree_conditioned_assignments,
    sample_degree_conditioned_assignment_dense_dp,
    sample_degree_conditioned_assignment_dp,
)


PLAN_PATH = ROOT / "validation" / "known_assignment_dense_dp_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "known_assignment_dense_dp_20260811"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify_locked_dependencies(plan: dict[str, object]) -> None:
    locked = plan["locked_dependencies"]
    assert isinstance(locked, dict)
    pairs = (
        ("mechanism", "mechanism_sha256"),
        ("retained_recursive_diagnostics", "retained_recursive_diagnostics_sha256"),
        ("person_coordinates", "person_coordinates_sha256"),
        ("base_plan", "base_plan_sha256"),
    )
    for path_key, hash_key in pairs:
        path = ROOT / str(locked[path_key])
        if _sha256(path).lower() != str(locked[hash_key]).lower():
            raise ValueError(f"Locked dependency hash mismatch: {path_key}")


def _small_contract() -> tuple[dict[str, int], dict[str, int], dict[str, float], dict[str, float]]:
    people = [f"p{index}" for index in range(6)]
    raters = [f"r{index}" for index in range(3)]
    return (
        {person: 2 for person in people},
        {rater: 4 for rater in raters},
        {person: float(index) for index, person in enumerate(people)},
        {rater: float(index) for index, rater in enumerate(raters)},
    )


def _degree_audit(
    edges: set[tuple[str, str]],
    person_degrees: dict[str, int],
    rater_degrees: dict[str, int],
) -> bool:
    return bool(
        all(sum(person == edge[0] for edge in edges) == degree for person, degree in person_degrees.items())
        and all(sum(rater == edge[1] for edge in edges) == degree for rater, degree in rater_degrees.items())
    )


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite frozen output: {OUTPUT_DIR}")
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    _verify_locked_dependencies(plan)

    small = plan["small_oracle"]
    small_gates = plan["registered_gates"]["small_oracle"]
    p_degree, r_degree, p_coordinate, r_coordinate = _small_contract()
    small_rows: list[dict[str, object]] = []
    for gamma, seed in zip(small["gammas"], small["seeds"], strict=True):
        gamma = float(gamma)
        oracle = enumerate_degree_conditioned_assignments(
            person_degrees=p_degree,
            rater_degrees=r_degree,
            person_coordinates=p_coordinate,
            rater_coordinates=r_coordinate,
            gamma=gamma,
            require_connected=True,
        )
        recursive = sample_degree_conditioned_assignment_dp(
            person_degrees=p_degree,
            rater_degrees=r_degree,
            person_coordinates=p_coordinate,
            rater_coordinates=r_coordinate,
            gamma=gamma,
            n_samples=1,
            seed=int(seed) + 100,
            require_connected=True,
        )
        dense = sample_degree_conditioned_assignment_dense_dp(
            person_degrees=p_degree,
            rater_degrees=r_degree,
            person_coordinates=p_coordinate,
            rater_coordinates=r_coordinate,
            gamma=gamma,
            n_samples=int(small["independent_draws_per_gamma"]),
            seed=int(seed),
            require_connected=True,
        )
        comparison = compare_chain_to_exact_oracle(dense, oracle)["summary"].iloc[0]
        dense_diagnostic = dense["diagnostics"].iloc[0]
        recursive_logz = float(recursive["diagnostics"].iloc[0]["DPLogNormalizerUnconditioned"])
        dense_logz = float(dense_diagnostic["DPLogNormalizerUnconditioned"])
        row = {
            "Gamma": gamma,
            "OracleStates": int(oracle["summary"].iloc[0]["States"]),
            "Available": bool(dense["available"]),
            "RecursiveLogNormalizer": recursive_logz,
            "DenseLogNormalizer": dense_logz,
            "AbsoluteLogNormalizerDifference": abs(recursive_logz - dense_logz),
            "ProbabilityOutsideOracle": float(comparison["ChainProbabilityOutsideOracle"]),
            "TotalVariationDistance": float(comparison["TotalVariationDistance"]),
            "MeanStatisticError": float(comparison["MeanStatisticError"]),
            "SampledStates": int(comparison["SampledStates"]),
            "MaximumConditionalLogResidual": float(dense_diagnostic["MaximumConditionalLogResidual"]),
            "MaximumStatisticUpdateResidual": float(dense_diagnostic["MaximumStatisticUpdateResidual"]),
        }
        row["Passed"] = bool(
            row["Available"] is small_gates["available"]
            and row["OracleStates"] == small_gates["states_each_gamma"]
            and row["AbsoluteLogNormalizerDifference"] <= small_gates["recursive_dense_log_normalizer_absolute_difference_max"]
            and row["ProbabilityOutsideOracle"] <= small_gates["probability_outside_oracle_max"]
            and row["TotalVariationDistance"] <= small_gates["total_variation_distance_max"]
            and abs(row["MeanStatisticError"]) <= small_gates["absolute_mean_statistic_error_max"]
            and row["SampledStates"] >= small_gates["unique_sampled_states_min"]
            and row["MaximumConditionalLogResidual"] <= small_gates["maximum_conditional_log_residual_max"]
            and row["MaximumStatisticUpdateResidual"] <= small_gates["maximum_statistic_update_residual_max"]
        )
        small_rows.append(row)

    large = plan["large_design"]
    large_gates = plan["registered_gates"]["large_design"]
    locked = plan["locked_dependencies"]
    base = json.loads((ROOT / locked["base_plan"]).read_text(encoding="utf-8"))
    coordinates = pd.read_csv(ROOT / locked["person_coordinates"])
    person_coordinates = dict(coordinates[["Person", "Theta"]].itertuples(index=False, name=None))
    people = sorted(person_coordinates)
    raters = [f"R{index:02d}" for index in range(1, 5)]
    rater_coordinates = dict(zip(raters, map(float, base["design"]["rater_coordinates"]), strict=True))
    person_degrees = {person: int(large["person_degree"]) for person in people}
    rater_degrees = {rater: int(large["rater_degree"]) for rater in raters}
    retained = pd.read_csv(ROOT / locked["retained_recursive_diagnostics"])
    retained_logz = {float(row.Gamma): float(row.DPLogNormalizerUnconditioned) for row in retained.itertuples()}
    registered_logz = {float(key): float(value) for key, value in large["retained_recursive_log_normalizers"].items()}
    if any(not np.isclose(retained_logz[key], value, atol=0.0, rtol=0.0) for key, value in registered_logz.items()):
        raise ValueError("Registered recursive log normalizers do not match the retained artifact.")

    large_rows: list[dict[str, object]] = []
    for gamma_index, gamma in enumerate(large["gammas"]):
        gamma = float(gamma)
        dense = sample_degree_conditioned_assignment_dense_dp(
            person_degrees=person_degrees,
            rater_degrees=rater_degrees,
            person_coordinates=person_coordinates,
            rater_coordinates=rater_coordinates,
            gamma=gamma,
            n_samples=int(large["independent_draws_per_gamma"]),
            seed=int(large["seed_base"]) + gamma_index,
            require_connected=True,
            max_dense_dp_cells=int(large["max_dense_dp_cells"]),
            retain_edge_signatures=False,
        )
        diagnostic = dense["diagnostics"].iloc[0]
        dense_logz = float(diagnostic["DPLogNormalizerUnconditioned"])
        partition_seconds = float(diagnostic["DenseDPPartitionSeconds"])
        row = {
            "Gamma": gamma,
            "Available": bool(dense["available"]),
            "RecursiveLogNormalizer": retained_logz[gamma],
            "DenseLogNormalizer": dense_logz,
            "AbsoluteLogNormalizerDifference": abs(retained_logz[gamma] - dense_logz),
            "DenseDPShape": str(diagnostic["DenseDPShape"]),
            "DenseDPCells": int(diagnostic["DenseDPCells"]),
            "DenseDPMemoryBytes": int(diagnostic["DenseDPMemoryBytes"]),
            "DenseDPPartitionSeconds": partition_seconds,
            "MaximumConditionalLogResidual": float(diagnostic["MaximumConditionalLogResidual"]),
            "MaximumStatisticUpdateResidual": float(diagnostic["MaximumStatisticUpdateResidual"]),
            "ConnectedRejections": int(diagnostic["ConnectedRejections"]),
            "FinalOverlapComponents": int(diagnostic["FinalOverlapComponents"]),
            "DegreeMarginsPreserved": bool(diagnostic["DegreeMarginsPreserved"] and _degree_audit(dense["final_edges"], person_degrees, rater_degrees)),
            "SpeedupVersusRetainedRecursiveAtPositive0p8": (
                float(large["retained_recursive_gamma_positive_0p8_partition_seconds"]) / partition_seconds
                if gamma == 0.8 else np.nan
            ),
        }
        row["Passed"] = bool(
            row["Available"] is large_gates["all_runs_available"]
            and row["AbsoluteLogNormalizerDifference"] <= large_gates["recursive_dense_log_normalizer_absolute_difference_max"]
            and row["MaximumConditionalLogResidual"] <= large_gates["maximum_conditional_log_residual_max"]
            and row["MaximumStatisticUpdateResidual"] <= large_gates["maximum_statistic_update_residual_max"]
            and row["ConnectedRejections"] <= large_gates["connected_rejections_max"]
            and row["FinalOverlapComponents"] == large_gates["final_overlap_components"]
            and row["DegreeMarginsPreserved"] is large_gates["degree_margins_preserved"]
            and row["DenseDPPartitionSeconds"] <= large_gates["maximum_partition_seconds_each_gamma"]
            and row["DenseDPShape"] == large["expected_dense_shape"]
            and row["DenseDPCells"] == large["expected_dense_cells"]
            and row["DenseDPMemoryBytes"] == large["expected_dense_memory_bytes"]
            and (
                gamma != 0.8
                or row["SpeedupVersusRetainedRecursiveAtPositive0p8"] >= large_gates["minimum_speedup_at_gamma_positive_0p8"]
            )
        )
        large_rows.append(row)

    small_table = pd.DataFrame(small_rows)
    large_table = pd.DataFrame(large_rows)
    gates = {
        "locked_dependency_hashes_match": True,
        "all_small_oracle_rows_pass": bool(small_table["Passed"].all()),
        "all_large_design_rows_pass": bool(large_table["Passed"].all()),
        "large_log_normalizer_symmetry": bool(
            np.isclose(large_table.loc[large_table["Gamma"].eq(-0.8), "DenseLogNormalizer"].iloc[0], large_table.loc[large_table["Gamma"].eq(0.8), "DenseLogNormalizer"].iloc[0], atol=1e-12)
            and np.isclose(large_table.loc[large_table["Gamma"].eq(-0.4), "DenseLogNormalizer"].iloc[0], large_table.loc[large_table["Gamma"].eq(0.4), "DenseLogNormalizer"].iloc[0], atol=1e-12)
        ),
    }
    assessment = {
        "schema_version": "known_assignment_dense_dp_assessment_v1",
        "qualification_pass": bool(all(gates.values())),
        "gates": gates,
        "maximum_small_log_normalizer_difference": float(small_table["AbsoluteLogNormalizerDifference"].max()),
        "maximum_small_total_variation_distance": float(small_table["TotalVariationDistance"].max()),
        "maximum_large_log_normalizer_difference": float(large_table["AbsoluteLogNormalizerDifference"].max()),
        "maximum_large_partition_seconds": float(large_table["DenseDPPartitionSeconds"].max()),
        "positive_0p8_speedup": float(large_table.loc[large_table["Gamma"].eq(0.8), "SpeedupVersusRetainedRecursiveAtPositive0p8"].iloc[0]),
        "claim_boundary": plan["claim_boundary"],
    }

    OUTPUT_DIR.mkdir(parents=False, exist_ok=False)
    small_table.to_csv(OUTPUT_DIR / "small_oracle_equivalence.csv", index=False)
    large_table.to_csv(OUTPUT_DIR / "large_design_equivalence.csv", index=False)
    (OUTPUT_DIR / "assessment.json").write_text(json.dumps(assessment, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    artifacts = ["small_oracle_equivalence.csv", "large_design_equivalence.csv", "assessment.json"]
    dependencies = {
        "plan": PLAN_PATH,
        "runner": Path(__file__).resolve(),
        "mechanism": ROOT / "mfrm_app" / "assignment_mechanism.py",
        "retained_recursive_diagnostics": ROOT / locked["retained_recursive_diagnostics"],
        "person_coordinates": ROOT / locked["person_coordinates"],
        "base_plan": ROOT / locked["base_plan"],
    }
    manifest = {
        "schema_version": "known_assignment_dense_dp_manifest_v1",
        "mechanism_schema_version": MECHANISM_SCHEMA_VERSION,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "dependency_sha256": {name: _sha256(path) for name, path in dependencies.items()},
        "artifact_sha256": {name: _sha256(OUTPUT_DIR / name) for name in artifacts},
    }
    (OUTPUT_DIR / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(assessment, ensure_ascii=False))


if __name__ == "__main__":
    main()
