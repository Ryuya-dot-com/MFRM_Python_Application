"""Frozen small-state qualification of the known assignment mechanism."""

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
    materialize_exchangeable_assignment_design,
    sample_degree_conditioned_assignment_chain,
)


PLAN_PATH = ROOT / "validation" / "known_assignment_mechanism_pilot_plan_20260811.json"
OUTPUT_DIR = ROOT / "validation" / "known_assignment_mechanism_pilot_20260811"

PEOPLE = [f"p{index}" for index in range(6)]
RATERS = [f"r{index}" for index in range(3)]
PERSON_DEGREES = {person: 2 for person in PEOPLE}
RATER_DEGREES = {rater: 4 for rater in RATERS}
PERSON_COORDINATES = {person: float(index) for index, person in enumerate(PEOPLE)}
RATER_COORDINATES = {rater: float(index) for index, rater in enumerate(RATERS)}
INITIAL_EDGES = {
    ("p0", "r0"), ("p0", "r1"),
    ("p1", "r0"), ("p1", "r2"),
    ("p2", "r0"), ("p2", "r1"),
    ("p3", "r1"), ("p3", "r2"),
    ("p4", "r0"), ("p4", "r2"),
    ("p5", "r1"), ("p5", "r2"),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _exchangeable_rows() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for person, rater in sorted(INITIAL_EDGES):
        for task in ("t1", "t2"):
            for criterion in ("c1", "c2"):
                rows.append(
                    {
                        "Person": person,
                        "Rater": rater,
                        "Task": task,
                        "Criterion": criterion,
                        "Score": (int(person[1:]) + int(rater[1:])) % 4,
                        "Weight": 1.0,
                    }
                )
    return pd.DataFrame(rows)


def _gamma_label(gamma: float) -> str:
    return f"{gamma:+.1f}"


def main() -> None:
    if OUTPUT_DIR.exists():
        raise FileExistsError(
            f"Frozen output already exists; refusing to overwrite: {OUTPUT_DIR}"
        )
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    OUTPUT_DIR.mkdir(parents=False, exist_ok=False)

    oracle_by_gamma: dict[float, dict[str, object]] = {}
    oracle_summaries: list[pd.DataFrame] = []
    oracle_states: list[pd.DataFrame] = []
    chain_diagnostics: list[pd.DataFrame] = []
    chain_comparisons: list[pd.DataFrame] = []
    state_frequencies: list[pd.DataFrame] = []
    pooled_comparisons: list[pd.DataFrame] = []
    materialization_audits: list[pd.DataFrame] = []
    samples_by_gamma: dict[float, list[pd.DataFrame]] = {}

    for gamma in [float(value) for value in plan["mechanism"]["gammas"]]:
        oracle = enumerate_degree_conditioned_assignments(
            person_degrees=PERSON_DEGREES,
            rater_degrees=RATER_DEGREES,
            person_coordinates=PERSON_COORDINATES,
            rater_coordinates=RATER_COORDINATES,
            gamma=gamma,
            require_connected=True,
        )
        oracle_by_gamma[gamma] = oracle
        oracle_summary = oracle["summary"].copy()
        oracle_summary.insert(0, "GammaLabel", _gamma_label(gamma))
        oracle_summaries.append(oracle_summary)
        states = oracle["states"].copy()
        states.insert(0, "Gamma", gamma)
        oracle_states.append(states)
        samples_by_gamma[gamma] = []

        seeds = plan["chains"]["seeds"][str(gamma)]
        for chain_number, seed in enumerate(seeds, start=1):
            chain_id = f"gamma_{_gamma_label(gamma)}_chain_{chain_number}"
            chain = sample_degree_conditioned_assignment_chain(
                initial_edges=INITIAL_EDGES,
                person_coordinates=PERSON_COORDINATES,
                rater_coordinates=RATER_COORDINATES,
                gamma=gamma,
                n_samples=int(plan["chains"]["samples_per_chain"]),
                burnin=int(plan["chains"]["burnin"]),
                thin=int(plan["chains"]["thin"]),
                seed=int(seed),
                require_connected=True,
            )
            samples = chain["samples"].copy()
            samples_by_gamma[gamma].append(samples)
            diagnostics = chain["diagnostics"].copy()
            diagnostics.insert(0, "ChainId", chain_id)
            diagnostics["Available"] = bool(chain["available"])
            chain_diagnostics.append(diagnostics)
            comparison = compare_chain_to_exact_oracle(chain, oracle)["summary"].copy()
            comparison.insert(0, "ChainId", chain_id)
            comparison.insert(1, "Gamma", gamma)
            chain_comparisons.append(comparison)
            frequencies = (
                samples["EdgeSignature"]
                .value_counts()
                .rename_axis("EdgeSignature")
                .reset_index(name="Count")
            )
            frequencies.insert(0, "ChainId", chain_id)
            frequencies.insert(1, "Gamma", gamma)
            state_frequencies.append(frequencies)

            materialized = materialize_exchangeable_assignment_design(
                _exchangeable_rows(),
                chain["final_edges"],
                facet_cols=["Rater", "Task", "Criterion"],
                context_cols=["Task", "Criterion"],
                rater_mapping_confirmed=True,
            )
            audit = materialized["invariants"].copy()
            audit.insert(0, "ChainId", chain_id)
            audit["ScoreColumnAbsent"] = "Score" not in materialized["design"].columns
            materialization_audits.append(audit)

    for gamma, sample_frames in samples_by_gamma.items():
        pooled_samples = pd.concat(sample_frames, ignore_index=True)
        pooled = compare_chain_to_exact_oracle(
            {"samples": pooled_samples}, oracle_by_gamma[gamma]
        )["summary"].copy()
        pooled.insert(0, "Gamma", gamma)
        pooled_comparisons.append(pooled)

    oracle_summary_table = pd.concat(oracle_summaries, ignore_index=True)
    oracle_state_table = pd.concat(oracle_states, ignore_index=True)
    diagnostics_table = pd.concat(chain_diagnostics, ignore_index=True)
    comparison_table = pd.concat(chain_comparisons, ignore_index=True)
    frequency_table = pd.concat(state_frequencies, ignore_index=True)
    pooled_table = pd.concat(pooled_comparisons, ignore_index=True)
    materialization_table = pd.concat(materialization_audits, ignore_index=True)

    # Outcome-blind negative control: arbitrary Score changes cannot alter the
    # score-free row materialization for a fixed target graph.
    target_edges = next(iter(oracle_by_gamma.values()))["states"].iloc[-1]["EdgeSignature"]
    target = {
        tuple(edge.split("::", maxsplit=1))
        for edge in str(target_edges).split("|")
    }
    source = _exchangeable_rows()
    changed = source.copy()
    changed["Score"] = np.arange(len(changed), dtype=int) * 101
    materialize_kwargs = dict(
        target_edges=target,
        facet_cols=["Rater", "Task", "Criterion"],
        context_cols=["Task", "Criterion"],
        rater_mapping_confirmed=True,
    )
    negative_original = materialize_exchangeable_assignment_design(
        source, **materialize_kwargs
    )
    negative_changed = materialize_exchangeable_assignment_design(
        changed, **materialize_kwargs
    )
    score_negative_control = bool(
        negative_original["available"]
        and negative_changed["available"]
        and negative_original["design"].equals(negative_changed["design"])
        and negative_original["assignment_map"].equals(
            negative_changed["assignment_map"]
        )
    )

    gates = plan["registered_gates"]
    oracle_gate = gates["oracle"]
    ordered = oracle_summary_table.sort_values("Gamma")
    oracle_checks = {
        "states_each_gamma": bool(
            oracle_summary_table["States"].eq(oracle_gate["states_each_gamma"]).all()
        ),
        "probability_normalization": bool(
            (oracle_summary_table["ProbabilitySum"] - 1.0).abs().le(
                oracle_gate["absolute_probability_sum_error_max"]
            ).all()
        ),
        "gamma_zero_uniform": bool(
            (
                oracle_state_table.loc[oracle_state_table["Gamma"].eq(0.0), "Probability"]
                - 1.0 / oracle_gate["states_each_gamma"]
            ).abs().max()
            <= oracle_gate["gamma_zero_max_uniform_probability_error"]
        ),
        "expected_statistic_direction": bool(
            np.all(np.diff(ordered["ExpectedStatistic"].to_numpy(dtype=float)) > 0)
        ),
        "expected_correlation_direction": bool(
            np.all(
                np.diff(
                    ordered["ExpectedAssignmentCorrelation"].to_numpy(dtype=float)
                ) > 0
            )
        ),
    }
    chain_gate = gates["each_chain"]
    merged_chain = diagnostics_table.merge(
        comparison_table, on=["ChainId", "Gamma"], validate="one_to_one"
    )
    chain_checks = {
        "available": bool(merged_chain["Available"].all()),
        "outside_oracle": bool(
            merged_chain["ChainProbabilityOutsideOracle"].le(
                chain_gate["probability_outside_oracle_max"]
            ).all()
        ),
        "total_variation": bool(
            merged_chain["TotalVariationDistance"].le(
                chain_gate["total_variation_distance_max"]
            ).all()
        ),
        "mean_statistic": bool(
            merged_chain["MeanStatisticError"].abs().le(
                chain_gate["absolute_mean_statistic_error_max"]
            ).all()
        ),
        "ess": bool(
            merged_chain["StatisticESS"].ge(chain_gate["statistic_ess_min"]).all()
        ),
        "state_coverage": bool(
            merged_chain["UniqueSampledStates"].ge(
                chain_gate["unique_sampled_states_min"]
            ).all()
        ),
        "numerical_residual": bool(
            merged_chain["MaximumStatisticUpdateResidual"].le(
                chain_gate["maximum_statistic_update_residual_max"]
            ).all()
        ),
        "connectivity": bool(
            merged_chain["FinalOverlapComponents"].eq(
                chain_gate["final_overlap_components"]
            ).all()
        ),
        "degree_margins": bool(merged_chain["DegreeMarginsPreserved"].all()),
    }
    pooled_gate = gates["pooled_two_chains_each_gamma"]
    pooled_checks = {
        "total_variation": bool(
            pooled_table["TotalVariationDistance"].le(
                pooled_gate["total_variation_distance_max"]
            ).all()
        ),
        "mean_statistic": bool(
            pooled_table["MeanStatisticError"].abs().le(
                pooled_gate["absolute_mean_statistic_error_max"]
            ).all()
        ),
        "outside_oracle": bool(
            pooled_table["ChainProbabilityOutsideOracle"].le(
                pooled_gate["probability_outside_oracle_max"]
            ).all()
        ),
    }
    materialization_checks = {
        "all_invariants": bool(materialization_table["Passed"].all()),
        "score_absent": bool(materialization_table["ScoreColumnAbsent"].all()),
        "score_negative_control": score_negative_control,
    }
    all_checks = {
        **{f"oracle.{key}": value for key, value in oracle_checks.items()},
        **{f"chain.{key}": value for key, value in chain_checks.items()},
        **{f"pooled.{key}": value for key, value in pooled_checks.items()},
        **{
            f"materialization.{key}": value
            for key, value in materialization_checks.items()
        },
    }
    assessment = {
        "schema_version": "known_assignment_mechanism_pilot_assessment_v1",
        "mechanism_schema_version": MECHANISM_SCHEMA_VERSION,
        "decision": "pass" if all(all_checks.values()) else "fail",
        "checks": all_checks,
        "min_statistic_ess": float(merged_chain["StatisticESS"].min()),
        "max_chain_total_variation": float(
            merged_chain["TotalVariationDistance"].max()
        ),
        "max_pooled_total_variation": float(
            pooled_table["TotalVariationDistance"].max()
        ),
        "max_absolute_chain_mean_error": float(
            merged_chain["MeanStatisticError"].abs().max()
        ),
        "max_absolute_pooled_mean_error": float(
            pooled_table["MeanStatisticError"].abs().max()
        ),
        "claim_boundary": plan["claim_boundary"],
    }

    artifacts = {
        "oracle_summary.csv": oracle_summary_table,
        "oracle_states.csv": oracle_state_table,
        "chain_diagnostics.csv": diagnostics_table,
        "chain_oracle_comparison.csv": comparison_table,
        "chain_state_frequencies.csv": frequency_table,
        "pooled_oracle_comparison.csv": pooled_table,
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
    artifact_hashes = {
        name: _sha256(OUTPUT_DIR / name)
        for name in [*artifacts, assessment_path.name]
    }
    manifest = {
        "schema_version": "known_assignment_mechanism_pilot_manifest_v1",
        "plan": str(PLAN_PATH.relative_to(ROOT)),
        "plan_sha256": _sha256(PLAN_PATH),
        "runner": str(script_path.relative_to(ROOT)),
        "runner_sha256": _sha256(script_path),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "artifact_sha256": artifact_hashes,
    }
    (OUTPUT_DIR / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(assessment, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
