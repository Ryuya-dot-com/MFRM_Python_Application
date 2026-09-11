"""Execute one immutable gamma shard of the frozen exact-DP calibration."""

from __future__ import annotations

import argparse
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
)


PLAN_PATH = ROOT / "validation" / "known_assignment_large_design_dp_plan_20260811.json"
AMENDMENT_PATH = ROOT / "validation" / "known_assignment_large_design_dp_execution_amendment_20260811.json"
SHARD_ROOT = ROOT / "validation" / "known_assignment_large_design_dp_shards_20260811"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _slug(gamma: float) -> str:
    return f"gamma_{gamma:+.1f}".replace("+", "p").replace("-", "m").replace(".", "p")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gamma", required=True, type=float)
    args = parser.parse_args()
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    gammas = [float(value) for value in plan["locked_inputs"]["candidate_gammas"]]
    if args.gamma not in gammas:
        raise ValueError(f"gamma must be one of {gammas}")
    gamma_index = gammas.index(args.gamma)
    output = SHARD_ROOT / _slug(args.gamma)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite gamma shard: {output}")

    base_path = ROOT / plan["locked_inputs"]["base_plan"]
    if _sha256(base_path).lower() != plan["locked_inputs"]["base_plan_sha256"]:
        raise ValueError("Base plan hash mismatch")
    base = json.loads(base_path.read_text(encoding="utf-8"))
    coordinates = pd.read_csv(ROOT / plan["locked_inputs"]["person_coordinates"])
    person_coordinates = dict(
        coordinates[["Person", "Theta"]].itertuples(index=False, name=None)
    )
    people = sorted(person_coordinates)
    raters = [f"R{index:02d}" for index in range(1, 5)]
    rater_coordinates = dict(
        zip(raters, map(float, base["design"]["rater_coordinates"]), strict=True)
    )
    exact = plan["exact_sampling"]
    batches = int(exact["batches_per_gamma"])
    batch_size = int(exact["independent_draws_per_batch"])
    sampled = sample_degree_conditioned_assignment_dp(
        person_degrees={person: 2 for person in people},
        rater_degrees={rater: 40 for rater in raters},
        person_coordinates=person_coordinates,
        rater_coordinates=rater_coordinates,
        gamma=args.gamma,
        n_samples=batches * batch_size,
        seed=int(exact["seed_base"]) + gamma_index,
        require_connected=True,
        max_dp_states=int(exact["max_dp_states"]),
        max_total_rejections=int(exact["max_total_connected_rejections"]),
        retain_edge_signatures=bool(exact["edge_signatures_retained"]),
    )
    diagnostics = sampled["diagnostics"].copy()
    diagnostics.insert(0, "RunId", f"gamma_{args.gamma:+.1f}")
    diagnostics["Available"] = bool(sampled["available"])
    trace = sampled["samples"].copy()
    trace.insert(0, "Gamma", args.gamma)
    trace["Batch"] = (np.arange(len(trace), dtype=int) // batch_size) + 1
    trace["BatchSample"] = (np.arange(len(trace), dtype=int) % batch_size) + 1
    batch_rows: list[dict[str, object]] = []
    for batch, frame in trace.groupby("Batch", sort=True):
        values = frame["Statistic"].to_numpy(dtype=float)
        batch_rows.append(
            {
                "Gamma": args.gamma,
                "Batch": int(batch),
                "Samples": len(frame),
                "MeanStatistic": float(np.mean(values)),
                "MeanAssignmentCorrelation": float(frame["AssignmentCorrelation"].mean()),
                "StatisticESSDiagnostic": _initial_positive_sequence_ess(values),
            }
        )
    source_edges = _planned_edges(people, raters)
    materialized = materialize_exchangeable_assignment_design(
        _exchangeable_rows(source_edges),
        sampled["final_edges"],
        facet_cols=["Rater", "Task", "Criterion"],
        context_cols=["Task", "Criterion"],
        rater_mapping_confirmed=True,
    )
    audit = materialized["invariants"].copy()
    audit.insert(0, "RunId", f"gamma_{args.gamma:+.1f}")
    audit["ScoreColumnAbsent"] = "Score" not in materialized["design"].columns
    if not sampled["available"] or not audit["Passed"].all() or not audit["ScoreColumnAbsent"].all():
        raise RuntimeError("Gamma shard failed structural qualification")

    output.mkdir(parents=True, exist_ok=False)
    artifacts = {
        "run_diagnostics.csv": diagnostics,
        "independent_trace.csv": trace,
        "batch_diagnostics.csv": pd.DataFrame(batch_rows),
        "materialization_audit.csv": audit,
    }
    for name, table in artifacts.items():
        table.to_csv(output / name, index=False)
    dependencies = {
        "plan": PLAN_PATH,
        "amendment": AMENDMENT_PATH,
        "runner": Path(__file__).resolve(),
        "mechanism": ROOT / "mfrm_app" / "assignment_mechanism.py",
    }
    manifest = {
        "schema_version": "known_assignment_large_design_dp_shard_manifest_v1",
        "mechanism_schema_version": MECHANISM_SCHEMA_VERSION,
        "gamma": args.gamma,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "dependency_sha256": {name: _sha256(path) for name, path in dependencies.items()},
        "artifact_sha256": {name: _sha256(output / name) for name in artifacts},
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"gamma": args.gamma, "output": str(output), "diagnostics": diagnostics.iloc[0].to_dict()}, default=str))


if __name__ == "__main__":
    main()
