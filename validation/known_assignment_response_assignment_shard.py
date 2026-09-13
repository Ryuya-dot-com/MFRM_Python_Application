"""Generate one frozen gamma shard of exact assignment graphs for response screening."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pandas as pd

from mfrm_app.assignment_mechanism import (
    materialize_exchangeable_assignment_design,
    sample_degree_conditioned_assignment_dp,
)
from validation.known_assignment_large_design_calibration import (
    _exchangeable_rows,
    _planned_edges,
)
from validation.known_assignment_large_design_dp_shard import _sha256, _slug


PLAN_PATH = ROOT / "validation" / "known_assignment_response_screening_plan_20260811.json"
OUTPUT_ROOT = ROOT / "validation" / "known_assignment_response_assignment_shards_20260811"


def _parse_edges(signature: str) -> set[tuple[str, str]]:
    return {tuple(value.split("::", maxsplit=1)) for value in signature.split("|")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gamma", required=True, type=float)
    args = parser.parse_args()
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    gammas = [float(value) for value in plan["assignment"]["gammas"]]
    if args.gamma not in gammas:
        raise ValueError(f"gamma must be one of {gammas}")
    output = OUTPUT_ROOT / _slug(args.gamma)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite assignment shard: {output}")
    assessment_path = ROOT / plan["assignment"]["qualification_assessment"]
    if _sha256(assessment_path).lower() != plan["assignment"]["qualification_assessment_sha256"]:
        raise ValueError("Assignment qualification assessment hash mismatch")
    coordinate_path = ROOT / plan["data_generating_process"]["fixed_person_coordinate_source"]
    if _sha256(coordinate_path).lower() != plan["data_generating_process"]["fixed_person_coordinate_sha256"]:
        raise ValueError("Fixed Person coordinate hash mismatch")
    coordinates = pd.read_csv(coordinate_path)
    person_coordinates = dict(coordinates[["Person", "Theta"]].itertuples(index=False, name=None))
    people = sorted(person_coordinates)
    rater_coordinates = {
        str(key): float(value)
        for key, value in plan["data_generating_process"]["raters"].items()
    }
    sampled = sample_degree_conditioned_assignment_dp(
        person_degrees={person: int(plan["assignment"]["person_degree"]) for person in people},
        rater_degrees={rater: int(plan["assignment"]["rater_degree"]) for rater in rater_coordinates},
        person_coordinates=person_coordinates,
        rater_coordinates=rater_coordinates,
        gamma=args.gamma,
        n_samples=int(plan["evidence_status"]["replicates"]),
        seed=int(plan["assignment"]["gamma_seeds"][str(args.gamma)]),
        require_connected=True,
        retain_edge_signatures=True,
    )
    source_edges = _planned_edges(people, sorted(rater_coordinates))
    source_rows = _exchangeable_rows(source_edges)
    edge_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    for sample in sampled["samples"].itertuples(index=False):
        replicate = int(sample.Sample)
        edges = _parse_edges(str(sample.EdgeSignature))
        materialized = materialize_exchangeable_assignment_design(
            source_rows,
            edges,
            facet_cols=["Rater", "Task", "Criterion"],
            context_cols=["Task", "Criterion"],
            rater_mapping_confirmed=True,
        )
        if not materialized["available"] or not materialized["invariants"]["Passed"].all():
            raise RuntimeError(f"Assignment materialization failed: gamma={args.gamma}, replicate={replicate}")
        for person, rater in sorted(edges):
            edge_rows.append({"Gamma": args.gamma, "Replicate": replicate, "Person": person, "Rater": rater})
        for row in materialized["invariants"].itertuples(index=False):
            audit_rows.append(
                {
                    "Gamma": args.gamma,
                    "Replicate": replicate,
                    "Invariant": row.Invariant,
                    "Passed": bool(row.Passed),
                    "Evidence": row.Evidence,
                }
            )
    output.mkdir(parents=True, exist_ok=False)
    artifacts = {
        "assignment_edges.csv": pd.DataFrame(edge_rows),
        "assignment_draws.csv": sampled["samples"],
        "assignment_diagnostics.csv": sampled["diagnostics"],
        "materialization_audit.csv": pd.DataFrame(audit_rows),
    }
    for name, table in artifacts.items():
        table.to_csv(output / name, index=False)
    dependencies = {
        "plan": PLAN_PATH,
        "runner": Path(__file__).resolve(),
        "mechanism": ROOT / "mfrm_app" / "assignment_mechanism.py",
    }
    manifest = {
        "schema_version": "known_assignment_response_assignment_shard_manifest_v1",
        "gamma": args.gamma,
        "dependency_sha256": {name: _sha256(path) for name, path in dependencies.items()},
        "artifact_sha256": {name: _sha256(output / name) for name in artifacts},
    }
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"gamma": args.gamma, "graphs": int(sampled["samples"].shape[0]), "diagnostics": sampled["diagnostics"].iloc[0].to_dict()}, default=str))


if __name__ == "__main__":
    main()
