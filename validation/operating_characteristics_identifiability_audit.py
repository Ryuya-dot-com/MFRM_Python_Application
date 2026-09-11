#!/usr/bin/env python3
"""Audit exact-free-coordinate JMLE structural identifiability on 160 RunIds."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import streamlit_app as app  # noqa: E402
from validation import operating_characteristics_strict_jmle as stage_a  # noqa: E402
from validation import operating_characteristics_jmle_movement_audit as movement  # noqa: E402


PLAN_SCHEMA = "mfrm-jmle-structural-identifiability-plan-v1"
DEFAULT_INPUT = REPO / "validation" / "operating_characteristics_pilot20_20260809"
DEFAULT_MOVEMENT = REPO / "validation" / "operating_characteristics_jmle_movement_20260809"
DEFAULT_PLAN = REPO / "validation" / "operating_characteristics_identifiability_plan_20260809.json"
DEFAULT_OUTPUT = REPO / "validation" / "operating_characteristics_identifiability_20260809"


def load_plan(path: Path) -> dict:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != PLAN_SCHEMA:
        raise ValueError(f"identifiability plan schema must be {PLAN_SCHEMA}")
    audit = plan.get("audit", {})
    if int(audit.get("attempted_runids", 0)) != 160:
        raise ValueError("identifiability audit must retain all 160 RunIds")
    if audit.get("linear_predictor_blocks") != ["theta", "Rater", "Task", "Criterion"]:
        raise ValueError("identifiability block order changed")
    if audit.get("rank_method") != "singular value decomposition":
        raise ValueError("identifiability rank method changed")
    boundaries = plan.get("scope_boundaries", {})
    if boundaries.get("public_surface_status") != "Withheld":
        raise ValueError("identifiability public surface must remain Withheld")
    if bool(boundaries.get("application_core_change_authorized", True)):
        raise ValueError("identifiability audit cannot authorize a core change")
    return plan


def validate_input_identity(input_dir: Path, movement_dir: Path, plan: dict) -> pd.DataFrame:
    expected = plan["input_identity"]
    paths = {
        "streamlit_app_sha256": REPO / "streamlit_app.py",
        "movement_plan_sha256": REPO / "validation" / "operating_characteristics_jmle_movement_plan_20260809.json",
        "movement_adapter_sha256": Path(movement.__file__).resolve(),
        "movement_run_summary_sha256": movement_dir / "jmle_movement_run_summary.csv",
        "movement_identity_sha256": movement_dir / "jmle_movement_identity.json",
        "pilot_manifest_file_sha256": input_dir / "manifest.csv",
        "pilot_generated_ratings_sha256": input_dir / "generated_ratings.csv",
        "pilot_generated_anchors_sha256": input_dir / "generated_anchors.csv",
    }
    rows = []
    for key, path in paths.items():
        actual = stage_a.sha256_file(path)
        target = str(expected[key])
        rows.append({
            "Check": key,
            "Passed": actual == target,
            "Evidence": f"actual={actual}; expected={target}",
        })
    identity = json.loads((movement_dir / "jmle_movement_identity.json").read_text(encoding="utf-8"))
    rows.append({
        "Check": "movement_audit_scope",
        "Passed": bool(
            identity.get("application_core_change_authorized") is False
            and identity.get("public_surface_enabled") is False
        ),
        "Evidence": (
            f"core={identity.get('application_core_change_authorized')}; "
            f"public={identity.get('public_surface_enabled')}; disposition={identity.get('disposition')}"
        ),
    })
    checks = pd.DataFrame(rows)
    if not checks["Passed"].all():
        failed = checks.loc[~checks["Passed"], "Check"].astype(str).tolist()
        raise ValueError(f"identifiability input identity failed: {failed}")
    return checks


def expansion_matrix(spec: dict) -> np.ndarray:
    size = int(spec["n_params"])
    baseline = app.expand_facet_with_constraints(np.zeros(size, dtype=float), spec)
    matrix = np.zeros((len(spec["levels"]), size), dtype=float)
    for coordinate in range(size):
        basis = np.zeros(size, dtype=float)
        basis[coordinate] = 1.0
        matrix[:, coordinate] = app.expand_facet_with_constraints(basis, spec) - baseline
    return matrix


def prepare_design(generated, manifest_row: pd.Series) -> tuple[dict, dict, dict, dict]:
    prep = app.prepare_mfrm_data(
        generated.data,
        person_col="Person",
        facet_cols=["Rater", "Task", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=int(manifest_row["Categories"]) - 1,
        keep_original=True,
    )
    specs = app.prepare_constraint_specs(
        prep,
        anchor_df=generated.anchors if not generated.anchors.empty else None,
        noncenter_facet="Person",
    )
    config = {
        "model": "RSM",
        "method": "JMLE",
        "n_cat": int(manifest_row["Categories"]),
        "n_person": len(prep["levels"]["Person"]),
        "facet_names": prep["facet_names"],
        "facet_levels": {facet: prep["levels"][facet] for facet in prep["facet_names"]},
        "facet_signs": {"Rater": -1, "Task": -1, "Criterion": -1},
        "theta_spec": specs["theta_spec"],
        "facet_specs": specs["facet_specs"],
        "population_model": {"enabled": False, "n_params": 0},
    }
    sizes = app.build_param_sizes(config)
    indices = app.build_indices(prep)
    return prep, config, sizes, indices


def eta_design_matrix(config: dict, sizes: dict, indices: dict) -> tuple[np.ndarray, pd.DataFrame]:
    pieces: list[np.ndarray] = []
    metadata: list[dict[str, object]] = []
    global_index = 0
    theta_matrix = expansion_matrix(config["theta_spec"])
    if theta_matrix.shape[1]:
        pieces.append(theta_matrix[indices["person"]])
        labels = movement.constrained_free_labels(config["theta_spec"])
        for local_index, label in enumerate(labels):
            metadata.append({
                "GlobalCoordinate": global_index,
                "Block": "theta",
                "BlockCoordinate": local_index,
                **label,
            })
            global_index += 1
    for facet in config["facet_names"]:
        spec = config["facet_specs"][facet]
        matrix = expansion_matrix(spec)
        if not matrix.shape[1]:
            continue
        sign = float(config["facet_signs"][facet])
        pieces.append(sign * matrix[indices["facets"][facet]])
        labels = movement.constrained_free_labels(spec)
        for local_index, label in enumerate(labels):
            metadata.append({
                "GlobalCoordinate": global_index,
                "Block": facet,
                "BlockCoordinate": local_index,
                **label,
            })
            global_index += 1
    design = np.column_stack(pieces) if pieces else np.zeros((len(indices["person"]), 0), dtype=float)
    metadata_frame = pd.DataFrame(metadata)
    if design.shape[1] != len(metadata_frame):
        raise RuntimeError("eta design columns do not match exact-free-coordinate metadata")
    return design, metadata_frame


def graph_components(data: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    adjacency: dict[str, set[str]] = {}
    for row in data[["Person", "Rater"]].drop_duplicates().itertuples(index=False):
        person = f"P::{row.Person}"
        rater = f"R::{row.Rater}"
        adjacency.setdefault(person, set()).add(rater)
        adjacency.setdefault(rater, set()).add(person)
    seen: set[str] = set()
    rows = []
    component_id = 0
    for node in sorted(adjacency):
        if node in seen:
            continue
        component_id += 1
        stack = [node]
        component: set[str] = set()
        while stack:
            current = stack.pop()
            if current in component:
                continue
            component.add(current)
            seen.add(current)
            stack.extend(adjacency[current] - component)
        persons = sorted(value[3:] for value in component if value.startswith("P::"))
        raters = sorted(value[3:] for value in component if value.startswith("R::"))
        rows.append({
            "Component": component_id,
            "Persons": len(persons),
            "Raters": len(raters),
            "PersonLevels": "|".join(persons),
            "RaterLevels": "|".join(raters),
        })
    return pd.DataFrame(rows), component_id


def audit_run(
    manifest_row: pd.Series,
    generated,
    movement_row: pd.Series,
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prep, config, sizes, indices = prepare_design(generated, manifest_row)
    design, metadata = eta_design_matrix(config, sizes, indices)
    _, singular_values, vh = np.linalg.svd(design, full_matrices=True)
    largest = float(singular_values[0]) if len(singular_values) else 0.0
    tolerance = max(design.shape) * np.finfo(float).eps * largest
    rank = int(np.sum(singular_values > tolerance))
    nullity = int(design.shape[1] - rank)
    smallest_nonzero = float(singular_values[rank - 1]) if rank > 0 else np.nan
    smallest = float(singular_values[-1]) if len(singular_values) else np.nan
    condition_number = largest / smallest_nonzero if rank > 0 and smallest_nonzero > 0 else np.inf
    null_vectors = vh[rank:, :] if nullity > 0 else np.zeros((0, design.shape[1]), dtype=float)
    null_weight = np.sum(null_vectors ** 2, axis=0) if nullity > 0 else np.zeros(design.shape[1], dtype=float)
    coordinate_null = metadata.copy()
    coordinate_null["NullProjectionWeight"] = null_weight
    coordinate_null["NullVulnerable"] = coordinate_null["NullProjectionWeight"].gt(1e-10)

    energy_rows: list[dict[str, object]] = []
    for vector_index, vector in enumerate(null_vectors):
        denominator = float(np.sum(vector ** 2))
        for block, block_frame in metadata.groupby("Block", sort=False):
            coordinates = block_frame["GlobalCoordinate"].to_numpy(dtype=int)
            energy_rows.append({
                "NullVector": vector_index + 1,
                "Block": block,
                "SquaredEnergyShare": float(np.sum(vector[coordinates] ** 2) / denominator) if denominator > 0 else np.nan,
                "MaximumAbsLoading": float(np.max(np.abs(vector[coordinates]))) if len(coordinates) else 0.0,
            })
    null_energy = pd.DataFrame(energy_rows)
    components, component_count = graph_components(prep["data"])
    identity = {
        "RunId": str(manifest_row["RunId"]),
        "ConditionId": str(manifest_row["ConditionId"]),
        "Design": str(manifest_row["Design"]),
        "TruthBias": float(manifest_row["TruthBias"]),
        "Replicate": int(manifest_row["Replicate"]),
    }
    for frame in (coordinate_null, null_energy, components):
        if frame.empty:
            continue
        for column, value in reversed(identity.items()):
            frame.insert(0, column, value)
    raters_per_person = prep["data"].groupby("Person", observed=True)["Rater"].nunique()
    persons_per_rater = prep["data"].groupby("Rater", observed=True)["Person"].nunique()
    row = {
        **identity,
        "Rows": int(design.shape[0]),
        "EtaFreeCoordinates": int(design.shape[1]),
        "EtaDesignRank": rank,
        "StructuralNullity": nullity,
        "RankTolerance": tolerance,
        "LargestSingularValue": largest,
        "SmallestSingularValue": smallest,
        "SmallestNonzeroSingularValue": smallest_nonzero,
        "NonzeroConditionNumber": condition_number,
        "PersonRaterComponents": component_count,
        "MinRatersPerPerson": int(raters_per_person.min()),
        "MedianRatersPerPerson": float(raters_per_person.median()),
        "MaxRatersPerPerson": int(raters_per_person.max()),
        "MinPersonsPerRater": int(persons_per_rater.min()),
        "MaxPersonsPerRater": int(persons_per_rater.max()),
        "FixedRaterAnchors": int(len(generated.anchors.loc[generated.anchors["Facet"].astype(str).eq("Rater")])),
        "StructurallyIdentified": nullity == 0,
        "MovementLocalizationClass": str(movement_row["LocalizationClass"]),
        "LargeNonThetaCoordinates": int(movement_row["LargeNonThetaCoordinates"]),
        "MaximumAbsFreeCoordinateMovement": float(movement_row["MaximumAbsFreeCoordinateMovement"]),
    }
    return row, coordinate_null, null_energy, components


def condition_summary(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, frame in runs.groupby(["ConditionId", "Design", "TruthBias"], sort=False):
        rows.append({
            "ConditionId": keys[0],
            "Design": keys[1],
            "TruthBias": keys[2],
            "Runs": len(frame),
            "StructurallyIdentifiedRuns": int(frame["StructurallyIdentified"].astype(bool).sum()),
            "StructuralNullityMin": int(frame["StructuralNullity"].min()),
            "StructuralNullityMax": int(frame["StructuralNullity"].max()),
            "PersonRaterComponentsMin": int(frame["PersonRaterComponents"].min()),
            "PersonRaterComponentsMax": int(frame["PersonRaterComponents"].max()),
            "RatersPerPersonMin": int(frame["MinRatersPerPerson"].min()),
            "RatersPerPersonMax": int(frame["MaxRatersPerPerson"].max()),
            "FixedRaterAnchors": int(frame["FixedRaterAnchors"].max()),
            "RunsWithLargeNonThetaMovement": int(frame["LargeNonThetaCoordinates"].gt(0).sum()),
            "MaximumAbsFreeCoordinateMovement": float(frame["MaximumAbsFreeCoordinateMovement"].max()),
            "IntegrationDisposition": (
                "eligible_for_separate_downstream_reconstruction_experiment"
                if frame["StructurallyIdentified"].all()
                else "blocked_structural_rank_deficiency"
            ),
        })
    return pd.DataFrame(rows)


def audit_gates(runs: pd.DataFrame, input_checks: pd.DataFrame, plan: dict) -> pd.DataFrame:
    expected = int(plan["audit"]["attempted_runids"])
    definitions = [
        ("input_identity", len(input_checks), int(input_checks["Passed"].sum())),
        ("runids_audited", expected, len(runs)),
        ("finite_rank_results", expected, int(runs["EtaDesignRank"].notna().sum())),
        ("graph_components_recorded", expected, int(runs["PersonRaterComponents"].ge(1).sum())),
    ]
    return pd.DataFrame([
        {"Gate": gate, "Required": required, "Passed": passed, "GatePassed": required == passed}
        for gate, required, passed in definitions
    ])


def decision_table(runs: pd.DataFrame, gates: pd.DataFrame) -> pd.DataFrame:
    complete = bool(gates["GatePassed"].all())
    deficient = int(runs["StructuralNullity"].gt(0).sum())
    affected_conditions = int(runs.loc[runs["StructuralNullity"].gt(0), "ConditionId"].nunique())
    disposition = (
        "blocked_incomplete_audit" if not complete else
        "blocked_unconditional_jmle_integration" if deficient else
        "eligible_for_downstream_reconstruction_experiment"
    )
    return pd.DataFrame([
        {
            "Decision": "StructuralIdentifiability",
            "Disposition": disposition,
            "Evidence": f"complete={complete}; rank_deficient_runs={deficient}; affected_conditions={affected_conditions}",
            "NextAction": (
                "Add a pre-fit rank/connectivity gate and redesign or reroute structurally deficient designs before downstream reconstruction."
                if deficient else
                "Proceed to a separate downstream reconstruction experiment."
            ),
            "ApplicationCoreChangeAuthorized": False,
            "PublicSurfaceEnabled": False,
        }
    ])


def plot_condition_nullity(conditions: pd.DataFrame, output: Path) -> None:
    labels = conditions["ConditionId"].str.replace("__", "\n", regex=False)
    values = conditions["StructuralNullityMax"]
    colors = np.where(values.gt(0), "#E15759", "#59A14F")
    fig, ax = plt.subplots(figsize=(12.5, 6.2))
    ax.bar(np.arange(len(conditions)), values, color=colors)
    ax.set_xticks(np.arange(len(conditions)), labels, rotation=35, ha="right")
    ax.set_ylabel("Exact-free-coordinate structural nullity")
    ax.set_title("JMLE linear-predictor structural identifiability by frozen condition")
    ax.grid(axis="y", alpha=0.22)
    fig.tight_layout()
    fig.savefig(output / "jmle_identifiability_nullity_by_condition.png", dpi=180)
    plt.close(fig)


def plot_nullity_movement(runs: pd.DataFrame, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    for design, frame in runs.groupby("Design", sort=False):
        jitter = np.linspace(-0.09, 0.09, len(frame))
        ax.scatter(
            frame["StructuralNullity"].to_numpy(dtype=float) + jitter,
            frame["MaximumAbsFreeCoordinateMovement"],
            alpha=0.65,
            s=28,
            label=str(design),
        )
    ax.set_yscale("log")
    ax.set_xlabel("Structural nullity (small horizontal jitter for visibility)")
    ax.set_ylabel("Maximum absolute Stage-B2 free-coordinate movement")
    ax.set_title("Numerical movement versus exact structural nullity")
    ax.grid(axis="y", which="both", alpha=0.22)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "jmle_identifiability_nullity_vs_movement.png", dpi=180)
    plt.close(fig)


def write_report(
    output: Path,
    runs: pd.DataFrame,
    conditions: pd.DataFrame,
    null_energy: pd.DataFrame,
    gates: pd.DataFrame,
    decision: pd.DataFrame,
) -> None:
    deficient = runs.loc[runs["StructuralNullity"].gt(0)]
    condition_lines = "\n".join(
        f"- `{row.ConditionId}`: nullity {int(row.StructuralNullityMin)}-{int(row.StructuralNullityMax)}, "
        f"person-rater components {int(row.PersonRaterComponentsMin)}-{int(row.PersonRaterComponentsMax)}, "
        f"disposition `{row.IntegrationDisposition}`."
        for row in conditions.itertuples(index=False)
    )
    if null_energy.empty:
        energy_text = "No null-space vectors were present."
    else:
        energy = (
            null_energy.groupby("Block", sort=False)["SquaredEnergyShare"]
            .agg(["mean", "max"])
            .reset_index()
        )
        energy_text = "\n".join(
            f"- `{row.Block}`: mean share {row.mean:.4f}, maximum share {row.max:.4f}."
            for row in energy.itertuples(index=False)
        )
    report = f"""# JMLE structural identifiability audit

## Decision

Disposition: `{decision.iloc[0]['Disposition']}`. The audit covered
{len(runs)}/160 RunIds; {len(deficient)} were structurally rank deficient.
No application-core change or public surface is authorized.

## Condition results

{condition_lines}

## Null-space localization

{energy_text}

## Interpretation

The audit uses the exact free-coordinate expansion already used by the Python
JMLE core. A positive nullity means at least one parameter direction leaves the
linear predictor unchanged exactly. In that case a small gradient and an
optimizer `success=True` cannot establish a unique estimate. Person-rater graph
components provide a design-level explanation, while the SVD is the decisive
algebraic check after anchors and centering constraints are applied.

For a one-rater-per-person design, person effects are nested within raters.
Without enough cross-rater overlap or identifying anchors, joint person/rater
location shifts remain possible. The remedy is a design/policy choice—add
connectivity, anchor appropriately, fit a model that explicitly handles the
nested/random structure, or withhold JMLE—not a looser floating-point threshold.

## Retained artifacts

- `jmle_identifiability_runs.csv`
- `jmle_identifiability_conditions.csv`
- `jmle_identifiability_graph_components.csv`
- `jmle_identifiability_coordinate_null_weight.csv`
- `jmle_identifiability_null_space_block_energy.csv`
- `jmle_identifiability_audit_gates.csv`
- `jmle_identifiability_decision.csv`
- `jmle_identifiability_input_checks.csv`
- `jmle_identifiability_first_read_summary.csv`
- `jmle_identifiability_identity.json`
- `jmle_identifiability_nullity_by_condition.png`
- `jmle_identifiability_nullity_vs_movement.png`
"""
    (output / "JMLE_IDENTIFIABILITY_AUDIT.md").write_text(report, encoding="utf-8")


def run(input_dir: Path, movement_dir: Path, output: Path, plan_path: Path) -> None:
    input_dir = input_dir.resolve()
    movement_dir = movement_dir.resolve()
    output = output.resolve()
    plan_path = plan_path.resolve()
    plan = load_plan(plan_path)
    input_checks = validate_input_identity(input_dir, movement_dir, plan)
    manifest = pd.read_csv(input_dir / "manifest.csv")
    ratings = pd.read_csv(input_dir / "generated_ratings.csv")
    truth = pd.read_csv(input_dir / "generated_facet_truth.csv")
    anchors = pd.read_csv(input_dir / "generated_anchors.csv")
    movement_summary = pd.read_csv(movement_dir / "jmle_movement_run_summary.csv").set_index("RunId", drop=False)
    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_jmle_identifiability_plan.json").write_bytes(plan_path.read_bytes())
    run_rows: list[dict[str, object]] = []
    coordinate_parts: list[pd.DataFrame] = []
    energy_parts: list[pd.DataFrame] = []
    component_parts: list[pd.DataFrame] = []
    started = time.perf_counter()
    for _, manifest_row in manifest.iterrows():
        run_id = str(manifest_row["RunId"])
        generated = stage_a.generated_design_for_run(manifest_row, ratings, truth, anchors)
        row, coordinates, energy, components = audit_run(
            manifest_row,
            generated,
            movement_summary.loc[run_id],
        )
        run_rows.append(row)
        coordinate_parts.append(coordinates)
        if not energy.empty:
            energy_parts.append(energy)
        component_parts.append(components)
    runs = pd.DataFrame(run_rows)
    coordinates = pd.concat(coordinate_parts, ignore_index=True)
    energy = pd.concat(energy_parts, ignore_index=True) if energy_parts else pd.DataFrame()
    components = pd.concat(component_parts, ignore_index=True)
    conditions = condition_summary(runs)
    gates = audit_gates(runs, input_checks, plan)
    decision = decision_table(runs, gates)
    deficient_conditions = conditions.loc[conditions["StructuralNullityMax"].gt(0), "ConditionId"].tolist()
    first_read = pd.DataFrame([
        {
            "Priority": 1,
            "Check": "Structural identifiability",
            "Status": str(decision.iloc[0]["Disposition"]),
            "Evidence": f"rank-deficient RunIds={int(runs['StructuralNullity'].gt(0).sum())}; conditions={deficient_conditions}",
            "NextAction": str(decision.iloc[0]["NextAction"]),
            "PublicSurfaceEnabled": False,
        },
        {
            "Priority": 2,
            "Check": "Numerical convergence",
            "Status": "Necessary, not sufficient",
            "Evidence": "Stage-B2 gradients passed even where the exact linear-predictor design was rank deficient.",
            "NextAction": "Place structural checks before optimizer and fit-threshold interpretation.",
            "PublicSurfaceEnabled": False,
        },
        {
            "Priority": 3,
            "Check": "Application core",
            "Status": "Withheld",
            "Evidence": "No unconditional JMLE integration is authorized.",
            "NextAction": "Design an automatic pre-fit route for connected, anchored, and structurally deficient inputs.",
            "PublicSurfaceEnabled": False,
        },
    ])
    runs.to_csv(output / "jmle_identifiability_runs.csv", index=False)
    conditions.to_csv(output / "jmle_identifiability_conditions.csv", index=False)
    components.to_csv(output / "jmle_identifiability_graph_components.csv", index=False)
    coordinates.to_csv(output / "jmle_identifiability_coordinate_null_weight.csv", index=False)
    energy.to_csv(output / "jmle_identifiability_null_space_block_energy.csv", index=False)
    gates.to_csv(output / "jmle_identifiability_audit_gates.csv", index=False)
    decision.to_csv(output / "jmle_identifiability_decision.csv", index=False)
    input_checks.to_csv(output / "jmle_identifiability_input_checks.csv", index=False)
    first_read.to_csv(output / "jmle_identifiability_first_read_summary.csv", index=False)
    plot_condition_nullity(conditions, output)
    plot_nullity_movement(runs, output)
    write_report(output, runs, conditions, energy, gates, decision)
    artifact_names = [
        "jmle_identifiability_runs.csv",
        "jmle_identifiability_conditions.csv",
        "jmle_identifiability_graph_components.csv",
        "jmle_identifiability_coordinate_null_weight.csv",
        "jmle_identifiability_null_space_block_energy.csv",
        "jmle_identifiability_audit_gates.csv",
        "jmle_identifiability_decision.csv",
        "jmle_identifiability_input_checks.csv",
        "jmle_identifiability_first_read_summary.csv",
        "jmle_identifiability_nullity_by_condition.png",
        "jmle_identifiability_nullity_vs_movement.png",
        "JMLE_IDENTIFIABILITY_AUDIT.md",
    ]
    identity = {
        "schema_version": PLAN_SCHEMA,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "wall_seconds": time.perf_counter() - started,
        "disposition": str(decision.iloc[0]["Disposition"]),
        "rank_deficient_runids": int(runs["StructuralNullity"].gt(0).sum()),
        "application_core_change_authorized": False,
        "public_surface_enabled": False,
        "identifiability_plan_sha256": stage_a.sha256_file(plan_path),
        "identifiability_adapter_sha256": stage_a.sha256_file(Path(__file__).resolve()),
        "streamlit_app_sha256": stage_a.sha256_file(REPO / "streamlit_app.py"),
        "artifacts": {name: stage_a.sha256_file(output / name) for name in artifact_names},
    }
    (output / "jmle_identifiability_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--movement", type=Path, default=DEFAULT_MOVEMENT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    args = parser.parse_args()
    run(args.input, args.movement, args.output, args.plan)


if __name__ == "__main__":
    main()
