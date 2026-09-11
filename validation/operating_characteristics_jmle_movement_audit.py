#!/usr/bin/env python3
"""Localize Stage-B2 JMLE polish movement by free and expanded parameter block."""

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
from scipy.optimize import minimize


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import streamlit_app as app  # noqa: E402
from validation import operating_characteristics_strict_jmle as stage_a  # noqa: E402
from validation import operating_characteristics_strict_jmle_a2 as stage_a2  # noqa: E402
from validation import operating_characteristics_strict_jmle_b2 as stage_b2  # noqa: E402


PLAN_SCHEMA = "mfrm-jmle-movement-localization-plan-v1"
DEFAULT_INPUT = REPO / "validation" / "operating_characteristics_pilot20_20260809"
DEFAULT_B2 = REPO / "validation" / "operating_characteristics_strict_jmle_b2_20260809"
DEFAULT_STAGE_A_PLAN = REPO / "validation" / "operating_characteristics_strict_jmle_plan_20260809.json"
DEFAULT_PLAN = REPO / "validation" / "operating_characteristics_jmle_movement_plan_20260809.json"
DEFAULT_OUTPUT = REPO / "validation" / "operating_characteristics_jmle_movement_20260809"


def load_plan(path: Path) -> dict:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != PLAN_SCHEMA:
        raise ValueError(f"movement plan schema must be {PLAN_SCHEMA}")
    reproduction = plan.get("reproduction", {})
    localization = plan.get("localization", {})
    if int(reproduction.get("attempted_runids", 0)) != 160:
        raise ValueError("movement audit must retain all 160 RunIds")
    if reproduction.get("polish_method") != "L-BFGS-B":
        raise ValueError("movement audit polish method changed")
    expected_options = {"maxiter": 1500, "gtol": 1e-7, "ftol": 1e-15, "maxls": 50, "maxcor": 20}
    if reproduction.get("polish_options") != expected_options:
        raise ValueError("movement audit polish options changed")
    if float(localization.get("large_movement_descriptive_threshold", np.nan)) != 1.0:
        raise ValueError("movement audit large-movement flag changed")
    if not bool(localization.get("thresholds_are_post_stage_b2_diagnostic_flags_not_qualification_gates", False)):
        raise ValueError("movement flags must remain descriptive")
    boundaries = plan.get("scope_boundaries", {})
    if boundaries.get("public_surface_status") != "Withheld":
        raise ValueError("movement public surface must remain Withheld")
    if bool(boundaries.get("application_core_change_authorized", True)):
        raise ValueError("movement audit cannot authorize a core change")
    return plan


def validate_input_identity(input_dir: Path, b2_dir: Path, plan: dict) -> pd.DataFrame:
    expected = plan["input_identity"]
    paths = {
        "streamlit_app_sha256": REPO / "streamlit_app.py",
        "stage_a2_adapter_sha256": Path(stage_a2.__file__).resolve(),
        "stage_b2_plan_sha256": REPO / "validation" / "operating_characteristics_strict_jmle_b2_plan_20260809.json",
        "stage_b2_adapter_sha256": Path(stage_b2.__file__).resolve(),
        "stage_b2_runs_sha256": b2_dir / "strict_jmle_b2_runs.csv",
        "stage_b2_baseline_sha256": b2_dir / "strict_jmle_b2_baseline.csv",
        "stage_b2_identity_sha256": b2_dir / "strict_jmle_b2_identity.json",
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
    identity = json.loads((b2_dir / "strict_jmle_b2_identity.json").read_text(encoding="utf-8"))
    rows.append({
        "Check": "stage_b2_decision_chain",
        "Passed": bool(
            identity.get("stage_b2_numerically_qualified") is True
            and identity.get("integration_experiment_authorized") is True
            and identity.get("application_core_change_authorized") is False
        ),
        "Evidence": (
            f"qualified={identity.get('stage_b2_numerically_qualified')}; "
            f"experiment={identity.get('integration_experiment_authorized')}; "
            f"core={identity.get('application_core_change_authorized')}"
        ),
    })
    checks = pd.DataFrame(rows)
    if not checks["Passed"].all():
        failed = checks.loc[~checks["Passed"], "Check"].astype(str).tolist()
        raise ValueError(f"movement input identity failed: {failed}")
    return checks


def constrained_free_labels(spec: dict) -> list[dict[str, object]]:
    anchors = np.asarray(spec["anchors"], dtype=float)
    groups = list(spec["groups"])
    levels = [str(value) for value in spec["levels"]]
    free_idx = np.where(np.isnan(anchors))[0].tolist()
    labels: list[dict[str, object]] = []
    grouped_free: set[int] = set()
    group_ids = sorted({
        groups[index]
        for index in free_idx
        if groups[index] not in (None, "", np.nan)
    })
    for group_id in group_ids:
        group_levels = [index for index, value in enumerate(groups) if value == group_id]
        free_in_group = [index for index in group_levels if np.isnan(anchors[index])]
        grouped_free.update(free_in_group)
        if len(free_in_group) > 1:
            derived = levels[free_in_group[-1]]
            for index in free_in_group[:-1]:
                labels.append({
                    "Level": levels[index],
                    "DerivedCounterpart": derived,
                    "CoordinateMeaning": f"group {group_id}: direct level with final level derived",
                })
    ungrouped = [
        index for index in free_idx
        if index not in grouped_free and groups[index] in (None, "", np.nan)
    ]
    if bool(spec["centered"]):
        if len(ungrouped) > 1:
            derived = levels[ungrouped[-1]]
            for index in ungrouped[:-1]:
                labels.append({
                    "Level": levels[index],
                    "DerivedCounterpart": derived,
                    "CoordinateMeaning": "direct level with final unanchored level derived by constraint",
                })
    else:
        for index in ungrouped:
            labels.append({
                "Level": levels[index],
                "DerivedCounterpart": "",
                "CoordinateMeaning": "direct unconstrained level",
            })
    return labels


def free_coordinate_metadata(sizes: dict, config: dict) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    global_index = 0
    for block, size_raw in sizes.items():
        size = int(size_raw)
        if size == 0:
            continue
        if block == "theta":
            labels = constrained_free_labels(config["theta_spec"])
        elif block in config["facet_names"]:
            labels = constrained_free_labels(config["facet_specs"][block])
        elif block == "steps":
            labels = [
                {
                    "Level": f"Step_{index + 1}",
                    "DerivedCounterpart": f"Step_{int(config['n_cat']) - 1}",
                    "CoordinateMeaning": "direct step with final step derived by exact sum-to-zero",
                }
                for index in range(size)
            ]
        else:
            labels = [
                {
                    "Level": f"Coordinate_{index + 1}",
                    "DerivedCounterpart": "",
                    "CoordinateMeaning": "free optimizer coordinate",
                }
                for index in range(size)
            ]
        if len(labels) != size:
            raise RuntimeError(f"free-coordinate labels for {block} have {len(labels)} rows; expected {size}")
        for local_index, label in enumerate(labels):
            rows.append({
                "GlobalCoordinate": global_index,
                "Block": str(block),
                "BlockCoordinate": local_index,
                **label,
            })
            global_index += 1
    return rows


def response_patterns(data: pd.DataFrame, rating_min: int, rating_max: int) -> pd.DataFrame:
    grouped = (
        data.assign(Person=data["Person"].astype(str))
        .groupby("Person", sort=False)["Score"]
        .agg(ObservedCount="count", ScoreMin="min", ScoreMax="max", ScoreMean="mean", ScoreTotal="sum")
        .reset_index()
    )
    grouped["MinimumPossibleTotal"] = int(rating_min) * grouped["ObservedCount"]
    grouped["MaximumPossibleTotal"] = int(rating_max) * grouped["ObservedCount"]
    grouped["ExtremeAllMinimum"] = grouped["ScoreMin"].eq(rating_min) & grouped["ScoreMax"].eq(rating_min)
    grouped["ExtremeAllMaximum"] = grouped["ScoreMin"].eq(rating_max) & grouped["ScoreMax"].eq(rating_max)
    grouped["ExtremeScorePattern"] = grouped["ExtremeAllMinimum"] | grouped["ExtremeAllMaximum"]
    span = max(int(rating_max) - int(rating_min), 1)
    mean_proportion = (grouped["ScoreMean"] - int(rating_min)) / span
    grouped["MeanScoreProportion"] = mean_proportion
    grouped["DistanceFromExtremeMeanProportion"] = np.minimum(mean_proportion, 1.0 - mean_proportion)
    return grouped


def exact_coordinate_rows(
    identity: dict[str, object],
    start: np.ndarray,
    polished: np.ndarray,
    sizes: dict,
    config: dict,
    patterns: pd.DataFrame,
    plan: dict,
) -> pd.DataFrame:
    metadata = pd.DataFrame(free_coordinate_metadata(sizes, config))
    if len(metadata) != len(start):
        raise RuntimeError(f"coordinate metadata has {len(metadata)} rows; expected {len(start)}")
    metadata["BaselineValue"] = np.asarray(start, dtype=float)
    metadata["PolishedValue"] = np.asarray(polished, dtype=float)
    metadata["Movement"] = metadata["PolishedValue"] - metadata["BaselineValue"]
    metadata["AbsMovement"] = metadata["Movement"].abs()
    metadata["MaterialMovement"] = metadata["AbsMovement"].ge(
        float(plan["localization"]["material_movement_descriptive_threshold"])
    )
    metadata["LargeMovement"] = metadata["AbsMovement"].ge(
        float(plan["localization"]["large_movement_descriptive_threshold"])
    )
    metadata["ExtremeScorePattern"] = False
    metadata["ObservedCount"] = np.nan
    metadata["ScoreMin"] = np.nan
    metadata["ScoreMax"] = np.nan
    metadata["ScoreMean"] = np.nan
    metadata["ScoreTotal"] = np.nan
    person_lookup = patterns.set_index("Person")
    person_mask = metadata["Block"].eq("theta")
    for index, level in metadata.loc[person_mask, "Level"].items():
        if str(level) not in person_lookup.index:
            continue
        pattern = person_lookup.loc[str(level)]
        for column in ("ExtremeScorePattern", "ObservedCount", "ScoreMin", "ScoreMax", "ScoreMean", "ScoreTotal"):
            metadata.loc[index, column] = pattern[column]
    for column, value in reversed(identity.items()):
        metadata.insert(0, column, value)
    return metadata


def expanded_parameter_rows(
    identity: dict[str, object],
    start: np.ndarray,
    polished: np.ndarray,
    sizes: dict,
    config: dict,
    patterns: pd.DataFrame,
    anchors: pd.DataFrame,
    plan: dict,
) -> pd.DataFrame:
    baseline = app.expand_params(start, sizes, config)
    final = app.expand_params(polished, sizes, config)
    rows: list[dict[str, object]] = []
    for level, before, after in zip(config["theta_spec"]["levels"], baseline["theta"], final["theta"]):
        rows.append({"Block": "theta", "Level": str(level), "BaselineValue": before, "PolishedValue": after})
    for facet in config["facet_names"]:
        for level, before, after in zip(config["facet_levels"][facet], baseline["facets"][facet], final["facets"][facet]):
            rows.append({"Block": facet, "Level": str(level), "BaselineValue": before, "PolishedValue": after})
    for index, (before, after) in enumerate(zip(baseline["steps"], final["steps"])):
        rows.append({"Block": "steps", "Level": f"Step_{index + 1}", "BaselineValue": before, "PolishedValue": after})
    frame = pd.DataFrame(rows)
    frame["Movement"] = frame["PolishedValue"] - frame["BaselineValue"]
    frame["AbsMovement"] = frame["Movement"].abs()
    frame["MaterialMovement"] = frame["AbsMovement"].ge(
        float(plan["localization"]["material_movement_descriptive_threshold"])
    )
    frame["LargeMovement"] = frame["AbsMovement"].ge(
        float(plan["localization"]["large_movement_descriptive_threshold"])
    )
    anchor_lookup = {
        (str(row.Facet), str(row.Level)): float(row.Anchor)
        for row in anchors.itertuples(index=False)
    }
    frame["Anchored"] = [
        (str(row.Block), str(row.Level)) in anchor_lookup
        for row in frame.itertuples(index=False)
    ]
    frame["AnchorValue"] = [
        anchor_lookup.get((str(row.Block), str(row.Level)), np.nan)
        for row in frame.itertuples(index=False)
    ]
    pattern_columns = [
        "Person", "ObservedCount", "ScoreMin", "ScoreMax", "ScoreMean", "ScoreTotal",
        "MinimumPossibleTotal", "MaximumPossibleTotal", "ExtremeAllMinimum",
        "ExtremeAllMaximum", "ExtremeScorePattern", "MeanScoreProportion",
        "DistanceFromExtremeMeanProportion",
    ]
    frame = frame.merge(
        patterns[pattern_columns],
        left_on="Level",
        right_on="Person",
        how="left",
    ).drop(columns="Person")
    frame.loc[~frame["Block"].eq("theta"), "ExtremeScorePattern"] = False
    for column, value in reversed(identity.items()):
        frame.insert(0, column, value)
    return frame


def reproduce_run(
    manifest_row: pd.Series,
    generated,
    strict_controls: dict,
    stored_baseline: pd.Series,
    stored_run: pd.Series,
    plan: dict,
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame]:
    result, _, _ = stage_a2.fit_stage_a_baseline(manifest_row, generated, strict_controls)
    start, idx, config, sizes = stage_a2.optimizer_context(result)
    baseline_objective, baseline_gradient = app.mfrm_loglik_jmle_value_grad(start, idx, config, sizes)
    optimizer = minimize(
        app.mfrm_loglik_jmle_value_grad,
        np.array(start, dtype=float, copy=True),
        args=(idx, config, sizes),
        jac=True,
        method=plan["reproduction"]["polish_method"],
        bounds=app.build_optimizer_bounds(sizes, config),
        options=dict(plan["reproduction"]["polish_options"]),
    )
    polished = np.asarray(optimizer.x, dtype=float)
    objective, gradient = app.mfrm_loglik_jmle_value_grad(polished, idx, config, sizes)
    baseline_gradient_sup = float(np.max(np.abs(baseline_gradient))) if len(baseline_gradient) else 0.0
    gradient_sup = float(np.max(np.abs(gradient))) if len(gradient) else 0.0
    movement_max = float(np.max(np.abs(polished - start))) if len(start) else 0.0
    reproduction = plan["reproduction"]
    comparisons = {
        "BaselineObjectiveDifference": abs(float(baseline_objective) - float(stored_baseline["StrictBaselineObjective"])),
        "BaselineGradientDifference": abs(baseline_gradient_sup - float(stored_baseline["StrictBaselineTerminalGradientSupNorm"])),
        "PolishedObjectiveDifference": abs(float(objective) - float(stored_run["Objective"])),
        "PolishedGradientDifference": abs(gradient_sup - float(stored_run["TerminalGradientSupNorm"])),
        "MaximumMovementDifference": abs(movement_max - float(stored_run["MaxAbsParameterMovement"])),
    }
    reproduced = bool(
        comparisons["BaselineObjectiveDifference"] <= float(reproduction["baseline_objective_max_abs_difference_tolerance"])
        and comparisons["BaselineGradientDifference"] <= float(reproduction["baseline_gradient_sup_norm_max_abs_difference_tolerance"])
        and comparisons["PolishedObjectiveDifference"] <= float(reproduction["polished_objective_max_abs_difference_tolerance"])
        and comparisons["PolishedGradientDifference"] <= float(reproduction["polished_gradient_sup_norm_max_abs_difference_tolerance"])
        and comparisons["MaximumMovementDifference"] <= float(reproduction["maximum_movement_max_abs_difference_tolerance"])
    )
    identity = {
        "RunId": str(manifest_row["RunId"]),
        "ConditionId": str(manifest_row["ConditionId"]),
        "Design": str(manifest_row["Design"]),
        "TruthBias": float(manifest_row["TruthBias"]),
        "Replicate": int(manifest_row["Replicate"]),
    }
    patterns = response_patterns(generated.data, 0, int(manifest_row["Categories"]) - 1)
    exact = exact_coordinate_rows(identity, start, polished, sizes, config, patterns, plan)
    expanded = expanded_parameter_rows(identity, start, polished, sizes, config, patterns, generated.anchors, plan)
    large = exact.loc[exact["LargeMovement"].astype(bool)]
    large_theta = large.loc[large["Block"].eq("theta")]
    large_non_theta = large.loc[~large["Block"].eq("theta")]
    large_nonextreme_theta = large_theta.loc[~large_theta["ExtremeScorePattern"].astype(bool)]
    if large.empty:
        localization_class = "no_large_movement"
    elif not large_non_theta.empty:
        localization_class = "non_theta_large_movement"
    elif not large_nonextreme_theta.empty:
        localization_class = "non_extreme_theta_large_movement"
    else:
        localization_class = "extreme_theta_only"
    non_theta = exact.loc[~exact["Block"].eq("theta"), "AbsMovement"]
    theta_extreme = exact.loc[exact["Block"].eq("theta") & exact["ExtremeScorePattern"].astype(bool), "AbsMovement"]
    theta_nonextreme = exact.loc[exact["Block"].eq("theta") & ~exact["ExtremeScorePattern"].astype(bool), "AbsMovement"]
    summary = {
        **identity,
        "Reproduced": reproduced,
        **comparisons,
        "OptimizerSuccess": bool(optimizer.success),
        "OptimizerMessage": str(optimizer.message),
        "BaselineGradientSupNorm": baseline_gradient_sup,
        "PolishedGradientSupNorm": gradient_sup,
        "LogLikGain": float(baseline_objective - objective),
        "MaximumAbsFreeCoordinateMovement": movement_max,
        "LargeMovementCoordinates": int(len(large)),
        "LargeThetaCoordinates": int(len(large_theta)),
        "LargeNonThetaCoordinates": int(len(large_non_theta)),
        "LargeNonExtremeThetaCoordinates": int(len(large_nonextreme_theta)),
        "ExtremePersons": int(patterns["ExtremeScorePattern"].astype(bool).sum()),
        "MaxAbsNonThetaMovement": float(non_theta.max()) if len(non_theta) else 0.0,
        "MaxAbsExtremeThetaMovement": float(theta_extreme.max()) if len(theta_extreme) else 0.0,
        "MaxAbsNonExtremeThetaMovement": float(theta_nonextreme.max()) if len(theta_nonextreme) else 0.0,
        "LocalizationClass": localization_class,
    }
    return summary, exact, expanded


def reproduction_gates(summary: pd.DataFrame, input_checks: pd.DataFrame, plan: dict) -> pd.DataFrame:
    expected = int(plan["reproduction"]["attempted_runids"])
    reproduction = plan["reproduction"]
    definitions = [
        ("input_identity", len(input_checks), int(input_checks["Passed"].sum())),
        ("runids_returned", expected, len(summary)),
        ("baseline_objective_reproduction", expected, int(summary["BaselineObjectiveDifference"].le(float(reproduction["baseline_objective_max_abs_difference_tolerance"])).sum())),
        ("baseline_gradient_reproduction", expected, int(summary["BaselineGradientDifference"].le(float(reproduction["baseline_gradient_sup_norm_max_abs_difference_tolerance"])).sum())),
        ("polished_objective_reproduction", expected, int(summary["PolishedObjectiveDifference"].le(float(reproduction["polished_objective_max_abs_difference_tolerance"])).sum())),
        ("polished_gradient_reproduction", expected, int(summary["PolishedGradientDifference"].le(float(reproduction["polished_gradient_sup_norm_max_abs_difference_tolerance"])).sum())),
        ("maximum_movement_reproduction", expected, int(summary["MaximumMovementDifference"].le(float(reproduction["maximum_movement_max_abs_difference_tolerance"])).sum())),
        ("complete_reproduction", expected, int(summary["Reproduced"].astype(bool).sum())),
    ]
    rows = [
        {"Gate": gate, "Required": required, "Passed": passed, "GatePassed": required == passed}
        for gate, required, passed in definitions
    ]
    return pd.DataFrame(rows)


def block_summary(exact: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for block, frame in exact.groupby("Block", sort=False):
        values = pd.to_numeric(frame["AbsMovement"], errors="coerce")
        rows.append({
            "Block": block,
            "Coordinates": len(frame),
            "MedianAbsMovement": float(values.median()),
            "P95AbsMovement": float(values.quantile(0.95)),
            "P99AbsMovement": float(values.quantile(0.99)),
            "MaximumAbsMovement": float(values.max()),
            "MaterialMovementCoordinates": int(frame["MaterialMovement"].astype(bool).sum()),
            "LargeMovementCoordinates": int(frame["LargeMovement"].astype(bool).sum()),
            "LargeExtremeScoreCoordinates": int(
                (frame["LargeMovement"].astype(bool) & frame["ExtremeScorePattern"].astype(bool)).sum()
            ),
        })
    return pd.DataFrame(rows)


def policy_decision(summary: pd.DataFrame, gates: pd.DataFrame) -> pd.DataFrame:
    reproduced = bool(gates["GatePassed"].all())
    classes = summary["LocalizationClass"].value_counts()
    non_theta = int(summary["LargeNonThetaCoordinates"].gt(0).sum())
    nonextreme = int(summary["LargeNonExtremeThetaCoordinates"].gt(0).sum())
    extreme_only = int(classes.get("extreme_theta_only", 0))
    if not reproduced:
        disposition = "blocked_reproduction_failure"
        action = "Repair audit reproducibility before interpreting movement localization."
    elif non_theta > 0:
        disposition = "blocked_non_theta_instability"
        action = "Investigate facet/step instability before downstream reconstruction."
    elif nonextreme > 0:
        disposition = "blocked_nonextreme_theta_instability"
        action = "Investigate sparse connectivity, quasi-separation, and information before integration."
    elif extreme_only > 0:
        disposition = "requires_explicit_extreme_score_policy"
        action = "Design and qualify an extreme-person policy; do not report arbitrary large finite JMLE person measures as stable."
    else:
        disposition = "movement_localization_clear"
        action = "Proceed to downstream reconstruction while retaining convergence and floating-point audits."
    return pd.DataFrame([
        {
            "Decision": "MovementLocalization",
            "Disposition": disposition,
            "Evidence": (
                f"reproduced={reproduced}; extreme_theta_only_runs={extreme_only}; "
                f"nonextreme_theta_runs={nonextreme}; non_theta_runs={non_theta}"
            ),
            "NextAction": action,
            "ApplicationCoreChangeAuthorized": False,
            "PublicSurfaceEnabled": False,
        }
    ])


def plot_person_boundary(expanded: pd.DataFrame, output: Path) -> None:
    person = expanded.loc[expanded["Block"].eq("theta")].copy()
    distance = pd.to_numeric(person["DistanceFromExtremeMeanProportion"], errors="coerce")
    movement = pd.to_numeric(person["AbsMovement"], errors="coerce")
    extreme = person["ExtremeScorePattern"].eq(True).to_numpy(dtype=bool)
    fig, ax = plt.subplots(figsize=(10.5, 6.3))
    ax.scatter(distance[~extreme], movement[~extreme], s=14, alpha=0.35, color="#4E79A7", label="non-extreme response pattern")
    ax.scatter(distance[extreme], movement[extreme], s=38, alpha=0.85, color="#E15759", label="all-minimum or all-maximum")
    ax.axhline(1.0, color="#222222", linestyle="--", linewidth=1, label="descriptive large-movement flag")
    ax.set_yscale("log")
    ax.set_xlabel("Distance of mean score proportion from 0 or 1")
    ax.set_ylabel("Absolute expanded person-measure movement")
    ax.set_title("JMLE polish movement localizes against response-score boundaries")
    ax.grid(axis="y", which="both", alpha=0.22)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "jmle_movement_person_boundary.png", dpi=180)
    plt.close(fig)


def plot_blocks(exact: pd.DataFrame, output: Path) -> None:
    order = exact["Block"].drop_duplicates().tolist()
    values = [pd.to_numeric(exact.loc[exact["Block"].eq(block), "AbsMovement"], errors="coerce").dropna() for block in order]
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    ax.boxplot(values, tick_labels=order, showmeans=True, showfliers=True)
    ax.axhline(1.0, color="#E15759", linestyle="--", linewidth=1, label="descriptive large-movement flag")
    ax.axhline(0.1, color="#F28E2B", linestyle=":", linewidth=1, label="descriptive material-movement flag")
    ax.set_yscale("log")
    ax.set_ylabel("Absolute exact-free-coordinate movement")
    ax.set_title("Movement localization by optimizer parameter block")
    ax.grid(axis="y", which="both", alpha=0.22)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output / "jmle_movement_by_block.png", dpi=180)
    plt.close(fig)


def write_report(
    output: Path,
    summary: pd.DataFrame,
    exact: pd.DataFrame,
    blocks: pd.DataFrame,
    gates: pd.DataFrame,
    policy: pd.DataFrame,
) -> None:
    classes = summary["LocalizationClass"].value_counts()
    non_theta_runs = int(summary["LargeNonThetaCoordinates"].gt(0).sum())
    nonextreme_theta_runs = int(summary["LargeNonExtremeThetaCoordinates"].gt(0).sum())
    large = exact.loc[exact["LargeMovement"].astype(bool)]
    top = exact.nlargest(10, "AbsMovement")
    top_lines = "\n".join(
        f"- `{row.RunId}` / `{row.Block}` / `{row.Level}`: {row.AbsMovement:.6g}; extreme={bool(row.ExtremeScorePattern)}"
        for row in top.itertuples(index=False)
    )
    failed = gates.loc[~gates["GatePassed"]]
    failed_text = "\n".join(
        f"- `{row.Gate}`: {int(row.Passed)}/{int(row.Required)}"
        for row in failed.itertuples(index=False)
    ) or "- None."
    disposition = str(policy.iloc[0]["Disposition"])
    report = f"""# JMLE Stage-B2 movement localization audit

## Decision

Disposition: `{disposition}`. This audit does not alter the completed Stage-B2
numerical result and does not authorize an application-core change or public
surface.

## Reproduction

- {int(summary['Reproduced'].astype(bool).sum())}/{len(summary)} Stage-B2 RunIds reproduced within every frozen tolerance.
- Failed reproduction gates:
{failed_text}

## Localization

- Runs with no coordinate movement >= 1 logit: {int(classes.get('no_large_movement', 0))}.
- Runs whose >=1-logit movements were confined to extreme-score theta coordinates: {int(classes.get('extreme_theta_only', 0))}.
- Runs with any >=1-logit non-extreme theta movement: {nonextreme_theta_runs}.
- Runs with any >=1-logit non-theta movement: {non_theta_runs}.
- Large exact-free-coordinate movements: {len(large)}; of these, {int((large['Block'].eq('theta') & large['ExtremeScorePattern'].astype(bool)).sum())} were extreme-score theta coordinates.
- Maximum non-theta movement across all runs: {summary['MaxAbsNonThetaMovement'].max():.6g}.
- Maximum non-extreme theta movement across all runs: {summary['MaxAbsNonExtremeThetaMovement'].max():.6g}.
- Maximum extreme-score theta movement across all runs: {summary['MaxAbsExtremeThetaMovement'].max():.6g}.

The 1.0 and 0.1 logit flags were registered after Stage B2 solely for
localization; they are not retroactive Stage-B2 qualification gates.

## Ten largest exact-coordinate movements

{top_lines}

## Interpretation boundary

JMLE person measures for all-minimum or all-maximum response patterns have no
finite maximum-likelihood estimate. A numerical optimizer can return a finite
value whose magnitude depends on stopping rules while the likelihood changes
only imperceptibly. Such values must not be described as stable merely because
the gradient is below a numerical threshold. Facet and step localization is
reported separately so an extreme-person policy is not used to conceal broader
parameter instability. Bias decisions are not recomputed here.

## Retained artifacts

- `jmle_movement_run_summary.csv`
- `jmle_movement_exact_coordinates.csv`
- `jmle_movement_expanded_parameters.csv`
- `jmle_movement_block_summary.csv`
- `jmle_movement_reproduction_gates.csv`
- `jmle_movement_policy_decision.csv`
- `jmle_movement_input_checks.csv`
- `jmle_movement_first_read_summary.csv`
- `jmle_movement_identity.json`
- `jmle_movement_person_boundary.png`
- `jmle_movement_by_block.png`
"""
    (output / "JMLE_MOVEMENT_AUDIT.md").write_text(report, encoding="utf-8")


def run(input_dir: Path, b2_dir: Path, stage_a_plan_path: Path, output: Path, plan_path: Path) -> None:
    input_dir = input_dir.resolve()
    b2_dir = b2_dir.resolve()
    stage_a_plan_path = stage_a_plan_path.resolve()
    output = output.resolve()
    plan_path = plan_path.resolve()
    plan = load_plan(plan_path)
    input_checks = validate_input_identity(input_dir, b2_dir, plan)
    strict_plan = stage_a.load_plan(stage_a_plan_path)
    manifest = pd.read_csv(input_dir / "manifest.csv")
    ratings = pd.read_csv(input_dir / "generated_ratings.csv")
    truth = pd.read_csv(input_dir / "generated_facet_truth.csv")
    anchors = pd.read_csv(input_dir / "generated_anchors.csv")
    stored_baseline = pd.read_csv(b2_dir / "strict_jmle_b2_baseline.csv").set_index("RunId", drop=False)
    stored_runs = pd.read_csv(b2_dir / "strict_jmle_b2_runs.csv").set_index("RunId", drop=False)
    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_jmle_movement_plan.json").write_bytes(plan_path.read_bytes())
    summary_rows: list[dict[str, object]] = []
    exact_parts: list[pd.DataFrame] = []
    expanded_parts: list[pd.DataFrame] = []
    started = time.perf_counter()
    for _, manifest_row in manifest.iterrows():
        run_id = str(manifest_row["RunId"])
        generated = stage_a.generated_design_for_run(manifest_row, ratings, truth, anchors)
        summary, exact, expanded = reproduce_run(
            manifest_row,
            generated,
            strict_plan["numerical_controls"],
            stored_baseline.loc[run_id],
            stored_runs.loc[run_id],
            plan,
        )
        summary_rows.append(summary)
        exact_parts.append(exact)
        expanded_parts.append(expanded)
    summary = pd.DataFrame(summary_rows)
    exact = pd.concat(exact_parts, ignore_index=True)
    expanded = pd.concat(expanded_parts, ignore_index=True)
    blocks = block_summary(exact)
    gates = reproduction_gates(summary, input_checks, plan)
    policy = policy_decision(summary, gates)
    first_read = pd.DataFrame([
        {
            "Priority": 1,
            "Check": "Movement localization",
            "Status": str(policy.iloc[0]["Disposition"]),
            "Evidence": str(policy.iloc[0]["Evidence"]),
            "NextAction": str(policy.iloc[0]["NextAction"]),
            "PublicSurfaceEnabled": False,
        },
        {
            "Priority": 2,
            "Check": "Facet/step stability",
            "Status": "Descriptive",
            "Evidence": f"maximum non-theta movement={summary['MaxAbsNonThetaMovement'].max():.6g}",
            "NextAction": "Keep separate from extreme-person handling and verify before bias reconstruction.",
            "PublicSurfaceEnabled": False,
        },
        {
            "Priority": 3,
            "Check": "Application core",
            "Status": "Withheld",
            "Evidence": "This audit authorizes no core change.",
            "NextAction": "Require an explicit extreme-score policy and downstream reconstruction plan.",
            "PublicSurfaceEnabled": False,
        },
    ])
    summary.to_csv(output / "jmle_movement_run_summary.csv", index=False)
    exact.to_csv(output / "jmle_movement_exact_coordinates.csv", index=False)
    expanded.to_csv(output / "jmle_movement_expanded_parameters.csv", index=False)
    blocks.to_csv(output / "jmle_movement_block_summary.csv", index=False)
    gates.to_csv(output / "jmle_movement_reproduction_gates.csv", index=False)
    policy.to_csv(output / "jmle_movement_policy_decision.csv", index=False)
    input_checks.to_csv(output / "jmle_movement_input_checks.csv", index=False)
    first_read.to_csv(output / "jmle_movement_first_read_summary.csv", index=False)
    plot_person_boundary(expanded, output)
    plot_blocks(exact, output)
    write_report(output, summary, exact, blocks, gates, policy)
    artifact_names = [
        "jmle_movement_run_summary.csv",
        "jmle_movement_exact_coordinates.csv",
        "jmle_movement_expanded_parameters.csv",
        "jmle_movement_block_summary.csv",
        "jmle_movement_reproduction_gates.csv",
        "jmle_movement_policy_decision.csv",
        "jmle_movement_input_checks.csv",
        "jmle_movement_first_read_summary.csv",
        "jmle_movement_person_boundary.png",
        "jmle_movement_by_block.png",
        "JMLE_MOVEMENT_AUDIT.md",
    ]
    identity = {
        "schema_version": PLAN_SCHEMA,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "wall_seconds": time.perf_counter() - started,
        "disposition": str(policy.iloc[0]["Disposition"]),
        "application_core_change_authorized": False,
        "public_surface_enabled": False,
        "movement_plan_sha256": stage_a.sha256_file(plan_path),
        "movement_adapter_sha256": stage_a.sha256_file(Path(__file__).resolve()),
        "streamlit_app_sha256": stage_a.sha256_file(REPO / "streamlit_app.py"),
        "artifacts": {name: stage_a.sha256_file(output / name) for name in artifact_names},
    }
    (output / "jmle_movement_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--stage-b2", type=Path, default=DEFAULT_B2)
    parser.add_argument("--stage-a-plan", type=Path, default=DEFAULT_STAGE_A_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    args = parser.parse_args()
    run(args.input, args.stage_b2, args.stage_a_plan, args.output, args.plan)


if __name__ == "__main__":
    main()
