#!/usr/bin/env python3
"""Run the frozen CMLE-WLE Person MnSq cross-engine validation."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle import fit_cmle  # noqa: E402
from mfrm_app.cmle_wle_fit import compute_cmle_wle_person_fit  # noqa: E402


DEFAULT_PLAN = ROOT / "validation/cmle_wle_person_fit_plan_20260810.json"
DEFAULT_FIXTURE = (
    ROOT / "validation/cmle_wle_bootstrap_pilot_20260810/bootstrap_fixture.csv"
)
DEFAULT_OUTPUT = ROOT / "validation/cmle_wle_person_fit_reproduction_current"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_plan(plan_path: Path) -> dict[str, object]:
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != "mfrm-cmle-wle-person-fit-plan-v1":
        raise ValueError("Unexpected CMLE-WLE Person-fit plan schema.")
    paths = {
        "cmle_core_sha256": ROOT / "mfrm_app/cmle.py",
        "fixed_calibration_person_scoring_sha256": ROOT / "mfrm_app/person_scoring.py",
        "cmle_wle_bridge_sha256": ROOT / "mfrm_app/cmle_person_scoring.py",
        "decision_stability_sha256": ROOT / "mfrm_app/decision_stability.py",
        "cmle_wle_bridge_decision_sha256": (
            ROOT / "validation/cmle_wle_bridge_20260810/cmle_wle_bridge_decision.json"
        ),
        "fixed_calibration_wle_parity_decision_sha256": (
            ROOT / "validation/fixed_calibration_wle_20260809/parity_decision.json"
        ),
        "bootstrap_pilot_decision_sha256": (
            ROOT
            / "validation/cmle_wle_bootstrap_pilot_20260810/bootstrap_pilot_decision.json"
        ),
        "bootstrap_design_sha256": ROOT / "docs/cmle_wle_bootstrap_design.md",
        "bootstrap_plan_sha256": ROOT / "validation/cmle_wle_bootstrap_plan_20260810.json",
    }
    failed = [
        key
        for key, path in paths.items()
        if not path.exists() or sha256_file(path) != str(plan["input_identity"][key])
    ]
    original_bootstrap_decision = json.loads(
        (
            ROOT
            / "validation/cmle_wle_bootstrap_pilot_20260810/bootstrap_pilot_decision.json"
        ).read_text(encoding="utf-8")
    )
    if str(original_bootstrap_decision.get("bootstrap_core_sha256")) != str(
        plan["input_identity"]["bootstrap_core_pre_fit_extension_sha256"]
    ):
        failed.append("bootstrap_core_pre_fit_extension_sha256")
    if failed:
        raise ValueError(f"CMLE-WLE Person-fit input identity failed: {failed}")
    return plan


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(
        path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )


def markdown_table(frame: pd.DataFrame) -> str:
    rendered = frame.copy()
    for column in rendered.columns:
        if pd.api.types.is_float_dtype(rendered[column]):
            rendered[column] = rendered[column].map(
                lambda value: "" if pd.isna(value) else f"{value:.10g}"
            )
        elif pd.api.types.is_bool_dtype(rendered[column]):
            rendered[column] = rendered[column].map(
                lambda value: "TRUE" if value else "FALSE"
            )
    header = "| " + " | ".join(rendered.columns.astype(str)) + " |"
    separator = "| " + " | ".join(["---"] * len(rendered.columns)) + " |"
    rows = [
        "| " + " | ".join(str(value) for value in row) + " |"
        for row in rendered.itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *rows])


def fit_fixture(frame: pd.DataFrame, model: str) -> dict[str, object]:
    return fit_cmle(
        frame,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model=model,
        step_facet="Criterion" if model == "PCM" else None,
        gtol=1e-8,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    plan = validate_plan(args.plan)
    if not args.fixture.exists():
        raise FileNotFoundError(args.fixture)
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    fixture = pd.read_csv(args.fixture)
    required = {"Person", "Rater", "Criterion", "Score"}
    if set(fixture.columns) != required or fixture.empty:
        raise ValueError("Fixture must contain exactly Person, Rater, Criterion, Score.")

    calibration_frames: list[pd.DataFrame] = []
    response_frames: list[pd.DataFrame] = []
    theta_frames: list[pd.DataFrame] = []
    python_person_frames: list[pd.DataFrame] = []
    python_observation_frames: list[pd.DataFrame] = []
    model_rows: list[dict[str, object]] = []
    for model in plan["validation_matrix"]["models"]:
        fit = fit_fixture(fixture, str(model))
        if not bool(fit["summary"].iloc[0]["InferenceReady"]):
            raise RuntimeError(f"Frozen {model} fixture is not inference-ready.")
        result = compute_cmle_wle_person_fit(fit)
        persons = result["persons"].copy()
        observations = result["observations"].copy()
        if not bool(persons["PersonFitReady"].all()):
            raise RuntimeError(f"Frozen {model} fixture has unavailable Person fit.")
        persons.insert(0, "Model", str(model))
        observations.insert(0, "Model", str(model))
        python_person_frames.append(persons)
        python_observation_frames.append(observations)

        surfaces = fit["surfaces"].loc[
            fit["surfaces"]["InternalCategory"].gt(0),
            ["VirtualUnit", "InternalCategory", "CumulativeDifficulty"],
        ].copy()
        surfaces.insert(0, "Model", str(model))
        calibration_frames.append(surfaces)
        response_frames.append(
            observations[
                ["Model", "Person", "VirtualUnit", "ObservedInternalCategory"]
            ].copy()
        )
        theta_frames.append(
            persons[["Model", "Person", "WLEEstimate"]].copy()
        )
        model_rows.append(
            {
                "Model": str(model),
                "Persons": int(len(persons)),
                "Observations": int(len(observations)),
                "ExactExtremePersons": int(persons["ExtremeScorePattern"].sum()),
                "MissingCellsFromRectangularDesign": int(
                    len(persons) * fixture[["Rater", "Criterion"]].drop_duplicates().shape[0]
                    - len(observations)
                ),
                "CMLEInferenceReady": True,
                "PythonPersonFitReady": bool(persons["PersonFitReady"].all()),
            }
        )

    calibration = pd.concat(calibration_frames, ignore_index=True)
    responses = pd.concat(response_frames, ignore_index=True)
    theta = pd.concat(theta_frames, ignore_index=True)
    python_persons = pd.concat(python_person_frames, ignore_index=True)
    python_observations = pd.concat(python_observation_frames, ignore_index=True)
    model_summary = pd.DataFrame(model_rows)
    write_csv(calibration, output / "sirt_calibration.csv")
    write_csv(responses, output / "sirt_responses.csv")
    write_csv(theta, output / "sirt_theta.csv")
    write_csv(python_persons, output / "python_person_fit.csv")
    write_csv(python_observations, output / "python_observation_moments.csv")
    write_csv(model_summary, output / "fixture_model_summary.csv")

    manifest_files = {}
    for filename in ("sirt_calibration.csv", "sirt_responses.csv", "sirt_theta.csv"):
        path = output / filename
        manifest_files[filename] = {
            "path": str(path.relative_to(ROOT)),
            "sha256": sha256_file(path),
        }
    manifest = {
        "schema_version": "mfrm-cmle-wle-person-fit-fixture-v1",
        "models": [str(value) for value in plan["validation_matrix"]["models"]],
        "source_fixture": str(args.fixture.relative_to(ROOT)),
        "source_fixture_sha256": sha256_file(args.fixture),
        "files": manifest_files,
    }
    (output / "fixture_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (output / "registered_person_fit_plan.json").write_text(
        args.plan.read_text(encoding="utf-8"), encoding="utf-8"
    )

    r_script = ROOT / "validation/cmle_wle_person_fit_sirt.R"
    completed = subprocess.run(
        [
            "Rscript",
            str(r_script),
            "--input",
            str(output),
            "--output",
            str(output),
            "--repo",
            str(ROOT),
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    (output / "sirt_adapter_stdout.txt").write_text(
        completed.stdout, encoding="utf-8"
    )
    (output / "sirt_adapter_stderr.txt").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "sirt Person-fit adapter failed: "
            + (completed.stderr.strip() or completed.stdout.strip())
        )

    sirt_persons = pd.read_csv(output / "sirt_person_fit.csv")
    sirt_observations = pd.read_csv(output / "sirt_observation_moments.csv")
    reference_identity = pd.read_csv(output / "reference_identity.csv")
    person_pairs = python_persons.merge(
        sirt_persons,
        on=["Model", "Person"],
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    person_pairs["BothPresent"] = person_pairs["_merge"].eq("both")
    person_pairs["AbsInfitDifference"] = (
        person_pairs["Infit"] - person_pairs["SirtInfit"]
    ).abs()
    person_pairs["AbsOutfitDifference"] = (
        person_pairs["Outfit"] - person_pairs["SirtOutfit"]
    ).abs()
    infit_tolerance = float(
        plan["validation_matrix"]["python_sirt_max_abs_infit_difference"]
    )
    outfit_tolerance = float(
        plan["validation_matrix"]["python_sirt_max_abs_outfit_difference"]
    )
    person_pairs["InfitPassed"] = person_pairs["AbsInfitDifference"].le(
        infit_tolerance
    )
    person_pairs["OutfitPassed"] = person_pairs["AbsOutfitDifference"].le(
        outfit_tolerance
    )
    person_pairs["PairPassed"] = (
        person_pairs["BothPresent"]
        & person_pairs["PersonFitReady"].fillna(False)
        & person_pairs["InfitPassed"]
        & person_pairs["OutfitPassed"]
    )

    observation_keys = [
        "Model",
        "Person",
        "VirtualUnit",
        "ObservedInternalCategory",
    ]
    observation_pairs = python_observations.merge(
        sirt_observations,
        on=observation_keys,
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    moment_pairs = {
        "Expected": ("ExpectedInternalCategory", "SirtExpectedInternalCategory"),
        "Variance": ("Variance", "SirtVariance"),
        "Fourth": ("FourthCentralMoment", "SirtFourthCentralMoment"),
        "Residual": ("Residual", "SirtResidual"),
        "StdSq": (
            "StandardizedSquaredResidual",
            "SirtStandardizedSquaredResidual",
        ),
        "ObservedProbability": ("ObservedProbability", "SirtObservedProbability"),
    }
    for label, (python_column, r_column) in moment_pairs.items():
        observation_pairs[f"Abs{label}Difference"] = (
            observation_pairs[python_column] - observation_pairs[r_column]
        ).abs()
    expected_tolerance = float(
        plan["validation_matrix"]["python_sirt_max_abs_expected_difference"]
    )
    variance_tolerance = float(
        plan["validation_matrix"]["python_sirt_max_abs_variance_difference"]
    )
    observation_pairs["PairPassed"] = (
        observation_pairs["_merge"].eq("both")
        & observation_pairs["AbsExpectedDifference"].le(expected_tolerance)
        & observation_pairs["AbsVarianceDifference"].le(variance_tolerance)
    )

    comparison_rows: list[dict[str, object]] = []
    for model in plan["validation_matrix"]["models"]:
        p = person_pairs.loc[person_pairs["Model"].eq(str(model))]
        o = observation_pairs.loc[observation_pairs["Model"].eq(str(model))]
        comparison_rows.append(
            {
                "Model": str(model),
                "PersonPairs": int(len(p)),
                "ObservationPairs": int(len(o)),
                "MaxAbsInfitDifference": float(p["AbsInfitDifference"].max()),
                "MaxAbsOutfitDifference": float(p["AbsOutfitDifference"].max()),
                "MaxAbsExpectedDifference": float(o["AbsExpectedDifference"].max()),
                "MaxAbsVarianceDifference": float(o["AbsVarianceDifference"].max()),
                "MaxAbsFourthMomentDifference": float(o["AbsFourthDifference"].max()),
                "MaxAbsStandardizedSquaredResidualDifference": float(
                    o["AbsStdSqDifference"].max()
                ),
                "AllPersonPairsPassed": bool(p["PairPassed"].all()),
                "AllObservationPairsPassed": bool(o["PairPassed"].all()),
            }
        )
    comparison = pd.DataFrame(comparison_rows)
    reference_passed = bool(reference_identity["IdentityPassed"].all())
    decision_audit = pd.concat(
        [
            compute_cmle_wle_person_fit(fit_fixture(fixture, str(model)))[
                "decision_audit"
            ].assign(Model=str(model))
            for model in plan["validation_matrix"]["models"]
        ],
        ignore_index=True,
    )
    finite_audit = decision_audit.loc[decision_audit["BoundaryStatus"].ne("unavailable")]
    display_mismatches = int((~finite_audit["DisplayDecisionConsistent"].astype(bool)).sum())
    boundary_statistics = int(
        finite_audit["BoundaryStatus"].isin(
            ["numerical_boundary", "display_rounding_boundary"]
        ).sum()
    )
    all_passed = bool(
        reference_passed
        and person_pairs["PairPassed"].all()
        and observation_pairs["PairPassed"].all()
        and len(person_pairs)
        == int(sum(row["Persons"] for row in model_rows))
        and len(observation_pairs)
        == int(sum(row["Observations"] for row in model_rows))
    )

    write_csv(person_pairs.drop(columns="_merge"), output / "python_sirt_person_pairs.csv")
    write_csv(
        observation_pairs.drop(columns="_merge"),
        output / "python_sirt_observation_pairs.csv",
    )
    write_csv(comparison, output / "python_sirt_fit_summary.csv")
    write_csv(decision_audit, output / "person_fit_decision_stability_audit.csv")
    decision = {
        "schema_version": "mfrm-cmle-wle-person-fit-parity-v1",
        "overall_passed": all_passed,
        "status": "research_parity_ready_bootstrap_extension_withheld"
        if all_passed
        else "parity_failed",
        "models": [str(value) for value in plan["validation_matrix"]["models"]],
        "persons_compared": int(len(person_pairs)),
        "observations_compared": int(len(observation_pairs)),
        "exact_extreme_persons_compared": int(
            python_persons["ExtremeScorePattern"].sum()
        ),
        "reference_identity_passed": reference_passed,
        "maximum_absolute_infit_difference": float(
            person_pairs["AbsInfitDifference"].max()
        ),
        "maximum_absolute_outfit_difference": float(
            person_pairs["AbsOutfitDifference"].max()
        ),
        "maximum_absolute_expected_difference": float(
            observation_pairs["AbsExpectedDifference"].max()
        ),
        "maximum_absolute_variance_difference": float(
            observation_pairs["AbsVarianceDifference"].max()
        ),
        "maximum_absolute_fourth_moment_difference": float(
            observation_pairs["AbsFourthDifference"].max()
        ),
        "maximum_absolute_standardized_squared_residual_difference": float(
            observation_pairs["AbsStdSqDifference"].max()
        ),
        "finite_fit_statistics_audited": int(len(finite_audit)),
        "fit_boundary_statistics": boundary_statistics,
        "raw_display_decision_mismatches": display_mismatches,
        "tolerances": {
            "infit": infit_tolerance,
            "outfit": outfit_tolerance,
            "expected": expected_tolerance,
            "variance": variance_tolerance,
            "decision_basis": "unrounded retained floating-point values",
        },
        "TAM_formula_audit": {
            "same_untrimmed_infit_formula": True,
            "same_pre_trim_outfit_contributions": True,
            "default_outfit_estimand_identical": False,
            "reason": "TAM::tam.jml.fit trims unusually large squared standardized residuals before Outfit aggregation.",
        },
        "withheld": {
            "ZSTD": True,
            "p_value": True,
            "confidence_interval": True,
            "total_standard_error": True,
            "bootstrap_fit_zone_transitions": True,
            "streamlit_integration": True,
        },
        "registered_plan_sha256": sha256_file(args.plan),
        "person_fit_core_sha256": sha256_file(ROOT / "mfrm_app/cmle_wle_fit.py"),
        "python_runner_sha256": sha256_file(Path(__file__)),
        "r_adapter_sha256": sha256_file(r_script),
        "source_fixture_sha256": sha256_file(args.fixture),
    }
    (output / "person_fit_parity_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    status = "PASS" if all_passed else "FAIL"
    report = f"""# CMLE-WLE Person MnSq: Python–sirt validation

## Decision

**{status}.** Python fixed-calibration Person Infit/Outfit matches `sirt::pcm.fit` at the same exact-CMLE category kernels and the same Warm WLE Person measures. Decisions use retained unrounded values. This is research parity evidence, not authorization for bootstrap fit transitions or Streamlit output.

{markdown_table(comparison)}

## Scope and numerical result

- Compared {len(person_pairs)} Person rows and {len(observation_pairs)} administered observation rows across RSM and PCM.
- Included {int(python_persons['ExtremeScorePattern'].sum())} exact minimum/maximum Person patterns and non-rectangular missingness.
- Maximum absolute Infit difference: {person_pairs['AbsInfitDifference'].max():.12g} (gate ≤ {infit_tolerance:.1g}).
- Maximum absolute Outfit difference: {person_pairs['AbsOutfitDifference'].max():.12g} (gate ≤ {outfit_tolerance:.1g}).
- Maximum expected-category difference: {observation_pairs['AbsExpectedDifference'].max():.12g}; variance difference: {observation_pairs['AbsVarianceDifference'].max():.12g}.
- Reference source identity passed for TAM 4.3.25 and sirt 4.2.133: {reference_passed}.
- Unrounded/display audit: {len(finite_audit)} finite MnSq values, {boundary_statistics} boundary-sensitive values, {display_mismatches} raw/display classification mismatches in this fixture.

## Critical interpretation

`sirt::pcm.fit` is the direct parity target because it evaluates untrimmed Person MnSq at fixed theta. TAM uses the same Person Infit formula and the same pre-trim squared-standardized contributions, but `TAM::tam.jml.fit` trims unusually large contributions before Outfit. Therefore the default TAM Outfit and this untrimmed CMLE-WLE Outfit must not be labelled universally identical.

ZSTD and p-values remain unavailable, rather than zero: their finite-sample reference distribution, effective degrees of freedom, trimming convention, and exact-extreme behavior have not been validated for this two-stage estimand. Confidence intervals, total standard errors, bootstrap fit-zone transitions, and public UI output also remain withheld.
"""
    (output / "CMLE_WLE_PERSON_FIT_RESULTS.md").write_text(report, encoding="utf-8")
    if not all_passed:
        raise SystemExit("CMLE-WLE Person-fit parity gates failed.")


if __name__ == "__main__":
    main()
