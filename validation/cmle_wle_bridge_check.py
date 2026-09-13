#!/usr/bin/env python3
"""Retain deterministic evidence for the repository-only CMLE-WLE bridge."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle import fit_cmle  # noqa: E402
from mfrm_app.cmle_person_scoring import score_cmle_persons_wle  # noqa: E402
from mfrm_app.person_scoring import score_fixed_calibration_persons  # noqa: E402


DEFAULT_PLAN = ROOT / "validation/cmle_wle_bridge_plan_20260810.json"
DEFAULT_OUTPUT = ROOT / "validation/cmle_wle_bridge_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fixture(include_extremes: bool) -> pd.DataFrame:
    scores = {
        "P1": [0, 1, 1, 2],
        "P2": [1, 0, 2, 1],
        "P3": [2, 1, 1, 0],
        "P4": [1, 2, 0, 1],
        "P5": [0, 2, 1, 1],
        "P6": [2, 0, 1, 1],
    }
    if include_extremes:
        scores.update({"P7": [0, 0, 0, 0], "P8": [2, 2, 2, 2]})
    units = [("R1", "C1"), ("R1", "C2"), ("R2", "C1"), ("R2", "C2")]
    return pd.DataFrame(
        [
            (person, rater, criterion, score)
            for person, values in scores.items()
            for (rater, criterion), score in zip(units, values)
        ],
        columns=["Person", "Rater", "Criterion", "Score"],
    )


def validate_plan(plan_path: Path) -> dict:
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != "mfrm-cmle-wle-bridge-plan-v1":
        raise ValueError("Unexpected CMLE-WLE bridge plan schema.")
    paths = {
        "cmle_core_sha256": ROOT / "mfrm_app/cmle.py",
        "fixed_calibration_person_scoring_sha256": ROOT / "mfrm_app/person_scoring.py",
        "wle_parity_decision_sha256": ROOT / "validation/fixed_calibration_wle_20260809/parity_decision.json",
        "wle_downstream_decision_sha256": ROOT / "validation/fixed_calibration_wle_downstream_20260809/wle_downstream_decision.json",
        "cmle_phase0_results_sha256": ROOT / "validation/cmle_phase0_20260809/RESULTS.md",
        "cmle_scaling_results_sha256": ROOT / "validation/cmle_scaling_20260809/RESULTS.md",
    }
    failed = [
        key for key, path in paths.items()
        if sha256_file(path) != str(plan["input_identity"][key])
    ]
    if failed:
        raise ValueError(f"CMLE-WLE bridge input identity failed: {failed}")
    return plan


def run(plan_path: Path, output: Path) -> None:
    plan_path = plan_path.resolve()
    output = output.resolve()
    plan = validate_plan(plan_path)
    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_cmle_wle_bridge_plan.json").write_bytes(plan_path.read_bytes())

    case_rows: list[dict[str, object]] = []
    person_parts: list[pd.DataFrame] = []
    for model, include_extremes in (("RSM", False), ("PCM", True)):
        fit = fit_cmle(
            fixture(include_extremes),
            person_col="Person",
            facet_cols=["Rater", "Criterion"],
            score_col="Score",
            rating_min=0,
            rating_max=2,
            model=model,
            step_facet="Criterion" if model == "PCM" else None,
            gtol=1e-8,
        )
        design = fit["design"]
        coefficients = fit["coefficients"].set_index("Parameter")["Estimate"]
        parameters = coefficients.reindex(design.parameter_names).to_numpy(dtype=float)
        intercepts = np.einsum("rkp,p->rk", design.row_design, parameters, optimize=True)
        slopes = np.tile(np.arange(design.n_categories, dtype=float), (len(design.data), 1))
        bridge = score_cmle_persons_wle(fit)
        direct = score_fixed_calibration_persons(
            design.data[design.person_col].astype(str),
            design.data[design.score_col].to_numpy(dtype=int) - design.rating_min,
            intercepts,
            slopes,
            person_levels=fit["person_status"]["Person"].astype(str).tolist(),
        )
        direct_map = direct.set_index("Person")
        bridge["DirectGenericEstimate"] = bridge["Person"].map(direct_map["Estimate"])
        bridge["DirectGenericSE"] = bridge["Person"].map(direct_map["StandardError"])
        bridge["EstimateBridgeMinusDirect"] = bridge["Estimate"] - bridge["DirectGenericEstimate"]
        bridge["SEBridgeMinusDirect"] = bridge["StandardError"] - bridge["DirectGenericSE"]
        bridge.insert(0, "Model", model)
        person_parts.append(bridge)

        surface_differences: list[float] = []
        seen: set[tuple[str, ...]] = set()
        surface = fit["surfaces"]
        for row_number, row in design.data.iterrows():
            key = tuple(str(row[facet]) for facet in design.facet_cols)
            if key in seen:
                continue
            seen.add(key)
            selected = np.ones(len(surface), dtype=bool)
            for facet, level in zip(design.facet_cols, key):
                selected &= surface[facet].astype(str).eq(level).to_numpy()
            values = surface.loc[selected].sort_values("InternalCategory")["LogKernelWithoutPerson"].to_numpy(dtype=float)
            surface_differences.extend((values - intercepts[row_number]).tolist())

        shuffled = copy.deepcopy(fit)
        shuffled["coefficients"] = shuffled["coefficients"].sample(
            frac=1.0, random_state=31
        ).reset_index(drop=True)
        shuffled_scores = score_cmle_persons_wle(shuffled)
        permutation_difference = float(
            np.max(np.abs(shuffled_scores["Estimate"] - bridge["Estimate"]))
        )
        missing_rejected = False
        missing = copy.deepcopy(fit)
        missing["coefficients"] = missing["coefficients"].iloc[:-1].copy()
        try:
            score_cmle_persons_wle(missing)
        except ValueError:
            missing_rejected = True
        extreme = bridge["ExtremeScorePattern"].astype(bool)
        case_rows.append(
            {
                "Model": model,
                "CalibrationInferenceReady": bool(fit["summary"].iloc[0]["InferenceReady"]),
                "Persons": len(bridge),
                "ExtremePersons": int(extreme.sum()),
                "AllWLEStatusOK": bool(bridge["Status"].eq("ok").all()),
                "AllExtremeWLEFinite": bool(np.isfinite(bridge.loc[extreme, "Estimate"]).all()),
                "MaxAbsBridgeDirectEstimateDifference": float(bridge["EstimateBridgeMinusDirect"].abs().max()),
                "MaxAbsBridgeDirectSEDifference": float(bridge["SEBridgeMinusDirect"].abs().max()),
                "MaxAbsSurfaceInterceptDifference": float(np.max(np.abs(surface_differences))),
                "MaxAdjustedScoreResidual": float(bridge["AdjustedScoreResidual"].max()),
                "CoefficientPermutationMaxAbsDifference": permutation_difference,
                "MissingCoefficientRejected": missing_rejected,
            }
        )

    cases = pd.DataFrame(case_rows)
    persons = pd.concat(person_parts, ignore_index=True)
    tolerance = float(plan["validation_contract"]["generic_bridge_equivalence_max_abs"])
    surface_tolerance = float(plan["validation_contract"]["cmle_surface_intercept_max_abs"])
    residual_tolerance = float(plan["validation_contract"]["all_person_adjusted_score_residuals_max_abs"])
    passed = bool(
        cases["CalibrationInferenceReady"].all()
        and cases["AllWLEStatusOK"].all()
        and cases["AllExtremeWLEFinite"].all()
        and cases["MaxAbsBridgeDirectEstimateDifference"].le(tolerance).all()
        and cases["MaxAbsBridgeDirectSEDifference"].le(tolerance).all()
        and cases["MaxAbsSurfaceInterceptDifference"].le(surface_tolerance).all()
        and cases["MaxAdjustedScoreResidual"].le(residual_tolerance).all()
        and cases["CoefficientPermutationMaxAbsDifference"].le(tolerance).all()
        and cases["MissingCoefficientRejected"].all()
    )
    cases.to_csv(output / "cmle_wle_bridge_cases.csv", index=False, float_format="%.17g")
    persons.to_csv(output / "cmle_wle_person_scores.csv", index=False, float_format="%.17g")
    decision = {
        "schema_version": "mfrm-cmle-wle-bridge-result-v1",
        "overall_passed": passed,
        "models": cases["Model"].tolist(),
        "persons": int(len(persons)),
        "extreme_persons": int(persons["ExtremeScorePattern"].sum()),
        "maximum_bridge_direct_estimate_difference": float(cases["MaxAbsBridgeDirectEstimateDifference"].max()),
        "maximum_bridge_direct_se_difference": float(cases["MaxAbsBridgeDirectSEDifference"].max()),
        "maximum_surface_intercept_difference": float(cases["MaxAbsSurfaceInterceptDifference"].max()),
        "maximum_adjusted_score_residual": float(cases["MaxAdjustedScoreResidual"].max()),
        "calibration_uncertainty_propagated": False,
        "streamlit_integration_authorized": False,
        "bridge_source_sha256": sha256_file(ROOT / "mfrm_app/cmle_person_scoring.py"),
        "validation_source_sha256": sha256_file(Path(__file__).resolve()),
    }
    (output / "cmle_wle_bridge_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    report = f"""# Exact CMLE structural calibration plus fixed-calibration Warm WLE

## Decision

**{'PASS' if passed else 'FAIL'}.** The repository-only bridge passed its frozen RSM/PCM numerical and identity gates. Streamlit integration remains withheld.

## Outcome

- Models: {', '.join(cases['Model'])}; Persons: {len(persons)}; exact-extreme Persons: {int(persons['ExtremeScorePattern'].sum())}.
- Maximum bridge-versus-direct generic theta difference: {cases['MaxAbsBridgeDirectEstimateDifference'].max():.12g}.
- Maximum bridge-versus-direct generic SE difference: {cases['MaxAbsBridgeDirectSEDifference'].max():.12g}.
- Maximum CMLE category-surface/intercept difference: {cases['MaxAbsSurfaceInterceptDifference'].max():.12g}.
- Maximum WLE adjusted-score residual: {cases['MaxAdjustedScoreResidual'].max():.12g}.
- Coefficient ordering was permutation-invariant; a missing coefficient was rejected in both model cases.

## Interpretation boundary

This validates a two-stage computational handoff: exact CMLE estimates structural RSM/PCM calibration, and Warm WLE scores fit-sample Persons with that calibration treated as fixed. WLE is not part of the conditional likelihood and is not a finite JMLE MLE. Its reported SE excludes CMLE calibration uncertainty. New-Person/unseen-unit scoring, anchors, downstream CMLE-WLE diagnostics, bootstrap uncertainty, GPCM, automatic fallback, and public UI claims remain out of scope.
"""
    (output / "CMLE_WLE_BRIDGE_RESULTS.md").write_text(report, encoding="utf-8")
    if not passed:
        raise SystemExit("CMLE-WLE bridge gates failed.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.plan, args.output)


if __name__ == "__main__":
    main()
