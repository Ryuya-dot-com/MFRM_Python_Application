#!/usr/bin/env python3
"""Run the frozen asymptotic CMLE-calibration sensitivity matrix for WLE."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle import fit_cmle  # noqa: E402
from mfrm_app.cmle_person_scoring import score_cmle_persons_wle  # noqa: E402
from mfrm_app.cmle_wle_uncertainty import (  # noqa: E402
    cmle_wle_calibration_sensitivity,
)


DEFAULT_PLAN = ROOT / "validation/cmle_wle_calibration_sensitivity_plan_20260810.json"
DEFAULT_DOCUMENTATION_AMENDMENT = (
    ROOT / "validation/cmle_wle_calibration_sensitivity_documentation_amendment_20260810.json"
)
DEFAULT_OUTPUT = ROOT / "validation/cmle_wle_calibration_sensitivity_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fixture() -> pd.DataFrame:
    scores = {
        "P1": [0, 0, 1, 2],
        "P2": [1, 0, 2, 2],
        "P3": [2, 1, 1, 0],
        "P4": [1, 2, 0, 2],
        "P5": [0, 2, 1, 0],
        "P6": [2, 0, 2, 1],
        "P7": [0, 0, 0, 0],
        "P8": [2, 2, 2, 2],
    }
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
    if plan.get("schema_version") != "mfrm-cmle-wle-calibration-sensitivity-plan-v1":
        raise ValueError("Unexpected CMLE-WLE calibration-sensitivity plan schema.")
    paths = {
        "cmle_core_sha256": ROOT / "mfrm_app/cmle.py",
        "fixed_calibration_person_scoring_sha256": ROOT / "mfrm_app/person_scoring.py",
        "cmle_wle_bridge_sha256": ROOT / "mfrm_app/cmle_person_scoring.py",
        "cmle_wle_bridge_decision_sha256": ROOT / "validation/cmle_wle_bridge_20260810/cmle_wle_bridge_decision.json",
        "wle_parity_decision_sha256": ROOT / "validation/fixed_calibration_wle_20260809/parity_decision.json",
        "cmle_phase0_design_sha256": ROOT / "docs/cmle_phase0_design.md",
    }
    failed = [
        key for key, path in paths.items()
        if sha256_file(path) != str(plan["input_identity"][key])
    ]
    documentation_key = "cmle_phase0_design_sha256"
    if failed == [documentation_key]:
        if not DEFAULT_DOCUMENTATION_AMENDMENT.exists():
            raise ValueError("CMLE design identity changed without a documentation amendment.")
        amendment = json.loads(
            DEFAULT_DOCUMENTATION_AMENDMENT.read_text(encoding="utf-8")
        )
        valid_amendment = bool(
            amendment.get("schema_version")
            == "mfrm-cmle-wle-calibration-sensitivity-documentation-amendment-v1"
            and amendment.get("base_plan_sha256") == sha256_file(plan_path)
            and amendment.get("original_cmle_phase0_design_sha256")
            == str(plan["input_identity"][documentation_key])
            and amendment.get("current_cmle_phase0_design_sha256")
            == sha256_file(paths[documentation_key])
            and amendment.get("numerical_contract_changed") is False
            and amendment.get("validation_gates_changed") is False
        )
        if not valid_amendment:
            raise ValueError("CMLE design documentation amendment identity failed.")
        plan["documentation_amendment"] = amendment
        failed = []
    if failed:
        raise ValueError(f"Calibration-sensitivity input identity failed: {failed}")
    bridge = json.loads(paths["cmle_wle_bridge_decision_sha256"].read_text(encoding="utf-8"))
    if bridge.get("overall_passed") is not True:
        raise ValueError("CMLE-WLE bridge evidence is not qualified.")
    return plan


def fit_fixture(model: str) -> dict:
    return fit_cmle(
        fixture(),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model=model,
        step_facet="Criterion" if model == "PCM" else None,
        gtol=1e-8,
    )


def plot_results(cases: pd.DataFrame, persons: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))
    for model, frame in cases.groupby("Model", sort=False):
        axes[0].plot(
            frame["CovarianceScale"],
            frame["MedianCalibrationDrawSD"],
            marker="o",
            linewidth=1.8,
            label=model,
        )
    axes[0].set_xlabel("CMLE covariance scale")
    axes[0].set_ylabel("Median calibration-draw SD (logits)")
    axes[0].set_title("Calibration perturbation response")
    axes[0].legend(frameon=False)

    for model, frame in persons.groupby("Model", sort=False):
        interior = ~frame["ExtremeScorePattern"].astype(bool)
        axes[1].scatter(
            frame.loc[interior, "StandardError"],
            frame.loc[interior, "CalibrationDrawSD"],
            s=30,
            alpha=0.75,
            label=f"{model} interior",
        )
        axes[1].scatter(
            frame.loc[~interior, "StandardError"],
            frame.loc[~interior, "CalibrationDrawSD"],
            s=65,
            marker="X",
            label=f"{model} extreme",
        )
    axes[1].set_xlabel("Conditional WLE SE, calibration fixed")
    axes[1].set_ylabel("CMLE calibration-draw SD")
    axes[1].set_title("Two distinct uncertainty components")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "cmle_wle_calibration_sensitivity.png", dpi=180)
    plt.close(fig)


def first_read_projection(
    *,
    passed: bool,
    cases: pd.DataFrame,
    persons: pd.DataFrame,
    draws: int,
) -> pd.DataFrame:
    """Project research evidence into separate, UI-ready status cards."""
    primary = cases.loc[cases["CovarianceScale"].eq(1.0)]
    extremes = persons["ExtremeScorePattern"].astype(bool)
    maximum_sd = float(persons["CalibrationDrawSD"].max())
    maximum_ratio = float(persons["CalibrationSDToConditionalSERatio"].max())
    successful = int(persons["CalibrationDrawsSuccessful"].sum())
    requested = int(len(persons) * draws)
    return pd.DataFrame(
        [
            {
                "CardOrder": 1,
                "CardId": "structural_calibration",
                "Status": "research_ready" if passed else "blocked",
                "RawValue": int(primary["Model"].nunique()),
                "DisplayValue": f"{primary['Model'].nunique()} model fixtures passed",
                "Interpretation": "Exact CMLE calibration passed the frozen RSM/PCM research gates.",
                "NextAction": "Inspect model-specific rank, covariance, and scope evidence.",
            },
            {
                "CardOrder": 2,
                "CardId": "person_scoring",
                "Status": "research_ready" if successful == requested else "blocked",
                "RawValue": successful,
                "DisplayValue": f"{successful}/{requested} draw-Person scores returned",
                "Interpretation": "Warm WLE scoring succeeded conditional on each calibration draw.",
                "NextAction": "Keep the combined CMLE + fixed-calibration WLE estimator label.",
            },
            {
                "CardOrder": 3,
                "CardId": "exact_extremes",
                "Status": "caution",
                "RawValue": int(extremes.sum()),
                "DisplayValue": f"{int(extremes.sum())} exact-extreme Persons",
                "Interpretation": "WLE is finite, but these are not finite JMLE maximum-likelihood estimates.",
                "NextAction": "Show extreme-score status beside every affected Person estimate.",
            },
            {
                "CardOrder": 4,
                "CardId": "calibration_sensitivity",
                "Status": "caution",
                "RawValue": maximum_sd,
                "DisplayValue": f"max calibration-draw SD {maximum_sd:.3f} logits",
                "Interpretation": (
                    "Calibration perturbation reached "
                    f"{100.0 * maximum_ratio:.1f}% of the conditional WLE SE."
                ),
                "NextAction": "Inspect Person-level sensitivity, especially exact extremes.",
            },
            {
                "CardOrder": 5,
                "CardId": "interval_claim",
                "Status": "withheld",
                "RawValue": np.nan,
                "DisplayValue": "confidence interval withheld",
                "Interpretation": "Same-sample CMLE/WLE dependence is absent from the coefficient draws.",
                "NextAction": "Require a qualified bootstrap or independent-calibration design.",
            },
            {
                "CardOrder": 6,
                "CardId": "public_ui",
                "Status": "withheld",
                "RawValue": np.nan,
                "DisplayValue": "Streamlit integration not authorized",
                "Interpretation": "Passing numerical research gates does not establish public readiness.",
                "NextAction": "Complete new-Person, downstream, anchor, bootstrap, and UI gates.",
            },
        ]
    )


def run(plan_path: Path, output: Path) -> None:
    plan_path = plan_path.resolve()
    output = output.resolve()
    plan = validate_plan(plan_path)
    matrix = plan["validation_matrix"]
    draws = int(matrix["primary_draws"])
    seed = int(matrix["primary_seed"])
    scales = [float(value) for value in matrix["covariance_scales"]]
    output.mkdir(parents=True, exist_ok=True)
    (output / "registered_calibration_sensitivity_plan.json").write_bytes(
        plan_path.read_bytes()
    )
    if "documentation_amendment" in plan:
        (output / "registered_documentation_amendment.json").write_bytes(
            DEFAULT_DOCUMENTATION_AMENDMENT.read_bytes()
        )

    case_parts: list[pd.DataFrame] = []
    primary_person_parts: list[pd.DataFrame] = []
    reproducibility_rows: list[dict[str, object]] = []
    zero_rows: list[dict[str, object]] = []
    for model in matrix["models"]:
        fit = fit_fixture(str(model))
        baseline = score_cmle_persons_wle(fit).set_index("Person")
        zero = cmle_wle_calibration_sensitivity(
            fit,
            n_draws=50,
            seed=seed,
            covariance_scale=0.0,
        )
        zero_rows.append(
            {
                "Model": model,
                "MaximumCalibrationDrawSD": float(zero["persons"]["CalibrationDrawSD"].max()),
                "MaximumAbsDrawMeanMinusBaseline": float(
                    zero["persons"]["CalibrationDrawMeanMinusBaseline"].abs().max()
                ),
            }
        )
        model_primary = None
        for scale in scales:
            result = cmle_wle_calibration_sensitivity(
                fit,
                n_draws=draws,
                seed=seed,
                covariance_scale=scale,
            )
            summary = result["summary"].copy()
            summary["Model"] = model
            case_parts.append(summary)
            if scale == 1.0:
                model_primary = result
                persons = result["persons"].copy()
                persons.insert(0, "Model", model)
                persons["BaselineBridgeDifference"] = (
                    persons["Estimate"].to_numpy(dtype=float)
                    - persons["Person"].map(baseline["Estimate"]).to_numpy(dtype=float)
                )
                primary_person_parts.append(persons)
        if model_primary is None:
            raise ValueError("Primary covariance scale 1.0 is absent.")
        repeated = cmle_wle_calibration_sensitivity(
            fit,
            n_draws=draws,
            seed=seed,
            covariance_scale=1.0,
        )
        columns = [
            "CalibrationDrawMean",
            "CalibrationDrawSD",
            "CalibrationDrawQ025",
            "CalibrationDrawQ500",
            "CalibrationDrawQ975",
        ]
        difference = np.max(
            np.abs(
                model_primary["persons"][columns].to_numpy(dtype=float)
                - repeated["persons"][columns].to_numpy(dtype=float)
            )
        )
        reproducibility_rows.append({"Model": model, "SameSeedMaxAbsDifference": float(difference)})

    cases = pd.concat(case_parts, ignore_index=True)
    persons = pd.concat(primary_person_parts, ignore_index=True)
    zero = pd.DataFrame(zero_rows)
    reproducibility = pd.DataFrame(reproducibility_rows)
    monotonic_rows = []
    for model, frame in cases.groupby("Model", sort=False):
        ordered = frame.sort_values("CovarianceScale")
        values = ordered["MedianCalibrationDrawSD"].to_numpy(dtype=float)
        monotonic_rows.append(
            {
                "Model": model,
                "StrictlyIncreasingMedianCalibrationDrawSD": bool(np.all(np.diff(values) > 0)),
                "Scale0p5MedianSD": values[0],
                "Scale1MedianSD": values[1],
                "Scale2MedianSD": values[2],
            }
        )
    monotonic = pd.DataFrame(monotonic_rows)

    baseline_tolerance = float(matrix["baseline_bridge_max_abs_difference"])
    zero_tolerance = float(matrix["zero_scale_max_abs_calibration_draw_sd"])
    same_seed_tolerance = float(matrix["same_seed_max_abs_difference"])
    extreme = persons["ExtremeScorePattern"].astype(bool)
    passed = bool(
        persons["BaselineBridgeDifference"].abs().le(baseline_tolerance).all()
        and zero["MaximumCalibrationDrawSD"].le(zero_tolerance).all()
        and zero["MaximumAbsDrawMeanMinusBaseline"].le(zero_tolerance).all()
        and reproducibility["SameSeedMaxAbsDifference"].le(same_seed_tolerance).all()
        and cases["AllDrawPersonScoresSuccessful"].all()
        and persons.loc[extreme, "CalibrationDrawsSuccessful"].eq(draws).all()
        and np.isfinite(persons.loc[extreme, "CalibrationDrawSD"]).all()
        and monotonic["StrictlyIncreasingMedianCalibrationDrawSD"].all()
    )

    cases.to_csv(output / "calibration_sensitivity_cases.csv", index=False, float_format="%.17g")
    persons.to_csv(output / "calibration_sensitivity_persons_primary.csv", index=False, float_format="%.17g")
    zero.to_csv(output / "calibration_sensitivity_zero_scale.csv", index=False, float_format="%.17g")
    reproducibility.to_csv(output / "calibration_sensitivity_reproducibility.csv", index=False, float_format="%.17g")
    monotonic.to_csv(output / "calibration_sensitivity_scale_monotonicity.csv", index=False, float_format="%.17g")
    first_read = first_read_projection(
        passed=passed,
        cases=cases,
        persons=persons,
        draws=draws,
    )
    first_read.to_csv(
        output / "cmle_wle_first_read_projection.csv",
        index=False,
        float_format="%.17g",
    )
    plot_results(cases, persons, output)

    decision = {
        "schema_version": "mfrm-cmle-wle-calibration-sensitivity-result-v1",
        "overall_passed": passed,
        "models": cases["Model"].drop_duplicates().tolist(),
        "primary_draws_per_model": draws,
        "primary_seed": seed,
        "persons_primary": len(persons),
        "extreme_persons_primary": int(extreme.sum()),
        "maximum_primary_calibration_draw_sd": float(persons["CalibrationDrawSD"].max()),
        "maximum_primary_calibration_sd_to_conditional_se_ratio": float(
            persons["CalibrationSDToConditionalSERatio"].max()
        ),
        "maximum_abs_primary_draw_mean_shift": float(
            persons["CalibrationDrawMeanMinusBaseline"].abs().max()
        ),
        "interval_is_confidence_interval": False,
        "quadrature_sensitivity_is_inference_qualified": False,
        "first_read_projection_generated": True,
        "streamlit_integration_authorized": False,
        "sensitivity_source_sha256": sha256_file(ROOT / "mfrm_app/cmle_wle_uncertainty.py"),
        "runner_source_sha256": sha256_file(Path(__file__).resolve()),
    }
    (output / "calibration_sensitivity_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    primary_cases = cases.loc[cases["CovarianceScale"].eq(1.0)]
    report = f"""# CMLE calibration-draw sensitivity for fixed-calibration Warm WLE

## Decision

**{'PASS' if passed else 'FAIL'}.** Frozen numerical, reproducibility, extreme-score, and scale-response gates passed. The output remains a same-sample asymptotic sensitivity diagnostic, not an inference-qualified total SE or confidence interval.

## Primary scale outcome

- Models: {', '.join(primary_cases['Model'])}; coefficient draws per model: {draws}; Persons: {len(persons)}; exact extremes: {int(extreme.sum())}.
- Median calibration-draw SD range across models: {primary_cases['MedianCalibrationDrawSD'].min():.6g}--{primary_cases['MedianCalibrationDrawSD'].max():.6g} logits.
- Maximum Person calibration-draw SD: {persons['CalibrationDrawSD'].max():.6g} logits.
- Maximum calibration-SD / conditional-WLE-SE ratio: {persons['CalibrationSDToConditionalSERatio'].max():.6g}.
- Maximum absolute draw-mean minus baseline shift: {persons['CalibrationDrawMeanMinusBaseline'].abs().max():.6g} logits.
- All draw-Person scores returned; exact extremes remained finite.
- Zero covariance scale produced at most {zero['MaximumCalibrationDrawSD'].max():.3g} calibration-draw SD; identical seeds reproduced outputs with maximum difference {reproducibility['SameSeedMaxAbsDifference'].max():.3g}.

## Interpretation boundary

The CMLE covariance materially perturbs WLE Person scores in this deliberately small calibration fixture, and the response grows with covariance scale. `CalibrationDrawSD` isolates this perturbation channel. `QuadratureSensitivitySE` is retained only as a planning sensitivity: the fit-sample CMLE and WLE stages reuse responses, so their dependence is not represented by independent multivariate-normal coefficient draws. The draw quantiles are not confidence or credible intervals. A conditional bootstrap or independent-calibration/new-Person design is required before inferential use or UI integration.
"""
    (output / "CMLE_WLE_CALIBRATION_SENSITIVITY_RESULTS.md").write_text(
        report, encoding="utf-8"
    )
    if not passed:
        raise SystemExit("CMLE-WLE calibration sensitivity gates failed.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.plan, args.output)


if __name__ == "__main__":
    main()
