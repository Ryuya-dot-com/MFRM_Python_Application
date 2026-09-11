#!/usr/bin/env python3
"""Compare frozen Python fixed-calibration WLE output with TAM."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
THETA_TOLERANCE = 5e-7
SE_TOLERANCE = 5e-7
RESIDUAL_TOLERANCE = 1e-8


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def markdown_table(frame: pd.DataFrame) -> str:
    rendered = frame.copy()
    for column in rendered.columns:
        if pd.api.types.is_float_dtype(rendered[column]):
            rendered[column] = rendered[column].map(
                lambda value: "" if pd.isna(value) else f"{value:.9g}"
            )
        elif pd.api.types.is_bool_dtype(rendered[column]):
            rendered[column] = rendered[column].map(lambda value: "TRUE" if value else "FALSE")
    header = "| " + " | ".join(rendered.columns.astype(str)) + " |"
    separator = "| " + " | ".join(["---"] * len(rendered.columns)) + " |"
    rows = [
        "| " + " | ".join(str(value) for value in row) + " |"
        for row in rendered.itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *rows])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    manifest_path = args.input / "fixture_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for entry in manifest["files"].values():
        path = ROOT / entry["path"]
        if not path.exists() or sha256_file(path) != entry["sha256"]:
            raise RuntimeError(f"Frozen fixture input changed: {entry['path']}")

    python = pd.read_csv(args.input / "python_wle.csv")
    tam = pd.read_csv(args.input / "tam_wle.csv")
    identity = pd.read_csv(args.input / "tam_identity.csv")
    pairs = python.merge(
        tam,
        on=["Case", "Person"],
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    pairs["BothPresent"] = pairs["_merge"].eq("both")
    pairs["ThetaDifferencePythonMinusTAM"] = pairs["Estimate"] - pairs["EstimateTAM"]
    pairs["AbsThetaDifference"] = pairs["ThetaDifferencePythonMinusTAM"].abs()
    pairs["SEDifferencePythonMinusTAM"] = (
        pairs["StandardError"] - pairs["StandardErrorTAM"]
    )
    pairs["AbsSEDifference"] = pairs["SEDifferencePythonMinusTAM"].abs()
    pairs["ObservedRowsMatch"] = pairs["ObservedRows"].eq(pairs["ObservedRowsTAM"])
    pairs["ThetaPassed"] = pairs["AbsThetaDifference"].le(THETA_TOLERANCE)
    pairs["SEPassed"] = pairs["AbsSEDifference"].le(SE_TOLERANCE)
    pairs["ResidualPassed"] = pairs["AdjustedScoreResidual"].le(RESIDUAL_TOLERANCE)
    pairs["PairPassed"] = (
        pairs["BothPresent"]
        & pairs["Status"].eq("ok")
        & pairs["ObservedRowsMatch"]
        & pairs["ThetaPassed"]
        & pairs["SEPassed"]
        & pairs["ResidualPassed"]
    )

    summary_rows: list[dict[str, object]] = []
    for case_name in manifest["case_names"]:
        case = pairs.loc[pairs["Case"].eq(case_name)]
        summary_rows.append(
            {
                "Case": case_name,
                "Pairs": len(case),
                "ExtremePatterns": int(case["ExtremeScorePattern"].fillna(False).sum()),
                "MaxAbsThetaDifference": float(case["AbsThetaDifference"].max()),
                "MaxAbsSEDifference": float(case["AbsSEDifference"].max()),
                "MaxAdjustedScoreResidual": float(case["AdjustedScoreResidual"].max()),
                "AllPassed": bool(case["PairPassed"].all()),
            }
        )
    summary = pd.DataFrame(summary_rows)

    identity_passed = (
        len(identity) == 1
        and bool(identity.loc[0, "FrozenAmendmentIdentityPassed"])
        and str(identity.loc[0, "Version"]) == "4.3.25"
    )
    case_order_passed = summary["Case"].tolist() == manifest["case_names"]
    all_extremes_finite = bool(
        pairs.loc[pairs["ExtremeScorePattern"].fillna(False), ["Estimate", "EstimateTAM"]]
        .apply(np.isfinite)
        .all()
        .all()
    )
    overall_passed = bool(
        identity_passed
        and case_order_passed
        and pairs["PairPassed"].all()
        and all_extremes_finite
        and len(pairs) == int(manifest["persons"])
    )

    pairs.drop(columns="_merge").to_csv(
        args.output / "python_tam_wle_pairs.csv",
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )
    summary.to_csv(
        args.output / "python_tam_wle_summary.csv",
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )
    decision = {
        "schema_version": "mfrm-fixed-calibration-wle-parity-v1",
        "overall_passed": overall_passed,
        "persons_compared": int(len(pairs)),
        "cases_compared": int(len(summary)),
        "identity_passed_against_pre_result_amendment": identity_passed,
        "original_plan_hash_reproducible": False,
        "all_extreme_pattern_wles_finite": all_extremes_finite,
        "maximum_absolute_theta_difference": float(pairs["AbsThetaDifference"].max()),
        "maximum_absolute_standard_error_difference": float(pairs["AbsSEDifference"].max()),
        "maximum_python_adjusted_score_residual": float(
            pairs["AdjustedScoreResidual"].max()
        ),
        "thresholds": {
            "theta": THETA_TOLERANCE,
            "standard_error": SE_TOLERANCE,
            "python_adjusted_score_residual": RESIDUAL_TOLERANCE,
            "basis": "unrounded retained floating-point values",
        },
        "scope": {
            "fixed_calibration_only": True,
            "streamlit_integration": False,
            "jmle_replacement": False,
            "facet_calibration_uncertainty": False,
            "downstream_fit_bias_visualization_rebuilt": False,
        },
    }
    (args.output / "parity_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    status = "PASS" if overall_passed else "FAIL"
    report = f"""# Fixed-calibration Warm WLE: Python–TAM parity

## Decision

**{status}.** The decision uses unrounded retained values. Python and TAM are compared as the same one-dimensional fixed-calibration Warm WLE estimand; this is not JMLE/MML parity and does not authorize application integration.

{markdown_table(summary)}

## Frozen gates

- Persons compared: {len(pairs)} across {len(summary)} cases.
- Maximum absolute theta difference: {pairs['AbsThetaDifference'].max():.12g} (gate ≤ {THETA_TOLERANCE:.1g}).
- Maximum absolute standard-error difference: {pairs['AbsSEDifference'].max():.12g} (gate ≤ {SE_TOLERANCE:.1g}).
- Maximum Python adjusted-score residual: {pairs['AdjustedScoreResidual'].max():.12g} (gate ≤ {RESIDUAL_TOLERANCE:.1g}).
- All exact extreme patterns returned finite WLEs in both engines: {all_extremes_finite}.
- Loaded TAM identity matched the pre-result amendment: {identity_passed}.

## Identity audit

The original plan's TAM function hash was not reproducible. It was not overwritten. Before fixtures or TAM results were created, `fixed_calibration_wle_plan_amendment_20260809.json` froze the explicit deparse calculation and a secondary formals/body digest; the loaded TAM 4.3.25 code matched both amended values.

## Interpretation boundary

This evidence validates the numerical translation of the fixed-calibration Warm adjusted score and its conditional information standard error for the frozen RSM/PCM/GPCM fixtures, including missing responses and exact minimum/maximum patterns. It does not propagate calibration uncertainty, estimate facets jointly, qualify Person reliability, replace finite JMLE MLEs, or establish downstream fit/bias/visualization behavior. Streamlit behavior remains unchanged.
"""
    (args.output / "FIXED_CALIBRATION_WLE_RESULTS.md").write_text(
        report,
        encoding="utf-8",
    )
    if not overall_passed:
        raise SystemExit("Fixed-calibration Python–TAM WLE parity gates failed.")


if __name__ == "__main__":
    main()
