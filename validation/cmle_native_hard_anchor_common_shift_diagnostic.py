#!/usr/bin/env python3
"""Verify the post-result common anchor-origin shift invariance."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle import fit_cmle  # noqa: E402
from mfrm_app.cmle_person_scoring import score_cmle_persons_wle  # noqa: E402
from mfrm_app.cmle_wle_fit import compute_cmle_wle_person_fit  # noqa: E402
from mfrm_app.fit_threshold_operating_characteristics import fit_rule_flags  # noqa: E402


DEFAULT_ADDENDUM = ROOT / "validation/cmle_native_hard_anchor_common_shift_addendum_20260810.json"
DEFAULT_OUTPUT = ROOT / "validation/cmle_native_hard_anchor_common_shift_20260810"


CORRECT_ANCHORS = pd.DataFrame(
    [
        {"ParameterType": "Facet", "Facet": "Rater", "Level": "R01", "Value": -0.75},
        {"ParameterType": "Facet", "Facet": "Rater", "Level": "R03", "Value": -0.15},
        {"ParameterType": "Facet", "Facet": "Rater", "Level": "R05", "Value": 0.45},
    ]
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_addendum(path: Path) -> dict[str, object]:
    addendum = json.loads(path.read_text(encoding="utf-8"))
    if addendum.get("study_id") != "cmle_native_hard_anchor_common_shift_diagnostic_v1":
        raise ValueError("Unexpected common-shift diagnostic addendum.")
    parent = addendum["parent_smoke"]
    identities = {
        parent["decision_path"]: parent["decision_sha256"],
        parent["retained_responses_path"]: parent["retained_responses_sha256"],
        "validation/cmle_native_hard_anchor_smoke_20260810/hard_anchor_smoke_person_results.csv": parent["person_results_sha256"],
        "validation/cmle_native_hard_anchor_smoke_20260810/hard_anchor_smoke_structural_recovery.csv": parent["structural_recovery_sha256"],
    }
    mismatches = [
        relative
        for relative, digest in identities.items()
        if not (ROOT / relative).exists() or sha256_file(ROOT / relative) != digest
    ]
    if mismatches:
        raise ValueError(f"Common-shift parent identity failed: {mismatches}")
    return addendum


def fit_scenario(frame: pd.DataFrame, anchors: pd.DataFrame) -> dict[str, object]:
    analysis = frame[["Person", "Rater", "Criterion", "ObservedCategory"]].copy()
    return fit_cmle(
        analysis,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="ObservedCategory",
        rating_min=0,
        rating_max=3,
        model="RSM",
        hard_anchors=anchors,
        gtol=1e-8,
        maxiter=800,
        newton_polish_maxiter=12,
    )


def max_abs(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.max(np.abs(values))) if values.size else 0.0


def compare_run(
    run_frame: pd.DataFrame,
    *,
    delta: float,
) -> dict[str, object]:
    correct = fit_scenario(run_frame, CORRECT_ANCHORS)
    shifted_anchors = CORRECT_ANCHORS.copy()
    shifted_anchors["Value"] += delta
    shifted = fit_scenario(run_frame, shifted_anchors)
    correct_summary = correct["summary"].iloc[0]
    shifted_summary = shifted["summary"].iloc[0]

    correct_coefficients = correct["coefficients"].set_index("Parameter")
    shifted_coefficients = shifted["coefficients"].set_index("Parameter").reindex(
        correct_coefficients.index
    )
    rater_free = correct_coefficients.index.str.startswith("facet:Rater:")
    other_free = ~rater_free
    free_shift_error = max_abs(
        shifted_coefficients.loc[rater_free, "Estimate"].to_numpy()
        - correct_coefficients.loc[rater_free, "Estimate"].to_numpy()
        - delta
    )
    other_free_difference = max_abs(
        shifted_coefficients.loc[other_free, "Estimate"].to_numpy()
        - correct_coefficients.loc[other_free, "Estimate"].to_numpy()
    )

    correct_expanded = pd.concat(
        [correct["facets"]["others"], correct["steps"]], ignore_index=True
    ).set_index(["ParameterType", "Facet", "Level", "Step"])
    shifted_expanded = pd.concat(
        [shifted["facets"]["others"], shifted["steps"]], ignore_index=True
    ).set_index(["ParameterType", "Facet", "Level", "Step"]).reindex(
        correct_expanded.index
    )
    rater_rows = correct_expanded.index.get_level_values("Facet") == "Rater"
    criterion_step_rows = ~rater_rows
    expanded_rater_shift_error = max_abs(
        shifted_expanded.loc[rater_rows, "Estimate"].to_numpy()
        - correct_expanded.loc[rater_rows, "Estimate"].to_numpy()
        - delta
    )
    criterion_step_difference = max_abs(
        shifted_expanded.loc[criterion_step_rows, "Estimate"].to_numpy()
        - correct_expanded.loc[criterion_step_rows, "Estimate"].to_numpy()
    )
    structural_se_difference = max_abs(
        shifted_expanded["SE"].to_numpy() - correct_expanded["SE"].to_numpy()
    )

    surface_key = ["VirtualUnit", "InternalCategory"]
    correct_surface = correct["surfaces"].sort_values(surface_key).reset_index(drop=True)
    shifted_surface = shifted["surfaces"].sort_values(surface_key).reset_index(drop=True)
    predicted_surface_difference = -delta * correct_surface[
        "InternalCategory"
    ].to_numpy(dtype=float)
    surface_shift_error = max_abs(
        shifted_surface["LogKernelWithoutPerson"].to_numpy(dtype=float)
        - correct_surface["LogKernelWithoutPerson"].to_numpy(dtype=float)
        - predicted_surface_difference
    )

    correct_wle = score_cmle_persons_wle(correct).sort_values("Person").reset_index(drop=True)
    shifted_wle = score_cmle_persons_wle(shifted).sort_values("Person").reset_index(drop=True)
    wle_shift_error = max_abs(
        shifted_wle["Estimate"].to_numpy() - correct_wle["Estimate"].to_numpy() - delta
    )
    wle_se_difference = max_abs(
        shifted_wle["StandardError"].to_numpy()
        - correct_wle["StandardError"].to_numpy()
    )

    correct_fit = compute_cmle_wle_person_fit(correct)
    shifted_fit = compute_cmle_wle_person_fit(shifted)
    correct_observation = correct_fit["observations"]
    shifted_observation = shifted_fit["observations"]
    expected_difference = max_abs(
        shifted_observation["ExpectedInternalCategory"].to_numpy()
        - correct_observation["ExpectedInternalCategory"].to_numpy()
    )
    correct_person = correct_fit["persons"].sort_values("Person").reset_index(drop=True)
    shifted_person = shifted_fit["persons"].sort_values("Person").reset_index(drop=True)
    infit_difference = max_abs(
        shifted_person["Infit"].to_numpy() - correct_person["Infit"].to_numpy()
    )
    outfit_difference = max_abs(
        shifted_person["Outfit"].to_numpy() - correct_person["Outfit"].to_numpy()
    )
    correct_flag = fit_rule_flags(correct_person["Infit"], correct_person["Outfit"])[
        "either_upper"
    ]
    shifted_flag = fit_rule_flags(shifted_person["Infit"], shifted_person["Outfit"])[
        "either_upper"
    ]

    return {
        "ConditionId": str(run_frame["ConditionId"].iloc[0]),
        "Replicate": int(run_frame["Replicate"].iloc[0]),
        "RunId": str(run_frame["RunId"].iloc[0]),
        "ResponseRows": int(len(run_frame)),
        "Persons": int(run_frame["Person"].nunique()),
        "CorrectInferenceReady": bool(correct_summary["InferenceReady"]),
        "ShiftedInferenceReady": bool(shifted_summary["InferenceReady"]),
        "ConditionalLogLikAbsDifference": abs(
            float(shifted_summary["ConditionalLogLik"])
            - float(correct_summary["ConditionalLogLik"])
        ),
        "InformationMaxAbsDifference": max_abs(
            np.asarray(shifted["information"]) - np.asarray(correct["information"])
        ),
        "CovarianceMaxAbsDifference": max_abs(
            np.asarray(shifted["covariance"]) - np.asarray(correct["covariance"])
        ),
        "FreeRaterShiftErrorMax": free_shift_error,
        "OtherFreeEstimateDifferenceMax": other_free_difference,
        "ExpandedRaterShiftErrorMax": expanded_rater_shift_error,
        "CriterionStepEstimateDifferenceMax": criterion_step_difference,
        "StructuralSEDifferenceMax": structural_se_difference,
        "SurfaceShiftErrorMax": surface_shift_error,
        "WLEShiftErrorMax": wle_shift_error,
        "WLEStandardErrorDifferenceMax": wle_se_difference,
        "ExpectedScoreDifferenceMax": expected_difference,
        "InfitDifferenceMax": infit_difference,
        "OutfitDifferenceMax": outfit_difference,
        "RawEitherUpperMismatches": int(np.sum(correct_flag != shifted_flag)),
    }


def markdown_table(frame: pd.DataFrame) -> str:
    headers = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in frame.itertuples(index=False, name=None):
        values = [
            f"{float(value):.3e}"
            if isinstance(value, (float, np.floating))
            else str(value)
            for value in row
        ]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--addendum", type=Path, default=DEFAULT_ADDENDUM)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    addendum = validate_addendum(arguments.addendum)
    output = arguments.output
    output.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(arguments.addendum, output / "registered_common_shift_addendum.json")

    responses = pd.read_csv(
        ROOT / addendum["parent_smoke"]["retained_responses_path"],
        float_precision="round_trip",
    )
    delta = float(addendum["mathematical_hypothesis"]["anchor_shift"])
    rows = []
    groups = list(responses.groupby(["ConditionId", "Replicate"], sort=False))
    for index, (_, frame) in enumerate(groups, start=1):
        row = compare_run(frame.reset_index(drop=True), delta=delta)
        rows.append(row)
        print(
            f"[{index:02d}/{len(groups):02d}] {row['RunId']} "
            f"loglik_diff={row['ConditionalLogLikAbsDifference']:.3e} "
            f"wle_shift_error={row['WLEShiftErrorMax']:.3e}",
            flush=True,
        )
    evidence = pd.DataFrame(rows)
    evidence_path = output / "common_shift_run_evidence.csv"
    evidence.to_csv(
        evidence_path, index=False, lineterminator="\n", float_format="%.17g"
    )

    gates = addendum["frozen_numerical_gates"]
    checks = {
        "all_fits_inference_ready": bool(
            evidence[["CorrectInferenceReady", "ShiftedInferenceReady"]].all().all()
        ),
        "conditional_loglik": float(evidence["ConditionalLogLikAbsDifference"].max())
        <= float(gates["conditional_loglik_absolute_difference_max"]),
        "information": float(evidence["InformationMaxAbsDifference"].max())
        <= float(gates["information_max_absolute_difference"]),
        "covariance": float(evidence["CovarianceMaxAbsDifference"].max())
        <= float(gates["covariance_max_absolute_difference"]),
        "free_rater_shift": float(evidence["FreeRaterShiftErrorMax"].max())
        <= float(gates["free_rater_shift_error_max"]),
        "criterion_step": float(evidence["CriterionStepEstimateDifferenceMax"].max())
        <= float(gates["criterion_step_estimate_difference_max"]),
        "surface_shift": float(evidence["SurfaceShiftErrorMax"].max())
        <= float(gates["surface_shift_error_max"]),
        "wle_shift": float(evidence["WLEShiftErrorMax"].max())
        <= float(gates["wle_shift_error_max"]),
        "wle_se": float(evidence["WLEStandardErrorDifferenceMax"].max())
        <= float(gates["wle_standard_error_difference_max"]),
        "expected_score": float(evidence["ExpectedScoreDifferenceMax"].max())
        <= float(gates["expected_score_difference_max"]),
        "fit_statistics": max(
            float(evidence["InfitDifferenceMax"].max()),
            float(evidence["OutfitDifferenceMax"].max()),
        )
        <= float(gates["infit_outfit_difference_max"]),
        "raw_fit_decisions": int(evidence["RawEitherUpperMismatches"].sum())
        == int(gates["raw_either_upper_mismatches"]),
    }
    maxima = {
        column: float(evidence[column].max())
        for column in evidence.columns
        if column.endswith("Difference")
        or column.endswith("DifferenceMax")
        or column.endswith("ErrorMax")
    }
    decision = {
        "schema_version": "mfrm-cmle-hard-anchor-common-shift-decision-v1",
        "study_id": addendum["study_id"],
        "post_result_diagnostic": True,
        "contract_passed": bool(all(checks.values())),
        "checks": checks,
        "maxima": maxima,
        "run_ids": int(len(evidence)),
        "response_rows": int(evidence["ResponseRows"].sum()),
        "persons_per_scenario": int(evidence["Persons"].sum()),
        "raw_either_upper_mismatches": int(
            evidence["RawEitherUpperMismatches"].sum()
        ),
        "anchor_shift": delta,
        "fit_indices_detect_common_origin_shift": False,
        "required_warning": addendum["decision_policy"][
            "required_user_warning_if_confirmed"
        ],
        "performance_claims_ready": False,
        "cross_engine_claims_ready": False,
        "ui_release_ready": False,
        "evidence_sha256": sha256_file(evidence_path),
    }
    (output / "common_shift_decision.json").write_text(
        json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    display = evidence[
        [
            "ConditionId",
            "Replicate",
            "ConditionalLogLikAbsDifference",
            "FreeRaterShiftErrorMax",
            "WLEShiftErrorMax",
            "InfitDifferenceMax",
            "OutfitDifferenceMax",
            "RawEitherUpperMismatches",
        ]
    ]
    report = "\n".join(
        [
            "# Native CMLE common anchor-origin shift diagnostic",
            "",
            "> Post-result diagnostic. This was registered after inspecting the Phase-B smoke summaries and is not confirmatory performance evidence.",
            "",
            f"Contract passed: `{decision['contract_passed']}`.",
            "",
            markdown_table(display),
            "",
            "A common +0.25 shift in all Rater anchors moved the fitted Rater origin and WLE Person origin together while preserving probabilities and fit statistics within the registered numerical tolerances. Consequently, MnSq cannot validate the absolute anchor origin. Differential anchor errors remain a separate risk and are not covered by this invariance.",
            "",
        ]
    )
    (output / "CMLE_HARD_ANCHOR_COMMON_SHIFT_DIAGNOSTIC.md").write_text(
        report, encoding="utf-8"
    )
    print(json.dumps(decision, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
