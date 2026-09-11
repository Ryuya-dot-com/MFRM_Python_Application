#!/usr/bin/env python3
"""Run the registered repository-only CMLE one-click result contract."""

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

from mfrm_app.cmle_one_click import (  # noqa: E402
    CMLE_ONE_CLICK_CARDS,
    run_cmle_one_click_analysis,
    summarize_person_fit_for_cards,
)
from mfrm_app.cmle_wle_fit import (  # noqa: E402
    compute_cmle_wle_person_fit,
    compute_fixed_calibration_person_fit,
)
from mfrm_app.decision_stability import evaluate_fit_mnsq  # noqa: E402


PLAN = ROOT / "validation/cmle_one_click_result_contract_plan_20260810.json"
OUTPUT = ROOT / "validation/cmle_one_click_result_contract_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, float_format="%.17g")


def validate_plan() -> dict[str, object]:
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    if plan.get("study_id") != "cmle_one_click_result_contract_v1":
        raise ValueError("Unexpected one-click result plan.")
    mismatches = [
        relative
        for relative, expected in plan["parent_identity"].items()
        if not (ROOT / relative).is_file()
        or sha256_file(ROOT / relative) != expected
    ]
    if mismatches:
        raise ValueError(f"One-click parent identity failed: {mismatches}")
    return plan


def interior_frame(*, extremes: bool) -> pd.DataFrame:
    scores = {
        "P1": [0, 1, 1, 2],
        "P2": [1, 0, 2, 1],
        "P3": [2, 1, 1, 0],
        "P4": [1, 2, 0, 1],
        "P5": [0, 2, 1, 1],
        "P6": [2, 0, 1, 1],
    }
    if extremes:
        scores.update({"P7": [0, 0, 0, 0], "P8": [2, 2, 2, 2]})
    units = [("R1", "C1"), ("R1", "C2"), ("R2", "C1"), ("R2", "C2")]
    return pd.DataFrame(
        [
            (person, rater, criterion, score)
            for person, values in scores.items()
            for (rater, criterion), score in zip(units, values, strict=True)
        ],
        columns=["Person", "Rater", "Criterion", "Score"],
    )


def boundary_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            (f"P{index:03d}", rater, score)
            for index in range(20)
            for rater, score in (("R1", 1), ("R2", 0))
        ],
        columns=["Person", "Rater", "Score"],
    )


def structural_frame() -> pd.DataFrame:
    patterns = ([0, 1, 2], [1, 2, 0], [2, 0, 1], [0, 2, 1])
    rows = []
    for person_index in range(12):
        rater = "R1" if person_index < 6 else "R2"
        rows.extend(
            (f"P{person_index}", rater, criterion, score)
            for criterion, score in zip(
                ["C1", "C2", "C3"], patterns[person_index % 4], strict=True
            )
        )
    return pd.DataFrame(rows, columns=["Person", "Rater", "Criterion", "Score"])


def case_inputs(plan: dict[str, object]):
    registered = {row["case_id"]: row for row in plan["analytical_cases"]}
    anchors = pd.DataFrame(
        [
            {
                "ParameterType": "Facet",
                "Facet": row["Facet"],
                "Level": row["Level"],
                "Value": row["Value"],
            }
            for row in registered["pcm_interior_differential_hard_anchors"][
                "anchors"
            ]
        ]
    )
    common = {
        "person_col": "Person",
        "score_col": "Score",
        "rating_min": 0,
        "gtol": 1e-8,
        "maxiter": 800,
        "display_decimals": 3,
    }
    return [
        (
            registered["rsm_interior_fit_sample_extremes"],
            interior_frame(extremes=True),
            {
                **common,
                "facet_cols": ["Rater", "Criterion"],
                "rating_max": 2,
                "model": "RSM",
            },
        ),
        (
            registered["pcm_interior_differential_hard_anchors"],
            interior_frame(extremes=False),
            {
                **common,
                "facet_cols": ["Rater", "Criterion"],
                "rating_max": 2,
                "model": "PCM",
                "step_facet": "Criterion",
                "hard_anchors": anchors,
            },
        ),
        (
            registered["binary_conditional_support_boundary"],
            boundary_frame(),
            {
                **common,
                "facet_cols": ["Rater"],
                "rating_max": 1,
                "model": "RSM",
            },
        ),
        (
            registered["rsm_structurally_disconnected"],
            structural_frame(),
            {
                **common,
                "facet_cols": ["Rater", "Criterion"],
                "rating_max": 2,
                "model": "RSM",
            },
        ),
        (
            registered["invalid_missing_score_column"],
            interior_frame(extremes=False).drop(columns=["Score"]),
            {
                **common,
                "facet_cols": ["Rater", "Criterion"],
                "rating_max": 2,
                "model": "RSM",
            },
        ),
    ]


def run_cases(plan: dict[str, object], output: Path):
    ledger_rows = []
    card_frames = []
    availability_frames = []
    handoff_rows = []
    anchor_rows = []
    for registered, frame, kwargs in case_inputs(plan):
        case_id = registered["case_id"]
        result = run_cmle_one_click_analysis(frame, **kwargs)
        summary = result["summary"].iloc[0]
        cards = result["cards"].copy()
        cards.insert(0, "CaseId", case_id)
        card_frames.append(cards)
        availability = result["availability"].copy()
        availability.insert(0, "CaseId", case_id)
        availability_frames.append(availability)
        ledger_rows.append(
            {
                "CaseId": case_id,
                "Model": registered["model"],
                "Mechanism": registered["mechanism"],
                "InputRows": len(frame),
                "ExpectedTerminalStatus": registered["expected_terminal_status"],
                "ObservedTerminalStatus": summary["TerminalStatus"],
                "TerminalMatch": (
                    summary["TerminalStatus"]
                    == registered["expected_terminal_status"]
                ),
                "ExpectedPersonScoringAttempted": registered[
                    "expected_person_scoring_attempted"
                ],
                "PersonScoringAttempted": summary["PersonScoringAttempted"],
                "PersonScoringAttemptMatch": bool(
                    summary["PersonScoringAttempted"]
                    == registered["expected_person_scoring_attempted"]
                ),
                "CalibrationReady": summary["CalibrationReady"],
                "PersonScoringReady": summary["PersonScoringReady"],
                "PersonCount": summary["PersonCount"],
                "ExtremePersonCount": summary["ExtremePersonCount"],
                "FitReviewPersonCount": summary["FitReviewPersonCount"],
                "RawDisplayMismatchCount": summary["RawDisplayMismatchCount"],
                "AnchoredMainEffectLevelCount": summary[
                    "AnchoredMainEffectLevelCount"
                ],
                "CardCount": len(cards),
                "CardOrderMatch": cards["Card"].tolist()
                == list(CMLE_ONE_CLICK_CARDS),
                "BilingualComplete": bool(
                    cards[["HeadlineEn", "HeadlineJa", "DetailEn", "DetailJa"]]
                    .astype(str)
                    .apply(lambda column: column.str.len().gt(0))
                    .all()
                    .all()
                ),
                "PublicSurfaceEnabled": summary["PublicSurfaceEnabled"],
                "Error": summary["Error"],
            }
        )
        if bool(summary["PersonScoringReady"]):
            observed = result["person_fit"]["persons"].copy()
            direct = compute_cmle_wle_person_fit(
                result["calibration"]["fit"], display_decimals=3
            )["persons"]
            paired = observed.merge(
                direct,
                on="Person",
                suffixes=("Observed", "Direct"),
                validate="one_to_one",
            )
            for _, row in paired.iterrows():
                handoff_rows.append(
                    {
                        "CaseId": case_id,
                        "Person": row["Person"],
                        **{
                            f"{column}AbsoluteDifference": abs(
                                float(row[f"{column}Observed"])
                                - float(row[f"{column}Direct"])
                            )
                            for column in (
                                "WLEEstimate",
                                "ConditionalWLEStandardError",
                                "Infit",
                                "Outfit",
                            )
                        },
                    }
                )
            facets = result["calibration"]["fit"]["facets"]["others"]
            for _, row in facets.loc[facets["Anchored"]].iterrows():
                expected = next(
                    item["Value"]
                    for item in registered["anchors"]
                    if item["Facet"] == row["Facet"]
                    and item["Level"] == row["Level"]
                )
                anchor_rows.append(
                    {
                        "CaseId": case_id,
                        "Facet": row["Facet"],
                        "Level": row["Level"],
                        "ExpectedValue": expected,
                        "Estimate": row["Estimate"],
                        "SE": row["SE"],
                        "EstimateExact": float(row["Estimate"]) == float(expected),
                        "SEExactZero": float(row["SE"]) == 0.0,
                    }
                )
    ledger = pd.DataFrame(ledger_rows)
    cards = pd.concat(card_frames, ignore_index=True)
    availability = pd.concat(availability_frames, ignore_index=True)
    handoff = pd.DataFrame(handoff_rows)
    anchors = pd.DataFrame(anchor_rows)
    write_csv(ledger, output / "case_ledger.csv")
    write_csv(cards, output / "result_cards.csv")
    write_csv(availability, output / "availability.csv")
    write_csv(handoff, output / "direct_handoff_comparison.csv")
    write_csv(anchors, output / "anchor_checks.csv")
    return ledger, cards, availability, handoff, anchors


def rounding_probe(output: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    raw = 1.5004
    persons = compute_fixed_calibration_person_fit(
        ["P1"],
        [0],
        np.array([[0.0, np.log(raw)]]),
        [0.0],
        display_decimals=3,
    )["persons"]
    persons["ExtremeScorePattern"] = False
    metrics = summarize_person_fit_for_cards(persons)
    evaluation = evaluate_fit_mnsq(raw, display_decimals=3)
    probe = pd.DataFrame(
        [
            {
                "Statistic": "Infit",
                "RawValue": evaluation["RawValue"],
                "DisplayValue": evaluation["DisplayValue"],
                "RawDecision": evaluation["RawDecision"],
                "DisplayDecision": evaluation["DisplayDecision"],
                "DisplayDecisionConsistent": evaluation[
                    "DisplayDecisionConsistent"
                ],
                "BoundaryStatus": evaluation["BoundaryStatus"],
                "DecisionInput": "finite_unrounded_mnsq",
                "AllPersonMNSQMismatchCount": metrics[
                    "RawDisplayMismatchCount"
                ],
            }
        ]
    )
    write_csv(probe, output / "rounding_boundary_probe.csv")
    return probe, metrics


def run_tests(output: Path) -> dict[str, object]:
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-q",
        "tests/test_cmle_one_click.py",
        "tests/test_cmle_workflow.py",
        "tests/test_cmle_person_scoring.py",
        "tests/test_cmle_wle_fit.py",
        "tests/test_cmle_hard_anchors.py",
    ]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    (output / "selected_tests_stdout.txt").write_text(
        completed.stdout, encoding="utf-8"
    )
    (output / "selected_tests_stderr.txt").write_text(
        completed.stderr, encoding="utf-8"
    )
    return {"passed": completed.returncode == 0, "returncode": completed.returncode}


def write_review(output: Path, results: dict[str, object]) -> None:
    text = f"""# CMLE One-Click Result Contract Critical Review

## Decision

The registered Streamlit-free contract **{'passed' if results['contract_passed'] else 'failed'}**. This is an orchestration and presentation-data result, not usability evidence and not a public estimator release.

## Guarded sequence

Five analytical cases were retained. Exactly {results['person_scoring_attempted']} reached fixed-calibration WLE/MnSq, while the input, structural-identification, and finite-existence stops did not run downstream scoring. The two ready handoffs reproduced direct WLE estimate, conditional SE, Infit, and Outfit values with maximum absolute difference `{results['handoff_max_absolute_difference']:.17g}`.

## Full precision and fit decisions

The registered `1.5004` Infit probe remained `noisy` even though its three-decimal display is `1.500` and would imply `acceptable`. Raw values drive all card counts. The probe is a constructed boundary demonstration, not an empirical rate or a guarantee that three or six display decimals are universally sufficient.

## Anchors

The PCM differential-anchor case retained two supplied main-effect anchors exactly with SE zero. Its structural card warns that exact enforcement cannot validate anchor provenance, bias, or common-scale adequacy. Criterion main-effect anchoring does not anchor Criterion-specific PCM transitions.

## Remaining risks

- The last uncertainty card is always withheld: conditional WLE SE fixes the fitted CMLE calibration, while ZSTD, p-values, confidence intervals, and total SE are not qualified.
- The contract has five small deterministic analytical cases. It is not high-dimensional sparse-PCM, drift, power, coverage, or task-comprehension evidence.
- New Persons/unseen units, group/step anchors, GPCM, automatic fallback, saved-result identity, and Streamlit wiring remain outside scope.
- A single function call reduces user actions but does not show that users understand the staged warnings.

## Product implication

The computation can now be triggered through one guarded function while returning six stable bilingual cards. Public exposure should wait for task-based comprehension testing and the remaining estimator/uncertainty gates; a button must never turn blocked or withheld stages green.
"""
    (output / "CMLE_ONE_CLICK_RESULT_CONTRACT_CRITICAL_REVIEW.md").write_text(
        text, encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    plan = validate_plan()
    if args.output.exists():
        raise FileExistsError(f"Evidence output exists: {args.output}")
    args.output.mkdir(parents=True)
    ledger, cards, availability, handoff, anchors = run_cases(plan, args.output)
    probe, probe_metrics = rounding_probe(args.output)
    tests = run_tests(args.output)

    diff_columns = [column for column in handoff if column.endswith("AbsoluteDifference")]
    handoff_max = float(handoff[diff_columns].to_numpy(dtype=float).max())
    expected_counts = {
        row["expected_terminal_status"]: sum(
            case["expected_terminal_status"] == row["expected_terminal_status"]
            for case in plan["analytical_cases"]
        )
        for row in plan["analytical_cases"]
    }
    observed_counts = ledger["ObservedTerminalStatus"].value_counts().to_dict()
    gates = {
        "identity_passed": True,
        "status_accounting_passed": bool(
            len(ledger) == 5
            and ledger["TerminalMatch"].all()
            and observed_counts == expected_counts
        ),
        "no_forbidden_scoring_passed": bool(
            ledger["PersonScoringAttemptMatch"].all()
            and int(ledger["PersonScoringAttempted"].sum()) == 2
        ),
        "handoff_identity_passed": bool(
            handoff["CaseId"].nunique() == 2
            and handoff_max <= 1e-12
            and len(handoff) == 14
            and handoff.groupby("CaseId").size().to_dict()
            == {
                "pcm_interior_differential_hard_anchors": 6,
                "rsm_interior_fit_sample_extremes": 8,
            }
        ),
        "anchor_contract_passed": bool(
            len(anchors) == 2
            and anchors["EstimateExact"].all()
            and anchors["SEExactZero"].all()
            and int(ledger["AnchoredMainEffectLevelCount"].sum()) == 2
        ),
        "card_contract_passed": bool(
            len(cards) == 30
            and ledger["CardOrderMatch"].all()
            and ledger["BilingualComplete"].all()
            and not ledger["PublicSurfaceEnabled"].any()
            and cards.groupby("CaseId")["Card"].apply(list).apply(
                lambda names: names == list(CMLE_ONE_CLICK_CARDS)
            ).all()
        ),
        "rounding_contract_passed": bool(
            len(probe) == 1
            and probe.iloc[0]["RawDecision"] == "noisy"
            and probe.iloc[0]["DisplayDecision"] == "acceptable"
            and not bool(probe.iloc[0]["DisplayDecisionConsistent"])
            and int(probe_metrics["RawDisplayMismatchCount"]) == 2
            and bool(probe_metrics["DecisionsUseUnroundedMNSQ"])
        ),
        "scope_passed": bool(
            cards.groupby("CaseId").tail(1)["Status"].eq("withheld").all()
            and not availability["PublicSurfaceEnabled"].any()
        ),
        "tests_passed": bool(tests["passed"]),
    }
    contract_passed = all(gates.values())
    results = {
        **gates,
        "contract_passed": contract_passed,
        "person_scoring_attempted": int(ledger["PersonScoringAttempted"].sum()),
        "handoff_max_absolute_difference": handoff_max,
        "terminal_status_counts": observed_counts,
        "selected_tests": tests,
    }
    write_review(args.output, results)
    output_files = sorted(
        path.relative_to(args.output).as_posix()
        for path in args.output.rglob("*")
        if path.is_file() and path.name != "decision.json"
    )
    decision = json_safe(
        {
            "study_id": plan["study_id"],
            "plan_sha256": sha256_file(PLAN),
            "contract_passed": contract_passed,
            "implementation_sha256": {
                "mfrm_app/cmle_one_click.py": sha256_file(
                    ROOT / "mfrm_app/cmle_one_click.py"
                ),
                "tests/test_cmle_one_click.py": sha256_file(
                    ROOT / "tests/test_cmle_one_click.py"
                ),
                "validation/cmle_one_click_result_contract.py": sha256_file(
                    Path(__file__)
                ),
            },
            "results": results,
            "output_sha256": {
                name: sha256_file(args.output / name) for name in output_files
            },
            "interpretation": {
                "one_click": "orchestration contract, not usability evidence",
                "decision_input": "finite unrounded MnSq",
                "public_ui": "withheld",
            },
        }
    )
    (args.output / "decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(decision, indent=2, ensure_ascii=False, allow_nan=False))
    if not contract_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
