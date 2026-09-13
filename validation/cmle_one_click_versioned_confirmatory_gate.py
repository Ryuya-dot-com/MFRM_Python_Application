#!/usr/bin/env python3
"""Validate instrument identity and directional Wilson gate mechanics."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist
import subprocess
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle_one_click_cognitive_interview import INSTRUMENT_VERSION  # noqa: E402
from mfrm_app.cmle_one_click_confirmatory_gate import (  # noqa: E402
    ANCHOR_PRIMARY_CASE,
    CONFIRMATORY_RESPONSE_COLUMNS,
    DIRECTIONAL_RULES,
    build_confirmatory_assignment_schedule,
    build_instrument_version_record,
    build_wilson_planning_table,
    compute_instrument_content_identity,
    confirmatory_human_gate_status,
    evaluate_directional_confirmatory_gate,
    validate_instrument_version_ledger,
    wilson_upper_bound,
)


PLAN = ROOT / "validation/cmle_one_click_versioned_confirmatory_gate_plan_20260810.json"
SOURCE = ROOT / "validation/cmle_one_click_comprehension_readiness_20260810"
PARTICIPANT = SOURCE / "participant_task_packet.csv"
KEY_PATH = SOURCE / "researcher_scoring_key.csv"
PREVIEWS = SOURCE / "previews"
OUTPUT = ROOT / "validation/cmle_one_click_versioned_confirmatory_gate_20260810"


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
    if plan.get("study_id") != "cmle_one_click_versioned_confirmatory_gate_v1":
        raise ValueError("Unexpected versioned confirmatory gate plan.")
    mismatches = [
        relative
        for relative, expected in plan["parent_identity"].items()
        if not (ROOT / relative).is_file()
        or sha256_file(ROOT / relative) != expected
    ]
    if mismatches:
        raise ValueError(f"Versioned confirmatory parent identity failed: {mismatches}")
    return plan


def directional_rules_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                **{key: value for key, value in rule.items() if key != "DangerousResponseCodes"},
                "DangerousResponseCodesJSON": json.dumps(
                    list(rule["DangerousResponseCodes"]), separators=(",", ":")
                ),
                "AnyIncorrectIsDangerous": False,
                "PrimaryUnit": "one_registered_case_role_per_slot",
            }
            for rule in DIRECTIONAL_RULES
        ]
    )


def identity_audit(
    ledger: pd.DataFrame, preview_manifest: pd.DataFrame
) -> pd.DataFrame:
    baseline = str(ledger.iloc[0]["InstrumentContentSHA256"])
    participant_bytes = PARTICIPANT.read_bytes()
    key_bytes = KEY_PATH.read_bytes()
    cases = []

    def identity(
        participant_value: bytes = participant_bytes,
        key_value: bytes = key_bytes,
        manifest_value: pd.DataFrame = preview_manifest,
        rules_value=DIRECTIONAL_RULES,
    ) -> str:
        return compute_instrument_content_identity(
            participant_value,
            key_value,
            manifest_value,
            directional_rules=rules_value,
        )["InstrumentContentSHA256"]

    cases.append(("identical_replay", identity(), True))
    cases.append(("participant_task_byte_mutation", identity(participant_value=participant_bytes + b"x"), False))
    cases.append(("researcher_key_byte_mutation", identity(key_value=key_bytes + b"x"), False))
    preview_mutation = preview_manifest.copy()
    preview_mutation.loc[0, "SHA256"] = "0" * 64
    cases.append(("preview_manifest_mutation", identity(manifest_value=preview_mutation), False))
    rule_mutation = copy.deepcopy(list(DIRECTIONAL_RULES))
    rule_mutation[0] = {
        **rule_mutation[0],
        "DangerousResponseCodes": ("ready", "synthetic_mutation"),
    }
    cases.append(("directional_rule_mutation", identity(rules_value=rule_mutation), False))
    return pd.DataFrame(
        [
            {
                "IdentityCase": case_id,
                "BaselineInstrumentContentSHA256": baseline,
                "ObservedInstrumentContentSHA256": observed,
                "ExpectedSameIdentity": expected_same,
                "ObservedSameIdentity": observed == baseline,
                "ContractMatch": (observed == baseline) == expected_same,
            }
            for case_id, observed, expected_same in cases
        ]
    )


def ledger_mutation_audit(ledger: pd.DataFrame) -> pd.DataFrame:
    cases = []
    cases.append(("valid_frozen_ledger", ledger, True, set()))
    alias = pd.concat(
        [ledger, ledger.assign(InstrumentVersion="content_alias")], ignore_index=True
    )
    cases.append(("duplicate_content_alias", alias, False, {"duplicate_content_identity"}))
    cases.append(
        (
            "pooling_permission",
            ledger.assign(PoolingWithDifferentContentAllowed=True),
            False,
            {"cross_version_pooling_allowed"},
        )
    )
    cases.append(
        (
            "responses_before_freeze",
            ledger.assign(HumanResponsesBeforeFreeze=1),
            False,
            {"responses_before_freeze"},
        )
    )
    rows = []
    for case_id, value, expected_valid, expected_codes in cases:
        result = validate_instrument_version_ledger(value)
        observed_codes = set(result["failure_codes"])
        rows.append(
            {
                "LedgerCase": case_id,
                "ExpectedValid": expected_valid,
                "ObservedValid": result["valid"],
                "ExpectedCodes": ";".join(sorted(expected_codes)),
                "ObservedCodes": ";".join(result["failure_codes"]),
                "ContractMatch": result["valid"] == expected_valid
                and expected_codes.issubset(observed_codes),
            }
        )
    return pd.DataFrame(rows)


def correct_responses(assignment: pd.DataFrame, key: pd.DataFrame) -> pd.DataFrame:
    merged = assignment.merge(
        key[["CaseId", "Language", "ItemId", "CorrectOption"]],
        on=["CaseId", "Language"],
        validate="many_to_many",
    )
    return merged[
        ["InstrumentVersion", "ConfirmatorySlotId", "Language", "CaseId", "ItemId"]
    ].assign(ResponseCode=merged["CorrectOption"].to_numpy())[
        list(CONFIRMATORY_RESPONSE_COLUMNS)
    ]


def run_gate_case(
    case_id: str,
    responses: pd.DataFrame,
    key: pd.DataFrame,
    assignment: pd.DataFrame,
    *,
    minimum_n: int,
    expected_overall: bool,
    expected_passing_cells: int,
    expected_invalid: int,
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame]:
    result = evaluate_directional_confirmatory_gate(
        responses,
        key,
        assignment,
        registered_instrument_version=INSTRUMENT_VERSION,
        minimum_valid_per_cell=minimum_n,
    )
    summary = result["study_summary"].iloc[0]
    row = {
        "SyntheticGateCase": case_id,
        "ExpectedOverallPrimaryGatePassed": expected_overall,
        "ObservedOverallPrimaryGatePassed": bool(summary["OverallPrimaryGatePassed"]),
        "ExpectedPassingCells": expected_passing_cells,
        "ObservedPassingCells": int(summary["PrimaryPassingCells"]),
        "ExpectedInvalidAttemptedSlots": expected_invalid,
        "ObservedInvalidAttemptedSlots": int(summary["InvalidAttemptedSlots"]),
        "ContractMatch": bool(
            bool(summary["OverallPrimaryGatePassed"]) == expected_overall
            and int(summary["PrimaryPassingCells"]) == expected_passing_cells
            and int(summary["InvalidAttemptedSlots"]) == expected_invalid
        ),
        "Interpretation": "synthetic_gate_logic_only_no_human_outcome",
    }
    cells = result["cell_summary"].copy()
    cells.insert(0, "SyntheticGateCase", case_id)
    slots = result["slot_audit"].copy()
    slots.insert(0, "SyntheticGateCase", case_id)
    return row, cells, slots


def synthetic_gate_audit(key: pd.DataFrame):
    assignment24 = build_confirmatory_assignment_schedule(24)
    assignment25 = build_confirmatory_assignment_schedule(25)
    correct24 = correct_responses(assignment24, key)
    correct25 = correct_responses(assignment25, key)
    cases = []

    cases.append(
        run_gate_case(
            "zero_dangerous_n24",
            correct24,
            key,
            assignment24,
            minimum_n=24,
            expected_overall=False,
            expected_passing_cells=0,
            expected_invalid=0,
        )
    )
    cases.append(
        run_gate_case(
            "zero_dangerous_n25",
            correct25,
            key,
            assignment25,
            minimum_n=25,
            expected_overall=True,
            expected_passing_cells=12,
            expected_invalid=0,
        )
    )

    directional = correct25.copy()
    ja_anchor_slot = assignment25.loc[
        assignment25["Language"].eq("ja")
        & assignment25["CaseRole"].eq("anchored_primary")
    ].iloc[0]["ConfirmatorySlotId"]
    public_mask = (
        directional["ConfirmatorySlotId"].eq(ja_anchor_slot)
        & directional["CaseId"].eq(ANCHOR_PRIMARY_CASE)
        & directional["ItemId"].eq("public_sharing")
    )
    directional.loc[public_mask, "ResponseCode"] = "permitted"
    cases.append(
        run_gate_case(
            "one_ja_directional_public_error",
            directional,
            key,
            assignment25,
            minimum_n=25,
            expected_overall=False,
            expected_passing_cells=11,
            expected_invalid=0,
        )
    )

    conservative = correct25.copy()
    en_anchor_slot = assignment25.loc[
        assignment25["Language"].eq("en")
        & assignment25["CaseRole"].eq("anchored_primary")
    ].iloc[0]["ConfirmatorySlotId"]
    rounding_mask = (
        conservative["ConfirmatorySlotId"].eq(en_anchor_slot)
        & conservative["CaseId"].eq(ANCHOR_PRIMARY_CASE)
        & conservative["ItemId"].eq("rounding_vignette")
    )
    conservative.loc[rounding_mask, "ResponseCode"] = "cannot_decide"
    cases.append(
        run_gate_case(
            "conservative_rounding_error_not_dangerous",
            conservative,
            key,
            assignment25,
            minimum_n=25,
            expected_overall=True,
            expected_passing_cells=12,
            expected_invalid=0,
        )
    )

    display_error = correct25.copy()
    display_error.loc[rounding_mask, "ResponseCode"] = "display_acceptable"
    cases.append(
        run_gate_case(
            "display_rounding_reclassification_dangerous",
            display_error,
            key,
            assignment25,
            minimum_n=25,
            expected_overall=False,
            expected_passing_cells=11,
            expected_invalid=0,
        )
    )

    en_slot = correct25.loc[
        correct25["Language"].eq("en"), "ConfirmatorySlotId"
    ].iloc[0]
    missing = correct25.drop(
        correct25.loc[
            correct25["ConfirmatorySlotId"].eq(en_slot)
        ].index[0]
    )
    cases.append(
        run_gate_case(
            "missing_response_slot_retained",
            missing,
            key,
            assignment25,
            minimum_n=25,
            expected_overall=False,
            expected_passing_cells=6,
            expected_invalid=1,
        )
    )
    duplicate = pd.concat(
        [correct25, correct25.loc[correct25["ConfirmatorySlotId"].eq(en_slot)].iloc[[0]]],
        ignore_index=True,
    )
    cases.append(
        run_gate_case(
            "duplicate_response_slot_retained",
            duplicate,
            key,
            assignment25,
            minimum_n=25,
            expected_overall=False,
            expected_passing_cells=6,
            expected_invalid=1,
        )
    )
    unknown = correct25.copy()
    unknown.loc[unknown["ConfirmatorySlotId"].eq(en_slot).idxmax(), "ResponseCode"] = "unknown_option"
    cases.append(
        run_gate_case(
            "unknown_option_slot_retained",
            unknown,
            key,
            assignment25,
            minimum_n=25,
            expected_overall=False,
            expected_passing_cells=6,
            expected_invalid=1,
        )
    )
    drift = correct25.copy()
    drift.loc[drift["ConfirmatorySlotId"].eq(en_slot).idxmax(), "CaseId"] = "assignment_drift"
    cases.append(
        run_gate_case(
            "assignment_drift_slot_retained",
            drift,
            key,
            assignment25,
            minimum_n=25,
            expected_overall=False,
            expected_passing_cells=6,
            expected_invalid=1,
        )
    )
    mixed = correct25.copy()
    mixed.loc[mixed["ConfirmatorySlotId"].eq(en_slot), "InstrumentVersion"] = "unregistered_version"
    cases.append(
        run_gate_case(
            "mixed_version_slot_retained",
            mixed,
            key,
            assignment25,
            minimum_n=25,
            expected_overall=False,
            expected_passing_cells=6,
            expected_invalid=1,
        )
    )

    ledger = pd.DataFrame([case[0] for case in cases])
    cells = pd.concat([case[1] for case in cases], ignore_index=True)
    slots = pd.concat([case[2] for case in cases], ignore_index=True)
    return ledger, cells, slots, assignment24, assignment25


def wilson_formula_audit(planning: pd.DataFrame) -> pd.DataFrame:
    rows = []
    z = NormalDist().inv_cdf(0.95)
    for n in (24, 25, 50, 100, 200):
        direct_zero = z * z / (n + z * z)
        observed = wilson_upper_bound(0, n, confidence=0.95)
        rows.append(
            {
                "Audit": "zero_error_closed_form",
                "Trials": n,
                "Errors": 0,
                "DirectValue": direct_zero,
                "FunctionValue": observed,
                "AbsoluteDifference": abs(direct_zero - observed),
                "Passed": abs(direct_zero - observed) <= 2e-16,
            }
        )
    for n in (25, 50, 100):
        values = [wilson_upper_bound(errors, n) for errors in range(n + 1)]
        rows.append(
            {
                "Audit": "monotone_in_error_count",
                "Trials": n,
                "Errors": -1,
                "DirectValue": math.nan,
                "FunctionValue": math.nan,
                "AbsoluteDifference": math.nan,
                "Passed": all(right >= left for left, right in zip(values, values[1:])),
            }
        )
    row24 = planning.set_index("ValidEligibleSlots").loc[24]
    row25 = planning.set_index("ValidEligibleSlots").loc[25]
    rows.append(
        {
            "Audit": "strict_n24_n25_boundary",
            "Trials": -1,
            "Errors": 0,
            "DirectValue": 0.10,
            "FunctionValue": math.nan,
            "AbsoluteDifference": math.nan,
            "Passed": bool(
                row24["PrimaryZeroErrorUpper"] >= 0.10
                and row25["PrimaryZeroErrorUpper"] < 0.10
                and row24["PrimaryMaximumPassingErrors"] == -1
                and row25["PrimaryMaximumPassingErrors"] == 0
            ),
        }
    )
    return pd.DataFrame(rows)


def run_tests(output: Path) -> dict[str, object]:
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-q",
        "tests/test_cmle_one_click_confirmatory_gate.py",
        "tests/test_cmle_one_click_cognitive_interview.py",
        "tests/test_cmle_one_click_comprehension.py",
        "tests/test_cmle_one_click.py",
        "tests/test_cmle_one_click_archive.py",
        "tests/test_decision_stability.py",
        "tests/test_threshold_decision_integration.py",
    ]
    completed = subprocess.run(
        command, cwd=ROOT, text=True, capture_output=True, check=False
    )
    (output / "selected_tests_stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output / "selected_tests_stderr.txt").write_text(completed.stderr, encoding="utf-8")
    return {"passed": completed.returncode == 0, "returncode": completed.returncode, "command": command}


def write_documents(output: Path, results: dict[str, object]) -> None:
    specification = """# Directional dangerous-error gate specification

## Separate estimands

The existing non-compensatory comprehension rule requires all six critical items correct. The confirmatory safety gate has a narrower numerator: only a registered unsafe response direction in its eligible primary analytical state. A conservative or uncertain wrong answer remains a comprehension error but is not called false-green, invented scoring, or another dangerous misconception.

## One primary observation per slot and domain

Each future slot receives one anchored PCM case and one rotating blocked case. False-green calibration and invented Person scores use the blocked case. Rounding, calibration uncertainty, anchor validity, and public sharing use the anchored PCM case. A slot therefore contributes at most once to each domain, avoiding repeated-response inflation.

## Primary rule

For each of six domains separately in English and Japanese, compute the one-sided 95% Wilson upper bound for the dangerous-error proportion. A cell passes only when its future prospectively registered minimum valid n is met and the unrounded bound is strictly below 0.10. All 12 cells must pass. No default minimum n exists in the implementation.

## Sensitivity and limits

A Bonferroni one-sided `1 - 0.05/12` Wilson upper bound is reported because 12 separate 95% cell bounds are not a simultaneous 95% familywise statement. The sensitivity does not silently replace the registered primary rule. The planning table is arithmetic, not a power analysis or sample-size recommendation.
"""
    (output / "DIRECTIONAL_GATE_SPECIFICATION.md").write_text(specification, encoding="utf-8")

    review = f"""# Versioned confirmatory gate critical review

## Decision

The no-human-data mechanics contract **{'passed' if results['contract_passed'] else 'failed'}**. Frozen instrument identity, directional error definitions, two-case slot assignment, Wilson arithmetic, invalid-record retention, and language-separated gate behavior are reproducible. No pilot or confirmatory result exists.

## The n=24 boundary matters

With zero dangerous errors, a one-sided 95% Wilson upper bound at `n=24` is still not strictly below 0.10; at `n=25` it becomes arithmetically below 0.10. This confirms that the planned at-least-24 pilot cannot promote the UI under the confirmatory rule even in a zero-error cell. It does not make 25 an adequate sample-size recommendation. The 12-cell Bonferroni sensitivity remains above 0.10 at `n=25` with zero errors.

## Directionality matters

`display_acceptable` for raw 1.5004 versus displayed 1.500 is a dangerous rounding reclassification. `cannot_decide` is incorrect but conservative and is retained in the any-comprehension-error count without entering the dangerous numerator. Likewise, a ready result mistakenly called blocked is not relabelled false-green. This prevents a broad error rate from masquerading as a directional safety rate.

## Invalidity and versioning

Missing, duplicate, unknown-option, assignment-drift, and mixed-version slots remain in the slot audit and lose eligibility rather than being scored correct. They reduce valid n and can therefore block cells. Different task/key/preview/rule content changes the composite instrument identity; content aliases and cross-version pooling permission are rejected.

## Remaining threats

- A future protocol must set minimum n, anticipated error rates, power/precision rationale, attrition handling, recruitment strata, moderator/site effects, and any design-effect correction before responses.
- Per-cell 95% Wilson bounds are not a familywise 95% guarantee; both primary and Bonferroni sensitivity should be reported without switching rules after results.
- Two cases per slot can create learning/order effects despite alternation. The primary-domain mapping prevents double counting but not carryover.
- Separate Japanese and English passes do not establish equivalence.
- Instrument hashes establish content identity, not translation quality, accessibility, ethics approval, privacy compliance, or human comprehension.

Human participants remain zero, no minimum confirmatory n is registered, and the public UI remains withheld.
"""
    (output / "VERSIONED_CONFIRMATORY_GATE_CRITICAL_REVIEW.md").write_text(review, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    plan = validate_plan()
    if args.output.exists():
        raise FileExistsError(f"Evidence output exists: {args.output}")
    args.output.mkdir(parents=True)

    key = pd.read_csv(KEY_PATH)
    ledger, preview_manifest = build_instrument_version_record(
        PARTICIPANT, KEY_PATH, PREVIEWS
    )
    rules = directional_rules_frame()
    identity_cases = identity_audit(ledger, preview_manifest)
    ledger_cases = ledger_mutation_audit(ledger)
    planning = build_wilson_planning_table()
    formula_audit = wilson_formula_audit(planning)
    gate_ledger, gate_cells, gate_slots, assignment24, assignment25 = synthetic_gate_audit(key)
    human_gate = confirmatory_human_gate_status()

    registered_rules = [
        {
            "ItemId": row["item_id"],
            "Domain": row["domain"],
            "PrimaryCaseRole": row["primary_case_role"],
            "EligibleCorrectOption": row["eligible_correct_option"],
            "DangerousResponseCodes": list(row["dangerous_response_codes"]),
        }
        for row in plan["directional_cells"]
    ]
    implemented_rules = [
        {
            "ItemId": row["ItemId"],
            "Domain": row["Domain"],
            "PrimaryCaseRole": row["PrimaryCaseRole"],
            "EligibleCorrectOption": row["EligibleCorrectOption"],
            "DangerousResponseCodes": list(row["DangerousResponseCodes"]),
        }
        for row in DIRECTIONAL_RULES
    ]
    rule_registration_match = implemented_rules == registered_rules

    write_csv(ledger, args.output / "instrument_version_ledger.csv")
    write_csv(preview_manifest, args.output / "preview_identity_manifest.csv")
    write_csv(rules, args.output / "directional_dangerous_error_rules.csv")
    write_csv(identity_cases, args.output / "instrument_identity_mutation_audit.csv")
    write_csv(ledger_cases, args.output / "version_ledger_mutation_audit.csv")
    write_csv(planning, args.output / "wilson_candidate_n_arithmetic.csv")
    write_csv(formula_audit, args.output / "wilson_formula_audit.csv")
    write_csv(assignment24, args.output / "synthetic_assignment_n24.csv")
    write_csv(assignment25, args.output / "synthetic_assignment_n25.csv")
    write_csv(gate_ledger, args.output / "synthetic_gate_case_ledger.csv")
    write_csv(gate_cells, args.output / "synthetic_gate_cell_summary.csv")
    write_csv(gate_slots, args.output / "synthetic_gate_slot_audit.csv")
    write_csv(human_gate, args.output / "human_gate_status.csv")
    tests = run_tests(args.output)

    n24 = planning.set_index("ValidEligibleSlots").loc[24]
    n25 = planning.set_index("ValidEligibleSlots").loc[25]
    gates = {
        "identity_passed": True,
        "version_identity_passed": bool(
            len(ledger) == 1
            and validate_instrument_version_ledger(ledger)["valid"]
            and len(preview_manifest) == 10
            and identity_cases["ContractMatch"].all()
            and ledger_cases["ContractMatch"].all()
        ),
        "directional_rules_passed": bool(
            len(rules) == 6
            and rule_registration_match
            and rules["Domain"].nunique() == 6
            and rules["ItemId"].nunique() == 6
            and not rules["AnyIncorrectIsDangerous"].any()
            and set(rules["PrimaryCaseRole"]) == {"anchored_primary", "blocked_primary"}
        ),
        "assignment_passed": bool(
            len(assignment24) == 96
            and len(assignment25) == 100
            and assignment25.groupby("ConfirmatorySlotId").size().eq(2).all()
            and assignment25.groupby(["ConfirmatorySlotId", "CaseRole"]).size().eq(1).all()
        ),
        "wilson_arithmetic_passed": bool(
            formula_audit["Passed"].all()
            and planning["ValidEligibleSlots"].tolist()
            == plan["planning_table_contract"]["candidate_n"]
            and n24["PrimaryZeroErrorUpper"] >= 0.10
            and n25["PrimaryZeroErrorUpper"] < 0.10
            and n24["PrimaryMaximumPassingErrors"] == -1
            and n25["PrimaryMaximumPassingErrors"] == 0
            and n25["FamilywiseZeroErrorUpper"] > 0.10
        ),
        "synthetic_gate_passed": bool(
            len(gate_ledger) == 10
            and gate_ledger["ContractMatch"].all()
            and len(gate_cells) == 120
            and gate_cells["DecisionUsesRawBound"].all()
        ),
        "human_and_public_withheld_passed": bool(
            len(human_gate) == 1
            and int(human_gate.iloc[0]["HumanParticipants"]) == 0
            and not bool(human_gate.iloc[0]["PilotResultAvailable"])
            and not bool(human_gate.iloc[0]["ConfirmatoryResultAvailable"])
            and not bool(human_gate.iloc[0]["MinimumValidPerCellRegistered"])
            and not bool(human_gate.iloc[0]["LanguageEquivalenceAvailable"])
            and not bool(human_gate.iloc[0]["PublicSurfaceEnabled"])
        ),
        "tests_passed": bool(tests["passed"]),
    }
    contract_passed = all(gates.values())
    results = {
        **gates,
        "contract_passed": contract_passed,
        "instrument_content_sha256": ledger.iloc[0]["InstrumentContentSHA256"],
        "directional_domains": len(rules),
        "primary_language_item_cells": 12,
        "n24_zero_error_primary_upper": float(n24["PrimaryZeroErrorUpper"]),
        "n25_zero_error_primary_upper": float(n25["PrimaryZeroErrorUpper"]),
        "n25_zero_error_familywise_sensitivity_upper": float(n25["FamilywiseZeroErrorUpper"]),
        "synthetic_gate_cases": len(gate_ledger),
        "human_participants": 0,
        "pilot_result_available": False,
        "confirmatory_result_available": False,
        "minimum_valid_per_cell_registered": False,
        "public_surface_enabled": False,
        "selected_tests": tests,
    }
    write_documents(args.output, results)
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
            "contract_interpretation": "no_human_data_version_and_gate_mechanics_only",
            "implementation_sha256": {
                "mfrm_app/cmle_one_click_confirmatory_gate.py": sha256_file(
                    ROOT / "mfrm_app/cmle_one_click_confirmatory_gate.py"
                ),
                "tests/test_cmle_one_click_confirmatory_gate.py": sha256_file(
                    ROOT / "tests/test_cmle_one_click_confirmatory_gate.py"
                ),
                "validation/cmle_one_click_versioned_confirmatory_gate.py": sha256_file(Path(__file__)),
            },
            "results": results,
            "output_sha256": {
                name: sha256_file(args.output / name) for name in output_files
            },
            "interpretation": {
                "n24_n25": "arithmetic boundary only, not sample-size recommendation",
                "dangerous_error": "directional eligible-state numerator, distinct from any comprehension error",
                "familywise": "Bonferroni sensitivity reported; primary remains per-cell one-sided 95% Wilson",
                "human_data": "none",
                "language_equivalence": "not assessed",
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
