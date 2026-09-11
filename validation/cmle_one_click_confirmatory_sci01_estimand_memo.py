#!/usr/bin/env python3
"""Generate the prospectively registered private SCI-01 estimand memo."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle_one_click_confirmatory_estimand import (  # noqa: E402
    SCI01_BUNDLE_MEMBERS,
    SCI01_OPTION_IDS,
    build_private_sci01_bundle,
    build_sci01_selection_template,
    sci01_human_gate_status,
    validate_private_sci01_bundle,
    validate_sci01_html,
    validate_sci01_selection,
)

PLAN = ROOT / "validation/cmle_one_click_confirmatory_sci01_estimand_memo_plan_20260810.json"
SURFACE = ROOT / "validation/cmle_one_click_confirmatory_protocol_preflight_20260810/partial_identification_sensitivity.csv"
OUTPUT = ROOT / "validation/cmle_one_click_confirmatory_sci01_estimand_memo_20260810"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.write_text(frame.to_csv(index=False, float_format="%.17g"), encoding="utf-8")


def validate_registration() -> dict[str, object]:
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    failures = []
    for relative, expected in plan["parent_evidence_sha256"].items():
        path = ROOT / relative
        if not path.is_file() or sha256_file(path) != expected:
            failures.append(relative)
    if plan["registered_options_in_display_order_only"] != list(SCI01_OPTION_IDS):
        failures.append("registered_options")
    if plan["repository_recommendation_allowed"] or plan["repository_default_allowed"] or plan["automatic_selection_allowed"]:
        failures.append("selection_policy")
    if failures:
        raise RuntimeError(f"Registration identity/contract failed: {failures}")
    return plan


def selection_attack_audit() -> pd.DataFrame:
    blank = build_sci01_selection_template()
    hidden = blank.copy(); hidden.loc[0, "SelectedOptionId"] = SCI01_OPTION_IDS[0]
    incomplete = blank.copy(); incomplete.loc[0, "DecisionStatus"] = "resolved_prospectively"; incomplete.loc[0, "SelectedOptionId"] = SCI01_OPTION_IDS[2]
    repository = incomplete.copy()
    for col in ("AssumptionSetReference", "RationaleReference", "OwnerAttestationReference", "InvalidityCodeMapReference"):
        repository.loc[0, col] = "SYNTHETIC"
    repository.loc[0, "ScientificOwnerRole"] = "repository agent"
    repository.loc[0, "DecisionRecordedBeforeOutcomes"] = True
    late = repository.copy(); late.loc[0, "ScientificOwnerRole"] = "external scientific lead"; late.loc[0, "HumanOutcomesInspectedBeforeDecision"] = True
    complete = late.copy(); complete.loc[0, "HumanOutcomesInspectedBeforeDecision"] = False
    cases = (
        ("blank_valid_incomplete", blank, True, False, ""),
        ("hidden_selection", hidden, False, False, "unresolved_row_contains_selection_material"),
        ("missing_evidence", incomplete, False, False, "missing_owner_or_evidence_reference"),
        ("repository_owner", repository, False, False, "repository_cannot_be_scientific_owner"),
        ("outcomes_inspected", late, False, False, "decision_not_prospectively_locked"),
        ("synthetic_external_complete_still_blocked", complete, True, True, ""),
    )
    rows = []
    for name, frame, expected_valid, expected_complete, required in cases:
        result = validate_sci01_selection(frame)
        rows.append({"Case": name, "ExpectedValid": expected_valid, "ActualValid": result["valid"],
            "ExpectedComplete": expected_complete, "ActualComplete": result["complete"],
            "RequiredFailureCode": required, "ActualFailureCodes": "|".join(result["failure_codes"]),
            "RecruitmentReady": result["RecruitmentReady"],
            "Passed": bool(result["valid"] == expected_valid and result["complete"] == expected_complete
                and not result["RecruitmentReady"] and (not required or required in result["failure_codes"]))})
    return pd.DataFrame(rows)


def same_data_divergence_audit(projection: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source_row, group in projection.groupby("SourceScenarioRow"):
        labels = set(group["ProjectedDecision"])
        if labels == {"pass", "proxy_fail", "observed_pass_not_robust"}:
            first = group.iloc[0]
            rows.append({"SourceScenarioRow": source_row, "DangerousErrors": int(first["DangerousErrors"]),
                "ValidEligibleN": int(first["ValidEligibleN"]), "InvalidEligibleN": int(first["InvalidEligibleN"]),
                "ObservedValidPointRaw": float(group.loc[group["OptionId"].eq(SCI01_OPTION_IDS[0]), "ReportedPointRaw"].iloc[0]),
                "AllInvalidDangerousProxyRaw": float(group.loc[group["OptionId"].eq(SCI01_OPTION_IDS[1]), "ReportedPointRaw"].iloc[0]),
                "ProjectedLabels": "|".join(sorted(labels)), "DecisionUsesDisplayedRounding": False})
    return pd.DataFrame(rows)


def render_projection_figure(projection: pd.DataFrame, output: Path) -> None:
    frame = projection.loc[
        projection["OptionId"].eq(SCI01_OPTION_IDS[2])
        & projection["ValidEligibleN"].eq(500)
        & projection["DangerousErrors"].eq(25)
    ].copy().sort_values("InvalidEligibleN")
    x = frame["InvalidEligibleN"].to_numpy(dtype=float) / frame["ValidEligibleN"].to_numpy(dtype=float)
    lower = frame["ReportedLowerRaw"].to_numpy(dtype=float)
    upper = frame["ReportedUpperRaw"].to_numpy(dtype=float)
    observed = frame["ReportedPointRaw"].to_numpy(dtype=float)
    worst_wilson = frame["GoverningWilsonUpper95Raw"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(8.6, 5.5))
    ax.fill_between(x, lower, upper, color="#93c5fd", alpha=.55, label="Partial-identification interval")
    ax.plot(x, observed, "o-", color="#0f766e", label="Observed-valid D/V")
    ax.plot(x, upper, "s-", color="#b45309", label="All-invalid-dangerous proxy")
    ax.plot(x, worst_wilson, "^-", color="#7f1d1d", label="Worst-case Wilson upper 95%")
    ax.axhline(.10, color="black", linestyle="--", linewidth=1.2, label="Strict threshold 0.10")
    ax.set(xlabel="Invalid / valid ratio (exact counts)", ylabel="Risk or upper confidence bound",
        title="Same data, different SCI-01 interpretation (V=500, D=25)")
    ax.legend(frameon=False, fontsize=8.5); ax.grid(alpha=.2); fig.tight_layout()
    fig.savefig(output / "sci01_estimand_projection.png", dpi=180); plt.close(fig)


def run_tests(output: Path) -> dict[str, object]:
    tests = [
        "tests/test_cmle_one_click_confirmatory_estimand.py",
        "tests/test_cmle_one_click_confirmatory_workbench.py",
        "tests/test_cmle_one_click_confirmatory_protocol.py",
        "tests/test_cmle_one_click_confirmatory_dependence.py",
        "tests/test_cmle_one_click_confirmatory_planning.py",
        "tests/test_cmle_one_click_confirmatory_gate.py",
        "tests/test_cmle_one_click_cognitive_interview.py",
        "tests/test_cmle_one_click_comprehension.py", "tests/test_cmle_one_click.py",
        "tests/test_cmle_one_click_archive.py", "tests/test_decision_stability.py",
        "tests/test_threshold_decision_integration.py",
    ]
    command = [sys.executable, "-m", "pytest", "-q", *tests]
    completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
    (output / "selected_tests_stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output / "selected_tests_stderr.txt").write_text(completed.stderr, encoding="utf-8")
    return {"passed": completed.returncode == 0, "returncode": completed.returncode,
        "command": command, "stdout_last_line": completed.stdout.strip().splitlines()[-1] if completed.stdout.strip() else ""}


def write_review(output: Path, result: dict[str, object]) -> None:
    (output / "SCI01_ESTIMAND_MEMO_CRITICAL_REVIEW.md").write_text(f"""# SCI-01 estimand memo — critical review

## Decision

The private software contract **{'passed' if result['contract_passed'] else 'failed'}**. No estimand was selected, substantive evidence remains unverified, and recruitment remains blocked.

## What changed

The three registered SCI-01 alternatives are now defined by target population, numerator/denominator, invalid-outcome treatment, identification requirements, allowed and forbidden claims, and downstream dependencies. The same 72 integer-count scenarios produce 216 explicitly labeled projections. There are {result['divergent_same_data_scenarios']} scenarios where the observed-valid gate passes while the all-invalid-dangerous proxy fails and the dual analysis is not robust.

## Critical limitation discovered

The frozen partial-identification surface does not contain a complete scheduled-slot denominator or a prospective invalidity-code-to-composite map. Therefore Option B is not directly calculable from that surface. Its displayed `(D+I)/(V+I)` value is intentionally labeled a diagnostic all-invalid-dangerous proxy, not the scheduled-slot estimand. Treating it as the final estimand would be denominator substitution.

## Decision integrity

The display order is not rank order. The repository supplies no recommendation, default, or automatic selection. Blank selection is valid but incomplete; hidden, unsupported, repository-owned, and post-outcome selections fail closed. Even a valid SCI-01 selection cannot authorize recruitment because 12 downstream decisions, substantive review, sample size, ethics, consent, and governance remain unresolved.

The seven-member offline ZIP reproduced byte-for-byte and passed manifest validation. Bundle SHA-256: `{result['bundle_sha256']}`. Human participants, confirmatory outcomes, selected sample size, and public UI remain zero/false.
""", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--output", type=Path, default=OUTPUT); args = parser.parse_args()
    plan = validate_registration()
    if args.output.exists():
        raise FileExistsError(f"Evidence output exists: {args.output}")
    args.output.mkdir(parents=True)
    surface = pd.read_csv(SURFACE)
    first = build_private_sci01_bundle(surface); second = build_private_sci01_bundle(surface)
    definitions, projection, impacts, selection = (first[k] for k in ("definitions", "projection", "impacts", "selection"))
    selection_audit = selection_attack_audit(); divergence = same_data_divergence_audit(projection)
    human = sci01_human_gate_status(first["selection_status"])
    html_result = validate_sci01_html(first["html"], first["content_sha256"])
    bundle_result = validate_private_sci01_bundle(first["bundle_bytes"])
    bundle_audit = pd.DataFrame([{"BundleSHA256": sha256_bytes(first["bundle_bytes"]),
        "RepeatedSHA256": sha256_bytes(second["bundle_bytes"]),
        "ByteIdentical": first["bundle_bytes"] == second["bundle_bytes"],
        "RegisteredMembers": len(SCI01_BUNDLE_MEMBERS), "ManifestRows": len(first["manifest"]),
        "Valid": bundle_result["valid"], "FailureCodes": "|".join(bundle_result["failure_codes"]),
        "Passed": bool(first["bundle_bytes"] == second["bundle_bytes"] and bundle_result["valid"]) }])
    parent_audit = pd.DataFrame([{"Path": relative, "ExpectedSHA256": expected,
        "ActualSHA256": sha256_file(ROOT / relative), "Passed": sha256_file(ROOT / relative) == expected}
        for relative, expected in plan["parent_evidence_sha256"].items()])
    formula_audit = pd.DataFrame([{"SourceRows": len(surface), "ProjectionRows": len(projection),
        "OptionsPerSourceMin": int(projection.groupby("SourceScenarioRow").size().min()),
        "OptionsPerSourceMax": int(projection.groupby("SourceScenarioRow").size().max()),
        "DisplayedRoundingUsed": bool(projection["DecisionUsesDisplayedRounding"].any()),
        "ScheduledSlotEstimandFalselyClaimedCalculated": bool(projection["FinalScheduledSlotEstimandCalculated"].any()),
        "Passed": bool(len(surface) == 72 and len(projection) == 216
            and projection.groupby("SourceScenarioRow").size().eq(3).all()
            and not projection["DecisionUsesDisplayedRounding"].any()
            and not projection["FinalScheduledSlotEstimandCalculated"].any())}])
    tables = {"sci01_estimand_definitions.csv": definitions, "sci01_same_data_projection.csv": projection,
        "sci01_downstream_decision_impacts.csv": impacts, "sci01_selection_template.csv": selection,
        "sci01_selection_validator_audit.csv": selection_audit, "sci01_same_data_divergence_audit.csv": divergence,
        "sci01_parent_identity_audit.csv": parent_audit, "sci01_formula_contract_audit.csv": formula_audit,
        "sci01_bundle_contract_audit.csv": bundle_audit, "sci01_human_gate_status.csv": human,
        "sci01_bundle_manifest.csv": first["manifest"]}
    for name, frame in tables.items(): write_csv(frame, args.output / name)
    (args.output / "sci01_estimand_comparison.html").write_text(first["html"], encoding="utf-8")
    (args.output / "private_sci01_estimand_memo.zip").write_bytes(first["bundle_bytes"])
    render_projection_figure(projection, args.output)
    tests = run_tests(args.output)
    human_row = human.iloc[0]
    gates = {"parent_identity_passed": bool(parent_audit["Passed"].all()),
        "definition_contract_passed": bool(len(definitions) == 3 and not definitions["DisplayOrderIsRank"].any()
            and definitions["RepositoryRecommendation"].eq("none").all()),
        "formula_contract_passed": bool(formula_audit["Passed"].all()),
        "same_data_divergence_demonstrated": bool(len(divergence) > 0),
        "downstream_map_passed": bool(len(impacts) == 36 and not impacts["DecisionResolved"].any()),
        "selection_fail_closed_passed": bool(selection_audit["Passed"].all()),
        "html_passed": bool(html_result["valid"]), "bundle_passed": bool(bundle_audit["Passed"].all()),
        "human_gate_passed": bool(int(human_row["HumanParticipants"]) == 0 and not bool(human_row["SCI01Complete"])
            and not bool(human_row["RecruitmentReady"]) and not bool(human_row["SampleSizeSelected"])
            and not bool(human_row["ConfirmatoryOutcomesAvailable"]) and not bool(human_row["PublicSurfaceEnabled"])),
        "tests_passed": bool(tests["passed"])}
    contract = all(gates.values())
    result = {**gates, "contract_passed": contract, "source_scenarios": len(surface),
        "projection_rows": len(projection), "divergent_same_data_scenarios": len(divergence),
        "options_selected": 0, "sci01_complete": False, "substantive_evidence_verified": False,
        "recruitment_ready": False, "sample_size_selected": False, "human_participants": 0,
        "confirmatory_outcomes_available": False, "public_surface_enabled": False,
        "bundle_members": len(SCI01_BUNDLE_MEMBERS), "bundle_sha256": sha256_bytes(first["bundle_bytes"]),
        "content_sha256": first["content_sha256"], "selected_tests": tests}
    write_review(args.output, result)
    files = sorted(path.relative_to(args.output).as_posix() for path in args.output.rglob("*") if path.is_file() and path.name != "decision.json")
    decision = {"study_id": plan["study_id"], "plan_sha256": sha256_file(PLAN),
        "contract_passed": contract, "contract_interpretation": "private_nonrecommending_SCI01_support_selection_blocked",
        "implementation_sha256": {
            "mfrm_app/cmle_one_click_confirmatory_estimand.py": sha256_file(ROOT / "mfrm_app/cmle_one_click_confirmatory_estimand.py"),
            "tests/test_cmle_one_click_confirmatory_estimand.py": sha256_file(ROOT / "tests/test_cmle_one_click_confirmatory_estimand.py"),
            "validation/cmle_one_click_confirmatory_sci01_estimand_memo.py": sha256_file(Path(__file__))},
        "results": result, "output_sha256": {name: sha256_file(args.output / name) for name in files},
        "interpretation": {"estimand_selection": "none", "option_B": "diagnostic_proxy_not_final_scheduled_slot_estimand",
            "rounding": "never_drives_decision", "human_study": "not_ready", "public_ui": "withheld"}}
    (args.output / "decision.json").write_text(json.dumps(decision, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(decision, indent=2, ensure_ascii=False, allow_nan=False))
    if not contract: raise SystemExit(1)


if __name__ == "__main__":
    main()
