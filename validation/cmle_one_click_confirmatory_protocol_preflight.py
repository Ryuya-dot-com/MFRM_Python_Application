#!/usr/bin/env python3
"""Build the registered fail-closed private confirmatory protocol preflight."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

_MPL_CONFIG = Path(tempfile.gettempdir()) / "mfrm_app_matplotlib"
_MPL_CONFIG.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CONFIG))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle_one_click_confirmatory_gate import (  # noqa: E402
    BLOCKED_PRIMARY_CASES,
    INSTRUMENT_VERSION,
)
from mfrm_app.cmle_one_click_confirmatory_protocol import (  # noqa: E402
    LEDGER_COLUMNS,
    PRIVATE_BUNDLE_MEMBERS,
    REGISTERED_INVALIDITY_CODES,
    build_private_protocol_preflight_bundle,
    compute_protocol_content_identity,
    invalidity_partial_identification,
    protocol_preflight_human_gate_status,
    validate_attempt_denominator_ledger,
    validate_private_protocol_preflight_bundle,
)


PLAN = ROOT / "validation/cmle_one_click_confirmatory_protocol_preflight_plan_20260810.json"
OUTPUT = ROOT / "validation/cmle_one_click_confirmatory_protocol_preflight_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


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


def validate_registration() -> dict[str, object]:
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    if plan.get("study_id") != "cmle_one_click_confirmatory_protocol_preflight_v1":
        raise ValueError("Unexpected protocol-preflight study identity.")
    mismatches = [
        relative
        for relative, expected in plan["parent_identity"].items()
        if not (ROOT / relative).is_file()
        or sha256_file(ROOT / relative) != expected
    ]
    if mismatches:
        raise ValueError(f"Protocol-preflight parent identity failed: {mismatches}")
    sensitivity = plan["partial_identification_sensitivity"]
    if sensitivity["valid_eligible_n"] != [100, 200, 300, 500]:
        raise ValueError("Unexpected registered valid-n surface.")
    if sensitivity["observed_dangerous_error_rates"] != [0.0, 0.025, 0.05]:
        raise ValueError("Unexpected registered observed-error surface.")
    if sensitivity["invalid_to_valid_ratios"] != [0.0, 0.01, 0.025, 0.05, 0.1, 0.2]:
        raise ValueError("Unexpected registered invalidity surface.")
    if tuple(sorted(plan["private_bundle_contract"]["members"])) != PRIVATE_BUNDLE_MEMBERS:
        raise ValueError("Registered private-bundle membership changed.")
    if tuple(plan["invalidity_reason_contract"]["registered_codes"]) != REGISTERED_INVALIDITY_CODES:
        raise ValueError("Registered invalidity vocabulary changed.")
    return plan


def build_registration_audit(plan: dict[str, object]) -> pd.DataFrame:
    rows = []
    for relative, expected in plan["parent_identity"].items():
        path = ROOT / relative
        actual = sha256_file(path) if path.is_file() else ""
        rows.append(
            {
                "Artifact": relative,
                "ExpectedSHA256": expected,
                "ActualSHA256": actual,
                "Passed": bool(actual == expected),
            }
        )
    return pd.DataFrame(rows)


def build_register_audit(
    register: pd.DataFrame, readiness: dict[str, object]
) -> pd.DataFrame:
    fixed = register.loc[~register["Blocking"]]
    blocking = register.loc[register["Blocking"]]
    return pd.DataFrame(
        [
            {
                "DecisionRows": len(register),
                "UniqueDecisionIds": register["DecisionId"].nunique(),
                "FixedRows": len(fixed),
                "BlockingRows": len(blocking),
                "UnresolvedBlockingRows": int(
                    blocking["Status"].eq("blocked_unresolved").sum()
                ),
                "BlankBlockingEvidenceRows": int(
                    blocking["EvidenceReference"].astype(str).str.strip().eq("").sum()
                ),
                "RepositorySelfApprovalRows": int(
                    register["RepositoryMaySelfApprove"].sum()
                ),
                "RecruitmentReady": bool(readiness["RecruitmentReady"]),
                "ReadinessStatus": readiness["Status"],
                "Passed": bool(
                    len(register) == 18
                    and register["DecisionId"].nunique() == 18
                    and len(fixed) == 5
                    and len(blocking) == 13
                    and blocking["Status"].eq("blocked_unresolved").all()
                    and blocking["EvidenceReference"].astype(str).str.strip().eq("").all()
                    and not register["RepositoryMaySelfApprove"].any()
                    and not bool(readiness["RecruitmentReady"])
                    and readiness["Status"] == "blocked_pre_recruitment"
                ),
            }
        ]
    )


def _valid_ledger_row(slot: str = "CF-EN-0001") -> dict[str, object]:
    return {
        "StudyRecordId": "SYNTH-REC-0001",
        "ConfirmatorySlotId": slot,
        "InstrumentVersion": INSTRUMENT_VERSION,
        "InstrumentContentSHA256": "a" * 64,
        "Language": "en",
        "SiteCode": "SITE-01",
        "ModeratorCode": "MOD-01",
        "CollectionWave": "WAVE-01",
        "BlockedMechanism": BLOCKED_PRIMARY_CASES[0],
        "ConsentDisposition": "consented",
        "EligibilityDisposition": "eligible",
        "SessionStarted": True,
        "ResponseLockCompleted": True,
        "CompletionDisposition": "complete",
        "InvalidityReasonCode": "none",
        "PrimaryAnalysisDisposition": "included",
        "ExclusionDecisionStage": "not_applicable",
        "AccessibilityStratumCode": "ACCESS-01",
        "DeviceStratumCode": "DEVICE-01",
        "ProtocolDeviationCode": "none",
        "PIIReviewed": True,
        "AuditReviewed": True,
    }


def build_ledger_validator_audit() -> pd.DataFrame:
    empty = pd.DataFrame(columns=LEDGER_COLUMNS)
    valid = pd.DataFrame([_valid_ledger_row()], columns=LEDGER_COLUMNS)
    pii = valid.copy()
    pii["ParticipantEmail"] = "synthetic@example.invalid"
    unknown_reason = valid.copy()
    unknown_reason.loc[0, "InvalidityReasonCode"] = "researcher_discretion"
    post_outcome = valid.copy()
    post_outcome.loc[0, "PrimaryAnalysisDisposition"] = "excluded"
    post_outcome.loc[0, "InvalidityReasonCode"] = "incomplete_primary_items"
    post_outcome.loc[0, "ExclusionDecisionStage"] = "post_outcome"
    duplicate = pd.concat([valid, valid], ignore_index=True)
    unreviewed = valid.copy()
    unreviewed.loc[0, "AuditReviewed"] = False
    cases = (
        ("zero_row_template", empty, {}, True, None),
        ("expected_slot_missing", empty, {"expected_slot_ids": ["CF-EN-0001"]}, False, "missing_expected_slots"),
        (
            "valid_synthetic_row",
            valid,
            {
                "expected_slot_ids": ["CF-EN-0001"],
                "instrument_version": INSTRUMENT_VERSION,
                "instrument_content_sha256": "a" * 64,
            },
            True,
            None,
        ),
        ("banned_email_column", pii, {}, False, "banned_identifier_or_free_text_column"),
        ("unregistered_reason", unknown_reason, {}, False, "unregistered_invalidity_reason"),
        ("post_outcome_exclusion", post_outcome, {}, False, "post_outcome_or_unknown_exclusion_stage"),
        ("duplicate_records", duplicate, {}, False, "duplicate_study_record_id"),
        ("audit_not_reviewed", unreviewed, {}, False, "audit_review_incomplete"),
    )
    rows = []
    for case, frame, kwargs, expected_valid, required_failure in cases:
        result = validate_attempt_denominator_ledger(frame, **kwargs)
        required_present = required_failure is None or required_failure in result[
            "failure_codes"
        ]
        rows.append(
            {
                "Case": case,
                "SyntheticRows": len(frame),
                "ExpectedValid": expected_valid,
                "ActualValid": bool(result["valid"]),
                "RequiredFailureCode": required_failure or "",
                "ActualFailureCodes": "|".join(result["failure_codes"]),
                "Passed": bool(result["valid"] == expected_valid and required_present),
            }
        )
    return pd.DataFrame(rows)


def build_partial_identification_audit(
    sensitivity: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for index, source in sensitivity.iterrows():
        valid_n = int(source["ValidEligibleN"])
        rate = float(source["ObservedDangerousErrorRateRegistered"])
        ratio = float(source["InvalidToValidRatioRegistered"])
        expected_errors = int(math.floor(valid_n * rate + 0.5))
        expected_invalid = int(math.ceil(valid_n * ratio))
        recomputed = invalidity_partial_identification(
            expected_errors, valid_n, expected_invalid
        )
        count_passed = bool(
            int(source["DangerousErrors"]) == expected_errors
            and int(source["InvalidEligibleN"]) == expected_invalid
        )
        numeric_columns = (
            "ObservedValidDangerousProportionRaw",
            "AllInvalidSafeLowerRiskRaw",
            "AllInvalidDangerousUpperRiskRaw",
            "ObservedValidWilsonUpper95Raw",
            "AllInvalidSafeWilsonUpper95Raw",
            "AllInvalidDangerousWilsonUpper95Raw",
        )
        numeric_passed = all(
            math.isclose(
                float(source[column]),
                float(recomputed[column]),
                rel_tol=0.0,
                abs_tol=2e-15,
            )
            for column in numeric_columns
        )
        decision_passed = bool(
            bool(source["ObservedValidWilsonPass"])
            == bool(recomputed["ObservedValidWilsonPass"])
            and bool(source["RobustToInvalidOutcomes"])
            == bool(recomputed["RobustToInvalidOutcomes"])
            and not bool(source["DecisionUsesDisplayedRounding"])
            and not bool(source["SampleSizeSelected"])
        )
        rows.append(
            {
                "SourceRow": index,
                "ValidEligibleN": valid_n,
                "ObservedDangerousErrorRateRegistered": rate,
                "InvalidToValidRatioRegistered": ratio,
                "CountFormulaPassed": count_passed,
                "RawNumericFormulaPassed": numeric_passed,
                "StrictDecisionPassed": decision_passed,
                "Passed": count_passed and numeric_passed and decision_passed,
            }
        )
    formula = pd.DataFrame(rows)
    monotone_rows = []
    groups = sensitivity.groupby(
        ["ValidEligibleN", "ObservedDangerousErrorRateRegistered"], sort=False
    )
    for (valid_n, rate), frame in groups:
        ordered = frame.sort_values("InvalidEligibleN")
        risk_monotone = bool(
            ordered["AllInvalidDangerousUpperRiskRaw"].is_monotonic_increasing
        )
        wilson_monotone = bool(
            ordered["AllInvalidDangerousWilsonUpper95Raw"].is_monotonic_increasing
        )
        robust_implies_observed = bool(
            (
                ~ordered["RobustToInvalidOutcomes"]
                | ordered["ObservedValidWilsonPass"]
            ).all()
        )
        monotone_rows.append(
            {
                "ValidEligibleN": valid_n,
                "ObservedDangerousErrorRateRegistered": rate,
                "Rows": len(ordered),
                "WorstCaseRiskMonotone": risk_monotone,
                "WorstCaseWilsonMonotone": wilson_monotone,
                "RobustImpliesObservedPass": robust_implies_observed,
                "Passed": risk_monotone and wilson_monotone and robust_implies_observed,
            }
        )
    return formula, pd.DataFrame(monotone_rows)


def _tamper_start_here(bundle_bytes: bytes) -> bytes:
    with zipfile.ZipFile(io.BytesIO(bundle_bytes), "r") as source:
        payloads = {name: source.read(name) for name in source.namelist()}
    payloads["START_HERE.md"] = payloads["START_HERE.md"].replace(
        b"BLOCKED", b"READY  ", 1
    )
    buffer = io.BytesIO()
    with zipfile.ZipFile(
        buffer, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as archive:
        for name in sorted(payloads):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            info.create_system = 3
            archive.writestr(info, payloads[name])
    return buffer.getvalue()


def build_bundle_audit(
    first: dict[str, object], second: dict[str, object]
) -> pd.DataFrame:
    first_bytes = first["bundle_bytes"]
    second_bytes = second["bundle_bytes"]
    first_validation = validate_private_protocol_preflight_bundle(first_bytes)
    tampered_validation = validate_private_protocol_preflight_bundle(
        _tamper_start_here(first_bytes)
    )
    return pd.DataFrame(
        [
            {
                "BundleSHA256": sha256_bytes(first_bytes),
                "RepeatedBundleSHA256": sha256_bytes(second_bytes),
                "ByteIdenticalOnRepeat": first_bytes == second_bytes,
                "RegisteredMemberCount": len(PRIVATE_BUNDLE_MEMBERS),
                "ManifestMemberCount": len(first["bundle_manifest"]),
                "HumanRowsInManifest": int(first["bundle_manifest"]["HumanRows"].sum()),
                "DirectIdentifierMembers": int(
                    first["bundle_manifest"]["ContainsDirectIdentifiers"].sum()
                ),
                "OriginalValid": bool(first_validation["valid"]),
                "TamperedValid": bool(tampered_validation["valid"]),
                "TamperedFailureCodes": "|".join(
                    tampered_validation["failure_codes"]
                ),
                "Passed": bool(
                    first_bytes == second_bytes
                    and first_validation["valid"]
                    and not tampered_validation["valid"]
                    and "manifest_hash_mismatch"
                    in tampered_validation["failure_codes"]
                    and "first_read_not_fail_closed"
                    in tampered_validation["failure_codes"]
                    and len(first["bundle_manifest"])
                    == len(PRIVATE_BUNDLE_MEMBERS) - 1
                    and int(first["bundle_manifest"]["HumanRows"].sum()) == 0
                    and not first["bundle_manifest"][
                        "ContainsDirectIdentifiers"
                    ].any()
                ),
            }
        ]
    )


def build_content_identity_mutation_audit(
    artifacts: dict[str, object]
) -> pd.DataFrame:
    changed = artifacts["decision_register"].copy()
    changed.loc[
        changed["DecisionId"].eq("SCI-01"), "Decision"
    ] += " synthetic mutation"
    changed_sha = compute_protocol_content_identity(
        changed,
        artifacts["invalidity_reason_codebook"],
        artifacts["ledger_data_dictionary"],
        artifacts["partial_identification_sensitivity"],
    )
    original = str(artifacts["protocol_content_sha256"])
    return pd.DataFrame(
        [
            {
                "OriginalProtocolContentSHA256": original,
                "MutatedProtocolContentSHA256": changed_sha,
                "MutationChangedIdentity": original != changed_sha,
                "Passed": original != changed_sha,
            }
        ]
    )


def render_figures(
    sensitivity: pd.DataFrame, preflight: pd.DataFrame, output: Path
) -> None:
    figure_data = sensitivity.loc[
        sensitivity["ObservedDangerousErrorRateRegistered"].eq(0.05)
    ]
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    for valid_n, frame in figure_data.groupby("ValidEligibleN", sort=True):
        ax.plot(
            frame["ActualInvalidToValidRatioRaw"],
            frame["AllInvalidDangerousWilsonUpper95Raw"],
            marker="o",
            linewidth=1.5,
            label=f"valid n={valid_n}",
        )
    ax.axhline(
        0.10,
        color="#b91c1c",
        linewidth=1.2,
        linestyle="--",
        label="strict decision threshold 0.10",
    )
    ax.set(
        xlabel="Invalid eligible / valid eligible ratio",
        ylabel="All-invalid-dangerous one-sided 95% Wilson upper bound",
        title="Observed 5% errors may not be robust to invalid outcomes",
        ylim=(0, max(0.24, float(figure_data["AllInvalidDangerousWilsonUpper95Raw"].max()) + 0.01)),
    )
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "partial_identification_invalidity_sensitivity.png", dpi=180)
    plt.close(fig)

    colors = {
        "blocked": "#b91c1c",
        "withheld": "#7c3aed",
        "fixed": "#2563eb",
        "template_ready": "#0f766e",
        "sensitivity_ready": "#0f766e",
    }
    ordered = preflight.sort_values("CardOrder", ascending=False)
    fig, ax = plt.subplots(figsize=(10, 5.6))
    y = np.arange(len(ordered))
    ax.barh(
        y,
        np.ones(len(ordered)),
        color=[colors.get(str(value), "#6b7280") for value in ordered["Status"]],
        alpha=0.9,
    )
    ax.set_yticks(y, ordered["CardId"])
    for index, (_, row) in enumerate(ordered.iterrows()):
        ax.text(
            0.02,
            index,
            f"{row['Status']}: {row['FirstRead']}",
            va="center",
            ha="left",
            color="white",
            fontsize=8.5,
        )
    ax.set(
        xlim=(0, 1),
        xticks=[],
        title="Private pre-recruitment first read: overall status remains blocked",
    )
    for spine in ("top", "right", "bottom"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(output / "pre_recruitment_readiness_status.png", dpi=180)
    plt.close(fig)


def run_tests(output: Path) -> dict[str, object]:
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-q",
        "tests/test_cmle_one_click_confirmatory_protocol.py",
        "tests/test_cmle_one_click_confirmatory_dependence.py",
        "tests/test_cmle_one_click_confirmatory_planning.py",
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
    (output / "selected_tests_stdout.txt").write_text(
        completed.stdout, encoding="utf-8"
    )
    (output / "selected_tests_stderr.txt").write_text(
        completed.stderr, encoding="utf-8"
    )
    return {
        "passed": completed.returncode == 0,
        "returncode": completed.returncode,
        "command": command,
    }


def write_review(output: Path, results: dict[str, object]) -> None:
    text = f"""# Private confirmatory protocol preflight critical review

## Decision

The fail-closed software contract **{'passed' if results['contract_passed'] else 'failed'}**. The human study is **not ready**: 13 blocking protocol/governance decisions remain unresolved, ethics approval is unavailable, recruitment is unauthorized, and no sample size is selected.

## Invalid outcomes change what can be claimed

The retained 72-row partial-identification surface separates the observed-valid error rate from the range obtained when invalid eligible outcomes are all safe versus all dangerous. With valid n=500 and an observed 5% dangerous-error rate, the registered grid remains worst-case robust at invalid/valid ratio `{results['n500_p05_max_robust_ratio']:.3f}`, but not at `{results['n500_p05_first_nonrobust_ratio']:.3f}`. At the latter ratio the all-invalid-dangerous Wilson upper bound is `{results['n500_p05_first_nonrobust_upper']:.6f}`. This does not identify MNAR; it shows exactly where an observed-valid pass ceases to survive the registered worst case.

## Denominators and exclusions fail closed

The exact-schema zero-row ledger retains one row per frozen slot after collection begins and rejects direct identifier/free-text columns, missing or duplicate slots, unknown reason codes, post-outcome exclusion, and incomplete PII/audit attestations. The categorical codebook retains withdrawn and invalid records in the attempted denominator. These checks do not prove de-identification, regulatory compliance, or unbiased missingness.

## Private one-click bundle

The deterministic ZIP has nine exact members, zero human rows, fixed metadata, and SHA-256 validation. Repeated generation produced identical bytes; a synthetic change from `BLOCKED` to `READY` failed both the first-read and manifest-hash checks. The archive SHA-256 is `{results['private_bundle_sha256']}`. It remains private and is not connected to Streamlit.

## What external owners must still decide

- target-population estimand and handling of invalid/withdrawn outcomes;
- minimum valid n, recruitment allowance, multiplicity, and cluster-aware analysis;
- mechanism-specific protection, language objective, and accessibility/device/proficiency strata;
- stopping, blinded exclusions, ethics, consent, privacy, retention, and incident response.

Passing this contract means the repository refuses to manufacture readiness. Human participants remain zero, recruitment remains blocked, and the public UI stays withheld.
"""
    (output / "CONFIRMATORY_PROTOCOL_PREFLIGHT_CRITICAL_REVIEW.md").write_text(
        text, encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    plan = validate_registration()
    if args.output.exists():
        raise FileExistsError(f"Evidence output exists: {args.output}")
    args.output.mkdir(parents=True)

    artifacts = build_private_protocol_preflight_bundle()
    repeated = build_private_protocol_preflight_bundle()
    register = artifacts["decision_register"]
    codebook = artifacts["invalidity_reason_codebook"]
    ledger = artifacts["ledger_template"]
    dictionary = artifacts["ledger_data_dictionary"]
    sensitivity = artifacts["partial_identification_sensitivity"]
    ethics = artifacts["ethics_template"]
    readiness = artifacts["readiness"]
    preflight = artifacts["preflight_summary"]
    manifest = artifacts["bundle_manifest"]
    human_gate = protocol_preflight_human_gate_status(readiness)

    registration_audit = build_registration_audit(plan)
    register_audit = build_register_audit(register, readiness)
    ledger_audit = build_ledger_validator_audit()
    formula_audit, monotonicity_audit = build_partial_identification_audit(
        sensitivity
    )
    bundle_audit = build_bundle_audit(artifacts, repeated)
    identity_mutation_audit = build_content_identity_mutation_audit(artifacts)

    outputs = {
        "protocol_decision_register.csv": register,
        "pre_recruitment_readiness_summary.csv": preflight,
        "attempt_denominator_ledger_template.csv": ledger,
        "attempt_denominator_data_dictionary.csv": dictionary,
        "invalidity_reason_codebook.csv": codebook,
        "ethics_governance_evidence_template.csv": ethics,
        "partial_identification_sensitivity.csv": sensitivity,
        "bundle_manifest.csv": manifest,
        "human_gate_status.csv": human_gate,
        "registration_identity_audit.csv": registration_audit,
        "protocol_decision_register_audit.csv": register_audit,
        "ledger_validator_synthetic_audit.csv": ledger_audit,
        "partial_identification_formula_audit.csv": formula_audit,
        "partial_identification_monotonicity_audit.csv": monotonicity_audit,
        "private_bundle_contract_audit.csv": bundle_audit,
        "protocol_content_identity_mutation_audit.csv": identity_mutation_audit,
    }
    for name, frame in outputs.items():
        write_csv(frame, args.output / name)
    (args.output / "private_confirmatory_protocol_preflight.zip").write_bytes(
        artifacts["bundle_bytes"]
    )
    render_figures(sensitivity, preflight, args.output)
    tests = run_tests(args.output)

    human = human_gate.iloc[0]
    gates = {
        "identity_passed": bool(registration_audit["Passed"].all()),
        "register_fail_closed_passed": bool(register_audit["Passed"].all()),
        "ledger_contract_passed": bool(
            len(ledger) == 0
            and tuple(ledger.columns) == LEDGER_COLUMNS
            and len(dictionary) == len(LEDGER_COLUMNS)
            and len(codebook) == len(REGISTERED_INVALIDITY_CODES)
            and ledger_audit["Passed"].all()
        ),
        "partial_identification_passed": bool(
            len(sensitivity) == 4 * 3 * 6
            and formula_audit["Passed"].all()
            and monotonicity_audit["Passed"].all()
            and not sensitivity["SampleSizeSelected"].any()
            and not sensitivity["MinimumValidNRegistered"].any()
            and not sensitivity["PlannedRecruitmentNSelected"].any()
        ),
        "bundle_contract_passed": bool(bundle_audit["Passed"].all()),
        "content_identity_passed": bool(identity_mutation_audit["Passed"].all()),
        "human_gate_passed": bool(
            int(human["HumanParticipants"]) == 0
            and not bool(human["EthicsApprovalAvailable"])
            and not bool(human["RecruitmentAuthorized"])
            and not bool(human["RecruitmentReady"])
            and not bool(human["MinimumValidPerCellRegistered"])
            and not bool(human["PlannedRecruitmentNSelected"])
            and not bool(human["ConfirmatoryResultAvailable"])
            and not bool(human["LanguageEquivalenceEstablished"])
            and not bool(human["PublicSurfaceEnabled"])
        ),
        "tests_passed": bool(tests["passed"]),
    }
    contract_passed = all(gates.values())

    n500_p05 = sensitivity.loc[
        sensitivity["ValidEligibleN"].eq(500)
        & sensitivity["ObservedDangerousErrorRateRegistered"].eq(0.05)
    ].sort_values("InvalidToValidRatioRegistered")
    robust = n500_p05.loc[n500_p05["RobustToInvalidOutcomes"]]
    nonrobust = n500_p05.loc[~n500_p05["RobustToInvalidOutcomes"]]
    first_nonrobust = nonrobust.iloc[0]
    results = {
        **gates,
        "contract_passed": contract_passed,
        "decision_rows": len(register),
        "blocking_decision_rows": int(register["Blocking"].sum()),
        "partial_identification_rows": len(sensitivity),
        "ledger_validator_cases": len(ledger_audit),
        "private_bundle_members": len(PRIVATE_BUNDLE_MEMBERS),
        "private_bundle_sha256": sha256_bytes(artifacts["bundle_bytes"]),
        "protocol_content_sha256": artifacts["protocol_content_sha256"],
        "n500_p05_max_robust_ratio": float(
            robust["InvalidToValidRatioRegistered"].max()
        ),
        "n500_p05_first_nonrobust_ratio": float(
            first_nonrobust["InvalidToValidRatioRegistered"]
        ),
        "n500_p05_first_nonrobust_upper": float(
            first_nonrobust["AllInvalidDangerousWilsonUpper95Raw"]
        ),
        "human_participants": 0,
        "ethics_approval_available": False,
        "recruitment_authorized": False,
        "recruitment_ready": False,
        "minimum_valid_per_cell_registered": False,
        "planned_recruitment_n_selected": False,
        "confirmatory_result_available": False,
        "public_surface_enabled": False,
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
            "contract_interpretation": "fail_closed_private_pre_recruitment_bundle_human_study_not_ready",
            "implementation_sha256": {
                "mfrm_app/cmle_one_click_confirmatory_protocol.py": sha256_file(
                    ROOT / "mfrm_app/cmle_one_click_confirmatory_protocol.py"
                ),
                "tests/test_cmle_one_click_confirmatory_protocol.py": sha256_file(
                    ROOT / "tests/test_cmle_one_click_confirmatory_protocol.py"
                ),
                "validation/cmle_one_click_confirmatory_protocol_preflight.py": sha256_file(
                    Path(__file__)
                ),
            },
            "results": results,
            "output_sha256": {
                name: sha256_file(args.output / name) for name in output_files
            },
            "interpretation": {
                "software_contract": "passed means fail-closed behavior is correct",
                "human_study": "not ready",
                "partial_identification": "worst-case sensitivity, not MNAR identification",
                "privacy": "schema checks do not prove de-identification or compliance",
                "ethics": "external independent approval absent",
                "sample_size": "not selected",
                "human_data": "none",
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
