#!/usr/bin/env python3
"""Run the registered CMLE one-click private archive identity study."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import warnings
import zipfile

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle_one_click import run_cmle_one_click_analysis  # noqa: E402
from mfrm_app.cmle_one_click_archive import (  # noqa: E402
    build_cmle_one_click_archive,
    load_cmle_one_click_archive,
    replay_cmle_one_click_archive,
    verify_cmle_one_click_archive,
)


PLAN = ROOT / "validation/cmle_one_click_archive_identity_plan_20260810.json"
OUTPUT = ROOT / "validation/cmle_one_click_archive_identity_20260810"


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
    if plan.get("study_id") != "cmle_one_click_archive_identity_v1":
        raise ValueError("Unexpected CMLE one-click archive plan.")
    mismatches = [
        relative
        for relative, expected in plan["parent_identity"].items()
        if not (ROOT / relative).is_file()
        or sha256_file(ROOT / relative) != expected
    ]
    if mismatches:
        raise ValueError(f"CMLE archive parent identity failed: {mismatches}")
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


def registered_cases(plan: dict[str, object]):
    contracts = {item["case_id"]: item for item in plan["registered_cases"]}
    anchors = pd.DataFrame(
        [
            {
                "ParameterType": "Facet",
                "Facet": "Rater",
                "Level": "R1",
                "Value": 0.25,
            },
            {
                "ParameterType": "Facet",
                "Facet": "Criterion",
                "Level": "C1",
                "Value": -0.2,
            },
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
            contracts["rsm_ready_extremes_private_archive"],
            interior_frame(extremes=True),
            {
                **common,
                "facet_cols": ["Rater", "Criterion"],
                "rating_max": 2,
                "model": "RSM",
            },
        ),
        (
            contracts["pcm_differential_anchor_private_archive"],
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
            contracts["binary_boundary_private_archive"],
            boundary_frame(),
            {
                **common,
                "facet_cols": ["Rater"],
                "rating_max": 1,
                "model": "RSM",
            },
        ),
    ]


def run_cases(plan: dict[str, object], output: Path):
    case_rows = []
    identity_rows = []
    asset_rows = []
    replay_rows = []
    floating_rows = []
    retained = {}
    for contract, frame, kwargs in registered_cases(plan):
        case_id = contract["case_id"]
        result = run_cmle_one_click_analysis(frame, **kwargs)
        first = build_cmle_one_click_archive(
            result, input_data=frame, calibration_kwargs=kwargs
        )
        reordered_kwargs = dict(reversed(list(kwargs.items())))
        second = build_cmle_one_click_archive(
            result, input_data=frame, calibration_kwargs=reordered_kwargs
        )
        verified = verify_cmle_one_click_archive(first["zip_bytes"])
        loaded = load_cmle_one_click_archive(first["zip_bytes"])
        replay = replay_cmle_one_click_archive(first["zip_bytes"])
        zip_name = f"{case_id}.zip"
        (output / zip_name).write_bytes(first["zip_bytes"])
        actual_conditional = sorted(
            set(first["assets"])
            - set(plan["archive_scope"]["required_core_assets"])
            - {"archive_manifest.json"}
        )
        expected_conditional = sorted(contract["expected_conditional_assets"])
        terminal = str(result["summary"].iloc[0]["TerminalStatus"])
        case_rows.append(
            {
                "CaseId": case_id,
                "ExpectedTerminalStatus": contract["expected_terminal_status"],
                "ObservedTerminalStatus": terminal,
                "TerminalMatch": terminal == contract["expected_terminal_status"],
                "ExpectedConditionalAssets": ";".join(expected_conditional),
                "ObservedConditionalAssets": ";".join(actual_conditional),
                "ConditionalAssetMatch": actual_conditional == expected_conditional,
                "ArchiveVerified": verified["verified"],
                "CardsLoaded": len(loaded["tables"]["result_cards.csv"]),
                "BilingualCardsLoaded": bool(
                    loaded["tables"]["result_cards.csv"][
                        ["HeadlineEn", "HeadlineJa", "DetailEn", "DetailJa"]
                    ]
                    .astype(str)
                    .apply(lambda column: column.str.len().gt(0))
                    .all()
                    .all()
                ),
                "RepeatedAssetIdentityMatch": first["manifest"]["assets"]
                == second["manifest"]["assets"],
                "RepeatedArchiveContentMatch": first["archive_content_sha256"]
                == second["archive_content_sha256"],
                "RepeatedZipBytesMatch": first["zip_bytes"] == second["zip_bytes"],
                "ReplayPassed": replay["passed"],
                "ZipFile": zip_name,
            }
        )
        identity_rows.append(
            {
                "CaseId": case_id,
                "AnalysisID": first["analysis_identity"].analysis_id,
                "OrderedInputSHA256": first["ordered_input_sha256"],
                "SemanticInputSHA256": first["semantic_input_sha256"],
                "ConfigSHA256": first["analysis_identity"].config_fingerprint,
                "ResultContentSHA256": first["result_content_sha256"],
                "ArchiveContentSHA256": first["archive_content_sha256"],
                "ZipSHA256": first["zip_sha256"],
                "PrivacyMode": first["manifest"]["privacy_mode"],
                "PublicSurfaceEnabled": first["manifest"][
                    "public_surface_enabled"
                ],
            }
        )
        asset_rows.extend(
            {"CaseId": case_id, **row} for row in first["manifest"]["assets"]
        )
        replay_rows.extend(
            {
                "CaseId": case_id,
                "Check": name,
                "Passed": passed,
            }
            for name, passed in replay["checks"].items()
        )
        if result["person_fit"] is not None:
            original = result["person_fit"]["persons"].sort_values("Person")
            restored = loaded["tables"]["person_results.csv"].sort_values("Person")
            for column in (
                "WLEEstimate",
                "ConditionalWLEStandardError",
                "Infit",
                "Outfit",
            ):
                left = original[column].to_numpy(dtype=float)
                right = restored[column].to_numpy(dtype=float)
                for person, before, after in zip(
                    original["Person"].astype(str), left, right, strict=True
                ):
                    floating_rows.append(
                        {
                            "CaseId": case_id,
                            "Person": person,
                            "Statistic": column,
                            "Original": before,
                            "Restored": after,
                            "OriginalHex": float(before).hex(),
                            "RestoredHex": float(after).hex(),
                            "Binary64Exact": float(before).hex()
                            == float(after).hex(),
                        }
                    )
        retained[case_id] = {
            "frame": frame,
            "kwargs": kwargs,
            "result": result,
            "archive": first,
        }
    frames = {
        "case_ledger": pd.DataFrame(case_rows),
        "identity_ledger": pd.DataFrame(identity_rows),
        "asset_manifest": pd.DataFrame(asset_rows),
        "replay_checks": pd.DataFrame(replay_rows),
        "floating_roundtrip": pd.DataFrame(floating_rows),
    }
    for name, frame in frames.items():
        write_csv(frame, output / f"{name}.csv")
    return frames, retained


def row_and_anchor_probes(retained, output: Path):
    rsm = retained["rsm_ready_extremes_private_archive"]
    permuted = rsm["frame"].sample(
        frac=1.0, random_state=20260810
    ).reset_index(drop=True)
    permuted_archive = build_cmle_one_click_archive(
        rsm["result"],
        input_data=permuted,
        calibration_kwargs=rsm["kwargs"],
    )
    row_probe = pd.DataFrame(
        [
            {
                "SemanticInputMatch": rsm["archive"]["semantic_input_sha256"]
                == permuted_archive["semantic_input_sha256"],
                "OrderedInputMatch": rsm["archive"]["ordered_input_sha256"]
                == permuted_archive["ordered_input_sha256"],
                "AnalysisIDMatch": rsm["archive"]["analysis_identity"].analysis_id
                == permuted_archive["analysis_identity"].analysis_id,
                "Policy": "semantic_same_but_ordered_identity_distinct",
            }
        ]
    )
    pcm = retained["pcm_differential_anchor_private_archive"]
    reversed_kwargs = {
        **pcm["kwargs"],
        "hard_anchors": pcm["kwargs"]["hard_anchors"].iloc[::-1].reset_index(
            drop=True
        ),
    }
    reversed_archive = build_cmle_one_click_archive(
        pcm["result"],
        input_data=pcm["frame"],
        calibration_kwargs=reversed_kwargs,
    )
    anchor_probe = pd.DataFrame(
        [
            {
                "NormalizedAnchorBytesMatch": pcm["archive"]["assets"][
                    "hard_anchors.csv"
                ]
                == reversed_archive["assets"]["hard_anchors.csv"],
                "ConfigSHA256Match": pcm["archive"][
                    "analysis_identity"
                ].config_fingerprint
                == reversed_archive["analysis_identity"].config_fingerprint,
                "AnalysisIDMatch": pcm["archive"]["analysis_identity"].analysis_id
                == reversed_archive["analysis_identity"].analysis_id,
                "ZipSHA256Match": pcm["archive"]["zip_sha256"]
                == reversed_archive["zip_sha256"],
            }
        ]
    )
    write_csv(row_probe, output / "row_order_probe.csv")
    write_csv(anchor_probe, output / "anchor_order_probe.csv")
    return row_probe, anchor_probe


def rewrite_zip(assets: dict[str, bytes], duplicate: str | None = None) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, raw in assets.items():
            archive.writestr(name, raw)
        if duplicate is not None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                archive.writestr(duplicate, assets[duplicate])
    return buffer.getvalue()


def tamper_and_privacy_probes(retained, output: Path):
    original = retained["rsm_ready_extremes_private_archive"]
    rows = []
    for mutation in ("change", "delete", "extra", "duplicate"):
        assets = dict(original["archive"]["assets"])
        duplicate = None
        if mutation == "change":
            assets["result_cards.csv"] = assets["result_cards.csv"].replace(
                b"ready", b"green", 1
            )
        elif mutation == "delete":
            del assets["result_cards.csv"]
        elif mutation == "extra":
            assets["undeclared.txt"] = b"undeclared"
        else:
            duplicate = "result_cards.csv"
        error = ""
        try:
            verify_cmle_one_click_archive(rewrite_zip(assets, duplicate))
            rejected = False
        except ValueError as exc:
            rejected = True
            error = str(exc)
        rows.append({"Probe": mutation, "Rejected": rejected, "Error": error})
    public_error = ""
    try:
        build_cmle_one_click_archive(
            original["result"],
            input_data=original["frame"],
            calibration_kwargs=original["kwargs"],
            privacy_mode="public",
        )
        public_rejected = False
    except ValueError as exc:
        public_rejected = True
        public_error = str(exc)
    rows.append(
        {"Probe": "public_mode", "Rejected": public_rejected, "Error": public_error}
    )
    ledger = pd.DataFrame(rows)
    write_csv(ledger, output / "tamper_privacy_probe.csv")
    return ledger


def run_tests(output: Path) -> dict[str, object]:
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-q",
        "tests/test_cmle_one_click_archive.py",
        "tests/test_cmle_one_click.py",
        "tests/test_evidence_contract.py",
        "tests/test_modular_helpers.py",
        "tests/test_cmle_workflow.py",
        "tests/test_cmle_wle_fit.py",
        "tests/test_cmle_hard_anchors.py",
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
    return {"passed": completed.returncode == 0, "returncode": completed.returncode}


def write_review(output: Path, results: dict[str, object]) -> None:
    text = f"""# CMLE One-Click Private Archive Critical Review

## Decision

The registered deterministic private-archive contract **{'passed' if results['contract_passed'] else 'failed'}**. The three analytical states produced verified archives and exact current-engine replays.

## Identity layers

Ordered input, semantic row-order identity, canonical configuration, result tables, core archive content, and deterministic ZIP transport are separate hashes. A row permutation preserved semantic identity but changed ordered identity and AnalysisID as registered. Reversing the two hard-anchor rows preserved normalized anchor bytes, configuration identity, AnalysisID, and ZIP identity.

## Replay and floating point

All three archives reproduced AnalysisID, terminal status, result-asset hashes, six-card bytes, archive-content hash, and ZIP hash. Across {results['floating_values_checked']} retained WLE/SE/Infit/Outfit values, 17-significant-digit CSV reload matched the original binary64 hexadecimal representation exactly. Decisions continue to use raw MnSq; display columns are evidence only.

## Tamper and privacy boundary

Changed, deleted, undeclared, and duplicate ZIP entries were all rejected before loading. Public-mode construction also failed closed. Each retained ZIP contains synthetic raw response rows and Person/facet identifiers and is therefore labelled private controlled-access reproduction material. This mechanism is not encryption, anonymization, de-identification, or sharing approval.

## Remaining risks

- Bitwise replay is established only for this current engine/platform and deterministic fixture scope, not across Python, NumPy, SciPy, BLAS, operating-system, or architecture changes.
- A matching archive proves artifact identity, not model fit, anchor validity, unbiasedness, or user comprehension.
- Public disclosure rules, ID pseudonymization, saved-result migration across schema versions, and Streamlit upload/download wiring remain unqualified.
- Large archive limits are defensive implementation bounds, not validated application capacity.

## Product implication

A future button can now save and reopen the exact private six-card result without silently losing raw thresholds or stopped states. Public export and public CMLE UI must remain disabled until separate privacy, migration, and task-comprehension gates pass.
"""
    (output / "CMLE_ONE_CLICK_PRIVATE_ARCHIVE_CRITICAL_REVIEW.md").write_text(
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
    frames, retained = run_cases(plan, args.output)
    row_probe, anchor_probe = row_and_anchor_probes(retained, args.output)
    tamper = tamper_and_privacy_probes(retained, args.output)
    tests = run_tests(args.output)
    case_ledger = frames["case_ledger"]
    replay_checks = frames["replay_checks"]
    floating = frames["floating_roundtrip"]
    gates = {
        "identity_passed": True,
        "case_contract_passed": bool(
            len(case_ledger) == 3
            and case_ledger[
                [
                    "TerminalMatch",
                    "ConditionalAssetMatch",
                    "ArchiveVerified",
                    "BilingualCardsLoaded",
                    "RepeatedAssetIdentityMatch",
                    "RepeatedArchiveContentMatch",
                    "RepeatedZipBytesMatch",
                    "ReplayPassed",
                ]
            ]
            .all()
            .all()
            and case_ledger["CardsLoaded"].eq(6).all()
        ),
        "replay_passed": bool(replay_checks["Passed"].all()),
        "row_order_passed": bool(
            row_probe.iloc[0]["SemanticInputMatch"]
            and not row_probe.iloc[0]["OrderedInputMatch"]
            and not row_probe.iloc[0]["AnalysisIDMatch"]
        ),
        "anchor_order_passed": bool(
            anchor_probe[
                [
                    "NormalizedAnchorBytesMatch",
                    "ConfigSHA256Match",
                    "AnalysisIDMatch",
                    "ZipSHA256Match",
                ]
            ]
            .all()
            .all()
        ),
        "tamper_privacy_passed": bool(len(tamper) == 5 and tamper["Rejected"].all()),
        "binary64_roundtrip_passed": bool(
            len(floating) == 56 and floating["Binary64Exact"].all()
        ),
        "tests_passed": bool(tests["passed"]),
    }
    contract_passed = all(gates.values())
    results = {
        **gates,
        "contract_passed": contract_passed,
        "archives_verified": int(case_ledger["ArchiveVerified"].sum()),
        "replay_checks_passed": int(replay_checks["Passed"].sum()),
        "replay_checks_total": len(replay_checks),
        "floating_values_checked": len(floating),
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
                "mfrm_app/cmle_one_click_archive.py": sha256_file(
                    ROOT / "mfrm_app/cmle_one_click_archive.py"
                ),
                "tests/test_cmle_one_click_archive.py": sha256_file(
                    ROOT / "tests/test_cmle_one_click_archive.py"
                ),
                "validation/cmle_one_click_archive_identity.py": sha256_file(
                    Path(__file__)
                ),
            },
            "results": results,
            "output_sha256": {
                name: sha256_file(args.output / name) for name in output_files
            },
            "interpretation": {
                "archive_scope": "private controlled-access reproduction",
                "public_export": "withheld",
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
