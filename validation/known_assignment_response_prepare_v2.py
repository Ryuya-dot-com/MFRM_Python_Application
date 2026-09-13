"""Create an immutable Seed-column remediation of the failed response input bundle."""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pandas as pd

from validation.known_assignment_large_design_dp_shard import _sha256
from validation.known_assignment_response_prepare import INPUT_FILES, STUDY_DIR as PARENT_STUDY


AMENDMENT_PATH = ROOT / "validation" / "known_assignment_response_seed_remediation_20260811.json"
STUDY_DIR = ROOT / "validation" / "known_assignment_response_screening_v2_20260811"


def prepare() -> dict[str, object]:
    if STUDY_DIR.exists():
        raise FileExistsError(f"Refusing to overwrite remediated study: {STUDY_DIR}")
    parent_identity_path = PARENT_STUDY / "input_identity.json"
    parent_identity = json.loads(parent_identity_path.read_text(encoding="utf-8"))
    for name, expected in parent_identity["retained_input_sha256"].items():
        if _sha256(PARENT_STUDY / "retained_input" / name) != expected:
            raise ValueError(f"Parent retained input changed: {name}")
    STUDY_DIR.mkdir(parents=False, exist_ok=False)
    input_dir = STUDY_DIR / "retained_input"
    input_dir.mkdir()
    for name in INPUT_FILES:
        if name in {"manifest.csv", "attempt_manifest.csv"}:
            continue
        shutil.copy2(PARENT_STUDY / "retained_input" / name, input_dir / name)
    manifest = pd.read_csv(PARENT_STUDY / "retained_input" / "manifest.csv")
    if "Seed" in manifest.columns:
        raise ValueError("Parent manifest unexpectedly already contains Seed")
    location = manifest.columns.get_loc("UniformSeed")
    manifest.insert(location, "Seed", manifest["UniformSeed"].astype("int64"))
    manifest.to_csv(input_dir / "manifest.csv", index=False, lineterminator="\n")
    attempts = pd.read_csv(PARENT_STUDY / "retained_input" / "attempt_manifest.csv")
    amendment_hash = _sha256(AMENDMENT_PATH)
    manifest_hash = _sha256(input_dir / "manifest.csv")
    attempts["ParentAttemptFingerprint"] = attempts["AttemptFingerprint"]
    attempts["AttemptFingerprint"] = [
        hashlib.sha256(
            f"{fingerprint}|{amendment_hash}|{manifest_hash}".encode("utf-8")
        ).hexdigest()
        for fingerprint in attempts["ParentAttemptFingerprint"]
    ]
    attempts.to_csv(input_dir / "attempt_manifest.csv", index=False, lineterminator="\n")
    hashes = {name: _sha256(input_dir / name) for name in INPUT_FILES}
    unchanged = {
        name: hashes[name] == parent_identity["retained_input_sha256"][name]
        for name in INPUT_FILES
        if name not in {"manifest.csv", "attempt_manifest.csv"}
    }
    checks = {
        "rows_unchanged": bool(all(unchanged.values())),
        "manifest_rows_30": len(manifest) == 30,
        "seed_present": "Seed" in manifest.columns,
        "seed_equals_uniform_seed": manifest["Seed"].equals(manifest["UniformSeed"]),
        "attempts_120": len(attempts) == 120,
        "attempt_fingerprints_rebound": bool(
            attempts["AttemptFingerprint"].ne(attempts["ParentAttemptFingerprint"]).all()
        ),
    }
    identity = {
        "schema_version": "known_assignment_response_screening_input_identity_v2",
        "parent_input_identity_sha256": _sha256(parent_identity_path),
        "remediation_sha256": amendment_hash,
        "plan_sha256": parent_identity["plan_sha256"],
        "retained_input_sha256": hashes,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    (STUDY_DIR / "input_identity.json").write_text(
        json.dumps(identity, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    if not identity["all_checks_pass"]:
        raise RuntimeError(f"Seed remediation audit failed: {checks}")
    print(json.dumps(identity, ensure_ascii=False, indent=2))
    return identity


if __name__ == "__main__":
    prepare()
