"""Seed-remediated wrapper around the immutable registered response runner."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validation import known_assignment_response_run as engine
from validation.known_assignment_large_design_dp_shard import _sha256


STUDY_DIR = ROOT / "validation" / "known_assignment_response_screening_v2_20260811"
REGISTRATION_PATH = ROOT / "validation" / "known_assignment_response_execution_registration_v3_20260811.json"
AMENDMENT_PATH = ROOT / "validation" / "known_assignment_response_seed_remediation_20260811.json"
IMPORT_AMENDMENT_PATH = ROOT / "validation" / "known_assignment_response_wrapper_import_amendment_20260811.json"
PREPARE_V2_PATH = ROOT / "validation" / "known_assignment_response_prepare_v2.py"
_BASE_BUILD_DEPENDENCIES = engine.build_dependency_manifest


def _build_dependencies_v2() -> dict[str, str]:
    manifest = _BASE_BUILD_DEPENDENCIES()
    for name, path in {
        "validation/known_assignment_response_run_v2.py": Path(__file__).resolve(),
        "validation/known_assignment_response_seed_remediation_20260811.json": AMENDMENT_PATH,
        "validation/known_assignment_response_wrapper_import_amendment_20260811.json": IMPORT_AMENDMENT_PATH,
        "validation/known_assignment_response_prepare_v2.py": PREPARE_V2_PATH,
    }.items():
        manifest[name] = _sha256(path)
    return dict(sorted(manifest.items()))


def configure_engine() -> None:
    engine.STUDY_DIR = STUDY_DIR
    engine.REGISTRATION_PATH = REGISTRATION_PATH
    engine.INPUT_IDENTITY_PATH = STUDY_DIR / "input_identity.json"
    engine.EXECUTION_IDENTITY_PATH = STUDY_DIR / "execution_identity.json"
    engine.build_dependency_manifest = _build_dependencies_v2


def validate_v2_registration() -> dict[str, object]:
    configure_engine()
    value = json.loads(REGISTRATION_PATH.read_text(encoding="utf-8"))
    expected = {
        "wrapper_sha256": _sha256(Path(__file__).resolve()),
        "seed_remediation_sha256": _sha256(AMENDMENT_PATH),
        "prepare_v2_sha256": _sha256(PREPARE_V2_PATH),
    }
    for key, digest in expected.items():
        if str(value.get(key, "")).lower() != digest.lower():
            raise ValueError(f"V2 execution registration mismatch: {key}")
    engine.validate_registration()
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--facets-exe", type=Path, default=engine.FACETS_EXE_DEFAULT)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--timeout-seconds", type=float, default=120.0)
    args = parser.parse_args()
    validate_v2_registration()
    engine.run_shard(
        facets_exe=args.facets_exe,
        shard_index=args.shard_index,
        shard_count=args.shard_count,
        resume=args.resume,
        timeout_seconds=args.timeout_seconds,
    )


if __name__ == "__main__":
    main()
