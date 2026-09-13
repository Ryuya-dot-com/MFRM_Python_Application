from __future__ import annotations

import json

import pandas as pd
import pytest

pytestmark = pytest.mark.retained_evidence

from validation import known_assignment_response_run as engine
from validation import known_assignment_response_run_v2 as wrapper
from validation.known_assignment_response_prepare import INPUT_FILES, STUDY_DIR as V1_STUDY


def test_seed_remediation_changes_only_manifest_and_attempt_identity_files():
    wrapper.configure_engine()
    identity = engine.validate_input_identity()
    parent = json.loads((V1_STUDY / "input_identity.json").read_text(encoding="utf-8"))
    manifest = pd.read_csv(wrapper.STUDY_DIR / "retained_input" / "manifest.csv")

    assert identity["all_checks_pass"] is True
    assert manifest["Seed"].equals(manifest["UniformSeed"])
    for name in INPUT_FILES:
        if name not in {"manifest.csv", "attempt_manifest.csv"}:
            assert identity["retained_input_sha256"][name] == parent["retained_input_sha256"][name]


def test_v2_attempt_fingerprints_are_rebound_without_denominator_change():
    attempts = pd.read_csv(wrapper.STUDY_DIR / "retained_input" / "attempt_manifest.csv")

    assert len(attempts) == 120
    assert attempts["AttemptFingerprint"].nunique() == 120
    assert attempts["ParentAttemptFingerprint"].nunique() == 120
    assert attempts["AttemptFingerprint"].ne(attempts["ParentAttemptFingerprint"]).all()
