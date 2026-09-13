from __future__ import annotations

import pandas as pd
import pytest

pytestmark = pytest.mark.retained_evidence

from validation import known_assignment_response_prepare as prepare
from validation import known_assignment_response_run as runner


def test_prepared_known_assignment_response_identity_and_denominators_are_frozen():
    identity = runner.validate_input_identity()
    manifest = pd.read_csv(prepare.STUDY_DIR / "retained_input" / "manifest.csv")
    attempts = pd.read_csv(prepare.STUDY_DIR / "retained_input" / "attempt_manifest.csv")

    assert identity["all_checks_pass"] is True
    assert len(manifest) == 30
    assert len(attempts) == 120
    assert manifest.groupby("Gamma").size().to_dict() == {-0.8: 10, 0.0: 10, 0.8: 10}
    assert attempts.groupby("AttemptType").size().eq(30).all()
    assert set(attempts["AttemptType"]) == set(prepare.ATTEMPT_TYPES)


def test_assignment_correlation_is_ordered_and_complete_response_is_shared():
    manifest = pd.read_csv(prepare.STUDY_DIR / "retained_input" / "manifest.csv")
    means = manifest.groupby("Gamma")["AssignmentCorrelation"].mean().sort_index()

    assert means.loc[-0.8] < means.loc[0.0] < means.loc[0.8]
    assert manifest.groupby("Replicate")["CompleteResponseSHA256"].nunique().eq(1).all()
    assert manifest["ExpectedStructuralNullity"].eq(0).all()
    assert manifest["CompleteCriterionCategorySupport"].all()


def test_attempt_shards_are_disjoint_complete_and_deterministic():
    attempts = pd.read_csv(prepare.STUDY_DIR / "retained_input" / "attempt_manifest.csv")
    shards = [
        runner.select_shard_attempts(attempts, shard_index=index, shard_count=7)
        for index in range(7)
    ]
    ordinals = pd.concat(shards)["AttemptOrdinal"].astype(int)

    assert len(ordinals) == 120
    assert ordinals.nunique() == 120
    assert set(ordinals) == set(range(120))
    assert runner.select_shard_attempts(
        attempts, shard_index=3, shard_count=7
    ).equals(shards[3])
