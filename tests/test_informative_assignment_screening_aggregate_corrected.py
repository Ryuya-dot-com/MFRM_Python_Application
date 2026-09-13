from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from validation.informative_assignment_screening_aggregate_corrected import (
    enrich_all_ledgers,
    validate_original_aggregate,
)


STUDY = Path("validation/informative_assignment_screening20_20260811")


def _manifest() -> pd.DataFrame:
    return pd.DataFrame({
        "RunId": ["run-1"],
        "ConditionId": ["aligned__normal"],
        "Design": ["ability_severity_aligned_connected"],
        "PersonDistribution": ["normal"],
        "Replicate": [221],
    })


def test_metadata_enrichment_fills_missing_threshold_fields_only_by_run_id() -> None:
    base = pd.DataFrame({"RunId": ["run-1"], "Value": [1.0]})
    thresholds = base.assign(
        ConditionId=pd.NA,
        Design=pd.NA,
        PersonDistribution=pd.NA,
        Replicate=pd.NA,
    )
    runs, recovery, enriched, constraints, audit = enrich_all_ledgers(
        base, base, thresholds, base, _manifest()
    )
    assert len(runs) == len(recovery) == len(enriched) == len(constraints) == 1
    assert enriched.iloc[0]["Design"] == "ability_severity_aligned_connected"
    assert int(enriched.iloc[0]["Replicate"]) == 221
    assert audit["missing_before"]["thresholds"]["Design"] == 1
    assert audit["missing_after"]["thresholds"]["Design"] == 0
    assert audit["all_registered_metadata_complete"]


def test_metadata_enrichment_rejects_non_null_conflict() -> None:
    conflicting = pd.DataFrame({"RunId": ["run-1"], "Design": ["complete"]})
    with pytest.raises(ValueError, match="conflicts"):
        enrich_all_ledgers(conflicting, conflicting, conflicting, conflicting, _manifest())


@pytest.mark.retained_evidence
def test_registered_original_aggregate_is_hash_valid() -> None:
    identity = validate_original_aggregate(STUDY)
    assert "screening_metrics.json" in identity["artifact_sha256"]

