"""Evidence-contract adapter for the legacy final-readiness table.

The Streamlit UI still exposes its established five human-facing columns.
This module appends a validated :class:`~mfrm_app.evidence.EvidenceRecord`
ledger without guessing what the ambiguous legacy ``Review`` label means.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

import pandas as pd

from . import evidence


LEGACY_READINESS_COLUMNS = (
    "Check",
    "Status",
    "Evidence",
    "ActionBeforeFinalReport",
    "Required",
)


@dataclass(frozen=True)
class ReadinessEvidenceSpec:
    """Explicit machine semantics for one legacy final-readiness row."""

    evidence_key: str
    domain: str
    question: str
    computation_state: evidence.ComputationState | str
    interpretation_boundary: str
    reason_code: evidence.ReasonCode | str | None = None
    stability_state: evidence.StabilityState | str = evidence.StabilityState.NOT_ASSESSED
    reason_args: Mapping[str, object] = field(default_factory=dict)
    scope: Mapping[str, object] = field(default_factory=dict)
    prerequisites: tuple[str, ...] = ()
    observed: Mapping[str, object] = field(default_factory=dict)
    uncertainty: Mapping[str, object] = field(default_factory=dict)
    next_inspection: str = ""
    source_artifacts: tuple[str, ...] = ()
    sensitivity_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class FinalReadinessEvidenceBundle:
    """Legacy-compatible frame plus its identity and validated records."""

    identity: evidence.AnalysisIdentity
    records: tuple[evidence.EvidenceRecord, ...]
    frame: pd.DataFrame


def attach_final_readiness_evidence(
    legacy_frame: pd.DataFrame,
    *,
    identity: evidence.AnalysisIdentity,
    row_specs: Mapping[str, ReadinessEvidenceSpec],
) -> FinalReadinessEvidenceBundle:
    """Append canonical evidence columns while preserving legacy columns.

    Every row must have an explicit specification.  In particular, callers
    must decide whether ``Review`` means computed-with-caution, held, or not
    assessable; the adapter never infers that distinction from prose.
    """
    if not isinstance(legacy_frame, pd.DataFrame):
        raise evidence.ContractValidationError("legacy_frame must be a DataFrame")
    missing_columns = [
        column for column in LEGACY_READINESS_COLUMNS if column not in legacy_frame.columns
    ]
    if missing_columns:
        raise evidence.ContractValidationError(
            f"Final-readiness frame is missing legacy columns: {missing_columns!r}"
        )
    legacy = legacy_frame.loc[:, LEGACY_READINESS_COLUMNS].copy().reset_index(drop=True)
    if legacy["Check"].astype(str).duplicated().any():
        duplicates = sorted(
            legacy.loc[legacy["Check"].astype(str).duplicated(False), "Check"]
            .astype(str)
            .unique()
        )
        raise evidence.ContractValidationError(
            f"Final-readiness Check values must be unique: {duplicates!r}"
        )
    checks = tuple(legacy["Check"].astype(str))
    missing_specs = sorted(set(checks).difference(row_specs))
    unknown_specs = sorted(set(row_specs).difference(checks))
    if missing_specs or unknown_specs:
        raise evidence.ContractValidationError(
            "Final-readiness evidence specification mismatch; "
            f"missing={missing_specs!r}, unknown={unknown_specs!r}"
        )

    records: list[evidence.EvidenceRecord] = []
    for _, row in legacy.iterrows():
        check = str(row["Check"])
        status = str(row["Status"])
        required_label = str(row["Required"])
        if required_label not in {"Yes", "No"}:
            raise evidence.ContractValidationError(
                f"Final-readiness Required must be 'Yes' or 'No', got {required_label!r}"
            )
        spec = row_specs[check]
        computation_state = evidence.normalize_legacy_status(
            status,
            vocabulary=evidence.LegacyVocabulary.FINAL_READINESS,
            explicit_state=spec.computation_state,
        )
        scope = {
            **dict(spec.scope),
            "legacy_status": status,
            "readiness_check": check,
        }
        observed = {
            **dict(spec.observed),
            "legacy_evidence": str(row["Evidence"]),
        }
        records.append(
            evidence.make_evidence_record(
                identity,
                evidence_key=spec.evidence_key,
                domain=spec.domain,
                question=spec.question,
                computation_state=computation_state,
                stability_state=spec.stability_state,
                required=required_label == "Yes",
                reason_code=spec.reason_code,
                reason_args=spec.reason_args,
                summary=str(row["Evidence"]),
                scope=scope,
                prerequisites=spec.prerequisites,
                observed=observed,
                uncertainty=spec.uncertainty,
                interpretation_boundary=spec.interpretation_boundary,
                recommended_action=str(row["ActionBeforeFinalReport"]),
                next_inspection=spec.next_inspection or str(row["ActionBeforeFinalReport"]),
                source_artifacts=spec.source_artifacts,
                sensitivity_ids=spec.sensitivity_ids,
            )
        )

    normalized_records = evidence.validate_evidence_records(records, identity=identity)
    contract_frame = evidence.evidence_records_to_frame(
        normalized_records,
        identity=identity,
    ).rename(columns={"Required": "ContractRequired"})
    combined = pd.concat([legacy, contract_frame], axis=1)
    if not combined.loc[:, LEGACY_READINESS_COLUMNS].equals(legacy):
        raise evidence.ContractValidationError(
            "Attaching the evidence ledger changed legacy final-readiness columns"
        )
    combined.attrs["analysis_identity"] = identity.to_payload()
    combined.attrs["evidence_records"] = [record.to_payload() for record in normalized_records]
    return FinalReadinessEvidenceBundle(
        identity=identity,
        records=normalized_records,
        frame=combined,
    )
