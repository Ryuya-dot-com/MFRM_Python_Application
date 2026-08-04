"""Fail-closed qualification contracts for statistically sensitive outputs.

The numerical engine can retain intermediate values for debugging and migration
while the public application withholds conclusions that have not cleared their
validation gate.  This module keeps that distinction machine-readable and free
of Streamlit dependencies.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from typing import Iterable


OUTPUT_QUALIFICATION_SCHEMA_VERSION = "mfrm_output_qualification_v1"


class OutputUse(str, Enum):
    """Allowed use of an output under the active remediation baseline."""

    READY = "READY"
    CONDITIONAL = "CONDITIONAL"
    SCREENING_ONLY = "SCREENING_ONLY"
    TECHNICAL_ONLY = "TECHNICAL_ONLY"
    WITHHELD = "WITHHELD"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class QualificationReason(str, Enum):
    """Stable reason codes linked to the statistical remediation ledger."""

    PARAMETER_COUNT_UNIDENTIFIED = "stat.stat_002.parameter_count_unidentified"
    MODEL_SELECTION_SCOPE_UNVALIDATED = (
        "model.model_001.selection_scope_unvalidated"
    )
    JMLE_AUTO_RECOMMENDATION_WITHHELD = (
        "model.model_001.jmle_auto_recommendation_withheld"
    )
    RSM_PCM_LRT_REGULARITY_UNVALIDATED = (
        "model.model_001.rsm_pcm_lrt_regularity_unvalidated"
    )
    GPCM_LRT_UNCALIBRATED = "model.model_002.gpcm_lrt_uncalibrated"
    COVARIANCE_RANK_DEFICIENT = "stat.stat_003.covariance_rank_deficient"
    COVARIANCE_COVERAGE_UNVALIDATED = (
        "stat.stat_003.covariance_coverage_unvalidated"
    )
    COVARIANCE_UNAVAILABLE = "stat.stat_003.covariance_unavailable"
    COVARIANCE_NOT_APPLICABLE = "stat.stat_003.covariance_not_applicable"
    BIAS_PAIRWISE_UNVERIFIED = "bias.bias_001_002.pairwise_unverified"


@dataclass(frozen=True)
class OutputQualification:
    """One exportable qualification decision for a named output."""

    output_id: str
    status: OutputUse
    reason_code: QualificationReason
    roadmap_issue_id: str
    next_gate: str
    raw_technical_export_allowed: bool
    public_conclusion_allowed: bool
    user_action: str
    schema_version: str = OUTPUT_QUALIFICATION_SCHEMA_VERSION

    def to_record(self) -> dict[str, object]:
        record = asdict(self)
        record["status"] = self.status.value
        record["reason_code"] = self.reason_code.value
        return {
            "SchemaVersion": record["schema_version"],
            "OutputID": record["output_id"],
            "QualificationStatus": record["status"],
            "ReasonCode": record["reason_code"],
            "RoadmapIssueID": record["roadmap_issue_id"],
            "NextGate": record["next_gate"],
            "RawTechnicalExportAllowed": record["raw_technical_export_allowed"],
            "PublicConclusionAllowed": record["public_conclusion_allowed"],
            "UserAction": record["user_action"],
        }


def _qualification(
    output_id: str,
    status: OutputUse,
    reason: QualificationReason,
    issue: str,
    gate: str,
    action: str,
    *,
    raw_export: bool = True,
    public_conclusion: bool = False,
) -> OutputQualification:
    return OutputQualification(
        output_id=output_id,
        status=status,
        reason_code=reason,
        roadmap_issue_id=issue,
        next_gate=gate,
        raw_technical_export_allowed=raw_export,
        public_conclusion_allowed=public_conclusion,
        user_action=action,
    )


def model_choice_qualifications(method: object) -> tuple[OutputQualification, ...]:
    """Qualification matrix for IC, recommendations, and nested LR outputs."""
    method_label = str(method or "").upper()
    recommendation_reason = (
        QualificationReason.JMLE_AUTO_RECOMMENDATION_WITHHELD
        if method_label == "JMLE"
        else QualificationReason.MODEL_SELECTION_SCOPE_UNVALIDATED
    )
    recommendation_issue = "MODEL-001"
    return (
        _qualification(
            "model_choice.information_criteria",
            OutputUse.TECHNICAL_ONLY,
            QualificationReason.MODEL_SELECTION_SCOPE_UNVALIDATED,
            "MODEL-001",
            "G2",
            "Parameter counts are identified; wait for estimator/comparison-scope validation before ranking models.",
        ),
        _qualification(
            "model_choice.automatic_recommendation",
            OutputUse.WITHHELD,
            recommendation_reason,
            recommendation_issue,
            "G2",
            "Do not generate or report an automatic preferred-model conclusion.",
        ),
        _qualification(
            "model_choice.lrt.rsm_pcm",
            OutputUse.WITHHELD,
            QualificationReason.RSM_PCM_LRT_REGULARITY_UNVALIDATED,
            "MODEL-001",
            "G2",
            "Do not report a chi-square p-value or retain/reject decision.",
        ),
        _qualification(
            "model_choice.lrt.pcm_gpcm",
            OutputUse.WITHHELD,
            QualificationReason.GPCM_LRT_UNCALIBRATED,
            "MODEL-002",
            "G2",
            "Do not use an ordinary chi-square reference distribution.",
        ),
        _qualification(
            "model_choice.lrt.rsm_gpcm",
            OutputUse.WITHHELD,
            QualificationReason.GPCM_LRT_UNCALIBRATED,
            "MODEL-002",
            "G2",
            "Do not use an ordinary chi-square reference distribution.",
        ),
    )


def covariance_qualification(
    *,
    status: object,
    rank: object = None,
    param_count: object = None,
) -> OutputQualification:
    """Qualify MML structural covariance without changing retained numerics."""
    status_label = str(status or "not_available").lower()
    try:
        rank_value = int(rank)
        param_value = int(param_count)
    except (TypeError, ValueError, OverflowError):
        rank_value = 0
        param_value = 0

    if status_label == "not_applicable":
        return _qualification(
            "mml.structural_covariance",
            OutputUse.NOT_APPLICABLE,
            QualificationReason.COVARIANCE_NOT_APPLICABLE,
            "STAT-003",
            "G2",
            "Use the estimation-method-specific uncertainty status instead.",
            raw_export=False,
        )
    if (
        status_label == "regularized"
        or (param_value > 0 and rank_value < param_value)
    ):
        return _qualification(
            "mml.structural_covariance",
            OutputUse.WITHHELD,
            QualificationReason.COVARIANCE_RANK_DEFICIENT,
            "STAT-003",
            "G1",
            "Do not use the regularized covariance for confirmatory SE or CI claims.",
        )
    if status_label == "ok" and param_value > 0 and rank_value == param_value:
        return _qualification(
            "mml.structural_covariance",
            OutputUse.TECHNICAL_ONLY,
            QualificationReason.COVARIANCE_COVERAGE_UNVALIDATED,
            "STAT-003",
            "G2",
            "Retain for audit; wait for coverage validation before confirmatory use.",
        )
    return _qualification(
        "mml.structural_covariance",
        OutputUse.WITHHELD,
        QualificationReason.COVARIANCE_UNAVAILABLE,
        "STAT-003",
        "G2",
        "Do not claim structural covariance-based SE or CI.",
        raw_export=False,
    )


def bias_pairwise_qualification() -> OutputQualification:
    """Hold pairwise local measures until direction and contrast SE are fixed."""
    return _qualification(
        "bias.pairwise_local_measure",
        OutputUse.WITHHELD,
        QualificationReason.BIAS_PAIRWISE_UNVERIFIED,
        "BIAS-001;BIAS-002",
        "G3",
        "Use cell-level bias screening only; do not interpret pairwise local measures.",
    )


def records(
    qualifications: Iterable[OutputQualification],
) -> list[dict[str, object]]:
    """Convert qualifications to a stable table-ready record sequence."""
    return [qualification.to_record() for qualification in qualifications]
