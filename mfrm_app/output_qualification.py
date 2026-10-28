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

import numpy as np


OUTPUT_QUALIFICATION_SCHEMA_VERSION = "mfrm_output_qualification_v1"
FREE_SD_MML_INFERENCE_HOLD_REASON = "stat.mml_free_sd.numerical_qualification_pending"
FIXED_SD_MML_INFERENCE_HOLD_REASON = "stat.mml_fixed_sd.numerical_qualification_pending"
OPTIMIZER_TERMINATION_HOLD_REASON = "stat.estimation.optimizer_termination_unverified"


def mml_inference_hold_reason(config) -> str:
    if not isinstance(config, dict) or str(config.get("method", "")).upper() != "MML":
        return ""
    return (FREE_SD_MML_INFERENCE_HOLD_REASON if config.get("estimate_population_sd", False)
            else FIXED_SD_MML_INFERENCE_HOLD_REASON)


def mml_inference_hold_message(config) -> str:
    if not mml_inference_hold_reason(config):
        return ""
    scope = "Free-population-SD" if config.get("estimate_population_sd", False) else "Fixed-population-SD"
    return f"{scope} MML numerical accuracy and sampling uncertainty remain unqualified."


def jmle_inference_hold_message(result) -> str:
    config = result.get("config", {})
    if (isinstance(config, dict) and str(config.get("method", "")).upper() == "JMLE"
            and config.get("jmle_facet_boundary_review", False)):
        from mfrm_app.jmle_boundary import HOLD_MESSAGE
        return HOLD_MESSAGE
    audit = result.get("identifiability")
    status = audit.get("status") if isinstance(audit, dict) else None
    if status in {"eta_rank_deficient", "audit_size_limit"}:
        return "JMLE person/facet identification has not been established; review the design before inference."
    return ""


def optimizer_termination_established(result) -> bool:
    """Require an explicitly recorded Boolean; truthy legacy values are not evidence."""
    stopped = getattr((result or {}).get("opt"), "success", None)
    return isinstance(stopped, (bool, np.bool_)) and bool(stopped)


def calibration_inference_hold_reason(result) -> str:
    """Requalify current and saved results, preserving method-specific reasons."""
    result = result or {}
    config = result.get("config", {})
    if reason := mml_inference_hold_reason(config):
        return reason
    if jmle_inference_hold_message(result):
        if (isinstance(config, dict) and str(config.get("method", "")).upper() == "JMLE"
                and config.get("jmle_facet_boundary_review", False)):
            from mfrm_app.jmle_boundary import HOLD_REASON
            return HOLD_REASON
        return "stat.jmle.identification_unqualified"
    return "" if optimizer_termination_established(result) else OPTIMIZER_TERMINATION_HOLD_REASON


def calibration_inference_hold_message(result) -> str:
    result = result or {}
    hold = mml_inference_hold_message(result.get("config")) or jmle_inference_hold_message(result)
    if hold:
        return hold
    if not optimizer_termination_established(result):
        return "Successful optimizer termination has not been established; review the numerical diagnostics."
    return ""


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
    POPULATION_SD_UNCERTAINTY_UNQUALIFIED = (
        "stat.stat_003.population_sd_uncertainty_unqualified"
    )
    BIAS_CELL_UNCALIBRATED = "bias.bias_003.cell_reference_uncalibrated"
    BIAS_PAIRWISE_UNVERIFIED = "bias.bias_001_002.pairwise_unverified"
    FACET_EQUIVALENCE_UNVALIDATED = "stat.stat_003.facet_equivalence_unvalidated"


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
            "Exact constraint coordinates are counted; observed-design identification and estimator/comparison scope still require validation before ranking models.",
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
    # A missing/skipped Hessian does not establish a deficient rank.
    if (
        status_label == "regularized"
        or (status_label == "ok" and param_value > 0 and rank_value < param_value)
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


def population_sd_uncertainty_qualification() -> OutputQualification:
    """Withhold free-SD inference, including legacy fixed-nuisance SE/CI."""
    return _qualification(
        "mml.population_sd_uncertainty",
        OutputUse.WITHHELD,
        QualificationReason.POPULATION_SD_UNCERTAINTY_UNQUALIFIED,
        "STAT-003",
        "G2",
        "Do not report population-SD SE/CI until joint stationarity, "
        "nuisance-adjusted uncertainty, and coverage are qualified. "
        "Fixed-nuisance curvature is a technical diagnostic, not profile information.",
    )


def bias_cell_qualification() -> OutputQualification:
    """Retain scalar likelihood diagnostics without promoting uncalibrated claims."""
    return _qualification(
        "bias.cell_inference", OutputUse.TECHNICAL_ONLY,
        QualificationReason.BIAS_CELL_UNCALIBRATED, "BIAS-003", "G3",
        "Conditional LR and cutoff intervals hold calibration fixed. "
        "Wait for error-rate, power, coverage, and calibration-uncertainty validation "
        "before reporting bias or no-bias conclusions.",
    )


def bias_pairwise_qualification() -> OutputQualification:
    """Keep conditional local contrasts technical until G3 inference is qualified."""
    return _qualification(
        "bias.pairwise_local_measure",
        OutputUse.WITHHELD,
        QualificationReason.BIAS_PAIRWISE_UNVERIFIED,
        "BIAS-001;BIAS-002",
        "G3",
        "Use conditional pairwise contrasts for technical audit only. "
        "Boundary behavior, joint calibration uncertainty, and reference-distribution "
        "and error-rate validation remain unqualified.",
    )


def facet_equivalence_qualification() -> OutputQualification:
    """Withhold equivalence claims until their uncertainty and decision rules are validated."""
    return _qualification(
        "facet.equivalence", OutputUse.WITHHELD,
        QualificationReason.FACET_EQUIVALENCE_UNVALIDATED, "STAT-003", "G2",
        "Specify the contrast, substantive equivalence margin and claim family. "
        "Validate joint uncertainty and the decision rule before reporting equivalence; "
        "nonsignificance and agreement among unvalidated methods do not establish it.",
        raw_export=False,
    )


def records(
    qualifications: Iterable[OutputQualification],
) -> list[dict[str, object]]:
    """Convert qualifications to a stable table-ready record sequence."""
    return [qualification.to_record() for qualification in qualifications]
