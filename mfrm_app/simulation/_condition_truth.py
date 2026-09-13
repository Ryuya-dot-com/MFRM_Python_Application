"""RSM truth and rater-effect scientific-condition contracts."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math

from ._condition_shared import (
    MAX_ABS_ADJACENT_THRESHOLD,
    MAX_CENTRAL_TENDENCY_LOG_SCALE_MEAN,
    MAX_GENERATING_STANDARD_DEVIATION,
    MAX_RATING_MIN_ABS,
    MAX_THRESHOLD_SPAN,
    MODEL_RSM,
    SimulationConditionValidationError,
    _fingerprint,
    _json_safe,
    _require_int,
    _require_literal,
    _require_real,
    _strict_mapping,
    _tuple_payload_field,
)


RATER_EFFECT_CONDITION_VERSION = "mfrm_rater_effect_condition_v1"
RSM_TRUTH_SPEC_VERSION = "mfrm_rsm_truth_spec_v1"


@dataclass(frozen=True, slots=True, kw_only=True)
class RaterEffectConditionV1:
    """Rater-side data-generating effects for one RSM stress condition.

    The adjacent-category equation is documented in
    ``docs/sample_size_design_mvp.md``.  Positive severity and local bias
    lower the expected score.  Positive central-tendency log scale expands
    the centered threshold spacing and therefore favors interior categories.
    """

    rater_severity_sd: float = 0.35
    central_tendency_log_scale_mean: float = 0.0
    central_tendency_log_scale_sd: float = 0.0
    rater_person_local_bias_sd: float = 0.0
    effect_distribution: str = "normal"
    severity_sd_basis: str = "pre_centering_superpopulation_draw_sd"
    severity_constraint: str = "sample_mean_zero"
    local_bias_sd_basis: str = "superpopulation_draw_sd"
    local_bias_centering: str = "population_mean_zero"
    cross_effect_dependence: str = "independent"
    schema_version: str = RATER_EFFECT_CONDITION_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != RATER_EFFECT_CONDITION_VERSION:
            raise SimulationConditionValidationError(
                f"Unsupported rater-effect version: {self.schema_version!r}"
            )
        for name in (
            "rater_severity_sd",
            "central_tendency_log_scale_sd",
            "rater_person_local_bias_sd",
        ):
            _require_real(
                name,
                getattr(self, name),
                minimum=0.0,
                maximum=MAX_GENERATING_STANDARD_DEVIATION,
            )
        _require_real(
            "central_tendency_log_scale_mean",
            self.central_tendency_log_scale_mean,
            minimum=0.0,
            maximum=MAX_CENTRAL_TENDENCY_LOG_SCALE_MEAN,
        )
        _require_literal("effect_distribution", self.effect_distribution, "normal")
        _require_literal(
            "severity_sd_basis",
            self.severity_sd_basis,
            "pre_centering_superpopulation_draw_sd",
        )
        _require_literal(
            "severity_constraint", self.severity_constraint, "sample_mean_zero"
        )
        _require_literal(
            "local_bias_centering", self.local_bias_centering, "population_mean_zero"
        )
        _require_literal(
            "local_bias_sd_basis",
            self.local_bias_sd_basis,
            "superpopulation_draw_sd",
        )
        _require_literal(
            "cross_effect_dependence", self.cross_effect_dependence, "independent"
        )

    def to_dict(self) -> dict:
        return _json_safe({
            "rater_severity_sd": float(self.rater_severity_sd),
            "central_tendency_log_scale_mean": float(
                self.central_tendency_log_scale_mean
            ),
            "central_tendency_log_scale_sd": float(
                self.central_tendency_log_scale_sd
            ),
            "rater_person_local_bias_sd": float(self.rater_person_local_bias_sd),
            "effect_distribution": self.effect_distribution,
            "severity_sd_basis": self.severity_sd_basis,
            "severity_constraint": self.severity_constraint,
            "local_bias_sd_basis": self.local_bias_sd_basis,
            "local_bias_centering": self.local_bias_centering,
            "cross_effect_dependence": self.cross_effect_dependence,
            "schema_version": self.schema_version,
        })


def normalize_rater_effect_condition(
    value: RaterEffectConditionV1 | Mapping,
) -> RaterEffectConditionV1:
    if type(value) is RaterEffectConditionV1:
        return value
    payload = _strict_mapping(value, RaterEffectConditionV1, label="rater effect condition")
    return RaterEffectConditionV1(**payload)


def rater_effect_condition_fingerprint(
    value: RaterEffectConditionV1 | Mapping,
    *,
    length: int = 16,
) -> str:
    return _fingerprint(normalize_rater_effect_condition(value).to_dict(), length=length)


@dataclass(frozen=True, slots=True, kw_only=True)
class RsmTruthSpecV1:
    """Explicit v1 RSM truth, excluding design-owned level counts."""

    rating_min: int
    n_categories: int
    adjacent_thresholds: tuple[float, ...]
    person_mean: float = 0.0
    person_sd: float = 1.0
    criterion_difficulty_sd: float = 0.25
    rater_effects: RaterEffectConditionV1 = RaterEffectConditionV1()
    model: str = MODEL_RSM
    person_distribution: str = "normal"
    person_centering: str = "population_mean_no_sample_centering"
    criterion_difficulty_distribution: str = "normal"
    criterion_difficulty_sd_basis: str = "pre_centering_superpopulation_draw_sd"
    criterion_difficulty_constraint: str = "sample_mean_zero"
    score_record_role: str = "criterion"
    artifact_effect: str = "none"
    missingness: str = "none"
    common_anchor_invariance: str = "exact_by_construction"
    truth_sampling_scope: str = "all_levels_resampled_each_replicate"
    conditional_response_dependence: str = (
        "within_arm_unique_score_records_independent_given_truth"
    )
    response_probability_normalization: str = "logsumexp_stable"
    dgp_algorithm: str = "adjacent_category_rsm_rater_style_v1"
    schema_version: str = RSM_TRUTH_SPEC_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != RSM_TRUTH_SPEC_VERSION:
            raise SimulationConditionValidationError(
                f"Unsupported RSM truth version: {self.schema_version!r}"
            )
        _require_int(
            "rating_min",
            self.rating_min,
            minimum=-MAX_RATING_MIN_ABS,
            maximum=MAX_RATING_MIN_ABS,
        )
        _require_int("n_categories", self.n_categories, minimum=2, maximum=100)
        if not isinstance(self.adjacent_thresholds, tuple):
            raise SimulationConditionValidationError(
                "adjacent_thresholds must be a tuple"
            )
        if len(self.adjacent_thresholds) != self.n_categories - 1:
            raise SimulationConditionValidationError(
                "adjacent_thresholds must have length n_categories - 1"
            )
        thresholds = tuple(
            _require_real(
                f"adjacent_thresholds[{index}]",
                value,
                minimum=-MAX_ABS_ADJACENT_THRESHOLD,
                maximum=MAX_ABS_ADJACENT_THRESHOLD,
            )
            for index, value in enumerate(self.adjacent_thresholds)
        )
        if any(right <= left for left, right in zip(thresholds, thresholds[1:])):
            raise SimulationConditionValidationError(
                "adjacent_thresholds must be strictly increasing"
            )
        threshold_mean = math.fsum(thresholds) / len(thresholds)
        if not math.isclose(threshold_mean, 0.0, rel_tol=0.0, abs_tol=1e-12):
            raise SimulationConditionValidationError(
                "adjacent_thresholds must have arithmetic mean zero"
            )
        person_mean = _require_real("person_mean", self.person_mean)
        if person_mean != 0.0:
            raise SimulationConditionValidationError(
                "person_mean must remain 0.0 in the v1 truth contract"
            )
        _require_real(
            "person_sd",
            self.person_sd,
            maximum=MAX_GENERATING_STANDARD_DEVIATION,
            strictly_positive=True,
        )
        _require_real(
            "criterion_difficulty_sd",
            self.criterion_difficulty_sd,
            minimum=0.0,
            maximum=MAX_GENERATING_STANDARD_DEVIATION,
        )
        if type(self.rater_effects) is not RaterEffectConditionV1:
            raise TypeError("rater_effects must be a RaterEffectConditionV1")
        has_central_tendency = (
            self.rater_effects.central_tendency_log_scale_mean > 0
            or self.rater_effects.central_tendency_log_scale_sd > 0
        )
        if has_central_tendency and self.n_categories < 3:
            raise SimulationConditionValidationError(
                "central tendency requires at least three response categories"
            )
        _require_literal("model", self.model, MODEL_RSM)
        _require_literal("person_distribution", self.person_distribution, "normal")
        _require_literal(
            "person_centering",
            self.person_centering,
            "population_mean_no_sample_centering",
        )
        _require_literal(
            "criterion_difficulty_distribution",
            self.criterion_difficulty_distribution,
            "normal",
        )
        _require_literal(
            "criterion_difficulty_sd_basis",
            self.criterion_difficulty_sd_basis,
            "pre_centering_superpopulation_draw_sd",
        )
        _require_literal(
            "criterion_difficulty_constraint",
            self.criterion_difficulty_constraint,
            "sample_mean_zero",
        )
        _require_literal("score_record_role", self.score_record_role, "criterion")
        _require_literal("artifact_effect", self.artifact_effect, "none")
        _require_literal("missingness", self.missingness, "none")
        _require_literal(
            "common_anchor_invariance",
            self.common_anchor_invariance,
            "exact_by_construction",
        )
        _require_literal(
            "truth_sampling_scope",
            self.truth_sampling_scope,
            "all_levels_resampled_each_replicate",
        )
        _require_literal(
            "conditional_response_dependence",
            self.conditional_response_dependence,
            "within_arm_unique_score_records_independent_given_truth",
        )
        _require_literal(
            "response_probability_normalization",
            self.response_probability_normalization,
            "logsumexp_stable",
        )
        _require_literal(
            "dgp_algorithm", self.dgp_algorithm, "adjacent_category_rsm_rater_style_v1"
        )

    @property
    def rating_max(self) -> int:
        return self.rating_min + self.n_categories - 1

    @property
    def is_misspecification_stress(self) -> bool:
        return bool(
            self.rater_effects.central_tendency_log_scale_mean > 0
            or self.rater_effects.central_tendency_log_scale_sd > 0
            or self.rater_effects.rater_person_local_bias_sd > 0
        )

    def to_dict(self) -> dict:
        return _json_safe({
            "rating_min": self.rating_min,
            "n_categories": self.n_categories,
            "adjacent_thresholds": [float(value) for value in self.adjacent_thresholds],
            "person_mean": float(self.person_mean),
            "person_sd": float(self.person_sd),
            "criterion_difficulty_sd": float(self.criterion_difficulty_sd),
            "rater_effects": self.rater_effects.to_dict(),
            "model": self.model,
            "person_distribution": self.person_distribution,
            "person_centering": self.person_centering,
            "criterion_difficulty_distribution": (
                self.criterion_difficulty_distribution
            ),
            "criterion_difficulty_sd_basis": self.criterion_difficulty_sd_basis,
            "criterion_difficulty_constraint": self.criterion_difficulty_constraint,
            "score_record_role": self.score_record_role,
            "artifact_effect": self.artifact_effect,
            "missingness": self.missingness,
            "common_anchor_invariance": self.common_anchor_invariance,
            "truth_sampling_scope": self.truth_sampling_scope,
            "conditional_response_dependence": self.conditional_response_dependence,
            "response_probability_normalization": (
                self.response_probability_normalization
            ),
            "dgp_algorithm": self.dgp_algorithm,
            "schema_version": self.schema_version,
        })


def symmetric_adjacent_thresholds(
    n_categories: int,
    *,
    span: float = 2.0,
) -> tuple[float, ...]:
    """Return equally spaced, mean-zero adjacent thresholds."""
    count = _require_int("n_categories", n_categories, minimum=2, maximum=100) - 1
    width = _require_real(
        "span", span, maximum=MAX_THRESHOLD_SPAN, strictly_positive=True
    )
    if count == 1:
        return (0.0,)
    half = count // 2
    if count % 2:
        positive = tuple(
            width * index / (count - 1)
            for index in range(1, half + 1)
        )
        centered = tuple(-value for value in reversed(positive)) + (0.0,) + positive
    else:
        positive = tuple(
            width * (2 * index + 1) / (2 * (count - 1))
            for index in range(half)
        )
        centered = tuple(-value for value in reversed(positive)) + positive
    if any(right <= left for left, right in zip(centered, centered[1:])):
        raise SimulationConditionValidationError(
            "span is too small to produce strictly increasing binary64 thresholds"
        )
    if not math.isclose(
        math.fsum(centered) / count, 0.0, rel_tol=0.0, abs_tol=1e-12
    ):
        raise SimulationConditionValidationError(
            "span cannot produce mean-zero binary64 thresholds"
        )
    return centered


def normalize_rsm_truth_spec(value: RsmTruthSpecV1 | Mapping) -> RsmTruthSpecV1:
    if type(value) is RsmTruthSpecV1:
        return value
    payload = _strict_mapping(value, RsmTruthSpecV1, label="RSM truth spec")
    payload["adjacent_thresholds"] = _tuple_payload_field(
        payload, "adjacent_thresholds", label="RSM truth spec"
    )
    payload["rater_effects"] = normalize_rater_effect_condition(payload["rater_effects"])
    return RsmTruthSpecV1(**payload)


def rsm_truth_spec_fingerprint(
    value: RsmTruthSpecV1 | Mapping,
    *,
    length: int = 16,
) -> str:
    return _fingerprint(normalize_rsm_truth_spec(value).to_dict(), length=length)


__all__ = [
    "RATER_EFFECT_CONDITION_VERSION",
    "RSM_TRUTH_SPEC_VERSION",
    "RaterEffectConditionV1",
    "RsmTruthSpecV1",
    "normalize_rater_effect_condition",
    "normalize_rsm_truth_spec",
    "rater_effect_condition_fingerprint",
    "rsm_truth_spec_fingerprint",
    "symmetric_adjacent_thresholds",
]
