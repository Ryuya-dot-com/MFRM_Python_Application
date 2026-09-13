"""Structured, Streamlit-free orchestration for repository-only exact CMLE.

The direct :func:`mfrm_app.cmle.fit_cmle` research route retains technical
optimizer output on a conditional-support boundary for numerical auditability.
This module is the complementary user-workflow route: it stops before
optimization when no finite unpenalized CMLE is supported and always returns a
stable, bilingual, machine-readable result rather than a partially populated
fit object.
"""

from __future__ import annotations

import time
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from mfrm_app.cmle import fit_cmle, prepare_cmle_design
from mfrm_app.cmle_existence import audit_cmle_finite_mle_oracle


CMLE_WORKFLOW_SCHEMA_VERSION = "native_exact_cmle_structured_workflow_v1"
_STAGES = (
    "input_and_design",
    "finite_mle_existence",
    "structural_optimization",
    "downstream_person_scoring",
)


_TERMINAL_MESSAGES = {
    "input_invalid": {
        "headline_en": "The input could not be prepared for exact CMLE.",
        "headline_ja": "入力をexact CMLE用に準備できませんでした。",
        "detail_en": "A declared-support, column, value, duplicate-unit, or option check failed before model geometry was evaluated.",
        "detail_ja": "モデルの幾何を点検する前に、得点範囲、列、値、重複単位、または設定の検査で停止しました。",
        "next_en": "Correct the reported input error and run the same analysis again.",
        "next_ja": "表示された入力エラーを修正し、同じ解析を再実行してください。",
    },
    "design_not_identified": {
        "headline_en": "The conditional design is not identified.",
        "headline_ja": "条件付きデザインを同定できません。",
        "detail_en": "The exact prefit rank audit found an unsupported free direction or another blocking design issue. Optimization was not attempted.",
        "detail_ja": "exact事前ランク監査で、支持されない自由方向または別の重大なデザイン問題が見つかりました。最適化は実行していません。",
        "next_en": "Review informative overlap, disconnected components, category use, and the provenance of any proposed anchors.",
        "next_ja": "情報をもつ重なり、非連結成分、カテゴリ使用状況、および候補アンカーの根拠を確認してください。",
    },
    "finite_mle_boundary": {
        "headline_en": "The unpenalized conditional likelihood has no finite maximizer.",
        "headline_ja": "無罰則の条件付き尤度に有限の最大化解がありません。",
        "detail_en": "The observed sufficient statistic is on the conditional convex-support boundary. Large finite optimizer values would be technical approximations to a direction at infinity, so optimization was not attempted.",
        "detail_ja": "観測十分統計量が条件付き凸支持の境界上にあります。大きな有限の最適化値は無限方向の技術的近似にすぎないため、最適化は実行していません。",
        "next_en": "Inspect the implicated sparse response/anchor structure. A prior or penalty would define a different estimator and must be selected explicitly.",
        "next_ja": "関連するスパースな応答構造とアンカー構造を確認してください。事前分布や罰則を使う場合は別の推定量として明示的に選択する必要があります。",
    },
    "finite_mle_unavailable": {
        "headline_en": "Finite-CMLE existence could not be certified.",
        "headline_ja": "有限CMLEの存在を確認できませんでした。",
        "detail_en": "The support-oracle/LP audit was unavailable, unstable across tolerances, or exceeded a declared guard. Uncertainty was not converted into a finite result.",
        "detail_ja": "支持オラクル／LP監査が利用不能、許容誤差間で不安定、または規定上限超過でした。不確実な状態を有限解として扱っていません。",
        "next_en": "Export the audit telemetry and review numerical tolerances, work limits, and design complexity before retrying.",
        "next_ja": "監査情報を出力し、数値許容誤差、計算上限、デザインの複雑さを確認してから再試行してください。",
    },
    "optimization_not_ready": {
        "headline_en": "A finite CMLE is supported, but the fitted calibration is not ready.",
        "headline_ja": "有限CMLEは支持されますが、推定された較正は利用準備未完了です。",
        "detail_en": "The optimizer, stationarity, final information rank, positive-definiteness, or covariance gate did not pass.",
        "detail_ja": "最適化、停留性、最終情報ランク、正定値性、または共分散の判定を通過しませんでした。",
        "next_en": "Review the fit telemetry; do not interpret structural estimates or compute downstream Person results.",
        "next_ja": "fitの監査情報を確認し、構造推定値の解釈や下流のPerson結果の計算を行わないでください。",
    },
    "calibration_ready": {
        "headline_en": "Finite exact-CMLE structural calibration is ready for repository research use.",
        "headline_ja": "有限のexact CMLE構造較正がリポジトリ内研究利用の準備を完了しました。",
        "detail_en": "Identification, finite-existence, optimizer, stationarity, final-rank, and covariance gates passed. This does not validate fit thresholds or anchors.",
        "detail_ja": "同定、有限存在、最適化、停留性、最終ランク、共分散の各判定を通過しました。fit閾値やアンカーの妥当性を保証するものではありません。",
        "next_en": "Downstream fixed-calibration Person scoring may be run with its own availability, uncertainty, fit, and rounding diagnostics.",
        "next_ja": "下流の固定較正Person scoringは、利用可能性、不確実性、fit、丸め感度を別に点検した上で実行できます。",
    },
}


def _stage(
    order: int,
    stage: str,
    status: str,
    reached: bool,
    *,
    headline_en: str,
    headline_ja: str,
    detail_en: str,
    detail_ja: str,
    next_en: str = "",
    next_ja: str = "",
    evidence_code: str = "",
) -> dict[str, object]:
    return {
        "Order": int(order),
        "Stage": stage,
        "Status": status,
        "Reached": bool(reached),
        "HeadlineEn": headline_en,
        "HeadlineJa": headline_ja,
        "DetailEn": detail_en,
        "DetailJa": detail_ja,
        "NextActionEn": next_en,
        "NextActionJa": next_ja,
        "EvidenceCode": evidence_code,
    }


def _not_run_stage(order: int, stage: str) -> dict[str, object]:
    return _stage(
        order,
        stage,
        "not_run",
        False,
        headline_en="This stage was not run.",
        headline_ja="この段階は実行していません。",
        detail_en="An earlier workflow gate stopped the analysis.",
        detail_ja="前段階のworkflow判定により解析を停止しました。",
        evidence_code="earlier_stage_stopped",
    )


def _availability(
    *,
    terminal_status: str,
    existence_direction: bool,
    fit_returned: bool,
    fit_ready: bool,
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Artifact": "workflow_diagnostics",
                "Available": True,
                "Performed": True,
                "Reason": terminal_status,
            },
            {
                "Artifact": "boundary_direction",
                "Available": bool(existence_direction),
                "Performed": bool(existence_direction),
                "Reason": (
                    "reported_boundary_direction"
                    if existence_direction
                    else "no_reported_boundary_direction"
                ),
            },
            {
                "Artifact": "structural_calibration",
                "Available": bool(fit_returned),
                "Performed": bool(fit_returned),
                "Reason": "fit_returned" if fit_returned else terminal_status,
            },
            {
                "Artifact": "downstream_person_scoring",
                "Available": bool(fit_ready),
                "Performed": False,
                "Reason": (
                    "available_not_performed_in_v1"
                    if fit_ready
                    else "structural_calibration_not_ready"
                ),
            },
        ]
    )


def _finalize(
    *,
    terminal_status: str,
    stages: list[dict[str, object]],
    started: float,
    input_rows: int,
    design_audit: dict[str, object] | None,
    existence_audit: dict[str, object] | None,
    fit: dict[str, object] | None,
    error: str,
) -> dict[str, object]:
    message = _TERMINAL_MESSAGES[terminal_status]
    existence_summary = (
        existence_audit["summary"].iloc[0]
        if existence_audit is not None
        and isinstance(existence_audit.get("summary"), pd.DataFrame)
        and not existence_audit["summary"].empty
        else None
    )
    fit_summary = (
        fit["summary"].iloc[0]
        if fit is not None
        and isinstance(fit.get("summary"), pd.DataFrame)
        and not fit["summary"].empty
        else None
    )
    fit_returned = fit is not None
    optimization_attempted = terminal_status in {
        "optimization_not_ready",
        "calibration_ready",
    }
    fit_ready = bool(
        fit_summary is not None and fit_summary.get("InferenceReady", False)
    )
    direction = (
        existence_audit.get("direction")
        if existence_audit is not None
        else None
    )
    elapsed = float(time.perf_counter() - started)
    summary = pd.DataFrame(
        [
            {
                "SchemaVersion": CMLE_WORKFLOW_SCHEMA_VERSION,
                "WorkflowStatus": terminal_status,
                "Ready": terminal_status == "calibration_ready",
                "StoppedBeforeOptimization": not optimization_attempted,
                "OptimizationAttempted": optimization_attempted,
                "PersonScoringAvailable": fit_ready,
                "PersonScoringPerformed": False,
                "HeadlineEn": message["headline_en"],
                "HeadlineJa": message["headline_ja"],
                "DetailEn": message["detail_en"],
                "DetailJa": message["detail_ja"],
                "NextActionEn": message["next_en"],
                "NextActionJa": message["next_ja"],
                "InputRows": int(input_rows),
                "PrefitEligible": (
                    bool(design_audit.get("eligible", False))
                    if design_audit is not None
                    else False
                ),
                "PrefitRank": (
                    design_audit.get("conditional_rank", np.nan)
                    if design_audit is not None
                    else np.nan
                ),
                "PrefitNullity": (
                    design_audit.get("conditional_nullity", np.nan)
                    if design_audit is not None
                    else np.nan
                ),
                "ExistenceStatus": (
                    str(existence_summary.get("Status", "not_run"))
                    if existence_summary is not None
                    else "not_run"
                ),
                "ExistenceQualified": bool(
                    existence_summary is not None
                    and existence_summary.get("ExistenceQualified", False)
                ),
                "FitReturned": fit_returned,
                "FitInferenceReady": fit_ready,
                "Error": error,
                "ElapsedSeconds": elapsed,
            }
        ]
    )
    return {
        "schema_version": CMLE_WORKFLOW_SCHEMA_VERSION,
        "summary": summary,
        "stages": pd.DataFrame(stages),
        "availability": _availability(
            terminal_status=terminal_status,
            existence_direction=direction is not None,
            fit_returned=fit_returned,
            fit_ready=fit_ready,
        ),
        "design_audit": design_audit,
        "existence_audit": existence_audit,
        "fit": fit,
        "error": error,
    }


def run_cmle_calibration_workflow(
    data: pd.DataFrame,
    *,
    person_col: str,
    facet_cols: Sequence[str],
    score_col: str,
    rating_min: int,
    rating_max: int,
    model: str = "RSM",
    step_facet: str | None = None,
    weight_col: str | None = None,
    response_unit_col: str | None = None,
    positive_facets: Iterable[str] | None = None,
    hard_anchors: pd.DataFrame | Sequence[Mapping[str, object]] | None = None,
    allow_duplicate_units: bool = False,
    rank_audit_max_parameters: int = 200,
    rank_audit_max_work: int = 50_000_000,
    rank_audit_max_bytes: int = 512 * 1024 * 1024,
    maxiter: int = 500,
    gtol: float = 1e-7,
    newton_polish_maxiter: int = 8,
) -> dict[str, object]:
    """Run the staged exact-CMLE calibration workflow without silent fallback."""

    started = time.perf_counter()
    try:
        input_rows = len(data)
    except Exception:
        input_rows = 0
    stages = [_not_run_stage(index + 1, stage) for index, stage in enumerate(_STAGES)]
    prepare_kwargs = {
        "person_col": person_col,
        "facet_cols": facet_cols,
        "score_col": score_col,
        "rating_min": rating_min,
        "rating_max": rating_max,
        "model": model,
        "step_facet": step_facet,
        "weight_col": weight_col,
        "response_unit_col": response_unit_col,
        "positive_facets": positive_facets,
        "hard_anchors": hard_anchors,
        "allow_duplicate_units": allow_duplicate_units,
        "rank_audit_max_parameters": rank_audit_max_parameters,
        "rank_audit_max_work": rank_audit_max_work,
        "rank_audit_max_bytes": rank_audit_max_bytes,
    }
    try:
        design = prepare_cmle_design(data, **prepare_kwargs)
    except Exception as exc:
        error = f"{type(exc).__name__}: {str(exc)[:1000]}"
        stages[0] = _stage(
            1,
            _STAGES[0],
            "block",
            True,
            headline_en="Input preparation failed.",
            headline_ja="入力準備に失敗しました。",
            detail_en=error,
            detail_ja=error,
            next_en=_TERMINAL_MESSAGES["input_invalid"]["next_en"],
            next_ja=_TERMINAL_MESSAGES["input_invalid"]["next_ja"],
            evidence_code="input_exception",
        )
        return _finalize(
            terminal_status="input_invalid",
            stages=stages,
            started=started,
            input_rows=input_rows,
            design_audit=None,
            existence_audit=None,
            fit=None,
            error=error,
        )

    design_audit = design.audit
    if not bool(design_audit.get("eligible", False)):
        issue_codes = ";".join(
            design_audit.get("issues_table", pd.DataFrame())
            .get("Code", pd.Series(dtype=str))
            .astype(str)
        )
        stages[0] = _stage(
            1,
            _STAGES[0],
            "block",
            True,
            headline_en="The exact prefit design audit blocked estimation.",
            headline_ja="exact事前デザイン監査が推定を停止しました。",
            detail_en=f"Blocking issue codes: {issue_codes or 'unspecified'}.",
            detail_ja=f"重大な問題コード: {issue_codes or 'unspecified'}。",
            next_en=_TERMINAL_MESSAGES["design_not_identified"]["next_en"],
            next_ja=_TERMINAL_MESSAGES["design_not_identified"]["next_ja"],
            evidence_code=issue_codes or "prefit_ineligible",
        )
        return _finalize(
            terminal_status="design_not_identified",
            stages=stages,
            started=started,
            input_rows=input_rows,
            design_audit=design_audit,
            existence_audit=None,
            fit=None,
            error="",
        )
    stages[0] = _stage(
        1,
        _STAGES[0],
        "pass",
        True,
        headline_en="Input and exact prefit design audit passed.",
        headline_ja="入力とexact事前デザイン監査を通過しました。",
        detail_en=(
            f"Conditional rank/nullity: {design_audit['conditional_rank']}/"
            f"{design_audit['conditional_nullity']}."
        ),
        detail_ja=(
            f"条件付きランク／nullity: {design_audit['conditional_rank']}/"
            f"{design_audit['conditional_nullity']}。"
        ),
        evidence_code="prefit_eligible",
    )

    try:
        existence_audit = audit_cmle_finite_mle_oracle(design)
        existence_summary = existence_audit["summary"].iloc[0]
    except Exception as exc:
        error = f"{type(exc).__name__}: {str(exc)[:1000]}"
        stages[1] = _stage(
            2,
            _STAGES[1],
            "review",
            True,
            headline_en="The finite-MLE audit raised an exception.",
            headline_ja="有限MLE監査で例外が発生しました。",
            detail_en=error,
            detail_ja=error,
            next_en=_TERMINAL_MESSAGES["finite_mle_unavailable"]["next_en"],
            next_ja=_TERMINAL_MESSAGES["finite_mle_unavailable"]["next_ja"],
            evidence_code="existence_exception",
        )
        return _finalize(
            terminal_status="finite_mle_unavailable",
            stages=stages,
            started=started,
            input_rows=input_rows,
            design_audit=design_audit,
            existence_audit=None,
            fit=None,
            error=error,
        )

    existence_status = str(existence_summary["Status"])
    if existence_status == "boundary_no_finite_cmle":
        stages[1] = _stage(
            2,
            _STAGES[1],
            "block",
            True,
            headline_en="A conditional-support boundary was certified.",
            headline_ja="条件付き支持の境界が確認されました。",
            detail_en=str(existence_summary["Reason"]),
            detail_ja="観測十分統計量が条件付き凸支持の境界上にあります。",
            next_en=_TERMINAL_MESSAGES["finite_mle_boundary"]["next_en"],
            next_ja=_TERMINAL_MESSAGES["finite_mle_boundary"]["next_ja"],
            evidence_code=existence_status,
        )
        return _finalize(
            terminal_status="finite_mle_boundary",
            stages=stages,
            started=started,
            input_rows=input_rows,
            design_audit=design_audit,
            existence_audit=existence_audit,
            fit=None,
            error="",
        )
    if not bool(existence_summary.get("ExistenceQualified", False)):
        stages[1] = _stage(
            2,
            _STAGES[1],
            "review",
            True,
            headline_en="Finite-CMLE existence was not certified.",
            headline_ja="有限CMLEの存在を確認できませんでした。",
            detail_en=(
                f"Status/reason: {existence_status} / "
                f"{existence_summary.get('Reason', '')}."
            ),
            detail_ja=(
                f"状態／理由: {existence_status} / "
                f"{existence_summary.get('Reason', '')}。"
            ),
            next_en=_TERMINAL_MESSAGES["finite_mle_unavailable"]["next_en"],
            next_ja=_TERMINAL_MESSAGES["finite_mle_unavailable"]["next_ja"],
            evidence_code=existence_status,
        )
        return _finalize(
            terminal_status="finite_mle_unavailable",
            stages=stages,
            started=started,
            input_rows=input_rows,
            design_audit=design_audit,
            existence_audit=existence_audit,
            fit=None,
            error="",
        )
    stages[1] = _stage(
        2,
        _STAGES[1],
        "pass",
        True,
        headline_en="Finite-CMLE existence is supported.",
        headline_ja="有限CMLEの存在が支持されました。",
        detail_en="No nonzero supporting direction was found across the tolerance grid.",
        detail_ja="許容誤差グリッド全体で非零の支持方向は見つかりませんでした。",
        evidence_code=existence_status,
    )

    fit_kwargs = {
        **prepare_kwargs,
        "maxiter": maxiter,
        "gtol": gtol,
        "newton_polish_maxiter": newton_polish_maxiter,
    }
    try:
        fit = fit_cmle(data, **fit_kwargs)
    except Exception as exc:
        error = f"{type(exc).__name__}: {str(exc)[:1000]}"
        stages[2] = _stage(
            3,
            _STAGES[2],
            "block",
            True,
            headline_en="Structural optimization raised an exception.",
            headline_ja="構造最適化で例外が発生しました。",
            detail_en=error,
            detail_ja=error,
            next_en=_TERMINAL_MESSAGES["optimization_not_ready"]["next_en"],
            next_ja=_TERMINAL_MESSAGES["optimization_not_ready"]["next_ja"],
            evidence_code="fit_exception",
        )
        return _finalize(
            terminal_status="optimization_not_ready",
            stages=stages,
            started=started,
            input_rows=input_rows,
            design_audit=design_audit,
            existence_audit=existence_audit,
            fit=None,
            error=error,
        )

    fit_summary = fit["summary"].iloc[0]
    if not bool(fit_summary["InferenceReady"]):
        stages[2] = _stage(
            3,
            _STAGES[2],
            "block",
            True,
            headline_en="The structural fit did not pass readiness.",
            headline_ja="構造fitがreadiness判定を通過しませんでした。",
            detail_en=str(fit_summary.get("ReadinessReasons", "")),
            detail_ja=f"readiness理由: {fit_summary.get('ReadinessReasons', '')}。",
            next_en=_TERMINAL_MESSAGES["optimization_not_ready"]["next_en"],
            next_ja=_TERMINAL_MESSAGES["optimization_not_ready"]["next_ja"],
            evidence_code="fit_not_ready",
        )
        return _finalize(
            terminal_status="optimization_not_ready",
            stages=stages,
            started=started,
            input_rows=input_rows,
            design_audit=design_audit,
            existence_audit=existence_audit,
            fit=fit,
            error="",
        )
    stages[2] = _stage(
        3,
        _STAGES[2],
        "pass",
        True,
        headline_en="Structural optimization and readiness gates passed.",
        headline_ja="構造最適化とreadiness判定を通過しました。",
        detail_en=(
            f"Final rank/nullity: {fit_summary['InformationRank']}/"
            f"{fit_summary['InformationNullity']}; gradient sup norm "
            f"{float(fit_summary['GradientSupNorm']):.6g}."
        ),
        detail_ja=(
            f"最終ランク／nullity: {fit_summary['InformationRank']}/"
            f"{fit_summary['InformationNullity']}、gradient sup norm "
            f"{float(fit_summary['GradientSupNorm']):.6g}。"
        ),
        evidence_code="fit_inference_ready",
    )
    stages[3] = _stage(
        4,
        _STAGES[3],
        "available",
        True,
        headline_en="Downstream Person scoring is available but was not run.",
        headline_ja="下流のPerson scoringは利用可能ですが、まだ実行していません。",
        detail_en="Workflow v1 stops after structural calibration and records availability explicitly.",
        detail_ja="workflow v1は構造較正後に停止し、利用可能性のみを明示します。",
        next_en=_TERMINAL_MESSAGES["calibration_ready"]["next_en"],
        next_ja=_TERMINAL_MESSAGES["calibration_ready"]["next_ja"],
        evidence_code="person_scoring_available_not_performed",
    )
    return _finalize(
        terminal_status="calibration_ready",
        stages=stages,
        started=started,
        input_rows=input_rows,
        design_audit=design_audit,
        existence_audit=existence_audit,
        fit=fit,
        error="",
    )


__all__ = [
    "CMLE_WORKFLOW_SCHEMA_VERSION",
    "run_cmle_calibration_workflow",
]
