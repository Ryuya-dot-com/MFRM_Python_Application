"""Streamlit-free one-click result contract for repository-only exact CMLE.

This module joins already-qualified research components without changing their
estimands or numerical gates.  Exact structural CMLE must pass the staged
workflow before fixed-calibration Warm WLE and untrimmed Person MnSq are run.
The returned cards are a view model, not a public UI or comprehension claim.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from mfrm_app.cmle_wle_fit import compute_cmle_wle_person_fit
from mfrm_app.cmle_workflow import run_cmle_calibration_workflow
from mfrm_app.decision_stability import FIT_DISPLAY_DECIMALS


CMLE_ONE_CLICK_SCHEMA_VERSION = "native_exact_cmle_one_click_result_v1"
CMLE_ONE_CLICK_CARDS = (
    "input_and_design",
    "finite_mle_existence",
    "structural_calibration",
    "person_scoring",
    "person_fit_and_rounding",
    "uncertainty_and_scope",
)

_CARD_COLUMNS = (
    "Order",
    "Card",
    "Status",
    "HeadlineEn",
    "HeadlineJa",
    "DetailEn",
    "DetailJa",
    "NextActionEn",
    "NextActionJa",
    "EvidenceCode",
)


def _card(
    order: int,
    card: str,
    status: str,
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
        "Card": card,
        "Status": status,
        "HeadlineEn": headline_en,
        "HeadlineJa": headline_ja,
        "DetailEn": detail_en,
        "DetailJa": detail_ja,
        "NextActionEn": next_en,
        "NextActionJa": next_ja,
        "EvidenceCode": evidence_code,
    }


def summarize_person_fit_for_cards(persons: pd.DataFrame) -> dict[str, object]:
    """Summarize Person MnSq using raw classifications, never display values."""

    required = {
        "Person",
        "PersonFitReady",
        "WorstFitClass",
        "Infit",
        "Outfit",
        "InfitBoundaryStatus",
        "OutfitBoundaryStatus",
        "InfitDisplayDecisionConsistent",
        "OutfitDisplayDecisionConsistent",
        "ExtremeScorePattern",
    }
    if not isinstance(persons, pd.DataFrame) or not required.issubset(persons.columns):
        missing = sorted(required - set(getattr(persons, "columns", [])))
        raise ValueError(f"Person-fit table lacks required card fields: {missing}")
    if persons["Person"].astype(str).duplicated().any():
        raise ValueError("Person-fit card input must contain one row per Person.")

    ready = persons["PersonFitReady"].fillna(False).astype(bool)
    raw_class = persons["WorstFitClass"].fillna("unavailable").astype(str)
    review = ready & ~raw_class.eq("acceptable")
    infit_consistent = (
        persons["InfitDisplayDecisionConsistent"].fillna(False).astype(bool)
    )
    outfit_consistent = (
        persons["OutfitDisplayDecisionConsistent"].fillna(False).astype(bool)
    )
    boundary_values = pd.concat(
        [
            persons["InfitBoundaryStatus"].astype(str),
            persons["OutfitBoundaryStatus"].astype(str),
        ],
        ignore_index=True,
    )
    finite_raw = np.isfinite(
        persons[["Infit", "Outfit"]].apply(pd.to_numeric, errors="coerce")
    )
    return {
        "PersonCount": int(len(persons)),
        "PersonFitReadyCount": int(ready.sum()),
        "PersonFitUnavailableCount": int((~ready).sum()),
        "FitReviewPersonCount": int(review.sum()),
        "ExtremePersonCount": int(
            persons["ExtremeScorePattern"].fillna(False).astype(bool).sum()
        ),
        "RawDisplayMismatchCount": int(
            (~infit_consistent).sum() + (~outfit_consistent).sum()
        ),
        "NumericalBoundaryStatisticCount": int(
            boundary_values.eq("numerical_boundary").sum()
        ),
        "DisplayBoundaryStatisticCount": int(
            boundary_values.eq("display_rounding_boundary").sum()
        ),
        "FiniteRawMNSQCount": int(finite_raw.to_numpy().sum()),
        "DecisionsUseUnroundedMNSQ": True,
    }


def _workflow_card(stage: pd.Series, order: int, card_name: str) -> dict[str, object]:
    status_map = {
        "pass": "ready",
        "block": "blocked",
        "review": "review",
        "not_run": "not_run",
        "available": "available",
    }
    return _card(
        order,
        card_name,
        status_map.get(str(stage["Status"]), "review"),
        headline_en=str(stage["HeadlineEn"]),
        headline_ja=str(stage["HeadlineJa"]),
        detail_en=str(stage["DetailEn"]),
        detail_ja=str(stage["DetailJa"]),
        next_en=str(stage.get("NextActionEn", "")),
        next_ja=str(stage.get("NextActionJa", "")),
        evidence_code=str(stage.get("EvidenceCode", "")),
    )


def _fit_parameter_counts(calibration: dict[str, object]) -> tuple[int, int]:
    fit = calibration.get("fit")
    if not isinstance(fit, dict):
        return 0, 0
    coefficients = fit.get("coefficients")
    free_count = len(coefficients) if isinstance(coefficients, pd.DataFrame) else 0
    facets = fit.get("facets")
    others = facets.get("others") if isinstance(facets, dict) else None
    anchored_count = 0
    if isinstance(others, pd.DataFrame) and "Anchored" in others.columns:
        anchored_count = int(others["Anchored"].fillna(False).astype(bool).sum())
    return anchored_count, int(free_count)


def build_cmle_one_click_result(
    calibration: dict[str, object],
    person_fit: dict[str, object] | None,
    *,
    display_decimals: int = FIT_DISPLAY_DECIMALS,
    scoring_error: str = "",
) -> dict[str, object]:
    """Build the stable six-card result from retained workflow artifacts."""

    if not isinstance(calibration, dict):
        raise ValueError("calibration must be a structured CMLE workflow result.")
    workflow_summary = calibration.get("summary")
    workflow_stages = calibration.get("stages")
    if (
        not isinstance(workflow_summary, pd.DataFrame)
        or len(workflow_summary) != 1
        or not isinstance(workflow_stages, pd.DataFrame)
        or len(workflow_stages) != 4
    ):
        raise ValueError("calibration lacks the stable structured-workflow schema.")

    decimals = int(display_decimals)
    if decimals < 0:
        raise ValueError("display_decimals must be non-negative.")
    workflow_status = str(workflow_summary.iloc[0]["WorkflowStatus"])
    calibration_ready = workflow_status == "calibration_ready"
    cards = [
        _workflow_card(workflow_stages.iloc[index], index + 1, CMLE_ONE_CLICK_CARDS[index])
        for index in range(3)
    ]
    fit_metrics = {
        "PersonCount": 0,
        "PersonFitReadyCount": 0,
        "PersonFitUnavailableCount": 0,
        "FitReviewPersonCount": 0,
        "ExtremePersonCount": 0,
        "RawDisplayMismatchCount": 0,
        "NumericalBoundaryStatisticCount": 0,
        "DisplayBoundaryStatisticCount": 0,
        "FiniteRawMNSQCount": 0,
        "DecisionsUseUnroundedMNSQ": True,
    }
    person_scoring_attempted = bool(calibration_ready)
    person_scoring_ready = False

    if not calibration_ready:
        cards.append(
            _card(
                4,
                CMLE_ONE_CLICK_CARDS[3],
                "not_run",
                headline_en="Person scoring was not run.",
                headline_ja="Person scoringは実行していません。",
                detail_en="A calibration-stage gate stopped the analysis.",
                detail_ja="較正段階の判定により解析を停止しました。",
                next_en=str(workflow_summary.iloc[0].get("NextActionEn", "")),
                next_ja=str(workflow_summary.iloc[0].get("NextActionJa", "")),
                evidence_code="calibration_not_ready",
            )
        )
        cards.append(
            _card(
                5,
                CMLE_ONE_CLICK_CARDS[4],
                "not_run",
                headline_en="Person fit and rounding were not evaluated.",
                headline_ja="Person fitと丸め感度は評価していません。",
                detail_en="No reportable fixed-calibration Person result was available.",
                detail_ja="報告可能な固定較正Person結果がありません。",
                evidence_code="person_scoring_not_run",
            )
        )
        terminal_status = workflow_status
        summary_headline_en = str(workflow_summary.iloc[0]["HeadlineEn"])
        summary_headline_ja = str(workflow_summary.iloc[0]["HeadlineJa"])
        summary_next_en = str(workflow_summary.iloc[0].get("NextActionEn", ""))
        summary_next_ja = str(workflow_summary.iloc[0].get("NextActionJa", ""))
    elif scoring_error:
        cards.append(
            _card(
                4,
                CMLE_ONE_CLICK_CARDS[3],
                "blocked",
                headline_en="Structural calibration is ready, but Person scoring failed.",
                headline_ja="構造較正はreadyですが、Person scoringに失敗しました。",
                detail_en=scoring_error,
                detail_ja=scoring_error,
                next_en="Inspect the downstream scoring error; do not substitute another estimator automatically.",
                next_ja="下流scoringのエラーを確認し、別の推定量へ自動的に切り替えないでください。",
                evidence_code="person_scoring_exception",
            )
        )
        cards.append(
            _card(
                5,
                CMLE_ONE_CLICK_CARDS[4],
                "not_run",
                headline_en="Person fit and rounding were not evaluated.",
                headline_ja="Person fitと丸め感度は評価していません。",
                detail_en="The downstream scoring failure prevented Person-fit computation.",
                detail_ja="下流scoringの失敗によりPerson-fit計算を実行できませんでした。",
                evidence_code="person_scoring_failed",
            )
        )
        terminal_status = "downstream_person_scoring_not_ready"
        summary_headline_en = str(cards[3]["HeadlineEn"])
        summary_headline_ja = str(cards[3]["HeadlineJa"])
        summary_next_en = str(cards[3]["NextActionEn"])
        summary_next_ja = str(cards[3]["NextActionJa"])
    else:
        if not isinstance(person_fit, dict) or not isinstance(
            person_fit.get("persons"), pd.DataFrame
        ):
            raise ValueError("A ready calibration requires a retained Person-fit result.")
        persons = person_fit["persons"]
        fit_metrics = summarize_person_fit_for_cards(persons)
        person_scoring_ready = bool(
            fit_metrics["PersonCount"] > 0
            and fit_metrics["PersonFitUnavailableCount"] == 0
        )
        cards.append(
            _card(
                4,
                CMLE_ONE_CLICK_CARDS[3],
                "ready" if person_scoring_ready else "review",
                headline_en=(
                    "Fixed-calibration Warm WLE Person scoring is ready."
                    if person_scoring_ready
                    else "Some fixed-calibration Person scores are unavailable."
                ),
                headline_ja=(
                    "固定較正Warm WLE Person scoringはreadyです。"
                    if person_scoring_ready
                    else "一部の固定較正Person得点を利用できません。"
                ),
                detail_en=(
                    f"Ready Persons: {fit_metrics['PersonFitReadyCount']}/"
                    f"{fit_metrics['PersonCount']}; exact extremes: "
                    f"{fit_metrics['ExtremePersonCount']}."
                ),
                detail_ja=(
                    f"ready Person: {fit_metrics['PersonFitReadyCount']}/"
                    f"{fit_metrics['PersonCount']}、exact extreme: "
                    f"{fit_metrics['ExtremePersonCount']}。"
                ),
                next_en="Interpret these as WLE scores conditional on the fitted CMLE calibration.",
                next_ja="推定CMLE較正を固定した条件付きWLE得点として解釈してください。",
                evidence_code="fixed_calibration_wle",
            )
        )
        fit_caution = bool(
            fit_metrics["PersonFitUnavailableCount"]
            or fit_metrics["FitReviewPersonCount"]
            or fit_metrics["ExtremePersonCount"]
            or fit_metrics["RawDisplayMismatchCount"]
            or fit_metrics["NumericalBoundaryStatisticCount"]
            or fit_metrics["DisplayBoundaryStatisticCount"]
        )
        cards.append(
            _card(
                5,
                CMLE_ONE_CLICK_CARDS[4],
                "caution" if fit_caution else "ready",
                headline_en="Person fit uses unrounded Infit and Outfit.",
                headline_ja="Person fitは丸め前のInfitとOutfitで判定します。",
                detail_en=(
                    f"Raw fit-review Persons: {fit_metrics['FitReviewPersonCount']}; "
                    f"raw/display mismatches: {fit_metrics['RawDisplayMismatchCount']}; "
                    f"display-boundary statistics: "
                    f"{fit_metrics['DisplayBoundaryStatisticCount']}."
                ),
                detail_ja=(
                    f"生値でfit要確認のPerson: {fit_metrics['FitReviewPersonCount']}、"
                    f"生値／表示値の不一致: {fit_metrics['RawDisplayMismatchCount']}、"
                    f"表示丸め境界の統計量: "
                    f"{fit_metrics['DisplayBoundaryStatisticCount']}。"
                ),
                next_en="Inspect raw values and boundary rows; displayed values never drive a classification.",
                next_ja="生値と境界行を確認してください。表示値を分類には使用しません。",
                evidence_code="finite_unrounded_mnsq",
            )
        )
        terminal_status = (
            "analysis_ready_with_cautions"
            if person_scoring_ready
            else "downstream_person_scoring_not_ready"
        )
        summary_headline_en = (
            "Calibration and fixed-calibration Person results are ready with cautions."
            if person_scoring_ready
            else "Calibration is ready, but some downstream Person results are unavailable."
        )
        summary_headline_ja = (
            "較正と固定較正Person結果は注意事項付きでreadyです。"
            if person_scoring_ready
            else "較正はreadyですが、一部の下流Person結果を利用できません。"
        )
        summary_next_en = (
            "Review raw Person-fit flags, exact extremes, anchors, and the withheld uncertainty scope."
        )
        summary_next_ja = (
            "生値のPerson-fit判定、exact extreme、アンカー、保留中の不確実性範囲を確認してください。"
        )

    anchored_count, free_count = _fit_parameter_counts(calibration)
    anchor_text_en = (
        f" {anchored_count} anchored main-effect level(s) are fixed exactly; "
        "this does not validate their values."
        if anchored_count
        else " No hard main-effect anchor is active."
    )
    anchor_text_ja = (
        f" 主効果{anchored_count}水準をexactに固定していますが、その値の妥当性は"
        "保証しません。"
        if anchored_count
        else " 主効果のhard anchorは使用していません。"
    )
    cards[2]["DetailEn"] = str(cards[2]["DetailEn"]) + anchor_text_en
    cards[2]["DetailJa"] = str(cards[2]["DetailJa"]) + anchor_text_ja
    cards.append(
        _card(
            6,
            CMLE_ONE_CLICK_CARDS[5],
            "withheld",
            headline_en="Calibration uncertainty and public use remain withheld.",
            headline_ja="較正不確実性と公開利用は保留中です。",
            detail_en=(
                "Conditional WLE SE treats fitted CMLE calibration as fixed. "
                "ZSTD, p-values, confidence intervals, total SE, and public UI "
                "are not qualified."
            ),
            detail_ja=(
                "条件付きWLE SEは推定CMLE較正を固定値として扱います。ZSTD、p値、"
                "信頼区間、total SE、公開UIは未検証です。"
            ),
            next_en="Use the technical evidence only within the documented repository-research scope.",
            next_ja="文書化されたリポジトリ内研究の範囲でのみ技術的証拠を使用してください。",
            evidence_code="uncertainty_and_public_surface_withheld",
        )
    )
    cards_frame = pd.DataFrame(cards, columns=_CARD_COLUMNS)
    if cards_frame["Card"].tolist() != list(CMLE_ONE_CLICK_CARDS):
        raise RuntimeError("One-click card order changed unexpectedly.")

    error = scoring_error or str(calibration.get("error", ""))
    summary = pd.DataFrame(
        [
            {
                "SchemaVersion": CMLE_ONE_CLICK_SCHEMA_VERSION,
                "TerminalStatus": terminal_status,
                "CalibrationWorkflowStatus": workflow_status,
                "CalibrationReady": calibration_ready,
                "PersonScoringAttempted": person_scoring_attempted,
                "PersonScoringReady": person_scoring_ready,
                **fit_metrics,
                "AnchoredMainEffectLevelCount": anchored_count,
                "FreeCoordinateCount": free_count,
                "DisplayDecimals": decimals,
                "CalibrationUncertaintyPropagated": False,
                "AutomaticEstimatorFallback": False,
                "PublicSurfaceEnabled": False,
                "HeadlineEn": summary_headline_en,
                "HeadlineJa": summary_headline_ja,
                "NextActionEn": summary_next_en,
                "NextActionJa": summary_next_ja,
                "Error": error,
            }
        ]
    )
    availability = pd.DataFrame(
        [
            {
                "Artifact": "first_read_summary",
                "Available": True,
                "Performed": True,
                "PublicSurfaceEnabled": False,
            },
            {
                "Artifact": "calibration_workflow_diagnostics",
                "Available": True,
                "Performed": True,
                "PublicSurfaceEnabled": False,
            },
            {
                "Artifact": "structural_calibration",
                "Available": calibration.get("fit") is not None,
                "Performed": calibration.get("fit") is not None,
                "PublicSurfaceEnabled": False,
            },
            {
                "Artifact": "fixed_calibration_person_scoring",
                "Available": calibration_ready,
                "Performed": person_fit is not None,
                "PublicSurfaceEnabled": False,
            },
            {
                "Artifact": "person_fit_and_rounding_audit",
                "Available": person_fit is not None,
                "Performed": person_fit is not None,
                "PublicSurfaceEnabled": False,
            },
            {
                "Artifact": "calibration_uncertainty_total_se",
                "Available": False,
                "Performed": False,
                "PublicSurfaceEnabled": False,
            },
        ]
    )
    return {
        "schema_version": CMLE_ONE_CLICK_SCHEMA_VERSION,
        "summary": summary,
        "cards": cards_frame,
        "availability": availability,
        "calibration": calibration,
        "person_fit": person_fit,
        "error": error,
    }


def run_cmle_one_click_analysis(
    data: pd.DataFrame,
    *,
    display_decimals: int = FIT_DISPLAY_DECIMALS,
    **calibration_kwargs: Any,
) -> dict[str, object]:
    """Run one guarded repository-only analysis action without fallback."""

    calibration = run_cmle_calibration_workflow(data, **calibration_kwargs)
    person_fit: dict[str, object] | None = None
    scoring_error = ""
    workflow_status = str(calibration["summary"].iloc[0]["WorkflowStatus"])
    if workflow_status == "calibration_ready":
        try:
            person_fit = compute_cmle_wle_person_fit(
                calibration["fit"], display_decimals=int(display_decimals)
            )
        except Exception as exc:
            scoring_error = f"{type(exc).__name__}: {str(exc)[:1000]}"
    return build_cmle_one_click_result(
        calibration,
        person_fit,
        display_decimals=display_decimals,
        scoring_error=scoring_error,
    )


__all__ = [
    "CMLE_ONE_CLICK_CARDS",
    "CMLE_ONE_CLICK_SCHEMA_VERSION",
    "build_cmle_one_click_result",
    "run_cmle_one_click_analysis",
    "summarize_person_fit_for_cards",
]
