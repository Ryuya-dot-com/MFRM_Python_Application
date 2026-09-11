"""Contracts for the repository-only CMLE one-click result view model."""

from __future__ import annotations

import numpy as np
import pandas as pd

import mfrm_app.cmle_one_click as one_click
from mfrm_app.cmle_wle_fit import (
    compute_cmle_wle_person_fit,
    compute_fixed_calibration_person_fit,
)


def _interior_frame(*, extremes: bool = True) -> pd.DataFrame:
    scores = {
        "P1": [0, 1, 1, 2],
        "P2": [1, 0, 2, 1],
        "P3": [2, 1, 1, 0],
        "P4": [1, 2, 0, 1],
        "P5": [0, 2, 1, 1],
        "P6": [2, 0, 1, 1],
    }
    if extremes:
        scores.update({"P7": [0, 0, 0, 0], "P8": [2, 2, 2, 2]})
    units = [("R1", "C1"), ("R1", "C2"), ("R2", "C1"), ("R2", "C2")]
    return pd.DataFrame(
        [
            (person, rater, criterion, score)
            for person, values in scores.items()
            for (rater, criterion), score in zip(units, values, strict=True)
        ],
        columns=["Person", "Rater", "Criterion", "Score"],
    )


def _boundary_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            (f"P{index:03d}", rater, score)
            for index in range(20)
            for rater, score in (("R1", 1), ("R2", 0))
        ],
        columns=["Person", "Rater", "Score"],
    )


def _structural_frame() -> pd.DataFrame:
    rows = []
    patterns = ([0, 1, 2], [1, 2, 0], [2, 0, 1], [0, 2, 1])
    for person_index in range(12):
        rater = "R1" if person_index < 6 else "R2"
        for criterion, score in zip(
            ["C1", "C2", "C3"], patterns[person_index % 4], strict=True
        ):
            rows.append((f"P{person_index}", rater, criterion, score))
    return pd.DataFrame(rows, columns=["Person", "Rater", "Criterion", "Score"])


def _run(
    frame: pd.DataFrame,
    *,
    model: str = "RSM",
    binary: bool = False,
    anchors: pd.DataFrame | None = None,
) -> dict[str, object]:
    return one_click.run_cmle_one_click_analysis(
        frame,
        person_col="Person",
        facet_cols=["Rater"] if binary else ["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=1 if binary else 2,
        model=model,
        step_facet="Criterion" if model == "PCM" else None,
        hard_anchors=anchors,
        gtol=1e-8,
        maxiter=800,
        display_decimals=3,
    )


def _assert_cards(result: dict[str, object]) -> None:
    assert result["schema_version"] == one_click.CMLE_ONE_CLICK_SCHEMA_VERSION
    cards = result["cards"]
    assert cards["Card"].tolist() == list(one_click.CMLE_ONE_CLICK_CARDS)
    assert cards["Order"].tolist() == list(range(1, 7))
    for language_column in (
        "HeadlineEn",
        "HeadlineJa",
        "DetailEn",
        "DetailJa",
    ):
        assert cards[language_column].astype(str).str.len().gt(0).all()
    assert cards.iloc[-1]["Status"] == "withheld"
    assert not result["availability"]["PublicSurfaceEnabled"].any()


def test_ready_rsm_runs_wle_fit_and_matches_direct_handoff() -> None:
    result = _run(_interior_frame())
    _assert_cards(result)
    summary = result["summary"].iloc[0]
    assert summary["TerminalStatus"] == "analysis_ready_with_cautions"
    assert bool(summary["CalibrationReady"])
    assert bool(summary["PersonScoringAttempted"])
    assert bool(summary["PersonScoringReady"])
    assert summary["ExtremePersonCount"] == 2
    assert bool(summary["DecisionsUseUnroundedMNSQ"])
    assert "ready with cautions" in summary["HeadlineEn"]
    assert "注意事項付き" in summary["HeadlineJa"]
    direct = compute_cmle_wle_person_fit(result["calibration"]["fit"])["persons"]
    observed = result["person_fit"]["persons"]
    paired = observed.merge(direct, on="Person", suffixes=("Observed", "Direct"))
    for column in ("WLEEstimate", "ConditionalWLEStandardError", "Infit", "Outfit"):
        assert np.max(
            np.abs(paired[f"{column}Observed"] - paired[f"{column}Direct"])
        ) < 1e-12


def test_ready_pcm_reports_exact_differential_anchors_and_warning() -> None:
    anchors = pd.DataFrame(
        [
            {"ParameterType": "Facet", "Facet": "Rater", "Level": "R1", "Value": 0.25},
            {"ParameterType": "Facet", "Facet": "Criterion", "Level": "C1", "Value": -0.2},
        ]
    )
    result = _run(_interior_frame(extremes=False), model="PCM", anchors=anchors)
    _assert_cards(result)
    summary = result["summary"].iloc[0]
    assert summary["TerminalStatus"] == "analysis_ready_with_cautions"
    assert summary["AnchoredMainEffectLevelCount"] == 2
    facets = result["calibration"]["fit"]["facets"]["others"]
    anchored = facets.loc[facets["Anchored"]].sort_values(["Facet", "Level"])
    assert anchored["Estimate"].tolist() == [-0.2, 0.25]
    assert anchored["SE"].tolist() == [0.0, 0.0]
    assert "does not validate" in result["cards"].iloc[2]["DetailEn"]
    assert "妥当性は保証しません" in result["cards"].iloc[2]["DetailJa"]


def test_blocked_calibration_never_attempts_person_scoring(monkeypatch) -> None:
    def forbidden(*_args, **_kwargs):
        raise AssertionError("Person scoring must not run after a calibration block")

    monkeypatch.setattr(one_click, "compute_cmle_wle_person_fit", forbidden)
    for result, expected in (
        (_run(_boundary_frame(), binary=True), "finite_mle_boundary"),
        (_run(_structural_frame()), "design_not_identified"),
        (_run(_interior_frame().drop(columns=["Score"])), "input_invalid"),
    ):
        _assert_cards(result)
        summary = result["summary"].iloc[0]
        assert summary["TerminalStatus"] == expected
        assert summary["HeadlineEn"] == result["calibration"]["summary"].iloc[0][
            "HeadlineEn"
        ]
        assert not bool(summary["PersonScoringAttempted"])
        assert result["person_fit"] is None
        assert result["cards"].iloc[3]["Status"] == "not_run"
        assert result["cards"].iloc[4]["Status"] == "not_run"


def test_downstream_exception_is_typed_without_estimator_fallback(monkeypatch) -> None:
    def failed(*_args, **_kwargs):
        raise FloatingPointError("injected downstream failure")

    monkeypatch.setattr(one_click, "compute_cmle_wle_person_fit", failed)
    result = _run(_interior_frame())
    _assert_cards(result)
    summary = result["summary"].iloc[0]
    assert summary["TerminalStatus"] == "downstream_person_scoring_not_ready"
    assert bool(summary["CalibrationReady"])
    assert bool(summary["PersonScoringAttempted"])
    assert not bool(summary["PersonScoringReady"])
    assert not bool(summary["AutomaticEstimatorFallback"])
    assert result["calibration"]["fit"] is not None
    assert "injected downstream failure" in result["error"]
    assert "Person scoring failed" in summary["HeadlineEn"]


def test_rounding_probe_counts_raw_display_mismatch_without_reclassification() -> None:
    raw_infit = 1.5004
    probe = compute_fixed_calibration_person_fit(
        ["P1"],
        [0],
        np.array([[0.0, np.log(raw_infit)]]),
        [0.0],
        display_decimals=3,
    )["persons"]
    probe["ExtremeScorePattern"] = False
    metrics = one_click.summarize_person_fit_for_cards(probe)
    row = probe.iloc[0]
    assert abs(row["Infit"] - raw_infit) < 1e-12
    assert row["InfitClass"] == "noisy"
    assert row["InfitDisplay"] == 1.5
    assert not bool(row["InfitDisplayDecisionConsistent"])
    assert metrics["RawDisplayMismatchCount"] == 2
    assert metrics["FitReviewPersonCount"] == 1
    assert bool(metrics["DecisionsUseUnroundedMNSQ"])


def test_person_card_summary_rejects_duplicate_identity() -> None:
    probe = compute_fixed_calibration_person_fit(
        ["P1"], [0], np.array([[0.0, 0.0]]), [0.0]
    )["persons"]
    probe["ExtremeScorePattern"] = False
    duplicate = pd.concat([probe, probe], ignore_index=True)
    try:
        one_click.summarize_person_fit_for_cards(duplicate)
    except ValueError as exc:
        assert "one row per Person" in str(exc)
    else:
        raise AssertionError("Duplicate Person identity must fail closed")
