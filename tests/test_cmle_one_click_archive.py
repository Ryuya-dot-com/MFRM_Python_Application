"""Deterministic private archive contracts for the CMLE one-click result."""

from __future__ import annotations

import io
import json
import warnings
import zipfile

import numpy as np
import pandas as pd
import pytest

from mfrm_app.cmle_one_click import run_cmle_one_click_analysis
from mfrm_app.cmle_one_click_archive import (
    PRIVATE_REPRODUCTION_MODE,
    build_cmle_one_click_archive,
    load_cmle_one_click_archive,
    replay_cmle_one_click_archive,
    verify_cmle_one_click_archive,
)


def _interior_frame() -> pd.DataFrame:
    scores = {
        "P1": [0, 1, 1, 2],
        "P2": [1, 0, 2, 1],
        "P3": [2, 1, 1, 0],
        "P4": [1, 2, 0, 1],
        "P5": [0, 2, 1, 1],
        "P6": [2, 0, 1, 1],
        "P7": [0, 0, 0, 0],
        "P8": [2, 2, 2, 2],
    }
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


def _kwargs(*, model: str = "RSM", binary: bool = False, anchors=None):
    values = {
        "person_col": "Person",
        "facet_cols": ["Rater"] if binary else ["Rater", "Criterion"],
        "score_col": "Score",
        "rating_min": 0,
        "rating_max": 1 if binary else 2,
        "model": model,
        "gtol": 1e-8,
        "maxiter": 800,
        "display_decimals": 3,
    }
    if model == "PCM":
        values["step_facet"] = "Criterion"
    if anchors is not None:
        values["hard_anchors"] = anchors
    return values


def _ready_archive(*, anchors=None, model="RSM"):
    frame = _interior_frame()
    kwargs = _kwargs(model=model, anchors=anchors)
    result = run_cmle_one_click_analysis(frame, **kwargs)
    return frame, kwargs, result, build_cmle_one_click_archive(
        result, input_data=frame, calibration_kwargs=kwargs
    )


def _rewrite_zip(assets: dict[str, bytes], *, duplicate: str | None = None) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, raw in assets.items():
            archive.writestr(name, raw)
        if duplicate is not None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                archive.writestr(duplicate, assets[duplicate])
    return buffer.getvalue()


def test_private_archive_is_byte_deterministic_and_round_trips_cards() -> None:
    frame, kwargs, result, first = _ready_archive()
    second = build_cmle_one_click_archive(
        result, input_data=frame, calibration_kwargs=dict(reversed(list(kwargs.items())))
    )
    assert first["zip_bytes"] == second["zip_bytes"]
    assert first["zip_sha256"] == second["zip_sha256"]
    assert first["analysis_identity"].analysis_id == second["analysis_identity"].analysis_id
    loaded = load_cmle_one_click_archive(first["zip_bytes"])
    cards = loaded["tables"]["result_cards.csv"]
    assert loaded["verified"]
    assert len(cards) == 6
    assert cards["HeadlineEn"].astype(str).str.len().gt(0).all()
    assert cards["HeadlineJa"].astype(str).str.len().gt(0).all()
    assert loaded["result_identity"]["terminal_status"] == "analysis_ready_with_cautions"
    loaded_persons = loaded["tables"]["person_results.csv"].sort_values("Person")
    original_persons = result["person_fit"]["persons"].sort_values("Person")
    for column in ("WLEEstimate", "ConditionalWLEStandardError", "Infit", "Outfit"):
        assert np.array_equal(
            loaded_persons[column].to_numpy(dtype=float),
            original_persons[column].to_numpy(dtype=float),
        )


def test_private_archive_replay_matches_all_retained_identities() -> None:
    _, _, _, archive = _ready_archive()
    replay = replay_cmle_one_click_archive(archive["zip_bytes"])
    assert replay["passed"]
    assert all(replay["checks"].values())


def test_boundary_archive_has_no_downstream_assets_and_replays() -> None:
    frame = _boundary_frame()
    kwargs = _kwargs(binary=True)
    result = run_cmle_one_click_analysis(frame, **kwargs)
    archive = build_cmle_one_click_archive(
        result, input_data=frame, calibration_kwargs=kwargs
    )
    verified = verify_cmle_one_click_archive(archive["zip_bytes"])
    assert verified["result_identity"]["terminal_status"] == "finite_mle_boundary"
    for name in (
        "structural_coefficients.csv",
        "facet_parameters.csv",
        "step_parameters.csv",
        "person_results.csv",
        "person_fit_decision_audit.csv",
    ):
        assert name not in verified["assets"]
    assert replay_cmle_one_click_archive(archive["zip_bytes"])["passed"]


def test_anchor_row_order_is_normalized_before_config_identity() -> None:
    anchors = pd.DataFrame(
        [
            {"ParameterType": "Facet", "Facet": "Rater", "Level": "R1", "Value": 0.25},
            {"ParameterType": "Facet", "Facet": "Criterion", "Level": "C1", "Value": -0.2},
        ]
    )
    frame, kwargs, result, first = _ready_archive(anchors=anchors, model="PCM")
    reversed_kwargs = {**kwargs, "hard_anchors": anchors.iloc[::-1].reset_index(drop=True)}
    second = build_cmle_one_click_archive(
        result, input_data=frame, calibration_kwargs=reversed_kwargs
    )
    assert first["assets"]["hard_anchors.csv"] == second["assets"]["hard_anchors.csv"]
    assert first["analysis_identity"].analysis_id == second["analysis_identity"].analysis_id
    assert first["zip_sha256"] == second["zip_sha256"]
    assert replay_cmle_one_click_archive(first["zip_bytes"])["passed"]


def test_row_permutation_preserves_semantic_but_changes_ordered_identity() -> None:
    frame, kwargs, result, first = _ready_archive()
    permuted = frame.sample(frac=1.0, random_state=20260810).reset_index(drop=True)
    second = build_cmle_one_click_archive(
        result, input_data=permuted, calibration_kwargs=kwargs
    )
    assert first["semantic_input_sha256"] == second["semantic_input_sha256"]
    assert first["ordered_input_sha256"] != second["ordered_input_sha256"]
    assert first["analysis_identity"].analysis_id != second["analysis_identity"].analysis_id


@pytest.mark.parametrize("mutation", ["change", "delete", "extra", "duplicate"])
def test_archive_tampering_fails_before_loading(mutation: str) -> None:
    _, _, _, archive = _ready_archive()
    assets = dict(archive["assets"])
    duplicate = None
    if mutation == "change":
        assets["result_cards.csv"] = assets["result_cards.csv"].replace(b"ready", b"green", 1)
    elif mutation == "delete":
        del assets["result_cards.csv"]
    elif mutation == "extra":
        assets["undeclared.txt"] = b"not declared"
    else:
        duplicate = "result_cards.csv"
    payload = _rewrite_zip(assets, duplicate=duplicate)
    with pytest.raises(ValueError):
        verify_cmle_one_click_archive(payload)


def test_public_archive_mode_fails_closed_and_private_readme_warns() -> None:
    frame, kwargs, result, archive = _ready_archive()
    with pytest.raises(ValueError, match="Public CMLE one-click export is withheld"):
        build_cmle_one_click_archive(
            result,
            input_data=frame,
            calibration_kwargs=kwargs,
            privacy_mode="public",
        )
    readme = archive["assets"]["README_FIRST.md"].decode("utf-8")
    manifest = json.loads(archive["assets"]["archive_manifest.json"])
    assert "private controlled-access reproduction artifact" in readme
    assert manifest["privacy_mode"] == PRIVATE_REPRODUCTION_MODE
    assert not manifest["public_surface_enabled"]
