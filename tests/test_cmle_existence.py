"""Finite-CMLE existence and conditional-separation research contracts."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mfrm_app.cmle import prepare_cmle_design
from mfrm_app.cmle_existence import (
    audit_cmle_finite_mle,
    audit_cmle_finite_mle_oracle,
    cmle_score_stratum_support_maximum,
)


def _binary_two_rater(
    patterns: list[tuple[int, int]], *, anchor_r1: bool = False
):
    rows = []
    for index, (first, second) in enumerate(patterns):
        rows.extend(
            [
                (f"P{index:03d}", "R1", first),
                (f"P{index:03d}", "R2", second),
            ]
        )
    anchors = (
        pd.DataFrame(
            [
                {
                    "ParameterType": "Facet",
                    "Facet": "Rater",
                    "Level": "R1",
                    "Value": 0.0,
                }
            ]
        )
        if anchor_r1
        else None
    )
    return prepare_cmle_design(
        pd.DataFrame(rows, columns=["Person", "Rater", "Score"]),
        person_col="Person",
        facet_cols=["Rater"],
        score_col="Score",
        rating_min=0,
        rating_max=1,
        hard_anchors=anchors,
    )


@pytest.mark.parametrize("pattern", [(1, 0), (0, 1)])
@pytest.mark.parametrize("persons", [2, 5, 10, 20, 100, 300])
@pytest.mark.parametrize("anchor_r1", [False, True])
def test_complete_binary_separation_is_a_stable_support_boundary(
    pattern, persons, anchor_r1
):
    audit = audit_cmle_finite_mle(
        _binary_two_rater([pattern] * persons, anchor_r1=anchor_r1)
    )
    summary = audit["summary"].iloc[0]
    assert summary["Status"] == "boundary_no_finite_cmle"
    assert bool(summary["BoundaryDetected"])
    assert not bool(summary["ExistenceQualified"])
    assert audit["lp_tolerances"]["BoundaryDetected"].all()
    assert np.max(audit["inequalities"] @ audit["direction"]) <= 1e-9
    trace = audit["directional_trace"]
    assert trace["Finite"].all()
    assert bool(summary["DirectionalObjectiveNonincreasing"])
    assert trace["Objective"].iloc[-1] <= trace["Objective"].iloc[0]


@pytest.mark.parametrize(
    "patterns",
    [
        [(1, 0)] * 5 + [(0, 1)] * 5,
        [(1, 0)] * 9 + [(0, 1)],
        [(1, 0)] + [(0, 1)] * 9,
    ],
)
@pytest.mark.parametrize("anchor_r1", [False, True])
def test_minority_discordant_direction_retains_finite_interior(
    patterns, anchor_r1
):
    audit = audit_cmle_finite_mle(
        _binary_two_rater(patterns, anchor_r1=anchor_r1)
    )
    summary = audit["summary"].iloc[0]
    assert summary["Status"] == "interior_finite_cmle_supported"
    assert not bool(summary["BoundaryDetected"])
    assert bool(summary["FiniteMLESupported"])
    assert bool(summary["ExistenceQualified"])
    assert not audit["lp_tolerances"]["BoundaryDetected"].any()
    assert audit["direction"] is None


def test_structural_nonidentification_precedes_existence_classification():
    rows = []
    score_patterns = ([0, 1, 2], [1, 2, 0], [2, 0, 1], [0, 2, 1])
    for person_index in range(12):
        rater = "R1" if person_index < 6 else "R2"
        for criterion, score in zip(
            ["C1", "C2", "C3"], score_patterns[person_index % 4], strict=True
        ):
            rows.append((f"P{person_index}", rater, criterion, score))
    design = prepare_cmle_design(
        pd.DataFrame(rows, columns=["Person", "Rater", "Criterion", "Score"]),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
    )
    audit = audit_cmle_finite_mle(design)
    summary = audit["summary"].iloc[0]
    assert summary["Status"] == "structural_nonidentification"
    assert summary["Reason"] == "exact_conditional_information_rank_deficient"
    assert pd.isna(summary["BoundaryDetected"])
    assert audit["direction"] is None


@pytest.mark.parametrize(
    ("keyword", "value", "reason"),
    [
        ("max_configurations", 1, "configuration_cap_exceeded"),
        ("max_constraints", 1, "constraint_cap_exceeded"),
        ("max_bytes", 1, "memory_cap_exceeded"),
    ],
)
def test_existence_work_caps_fail_closed(keyword, value, reason):
    audit = audit_cmle_finite_mle(
        _binary_two_rater([(1, 0)] * 10), **{keyword: value}
    )
    summary = audit["summary"].iloc[0]
    assert summary["Status"] == "unavailable"
    assert summary["Reason"] == reason
    assert not bool(summary["ExistenceQualified"])


def test_existence_input_contract_rejects_invalid_tolerances_and_radii():
    design = _binary_two_rater([(1, 0), (0, 1)])
    with pytest.raises(ValueError, match="tolerances"):
        audit_cmle_finite_mle(design, tolerances=[1e-9, 1e-9])
    with pytest.raises(ValueError, match="directional_radii"):
        audit_cmle_finite_mle(design, directional_radii=[0.0, 2.0, 1.0])


@pytest.mark.parametrize(
    ("patterns", "expected"),
    [
        ([(1, 0)] * 20, "boundary_no_finite_cmle"),
        ([(0, 1)] * 20, "boundary_no_finite_cmle"),
        (
            [(1, 0)] * 10 + [(0, 1)] * 10,
            "interior_finite_cmle_supported",
        ),
        (
            [(1, 0)] * 19 + [(0, 1)],
            "interior_finite_cmle_supported",
        ),
    ],
)
@pytest.mark.parametrize("anchor_r1", [False, True])
def test_support_oracle_matches_exhaustive_binary_status(
    patterns, expected, anchor_r1
):
    design = _binary_two_rater(patterns, anchor_r1=anchor_r1)
    exhaustive = audit_cmle_finite_mle(design)
    oracle = audit_cmle_finite_mle_oracle(design)
    exhaustive_summary = exhaustive["summary"].iloc[0]
    oracle_summary = oracle["summary"].iloc[0]
    assert exhaustive_summary["Status"] == expected
    assert oracle_summary["Status"] == expected
    assert bool(oracle_summary["OracleComplete"])
    assert bool(oracle_summary["LPComplete"])
    assert (oracle["lp_tolerances"]["MaxNormalizedOracleViolation"] <= 1e-9).all()
    if expected == "boundary_no_finite_cmle":
        assert bool(oracle_summary["DirectionalObjectiveNonincreasing"])
        assert oracle["direction"] is not None
    else:
        assert bool(oracle_summary["ExistenceQualified"])
        assert oracle["direction"] is None


def test_fixed_score_support_oracle_matches_direct_enumeration():
    design = _binary_two_rater([(1, 0), (0, 1)])
    pattern = design.patterns[0]
    for direction in (
        np.ones(design.n_parameters),
        -np.ones(design.n_parameters),
        np.arange(1, design.n_parameters + 1, dtype=float),
    ):
        observed = cmle_score_stratum_support_maximum(
            design,
            pattern_index=0,
            raw_score=1,
            direction=direction,
        )
        direct = []
        for first in range(design.n_categories):
            second = 1 - first
            if not 0 <= second < design.n_categories:
                continue
            statistic = pattern.design[0, first] + pattern.design[1, second]
            direct.append(float(statistic @ direction))
        assert observed["Maximum"] == pytest.approx(max(direct), abs=1e-12)
        assert int(np.sum(observed["Categories"])) == 1


@pytest.mark.parametrize(
    ("keyword", "value", "reason"),
    [
        ("max_oracle_state_cells", 1, "oracle_state_cell_cap_exceeded"),
        ("max_bytes", 1, "oracle_memory_cap_exceeded"),
    ],
)
def test_support_oracle_preparation_caps_fail_closed(keyword, value, reason):
    audit = audit_cmle_finite_mle_oracle(
        _binary_two_rater([(1, 0)] * 10), **{keyword: value}
    )
    summary = audit["summary"].iloc[0]
    assert summary["Status"] == "unavailable"
    assert summary["Reason"] == reason
    assert not bool(summary["ExistenceQualified"])


def test_support_oracle_round_cap_fails_closed():
    audit = audit_cmle_finite_mle_oracle(
        _binary_two_rater([(1, 0)] * 10), max_cutting_plane_rounds=1
    )
    summary = audit["summary"].iloc[0]
    assert summary["Status"] == "oracle_unavailable"
    assert summary["Reason"] == "oracle_round_cap_exceeded"
    assert not bool(summary["ExistenceQualified"])
