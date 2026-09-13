from __future__ import annotations

import pandas as pd
import pandas.testing as pdt

from mfrm_app.assignment_sensitivity import (
    build_degree_preserving_assignment_perturbation,
    build_perturbation_path_snapshots,
    evaluate_fixed_density_perturbation_feasibility,
)


def _exchangeable_sparse_design() -> pd.DataFrame:
    assignment = {
        "p0": ("r2", "r3"),
        "p1": ("r2", "r3"),
        "p2": ("r1", "r3"),
        "p3": ("r1", "r2"),
        "p4": ("r0", "r3"),
        "p5": ("r0", "r2"),
        "p6": ("r0", "r1"),
        "p7": ("r0", "r1"),
    }
    rows = []
    for person, raters in assignment.items():
        for rater in raters:
            for task in ("t1", "t2"):
                rows.append(
                    {
                        "Person": person,
                        "Rater": rater,
                        "Task": task,
                        "Score": (int(person[1:]) + int(rater[1:])) % 4,
                        "Weight": 1.0,
                    }
                )
    return pd.DataFrame(rows)


PERSON_SCORES = {f"p{idx}": float(idx) for idx in range(8)}
RATER_SCORES = {f"r{idx}": float(idx) for idx in range(4)}


def test_exchangeable_sparse_design_passes_strict_fixed_density_gate():
    audit = evaluate_fixed_density_perturbation_feasibility(
        _exchangeable_sparse_design(),
        person_col="Person",
        rater_col="Rater",
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )

    assert audit["available"] is True
    assert audit["gates"]["Passed"].all()
    assert audit["block_profiles"]["Rows"].eq(2).all()
    assert audit["block_profiles"]["ContextSignature"].nunique() == 1


def test_aligned_perturbation_preserves_degrees_rows_exposure_and_connectivity():
    original = _exchangeable_sparse_design()
    bundle = build_degree_preserving_assignment_perturbation(
        original,
        person_scores=PERSON_SCORES,
        rater_scores=RATER_SCORES,
        person_col="Person",
        rater_col="Rater",
        facet_cols=["Rater", "Task"],
        context_cols=["Task"],
        rater_mapping_confirmed=True,
        direction="aligned",
    )

    assert bundle["available"] is True
    assert bundle["invariants"]["Passed"].all()
    assert bundle["assignment_map"]["Changed"].any()
    assert "Score" not in bundle["design"].columns
    assert bundle["design"]["_SourceRow"].tolist() == list(range(len(original)))
    trajectory = bundle["trajectory"]
    assert trajectory.iloc[-1]["Objective"] > trajectory.iloc[0]["Objective"]
    assert trajectory.iloc[-1]["AssignmentRankCorrelation"] > trajectory.iloc[0]["AssignmentRankCorrelation"]

    before_edges = original[["Person", "Rater"]].drop_duplicates()
    after_edges = bundle["design"][["Person", "Rater"]].drop_duplicates()
    pdt.assert_series_equal(
        before_edges.groupby("Person")["Rater"].nunique().sort_index(),
        after_edges.groupby("Person")["Rater"].nunique().sort_index(),
    )
    pdt.assert_series_equal(
        original.groupby("Rater").size().sort_index(),
        bundle["design"].groupby("Rater").size().sort_index(),
    )


def test_perturbation_is_deterministic_and_does_not_depend_on_observed_scores():
    original = _exchangeable_sparse_design()
    changed_scores = original.copy()
    changed_scores["Score"] = range(len(changed_scores))
    kwargs = dict(
        person_scores=PERSON_SCORES,
        rater_scores=RATER_SCORES,
        facet_cols=["Rater", "Task"],
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )
    first = build_degree_preserving_assignment_perturbation(original, **kwargs)
    second = build_degree_preserving_assignment_perturbation(changed_scores, **kwargs)

    pdt.assert_frame_equal(first["design"], second["design"])
    pdt.assert_frame_equal(first["assignment_map"], second["assignment_map"])
    pdt.assert_frame_equal(first["switch_ledger"], second["switch_ledger"])


def test_anti_aligned_direction_reduces_assignment_rank_correlation():
    original = _exchangeable_sparse_design()
    # Reverse the Rater coordinate so the observed design has room to move in
    # the requested anti-aligned direction.
    bundle = build_degree_preserving_assignment_perturbation(
        original,
        person_scores=PERSON_SCORES,
        rater_scores={key: -value for key, value in RATER_SCORES.items()},
        facet_cols=["Rater", "Task"],
        context_cols=["Task"],
        rater_mapping_confirmed=True,
        direction="anti_aligned",
    )
    assert bundle["available"] is True
    trajectory = bundle["trajectory"]
    assert trajectory.iloc[-1]["AssignmentRankCorrelation"] < trajectory.iloc[0]["AssignmentRankCorrelation"]


def test_complete_design_is_blocked_because_no_degree_preserving_contrast_exists():
    rows = [
        {"Person": person, "Rater": rater, "Task": "t1", "Score": 1}
        for person in ("p1", "p2", "p3")
        for rater in ("r1", "r2", "r3")
    ]
    audit = evaluate_fixed_density_perturbation_feasibility(
        pd.DataFrame(rows),
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )

    assert audit["available"] is False
    switch_gate = audit["gates"].loc[
        audit["gates"]["Gate"].eq("At least one connected degree-preserving 2-switch")
    ].iloc[0]
    assert not switch_gate["Passed"]


def test_unequal_context_blocks_and_frequency_weights_fail_closed():
    unequal = _exchangeable_sparse_design().copy()
    unequal = unequal.loc[
        ~(
            unequal["Person"].eq("p0")
            & unequal["Rater"].eq("r2")
            & unequal["Task"].eq("t2")
        )
    ].copy()
    audit_unequal = evaluate_fixed_density_perturbation_feasibility(
        unequal,
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )
    assert audit_unequal["available"] is False
    assert "Common context signature" in audit_unequal["reason"]

    weighted = _exchangeable_sparse_design().copy()
    weighted.loc[0, "Weight"] = 2.0
    audit_weighted = evaluate_fixed_density_perturbation_feasibility(
        weighted,
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )
    assert audit_weighted["available"] is False
    assert "Uncompressed row topology" in audit_weighted["reason"]


def test_unconfirmed_rater_mapping_fails_closed_even_when_topology_is_eligible():
    audit = evaluate_fixed_density_perturbation_feasibility(
        _exchangeable_sparse_design(),
        context_cols=["Task"],
        rater_mapping_confirmed=False,
    )
    assert audit["available"] is False
    assert "Confirmed Rater role" in audit["reason"]


def test_perturbation_core_retains_required_rater_and_context_without_facet_list():
    bundle = build_degree_preserving_assignment_perturbation(
        _exchangeable_sparse_design(),
        person_scores=PERSON_SCORES,
        rater_scores=RATER_SCORES,
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )

    assert bundle["available"] is True
    assert {"Person", "Rater", "Task", "_SourceRow"}.issubset(bundle["design"].columns)
    assert bundle["invariants"]["Passed"].all()


def test_path_snapshots_are_score_free_monotone_and_preserve_required_invariants():
    original = _exchangeable_sparse_design()
    perturbation = build_degree_preserving_assignment_perturbation(
        original,
        person_scores=PERSON_SCORES,
        rater_scores=RATER_SCORES,
        facet_cols=["Rater", "Task"],
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )
    path = build_perturbation_path_snapshots(
        original,
        perturbation,
        requested_doses=[0.0, 0.25, 0.5, 0.75, 1.0],
        facet_cols=["Rater", "Task"],
        context_cols=["Task"],
    )

    assert path["available"] is True
    dose = path["dose_table"]
    assert dose.iloc[0]["Scenario"] == "observed_assignment"
    assert dose.iloc[-1]["Scenario"] == "aligned_counterfactual"
    assert dose.iloc[0]["AchievedAlignmentDose"] == 0.0
    assert dose.iloc[-1]["AchievedAlignmentDose"] == 1.0
    assert dose["AchievedAlignmentDose"].is_monotonic_increasing
    assert dose["SwitchesApplied"].is_monotonic_increasing
    assert dose["RequestedDoseCount"].sum() == 5
    assert set(path["designs"]) == set(dose["Scenario"])
    assert all("Score" not in design.columns for design in path["designs"].values())
    required = path["path_invariants"].loc[
        path["path_invariants"]["RequiredAtDose"]
    ]
    assert required["PassedForDose"].all()
    pdt.assert_frame_equal(
        path["designs"]["aligned_counterfactual"].reset_index(drop=True),
        perturbation["design"].reset_index(drop=True),
    )
    assert {"Scenario", "Person", "SourceRater", "CounterfactualRater"}.issubset(
        path["path_assignment_map"].columns
    )
