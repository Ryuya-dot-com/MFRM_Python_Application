from __future__ import annotations

import pandas as pd
import pandas.testing as pdt

from mfrm_app.assignment_context_milp import (
    build_context_margin_assignment_perturbation,
    build_context_margin_endpoint_pair,
    evaluate_context_margin_milp_feasibility,
)
from mfrm_app.assignment_sensitivity import (
    evaluate_fixed_density_perturbation_feasibility,
)


def _unequal_context_design() -> pd.DataFrame:
    # Each Person has one one-row block and one two-row block.  Every Rater has
    # the same aggregate Task margins, but Person-Rater block signatures differ.
    assignment = {
        "p0": (("r0", "A"), ("r1", "B")),
        "p1": (("r0", "A"), ("r2", "B")),
        "p2": (("r1", "A"), ("r0", "B")),
        "p3": (("r1", "A"), ("r2", "B")),
        "p4": (("r2", "A"), ("r0", "B")),
        "p5": (("r2", "A"), ("r1", "B")),
    }
    rows = []
    for person, blocks in assignment.items():
        for rater, signature in blocks:
            rows.append(
                {"Person": person, "Rater": rater, "Task": "t1", "Score": 1, "Weight": 1.0}
            )
            if signature == "B":
                rows.append(
                    {"Person": person, "Rater": rater, "Task": "t2", "Score": 2, "Weight": 1.0}
                )
    return pd.DataFrame(rows)


PERSON_SCORES = {f"p{index}": float(index) for index in range(6)}
RATER_SCORES = {f"r{index}": float(index) for index in range(3)}


def test_generalized_gate_accepts_design_rejected_by_equal_signature_runner():
    data = _unequal_context_design()
    strict = evaluate_fixed_density_perturbation_feasibility(
        data,
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )
    generalized = evaluate_context_margin_milp_feasibility(
        data,
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )

    assert strict["available"] is False
    assert "Common context signature" in strict["reason"]
    assert generalized["available"] is True
    assert generalized["gates"]["Passed"].all()
    assert len(generalized["witnesses"]) == 2
    assert generalized["block_profiles"]["ContextCountSignature"].nunique() == 2


def test_milp_endpoint_preserves_exact_context_margins_and_connectivity():
    original = _unequal_context_design()
    bundle = build_context_margin_assignment_perturbation(
        original,
        person_scores=PERSON_SCORES,
        rater_scores=RATER_SCORES,
        facet_cols=["Rater", "Task"],
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )

    assert bundle["available"] is True
    assert bundle["solver_audit"].iloc[0]["ObservedBaselineFeasible"]
    assert bundle["solver_audit"].iloc[0]["GlobalOptimumCertified"]
    assert bundle["solver_audit"].iloc[0]["MIPGap"] == 0.0
    assert bundle["invariants"]["Passed"].all()
    assert bundle["context_margin_audit"]["ExactMatch"].all()
    assert bundle["witnesses"]["LockSatisfied"].all()
    assert bundle["assignment_map"]["Changed"].any()
    assert "Score" not in bundle["design"].columns
    assert bundle["design"]["_SourceRow"].tolist() == list(range(len(original)))
    trajectory = bundle["trajectory"]
    assert trajectory.iloc[-1]["DirectionAdjustedObjective"] > trajectory.iloc[0]["DirectionAdjustedObjective"]

    pdt.assert_series_equal(
        original.groupby(["Rater", "Task"]).size().sort_index(),
        bundle["design"].groupby(["Rater", "Task"]).size().sort_index(),
    )


def test_milp_design_is_independent_of_observed_scores():
    original = _unequal_context_design()
    changed = original.copy()
    changed["Score"] = range(len(changed))
    kwargs = dict(
        person_scores=PERSON_SCORES,
        rater_scores=RATER_SCORES,
        facet_cols=["Rater", "Task"],
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )
    first = build_context_margin_assignment_perturbation(original, **kwargs)
    second = build_context_margin_assignment_perturbation(changed, **kwargs)

    pdt.assert_frame_equal(first["design"], second["design"])
    pdt.assert_frame_equal(first["assignment_map"], second["assignment_map"])
    pdt.assert_frame_equal(first["trajectory"], second["trajectory"])


def test_endpoint_pair_exposes_only_observed_and_optimized_designs():
    original = _unequal_context_design()
    perturbation = build_context_margin_assignment_perturbation(
        original,
        person_scores=PERSON_SCORES,
        rater_scores=RATER_SCORES,
        facet_cols=["Rater", "Task"],
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )
    pair = build_context_margin_endpoint_pair(
        original,
        perturbation,
        facet_cols=["Rater", "Task"],
        context_cols=["Task"],
    )

    assert pair["available"] is True
    assert list(pair["designs"]) == ["observed_assignment", "aligned_counterfactual"]
    assert pair["dose_table"]["AchievedAlignmentDose"].tolist() == [0.0, 1.0]
    assert pair["dose_table"]["SwitchesApplied"].eq(0).all()
    assert pair["dose_table"].iloc[-1]["OptimizationEndpoint"]
    assert all("Score" not in design.columns for design in pair["designs"].values())
    required = pair["path_invariants"].loc[pair["path_invariants"]["RequiredAtDose"]]
    assert required["PassedForDose"].all()


def test_frequency_weights_and_disconnected_overlap_fail_closed():
    weighted = _unequal_context_design()
    weighted.loc[0, "Weight"] = 2.0
    weighted_audit = evaluate_context_margin_milp_feasibility(
        weighted,
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )
    assert weighted_audit["available"] is False
    assert "Uncompressed row topology" in weighted_audit["reason"]

    disconnected = pd.DataFrame(
        [
            {"Person": "p0", "Rater": "r0", "Task": "t1"},
            {"Person": "p1", "Rater": "r0", "Task": "t1"},
            {"Person": "p2", "Rater": "r1", "Task": "t1"},
            {"Person": "p3", "Rater": "r1", "Task": "t1"},
        ]
    )
    disconnected_audit = evaluate_context_margin_milp_feasibility(
        disconnected,
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )
    assert disconnected_audit["available"] is False
    assert "Connected direct Rater overlap" in disconnected_audit["reason"]


def test_structurally_fixed_optimum_is_not_exposed_as_a_sensitivity_contrast():
    complete = pd.DataFrame(
        [
            {"Person": person, "Rater": rater, "Task": "t1", "Score": 1}
            for person in ("p0", "p1", "p2")
            for rater in ("r0", "r1", "r2")
        ]
    )
    bundle = build_context_margin_assignment_perturbation(
        complete,
        person_scores={"p0": 0.0, "p1": 1.0, "p2": 2.0},
        rater_scores={"r0": 0.0, "r1": 1.0, "r2": 2.0},
        context_cols=["Task"],
        rater_mapping_confirmed=True,
    )

    assert bundle["available"] is False
    failed = bundle["invariants"].loc[~bundle["invariants"]["Passed"], "Invariant"].tolist()
    assert "Assignment changed" in failed
    assert "Direction-adjusted objective improved" in failed
