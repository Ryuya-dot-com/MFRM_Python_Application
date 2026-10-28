import json
from pathlib import Path
import numpy as np
import pytest

from validation import robustness_writing_check as study

PLAN = json.loads((Path(study.ROOT)/"validation/robustness_writing_protocol.json").read_text())


@pytest.mark.parametrize("design", PLAN["designs"])
def test_same_mixture_pool_and_scores_with_only_assignment_changed(design):
    independent, a = study.generate(design,"mixture_independent",0,PLAN["development_seed_entropy"])
    selected, b = study.generate(design,"mixture_selected",0,PLAN["development_seed_entropy"])
    for field in ("person_effect","z","task_shock","group","permutation"):
        np.testing.assert_array_equal(a[field],b[field])
    np.testing.assert_array_equal(independent.Score,selected.Score)
    assert not independent.Rater.equals(selected.Rater)
    expected = np.full(10,20) if design=="rotating" else np.r_[np.full(5,38),np.full(5,2)]
    np.testing.assert_array_equal(b["positive_group_persons_per_rater"],expected)
    for data in (independent,selected):
        assert len(data)==1200
        assert data.groupby(["Rater","Task"]).size().eq(40).all()
        assert data.groupby(["Person","Task"]).size().eq(2).all()


def test_shared_task_condition_reuses_draws_and_keeps_task_effect_shared_by_raters():
    baseline,a = study.generate("weak_bridge","normal_independent",0,PLAN["development_seed_entropy"])
    shared,b = study.generate("weak_bridge","shared_task_total_variance",0,PLAN["development_seed_entropy"])
    np.testing.assert_array_equal(a["z"],b["z"])
    np.testing.assert_array_equal(a["task_shock"],b["task_shock"])
    assert baseline.Rater.equals(shared.Rater)
    np.testing.assert_allclose(b["person_effect"],.8*a["z"])
    np.testing.assert_allclose(b["person_task_effect"],.6*a["task_shock"])
    np.testing.assert_allclose(b["generating_latent_predictor"],.8*a["z"][:,None]+.6*a["task_shock"]-[-.3,0.,.3])
    assert np.count_nonzero(b["group"]==1)==100
    assert abs(np.mean(a["z"]))>1e-6
    _,other = study.generate("weak_bridge","normal_independent",1,PLAN["development_seed_entropy"])
    assert not np.array_equal(a["z"],other["z"])


def test_task_contrast_uses_global_task_coordinates_and_covariance():
    point=np.zeros(14);point[9],point[11]=-.3,.3
    covariance=np.eye(14);covariance[9,11]=covariance[11,9]=.5
    value=study.task_contrast(point,covariance)
    assert value["estimate"]==pytest.approx(.6)
    assert value["full"]["se"]==pytest.approx(1.)
    assert value["ignored_covariance"]["se"]==pytest.approx(np.sqrt(2))


def test_analysis_preserves_unresolved_rows_and_paired_replication_unit(tmp_path):
    plan=dict(PLAN,replicates_per_cell=2)
    contrasts={target:dict(estimate=spec["truth"],truth=spec["truth"],
        full=dict(se=.1,covered=True,excludes_zero=bool(spec["truth"])),
        ignored_covariance=dict(se=.2,covered=True,excludes_zero=bool(spec["truth"])))
        for target,spec in plan["targets"].items()}
    for mechanism in ("mixture_independent","mixture_selected"):
        path=tmp_path/"cases"/study.case_id("weak_bridge",mechanism)/"rep_000/q81_grouped.json"
        study.base.save(path,dict(status="finite_rule_stationary",numerical_pass=True,contrasts=contrasts))
    result=study.analyze(tmp_path,plan)
    row=next(r for r in result["cells"] if r["design"]=="weak_bridge" and r["mechanism"]=="mixture_selected" and r["target"]=="primary")
    assert (row["planned"],row["eligible"])==(2,1)
    assert row["full"]["generating_truth_containment"]["planned_lower"]==.5
    pair=next(r for r in result["paired_selected_minus_independent"] if r["design"]=="weak_bridge" and r["target"]=="primary")
    assert pair["both_eligible"]==1 and len(pair["rows"])==2
    assert pair["rows"][1]["estimate_difference"] is None
    task=next(r for r in result["cells"] if r["target"]=="task")
    assert "generating_null_zero_exclusion" not in task["full"]
