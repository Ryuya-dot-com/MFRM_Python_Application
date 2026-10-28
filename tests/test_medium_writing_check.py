"""Deterministic checks for the separately frozen exploratory study driver."""
import json
from pathlib import Path

import numpy as np
import pytest

from validation import medium_writing_check as study


PLAN = json.loads((Path(study.ROOT) / "validation/medium_writing_protocol.json").read_text())


@pytest.mark.parametrize("n", [200, 500])
def test_assignments_change_overlap_without_changing_workload(n):
    workloads = []
    for design, expected_bridge in [("rotating", n), ("weak_bridge", n // 10)]:
        pairs = study.assignment(n, design)
        assert np.all(pairs[:, 0] != pairs[:, 1])
        assert np.sum(pairs[:, 0] // 5 != pairs[:, 1] // 5) == expected_bridge
        data, truth = study.generate(n, design, .4, 0, PLAN["development_seed_entropy"])
        assert len(data) == 6 * n
        assert data.groupby("Person").size().eq(6).all()
        assert data.groupby(["Person", "Task"]).size().eq(2).all()
        assert set(data.Score) == set(range(4))
        workloads.append(data.groupby(["Rater", "Task"]).size().to_numpy())
        assert abs(np.mean(truth["theta"])) > 1e-6  # No realized-sample recentering.
    np.testing.assert_array_equal(workloads[0], workloads[1])
    np.testing.assert_array_equal(workloads[0], np.full(30, 2 * n // 10))


def test_generation_is_reproducible_and_development_seed_is_separate():
    args = (200, "weak_bridge", .4, 0)
    data, truth = study.generate(*args, PLAN["development_seed_entropy"])
    repeated, repeated_truth = study.generate(*args, PLAN["development_seed_entropy"])
    assert data.equals(repeated)
    np.testing.assert_array_equal(truth["theta"], repeated_truth["theta"])
    different, different_truth = study.generate(*args, PLAN["development_seed_entropy"] + 1)
    assert not data.Score.equals(different.Score)
    assert not np.array_equal(truth["theta"], different_truth["theta"])


def test_contrasts_expand_the_constrained_rater_and_keep_covariance():
    point = np.r_[np.full(5, .2), np.full(4, -.2), np.zeros(5)]
    covariance = np.eye(14)
    result = study.contrast_results(point, covariance, {"primary": .4, "secondary": 0.})
    assert result["primary"]["estimate"] == pytest.approx(.4)
    assert result["secondary"]["estimate"] == 0
    assert result["primary"]["full"]["se"] ** 2 == pytest.approx(.8)
    assert result["primary"]["ignored_covariance"]["se"] ** 2 == pytest.approx(.72)
    assert result["secondary"]["full"]["se"] ** 2 == pytest.approx(2.)


def test_unresolved_results_remain_in_planned_denominator():
    summary = study.binary_summary(8, 10, 24)
    assert summary["proportion"] == .8
    assert summary["planned_lower"] == pytest.approx(8 / 24)
    assert summary["planned_upper"] == pytest.approx(22 / 24)
    assert summary["mcse"] == pytest.approx(np.sqrt(.8 * .2 / 10))
    empty = study.binary_summary(0, 0, 24)
    assert empty["proportion"] is None
    assert (empty["planned_lower"], empty["planned_upper"]) == (0., 1.)
    complete = study.binary_summary(24, 24, 24)
    assert complete["wilson_lower"] < 1 and complete["wilson_upper"] == pytest.approx(1.)


def test_native_grouped_parity_at_prespecified_start_and_perturbed_point():
    data, _ = study.generate(200, "weak_bridge", .4, 0, PLAN["development_seed_entropy"])
    prep, problem = study.prepare(data, 31)
    grouped = study.GroupedMmlLikelihood(prep["data"], problem.config)
    nodes, weights = study.normal_rule(31)
    for point in [problem.structural_start, problem.structural_start + np.linspace(-.03, .03, 14)]:
        evaluated = grouped.evaluate(np.r_[point, 0.], nodes, weights, information=True)
        native, gradient = study.app.mfrm_loglik_mml_value_grad(point, problem.idx,
            problem.config, problem.sizes, dict(nodes=nodes, weights=weights, sd=1.))
        assert native == pytest.approx(evaluated["nll"], abs=1e-9)
        np.testing.assert_allclose(gradient, evaluated["gradient"][:-1], atol=1e-10, rtol=1e-10)
        # Fixed-SD covariance must use the structural block, before inversion.
        assert evaluated["information"][:-1, :-1].shape == (14, 14)


def test_interrupted_and_timed_out_attempts_are_retained_without_retry(tmp_path):
    path = tmp_path / "interrupted.json"
    study.save(path.with_suffix(".started.json"), {"started_unix": 0})
    def must_not_run():
        raise AssertionError("An interrupted attempt must not restart")
    result = study.attempt(path, 1., must_not_run)
    assert result["status"] == "interrupted_attempt_not_retried"
    assert json.loads(path.read_text()) == result
    def timeout():
        raise study.StudyTimeout("test timeout")
    result = study.attempt(tmp_path / "timeout.json", 1., timeout)
    assert result["status"] == "timeout"
    assert result["seconds"] >= 0
    assert result["reason_code"] == "stat.medium_writing.timeout"
