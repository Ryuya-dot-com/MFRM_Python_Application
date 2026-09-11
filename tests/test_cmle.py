"""Exact CMLE core, eligibility, and estimator-contract tests."""

from __future__ import annotations

from itertools import product
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy.special import logsumexp

from mfrm_app.cmle import (
    CMLEEligibilityError,
    audit_cmle_eligibility,
    cmle_objective_value_grad,
    cmle_objective_value_grad_hessian,
    fit_cmle,
    prepare_cmle_design,
)


def _small_complete_frame() -> pd.DataFrame:
    scores = {
        "P1": [0, 1, 1, 2],
        "P2": [1, 0, 2, 1],
        "P3": [2, 1, 1, 0],
        "P4": [1, 2, 0, 1],
        "P5": [0, 2, 1, 1],
        "P6": [2, 0, 1, 1],
    }
    units = [("R1", "C1"), ("R1", "C2"), ("R2", "C1"), ("R2", "C2")]
    rows = []
    for person, values in scores.items():
        rows.extend(
            (person, rater, criterion, value)
            for (rater, criterion), value in zip(units, values)
        )
    return pd.DataFrame(rows, columns=["Person", "Rater", "Criterion", "Score"])


def _binary_separation_frame(persons: int = 20) -> pd.DataFrame:
    rows = []
    for person_index in range(persons):
        rows.extend(
            [
                (f"P{person_index:03d}", "R1", 1),
                (f"P{person_index:03d}", "R2", 0),
            ]
        )
    return pd.DataFrame(rows, columns=["Person", "Rater", "Score"])


def _manual_conditional_nll(parameters: np.ndarray, design) -> float:
    log_likelihood = 0.0
    max_category = design.n_categories - 1
    for pattern in design.patterns:
        log_kernel = np.einsum("ukp,p->uk", pattern.design, parameters)
        log_likelihood += float(pattern.observed_statistics @ parameters)
        for score in np.flatnonzero(pattern.score_frequencies > 0):
            values = []
            for response in product(range(design.n_categories), repeat=len(pattern.signature)):
                if sum(response) == int(score):
                    values.append(sum(log_kernel[u, category] for u, category in enumerate(response)))
            log_likelihood -= float(pattern.score_frequencies[score]) * float(logsumexp(values))
    return -log_likelihood


def _finite_difference_gradient(function, point: np.ndarray, step: float = 1e-6) -> np.ndarray:
    out = np.zeros_like(point, dtype=float)
    for index in range(len(point)):
        plus = point.copy()
        minus = point.copy()
        plus[index] += step
        minus[index] -= step
        out[index] = (function(plus) - function(minus)) / (2.0 * step)
    return out


def _manual_conditional_information(parameters: np.ndarray, design) -> np.ndarray:
    information = np.zeros((design.n_parameters, design.n_parameters), dtype=float)
    for pattern in design.patterns:
        log_kernel = np.einsum("ukp,p->uk", pattern.design, parameters)
        for score in np.flatnonzero(pattern.score_frequencies > 0):
            statistics = []
            log_weights = []
            for response in product(
                range(design.n_categories), repeat=len(pattern.signature)
            ):
                if sum(response) != int(score):
                    continue
                statistics.append(
                    np.sum(
                        pattern.design[np.arange(len(response)), response, :], axis=0
                    )
                )
                log_weights.append(
                    sum(
                        log_kernel[unit, category]
                        for unit, category in enumerate(response)
                    )
                )
            statistic_matrix = np.asarray(statistics, dtype=float)
            probabilities = np.exp(log_weights - logsumexp(log_weights))
            mean = probabilities @ statistic_matrix
            centered = statistic_matrix - mean
            covariance = (centered * probabilities[:, None]).T @ centered
            information += float(pattern.score_frequencies[score]) * covariance
    return information


def test_exact_objective_matches_full_response_enumeration_and_gradient():
    design = prepare_cmle_design(
        _small_complete_frame(),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="RSM",
    )
    parameters = np.array([0.31, -0.22, 0.47])
    value, gradient = cmle_objective_value_grad(parameters, design)
    manual = _manual_conditional_nll(parameters, design)
    numeric = _finite_difference_gradient(
        lambda value: cmle_objective_value_grad(value, design)[0], parameters
    )
    assert value == pytest.approx(manual, abs=1e-12)
    assert np.max(np.abs(gradient - numeric)) < 2e-8


@pytest.mark.parametrize(
    ("model", "step_facet", "parameters"),
    [
        ("RSM", None, np.array([0.31, -0.22, 0.47])),
        ("PCM", "Criterion", np.array([0.31, -0.22, 0.47, -0.18])),
    ],
)
def test_analytical_information_matches_enumeration_and_gradient_jacobian(
    model, step_facet, parameters
):
    design = prepare_cmle_design(
        _small_complete_frame(),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model=model,
        step_facet=step_facet,
    )
    value, gradient, information = cmle_objective_value_grad_hessian(
        parameters, design
    )
    plain_value, plain_gradient = cmle_objective_value_grad(parameters, design)
    enumerated = _manual_conditional_information(parameters, design)
    finite_difference = np.column_stack(
        [
            _finite_difference_gradient(
                lambda point: cmle_objective_value_grad(point, design)[1][column],
                parameters,
                step=1e-5,
            )
            for column in range(design.n_parameters)
        ]
    ).T

    assert value == pytest.approx(plain_value, abs=1e-12)
    assert np.max(np.abs(gradient - plain_gradient)) < 1e-12
    assert np.max(np.abs(information - enumerated)) < 2e-12
    assert np.max(np.abs(information - finite_difference)) < 2e-9
    assert np.max(np.abs(information - information.T)) < 1e-14
    assert float(np.min(np.linalg.eigvalsh(information))) > 0.0


def test_rsm_fit_returns_ready_structural_fit_and_no_cmle_person_estimates():
    fit = fit_cmle(
        _small_complete_frame(),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="RSM",
        gtol=1e-8,
    )
    summary = fit["summary"].iloc[0]
    assert bool(summary["Converged"])
    assert bool(summary["InferenceReady"])
    assert bool(summary["FiniteMLEGateEnabled"])
    assert summary["FiniteMLEStatus"] == "interior_finite_cmle_supported"
    assert bool(summary["FiniteMLEExistenceQualified"])
    assert fit["finite_mle_audit"] is not None
    assert summary["Method"] == "CMLE"
    assert summary["ConditionalAICScope"] == "compare only matched conditional likelihoods"
    assert fit["facets"]["person"].empty
    assert fit["config"]["person_scoring"] == "not_part_of_cmle_phase0"
    facet_sums = fit["facets"]["others"].groupby("Facet")["Estimate"].sum()
    assert np.max(np.abs(facet_sums.to_numpy())) < 1e-12
    assert abs(float(fit["steps"]["Estimate"].sum())) < 1e-12


def test_default_finite_mle_gate_withholds_optimizer_ready_boundary_fit():
    fit = fit_cmle(
        _binary_separation_frame(),
        person_col="Person",
        facet_cols=["Rater"],
        score_col="Score",
        rating_min=0,
        rating_max=1,
        model="RSM",
        maxiter=800,
        gtol=1e-8,
    )
    summary = fit["summary"].iloc[0]
    assert bool(summary["Converged"])
    assert int(summary["InformationNullity"]) == 0
    assert summary["FiniteMLEStatus"] == "boundary_no_finite_cmle"
    assert bool(summary["FiniteMLEBoundaryDetected"])
    assert not bool(summary["FiniteMLEExistenceQualified"])
    assert not bool(summary["InferenceReady"])
    assert (
        "finite_mle_boundary_no_finite_cmle"
        in str(summary["ReadinessReasons"]).split(";")
    )


def test_explicit_gate_disabled_path_reproduces_technical_boundary_readiness():
    fit = fit_cmle(
        _binary_separation_frame(),
        person_col="Person",
        facet_cols=["Rater"],
        score_col="Score",
        rating_min=0,
        rating_max=1,
        model="RSM",
        maxiter=800,
        gtol=1e-8,
        finite_mle_gate=False,
    )
    summary = fit["summary"].iloc[0]
    assert not bool(summary["FiniteMLEGateEnabled"])
    assert summary["FiniteMLEStatus"] == "not_evaluated"
    assert bool(summary["InferenceReady"])
    assert fit["finite_mle_audit"] is None
    assert fit["config"]["finite_mle_gate"] is False


def test_finite_mle_gate_does_not_move_interior_numeric_outputs():
    kwargs = dict(
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="RSM",
        gtol=1e-8,
    )
    gated = fit_cmle(_small_complete_frame(), **kwargs)
    comparison = fit_cmle(
        _small_complete_frame(), finite_mle_gate=False, **kwargs
    )
    assert bool(gated["summary"].iloc[0]["InferenceReady"])
    assert bool(comparison["summary"].iloc[0]["InferenceReady"])
    assert np.array_equal(
        gated["coefficients"]["Estimate"].to_numpy(),
        comparison["coefficients"]["Estimate"].to_numpy(),
    )
    assert np.array_equal(
        gated["coefficients"]["SE"].to_numpy(),
        comparison["coefficients"]["SE"].to_numpy(),
    )
    assert (
        gated["summary"].iloc[0]["ConditionalLogLik"]
        == comparison["summary"].iloc[0]["ConditionalLogLik"]
    )


def test_finite_mle_unavailable_result_fails_readiness_closed(monkeypatch):
    import mfrm_app.cmle_existence as existence

    def unavailable(_design):
        return {
            "summary": pd.DataFrame(
                [
                    {
                        "Status": "tolerance_unstable",
                        "Reason": "injected_test_instability",
                        "BoundaryDetected": pd.NA,
                        "ExistenceQualified": False,
                        "ToleranceGrid": "1e-10;1e-9;1e-8",
                    }
                ]
            )
        }

    monkeypatch.setattr(existence, "audit_cmle_finite_mle_oracle", unavailable)
    fit = fit_cmle(
        _small_complete_frame(),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="RSM",
        gtol=1e-8,
    )
    summary = fit["summary"].iloc[0]
    assert bool(summary["Converged"])
    assert summary["FiniteMLEStatus"] == "tolerance_unstable"
    assert not bool(summary["InferenceReady"])
    assert "finite_mle_tolerance_unstable" in summary["ReadinessReasons"]


def test_finite_mle_gate_requires_boolean():
    with pytest.raises(ValueError, match="finite_mle_gate"):
        fit_cmle(
            _small_complete_frame(),
            person_col="Person",
            facet_cols=["Rater", "Criterion"],
            score_col="Score",
            rating_min=0,
            rating_max=2,
            finite_mle_gate="yes",
        )


def test_exact_newton_polish_can_finish_an_iteration_limited_bfgs_fit():
    fit = fit_cmle(
        _small_complete_frame(),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="RSM",
        maxiter=0,
        gtol=1e-8,
        newton_polish_maxiter=12,
    )
    summary = fit["summary"].iloc[0]
    assert bool(summary["OptimizerSuccess"]) is False
    assert int(summary["NewtonPolishAcceptedSteps"]) >= 1
    assert summary["NewtonPolishReason"] == "stationarity_reached"
    assert bool(summary["Converged"])
    assert bool(summary["InferenceReady"])
    assert float(summary["GradientSupNorm"]) <= float(
        summary["StationarityTolerance"]
    )
    assert len(fit["newton_polish_trace"]) == int(summary["NewtonPolishSteps"])
    assert fit["objective_trace"][-1] == pytest.approx(
        float(fit["newton_polish_trace"].iloc[-1]["ObjectiveAfter"])
    )


def _simulate_pcm(seed: int = 1042, n_persons: int = 80) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rater_effect = {"R1": -0.35, "R2": 0.35}
    criterion_effect = {"C1": -0.2, "C2": 0.2}
    thresholds = {
        "C1": np.array([-1.0, 0.1, 0.9]),
        "C2": np.array([-0.7, -0.1, 0.8]),
    }
    rows = []
    for person_index in range(n_persons):
        theta = rng.normal()
        for rater, rater_value in rater_effect.items():
            for criterion, criterion_value in criterion_effect.items():
                cumulative = np.r_[0.0, np.cumsum(thresholds[criterion])]
                categories = np.arange(4, dtype=float)
                logits = categories * (theta - rater_value - criterion_value) - cumulative
                probabilities = np.exp(logits - logsumexp(logits))
                score = int(rng.choice(4, p=probabilities))
                rows.append((f"P{person_index:03d}", rater, criterion, score))
    return pd.DataFrame(rows, columns=["Person", "Rater", "Criterion", "Score"])


def test_pcm_fit_preserves_per_step_facet_threshold_constraints():
    fit = fit_cmle(
        _simulate_pcm(),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=3,
        model="PCM",
        step_facet="Criterion",
        maxiter=800,
    )
    assert bool(fit["summary"].iloc[0]["InferenceReady"])
    threshold_sums = fit["steps"].groupby("Level")["Estimate"].sum()
    assert np.max(np.abs(threshold_sums.to_numpy())) < 1e-12
    assert set(fit["steps"]["Level"]) == {"C1", "C2"}


def test_conditional_rank_blocks_facet_constant_within_every_person():
    rows = []
    scores = ([0, 1, 2], [1, 2, 0], [2, 0, 1], [0, 2, 1])
    for person_index in range(12):
        rater = "R1" if person_index < 6 else "R2"
        for criterion, score in zip(["C1", "C2", "C3"], scores[person_index % 4]):
            rows.append((f"P{person_index}", rater, criterion, score))
    frame = pd.DataFrame(rows, columns=["Person", "Rater", "Criterion", "Score"])
    audit = audit_cmle_eligibility(
        frame,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="RSM",
    )
    assert audit["eligible"] is False
    assert audit["conditional_nullity"] >= 1
    assert "conditional_information_rank_deficient" in set(audit["issues_table"]["Code"])
    with pytest.raises(CMLEEligibilityError, match="Conditional information rank"):
        fit_cmle(
            frame,
            person_col="Person",
            facet_cols=["Rater", "Criterion"],
            score_col="Score",
            rating_min=0,
            rating_max=2,
        )


def test_primary_optimizer_rejects_nonfinite_trials_and_returns_fail_closed():
    response_blocks = {
        "P057": {"R01": [2, 1, 1, 0], "R02": [1, 1, 0, 0]},
        "P114": {"R03": [3, 3, 3, 3], "R04": [2, 2, 2, 3]},
        "P178": {"R02": [3, 3, 3, 3], "R03": [3, 1, 2, 3]},
        "P269": {"R05": [1, 2, 2, 0], "R06": [2, 2, 0, 0]},
        "P285": {"R04": [2, 0, 0, 0], "R05": [0, 0, 0, 0]},
        "P001": {"R03": [2, 3, 2, 2]},
        "P002": {"R01": [2, 2, 3, 1]},
        "P003": {"R05": [3, 1, 1, 0]},
        "P004": {"R04": [1, 0, 0, 0]},
        "P005": {"R05": [3, 2, 0, 1]},
        "P006": {"R05": [1, 2, 0, 0]},
    }
    rows = []
    criteria = ["C01", "C02", "C03", "C04"]
    for person, rater_blocks in response_blocks.items():
        for rater, scores in rater_blocks.items():
            rows.extend(
                (person, rater, criterion, score)
                for criterion, score in zip(criteria, scores, strict=True)
            )
    frame = pd.DataFrame(
        rows, columns=["Person", "Rater", "Criterion", "Score"]
    )
    fit = fit_cmle(
        frame,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=3,
        gtol=1e-8,
        maxiter=800,
        newton_polish_maxiter=12,
    )
    summary = fit["summary"].iloc[0]
    assert int(summary["OptimizerInvalidEvaluations"]) > 0
    assert int(summary["OptimizerFiniteEvaluations"]) > 0
    assert not bool(summary["InferenceReady"])
    assert int(summary["InformationNullity"]) > 0

    design = fit["design"]
    with np.errstate(over="ignore", invalid="ignore"):
        with pytest.raises(FloatingPointError):
            cmle_objective_value_grad(
                np.full(design.n_parameters, np.finfo(float).max), design
            )


def test_primary_optimizer_falls_back_to_best_finite_point(monkeypatch):
    def terminate_at_invalid_point(
        fun, x0, args=(), method=None, jac=None, callback=None, options=None
    ):
        start = np.asarray(x0, dtype=float)
        start_value, _ = fun(start, *args)
        assert np.isfinite(start_value)
        invalid = np.full_like(start, np.finfo(float).max)
        with np.errstate(over="ignore", invalid="ignore"):
            invalid_value, invalid_gradient = fun(invalid, *args)
        assert np.isinf(invalid_value)
        assert np.all(invalid_gradient == 0.0)
        return SimpleNamespace(
            x=invalid,
            success=True,
            nit=1,
            nfev=2,
            message="deliberate non-finite terminal trial",
        )

    monkeypatch.setattr("mfrm_app.cmle.minimize", terminate_at_invalid_point)
    fit = fit_cmle(
        _small_complete_frame(),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        gtol=1e-8,
        maxiter=20,
        newton_polish_maxiter=12,
    )
    summary = fit["summary"].iloc[0]
    assert int(summary["OptimizerInvalidEvaluations"]) == 1
    assert bool(summary["OptimizerBestFiniteFallbackUsed"])
    assert not bool(summary["OptimizerSuccess"])
    assert "finite" in str(summary["OptimizerLastInvalidReason"])


def test_duplicate_response_units_fail_closed_without_event_identity():
    frame = pd.concat([_small_complete_frame(), _small_complete_frame().iloc[[0]]], ignore_index=True)
    audit = audit_cmle_eligibility(
        frame,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
    )
    assert audit["eligible"] is False
    assert "duplicate_response_units" in set(audit["issues_table"]["Code"])


def test_nonunit_row_weights_and_gpcm_fail_closed():
    frame = _small_complete_frame()
    frame["Weight"] = 1.0
    frame.loc[0, "Weight"] = 0.5
    weighted = audit_cmle_eligibility(
        frame,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        weight_col="Weight",
    )
    assert weighted["eligible"] is False
    assert "nonunit_row_weights_unsupported" in set(weighted["issues_table"]["Code"])

    gpcm = audit_cmle_eligibility(
        frame.drop(columns="Weight"),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="GPCM",
    )
    assert gpcm["eligible"] is False
    assert "GPCM" in gpcm["issues_table"].iloc[0]["Message"]


def test_phase0_rank_work_caps_fail_closed_before_optimization():
    parameter_audit = audit_cmle_eligibility(
        _small_complete_frame(),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        rank_audit_max_parameters=2,
    )
    assert parameter_audit["eligible"] is False
    assert "phase0_rank_audit_parameter_cap" in set(
        parameter_audit["issues_table"]["Code"]
    )

    work_audit = audit_cmle_eligibility(
        _small_complete_frame(),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        rank_audit_max_work=1,
    )
    assert work_audit["eligible"] is False
    assert work_audit["rank_audit_method"] == "exact_conditional_second_moment_dp"
    assert "phase0_rank_audit_work_cap" in set(work_audit["issues_table"]["Code"])

    memory_audit = audit_cmle_eligibility(
        _small_complete_frame(),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        rank_audit_max_bytes=1,
    )
    assert memory_audit["eligible"] is False
    assert memory_audit["rank_audit_peak_bytes_proxy"] > 1
    assert "phase0_rank_audit_memory_cap" in set(
        memory_audit["issues_table"]["Code"]
    )


def test_declared_category_support_is_not_truncated_by_virtual_unit_maximum():
    frame = _simulate_pcm(n_persons=50)
    mask = (frame["Rater"] == "R2") & (frame["Criterion"] == "C2") & (frame["Score"] == 3)
    frame.loc[mask, "Score"] = 2
    # Keep the global top category present in other virtual units.
    assert 3 in set(frame["Score"])
    fit = fit_cmle(
        frame,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=3,
        model="RSM",
    )
    target = fit["surfaces"].query("Rater == 'R2' and Criterion == 'C2'")
    assert set(target["Category"]) == {0, 1, 2, 3}


def _immer_parity_frame(seed: int = 9981, n_persons: int = 70) -> pd.DataFrame:
    frame = _simulate_pcm(seed=seed, n_persons=n_persons)
    # Force every virtual item to expose all four categories so immer 1.5-13's
    # observed per-column maxK equals the declared Python support.
    for person_index, category in enumerate(range(4)):
        frame.loc[frame["Person"] == f"P{person_index:03d}", "Score"] = category
    return frame


def _r_has_immer() -> bool:
    if shutil.which("Rscript") is None:
        return False
    probe = subprocess.run(
        ["Rscript", "-e", "quit(status=ifelse(requireNamespace('immer', quietly=TRUE),0,1))"],
        capture_output=True,
        text=True,
        check=False,
    )
    return probe.returncode == 0


@pytest.mark.skipif(not _r_has_immer(), reason="R package immer is unavailable")
@pytest.mark.parametrize(
    ("model", "planned_missing"),
    [("RSM", False), ("PCM", True)],
)
def test_free_coordinates_and_conditional_loglik_match_immer(
    tmp_path, model, planned_missing
):
    frame = _immer_parity_frame()
    if planned_missing:
        person_number = frame["Person"].str.removeprefix("P").astype(int)
        unit_number = frame["Rater"].map({"R1": 0, "R2": 2}) + frame["Criterion"].map(
            {"C1": 0, "C2": 1}
        )
        frame = frame.loc[~((person_number >= 4) & (unit_number == person_number % 4))].copy()
    fit = fit_cmle(
        frame,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=3,
        model=model,
        step_facet="Criterion" if model == "PCM" else None,
        maxiter=1000,
        gtol=1e-9,
    )
    frame = frame.copy()
    frame["VirtualUnit"] = frame["Rater"] + "__" + frame["Criterion"]
    wide = frame.pivot(index="Person", columns="VirtualUnit", values="Score").reset_index()
    mapping = frame[["VirtualUnit", "Rater", "Criterion"]].drop_duplicates().sort_values("VirtualUnit")
    wide_path = tmp_path / "wide.csv"
    mapping_path = tmp_path / "mapping.csv"
    output_path = tmp_path / "immer.csv"
    wide.to_csv(wide_path, index=False)
    mapping.to_csv(mapping_path, index=False)
    runner = Path(__file__).parent / "data" / "run_immer_cmle_parity.R"
    completed = subprocess.run(
        [
            "Rscript",
            str(runner),
            str(wide_path),
            str(mapping_path),
            str(output_path),
            model,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    r_result = pd.read_csv(output_path)
    py_result = fit["coefficients"][["Parameter", "Estimate"]]
    joined = py_result.merge(r_result[["Parameter", "Estimate"]], on="Parameter", suffixes=("_Python", "_immer"))
    assert len(joined) == len(py_result) == len(r_result)
    assert np.max(np.abs(joined["Estimate_Python"] - joined["Estimate_immer"])) < 2e-5
    py_loglik = float(fit["summary"].iloc[0]["ConditionalLogLik"])
    r_loglik = float(r_result["ConditionalLogLik"].iloc[0])
    assert py_loglik == pytest.approx(r_loglik, abs=2e-7)
