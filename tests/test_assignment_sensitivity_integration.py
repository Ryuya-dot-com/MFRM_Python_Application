from __future__ import annotations

import numpy as np
import pandas as pd

import streamlit_app as app


def _exchangeable_sparse_responses() -> pd.DataFrame:
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


def _fitted_jmle() -> dict:
    return app.mfrm_estimate(
        _exchangeable_sparse_responses(),
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=0,
        rating_max=3,
        weight_col="Weight",
        keep_original=True,
        model="RSM",
        method="JMLE",
        maxit=120,
        reltol=1e-5,
    )


def _fitted_mml(*, free_sd: bool = False) -> dict:
    return app.mfrm_estimate(
        _exchangeable_sparse_responses(),
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=0,
        rating_max=3,
        weight_col="Weight",
        keep_original=True,
        model="RSM",
        method="MML",
        quad_points=11,
        population_prior_sd=1.0,
        estimate_population_sd=free_sd,
        maxit=300,
        reltol=1e-5,
    )


def _unequal_context_responses() -> pd.DataFrame:
    base = [
        (("r0", "A"), ("r1", "B")),
        (("r0", "A"), ("r2", "B")),
        (("r1", "A"), ("r0", "B")),
        (("r1", "A"), ("r2", "B")),
        (("r2", "A"), ("r0", "B")),
        (("r2", "A"), ("r1", "B")),
    ]
    rows = []
    for person_index in range(12):
        for rater, signature in base[person_index % len(base)]:
            tasks = ["t1"] if signature == "A" else ["t1", "t2"]
            for task in tasks:
                for criterion_index in range(5):
                    latent = (
                        1.5
                        + 0.22 * (person_index - 5.5)
                        - 0.25 * int(rater[1])
                        + 0.35 * (criterion_index - 2)
                        + (0.15 if task == "t2" else -0.15)
                    )
                    rows.append({
                        "Person": f"p{person_index:02d}",
                        "Rater": rater,
                        "Task": task,
                        "Criterion": f"c{criterion_index}",
                        "Score": int(np.clip(np.rint(latent), 0, 3)),
                        "Weight": 1.0,
                    })
    return pd.DataFrame(rows)


def _unequal_context_jmle() -> dict:
    return app.mfrm_estimate(
        _unequal_context_responses(),
        person_col="Person",
        facet_cols=["Rater", "Task", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=3,
        weight_col="Weight",
        keep_original=True,
        model="RSM",
        method="JMLE",
        maxit=250,
        reltol=1e-5,
    )


def test_real_jmle_fit_runs_paired_fixed_density_sensitivity_end_to_end():
    fitted = _fitted_jmle()
    preflight = app.build_assignment_sensitivity_preflight(fitted)

    assert preflight["available"] is True
    assert preflight["gates"]["Passed"].all()
    bundle = app.simulate_fixed_density_assignment_sensitivity(
        fitted,
        n_replicates=1,
        seed=260811,
        refit_maxit=100,
        refit_reltol=1e-4,
    )

    assert bundle["available"] is True
    assert bundle["invariants"]["Passed"].all()
    assert bundle["summary"]["Method"].eq("JMLE").all()
    assert bundle["summary"]["Scenario"].tolist() == [
        "observed_assignment",
        "aligned_counterfactual",
    ]
    assert bundle["summary"]["Completed"].all()
    assert bundle["contrasts"]["PairComplete"].all()
    assert not bundle["contrasts"]["LikelihoodComparisonPerformed"].any()
    assert not bundle["completion"].iloc[0]["CrossBasisLikelihoodCompared"]
    assert bundle["completion"].iloc[0]["ScenarioFitsExpected"] == 2
    assert bundle["completion"].iloc[0]["ScenarioFitsCompleted"] == 2
    assert bundle["contrast_summary"]["InferenceReadyPairs"].eq(1).all()
    assert not bundle["contrast_summary"]["LikelihoodComparisonPerformed"].any()
    assert "Score" not in bundle["assignment_map"].columns
    assert bundle["trajectory"].iloc[-1]["AssignmentRankCorrelation"] > bundle["trajectory"].iloc[0]["AssignmentRankCorrelation"]
    assert np.isfinite(bundle["summary"]["RaterReferenceCenteredRMSE"]).all()
    assert not bundle["summary"]["ReferenceIsKnownTruth"].any()
    assert not bundle["reference_contract"].iloc[0]["ReferenceKnownTruth"]

    private_frames = app.assignment_sensitivity_bundle_frames(bundle)
    public_frames = app.prepare_download_frames_for_privacy(
        private_frames,
        public_export_mode=True,
    )
    assert "fixed_density_assignment_refit_summary" in public_frames
    assert "fixed_density_assignment_map" not in public_frames
    assert "fixed_density_assignment_path_map" not in public_frames
    assert "fixed_density_assignment_block_profiles" not in public_frames
    assert "fixed_density_assignment_switch_ledger" not in public_frames
    assert "export_privacy_manifest" in public_frames


def test_unequal_context_fit_routes_to_audited_milp_endpoint_and_refits():
    fitted = _unequal_context_jmle()
    strict = app.build_assignment_sensitivity_preflight(fitted)
    selected = app.select_assignment_sensitivity_preflight(fitted)

    assert strict["available"] is False
    assert "Common context signature" in strict["reason"]
    assert selected["available"] is True
    assert selected["design_engine"] == "context_margin_milp"
    bundle = app.simulate_fixed_density_assignment_sensitivity(
        fitted,
        n_replicates=1,
        seed=260816,
        refit_maxit=100,
        alignment_doses=[0.0, 0.5, 1.0],
        design_engine="context_margin_milp",
    )

    assert bundle["available"] is True
    assert bundle["design_engine"] == "context_margin_milp"
    assert bundle["dose_table"]["AchievedAlignmentDose"].tolist() == [0.0, 1.0]
    assert bundle["completion"].iloc[0]["DesignEngine"] == "context_margin_milp"
    assert bundle["completion"].iloc[0]["DoseContrastsComplete"] == 1
    assert bundle["solver_audit"].iloc[0]["ObservedBaselineFeasible"]
    assert bundle["solver_audit"].iloc[0]["GlobalOptimumCertified"]
    assert bundle["context_margin_audit"]["ExactMatch"].all()
    assert bundle["invariants"]["Passed"].all()
    assert len(bundle["witnesses"]) == 2
    assert not bundle["reference_contract"].iloc[0]["IntermediateDoseQualified"]
    assert not bundle["contrasts"]["LikelihoodComparisonPerformed"].any()

    private_frames = app.assignment_sensitivity_bundle_frames(bundle)
    public_frames = app.prepare_download_frames_for_privacy(
        private_frames,
        public_export_mode=True,
    )
    assert "assignment_context_margin_audit" in public_frames
    assert "assignment_context_milp_solver_audit" in public_frames
    assert "assignment_context_connectivity_witnesses" not in public_frames


def test_unequal_context_fixed_sd_mml_uses_same_milp_endpoint_without_cross_basis_ranking():
    source = _unequal_context_responses()
    fitted = app.mfrm_estimate(
        source,
        person_col="Person",
        facet_cols=["Rater", "Task", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=3,
        weight_col="Weight",
        keep_original=True,
        model="RSM",
        method="MML",
        quad_points=11,
        population_prior_sd=1.0,
        estimate_population_sd=False,
        maxit=300,
        reltol=1e-5,
    )
    selected = app.select_assignment_sensitivity_preflight(fitted)
    assert selected["available"] is True
    assert selected["design_engine"] == "context_margin_milp"

    bundle = app.simulate_fixed_density_assignment_sensitivity(
        fitted,
        n_replicates=1,
        seed=260817,
        refit_maxit=120,
        design_engine="context_margin_milp",
        person_generation_mode="mml_population_rank_preserving",
    )

    assert bundle["available"] is True
    assert bundle["person_generation_mode"] == "mml_population_rank_preserving"
    assert bundle["person_generation_summary"].iloc[0]["RankPreservedExactly"]
    assert bundle["summary"]["Method"].eq("MML").all()
    assert bundle["summary"]["PopulationSDUsed"].eq(1.0).all()
    assert bundle["completion"].iloc[0]["DoseContrastsComplete"] == 1
    assert bundle["solver_audit"].iloc[0]["GlobalOptimumCertified"]
    assert not bundle["contrasts"]["LikelihoodComparisonPerformed"].any()
    assert not bundle["completion"].iloc[0]["CrossBasisLikelihoodCompared"]


def test_preflight_fails_closed_for_estimators_outside_v1_scope():
    fitted = _fitted_jmle()
    fitted["config"]["method"] = "EXACT_CMLE"
    preflight = app.build_assignment_sensitivity_preflight(fitted)

    assert preflight["available"] is False
    gate = preflight["gates"].loc[
        preflight["gates"]["Gate"].eq("Current public estimator basis")
    ].iloc[0]
    assert not gate["Passed"]


def test_jmle_rejects_mml_population_person_generator_before_refits():
    fitted = _fitted_jmle()
    bundle = app.simulate_fixed_density_assignment_sensitivity(
        fitted,
        n_replicates=1,
        person_generation_mode="mml_population_rank_preserving",
    )

    assert bundle["available"] is False
    assert "Generator-estimator compatibility" in bundle["reason"]
    assert bundle["completion"].empty
    assert bundle["summary"].empty
    assert not bundle["person_generator_gates"].loc[
        bundle["person_generator_gates"]["Gate"].eq(
            "Generator-estimator compatibility"
        ),
        "Passed",
    ].iloc[0]


def test_preflight_blocks_degenerate_fitted_alignment_coordinates():
    fitted = _fitted_jmle()
    fitted["facets"]["person"]["Estimate"] = 0.0
    preflight = app.build_assignment_sensitivity_preflight(fitted)

    assert preflight["available"] is False
    gate = preflight["gates"].loc[
        preflight["gates"]["Gate"].eq("Nondegenerate fitted Person ordering")
    ].iloc[0]
    assert not gate["Passed"]


def test_fixed_sd_mml_runner_locks_marginal_method_and_sd_contract():
    fitted = _fitted_mml(free_sd=False)
    bundle = app.simulate_fixed_density_assignment_sensitivity(
        fitted,
        n_replicates=1,
        seed=260812,
        refit_maxit=120,
        refit_reltol=1e-4,
        alignment_doses=[0.0, 0.5, 1.0],
    )

    assert bundle["available"] is True
    assert bundle["summary"]["Method"].eq("MML").all()
    assert bundle["summary"]["PopulationSDUsed"].eq(1.0).all()
    assert not bundle["summary"]["PopulationSDFreeEstimated"].any()
    assert len(bundle["dose_table"]) >= 2
    pop_contrast = bundle["contrast_summary"].loc[
        bundle["contrast_summary"]["Metric"].eq("PopulationSDUsed")
    ].iloc[0]
    assert pop_contrast["MeanDifference"] == 0.0
    assert pop_contrast["FiniteDifferences"] == 1
    pop_curve = bundle["dose_contrast_summary"].loc[
        bundle["dose_contrast_summary"]["Metric"].eq("PopulationSDUsed")
    ]
    assert np.allclose(pop_curve["MeanDifference"], 0.0)
    assert not bundle["contrasts"]["LikelihoodComparisonPerformed"].any()


def test_free_sd_mml_runner_withholds_even_a_legacy_inference_ready_flag(monkeypatch):
    fitted = _fitted_mml(free_sd=True)
    assert not bool(fitted["summary"].iloc[0]["InferenceReady"])
    fitted["summary"]["InferenceReady"] = True  # simulate a pre-audit saved result
    calls = []

    def forbidden_refit(*args, **kwargs):
        calls.append(True)
        raise AssertionError("Unqualified source must not launch sensitivity refits")

    monkeypatch.setattr(app, "_refit_assignment_sensitivity_dataset", forbidden_refit)
    bundle = app.simulate_fixed_density_assignment_sensitivity(
        fitted, n_replicates=1, seed=260813, refit_maxit=100,
        refit_reltol=1e-4, person_generation_mode="mml_population_rank_preserving",
    )
    assert bundle["available"] is False
    assert bundle["summary"].empty
    gate = bundle["gates"].set_index("Gate").loc["Source fit inference readiness"]
    assert not bool(gate["Passed"])
    assert app.FREE_SD_MML_INFERENCE_HOLD_REASON in gate["Evidence"]
    assert not calls


def test_fixed_sd_mml_rank_preserving_population_generator_uses_one_draw_per_pair():
    fitted = _fitted_mml(free_sd=False)
    bundle = app.simulate_fixed_density_assignment_sensitivity(
        fitted,
        n_replicates=2,
        seed=260818,
        refit_maxit=120,
        alignment_doses=[0.0, 0.5, 1.0],
        person_generation_mode="mml_population_rank_preserving",
    )

    assert bundle["available"] is True
    assert bundle["person_generation_mode"] == "mml_population_rank_preserving"
    assert len(bundle["person_generation_summary"]) == 2
    assert bundle["person_generation_summary"]["RankPreservedExactly"].all()
    assert bundle["person_generation_summary"]["GeneratorPopulationSD"].eq(1.0).all()
    assert bundle["summary"]["GeneratedPersonTruthKnownWithinReplicate"].all()
    assert bundle["summary"]["GeneratorPredictionScope"].eq(
        "known_person_conditional_on_supplied_measure"
    ).all()
    for _, replicate_rows in bundle["summary"].groupby("Replicate"):
        assert replicate_rows["GeneratedPersonMean"].nunique() == 1
        assert replicate_rows["GeneratedPersonSD"].nunique() == 1
    assert len(bundle["person_generation_draws"]) == 2 * len(
        fitted["facets"]["person"]
    )
    assert not bundle["reference_contract"].iloc[0]["UnconditionalNewSample"]
    assert bundle["reference_contract"].iloc[0]["RankRelationPreserved"]
    assert not bundle["contrasts"]["LikelihoodComparisonPerformed"].any()

    private_frames = app.assignment_sensitivity_bundle_frames(bundle)
    public_frames = app.prepare_download_frames_for_privacy(
        private_frames, public_export_mode=True
    )
    assert "assignment_generation_summary" in public_frames
    assert "assignment_generation_draws" not in public_frames


def test_prediction_person_override_is_explicit_and_rejects_invalid_people():
    fitted = _fitted_mml(free_sd=False)
    design = fitted["prep"]["data"].drop(columns=["Score"]).head(4).copy()
    persons = design["Person"].astype(str).drop_duplicates().tolist()
    overrides = {person: 0.25 + index for index, person in enumerate(persons)}
    prediction = app.predict_mfrm_design(
        fitted,
        design,
        person_measure_overrides=overrides,
    )

    assert prediction["available"] is True
    table = prediction["table"]
    assert table["PredictionScope"].eq(
        "known_person_conditional_on_supplied_measure"
    ).all()
    assert np.allclose(
        table["ThetaUsed"], table["Person"].astype(str).map(overrides)
    )

    with np.testing.assert_raises(ValueError):
        app.predict_mfrm_design(
            fitted,
            design,
            person_measure_overrides={"not_in_fit": 0.0},
        )


def test_nonready_refits_are_counted_but_never_enter_metric_contrasts(monkeypatch):
    fitted = _fitted_jmle()

    def nonready_refit(*args, **kwargs):
        return {
            "summary": pd.DataFrame([{
                "Converged": False,
                "InferenceReady": False,
            }]),
            "facets": {},
            "config": {},
        }

    monkeypatch.setattr(app, "_refit_assignment_sensitivity_dataset", nonready_refit)
    bundle = app.simulate_fixed_density_assignment_sensitivity(
        fitted,
        n_replicates=1,
        seed=260814,
    )

    completion = bundle["completion"].iloc[0]
    assert bundle["available"] is True
    assert bundle["evidence_ready"] is False
    assert completion["ScenarioFitsCompleted"] == 2
    assert completion["ScenarioFitsInferenceReady"] == 0
    assert completion["ScenarioFitsMetricsAvailable"] == 0
    assert completion["PairedContrastsComplete"] == 0
    assert not bundle["contrasts"]["PairComplete"].any()
    contrast_columns = [
        column for column in bundle["contrasts"].columns
        if column.startswith("AlignedMinusObserved_")
    ]
    assert bundle["contrasts"][contrast_columns].isna().all().all()


def test_multi_dose_jmle_curve_refits_each_achievable_snapshot_once():
    fitted = _fitted_jmle()
    bundle = app.simulate_fixed_density_assignment_sensitivity(
        fitted,
        n_replicates=1,
        seed=260815,
        alignment_doses=[0.0, 0.25, 0.5, 0.75, 1.0],
    )

    dose_table = bundle["dose_table"]
    scenario_count = len(dose_table)
    assert bundle["available"] is True
    assert 2 <= scenario_count <= 5
    assert dose_table["AchievedAlignmentDose"].is_monotonic_increasing
    assert dose_table.iloc[0]["AchievedAlignmentDose"] == 0.0
    assert dose_table.iloc[-1]["AchievedAlignmentDose"] == 1.0
    assert len(bundle["summary"]) == scenario_count
    assert len(bundle["dose_contrasts"]) == scenario_count - 1
    assert len(bundle["dose_contrast_summary"]) == scenario_count * 5
    completion = bundle["completion"].iloc[0]
    assert completion["AssignmentDoseScenarios"] == scenario_count
    assert completion["ScenarioFitsExpected"] == scenario_count
    assert completion["DoseContrastsExpected"] == scenario_count - 1
    assert completion["DoseContrastsComplete"] == scenario_count - 1
    assert not bundle["dose_contrast_summary"]["LikelihoodComparisonPerformed"].any()
    assert "Score" not in bundle["path_assignment_map"].columns
    required = bundle["path_invariants"].loc[
        bundle["path_invariants"]["RequiredAtDose"]
    ]
    assert required["PassedForDose"].all()
