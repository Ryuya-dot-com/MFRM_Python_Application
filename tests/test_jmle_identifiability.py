"""Contracts for exact-free-coordinate JMLE structural-identifiability guards."""

from __future__ import annotations

import pandas as pd

import streamlit_app as app
from mfrm_app.operating_characteristics import build_replicate_manifest
from validation import operating_characteristics_pilot as pilot


def _generated(design: str, truth_bias: float = 0.0):
    manifest = build_replicate_manifest(
        pilot.study_conditions(),
        replicates=1,
        base_seed=20260809,
    )
    row = manifest.loc[
        manifest["Design"].eq(design)
        & manifest["TruthBias"].eq(float(truth_bias))
    ].iloc[0]
    return row, pilot.generate_condition_data(row)


def _audit(row: pd.Series, generated, *, method: str = "JMLE", max_dense_cells: int = 5_000_000):
    prep = app.prepare_mfrm_data(
        generated.data,
        person_col="Person",
        facet_cols=["Rater", "Task", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=int(row["Categories"]) - 1,
        keep_original=True,
    )
    specs = app.prepare_constraint_specs(
        prep,
        anchor_df=generated.anchors if not generated.anchors.empty else None,
        noncenter_facet="Person",
    )
    config = {
        "model": "RSM",
        "method": method,
        "n_cat": int(row["Categories"]),
        "n_person": len(prep["levels"]["Person"]),
        "facet_names": prep["facet_names"],
        "facet_levels": {facet: prep["levels"][facet] for facet in prep["facet_names"]},
        "facet_signs": {"Rater": -1, "Task": -1, "Criterion": -1},
        "theta_spec": specs["theta_spec"],
        "facet_specs": specs["facet_specs"],
        "population_model": {"enabled": False, "n_params": 0},
    }
    return app.build_jmle_eta_identifiability_audit(
        prep,
        config,
        max_dense_cells=max_dense_cells,
    )


def test_balanced_connected_design_has_full_eta_rank():
    row, generated = _generated("balanced_small")

    audit = _audit(row, generated)
    summary = audit["summary"].iloc[0]
    rater = audit["connectivity"].loc[
        audit["connectivity"]["Facet"].eq("Rater")
    ].iloc[0]

    assert audit["scope"] == "eta_person_facet_blocks_only"
    assert audit["status"] == "eta_identified"
    assert bool(audit["eta_structurally_identified"])
    assert int(summary["EtaStructuralNullity"]) == 0
    assert int(rater["Components"]) == 1
    assert int(rater["MinimumRatersOrLevelsPerPerson"]) == 4


def test_one_rater_per_person_design_has_seven_person_rater_null_directions():
    row, generated = _generated("sparse_missing")

    audit = _audit(row, generated)
    summary = audit["summary"].iloc[0]
    rater = audit["connectivity"].loc[
        audit["connectivity"]["Facet"].eq("Rater")
    ].iloc[0]
    energy = audit["null_space_block_energy"]

    assert audit["status"] == "eta_rank_deficient"
    assert not bool(audit["inference_ready"])
    assert int(summary["EtaStructuralNullity"]) == 7
    assert int(rater["Components"]) == 8
    assert int(rater["MinimumRatersOrLevelsPerPerson"]) == 1
    assert set(energy["Block"]) == {"theta", "Rater", "Task", "Criterion"}
    mean_energy = energy.groupby("Block")["SquaredEnergyShare"].mean()
    assert mean_energy["theta"] > 0.65
    assert mean_energy["Rater"] > 0.08
    assert mean_energy["Task"] < 1e-20
    assert mean_energy["Criterion"] < 1e-20


def test_mml_is_explicitly_outside_jmle_eta_rank_gate():
    row, generated = _generated("balanced_small")

    audit = _audit(row, generated, method="MML")

    assert audit["status"] == "not_applicable_mml"
    assert audit["eta_structurally_identified"] is None
    assert audit["summary"].iloc[0]["Scope"] == "eta_person_facet_blocks_only"


def test_large_design_limit_fails_closed_instead_of_claiming_identification():
    row, generated = _generated("balanced_small")

    audit = _audit(row, generated, max_dense_cells=1)

    assert audit["status"] == "audit_size_limit"
    assert not bool(audit["inference_ready"])
    assert not bool(audit["summary"].iloc[0]["InferenceReady"])


def test_wide_design_uses_rank_gate_without_claiming_null_space_localization():
    frame = pd.DataFrame({
        "Person": ["P1", "P1"],
        "Rater": ["R1", "R1"],
        "Task": ["T1", "T1"],
        "Criterion": ["C1", "C1"],
        "Score": [0, 1],
    })
    prep = app.prepare_mfrm_data(
        frame,
        person_col="Person",
        facet_cols=["Rater", "Task", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=1,
        keep_original=True,
    )
    prep["levels"]["Person"] = ["P1", "P2", "P3"]
    specs = app.prepare_constraint_specs(prep, noncenter_facet="Person")
    config = {
        "model": "RSM",
        "method": "JMLE",
        "n_cat": 2,
        "n_person": 3,
        "facet_names": prep["facet_names"],
        "facet_levels": {facet: prep["levels"][facet] for facet in prep["facet_names"]},
        "facet_signs": {"Rater": -1, "Task": -1, "Criterion": -1},
        "theta_spec": specs["theta_spec"],
        "facet_specs": specs["facet_specs"],
        "population_model": {"enabled": False, "n_params": 0},
    }

    audit = app.build_jmle_eta_identifiability_audit(prep, config)
    summary = audit["summary"].iloc[0]

    assert audit["status"] == "eta_rank_deficient"
    assert int(summary["EtaStructuralNullity"]) == 2
    assert not bool(summary["NullSpaceBasisComplete"])
    assert audit["coordinate_null_weight"]["NullProjectionWeight"].isna().all()
    assert audit["coordinate_null_weight"]["NullVulnerable"].isna().all()
    assert audit["null_space_block_energy"].empty


def test_rank_deficient_fit_is_not_inference_ready_and_bias_is_withheld():
    row, generated = _generated("sparse_missing")
    result = app.mfrm_estimate(
        generated.data,
        person_col="Person",
        facet_cols=["Rater", "Task", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=int(row["Categories"]) - 1,
        model="RSM",
        method="JMLE",
        noncenter_facet="Person",
        maxit=80,
        reltol=1e-6,
        keep_original=True,
    )

    summary = result["summary"].iloc[0]
    skipped = app.estimate_bias_interaction(result, {}, "Rater", "Task")
    quick_frames = app.build_result_bundle_frames(result, {})
    apa_tables = app._collect_apa_exportable_tables(result, {})

    assert result["identifiability"]["status"] == "eta_rank_deficient"
    assert int(summary["EtaStructuralNullity"]) == 7
    assert not bool(summary["InferenceReady"])
    assert "structurally rank deficient" in skipped["_skip_reason"]
    assert {
        "identifiability_summary",
        "identifiability_connectivity",
        "identifiability_null_space_block_energy",
        "identifiability_coordinate_null_weight",
    }.issubset(quick_frames)
    assert "Structural identifiability summary" in apa_tables
    assert "Person/facet connectivity" in apa_tables


def test_rank_full_fit_retains_optimizer_parameter_vector_exactly():
    row, generated = _generated("balanced_small")
    result = app.mfrm_estimate(
        generated.data,
        person_col="Person",
        facet_cols=["Rater", "Task", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=int(row["Categories"]) - 1,
        model="RSM",
        method="JMLE",
        noncenter_facet="Person",
        maxit=80,
        reltol=1e-6,
        keep_original=True,
    )
    sizes = app.build_param_sizes(result["config"])
    expanded = app.expand_params(result["opt"].x, sizes, result["config"])

    assert result["identifiability"]["status"] == "eta_identified"
    assert bool(result["summary"].iloc[0]["InferenceReady"]) == bool(result["opt"].success)
    assert (expanded["theta"] == result["params"]["theta"]).all()
    for facet in result["config"]["facet_names"]:
        assert (expanded["facets"][facet] == result["params"]["facets"][facet]).all()
