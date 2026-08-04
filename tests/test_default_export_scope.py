"""Exact-name boundary checks for the standalone default export surface."""

from __future__ import annotations

import inspect
import io
import zipfile
from types import SimpleNamespace

import pandas as pd

import streamlit_app as app


FORBIDDEN_DEFAULT_FRAMES = frozenset({
    "external_simulation_reference_inventory",
    "external_simulation_template_inventory",
    "reproducibility_script_matrix",
    "bayesian_mfrm_stan_refinement_plan",
    "stan_reproducibility_archive_contract",
    "stan_posterior_reproducibility_route",
    "stan_posterior_handoff_checklist",
    "stan_run_manifest_template",
    "mfrm_stan_data_manifest",
    "mfrm_stan_data_dictionary",
    "mfrm_stan_id_index_map",
    "mfrm_stan_prior_guidance",
    "mfrm_stan_prior_sensitivity_grid",
    "mfrm_stan_prior_decision_log_template",
    "mfrm_complete_stan_reproducibility_manifest",
    "mfrm_uto_bayesian_mfrm_design_audit",
    "mfrm_uto_bayesian_mfrm_claim_wording",
    "mfrm_uto_bayesian_mfrm_mapping_manifest",
    "mfrm_uto_bayesian_mfrm_data_dictionary",
    "mfrm_uto_bayesian_mfrm_id_index_map",
    "mfrmr_015_migration_coverage",
    "mfrmr_016_migration_coverage",
    "mfrmr_020_migration_coverage",
})

FORBIDDEN_DEFAULT_FILES = frozenset({
    "README_external_simulation_templates.md",
    "simulation_validation_python_template.py",
    "simulation_validation_r_template.R",
    "simulation_validation_julia_template.jl",
    "mfrm_analysis.R",
    "mfrm_yardstick_geom_text.R",
    "mfrm_yardstick_makie.jl",
    "apply_rating_scale_recodes.R",
    "apply_rating_scale_recodes.jl",
    "MFRM_Cross_Engine_Validation_Bundle.zip",
    "MFRM_External_Simulation_Validation_Templates.zip",
    "MFRM_Complete_Stan_Reproducibility_Package.zip",
    "MFRM_Bayesian_Stan_Runners.zip",
    "MFRM_Posterior_Viewer_Example.zip",
    "MFRM_Stan_Data_Package.zip",
    "README_bayesian_mfrm_stan_runners.md",
    "run_bayesian_mfrm_cmdstan_cli.jl",
    "run_bayesian_mfrm_cmdstanpy.py",
    "run_bayesian_mfrm_cmdstanr.R",
    "stan_posterior_handoff_checklist.csv",
    "stan_posterior_reproducibility_handoff.md",
    "stan_posterior_reproducibility_route.csv",
    "stan_reproducibility_archive_contract.csv",
    "stan_run_manifest_template.csv",
    "mfrm_stan_data.json",
    "mfrm_stan_id_index_map.csv",
})

FORBIDDEN_DEFAULT_ARCHIVE_NAMES = frozenset({
    *FORBIDDEN_DEFAULT_FILES,
    *(f"{frame_name}.csv" for frame_name in FORBIDDEN_DEFAULT_FRAMES),
})

REQUIRED_PYTHON_NATIVE_FILES = frozenset({
    "mfrm_app_engine_runner.py",
    "requirements.txt",
    "mfrm_local_batch_workflow.md",
    "mfrm_jmle_self_contained.py",
    "mfrm_yardstick_plotly.py",
    "apply_rating_scale_recodes.py",
})

OBSOLETE_DEMO_ARTIFACT_REFERENCES = frozenset({
    "MFRM_Tables.zip",
    "MFRM_Report.xlsx",
    "MFRM_Report.html",
    "MFRM_Publication_Figures.zip",
    "MFRM_Manuscript_Binder.zip",
    "MFRM_OSF_Package.zip",
    "mfrm_method_appendix.md",
    "mfrm_manuscript_template.md",
    "mfrm_manuscript_handoff.md",
    "analysis_identity.json",
    "evidence_records.json",
})


def _minimal_result_and_diagnostics():
    result = {
        "config": {
            "app_version": app.APP_VERSION,
            "model": "RSM",
            "method": "JMLE",
            "facet_names": ["Rater"],
            "n_cat": 5,
        },
        "prep": {
            "n_obs": 12,
            "n_person": 4,
            "rating_min": 0,
            "rating_max": 4,
            "score_map": pd.DataFrame({"RawScore": [0, 1], "MappedScore": [0, 1]}),
            "audit_summary": pd.DataFrame({"Check": ["rows"], "Status": ["ok"]}),
        },
        "summary": pd.DataFrame([{"Model": "RSM", "Method": "JMLE", "Converged": True}]),
        "convergence": pd.DataFrame([{"Converged": True}]),
        "steps": pd.DataFrame({"Step": [1, 2], "Estimate": [-0.5, 0.5]}),
        "facets": {
            "person": pd.DataFrame({"Facet": ["Person"], "Level": ["P1"], "Estimate": [0.0]}),
            "others": pd.DataFrame({"Facet": ["Rater"], "Level": ["R1"], "Estimate": [0.0]}),
        },
        "opt": SimpleNamespace(success=True, message="ok"),
    }
    diagnostics = {
        "measures": pd.DataFrame({
            "Facet": ["Rater"],
            "Level": ["R1"],
            "Estimate": [0.0],
            "SE": [0.1],
        }),
        "reliability": pd.DataFrame({"Facet": ["Rater"], "Reliability": [0.8]}),
        "fit": pd.DataFrame({
            "Facet": ["Rater"],
            "Level": ["R1"],
            "Infit": [1.0],
            "Outfit": [1.0],
        }),
        "obs": pd.DataFrame(),
    }
    data = pd.DataFrame({"Person": ["P1"], "Rater": ["R1"], "Score": [1]})
    return result, diagnostics, data


def test_default_download_frames_exclude_exact_legacy_names():
    result, diagnostics, _ = _minimal_result_and_diagnostics()
    frames, context = app.collect_download_frames(
        result,
        diagnostics,
        {},
        pd.DataFrame(),
        pd.DataFrame(),
        public_export_mode=True,
    )

    assert FORBIDDEN_DEFAULT_FRAMES.isdisjoint(frames)
    assert {"summary", "measures", "public_release_readiness"}.issubset(frames)
    assert "generic_stan_data_dl" not in context
    assert "uto_stan_data_dl" not in context
    method_areas = set(frames["method_reference_audit"]["MethodArea"].astype(str))
    assert "External R ecosystem and cross-package reference checks" not in method_areas
    assert "Recent MFRM extensions and Bayesian/rater drift context" not in method_areas


def test_demo_archive_keeps_python_assets_and_excludes_exact_legacy_names():
    result, diagnostics, data = _minimal_result_and_diagnostics()
    frames = app.build_demo_report_frames(result, diagnostics, data)
    readiness = frames.get("final_report_readiness", pd.DataFrame())
    text_assets = {
        **app.build_evidence_contract_text_assets(result, readiness),
        **app.python_native_reproducibility_assets(result, diagnostics),
    }
    archive = app.build_osf_zip(
        frames,
        title="MFRM_Demo_Report",
        text_assets=text_assets,
    )
    with zipfile.ZipFile(io.BytesIO(archive), "r") as bundle:
        archive_names = set(bundle.namelist())

    assert FORBIDDEN_DEFAULT_FRAMES.isdisjoint(frames)
    assert FORBIDDEN_DEFAULT_ARCHIVE_NAMES.isdisjoint(archive_names)
    assert REQUIRED_PYTHON_NATIVE_FILES.issubset(archive_names)
    assert {
        "final_report_readiness_evidence_contract.json",
        "final_report_readiness_evidence_records.csv",
    }.issubset(archive_names)
    handoff_files = " ".join(
        frames["manuscript_handoff_checklist"]["DownloadFile"].astype(str)
    )
    assert "MFRM_Demo_Report.zip" in handoff_files
    assert "MFRM_Demo_Publication_Figures.zip" in handoff_files
    assert "method_appendix.md" in handoff_files
    assert "MFRM_Tables.zip" not in handoff_files
    assert "mfrm_method_appendix.md" not in handoff_files


def test_demo_artifact_profile_matches_real_files_and_is_idempotent():
    result, diagnostics, data = _minimal_result_and_diagnostics()
    frames = app.build_demo_report_frames(result, diagnostics, data)
    frame_text = "\n".join(
        frame.to_csv(index=False)
        for frame in frames.values()
        if isinstance(frame, pd.DataFrame)
    )

    assert not {
        match.group(0)
        for match in app._DEMO_ARTIFACT_REFERENCE_PATTERN.finditer(frame_text)
    }
    expected_references = {
        "MFRM_Demo_Report.zip",
        "MFRM_Demo_Report.xlsx",
        "MFRM_Demo_Report.html",
        "MFRM_Demo_Publication_Figures.zip",
        "MFRM_Demo_Manuscript_Binder.zip",
        "method_appendix.md",
        "final_report_readiness_analysis_identity.json",
        "final_report_readiness_evidence_contract.json",
    }
    assert all(reference in frame_text for reference in expected_references)

    source_text = "; ".join(sorted(OBSOLETE_DEMO_ARTIFACT_REFERENCES))
    profiled_once = app._demo_artifact_reference_text(source_text)
    profiled_twice = app._demo_artifact_reference_text(profiled_once)
    assert profiled_once == profiled_twice
    assert not {
        match.group(0)
        for match in app._DEMO_ARTIFACT_REFERENCE_PATTERN.finditer(str(profiled_once))
    }

    manifest = frames["export_privacy_manifest"]
    assert "sample_data" in manifest["Frame"].astype(str).tolist()
    assert "export_privacy_manifest" in manifest["Frame"].astype(str).tolist()
    assert set(manifest["Status"].astype(str)) == {
        "included_synthetic_demo",
        "included",
    }
    assert manifest.iloc[-1]["Rows"] == len(manifest)

    binder_assets = app._demo_artifact_reference_assets(
        app.build_manuscript_binder_assets(
            frames,
            {"manuscript_handoff.md": "MFRM_OSF_Package.zip"},
            public_export_mode=True,
        )
    )
    binder_text = "\n".join(
        value for value in binder_assets.values() if isinstance(value, str)
    )
    assert not {
        match.group(0)
        for match in app._DEMO_ARTIFACT_REFERENCE_PATTERN.finditer(binder_text)
    }
    assert "MFRM_Demo_Report.zip" in binder_text


def test_demo_export_profiles_late_figure_and_binder_assets():
    export_source = inspect.getsource(app.export_demo_report)
    assert "_demo_artifact_reference_frame" in export_source
    assert "_demo_artifact_reference_assets" in export_source
    assert "_refresh_demo_export_privacy_manifest" in export_source


def test_download_ui_does_not_offer_exact_legacy_archive_names():
    render_source = inspect.getsource(app._render_downloads)
    offered_legacy_names = {
        name for name in FORBIDDEN_DEFAULT_FILES if name in render_source
    }
    assert offered_legacy_names == set()
