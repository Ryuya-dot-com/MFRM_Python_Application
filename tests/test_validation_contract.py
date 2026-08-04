import io
import inspect
import json
import zipfile

import pandas as pd
import pytest

import streamlit_app as app


@pytest.mark.legacy_compat
def test_cross_package_validation_plan_contract():
    plan = app.cross_package_validation_plan()
    assert not plan.empty
    assert {"TAM", "mirt", "sirt"}.issubset(set(plan["ExternalReference"].str.split().str[0]))
    assert "mml_rsm_latent_regression" in " ".join(plan["PythonScenario"].astype(str))
    assert "mml_gpcm_task_steps" in " ".join(plan["PythonScenario"].astype(str))

    notes = app.cross_package_parameterization_notes()
    assert "Latent variance" in notes["Topic"].tolist()
    assert "GPCM slope identification" in notes["Topic"].tolist()

    tolerance = app.cross_package_tolerance_policy()
    assert "Log-likelihood" in tolerance["EvidenceType"].tolist()
    assert "Plausible-value draws" in tolerance["EvidenceType"].tolist()

    docs = app.external_reference_documentation_table()
    assert {"TAM tam.mml.mfr", "mirt mirt", "sirt rm.facets", "mfrmr 0.1.6 package source"}.issubset(set(docs["Reference"]))

    simulation_inventory = app.external_simulation_reference_inventory()
    assert {"Main engine refit sweep", "Validation input replicates"}.issubset(set(simulation_inventory["ReferenceSet"]))
    assert simulation_inventory["PublicHandling"].str.contains("Do not", case=False, na=False).any()

    simulation_templates = app.external_simulation_template_inventory()
    assert {"Python", "R", "Julia"}.issubset(set(simulation_templates["Runtime"]))
    script_bundle = app.external_simulation_template_scripts()
    assert {
        "simulation_validation_python_template.py",
        "simulation_validation_r_template.R",
        "simulation_validation_julia_template.jl",
    }.issubset(set(script_bundle))
    joined = "\n".join(script_bundle.values())
    assert "MFRM_INPUT_CSV" in joined
    assert "/Users/" not in joined
    assert "C:/Users/" not in joined

    script_matrix = app.reproducibility_script_export_matrix()
    assert not script_matrix.empty
    assert {
        "Python",
        "R",
        "Julia",
        "Python + CmdStanPy",
        "R + CmdStanR",
        "Julia + CmdStan CLI",
    }.issubset(set(script_matrix["Runtime"]))
    assert "mfrm_app_engine_runner.py" in script_matrix["Artifact"].tolist()
    assert "run_bayesian_mfrm_cmdstan_cli.jl" in script_matrix["Artifact"].tolist()
    assert "MFRM_Complete_Stan_Reproducibility_Package.zip" in script_matrix["Artifact"].tolist()

    stan_plan = app.bayesian_mfrm_stan_refinement_plan()
    stan_plan_text = " ".join(stan_plan.astype(str).to_numpy().ravel().tolist())
    assert "Uto_Ueno_2020" in stan_plan_text
    assert "Uto_2021" in stan_plan_text
    assert "Uto_2023" in stan_plan_text
    assert "rater severity drift" in stan_plan_text.lower()

    stan_runners = app.bayesian_stan_runner_templates()
    assert {
        "stan_reproducibility_archive_contract.csv",
        "stan_posterior_reproducibility_handoff.md",
        "stan_posterior_reproducibility_route.csv",
        "stan_posterior_handoff_checklist.csv",
        "stan_run_manifest_template.csv",
        "run_bayesian_mfrm_cmdstanpy.py",
        "run_bayesian_mfrm_cmdstanr.R",
        "run_bayesian_mfrm_cmdstan_cli.jl",
    }.issubset(set(stan_runners))
    stan_runner_text = "\n".join(stan_runners.values())
    assert "MFRM_STAN_FILE" in stan_runner_text
    assert "MFRM_STAN_DATA_JSON" in stan_runner_text
    assert "mfrm_stan_prior_guidance.csv" in stan_runner_text
    assert "mfrm_stan_prior_decision_log_template.csv" in stan_runner_text
    assert "mfrm_stan_prior_sensitivity_grid.csv" in stan_runner_text
    assert "MFRM_Complete_Stan_Reproducibility_Package.zip" in stan_runner_text
    assert "stan_run_manifest.json" in stan_runner_text
    assert "stan_run_manifest.csv" in stan_runner_text
    assert "Stan Posterior Reproducibility Handoff" in stan_runner_text
    assert "Prior Setting And Sensitivity" in stan_runner_text
    assert "Posterior Viewer" in stan_runner_text
    assert "stan_file_sha256" in stan_runner_text
    assert "data_json_sha256" in stan_runner_text
    assert "/Users/" not in stan_runner_text
    assert "C:/Users/" not in stan_runner_text
    stan_runner_zip = app.cached_mixed_asset_zip(
        stan_runners,
        app.bytes_mapping_fingerprint(stan_runners),
    )
    with zipfile.ZipFile(io.BytesIO(stan_runner_zip), "r") as zf:
        runner_zip_names = set(zf.namelist())
    assert {
        "README_bayesian_mfrm_stan_runners.md",
        "stan_reproducibility_archive_contract.csv",
        "stan_posterior_reproducibility_handoff.md",
        "stan_posterior_reproducibility_route.csv",
        "stan_posterior_handoff_checklist.csv",
        "stan_run_manifest_template.csv",
        "run_bayesian_mfrm_cmdstanpy.py",
        "run_bayesian_mfrm_cmdstanr.R",
        "run_bayesian_mfrm_cmdstan_cli.jl",
    }.issubset(runner_zip_names)

    artifacts = app.external_validation_artifact_checklist()
    assert "External package versions" in artifacts["Artifact"].tolist()
    assert "Simulation reference inventory" in artifacts["Artifact"].tolist()
    assert "Simulation validation template scripts" in artifacts["Artifact"].tolist()
    assert "Bayesian Stan runner templates" in artifacts["Artifact"].tolist()
    assert "Stan data JSON and coding maps" in artifacts["Artifact"].tolist()
    assert "Bayesian prior-setting rationale" in artifacts["Artifact"].tolist()
    assert "Stan posterior handoff manifest" in artifacts["Artifact"].tolist()
    assert "Uto-family Bayesian MFRM refinement plan" in artifacts["Artifact"].tolist()
    stan_runner_artifact = artifacts.loc[artifacts["Artifact"] == "Bayesian Stan runner templates"].iloc[0]
    assert "stan_reproducibility_archive_contract.csv" in stan_runner_artifact["ExpectedFileOrEvidence"]
    assert "stan_posterior_reproducibility_handoff.md" in stan_runner_artifact["ExpectedFileOrEvidence"]
    assert "stan_posterior_reproducibility_route.csv" in stan_runner_artifact["ExpectedFileOrEvidence"]
    stan_data_artifact = artifacts.loc[artifacts["Artifact"] == "Stan data JSON and coding maps"].iloc[0]
    assert "mfrm_stan_prior_decision_log_template.csv" in stan_data_artifact["ExpectedFileOrEvidence"]
    assert "mfrm_complete_stan_reproducibility_manifest.csv" in stan_data_artifact["ExpectedFileOrEvidence"]
    stan_handoff_artifact = artifacts.loc[artifacts["Artifact"] == "Stan posterior handoff manifest"].iloc[0]
    assert "stan_reproducibility_archive_contract.csv" in stan_handoff_artifact["ExpectedFileOrEvidence"]
    assert "stan_posterior_reproducibility_handoff.md" in stan_handoff_artifact["ExpectedFileOrEvidence"]
    assert "stan_posterior_reproducibility_route.csv" in stan_handoff_artifact["ExpectedFileOrEvidence"]

    template = app.external_validation_report_template()
    assert "mfrmr 0.1.6 migration" in template["ClaimArea"].tolist()
    assert "mfrmr 0.2.0 migration" in template["ClaimArea"].tolist()
    assert "External Simulation numerical validation" in template["ClaimArea"].tolist()
    assert "Uto-family Bayesian MFRM Stan extension" in template["ClaimArea"].tolist()

    coverage = app.mfrmr_015_migration_coverage_table()
    assert "Bounded GPCM" in coverage["mfrmr015Area"].tolist()
    assert "Latent regression / population_formula" in coverage["mfrmr015Area"].tolist()

    coverage_020 = app.mfrmr_020_migration_coverage_table()
    # The 0.2.0 table inherits every 0.1.5 / 0.1.6 row and renames the area
    # column. Inherited rows must still be findable on the new column.
    assert "Bounded GPCM" in coverage_020["mfrmr020Area"].tolist()
    assert "Empirical-Bayes facet shrinkage" in coverage_020["mfrmr020Area"].tolist()
    # The slope-aware fair-average / structural SE work shipped in this
    # release must show up as a Ready row.
    assert "GPCM Linacre fair-average and structural SE" in coverage_020["mfrmr020Area"].tolist()
    assert "MML observed-information covariance" in coverage_020["mfrmr020Area"].tolist()
    # The slope-aware bias inference shipped: the row must say so. The
    # ``Ready`` status, the Muraki (1993) / Wilks (1938) / Cox (1975)
    # citations in the evidence string, and the references to the
    # publication-document integration are all part of the manuscript
    # claim surface the public README points at.
    bias_row = coverage_020.loc[
        coverage_020["mfrmr020Area"] == "GPCM bias inference - slope-aware"
    ]
    assert not bias_row.empty
    assert bias_row.iloc[0]["PythonStatus"] == "Ready"
    bias_evidence = str(bias_row.iloc[0]["PythonEvidence"])
    assert "Muraki, 1993" in bias_evidence
    assert "Wilks, 1938" in bias_evidence
    assert "Cox, 1975" in bias_evidence
    assert "LR ChiSq" in bias_evidence
    assert "profile-likelihood" in bias_evidence.lower() or "Profile CI" in bias_evidence
    # The bounded-GPCM row must have been overridden to reflect the
    # slope-aware fair-average / structural SE delivery; the fair-average
    # half of the historic limitation is gone.
    gpcm_row = coverage_020.loc[coverage_020["mfrmr020Area"] == "Bounded GPCM"]
    assert not gpcm_row.empty
    assert "slope-aware" in gpcm_row.iloc[0]["PythonStatus"].lower()


def test_legacy_extensions_are_excluded_from_native_claim_surfaces():
    result = {
        "config": {"model": "RSM", "method": "JMLE", "facet_names": ["Rater", "Task"], "n_cat": 5},
        "prep": {"n_obs": 10, "n_person": 3, "rating_min": 0, "rating_max": 4},
    }
    diagnostics = {"obs": app.pd.DataFrame({"StdResidual": [0.0, 0.5]})}

    guide = app.build_manuscript_claim_guide(result, diagnostics, {})
    areas = set(guide["ManuscriptArea"].astype(str))
    assert "Bayesian / Uto-family extension" not in areas
    assert "External package comparison" not in areas

    matrix = app.build_claim_to_evidence_matrix(result, diagnostics, {})
    matrix_areas = set(matrix["ManuscriptArea"].astype(str))
    assert "Bayesian / Uto-family extension" not in matrix_areas
    assert "External package comparison" not in matrix_areas


@pytest.mark.legacy_compat
def test_stan_reproducibility_archive_contract_tracks_package_boundaries():
    public_contract = app.stan_reproducibility_archive_contract_table(public_export_mode=True)
    private_contract = app.stan_reproducibility_archive_contract_table(public_export_mode=False)

    assert {
        "PackageOrSurface",
        "ExpectedFiles",
        "Role",
        "CheckBeforeSharing",
        "PrivacyBoundary",
    }.issubset(public_contract.columns)
    assert {
        "MFRM_Bayesian_Stan_Runners.zip",
        "MFRM_Stan_Data_Package.zip",
        "MFRM_OSF_Package.zip",
        "MFRM_Manuscript_Binder.zip",
        "Posterior Viewer return path",
    }.issubset(set(public_contract["PackageOrSurface"]))
    public_data_row = public_contract.loc[
        public_contract["PackageOrSurface"] == "MFRM_Stan_Data_Package.zip"
    ].iloc[0]
    private_data_row = private_contract.loc[
        private_contract["PackageOrSurface"] == "MFRM_Stan_Data_Package.zip"
    ].iloc[0]
    assert "stan_posterior_reproducibility_handoff.md" in public_data_row["ExpectedFiles"]
    assert "stan_reproducibility_archive_contract.csv" in public_data_row["ExpectedFiles"]
    assert "mfrm_stan_data.json" not in public_data_row["ExpectedFiles"]
    assert "mfrm_stan_data.json" in private_data_row["ExpectedFiles"]
    assert "private" in private_data_row["PrivacyBoundary"].lower()


def _stan_export_fixture_result(model: str = "RSM", *, step_facet: str | None = None) -> dict:
    data = app.pd.DataFrame(
        {
            "Person": ["P1", "P1", "P2", "P2", "P3", "P3"],
            "Rater": ["R1", "R2", "R1", "R2", "R1", "R2"],
            "Task": ["T1", "T1", "T2", "T2", "T1", "T2"],
            "Criterion": ["C1", "C2", "C1", "C2", "C1", "C2"],
            "TimeBlock": ["B1", "B1", "B2", "B2", "B1", "B2"],
            "Score": [0, 1, 2, 3, 4, 2],
        }
    )
    facets = ["Rater", "Task", "Criterion", "TimeBlock"]
    prep = app.prepare_mfrm_data(
        data,
        person_col="Person",
        facet_cols=facets,
        score_col="Score",
        rating_min=0,
        rating_max=4,
    )
    return {
        "config": {
            "model": model,
            "method": "JMLE",
            "facet_names": facets,
            "n_cat": 5,
            "step_facet": step_facet,
        },
        "prep": prep,
    }


@pytest.mark.legacy_compat
def test_generic_stan_data_export_is_cmdstan_json_with_private_maps():
    result = _stan_export_fixture_result(model="PCM", step_facet="Task")

    export = app.build_generic_mfrm_stan_data_export(result)

    assert export["available"]
    payload = json.loads(export["json_text"])
    assert payload["N"] == 6
    assert payload["J"] == 3
    assert payload["C"] == 5
    assert {"K_rater", "K_task", "K_criterion", "K_timeblock"}.issubset(payload)
    for array_name, upper_name in {
        "person": "J",
        "rater": "K_rater",
        "task": "K_task",
        "criterion": "K_criterion",
        "timeblock": "K_timeblock",
        "y": "C",
    }.items():
        values = payload[array_name]
        assert len(values) == payload["N"]
        assert min(values) >= 1
        assert max(values) <= payload[upper_name]

    id_map = export["id_map"]
    assert {"OriginalLabel", "StanIndex", "StanVariable"}.issubset(id_map.columns)
    assert id_map.loc[id_map["StanVariable"] == "person", "StanIndex"].min() == 1
    assert "row-level coded data" in " ".join(export["data_dictionary"]["PrivacyBoundary"].astype(str)).lower()
    pcm_row = export["manifest"].loc[export["manifest"]["Check"] == "PCM step facet"].iloc[0]
    assert "`Task` / `task`" in pcm_row["Value"]
    prior_row = export["manifest"].loc[export["manifest"]["Check"] == "Prior settings"].iloc[0]
    assert "mfrm_stan_prior_guidance.csv" in prior_row["Action"]
    assert {"PriorTarget", "TemplatePrior", "HowToSet", "SensitivityCheck"}.issubset(export["prior_guidance"].columns)
    assert "sigma_theta_prior_scale" in export["prior_sensitivity_grid"].columns


@pytest.mark.legacy_compat
def test_stan_data_assets_respect_public_row_level_boundary():
    result = _stan_export_fixture_result()

    private_assets = app.stan_data_export_assets(result, include_row_level=True)
    public_assets = app.stan_data_export_assets(result, include_row_level=False)

    assert "mfrm_stan_data.json" in private_assets
    assert "mfrm_stan_id_index_map.csv" in private_assets
    assert "mfrm_stan_prior_decision_log_template.csv" in private_assets
    assert "mfrm_uto_bayesian_mfrm_data_template.json" in private_assets
    assert "mfrm_uto_bayesian_mfrm_design_audit.csv" in private_assets
    assert "mfrm_uto_bayesian_mfrm_claim_wording.csv" in private_assets
    assert "mfrm_stan_data_dictionary.csv" in public_assets
    assert "mfrm_stan_data_manifest.csv" in public_assets
    assert "mfrm_stan_prior_guidance.csv" in public_assets
    assert "mfrm_stan_prior_sensitivity_grid.csv" in public_assets
    assert "mfrm_stan_prior_decision_log_template.csv" in public_assets
    assert "mfrm_uto_bayesian_mfrm_design_audit.csv" in public_assets
    assert "mfrm_uto_bayesian_mfrm_claim_wording.csv" in public_assets
    assert "stan_reproducibility_archive_contract.csv" in public_assets
    assert "stan_posterior_reproducibility_handoff.md" in public_assets
    assert "stan_posterior_reproducibility_route.csv" in public_assets
    assert "stan_posterior_handoff_checklist.csv" in public_assets
    assert "stan_run_manifest_template.csv" in public_assets
    assert "mfrm_stan_data.json" not in public_assets
    assert "mfrm_stan_id_index_map.csv" not in public_assets
    assert "mfrm_uto_bayesian_mfrm_data_template.json" not in public_assets

    public_zip = app.cached_mixed_asset_zip(
        public_assets,
        app.bytes_mapping_fingerprint(public_assets),
    )
    with zipfile.ZipFile(io.BytesIO(public_zip), "r") as zf:
        public_zip_names = set(zf.namelist())
    assert "stan_reproducibility_archive_contract.csv" in public_zip_names
    assert "stan_posterior_reproducibility_handoff.md" in public_zip_names
    assert "mfrm_stan_data.json" not in public_zip_names
    assert "mfrm_stan_id_index_map.csv" not in public_zip_names


@pytest.mark.legacy_compat
def test_complete_stan_reproducibility_package_includes_models_priors_and_privacy_manifest():
    result = _stan_export_fixture_result(model="PCM", step_facet="Task")

    generic_code = app.build_generic_mfrm_stan_code(result)
    assert generic_code["available"]
    assert "matrix[K_task, C-1] tau" in generic_code["stan_code"]
    assert "Sensitivity notice" in generic_code["stan_code"]

    private_assets = app.stan_reproducibility_package_assets(result, include_row_level=True)
    assert {
        "README_complete_stan_reproducibility_package.md",
        "mfrm_model.stan",
        "mfrm_uto_bayesian_mfrm.stan",
        "mfrm_stan_data.json",
        "mfrm_stan_id_index_map.csv",
        "mfrm_stan_prior_guidance.csv",
        "mfrm_stan_prior_decision_log_template.csv",
        "mfrm_stan_prior_sensitivity_grid.csv",
        "mfrm_uto_bayesian_mfrm_design_audit.csv",
        "mfrm_uto_bayesian_mfrm_claim_wording.csv",
        "run_bayesian_mfrm_cmdstanpy.py",
        "run_bayesian_mfrm_cmdstanr.R",
        "run_bayesian_mfrm_cmdstan_cli.jl",
        "mfrm_complete_stan_reproducibility_manifest.csv",
    }.issubset(private_assets)
    assert "prior decision log" in private_assets["README_complete_stan_reproducibility_package.md"].lower()
    prior_log = pd.read_csv(io.StringIO(private_assets["mfrm_stan_prior_decision_log_template.csv"]))
    assert {"PriorTarget", "ChosenValue", "PriorDecisionRationale", "DiagnosticsToAttach"}.issubset(prior_log.columns)
    manifest = pd.read_csv(io.StringIO(private_assets["mfrm_complete_stan_reproducibility_manifest.csv"]))
    assert "mfrm_model.stan" in manifest["PackageFile"].tolist()
    private_row = manifest.loc[manifest["PackageFile"] == "mfrm_stan_data.json"].iloc[0]
    assert bool(private_row["ContainsRowLevelOrPrivateMap"])

    public_assets = app.stan_reproducibility_package_assets(result, include_row_level=False)
    assert "mfrm_stan_data.json" not in public_assets
    assert "mfrm_stan_id_index_map.csv" not in public_assets
    assert "mfrm_model.stan" in public_assets
    public_zip = app.cached_mixed_asset_zip(public_assets, app.bytes_mapping_fingerprint(public_assets))
    with zipfile.ZipFile(io.BytesIO(public_zip), "r") as zf:
        public_zip_names = set(zf.namelist())
    assert "mfrm_complete_stan_reproducibility_manifest.csv" in public_zip_names
    assert "mfrm_stan_data.json" not in public_zip_names


@pytest.mark.legacy_compat
def test_stan_prior_guidance_connects_defaults_to_sensitivity_variants():
    guidance = app.stan_prior_setting_guidance(
        sigma_theta_prior_scale=2.5,
        facet_prior_scale=2.0,
        step_prior_scale=5.0,
    )
    grid = app.stan_prior_sensitivity_grid(
        sigma_theta_prior_scale=2.5,
        facet_prior_scale=2.0,
        step_prior_scale=5.0,
    )

    assert not guidance.empty
    assert not grid.empty
    assert "Person ability SD sigma_theta" in guidance["PriorTarget"].tolist()
    assert "RSM/PCM thresholds tau" in guidance["PriorTarget"].tolist()
    assert {"current_template", "regularizing", "wide_scale", "threshold_regularized"}.issubset(
        set(grid["Variant"])
    )
    decision_log = app.stan_prior_decision_log_template(
        sigma_theta_prior_scale=2.5,
        facet_prior_scale=2.0,
        step_prior_scale=5.0,
    )
    assert "PriorDecisionRationale" in decision_log.columns
    assert decision_log["SensitivityVariantsToRun"].str.contains("regularizing", case=False, na=False).any()
    regularizing = grid.loc[grid["Variant"] == "regularizing"].iloc[0]
    wide = grid.loc[grid["Variant"] == "wide_scale"].iloc[0]
    assert regularizing["facet_prior_scale"] < 2.0
    assert wide["sigma_theta_prior_scale"] > 2.5


@pytest.mark.legacy_compat
def test_stan_manifest_handoff_tables_and_parser_validate_run_identity():
    checklist = app.stan_posterior_handoff_checklist()
    template = app.stan_run_manifest_template()

    assert "stan_run_manifest.json or stan_run_manifest.csv" in checklist["Artifact"].tolist()
    assert {"Field", "Required", "Meaning"}.issubset(template.columns)
    assert {"stan_file_sha256", "data_json_sha256", "posterior_csv_count", "posterior_csv_sha256"}.issubset(set(template["Field"]))

    manifest_payload = {
        "schema_version": "mfrm_stan_run_manifest_v1",
        "producer": "MFRM Streamlit Bayesian Stan runner",
        "runtime": "Python + CmdStanPy",
        "stan_file_name": "mfrm_model.stan",
        "stan_file_sha256": "a" * 64,
        "data_json_name": "mfrm_stan_data.json",
        "data_json_sha256": "b" * 64,
        "prior_scales": {
            "sigma_theta_prior_scale": 2.5,
            "facet_prior_scale": 2.0,
            "step_prior_scale": 5.0,
        },
        "seed": 42,
        "chains": 4,
        "iter_warmup": 1000,
        "iter_sampling": 2000,
        "adapt_delta": 0.95,
        "max_treedepth": 12,
        "cmdstan_version": "2.35.0",
        "posterior_csv_files": ["chain_1.csv", "chain_2.csv", "chain_3.csv", "chain_4.csv"],
        "posterior_csv_count": 4,
        "posterior_csv_sha256": ";".join(
            f"chain_{i}.csv={'%064x' % i}" for i in range(1, 5)
        ),
        "summary_file": "summary.csv",
        "diagnose_file": "diagnose.txt",
    }
    upload = io.BytesIO(json.dumps(manifest_payload).encode("utf-8"))
    upload.name = "stan_run_manifest.json"

    flat, errors = app.parse_stan_run_manifest_upload(upload)
    assert not errors
    assert flat["prior_scales.sigma_theta_prior_scale"] == 2.5

    loaded_payload = {
        "n_chains": 4,
        "uploaded_file_count": 4,
        "uploaded_file_names": ["chain_1.csv", "chain_2.csv", "chain_3.csv", "chain_4.csv"],
        "uploaded_file_sha256": {
            f"chain_{i}.csv": "%064x" % i for i in range(1, 5)
        },
        "diagnostics": {
            "divergent": [[0, 0], [0, 0], [0, 0], [0, 0]],
            "treedepth": [[5, 6], [5, 6], [5, 6], [5, 6]],
            "energy": [[1.0, 1.4, 1.1], [1.2, 1.5, 1.0], [0.8, 1.2, 1.1], [1.1, 1.6, 1.3]],
        },
    }
    checks = app.stan_run_manifest_check_table(flat, loaded_payload)
    indexed = checks.set_index("Check")
    assert indexed.loc["Manifest uploaded", "Status"] == "Ready"
    assert indexed.loc["Stan source hash", "Status"] == "Ready"
    assert indexed.loc["Data JSON hash", "Status"] == "Ready"
    assert indexed.loc["Prior scales", "Status"] == "Ready"
    assert indexed.loc["Sampler controls", "Status"] == "Ready"
    assert indexed.loc["Posterior chain count", "Status"] == "Ready"
    assert indexed.loc["Posterior CSV filenames", "Status"] == "Ready"
    assert indexed.loc["Posterior CSV hashes", "Status"] == "Ready"
    assert indexed.loc["External diagnostic artifacts", "Status"] == "Ready"
    assert indexed.loc["Transition diagnostics loaded", "Status"] == "Ready"

    mismatch = dict(flat)
    mismatch["posterior_csv_count"] = 3
    mismatch["posterior_csv_sha256"] = "chain_1.csv=" + ("f" * 64)
    mismatch_checks = app.stan_run_manifest_check_table(mismatch, loaded_payload).set_index("Check")
    assert mismatch_checks.loc["Posterior chain count", "Status"] == "Review"
    assert mismatch_checks.loc["Posterior CSV hashes", "Status"] == "Review"


@pytest.mark.legacy_compat
def test_stan_manifest_csv_upload_accepts_field_value_rows():
    csv_upload = io.BytesIO(
        "\n".join([
            "Field,Value",
            "schema_version,mfrm_stan_run_manifest_v1",
            "chains,4",
            "posterior_csv_count,4",
        ]).encode("utf-8")
    )
    csv_upload.name = "stan_run_manifest.csv"

    flat, errors = app.parse_stan_run_manifest_upload(csv_upload)

    assert not errors
    assert flat["schema_version"] == "mfrm_stan_run_manifest_v1"
    assert int(flat["chains"]) == 4


@pytest.mark.legacy_compat
def test_posterior_viewer_example_package_matches_manifest_contract():
    assets = app.posterior_viewer_example_package_assets()

    required = {
        "README_posterior_viewer_example.md",
        "example_mfrm_model.stan",
        "example_mfrm_stan_data.json",
        "example_chain_1.csv",
        "example_chain_2.csv",
        "stan_run_manifest.json",
        "stan_run_manifest.csv",
    }
    assert required.issubset(set(assets))
    assert "/Users/" not in "\n".join(assets.values())
    assert "C:/Users/" not in "\n".join(assets.values())

    manifest = json.loads(assets["stan_run_manifest.json"])
    assert manifest["schema_version"] == "mfrm_stan_run_manifest_v1"
    assert manifest["posterior_csv_count"] == 2
    assert manifest["chains"] == 2
    assert manifest["stan_file_sha256"] == app.hashlib.sha256(
        assets["example_mfrm_model.stan"].encode("utf-8")
    ).hexdigest()
    assert manifest["data_json_sha256"] == app.hashlib.sha256(
        assets["example_mfrm_stan_data.json"].encode("utf-8")
    ).hexdigest()

    upload = io.BytesIO(assets["stan_run_manifest.json"].encode("utf-8"))
    upload.name = "stan_run_manifest.json"
    flat, errors = app.parse_stan_run_manifest_upload(upload)
    assert not errors
    checks = app.stan_run_manifest_check_table(flat, {"n_chains": 2}).set_index("Check")
    assert checks.loc["Manifest uploaded", "Status"] == "Ready"
    assert checks.loc["Stan source hash", "Status"] == "Ready"
    assert checks.loc["Data JSON hash", "Status"] == "Ready"
    assert checks.loc["Prior scales", "Status"] == "Ready"
    assert checks.loc["Sampler controls", "Status"] == "Ready"
    assert checks.loc["Posterior chain count", "Status"] == "Ready"
    assert "posterior_csv_sha256" in manifest


def test_facets_yardstick_uses_text_first_columns_and_thresholds():
    person_tbl = pd.DataFrame({"Estimate": [-1.0, -0.4, 0.1, 0.6, 1.2]})
    facet_tbl = pd.DataFrame({
        "Facet": ["Senior scientists", "Senior scientists", "Junior scientists", "Traits"],
        "Level": ["Brahe", "Avogadro", "Betty", "Enthusiasm"],
        "Estimate": [0.15, -0.05, 0.55, 0.75],
    })
    step_tbl = pd.DataFrame({
        "Step": ["1|2", "2|3", "3|4"],
        "Estimate": [-0.75, 0.0, 0.72],
    })

    fig = app._make_yardstick_figure(
        person_tbl,
        facet_tbl,
        step_tbl,
        show_direct_labels=True,
    )

    assert fig is not None
    annotation_text = " ".join(str(item.text) for item in fig.layout.annotations)
    assert "Senior scientists" in annotation_text
    assert "Thresholds" in annotation_text
    trace_modes = [getattr(trace, "mode", "") for trace in fig.data]
    assert "text" in trace_modes
    assert "markers+text" not in trace_modes
    threshold_line_traces = [
        trace for trace in fig.data
        if getattr(trace, "name", "") == "Threshold boundary lines"
    ]
    assert threshold_line_traces
    line_trace = threshold_line_traces[0]
    assert list(line_trace.x)[:3] == [0.08, 0.20, None]
    assert list(line_trace.y)[:2] == [-0.75, -0.75]
    trace_text_parts = []
    for trace in fig.data:
        text_values = getattr(trace, "text", None)
        if text_values is None:
            continue
        if isinstance(text_values, str):
            trace_text_parts.append(text_values)
        else:
            trace_text_parts.extend(str(value) for value in list(text_values))
    trace_text = " ".join(trace_text_parts)
    assert "Betty" in trace_text
    assert "Enthusiasm" in trace_text
    assert "1|2" in trace_text
    tick_traces = [
        trace for trace in fig.data
        if getattr(trace, "mode", "") == "lines"
        and len(getattr(trace, "x", [])) >= 3
        and list(trace.x)[:3] == [0.02, 0.1, None]
    ]
    assert tick_traces


def test_plot_label_lookup_supports_plotnumber_style_ids():
    plot_rows = pd.DataFrame({
        "Role": ["FacetElement", "FacetElement", "FacetElement", "Threshold"],
        "Facet": ["Rater", "Rater", "Criterion", "Common"],
        "RawLabel": ["Rater with a very long display name", "R2", "Grammar", "Step_1"],
        "DisplayLabel": ["Rater with a very long display name", "R2", "Grammar", "0|1"],
        "PlotColumn": ["Rater", "Rater", "Criterion", "Thresholds"],
    })

    lookup = app.build_plot_label_lookup(plot_rows, label_mode="Number", max_chars=12)

    assert lookup["PlotLabelID"].tolist() == ["R001", "R002", "C001", "TH001"]
    assert lookup["DisplayLabel"].tolist() == ["R001", "R002", "C001", "TH001"]
    assert lookup.loc[0, "PlotLabelShort"].endswith("...")
    assert lookup.loc[3, "PlotLabelFull"] == "0|1"


def test_krippendorff_alpha_is_one_for_perfect_shared_scores():
    obs = pd.DataFrame({
        "Person": ["P1", "P1", "P2", "P2", "P3", "P3"],
        "Task": ["T1", "T1", "T1", "T1", "T1", "T1"],
        "Rater": ["A", "B", "A", "B", "A", "B"],
        "Observed": [1, 1, 2, 2, 3, 3],
    })

    out = app.calc_krippendorff_alpha(
        obs,
        ["Person", "Task", "Rater"],
        "Rater",
        metric="nominal",
    )
    summary = out["summary"].iloc[0]

    assert abs(summary["KrippendorffAlpha"] - 1.0) < 1e-12
    assert summary["Status"] == "Ready"
    assert summary["PairableUnits"] == 3
    assert not out["coverage"].empty


def test_krippendorff_alpha_flags_systematic_disagreement():
    obs = pd.DataFrame({
        "Person": ["P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4"],
        "Task": ["T1"] * 8,
        "Rater": ["A", "B", "A", "B", "A", "B", "A", "B"],
        "Observed": [1, 3, 3, 1, 1, 3, 3, 1],
    })

    out = app.calc_krippendorff_alpha(
        obs,
        ["Person", "Task", "Rater"],
        "Rater",
        metric="interval",
    )
    alpha = float(out["summary"].iloc[0]["KrippendorffAlpha"])

    assert alpha < 0.0
    assert out["summary"].iloc[0]["Status"] == "Review"


def test_krippendorff_alpha_handles_missing_and_by_facet_exports():
    obs = pd.DataFrame({
        "Person": ["P1", "P1", "P2", "P3", "P3"],
        "Task": ["T1", "T1", "T1", "T1", "T1"],
        "Rater": ["A", "B", "A", "A", "B"],
        "Observed": [1, 1, 2, 1, 2],
    })

    out = app.calc_krippendorff_alpha(
        obs,
        ["Person", "Task", "Rater"],
        "Rater",
        metric="ordinal",
    )
    summary = out["summary"].iloc[0]
    assert summary["PairableUnits"] == 2
    assert summary["UnitsWithSingleRating"] == 1
    assert "Caveat" in out["summary"].columns

    all_facets = app.calc_krippendorff_alpha_by_facet(obs, ["Person", "Task", "Rater"])
    assert not all_facets.empty
    assert {"ordinal", "nominal", "interval"}.issubset(set(all_facets["DistanceMetric"]))


def _scoring_qc_fixture(*, low_agreement: bool = False):
    if low_agreement:
        observed = [1, 3, 3, 1, 1, 3, 3, 1]
    else:
        observed = [1, 1, 2, 2, 3, 3, 2, 2]
    obs = pd.DataFrame({
        "Person": ["P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4"],
        "Rater": ["A", "B", "A", "B", "A", "B", "A", "B"],
        "Criterion": ["C1"] * 8,
        "Observed": observed,
    })
    result = {
        "config": {
            "facet_names": ["Rater", "Criterion"],
            "model": "RSM",
            "method": "JMLE",
        },
        "prep": {
            "rating_min": 1,
            "rating_max": 3,
        },
    }
    diagnostics = {
        "obs": obs,
        "fit": pd.DataFrame({
            "Facet": ["Rater", "Rater"],
            "Level": ["A", "B"],
            "Infit": [1.0, 1.0],
            "Outfit": [1.0, 1.0],
            "InfitZSTD": [0.0, 0.0],
            "OutfitZSTD": [0.0, 0.0],
        }),
    }
    diagnostics["krippendorff_alpha"] = app.calc_krippendorff_alpha_by_facet(
        obs,
        ["Person", "Rater", "Criterion"],
    )
    return result, diagnostics


def test_scoring_consistency_decision_contract_and_bundle_exports():
    result, diagnostics = _scoring_qc_fixture()

    decision = app.build_scoring_consistency_decision(
        result,
        diagnostics,
        ["Rater", "Criterion"],
        rater_facet="Rater",
    )

    assert not decision.empty
    assert {
        "DecisionArea",
        "Status",
        "PrimaryEvidence",
        "RecommendedAction",
        "APAReporting",
        "EvidenceFile",
        "DoNotClaim",
        "ReportBlocking",
        "DecisionRule",
        "PlainLanguageSummary",
        "NextCheckpoint",
        "FACETSCrosswalk",
        "PrimaryAudience",
    }.issubset(decision.columns)
    assert "Overall scoring consistency" in decision["DecisionArea"].tolist()
    assert "Krippendorff alpha" in decision["DecisionArea"].tolist()

    bundle = app.build_result_bundle_frames(result, diagnostics)
    assert "scoring_consistency_decision" in bundle
    assert "quality_control_recommendations" in bundle
    assert "quality_control_todo_checklist" in bundle
    assert "scoring_quality_first_read_summary" in bundle
    assert "scoring_quality_term_guide" in bundle
    assert "help_reference_coverage" in bundle
    assert "krippendorff_alpha_interpretation" in bundle
    assert "report_ready_summary_panel" in bundle
    assert "role_based_action_memos" in bundle
    assert "role_based_action_checklist" in bundle
    assert "critical_final_review_panel" in bundle
    assert "status_rationale_drilldown" in bundle
    assert "result_reading_path" in bundle
    assert "operational_decision_board" in bundle
    assert "reporting_action_bridge" in bundle
    assert "category_collapse_action_center" in bundle
    assert "rating_scale_recode_candidates" in bundle
    assert "rating_scale_decision_support" in bundle


def test_quality_control_recommendations_prioritize_rater_training_on_low_alpha():
    result, diagnostics = _scoring_qc_fixture(low_agreement=True)

    qc = app.build_quality_control_recommendations(
        result,
        diagnostics,
        ["Rater", "Criterion"],
        rater_facet="Rater",
    )
    rater_training = qc.loc[qc["QCDecision"] == "Rater training"].iloc[0]

    assert rater_training["Recommendation"] == "Prioritize"
    assert rater_training["Priority"] == "High"
    assert rater_training["OwnerRole"] == "Scoring manager"
    assert "CompletionCriterion" in qc.columns
    assert "Krippendorff" in " ".join(
        app.build_krippendorff_alpha_interpretation(result, diagnostics)["ReportWording"].astype(str)
    )


def test_quality_control_todo_checklist_is_checkbox_ready_and_blocks_low_evidence():
    result, diagnostics = _scoring_qc_fixture(low_agreement=True)

    todo = app.build_quality_control_todo_checklist(
        result,
        diagnostics,
        ["Rater", "Criterion"],
        rater_facet="Rater",
    )

    assert not todo.empty
    assert {
        "Done",
        "Step",
        "ToDoItem",
        "BlockingForReport",
        "NeedsCaveat",
        "CompletionCriterion",
        "EvidenceToOpen",
        "ClaimUnlocked",
        "FACETSCrosswalk",
    }.issubset(todo.columns)
    assert todo["Done"].eq(False).all()
    assert todo["BlockingForReport"].astype(bool).any()
    assert "Finalize APA scoring-quality wording" in todo["ToDoItem"].tolist()


def test_scoring_quality_first_read_summary_points_to_next_action():
    result, diagnostics = _scoring_qc_fixture(low_agreement=True)

    summary = app.build_scoring_quality_first_read_summary(
        result,
        diagnostics,
        ["Rater", "Criterion"],
        rater_facet="Rater",
    )

    assert not summary.empty
    row = summary.iloc[0]
    assert row["ReadinessLabel"] == "Do not report yet"
    assert row["CanUseApaDraft"] == "No"
    assert int(row["BlockingItems"]) > 0
    assert row["NextAction"]
    assert row["FirstEvidenceToOpen"]


def test_scoring_quality_term_guide_covers_agreement_fit_and_category_terms():
    guide = app.scoring_quality_term_guide()

    assert not guide.empty
    terms = set(guide["Term"])
    assert {"Exact agreement", "Krippendorff alpha", "Rater fit/severity", "Category collapse"}.issubset(terms)
    assert {"PlainMeaning", "WhyItMatters", "LookAt", "DoNotConfuseWith"}.issubset(guide.columns)


def test_apa_scoring_quality_draft_is_conservative_and_report_ready():
    result, diagnostics = _scoring_qc_fixture()
    decision = app.build_scoring_consistency_decision(result, diagnostics, ["Rater", "Criterion"], rater_facet="Rater")
    qc = app.build_quality_control_recommendations(result, diagnostics, ["Rater", "Criterion"], rater_facet="Rater")

    draft = app.generate_apa_scoring_quality_draft(result, diagnostics, decision, qc)

    assert "# APA Scoring Quality Draft" in draft
    assert "Krippendorff" in draft
    assert "MFRM" in draft
    assert "Category collapse" in draft
    assert "Quality-control gate" in draft
    assert "First-read summary" in draft
    assert "FACETS" in draft
    assert "help_reference_coverage.csv" in draft
    assert "Krippendorff, 2004" in draft
    assert "Do not overclaim" in draft


def test_yardstick_panel_keeps_direct_text_labels_by_default():
    source = inspect.getsource(app._draw_yardstick)
    assert "st.checkbox" not in source
    assert "show_direct_labels = True" in source
    assert "dense_yardstick_caption" in source


def test_plot_label_mode_defaults_to_full_label():
    assert app.PLOT_LABEL_MODE_DEFAULT == "Full label"
    assert app._normalize_plot_label_mode(None) == "Full label"
    prefs = app.get_visualization_preferences()
    assert prefs["plot_label_mode"] == "Full label"


def test_facets_yardstick_help_table_matches_current_renderer_contract():
    table = app.facets_yardstick_help_table()
    joined = " ".join(table.astype(str).to_numpy().ravel())
    assert "direct text labels" in joined
    assert "PlotY" in joined
    assert "TextLane" in joined
    assert "mfrm_yardstick_map.csv" in joined
    source = inspect.getsource(app.show_help_section)
    assert "facets_yardstick_help_table()" in source


def test_yardstick_threshold_labels_use_rating_scale_boundaries():
    step_tbl = pd.DataFrame({
        "Step": ["Step_1", "Step_2", "Step_5"],
        "Estimate": [-2.0, -1.0, 2.0],
    })

    thresholds = app._yardstick_threshold_frame(step_tbl, rating_min=0)

    assert thresholds["DisplayLabel"].tolist() == ["0|1", "1|2", "4|5"]
    assert thresholds["BoundaryLower"].tolist() == [0.0, 1.0, 4.0]
    assert thresholds["BoundaryUpper"].tolist() == [1.0, 2.0, 5.0]


def test_yardstick_label_layout_dodges_dense_direct_text():
    estimates = pd.Series([0.10, 0.11, 0.12, 0.13])

    text_x, plot_y, lanes = app._yardstick_label_layout(estimates)

    assert len(text_x) == len(estimates)
    assert len(plot_y) == len(estimates)
    assert len(set(lanes)) > 1
    assert max(plot_y) - min(plot_y) > max(estimates) - min(estimates)


def test_yardstick_export_and_reproduction_scripts_contract():
    person_tbl = pd.DataFrame({"Person": ["P1", "P2"], "Estimate": [-0.5, 0.4]})
    facet_tbl = pd.DataFrame({
        "Facet": ["Rater", "Rater", "Criterion"],
        "Level": ["R1", "R2", "Mechanics"],
        "Estimate": [0.1, 0.12, 0.7],
    })
    step_tbl = pd.DataFrame({"Step": ["Step_1", "Step_2"], "Estimate": [-1.0, 1.0]})

    yardstick = app.make_yardstick_export_table(
        person_tbl,
        facet_tbl,
        step_tbl,
        rating_min=0,
    )
    assert {
        "Role",
        "Facet",
        "DisplayLabel",
        "PlotLabelID",
        "PlotLabelShort",
        "PlotLabelFull",
        "PlotLabelMode",
        "Estimate",
        "PlotY",
        "TextLane",
        "LineXStart",
        "LineXEnd",
        "PlotColumn",
    }.issubset(yardstick.columns)
    assert {"Person", "FacetElement", "Threshold"}.issubset(set(yardstick["Role"]))
    assert "0|1" in yardstick["DisplayLabel"].astype(str).tolist()
    assert "Step_1" in yardstick["RawLabel"].astype(str).tolist()
    direct_text = yardstick[yardstick["PlotKind"] == "direct_text"]
    assert direct_text["PlotY"].notna().all()
    assert direct_text["TextLane"].notna().all()
    threshold_rows = yardstick[yardstick["Role"] == "Threshold"]
    assert threshold_rows["LineXStart"].notna().all()
    assert threshold_rows["LineXEnd"].notna().all()
    assert (threshold_rows["LineXStart"] < threshold_rows["LineXEnd"]).all()

    numbered_yardstick = app.make_yardstick_export_table(
        person_tbl,
        facet_tbl,
        step_tbl,
        rating_min=0,
        plot_label_mode="Number",
    )
    assert {"R001", "R002", "C001", "TH001", "TH002"}.issubset(
        set(numbered_yardstick["DisplayLabel"].astype(str))
    )
    assert "Step_1" in numbered_yardstick["RawLabel"].astype(str).tolist()
    assert numbered_yardstick["PlotLabelFull"].astype(str).str.len().gt(0).all()

    scripts = app.python_yardstick_reproducibility_assets()
    assert {
        "README_facets_yardstick_reproduction.md",
        "mfrm_yardstick_plotly.py",
    }.issubset(set(scripts))
    assert not any(name.endswith((".R", ".jl")) for name in scripts)
    joined = "\n".join(scripts.values())
    assert "mfrm_yardstick_map.csv" in joined
    assert "PlotY" in joined
    assert "TextLane" in joined
    assert "LineXStart" in joined
    assert "LineXEnd" in joined
    assert "/Users/" not in joined
    assert "C:/Users/" not in joined


@pytest.mark.legacy_compat
def test_uto_family_stan_data_template_maps_current_design_fields():
    result = _stan_export_fixture_result()

    export = app.build_uto_bayesian_mfrm_stan_data_export(result)

    assert export["available"]
    payload = json.loads(export["json_text"])
    assert {
        "N",
        "J",
        "R",
        "I",
        "C_dim",
        "S",
        "K",
        "person",
        "rater",
        "task",
        "criterion",
        "time_block",
        "y",
    }.issubset(payload)
    assert payload["R"] == 2
    assert payload["I"] == 2
    assert payload["C_dim"] == 2
    assert payload["S"] == 2
    assert max(payload["time_block"]) == payload["S"]
    assert export["mapping"]["rater_facet"] == "Rater"
    assert export["mapping"]["criterion_facet"] == "Criterion"
    assert export["mapping"]["time_facet"] == "TimeBlock"
    audit = export["design_audit"]
    assert {
        "DesignRole",
        "MappedFacet",
        "StanVariable",
        "Status",
        "ClaimReadiness",
        "AllowedClaimBoundary",
    }.issubset(audit.columns)
    assert {
        "Rater severity facet",
        "Rubric criterion / latent dimension",
        "Time/order block for rater drift",
    }.issubset(set(audit["DesignRole"]))
    time_row = audit.loc[audit["DesignRole"] == "Time/order block for rater drift"].iloc[0]
    assert time_row["Status"] == "Ready"
    assert "drift" in time_row["ClaimReadiness"].lower()
    criterion_row = audit.loc[audit["DesignRole"] == "Rubric criterion / latent dimension"].iloc[0]
    assert "multidimensional" in criterion_row["ClaimReadiness"].lower()
    notice = app.uto_design_audit_claim_notice(audit)
    assert notice["level"] == "info"
    assert "privacy/archive review" in notice["message"]
    wording = export["claim_wording"]
    assert {
        "DesignRole",
        "Status",
        "SafeReportWording",
        "DoNotWrite",
        "EvidenceToAttach",
        "NextAction",
        "ArchiveFile",
    }.issubset(wording.columns)
    assert len(wording) == len(audit)
    time_wording = wording.loc[wording["DesignRole"] == "Time/order block for rater drift"].iloc[0]
    assert time_wording["Status"] == "Ready"
    assert "documented block order" in time_wording["SafeReportWording"]
    assert "mfrm_uto_bayesian_mfrm_claim_wording.csv" in set(wording["ArchiveFile"])


@pytest.mark.legacy_compat
def test_uto_family_design_audit_blocks_unmapped_drift_and_multidimensional_claims():
    result = _stan_export_fixture_result()
    data = result["prep"]["data"][["Person", "Rater", "Task", "Score"]].copy()
    facets = ["Rater", "Task"]
    prep = app.prepare_mfrm_data(
        data,
        person_col="Person",
        facet_cols=facets,
        score_col="Score",
        rating_min=0,
        rating_max=4,
    )
    compact_result = {
        "config": {
            "model": "RSM",
            "method": "JMLE",
            "facet_names": facets,
            "n_cat": 5,
        },
        "prep": prep,
    }

    export = app.build_uto_bayesian_mfrm_stan_data_export(compact_result)

    assert export["available"]
    payload = json.loads(export["json_text"])
    assert payload["C_dim"] == 1
    assert payload["S"] == 1
    audit = export["design_audit"].set_index("DesignRole")
    criterion_row = audit.loc["Rubric criterion / latent dimension"]
    time_row = audit.loc["Time/order block for rater drift"]
    assert criterion_row["Status"] == "Review"
    assert "do not make multidimensional" in criterion_row["ClaimReadiness"].lower()
    assert time_row["Status"] == "Review"
    assert "arbitrary row order is not a time design" in time_row["ClaimReadiness"].lower()
    notice = app.uto_design_audit_claim_notice(export["design_audit"])
    assert notice["level"] == "warning"
    assert "Time/order block for rater drift" in notice["message"]
    assert "multidimensional rubric criteria" in notice["message"]
    wording = export["claim_wording"].set_index("DesignRole")
    criterion_wording = wording.loc["Rubric criterion / latent dimension"]
    time_wording = wording.loc["Time/order block for rater drift"]
    assert criterion_wording["Status"] == "Review"
    assert "do not write" not in criterion_wording["SafeReportWording"].lower()
    assert "do not write" not in time_wording["SafeReportWording"].lower()
    assert "multidimensional" in criterion_wording["DoNotWrite"].lower()
    assert "rater severity drift" in time_wording["DoNotWrite"].lower()


@pytest.mark.legacy_compat
def test_stan_facet_variable_names_are_valid_and_unique():
    names = app.stan_facet_variable_names(["Task Type", "Task-Type", "123", "data"])

    assert names == ["task_type", "task_type_2", "facet_3_123", "data_facet_4"]
    assert len(names) == len(set(names))
    assert all(name[0].isalpha() for name in names)


def test_final_readiness_uses_five_percent_residual_benchmark():
    class Opt:
        success = True
        message = "ok"

    result = {
        "opt": Opt(),
        "config": {},
        "prep": {},
        "facets": {"person": app.pd.DataFrame(), "others": app.pd.DataFrame()},
    }
    diagnostics_ready = {
        "obs": app.pd.DataFrame({"StdResidual": [2.1] * 5 + [0.0] * 95}),
        "reliability": app.pd.DataFrame(),
        "pca_enabled": False,
    }
    diagnostics_review = {
        "obs": app.pd.DataFrame({"StdResidual": [2.1] * 6 + [0.0] * 94}),
        "reliability": app.pd.DataFrame(),
        "pca_enabled": False,
    }

    ready = app.build_final_report_readiness(result, diagnostics_ready, all_bias_results={})
    review = app.build_final_report_readiness(result, diagnostics_review, all_bias_results={})
    ready_status = ready.loc[ready["Check"] == "Global residual fit", "Status"].iloc[0]
    review_status = review.loc[review["Check"] == "Global residual fit", "Status"].iloc[0]

    assert app.FINAL_RESIDUAL_PCT_GE2_READY == 5.0
    assert ready_status == "Ready"
    assert review_status == "Review"


def test_pca_bundle_carries_stability_audit_for_sparse_residual_matrix():
    mat = app.pd.DataFrame(
        {
            "A": [0.1, 0.2, app.np.nan, app.np.nan, app.np.nan, app.np.nan],
            "B": [0.0, 0.3, 0.4, app.np.nan, app.np.nan, app.np.nan],
            "C": [app.np.nan, 0.1, app.np.nan, 0.2, 0.3, app.np.nan],
        },
        index=[f"P{i}" for i in range(6)],
    )
    bundle = app.compute_pca_bundle(mat)

    assert bundle is not None
    assert "stability_table" in bundle
    stability = bundle["stability_table"]
    assert stability["PCAStabilityStatus"].iloc[0] == "Review"
    assert "missing" in stability["Caution"].iloc[0] or "overlap" in stability["Caution"].iloc[0]


def test_mml_prior_sensitivity_plan_and_assumption_audit_contract():
    result = {
        "config": {
            "method": "MML",
            "population_prior_sd": 1.2,
            "population_model": {"enabled": True},
        }
    }
    diagnostics = {
        "uncertainty": {
            "summary": app.pd.DataFrame(
                [{
                    "Area": "Measure SE/CI",
                    "Status": "ok",
                    "Rows": 3,
                    "Method": "MML observed-information delta method",
                    "Caution": "large-sample Wald CI",
                }]
            )
        },
        "reliability": app.pd.DataFrame(),
    }

    plan = app.build_mml_prior_sensitivity_plan(result)
    assert not plan.empty
    assert 1.2 in plan["PopulationPriorSD"].round(8).tolist()
    assert plan["Boundary"].str.contains("not estimated", case=False).any()
    assert {"PrimaryComparisons", "DecisionRule", "RunStatus"}.issubset(plan.columns)

    audit = app.build_statistical_assumption_audit(result, diagnostics, all_bias_results={})
    assert {"Measure SE/CI", "MML fixed population prior SD"}.issubset(set(audit["Area"]))
    prior_row = audit[audit["Area"] == "MML fixed population prior SD"].iloc[0]
    assert prior_row["Status"] == "Review"
    assert "fixed" in prior_row["Implication"].lower()


def test_manuscript_claim_guide_contract():
    class Opt:
        success = True
        message = "ok"

    result = {
        "opt": Opt(),
        "config": {
            "app_version": app.APP_VERSION,
            "model": "RSM",
            "method": "JMLE",
            "facet_names": ["Rater", "Task"],
            "anchor_audit": {"overall_status": "ok", "message": "No anchor issues."},
        },
        "prep": {"n_obs": 100, "n_person": 20, "rating_min": 1, "rating_max": 5},
        "facets": {
            "person": app.pd.DataFrame({"Estimate": [-0.5, 0.5]}),
            "others": app.pd.DataFrame({"Estimate": [-0.2, 0.2]}),
        },
    }
    diagnostics = {
        "obs": app.pd.DataFrame({
            "Observed": [1, 2, 3, 4, 5] * 20,
            "StdResidual": [2.1] * 5 + [0.0] * 95,
        }),
        "reliability": app.pd.DataFrame({"Facet": ["Person"], "Reliability": [0.85]}),
        "pca_enabled": False,
    }

    guide = app.build_manuscript_claim_guide(result, diagnostics, all_bias_results={})
    assert not guide.empty
    assert {
        "ManuscriptArea",
        "ClaimStatus",
        "SafeManuscriptWording",
        "EvidenceToReport",
        "DoNotClaim",
        "NextAction",
    }.issubset(guide.columns)
    assert "External package comparison" not in guide["ManuscriptArea"].tolist()
    assert "Bayesian / Uto-family extension" not in guide["ManuscriptArea"].tolist()
    bias = guide.loc[guide["ManuscriptArea"] == "Bias / local interaction"].iloc[0]
    assert bias["ClaimStatus"] == "Do not claim"
    fit = guide.loc[guide["ManuscriptArea"] == "Fit and dimensionality"].iloc[0]
    assert fit["ClaimStatus"] == "Report with caveat"

    template = app.generate_manuscript_reporting_template(result, diagnostics, all_bias_results={})
    assert "## Methods Draft" in template
    assert "## Results Draft" in template
    assert "standalone Python" in template
    assert "Claims Requiring Caution" in template
    assert "Keep every table and figure tied to the same AnalysisID" in template

    gate = app.build_publication_gate_summary(result, diagnostics, all_bias_results={})
    assert not gate.empty
    assert {"GateArea", "GateStatus", "Evidence", "ManuscriptAction"}.issubset(gate.columns)
    overall = gate.loc[gate["GateArea"] == "Overall manuscript gate"].iloc[0]
    assert overall["GateStatus"] == "Report with caveat"
    assert "Fit and dimensionality" in gate["GateArea"].astype(str).tolist()
    assert "Bias / local interaction" in gate["GateArea"].astype(str).tolist()

    cases = app.build_case_interpretation_guidance(result, diagnostics, all_bias_results={})
    assert not cases.empty
    assert {
        "Priority",
        "Case",
        "Status",
        "Evidence",
        "InterpretationNote",
        "ManuscriptGuardrail",
        "NextAction",
        "WhereToInspect",
        "AvoidWording",
        "SaferWording",
    }.issubset(cases.columns)
    assert "Dimensionality caveat" in cases["Case"].astype(str).tolist()
    assert "Bias/local interaction screen" in cases["Case"].astype(str).tolist()
    assert cases["AvoidWording"].astype(str).str.len().gt(0).all()
    assert cases["SaferWording"].astype(str).str.len().gt(0).all()

    action_plan = app.build_submission_action_plan(result, diagnostics, all_bias_results={})
    assert not action_plan.empty
    assert {
        "Priority",
        "Source",
        "Topic",
        "Status",
        "Evidence",
        "ImmediateAction",
        "BeforeCopying",
        "AvoidWording",
        "SaferWording",
        "WhereToInspect",
    }.issubset(action_plan.columns)
    assert action_plan["Priority"].tolist() == sorted(action_plan["Priority"].tolist())
    assert "Overall manuscript gate" in action_plan["Topic"].astype(str).tolist()
    assert "Do not claim" in action_plan["Status"].astype(str).tolist()

    claim_evidence = app.build_claim_to_evidence_matrix(result, diagnostics, all_bias_results={})
    assert not claim_evidence.empty
    assert {
        "Priority",
        "ManuscriptArea",
        "ManuscriptSection",
        "PrimaryClaim",
        "ClaimStatus",
        "GateStatus",
        "EvidenceToReport",
        "PrimaryTables",
        "PrimaryFigures",
        "ReadinessChecks",
        "ReviewerQuestion",
        "ArchiveFiles",
        "AppLocation",
    }.issubset(claim_evidence.columns)
    assert "Bias / local interaction" in claim_evidence["ManuscriptArea"].astype(str).tolist()
    assert claim_evidence["ArchiveFiles"].astype(str).str.len().gt(0).all()

    report_ready = app.build_report_ready_summary_panel(result, diagnostics, all_bias_results={})
    assert not report_ready.empty
    assert {
        "Priority",
        "ReportArea",
        "ReportStatus",
        "PrimaryUsers",
        "UseNow",
        "OpenFirst",
        "EvidenceFiles",
        "ActionBeforeWriting",
        "SafeWording",
        "DoNotWrite",
        "Details",
    }.issubset(report_ready.columns)
    assert {
        "Overall manuscript conclusion",
        "APA sentence draft",
        "Claim boundary and evidence",
        "Scoring quality and QC decision",
        "Rubric categories and collapse sensitivity",
    }.issubset(set(report_ready["ReportArea"].astype(str)))
    assert set(report_ready["ReportStatus"].astype(str)).issubset(
        {"Ready", "Needs caveat", "Needs rerun", "Do not report yet"}
    )
    joined_report_ready = " ".join(report_ready.astype(str).to_numpy().ravel().tolist())
    assert "publication_gate_summary.csv" in joined_report_ready
    assert "apa_report_sentence_audit.csv" in joined_report_ready
    assert "claim_to_evidence_matrix.csv" in joined_report_ready
    assert "category_collapse_action_center.csv" in joined_report_ready

    decision_brief = app.generate_report_ready_decision_brief(result, diagnostics, all_bias_results={})
    assert "# Report-Ready Decision Brief" in decision_brief
    assert "Overall status:" in decision_brief
    assert "publication_gate_summary.csv" in decision_brief
    assert "report_ready_summary_panel.csv" in decision_brief or "apa_report_sentence_audit.csv" in decision_brief

    apa_results_draft = app.generate_report_ready_apa_results_draft(result, diagnostics, all_bias_results={})
    assert "# APA Results Paragraph Draft" in apa_results_draft
    assert "Overall report-ready status:" in apa_results_draft
    assert "## Final-Output Blockers" in apa_results_draft
    assert "## Critical Review Caveats" in apa_results_draft
    assert "## Bridge-Constrained Use Conditions" in apa_results_draft
    assert "## Draft Paragraph" in apa_results_draft
    assert "## Bridge Revision Plan" in apa_results_draft
    assert "## Do Not Use Yet" in apa_results_draft
    assert "apa_report_sentence_audit.csv" in apa_results_draft or "claim_to_evidence_matrix.csv" in apa_results_draft
    assert "critical_final_review_panel.csv" in apa_results_draft
    assert "reporting_action_bridge.csv" in apa_results_draft

    role_memos = app.build_role_based_action_memos(result, diagnostics, all_bias_results={})
    assert not role_memos.empty
    assert {
        "Role",
        "Priority",
        "DecisionFocus",
        "ReportStatus",
        "DecisionQuestion",
        "WhyThisMatters",
        "ImmediateAction",
        "EvidenceToOpen",
        "SafeOutput",
        "DoNotDo",
        "AppLocation",
        "DownloadFile",
        "Details",
    }.issubset(role_memos.columns)
    assert {
        "Paper author",
        "Scoring manager",
        "Measurement practitioner",
    }.issubset(set(role_memos["Role"].astype(str)))
    assert set(role_memos["ReportStatus"].astype(str)).issubset(
        {"Ready", "Needs caveat", "Needs rerun", "Do not report yet"}
    )
    joined_role_memos = " ".join(role_memos.astype(str).to_numpy().ravel().tolist())
    assert "apa_results_paragraph_draft.md" in joined_role_memos
    assert "quality_control_todo_checklist.csv" in joined_role_memos
    assert "final_report_readiness.csv" in joined_role_memos

    role_memos_md = app.generate_role_based_action_memos_markdown(result, diagnostics, all_bias_results={})
    assert "# Role-Based Action Memos" in role_memos_md
    assert "## Paper author" in role_memos_md
    assert "## Scoring manager" in role_memos_md
    assert "## Measurement practitioner" in role_memos_md
    assert "role_based_action_memos.csv" in role_memos_md or "report_ready_summary_panel.csv" in role_memos_md

    role_checklist = app.build_role_based_action_checklist(result, diagnostics, all_bias_results={})
    assert not role_checklist.empty
    assert {
        "Done",
        "ChecklistId",
        "Role",
        "Step",
        "ToDoItem",
        "ReportStatus",
        "BlockingForFinalOutput",
        "NeedsCaveat",
        "Action",
        "EvidenceToOpen",
        "CompletionCriterion",
        "OutputArtifact",
        "SafeOutput",
        "DoNotDo",
        "AppLocation",
        "DownloadFile",
    }.issubset(role_checklist.columns)
    assert role_checklist["Done"].eq(False).all()
    assert {
        "Paper author",
        "Scoring manager",
        "Measurement practitioner",
    }.issubset(set(role_checklist["Role"].astype(str)))
    assert role_checklist["CompletionCriterion"].astype(str).str.len().gt(0).all()
    joined_role_checklist = " ".join(role_checklist.astype(str).to_numpy().ravel().tolist())
    assert "apa_results_paragraph_draft.md" in joined_role_checklist
    assert "quality_control_todo_checklist.csv" in joined_role_checklist

    role_checklist_md = app.generate_role_based_action_checklist_markdown(result, diagnostics, all_bias_results={})
    assert "# Role-Based Action Checklist" in role_checklist_md
    assert "- [ ]" in role_checklist_md
    assert "role_based_action_checklist.csv" in role_checklist_md

    reanalysis_checklist = app.build_report_ready_reanalysis_checklist(
        result,
        diagnostics,
        all_bias_results={},
        report_ready_summary=report_ready,
        role_action_checklist=role_checklist,
    )
    assert not reanalysis_checklist.empty
    assert {
        "Done",
        "ChecklistId",
        "Priority",
        "Phase",
        "Audience",
        "Task",
        "ReportStatus",
        "BlocksReportReady",
        "NeedsCaveat",
        "Action",
        "EvidenceToOpen",
        "CompletionCriterion",
        "OutputArtifact",
        "CaveatToCarry",
        "DoNotDo",
        "AppLocation",
        "DownloadFile",
    }.issubset(reanalysis_checklist.columns)
    assert reanalysis_checklist["Done"].eq(False).all()
    assert reanalysis_checklist["ChecklistId"].astype(str).str.len().gt(0).all()
    assert {
        "Study result handoff",
        "Reanalysis reproducibility",
        "Model diagnostic rerun",
        "Scoring quality control",
        "Rubric/category sensitivity",
        "APA/manuscript wording",
    }.issubset(set(reanalysis_checklist["Phase"].astype(str)))
    joined_reanalysis = " ".join(reanalysis_checklist.astype(str).to_numpy().ravel().tolist())
    assert "report_ready_summary_panel.csv" in joined_reanalysis
    assert "MFRM_Manuscript_Binder.zip" in joined_reanalysis

    reanalysis_md = app.generate_report_ready_reanalysis_checklist_markdown(reanalysis_checklist)
    assert "# Report-Ready Simulation / Reanalysis Checklist" in reanalysis_md
    assert "- [ ]" in reanalysis_md
    assert "report_ready_reanalysis_checklist.csv" in reanalysis_md

    first_checklist_id = str(role_checklist.iloc[0]["ChecklistId"])
    persisted_role_checklist = app.apply_role_based_action_checklist_done_state(
        role_checklist,
        {first_checklist_id: True},
    )
    assert persisted_role_checklist.loc[
        persisted_role_checklist["ChecklistId"].astype(str).eq(first_checklist_id),
        "Done",
    ].eq(True).all()
    role_checklist_progress = app.role_based_action_checklist_progress_summary(persisted_role_checklist)
    assert not role_checklist_progress.empty
    assert {
        "Role",
        "DoneItems",
        "TotalItems",
        "RemainingItems",
        "RemainingBlockingItems",
        "RemainingCaveatItems",
        "CompletionPercent",
        "State",
        "NextAction",
        "NextEvidenceToOpen",
    }.issubset(role_checklist_progress.columns)
    assert role_checklist_progress["TotalItems"].sum() == len(role_checklist)
    assert role_checklist_progress["DoneItems"].sum() == 1
    role_checklist_done_md = app.generate_role_based_action_checklist_markdown(
        result,
        diagnostics,
        all_bias_results={},
        checklist=persisted_role_checklist,
    )
    assert "- [x]" in role_checklist_done_md
    assert "Progress Snapshot" in role_checklist_done_md

    critical_review = app.build_critical_final_review_panel(
        result,
        diagnostics,
        all_bias_results={},
        role_action_checklist=persisted_role_checklist,
    )
    assert not critical_review.empty
    assert {
        "Priority",
        "Lens",
        "OwnerRole",
        "CriticalQuestion",
        "CurrentStatus",
        "Decision",
        "StopBeforeFinalOutput",
        "EvidenceBlocksDirectClaim",
        "UnresolvedWorkflowItems",
        "UnresolvedWorkflowBlockers",
        "NeedsCaveat",
        "WhatCouldGoWrong",
        "ActionBeforeContinuing",
        "EvidenceToOpen",
        "DoNotClaim",
        "SafeNextOutput",
        "DownloadFile",
    }.issubset(critical_review.columns)
    assert {
        "Manuscript conclusion",
        "APA wording",
        "Claim scope",
        "Scoring operations",
        "Rubric/category decision",
        "Measurement diagnostics",
        "Reproducibility and archive",
        "Sharing and privacy",
    }.issubset(set(critical_review["Lens"].astype(str)))
    assert critical_review["WhatCouldGoWrong"].astype(str).str.len().gt(0).all()
    assert critical_review["EvidenceToOpen"].astype(str).str.contains(".csv|.zip|.json|.md", regex=True).any()
    critical_review_md = app.generate_critical_final_review_markdown(
        result,
        diagnostics,
        all_bias_results={},
        role_action_checklist=persisted_role_checklist,
    )
    assert "# Critical Final Review" in critical_review_md
    assert "Stop-before-final-output rows" in critical_review_md
    assert "critical_final_review_panel.csv" in critical_review_md

    status_rationale = app.build_status_rationale_drilldown(
        result,
        diagnostics,
        all_bias_results={},
        role_action_checklist=persisted_role_checklist,
    )
    assert not status_rationale.empty
    assert {
        "Priority",
        "Audience",
        "StatusSource",
        "CurrentStatus",
        "Decision",
        "DecisionLogic",
        "TriggerDiagnostics",
        "MinimumEvidenceNeeded",
        "ClaimBoundary",
        "NextInspection",
        "SafeUse",
        "DownloadFile",
    }.issubset(status_rationale.columns)
    assert status_rationale["Audience"].astype(str).str.contains("Psychometric").any()
    assert status_rationale["DecisionLogic"].astype(str).str.contains("stop_before_final_output").all()

    reading_path = app.build_result_reading_path(result, diagnostics, all_bias_results={})
    assert not reading_path.empty
    assert {
        "Step",
        "Audience",
        "ReadingGoal",
        "GuidingQuestion",
        "CurrentStatus",
        "OpenFirst",
        "EvidenceFile",
        "HowToRead",
        "WhatCanBeWritten",
        "DoNotWriteYet",
        "NextMove",
    }.issubset(reading_path.columns)
    assert reading_path["Step"].tolist() == sorted(reading_path["Step"].tolist())
    assert "Draft only after evidence mapping" in reading_path["ReadingGoal"].astype(str).tolist()
    assert reading_path["EvidenceFile"].astype(str).str.contains("critical_final_review_panel.csv").any()

    operational_board = app.build_operational_decision_board(result, diagnostics, all_bias_results={})
    assert not operational_board.empty
    assert {
        "Priority",
        "Audience",
        "Decision",
        "OwnerRole",
        "CurrentSignal",
        "EvidenceToOpen",
        "RecommendedAction",
        "RiskIfIgnored",
        "BeforeNextAdministration",
        "DoNotDo",
    }.issubset(operational_board.columns)
    assert operational_board["Decision"].astype(str).str.contains("Rater training|Rubric revision|Category", regex=True).any()
    assert operational_board["EvidenceToOpen"].astype(str).str.contains(".csv", regex=False).any()

    reporting_bridge = app.build_reporting_action_bridge(
        result,
        diagnostics,
        all_bias_results={},
        role_action_checklist=persisted_role_checklist,
        critical_final_review=critical_review,
        status_rationale=status_rationale,
        result_reading_path=reading_path,
        operational_decision_board=operational_board,
    )
    assert not reporting_bridge.empty
    assert {
        "Priority",
        "BridgeTrack",
        "PrimaryAudience",
        "SourceBoard",
        "CurrentStatus",
        "OutputToProduce",
        "DraftOrDecisionMove",
        "RequiredCheck",
        "CaveatToCarry",
        "StopRule",
        "EvidenceFiles",
        "AppLocation",
        "DownloadFile",
    }.issubset(reporting_bridge.columns)
    joined_bridge = " ".join(reporting_bridge.astype(str).to_numpy().ravel().tolist())
    assert "APA final wording gate" in joined_bridge
    assert "Rater training action memo" in joined_bridge
    assert "Rubric revision action memo" in joined_bridge
    assert "Category collapse sensitivity memo" in joined_bridge
    assert "reporting_action_bridge.csv" in joined_bridge
    assert "apa_results_paragraph_draft.md" in joined_bridge
    reporting_bridge_md = app.generate_reporting_action_bridge_markdown(
        result,
        diagnostics,
        all_bias_results={},
        bridge=reporting_bridge,
    )
    assert "# Reporting Action Bridge" in reporting_bridge_md
    assert "APA wording" in reporting_bridge_md
    assert "operational_decision_board.csv" in reporting_bridge_md

    apa_bridge_plan = app.build_apa_bridge_revision_plan(
        result,
        diagnostics,
        all_bias_results={},
        critical_review=critical_review,
        reporting_action_bridge=reporting_bridge,
    )
    assert not apa_bridge_plan.empty
    assert {
        "Priority",
        "APAPlacement",
        "SourceTrack",
        "CurrentStatus",
        "BlocksFinalCopy",
        "RevisionInstruction",
        "RequiredCheck",
        "CaveatToCarry",
        "StopRule",
        "EvidenceFiles",
        "DownloadFile",
    }.issubset(apa_bridge_plan.columns)
    assert apa_bridge_plan["SourceTrack"].astype(str).str.contains("APA|Claim|Guided", regex=True).any()
    assert apa_bridge_plan["DownloadFile"].astype(str).str.contains("reporting_action_bridge.csv").any()

    apa_results_draft_with_bridge = app.generate_report_ready_apa_results_draft(
        result,
        diagnostics,
        all_bias_results={},
        critical_review=critical_review,
        reporting_action_bridge=reporting_bridge,
        apa_bridge_revision_plan=apa_bridge_plan,
    )
    assert "Bridge-Constrained Use Conditions" in apa_results_draft_with_bridge
    assert "Bridge Revision Plan" in apa_results_draft_with_bridge
    assert "Caveat-carrying draft" in apa_results_draft_with_bridge or "Editing worksheet only" in apa_results_draft_with_bridge

    handoff = app.build_manuscript_handoff_checklist(
        result,
        diagnostics,
        all_bias_results={},
        public_export_mode=True,
    )
    assert not handoff.empty
    assert {
        "Step",
        "Phase",
        "Status",
        "Task",
        "AppLocation",
        "DownloadFile",
        "Action",
        "WhyItMatters",
    }.issubset(handoff.columns)
    assert handoff["Step"].tolist() == sorted(handoff["Step"].tolist())
    assert "MFRM_OSF_Package.zip; mfrm_manuscript_handoff.md" in handoff["DownloadFile"].astype(str).tolist()
    assert handoff["DownloadFile"].astype(str).str.contains("report_ready_decision_brief.md").any()
    assert handoff["DownloadFile"].astype(str).str.contains("role_based_action_memos.md").any()
    assert handoff["DownloadFile"].astype(str).str.contains("role_based_action_checklist.csv").any()
    assert handoff["DownloadFile"].astype(str).str.contains("critical_final_review_panel.csv").any()
    assert handoff["DownloadFile"].astype(str).str.contains("status_rationale_drilldown.csv").any()
    assert handoff["DownloadFile"].astype(str).str.contains("result_reading_path.csv").any()
    assert handoff["DownloadFile"].astype(str).str.contains("operational_decision_board.csv").any()
    assert handoff["DownloadFile"].astype(str).str.contains("reporting_action_bridge.csv").any()
    assert handoff["DownloadFile"].astype(str).str.contains("reporting_action_bridge.md").any()

    handoff_md = app.generate_manuscript_handoff_markdown(
        result,
        diagnostics,
        all_bias_results={},
        public_export_mode=True,
    )
    assert "Final Results and Manuscript Handoff" in handoff_md
    assert "## Download Package" in handoff_md
    assert "## Before Manuscript Submission" in handoff_md
    assert "claim_to_evidence_matrix.csv" in handoff_md
    assert "report_ready_decision_brief.md" in handoff_md
    assert "role_based_action_memos.md" in handoff_md
    assert "role_based_action_checklist.csv" in handoff_md
    assert "critical_final_review_panel.csv" in handoff_md
    assert "status_rationale_drilldown.csv" in handoff_md
    assert "result_reading_path.csv" in handoff_md
    assert "operational_decision_board.csv" in handoff_md
    assert "reporting_action_bridge.csv" in handoff_md
    assert "reporting_action_bridge.md" in handoff_md
    assert "apa_results_paragraph_draft.md" in handoff_md
    assert "MFRM_OSF_Package.zip" in handoff_md

    binder_assets = app.build_manuscript_binder_assets(
        {
            "report_ready_summary_panel": report_ready,
            "role_based_action_memos": role_memos,
            "role_based_action_checklist": role_checklist,
            "report_ready_reanalysis_checklist": reanalysis_checklist,
            "critical_final_review_panel": critical_review,
            "status_rationale_drilldown": status_rationale,
            "result_reading_path": reading_path,
            "operational_decision_board": operational_board,
            "reporting_action_bridge": reporting_bridge,
            "claim_to_evidence_matrix": claim_evidence,
            "submission_action_plan": action_plan,
            "manuscript_handoff_checklist": handoff,
        },
        {
            "manuscript_handoff.md": handoff_md,
            "report_ready_decision_brief.md": decision_brief,
            "apa_results_paragraph_draft.md": apa_results_draft,
            "role_based_action_memos.md": role_memos_md,
            "role_based_action_checklist.md": role_checklist_md,
            "report_ready_reanalysis_checklist.md": reanalysis_md,
            "critical_final_review.md": critical_review_md,
            "reporting_action_bridge.md": reporting_bridge_md,
            "method_appendix.md": "method",
            "manuscript_template.md": template,
            "mfrm_app_engine_runner.py": "print('runner')",
            "requirements.txt": "pandas\n",
            "mfrm_local_batch_workflow.md": "# Local batch",
            "mfrm_jmle_self_contained.py": "print('jmle')",
        },
        public_export_mode=True,
    )
    assert "README_first.md" in binder_assets
    assert "report_ready_summary_panel.csv" in binder_assets
    assert "report_ready_decision_brief.md" in binder_assets
    assert "apa_results_paragraph_draft.md" in binder_assets
    assert "role_based_action_memos.csv" in binder_assets
    assert "role_based_action_memos.md" in binder_assets
    assert "role_based_action_checklist.csv" in binder_assets
    assert "role_based_action_checklist.md" in binder_assets
    assert "report_ready_reanalysis_checklist.csv" in binder_assets
    assert "report_ready_reanalysis_checklist.md" in binder_assets
    assert "critical_final_review_panel.csv" in binder_assets
    assert "critical_final_review.md" in binder_assets
    assert "status_rationale_drilldown.csv" in binder_assets
    assert "result_reading_path.csv" in binder_assets
    assert "operational_decision_board.csv" in binder_assets
    assert "reporting_action_bridge.csv" in binder_assets
    assert "reporting_action_bridge.md" in binder_assets
    assert "claim_to_evidence_matrix.csv" in binder_assets
    assert "manuscript_handoff.md" in binder_assets
    assert "mfrm_app_engine_runner.py" in binder_assets
    assert "requirements.txt" in binder_assets
    assert "mfrm_local_batch_workflow.md" in binder_assets
    assert "mfrm_jmle_self_contained.py" in binder_assets
    assert "stan_reproducibility_archive_contract.csv" not in binder_assets
    assert "stan_posterior_reproducibility_handoff.md" not in binder_assets
    assert "stan_posterior_reproducibility_route.csv" not in binder_assets
    assert "stan_posterior_handoff_checklist.csv" not in binder_assets
    assert "visualization_settings.json" in binder_assets
    assert "visualization_settings.csv" in binder_assets
    binder_zip = app.cached_mixed_asset_zip(
        binder_assets,
        app.bytes_mapping_fingerprint(binder_assets),
    )
    with zipfile.ZipFile(io.BytesIO(binder_zip), "r") as zf:
        binder_zip_names = set(zf.namelist())
    assert "README_first.md" in binder_zip_names
    assert "mfrm_app_engine_runner.py" in binder_zip_names
    assert "mfrm_jmle_self_contained.py" in binder_zip_names
    assert "stan_reproducibility_archive_contract.csv" not in binder_zip_names
    assert "stan_posterior_reproducibility_handoff.md" not in binder_zip_names

    osf_zip = app.build_osf_zip(
        {"claim_to_evidence_matrix": claim_evidence},
        text_assets={
            "manuscript_handoff.md": handoff_md,
            "mfrm_app_engine_runner.py": "print('runner')",
            "requirements.txt": "pandas\n",
        },
    )
    with zipfile.ZipFile(io.BytesIO(osf_zip), "r") as zf:
        osf_namelist = zf.namelist()
        osf_zip_names = set(osf_namelist)
    assert len(osf_namelist) == len(osf_zip_names)
    assert "claim_to_evidence_matrix.csv" in osf_zip_names
    assert "manuscript_handoff.md" in osf_zip_names
    assert "mfrm_app_engine_runner.py" in osf_zip_names
    assert "stan_reproducibility_archive_contract.csv" not in osf_zip_names
    assert "stan_posterior_reproducibility_handoff.md" not in osf_zip_names
    assert "stan_posterior_reproducibility_route.csv" not in osf_zip_names
    assert "mfrm_stan_data.json" not in osf_zip_names
    assert "mfrm_stan_id_index_map.csv" not in osf_zip_names

    settings_table = app.build_visualization_preferences_table()
    assert {"Setting", "Value", "Affects"}.issubset(settings_table.columns)

    figure_manifest = app._figure_export_manifest(
        {"wright_map": "<html></html>", "fit_scatter": "<html></html>"},
        {"wright_map": b"png"},
    )
    assert {
        "Theme",
        "BaseFontSize",
        "LabelPolicy",
        "PlotLabelMode",
        "LabelMaxChars",
        "CaptionDetail",
    }.issubset(figure_manifest.columns)
    visual_map = app.build_visual_evidence_map(figure_manifest, claim_evidence)
    assert not visual_map.empty
    assert {
        "FigureName",
        "LinkedManuscriptArea",
        "ManuscriptSection",
        "Theme",
        "BaseFontSize",
        "LabelPolicy",
        "PlotLabelMode",
        "CaptionDetail",
        "CaptionDraft",
        "EvidenceTables",
        "ReviewerQuestion",
        "ReportingCaution",
        "ArchiveFiles",
        "VisualGuardrailID",
        "SafeReportWording",
        "DoNotWrite",
        "RequiredEvidence",
        "VisualGuardrailArchive",
    }.issubset(visual_map.columns)
    assert "Wright map targeting and measure interpretation" in visual_map["LinkedManuscriptArea"].astype(str).tolist()
    assert visual_map["DoNotWrite"].astype(str).str.contains("alone|solely", case=False, regex=True).any()
    captions = app.generate_visual_caption_drafts(visual_map)
    assert "Figure 1" in captions
    assert "Reviewer question" in captions
    assert "Do not write" in captions
    visual_qa = app.build_visual_qa_preflight(figure_manifest, visual_map)
    assert not visual_qa.empty
    assert {"Check", "Status", "Evidence", "Action"}.issubset(visual_qa.columns)

    visual_assets = app.build_visual_evidence_binder_assets(
        {"wright_map": "<html></html>"},
        {"wright_map": b"png"},
        figure_manifest,
        visual_map,
        claim_evidence,
        public_export_mode=True,
    )
    assert "README_visual_evidence.md" in visual_assets
    assert "visual_evidence_map.csv" in visual_assets
    assert "visual_caption_drafts.md" in visual_assets
    assert "visualization_settings.json" in visual_assets
    assert "visualization_settings.csv" in visual_assets
    assert "visual_qa_preflight.csv" in visual_assets
    assert "visual_claim_guardrails.csv" in visual_assets
    assert "figures_png_wright_map.png" in visual_assets


def test_uploaded_file_fingerprint_uses_content_not_only_name_and_size():
    first = io.BytesIO(b"abc")
    second = io.BytesIO(b"abd")
    first.name = "ratings.csv"
    second.name = "ratings.csv"
    first.size = 3
    second.size = 3

    assert app.uploaded_file_fingerprint(first) != app.uploaded_file_fingerprint(second)


def test_population_covariate_type_summary_flags_integer_codes():
    person_data = app.pd.DataFrame({
        "Person": ["P1", "P2", "P3", "P4"],
        "GradeCode": [1, 2, 1, 2],
        "SES": [0.1, -0.4, 0.6, 0.2],
        "Group": ["A", "A", "B", "B"],
    })
    summary = app.summarize_population_covariate_types(
        person_data,
        "Person",
        "~ GradeCode + SES + Group",
        person_levels=["P1", "P2", "P3", "P4"],
    )
    grade_code = summary.loc[summary["Term"] == "GradeCode"].iloc[0]
    ses = summary.loc[summary["Term"] == "SES"].iloc[0]

    assert grade_code["InferredType"] == "numeric"
    assert bool(grade_code["ReviewFlag"])
    assert ses["InferredType"] == "numeric"
    assert not bool(ses["ReviewFlag"])

    forced = app.summarize_population_covariate_types(
        person_data,
        "Person",
        "~ GradeCode + SES + Group",
        categorical_terms=["GradeCode"],
        person_levels=["P1", "P2", "P3", "P4"],
    )
    forced_grade_code = forced.loc[forced["Term"] == "GradeCode"].iloc[0]
    assert forced_grade_code["InferredType"] == "categorical"
    assert forced_grade_code["Override"] == "categorical"
    assert not bool(forced_grade_code["ReviewFlag"])


def test_public_beta_release_contract_tables():
    limitations = app.standalone_release_limitations_table()
    readiness = app.public_release_readiness_table()

    assert not limitations.empty
    assert {"Area", "PublicBetaStatus", "SupportedNow", "Boundary", "UserAction"}.issubset(limitations.columns)
    joined = " ".join(limitations["Area"].astype(str).tolist() + limitations["Boundary"].astype(str).tolist())
    assert "GPCM" in joined
    assert "privacy" in joined
    assert "Cross-package" not in joined
    assert "external estimation-engine" in joined

    assert not readiness.empty
    assert {"Check", "Status", "Evidence", "Action"}.issubset(readiness.columns)
    assert "License" in readiness["Check"].astype(str).tolist()
    assert not (readiness["Status"].astype(str) == "Blocker").any()
