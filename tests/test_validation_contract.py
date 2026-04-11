import streamlit_app as app


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
