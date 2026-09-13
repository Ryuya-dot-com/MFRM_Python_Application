"""Negative boundaries for the standalone Python release-check surface."""

from __future__ import annotations

import json

import streamlit_app as app


EXCLUDED_RELEASE_AREAS = {
    "Cross-package equivalence",
    "Advanced models (Stan, download-only)",
    "Posterior Viewer (upload)",
}

EXCLUDED_READINESS_CHECKS = {
    "External validation stance",
    "Archived validation artifact inventory",
    "Simulation validation templates",
    "Cross-language reproducibility scripts",
    "Uto-family Bayesian MFRM Stan roadmap",
    "mfrmr 0.2.0 migration coverage",
}


def test_standalone_release_limitations_are_python_native():
    limitations = app.standalone_release_limitations_table()
    areas = set(limitations["Area"].astype(str))

    assert areas.isdisjoint(EXCLUDED_RELEASE_AREAS)
    portable = limitations.loc[limitations["Area"].eq("Portable scripts")]
    assert len(portable) == 1
    wording = " ".join(
        portable.iloc[0][column]
        for column in ("SupportedNow", "Boundary", "UserAction")
    )
    assert "Python" in wording
    assert "JMLE/R workflows" not in wording


def test_public_release_readiness_excludes_external_handoff_checks():
    readiness = app.public_release_readiness_table()
    checks = set(readiness["Check"].astype(str))

    assert checks.isdisjoint(EXCLUDED_READINESS_CHECKS)


def test_release_check_json_has_exact_native_payload(capsys, monkeypatch):
    def fail_if_called(*_args, **_kwargs):
        raise AssertionError("legacy external release generator was called")

    for function_name in (
        "external_simulation_reference_inventory",
        "external_simulation_template_inventory",
        "reproducibility_script_export_matrix",
        "bayesian_mfrm_stan_refinement_plan",
        "mfrmr_015_migration_coverage_table",
        "mfrmr_016_migration_coverage_table",
        "mfrmr_020_migration_coverage_table",
    ):
        monkeypatch.setattr(app, function_name, fail_if_called)

    assert app.run_release_check(json_output=True) == 0
    payload = json.loads(capsys.readouterr().out)

    assert set(payload) == {
        "release_status",
        "app_version",
        "release_label",
        "readiness",
        "limitations",
    }
    assert {
        row["Check"] for row in payload["readiness"]
    }.isdisjoint(EXCLUDED_READINESS_CHECKS)
    assert {
        row["Area"] for row in payload["limitations"]
    }.isdisjoint(EXCLUDED_RELEASE_AREAS)
