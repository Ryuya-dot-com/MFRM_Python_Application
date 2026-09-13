"""Contract tests for the cross-engine validation bundle.

`build_cross_engine_validation_bundle` ships the exact fitted data, the run
settings, the app's own estimates, an R refit script (TAM/sirt/mirt), and a
README of expected differences, so a researcher can cross-check the main
MML/GPCM estimates against established R packages. These tests pin the bundle's
contract: the right files, the right CSV columns (no internal `score_k` leak),
the R tokens, the equivalence-boundary wording, GPCM's extra slope file, and the
empty-result guard. They do not assert numerical parity — that is, by design,
not claimed.
"""

from __future__ import annotations

import io

import numpy as np
import pandas as pd

import streamlit_app as app


def _sim_rsm(seed: int = 4242, n_person: int = 80):
    params = {
        "persons": [f"P{i:03d}" for i in range(n_person)],
        "raters": ["R1", "R2", "R3"],
        "tasks": ["T1", "T2"],
        "criteria": ["C1", "C2"],
        "theta_sd": 1.3,
        "rater_severities": np.array([-0.3, 0.0, 0.3]),
        "task_difficulties": np.array([-0.2, 0.2]),
        "criterion_difficulties": np.array([-0.1, 0.1]),
        "tau": np.array([-1.0, 0.0, 1.0]),
    }
    return app._generate_mfrm_rsm_from_params(params, seed=seed)


def _fit(df, model="RSM", **kw):
    common = dict(
        data=df, person_col="Person", facet_cols=["Rater", "Task", "Criterion"],
        score_col="Score", model=model, method="MML", mml_engine="EM",
        maxit=150, reltol=1e-6, quad_points=13, population_prior_sd=1.0,
    )
    common.update(kw)
    return app.mfrm_estimate(**common)


def _read_csv(byts: bytes) -> pd.DataFrame:
    return pd.read_csv(io.BytesIO(byts))


def test_rsm_bundle_has_expected_files_and_columns():
    res = _fit(_sim_rsm(), model="RSM")
    bundle = app.build_cross_engine_validation_bundle(res)
    for key in ("data.csv", "settings.csv", "app_person.csv", "app_facets.csv",
                "app_steps.csv", "run_tam_mirt_crosscheck.R", "README_cross_engine.md"):
        assert key in bundle, f"missing {key}"
    # RSM has no per-item slope facet, so no slopes file.
    assert "app_slopes.csv" not in bundle

    data = _read_csv(bundle["data.csv"])
    assert "Person" in data.columns and "Score" in data.columns
    for facet in ("Rater", "Task", "Criterion"):
        assert facet in data.columns
    assert "score_k" not in data.columns, "internal score_k column leaked into data.csv"


def test_settings_records_model_and_identification():
    res = _fit(_sim_rsm(), model="RSM")
    settings = _read_csv(app.build_cross_engine_validation_bundle(res)["settings.csv"])
    row = settings.iloc[0]
    assert row["model"] == "RSM"
    assert row["method"] == "MML"
    assert str(row["facet_names"]) == "Rater;Task;Criterion"
    assert "noncenter_facet" in settings.columns
    assert "population_prior_sd" in settings.columns
    assert "quad_points" in settings.columns


def test_gpcm_bundle_adds_slopes_file():
    res = _fit(_sim_rsm(seed=77), model="GPCM", step_facet="Criterion", slope_facet="Criterion")
    bundle = app.build_cross_engine_validation_bundle(res)
    assert "app_slopes.csv" in bundle
    slopes = _read_csv(bundle["app_slopes.csv"])
    assert not slopes.empty
    settings = _read_csv(bundle["settings.csv"])
    assert settings.iloc[0]["model"] == "GPCM"
    assert settings.iloc[0]["step_facet"] == "Criterion"
    assert settings.iloc[0]["slope_facet"] == "Criterion"


def test_r_script_targets_three_engines():
    res = _fit(_sim_rsm(), model="RSM")
    r_txt = app.build_cross_engine_validation_bundle(res)["run_tam_mirt_crosscheck.R"]
    for token in ("tam.mml.mfr", "rm.facets", "mirt", "read.csv", "settings.csv", "data.csv"):
        assert token in r_txt, f"R script missing {token}"


def test_readme_states_equivalence_boundary():
    res = _fit(_sim_rsm(), model="RSM")
    readme = app.build_cross_engine_validation_bundle(res)["README_cross_engine.md"]
    low = readme.lower()
    assert "not claimed" in low
    assert "rank" in low
    assert "expected" in low
    # Names the three packages so the artifact is self-describing.
    for pkg in ("TAM", "sirt", "mirt"):
        assert pkg in readme


def test_bundle_round_trips_through_zip():
    res = _fit(_sim_rsm(), model="RSM")
    bundle = app.build_cross_engine_validation_bundle(res)
    zbytes = app.cached_mixed_asset_zip(bundle, "test_cross_engine")
    assert isinstance(zbytes, (bytes, bytearray)) and len(zbytes) > 0
    import zipfile
    with zipfile.ZipFile(io.BytesIO(zbytes)) as zf:
        names = set(zf.namelist())
    assert {"data.csv", "settings.csv", "run_tam_mirt_crosscheck.R",
            "README_cross_engine.md"}.issubset(names)


def test_empty_or_invalid_result_yields_empty_bundle():
    assert app.build_cross_engine_validation_bundle({}) == {}
    assert app.build_cross_engine_validation_bundle(None) == {}
    assert app.build_cross_engine_validation_bundle({"prep": {"data": pd.DataFrame()}}) == {}


def test_bundle_registered_in_export_matrix():
    matrix = app.reproducibility_script_export_matrix()
    assert matrix["Artifact"].astype(str).str.contains("Cross_Engine_Validation").any()
