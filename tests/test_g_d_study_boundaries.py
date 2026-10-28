"""Independent ANOVA controls, rejected designs and public G/D output limits."""
import io
import itertools

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app


def ratings(two_facets=True):
    rows = []
    for p, r, i in itertools.product((-1, 1), (-1, 1), (-1, 1) if two_facets else (0,)):
        score = (20 + 4*p + 3*r + 2*i + 1.5*p*r + p*i + .75*r*i + .25*p*r*i
                 if two_facets else 20 + 4*p + 2*r + p*r)
        row = {"Person": str(p), "Rater": str(r), "Score": score}
        if two_facets:
            row["Item"] = str(i)
        rows.append(row)
    return pd.DataFrame(rows)


def result(data):
    return {"prep": {"data": data}, "config": {
        "facet_names": [c for c in data.columns if c not in ("Score", "Weight")]}}


@pytest.mark.parametrize("two_facets", [False, True])
def test_independent_orthogonal_contrast_reference(two_facets):
    """Expected components follow hand-expanded orthogonal contrasts, not app output."""
    data = ratings(two_facets)
    out = app.compute_generalizability_study(result(data))
    assert out["available"]
    vc = out["variance_components"].set_index("Source")
    expected = ({"Person": 25.625, "Rater": 12.5, "Item": 5, "Person:Rater": 8.75,
                 "Person:Item": 3.75, "Rater:Item": 2, "Residual": .5}
                if two_facets else {"Person": 30, "Rater": 6, "Residual": 4})
    for source, value in expected.items():
        assert vc.loc[source, "RawVariance"] == pytest.approx(value, abs=1e-12)
        assert vc.loc[source, "Variance"] == pytest.approx(value, abs=1e-12)
    assert vc["df"].sum() == len(data)-1
    assert vc.SumSquares.sum() == pytest.approx(((data.Score-data.Score.mean())**2).sum())
    expected_g, expected_phi = (25.625/32, 25.625/41.25) if two_facets else (30/32, 30/35)
    assert out["observed_coefficients"]["G"] == pytest.approx(expected_g, abs=1e-12)
    assert out["observed_coefficients"]["Phi"] == pytest.approx(expected_phi, abs=1e-12)


@pytest.mark.parametrize("case,reason", [
    ("missing_cell", "missing_cells"), ("duplicate", "duplicate_cells"),
    ("single_level", "insufficient_levels"), ("nan", "invalid_scores"),
    ("infinite", "invalid_scores"), ("blank", "invalid_labels"),
    ("missing_label", "invalid_labels"), ("colliding_labels", "invalid_labels"),
    ("weights", "unsupported_weights"), ("third_facet", "unsupported_facets"),
    ("reserved_name", "invalid_columns"),
])
def test_unsupported_designs_never_fall_back(case, reason):
    data = ratings()
    if case == "missing_cell":
        data = data.iloc[:-1].copy()
    elif case == "duplicate":
        data = pd.concat([data, data.iloc[[0]]], ignore_index=True)
    elif case == "single_level":
        data = data[data.Rater == "1"].copy()
    elif case in ("nan", "infinite"):
        data.loc[0, "Score"] = np.nan if case == "nan" else np.inf
    elif case in ("blank", "missing_label", "colliding_labels"):
        data["Person"] = data["Person"].astype(object)
        data.loc[0, "Person"] = {"blank": " ", "missing_label": None,
                                 "colliding_labels": 1}[case]
    elif case == "weights":
        data["Weight"] = 1
        data.loc[0, "Weight"] = 2
    elif case == "third_facet":
        data["Occasion"] = "O1"
    elif case == "reserved_name":
        data = data.rename(columns={"Rater": "Residual"})
    original = data.copy(deep=True)
    fit = result(data)
    out = app.compute_generalizability_study(fit)
    assert not out["available"] and out["reason_code"] == reason
    assert out["d_study"].empty and out["variance_components"].empty
    assert out["design"]["n_input_rows"] == len(data)
    assert out["design"]["n_observations"] == 0
    assert not out["bootstrap_ci"]["inference_allowed"]
    assert app._crossed_anova_variance_components(
        data, object_facet="Person",
        random_facets=[c for c in fit["config"]["facet_names"] if c != "Person"],
    ).empty
    pd.testing.assert_frame_equal(data, original)


def test_missing_cells_detected_with_strings_and_unused_categorical_levels():
    data = ratings()
    data["Rater"] = pd.Categorical(data.Rater, categories=["-1", "1", "unused"])
    assert app.compute_generalizability_study(result(data))["available"]
    assert not app._is_balanced_one_obs_per_cell(data.iloc[:-1], ["Person", "Rater", "Item"])


@pytest.mark.parametrize("name", ["G", "Phi", "IsObservedDesign", "_section", "Rater:Item"])
def test_facet_names_cannot_overwrite_coefficients_or_interaction_labels(name):
    out = app.compute_generalizability_study(result(ratings().rename(columns={"Rater": name})))
    assert not out["available"] and out["reason_code"] == "invalid_columns"


@pytest.mark.parametrize("limit", [2, 10])
def test_observed_seven_raters_preserved_outside_planned_grid(limit):
    data = pd.DataFrame([
        {"Person": str(p), "Rater": str(r), "Item": str(i), "Score": 20+4*p+r+i+p*r/10}
        for p, r, i in itertools.product(range(3), range(7), range(2))])
    out = app.compute_generalizability_study(result(data), n_max_per_facet=limit)
    grid = out["d_study"]
    observed = grid[grid.IsObservedDesign]
    assert len(observed) == 1
    assert observed.iloc[0]["Rater"] == 7 and observed.iloc[0]["Item"] == 2
    assert observed.iloc[0]["G"] == out["observed_coefficients"]["G"]
    assert observed.iloc[0]["Phi"] == out["observed_coefficients"]["Phi"]
    assert (grid.loc[~grid.IsObservedDesign, ["Rater", "Item"]] <= limit).all().all()


@pytest.mark.parametrize("limit", [0, -1, True, 2.5])
def test_invalid_planning_count_is_explained(limit):
    out = app.compute_generalizability_study(result(ratings()), n_max_per_facet=limit)
    assert not out["available"] and out["reason_code"] == "invalid_grid"


def test_raw_negative_component_and_undefined_object_variance():
    data = ratings(False)
    data["Score"] -= 2*data.Rater.astype(int)
    out = app.compute_generalizability_study(result(data))
    assert out["available"]
    vc = out["variance_components"].set_index("Source")
    assert vc.loc["Rater", "RawVariance"] == -2
    assert vc.loc["Rater", "Variance"] == 0
    data["Score"] = 20+3*data.Rater.astype(int)
    out = app.compute_generalizability_study(result(data))
    assert not out["available"] and out["reason_code"] == "undefined_coefficients"
    assert not out["variance_components"].empty
    assert out["d_study"].empty


def test_source_and_settings_identity_without_mutation():
    data = ratings()
    first = app.compute_generalizability_study(result(data))
    same = app.compute_generalizability_study(result(data.copy()))
    assert first["design"] == same["design"]
    data.loc[0, "Score"] += .5
    changed = app.compute_generalizability_study(result(data))
    assert first["design"]["data_fingerprint"] != changed["design"]["data_fingerprint"]
    settings = app.compute_generalizability_study(result(data), n_max_per_facet=2)
    assert changed["design"]["analysis_id"] != settings["design"]["analysis_id"]
    assert changed["design"]["data_fingerprint"] == settings["design"]["data_fingerprint"]


def test_bootstrap_failures_keep_coefficient_specific_denominators(monkeypatch):
    real = app._crossed_anova_variance_components
    attempts = []
    def sometimes_fails(*args, **kwargs):
        attempts.append(1)
        if len(attempts) == 1:
            raise ValueError("controlled failure")
        if len(attempts) == 2:
            return pd.DataFrame()
        return real(*args, **kwargs)
    monkeypatch.setattr(app, "_crossed_anova_variance_components", sometimes_fails)
    monkeypatch.setattr(app, "_g_phi_from_variance_components", lambda *a, **k: {"G": .5, "Phi": np.nan})
    ci = app._bootstrap_g_phi_ci(
        ratings(), object_facet="Person", random_facets=["Rater", "Item"],
        score_col="Score", facet_observations={"Rater": 2, "Item": 2},
        n_bootstrap=4, extra_designs=[{"Rater": 3, "Item": 2}])
    assert ci["available"] and not ci["inference_allowed"]
    for entry in [ci]+ci["extra_design_replicates"]:
        assert entry["n_attempted"] == 4
        assert entry["n_success"] == 2 and entry["n_success_phi"] == 0
        assert entry["failure_counts"]["G"] == {"decomposition_exception:ValueError": 1, "invalid_components": 1}
        for coefficient, n_success in (("G", entry["n_success"]), ("Phi", entry["n_success_phi"])):
            assert n_success+sum(entry["failure_counts"][coefficient].values()) == 4


@pytest.mark.parametrize("valid", [True, False])
def test_csv_carries_scope_and_reason_without_unqualified_intervals(valid):
    data = ratings()
    if not valid:
        data = data.iloc[:-1]
    out = app.compute_generalizability_study(result(data))
    # A stale saved claim must not re-enable public intervals.
    out["bootstrap_ci"].update(available=True, inference_allowed=True, G_lower=.2, G_upper=.9)
    out["d_study"]["G_lower"] = .2
    out["d_study"]["G_upper"] = .9
    csv = app.to_csv_bytes(app.g_d_study_export_frame(out))
    frame = pd.read_csv(io.BytesIO(csv), keep_default_na=False)
    meta = frame[frame._section == "metadata"].set_index("Field").Value.to_dict()
    assert meta["calculation_version"] == "balanced_anova_v2"
    assert meta["interval_inference_allowed"] == "False"
    assert meta["interval_hold_reason"]
    assert meta["bootstrap_n_attempted"] == "0" and meta["bootstrap_n_requested"] == "0"
    assert meta["data_fingerprint"] and meta["n_input_rows"] == str(len(data))
    assert meta["weight_basis"] == "unit_observation_weights"
    assert "clipped" in meta["variance_policy"]
    assert "G_lower" not in frame and "G_upper" not in frame
    if valid:
        assert meta["decomposition"] == "full_3way_brennan_eq_3_18"
        observed = frame[frame._section == "observed_coefficients"]
        assert float(observed.iloc[0].G) == out["observed_coefficients"]["G"]
    else:
        assert meta["reason_code"] == "missing_cells"
        assert set(frame._section) == {"metadata"}


@pytest.mark.parametrize("lang", ["en", "ja"])
@pytest.mark.parametrize("valid", [True, False])
def test_panel_language_and_rerender_do_not_run_bootstrap(monkeypatch, lang, valid):
    def forbidden(*args, **kwargs):
        raise AssertionError("The public panel must not run bootstrap or refit MFRM")
    monkeypatch.setattr(app, "_bootstrap_g_phi_ci", forbidden)
    monkeypatch.setattr(app, "mfrm_estimate", forbidden)
    def view(language, valid_data):
        import itertools
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        st.session_state["lang"] = language
        data = pd.DataFrame([
            {"Person": str(p), "Rater": str(r), "Score": 20+4*p+2*r+p*r}
            for p, r in itertools.product((-1, 1), repeat=2)])
        if not valid_data:
            data = data.iloc[:-1]
        fit = {"prep": {"data": data}, "config": {"facet_names": ["Person", "Rater"]}}
        # Old widget state must not turn the removed bootstrap option back on.
        st.session_state[f"g_d_study_bootstrap_ci::{id(fit)}"] = True
        app.render_g_d_study_section(fit)
    at = AppTest.from_function(view, args=(lang, valid)).run(timeout=30)
    assert not at.exception
    assert not at.checkbox
    text = " ".join(c.value for c in at.caption)
    assert ("区間は保留" if lang == "ja" else "intervals are withheld") in text
    if valid:
        assert [m.value for m in at.metric] == ["0.938", "0.857"]
        assert next(d.value for d in at.dataframe if "IsObservedDesign" in d.value).IsObservedDesign.sum() == 1
    else:
        assert not at.metric
        assert ("組合せ" if lang == "ja" else "combinations") in at.info[0].value
    at.run(timeout=30)
    assert not at.exception
