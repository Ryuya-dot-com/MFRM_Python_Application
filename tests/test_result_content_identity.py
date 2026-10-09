"""Current scientific inputs and estimates must own their saved evidence."""
import ast
import builtins
from contextlib import nullcontext
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app
from tests.test_mml_prior_sensitivity import _small_mml_fit
from tests.test_app_joint_mml import arguments, sample


def test_complete_arrays_and_display_independence():
    left = np.full(1501, np.nan)
    right = left.copy()
    right[777] = 0.0
    assert str(left) == str(right)  # The original default=str defect.
    fingerprint = lambda x: app.config_fingerprint({"theta_spec": {"anchors": x}}, 64)
    assert fingerprint(left) != fingerprint(right)
    baseline = fingerprint(right)
    with np.printoptions(threshold=5, precision=1):
        assert fingerprint(right) == baseline
    assert fingerprint(right.copy()) == baseline
    assert fingerprint(right[:, None]) != baseline
    assert fingerprint(right.astype(np.float32)) != baseline
    matrix = np.arange(120).reshape(12, 10)
    assert fingerprint(matrix) == fingerprint(np.asfortranarray(matrix))


@pytest.mark.parametrize("left,right", [
    (pd.DataFrame({"x": [None]}, dtype=object), pd.DataFrame({"x": ["<NA>"]}, dtype=object)),
    (pd.DataFrame({"x": [1]}, dtype=object), pd.DataFrame({"x": ["1"]}, dtype=object)),
    (pd.DataFrame({1: [1]}), pd.DataFrame({"1": [1]})),
    (pd.DataFrame({"x": [1]}, index=[1]), pd.DataFrame({"x": [1]}, index=["1"])),
    (pd.DataFrame({"x": [1]}), pd.DataFrame({"x": [1]}).rename_axis("row")),
    (pd.DataFrame({"x": pd.Categorical(["a"], categories=["a", "b"])}),
     pd.DataFrame({"x": pd.Categorical(["a"], categories=["a", "c"])})),
    (pd.DataFrame({"x": pd.Categorical(["a"], ordered=False)}),
     pd.DataFrame({"x": pd.Categorical(["a"], ordered=True)})),
])
def test_table_content_types_missingness_and_category_definitions(left, right):
    assert app.dataframe_fingerprint(left, 64) != app.dataframe_fingerprint(right, 64)
    assert app.dataframe_fingerprint(left, 64) == app.dataframe_fingerprint(left.copy(), 64)


def test_container_types_order_missing_markers_and_unsupported_objects():
    fp = app.stable_json_fingerprint
    assert fp({"a": 1, "b": 2}) == fp({"b": 2, "a": 1})
    assert fp({1: "a"}) != fp({"1": "a"})
    assert fp([1, 2]) != fp((1, 2))
    assert fp({2, 1}) == fp({1, 2})
    markers = [None, pd.NA, pd.NaT, np.nan, "nan", "<NA>", ["float", "nan"]]
    assert len({fp(value) for value in markers}) == len(markers)
    with pytest.raises(TypeError, match="Unsupported fingerprint"):
        fp(object())
    with pytest.raises(TypeError, match="requires a DataFrame"):
        app.dataframe_fingerprint(object())
    app._self_test_reproducibility_fingerprints()


@pytest.fixture(scope="module")
def sensitivity_result():
    result = _small_mml_fit()
    result["config"].update(analysis_config_fingerprint="cached-config", run_fingerprint="cached-run",
                            input_data_fingerprint="recorded-upload")
    bundle = app.evaluate_mml_prior_sd_sensitivity(result, multipliers=(1.0, 1.25),
                                                 maxit=10, reltol=1e-3)
    assert bundle["available"]
    return result


def diagnostics():
    return {"obs": pd.DataFrame({"StdResidual": [0.0, 0.1]}),
            "reliability": pd.DataFrame(), "pca_enabled": False}


@pytest.mark.parametrize("mutation", ["anchor", "response", "point", "params", "reported_estimate", "summary"])
def test_stale_evidence_is_rejected_in_readiness_tables_and_json(sensitivity_result, mutation):
    before = sensitivity_result
    changed = deepcopy(before)
    if mutation == "anchor":
        changed["config"]["theta_spec"]["anchors"][2] = 0.125
    elif mutation == "response":
        changed["prep"]["data"].loc[0, "Score"] += 1
    elif mutation == "point":
        changed["opt"].x[0] += 0.01
    elif mutation == "params":
        changed["params"]["facets"]["Rater"][0] += 0.01
    elif mutation == "reported_estimate":
        changed["facets"]["person"].loc[0, "Estimate"] += 0.01
    else:
        changed["summary"].loc[0, "Converged"] = not bool(changed["summary"].loc[0, "Converged"])
    assert app.build_result_analysis_identity(before).analysis_id != app.build_result_analysis_identity(changed).analysis_id
    readiness = app.build_final_report_readiness(changed, diagnostics())
    row = readiness.set_index("Check").loc["MML fixed population prior SD"]
    assert row.StabilityState == "NOT_ASSESSED"
    assert "does not match" in row.Evidence
    frames = {}
    app._add_mml_prior_sensitivity_export_frames(frames, changed)
    assert frames == {}
    assets = app.build_evidence_contract_text_assets(changed, readiness)
    assert "mml_prior_sd_sensitivity_evidence_contract.json" not in assets
    previous_readiness = app.build_final_report_readiness(before, diagnostics())
    with pytest.raises(app._evidence.ContractValidationError, match="does not match"):
        app.build_evidence_contract_text_assets(changed, previous_readiness)


def test_current_evidence_survives_diagnostics_but_v1_binding_does_not(sensitivity_result, monkeypatch):
    result = deepcopy(sensitivity_result)
    identity = app.build_result_analysis_identity(result)
    app.mfrm_diagnostics(result, compute_pca=False)
    assert app.build_result_analysis_identity(result) == identity
    frames = {}
    app._add_mml_prior_sensitivity_export_frames(frames, result)
    assert "mml_prior_sd_sensitivity_report" in frames
    assets = app.build_evidence_contract_text_assets(result, app.build_final_report_readiness(result, diagnostics()))
    assert "mml_prior_sd_sensitivity_evidence_contract.json" in assets
    current_settings = app._result_identity_settings
    def old_settings(run):
        settings = current_settings(run)
        settings["schema_version"] = "mfrm_result_identity_settings_v1"
        settings.pop("normalized_response_fingerprint")
        settings.pop("source_fitted_state_fingerprint")
        return settings
    with monkeypatch.context() as patch:
        patch.setattr(app, "_result_identity_settings", old_settings)
        old = app.evaluate_mml_prior_sd_sensitivity(result, multipliers=(1.0, 1.25), maxit=10, reltol=1e-3)
        assert old["available"]
    # A structurally valid old contract still cannot authorize the current result.
    app._mml_prior_contract.validate_contract_bundle(old["evidence_contract"])
    with pytest.raises(app._evidence.ContractValidationError, match="does not match"):
        app._validate_mml_prior_contract_for_result(old["evidence_contract"], result)


def test_sensitivity_screen_never_displays_or_downloads_stale_cached_evidence(sensitivity_result, monkeypatch):
    result = deepcopy(sensitivity_result)
    bundle = result["evidence_bundles"]["mml_prior_sd"]
    class UI:
        def __init__(self):
            self.session_state = {}
            self.downloads = []
            self.messages = []
        def subheader(self, *a, **k): pass
        def caption(self, *a, **k): pass
        def info(self, message): self.messages.append(message)
        def success(self, message): self.messages.append(message)
        def warning(self, message): self.messages.append(message)
        def multiselect(self, *a, **k): return [1.0, 1.25]
        def number_input(self, *a, **k): return 10
        def selectbox(self, *a, **k): return 1e-3
        def columns(self, *a, **k): return [nullcontext(), nullcontext()]
        def button(self, *a, **k): return False
        def dataframe(self, *a, **k): pass
        def expander(self, *a, **k): return nullcontext()
        def download_button(self, *a, **k): self.downloads.append((a, k))
    ui = UI()
    monkeypatch.setattr(app, "st", ui)
    monkeypatch.setattr(app, "_render_compact_dataframe", lambda *a, **k: None)
    def key():
        digest = app.stable_json_fingerprint({
            "run": app.build_result_analysis_identity(result).analysis_id,
            "multipliers": [1.0, 1.25], "maxit": 10, "reltol": 1e-3,
            "base_sd": result["config"]["population_prior_sd"],
        })
        return "mml_prior_sd_sensitivity::" + digest
    ui.session_state[key()] = bundle
    app.show_mml_prior_sd_sensitivity_section(result)
    assert ui.downloads
    ui.downloads.clear()
    result["prep"]["data"].loc[0, "Score"] += 1
    app.show_mml_prior_sd_sensitivity_section(result)
    assert not ui.downloads
    ui.session_state[key()] = bundle  # Even forced reuse must pass the common guard.
    app.show_mml_prior_sd_sensitivity_section(result)
    assert not ui.downloads
    assert "failed validation" in ui.messages[-1]


def ui_metadata_update(result, initial_sd):
    """Execute the actual UI snapshot statement without starting a browser."""
    tree = ast.parse(Path(app.__file__).read_text())
    statements = [node for node in ast.walk(tree) if isinstance(node, ast.Expr)
                  and isinstance(node.value, ast.Call) and node.value.args
                  and isinstance(node.value.args[0], ast.Dict)
                  and "ui_mml_engine" in [key.value for key in node.value.args[0].keys
                                          if isinstance(key, ast.Constant)]]
    assert len(statements) == 1
    statement = statements[0]
    scope = dict(vars(builtins), **vars(app))
    for node in ast.walk(statement):
        if isinstance(node, ast.Name) and node.id not in scope:
            scope[node.id] = 0
    scope.update(result=result, est_method=result["config"]["method"], population_prior_sd=initial_sd,
                 selected_bias_pair=None, anchor_df=pd.DataFrame(), group_anchor_df=pd.DataFrame(),
                 person_data_df=pd.DataFrame(), quad_points=result["config"]["quad_points"])
    exec(compile(ast.Module(body=[statement], type_ignores=[]), "ui_snapshot", "exec"), scope)


@pytest.mark.parametrize("model,engine,free", [("RSM", "EM", True), ("RSM", "joint", True),
                                               ("PCM", "joint", True), ("GPCM", "joint", True),
                                               ("RSM", "direct", False)])
def test_ui_snapshot_preserves_fitted_scale_and_downstream_quadrature(model, engine, free):
    kwargs = arguments(sample(model), model, population_prior_sd=1.2)
    kwargs.update(mml_engine=engine, estimate_population_sd=free)
    result = app.mfrm_estimate(**kwargs)
    fitted = result["config"]["population_prior_sd"]
    if free:
        assert fitted == result["config"]["estimated_population_sd"]
        assert abs(fitted - 1.2) > 1e-5
    quad_before = app.make_mml_quadrature(result["config"], 31)
    scores_before = app.compute_scorefile(result)
    ui_metadata_update(result, 1.2)
    assert result["config"]["population_prior_sd"] == fitted
    assert result["config"]["ui_population_prior_sd"] == 1.2
    assert result["config"].get("population_prior_sd_input", 1.2) == 1.2
    np.testing.assert_array_equal(app.make_mml_quadrature(result["config"], 31)["nodes"], quad_before["nodes"])
    pd.testing.assert_frame_equal(app.compute_scorefile(result), scores_before)
    if engine == "joint":
        assert app.build_joint_mml_audit_frames(result)
