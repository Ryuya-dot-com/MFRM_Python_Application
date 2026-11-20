"""Installed public workflows; runnable independently of the application checkout."""
import importlib.util
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

import mfrmpy


def ratings():
    return pd.DataFrame([
        dict(Person=f"{p:03d}", Rater=f"R{r}", Task=f"T{t}", Score=(p + 2*r + t) % 3)
        for p in range(30) for r in range(3) for t in range(2)
    ])


def settings(model="RSM", method="JMLE"):
    return dict(person_col="Person", facet_cols=["Rater", "Task"], score_col="Score",
                rating_min=0, rating_max=2, model=model, method=method,
                step_facet="Task" if model == "PCM" else None,
                quad_points=61 if method == "MML" else 15, maxit=100)


def test_import_keeps_ui_and_plotting_optional(tmp_path):
    # Isolated Python cannot accidentally import modules from the checkout.
    script = """
import importlib.abc, sys
class NoUI(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'streamlit', 'streamlit_app', 'tkinter', 'matplotlib', 'plotly'}:
            raise AssertionError('Core imported optional UI: ' + fullname)
sys.meta_path.insert(0, NoUI())
import mfrmpy
data = mfrmpy.load_example()
chk = mfrmpy.check_design(data, person='Person', facets=['Rater', 'Task', 'Criterion'],
                          score='Score', rating_min=data['Score'].min(), rating_max=data['Score'].max())
assert chk['ok']
analysis = mfrmpy.analyze(data, person='Person', facets=['Rater', 'Task', 'Criterion'],
                         score='Score', rating_min=data['Score'].min(), rating_max=data['Score'].max())
assert mfrmpy.tables(analysis)
assert len(mfrmpy.tables(analysis, 'measures')) == 1
assert not mfrmpy.summary(analysis).empty
"""
    subprocess.run([sys.executable, "-I", "-c", script], cwd=tmp_path, check=True)


def test_example_is_bundled_and_not_shared_mutable_state():
    data = mfrmpy.load_example()
    assert set(data) == {"Person", "Rater", "Task", "Criterion", "Score"}
    assert not data.empty and data.Score.between(0, 5).all()
    data.iloc[0, data.columns.get_loc("Score")] = -100
    assert mfrmpy.load_example().Score.between(0, 5).all()


@pytest.fixture(scope="module", params=[("RSM", "JMLE"), ("PCM", "JMLE"),
                                       ("RSM", "MML"), ("PCM", "MML")])
def fitted(request):
    model, method = request.param
    return mfrmpy.analyze(ratings(), **settings(model, method),
        bias_requests=[dict(facet_a="Rater", facet_b="Task", omit_extreme=False)],
        interval_requests=[dict(facet="Rater", difference=["R0", "R2"], level=.9)]
            if method == "MML" else None)


def forbid_refit(monkeypatch):
    from mfrm_app import engine, bias_analysis, generalizability, structural_intervals, joint_interaction
    def forbidden(*args, **kwargs):
        pytest.fail("Reading retained public results must not refit")
    for module, names in (
        (engine, ["mfrm_estimate", "mfrm_diagnostics"]),
        (bias_analysis, ["estimate_bias_interaction"]),
        (generalizability, ["compute_generalizability_from_data"]),
        (structural_intervals, ["compute_structural_intervals"]),
        (joint_interaction, ["fit_joint_interaction"]),
    ):
        for name in names:
            monkeypatch.setattr(module, name, forbidden)


def test_model_routes_tables_and_zero_refit_roundtrip(fitted, tmp_path, monkeypatch):
    before = mfrmpy.tables(fitted)
    forbid_refit(monkeypatch)
    path = tmp_path / "保存.mfrm"
    mfrmpy.save_project(path, fitted)
    restored = mfrmpy.load_project(path)
    assert restored["analysis_id"] == fitted["analysis_id"]
    assert restored["content_fingerprint"] == fitted["content_fingerprint"]
    pd.testing.assert_frame_equal(restored["input_data"], ratings())
    assert restored["input_data"].Person.iloc[0] == "000"
    after = mfrmpy.tables(restored)
    assert len(before) == len(after)
    for (title, frame, note), (other_title, other, other_note) in zip(before, after):
        assert (title, note) == (other_title, other_note)
        pd.testing.assert_frame_equal(frame, other)
        if title.startswith("Bias:"):
            assert frame.Use.isin(["screening_only", "withheld"]).all()
            assert frame.SupportsFormalInference.eq(False).all()
            assert {"ObsScore", "ExpScore", "ObsExpAvg", "BiasSize", "SE", "Z", "LR_ChiSq", "LR_Prob", "StatusReason"}.issubset(frame.columns)
        if title.startswith("Specified intervals:"):
            assert frame.Use.eq("Approximate").all()
            assert frame.SupportsFormalInference.eq(False).all()
            assert "Reason" in frame.columns and frame.Reason.eq("").all()
            assert pd.to_numeric(frame.SE).gt(0).all()
            assert "finite-sample" in note
    assert {"Step Thresholds", "Separation and Reliability by Facet"}.issubset({title for title, _, _ in after})


def test_joint_interaction_and_its_qualified_table(tmp_path, monkeypatch):
    analysis = mfrmpy.analyze(ratings(), **settings(method="MML"),
        joint_settings=dict(facet_a="Rater", facet_b="Task"))
    assert analysis["joint_analysis"]["status"] == "computed"
    forbid_refit(monkeypatch)
    path = tmp_path / "joint.mfrm"
    mfrmpy.save_project(path, analysis)
    entries = mfrmpy.tables(mfrmpy.load_project(path))
    title, frame, note = next(entry for entry in entries if entry[0].startswith("Joint interaction:"))
    assert not frame.empty and note
    assert "CIEligible" in frame or "Reason" in frame


@pytest.mark.parametrize("two_facets", [False, True])
def test_gd_planning_and_unavailable_design(two_facets, tmp_path, monkeypatch):
    data = pd.DataFrame([
        dict(Person=f"{p:03d}", Rater=f"R{r}", Task=f"T{t}",
             Score=.4*p + .3*r + .2*t + .1*((p+r+t) % 3))
        for p in range(8) for r in range(3) for t in range(2)
    ])
    facets = ["Rater", "Task"] if two_facets else ["Rater"]
    if not two_facets:
        data = data[data.Task.eq("T0")].drop(columns="Task")
    gd_settings = dict(object_facet="Person", random_facets=facets,
                       selected_plan={f: 4 for f in facets}, target_population="Synthetic population")
    analysis = mfrmpy.analyze(data, workflow="gd", gd_settings=gd_settings)
    missing = mfrmpy.analyze(data.iloc[1:], workflow="gd", gd_settings=gd_settings)
    assert analysis["gd_analysis"]["result"]["available"]
    assert missing["gd_analysis"]["result"]["reason_code"] == "missing_cells"
    forbid_refit(monkeypatch)
    for index, value in enumerate((analysis, missing)):
        path = tmp_path / f"gd-{index}.mfrm"
        mfrmpy.save_project(path, value)
        restored = mfrmpy.load_project(path)
        assert restored["content_fingerprint"] == value["content_fingerprint"]
        assert mfrmpy.tables(restored)


@pytest.mark.parametrize("changed", [dict(model="GPCM"), dict(estimate_population_sd=True),
                                     dict(rating_min=0.5), dict(facet_cols=["Unknown"])])
def test_unsupported_or_invalid_requests_are_rejected_before_estimation(changed, monkeypatch):
    from mfrm_app import engine
    def forbidden(*args, **kwargs):
        pytest.fail("Invalid request reached estimation")
    monkeypatch.setattr(engine, "mfrm_estimate", forbidden)
    with pytest.raises(ValueError):
        mfrmpy.analyze(ratings(), **(settings() | changed))


def test_corrupt_project_is_not_accepted(tmp_path):
    path = tmp_path / "bad.mfrm"
    path.write_text("not a project", encoding="utf-8")
    with pytest.raises(mfrmpy.ProjectError):
        mfrmpy.load_project(path)


@pytest.mark.parametrize("storage", ["python", "pyarrow"])
@pytest.mark.parametrize("missing", [pd.NA, np.nan], ids=["NA", "NaN"])
def test_project_retains_string_storage_missing_values_and_fingerprints(storage, missing, tmp_path, monkeypatch):
    dtype = pd.StringDtype(storage=storage, na_value=missing)
    data = ratings().astype({name: dtype for name in ("Person", "Rater", "Task")})
    data["Label"] = pd.Series(["001", None, "nan"] * (len(data) // 3), dtype=dtype)
    data.columns = data.columns.astype(dtype)
    analysis = mfrmpy.analyze(data, **settings())
    forbid_refit(monkeypatch)
    path = tmp_path / "typed-strings.mfrm"
    mfrmpy.save_project(path, analysis)
    restored = mfrmpy.load_project(path)
    assert restored["content_fingerprint"] == analysis["content_fingerprint"]
    assert restored["analysis_id"] == analysis["analysis_id"]
    pd.testing.assert_frame_equal(restored["input_data"], data)
    assert restored["input_data"].Label.dtype.storage == storage
    assert restored["input_data"].Label.iloc[0] == "001"
    assert restored["input_data"].Label.iloc[2] == "nan"
    assert mfrmpy.tables(restored)


@pytest.mark.skipif(importlib.util.find_spec("matplotlib") is None, reason="Optional figures extra")
@pytest.mark.figures
def test_optional_figures_and_report_reuse_retained_results(fitted, monkeypatch):
    forbid_refit(monkeypatch)
    for kind in ("wright", "pathway", "category"):
        level = "T0" if fitted["fit_settings"]["model"] == "PCM" else None
        figure = mfrmpy.plot(fitted, kind, level=level)
        assert figure.axes
        import io
        buffer = io.BytesIO()
        figure.savefig(buffer, format="png", dpi=72)
        assert buffer.getvalue().startswith(b"\x89PNG")
    report = mfrmpy.report_html(fitted)
    assert "<html" in report.lower() and "<svg" in report.lower()


def varied_ratings():
    rng = np.random.default_rng(42)
    rows = [
        dict(Person=f"P{p:02d}", Rater=f"R{r}", Task=f"T{t}",
             Score=int(np.clip(np.round(1.0 + 0.08 * (p - 7.5) - 0.3 * (r - 1) - 0.2 * (t - 0.5) + rng.normal(0, 0.7)), 0, 2)))
        for p in range(16) for r in range(3) for t in range(2)
    ]
    df = pd.DataFrame(rows)
    for _, grp in df.groupby("Person"):
        if grp["Score"].nunique() == 1:
            df.loc[grp.index[0], "Score"] = 1
    return df


def test_screening_only_se_intervals_and_four_series_reliability(fitted, tmp_path, monkeypatch):
    entries = {title: (frame, note) for title, frame, note in mfrmpy.tables(fitted)}
    measures, m_note = entries["Measures and Fit"]
    assert {"ModelSE", "RealSE", "Lower", "Upper", "Basis", "SEUse", "CIUse", "SupportsFormalInference", "FitStatus", "FitReason"}.issubset(measures.columns)
    assert measures["SEUse"].eq("screening_only").all()
    assert measures["CIUse"].eq("screening_only").all()
    assert measures["SupportsFormalInference"].eq(False).all()
    model_se = pd.to_numeric(measures["ModelSE"], errors="coerce")
    real_se = pd.to_numeric(measures["RealSE"], errors="coerce")
    est = pd.to_numeric(measures["Estimate"], errors="coerce")
    lo = pd.to_numeric(measures["Lower"], errors="coerce")
    hi = pd.to_numeric(measures["Upper"], errors="coerce")
    assert model_se.gt(0).all() and (real_se + 1e-12 >= model_se).all()
    assert (lo < est).all() and (est < hi).all()
    assert "screening_only" in m_note and "SupportsFormalInference = FALSE" in m_note

    # On varied data, all 4 series (Model/Real × Sample/Population) are finite and ordered.
    model_name = fitted["fit_settings"]["model"]
    method_name = fitted["fit_settings"]["method"]
    varied_fit = mfrmpy.analyze(varied_ratings(), **settings(model_name, method_name))
    rel, r_note = dict((t, (f, n)) for t, f, n in mfrmpy.tables(varied_fit))["Separation and Reliability by Facet"]
    assert set(zip(rel["SEMode"], rel["DistributionBasis"])) == {
        ("model", "sample"), ("model", "population"),
        ("fit_adjusted", "sample"), ("fit_adjusted", "population"),
    }
    assert len(rel) == 4 * 3  # Person, Rater, Task × 4 series
    assert rel["Use"].eq("screening_only").all()
    assert rel["SupportsFormalInference"].eq(False).all()
    assert "SupportsFormalInference = FALSE" in r_note

    for facet in ("Person", "Rater", "Task"):
        sub = rel.loc[rel["Facet"].eq(facet)].set_index(["SEMode", "DistributionBasis"])
        assert (sub["EstimateAvailable"] == sub["Levels"]).all()
        assert (sub["ExcludedLevels"] == 0).all()
        assert float(sub.loc[("model", "sample"), "SD"]) > float(sub.loc[("model", "population"), "SD"])
        assert float(sub.loc[("fit_adjusted", "sample"), "RMSE"]) >= float(sub.loc[("model", "sample"), "RMSE"])
        assert float(sub.loc[("fit_adjusted", "sample"), "Reliability"]) <= float(sub.loc[("model", "sample"), "Reliability"]) + 1e-6

    forbid_refit(monkeypatch)
    written = mfrmpy.export_csv(varied_fit, tmp_path / "csv_out")
    assert written[0].name == "00_manifest.csv"
    manifest = pd.read_csv(written[0])
    assert manifest["SupportsFormalInference"].eq(False).all()
    assert len(written) == 1 + len(mfrmpy.tables(varied_fit))


def test_audit_size_limit_rank_deficiency_and_non_convergence_are_separated(tmp_path, monkeypatch):
    from copy import deepcopy
    from mfrm_app import basic_analysis, engine, output_qualification, project

    # 1) Connected JMLE where dense SVD hit audit_size_limit: unverified dense rank,
    #    yet screening_only conditional SEs remain visible with SupportsFormalInference=False.
    base = mfrmpy.analyze(ratings(), **settings("RSM", "JMLE"))
    size_limited = deepcopy(base)
    size_limited["result"]["identifiability"]["status"] = "audit_size_limit"
    size_limited["qualification_at_save"] = basic_analysis.qualification_snapshot(size_limited["result"])
    size_limited["content_fingerprint"] = project.content_fingerprint(size_limited)

    # 2) Actually unidentifiable (disconnected subsets without anchors): SEs/intervals withheld.
    disc_data = pd.DataFrame([
        dict(Person=f"P{p}", Rater=f"R{r}", Task="T0" if p < 4 else "T1", Score=(p + r) % 3)
        for p in range(8)
        for r in ((0, 1) if p < 4 else (2, 3))
    ])
    rank_def = mfrmpy.analyze(disc_data, **settings("RSM", "JMLE"))

    # 3) Non-converged (maxit=1): SEs/intervals withheld and reason is non-convergence,
    #    even when audit_size_limit or fixed-SD MML is also active.
    non_conv_jmle = mfrmpy.analyze(ratings(), **(settings("RSM", "JMLE") | {"maxit": 1}))
    non_conv_mml = mfrmpy.analyze(ratings(), **(settings("RSM", "MML") | {"maxit": 1}))

    st_limit = output_qualification.run_qualification_state(size_limited["result"])
    st_rank = output_qualification.run_qualification_state(rank_def["result"])
    st_nc_j = output_qualification.run_qualification_state(non_conv_jmle["result"])
    st_nc_m = output_qualification.run_qualification_state(non_conv_mml["result"])

    assert st_limit["state"] == "audit_size_limit_unverified"
    assert st_limit["formal_hold_reason"] == "stat.jmle.identifiability_audit_size_limit"
    assert st_limit["screening_block_reason"] == ""

    assert st_rank["state"] == "unidentifiable"
    assert st_rank["formal_hold_reason"] == "stat.jmle.eta_rank_deficient"
    assert st_rank["screening_block_reason"] == "stat.jmle.eta_rank_deficient"

    assert st_nc_j["state"] == st_nc_m["state"] == "non_converged"
    assert st_nc_j["formal_hold_reason"] == st_nc_m["formal_hold_reason"] == "stat.estimation.optimizer_termination_unverified"
    assert st_nc_j["screening_block_reason"] == st_nc_m["screening_block_reason"] == "stat.estimation.optimizer_termination_unverified"

    forbid_refit(monkeypatch)
    m_limit = dict((t, f) for t, f, _ in mfrmpy.tables(size_limited))["Measures and Fit"]
    assert m_limit["Basis"].eq("Conditional (audit size limit)").all()
    assert m_limit["SEUse"].eq("screening_only").all()
    assert m_limit["SupportsFormalInference"].eq(False).all()
    assert m_limit["StatusReason"].eq("stat.jmle.identifiability_audit_size_limit").all()
    assert pd.to_numeric(m_limit["SE / SD"], errors="coerce").gt(0).all()

    for label, blocked, expected_reason in (
        ("rank_def", rank_def, "stat.jmle.eta_rank_deficient"),
        ("non_conv_jmle", non_conv_jmle, "stat.estimation.optimizer_termination_unverified"),
        ("non_conv_mml", non_conv_mml, "stat.estimation.optimizer_termination_unverified"),
    ):
        tbls = dict((t, f) for t, f, _ in mfrmpy.tables(blocked))
        mf = tbls["Measures and Fit"]
        rel = tbls["Separation and Reliability by Facet"]
        assert mf["SEUse"].eq("withheld").all(), label
        assert mf["SupportsFormalInference"].eq(False).all(), label
        assert mf[["SE / SD", "Lower", "Upper"]].eq("").all().all(), label
        assert mf["StatusReason"].eq(expected_reason).all(), label
        assert rel["Use"].eq("withheld").all(), label
        assert rel["StatusReason"].eq(expected_reason).all(), label


@pytest.mark.figures
def test_boundary_exclusion_counts_low_df_fit_guards_and_legacy_reload(tmp_path, monkeypatch):
    from copy import deepcopy
    from mfrm_app import project

    data = varied_ratings()
    data.loc[data["Person"].eq("P00"), "Score"] = 0
    data.loc[data["Person"].eq("P15"), "Score"] = 2
    analysis = mfrmpy.analyze(data, **settings("RSM", "JMLE"))

    tbls = {t: f for t, f, _ in mfrmpy.tables(analysis)}
    mf = tbls["Measures and Fit"].set_index(["Facet", "Level"])
    for extreme_p in ("P00", "P15"):
        assert mf.loc[("Person", extreme_p), "Basis"] == "No finite JMLE"
        assert mf.loc[("Person", extreme_p), "FitStatus"] == "unavailable_extreme_score"
        assert mf.loc[("Person", extreme_p), ["Estimate", "SE / SD", "Infit", "InfitZSTD"]].eq("").all()
    # Rater rows include observations on extreme persons: MnSq stays descriptive while ZSTD is protected.
    assert mf.loc[("Rater", "R0"), "Infit"] != ""
    assert mf.loc[("Rater", "R0"), "InfitZSTD"] == ""
    assert mf.loc[("Rater", "R0"), "FitStatus"] == "zstd_unavailable_boundary_proxy"
    assert "JMLE optimizer terminal proxies" in mf.loc[("Rater", "R0"), "FitReason"]

    rel = tbls["Separation and Reliability by Facet"]
    person_rel = rel.loc[rel["Facet"].eq("Person")]
    assert len(person_rel) == 4
    assert (person_rel["Levels"] == 16).all()
    assert (person_rel["EstimateAvailable"] == 14).all()
    assert (person_rel["ExcludedLevels"] == 2).all()
    assert person_rel["AggregationTarget"].str.contains("2 extreme boundary excluded").all()
    assert person_rel["Use"].eq("screening_only").all()
    assert pd.to_numeric(person_rel["Reliability"], errors="coerce").between(0.0, 1.0).all()

    # Simulate an older 0.1.0 saved MML archive with withheld status and a pre-guard df < 1 ZSTD = 9.0.
    mml_legacy = deepcopy(mfrmpy.analyze(ratings(), **settings("RSM", "MML")))
    mml_legacy["qualification_at_save"] = {
        "schema": "mfrm_output_qualification_v1",
        "optimizer_terminated": True,
        "calibration_hold_reason": "stat.mml_fixed_sd.numerical_qualification_pending",
        "calibration_hold_message": "Fixed-population-SD MML numerical accuracy and sampling uncertainty remain unqualified.",
    }
    mml_legacy["diagnostics"]["measures"].loc[0, "DF_Infit"] = 0.25
    mml_legacy["diagnostics"]["measures"].loc[0, "DF_Infit_FACETS"] = 0.25
    mml_legacy["diagnostics"]["measures"].loc[0, "InfitZSTD"] = 9.0
    mml_legacy["diagnostics"]["measures"].loc[0, "InfitZSTD_FACETS"] = 9.0
    mml_legacy["diagnostics"]["measures"] = mml_legacy["diagnostics"]["measures"].drop(
        columns=[c for c in ("ModelSE", "RealSE", "FitStatus", "FitReason") if c in mml_legacy["diagnostics"]["measures"].columns]
    )
    mml_legacy["content_fingerprint"] = project.content_fingerprint(mml_legacy)

    archive_path = tmp_path / "legacy_010_mml.mfrm"
    project.save_project(archive_path, mml_legacy, compress=False)
    assert not archive_path.read_bytes().startswith(project.GZIP_MAGIC)

    gz_path = tmp_path / "compressed_020_mml.mfrm"
    mfrmpy.save_project(gz_path, mml_legacy)
    assert gz_path.read_bytes().startswith(project.GZIP_MAGIC)
    assert len(gz_path.read_bytes()) * 3 < len(archive_path.read_bytes())

    forbid_refit(monkeypatch)
    restored = mfrmpy.load_project(archive_path)
    restored_gz = mfrmpy.load_project(gz_path)
    assert restored_gz["content_fingerprint"] == restored["content_fingerprint"]
    assert restored["content_fingerprint"] == mml_legacy["content_fingerprint"]
    assert float(restored["diagnostics"]["measures"].loc[0, "InfitZSTD"]) == 9.0
    assert restored["qualification_at_save"]["calibration_hold_reason"] == "stat.mml_fixed_sd.numerical_qualification_pending"

    re_tables = {t: f for t, f, _ in mfrmpy.tables(restored)}
    re_mf = re_tables["Measures and Fit"]
    assert re_mf.loc[0, "InfitZSTD"] == ""
    assert float(re_mf.loc[0, "SavedInfitZSTD"]) == 9.0
    assert re_mf.loc[0, "FitStatus"] == "zstd_unavailable_low_df"
    assert "degrees of freedom are below 1" in re_mf.loc[0, "FitReason"]
    assert re_mf["SEUse"].eq("screening_only").all()
    assert re_mf["SupportsFormalInference"].eq(False).all()
    assert "SavedReasonCode" in re_mf.columns

    csv_files = mfrmpy.export_csv(restored, tmp_path / "legacy_csv")
    manifest = pd.read_csv(csv_files[0])
    assert manifest["SavedHoldReason"].eq("stat.mml_fixed_sd.numerical_qualification_pending").all()
    measures_csv = pd.read_csv(next(p for p in csv_files if "Measures_and_Fit" in p.name))
    assert pd.isna(measures_csv.loc[0, "InfitZSTD"])
    assert float(measures_csv.loc[0, "SavedInfitZSTD"]) == 9.0

    if importlib.util.find_spec("matplotlib") is not None:
        html_doc = mfrmpy.report_html(restored)
        assert "Saved qualification at archive creation (retained)" in html_doc
        assert "SupportsFormalInference=FALSE" in html_doc
        fig = mfrmpy.plot(restored, "pathway", statistic="InfitZSTD")
        assert fig.axes


@pytest.mark.figures
def test_diagnostic_plots_coordinates_exclusions_and_qualification(monkeypatch):
    pytest.importorskip("matplotlib")
    from copy import deepcopy

    data = varied_ratings()
    an = mfrmpy.analyze(
        data,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="RSM",
        method="JMLE",
        compute_pca=True,
        bias_requests=[{"facet_a": "Rater", "facet_b": "Task", "omit_extreme": False}],
    )
    pcm_an = mfrmpy.analyze(
        data,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="PCM",
        step_facet="Task",
        method="MML",
        quad_points=15,
        population_prior_sd=1.0,
    )
    forbid_refit(monkeypatch)

    # 1. Expected score & Information curves match independent closed-form Rasch math
    steps = np.asarray(an["result"]["params"]["steps"], dtype=float)
    fig_exp = mfrmpy.plot(an, "expected")
    eta = fig_exp.axes[0].lines[-1].get_xdata()
    weights = np.exp(np.outer(eta, np.arange(len(steps) + 1)) - np.r_[0.0, np.cumsum(steps)])
    probs = weights / weights.sum(axis=1, keepdims=True)
    ks = np.arange(len(steps) + 1, dtype=float)
    expected_score = probs @ ks
    np.testing.assert_allclose(fig_exp.axes[0].lines[-1].get_ydata(), expected_score, rtol=1e-12)
    assert mfrmpy.plot(an, "tcc").axes

    fig_info = mfrmpy.plot(an, "information", presentation="manuscript")
    info_expected = np.maximum((probs * (ks - expected_score[:, None]) ** 2).sum(axis=1), 1e-12)
    sem_expected = 1.0 / np.sqrt(info_expected)
    np.testing.assert_allclose(fig_info.axes[0].lines[0].get_ydata(), info_expected, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(fig_info.axes[1].lines[0].get_ydata(), sem_expected, rtol=1e-11, atol=1e-12)
    assert fig_info.axes[0].get_title() == "Category information"
    assert "SupportsFormalInference=FALSE" in " ".join(fig_info._mfrm_caption["notes"])

    # PCM multi-level overlay when level=None vs single-level when level="T0"
    assert len(mfrmpy.plot(pcm_an, "step").axes[0].lines) >= 3  # zero ref + 2 Task ladders
    assert len(mfrmpy.plot(pcm_an, "thresholds", level="T0").axes[0].lines) == 2

    # Disordered steps highlight crossing without mutating archive
    disordered = deepcopy(an)
    disordered["result"]["params"]["steps"] = np.array([0.6, -0.6])
    fig_dis = mfrmpy.plot(disordered, "step")
    assert "disordered" in fig_dis._mfrm_caption["title"].lower()

    # 2. Facet / rater_severity: screening intervals on converged, point-only on rank-deficient, withheld on non-converged
    fig_rater = mfrmpy.plot(an, "rater_severity")
    assert any(type(c).__name__ == "ErrorbarContainer" for c in fig_rater.axes[0].containers)
    assert "SupportsFormalInference=FALSE" in " ".join(fig_rater._mfrm_caption["notes"])

    rank_def = deepcopy(an)
    rank_def["result"]["identifiability"] = {"status": "eta_rank_deficient"}
    fig_rd = mfrmpy.plot(rank_def, "facet")
    assert not any(type(c).__name__ == "ErrorbarContainer" for c in fig_rd.axes[0].containers)

    non_conv = deepcopy(an)
    non_conv["result"]["opt"].success = False
    fig_nc = mfrmpy.plot(non_conv, "facet")
    assert any("Withheld:" in t.get_text() for t in fig_nc.axes[0].texts)

    # 3. Pathway & Bubble support Outfit / OutfitZSTD and enforce df < 1 ZSTD guard
    for stat in ("Infit", "Outfit", "InfitZSTD", "OutfitZSTD"):
        assert mfrmpy.plot(an, "pathway", statistic=stat).axes
        assert mfrmpy.plot(an, "bubble", statistic=stat).axes
    assert mfrmpy.plot(an, "fit_pathway", statistic="Outfit").axes

    low_df = deepcopy(pcm_an)
    low_df["diagnostics"]["measures"].loc[
        low_df["diagnostics"]["measures"]["Facet"].eq("Rater"), ["DF_Outfit", "DF_Outfit_FACETS"]
    ] = 0.4
    fig_bub_low = mfrmpy.plot(low_df, "bubble", facet="Rater", statistic="OutfitZSTD")
    assert any("degrees of freedom < 1" in " ".join(t.get_text().split()) for t in fig_bub_low.axes[0].texts)

    # 4. Fair average coordinates match report_tables
    fig_fair = mfrmpy.plot(an, "fair_average", facet="Rater")
    r_tbl = an["report_tables"]["Rater"]
    pts = fig_fair._mfrm_hover_points[0][0].get_offsets()
    np.testing.assert_allclose(pts[:, 0], r_tbl["ObservedAverage"].to_numpy(dtype=float))
    np.testing.assert_allclose(pts[:, 1], r_tbl["FairM"].to_numpy(dtype=float))

    # 5. Unexpected, Residual Q-Q (excluding extreme persons), PCA loadings, Bias ranked, and HTML report
    assert mfrmpy.plot(an, "unexpected").axes
    fig_qq = mfrmpy.plot(an, "residual_qq")
    qq_df = fig_qq._mfrm_hover_points[0][1]
    ext_persons = set(
        an["result"]["facets"]["person"]
        .loc[an["result"]["facets"]["person"]["ExtremeScorePattern"].fillna(False), "Person"]
        .astype(str)
    )
    assert not (set(qq_df["Person"].astype(str)) & ext_persons)

    fig_load = mfrmpy.plot(an, "pca_loadings")
    assert len(fig_load.axes[0].get_yticklabels()) == len(an["diagnostics"]["pca"]["loadings"])
    assert mfrmpy.plot(an, "residual_pca").axes

    fig_br = mfrmpy.plot(an, "bias_ranked", index=0)
    assert fig_br.axes
    report = mfrmpy.report_html(an)
    assert "S.E." not in report and "LR Prob." not in report


@pytest.mark.figures
def test_boundary_scenarios_all_ones_category_gaps_disconnected_mml_and_zero_df(tmp_path, monkeypatch):
    import warnings

    def table_map(analysis):
        return {title: (frame, note) for title, frame, note in mfrmpy.tables(analysis)}

    # 1) Endpoint constant rater (all-1 on 1..4): JMLE & opposite-endpoint MML withhold unbounded rater without RuntimeWarning.
    ep_rows = [
        dict(
            Person=f"P{p:02d}",
            Rater=r,
            Task=t,
            Score=(1 if r == "R_all1" else 1 + ((p + (0 if r == "R1" else 2) + (0 if t == "T1" else 1)) % 4)),
        )
        for p in range(12) for r in ("R1", "R2", "R_all1") for t in ("T1", "T2")
    ]
    ep_df = pd.DataFrame(ep_rows)
    ep_jmle = mfrmpy.analyze(
        ep_df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=1,
        rating_max=4,
        model="RSM",
        method="JMLE",
        bias_requests=[{"facet_a": "Rater", "facet_b": "Task"}],
    )
    ep_tbls = table_map(ep_jmle)
    ep_mf = ep_tbls["Measures and Fit"][0].set_index(["Facet", "Level"])
    assert ep_mf.loc[("Rater", "R_all1"), "Estimate"] == ""
    assert ep_mf.loc[("Rater", "R_all1"), "Infit"] == ""
    assert ep_mf.loc[("Rater", "R_all1"), "FitStatus"] == "unavailable_boundary"
    assert int(ep_mf.loc[("Rater", "R_all1"), "ObservedCategories"]) == 1
    assert pd.notna(pd.to_numeric(ep_mf.loc[("Rater", "R_all1"), "SavedEstimate"], errors="coerce"))
    ep_rel = ep_tbls["Separation and Reliability by Facet"][0]
    ep_r_rel = ep_rel.loc[ep_rel["Facet"].eq("Rater")].iloc[0]
    assert int(ep_r_rel["EstimateAvailable"]) == 2 and int(ep_r_rel["ExcludedLevels"]) == 1

    # Case B: Fixed-SD MML with noncenter_facet="Person" (centered Rater & Task) has a genuinely finite MLE for single endpoint rater.
    ep_mml_case_b = mfrmpy.analyze(
        ep_df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=1,
        rating_max=4,
        model="RSM",
        method="MML",
        quad_points=15,
    )
    mb_mf = table_map(ep_mml_case_b)["Measures and Fit"][0].set_index(["Facet", "Level"])
    assert mb_mf.loc[("Rater", "R_all1"), "Estimate"] != ""
    assert mb_mf.loc[("Rater", "R_all1"), "SEUse"] == "screening_only"

    # Opposite-endpoint free raters in centered MML (S8) diverge along a sum-preserving recession ray; must hold with zero RuntimeWarning.
    opp_rows = [
        dict(
            Person=f"P{p:02d}",
            Rater=r,
            Task=t,
            Score=(1 if r == "R_min" else 4 if r == "R_max" else 1 + ((p + (0 if t == "T1" else 1)) % 4)),
        )
        for p in range(12) for r in ("R1", "R_min", "R_max") for t in ("T1", "T2")
    ]
    with warnings.catch_warnings(record=True) as caught_warns:
        warnings.simplefilter("always")
        opp_mml = mfrmpy.analyze(
            pd.DataFrame(opp_rows),
            person_col="Person",
            facet_cols=["Rater", "Task"],
            score_col="Score",
            rating_min=1,
            rating_max=4,
            model="RSM",
            method="MML",
            quad_points=15,
            bias_requests=[{"facet_a": "Rater", "facet_b": "Task"}],
            interval_requests=[{"facet": "Rater"}],
            joint_settings={"facet_a": "Rater", "facet_b": "Task"},
        )
    runtime_divs = [w for w in caught_warns if issubclass(w.category, RuntimeWarning) and "divide" in str(w.message).lower()]
    assert not runtime_divs
    opp_tbls = table_map(opp_mml)
    opp_mf = opp_tbls["Measures and Fit"][0].set_index(["Facet", "Level"])
    for r_ext in ("R_min", "R_max"):
        assert opp_mf.loc[("Rater", r_ext), "Estimate"] == ""
        assert opp_mf.loc[("Rater", r_ext), "FitStatus"] == "unavailable_boundary"
    assert opp_mf["SEUse"].ne("screening_only").all()
    assert opp_tbls["Specified intervals: Rater"][0]["Use"].eq("Withheld").all()
    assert opp_tbls["Joint interaction: Rater × Task"][0]["Use"].eq("Withheld").all()
    assert opp_tbls["Bias: Rater × Task"][0]["Use"].eq("Unavailable").all()

    # 2) Interior constant rater (all-1 on 0..3): RSM agreement/G-D notes + PCM step_facet="Rater" hold & partner fit guard.
    int_rows = [
        dict(
            Person=f"P{p:02d}",
            Rater=r,
            Task=t,
            Score=(
                1
                if r == "R_all1"
                else (p // 4) + (1 if (p % 2 == 1 and t == "T2") else 0)
                if r == "R1"
                else min(3, (p // 4) + (1 if t == "T1" else 0) + (1 if p >= 8 else 0))
            ),
        )
        for p in range(12) for r in ("R1", "R2", "R_all1") for t in ("T1", "T2")
    ]
    int_df = pd.DataFrame(int_rows)
    int_rsm = mfrmpy.analyze(
        int_df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=0,
        rating_max=3,
        model="RSM",
        method="JMLE",
        gd_settings={
            "object_facet": "Person",
            "random_facets": ["Rater", "Task"],
            "score_col": "Score",
            "target_population": "Synthetic cohort",
        },
    )
    int_rsm_tbls = table_map(int_rsm)
    rsm_mf, rsm_mf_note = int_rsm_tbls["Measures and Fit"]
    assert int(rsm_mf.set_index(["Facet", "Level"]).loc[("Rater", "R_all1"), "ObservedCategories"]) == 1
    assert "ObservedCategories = 1" in rsm_mf_note
    _, agr_note = int_rsm_tbls["Inter-rater Agreement and Krippendorff's Alpha"]
    assert "zero score variance" in agr_note and "excluded from MeanCorr" in agr_note
    gd_meta = int_rsm_tbls["G/D: metadata"][0].set_index("Field")["Value"].to_dict()
    assert "R_all1" in gd_meta["zero_variance_levels"]
    assert "Zero-variance level(s) detected" in gd_meta["scope"]

    int_pcm = mfrmpy.analyze(
        int_df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=0,
        rating_max=3,
        model="PCM",
        step_facet="Rater",
        method="JMLE",
        bias_requests=[{"facet_a": "Rater", "facet_b": "Task"}],
    )
    pcm_tbls = table_map(int_pcm)
    pcm_mf = pcm_tbls["Measures and Fit"][0].set_index(["Facet", "Level"])
    assert pcm_mf.loc[("Rater", "R_all1"), "Estimate"] == ""
    assert pcm_mf.loc[("Rater", "R_all1"), "FitStatus"] == "unavailable_boundary"
    assert pcm_mf.loc[("Task", "T1"), "Outfit"] == ""  # Collapsed boundary-proxy rows do not distort Task Outfit
    assert (pcm_mf["StatusReason"] == "stat.category.unsupported_internal_step_contrast").any()
    assert pcm_tbls["Bias: Rater × Task"][0]["Use"].eq("Unavailable").all()

    # 3) Unobserved internal category (score 1 absent on 0..3) & unobserved boundary category (score 3 absent on 0..3).
    gap_df = int_df.copy()
    gap_df["Score"] = gap_df["Score"].replace({1: 2})
    gap_an = mfrmpy.analyze(
        gap_df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=0,
        rating_max=3,
        model="RSM",
        method="MML",
        quad_points=15,
        interval_requests=[{"facet": "Rater"}],
    )
    gap_tbls = table_map(gap_an)
    gap_mf = gap_tbls["Measures and Fit"][0]
    assert gap_mf["SEUse"].ne("screening_only").all()
    assert (gap_mf["StatusReason"] == "stat.category.unsupported_internal_step_contrast").any()
    assert gap_tbls["Specified intervals: Rater"][0]["Use"].eq("Withheld").all()

    b_gap_df = int_df.copy()
    b_gap_df["Score"] = b_gap_df["Score"].replace({3: 2})
    b_gap_an = mfrmpy.analyze(
        b_gap_df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=0,
        rating_max=3,
        model="RSM",
        method="JMLE",
    )
    b_gap_path = tmp_path / "boundary_gap.mfrm"
    mfrmpy.save_project(b_gap_path, b_gap_an)
    forbid_refit(monkeypatch)
    b_gap_re = mfrmpy.load_project(b_gap_path)
    monkeypatch.undo()
    b_re_mf = table_map(b_gap_re)["Measures and Fit"][0]
    assert b_re_mf["SEUse"].ne("screening_only").all()
    assert (b_re_mf["StatusReason"] == "stat.category.unobserved_boundary_category").any()
    if importlib.util.find_spec("matplotlib") is not None:
        fig_step_held = mfrmpy.plot(b_gap_re, "step")
        assert any("Withheld:" in t.get_text() for t in fig_step_held.axes[0].texts)

    # 4) Disconnected MML cohorts (Cohort A vs Cohort B sharing no persons, raters, or tasks) + JMLE SavedReliability hygiene.
    disc_rows = []
    for cohort, p_off, raters, tasks in ((0, 0, ("R1", "R2"), ("T1", "T2")), (1, 6, ("R3", "R4"), ("T3", "T4"))):
        for p in range(6):
            for r_i, r in enumerate(raters):
                for t_i, t in enumerate(tasks):
                    disc_rows.append(
                        dict(Person=f"P{p_off + p:02d}", Rater=r, Task=t, Score=(p + r_i + t_i + cohort) % 3)
                    )
    disc_df = pd.DataFrame(disc_rows)
    disc_mml = mfrmpy.analyze(
        disc_df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="RSM",
        method="MML",
        quad_points=15,
        bias_requests=[{"facet_a": "Rater", "facet_b": "Task"}],
        interval_requests=[{"facet": "Rater", "difference": ["R1", "R3"]}],
        joint_settings={"facet_a": "Rater", "facet_b": "Task"},
    )
    disc_tbls = table_map(disc_mml)
    disc_mf = disc_tbls["Measures and Fit"][0]
    assert disc_mf["SEUse"].ne("screening_only").all()
    assert (disc_mf["StatusReason"] == "stat.mml.population_assumption_linked").all()
    disc_iv, disc_iv_note = disc_tbls["Specified intervals: Rater"]
    assert disc_iv["Use"].eq("Withheld").all()
    assert "person_facet_connectivity" in disc_iv_note
    assert disc_tbls["Joint interaction: Rater × Task"][0]["Use"].eq("Withheld").all()
    assert disc_tbls["Bias: Rater × Task"][0]["Use"].eq("Unavailable").all()

    disc_jmle = mfrmpy.analyze(
        disc_df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="RSM",
        method="JMLE",
    )
    disc_j_rel = table_map(disc_jmle)["Separation and Reliability by Facet"][0]
    assert disc_j_rel["Use"].eq("withheld").all()
    assert "SavedReliability" not in disc_j_rel.columns

    # 5) Single-observation free JMLE parameter (N = 1) withholds 0-df fit mean squares instead of reporting 0.000.
    one_obs_df = pd.concat(
        [
            varied_ratings(),
            pd.DataFrame([{"Person": "P00", "Rater": "R_single", "Task": "T0", "Score": 1}]),
        ],
        ignore_index=True,
    )
    one_obs_an = mfrmpy.analyze(
        one_obs_df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="RSM",
        method="JMLE",
    )
    one_mf = table_map(one_obs_an)["Measures and Fit"][0].set_index(["Facet", "Level"])
    assert one_mf.loc[("Rater", "R_single"), "Infit"] == ""
    assert one_mf.loc[("Rater", "R_single"), "Outfit"] == ""
    assert one_mf.loc[("Rater", "R_single"), "FitStatus"] == "unavailable_zero_df"


@pytest.mark.figures
def test_api_consolidation_check_design_summary_table_kinds_and_aliases(monkeypatch):
    df = ratings()

    # 1) Optimization-free check_design for MFRM and G/D
    forbid_refit(monkeypatch)
    ok_chk = mfrmpy.check_design(
        df,
        person="Person",
        facets=["Rater", "Task"],
        score="Score",
        rating_min=df["Score"].min(),
        rating_max=df["Score"].max(),
    )
    assert ok_chk["ok"] is True
    assert ok_chk["connectivity_ok"] is True
    assert ok_chk["category_support_ok"] is True
    assert ok_chk["nullity"] == 0

    gap_chk = mfrmpy.check_design(
        df,
        person="Person",
        facets=["Rater", "Task"],
        score="Score",
        rating_min=0,
        rating_max=3,
    )
    assert gap_chk["ok"] is False
    assert gap_chk["connectivity_ok"] is True
    assert gap_chk["category_support_ok"] is False
    assert gap_chk["unobserved_boundary_categories"] == [3]

    gd_ok = mfrmpy.check_design(df, workflow="gd", person="Person", facets=["Rater", "Task"], score="Score")
    assert gd_ok["ok"] is True
    gd_bad = mfrmpy.check_design(df.iloc[1:], workflow="gd", person="Person", facets=["Rater", "Task"], score="Score")
    assert gd_bad["ok"] is False and gd_bad["reason_code"] == "missing_cells"
    monkeypatch.undo()

    # 2) Column role aliases + numpy.integer / numpy.bool_ scalars match canonical calls identically
    canonical = mfrmpy.analyze(
        df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        method="MML",
        quad_points=15,
        maxit=80,
    )
    aliased = mfrmpy.analyze(
        df,
        person="Person",
        facets=["Rater", "Task"],
        score="Score",
        rating_min=df["Score"].min(),
        rating_max=df["Score"].max(),
        method="MML",
        quad_points=np.int64(15),
        maxit=np.int64(80),
        compute_pca=np.bool_(False),
    )
    assert aliased["analysis_id"] == canonical["analysis_id"]
    assert aliased["fit_settings"] == canonical["fit_settings"]

    with pytest.raises(ValueError, match="person or person_col"):
        mfrmpy.analyze(df, person="Person", person_col="Person", facets=["Rater"], score="Score", rating_min=0, rating_max=2)
    with pytest.raises(ValueError):
        mfrmpy.analyze(df, person="Person", facets=["Rater", "Task"], score="Score", rating_min=True, rating_max=2)

    # 3) tables(analysis, kind=...) filtering + summary(analysis) without refitting
    gd_an = mfrmpy.analyze(
        df,
        workflow="gd",
        gd_settings={
            "object_facet": "Person",
            "random_facets": ["Rater", "Task"],
            "n_max_per_facet": np.int64(6),
            "selected_plan": {"Rater": np.int64(4), "Task": np.int64(2)},
        },
    )
    forbid_refit(monkeypatch)
    assert [t for t, _, _ in mfrmpy.tables(aliased, "measures")] == ["Measures and Fit"]
    assert [t for t, _, _ in mfrmpy.tables(aliased, "estimates")] == ["Measures and Fit"]
    assert [t for t, _, _ in mfrmpy.tables(aliased, "reliability")] == ["Separation and Reliability by Facet"]
    assert len(mfrmpy.tables(aliased, "agreement")) == 2
    assert len(mfrmpy.tables(aliased, "agreement_summary")) == 1
    assert len(mfrmpy.tables(aliased, "pairwise_agreement")) == 1
    assert [t for t, _, _ in mfrmpy.tables(aliased, ("design", "categories"))] == [
        "Analysis Design and Retained Ratings",
        "Observed Category Use",
    ]
    assert len(mfrmpy.tables(gd_an, "gd")) == 5
    assert [t for t, _, _ in mfrmpy.tables(gd_an, "gd_variance")] == ["G/D: variance_components"]
    with pytest.raises(ValueError, match="Unknown table kind"):
        mfrmpy.tables(aliased, "nonexistent_kind")

    s_mfrm = mfrmpy.summary(aliased).set_index("Field")["Value"].to_dict()
    assert s_mfrm["Workflow"] == "mfrm"
    assert s_mfrm["Model"] == "RSM" and s_mfrm["Method"] == "MML"
    assert "measures" in s_mfrm["AvailableTables"]
    assert "rater_severity" in s_mfrm["AvailablePlots"]

    s_gd = mfrmpy.summary(gd_an).set_index("Field")["Value"].to_dict()
    assert s_gd["Workflow"] == "gd" and s_gd["Available"] == "True"
    assert s_gd["AvailableTables"] == "gd" and s_gd["AvailablePlots"] == "gd"

    # 4) report_html accepts facet/statistic aliases matching plot()
    if importlib.util.find_spec("matplotlib") is not None:
        h1 = mfrmpy.report_html(aliased, pathway_facet="Rater", pathway_statistic="Outfit")
        h2 = mfrmpy.report_html(aliased, facet="Rater", statistic="Outfit")
        assert h1 == h2
        with pytest.raises(ValueError, match="facet or pathway_facet"):
            mfrmpy.report_html(aliased, facet="Rater", pathway_facet="Task")


def test_cross_engine_parity_with_mfrmr_030_and_tam_4325():
    """Pin numerical agreement with mfrmr 0.3.0 and TAM 4.3.25 across JMLE, MML, intervals and G/D."""
    from mfrm_app import engine

    rng = np.random.default_rng(20261118)
    persons = [f"P{i:02d}" for i in range(1, 31)]
    raters = ["R1", "R2", "R3", "R4"]
    tasks = ["T1", "T2", "T3"]
    theta = np.linspace(-1.45, 1.45, len(persons))
    theta -= theta.mean()
    b_rater = {"R1": -0.65, "R2": -0.15, "R3": 0.25, "R4": 0.55}
    b_task = {"T1": -0.45, "T2": 0.05, "T3": 0.40}
    cum_steps = np.r_[0.0, np.cumsum([-1.15, 0.05, 1.10])]
    p_map = dict(zip(persons, theta))
    rows = []
    for p in persons:
        for r in raters:
            for t in tasks:
                inter = 0.55 if (r == "R1" and t == "T1") else (-0.25 if (r == "R4" and t == "T3") else 0.0)
                probs = engine.category_prob_rsm(np.array([p_map[p] - b_rater[r] - b_task[t] + inter]), cum_steps)
                cdf = np.cumsum(probs, axis=1)
                score = int((rng.random(1)[:, None] >= cdf[:, :-1]).sum(axis=1)[0])
                rows.append({"Person": p, "Rater": r, "Task": t, "Score": score})
    df = pd.DataFrame(rows)
    for p, m in df.groupby("Person")["Score"].mean().items():
        if m <= 0.0:
            df.loc[df.index[df["Person"] == p][0], "Score"] = 1
        elif m >= 3.0:
            df.loc[df.index[df["Person"] == p][0], "Score"] = 2

    jmle = mfrmpy.analyze(
        df, person="Person", facets=["Rater", "Task"], score="Score",
        rating_min=0, rating_max=3, model="RSM", method="JMLE", maxit=800, reltol=1e-11,
        gd_settings={"object_facet": "Person", "random_facets": ["Rater", "Task"], "selected_plan": {"Rater": 4, "Task": 3}},
    )
    assert jmle["result"]["summary"]["LogLik"].iloc[0] == pytest.approx(-379.1866896972, abs=1e-6)
    assert jmle["result"]["params"]["steps"][0] == pytest.approx(-1.30619422, abs=2e-5)
    qj = engine.qualify_measures_for_reporting(jmle["diagnostics"]["measures"], jmle["result"], obs_df=jmle["diagnostics"]["obs"])
    r1_j = qj.loc[(qj["Facet"] == "Rater") & (qj["Level"] == "R1")].iloc[0]
    assert float(r1_j["Estimate"]) == pytest.approx(-0.62556917, abs=2e-5)
    assert float(r1_j["ModelSE"]) == pytest.approx(0.14364675, abs=1e-5)
    assert float(r1_j["RealSE"]) == pytest.approx(0.14534520, abs=1e-5)
    assert float(r1_j["Infit"]) == pytest.approx(1.02378737, abs=2e-5)
    assert float(r1_j["InfitZSTD"]) == pytest.approx(0.21292925, abs=5e-5)

    gd_sel = jmle["gd_analysis"]["result"]["selected_plan"]
    assert float(gd_sel["G"]) == pytest.approx(0.855446765862, abs=1e-11)
    assert float(gd_sel["Phi"]) == pytest.approx(0.784668693869, abs=1e-11)
    assert float(gd_sel["RelativeSEM"]) == pytest.approx(0.244675091464, abs=1e-11)
    assert float(gd_sel["AbsoluteSEM"]) == pytest.approx(0.311804782231, abs=1e-11)

    mml = mfrmpy.analyze(
        df, person="Person", facets=["Rater", "Task"], score="Score",
        rating_min=0, rating_max=3, model="RSM", method="MML", quad_points=81, maxit=1000, reltol=1e-12,
        interval_requests=[{"facet": "Rater"}],
    )
    # Matches both mfrmr 0.3.0 (-424.1575682751) and TAM 4.3.25 (-424.1575684631)
    assert mml["result"]["summary"]["LogLik"].iloc[0] == pytest.approx(-424.15756846, abs=1e-6)
    assert mml["result"]["params"]["steps"][0] == pytest.approx(-1.17827148, abs=2e-5)
    qm = engine.qualify_measures_for_reporting(mml["diagnostics"]["measures"], mml["result"], obs_df=mml["diagnostics"]["obs"])
    r1_m = qm.loc[(qm["Facet"] == "Rater") & (qm["Level"] == "R1")].iloc[0]
    assert float(r1_m["Estimate"]) == pytest.approx(-0.57767929, abs=2e-5)
    assert float(r1_m["ModelSE"]) == pytest.approx(0.12387120, abs=1e-5)
    assert float(r1_m["Infit"]) == pytest.approx(0.92304306, abs=2e-5)
    assert float(r1_m["InfitZSTD"]) == pytest.approx(-0.52573641, abs=5e-5)

    rel_tbl = mfrmpy.tables(mml, "reliability")[0][1]
    r_rel = rel_tbl.loc[(rel_tbl["Facet"] == "Rater") & (rel_tbl["SEMode"] == "model") & (rel_tbl["DistributionBasis"] == "sample")].iloc[0]
    p_rel = rel_tbl.loc[(rel_tbl["Facet"] == "Person") & (rel_tbl["SEMode"] == "model") & (rel_tbl["DistributionBasis"] == "population")].iloc[0]
    assert float(r_rel["Reliability"]) == pytest.approx(0.94417929, abs=5e-4)  # public table rounded to 3 dp: 0.944
    assert float(p_rel["Reliability"]) == pytest.approx(0.86098196, abs=5e-4)  # TAM 4.3.25 EAP.rel = 0.86098196

    int_rec = pd.DataFrame(mml["interval_analyses"][0]["record"]["table"]).set_index("Target").loc["R1"]
    assert float(int_rec["Lower"]) == pytest.approx(-0.82046239, abs=2e-5)
    assert float(int_rec["Upper"]) == pytest.approx(-0.33489620, abs=2e-5)


def test_plan_a_quad31_default_step_thresholds_scale_orientation_and_bias_screening(tmp_path, monkeypatch):
    """Verify Plan A refinements: MML default quad_points=31, Step Thresholds, noncenter/positive_facets, and Bias screening."""
    from mfrm_app import engine

    rng = np.random.default_rng(20261118)
    persons = [f"P{i:02d}" for i in range(1, 31)]
    raters = ["R1", "R2", "R3", "R4"]
    tasks = ["T1", "T2", "T3"]
    theta = np.linspace(-1.45, 1.45, len(persons))
    theta -= theta.mean()
    b_rater = {"R1": -0.65, "R2": -0.15, "R3": 0.25, "R4": 0.55}
    b_task = {"T1": -0.45, "T2": 0.05, "T3": 0.40}
    cum_steps = np.r_[0.0, np.cumsum([-1.15, 0.05, 1.10])]
    p_map = dict(zip(persons, theta))
    rows = []
    for p in persons:
        for r in raters:
            for t in tasks:
                inter = 0.55 if (r == "R1" and t == "T1") else (-0.25 if (r == "R4" and t == "T3") else 0.0)
                probs = engine.category_prob_rsm(np.array([p_map[p] - b_rater[r] - b_task[t] + inter]), cum_steps)
                cdf = np.cumsum(probs, axis=1)
                score = int((rng.random(1)[:, None] >= cdf[:, :-1]).sum(axis=1)[0])
                rows.append({"Person": p, "Rater": r, "Task": t, "Score": score})
    df = pd.DataFrame(rows)
    for p, m in df.groupby("Person")["Score"].mean().items():
        if m <= 0.0:
            df.loc[df.index[df["Person"] == p][0], "Score"] = 1
        elif m >= 3.0:
            df.loc[df.index[df["Person"] == p][0], "Score"] = 2

    # A1: Default MML quad_points=31 passes Specified intervals & Joint interaction without explicit quad_points;
    #     explicit quad_points=15 surfaces Reason in the public Specified intervals DataFrame.
    mml_default = mfrmpy.analyze(
        df, person="Person", facets=["Rater", "Task"], score="Score",
        rating_min=0, rating_max=3, model="RSM", method="MML",
        interval_requests=[{"facet": "Rater", "difference": ["R1", "R4"]}],
        joint_settings={"facet_a": "Rater", "facet_b": "Task"},
        bias_requests=[{"facet_a": "Rater", "facet_b": "Task"}],
    )
    assert mml_default["fit_settings"]["quad_points"] == 31
    assert mml_default["result"]["config"]["quad_points"] == 31
    iv_def = mfrmpy.tables(mml_default, "intervals")[0][1]
    assert iv_def["Use"].eq("Approximate").all()
    assert iv_def["SupportsFormalInference"].eq(False).all()
    assert iv_def["Reason"].eq("").all()

    mml_q15 = mfrmpy.analyze(
        df, person="Person", facets=["Rater", "Task"], score="Score",
        rating_min=0, rating_max=3, model="RSM", method="MML", quad_points=15,
        interval_requests=[{"facet": "Rater"}],
    )
    iv_q15 = mfrmpy.tables(mml_q15, "intervals")[0][1]
    assert iv_q15["Use"].eq("Withheld").all()
    assert iv_q15["SupportsFormalInference"].eq(False).all()
    assert iv_q15["Reason"].str.contains("quadrature_score_shift|full_covariance_change").all()

    with pytest.raises(ValueError, match="only to MML"):
        mfrmpy.analyze(df, person="Person", facets=["Rater", "Task"], score="Score", rating_min=0, rating_max=3, method="JMLE", quad_points=31)

    # A2: Step Thresholds table (RSM & PCM × JMLE & MML + disordered steps + category hold)
    st_mml = mfrmpy.tables(mml_default, "steps")[0][1]
    assert mfrmpy.tables(mml_default, "thresholds")[0][0] == "Step Thresholds"
    assert list(st_mml["Step"]) == ["Step_1", "Step_2", "Step_3"]
    assert list(st_mml["Transition"]) == ["0->1", "1->2", "2->3"]
    assert st_mml["StepFacet"].eq("Common").all()
    assert st_mml["SEUse"].eq("screening_only").all()
    assert st_mml["SupportsFormalInference"].eq(False).all()
    assert pd.to_numeric(st_mml["SE"]).gt(0).all()  # includes final sum-to-zero step SE!
    assert st_mml["Disordered"].eq(False).all()

    pcm_jmle = mfrmpy.analyze(
        df, person="Person", facets=["Rater", "Task"], score="Score",
        rating_min=0, rating_max=3, model="PCM", step_facet="Task", method="JMLE",
    )
    st_pcm = mfrmpy.tables(pcm_jmle, "steps")[0][1]
    assert set(st_pcm["StepFacet"]) == {"T1", "T2", "T3"}
    assert len(st_pcm) == 9
    assert st_pcm["SEUse"].eq("screening_only").all()
    assert pd.to_numeric(st_pcm["SE"]).gt(0).all()

    # Disordered intermediate category (rare score=1 produces tau_2 < tau_1 => Disordered=True)
    dis_df = df.copy()
    dis_df.loc[dis_df["Score"] == 1, "Score"] = 2
    dis_df.loc[dis_df.index[:3], "Score"] = 1
    dis_an = mfrmpy.analyze(dis_df, person="Person", facets=["Rater", "Task"], score="Score", rating_min=0, rating_max=3)
    st_dis = mfrmpy.tables(dis_an, "steps")[0][1].set_index("Step")
    assert bool(st_dis.loc["Step_2", "Disordered"]) is True
    assert float(st_dis.loc["Step_2", "Spacing"]) < 0.0

    # Unobserved internal category withholds step estimates and SEs
    gap_df = df.copy()
    gap_df.loc[gap_df["Score"] == 1, "Score"] = 2
    gap_an = mfrmpy.analyze(gap_df, person="Person", facets=["Rater", "Task"], score="Score", rating_min=0, rating_max=3)
    st_gap = mfrmpy.tables(gap_an, "steps")[0][1]
    assert st_gap["Estimate"].eq("").all()
    assert st_gap["SE"].eq("").all()
    assert st_gap["SEUse"].eq("withheld").all()
    assert st_gap["StatusReason"].eq("stat.category.unsupported_internal_step_contrast").all()

    # A3: noncenter_facet="Task" matches TAM 4.3.25 & ConQuest 5.47.5 cases-constraint logLik (-423.8544085);
    #     positive_facets=["Rater"] flips Rater orientation to +1 leniency without changing logLik or other parameters;
    #     and .mfrm round-trips with zero refit.
    chk_nc = mfrmpy.check_design(
        df, person="Person", facets=["Rater", "Task"], score="Score",
        rating_min=0, rating_max=3, noncenter_facet="Task", positive_facets=["Rater"],
    )
    assert chk_nc["ok"] is True
    assert chk_nc["noncenter_facet"] == "Task"
    assert chk_nc["positive_facets"] == ["Rater"]
    with pytest.raises(ValueError, match="noncenter_facet"):
        mfrmpy.check_design(df, person="Person", facets=["Rater", "Task"], score="Score", noncenter_facet="Unknown")

    mml_cases = mfrmpy.analyze(
        df, person="Person", facets=["Rater", "Task"], score="Score",
        rating_min=0, rating_max=3, model="RSM", method="MML",
        noncenter_facet="Task", positive_facets=["Rater"],
        quad_points=81, maxit=1000, reltol=1e-12,
        interval_requests=[{"facet": "Rater", "difference": ["R1", "R4"]}],
    )
    mml_cases_neg = mfrmpy.analyze(
        df, person="Person", facets=["Rater", "Task"], score="Score",
        rating_min=0, rating_max=3, model="RSM", method="MML",
        noncenter_facet="Task",
        quad_points=81, maxit=1000, reltol=1e-12,
    )
    assert mml_cases["result"]["summary"]["LogLik"].iloc[0] == pytest.approx(-423.8544085, abs=5e-5)
    assert mml_cases["result"]["summary"]["LogLik"].iloc[0] == pytest.approx(
        mml_cases_neg["result"]["summary"]["LogLik"].iloc[0], abs=1e-10
    )
    # Under positive_facets=["Rater"], lenient R1 has exact opposite sign (+0.5792) to subtractive severity (-0.5792)
    q_cases = engine.qualify_measures_for_reporting(mml_cases["diagnostics"]["measures"], mml_cases["result"], obs_df=mml_cases["diagnostics"]["obs"])
    q_neg = engine.qualify_measures_for_reporting(mml_cases_neg["diagnostics"]["measures"], mml_cases_neg["result"], obs_df=mml_cases_neg["diagnostics"]["obs"])
    r1_pos = float(q_cases.loc[(q_cases["Facet"] == "Rater") & (q_cases["Level"] == "R1"), "Estimate"].iloc[0])
    r1_neg = float(q_neg.loc[(q_neg["Facet"] == "Rater") & (q_neg["Level"] == "R1"), "Estimate"].iloc[0])
    assert r1_pos > 0.57 and r1_pos == pytest.approx(-r1_neg, abs=1e-8)
    # Rater remains sum-to-zero centered while non-centered Task absorbs the grand-mean offset (-0.1527 logits)
    r_mean = float(q_cases.loc[q_cases["Facet"] == "Rater", "Estimate"].mean())
    t_mean = float(q_cases.loc[q_cases["Facet"] == "Task", "Estimate"].mean())
    assert abs(r_mean) < 1e-10 and t_mean == pytest.approx(-0.152704, abs=1e-3)
    sum_cases = mfrmpy.summary(mml_cases).set_index("Field")["Value"].to_dict()
    assert sum_cases["NoncenterFacet"] == "Task"
    assert sum_cases["PositiveFacets"] == "Rater"
    assert "steps" in sum_cases["AvailableTables"]

    # A4: Conditional Bias screening table + withholding under calibration block + .mfrm zero-refit reload
    bias_tbl = mfrmpy.tables(mml_default, "bias")[0][1]
    assert bias_tbl["Use"].eq("screening_only").all()
    assert bias_tbl["SupportsFormalInference"].eq(False).all()
    raw_b = mml_default["bias_analyses"][0]["result"]["table"]
    np.testing.assert_allclose(pd.to_numeric(bias_tbl["BiasSize"]), raw_b["Bias Size"].round(3))
    np.testing.assert_allclose(pd.to_numeric(bias_tbl["SE"]), raw_b["S.E."].round(3))
    np.testing.assert_allclose(pd.to_numeric(bias_tbl["Z"]), raw_b["Standardized Bias"].round(3))
    np.testing.assert_allclose(pd.to_numeric(bias_tbl["LR_ChiSq"]), raw_b["LR ChiSq"].round(3))
    np.testing.assert_allclose(pd.to_numeric(bias_tbl["LR_Prob"]), raw_b["LR Prob."].round(4))

    nonconv_bias = mfrmpy.analyze(
        df, person="Person", facets=["Rater", "Task"], score="Score",
        rating_min=0, rating_max=3, maxit=1, bias_requests=[{"facet_a": "Rater", "facet_b": "Task"}],
    )
    b_held = mfrmpy.tables(nonconv_bias, "bias")[0][1]
    assert b_held["Use"].eq("withheld").all()
    assert b_held[["SE", "Z", "LR_ChiSq", "LR_Prob"]].eq("").all().all()
    assert b_held["StatusReason"].eq("stat.estimation.optimizer_termination_unverified").all()

    nc_path = tmp_path / "cases_pos.mfrm"
    mfrmpy.save_project(nc_path, mml_cases)
    forbid_refit(monkeypatch)
    nc_re = mfrmpy.load_project(nc_path)
    assert nc_re["fit_settings"]["noncenter_facet"] == "Task"
    assert nc_re["fit_settings"]["positive_facets"] == ["Rater"]
    pd.testing.assert_frame_equal(mfrmpy.tables(nc_re, "steps")[0][1], mfrmpy.tables(mml_cases, "steps")[0][1])
    pd.testing.assert_frame_equal(mfrmpy.tables(nc_re, "intervals")[0][1], mfrmpy.tables(mml_cases, "intervals")[0][1])


def test_public_score_preparation_and_preflight():
    data = ratings().astype({"Score": float})
    original = data.Score.copy()
    data.loc[data.Score.eq(1), "Score"] -= 1e-10
    data = data.rename(columns={"Rater": "score_k"})
    cols = dict(person="Person", facets=["score_k", "Task"], score="Score", rating_min=0, rating_max=2)
    fit = mfrmpy.analyze(data, **cols)
    prepared = fit["result"]["prep"]["data"]
    assert prepared.score_k.astype(str).tolist() == data.score_k.tolist()
    np.testing.assert_array_equal(prepared.Score, original)
    assert mfrmpy.check_design(pd.concat([data, data]), **cols)["density_pct"] == 100
    for invalid in [1.5, np.inf]:
        bad = data.copy()
        bad.iloc[0, bad.columns.get_loc("Score")] = invalid
        assert not mfrmpy.check_design(bad, **cols)["ok"]
        with pytest.raises(ValueError, match="row"):
            mfrmpy.analyze(bad, **cols)


def test_anchor_labels_and_unknown_rank():
    data = ratings()
    data.loc[data.Rater.eq("R1"), "Rater"] = " R1 "
    cols = dict(person="Person", facets=["Rater", "Task"], score="Score", rating_min=0, rating_max=2)
    anchor = [dict(Facet="Rater", Level="R1", Anchor=-0.6)]
    assert mfrmpy.check_design(data, **cols, anchors=anchor)["ok"]
    fit = mfrmpy.analyze(data, **cols, anchors=anchor)
    assert "R1" in fit["result"]["prep"]["levels"]["Rater"]
    data.loc[2, "Rater"] = "R1"
    assert not mfrmpy.check_design(data, **cols, anchors=anchor)["ok"]
    with pytest.raises(ValueError, match="R1"):
        mfrmpy.analyze(data, **cols, anchors=anchor)
    large = pd.DataFrame([(p,r,r,(p+r)%3) for p in range(1600) for r in range(2)], columns=["Person","Rater","Task","Score"])
    pre = mfrmpy.check_design(large, **cols)
    assert not pre["ok"] and pre["nullity"] is None


@pytest.mark.parametrize("scores,expected", [([0,2,2,2], (0,0)), ([0,2,0,2], (np.nan,0))])
def test_zero_object_variance_by_coefficient(scores, expected):
    data = pd.DataFrame(dict(Person=["A","A","B","B"], Rater=["X","Y","X","Y"], Score=scores))
    cfg = dict(object_facet="Person", random_facets=["Rater"])
    fit = mfrmpy.analyze(data, workflow="gd", gd_settings=cfg)
    result = fit["gd_analysis"]["result"]
    assert result["available"]
    for _, frame, _ in mfrmpy.tables(fit, ("gd_observed", "gd_selected")):
        np.testing.assert_allclose(frame[["G", "Phi"]].iloc[0], expected, equal_nan=True)
    for coeff in (result["observed_coefficients"], result["selected_plan"]):
        np.testing.assert_allclose([coeff["G"],coeff["Phi"]], expected, equal_nan=True)


def test_partial_anchor_links_only_estimable_levels(tmp_path):
    cols = dict(person="Person", facets=["Rater", "Task"], score="Score", rating_min=0, rating_max=2)
    data = ratings()
    data = data.loc[(data.Person.astype(int) % 3 == 0) == data.Rater.eq("R2")]
    anchor = dict(Facet="Rater", Level="R0", Anchor=-0.6)
    fit = mfrmpy.analyze(data, **cols, method="MML", anchors=[anchor])
    path = tmp_path / "partial.mfrm"
    mfrmpy.save_project(path, fit)
    for record in [fit, mfrmpy.load_project(path)]:
        table = mfrmpy.tables(record, "measures")[0][1].set_index(["Facet", "Level"])
        assert table.loc[("Rater", "R2"), "SEUse"] == "withheld"
        assert table.loc[("Rater", "R1"), "SEUse"] == "screening_only"
    # A mean across both cohorts supplies a single equation for two independent shifts.
    group = [dict(Facet="Rater", Level=level, Group="reference", GroupValue=0)
             for level in ("R0", "R2")]
    record = mfrmpy.analyze(data, **cols, method="MML", group_anchors=group)
    table = mfrmpy.tables(record, "measures")[0][1]
    assert table.loc[table.Facet.eq("Rater"), "SEUse"].eq("withheld").all()


def test_export_switches_complete_sets_and_preserves_user_files(tmp_path, monkeypatch):
    from mfrm_app import project
    opts = dict(facet_a="Rater", facet_b="Task", omit_extreme=False)
    a = mfrmpy.analyze(ratings(), **settings(), bias_requests=[opts])
    b = mfrmpy.analyze(ratings(), **settings())
    old = mfrmpy.export_csv(a, tmp_path)
    manifest = old[0].read_bytes()
    mine = tmp_path / "user.csv"
    mine.write_text("my work")
    write = project.write_bytes_atomic

    def fail_manifest(path, data):
        if path == old[0]:
            raise OSError("simulated interrupted write")
        write(path, data)

    with monkeypatch.context() as patch:
        patch.setattr(project, "write_bytes_atomic", fail_manifest)
        with pytest.raises(OSError, match="interrupted"):
            mfrmpy.export_csv(b, tmp_path)
    assert old[0].read_bytes() == manifest and all(p.exists() for p in old)
    latest = mfrmpy.export_csv(b, tmp_path)
    current = pd.read_csv(latest[0])
    assert all((tmp_path / name).is_file() for name in current.File)
    assert current.ContentFingerprint.eq(b["content_fingerprint"]).all()
    assert all(not p.exists() for p in old[1:]) and mine.read_text() == "my work"


@pytest.mark.figures
def test_rater_profile_respects_positive_direction_and_facet_role():
    pytest.importorskip("matplotlib")
    an = mfrmpy.analyze(ratings(), **settings(), positive_facets=["Rater"])
    fig = mfrmpy.plot(an, "rater_severity")
    assert "leniency" in fig._suptitle.get_text()
    assert "raise" in fig.axes[0].get_xlabel()
    assert "Rater" not in mfrmpy.plot(an, "rater_severity", facet="Task")._suptitle.get_text()
