"""Installed public workflows; runnable independently of the application checkout."""
import importlib.util
import subprocess
import sys

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
analysis = mfrmpy.analyze(data, person_col='Person', facet_cols=['Rater', 'Task', 'Criterion'],
                         score_col='Score', rating_min=0, rating_max=5)
assert mfrmpy.tables(analysis)
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
            assert frame.Use.eq("Withheld").all()
        if title.startswith("Specified intervals:"):
            assert frame.Use.eq("Approximate").all()
            assert pd.to_numeric(frame.SE).gt(0).all()
            assert "finite-sample" in note
    assert "Separation and Reliability by Facet" in {title for title, _, _ in after}


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


@pytest.mark.skipif(importlib.util.find_spec("matplotlib") is None, reason="Optional figures extra")
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
