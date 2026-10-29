"""Extended project contracts: stored requests, no refit, and qualified replay."""
from copy import deepcopy
from collections import OrderedDict
import json
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import pytest

from desktop_app import AnalysisWorker, DesktopApp
from mfrm_app import basic_analysis as basic, desktop_analysis as desktop, engine, project
from mfrm_app import bias_analysis as bias, generalizability as gd, structural_intervals as intervals
from tests.test_g_d_study_boundaries import ratings


def test_ordered_parameter_sizes_keep_type_order_and_reject_duplicate_keys():
    sizes = OrderedDict([("steps", 2), ("rater", 3)])
    restored = project._decode(json.loads(json.dumps(project._encode(sizes))))
    assert type(restored) is OrderedDict and list(restored.items()) == list(sizes.items())
    with pytest.raises(project.ProjectError, match="Duplicate"):
        project._decode(["ordered-dict", [["steps", 2], ["steps", 3]]])


def test_parameter_slices_are_typed_data_and_never_dynamic_objects():
    value = OrderedDict([("empty", slice(0, 0)), ("steps", slice(2, 6))])
    assert project._decode(project._encode(value)) == value
    for malformed in (["slice", "0", 1, None], ["slice", 0, 1, 0], ["slice", False, 1, None]):
        with pytest.raises(project.ProjectError, match="slice"):
            project._decode(malformed)


def test_legacy_schema_one_still_loads_and_new_writes_are_version_two(tmp_path, monkeypatch):
    from desktop_app import demonstration_data
    original = basic.run_basic_analysis(demonstration_data(), person_col="Person", facet_cols=["Rater", "Task"],
                                       score_col="Score", rating_min=0, rating_max=3)
    path = tmp_path / "schema1.mfrm"
    path.write_bytes(project._json_bytes(dict(format=project.FORMAT, schema_version=1, analysis=project._encode(original))))
    _forbid_calculation(monkeypatch)
    restored = project.load_project(path)
    assert restored["content_fingerprint"] == original["content_fingerprint"]
    project.save_project(tmp_path / "schema2.mfrm", restored)
    assert json.loads((tmp_path / "schema2.mfrm").read_text(encoding='utf-8'))["schema_version"] == 2
    assert json.loads(path.read_text(encoding='utf-8'))["schema_version"] == 1


@pytest.fixture(scope="module", params=[("RSM", "JMLE"), ("PCM", "JMLE"), ("RSM", "MML"), ("PCM", "MML")])
def analysis(request):
    model, method = request.param
    data = pd.DataFrame([dict(Person=f"{p:03d}", Rater=f"R{r}", Task=f"T{t}", Score=(p+2*r+t)%3)
                        for p in range(30) for r in range(3) for t in range(2)])
    return desktop.run_analysis(data, person_col="Person", facet_cols=["Rater", "Task"], score_col="Score",
        rating_min=0, rating_max=2, model=model, method=method, maxit=100,
        step_facet="Task" if model == "PCM" else None, quad_points=61 if method == "MML" else 15,
        bias_requests=[dict(facet_a="Rater", facet_b="Task", omit_extreme=False)],
        interval_requests=[dict(facet="Rater", difference=["R0", "R2"], level=.9)] if method == "MML" else [],
        gd_settings=dict(object_facet="Person", random_facets=["Rater", "Task"], target_population="Declared population"))


def _forbid_calculation(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("Restoring or displaying a project must not recompute its analysis")
    for module, names in ((engine, ["mfrm_estimate", "mfrm_diagnostics", "mfrm_loglik_mml_value_grad"]),
                          (bias, ["estimate_bias_interaction"]),
                          (gd, ["compute_generalizability_from_data", "_crossed_anova_variance_components"]),
                          (intervals, ["compute_structural_intervals"])):
        for name in names:
            monkeypatch.setattr(module, name, forbidden)


def test_all_model_method_projects_roundtrip_requests_results_and_current_guards(analysis, tmp_path, monkeypatch):
    path = tmp_path / "分析と条件.mfrm"
    before = desktop.current_additional_tables(analysis)
    _forbid_calculation(monkeypatch)
    project.save_project(path, analysis)
    assert json.loads(path.read_text(encoding='utf-8'))["schema_version"] == 2
    restored = project.load_project(path)
    assert restored["kind"] == "mfrm-analysis-v2"
    assert restored["fit_settings"] == analysis["fit_settings"]
    assert restored["qualification_at_save"] == analysis["qualification_at_save"]
    assert restored["content_fingerprint"] == analysis["content_fingerprint"]
    pd.testing.assert_frame_equal(restored["input_data"], analysis["input_data"])
    assert restored["input_data"].Person.iloc[0] == "000"
    np.testing.assert_array_equal(restored["result"]["opt"].x, analysis["result"]["opt"].x)
    assert "hess_inv" not in restored["result"]["opt"]
    for (name, table, note), (other_name, other_table, other_note) in zip(before, desktop.current_additional_tables(restored)):
        assert (name, note) == (other_name, other_note)
        pd.testing.assert_frame_equal(table, other_table)
        if name.startswith("Bias:"):
            assert table.Use.eq("Withheld").all()
            assert "LR Prob." not in table and "S.E." not in table
    for entry in restored["interval_analyses"]:
        assert entry["settings"] == dict(facet="Rater", difference=["R0", "R2"], level=.9)
        assert entry["status"] == "computed"
        table = intervals.read_structural_intervals(restored["result"], entry["record"])
        assert table.CIEligible.all() and not table.CoverageQualified.any()


def test_schema_and_foreign_binding_are_checked_even_with_recomputed_outer_hash(analysis, tmp_path):
    mutations = [lambda a: a["fit_settings"].update(model="GPCM"),
                 lambda a: a["bias_analyses"][0].update(analysis_id="different_calibration"),
                 lambda a: a["gd_analysis"]["settings"].update(target_population="Different population")]
    if analysis["interval_analyses"]:
        mutations += [lambda a: a["interval_analyses"][0].update(analysis_id="different_calibration"),
                      lambda a: a["interval_analyses"][0]["settings"].update(level=.95)]
    for mutate in mutations:
        other = deepcopy(analysis); mutate(other)
        other["content_fingerprint"] = project.content_fingerprint(other)
        with pytest.raises(project.ProjectError):
            project.save_project(tmp_path / "foreign.mfrm", other)
    path = tmp_path / "downgraded.mfrm"
    project.save_project(path, analysis)
    document = json.loads(path.read_text(encoding='utf-8')); document["schema_version"] = 1
    path.write_bytes(project._json_bytes(document))
    with pytest.raises(project.ProjectError, match="Schema 1"):
        project.load_project(path)


def test_interval_source_change_preserves_archive_but_withholds_current_output(analysis, tmp_path, monkeypatch):
    if not analysis["interval_analyses"]:
        return
    path = tmp_path / "old_source.mfrm"; project.save_project(path, analysis)
    monkeypatch.setattr(intervals, "_source_hash", lambda: "different-numerical-source")
    _forbid_calculation(monkeypatch)
    restored = project.load_project(path)
    selected = [table for name, table, _ in desktop.current_additional_tables(restored) if name.startswith("Specified intervals:")]
    assert selected and all(table.Use.eq("Withheld").all() for table in selected)
    assert all("SE" not in table for table in selected)
    assert restored["content_fingerprint"] == analysis["content_fingerprint"]


def test_interval_review_failure_preserves_the_fitted_result_and_reason(analysis, tmp_path, monkeypatch):
    if not analysis["interval_analyses"]:
        return
    monkeypatch.setattr(engine, "mfrm_estimate", lambda *a, **k: deepcopy(analysis["result"]))
    def unsupported(*a, **k):
        raise ValueError("Controlled numerical review limit")
    monkeypatch.setattr(intervals, "compute_structural_intervals", unsupported)
    out = desktop.run_analysis(analysis["input_data"], **analysis["fit_settings"],
                               interval_requests=[dict(facet="Rater")])
    assert out["analysis_id"] == analysis["analysis_id"]
    entry = out["interval_analyses"][0]
    assert entry["status"] == "unavailable" and entry["record"] is None
    path = tmp_path / "review_unavailable.mfrm"; project.save_project(path, out)
    table = desktop.current_additional_tables(project.load_project(path))[0][1]
    assert table.Use.eq("Withheld").all() and table.Reason.eq("Controlled numerical review limit").all()


def test_basic_report_cannot_mislabel_an_extended_project_as_rsm_jmle(analysis):
    with pytest.raises(ValueError, match="original RSM/JMLE"):
        basic.basic_report_html(analysis)


@pytest.mark.parametrize("two_facets", [False, True])
def test_standalone_gd_preserves_original_scores_selected_plan_and_reasons(two_facets, tmp_path, monkeypatch):
    data = ratings(two_facets)
    facets = ["Rater", "Item"] if two_facets else ["Rater"]
    settings = dict(object_facet="Person", random_facets=facets, selected_plan={f: 4 for f in facets},
                    target_population="Students and comparable ratings")
    monkeypatch.setattr(engine, "mfrm_estimate", lambda *a, **k: pytest.fail("G/D must not fit MFRM"))
    good = desktop.run_analysis(data, workflow="gd", gd_settings=settings)
    missing = desktop.run_analysis(data.iloc[1:], workflow="gd", gd_settings=settings)
    assert good["gd_analysis"]["result"]["available"]
    assert missing["gd_analysis"]["result"]["reason_code"] == "missing_cells"
    _forbid_calculation(monkeypatch)
    for index, original in enumerate((good, missing)):
        path = tmp_path / f"gd_{index}.mfrm"; project.save_project(path, original)
        restored = project.load_project(path)
        pd.testing.assert_frame_equal(restored["input_data"], original["input_data"])
        assert restored["gd_analysis"]["settings"]["selected_plan"] == settings["selected_plan"]
        assert restored["content_fingerprint"] == original["content_fingerprint"]
        assert desktop.current_additional_tables(restored)


def test_gd_reopens_in_fresh_process_with_no_web_or_calculation(tmp_path):
    path = tmp_path / "G研究.mfrm"
    project.save_project(path, desktop.run_analysis(ratings(), workflow="gd", gd_settings=dict(
        object_facet="Person", random_facets=["Rater", "Item"], selected_plan={"Rater": 4, "Item": 2})))
    script = '''
import importlib.abc, sys
class NoWeb(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'streamlit', 'streamlit_app'}:
            raise AssertionError('No Web imports')
sys.meta_path.insert(0, NoWeb())
from mfrm_app import project, desktop_analysis, engine, generalizability
def forbidden(*a, **k): raise AssertionError('No recalculation on reopen')
engine.mfrm_estimate = engine.mfrm_diagnostics = forbidden
generalizability.compute_generalizability_from_data = forbidden
generalizability._crossed_anova_variance_components = forbidden
a = project.load_project(sys.argv[1])
assert a['gd_analysis']['settings']['selected_plan'] == {'Rater': 4, 'Item': 2}
assert desktop_analysis.current_additional_tables(a)
'''
    result = subprocess.run([sys.executable, "-c", script, str(path)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("extra", [dict(model="PCM"), dict(model="GPCM"), dict(method="MML", quad_points=True),
    dict(method="MML", population_prior_sd=0), dict(estimate_population_sd=True),
    dict(bias_requests=[dict(facet_a="Rater", facet_b="Rater")]),
    dict(interval_requests=[dict(facet="Rater")])])
def test_invalid_requests_fail_before_estimation(extra, monkeypatch):
    monkeypatch.setattr(engine, "mfrm_estimate", lambda *a, **k: pytest.fail("Validate before fitting"))
    with pytest.raises(ValueError):
        desktop.run_analysis(ratings(False), person_col="Person", facet_cols=["Rater"], score_col="Score",
                             rating_min=0, rating_max=40, **extra)


def test_worker_cancellation_failure_and_gd_restart_use_the_existing_transport(monkeypatch):
    from pathlib import Path
    # Streamlit AppTest can leave __main__ pointing at its temporary view script.
    # Spawn from the real desktop entrypoint, as the application does.
    monkeypatch.setattr(sys.modules["__main__"], "__file__", str(Path(__file__).resolve().parents[1] / "desktop_app.py"))
    monkeypatch.setattr(sys.modules["__main__"], "__spec__", None)
    worker = AnalysisWorker()
    data = ratings()
    settings = dict(workflow="gd", gd_settings=dict(object_facet="Person", random_facets=["Rater", "Item"]))
    def drain():
        events = []
        while worker.process is not None:
            events.extend(worker.poll())
            if worker.process is not None:
                time.sleep(.02)
        return events[-1]
    try:
        worker.start(data, settings); path = worker.output_path
        assert worker.cancel()[0] == "cancelled" and not path.exists()
        worker.start(data, dict(workflow="unknown"))
        assert drain()[0] == "failed"
        worker.start(data, settings)
        status, output = drain()
        assert status == "complete", output
        assert output["kind"] == "gd-analysis-v2" and output["gd_analysis"]["result"]["available"]
        from desktop_app import demonstration_data
        worker.start(demonstration_data(), dict(person_col="Person", facet_cols=["Rater", "Task"],
            score_col="Score", rating_min=0, rating_max=3, model="PCM", step_facet="Task", method="MML",
            quad_points=5, population_prior_sd=1.25, maxit=3,
            bias_requests=[dict(facet_a="Rater", facet_b="Task")],
            interval_requests=[dict(facet="Rater")]))
        status, output = drain()
        assert status == "complete", output
        assert output["fit_settings"]["model"] == "PCM" and output["fit_settings"]["method"] == "MML"
        assert output["fit_settings"]["population_prior_sd"] == output["result"]["config"]["population_prior_sd"] == 1.25
        assert output["bias_analyses"][0]["analysis_id"] == output["analysis_id"]
        assert output["interval_analyses"] and desktop.current_additional_tables(output)
    finally:
        worker.close()


def test_malformed_project_does_not_replace_an_existing_project(tmp_path, monkeypatch):
    from tkinter import filedialog
    path = tmp_path / "extended.mfrm"
    path.write_text("not a project")
    ui = DesktopApp.__new__(DesktopApp)
    ui.busy = False
    ui.root, ui.analysis, ui.discard = None, {"existing": "retained"}, lambda: True
    errors = []; ui.error = errors.append
    ui.set_data = lambda *a: pytest.fail("A malformed project must not replace current input")
    monkeypatch.setattr(filedialog, "askopenfilename", lambda **k: str(path))
    ui.open_project()
    assert errors
    assert ui.analysis == {"existing": "retained"}
