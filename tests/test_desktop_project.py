"""Desktop/project acceptance: shared calculation, typed persistence, no refit."""
from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from desktop_app import demonstration_data
from mfrm_app import basic_analysis as basic, engine, project


@pytest.fixture(scope="module")
def analysis():
    return basic.run_basic_analysis(demonstration_data(), person_col="Person", facet_cols=["Rater", "Task"],
                                    score_col="Score", rating_min=0, rating_max=3)


def test_engine_is_shared_and_imports_without_ui():
    import streamlit_app as app
    for name in ("mfrm_estimate", "mfrm_diagnostics", "calc_facets_report_tbls", "build_result_analysis_identity"):
        assert getattr(app, name) is getattr(engine, name)
    check = subprocess.run([sys.executable, "-c", "import sys; from mfrm_app import basic_analysis; "
                            "assert 'streamlit' not in sys.modules; assert 'streamlit_app' not in sys.modules"],
                           capture_output=True, text=True)
    assert check.returncode == 0, check.stderr


def test_semantic_pca_reason_still_switches_language_in_web(monkeypatch):
    import streamlit_app as app
    from types import SimpleNamespace
    state = {"lang": "ja"}
    monkeypatch.setattr(app, "st", SimpleNamespace(session_state=state))
    payload = app._make_reason("dimensionality.skip_reason_too_few_cols_template",
                               n_cols=1, key_label=app._LazyKeyLabel(mode="overall"))
    japanese = app._render_reason(payload)
    assert japanese == app._render_reason(payload, lang="ja")
    state["lang"] = "en"
    english = app._render_reason(payload)
    assert english == app._render_reason(payload, lang="en")
    assert japanese != english and "item-combination" in english


def test_basic_workflow_matches_web_entrypoint(analysis):
    import streamlit_app as app
    direct = app.mfrm_estimate(data=analysis["input_data"], **analysis["fit_settings"])
    np.testing.assert_array_equal(direct["opt"].x, analysis["result"]["opt"].x)
    assert app.build_result_analysis_identity(direct).analysis_id == analysis["analysis_id"]
    diagnostics = app.mfrm_diagnostics(direct, compute_pca=False, compute_marginal=False, compute_eb_shrinkage=False)
    for name, frame in app.calc_facets_report_tbls(direct, diagnostics).items():
        pd.testing.assert_frame_equal(frame, analysis["report_tables"][name])
    publication, _ = app._publication_measure_table(analysis["result"], analysis["diagnostics"])
    pd.testing.assert_frame_equal(publication, basic.current_tables(analysis)[1][1])


def test_new_process_restore_keeps_identity_tables_report_and_never_fits(analysis, tmp_path):
    path = tmp_path / "日本語の案件.mfrm"
    project.save_project(path, analysis)
    before = basic.basic_report_html(analysis)
    script = """
from pathlib import Path
import sys
from mfrm_app import engine, project, basic_analysis as basic
def forbidden(*args, **kwargs):
    raise AssertionError('Restoration must not fit')
engine.mfrm_estimate = forbidden
engine.mfrm_diagnostics = forbidden
a = project.load_project(sys.argv[1])
assert 'hess_inv' not in a['result']['opt']
assert a['scope']['includes_input'] and not a['scope']['optimizer_checkpoint']
Path(sys.argv[2]).write_text(basic.basic_report_html(a), encoding='utf-8')
print(a['analysis_id'])
"""
    report = tmp_path / "復元後.html"
    check = subprocess.run([sys.executable, "-c", script, str(path), str(report)], capture_output=True, text=True)
    assert check.returncode == 0, check.stderr
    assert check.stdout.strip() == analysis["analysis_id"]
    assert before == report.read_text(encoding="utf-8")
    restored = project.load_project(path)
    pd.testing.assert_frame_equal(restored["input_data"], analysis["input_data"])
    pd.testing.assert_frame_equal(restored["result"]["prep"]["data"], analysis["result"]["prep"]["data"])
    assert restored["qualification_at_save"] == analysis["qualification_at_save"]


@pytest.mark.parametrize("value", [
    pd.DataFrame({"id": pd.Series(["001", 1, None, pd.NA, np.nan, "nan"], dtype=object)}),
    pd.DataFrame({"cat": pd.Categorical(["b", None], categories=["b", "a", "unused"], ordered=True)}),
    pd.DataFrame({"x": pd.Series([1, pd.NA], dtype="Int64"), "y": pd.Series([True, pd.NA], dtype="boolean")}),
    pd.DataFrame({1: pd.Series(["01", None], dtype=pd.StringDtype(storage="python", na_value=pd.NA))}),
    pd.DataFrame([[1., np.inf], [-np.inf, np.nan]], columns=pd.Index([1, "1"], dtype=object), index=pd.Index([1, "1"], dtype=object, name="識別子")),
    np.array([[1., np.nan], [np.inf, -np.inf]], dtype=np.float32),
    np.array(["001", None, 1, pd.NA], dtype=object),
    {1: "number", "1": "string", "tuple": (1, np.int64(2)), "missing": [None, pd.NA, pd.NaT, np.nan, "nan"]},
    pd.DataFrame({"x": [1, 2]}, index=pd.MultiIndex.from_tuples([("a", 1), ("b", 2)], names=["label", "number"])),
])
def test_typed_codec_preserves_identity(value):
    restored = project._decode(json.loads(json.dumps(project._encode(value), allow_nan=False)))
    assert engine.stable_json_fingerprint(value, 64) == engine.stable_json_fingerprint(restored, 64)
    if isinstance(value, pd.DataFrame):
        pd.testing.assert_frame_equal(value, restored)


@pytest.mark.parametrize("mutation", ["input", "coordinates", "table", "settings", "qualification"])
def test_changed_or_foreign_content_is_rejected(analysis, tmp_path, mutation):
    other = deepcopy(analysis)
    if mutation == "input":
        other["result"]["prep"]["data"].loc[0, "Score"] += 1
    elif mutation == "coordinates":
        other["result"]["opt"].x[0] += .1
    elif mutation == "table":
        other["report_tables"]["Rater"].loc[0, "Measure"] += 1
    elif mutation == "settings":
        other["fit_settings"]["rating_max"] += 1
    else:
        other["qualification_at_save"]["calibration_hold_reason"] = "forged-ready"
    with pytest.raises(project.ProjectError, match="content changed"):
        project.save_project(tmp_path / "bad.mfrm", other)


def test_future_corrupt_and_nonproject_files_fail_closed(analysis, tmp_path):
    path = tmp_path / "bad.mfrm"
    project.save_project(path, analysis)
    doc = json.loads(path.read_text())
    doc["schema_version"] = 999
    path.write_text(json.dumps(doc))
    with pytest.raises(project.ProjectError, match="schema"):
        project.load_project(path)
    path.write_bytes(b"\x80\x04pickle-is-not-a-project")
    with pytest.raises(project.ProjectError):
        project.load_project(path)
    path.write_text('{"format": "mfrm.basic-project", "schema_version": 1, "analysis": ["python-object", "os.system"]}')
    with pytest.raises(project.ProjectError, match="tag"):
        project.load_project(path)


def test_missing_column_mapping_is_rejected_before_opening_in_ui(analysis, tmp_path):
    changed = deepcopy(analysis)
    changed["fit_settings"] = {}
    changed["content_fingerprint"] = project.content_fingerprint(changed)
    with pytest.raises(project.ProjectError, match="complete column"):
        project.save_project(tmp_path / "incomplete.mfrm", changed)


def test_arrays_cannot_allocate_beyond_declared_content_or_limit():
    for value in [["array", ["numpy", "<f8"], [10**15], ""],
                  ["array", ["numpy", "<f8"], [1], ""],
                  ["array", ["numpy", "|O"], [2], ["list", [1]]]]:
        with pytest.raises(project.ProjectError):
            project._decode(value)
    with pytest.raises(project.ProjectError, match="Unsupported project value"):
        project._encode(object())


def test_failed_atomic_save_preserves_existing_file(analysis, tmp_path, monkeypatch):
    path = tmp_path / "existing.mfrm"
    path.write_text("original")
    def reject(*args):
        raise PermissionError("read-only destination")
    monkeypatch.setattr(project.os, "replace", reject)
    with pytest.raises(PermissionError):
        project.save_project(path, analysis)
    assert path.read_text() == "original"
    assert list(tmp_path.iterdir()) == [path]


def test_failed_termination_survives_save_load_and_current_policy_is_reapplied(analysis, tmp_path, monkeypatch):
    failed = deepcopy(analysis)
    failed["result"]["opt"].success = False
    failed["analysis_id"] = engine.build_result_analysis_identity(failed["result"]).analysis_id
    failed["qualification_at_save"] = basic.qualification_snapshot(failed["result"])
    failed["content_fingerprint"] = project.content_fingerprint(failed)
    path = tmp_path / "failed.mfrm"
    project.save_project(path, failed)
    loaded = project.load_project(path)
    assert loaded["result"]["opt"].success is False
    assert basic.qualification_snapshot(loaded["result"])["calibration_hold_reason"]
    assert not engine.qualified_estimation_summary(loaded["result"])["InferenceReady"].any()
    public = basic.current_tables(loaded)[1][1]
    assert public["Basis"].eq("Withheld").all()
    assert public[["SE / SD", "Lower", "Upper"]].eq("").all().all()
    saved = deepcopy(loaded["qualification_at_save"])
    monkeypatch.setattr(basic.output_qualification, "calibration_inference_hold_reason", lambda result: "new-policy")
    monkeypatch.setattr(basic.output_qualification, "calibration_inference_hold_message", lambda result: "New hold reason")
    report = basic.basic_report_html(loaded)
    assert "New hold reason" in report and "conditions for interpreting this result have changed" in report
    assert loaded["qualification_at_save"] == saved


def test_report_escapes_identifiers_and_is_offline(analysis):
    copy = deepcopy(analysis)
    copy["diagnostics"]["measures"].loc[0, "Level"] = '<script>alert("x")</script>'
    copy["content_fingerprint"] = project.content_fingerprint(copy)
    report = basic.basic_report_html(copy)
    assert "<script>" not in report
    assert "&lt;script&gt;" in report
    from mfrm_app.desktop_guidance import WEB_APP_URL
    assert report.count("<svg") == 3 and 'src="http' not in report
    assert report.count(f'href="{WEB_APP_URL}"') == 1
    assert 'href="http' not in report.replace(f'href="{WEB_APP_URL}"', '')


@pytest.mark.parametrize("settings", [{"facet_cols": []}, {"rating_min": 3}, {"person_col": "absent"}])
def test_invalid_request_is_rejected_before_fitting(settings, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Must reject before fitting")
    monkeypatch.setattr(engine, "mfrm_estimate", forbidden)
    arguments = dict(person_col="Person", score_col="Score", facet_cols=["Rater", "Task"], rating_min=0, rating_max=3)
    arguments.update(settings)
    with pytest.raises(ValueError):
        basic.run_basic_analysis(demonstration_data(), **arguments)
