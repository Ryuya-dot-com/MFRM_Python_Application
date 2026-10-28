"""History copy preserves restore/delete behavior and technical export holds."""
from copy import deepcopy
from io import BytesIO
import json

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_help_contextual_app_integration import _block_analysis_calls, _locale_value


@pytest.fixture(scope="module")
def history_output():
    result = app.mfrm_estimate(app._make_self_test_gradient_data(), "Person", ["Rater", "Task"], "Score",
        method="MML", mml_engine="EM", estimate_population_sd=True, quad_points=15, maxit=100)
    return {"result": result, "diagnostics": app.mfrm_diagnostics(result, compute_pca=False),
            "model_type": "RSM", "_est_method": "MML", "facet_cols": ["Rater", "Task"]}


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_bilingual_history_preserves_holds_restore_and_delete_confirmation(lang, history_output, monkeypatch):
    at = AppTest.from_string('''
import streamlit_app as app
if not app.st.session_state.get("_history_test_seeded"):
    app.st.session_state[app._RUN_HISTORY_STORE_SNAPSHOTS_KEY] = True
    app.record_run_in_history(output=app.st.session_state["seed_output"], elapsed_sec=1.2)
    app.st.session_state[app._RUN_HISTORY_STORE_SNAPSHOTS_KEY] = False
    app.record_run_in_history(output=app.st.session_state["seed_output"], elapsed_sec=1.3)
    app.st.session_state["_history_test_seeded"] = True
app.render_run_history_panel()
''', default_timeout=120)
    at.session_state["lang"] = lang
    at.session_state["seed_output"] = deepcopy(history_output)
    calls = _block_analysis_calls(monkeypatch)
    downloads = {}
    original = app.st.download_button

    def capture(label, data, *args, **kwargs):
        downloads[kwargs.get("file_name")] = data
        return original(label, data, *args, **kwargs)

    monkeypatch.setattr(app.st, "download_button", capture)

    def run():
        at.run(timeout=120)
        assert not at.exception

    run()
    full, compact = at.session_state[app._RUN_HISTORY_KEY]
    saved = deepcopy(full["output_snapshot"])
    identity = app.build_result_analysis_identity(saved["result"]).analysis_id
    assert at.expander[0].label == _locale_value(lang, "history.title").format(n=2)
    assert at.checkbox(key=f"{app._RUN_HISTORY_STORE_SNAPSHOTS_KEY}_widget").label == _locale_value(lang, "history.snapshot_label")
    assert _locale_value(lang, "history.criteria_hold") in [w.value for w in at.info]
    assert not at.button(key=f"restore_{full['run_id']}").disabled
    assert at.button(key=f"restore_{compact['run_id']}").disabled
    assert at.button(key=f"restore_{compact['run_id']}").proto.help == _locale_value(lang, "history.no_snapshot_help")
    assert json.loads(at.dataframe[0].proto.columns)["timestamp"]["label"] == _locale_value(lang, "history.columns.time")
    for table in at.dataframe:
        assert not {"QualificationStatus", "ReasonCode", "PreferredDirection"}.intersection(table.value.columns)
    visible = " ".join(w.value for kind in (at.caption, at.markdown, at.info, at.warning) for w in kind)
    assert "G2" not in visible and "TECHNICAL_ONLY" not in visible
    csv = pd.read_csv(BytesIO(downloads["mfrm_run_history_likelihood_criteria.csv"]))
    assert csv.QualificationStatus.eq("TECHNICAL_ONLY").all()
    assert csv.ReasonCode.eq("model.model_001.selection_scope_unvalidated").all()
    assert not csv.PublicConclusionAllowed.any()

    at.button(key=f"restore_{full['run_id']}").click()
    run()
    restored = at.session_state["facets_mode_output"]
    assert app.build_result_analysis_identity(restored["result"]).analysis_id == identity
    assert restored is not full["output_snapshot"]
    pd.testing.assert_frame_equal(restored["diagnostics"]["measures"], saved["diagnostics"]["measures"], check_exact=True)
    at.button(key="run_history_clear").click()
    run()
    assert _locale_value(lang, "history.clear_warning").format(n=2) in [w.value for w in at.warning]
    assert len(at.session_state[app._RUN_HISTORY_KEY]) == 2
    at.button(key="run_history_clear_no").click()
    run()
    assert len(at.session_state[app._RUN_HISTORY_KEY]) == 2
    pd.testing.assert_frame_equal(full["output_snapshot"]["diagnostics"]["measures"], saved["diagnostics"]["measures"], check_exact=True)
    at.button(key="run_history_clear").click()
    run()
    at.button(key="run_history_clear_yes").click()
    run()
    assert app._RUN_HISTORY_KEY not in at.session_state
    assert at.session_state["facets_mode_output"] is restored
    assert calls == {"estimate": 0, "refit": 0}


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_snapshot_preference_survives_navigation_and_opt_out(lang, history_output):
    at = AppTest.from_string('''
import streamlit_app as app
st = app.st
if not st.session_state.get("seeded"):
    app.record_run_in_history(output=st.session_state["seed_output"], elapsed_sec=1.)
    st.session_state["seeded"] = True
panel = st.radio("Panel", ["history", "other"], key="panel")
if st.button("Record", key="record"):
    app.record_run_in_history(output=st.session_state["seed_output"], elapsed_sec=1.)
if panel == "history":
    app.render_run_history_panel()
''', default_timeout=60)
    at.session_state["lang"] = lang
    at.session_state["seed_output"] = deepcopy(history_output)
    at.run()
    preference = next(c for c in at.checkbox if c.label == _locale_value(lang, "history.snapshot_label"))
    assert not preference.value
    preference.check().run()
    at.radio(key="panel").set_value("other").run()
    at.button(key="record").click().run()
    assert not at.exception
    entry = at.session_state[app._RUN_HISTORY_KEY][-1]
    assert entry["snapshot_available"]
    at.radio(key="panel").set_value("history").run()
    preference = next(c for c in at.checkbox if c.label == _locale_value(lang, "history.snapshot_label"))
    assert preference.value
    at.button(key=f"restore_{entry['run_id']}").click().run()
    assert not at.exception
    pd.testing.assert_frame_equal(at.session_state["facets_mode_output"]["diagnostics"]["measures"],
                                  history_output["diagnostics"]["measures"])
    next(c for c in at.checkbox if c.label == preference.label).uncheck().run()
    at.radio(key="panel").set_value("other").run()
    at.button(key="record").click().run()
    assert not at.session_state[app._RUN_HISTORY_KEY][-1]["snapshot_available"]
