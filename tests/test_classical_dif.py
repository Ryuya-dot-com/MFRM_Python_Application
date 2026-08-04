"""Tests for classical polytomous DIF screening (Mantel + ordinal logistic).

These pin the internal statistical contract — a planted DIF item is detected,
null data is not flagged, purification removes matching contamination, the
result is deterministic, and bad input is refused — plus a light
operating-characteristic sanity band. A full cross-implementation equivalence
check against difR / lordif is out of scope (those packages are not in this
environment); the contract here is internal correctness, not numeric parity.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import streamlit_app as app


def _two_group_data(dif_shift=None, n_person=240, seed=2024):
    persons = [f"P{i:03d}" for i in range(n_person)]
    groups = {p: ("B" if i % 2 else "A") for i, p in enumerate(persons)}
    params = dict(
        persons=persons, raters=["R1", "R2", "R3"], tasks=["T1"],
        criteria=["C1", "C2", "C3", "C4", "C5"], theta_sd=1.2,
        rater_severities=np.array([-0.3, 0.0, 0.3]), task_difficulties=np.array([0.0]),
        criterion_difficulties=np.array([-0.4, -0.2, 0.0, 0.2, 0.4]),
        tau=np.array([-1.5, -0.5, 0.5, 1.5]), group_assignment=groups,
    )
    if dif_shift:
        params["dif_shift"] = dif_shift
    df = app._generate_mfrm_rsm_from_params(params, seed=seed)
    return df, pd.Series(groups)


def test_planted_dif_item_is_flagged_and_directional():
    df, gl = _two_group_data({("C3", "B"): 0.9})
    res = app.analyze_classical_dif(df, "Criterion", gl, alpha=0.05, min_group_n=20)
    tbl = res["table"]
    c3 = tbl[tbl["Item"] == "C3"].iloc[0]
    assert bool(c3["Flag"])
    assert c3["EvidenceLevel"] == "Strong review"
    assert str(c3["Direction"]) == "favors reference"
    assert c3["MH_p"] < 0.01 and c3["Logit_Total_p"] < 0.01


def test_purification_removes_spurious_flags():
    df, gl = _two_group_data({("C3", "B"): 0.9})
    res = app.analyze_classical_dif(df, "Criterion", gl, alpha=0.05, min_group_n=20, purify=True)
    tbl = res["table"]
    assert res["settings"]["purified_items"] == ["C3"]
    # Only the planted item is flagged once contamination is purged.
    assert int(tbl[tbl["Item"] != "C3"]["Flag"].astype(bool).sum()) == 0


def test_null_data_not_flagged():
    df, gl = _two_group_data(None)
    res = app.analyze_classical_dif(df, "Criterion", gl, alpha=0.05, min_group_n=20)
    assert int(res["table"]["Flag"].astype(bool).sum()) == 0


def test_deterministic_same_seed():
    df, gl = _two_group_data({("C3", "B"): 0.8})
    a1 = app.analyze_classical_dif(df, "Criterion", gl, alpha=0.05, min_group_n=20)["table"]
    a2 = app.analyze_classical_dif(df, "Criterion", gl, alpha=0.05, min_group_n=20)["table"]
    pd.testing.assert_frame_equal(a1.reset_index(drop=True), a2.reset_index(drop=True))


def test_skip_contract_single_group_and_few_items():
    df, gl = _two_group_data(None)
    one = pd.Series({p: "A" for p in gl.index})
    assert "_skip_reason" in app.analyze_classical_dif(df, "Criterion", one, min_group_n=20)
    # Too few groups with enough persons
    tiny = pd.Series({p: ("A" if i < 5 else "B") for i, p in enumerate(gl.index)})
    assert "_skip_reason" in app.analyze_classical_dif(df, "Criterion", tiny, min_group_n=20)


def test_mh_only_when_logistic_disabled():
    df, gl = _two_group_data({("C3", "B"): 0.9})
    res = app.analyze_classical_dif(df, "Criterion", gl, alpha=0.05, min_group_n=20, fit_logistic=False)
    tbl = res["table"]
    assert tbl["LogisticStatus"].eq("not_fitted").all()
    # MH alone still flags the planted item.
    assert bool(tbl[tbl["Item"] == "C3"].iloc[0]["Flag"])


def test_proportional_odds_cutpoints_are_ordered():
    p = np.array([-1.2, 0.1, 0.3, -0.5])
    alpha = app._po_cutpoints(p)
    assert np.all(np.diff(alpha) > 0)


def test_operating_characteristic_sanity():
    """Light OC band: planted item has high power, null items low false-flag."""
    reps = 12
    power_hits = 0
    false_flags = 0
    false_trials = 0
    for r in range(reps):
        df, gl = _two_group_data({("C3", "B"): 0.8}, n_person=160, seed=100 + r)
        tbl = app.analyze_classical_dif(df, "Criterion", gl, alpha=0.05, min_group_n=20)["table"]
        if bool(tbl[tbl["Item"] == "C3"].iloc[0]["Flag"]):
            power_hits += 1
        for it in ("C1", "C2", "C4", "C5"):
            false_trials += 1
            if bool(tbl[tbl["Item"] == it].iloc[0]["Flag"]):
                false_flags += 1
    assert power_hits / reps >= 0.7, f"low power {power_hits}/{reps}"
    assert false_flags / max(false_trials, 1) <= 0.15, f"high false-flag {false_flags}/{false_trials}"


def test_sibtest_detects_planted_dif_and_agrees_with_mantel():
    df, gl = _two_group_data({("C3", "B"): 0.9})
    tbl = app.analyze_classical_dif(df, "Criterion", gl, alpha=0.05, min_group_n=20)["table"]
    for col in ("SIBTEST_Beta", "SIBTEST_p", "SIBTEST_class", "SIBTEST_Direction", "SIBTEST_Agreement"):
        assert col in tbl.columns
    c3 = tbl[tbl["Item"] == "C3"].iloc[0]
    assert c3["SIBTEST_p"] < 0.01
    assert c3["SIBTEST_class"] in ("moderate (B)", "large (C)")
    assert c3["SIBTEST_Direction"] == "favors reference"
    assert float(c3["SIBTEST_Beta"]) > 0  # reference-minus-focal > 0 when reference favored
    assert c3["SIBTEST_Agreement"] == "agrees"  # corroborates the Mantel direction


def test_sibtest_null_no_large_effect():
    df, gl = _two_group_data(None)
    tbl = app.analyze_classical_dif(df, "Criterion", gl, alpha=0.05, min_group_n=20)["table"]
    assert not tbl["SIBTEST_class"].astype(str).isin(["moderate (B)", "large (C)"]).any()


def test_validation_bundle_contents():
    df, gl = _two_group_data({("C3", "B"): 0.9})
    res = app.analyze_classical_dif(df, "Criterion", gl, alpha=0.05, min_group_n=20)
    assert isinstance(res.get("analysis_frame"), dict)
    bundle = app.build_dif_validation_bundle(res)
    for name in ("dif_analysis_frame.csv", "dif_analysis_settings.csv", "dif_app_results.csv",
                 "run_difR_crosscheck.R", "README_difR_crosscheck.md"):
        assert name in bundle
    frame = pd.read_csv(pd.io.common.BytesIO(bundle["dif_analysis_frame.csv"]))
    assert "Person" in frame.columns and "group" in frame.columns
    assert any(c.startswith("C") for c in frame.columns)  # item columns present
    r_text = bundle["run_difR_crosscheck.R"]
    assert "difGMH" in r_text and "difSIBTEST" in r_text and "lordif" in r_text
    assert "expected" in bundle["README_difR_crosscheck.md"].lower()
