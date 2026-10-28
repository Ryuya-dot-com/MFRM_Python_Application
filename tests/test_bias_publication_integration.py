"""Public bias exports keep the G3 hold; technical records remain unchanged.

Legacy GPCM tables may contain finite reference probabilities and intervals.
They must not be promoted through Word/PDF/HTML table helpers. Bibliography
registration is tested separately from statistical qualification.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


# -----------------------------------------------------------------------------
# Unit: _bias_table_for_publication
# -----------------------------------------------------------------------------


def _synthetic_bias_table(model: str) -> pd.DataFrame:
    """Build a tiny bias table in the same shape ``estimate_bias_interaction``
    produces, for use in unit tests of the publication helpers."""
    rows = [
        {
            "Sq": 1,
            "FacetA": "Rater",
            "FacetA_Level": "R1",
            "FacetB": "Criterion",
            "FacetB_Level": "C1",
            "Bias Size": 0.12,
            "S.E.": 0.18,
            "t": 0.67,
            "d.f.": 30.0,
            "Prob.": 0.51,
            "Infit": 1.0,
            "Outfit": 1.0,
            "ObsN": 30,
        },
        {
            "Sq": 2,
            "FacetA": "Rater",
            "FacetA_Level": "R2",
            "FacetB": "Criterion",
            "FacetB_Level": "C1",
            "Bias Size": -0.15,
            "S.E.": 0.20,
            "t": -0.75,
            "d.f.": 30.0,
            "Prob.": 0.46,
            "Infit": 1.0,
            "Outfit": 1.0,
            "ObsN": 30,
        },
    ]
    df = pd.DataFrame(rows)
    if model == "GPCM":
        df["LR ChiSq"] = [0.45, 0.58]
        df["LR Prob."] = [0.502, 0.446]
        df["Profile CI Lower"] = [-0.24, -0.55]
        df["Profile CI Upper"] = [0.48, 0.25]
        df["Profile CI Level"] = [0.95, 0.95]
        df["Profile CI Status"] = ["ok", "ok"]
        df["Likelihood Basis"] = ["conditional profile likelihood" * 2] * 2
        df["InferenceTier"] = ["screening", "screening"]
    return df


def test_bias_table_for_publication_returns_none_for_empty_input():
    assert app._bias_table_for_publication(None, "GPCM") is None
    assert app._bias_table_for_publication({}, "GPCM") is None


def test_bias_table_for_publication_returns_none_when_all_results_are_skip_reasons():
    bundle = {
        "Rater x Criterion": {"_skip_reason": "no usable cells"},
        "Rater x Task": {"_skip_reason": "extreme levels excluded everything"},
    }
    out = app._bias_table_for_publication(bundle, "GPCM")
    assert out is None


def test_bias_table_for_publication_withholds_unqualified_inference_even_when_finite():
    bundle = {"Rater x Criterion": {"table": _synthetic_bias_table("GPCM")}}
    original = bundle["Rater x Criterion"]["table"].copy(deep=True)
    out = app._bias_table_for_publication(bundle, "GPCM")
    assert {"Pair", "Level A", "Level B", "Bias Size", "Calculation", "Use"}.issubset(out.columns)
    assert not {"S.E.", "t", "d.f.", "Prob.", "LR Prob.", "Profile CI Lower"} & set(out.columns)
    assert out["Use"].eq("Withheld").all()
    assert out["Calculation"].eq("Not verified").all()
    assert not {"QualificationStatus", "PublicConclusionAllowed"} & set(out.columns)
    pd.testing.assert_frame_equal(original, bundle["Rater x Criterion"]["table"])


@pytest.mark.parametrize("status,label", [
    ("interior", "Interior solution"), ("excluded_extreme_level", "Extreme level excluded"),
    ("no_finite_mle", "No finite estimate"), ("unknown_status", "Needs review"),
])
def test_bias_publication_calculation_labels_never_release_hold(status, label):
    table = _synthetic_bias_table("RSM").assign(**{"Cell Status": status})
    out = app._bias_table_for_publication({"Rater x Criterion": {"table": table}}, "RSM")
    assert out["Calculation"].eq(label).all()
    assert out["Use"].eq("Withheld").all()
    np.testing.assert_array_equal(out["Bias Size"], table["Bias Size"])


def test_bias_table_for_publication_omits_gpcm_columns_when_model_is_not_gpcm():
    bundle = {"Rater x Criterion": {"table": _synthetic_bias_table("RSM")}}
    out = app._bias_table_for_publication(bundle, "RSM")
    assert out is not None
    for col in ("LR ChiSq", "LR Prob.", "Profile CI Lower", "Profile CI Upper", "Profile CI Status"):
        assert col not in out.columns


def test_bias_table_for_publication_omits_gpcm_columns_when_lr_chisq_is_all_nan():
    """Even under model == GPCM, if the LR column is uniformly NaN
    (e.g. estimation hit an error) the helper keeps the G3 hold."""
    df = _synthetic_bias_table("GPCM").copy()
    df["LR ChiSq"] = np.nan
    bundle = {"Rater x Criterion": {"table": df}}
    out = app._bias_table_for_publication(bundle, "GPCM")
    assert out is not None
    assert "LR ChiSq" not in out.columns


def test_bias_table_for_publication_concatenates_multiple_facet_pairs():
    bundle = {
        "Rater x Criterion": {"table": _synthetic_bias_table("GPCM")},
        "Rater x Task": {"table": _synthetic_bias_table("GPCM")},
    }
    out = app._bias_table_for_publication(bundle, "GPCM")
    assert out is not None
    assert set(out["Pair"]) == {"Rater x Criterion", "Rater x Task"}
    assert len(out) == 4  # 2 rows per pair, 2 pairs


# -----------------------------------------------------------------------------
# Unit: _bias_table_caption
# -----------------------------------------------------------------------------


@pytest.mark.parametrize("model", ["RSM", "PCM", "GPCM", "gpcm", "", None])
def test_bias_table_caption_explains_hold_for_every_model(model):
    caption = app._bias_table_caption(model, table_idx=3)
    assert "Table 3." in caption
    assert "Validation of tests and intervals is pending" in caption
    assert "G3" not in caption
    assert "calibration held fixed" in caption
    assert "technical exports" in caption
    assert "do not establish bias or absence of bias" in caption


# -----------------------------------------------------------------------------
# Unit: APA reference library registers the new citations
# -----------------------------------------------------------------------------


def test_apa_reference_library_includes_new_inference_citations():
    library = app._APA_REFERENCE_LIBRARY
    for key in ("Wilks_1938", "Cox_1975", "Louis_1982", "Cramer_1946", "Muraki_1993"):
        assert key in library, f"{key} missing from APA reference library"
        # Each entry must be a non-empty string with the year visible so
        # the rendered reference is human-readable.
        assert isinstance(library[key], str) and len(library[key]) > 20
    assert "1938" in library["Wilks_1938"]
    assert "1975" in library["Cox_1975"]
    assert "1982" in library["Louis_1982"]
    assert "1946" in library["Cramer_1946"]
    assert "1993" in library["Muraki_1993"]


def test_citation_map_resolves_new_inference_citations():
    citation_map = app._CITATION_TO_KEY
    for cite, key in (
        ("(Wilks, 1938)", "Wilks_1938"),
        ("(Cox, 1975)", "Cox_1975"),
        ("(Louis, 1982)", "Louis_1982"),
        ("(Cramer, 1946)", "Cramer_1946"),
        ("(Muraki, 1993)", "Muraki_1993"),
    ):
        assert citation_map.get(cite) == key


def test_build_apa_reference_list_emits_new_citations_when_text_mentions_them():
    """A Methods narrative that names the new citations in `(Author, Year)`
    form must produce a references list that includes their full entries."""
    text = (
        "Bias inference uses the likelihood-ratio test (Wilks, 1938) and "
        "the profile-likelihood confidence interval (Cox, 1975). Standard "
        "errors follow the delta method (Cramer, 1946) on the marginal "
        "likelihood's observed information (Louis, 1982). The slope-aware "
        "information identity follows (Muraki, 1993)."
    )
    refs = app.build_apa_reference_list(text)
    rendered = "\n".join(refs)
    for marker in ("Wilks", "Cox", "Cramer", "Louis", "Muraki"):
        assert marker in rendered, f"{marker} missing from rendered references"
    # APA 7 references list is alphabetised by first-author surname.
    surnames = [r.split(",")[0] for r in refs]
    assert surnames == sorted(surnames), (
        f"references are not alphabetised: {surnames!r}"
    )
