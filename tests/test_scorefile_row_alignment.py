"""Score probabilities must remain attached to the observations they describe."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from desktop_app import demonstration_data
from mfrm_app import basic_analysis, engine, project


@pytest.fixture(scope="module")
def fitted():
    return engine.mfrm_estimate(
        demonstration_data(), person_col="Person", facet_cols=["Rater", "Task"],
        score_col="Score", rating_min=0, rating_max=3, model="RSM", method="JMLE",
    )


@pytest.mark.parametrize("index_kind", ["permuted", "gapped", "duplicate", "text"])
def test_scorefile_probabilities_follow_positions_not_index_labels(fitted, index_kind):
    expected = engine.compute_scorefile(fitted).reset_index(drop=True)
    result = deepcopy(fitted)
    n = len(result["prep"]["data"])
    indices = {
        "permuted": np.arange(n)[::-1],
        "gapped": np.arange(n) * 2 + 100,
        "duplicate": np.zeros(n, dtype=int),
        "text": [f"source-{i}" for i in range(n)],
    }
    result["prep"]["data"].index = pd.Index(indices[index_kind], name="source_row")
    actual = engine.compute_scorefile(result)
    pd.testing.assert_index_equal(actual.index, result["prep"]["data"].index)
    pd.testing.assert_frame_equal(actual.reset_index(drop=True), expected)


def test_desktop_dropped_score_row_survives_save_restore(tmp_path):
    data = demonstration_data()
    # A missing score leaves a gap in the retained row labels, even when the
    # imported CSV originally had an ordinary RangeIndex.
    data.loc[1, "Score"] = np.nan
    analysis = basic_analysis.run_basic_analysis(
        data, person_col="Person", facet_cols=["Rater", "Task"], score_col="Score",
        rating_min=0, rating_max=3,
    )
    scorefile = analysis["scorefile"]
    assert len(scorefile) == len(data) - 1
    assert 1 not in scorefile.index
    probs = engine.compute_prob_matrix(analysis["result"])
    np.testing.assert_allclose(scorefile[["P_0", "P_1", "P_2", "P_3"]], probs)
    np.testing.assert_array_equal(scorefile["MostLikely"], np.argmax(probs, axis=1))
    expected = (probs * np.arange(4)).sum(axis=1)
    np.testing.assert_allclose(scorefile["Expected"], expected)
    path = tmp_path / "missing-score.mfrm"
    project.save_project(path, analysis)
    restored = project.load_project(path)
    pd.testing.assert_frame_equal(restored["scorefile"], scorefile)
