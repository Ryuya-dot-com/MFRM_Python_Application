import pandas as pd
import pytest

from validation.known_assignment_multivector_visible_join import (
    _hash_mapping,
    _numeric_frame_difference,
)


def test_numeric_frame_difference_is_keyed_and_nan_aware():
    left = pd.DataFrame({
        "RunId": ["b", "a"],
        "Facet": ["Rater", "Rater"],
        "Level": ["R02", "R01"],
        "Estimate": [float("nan"), 0.25],
    })
    right = pd.DataFrame({
        "RunId": ["a", "b"],
        "Facet": ["Rater", "Rater"],
        "Level": ["R01", "R02"],
        "Estimate": [0.2500000000005, float("nan")],
    })
    difference, rows = _numeric_frame_difference(
        left,
        right,
        keys=["RunId", "Facet", "Level"],
        values=["Estimate"],
    )
    assert rows == 2
    assert difference == pytest.approx(5e-13, rel=0, abs=1e-16)


def test_numeric_frame_difference_rejects_missing_keys():
    left = pd.DataFrame({"Key": ["a"], "Value": [1.0]})
    right = pd.DataFrame({"Key": ["b"], "Value": [1.0]})
    with pytest.raises(ValueError, match="keys differ"):
        _numeric_frame_difference(left, right, keys=["Key"], values=["Value"])


def test_marker_mapping_digest_is_order_independent():
    assert _hash_mapping({"b": "2", "a": "1"}) == _hash_mapping(
        {"a": "1", "b": "2"}
    )
