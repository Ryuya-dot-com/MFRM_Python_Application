"""Whole-vector and partitioned independent integrals retain the same target."""
import numpy as np
import pandas as pd
import pytest

from validation.mml_rasch_integral import RaschIntegral
from validation.mml_rasch_block_reference import BlockRaschIntegral


def problem():
    config = {"model": "RSM", "noncenter_facet": "Task", "facet_names": ["Rater", "Task"],
              "facet_levels": {"Rater": ["R0", "R1", "R2"], "Task": ["T0", "T1"]},
              "facet_signs": {"Rater": 1., "Task": -1.}, "n_cat": 4, "rating_min": 0}
    rows = [(f"P{p:02}", f"R{r}", f"T{t}", (p + 2*r + t) % 4)
            for p in range(13) for r in (p % 3, (p + 1) % 3) for t in range(2)]
    data = pd.DataFrame(rows, columns=["Person", "Rater", "Task", "Score"])
    point = np.array([.1, -.2, -.15, .2, -.8, .1, 0.])
    return data.iloc[::-1], config, point


@pytest.mark.parametrize("information", [False, True])
def test_partition_retains_global_coordinates_person_order_and_tail_bounds(information):
    data, config, point = problem()
    before = data.copy(deep=True)
    settings = dict(bound=12., rel_tol=1e-10, abs_tol=1e-12, information=information)
    whole = RaschIntegral(data, config).continuous(point, **settings)
    split = BlockRaschIntegral(data, config).continuous(point, **settings)
    for name in ("nll", "gradient", "eap", "sd", "density_log_sigma_score"):
        np.testing.assert_allclose(split[name], whole[name], atol=1e-10, rtol=1e-11)
    if information:
        np.testing.assert_allclose(split["information"], whole["information"], atol=1e-10, rtol=1e-11)
        assert split["information"].shape == (7, 7)  # Full log-SD coordinate remains present.
    assert split["integration_blocks"] == 2  # Includes a partial final block.
    assert split["person_labels"] == sorted(data.Person.unique())
    assert split["tail_gradient_abs_bound"] >= whole["tail_gradient_abs_bound"] * (1 - 1e-12)
    assert split["numeric_gradient_error_allowance"] >= 0
    pd.testing.assert_frame_equal(data, before)


def test_partition_preserves_reference_rejection_of_nonunit_weights_and_empty_data():
    data, config, point = problem()
    data["Weight"] = 2.
    with pytest.raises(ValueError, match="unit-weight"):
        BlockRaschIntegral(data, config)
    with pytest.raises(ValueError, match="nonempty"):
        BlockRaschIntegral(data.iloc[:0], config)
