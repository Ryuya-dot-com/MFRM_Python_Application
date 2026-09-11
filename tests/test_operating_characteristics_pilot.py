"""Design contracts for the retained Python operating-characteristics pilot."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mfrm_app.operating_characteristics import build_replicate_manifest
from validation import operating_characteristics_pilot as pilot


def _manifest():
    return build_replicate_manifest(
        pilot.study_conditions(),
        replicates=1,
        base_seed=20260809,
    )


def _write_generated_bundle(output):
    manifest = _manifest()
    manifest.to_csv(output / "manifest.csv", index=False)
    identities = []
    for index, (_, row) in enumerate(manifest.iterrows()):
        generated = pilot.generate_condition_data(row)
        identities.append(
            pilot._append_generated_bundle(
                output,
                row,
                generated,
                first=index == 0,
            )
        )
    identity = pd.DataFrame(identities)
    inventory = pilot._finalize_generated_bundle(output, identity)
    return manifest, identity, inventory


def test_condition_matrix_pairs_null_and_alternative_for_four_designs():
    conditions = pilot.study_conditions()

    assert len(conditions) == 8
    assert {row["Design"] for row in conditions} == {
        "balanced_small",
        "balanced_large_anchors",
        "sparse_missing",
        "anchor_drift",
    }
    for design in {row["Design"] for row in conditions}:
        pair = [row for row in conditions if row["Design"] == design]
        assert {row["TruthBias"] for row in pair} == {0.0, 0.6}
        assert len({row["SeedGroup"] for row in pair}) == 1
    clean = next(row for row in conditions if row["Design"] == "balanced_large_anchors")
    drifted = next(row for row in conditions if row["Design"] == "anchor_drift")
    assert clean["SeedGroup"] == drifted["SeedGroup"]


def test_common_random_numbers_keep_design_fixed_and_change_only_responses():
    manifest = _manifest()
    null_row = manifest.loc[
        manifest["ConditionId"].eq("balanced_small__bias_0p0")
    ].iloc[0]
    alt_row = manifest.loc[
        manifest["ConditionId"].eq("balanced_small__bias_0p6")
    ].iloc[0]
    null = pilot.generate_condition_data(null_row)
    alternative = pilot.generate_condition_data(alt_row)

    assert null.data.drop(columns="Score").equals(alternative.data.drop(columns="Score"))
    assert (null.data["Score"] != alternative.data["Score"]).any()
    assert null.facet_truth.equals(alternative.facet_truth)


def test_sparse_negative_control_creates_low_count_pair_cells():
    manifest = _manifest()
    row = manifest.loc[
        manifest["ConditionId"].eq("sparse_missing__bias_0p0")
    ].iloc[0]
    generated = pilot.generate_condition_data(row)
    bundle = {
        "data": generated.data,
        "meta": {
            "facet_names": ["Rater", "Task", "Criterion"],
            "facet_level_counts": [int(row["Raters"]), int(row["Tasks"]), int(row["Criteria"])],
        },
    }
    pair_cells = pilot.app.build_custom_simulation_sparse_pair_cells(bundle)

    assert generated.realized_missing_rate > 0.20
    assert int(pair_cells["LowCountCells"].sum()) > 0


def test_anchor_drift_condition_moves_only_supplied_anchor_targets():
    manifest = _manifest()
    clean = pilot.generate_condition_data(
        manifest.loc[
            manifest["ConditionId"].eq("balanced_large_anchors__bias_0p0")
        ].iloc[0]
    )
    drifted = pilot.generate_condition_data(
        manifest.loc[
            manifest["ConditionId"].eq("anchor_drift__bias_0p0")
        ].iloc[0]
    )

    assert len(clean.anchors) == len(drifted.anchors) == 2
    assert clean.data.equals(drifted.data)
    assert clean.facet_truth.equals(drifted.facet_truth)
    clean_offsets = clean.anchors.merge(
        clean.facet_truth,
        on=["Facet", "Level"],
    )
    drift_offsets = drifted.anchors.merge(
        drifted.facet_truth,
        on=["Facet", "Level"],
    )
    assert np.allclose(clean_offsets["Anchor"] - clean_offsets["Truth"], 0.0)
    assert np.allclose(drift_offsets["Anchor"] - drift_offsets["Truth"], 0.25)


def test_parameter_recovery_uses_absolute_scale_for_anchored_facet():
    manifest = _manifest()
    row = manifest.loc[
        manifest["ConditionId"].eq("anchor_drift__bias_0p0")
    ].iloc[0]
    generated = pilot.generate_condition_data(row)
    measures = generated.facet_truth.rename(columns={"Truth": "Estimate"}).copy()
    anchor_values = generated.anchors.set_index(["Facet", "Level"])["Anchor"]
    for index, measure_row in measures.iterrows():
        key = (str(measure_row["Facet"]), str(measure_row["Level"]))
        if key in anchor_values.index:
            measures.loc[index, "Estimate"] = float(anchor_values.loc[key])
    measures["SE"] = 0.1

    recovery = pilot._parameter_rows(  # contract-level test of the runner adapter
        row,
        {"measures": measures},
        generated,
        include=True,
    )
    anchored = recovery.loc[recovery["Anchored"]]
    unanchored_task = recovery.loc[recovery["Facet"].eq("Task")]

    assert set(anchored["ComparisonScale"]) == {"anchor_identified_absolute"}
    assert np.allclose(anchored["ErrorAligned"], 0.25)
    assert set(unanchored_task["ComparisonScale"]) == {"mean_aligned_location"}
    assert np.allclose(unanchored_task["ErrorAligned"], 0.0)


def test_generated_bundle_retains_exact_shared_data_and_distinct_anchor_inputs(tmp_path):
    manifest, identity, inventory = _write_generated_bundle(tmp_path)

    validation = pilot.validate_generated_bundle(tmp_path)
    assert validation["Passed"].all()
    assert set(inventory["File"]) == {
        "generated_ratings.csv",
        "generated_facet_truth.csv",
        "generated_anchors.csv",
        "generated_data_identity.csv",
    }
    assert set(identity["RunId"]) == set(manifest["RunId"])

    clean = identity.loc[
        identity["ConditionId"].eq("balanced_large_anchors__bias_0p0")
    ].iloc[0]
    drift = identity.loc[
        identity["ConditionId"].eq("anchor_drift__bias_0p0")
    ].iloc[0]
    assert clean["DataId"] == drift["DataId"]
    assert clean["RatingsFingerprint"] == drift["RatingsFingerprint"]
    assert clean["FacetTruthFingerprint"] == drift["FacetTruthFingerprint"]
    assert clean["AnchorFingerprint"] != drift["AnchorFingerprint"]
    assert clean["FitInputId"] != drift["FitInputId"]


def test_generated_bundle_validation_fails_after_byte_tampering(tmp_path):
    _write_generated_bundle(tmp_path)
    ratings_path = tmp_path / "generated_ratings.csv"
    ratings_path.write_bytes(ratings_path.read_bytes() + b"\n")

    with pytest.raises(ValueError, match="bundle_file_sha256"):
        pilot.validate_generated_bundle(tmp_path)
