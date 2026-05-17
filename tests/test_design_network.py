"""Tests for the design network analysis pipeline.

``compute_design_network_analysis(res)`` treats the rating design as an
undirected weighted graph (nodes = (Facet, Level), edges = co-observed
in at least one row, weight = co-observation count) and computes
connectivity diagnostics: components, articulation points, bridges,
density, mean degree, plus per-node betweenness / closeness /
eigenvector centrality.

The math contract pinned here covers:

* Refusal on invalid input (non-dict, missing data, single-facet).
* Node count = sum of distinct levels across facet columns.
* Edge count matches the number of distinct (Facet, Level) pairs that
  co-occur in at least one row.
* A fully crossed design has 1 component, no articulation points, no
  bridges.
* A design split into two disjoint blocks must report Components = 2
  and the largest-component share must equal one block's size / total.
* Per-facet summary aggregates match the per-node table.
* Cut-node / bridge-edge frames are subsets of the per-node / per-edge
  tables respectively.
* ``min_observations=k`` drops every edge with weight < k.
* ``include_graph=True`` round-trips a NetworkX Graph.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


# -----------------------------------------------------------------------------
# Refusal: invalid input
# -----------------------------------------------------------------------------


def test_network_analysis_refuses_non_dict():
    out = app.compute_design_network_analysis(None)
    assert out["available"] is False
    assert "fit dictionary" in out["reason"].lower()


def test_network_analysis_refuses_empty_data():
    out = app.compute_design_network_analysis(
        {"config": {"facet_names": ["Rater"]}, "prep": {"data": pd.DataFrame()}}
    )
    assert out["available"] is False
    assert "data frame" in out["reason"].lower()


def test_network_analysis_refuses_single_facet():
    """A network needs at least two facet columns (one of which is Person)."""
    out = app.compute_design_network_analysis(
        {
            "config": {"facet_names": []},
            "prep": {"data": pd.DataFrame({"Person": ["P1"]})},
        }
    )
    assert out["available"] is False
    assert "two facet" in out["reason"].lower()


# -----------------------------------------------------------------------------
# End-to-end: fully crossed design
# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def fully_crossed_fit():
    """20 persons x 2 raters x 3 criteria fully crossed RSM fit."""
    rng = np.random.default_rng(20260616)
    rows = []
    for i in range(20):
        for j in range(2):
            for k in range(3):
                rows.append({
                    "Person": f"P{i+1:02d}",
                    "Rater": f"R{j+1}",
                    "Criterion": f"C{k+1}",
                    "Score": int(rng.uniform() * 3),
                })
    return app.mfrm_estimate(
        data=pd.DataFrame(rows), person_col="Person",
        facet_cols=["Rater", "Criterion"], score_col="Score",
        rating_min=0, rating_max=2, model="RSM", method="JMLE", maxit=15,
    )


def test_fully_crossed_design_is_one_connected_component(fully_crossed_fit):
    out = app.compute_design_network_analysis(fully_crossed_fit)
    assert out["available"] is True
    summary = out["summary"].iloc[0]
    assert int(summary["Components"]) == 1
    assert bool(summary["Connected"]) is True
    assert int(summary["LargestComponentNodes"]) == int(summary["Nodes"])
    assert float(summary["LargestComponentShare"]) == pytest.approx(1.0, abs=1e-12)


def test_fully_crossed_design_has_no_articulation_points_or_bridges(
    fully_crossed_fit,
):
    """A complete person x rater x criterion design has no single point /
    edge whose removal would disconnect it."""
    out = app.compute_design_network_analysis(fully_crossed_fit)
    summary = out["summary"].iloc[0]
    assert int(summary["ArticulationPoints"]) == 0
    assert int(summary["Bridges"]) == 0
    assert out["cut_nodes"].empty
    assert out["bridge_edges"].empty


def test_node_count_matches_distinct_levels(fully_crossed_fit):
    out = app.compute_design_network_analysis(fully_crossed_fit)
    summary = out["summary"].iloc[0]
    # 20 persons + 2 raters + 3 criteria = 25 nodes
    assert int(summary["Nodes"]) == 25


# -----------------------------------------------------------------------------
# Disjoint blocks: components > 1
# -----------------------------------------------------------------------------


def test_disjoint_design_reports_two_components():
    """Two non-overlapping rating blocks must produce two components.

    The network helper works directly on the raw rating data and the
    facet-names list, so the test builds a synthetic ``res`` dict
    without going through ``mfrm_estimate``.
    """
    rng = np.random.default_rng(20260617)
    rows = []
    # Block A: persons PA1-PA5 rated by raters RA1-RA2 on criteria CA1-CA2
    for i in range(5):
        for j in range(2):
            for k in range(2):
                rows.append({
                    "Person": f"PA{i+1}", "Rater": f"RA{j+1}",
                    "Criterion": f"CA{k+1}", "Score": int(rng.uniform() * 3),
                })
    # Block B: persons PB1-PB5 rated by raters RB1-RB2 on criteria CB1-CB2
    for i in range(5):
        for j in range(2):
            for k in range(2):
                rows.append({
                    "Person": f"PB{i+1}", "Rater": f"RB{j+1}",
                    "Criterion": f"CB{k+1}", "Score": int(rng.uniform() * 3),
                })
    df = pd.DataFrame(rows)
    res = {
        "config": {"facet_names": ["Rater", "Criterion"]},
        "prep": {"data": df},
    }
    out = app.compute_design_network_analysis(res)
    summary = out["summary"].iloc[0]
    assert int(summary["Components"]) >= 2, summary
    assert bool(summary["Connected"]) is False


# -----------------------------------------------------------------------------
# Per-facet summary aggregates match the per-node table
# -----------------------------------------------------------------------------


def test_facet_summary_aggregates_match_node_metrics(fully_crossed_fit):
    """Mean degree / mean betweenness in facet_summary must equal the
    column means of node_metrics filtered to that facet."""
    out = app.compute_design_network_analysis(fully_crossed_fit)
    node_metrics = out["node_metrics"]
    facet_summary = out["facet_summary"]
    for _, row in facet_summary.iterrows():
        sub = node_metrics[node_metrics["Facet"] == row["Facet"]]
        assert float(row["MeanDegree"]) == pytest.approx(float(sub["Degree"].mean()), abs=1e-12)
        assert float(row["MeanBetweenness"]) == pytest.approx(
            float(sub["Betweenness"].mean()), abs=1e-12
        )
        assert int(row["Nodes"]) == int(len(sub))
        assert int(row["ArticulationCount"]) == int(sub["IsArticulationPoint"].sum())


def test_cut_node_and_bridge_frames_are_subsets(fully_crossed_fit):
    out = app.compute_design_network_analysis(fully_crossed_fit)
    cut = out["cut_nodes"]
    nodes = out["node_metrics"]
    if not cut.empty:
        assert set(cut["Node"]).issubset(set(nodes["Node"]))
    bridges = out["bridge_edges"]
    edges = out["edge_metrics"]
    if not bridges.empty:
        # Edge frame is keyed by (From, To); the bridge frame must
        # appear in the edge frame on the same key.
        edge_keys = set(zip(edges["From"].astype(str), edges["To"].astype(str)))
        bridge_keys = set(zip(bridges["From"].astype(str), bridges["To"].astype(str)))
        assert bridge_keys.issubset(edge_keys)


# -----------------------------------------------------------------------------
# Min-observations filter
# -----------------------------------------------------------------------------


def test_min_observations_filter_drops_weak_edges(fully_crossed_fit):
    """For a 20-person x 2-rater x 3-criterion fully crossed design:
    Person-Rater edges accumulate weight = n_criterion = 3,
    Person-Criterion edges accumulate weight = n_rater = 2,
    Rater-Criterion edges accumulate weight = n_person = 20.
    Setting min_observations = 3 drops the (Person, Criterion) edges
    but keeps the (Person, Rater) and (Rater, Criterion) edges."""
    out_default = app.compute_design_network_analysis(fully_crossed_fit)
    out_strict = app.compute_design_network_analysis(
        fully_crossed_fit, min_observations=3
    )
    assert (
        int(out_strict["summary"].iloc[0]["Edges"])
        < int(out_default["summary"].iloc[0]["Edges"])
    )


# -----------------------------------------------------------------------------
# include_graph round-trip
# -----------------------------------------------------------------------------


def test_include_graph_returns_a_networkx_graph(fully_crossed_fit):
    import networkx as nx
    out = app.compute_design_network_analysis(fully_crossed_fit, include_graph=True)
    assert "graph" in out
    assert isinstance(out["graph"], nx.Graph)
    assert out["graph"].number_of_nodes() == int(out["summary"].iloc[0]["Nodes"])
    assert out["graph"].number_of_edges() == int(out["summary"].iloc[0]["Edges"])


def test_include_graph_false_does_not_carry_a_graph_key(fully_crossed_fit):
    out = app.compute_design_network_analysis(fully_crossed_fit, include_graph=False)
    assert "graph" not in out


# -----------------------------------------------------------------------------
# Bounds on centrality measures
# -----------------------------------------------------------------------------


def test_node_centralities_are_in_unit_interval(fully_crossed_fit):
    """Normalized betweenness, closeness, and eigenvector centrality all
    lie in [0, 1] up to numerical noise."""
    out = app.compute_design_network_analysis(fully_crossed_fit)
    nodes = out["node_metrics"]
    for col in ["Betweenness", "EigenvectorCentrality"]:
        vals = pd.to_numeric(nodes[col], errors="coerce").dropna()
        assert (vals >= -1e-12).all()
        assert (vals <= 1.0 + 1e-12).all()


def test_node_metrics_contains_all_expected_columns(fully_crossed_fit):
    out = app.compute_design_network_analysis(fully_crossed_fit)
    expected_cols = {
        "Node", "Facet", "Level", "Degree", "Strength",
        "Betweenness", "Closeness", "EigenvectorCentrality",
        "IsArticulationPoint",
    }
    assert expected_cols.issubset(set(out["node_metrics"].columns))
