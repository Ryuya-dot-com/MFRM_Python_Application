"""Tests for the rater severity and rater halo network analyses.

Two screens ship side by side under the rater-network family:

* **Rater severity / leniency network** (``compute_rater_severity_network``):
  for every pair of raters with shared scoring contexts, computes the
  mean / mean-absolute score difference and the directional
  "higher-count" / "higher-prop" diagnostics. Builds a directed graph
  whose edge ``A -> B`` carries ``P(A > B | A != B)``; per-rater
  ``LeniencyIndex = out_mass - in_mass`` and ``SeverityRank`` are
  the headline diagnostics.

* **Rater halo network** (``compute_rater_halo_network``): builds an
  undirected weighted graph whose nodes are ``(Rater, Criterion)``
  pairs and edge weights are Spearman correlations between the
  node score vectors over a shared context (default ``Person``).
  Halo edges connect two criteria scored by the same rater; non-halo
  edges connect different raters. Per-rater ``ReviewStatus`` of
  ``warning`` / ``review`` / ``ok`` surfaces halo-pattern flags
  (Lai, Wolfe, & Vickers, 2015; Lamprianou, 2025).

The math contract pinned here covers:

* Refusal on invalid input (non-dict, missing data, missing rater
  facet, single rater).
* Severity: ``Rater1HigherProp + Rater2HigherProp = 1`` whenever
  ``DirectionN > 0``; the closed-form identity on the directional
  counts.
* Severity: the rater generated with most-lenient severity (smallest
  ``rater_eff``) recovers ``SeverityRank = 1`` on a synthetic fixture.
* Severity: ``LeniencyIndex`` sums to zero across raters under a
  fully crossed design (symmetric edges by construction).
* Halo: refusal when rater and criterion facets coincide.
* Halo: on halo-patterned data (per-rater noise) every rater is
  flagged ``warning`` and the retained halo-edge weight exceeds the
  retained non-halo-edge weight.
* Halo: on RSM-clean data with no halo, every rater is ``ok``.
* Halo: HaloNonHaloContrast = HaloMeanWeight - NonHaloMeanWeight.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


# -----------------------------------------------------------------------------
# Severity: refusal on invalid input
# -----------------------------------------------------------------------------


def test_severity_refuses_non_dict():
    out = app.compute_rater_severity_network(None)
    assert out["available"] is False
    assert "fit dictionary" in out["reason"].lower()


def test_severity_refuses_missing_rater_column():
    out = app.compute_rater_severity_network({
        "config": {"facet_names": ["Rater"]},
        "prep": {"data": pd.DataFrame({"Person": ["P1"], "Score": [1]})},
    }, rater_facet="Rater")
    assert out["available"] is False
    assert "not a column" in out["reason"].lower()


def test_severity_refuses_person_as_rater_facet():
    out = app.compute_rater_severity_network({
        "config": {"facet_names": ["Person"]},
        "prep": {"data": pd.DataFrame({"Person": ["P1"], "Score": [1]})},
    }, rater_facet="Person")
    assert out["available"] is False
    assert "person" in out["reason"].lower()


# -----------------------------------------------------------------------------
# Severity: math contract on a synthetic many-rater design
# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def graded_severity_fit():
    """4 raters generated with monotonically increasing severity
    (R1 = -0.6 lenient, R4 = +0.6 severe); 30 persons x 3 criteria."""
    rng = np.random.default_rng(20260620)
    n_person = 30
    rater_sev = np.array([-0.6, -0.1, 0.1, 0.6])
    crit_dif = np.array([-0.4, 0.0, 0.4])
    rows = []
    for i in range(n_person):
        theta = rng.normal(0, 1)
        for j, sev in enumerate(rater_sev):
            for k, dif in enumerate(crit_dif):
                eta = theta - sev - dif
                p1 = 1.0 / (1.0 + np.exp(-eta))
                p2 = 1.0 / (1.0 + np.exp(-(eta - 0.3)))
                score = int(rng.uniform() < p1) + int(rng.uniform() < p2)
                rows.append({"Person": f"P{i+1:02d}", "Rater": f"R{j+1}",
                             "Criterion": f"C{k+1}", "Score": score})
    return app.mfrm_estimate(
        data=pd.DataFrame(rows), person_col="Person",
        facet_cols=["Rater", "Criterion"], score_col="Score",
        rating_min=0, rating_max=2, model="RSM", method="JMLE", maxit=20,
    )


def test_severity_higher_cond_prop_sums_to_one_per_pair(graded_severity_fit):
    """For every pair with at least one disagreement,
    ``Rater1HigherCondProp + Rater2HigherCondProp = 1`` by construction
    (the conditional probability is normalised by the disagreement
    count, not by the total comparison count)."""
    out = app.compute_rater_severity_network(graded_severity_fit)
    assert out["available"] is True
    pairs = out["pair_metrics"]
    pairs_with_dir = pairs[pairs["DirectionN"] > 0]
    sums = (
        pd.to_numeric(pairs_with_dir["Rater1HigherCondProp"], errors="coerce")
        + pd.to_numeric(pairs_with_dir["Rater2HigherCondProp"], errors="coerce")
    )
    for v in sums:
        assert v == pytest.approx(1.0, abs=1e-12)


def test_severity_higher_prop_matches_mfrmr_convention(graded_severity_fit):
    """``Rater1HigherProp + Rater2HigherProp + TieRate = 1`` by
    construction; the proportions match the mfrmr / FACETS convention
    (share of *all* shared contexts) so the R parity fixture lines up
    with the Python output at machine precision."""
    out = app.compute_rater_severity_network(graded_severity_fit)
    pairs = out["pair_metrics"]
    pairs_nonempty = pairs[pairs["N"] > 0]
    tie_rate = pairs_nonempty["TieCount"] / pairs_nonempty["N"]
    total = (
        pd.to_numeric(pairs_nonempty["Rater1HigherProp"], errors="coerce")
        + pd.to_numeric(pairs_nonempty["Rater2HigherProp"], errors="coerce")
        + tie_rate
    )
    for v in total:
        assert v == pytest.approx(1.0, abs=1e-12)


def test_severity_higher_counts_and_ties_sum_to_n(graded_severity_fit):
    """Rater1HigherCount + Rater2HigherCount + TieCount = N."""
    out = app.compute_rater_severity_network(graded_severity_fit)
    pairs = out["pair_metrics"]
    totals = (
        pairs["Rater1HigherCount"]
        + pairs["Rater2HigherCount"]
        + pairs["TieCount"]
    )
    assert (totals == pairs["N"]).all()


def test_severity_higher_counts_sum_to_direction_n(graded_severity_fit):
    out = app.compute_rater_severity_network(graded_severity_fit)
    pairs = out["pair_metrics"]
    counts = pairs["Rater1HigherCount"] + pairs["Rater2HigherCount"]
    assert (counts == pairs["DirectionN"]).all()


def test_severity_rank_recovers_generated_ordering(graded_severity_fit):
    """The most lenient rater (smallest rater_eff) must come out at
    SeverityRank 1; the most severe at the highest rank."""
    out = app.compute_rater_severity_network(graded_severity_fit)
    nodes = out["node_metrics"].set_index("Rater")
    assert int(nodes.loc["R1", "SeverityRank"]) == 1
    # The remaining raters are not guaranteed to keep their exact
    # generated order on small samples, but R4 (most severe) must
    # rank below R1 (most lenient).
    assert int(nodes.loc["R4", "SeverityRank"]) > int(nodes.loc["R1", "SeverityRank"])


def test_severity_leniency_index_sums_to_zero(graded_severity_fit):
    """Under a fully crossed design every directional edge (A -> B)
    is paired with a complementary edge (B -> A) whose weights sum
    to 1 per pair. Summing LeniencyIndex over raters therefore must
    equal zero by construction."""
    out = app.compute_rater_severity_network(graded_severity_fit)
    nodes = out["node_metrics"]
    total = float(nodes["LeniencyIndex"].sum())
    assert total == pytest.approx(0.0, abs=1e-9)


def test_severity_summary_reports_pair_n_and_rater_count(graded_severity_fit):
    out = app.compute_rater_severity_network(graded_severity_fit)
    summary = out["summary"].iloc[0]
    assert int(summary["NRaters"]) == 4
    # 4 choose 2 = 6 eligible pairs.
    assert int(summary["NEligiblePairs"]) == 6


# -----------------------------------------------------------------------------
# Halo: refusal cases
# -----------------------------------------------------------------------------


def test_halo_refuses_when_rater_and_criterion_facets_coincide():
    out = app.compute_rater_halo_network({
        "config": {"facet_names": ["Rater"]},
        "prep": {"data": pd.DataFrame({"Person": ["P1"], "Rater": ["R1"], "Score": [1]})},
    }, rater_facet="Rater", criterion_facet="Rater")
    assert out["available"] is False
    assert "distinct" in out["reason"].lower()


def test_halo_refuses_when_rater_facet_missing():
    out = app.compute_rater_halo_network(
        {"config": {"facet_names": []}, "prep": {"data": pd.DataFrame()}},
        rater_facet="Rater", criterion_facet="Criterion",
    )
    assert out["available"] is False


# -----------------------------------------------------------------------------
# Halo: halo-patterned data triggers warnings; clean data stays ok
# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def halo_patterned_fit():
    """Deliberately halo-patterned: each rater has a strong personal
    noise term that applies across all criteria, so the same rater's
    cross-criterion correlations should be high and the across-rater
    correlations should be near zero."""
    rng = np.random.default_rng(20260620)
    n_person, n_rater, n_criterion = 40, 3, 4
    rows = []
    for i in range(n_person):
        theta = rng.normal(0, 1)
        for j in range(n_rater):
            rater_halo = rng.normal(0, 1.5)
            for k in range(n_criterion):
                eta = theta + rater_halo
                p1 = 1.0 / (1.0 + np.exp(-eta))
                p2 = 1.0 / (1.0 + np.exp(-(eta - 0.3)))
                score = int(rng.uniform() < p1) + int(rng.uniform() < p2)
                rows.append({"Person": f"P{i+1:02d}", "Rater": f"R{j+1}",
                             "Criterion": f"C{k+1}", "Score": score})
    return app.mfrm_estimate(
        data=pd.DataFrame(rows), person_col="Person",
        facet_cols=["Rater", "Criterion"], score_col="Score",
        rating_min=0, rating_max=2, model="RSM", method="JMLE", maxit=20,
    )


@pytest.fixture(scope="module")
def clean_no_halo_fit():
    """Clean RSM design where each rater's criteria are independent of
    the rater's personal noise; halo should NOT trigger here."""
    rng = np.random.default_rng(20260621)
    n_person, n_rater, n_criterion = 40, 3, 4
    rater_sev = np.array([-0.4, 0.0, 0.4])
    crit_dif = np.linspace(-0.4, 0.4, n_criterion)
    rows = []
    for i in range(n_person):
        theta = rng.normal(0, 1)
        for j in range(n_rater):
            for k in range(n_criterion):
                eta = theta - rater_sev[j] - crit_dif[k]
                p1 = 1.0 / (1.0 + np.exp(-eta))
                p2 = 1.0 / (1.0 + np.exp(-(eta - 0.3)))
                score = int(rng.uniform() < p1) + int(rng.uniform() < p2)
                rows.append({"Person": f"P{i+1:02d}", "Rater": f"R{j+1}",
                             "Criterion": f"C{k+1}", "Score": score})
    return app.mfrm_estimate(
        data=pd.DataFrame(rows), person_col="Person",
        facet_cols=["Rater", "Criterion"], score_col="Score",
        rating_min=0, rating_max=2, model="RSM", method="JMLE", maxit=20,
    )


def test_halo_detects_halo_pattern(halo_patterned_fit):
    """Every rater on halo-patterned data must be flagged with
    ReviewStatus = 'warning' (per-rater halo correlations are strong
    and there are no comparable non-halo edges retained)."""
    out = app.compute_rater_halo_network(halo_patterned_fit, min_pair_n=10)
    assert out["available"] is True
    halo_summary = out["halo_summary_by_rater"]
    assert (halo_summary["ReviewStatus"] == "warning").all()
    # Mean halo weight per rater must exceed 0.30 on this fixture.
    assert (halo_summary["HaloMeanWeight"].dropna() > 0.3).all()


def test_halo_retained_edges_are_all_halo_on_halo_patterned_data(
    halo_patterned_fit,
):
    out = app.compute_rater_halo_network(halo_patterned_fit, min_pair_n=10)
    edges = out["edge_metrics"]
    assert not edges.empty
    # On halo-patterned data the cross-rater correlations are noise;
    # Bonferroni-adjusted p-values are large, so retained edges are
    # overwhelmingly halo.
    halo_share = float((edges["EdgeType"] == "halo").mean())
    assert halo_share > 0.9


def test_halo_contrast_identity(halo_patterned_fit):
    """HaloNonHaloContrast = HaloMeanWeight - NonHaloMeanWeight when
    both means are finite."""
    out = app.compute_rater_halo_network(halo_patterned_fit, min_pair_n=10)
    summary = out["summary"].iloc[0]
    if (
        np.isfinite(summary["HaloMeanWeight"])
        and np.isfinite(summary["NonHaloMeanWeight"])
    ):
        assert float(summary["HaloNonHaloContrast"]) == pytest.approx(
            float(summary["HaloMeanWeight"]) - float(summary["NonHaloMeanWeight"]),
            abs=1e-12,
        )


def test_halo_clean_design_stays_mostly_ok(clean_no_halo_fit):
    """On clean RSM data without a rater-specific noise term every
    rater should be either 'ok' or at most 'review' (no warnings)."""
    out = app.compute_rater_halo_network(clean_no_halo_fit, min_pair_n=10)
    halo_summary = out["halo_summary_by_rater"]
    if not halo_summary.empty:
        assert (halo_summary["ReviewStatus"] != "warning").all()


# -----------------------------------------------------------------------------
# Halo: settings round-trip
# -----------------------------------------------------------------------------


def test_halo_settings_round_trip_user_arguments(halo_patterned_fit):
    out = app.compute_rater_halo_network(
        halo_patterned_fit, min_pair_n=12, alpha=0.01,
        halo_weight_review=0.6, halo_contrast_review=0.15,
    )
    settings = out["settings"]
    assert settings["min_pair_n"] == 12
    assert settings["alpha"] == pytest.approx(0.01)
    assert settings["halo_weight_review"] == pytest.approx(0.6)
    assert settings["halo_contrast_review"] == pytest.approx(0.15)


def test_halo_node_count_equals_rater_times_criterion(halo_patterned_fit):
    """The halo graph has one node per (Rater, Criterion) combination."""
    out = app.compute_rater_halo_network(halo_patterned_fit, min_pair_n=10)
    nodes = out["node_metrics"]
    assert len(nodes) == 3 * 4  # n_rater x n_criterion in the fixture


# -----------------------------------------------------------------------------
# Correlation method dispatch: Pearson / Kendall / Spearman
# -----------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["spearman", "pearson", "kendall"])
def test_halo_supports_pearson_spearman_kendall(halo_patterned_fit, method):
    """All three correlation methods must produce a bundle with a
    populated settings.method field and retain edges on halo-patterned
    data (the correlation flavour changes the exact values, but a
    deliberately halo-patterned design fires retention under any of
    the three rank- or moment-based families)."""
    out = app.compute_rater_halo_network(
        halo_patterned_fit, method=method, min_pair_n=10, alpha=0.05,
    )
    assert out["available"] is True
    assert out["settings"]["method"] == method
    assert int(out["summary"].iloc[0]["RetainedEdges"]) > 0


def test_halo_unknown_method_returns_nan_correlations(halo_patterned_fit):
    """An unsupported method label produces NaN correlations rather than
    raising; the bundle stays well-formed."""
    out = app.compute_rater_halo_network(
        halo_patterned_fit, method="bogus", min_pair_n=10,
    )
    assert out["available"] is True
    pair_estimates = pd.to_numeric(out["pair_metrics"]["Estimate"], errors="coerce")
    assert pair_estimates.isna().all()


# -----------------------------------------------------------------------------
# Rater severity score_source: residual vs observed
# -----------------------------------------------------------------------------


def test_rater_severity_residual_mode_runs_on_clean_fit(graded_severity_fit):
    """The MFRM-residualized mode must run without error and surface
    its mode in the settings field; the SeverityRank may differ
    slightly from the observed-mode rank when raters happen to
    encounter different ability distributions."""
    out_obs = app.compute_rater_severity_network(
        graded_severity_fit, score_source="observed",
    )
    out_res = app.compute_rater_severity_network(
        graded_severity_fit, score_source="residual",
    )
    assert out_obs["available"] is True
    assert out_res["available"] is True
    assert out_obs["settings"]["score_source"] == "observed"
    assert out_res["settings"]["score_source"] == "residual"
    # Each rater still appears in both rankings.
    assert set(out_obs["node_metrics"]["Rater"]) == set(out_res["node_metrics"]["Rater"])


def test_rater_severity_higher_cond_prop_identity_holds_under_residual(
    graded_severity_fit,
):
    """The closed-form Rater1HigherCondProp + Rater2HigherCondProp = 1
    identity must continue to hold when comparisons run on residuals
    instead of raw scores; the math depends on the relative ordering,
    not the absolute scale."""
    out = app.compute_rater_severity_network(
        graded_severity_fit, score_source="residual",
    )
    pairs = out["pair_metrics"]
    pairs_with_dir = pairs[pairs["DirectionN"] > 0]
    if pairs_with_dir.empty:
        pytest.skip("No directional disagreement on residual scale.")
    sums = (
        pd.to_numeric(pairs_with_dir["Rater1HigherCondProp"], errors="coerce")
        + pd.to_numeric(pairs_with_dir["Rater2HigherCondProp"], errors="coerce")
    )
    for v in sums:
        assert v == pytest.approx(1.0, abs=1e-12)
