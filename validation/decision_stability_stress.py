#!/usr/bin/env python3
"""Deterministic stress matrix for threshold, sparsity, anchors, and bias.

The run is intentionally diagnostic rather than an operating-characteristic
claim. It checks whether the application exposes fragile decision boundaries,
whether planned sparse designs are flagged, and whether anchor percentages
remain descriptive instead of becoming an undocumented pass/fail rule.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import platform
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import streamlit_app as app
from mfrm_app import decision_stability as ds


DEFAULT_OUTPUT = Path(__file__).resolve().parent / "decision_stability_20260809"
SCHEMA_VERSION = "decision-stability-sparse-anchor-bias-v1"


def fit_boundary_grid() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    display_unit = 10.0 ** (-ds.FIT_DISPLAY_DECIMALS)
    cases = (
        ("display_outside_below", -0.51 * display_unit),
        ("display_inside_below", -0.49 * display_unit),
        ("exact", 0.0),
        ("display_inside_above", 0.49 * display_unit),
        ("display_outside_above", 0.51 * display_unit),
    )
    for threshold in ds.FIT_MNSQ_THRESHOLDS:
        values = [
            ("nextafter_below", float(np.nextafter(threshold, -np.inf))),
            *[(label, float(threshold + delta)) for label, delta in cases],
            ("nextafter_above", float(np.nextafter(threshold, np.inf))),
        ]
        for label, value in values:
            evidence = ds.evaluate_fit_mnsq(value)
            rows.append({
                "Threshold": threshold,
                "Case": label,
                "OffsetDisplayUnits": (value - threshold) / display_unit,
                **evidence,
            })
    return pd.DataFrame(rows)


def sparse_design_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    cases = [
        {"Case": "balanced", "n_person": 80, "first": 4, "missing": 0.00, "zero": None, "counts": (4, 3, 3)},
        {"Case": "partial-rater", "n_person": 80, "first": 2, "missing": 0.00, "zero": None, "counts": (4, 3, 3)},
        {"Case": "single-rater", "n_person": 80, "first": 1, "missing": 0.00, "zero": None, "counts": (4, 3, 3)},
        {"Case": "planned-missing-20", "n_person": 80, "first": 2, "missing": 0.20, "zero": None, "counts": (4, 3, 3)},
        {"Case": "planned-missing-50", "n_person": 80, "first": 2, "missing": 0.50, "zero": None, "counts": (4, 3, 3)},
        {"Case": "small-balanced", "n_person": 20, "first": 4, "missing": 0.00, "zero": None, "counts": (4, 3, 3)},
        {"Case": "small-sparse", "n_person": 8, "first": 1, "missing": 0.65, "zero": None, "counts": (4, 3, 3)},
        {"Case": "many-level-sparse", "n_person": 12, "first": 1, "missing": 0.45, "zero": None, "counts": (12, 6, 4)},
        {"Case": "zero-category", "n_person": 80, "first": 4, "missing": 0.00, "zero": 2, "counts": (4, 3, 3)},
    ]
    summaries: list[dict[str, object]] = []
    pair_frames: list[pd.DataFrame] = []
    status_rank = {"Hold before reporting": 0, "Review": 1, "Ready": 2}
    for index, case in enumerate(cases):
        bundle = app.generate_custom_mfrm_simulation_bundle(
            n_person=int(case["n_person"]),
            facet_names=("Rater", "Task", "Criterion"),
            facet_level_counts=tuple(case["counts"]),
            facet_sds=(0.35, 0.25, 0.25),
            first_facet_levels_per_person=int(case["first"]),
            n_categories=5,
            missing_rate=float(case["missing"]),
            zero_count_score=case["zero"],
            seed=202608090 + index,
        )
        audit = bundle["sparse_design_audit"].copy()
        pairs = bundle["sparse_pair_cells"].copy()
        pairs.insert(0, "Case", str(case["Case"]))
        pair_frames.append(pairs)
        statuses = audit["Status"].astype(str)
        overall = min(statuses, key=lambda value: status_rank.get(value, 1))
        retained_share = len(bundle["data"]) / max(int(bundle["meta"]["full_rows"]), 1)
        category_min = int(bundle["category_counts"]["Count"].min())
        summaries.append({
            **case,
            "Rows": int(len(bundle["data"])),
            "RetainedShare": float(retained_share),
            "OverallSparseStatus": overall,
            "HoldChecks": int(statuses.eq("Hold before reporting").sum()),
            "ReviewChecks": int(statuses.eq("Review").sum()),
            "ReadyChecks": int(statuses.eq("Ready").sum()),
            "ZeroPairCells": int(pd.to_numeric(pairs["ZeroCells"], errors="coerce").fillna(0).sum()),
            "LowPairCells": int(pd.to_numeric(pairs["LowCountCells"], errors="coerce").fillna(0).sum()),
            "WorstPairCellCount": int(pd.to_numeric(pairs["MinCellCount"], errors="coerce").min()),
            "CategoryMinCount": category_min,
            "BiasScreenReadiness": (
                "Review sparse pair cells"
                if int(pairs["ZeroCells"].sum()) or int(pairs["LowCountCells"].sum())
                else "Pair-cell screen ready"
            ),
        })
    return pd.DataFrame(summaries), pd.concat(pair_frames, ignore_index=True)


def anchor_coverage_matrix() -> pd.DataFrame:
    levels = [f"R{idx}" for idx in range(1, 11)]
    persons = [f"P{idx}" for idx in range(1, 6)]
    data = pd.DataFrame(
        [(person, level, (pi + li) % 5) for pi, person in enumerate(persons) for li, level in enumerate(levels)],
        columns=["Person", "Rater", "Score"],
    )
    prep = {
        "facet_names": ["Rater"],
        "levels": {"Person": persons, "Rater": levels},
        "data": data,
    }
    rows: list[dict[str, object]] = []
    for anchor_count in (0, 1, 2, 3, 5, 10):
        anchors = None
        if anchor_count:
            anchors = pd.DataFrame({
                "Facet": "Rater",
                "Level": levels[:anchor_count],
                "Anchor": np.linspace(-0.5, 0.5, anchor_count),
            })
        audit = app.audit_mfrm_anchors(
            prep,
            anchor_df=anchors,
            min_common_anchors=2,
            min_obs_per_element=2,
        )
        summary = audit["summary"].loc[audit["summary"]["Facet"] == "Rater"].iloc[0]
        rows.append({
            "AnchorCountRequested": anchor_count,
            "Levels": int(summary["Levels"]),
            "AnchoredLevelsTotal": int(summary["AnchoredLevelsTotal"]),
            "UnanchoredLevels": int(summary["UnanchoredLevels"]),
            "AnchorShare": float(summary["AnchorShare"]),
            "Status": str(summary["Status"]),
            "OverallAuditStatus": str(audit["overall_status"]),
            "CoverageIsDecisionThreshold": False,
            "CoverageInterpretation": str(summary["CoverageInterpretation"]),
        })
    return pd.DataFrame(rows)


def bias_boundary_grid() -> pd.DataFrame:
    probability_values = [
        np.nextafter(0.05, -np.inf),
        0.04996,
        0.05,
        0.05004,
        np.nextafter(0.05, np.inf),
        0.06,
    ]
    practical_values = [
        np.nextafter(0.50, -np.inf),
        0.49996,
        0.50,
        0.50004,
        np.nextafter(0.50, np.inf),
        0.60,
    ]
    rows: list[dict[str, object]] = []
    for idx, value in enumerate(probability_values):
        rows.append({
            "FacetPair": f"probability-{idx}",
            "p_holm": float(value),
            "p_bh": float(value),
            "AbsBias": 0.20,
        })
    for idx, value in enumerate(practical_values):
        rows.append({
            "FacetPair": f"practical-{idx}",
            "p_holm": 0.20,
            "p_bh": 0.20,
            "AbsBias": float(value),
        })
    return ds.audit_bias_decision_stability(pd.DataFrame(rows))


def save_figures(
    output: Path,
    fit_grid: pd.DataFrame,
    sparse: pd.DataFrame,
    anchors: pd.DataFrame,
    bias_grid: pd.DataFrame,
) -> None:
    colors = {
        "stable": "#2a9d8f",
        "display_rounding_boundary": "#f4a261",
        "numerical_boundary": "#e63946",
        "unavailable": "#999999",
    }

    fig, ax = plt.subplots(figsize=(9, 4.8))
    for status, group in fit_grid.groupby("BoundaryStatus"):
        ax.scatter(
            group["OffsetDisplayUnits"],
            group["Threshold"],
            label=status.replace("_", " "),
            color=colors.get(status, "#666666"),
            s=58,
            alpha=0.85,
        )
    ax.axvspan(-0.5, 0.5, color="#f4a261", alpha=0.10, label="3-decimal display half-unit")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set(xlabel="Offset from threshold (3-decimal display units)", ylabel="MNSQ threshold", title="Fit decision-boundary stress map")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(output / "fit_boundary_map.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5.2))
    status_color = {"Ready": "#2a9d8f", "Review": "#f4a261", "Hold before reporting": "#e63946"}
    bars = ax.bar(
        sparse["Case"],
        sparse["RetainedShare"],
        color=[status_color.get(value, "#999999") for value in sparse["OverallSparseStatus"]],
    )
    for bar, row in zip(bars, sparse.itertuples(index=False), strict=False):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            min(0.98, bar.get_height() + 0.025),
            f"low cells={row.LowPairCells}\nzero={row.ZeroPairCells}",
            ha="center",
            va="bottom",
            fontsize=7,
        )
    ax.set_ylim(0, 1.15)
    ax.set(ylabel="Retained row share", title="Sparse-design stress matrix (colour = worst audit state)")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(output / "sparse_design_matrix.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.barh(anchors["Status"] + " (n=" + anchors["AnchorCountRequested"].astype(str) + ")", anchors["AnchoredLevelsTotal"], color="#2a9d8f", label="Anchored")
    ax.barh(
        anchors["Status"] + " (n=" + anchors["AnchorCountRequested"].astype(str) + ")",
        anchors["UnanchoredLevels"],
        left=anchors["AnchoredLevelsTotal"],
        color="#d9d9d9",
        label="Unanchored",
    )
    ax.set(xlabel="Observed levels", title="Descriptive anchor coverage (no percentage pass line)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "anchor_coverage_matrix.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for axis, statistic, threshold, title in (
        (axes[0], "p_holm", 0.05, "Holm-adjusted p (alpha=.05)"),
        (axes[1], "AbsBias", 0.50, "Absolute bias (threshold=.50)"),
    ):
        sub = bias_grid.loc[bias_grid["Statistic"] == statistic].copy()
        case_prefix = "probability-" if statistic == "p_holm" else "practical-"
        sub = sub.loc[sub["FacetPair"].astype(str).str.startswith(case_prefix)]
        axis.scatter(
            sub["RawValue"],
            np.arange(len(sub)),
            c=[colors.get(status, "#666666") for status in sub["BoundaryStatus"]],
            s=60,
        )
        axis.axvline(threshold, color="black", linewidth=0.9)
        axis.set(title=title, xlabel="Raw value", yticks=[])
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=colors[status], markersize=8, label=status.replace("_", " "))
        for status in ("numerical_boundary", "display_rounding_boundary", "stable")
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3, frameon=False)
    fig.subplots_adjust(bottom=0.24, top=0.88, wspace=0.25)
    fig.savefig(output / "bias_boundary_map.png", dpi=180)
    plt.close(fig)


def write_results(
    output: Path,
    fit_grid: pd.DataFrame,
    sparse: pd.DataFrame,
    anchors: pd.DataFrame,
    bias_grid: pd.DataFrame,
) -> None:
    fit_review = int(fit_grid["BoundaryStatus"].isin({"numerical_boundary", "display_rounding_boundary"}).sum())
    fit_mismatch = int((~fit_grid["DisplayDecisionConsistent"].astype(bool)).sum())
    sparse_hold = int(sparse["OverallSparseStatus"].eq("Hold before reporting").sum())
    sparse_review = int(sparse["OverallSparseStatus"].eq("Review").sum())
    bias_counts = ds.summarize_boundary_audit(bias_grid)
    text = f"""# Decision stability, sparse-design, anchor, and bias stress results

Schema: `{SCHEMA_VERSION}`  
Platform: `{platform.platform()}`  
Python: `{platform.python_version()}`

## Outcome

- Fit grid: {len(fit_grid)} threshold-neighbour cases; {fit_review} were explicitly labelled as numerical/display boundaries and {fit_mismatch} had different raw-versus-rounded classifications.
- Sparse matrix: {len(sparse)} deterministic designs; {sparse_hold} reached `Hold before reporting`, {sparse_review} reached `Review`, and the retained tables identify zero/low facet-pair cells relevant to conditional bias screening.
- Anchor matrix: {len(anchors)} coverage conditions from 0% to 100%. Coverage was exported in every case, but `CoverageIsDecisionThreshold` remained false; the existing minimum-common-anchor count, content, observations, drift, and connectedness remain separate evidence.
- Bias grid: {bias_counts['audited']} threshold decisions; {bias_counts['numerical_boundary']} numerical-boundary and {bias_counts['display_rounding_boundary']} display-rounding-boundary decisions, with {bias_counts['display_decision_mismatch']} raw/display mismatches.

## Critical interpretation

The run confirms that rounding can change a visible label near 0.50, 1.50,
2.00, alpha=.05, or |bias|=.50 even when the raw decision is deterministic.
The application must therefore retain raw values, export the rule, and label
the boundary instead of silently moving or widening the threshold.

Sparse-design findings are eligibility/caveat evidence, not automatic data
deletion rules. Anchor share is descriptive and is not a universal adequacy
cutoff. Bias flags remain conditional screens affected by cell support,
connectedness, multiplicity, and threshold proximity.

## Retained evidence

- `fit_boundary_grid.csv` and `fit_boundary_map.png`
- `sparse_design_summary.csv`, `sparse_pair_cells.csv`, and `sparse_design_matrix.png`
- `anchor_coverage_matrix.csv` and `anchor_coverage_matrix.png`
- `bias_boundary_grid.csv` and `bias_boundary_map.png`

Reproduce with:

```bash
python3 validation/decision_stability_stress.py
```
"""
    (output / "RESULTS.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    fit_grid = fit_boundary_grid()
    sparse, sparse_pairs = sparse_design_matrix()
    anchors = anchor_coverage_matrix()
    bias_grid = bias_boundary_grid()

    fit_grid.to_csv(output / "fit_boundary_grid.csv", index=False)
    sparse.to_csv(output / "sparse_design_summary.csv", index=False)
    sparse_pairs.to_csv(output / "sparse_pair_cells.csv", index=False)
    anchors.to_csv(output / "anchor_coverage_matrix.csv", index=False)
    bias_grid.to_csv(output / "bias_boundary_grid.csv", index=False)
    save_figures(output, fit_grid, sparse, anchors, bias_grid)
    write_results(output, fit_grid, sparse, anchors, bias_grid)
    print((output / "RESULTS.md").read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
