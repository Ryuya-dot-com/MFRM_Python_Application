import json
from pathlib import Path

import numpy as np
import pandas as pd

from validation import operating_characteristics_facets as facets


def _manifest_row() -> pd.Series:
    return pd.Series({
        "RunId": "balanced_small__bias_0p0::rep-00001",
        "ConditionId": "balanced_small__bias_0p0",
        "Design": "balanced_small",
        "TruthBias": 0.0,
        "Replicate": 1,
        "Seed": 123,
        "Categories": 4,
    })


def _ratings() -> pd.DataFrame:
    return pd.DataFrame({
        "Person": ["P001", "P001", "P002", "P002"],
        "Rater": ["R01", "R02", "R01", "R02"],
        "Task": ["T01", "T01", "T01", "T01"],
        "Criterion": ["C01", "C01", "C01", "C01"],
        "Score": [0, 1, 2, 3],
    })


def test_build_facets_spec_declares_matched_rsm_bias_and_partial_anchors(tmp_path):
    anchors = pd.DataFrame({"Facet": ["Rater"], "Level": ["R01"], "Anchor": [0.25]})
    spec, mappings = facets.build_facets_spec(
        _manifest_row(),
        _ratings(),
        anchors,
        score_base=tmp_path / "scores.txt",
    )

    assert "Noncenter=1" in spec
    assert "Positive=1" in spec
    assert "Umean=0,1,6" in spec
    assert "?,?,?,?,R3" in spec
    assert "?,?B,?B,?,R3" in spec
    assert "2,Rater,A" in spec
    assert "1=R01,0.25" in spec
    assert "2=R02\n" in spec
    assert mappings["Person"] == {"P001": 1, "P002": 2}
    assert spec.rstrip().endswith("2,2,1,1,3")


def test_parse_score_file_names_fit_as_displayed_and_adds_rounding_contract(tmp_path):
    score = tmp_path / "scores.2.txt"
    score.write_text(
        "2\tRater\n"
        "Measure\tS.E.\tInfitMS\tInfitZ\tOutfitMS\tOutfitZ\tStatus\t2\tRater\n"
        "0.123456\t0.045678\t1.50\t2.00\t1.49\t1.99\t-1\t1\tR01\n",
        encoding="utf-8",
    )

    parsed = facets.parse_score_file(score, facet_number=2, facet_name="Rater")

    assert parsed.loc[0, "Estimate"] == 0.123456
    assert parsed.loc[0, "SE"] == 0.045678
    assert parsed.loc[0, "InfitMSDisplayed"] == 1.50
    assert "InfitMS" not in parsed.columns
    assert np.isclose(parsed.loc[0, "InfitMSRawLowerBound"], 1.495)
    assert np.isclose(parsed.loc[0, "InfitMSRawUpperBound"], 1.505)
    assert parsed.loc[0, "InfitMSDecision0p5To1p5"] == "boundary_uncertain"
    assert parsed.loc[0, "OutfitMSDecision0p5To1p5"] == "pass"
    assert parsed.loc[0, "InfitZDecisionAbs2"] == "boundary_uncertain"
    assert parsed.loc[0, "FitPrecisionContract"] == "FACETS_display_2dp_not_raw"


def test_displayed_band_requires_entire_rounding_interval_on_one_side():
    displayed = pd.Series([0.49, 0.50, 0.51, 1.49, 1.50, 1.51, np.nan])
    classified = facets.classify_displayed_band(displayed, lower=0.5, upper=1.5)

    assert list(classified.iloc[:6]) == [
        "flag",
        "boundary_uncertain",
        "pass",
        "pass",
        "boundary_uncertain",
        "flag",
    ]
    assert pd.isna(classified.iloc[6])


def test_parse_iteration_report_uses_final_jmle_row(tmp_path):
    report = tmp_path / "report.out.txt"
    report.write_text(
        "Facets 64-bit (Many-Facet Rasch Measurement) 4.5.0 Copyright\n"
        "| JMLE  4  -1.5098  -8.6  .0000  -.2224  .0000 |\n"
        "| JMLE 37  -.0550   -.3  .0000   .0094  .0000 |\n"
        "Subset connection O.K.\n",
        encoding="utf-8",
    )

    parsed = facets.parse_iteration_report(report)

    assert parsed["FacetsVersion"] == "4.5.0"
    assert parsed["Iterations"] == 37
    assert parsed["MaxScoreResidual"] == -0.055
    assert parsed["MaxElementLogitChange"] == 0.0094
    assert parsed["Converged"] is True
    assert parsed["SubsetConnected"] is True


def test_alignment_does_not_recenter_an_anchored_facet():
    estimates = pd.DataFrame({
        "Facet": ["Rater", "Rater", "Task", "Task"],
        "Level": ["R01", "R02", "T01", "T02"],
        "Estimate": [0.30, -0.10, 0.50, 0.10],
        "SE": [0.0, 0.1, 0.1, 0.1],
        "Status": [1, -1, -1, -1],
    })
    truth = pd.DataFrame({
        "Facet": ["Rater", "Rater", "Task", "Task"],
        "Level": ["R01", "R02", "T01", "T02"],
        "Truth": [0.20, -0.20, 0.30, -0.30],
    })
    anchors = pd.DataFrame({"Facet": ["Rater"], "Level": ["R01"], "Anchor": [0.30]})

    aligned = facets.align_recovery(estimates, truth, anchors, _manifest_row())
    rater = aligned[aligned["Facet"].eq("Rater")]
    task = aligned[aligned["Facet"].eq("Task")]

    assert rater["ComparisonScale"].eq("anchor_identified_absolute").all()
    assert np.allclose(rater["ErrorAligned"], [0.10, 0.10])
    assert task["ComparisonScale"].eq("mean_aligned_location").all()
    assert np.allclose(task["ErrorAligned"], [-0.10, 0.10])
    assert np.isclose(task["ErrorAligned"].mean(), 0.0)


def test_facets_contract_files_are_machine_readable():
    repo = Path(__file__).resolve().parents[1]
    plan = json.loads((repo / "validation" / "facets_complementary_workbench_plan_20260810.json").read_text(encoding="utf-8"))
    matrix = json.loads((repo / "validation" / "facets_compatibility_matrix.json").read_text(encoding="utf-8"))

    assert plan["product_position"] == "FACETS-validated complementary MFRM workbench"
    assert plan["rounding_contract"]["fit"].startswith("FACETS 4.5.0")
    statuses = set(matrix["status_vocabulary"])
    assert statuses == {
        "externally_matched",
        "formula_aligned",
        "analogous",
        "workbench_addition",
        "unsupported",
        "unvalidated",
        "operationally_qualified",
    }


def test_facets_agreement_crosswalk_uses_table_7_not_table_8():
    repo = Path(__file__).resolve().parents[1]
    source = (repo / "streamlit_app.py").read_text(encoding="utf-8")

    assert "FACETS Table 7-style agreement screen." in source
    assert "Table 7-style shared-context screen" in source
    assert "FACETS Table 8-style agreement screen." not in source
    assert "Table 8-style shared-context screen" not in source


def test_anchored_runs_require_explicit_absolute_origin_declaration(tmp_path):
    common = [
        "--input-dir",
        str(tmp_path / "input"),
        "--output-dir",
        str(tmp_path / "output"),
    ]

    legacy = facets.parse_args(common)
    corrected = facets.parse_args(common + ["--python-anchor-constraint", "absolute_origin"])

    assert legacy.python_anchor_constraint == "legacy_free_levels_centered"
    assert corrected.python_anchor_constraint == "absolute_origin"


def test_facets_cli_separates_immutable_input_from_python_comparison_results(tmp_path):
    args = facets.parse_args([
        "--input-dir",
        str(tmp_path / "immutable-input"),
        "--comparison-dir",
        str(tmp_path / "current-python-results"),
        "--output-dir",
        str(tmp_path / "facets-output"),
    ])

    assert args.input_dir == tmp_path / "immutable-input"
    assert args.comparison_dir == tmp_path / "current-python-results"


def test_display_parser_rejects_truncated_table8_tokens():
    assert facets.parse_display_token("1.")["status"] == "ambiguous_truncated_display"
    assert facets.parse_display_token(".")["status"] == "ambiguous_truncated_display"
    parsed = facets.parse_display_token("-.93*")
    assert parsed["value"] == -0.93
    assert parsed["decimals"] == 2
    assert parsed["nonmonotonic_marker"] is True


def test_parse_table8_uses_auxiliary_precision_and_conservative_order(tmp_path):
    report = tmp_path / "table8.out.txt"
    report.write_text(
        "Table 8.1  Category Statistics.\n\n"
        " Model = ?,?,?,?,R3\n"
        "+-----------------------------------------------------------------------------------------------------------+\n"
        "|           DATA                 |   QUALITY CONTROL |RASCH-ANDRICH|  EXPECTATION  |  MOST  |  RASCH-  | Cat|\n"
        "|      Category Counts       Cum.|  Avge  Exp. OUTFIT| Thresholds  |  Measure at   |PROBABLE| THURSTONE|PEAK|\n"
        "|Score Total      Used    %    % |  Meas  Meas  MnSq |Measure  S.E.|Category  -0.5 |  from  |Thresholds|Prob|\n"
        "|--------------------------------+-------------------+-------------+---------------+--------+----------+----|\n"
        "|  0      20        20   20%  20%| -1.41  -1.41  1.0 |             |( -2.62)       |   low  |   low    |100%|\n"
        "|  1      30        30   30%  50%|  -.45   -.44   .9 | -1.40    .07|   -.79   -1.81|  -1.40 |  -1.59   | 49%|\n"
        "|  2      30        30   30%  80%|   .52    .48  1.0 |   .01    .06|    .81     .00|    .01 |    .00   | 49%|\n"
        "|  3      20        20   20% 100%|  1.30   1.33  1.5 |  1.39    .07|(  2.63)   1.81|   1.39 |   1.58   |100%|\n"
        "+---------------------------------------------------------------------(Mean)---------(Modal)--(Median)------+\n",
        encoding="utf-8",
    )

    table = facets.parse_table8_categories(report)

    assert len(table) == 4
    assert list(table["TotalCount"]) == [20.0, 30.0, 30.0, 20.0]
    assert table.loc[0, "CategoryOutfitMSDisplayDecimals"] == 1
    assert np.isclose(table.loc[0, "CategoryOutfitMSRawLowerBound"], 0.95)
    assert table.loc[3, "CategoryOutfitDecision0p5To1p5"] == "boundary_uncertain"
    assert list(table.loc[table["Category"].ge(2), "ThresholdOrderDecision"]) == ["ordered", "ordered"]
    assert table.loc[1:, "AverageMeasureOrderDecision"].eq("ordered").all()


def test_parse_table13_selects_n_arrangement_and_applies_holm(tmp_path):
    report = tmp_path / "table13.out.txt"
    report.write_text(
        "Table 13.2.1  Bias/Interaction Report (arranged by mN).\n"
        "Bias/Interaction: 2. Rater, 3. Task (higher score = higher bias measure)\n"
        "|   10       9.00    10        .10| .90000 .2000   4.50     9 .0001 |  1.0   1.0 |  1 1 R01  .100000 1 T01  .100|\n"
        "Table 13.2.3  Bias/Interaction Report (arranged by N).\n\n"
        "Bias/Interaction: 2. Rater, 3. Task (higher score = higher bias measure)\n"
        "|Observd  Expctd  Observd  Obs-Exp| Bias+  Model                    |Infit Outfit|    Rater        Task         |\n"
        "|  Score   Score    Count  Average|  Size   S.E.     t   d.f. Prob. | MnSq  MnSq | Sq N Rat measr- N Tas measr- |\n"
        "|---------------------------------+---------------------------------+------------+------------------------------|\n"
        "|   10       9.00    10        .10| .60000 .1759   3.41     9 .0010 |  1.3   1.3 |  1 1 R01  .100000 1 T01  .100|\n"
        "|   11      10.00    10        .10| .10000 .1759    .57     9 .2000 |  1.0   1.0 |  2 2 R02 -.100000 1 T01  .100|\n"
        "|    9       9.50    10       -.05|-.10000 .1759   -.57     9 .3000 |   .9    .9 |  3 1 R01  .100000 2 T02 -.100|\n"
        "|   10      10.50    10       -.05|-.50000 .1759  -2.84     9 .4000 |  1.1   1.1 |  4 2 R02 -.100000 2 T02 -.100|\n"
        "Fixed (all = 0) chi-squared: 1.0  d.f.: 4  significance (probability): .90\n",
        encoding="utf-8",
    )

    table = facets.parse_table13_bias(
        report,
        level_maps={"Rater": {"R01": 1, "R02": 2}, "Task": {"T01": 1, "T02": 2}},
    )

    assert len(table) == 4
    focal = table[(table["Rater"].eq("R01")) & (table["Task"].eq("T01"))].iloc[0]
    assert np.isclose(focal["HolmPDisplayed"], 0.004)
    assert focal["HolmDecision"] == "significant"
    assert focal["PracticalDecision"] == "practically_large"
    assert focal["StrongDecision"] == "flag"
    boundary = table[(table["Rater"].eq("R02")) & (table["Task"].eq("T02"))].iloc[0]
    assert boundary["PracticalDecision"] == "boundary_uncertain"
    assert boundary["StrongDecision"] == "no_flag"


def test_table8_table13_plan_is_machine_readable():
    repo = Path(__file__).resolve().parents[1]
    plan = json.loads(
        (repo / "validation" / "facets_table8_table13_parser_plan_20260811.json")
        .read_text(encoding="utf-8")
    )
    assert plan["table13_contract"]["canonical_arrangement"].startswith("Table 13")
    assert plan["table13_contract"]["practical_rule"] == "Absolute bias size >= 0.50 logits."
    assert "facets_table8_categories.csv" in plan["registered_outputs"]


def test_table13_expected_cell_audit_retains_unreported_cells():
    manifest = pd.DataFrame({
        "RunId": ["run-1"],
        "ConditionId": ["condition-1"],
        "Design": ["sparse"],
        "TruthBias": [0.0],
        "TruthPositive": [False],
        "BiasCell": ["R01 x T01"],
        "Categories": [4],
        "Replicate": [1],
        "Seed": [123],
    })
    ratings = pd.DataFrame({
        "RunId": ["run-1", "run-1"],
        "Rater": ["R01", "R02"],
        "Task": ["T01", "T02"],
        "Score": [3, 0],
    })
    reported = pd.DataFrame({
        "RunId": ["run-1"],
        "FacetA": ["Rater"],
        "FacetA_Level": ["R01"],
        "FacetB": ["Task"],
        "FacetB_Level": ["T01"],
        "ObservedCountDisplayed": [1.0],
        "ObservedScoreDisplayed": [3.0],
        "BiasSizeDisplayed": [0.1],
        "PDisplayed": [0.5],
        "StrongDecision": ["sparse_no_claim"],
    })

    audit = facets.table13_expected_cell_audit(manifest, ratings, reported)

    assert len(audit) == 4
    assert audit["FACETSReported"].sum() == 1
    assert (~audit["FACETSReported"]).sum() == 3
    assert audit.loc[~audit["FACETSReported"], "NegativeFindingAllowed"].eq(False).all()
    assert audit.loc[audit["FACETSReported"], "FACETSCountNotAboveGenerated"].iloc[0]


def test_bias_agreement_fails_closed_on_run_ineligibility(tmp_path):
    stability = pd.DataFrame({
        "RunId": ["run-1", "run-1"],
        "FacetA": ["Rater", "Rater"],
        "FacetA_Level": ["R01", "R01"],
        "FacetB": ["Task", "Task"],
        "FacetB_Level": ["T01", "T01"],
        "Statistic": ["p_holm", "AbsBias"],
        "RawValue": [0.001, 0.8],
    })
    stability.to_csv(tmp_path / "bias_decision_stability.csv", index=False)
    pd.DataFrame({
        "RunId": ["run-1"],
        "FitReturned": [True],
        "Converged": [False],
        "AnalysisEligible": [False],
    }).to_csv(tmp_path / "runs.csv", index=False)
    facets_bias = pd.DataFrame({
        "RunId": ["run-1"],
        "Design": ["balanced"],
        "TruthBias": [0.6],
        "BiasCell": ["R01 x T01"],
        "FacetA": ["Rater"],
        "FacetA_Level": ["R01"],
        "FacetB": ["Task"],
        "FacetB_Level": ["T01"],
        "BiasSizeDisplayed": [0.8],
        "StrongDecision": ["flag"],
        "SparseCell": [False],
        "DirectComparisonEligible": [True],
    })

    pairs, summary = facets.facets_python_bias_agreement(facets_bias, tmp_path)

    assert not pairs.loc[0, "ComparisonEligible"]
    assert summary.loc[0, "ComparisonEligibleCells"] == 0
