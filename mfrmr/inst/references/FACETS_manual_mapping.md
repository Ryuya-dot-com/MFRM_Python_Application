# FACETS Manual Mapping

This package mirrors FACETS-style workflows while keeping a native R implementation.

## Manual References Used
- Output table index: <https://www.winsteps.com/facetman/outputtableindex.htm>
- Table 5 measurable data summary: <https://www.winsteps.com/facetman/table5.htm>
- Table 7 measurement report: <https://www.winsteps.com/facetman/table7.htm>
- Table 7 reliability/chi-square notes: <https://www.winsteps.com/facetman/table7summarystatistics.htm>
- Table 7 agreement notes: <https://www.winsteps.com/facetman/table7agreementstatistics.htm>
- Table 8.1 rating scale report: <https://www.winsteps.com/facetman/table8_1ratingscale.htm>
- Table 8 bar-chart report: <https://www.winsteps.com/facetman/table8barchart.htm>
- Table 8 probability-curves report: <https://www.winsteps.com/facetman/table8curves.htm>
- Graph output (`Graphfile=`): <https://www.winsteps.com/facetman/graphoutputfile.htm>
- Table 9 bias-iteration report: <https://www.winsteps.com/facetman/table9.htm>
- Table 10 unexpected-after-bias report: <https://www.winsteps.com/facetman/table10.htm>
- Table 11 bias-calculation counts report: <https://www.winsteps.com/facetman/table11.htm>
- Table 12 DIF/bias summary report: <https://www.winsteps.com/facetman/table12.htm>
- Table 13 DIF/bias detail report: <https://www.winsteps.com/facetman/table13.htm>
- Table 14 pairwise bias report: <https://www.winsteps.com/facetman/table14.htm>

## Implemented (Direct or Close FACETS Equivalent)
- Core multifacet estimation (RSM/PCM, MML/JML): `fit_mfrm()` / `mfrm_estimate()`
- Diagnostics core bundle (obs, fit, reliability, interactions, subsets): `diagnose_mfrm()`
- Table 1-style specification summary: `specifications_report()`
- Table 2-style data summary report: `data_quality_report()`
- Table 3-style iteration report (replayed): `estimation_iteration_report()`
- Table 4-style unexpected responses: `unexpected_response_table()`, `plot_unexpected()`
- Table 5-style measurable summary bundle: `measurable_summary_table()`
- Table 6.0.0-style subset/disjoint listing: `subset_connectivity_report()`
- Table 6.2-style facet-statistics graphic summary: `facet_statistics_report()`
- Table 7-style facet/person measures and fit summary: `diagnose_mfrm()` + `summary.mfrm_fit()`
- Table 7 reliability + facet chi-square style summaries: `diagnose_mfrm()$reliability`, `diagnose_mfrm()$facets_chisq`, `facets_chisq_table()`, `plot_facets_chisq()`
- Table 7 agreement style summaries: `diagnose_mfrm()$interrater`, `interrater_agreement_table()`, `plot_interrater_agreement()`
- Table 8.1-style rating scale bundle: `rating_scale_table()`
- Table 8-style bar-chart and curves exporters: `category_structure_report()`, `category_curves_report()`
- Output-file emulation (`GRAPH=` / `SCORE=` style): `facets_output_file_bundle()`
- Table 12 fair-average style output bundle: `fair_average_table()`, `plot_fair_average()`
- Displacement diagnostics (FACETS-style anchor drift check): `displacement_table()`, `plot_displacement()`
- Bias re-estimation iteration (Table 9 workflow): `estimate_bias()` / `estimate_bias_interaction()`
- Table 10-style unexpected-after-bias output: `unexpected_after_bias_table()`
- Table 11-style bias-count report: `bias_count_table()`
- Table 12/13/14 style bias outputs:
  - summary/detail: `estimate_bias()` result (`summary`, `table`, `chi_sq`)
  - higher-order mode (3+ facets): `estimate_bias(..., interaction_facets = c(...))`
  - pairwise: `build_fixed_reports()` (`pairwise_table`, `pairwise_fixed`) for 2-way runs
  - Table 13 plot export: `bias_interaction_report()`, `plot_bias_interaction()` (including `facet_profile` mode)
- Fixed-width report generation: `build_fixed_reports()`
- APA narrative/table helpers: `build_apa_outputs()`, `apa_table()`
- Residual PCA checks (overall + by facet): `analyze_residual_pca()`, `plot_residual_pca()`
- QC dashboard (base graphics): `plot_qc_dashboard()`
- Anchor workflow: `audit_mfrm_anchors()`, `make_anchor_table()`
- Data packaging/loading helpers: `list_mfrmr_data()`, `load_mfrmr_data()`
- Automated FACETS parity checks (columns + core metrics): `facets_parity_report()`, `tests/testthat/test-facets-column-contract.R`, `tests/testthat/test-facets-metric-contract.R`, `inst/references/facets_column_contract.csv`

## Partial (Implemented Concept, Not Yet 1:1 FACETS Output)
- Design policy:
  - structured tables and visualization APIs are primary deliverables
  - fixed-width / line-printer text is optional and secondary (audit/log use)
  - exact FACETS line-printer emulation is intentionally out of scope
  - legacy numbered `table*` names are internal and not exported
- Table 1/2/3 reports:
  - current: `specifications_report()`, `data_quality_report()`, `estimation_iteration_report()` with structured output and optional fixed-width text
  - gap: FACETS fixed-width text layout and exact optimizer-internal iteration path are not yet 1:1
- Table 5 measurable data summary:
  - current: `measurable_summary_table()` and `describe_mfrm_data()` (including observed inter-rater agreement bundle)
  - gap: FACETS column-by-column textual layout parity is not exact
- Table 8.1 rating-scale report:
  - current: `rating_scale_table()` plus CCC/pathway visualization (`plot.mfrm_fit`, QC category panel)
  - gap: FACETS text layout and all legacy columns/order are not yet 1:1
- Table 8 bar-chart / probability-curves exporters:
  - current: `category_structure_report()` and `category_curves_report()` including Graphfile-style wide output and optional fixed-width text mirrors
  - gap: exact FACETS line-printer artwork/fixed-column parity is intentionally not targeted
- Table 6.2 graphical facet-statistics report:
  - current: `facet_statistics_report()` with fixed-width `M/S/Q/X` rulers
  - gap: FACETS native table layout and printer-graph formatting are not yet 1:1
- Output-file emulation (`GRAPH=` / `SCORE=`):
  - current: `facets_output_file_bundle()` with graph coordinates, observation-level modeled score export, optional fixed-width mirrors, and optional file writing
  - gap: FACETS command-level options and fixed-column file-writing parity are not yet 1:1
- Table 14 pairwise contrast report:
  - current: available for 2-way bias runs via `build_fixed_reports()`
  - gap: FACETS native layout and options are broader; higher-order runs intentionally omit pairwise section

## Not Yet Implemented from Output Index Scope
- No high-priority items in the current Table 5-14 scope; remaining gaps are mostly formatting/options parity.

## Anchoring Rules Encoded
- Direct anchors (`Facet`, `Level`, `Anchor`) are fixed.
- Group anchors (`Facet`, `Level`, `Group`, `GroupValue`) constrain group means.
- If a level appears in both tables, direct anchor takes precedence.
- Missing `GroupValue` is treated as 0.
- Default recommendation thresholds:
  - common anchors per linking facet: `>= 5`
  - observations per element: `>= 30`
  - observations per score category: `>= 10`

## Pre-release Status (Current)
- Core estimation and diagnostics: ready
- New inter-rater and facet-chi-square APIs: implemented and tested
- Remaining work for closer FACETS parity is mostly report-format completeness, not model-estimation correctness
