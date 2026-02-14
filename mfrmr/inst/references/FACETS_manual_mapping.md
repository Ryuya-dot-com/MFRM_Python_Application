# FACETS Manual Mapping

This package is being developed to mirror FACETS-style workflows while keeping a native R implementation.

Reference pages (Winsteps/FACETS manual):
- Theory overview: [FACETS Theory](http://www.winsteps.com/facetman/theory.htm)
- Estimation notes: [Estimation Considerations](http://www.winsteps.com/facetman/estimationconsiderations.htm)
- Measurement report structure: [Table 7](http://www.winsteps.com/facetman/table7.htm)
- Bias/interaction workflow: [Bias Analysis](http://www.winsteps.com/facetman/biasanalysis.htm)

Current implementation mapping:
- Flexible multifacet estimation (RSM/PCM, MML/JML): `fit_mfrm()` / `mfrm_estimate()`
- Bias re-estimation with fixed main effects: `estimate_bias()` / `estimate_bias_interaction()`
- FACETS-style fixed-width text reports: `build_fixed_reports()`
- APA-oriented automatic text output: `build_apa_outputs()`
- Visual warning and summary layers: `build_visual_summaries()`

Note:
- This mapping file is documentation scaffolding. As implementation choices are finalized,
  each section should be cross-checked against the corresponding FACETS manual behavior.
