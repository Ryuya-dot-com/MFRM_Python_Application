# R Cross-Check Status

Last local smoke run: 2026-04-12

This file is a frozen record of a local compatibility smoke check performed in
2026. It is not current release evidence and does not define a supported app
handoff or integration.

## Historical Command Status

The former `--export-parity-fixture` CLI route has been removed. The command
used for this archived run is intentionally not presented as a current workflow.
Any retained files are compatibility-only and are excluded from normal CI.

## Environment

| Package | Version | Installed |
| --- | --- | --- |
| R | 4.5.2 | TRUE |
| TAM | 4.3.25 | TRUE |
| sirt | 4.2.133 | TRUE |
| mirt | 1.46.1 | TRUE |
| tidyr | 1.3.2 | TRUE |
| dplyr | 1.2.1 | TRUE |

## Result

| Package | Check | Status | Detail |
| --- | --- | --- | --- |
| TAM | `tam.mml.mfr` | ok | Wrote `tam_xsi_facets.csv` and `tam_person.csv`. |
| mirt | `gpcm item-level` | ok | Wrote `mirt_gpcm_eap_scores.csv` and `mirt_gpcm_item_parameters.csv`. |
| sirt | `rm.facets rater fixture` | ok | Wrote `sirt_rater_facets_model.rds`, fit log, and summary text. |

## Interpretation

- The generated fixture can be handed off to TAM, mirt, and sirt in the tested
  local R environment.
- The R outputs should be used for directional validation and reviewer
  evidence, not as runtime dependencies.
- Exact equality is not expected because constraints, parameterization,
  quadrature, latent variance handling, optimizer details, and omitted
  likelihood constants differ across packages.
- Before making a public parity claim, archive the generated fixture folder,
  `r_crosscheck_status.csv`, `r_crosscheck_package_versions.csv`,
  `r_crosscheck_report.md`, and the filled
  `external_validation_report_template.csv`.
