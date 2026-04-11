# Simulation Reference Status

Last local inspection: 2026-04-12

This file records what the sibling/local Simulation directory can contribute to
numerical validation. It is not a data bundle and should not be treated as a
runtime dependency for the Streamlit app.

## Local Evidence Inventory

| Reference set | Artifacts | Local inspection result | Validation use |
| --- | --- | --- | --- |
| Observed long-form rating data | `data/writing_long.csv`, `data/speaking_long.csv` | writing: 111,995 rows x 8 columns; speaking: 111,993 rows x 8 columns | Realistic import, score support, missingness, and category-distribution checks |
| Main engine refit sweep | `results/engine_refit_manifest.csv` | 4,004 rows; 1,001 rows each for mfrmr JMLE, mfrmr MML, Python JMLE, and Julia JMLE; all inspected rows marked `ok` | Large repeated-refit regression evidence |
| Runtime benchmark sweep | `results/runtime_validation_summary.csv`, `results/runtime_validation_long.csv`, `results/runtime_fullref.csv` | 4 summary rows, 4,000 validation runtime rows, and 4 full-reference runtime rows | Performance and convergence smoke evidence |
| Full-reference backup outputs | `results_backup_20260407_pre_diag/*_full_ref/*.csv` | R/mfrmr JMLE and MML, Python JMLE, Julia JMLE, and FACETS-style summary/person/facet outputs are present | Small auditable cross-engine comparison set |
| Diagnostic spot checks | `results_diagcheck_all/`, `results_julia_diag/`, `results_mmlcheck/`, `results_backup_20260407_pre_diag/engine_refit_manifest.csv` | Inspected manifests had 3, 1, 2, and 16 `ok` rows respectively | Targeted checks for diagnostics, Julia output, and mfrmr MML/JMLE paths |
| Validation input replicates | `results/validation_inputs/rep_*_writing_long.csv` | Some large replicate CSVs are present; inspected examples included a non-empty `rep_010` with 83,996 lines and zero-line placeholders for `rep_001` and `rep_100` | Parser stress tests and empty-input handling |

## Public Reporting Rule

Do not copy raw observed data, generated replicate data, or historical manifests
with local absolute paths into the public repository. For public validation,
archive sanitized summary tables, package versions, parameterization notes, and
the tolerance policy before making any numerical-validation claim.

The app exports the same policy as
`external_simulation_reference_inventory.csv` through the parity fixture and
demo/report downloads.
