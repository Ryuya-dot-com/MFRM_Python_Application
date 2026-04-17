# Security and Data Privacy

## Supported Versions

This repository is currently a public-beta research preview. Security and privacy fixes should target the `main` branch first and be included in the next supported beta tag when applicable.

## Reporting

Do not include confidential rating data in public issues, pull requests, screenshots, logs, or exported files.

If this project is hosted on GitHub, prefer a private GitHub security advisory for security-sensitive reports. If no private reporting channel is available yet, contact the repository owner privately before sharing details.

## Data Handling Principles

- Uploaded or pasted rating data can include person IDs, rater IDs, institution labels, and subgroup variables.
- The app is intended to keep user-provided data in memory unless a user explicitly downloads or exports outputs.
- `st.cache_resource` is used only for the read-only core namespace.
- `st.cache_data` is used for built-in sample data, bundled templates/guidelines, and bounded short-lived export bytes keyed by fingerprints.
- Reproducibility fingerprints are not encryption and are not a privacy guarantee.
- Hosted deployments can process uploaded data on remote infrastructure. Confidential data should be analyzed locally on controlled machines.

## Posterior Viewer Uploads (v0.2.0 and later)

The Posterior Viewer app mode accepts externally-produced posterior draws in three formats: CmdStan CSV (one file per chain), Apache Parquet (any DataFrame with `chain`, `draw`, and one column per parameter), and ArviZ NetCDF (`.nc`). A few points that are distinct from the rating-data upload path:

- Posterior-draw files are treated as user-generated analysis artefacts rather than raw rating data. They still go through the `st.file_uploader` path and obey the same in-memory-only handling principle above.
- Do not upload posterior files that embed person / rater identifiers in their parameter names (e.g. `theta_R1234`). The Posterior Viewer renders parameter names verbatim in plots, summary tables, and downloads.
- The Posterior Viewer does not compile or sample Stan models; it only ingests draws produced elsewhere. Source your draws from trusted tooling (CmdStanPy, cmdstanr, rstan, PyMC, etc.) — uploading an unknown `.nc` file means importing whatever numerical content the producer chose.
- ArviZ and NetCDF4 are required runtime dependencies for this mode (≈ 80 MB combined). Refer to `requirements.txt` for pinned versions.

## Advanced-Model Stan Code Downloads (v0.2.0 and later)

The sidebar's "🧪 Advanced models (Stan, download only)" expander emits Stan programs for DINA, HRM, Testlet RI / Bifactor, Mixture Rasch, 2PL Binary, and Pairwise BTL. The app never compiles, samples, or executes these programs in-process — they are downloaded as `.stan` files and must be run locally with cmdstan / cmdstanpy / rstan. The Q-matrix uploader used by DINA treats the uploaded file as a configuration document, not as rating data, but users should still avoid embedding identifying labels in attribute / item column names.

## Before Sharing Reproducible Examples

- Replace person/rater/institution IDs with synthetic labels.
- Remove columns that are not needed to reproduce the issue.
- Prefer the built-in sample data or a small simulated dataset.
- Confirm that `validation/generated/` and exported ZIP files do not contain sensitive information before attaching them anywhere.
