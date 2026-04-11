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

## Before Sharing Reproducible Examples

- Replace person/rater/institution IDs with synthetic labels.
- Remove columns that are not needed to reproduce the issue.
- Prefer the built-in sample data or a small simulated dataset.
- Confirm that `validation/generated/` and exported ZIP files do not contain sensitive information before attaching them anywhere.
