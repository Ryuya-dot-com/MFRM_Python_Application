# Release Checklist

Use this checklist before tagging or pushing a public beta release.

## Scope

- [ ] Confirm this repository is being used as an independent repo, not accidentally committed as ordinary files inside a private parent workspace.
- [ ] If the parent repo should reference this project, decide between a Git submodule and a parent-level ignore rule.
- [ ] Confirm the release label remains `standalone Python beta` unless the validation scope has changed.
- [ ] Confirm `ROADMAP.md` remains the active product boundary and that no core
      feature invokes or requires an external estimation engine.
- [x] Confirm the repository-external Simulation engine loader and dead dispatch
      have been removed from `streamlit_app.py` and are covered by the
      entrypoint boundary test.

## Privacy

- [ ] Confirm no example or generated files contain confidential rating data.
- [ ] Confirm generated outputs under `validation/generated/` and retained source
      snapshots are excluded from the commit and deployment bundle; preserve the
      local research evidence.
- [ ] Confirm `.streamlit/secrets.toml` is not tracked.
- [ ] Confirm the in-app and README privacy warnings are still visible.
- [ ] Confirm README screenshots use only built-in, synthetic, or fully de-identified data.

## Deployment

- [ ] Review `DEPLOYMENT.md` before hosted deployment.
- [ ] Confirm hosted demos use synthetic, built-in, or fully de-identified data.
- [ ] Record the Streamlit Community Cloud Python version selected in Advanced settings.
- [ ] If a contextual Help adapter is marked `BROWSER_ACCEPTED`, confirm its
      stable evidence reference resolves to a reviewed record that passed
      `docs/help_browser_acceptance.md` against this exact deployed commit.

## Representative journeys

- [ ] Review the [current UX1 task matrix](docs/ux_adversarial_audit.md#next-acceptance-slice--2026-10-14-planned)
      against the candidate: ordinary and held JMLE, fixed/free-SD MML, saved
      results and input recovery, in both languages.
- [ ] Record output consistency, presentation-only state preservation, browser
      operation and actual-user understanding separately. Keep unrun cases
      explicitly open; test counts do not complete the human review.
- [ ] Confirm that publication blanks and technical CSV values follow the same
      current hold, uncertainty basis and analysis identity. Inspect affected
      rendered tables, not only whether document generation succeeds.
- [ ] Record cold/repeated timing and memory for affected journeys when their
      computation or layout changes, following the [J0 comparison conditions](ROADMAP.md#performance-priority-20261014).
      Compare the same requested work and numerical accuracy; distinguish
      first use, reuse and changed-input recomputation. Retain failures and
      timeouts, and check full-bundle time when deferring individual formats.

## Statistical Validation

- [ ] Run `make verify` or the equivalent commands below.
- [ ] Confirm `tests/test_standalone_core_boundary.py` passes.
- [ ] Confirm every new unavailable/held computation exposes a stable ReasonCode.
- [ ] Confirm saved decisions reproduce from their exact AnalysisID, evidence
      fingerprints, prespecified sensitivity plan, and source records.
- [ ] Confirm native deterministic fixtures cover the supported model and design conditions changed in this release.
- [ ] Confirm the README preview image still matches the current public-beta UI after material layout changes.

The release-check payload, default Downloads archives, deterministic demo
archives, CLI, Make targets, built-in self-test registry, and normal pytest
selection are Python-native. Dormant compatibility tests carry the
`legacy_compat` marker and are excluded from `make apptest`, `make verify`, and
GitHub CI. They must be explicitly retained, deprecated, or retired before the
product-boundary roadmap item is complete; they are not release evidence.

## Commands

```bash
python -m pip install -r requirements-dev.txt
python -m py_compile streamlit_app.py
python streamlit_app.py --doctor
python streamlit_app.py --release-check
python streamlit_app.py --self-test
make apptest
python streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv
make clean
```

## GitHub Release

- [ ] Confirm the version, release notes, license notice, and documented feature
      scope match the candidate commit.
- [ ] Record successful required CI checks for that commit on Python 3.11 and
      3.12; earlier checked boxes are not evidence for a new candidate.
- [ ] Use the agreed publication workflow for this release, preserving repository
      history. A past history-replacement approval does not authorize a new one.
- [ ] Create the tag and prerelease only after the candidate checks pass and the
      current release is authorized.
- [ ] Record the published commit and verify the deployed app against it.

## Historical initial publication record

The completed items below describe the v0.1.x publication. They are retained
for provenance, not as a procedure or authorization for the next release.
In particular, the history replacement was specific to that publication.

- [x] Choose a license before public distribution.
- [ ] Confirm `LICENSE_NOTICE.md` still states commercial use is permitted and the software is provided as-is without warranty.
- [x] Add a remote with `git remote add origin <repo-url>`.
- [x] Because `origin/main` already has a separate history, choose a publication path: branch/PR, new repository, or explicit replacement.
- [x] Explicit replacement of `origin/main` approved by the project owner.
- [x] Push with `git push --force-with-lease origin main`.
- [x] Confirm GitHub Actions pass on Python 3.11 and 3.12.
- [x] Add a tag only after CI passes.
- [x] Create GitHub prerelease `v0.1.0-beta`.
- [x] Prepare `v0.1.1-beta` notes after cross-package validation protocol changes.
- [x] Prepare `v0.1.2-beta` notes after publication figure export and Simulation validation-template changes.

## Retained research evidence

`retained_evidence` tests replay local, hash-bound study outputs. They remain
outside `make apptest` and CI because those output bundles are deliberately
not distributed. Restore the registered evidence before running
`python -m pytest -m retained_evidence tests`. Native unit tests generate
synthetic inputs in temporary directories and continue to run in CI.
A passing release gate does not certify retained research findings.
