# Deployment Guide

This app is safest for confidential rating data when run locally. Hosted deployment is appropriate for public demonstrations, teaching data, synthetic data, or non-confidential research prototypes.

## Local Deployment

Use local execution for confidential, regulated, or institutionally restricted data:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
streamlit run streamlit_app.py
```

Before sharing results, remove direct identifiers where possible and keep only the columns required for the analysis.

## Streamlit Community Cloud

Use these settings when deploying from GitHub:

- App: https://mfrmpythonapplication-xodghmcghngx8cz5zmhfeg.streamlit.app
- Repository: `Ryuya-dot-com/MFRM_Python_Application`
- Branch: `main`
- Main file path: `streamlit_app.py`
- Python dependencies: `requirements.txt`
- Development-only dependencies: `requirements-dev.txt` is for CI and local verification, not the hosted runtime.
- Secrets: none are required by the app. Do not add confidential data to `.streamlit/secrets.toml`.

Open **About MFRM / アプリ情報** to inspect the version and source commit when
available. Compare these with the adopted release commit, not an unrelated local
checkout. If they differ, inspect the deployment branch and Cloud logs before
changing code or rebooting. A missing commit display requires corroboration from
the deployment record; the version string alone does not identify the source.

The root `requirements.txt` is intentional because Streamlit Community Cloud installs app dependencies from a requirements file in the repository root or next to the entrypoint file. A `packages.txt` file is not currently needed for core analysis because the app falls back to interactive HTML figure exports if Chrome or Chromium is unavailable for Kaleido PNG export.

Set the Python version deliberately in the Streamlit Community Cloud Advanced settings. After deployment, changing the Python version requires deleting and redeploying the app, so record the chosen version with the app release notes.

## Privacy Gate Before Hosted Deployment

Do not use a hosted deployment for:

- confidential operational scoring data
- person-level records with direct identifiers
- regulated or contract-restricted assessment data
- data whose sharing policy has not been confirmed

For public demos, use synthetic data, the built-in example, or fully de-identified teaching data.

## Preflight Checks

Verify the complete candidate source in a separate directory and virtual
environment before updating the public branch. Include newly added modules and
resources as well as tracked changes: `streamlit_app.py` imports the shared
`mfrm_app/` package and reads `locales/`, bundled anchor guides, `CITATION.cff`
and selected files in `docs/`. Copying the entrypoint alone is insufficient.
The desktop application and its build requirements are not required to run the Web
app; their presence does not establish desktop distribution support.

Keep a source manifest and installed package versions with the local verification
record. Exclude `validation/generated/`, `validation/source_snapshots/`, backups,
virtual environments, secrets and local `AGENTS.md` instructions from the public
candidate. Use existing results only where their source, environment and checked
behavior still apply; rerun affected or missing checks. Human comprehension
evaluation is separate from these technical checks.

Inspect the isolated release checkout before deployment:

```bash
git status --short --branch
git branch --show-current
git rev-parse HEAD
```

Then confirm:

- The candidate and target `main` are identified separately. The candidate need
  not already be merged to be reviewed or tested.
- `git status --short` is clean.
- Local validation evidence in `validation/generated/` and
  `validation/source_snapshots/` is preserved and excluded from the deployment
  bundle. `make clean` removes development caches only.
- Development caches are excluded from the source archive.
- `.streamlit/secrets.toml` is not tracked.
- README and in-app privacy warnings still appear.

Run `make verify` when the complete required checks lack applicable evidence.
Do not repeat successful numerical studies or a complete CI matrix solely to
refresh progress text. For later changes, compare the source and dependency
manifests, run affected checks, and satisfy any required final-head checks.
`make clean` is optional checkout housekeeping, not a numerical qualification.

## Publication and hosted acceptance

Use the [hosted release checks](RELEASE_CHECKLIST.md#hosted-web-release) and
record their outcome for the specific candidate. Prepare the following before
an authorized public update:

1. Record the candidate commit, CI run, pre-update `main` commit, installed Python
   and package versions, hosted URL and deployment coordinates. Compare the
   final diff with the tested candidate. Keep local evidence outside the public
   source bundle.
2. For a new version, set `date-released` in `CITATION.cff` on its publication day; use that
   date/year consistently in the README, changelog and release body. Remove the
   unpublished-candidate wording for that version together. Check the generated
   citation and its saved-report assets. Do not change a published version’s
   date when correcting documentation, or describe a Web-only release as a
   desktop release.
3. Merge the reviewed candidate once. A `main` update is a deployment action:
   Community Cloud follows the configured repository and reinstalls dependencies
   when their files change. Record the actual merge commit and confirm the
   running source. See [Streamlit's update behavior](https://docs.streamlit.io/deploy/streamlit-community-cloud/manage-your-app#update-your-app).
4. Run the representative checks below from the public app using only synthetic
   data. Keep the input, source identity and downloaded files with the acceptance
   record. These checks address the hosted environment; they do not repeat the
   statistical validation studies.
5. After acceptance, tag the adopted commit and publish the GitHub prerelease
   with the actual date, CI link, tested scope and remaining holds. Explain that
   GitHub source archives require Python and are not desktop installers. Verify
   the public release page, source archive identity and app link. If acceptance
   fails, follow the rollback section before announcing availability.

| Hosted check | Required observation |
|---|---|
| About / citation | Adopted version and source; actual publication year in citation assets; Python and dependencies recorded from deployment evidence |
| Built-in writing example, guided RSM/JMLE | Input → analysis → results/figures → HTML/PDF report; interpretation and person-boundary holds remain visible; Japanese/English navigation retains the analysis |
| Synthetic fixed-SD RSM/PCM selected intervals | Use a recorded input/settings pair with an accepted numerical reference; obtain eligible output, save CSV/JSON/HTML, reload the matching JSON, reject the other model's record and restore opted-in history |
| Held output | Use a recorded interval counterexample that fails the applicable numerical checks and retains unavailable SEs/bounds; invalid/missing-cell G/D input explains its reason in the downloaded CSV |
| Native downloads | Read the resulting files, verify analysis identity, values/holds and readable report/figure content; saving leaves the selected result available |
| Primary controls and reflow | At desktop and 800px width, complete the sample's primary navigation and saving with the keyboard; confirm visible focus, reachable controls and retained warnings. This is not certification of every detailed screen or assistive technology |
| Static figure availability | Save a PNG when the host supports the renderer. If unavailable, show its explanation and successfully save/open interactive HTML; do not claim hosted PNG support based on a local success |

Reuse already accepted layouts and numerical reference values. A failed check
must be classified as source mismatch, application behavior, export rendering or
host environment; fix and repeat the affected route rather than restarting all
studies. Record each required check as passed, failed or not performed. A
successful push or a healthy landing page alone does not complete acceptance.

## Rollback

Before merging, retain the actual pre-update deployment commit, dependencies
and settings. Recheck the target branch immediately before publication, and
identify concurrent unrelated changes that must be preserved. Do not use
an old release’s rollback commit as the baseline for a new update.

If the deployed source is wrong, required analysis/saving fails, or holds are
lost, stop the release announcement. Revert the release merge through a new
reviewable PR, preserve history, and merge that revert to the deployment branch.
[GitHub's Revert workflow](https://docs.github.com/en/pull-requests/how-tos/merge-and-close-pull-requests/reverting-a-pull-request)
creates such a PR; conflicts or nonstandard merge history require a reviewed
manual revert in an isolated checkout. Do not reset or force-push public `main`.

Confirm the reverted application source/dependencies match the retained baseline
(the revert commit itself has a new SHA), wait for deployment, then check About,
the built-in analysis and an actual report download. Inspect Cloud logs if the
deployment remains stale. A revert cannot guarantee restoration of an old
dependency environment: compare the recorded installed versions, and do not
call recovery complete while the runtime or required behavior differs.

Keep failure logs and the rejected candidate for correction. Leave any published
tag intact; mark an already published prerelease as withdrawn/affected with its
reason and working version. Do not silently retarget its tag. Rollback restores
the earlier supported behavior; it does not validate the rejected feature.

## Official References

- Streamlit Community Cloud app dependencies: https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/app-dependencies
- Streamlit Community Cloud Python version management: https://docs.streamlit.io/deploy/streamlit-community-cloud/manage-your-app/upgrade-python
