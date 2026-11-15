# Release checklist

This is a reusable checklist for maintainers. It states acceptance conditions,
not the completion status of a particular candidate. Record the candidate,
applicable checks, results and evidence with that release. Published changes
belong in [CHANGELOG](CHANGELOG.md); current product direction is in
[ROADMAP](ROADMAP.md).

<a id="release-milestones-20261028"></a>

## Candidate and evidence

- [ ] Identify the source commit, target platform, dependency versions, build
      inputs and distributed artifact. Keep source, hosted app and desktop
      binary identities separate, even when their version labels match.
- [ ] Review the changes and dependencies against existing evidence. Reuse
      successful checks only within their source, environment and behavior
      scope; verify changed or previously untested routes.
- [ ] Satisfy required CI for the adopted candidate without duplicating runs
      solely to update documentation or status text. Do not repeat numerical
      studies when their inputs and computation have not changed.
- [ ] Record failures, skips and checks not performed. Unit tests, numerical
      studies, browser/GUI operation and human comprehension are different
      evidence; none substitutes for all the others.
- [ ] Keep private input data, secrets, local instructions, generated research
      outputs and unrelated work out of the public source and application
      bundle. Preserve the underlying validation evidence outside that bundle.

A known failure of a required operation is not resolved by calling the release
beta. Conversely, a bounded release need not qualify unrelated experimental
methods. Record what users can do and what each result can support.

<a id="milestone-m0"></a>
<a id="milestone-t0"></a>
<a id="milestone-x1"></a>

## Scope and shared behavior

- [ ] Match README, guides, Help, release notes and actual outputs to the
      candidate's model/estimator and platform scope. Do not advertise an
      unavailable download, untested OS version or unsupported inference.
- [ ] The initial desktop distribution includes RSM/PCM × JMLE/fixed-SD MML,
      Wright/Pathway/category figures, conditional bias, specified intervals,
      joint interactions and supported G/D in one application.
- [ ] Include CSV input and column/rubric mapping, cancellable background
      analysis, English/Japanese controls, explanatory/manuscript figures,
      project save/reopen and HTML/CSV/PNG/SVG exports.
- [ ] Web and desktop use the same calculation and output-qualification rules.
      Core workflows do not invoke an external estimation engine.
- [ ] Optional advanced Web features and historical compatibility surfaces
      are distinguished from the required desktop scope.

## Statistical and saved-result checks

Select ordinary cases and relevant counterexamples for each changed feature.
Do not replace statistical acceptance with a test count or an attractive plot.

- [ ] Match the estimand, fixed/estimated quantities, constraints, covariance
      and uncertainty method to the claimed interpretation. Use
      [analysis scope](docs/analysis_scope.md) and the
      [target-specific SE specification](docs/mfrm_reporting_requirements.md#se-reference-contract-20261028).
- [ ] Preserve protections for unsuccessful, missing or invalid optimizer
      termination, unidentified designs, extreme scores and parameter boundaries.
      Fresh, cached and reopened results apply the current output rules,
      including when old saved tables say Ready or Converged.
- [ ] All-minimum/all-maximum person scores do not acquire publishable point,
      fit or uncertainty cells through report generation or reopening. Preserve
      raw calculation records, explain any technical plotted locations and
      retain valid results for other facets. Person reliability must not be
      repaired by silently dropping the affected people.
- [ ] Specified fixed-SD intervals require their target-specific numerical
      checks and full covariance for differences. Failed, stale or mismatched
      interval records remain unavailable. Do not describe their eligibility
      as finite-sample coverage or simultaneous inference.
- [ ] Fixed-calibration bias diagnostics do not establish significance or
      fairness. Joint interaction remains a separate model with fixed N(0,1),
      unit weights and its supported facet-pair scope. Verify pointwise
      intervals, the recorded test family and Holm adjustment; the adjustment
      does not turn intervals into simultaneous intervals.
- [ ] G/D accepts the supported balanced, fully crossed one-/two-random-facet
      design and rejects missing/repeated cells or unsupported weights/designs.
      Distinguish signed variance estimates from zero-clipped values used in
      coefficients, and retain the observed/planned design and interval limits.
- [ ] General MML, person-score, free-SD and equivalence restrictions remain
      consistent in tables, figures, technical exports and prose. Zero or
      missing SEs, zero plug-in EB variance, stale labels and interval overlap
      do not bypass these restrictions or become equivalence/contrast tests.
- [ ] Reliability uses the intended matched estimate/uncertainty pairs and
      records its calculation basis. EAP posterior uncertainty is not silently
      treated as a structural-effect SE or used to select another formula.
- [ ] Model comparisons retain their likelihood basis, normalized input,
      constraints, population treatment, quadrature and termination context.
      Restricted IC/LR results do not produce automatic model recommendations.
- [ ] Saved analyses and decisions retain exact identity, data/settings,
      provenance and relevant sensitivity references. Viewing, language/style
      changes, exporting and restoration do not refit or relabel another run.
      Changed settings are clearly separated from the displayed prior result.
- [ ] Reopening unsupported or altered files fails clearly without overwriting
      the original. Do not weaken source/hash checks to make old intervals appear
      eligible after implementation changes.

## User journeys, reports and privacy

- [ ] Complete input → settings → analysis → figures/tables → saving for each
      affected route. Inspect actual saved files; generation success does not
      establish that the browser or file dialog saved them correctly.
- [ ] Inspect ordinary, held, unsupported and failed results in both languages.
      Messages identify the affected result, reason and useful action. A method
      restriction is not described as a correctable input error.
- [ ] Validate retained rows, rubric bounds, unused categories and column roles.
      Distinguish observation counts, distinct people and likelihood-weight sums.
      Check duplicates, excluded rows, missing IDs and malformed input.
- [ ] Navigation, optional guidance, basic/detailed views and output selection
      preserve results and relevant settings. Saved histories remain accessible
      after inputs change; restoration does not add a new fit.
- [ ] Check primary keyboard navigation, focus, narrow windows, display scaling
      and reachable controls. Retain warnings after navigation; visual or
      automated checks alone do not establish comprehension or full accessibility.
- [ ] Keep titles, units, table notes, analysis identity and interpretation
      limits consistent across requested formats. Inspect affected Word/PDF
      pagination and HTML/figure rendering rather than only file existence.
- [ ] Manuscript figure mode retains axes, legends, labels and unavailable-result
      messages; removed explanations remain in captions/report notes. Do not
      describe this switch as comprehensive APA compliance.
- [ ] Explanations do not equate Strata with actual groups, underfit with a
      decision to exclude someone, an interaction with unfairness, or descriptive
      fit conventions with universal acceptance thresholds.
- [ ] Examples and screenshots use synthetic or appropriately de-identified
      data and match the described interface. Before real-data upload or sharing,
      explain hosted/local processing and which saved files contain input data
      or identifiers. Analysis fingerprints are not anonymization.
- [ ] Errors and support artifacts do not disclose confidential rows or secrets.
      Missing or withheld computations retain their stable diagnostic reason.
- [ ] If contextual Help is labelled `BROWSER_ACCEPTED`, its reviewed evidence
      applies to the deployed source and the
      [Help acceptance protocol](docs/help_browser_acceptance.md). Help must
      return to the correct analysis rather than a stale context.

<a id="milestone-w1"></a>
<a id="milestone-w2"></a>
<a id="milestone-w3"></a>

## Python library distribution

- [ ] Match the `mfrmpy` public guide, metadata and API to the desktop analysis
      scope; distinguish package API version from the shared engine version.
- [ ] Build both sdist and wheel, with the wheel built from that sdist; pass
      `twine check` and inspect the included modules, license and resource files.
- [ ] Install the wheel in a fresh numerical-only environment, outside the
      checkout; run `tests/test_python_library.py` with isolated Python.
- [ ] Check RSM/PCM × JMLE/fixed-SD MML, qualified supplementary results,
      one-/two-facet G/D, invalid/unsupported requests and zero-refit .mfrm replay.
- [ ] Check the figures extra and HTML export using the same built wheel.
- [ ] Complete the declared Python/OS matrix and minimum-dependency checks;
      record which environments actually ran. Do not infer platform support
      from a pure-Python wheel tag alone.
- [ ] Confirm that affected existing application workflows still work. Separate
      pre-existing CI failures from regressions caused by the library changes.
- [ ] Configure the appropriate PyPI/TestPyPI Trusted Publisher and GitHub
      environment. Ordinary pushes and PRs must only build/test, never upload.
- [ ] Publish the reviewed candidate only by an explicit manual action. Install
      from the advertised index afterward and record the actual release URL.

Initial preparation is based on application commit
`5b2c20546eb94b8a846bd48a7da0cdb031377739` (2026-10-30). Its Python 3.11 and 3.12
CI each recorded 4,161 passing, 14 skipped, 91 deselected and 3 failing tests:
two outdated button-budget assertions and one AppTest export timeout.
The engine self-tests and other preceding checks passed; benchmark/export and
clean-checkout steps were skipped after the failures. These are baseline
records, not a claim that the whole candidate passed.
[Baseline run](https://github.com/Ryuya-dot-com/MFRM_Python_Application/actions/runs/41133816447).

Local preparation evidence (2026-11-15; macOS arm64, Python 3.12.13;
NumPy 2.5.4, pandas 3.1.1, SciPy 1.18.1):

- Wheel built from sdist; both artifacts passed `twine check`. Archive inspection
  confirmed the shared engine and resources match their source files.
- Installed core, isolated Python outside the checkout: 14 passed; 4 optional
  figure tests skipped because Matplotlib was absent. Import/use was also
  tested with UI and plotting imports explicitly blocked.
- Installed optional figures/reports (Matplotlib 3.11.2): all 4 model/estimator
  cases passed, with refitting blocked during plotting/report generation.
- Existing desktop project/workflow/output/packaging tests: 101 passed.
- Existing representative Web scenarios: 2 passed after the documented budget
  correction. The pre-existing Web/desktop archive integration test passed
  locally without extending its AppTest limits; remote timing remains to be
  checked in CI.

These are local package/application checks, not new sampling validation or
interactive desktop OS acceptance. Later CI and publication evidence follows.

Published library `mfrmpy` 0.1.0a1 (2026-11-15):

- Adopted source: `c2cfae046368cfd212fa04a53b4a5e7e4a1fe2ea`, merging
  [PR #7](https://github.com/Ryuya-dot-com/MFRM_Python_Application/pull/7).
  Its tree matches approved candidate `5b7b701796df5805149ddbb2c2e56257da343b05`.
- [Full application CI](https://github.com/Ryuya-dot-com/MFRM_Python_Application/actions/runs/43689546103)
  passed 4,185 tests on each Python version, with 14 skipped and 91 deselected;
  the one failure was the existing 45-second export test limit. After removing
  that limit, [remaining-checks CI](https://github.com/Ryuya-dot-com/MFRM_Python_Application/actions/runs/43691725429)
  passed archive interoperability, benchmark, demo report export and clean-tree
  checks on Linux / Python 3.11.18 and 3.12.16. The application and numerical
  sources were unchanged. This combines evidence from separate runs; the earlier
  failed run is retained, not relabelled as an entirely successful full run.
- [Windows desktop source checks](https://github.com/Ryuya-dot-com/MFRM_Python_Application/actions/runs/43689546117)
  passed all 252 tests. No desktop binary was replaced by this library release.
- [TestPyPI publication](https://github.com/Ryuya-dot-com/MFRM_Python_Application/actions/runs/43692848976)
  and [PyPI publication](https://github.com/Ryuya-dot-com/MFRM_Python_Application/actions/runs/43693328216)
  each passed all seven package environments: Linux/macOS/Windows × Python
  3.11/3.12, plus Linux 3.11 with NumPy 1.24.0, pandas 2.0.0, SciPy 1.10.0 and
  Matplotlib 3.8.0 retained throughout the optional-figures check. Each environment
  passed 14 core and four optional figure/report tests.
- Both indexes serve the same tested wheel and sdist. Fresh, isolated Python
  3.12 installations from each index matched the wheel's installed file contents
  and completed the 960-row example, qualified tables and zero-refit .mfrm
  restoration. The PyPI installation also produced a PNG and HTML report with
  the `figures` extra; dependency consistency checks passed.
- Public release: [PyPI](https://pypi.org/project/mfrmpy/0.1.0a1/);
  trial release: [TestPyPI](https://test.pypi.org/project/mfrmpy/0.1.0a1/).
  The main-branch-only GitHub environments use the registered Trusted Publishers;
  ordinary pushes and PRs still cannot publish a package.

Published SHA-256:

| Artifact | SHA-256 |
|---|---|
| `mfrmpy-0.1.0a1-py3-none-any.whl` | `de04bb2564a5853c937cffbb42270b876bc39b198d7dfc31943b6b88b094323e` |
| `mfrmpy-0.1.0a1.tar.gz` | `78364188f5774e1b83f2ac5b2fd0f84ffc91c86cc61fb21a9511ef185747e1d8` |

The subsequent hosted interoperability check found a gap in the first
package's dependency declaration: a real Web `.mfrm` contained Arrow-backed
strings and could not be opened without PyArrow. pandas 2 also lacks the
`StringDtype(na_value=...)` constructor used by the existing codec. The same-
environment round trips above did not cover that case. Version 0.1.0a2 declares
the needed dependencies, with no change to `mfrm_app`, its numerical algorithms,
content/analysis fingerprints or archive format.
All 96 shared `mfrm_app` files in the two wheels are byte-identical.

Local correction evidence (2026-11-15):

- Source `d0e8c2fae6cefac215729b59cfc63fea57fea632`; wheel built from sdist,
  both passing `twine check`. Fresh installed core: 18 passed, four optional
  tests skipped. Installed minimum dependencies (Python 3.11.14, NumPy 1.26.0,
  pandas 3.0.0, SciPy 1.10.0, PyArrow 16.0.0, Matplotlib 3.8.0): all 22 passed.
- Public Web source `c2cfae0`, synthetic 960-row writing analysis: completed
  analysis, downloaded a 500,827-byte HTML report retaining dimensionality and
  stopping notes, and a 2,001,781-byte `.mfrm` project. Both minimum and current
  installed environments restored and re-saved that actual archive with nine
  qualified tables, unchanged analysis ID/content fingerprint and no refitting.
  These checks cover the affected hosted analysis/report/archive route, not all
  advanced Web features or new desktop binary acceptance.

Published correction `mfrmpy` 0.1.0a2 (2026-11-15):

- [TestPyPI workflow](https://github.com/Ryuya-dot-com/MFRM_Python_Application/actions/runs/43694493879)
  and [PyPI workflow](https://github.com/Ryuya-dot-com/MFRM_Python_Application/actions/runs/43694850850)
  passed all seven environments, each with 18 core tests and four figure/report
  tests. The minimum-dependency environment retained the declared versions.
- [PyPI](https://pypi.org/project/mfrmpy/0.1.0a2/) and
  [TestPyPI](https://test.pypi.org/project/mfrmpy/0.1.0a2/) serve wheel and sdist
  identical to the local and CI-tested files. A fresh TestPyPI installation
  completed the sample analysis, zero-refit restoration and the actual hosted
  project's nine qualified tables, with matching installed wheel contents.
- A fresh production-PyPI install of `mfrmpy[figures]==0.1.0a2` also passed
  those checks, rendered PNG/HTML output, and passed `pip check`. The hosted
  download and both index-install records are retained under
  `validation/generated/pypi-preparation-20261115/` in the local workspace.

| Artifact | SHA-256 |
|---|---|
| `mfrmpy-0.1.0a2-py3-none-any.whl` | `4c7df8838578189a51396225afd771f099a585c678a26d2edd1a407c3532887b` |
| `mfrmpy-0.1.0a2.tar.gz` | `482735d96a074879f88a7b4e2337094710b5c037c80852057353b9b6227fa529` |

## Hosted Web release

- [ ] Follow [DEPLOYMENT](DEPLOYMENT.md#publication-and-hosted-acceptance).
      Include all shared modules, locale/assets and required runtime dependencies.
- [ ] Record the configured deployment repository/branch/entrypoint, installed
      Python/package versions and recoverable pre-update source/environment.
- [ ] Confirm the running version/source and citation match the adopted release.
      A version string, successful push or healthy landing page is insufficient.
- [ ] Use synthetic ordinary and relevant held examples for bilingual analysis,
      selected intervals, reports, actual downloads, matching-record reload,
      incompatible-record rejection and opted-in history restoration.
- [ ] Check the host's actual static-image capability. If PNG rendering is
      unavailable, explain it and verify the offered HTML figure opens correctly.
- [ ] Verify the recovery procedure and preserve downloaded results before
      planned restarts. Do not advertise the update while a required hosted
      operation or source match remains unresolved.

<a id="milestone-t1"></a>

## Desktop build and bundled resources

- [ ] Build on the target OS/CPU using the reviewed source and build environment.
      Match app version, platform metadata and the included build manifest to
      the executable, Python runtime, shared modules, translations and assets.
- [ ] Include required dependency/runtime license and notice files, including
      Python and Tcl/Tk as applicable. Check the actual bundle, not filenames
      alone or the existence of the application's own license.
- [ ] Include a self-contained first-use guide matching the packaged app.
      A changed START_HERE or other resource creates a different archive;
      update the artifact hash and distribution record accordingly.
- [ ] Determine the minimum OS from all bundled dependencies and actual tests.
      The executable's own minimum alone does not establish the bundle's support.
- [ ] Exercise the frozen executable's worker and supported analyses, exports
      and independent restoration without re-estimation. Keep source-test and
      frozen-app evidence distinct.

Build entrypoints are [scripts/build_desktop.py](scripts/build_desktop.py) and
[scripts/desktop.spec](scripts/desktop.spec). The
[Windows workflow](.github/workflows/desktop-windows.yml) runs source and frozen
checks; [check_windows_desktop.py](scripts/check_windows_desktop.py) prepares
and checks the Windows archive. Consult those inputs when changing packaging;
GitHub-generated source archives are not application installers.

<a id="milestone-t2"></a>
<a id="milestone-t3"></a>

## Desktop target-environment acceptance

- [ ] On each claimed OS/CPU, use native dialogs for CSV input and actual project,
      HTML, CSV, PNG and SVG saving. Include Japanese paths and identifiers,
      keyboard use, relevant window sizes and display scaling.
- [ ] Check cancellation and successful retry, worker failure, unsaved-work
      handling, input changes, read-only inputs and unwritable destinations.
      Preserve the prior successful result and existing files on failure.
- [ ] Quit and restart the application, then reopen a saved project without
      refitting. Test previous supported project formats, current extended
      projects and unsupported/altered files. Keep originals for rollback.
- [ ] Verify startup, analysis and saving in a target environment without a
      development Python installation. A modified PATH on a development machine
      is narrower evidence than a clean-environment test; record which was used.
- [ ] Verify normal analysis, figures, local HTML and projects offline. Explain
      that online Help requires a connection; do not make a broader offline claim.

<a id="milestone-t4"></a>

## macOS acquisition and startup

- [ ] Download the intended archive through the intended distribution route,
      verify its hash, extract and launch it with normal download quarantine
      behavior. Check analysis, saving and project reopening on that artifact.
- [ ] State architecture, tested OS range and actual signing/notarization status.
      Distinguish ad-hoc signing from Developer ID signing and notarization.
- [ ] Document only the applicable per-app first-launch procedure. Do not make
      disabling OS protection or bulk quarantine removal a normal installation
      step. An environment that cannot open the app is not accepted by a note.
- [ ] Check the supported Finder/document-opening route, including when the app
      is already running, alongside the in-app project-opening route.

<a id="milestone-t5"></a>

## Windows acquisition and startup

- [ ] Download and extract the complete application ZIP on Windows 11 x64;
      retain the EXE and required adjacent resources. Verify the hash and
      record any signing/security prompt and the actual launch result.
- [ ] Complete the applicable GUI and clean-environment checks above on Windows
      11. A Windows Server CI run or Mac source check is not Windows 11 GUI
      acceptance. Reuse earlier real-PC evidence only for the behavior and
      artifact scope it supports.
- [ ] Verify in-app project opening and explain whether file association is
      installed. Do not require Python installation for the bundled application.

<a id="milestone-t6"></a>

## Publication, updating and recovery

- [ ] Match version, citation/date, license, scope, tested environments, signing
      status, known limits, startup instructions and SHA-256 to the actual release.
      Do not predate an unpublished version or replace the date of an old version.
- [ ] Publish the reviewed source/tag and accepted OS artifacts through the
      authorized release process. GitHub can host the announcement with a
      version-specific external download link; verify that actual destination.
- [ ] Retrieve the publicly offered artifact, match its hash to the accepted
      one, and confirm the applicable startup/save/reopen route. A source ZIP
      or successful file upload is not evidence of a working app download.
- [ ] Keep the previous app and project originals available during an update.
      State that new extended projects cannot be opened by older basic-only
      apps; test supported forward migration without overwriting originals.
- [ ] Keep a recovery path to the previous supported behavior. Preserve public
      Git history and published tag identity; withdraw or correct release notes
      rather than silently retargeting a published tag.
- [ ] Verify the final public entry points, release body and downloads. Each
      OS can be released once its own requirements pass; simultaneous OS
      publication is not required.

## Verification commands and research records

Use [CONTRIBUTING](CONTRIBUTING.md#verification) for source verification commands
and the workflows for required platform checks. The static `--release-check`
checks prerequisites; Ready is not a release authorization or proof of numerical
validity, hosted operation or target-OS acceptance.

Normal tests and default exports are Python-native. `legacy_compat` tests cover
dormant external-engine compatibility and are excluded from normal verification;
they are not release evidence. `retained_evidence` tests require separately
retained, hash-bound study outputs, including timestamp-controller trust roots
and receipts where applicable. These assets are not distributed in the public
checkout. Restore their registered inputs before running that selection; do not
count those tests as a clean-checkout pass without the required environment.
A passing product release check does not certify every retained research claim.

Human usability/comprehension findings should identify the evaluated version,
tasks, assistance and observed errors. They can guide subsequent improvements;
no participant recruitment count is a gate for this release. This does not
excuse known technical errors in the required user journeys.
