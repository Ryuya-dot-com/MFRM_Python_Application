# Product roadmap

MFRM supports the path from rating data to analysis, interpretation and a
reproducible report. Development prioritizes usable, scientifically bounded
workflows over feature counts or a fixed release date.

For available versions and first-use instructions, see the
[English](README.md) or [Japanese](README.ja.md) introduction.
[Analysis scope](docs/analysis_scope.md) explains the current interpretation
limits. This roadmap describes direction; it is not a download announcement
or a claim that a planned feature has passed validation.

<a id="release-focus-20261028"></a>
<a id="現在の実行計画--2026-09-15"></a>

## Current delivery priority

The next distribution priority is a desktop application that includes the
complete supported workflow in one package:

| Area | Initial desktop scope |
|---|---|
| Models and estimation | RSM/PCM with JMLE or fixed-population-SD MML |
| Figures and diagnostics | Wright maps, Pathway maps, category curves, estimates and fit tables |
| Facet combinations and uncertainty | Conditional bias diagnostics, specified structural intervals and a separately fitted joint-interaction model, each with its own conditions |
| Assessment planning | Observed-score G/D studies for supported balanced, fully crossed one- or two-random-facet designs |
| Working with results | CSV input; column/rubric selection; cancellable analysis; projects that reopen without refitting; HTML, CSV, PNG and SVG output |
| Presentation | English by default, Japanese switching, explanatory and manuscript figures with retained captions |

The desktop package includes its Python runtime. Initial distribution targets
are Windows 11 x64 and macOS Apple Silicon, with acceptance and publication
performed separately for each OS. Tested OS versions, signing status and
first-launch instructions belong with each actual download. Source code,
compiled applications and a hosted Web deployment are distinct deliverables.

The initial desktop scope includes PCM, MML, bias, intervals and G/D; these are
not deferred to a reduced later edition. The desktop does not need to reproduce
every advanced Web screen to provide that scope. Native Word/PDF exports,
automatic updating and batch handling of multiple projects are outside this
initial desktop requirement.

Before distribution, the complete application must pass the applicable
[release checklist](RELEASE_CHECKLIST.md): matching source and packaged
resources, supported analysis routes, target-OS startup, saving/restoration,
and acquisition through the advertised download route. A build or CI success
alone does not establish desktop usability on a supported OS.

## Shared calculation and interpretation

Web and desktop views use the shared Python computation and output rules.
Changing language, figure style, navigating results or reopening a project
must not silently recompute the analysis or strengthen its conclusions.
Unsupported designs, failed fits, extreme scores and unidentified parameters
must remain distinguishable from valid results.

Priority corrections address incorrect numbers, mismatched saved results,
unsupported inferential claims, lost data, failed required operations and
misleading instructions. A restriction should identify the affected output,
its reason and a useful response. A method limitation must not be presented
as a data-entry error that repeated runs will resolve.

The current distinction between conditional bias diagnostics and a separately
estimated interaction remains essential. So does the distinction between
pointwise model-based intervals, finite-sample performance and simultaneous
inference. G/D uses observed scores under a declared design; it is not MFRM
person reliability. See [analysis scope](docs/analysis_scope.md) for the
conditions rather than treating all results as either unrestricted or unusable.

<a id="learning-reporting-milestones"></a>

## Guidance and reporting

Improve the connection between a research question, the input design, selected
outputs and what a report can conclude. The maintained writing example and
[reporting guide](docs/apa_manuscript_template.md) provide a concrete starting
point. Extend examples when they address a demonstrated use, rather than
requiring users to navigate every available analysis.

The criteria for these improvements are:

- The interface and saved outputs describe the same analysis, assumptions and
  limitations in English and Japanese.
- A reader can identify the relevant table or figure, its units and uncertainty
  method, and what cannot be concluded from it.
- Manuscript figures retain necessary axes and labels; removed explanatory text
  remains available as captions or report notes.
- Optional guidance can be left or revisited without losing the current work.
- Functional tests, browser/desktop checks and observed human understanding are
  reported separately. Opening a guide or completing a workflow does not prove
  comprehension or qualify a statistical claim.

Usability and comprehension reviews can inform later changes. They are not
represented as completed research, and no fixed recruitment target is a
prerequisite for the current desktop release.

<a id="performance-priority-20261014"></a>

## Performance and maintainability

Optimize a measured bottleneck in an actual input-to-report route. Compare the
same data, requested outputs, numerical accuracy and environment, separating
first use, repeated display and recomputation after changed inputs. Include
memory use, failures and full-export time where relevant.

Reuse stored results and avoid repeated calculations during presentation.
Cache identity must include inputs and settings that change the result;
faster rendering must not reuse an unrelated or stale analysis. Shared
computation belongs in `mfrm_app/`, with views responsible for presentation.

A new numerical backend or dependency needs evidence that it solves a relevant
problem without losing supported behavior, accuracy or portability. GPU/JAX,
Bayesian backends and broad architectural rewrites are not prerequisites for
shipping the current desktop workflow. Contributor setup and verification are
in [CONTRIBUTING](CONTRIBUTING.md).

<a id="tam-model-comparison-plan"></a>

## Research and possible extensions

Further uncertainty methods, multivariate G-theory, G/D intervals, broader
population models and additional assessment designs require their own
statistical and product decisions. They are not promised by the current
interface or by this roadmap.

For an extension, first specify the question, estimand, assumptions, existing
literature/manuals and unresolved implementation issue. Choose the evidence
appropriate to the claim: algebra or implementation checks, independent
reference calculations, and sampling studies address different questions.
Numerical agreement with another program is not by itself a guarantee of
coverage, robustness or estimator superiority. Comparisons must align models,
constraints, weights, population treatment, parameterization and uncertainty
methods; simulations must report failures and Monte Carlo uncertainty.

External packages can serve as methodological or numerical references. Core
application use does not require an external estimation engine. Historical
code generators, posterior viewers and compatibility fixtures do not define
new core functionality; changes to those surfaces need an explicit scope and
maintenance decision.

Research protocols and results remain in their dedicated technical records.
Their existence does not promote experimental functionality into the product
or make every unfinished study a release prerequisite. Release-specific
changes belong in [CHANGELOG](CHANGELOG.md); deployment and recovery procedures
belong in [DEPLOYMENT](DEPLOYMENT.md).
