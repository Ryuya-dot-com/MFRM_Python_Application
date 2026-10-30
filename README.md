<p align="center">
  <img src="docs/images/mfrm-desktop-icon.svg" width="92" height="92" alt="MFRM Desktop icon" />
</p>

<h1 align="center">MFRM — Fairer Rating Analysis, Visually Explained</h1>

<p align="center">
  <strong>Separate examinee ability from rater severity, task difficulty, and rubric behavior on one shared scale.</strong><br />
  <a href="https://ryuya-dot-com.github.io/MFRM_Python_Application/"><strong>Website &amp; Docs</strong></a> ·
  <a href="https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-macos-20261030"><strong>Download for macOS</strong></a> ·
  <a href="https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-windows-20261030"><strong>Download for Windows</strong></a> ·
  <a href="https://mfrmpythonapplication-xodghmcghngx8cz5zmhfeg.streamlit.app/"><strong>Open Web App</strong></a> ·
  <a href="docs/desktop_user_guide.en.md"><strong>User Guide</strong></a> ·
  <a href="README.ja.md"><strong>日本語</strong></a>
</p>

---

## Why Use MFRM? A Concrete Example

Imagine you are collecting essay, speaking, or clinical OSCE scores for a thesis or program evaluation on a **0–5 rubric**:

| Examinee | Assigned rater | Assigned prompt | Raw mean score | Who looks stronger from raw averages? |
|---|---|---|---|---|
| **Student A** | **Rater 4** *(severe)* | **Task 2** *(complex argumentative prompt)* | **3.2 / 5.0** | Looks weaker in Excel (`3.2`) |
| **Student B** | **Rater 1** *(lenient)* | **Task 1** *(familiar personal prompt)* | **3.6 / 5.0** | Looks stronger in Excel (`3.6`) |

Why can't raw averages—or even a high inter-rater correlation (`r = 0.90`)—solve this?

1. **Raw averages mix examinee ability with luck of the draw:** Student A was graded by a tougher rater on a harder prompt.
2. **Correlation ignores constant severity gaps:** Two raters can rank students in nearly the same order while one consistently scores a full point lower than the other.
3. **Rubric numbers are not automatically equal steps:** Moving from `1 → 2` on *Content* is rarely the same jump as moving from `4 → 5` on *Mechanics*.

**Many-Facet Rasch Measurement (MFRM)** estimates examinee ability, rater severity, task difficulty, and criterion difficulty together on a single **logit ruler** (`0.00` = average facet location), while **Generalizability / Decision (G/D) studies** show how changing the number of raters or tasks affects score reliability.

**MFRM Desktop** brings these workflows into one self-contained app (`Data → Settings → Results`)—with **Python bundled inside and no FACETS, ConQuest, R, or command-line setup required** in the packaged macOS and Windows previews.

![Wright map from the built-in synthetic writing sample](docs/images/desktop-wright-map.png)

> **Reading the built-in writing sample (`Data → Sample → Next → Analyze`) at a glance:**
> - **Left (`Persons · n = 30`):** Examinee ability estimates span roughly `-2.4` to `+2.4` logits (higher = higher ability).
> - **Right (`Rater`, `Task`, `Criterion`):** `R4` (`+0.47` logits) grades much more severely than `R1` (`-0.73` logits—nearly a **1.2-logit** severity gap!); `Task2` (`+0.31`) is harder than `Task1` (`-0.31`); and earning high marks on `Mechanics` (`+0.28`) is harder than on `Content` (`-0.56`).
> - **Dotted horizontal lines:** Common rating-scale step thresholds across the 0–5 rubric. Because student locations and facet/step locations overlap vertically, the test is **well targeted** to this group.

---

## Choose How to Run MFRM

| Platform | Best for | Current availability & how to start |
|---|---|---|
| **macOS (Apple Silicon) Desktop** | Offline local analysis, `.mfrm` project files, manuscript PNG/SVG exports on M1 or later Macs — **no Python install needed** | **Public preview available now.** [Download from GitHub Releases (`desktop-macos-20261030`)](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-macos-20261030) and follow the [Desktop Guide](docs/desktop_user_guide.en.md). Ad-hoc signed / unnotarized preview (`0.2.19-beta`, macOS 14.0+); Intel Macs can run from Python source. |
| **Windows 11 (x64) Desktop** | Same offline local `Data → Settings → Results` workflow on Windows — **no Python install needed** | **Public preview available now.** [Download from GitHub Releases (`desktop-windows-20261030`)](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-windows-20261030) (application ZIP hosted on Dropbox) and follow the [Desktop Guide](docs/desktop_user_guide.en.md). Unsigned preview (`0.2.19-beta`); Windows on ARM and other Windows editions are not yet verified. |
| **Web Browser (Any OS)** | Instant browser exploration, plus extra research modules (bounded GPCM, latent regression, anchors/linking, residual PCA, simulation) | **[Open the hosted Web app](https://mfrmpythonapplication-xodghmcghngx8cz5zmhfeg.streamlit.app/)** (`v0.2.19-beta`) with synthetic or appropriately de-identified data, or [install and run locally](#install) for confidential datasets. |
| **R (`mfrmr` on CRAN)** | Reproducible scripted MFRM pipelines in R, RStudio, Positron, or Quarto | Install via `install.packages("mfrmr")` from [CRAN](https://cran.r-project.org/package=mfrmr) and see the [`mfrmr` documentation](https://ryuya-dot-com.github.io/mfrmr/) (independent R implementation; align settings before comparing outputs). |
| **Contributors & Builders** | Running tests, checking build manifests, or packaging desktop prototypes | See [CONTRIBUTING.md](CONTRIBUTING.md) and [DEPLOYMENT.md](DEPLOYMENT.md). |

> [!NOTE]
> **Downloading the right file:** On GitHub Releases, GitHub's automatic **Source code (zip)** contains Python source code, not the standalone desktop app. Use the application ZIP on [`desktop-macos-20261030`](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-macos-20261030) (macOS Apple Silicon) or [`desktop-windows-20261030`](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-windows-20261030) (Windows 11 x64). The separate `v0.2.19-beta` tag is the Web/source release. In the Web app, check **About MFRM** for the live commit and version.

---

## Your First Analysis in 60 Seconds

MFRM Desktop guides you through three screens—**Data → Settings → Results**—with an instant **English / 日本語** selector in the top bar:

```mermaid
flowchart LR
    A["1. Data<br/>Load built-in Sample,<br/>your UTF-8 CSV,<br/>or a saved .mfrm project"] --> B["2. Settings<br/>Map Person / Score / Facets,<br/>set rubric min–max, and pick<br/>MFRM or G/D options"]
    B --> C["3. Results &amp; Save<br/>Read figures &amp; tables,<br/>switch Explanatory / Manuscript mode,<br/>save .mfrm, .html, .png/.svg, .csv"]
```

- **Desktop quick start:** Click **Data → Sample → Next → Analyze**. Explore **Figure** (`Wright map`, `Pathway map`, `Category probabilities`—click **Reading guide** next to any chart for plain-English help) and **Table**, then save a `.mfrm` project or `.html` report. Follow the **[Desktop User Guide (Windows & macOS)](docs/desktop_user_guide.en.md)** for a full walkthrough.
- **Browser quick start:** Open the [Web app](https://mfrmpythonapplication-xodghmcghngx8cz5zmhfeg.streamlit.app/), choose the sample tour or **Continue without the guide**, run **Built-in example → Writing essay**, read **First Read**, and export from **Report & Export**.

---

## What Results Look Like — and How to Interpret Them

Once an analysis finishes, here is how to connect what you see on screen to concrete assessment conclusions (and what **not** to overclaim in a thesis or manuscript):

| Output | If you see this pattern… | What it tells you in practice |
|---|---|---|
| **Wright map** | `Persons` (left) line up vertically with `Rater` / `Task` / `Criterion` points and dotted step lines (right) | **Good targeting:** tasks and rubric steps match the ability range of your examinees. If almost all persons sit far above or below the facets, the tasks were too easy (ceiling) or too hard (floor). |
| **Pathway map & Fit (`Infit MnSq`)** | A severe rater like `R4` (`+0.47` logits) sits near the vertical `Infit MnSq = 1.00` line (`Infit = 1.05`) | **Predictably severe vs. erratic:** Severity (vertical) is distinct from consistency (horizontal). Near `1.00` means consistent application of standards; far right (`> 1`, *underfit*) flags noisy, erratic ratings to inspect; far left (`< 1`, *overfit*) suggests muted/central-tendency scoring. Never drop a rater or examinee from a fit value alone. |
| **Category probability curves** | Each score (`0`, `1`, `2`, `3`, `4`, `5`) forms its own distinct hill from left to right | **Healthy rubric scale:** Raters use every category to distinguish ascending performance bands. If a middle curve stays flat or buried underneath neighbors, check **Observed Category Use** and rubric descriptors (or test **PCM** if criteria use steps differently). |
| **Conditional bias heatmap** | Most `Facet A × Facet B` cells stay near `0.00`, while one cell shifts noticeably (`+` or `−` logits) | **Exploratory diagnostic prompt:** Flags a rater-by-task or rater-by-criterion pairing to review during rater training, with main calibration held fixed. It does **not** prove statistically significant bias or unfairness. |
| **Specified intervals / Joint interaction** | Horizontal interval bars appear for eligible fixed-SD MML facet levels, `A − B` contrasts, or refitted interaction cells | **Scoped model-based uncertainty:** Shows pointwise approximate intervals where numerical checks pass (using full covariance for contrasts; Holm adjusts the recorded family of joint-interaction tests, not the intervals). Blank/withheld rows indicate unmet checks—never zero uncertainty. |
| **D-Study (`G` & `Phi` curves)** | Vertical jumps across `Task = 1..4` lines vs. the slope along `Rater count = 1..6` show how each design reaches higher `G` / `Phi` | **Actionable design trade-off:** Compare equal-effort designs (e.g., `2 tasks × 3 raters` vs. `1 task × 6 raters`) alongside the variance-component table to see whether more prompts or more raters reduce error faster (`G` = relative rank; `Phi` = absolute cutoff; intervals withheld). |

| Rubric Category Functioning (`RSM`) | Rater & Task Count Planning (`D Study`) |
|---|---|
| ![Category probability curves for a 0-5 writing rubric](docs/images/desktop-category-curves.png) | ![D-Study G and Phi projections across rater and task counts](docs/images/desktop-d-study.png) |
| *Built-in writing sample (`0`–`5`): six ordered peaks show each rubric level becoming most likely in turn along the logit scale.* | *Synthetic `GD_ANOVA.mfrm` example: projects observed-score `G` (solid) and `Phi` (dashed) across 1–6 raters and 1–4 tasks.* |

---

<a id="model-scope"></a>

## Model Scope

MFRM Desktop integrates core Rasch calibration, rubric diagnostics, interaction checks, selected fixed-SD intervals, and balanced G/D design planning in one application:

| Research or operational task | Available workflow |
|---|---|
| **Estimate examinee measures & facet effects** | **RSM** (shared rubric steps) or **PCM** (facet-specific steps), with **JMLE** or **fixed-population-SD MML** (EM) |
| **Inspect scale targeting & fit** | Wright maps, Pathway maps (`Infit`, `Outfit`, `Approx. Infit ZSTD`), category probability curves, estimate & fit tables |
| **Examine combinations of facet levels** | Fixed-calibration **Conditional bias** diagnostics; Desktop also includes a separately refitted fixed-`N(0,1)` MML **Joint interaction** model |
| **Quantify selected facet uncertainty** | **Specified pointwise intervals** for numerically eligible fixed-SD MML facet levels and prespecified `A − B` differences |
| **Compare rater & task count designs** | Observed-score **G/D studies** for completely crossed, balanced 1- or 2-random-facet designs |
| **Save, restore & publish** | Web reports/downloads; Desktop `.mfrm` projects (reopened without refitting), offline `.html` reports, `.csv`, `.png` (300 DPI), and `.svg` |

The Web interface adds bounded GPCM, latent regression, anchors/linking, residual PCA, and simulation. Read [Analyses, interpretation and reporting limits](docs/analysis_scope.md) before using any output to support a statistical claim.

---

## Preparing Your Data (`CSV` Format)

Arrange ratings in **long format**—**one row per single rating event** (one examinee rated by one rater on one task and criterion):

```csv
Person,Rater,Task,Criterion,Score
P01,R1,T1,Content,4
P01,R1,T1,Organization,3
P01,R2,T1,Content,5
P02,R1,T1,Content,2
```

*(Practice with the full 960-row [synthetic writing CSV](examples/desktop_writing.csv) or click **CSV format → Save example CSV** inside the desktop app.)*

- **Three essentials:** (1) Use integer ratings for MFRM rather than averaged decimals, and reshape wide `Rater1, Rater2` columns into rows; (2) enter the rubric's **full intended minimum and maximum** (e.g., `0` and `5`) even if an extreme category was never observed; (3) **never enter `0` for an unrated cell**—leave missing scores blank or omit the row (`0` means an actual score of zero).
- **Supported formats:** Desktop reads **UTF-8 CSV**. Web also accepts spreadsheet paste, TSV/TXT, Excel, Parquet, and JSON. See [Use your own CSV](docs/desktop_user_guide.en.md#use-your-own-csv) for column mapping details.

---

## Statistical Caveats & Scope Boundaries

A completed run does not make every displayed number suitable for statistical inference. Always check the message attached to the specific table or figure:

- **General MML tables, Person-score uncertainty, estimated population-SD uncertainty, and automatic equivalence claims** remain restricted. JMLE `Conditional` bounds hold other fitted parameters fixed as screening approximations; fixed-SD MML **Specified structural intervals** apply only to their own numerically checked targets and do not guarantee finite-sample coverage or simultaneous intervals.
- **Conditional bias diagnostics** hold main calibration fixed and do **not** establish statistical significance, absence of bias, or fairness. **Joint interaction** refits a separate fixed-`N(0,1)` unit-weight model; Holm adjusts its recorded family of requested tests, not its intervals.
- **G/D studies** require complete, balanced 1- or 2-random-facet designs with >= 2 levels per factor, 1 observation per cell, and unit weights. Multivariate G-theory and `G`/`Phi` confidence intervals are not implemented.
- **Blank `SE` or interval cells mean unavailable, not zero.** Input/design problems (such as disconnected rating groups) can be corrected in data setup, whereas method-scope limits cannot be removed by rerunning the same analysis. See the [interpretation guide](docs/analysis_scope.md).

---

## Guided Reporting, Manuscript Figures & Data Privacy

- **Saving & reopening:** On Desktop, **Save project** (`.mfrm`) preserves input data, settings, and results for instant reopening via **Data → Open project** without re-estimation *(older basic projects open cleanly, while older schema-1 apps cannot read extended projects saved by this preview—keep original files when updating)*. **Save report** writes a standalone offline `.html` report. See [Save and reopen](docs/desktop_user_guide.en.md#save-and-reopen). On Web, use **Report & Export** ([manuscript guide](docs/apa_manuscript_template.md)).
- **Manuscript figures:** Desktop **Manuscript** mode moves overall titles and long explanatory notes out of the image into an accompanying **`.caption.txt`** file while keeping axes, legends, labels, and withheld-result notices intact. On the hosted Web app, use the offered interactive HTML figure download whenever host PNG export is unavailable.
- **Data privacy:** Ratings and `.mfrm` projects can contain examinee/rater IDs and sensitive responses. Run confidential analyses locally on a controlled machine under your institution's rules. Uploaded Web data run on remote servers, and analysis IDs/fingerprints verify provenance rather than encrypting or anonymizing data.

---

<a id="faq"></a>

## Frequently Asked Questions (Q&A)

1. **Can I analyze confidential or IRB-restricted assessment data in MFRM Desktop?**
   **Yes.** MFRM Desktop runs **100% locally offline** on your Mac or Windows PC—CSV import, estimation, plotting, and `.mfrm` / `.html` / `.csv` saving never transmit ratings or IDs to external servers (only **Help → Web theory and help** opens online documentation in your browser, still without sending data). Because saved `.mfrm` projects and exported reports retain input ratings and examinee/rater identifiers for zero-re-estimation reopening, store them according to your institution's rules and de-identify records before sharing files externally or using the public Streamlit Web app.
2. **How many facets can I include in an analysis?**
   **Core MFRM calibration (`RSM` / `PCM`, plus `GPCM` on Web) has no fixed upper limit on the number of facets** (`Person` plus one or more non-person facets such as `Rater`, `Task`, `Criterion`, `Occasion`, etc., provided all levels are connected through shared ratings; Wright maps print direct labels for up to 4 non-person facets and show hover values beyond that). Supplementary modules have specific bounds: **Conditional bias** screens 2-facet pairs (`Facet A × Facet B`); fixed-SD MML **Specified intervals** and **Joint interaction** support up to `64` free structural coordinates within synchronous memory limits (with Joint interaction refitting one 2-facet pair); and **G / D reliability studies** support **1 or 2 random facets** (`p × r` or `p × t × r` complete balanced designs).
3. **Do I need Python or R installed, and can macOS and Windows collaborators share `.mfrm` files?**
   **No Python, R, or command-line setup is needed** for MFRM Desktop (Python is bundled inside both the macOS Apple Silicon and Windows 11 x64 preview ZIPs) or MFRM Web. Saved `.mfrm` project files reopen across macOS and Windows via **Data → Open project** without re-estimating *(keep original `.mfrm` backups when updating preview builds, as older apps cannot open newer extended projects)*. For scripted R workflows, use CRAN **[`mfrmr`](https://cran.r-project.org/package=mfrmr)**.

---

## Install

*(Skip this section if you are using the standalone [macOS](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-macos-20261030) or [Windows](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-windows-20261030) Desktop Preview ZIPs, which bundle Python.)*

To run the Web app from source, use **Python 3.11 or 3.12** in a virtual environment:

```bash
git clone https://github.com/Ryuya-dot-com/MFRM_Python_Application.git
cd MFRM_Python_Application
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, create the environment with `py -3.12 -m venv .venv` and invoke `.venv\Scripts\python.exe` in place of `python` for installation and startup commands (activation is optional).

## Run

From the repository root with your virtual environment:

```bash
python -m streamlit run streamlit_app.py
```

For Desktop source startup on macOS or Windows, see [CONTRIBUTING.md — Run from source](CONTRIBUTING.md#run-from-source). For hosted deployment, see [DEPLOYMENT.md](DEPLOYMENT.md).

---

## License & User Responsibility

Released under the **[MIT License](LICENSE)** ([License & Responsible Use Notice](LICENSE_NOTICE.md)).

- **Permitted uses (free of charge):** Academic research, thesis work, classroom teaching, institutional assessment, consulting, and **commercial use** are all permitted free of charge (retain copyright and license notices when redistributing; bundled third-party components keep their own licenses).
- **No warranty & user responsibility:** Because this software is provided **free of charge and "AS IS" without warranty of any kind**, **all responsibility rests with the user**—including data privacy and institutional compliance, verifying model assumptions and design eligibility, respecting withheld-inference boundaries, and any substantive, educational, clinical, or operational decisions made from the outputs. Neither the author (**Ryuya Komuro**) nor copyright holders assume liability for any claim, loss, or damages arising from its use.

<a id="citing-this-software"></a>

## Citing this software (APA 7th edition)

Cite **Ryuya Komuro** in APA 7th edition format (`(Komuro, 2026)` or `Komuro (2026)` in text) for the edition and version you actually used, alongside the underlying statistical methods:

- **MFRM Desktop (macOS / Windows):**
  > Komuro, R. (2026). *MFRM Desktop* (Version 0.2.19-beta) [Computer software]. https://github.com/Ryuya-dot-com/MFRM_Python_Application

  *(Also record the OS release tag—`desktop-macos-20261030` or `desktop-windows-20261030`—and archive filename in your Methods section; built-in [`CITATION.cff`](CITATION.cff) sidecars record the shared repository title `MFRM Streamlit`.)*
- **MFRM Web (`MFRM Streamlit`):**
  > Komuro, R. (2026). *MFRM Streamlit* (Version 0.2.19-beta) [Computer software]. https://github.com/Ryuya-dot-com/MFRM_Python_Application

  *(BibTeX and `CITATION.cff` downloads are available under **About MFRM** in the Web app.)*
- **R package (`mfrmr` on CRAN):**
  > Komuro, R. (2026). *mfrmr: Estimation and diagnostics for many-facet measurement models* (Version 0.2.3.1) [R package]. https://CRAN.R-project.org/package=mfrmr

  *(Replace `0.2.3.1` with your installed version via `packageVersion("mfrmr")` or `citation("mfrmr")` in R.)*

Citation is appreciated and is not an additional condition of the MIT License.

## Bug reports & contributions

If you encounter a bug, crash, or unexpected behavior, please report it by email to the author, **Ryuya Komuro**, at **[ryuya.komuro.c4@tohoku.ac.jp](mailto:ryuya.komuro.c4@tohoku.ac.jp)** with your OS and version, edition/app version (`0.2.19-beta`) and archive filename (or `mfrmr` version), steps to reproduce, and the exact error message or screenshot. **Never include confidential or identifiable examinee/rater data in emails or attachments**—reproduce with the built-in sample or anonymized synthetic data whenever possible. Contribution and build instructions are in [CONTRIBUTING.md](CONTRIBUTING.md).
