# MFRM — Many-Facet Rasch analysis

[日本語](README.ja.md) · [Open the Web app](https://mfrmpythonapplication-xodghmcghngx8cz5zmhfeg.streamlit.app/) · [Windows guide](docs/desktop_user_guide.en.md)

Analyze ratings from writing, speaking, practical examinations and other
assessments. MFRM helps you examine person measures, rater severity, task
difficulty and rating-category use, then save figures and reports.

Start with the built-in synthetic writing example before preparing your own
data. No FACETS, ConQuest or R installation is required.

## Choose how to use MFRM

| Use | Start here |
|---|---|
| Try the application in a browser | [Web app](https://mfrmpythonapplication-xodghmcghngx8cz5zmhfeg.streamlit.app/). Use synthetic or appropriately de-identified data on the hosted service. |
| Use the Windows desktop preview | [Download the Windows preview](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-windows-20261030) and follow the [Windows guide](docs/desktop_user_guide.en.md). Python is included; extract the whole folder before opening the EXE. |
| Run the Web application on your computer | [Install](#install) and [run](#run) from source. |
| Develop or build the application | [Contributing](CONTRIBUTING.md), including desktop source startup and verification. |

**Release availability:** a **Windows 11 x64 desktop preview** is available
from the release link above, with the application ZIP hosted on Dropbox.
It includes the complete desktop workflow described below, uses app version
**0.2.19-beta**, and is unsigned. See the release for the archive name, SHA-256,
tested environment and known limits. Windows on ARM and other Windows versions
have not been verified. The macOS application is not yet publicly released.
The separate **v0.2.19-beta** Web release provides the Web application and Python
source; it is not the Windows download.

GitHub's **Source code (zip)** is source code, not a desktop installer.
The Web app's **About MFRM** identifies its running version and source.
Features in a source branch may differ from the published release.

## Model scope

| Task | Available workflow |
|---|---|
| Estimate person measures and facet effects | RSM or PCM, with JMLE or fixed-population-SD MML |
| Examine the rating scale and fit | Wright maps, Pathway maps, category curves, estimates and fit tables |
| Investigate combinations of facet levels | Conditional bias diagnostics; the desktop preview also includes a separately fitted joint-interaction model |
| Examine selected uncertainty | Specified pointwise intervals for eligible fixed-SD MML facet levels and prespecified differences |
| Compare numbers of raters or tasks | G/D studies for supported complete, balanced designs |
| Save and report | Web reports and downloads; desktop projects, HTML, CSV, PNG and SVG |

The desktop preview brings these workflows together in one application,
with English/Japanese controls and explained/manuscript figure modes.
The Web interface has additional advanced workflows, including bounded GPCM,
latent regression, anchors/linking, residual PCA and simulation.
See [analyses and interpretation limits](docs/analysis_scope.md) before using
an output to support a statistical conclusion.

## Try the sample

**Web:** open the app, choose **Continue without the guide** or the sample
learning route, and use **Built-in example → Writing essay**. Run the analysis,
read **First Read**, then choose a relevant diagnostic and **Report & Export**.

**Windows preview:** choose **Data → Sample → Next → Analyze**. Inspect
**Figure** and **Table**, then use **Save project** and **Save report**.
The [Windows guide](docs/desktop_user_guide.en.md) explains each step.

## Data format

Use long-format ratings: one row per person, rater, task and criterion rating.

    Person,Rater,Task,Criterion,Score
    P01,R1,T1,Content,4
    P01,R1,T1,Organization,3
    P01,R2,T1,Content,5
    P02,R1,T1,Content,2

This is a format example, not a complete analysis dataset. A full
[synthetic writing CSV](examples/desktop_writing.csv) is included.

Map the person, score and facet columns in the application. For MFRM, enter
integer scores and the minimum and maximum allowed by the scoring rubric,
including categories absent from the observed data. Do not record an unscored
response as zero unless zero is its actual score.

The desktop entry reads CSV. The Web entry also supports spreadsheet paste,
TSV/TXT, Excel, Parquet and JSON formats. See the Windows guide for
[data preparation](docs/desktop_user_guide.en.md#use-your-own-csv);
use the Web app's input preview to check imported rows and column roles.

## Statistical caveats

A completed calculation does not by itself establish that its results support
inference. Read the message beside the specific output.

- General MML tables, Person-score uncertainty, estimated population-SD
  uncertainty and automatic equivalence conclusions retain restrictions.
  Eligible **specified structural intervals** have a separate, limited scope.
- Conditional bias diagnostics do not confirm the presence or absence of bias.
  Joint-interaction inference uses a separate model and explicit assumptions.
- G/D supports complete balanced designs with one or two random facets.
  Multivariate G-theory and G/Phi uncertainty intervals are not implemented.
- Blank SEs or intervals are unavailable, not zero. Some limitations concern
  the implemented method and cannot be removed by rerunning the same analysis.

The [interpretation guide](docs/analysis_scope.md) explains these distinctions,
the assumptions behind approximate intervals, and differences from other
software. It also distinguishes correctable input problems from unsupported
uses of a result.

## Guided reporting workflow

Save results you need outside the current session.

- **Desktop:** save a project to reopen the data, settings and results without
  re-estimation. Save an HTML report for reading and sharing, or export tables
  and figures separately. Follow the [saving guide](docs/desktop_user_guide.en.md#save-and-reopen).
- **Web:** use **Report & Export → Create a report**, **Save tables & figures**
  or **Start a paper** for the output you need. The
  [manuscript guide](docs/apa_manuscript_template.md) includes a worked writing
  example and reporting guidance in Japanese.

### Figure exports

Desktop manuscript mode moves the overall title and explanatory notes outside
the graphic while retaining axes, legends and data labels. Keep the accompanying
caption file when using the image in a paper.

Static image export in the Web app depends on its host. If PNG export is
unavailable, use the offered HTML figure download. The desktop preview saves
PNG/SVG locally. Read each output's notes when moving a figure or table into
another document.

## Data privacy

Ratings and saved projects can include person/rater IDs and other sensitive
data. For confidential work, use a controlled local environment and follow
your institution's data-handling requirements.

Uploaded data in the hosted Web app are processed on remote infrastructure.
Desktop analysis runs on your computer; opening the linked Web Help uses the
internet. Review project files, reports, tables and screenshots before sharing
them. Analysis IDs and reproducibility fingerprints are not encryption or
anonymization. A browser session is not a backup.

## Install

This section is for running the Web app from source. A supplied desktop
application ZIP does not require these steps.

Use Python 3.11 or 3.12 and a virtual environment:

    git clone https://github.com/Ryuya-dot-com/MFRM_Python_Application.git
    cd MFRM_Python_Application
    python3 -m venv .venv
    source .venv/bin/activate
    python -m pip install -r requirements.txt

On Windows, create the environment with **py -3.12 -m venv .venv** and use
**.venv\Scripts\python.exe** in place of **python** for the installation and
run commands; activation is optional.

## Run

From the source folder with the chosen environment:

    python -m streamlit run streamlit_app.py

For desktop source startup, use the
[developer instructions](CONTRIBUTING.md#run-from-source).
For hosted deployment, see [DEPLOYMENT](DEPLOYMENT.md).

## License

[MIT License](LICENSE) · [License notice](LICENSE_NOTICE.md).
Research, teaching and commercial use are permitted under the license.
Retain the required copyright and license notices when redistributing.
The software is provided without warranty; users remain responsible for
data handling, model assumptions and the interpretation and use of results.
Third-party components retain their own licenses.

## Citing this software

Cite the version actually used and record its source commit. GitHub's citation
menu uses [CITATION.cff](CITATION.cff); the Web app also provides citation
downloads under **About MFRM**.

For version 0.2.19-beta:

> Komuro, R. (2026). *MFRM Streamlit* (Version 0.2.19-beta) [Computer software]. https://github.com/Ryuya-dot-com/MFRM_Python_Application

Report the model, estimator, constraints and uncertainty method separately.
The software citation does not replace citations for statistical methods.
Do not substitute the current version for an older saved analysis. Citation is
appreciated, not an additional MIT license condition.

## Help and contributions

For questions or bugs, use
[GitHub Issues](https://github.com/Ryuya-dot-com/MFRM_Python_Application/issues).
Include the OS, version, steps and exact message. Remove identifying data
before attaching files or screenshots.

Development setup, verification and contribution guidance are in
[CONTRIBUTING](CONTRIBUTING.md).
