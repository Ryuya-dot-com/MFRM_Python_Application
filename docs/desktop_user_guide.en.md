# MFRM Desktop guide for Windows

[日本語](desktop_user_guide.ja.md) · [Introduction](../README.md) · [Get and open the app](#get-and-open-the-app) · [Try the sample](#try-the-sample) · [Use your own CSV](#use-your-own-csv) · [Save and reopen](#save-and-reopen)

MFRM Desktop analyzes ratings from writing, speaking, practical examinations
and similar assessments. This guide takes you from opening the application to
viewing and saving results.

The desktop preview includes RSM/PCM, JMLE/fixed-population-SD MML, bias
diagnostics, specified intervals, joint interactions and G/D studies in one
application. FACETS, ConQuest and R are not required. The supported designs and
interpretation limits differ by analysis; see [Choose an analysis](#choose-an-analysis).

The interface starts in English. The language selector switches to **日本語**
without translating your column names or level labels or changing calculations.

## Get and open the app

Download the application ZIP from the
[Windows desktop preview release](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-windows-20261030).
The release links to the Dropbox download and gives its SHA-256. Python is
included. GitHub's automatically generated source ZIP and the separate
v0.2.19-beta Web release are not the Windows application download.

The Windows preview targets **Windows 11 x64**. Windows on ARM has not been
verified. The release states the tested environment and remaining limitations;
this is a preview rather than a stable release for every Windows configuration.

To open the application ZIP:

1. Obtain the ZIP and its matching SHA-256 from the supplied location.
2. Use File Explorer's **Extract All**. Keep `MFRM-Desktop.exe` and the complete
   `_internal` folder together.
3. Open the EXE. Do not run it inside the ZIP or move the EXE alone.
4. Try the sample below before preparing your own data.

The preview is unsigned. If Windows blocks opening it, record the exact
message and include it in a support request. Do not disable Windows protection
to run it.

If you do not have an application ZIP, try the
[Web app](https://mfrmpythonapplication-xodghmcghngx8cz5zmhfeg.streamlit.app/)
without installing Python. For source startup only, use the
[developer instructions](../CONTRIBUTING.md#run-from-source); available features
depend on the source branch or version.

Normal desktop analysis, figures and saving run on your computer and require
no account. **Help → Web theory and help** needs an internet connection;
opening that link does not upload your desktop data or results.

## Try the sample

1. In **Data**, choose **Sample** to load synthetic writing ratings.
2. Choose **Next**.
3. In **Settings**, check the sample's automatically selected values:

| Setting | Sample value |
|---|---|
| Person | `Person` |
| Score | `Score` |
| Raters, tasks and other facets | `Rater`, `Task`, `Criterion` |
| Minimum score / Maximum score | `0` / `5` |

4. Choose **Analyze**. The app moves to **Results** when calculation finishes.
5. Use **Figure** to select Wright map, Pathway map or category probabilities.
   **Reading guide** explains the selected figure.
6. Use **Table** to inspect estimates, category counts and data summaries.
   Read the accompanying notes.
7. Choose **Save project** and save a practice file.

In Pathway map, the selector beside the figure changes the displayed facet,
such as rater, task or criterion.

## Use your own CSV

### Arrange the ratings

Use **one row for one person's rating by one rater on one task and criterion**.
Repeated ratings of the same person belong on separate rows.

```csv
Person,Rater,Task,Criterion,Score
E01,R1,Task1,Content,4
E01,R1,Task1,Organization,3
E01,R2,Task1,Content,5
E02,R1,Task1,Content,2
```

This is a format excerpt, not a complete analysis dataset. Use the full
[writing example CSV](../examples/desktop_writing.csv) for practice, or choose
**CSV format → Save example CSV** in the application.

| Column | Contents |
|---|---|
| `Person` | A consistent ID for each examinee or learner. |
| `Rater` | Rater ID. |
| `Task` | Task ID, if tasks are part of the design. |
| `Criterion` | Rating criterion, such as content or organization. |
| `Score` | The integer score for that rating, not an average or total. |

Japanese column names are supported. Include the facets present in your design.
When saving from Excel, use **CSV UTF-8 (Comma delimited)**. Do not turn an
unscored response into zero: omit that rating or leave its score blank. Zero
means an actual rating of zero.

If each rater has a separate score column, reshape the data into rater/score
rows as above. Check duplicate person/rater/task/criterion combinations. If
ratings are repeated on separate occasions, represent that distinction in the
data and choose an analysis appropriate to the resulting design.

### Import and map columns

1. Choose **Open CSV** and inspect the preview of the first 10 rows.
2. Choose **Next**, then select the **Person** and **Score** columns.
3. Select one or more columns in **Raters, tasks and other facets**. Their roles
   must not overlap the person or score column.
4. Enter the scoring rubric's minimum and maximum. For a 0–5 rubric, enter `0`
   even if no observed rating is zero.
5. Choose **Analyze**. Review any excluded rows or duplicate combinations in
   the confirmation message before continuing.

Blank or nonnumeric scores are excluded after review. Fractional, infinite or
out-of-range MFRM scores require correction. This desktop version treats each
rating equally; adding a weight column does not enable weighted analysis.

## Choose an analysis

In **Settings**, select **MFRM** or **G/D Study**. MFRM settings use the
**Columns / rubric**, **Model / estimation**, **Bias / intervals** and
**G/D plan** tabs.

| Purpose | Selection and limits |
|---|---|
| Share category steps across ratings | **RSM**, with JMLE or fixed-SD MML. |
| Allow steps to differ across levels of a facet | **PCM**, then select the **PCM step facet**. Both estimators are available. |
| Estimate marginally with a fixed population SD | **MML** uses EM. Set quadrature points and the fixed population SD; the app does not estimate that SD. |
| Investigate combinations of facet levels | In **Bias / intervals**, register a facet pair with **Add**. These conditional diagnostics hold calibration fixed and do not establish significance or presence/absence of bias. |
| Obtain intervals for a level or planned difference | Register **Specified intervals** with fixed-SD MML. Only numerically eligible targets receive pointwise approximate intervals; other targets retain a reason. |
| Fit an interaction together with main effects and steps | Select one pair and its cells/contrast in **Joint interaction**, with **MML**, fixed SD `1`, and **Add to next analysis** enabled. This refits a separate RSM/PCM model with a fixed N(0,1) person population. Eligible outputs include approximate Wald intervals/p values and Holm-adjusted p values. Holm does not adjust the intervals. |
| Examine reliability under alternative rater/task counts | In **G/D Study**, select the object, observed score and one or two random facets, and planned counts. G/D can also run alongside MFRM. |

For bias and specified intervals, entering values alone does not register a
request: choose **Add** and check the resulting list. Joint interaction keeps
the original main-effects model and fixed-calibration diagnostic separate.

G/D uses original observed scores, without an MFRM fit or integer rubric
requirement in its independent entry. It requires fully crossed, balanced
data, at least two levels of each factor, one observation per cell and unit
weights. G describes relative decisions and Phi absolute decisions. Target
populations are your declarations. Negative variance components remain visible
as signed values and are set to zero for coefficient calculation. Multivariate
G-theory and G/Phi intervals are not implemented.

Specified intervals and joint interactions assume independent persons, the
specified response model and normal population, ignorable assignment/missingness
and an identified interior solution, among other conditions. Numerical checks
do not guarantee finite-sample coverage. See the
[interpretation guide](analysis_scope.md) for these assumptions and SE definitions.

## Read the results

Check the analyzed data and any messages before interpreting a figure.

| Output | What to examine |
|---|---|
| Wright map | The distribution of person measures and facet locations. Under the usual sign convention, higher raters are more severe and higher tasks more difficult. |
| Pathway map | Measures vertically and fit horizontally; the initial fit measure is Infit MnSq. Interpret deviations from 1 alongside rating counts and the original ratings. |
| Category probabilities | Where each score category is likely; compare with category frequencies and the scoring rubric. |
| Estimates and fit tables | Person and facet results, together with the output's notes. |
| Category counts | How often each score category was observed. |

General MML tables and person-score uncertainty retain inferential limits.
Dedicated specified-interval and joint-interaction outputs have their own
eligibility conditions. A general warning does not mean every dedicated output
is unavailable; a drawn figure does not mean inference is valid either.
Read the message for the result you intend to use.

## Save and reopen

| Purpose | Action and format |
|---|---|
| Continue the analysis later | **Save project** → `.mfrm`: input data, settings and results. |
| Read or share figures and tables together | **Save report** → `.html`, then open it in a browser. Figures are included locally. |
| Use a figure in a document | Select the figure, then **Save figure…** → `.png` or `.svg`. |
| Use a table in another program | Select the table, then **Save CSV…**. Retain the output's notes. |

**Figure presentation** switches between **Explanatory** and **Manuscript**.
Manuscript mode moves overall titles and long explanatory notes outside the
graphic while retaining axes, legends and labels. **Title and note** provides
copyable text. A manuscript PNG/SVG also saves a `.caption.txt` file: use its
content in the paper's figure note. This option does not by itself establish
compliance with every journal or APA requirement.

To reopen without refitting, start the app and choose **Data → Open project**,
then select the `.mfrm` file. Images and HTML cannot restore an analysis.
The portable Windows app does not register a default application for `.mfrm`;
use the in-app button if Windows asks which program should open it.

This preview reads older basic-analysis projects, but older applications cannot
read the extended projects saved by this version. Keep the old app and original
projects. Extract an update into a separate folder and open a copy first.
Record the supplied ZIP name as well as the displayed version: preview builds
can share a version label.

Projects contain input ratings and person/rater identifiers. Review their
contents before sharing. Choose a writable folder for your saved work.

Changing settings does not recalculate displayed results. If they are marked
as belonging to the previous analysis, preserve any result you need, check the
new settings and choose **Analyze**.

## Troubleshooting

| Situation | Action |
|---|---|
| The EXE will not open | Confirm full extraction, the adjacent `_internal` folder and an x64 Windows system. Record the exact security message if blocked. |
| CSV import fails or text is garbled | Save as UTF-8 CSV; check for blank or duplicate column names. |
| Analyze leaves you in Settings | Follow the validation message for column roles, rubric bounds or score values. |
| Fewer rows are analyzed than expected | Review the excluded rows and reasons, including missing scores or required IDs. |
| An interval is unavailable because of data/design or numerical conditions | Check the particular reason: column selection, connections between rating groups, extreme scores, termination or other stated conditions. A plotted point alone does not establish a difference or ranking. |
| Inference is unavailable because the method/use is unsupported | General MML tables, person-score uncertainty and G/D intervals have restrictions that repeated runs or data edits do not resolve. Check the scope of any available dedicated output. |
| An SE or interval is blank | Read the notes; blank does not mean zero uncertainty. |
| A figure has not changed after editing settings | It belongs to the previous analysis. Preserve it if needed, then run the new settings. |
| You need to stop calculation | Choose **Cancel** while the analysis is running. |

For unresolved problems, use
[GitHub Issues](https://github.com/Ryuya-dot-com/MFRM_Python_Application/issues).
Include the OS, app version and ZIP name, steps, and exact message. Remove
identifying data from files and screenshots. Additional theory is available
through **Help → Web theory and help**; the public Web version may differ.
