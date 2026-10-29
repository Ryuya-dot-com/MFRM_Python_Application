# MFRM Desktop for Windows — getting started

This desktop preview is for Windows 11 x64. It is not yet a general desktop
release. Windows on ARM has not been verified. Python is included; installing
Python, R, FACETS, ConQuest or developer tools is not required.

## Open and try the application

1. Obtain the application ZIP and matching SHA-256 from the supplied location.
   GitHub's automatically generated source ZIP is not the application.
2. Use File Explorer's **Extract All**. Keep `MFRM-Desktop.exe` together with the
   complete `_internal` folder. Do not run the EXE inside the ZIP or copy it alone.
3. Open `MFRM-Desktop.exe`. The preview is unsigned. If Windows blocks it,
   record the exact message for support; do not disable Windows protection.
4. Select **Data → Sample → Next → Analyze** to try synthetic writing ratings.
   The sample sets the scoring range to 0–5. Inspect **Figure** and **Table**;
   use **Reading guide** and table notes to interpret the outputs.
5. Use **Save project** to save a `.mfrm` file in a writable folder. Use
   **Save report** for an HTML report that opens in a browser with local figures.
6. To return to the saved analysis, start MFRM and select **Data → Open project**.
   Results reopen without refitting. An HTML report or image cannot restore
   an analysis, so keep the project if you intend to continue working.

The app starts in English; the language selector also offers **日本語**.
The portable ZIP does not register a default application for `.mfrm` files.
Use **Open project** if Windows asks which program should open one. Open one
project at a time.

## Prepare your data and choose an analysis

Use **Open CSV** to import your data. **CSV format → Save example CSV** provides
a complete example. Each row is one rating, with person, score and the facet
IDs relevant to your design. Review the preview and choose each column's role.
For MFRM, use integer scores and enter the scoring rubric's full minimum and
maximum, even if some categories were not observed. Missing ratings are not
zero scores. Each rating has equal weight in this desktop version.

RSM and PCM support JMLE and fixed-population-SD MML. Outputs include Wright
and Pathway maps, category probabilities, conditional bias diagnostics,
specified intervals and joint interaction analysis. Register bias/interval
requests with **Add** before analysis. Joint interaction is a separate jointly
refitted model and requires MML with fixed SD 1.

G/D Study uses observed scores and supports fully crossed, balanced designs
with one or two random facets, at least two levels of each factor and one
observation per cell. Multivariate G-theory and G/Phi uncertainty intervals
are not implemented.

Read the message for each result. Conditional bias diagnostics do not establish
significance or fairness. General MML tables and person-score uncertainty retain
inferential limits. Dedicated intervals and joint interactions have their own
conditions and pointwise approximations; numerical checks do not guarantee
finite-sample coverage. Holm adjustment in joint interaction applies to the
recorded tests, not the intervals. Blank uncertainty values are unavailable,
not zero. Unsupported methods cannot be enabled by repeatedly rerunning an
analysis. **Settings help**, figure guides and table notes give further details.

## Save figures, update and share

Tables save as CSV, and figures as PNG or SVG. **Manuscript** figure mode moves
titles and long notes outside the image. Keep the accompanying `.caption.txt`
file and use it when writing a figure note. HTML reports retain figures and
notes for viewing without a plotting service.

Projects include input data and identifiers. Review projects, reports and
screenshots before sharing. Normal analysis and saving run locally without an
account or internet access. **Help → Web theory and help** opens an online page;
it does not upload desktop data or results.

For updates, keep an untouched copy of the previous app and projects. Extract
the new app into a separate folder and try a copy of a project first. This
preview reads older basic-analysis projects, but its extended saved projects
cannot be read by older applications. There is no automatic update.

For support, include the OS, app version, supplied ZIP name, steps and exact
message. Preview builds may share the in-app version 0.2.19-beta, so keep the
archive name and supplied SHA-256. Remove identifying data before attaching
files or screenshots.

Online project information and release announcements:
https://github.com/Ryuya-dot-com/MFRM_Python_Application

Questions and bug reports:
https://github.com/Ryuya-dot-com/MFRM_Python_Application/issues
