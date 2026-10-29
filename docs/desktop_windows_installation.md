# MFRM Desktop for Windows — validation candidate

This is a candidate for Windows 11 x64 acceptance, not an accepted public release.
The current in-app version is 0.2.19-beta; the build manifest and archive SHA-256
identify the exact candidate. Windows on ARM is not an accepted target.

## Open and try the application

1. Download the application ZIP and its SHA-256 from the supplied candidate link.
   GitHub's automatically generated source ZIP is not the application.
2. Use File Explorer's **Extract All**. Keep `MFRM-Desktop.exe` together with the
   complete `_internal` folder. Do not run the exe from inside the ZIP or copy it alone.
3. Open `MFRM-Desktop.exe`. Python, R and developer tools are not required.
   This candidate is unsigned. Record the exact Windows security message if
   opening is blocked; do not disable Windows protection to complete a test.
4. Select **Data → Sample → Next → Analyze**. Check both **Figure** and **Table**.
5. Save a project and an HTML report in a writable folder. Quit and restart MFRM,
   then select **Data → Open project** and reopen the saved `.mfrm` file.
   Open the HTML in a browser and check the figures as well as the text and tables.

The portable ZIP does not register a default application for `.mfrm` files or
change the Windows registry. Use **Open project**, or choose the extracted exe
using Windows **Open with**. Open one project at a time.

## Included analyses and outputs

The interface starts in English and can switch to Japanese. RSM and PCM support
JMLE and fixed-population-SD MML. Outputs include Wright and Pathway maps,
category curves, conditional bias diagnostics, eligible specified intervals,
and the supported joint interaction analysis. G/D studies support complete,
balanced designs with one or two random facets. Multivariate G-theory and G/D
uncertainty intervals are not included.

Read each result's availability and assumptions. Conditional bias diagnostics
do not establish significance or fairness. Joint interaction inference has its
own approximate pointwise intervals and multiplicity rules.

Projects include the input data and identifiers. HTML reports contain local
figures; viewing them does not require a plotting service. CSV, PNG and SVG can
also be saved. Manuscript figures have a separate `.caption.txt` file: retain
its explanation when preparing an external figure note.

## Update and acceptance

Keep an untouched copy of the earlier app and projects. Extract a newer version
into a separate folder and open a copy of a project before replacing anything.
Current projects use schema 2; older schema-1-only apps cannot read schema 2.

For this candidate, record the Windows edition/version, CPU type, app archive
hash, operation and exact error. Real-PC checks include Japanese paths, display
scaling, cancellation, output files, restart/restore and downloaded-app security
messages. Successful Windows Server build checks do not establish Windows 11
GUI acceptance. No automatic update or account is required.
