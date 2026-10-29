MFRM Desktop — Windows x64 validation candidate

Extract the complete ZIP, then open MFRM-Desktop/MFRM-Desktop.exe.
Keep the EXE and the entire _internal folder together. Python, R, FACETS and
ConQuest do not need to be installed separately.

This candidate rebuilds GitHub main d78fa706 with the corrected handling of
missing score rows. It is not the earlier 49eec16 integration candidate.
The current main desktop has a Japanese interface and supports RSM/JMLE,
Wright/Pathway/category-probability plots, HTML and PNG/SVG output, and .mfrm
project save/restore. Other models and analyses in the web application or
separate development branches are not included in this executable.

Start with Data (データ) → Sample (サンプル) → Next (次へ) → Analyze (分析).
Save project (プロジェクトを保存) stores input and results in a .mfrm file.
Restart the app and use Open project (プロジェクトを開く) inside the app.
Save report (レポートを保存) creates an HTML file that opens in a browser.

The candidate is unsigned. Record the text of any Windows startup warning.
Windows-runner checks cover the missing-score calculation, worker cancel and
restart, three report figures, and a separate-process restore without refitting.
This newly rebuilt executable has not yet completed Windows 11 GUI/download
acceptance. These checks do not establish statistical coverage guarantees.
Read the inference holds and interpretation limits shown by the application.

Extract each version into a separate folder and keep previous versions and
projects. Open a copy of an existing project when checking a replacement.
There is no installer, automatic update, or file association setup.
The bundled third_party_notices directory records dependency licenses.

日本語

ZIPを「すべて展開」して、MFRM-Desktop.exeを開いてください。
EXEと_internalは一緒に保持します。Python等の追加インストールは不要です。
GitHub mainの欠損得点行の修正を含む再ビルド版です。
現在のmainの画面・機能（日本語、RSM/JMLE）に対応します。
別ブランチのPCM・MML・G/D等を含む版とは異なります。
新しい実行物のWindows 11実機受入は未完了の検証用候補です。
既存プロジェクトのコピーを使い、旧版は別フォルダーに保持してください。
