"""Desktop UI translations; user column names and stored calculation settings stay intact."""
import re

EN = {
'プロジェクトは一度に1つ開いてください。': 'Open one project at a time.',
'分析の完了または中止後にプロジェクトを開いてください。': 'Wait for the analysis to finish or cancel it before opening a project.',
'データ':'Data','設定':'Settings','結果':'Results','ヘルプ':'Help','CSVを開く':'Open CSV','サンプル':'Sample','G/D 2相サンプル':'G/D 2-facet sample','プロジェクトを開く':'Open project','CSV形式':'CSV format','次へ':'Next','設定の説明':'Settings help','列・採点基準':'Columns / rubric','モデル・推定':'Model / estimation','バイアス・区間':'Bias / intervals','G/D-Study計画':'G/D-Study plan','列':'Columns','対象者':'Person','得点':'Score','評価者・課題など（相）':'Raters, tasks and other facets','採点基準':'Scoring rubric','最低点':'Minimum score','最高点':'Maximum score','戻る':'Back','分析':'Analyze','中止':'Cancel','プロジェクトを保存':'Save project','レポートを保存':'Save report','分析情報':'Analysis information','モデル':'Model','推定法':'Estimation method','PCMの段階相':'PCM step facet','MMLの求積点数':'MML quadrature points','MMLの固定母集団SD':'MML fixed population SD','最大反復数':'Maximum iterations','停止許容誤差':'Stopping tolerance','バイアス（一覧に登録した条件で計算）':'Bias (registered requests)','相A':'Facet A','相B':'Facet B','探索範囲の絶対値':'Absolute search bound','極端な水準を除外':'Omit extreme levels','追加':'Add','選択を更新':'Update selected','選択を削除':'Remove selected','有意性・バイアスの有無の結論は保留。':'Significance and presence/absence of bias remain withheld.','指定区間：固定SD MML（一覧に登録した条件で計算）':'Specified intervals: fixed-SD MML (registered requests)','相':'Facet','信頼水準':'Confidence level','対象':'Targets','全水準':'All levels','A − B（事前指定）':'A − B (prespecified)','水準A':'Level A','水準B':'Level B','保留対象は区間を描かず理由を保持します。':'Withheld intervals are not drawn; reasons remain available.','共同交互作用：主効果・段階も同時に再推定する別モデル':'Joint interaction: a separate model refitting main effects and steps','次回の分析に追加（MML・固定SD 1のみ）':'Add to next analysis (MML, fixed SD 1 only)','全セル':'All cells','全セル＋指定差':'All cells + specified contrast','差の左側':'Left cell','差の右側':'Right cell','：A水準':': A level','：B水準':': B level','MFRMと併せてG/D-Studyも計算（G/D-Study入口では常に実行）':'Add G/D-Study to MFRM (always computed in G/D-Study)','測定対象の列':'Object of measurement','観測得点の列':'Observed score column','ランダム相（1〜2列）':'Random facets (1–2 columns)','順序入替':'Swap order','固定相の扱い（3相以上のデータ）':'Fixed-facet handling (3+ facets)','指定しない（選択した相のみで均衡交差）':'None (selected facets already crossed)','1水準を抽出':'Condition on one level','全水準の合計得点（完全均衡を確認）':'Sum across all levels (verified balanced)','全水準の平均得点（完全均衡を確認）':'Mean across all levels (verified balanced)','固定相の列':'Fixed-facet column','抽出する水準':'Selected level','D-Study計画の指定':'D-Study plan mode','観測水準数を使う':'Use observed counts','計画数を指定':'Specify custom counts','計画表の上限':'Maximum planned count','対象母集団（宣言）':'Target population (declared)','横軸を入れ替え':'Swap X-axis','入力データにある水準を選択してください。':'Choose levels present in the input data.','更新する条件を一覧から選択してください。':'Select a registered request to update.','G/D-Studyの測定対象・得点・ランダム相1〜2列を、重複なく選択してください。':'Select distinct object, score and one or two random-facet columns for G/D-Study.','固定相の列と水準を選択してください。':'Select a fixed-facet column and level.','採点基準の最低点':'Rubric minimum','採点基準の最高点':'Rubric maximum','反復数・許容誤差には正の値を入力してください。':'Iterations and tolerance must be positive.','PCMの段階相を選択してください。':'Select a PCM step facet.','求積点数':'Quadrature points','固定母集団SD':'Fixed population SD','MMLの求積点数は3以上、固定SDは正の有限値にしてください。':'MML requires at least 3 quadrature points and a positive finite fixed SD.','共同交互作用にはMML・固定母集団SD 1を指定してください。':'Joint interaction requires MML with fixed population SD 1.','指定差には相A・Bの観測された水準を選択してください。':'Choose observed levels of facets A and B for the contrast.','共同交互作用の信頼水準':'Joint-interaction confidence level','設定は変更されています。表示中の結果は前回の分析です。':'Settings have changed. The displayed results are from the previous analysis.','データ内の得点：':'Observed scores: ','数値なし':'No numeric scores','\n\n分析ID: ':'\n\nAnalysis ID: ','\n\nプロジェクトファイルには入力データと識別子を含みます。':'\n\nProject files include input data and identifiers.','\n\n計算版・出典: ':'\n\nCalculation version / source: ','未保存の結果':'Unsaved results','現在の結果を保存せずに進みますか？':'Continue without saving the current results?','処理できませんでした。入力や保存先を確認してください。':'Operation failed. Check the input or destination.','処理できませんでした':'Operation failed','指定した場所へのアクセスが拒否されました。\n読込・保存先の権限を確認し、利用できる場所を選んでください。\n\n詳細: ':'Access to this location was denied.\nCheck permissions and choose an accessible location.\n\nDetails: ','閉じる':'Close','CSVの形式と見本':'CSV format and example','見本CSVを保存':'Save example CSV','作文_見本.csv':'writing_example.csv','見本は .csv 拡張子で保存してください。':'Save the example with a .csv extension.','使い方・理論':'Help and theory','Web版の理論・ヘルプ':'Web theory and help','ブラウザーを開けませんでした。ヘルプ欄のURLをコピーして開いてください。':'Could not open the browser. Copy the URL from Help and open it manually.','作文サンプル（人工データ）':'Writing sample (synthetic data)','G/D 2相サンプル（人工データ）':'G/D 2-facet sample (synthetic data)','先頭10行':'First 10 rows','入力データの確認':'Review input data','中止しました。':'Analysis cancelled.','推定しています…':'Estimating…','診断と利用条件を確認しています…':'Checking diagnostics and eligibility…','表を準備しています…':'Preparing tables…','バイアスを計算しています…':'Computing conditional bias…','指定区間の数値条件を確認しています…':'Checking specified intervals…','G/D-Studyを計算しています…':'Computing G/D-Study…','主効果・段階・交互作用を共同推定し、区間を確認しています…':'Jointly estimating main effects, steps and interaction; checking intervals…','計算に失敗しました。入力と設定を確認してやり直してください。':'Analysis failed. Check the input and settings before retrying.','図':'Figure','すべて':'All','近似Infit ZSTD':'Approx. Infit ZSTD','図を保存…':'Save figure…','読み方':'Reading guide','表':'Table','CSVを保存…':'Save CSV…','一般表の推論は保留です：':'General-table inference is withheld: ',' 保存時から利用条件が変わったため、現在の判定を適用しています。':' Eligibility has changed since saving; current rules apply.',' 固定較正のバイアス診断では有無の推論は保留です。':' Fixed-calibration bias diagnostics do not support presence/absence conclusions.',' 共同交互作用は別モデルの点ごとの近似推論です。利用可否は専用の表・図で確認してください。':' Joint interaction uses a separate model with pointwise approximate inference. Check its dedicated tables and figures for eligibility.',' G/D-Studyの区間推定は保留です。':' G/D-Study interval estimates are withheld.','.csv 拡張子で保存してください。':'Save with a .csv extension.','MFRMプロジェクト':'MFRM project','プロジェクトは .mfrm 拡張子で保存してください。':'Save projects with a .mfrm extension.','レポートは .html 拡張子で保存してください。':'Save reports with an .html extension.','説明付き':'Explanatory','論文掲載用':'Manuscript','図の形式':'Figure presentation','タイトル・図注':'Title and note','カテゴリー確率':'Category probabilities','バイアス':'Bias','指定区間':'Specified intervals','共同交互作用':'Joint interaction',
'RSM：共通の段階。PCM：選んだ相の水準別の段階。\nMMLは固定SD・EMです。指定区間は対象ごとに数値条件を確認します。':'RSM: common steps. PCM: separate steps for each level of the selected facet.\nMML uses EM with a fixed population SD. Specified intervals require numerical checks for each target.',
'数値条件を満たすセル・指定差に点ごとの区間と近似p値を表示します。\nHolm補正の対象は、この分析の全セルと指定差です。区間は同時区間ではありません。\n元の主効果モデルと固定較正の診断結果は保持します。公平性の判定にはなりません。':'Pointwise intervals and approximate p values are shown only for eligible cells and contrasts.\nHolm adjustment covers all cells and requested contrasts in this analysis. Intervals are not simultaneous.\nThe original main-effects model and fixed-calibration diagnostics are retained. This is not a fairness verdict.',
'元の観測得点を使用します。MFRM推定・採点範囲の指定は不要です。\n完全均衡交差、各セル1観測、単位重み。3相以上のデータは「固定相の扱い」で1水準抽出または完全均衡を確認した合計／平均集約を選べます。\nG／Phiは点推定です。区間は保留です。':'Uses original observed scores; no MFRM fit or rubric range is required.\nFully crossed, balanced data; one observation per cell; unit weights. For 3+ facets, choose single-level slice or verified balanced sum/mean composite above.\nG/Phi are point estimates. Intervals are withheld.',
'校正を固定した条件付きのセル推定です。主効果の厳しさとは別です。空白セルの状態と評定数は表で確認してください。色の大小からバイアスの有無・公平性を結論にしません。':'Conditional cell estimates with calibration held fixed, distinct from main-effect severity. Check cell status and counts in the table. Colors do not establish bias presence/absence or fairness.',
'主効果・段階・交互作用を同時に推定した別モデルです。区間は点ごとの近似区間です。表にはWald近似p値と、全セル・指定差を対象としたHolm補正を示します。空白や保留は「交互作用なし」を意味しません。':'A separate model jointly estimating main effects, steps and interaction. Intervals are pointwise approximations. Tables show Wald approximate p values and Holm adjustment across all cells and requested contrasts. Missing or withheld values do not mean no interaction.',
'指定した水準・差の点ごとの近似区間です。数値条件を満たした対象だけを描きます。保留の理由は表の注記を参照してください。有限標本被覆率・個人得点の区間は保証していません。':'Pointwise approximate intervals for specified levels and contrasts. Only numerically eligible targets are drawn. See table notes for withholding reasons. Finite-sample coverage and person-score intervals are not guaranteed.',
'元の観測得点のG（相対評価）とPhi（絶対評価）を、評定回数の計画別に示します。星は保存した計画です。完全均衡交差・1／2ランダム相・各セル1観測・単位重みの範囲です。負の分散成分は0として係数を算出し、符号付きの元の値も表に残します。区間は保留です。':'G (relative decisions) and Phi (absolute decisions) use original observed scores, projected for alternative counts. Stars mark the saved plan. Requires fully crossed, balanced data, one or two random facets, one observation per cell and unit weights. Negative variance components are set to zero for coefficients; signed raw values remain in the table. Intervals are withheld.',
'\n点線は選択した尺度の段階オフセットです。対象者・相の位置とは区別して読みます。PCMでは表示中の水準を確認してください。':'\nDotted lines show step offsets for the selected scale, distinct from person/facet locations. In PCM, check the displayed step-facet level.',
}

# Only application-generated messages use these patterns; data widgets are excluded.
TEMPLATES = [
('MFRM Desktop 試作 — {version}', 'MFRM Desktop prototype — {version}'),
('{label} · {rows}行 × {columns}列', '{label} · {rows} rows × {columns} columns'),
('計算中 · {rows}行', 'Computing · {rows} rows'),
('計画の相{index}の数', 'Planned count for facet {index}'),
('計画の {facet} の数', 'Planned {facet} count'),
('選択中：{count}列（1〜2列を選択）', 'Selected: {count} columns (select 1–2)'),
('選択中：1列（横軸: {x}）', 'Selected: 1 column (X-axis: {x})'),
('選択中：2列（横軸: {x} / 系列: {series}）', 'Selected: 2 columns (X-axis: {x} / Series: {series})'),
('{pair} · 範囲 ±{bound} · 極端水準除外', '{pair} · bound ±{bound} · omit extremes'),
('{pair} · 範囲 ±{bound} · 極端水準保持', '{pair} · bound ±{bound} · retain extremes'),
('{label}を整数で入力してください。', 'Enter an integer for {label}.'),
('{label}を数値で入力してください。', 'Enter a number for {label}.'),
('作文の人工データを保存しました：{path}', 'Synthetic writing data saved: {path}'),
('表と利用条件を保存しました：{path}', 'Table and notes saved: {path}'),
('表示中の対象・条件の図を保存しました：{path}', 'Figure saved with the displayed targets and settings: {path}'),
('プロジェクトを保存しました：{path}（入力データを含む）', 'Project saved: {path} (includes input data)'),
('レポートを保存しました：{path}', 'Report saved: {path}'),
('\n\n計算対象は{rows}行です。この内容で分析しますか？', '\n\nAnalyze the {rows} included rows?'),
]


def text(message, language='en'):
    if language not in ('en', 'ja'): raise ValueError('Unsupported desktop language')
    source = EN if language == 'en' else {v:k for k,v in EN.items()}
    if message in source: return source[message]
    for ja, en in TEMPLATES:
        old, new = (ja, en) if language == 'en' else (en, ja)
        pattern = re.sub(r'\\\{(\w+)\\\}', lambda m: '(?P<'+m[1]+'>.*?)', re.escape(old))
        match = re.fullmatch(pattern, message, re.S)
        if match: return new.format(**match.groupdict())
    return message

TEMPLATES += [
('バイアス · {pair}', 'Bias · {pair}'),
('指定区間 {target}', 'Specified intervals {target}'),
('共同交互作用 · {pair}', 'Joint interaction · {pair}'),
('カテゴリー確率{scope}', 'Category probabilities{scope}'),
('{facet} · 全水準 · {level}%', '{facet} · All levels · {level}%'),
('CSVの{row}行目の列数が1行目と違います。区切りや引用符を確認してください。', 'CSV row {row} has a different column count from the header. Check delimiters and quotation marks.'),
('得点に指定範囲（{minimum}〜{maximum}）の外の値があります。得点列と採点基準を確認してください。', 'Scores fall outside the specified range ({minimum}–{maximum}). Check the score column and rubric.'),
('空欄や数値でない得点などにより、{rows}行のうち{excluded}行が対象外になります。', '{excluded} of {rows} rows will be excluded because of missing or nonnumeric values.'),
('同じ対象者・相の組合せが重なる行が{rows}行あります。重複コピーや採点回の区別を確認してください。', '{rows} rows share person/facet combinations. Check duplicate copies and distinguish repeated rating occasions.'),
]
EN.update({
'CSVが読込上限（100 MiB）を超えています。':'CSV exceeds the 100 MiB input limit.',
'文字コードを読み取れません。Excelで「CSV UTF-8（コンマ区切り）」として保存し直してください。':'Cannot decode the file. In Excel, save it as CSV UTF-8 (comma delimited).',
'1行目に対象者・評価者等・得点の列名が必要です。CSVの形式と見本を確認してください。':'The first row needs column names for persons, facets and scores. Check the CSV example.',
'同じ列名が重複しています。それぞれ異なる名前にしてください。':'Column names are duplicated. Give each column a unique name.',
'CSVを読み取れません。空のファイルや引用符の閉じ忘れがないか確認してください。':'Cannot read the CSV. Check for an empty file or unclosed quotation marks.',
'列名だけでなく、評定データの行も必要です。':'The CSV needs rating rows as well as column headers.',
'対象者と得点の列、および評価者・課題等を1つ以上選択してください。':'Select person and score columns and at least one rater, task or other facet.',
'対象者・得点・相には、それぞれ別の列を選択してください。':'Select different columns for persons, scores and facets.',
'採点基準の最低点より大きい最高点を入力してください。':'The rubric maximum must exceed its minimum.',
'得点に無限大があります。有限の整数得点に修正してください。':'Scores include infinity. Use finite integer ratings.',
'得点に小数があります。平均点ではなく、1回ごとの整数の評定を入力してください。':'Scores contain fractional values. Enter individual integer ratings, not averages.',
'計算対象となる得点が2種類以上必要です。列の対応や空欄を確認してください。':'At least two different included scores are required. Check columns and missing values.',
'プロジェクトには入力データと識別子を含みます。表示・報告では再推定しません。':'Projects contain input data and identifiers. Viewing and reporting do not refit the model.',
'個人得点と一般表の利用条件を、指定区間の利用条件から判断しないでください。自動的な同等性の結論は保留です。':'Eligibility of specified intervals does not establish eligibility of person scores or general tables. Automatic equivalence conclusions remain withheld.',
'通常の表・Wright map等は元の主効果モデルです。共同交互作用は別モデルとして再推定した結果を専用の表・図に示します。':'General tables and Wright maps use the original main-effects model. Joint-interaction tables and figures show a separately refitted model.',
})
GUIDES_EN = {
'CSV_GUIDE': '''Use long-form UTF-8 CSV: one row per individual rating, with a header.
Person,Rater,Task,Criterion,Score
E01,R1,Task1,Content,4
E01,R1,Task1,Organization,3
E01,R2,Task1,Content,5

Person identifies the person being measured. Rater identifies the evaluator. Task and Criterion identify facets in the design. Score is one observed rating, not a mean or total. Column names may be in any language; select their roles after loading. Identifiers, including leading zeros, are preserved.

Save as CSV UTF-8 (comma delimited). Convert wide rater columns into rows. Do not encode an unobserved rating as zero; omit the row or leave its score blank. Check duplicate person/facet combinations; distinguish repeated rating occasions with the appropriate design. Disconnected raters or tasks may not be comparable.

For MFRM, enter the rubric minimum and maximum, including unobserved score categories. Missing or nonnumeric scores are excluded and reviewed before calculation. Fractional, infinite or out-of-range MFRM scores must be corrected. Independent G/D-Study uses original observed scores and requires complete balanced crossing with one observation per cell; for 3+ facets, you can explicitly condition on a single fixed level or compute a verified balanced sum/mean composite. All observations have unit weight; adding a weight column does not enable weighting.

The example CSV contains synthetic writing ratings, shared with the Web example. It is for practicing operations and interpretation, not real learner or rater data.''',
'SETTINGS_HELP': '''MFRM: choose RSM or PCM and JMLE or fixed-SD MML. Independent G/D-Study uses original observed scores and needs no MFRM fit or rubric range.

Columns and rubric: select distinct person, score and facet columns. Enter the rubric minimum and maximum, including categories absent in the data. Observed scores are context, not an automatically determined rubric.

Model: RSM shares steps; PCM estimates separate steps by the selected facet. MML uses EM with fixed population SD; quadrature and SD settings are retained.

Bias and intervals: use Add to register each request. Fixed-calibration cell estimates do not establish bias presence/absence. Specified structural intervals require fixed-SD MML and target-specific numerical checks. Joint interaction is a separate, jointly refitted model; eligible pointwise intervals and Wald approximate p values have explicitly stated conditions. Holm adjustment covers requested cells and contrasts, not simultaneous intervals.

G/D-Study: choose the object, score and one or two random facets. Requires full balanced crossing, one observation per cell and unit weights. For 3+ facets, Fixed-facet handling lets you extract one level or compute a sum/mean composite after verifying complete balanced elementary cells. Population descriptions are user declarations. Negative variance components are zeroed for coefficients while signed raw values are retained. G and Phi are point estimates; intervals are withheld.

Validation identifies input problems before calculation. Excluded rows and duplicate combinations require review.''',
'QUICK_HELP': '''Data → Settings → Results
Open CSV, Sample (3-facet) or G/D 2-facet sample, select columns and the rubric, then Analyze. Independent G/D-Study uses its own settings. Select figures and tables in Results. Reading guides and table notes describe interpretation.

English is the default interface language; 日本語 switches to Japanese without changing data or calculations. Figure labels remain in English for manuscript use. Explanatory figures include titles and notes. Manuscript mode moves overall titles and explanatory notes outside the graphic while retaining axes, legends, panel identities, data labels and unavailable-result messages. Title and note opens copyable text; manuscript PNG/SVG saves an accompanying .caption.txt file. Place the figure number/title above the image and necessary notes below it according to the journal's requirements; this switch is not a complete APA-compliance certification.

Projects include input data, identifiers, settings and results, and reopen without refitting. HTML reports contain all requested figures and tables. Tables also save as CSV with notes. New schema 2 projects cannot be opened by older applications; preserve originals when needed. Returning to Data or Settings keeps results, and changed settings are marked as belonging to a new request.

RSM/PCM, main effects, conditional bias, joint interaction and observed-score G/D-Study address different quantities. A drawn figure does not override an inference hold. G/D-Study intervals and multivariate G/D-Study are not available in this version.

The Web theory/help link opens the public Web entrypoint; it does not upload desktop data or results. Internet access is required and the public Web version may differ.

MIT license: research, teaching, professional and commercial use are allowed. Users are responsible for data handling, applicability and interpretation. No warranty is provided. Redistributors must retain copyright/license notices; third-party components have their own licenses.'''
}

EN.update({ja + suffix: en + ending for ja, en in [('差の左側', 'Left cell'), ('差の右側', 'Right cell')]
           for suffix, ending in [('：A水準', ': A level'), ('：B水準', ': B level')]})
