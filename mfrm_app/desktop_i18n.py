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
('サンプルデータを保存しました：{path}', 'Sample CSV saved: {path}'),
('アンカーを{count}件読み込みました。', 'Imported {count} anchor entries.'),
('横長表を縦長へ変換しました：{rows}行 × {columns}列（空欄除外 {dropped}セル）', 'Converted wide table to long format: {rows} rows × {columns} columns ({dropped} blank cells omitted)'),
('プレビュー：{rows}行 × {columns}列（空欄除外 {dropped}セル）', 'Preview: {rows} rows × {columns} columns ({dropped} blank cells omitted)'),
('{facet} : {level}  =  {value} logits', '{facet} : {level}  =  {value} logits'),
('{facet} [Group={group}, 平均={mean}] : {levels}', '{facet} [Group={group}, Mean={mean}] : {levels}'),
('バイアス · {pair}', 'Bias · {pair}'),
('指定区間 {target}', 'Specified intervals {target}'),
('共同交互作用 · {pair}', 'Joint interaction · {pair}'),
('カテゴリー確率{scope}', 'Category probabilities{scope}'),
('{facet} · 全水準 · {level}%', '{facet} · All levels · {level}%'),
('CSVの{row}行目の列数が1行目と違います。区切りや引用符を確認してください。', 'CSV row {row} has a different column count from the header. Check delimiters and quotation marks.'),
('得点に指定範囲（{minimum}〜{maximum}）の外の値があります。得点列と採点基準を確認してください。', 'Scores fall outside the specified range ({minimum}–{maximum}). Check the score column and rubric.'),
('空欄や数値でない得点などにより、{rows}行のうち{excluded}行が対象外になります。', '{excluded} of {rows} rows will be excluded because of missing or nonnumeric values.'),
('同じ対象者・相の組合せが重なる行が{rows}行あります。重複コピーや採点回の区別を確認してください。', '{rows} rows share person/facet combinations. Check duplicate copies and distinguish repeated rating occasions.'),
(' G/D-Studyは計算できませんでした：{reason}', ' G/D-Study could not be computed: {reason}'),
]
EN.update({
'G/D: 観測設計のG・Phiと誤差':'G/D: observed_coefficients',
'G/D: 選択した計画（D-Study）':'G/D: selected_d_study',
'G/D: 分散成分表':'G/D: variance_components',
'G/D: D-Study計画予測表':'G/D: d_study',
'G/D: 計算条件と記録':'G/D: metadata',
'元の観測得点における観測設計の点推定です（G：相対決定、Phi：絶対決定、誤差分散・SEM）。G／Phiの区間は保留です。':'Observed-design point estimates on original scores: G for relative decisions, Phi for absolute decisions, and error variances / SEMs. Intervals for G and Phi are withheld.',
'同一の母集団・完全交差観測を仮定した、選択計画でのG・Phi・誤差分散・SEMの予測値です。':'Projected G, Phi, error variances and SEMs at the selected D-Study counts assuming the same universe and fully crossed observations.',
'均衡ランダム効果ANOVAの分散成分です。RawVarianceは符号付きのモーメント法推定値、VarianceはG／Phi算出用に負値を0へ補正した値です。最高次の交互作用は残差と分離できません。':'Balanced random-effects ANOVA components. RawVariance retains signed method-of-moments estimates; Variance clips negative values at zero for G/Phi. Highest-order interaction is confounded with Residual.',
'均衡ランダム効果ANOVAの分散成分です。負のRawVariance推定値が含まれており（真の分散0付近の標本変動または設計不適合の可能性）、G／Phi算出では0に補正しています。最高次の交互作用は残差と分離できません。':'Balanced random-effects ANOVA components. One or more RawVariance estimates are negative (possible sampling variability near zero or design misspecification) and were clipped at zero for G/Phi. Highest-order interaction is confounded with Residual.',
'計画水準数ごとの閉形式D-Study点予測です（IsObservedDesign=Trueが実際の観測水準数）。':'Closed-form D-Study point projections across planned facet counts (IsObservedDesign=True marks the observed counts).',
'設計条件の検証結果・得点基準・宣言母集団・計算設定の監査記録です（interval_inference_allowed=False）。':'Audit record of design checks, score basis, declared population and calculation settings (interval_inference_allowed=False).',
'測定対象以外に1〜2個の異なるランダム相が必要です。':'One or two distinct non-object random facets are required for a G-study.',
'必要な列が不足・重複しているか、予約名を含んでいます。':'Required columns are missing, duplicated, or use reserved component names.',
'保持された得点はすべて有限値である必要があります（自動除外は行いません）。':'Every retained score must be finite; no rows were silently removed.',
'相のラベルに欠測・空欄または文字列変換後の重複があります。':'Facet labels must be nonmissing, nonblank and unambiguous after text conversion.',
'測定対象を含むすべての相に2水準以上が必要です。':'Each facet, including the object of measurement, needs at least two levels.',
'同一セルの反復観測には対応していません。分析前に観測単位を定義してください。':'Repeated cells are unsupported; define the observation unit before analysis.',
'完全交差設計に欠測セルがあります（欠測は0点とみなしません）。':'The full crossed design is incomplete; missing ratings are not zero scores.',
'単位重み（Weight=1）の観測のみ対応しています。':'This G-study supports unit observation weights only.',
'計画水準数の上限は正の整数にしてください。':'The maximum planned facet count must be a positive integer.',
'分散分解で有限の分散成分が得られませんでした。':'The variance decomposition did not produce finite components.',
'有効な測定対象分散が0のためGおよびPhiを定義できません。':'G and Phi are undefined because usable object variance is absent.',
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
EN.update({
'サンプル集（8種）…':'Sample gallery (8)…',
'横長→縦長変換…':'Wide → Long…',
'アンカー・等化':'Anchors / equating',
'データ未選択':'No dataset loaded',
'「CSVを開く」または「サンプル集（8種）…」からデータを読み込んでください。':'Load ratings from "Open CSV" or "Sample gallery (8)…" to begin.',
'オフライン実データ分析ワークステーション':'Offline Real-Data Psychometric Workstation',
'1. データ':'1. Data',
'2. 設定':'2. Settings',
'3. 結果':'3. Results',
'分析ワークフロー：':'Workflow:',
'サンプルデータ集（人工データ・8シナリオ）':'Sample Dataset Gallery (8 Synthetic Scenarios)',
'試したい研究・実務シナリオを選択すると、列・採点基準（アンカー例ではアンカー値も）を自動設定して読み込みます。':'Select a scenario to load synthetic ratings with pre-configured columns, rubric bounds, and (where applicable) anchors.',
'このサンプルを読み込む':'Load selected sample',
'選択サンプルをCSV保存…':'Save sample CSV…',
'横長（ワイド形式）CSVの縦長変換':'Convert Wide Spreadsheet CSV to Long Format',
'評価者や観点が列方向に並んだ表を、1行＝1評定の縦長形式へ変換します（空欄セルは自動除外）。':'Melt wide rater/item columns into one-row-per-rating long format (blank unrated cells are automatically omitted).',
'行に残す列（対象者・課題など）':'ID / fixed columns to keep',
'縦に展開する得点列（複数選択）':'Score columns to melt (multi-select)',
'新しい相の列名':'New facet column name',
'新しい得点の列名':'New score column name',
'変換して読み込む':'Convert and load',
'変換せずそのまま読み込む':'Load as-is without converting',
'行に残す列と縦に展開する得点列を重複なく1列以上ずつ選択してください。':'Select at least one ID column and one non-overlapping score column to melt.',
'新しい相の列名と得点の列名を入力してください（既存ID列と重複不可）。':'Enter distinct nonblank names for the new facet and score columns.',
'固定アンカー（水準のロジット値を固定：Status A）':'Fixed anchors (fix individual level logits: Status A)',
'グループアンカー（指定水準群の平均ロジットを固定：Status G）':'Group anchors (fix mean logit across levels: Status G)',
'水準':'Level',
'読み込む':'Load',
'アンカー値':'Anchor logit',
'CSVから一括読込…':'Import CSV…',
'全解除':'Clear all',
'グループ名':'Group name',
'グループ平均':'Group mean',
'グループに含める水準（2つ以上選択）':'Levels in group (select 2+)',
'グループを追加・更新':'Add / update group',
'アンカーCSV（Facet, Level, Anchor または Estimate 列）を選択':'Select anchor CSV (Facet, Level, Anchor or Estimate columns)',
'アンカーCSVには Facet, Level, Anchor（または Estimate）列が必要です。':'Anchor CSV requires Facet, Level, and Anchor (or Estimate) columns.',
'グループアンカーには2つ以上の水準を選択してください。':'Select at least 2 levels for a group anchor.',
'アンカー指定時は主効果・適合度・条件付きバイアスを等化尺度上で推定します（固定SD MMLの指定区間・共同交互作用は総和ゼロ制約専用のためアンカー未使用時のみ有効）。':'With anchors active, main effects, fit, and conditional bias are estimated on the anchored scale (fixed-SD MML specified intervals and joint interaction require unanchored sum-to-zero facets).',
'連結性・アンカー等化サマリ':'Connectivity and Scale Linking',
'保持された Person × 各相の連結成分数および尺度識別に用いた固定・グループアンカー制約の一覧です。アンカー等化は実施回間の不変性を仮定するため、適合度・ドリフトも併せて確認してください。':'Retained bipartite Person × Facet connectivity counts and explicit fixed/group anchor constraints used to identify the logit scale. Anchor equating assumes invariance across administrations; review fit and drift.',
'作文評定（3相・完全交差・960行）':'Essay Writing (3 Facets, Crossed, 960 rows)',
'臨床OSCE評価（4相・1〜4点・576行）':'Clinical OSCE (4 Facets, 1–4 Rubric, 576 rows)',
'2コホート尺度連結（アンカー自動設定・540行）':'2-Cohort Scale Equating (Preset Anchors, 540 rows)',
'欠損・空欄セル混在データ（900行・空欄12行）':'Ratings with Missing Rows & Blank Cells (900 rows)',
'評価者スパイラル割当（充填率25%・768行）':'Spiral Rater Assignment (25% Density, 768 rows)',
'午前・午後シフト分断例（未接続サブセット診断・128行）':'Unlinked AM/PM Shifts (Disconnected Subset Warning, 128 rows)',
'横長Excel形式シート（評価者列 R1〜R4 を縦長変換）':'Wide Spreadsheet Layout (Convert R1–R4 Columns to Long)',
'G/D-Study 2相完全均衡サンプル（240行）':'G/D-Study Balanced 2-Facet Sample (240 rows)',
})
GUIDES_EN = {
'CSV_GUIDE': '''Use long-form UTF-8 CSV: one row per individual rating, with a header.
Person,Rater,Task,Criterion,Score
E01,R1,Task1,Content,4
E01,R1,Task1,Organization,3
E01,R2,Task1,Content,5

Person identifies the person being measured. Rater identifies the evaluator. Task and Criterion identify facets in the design. Score is one observed rating, not a mean or total. Column names may be in any language; select their roles after loading. Identifiers, including leading zeros, are preserved.

Save as CSV UTF-8 (comma delimited). Wide tables with separate rater/item columns can be melted into long format inside the app via "Wide → Long…". Do not encode an unobserved rating as zero; omit the row or leave its score blank. Check duplicate person/facet combinations; distinguish repeated rating occasions with the appropriate design. Disconnected raters or tasks may not be comparable without shared observations or anchors.

For MFRM, enter the rubric minimum and maximum, including unobserved score categories. Missing or nonnumeric scores are excluded and reviewed before calculation. Fractional, infinite or out-of-range MFRM scores must be corrected. Independent G/D-Study uses original observed scores and requires complete balanced crossing with one observation per cell; for 3+ facets, you can explicitly condition on a single fixed level or compute a verified balanced sum/mean composite. All observations have unit weight; adding a weight column does not enable weighting.

The example CSV contains synthetic writing ratings, shared with the Web example. It is for practicing operations and interpretation, not real learner or rater data.''',
'SETTINGS_HELP': '''MFRM: choose RSM or PCM and JMLE or fixed-SD MML. Independent G/D-Study uses original observed scores and needs no MFRM fit or rubric range.

Columns and rubric: select distinct person, score and facet columns. Enter the rubric minimum and maximum, including categories absent in the data. Observed scores are context, not an automatically determined rubric.

Model and anchors: RSM shares steps; PCM estimates separate steps by the selected facet. MML uses EM with fixed population SD; quadrature and SD settings are retained. Under Anchors / equating, you can specify fixed element logits (Status A) or group-mean anchors (Status G) and verify live connectivity and scale identification.

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
