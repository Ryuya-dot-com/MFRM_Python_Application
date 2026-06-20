# FACETS比較ロードマップ

作成日: 2026-06-21

この文書は、MFRM StreamlitアプリをFACETSと比較し、強み、弱み、次の開発
優先度を整理するためのロードマップです。FACETSは、多相Rasch分析の実務で
広く使われてきた基準点です。このアプリは、Pythonだけで動く透明な
public beta / research previewとして、探索、教育、報告支援、研究
workflow prototypingを助ける位置づけです。

## エグゼクティブサマリー

- 現在の妥当な主張は「主要MFRM workflowのfunctional parityを目標にした
  beta」であり、FACETS等との数値一致や高リスク運用での妥当性は、個別の
  比較証拠なしには主張しない。
- 高 stakes な採点、認定、配置、雇用、制度判断では、FACETSまたは別の確立
  した手順で外部クロスチェックを行い、モデル仮定と差分の解釈を保存する。
- 直近の実装優先度は、結果への導線、claim boundary、FACETS出力crosswalk、
  比較fixture/parameterization noteの4点をそろえること。
- FACETS output crosswalkはP0 backlogで扱う。FACETSのtable/graph/file名を
  アプリ内のtab、table、figure、download、status、caveatに対応づける。

## 位置づけ

FACETSを優先すべき場面:

- 既存のFACETS制御ファイル、表、出力慣行に合わせる必要がある。
- 高 stakes な判断に直結する分析で、組織や査読者がFACETS出力を要求している。
- 係留、等化、bias/DIF、欠測の多い大規模運用を既存FACETS手順で継続している。

このアプリを優先しやすい場面:

- データ準備、列マッピング、事前チェック、図表確認をブラウザ上で誘導したい。
- 日本語/英語のhelp、用語集、tutorial、図の読み方を同じ画面で提供したい。
- 結果表、Wright Map、Category Probability Curves (CCC)、図ファイル、報告
  ドラフトをまとめて出力したい。
- 実装、テスト、検証方針、出力テンプレートをリポジトリ上で確認したい。

併用すべき場面:

- 論文、報告書、制度判断に使う場合は、アプリで探索・図表化・報告package化を
  行い、外部workflowの比較結果を保存する。
- 「FACETSと同じ」と読める表現を使う前に、Python出力、FACETS制御ファイル、
  FACETS出力、parameterization memo、tolerance policy、解釈メモを一式で
  残す。

## 比較の前提

FACETS公式資料では、受験者、項目、課題、評定者など複数facetの組み合わせ、
rating scale、partial credit、係留、fit評価、bias/DIF interaction、
weighting、欠測の多い設計、出力表、グラフ、Category Probability Curves、
Expected Score ICC/IRF、residual/score/anchor fileなどが扱われています。

したがって比較は「どちらが優れているか」ではなく、どのscenarioで、どの出力を、
どのparameterizationと許容差で比較できるかを明示する作業です。制約、
centering、quadrature、optimizer、潜在分散、extreme score処理、省略された
尤度定数が違えば、正当な理由で数値が変わることがあります。

## このアプリの強み

- 制御ファイルからではなく、sample data、貼り付け、upload、列マッピング、
  preflight、run、results、figures、downloadsへ進むguided UIがある。
- 日英localizationにより、MFRM用語に慣れていないユーザーにも説明を添えられる。
- 設定fingerprint、表、Quick result bundle、figure bundle、manuscript
  handoff、publication documentを出力できる。
- TAM、sirt、mirt、FACETS、`mfrmr`、R、`rpy2`を実行時依存にせず、外部検証の
  handoff artifactsを生成できる。
- RSM、PCM、bounded GPCM、JMLE、MML、latent regression、EAP scoring、
  plausible values、strict marginal diagnostics、residual PCA、
  bias/local interaction screening、anchor/linking review、simulation、
  posterior viewerなどを一つのPython workflowで扱う。

## 弱みとリスク

- FACETSには長年の運用実績、support、sample、出力慣行があり、このアプリは
  public beta / research previewである。
- FACETS経験者向けには、画面名、表名、download名とFACETS出力名の対応がまだ
  足りない。
- TAM/sirt/mirt handoff smokeは有用だが、FACETSを含む大規模・疎・係留・
  bias・disconnect・edge-category設計全体の検証programとしては未完成である。
- Streamlitは便利だが、大規模dataではbrowser memory、rerun、figure generation、
  ZIP export、hosted deploymentの制約を受ける。
- Bias/DFF/local interactionは、設計、精度、多重性、linking evidenceが揃わない
  限りscreeningとして扱う。Bounded GPCMとMML母集団scaleのclaimも、Rasch同等
  または安定推定として扱う前に感度分析が必要である。

## 直近の実装ステップ

1. 結果導線を固定する。
   Tutorial、Guided view、Helpで、sample data、run、Results、Figures、
   Wright Map、CCC、downloadへの経路を同じ表現にする。英日locale aliasと
   AppTestで主要labelを固定する。

2. Claim boundaryを画面内に置く。
   README、Help、Tutorial、Validation、Export周辺で、外部クロスチェックが必要な
   場面と、数値一致を保証しない理由を同じ短い表現で示す。

3. FACETS output crosswalkの仕様を先に作る。
   列は、`FACETS output name`、`user question`、`app location`、`download file`、
   `status`、`validation evidence`、`caveat`、`next action`を基本形にする。
   初期rowはmeasure table、fit statistics、Wright Map、CCC、expected score
   curves、residual file、score file、anchor file、bias/DFF interaction table。

4. 比較fixtureとparameterization noteを作る。
   balanced RSM、PCM、sparse ratings、extreme scores、missing categories、
   rater severity、task difficulty、anchors、group anchors、bias interactions、
   disconnected/subset riskを最小fixtureとして定義する。FACETS実行が手動でも、
   fixture生成とmetadataは決定的にする。

5. Review packetの形を決める。
   App reported、imported from FACETS、not yet implemented、not comparableを
   分け、比較reportなしで強いmethod claimが出ないようにする。

## ロードマップ

### Phase 0: 明瞭性とclaim boundary (0-1か月)

- 初回ユーザーが、主要結果、Wright Map、CCC、図downloadを迷わず見つけられる。
- FACETS経験者向けの「いつもの出力はどこか」入口をHelp/Tutorialに追加する。
- README、Help、Tutorial、Guided view、Validation、Export周辺のclaim boundaryを
  同じ語調にする。
- FACETS output crosswalkをP0 backlog itemとして仕様化する。

完了条件:

- Guided UIから主要結果と図に到達できることをAppTestで確認できる。
- Crosswalkの初期schemaと最初の対象出力リストがreview可能になる。

### Phase 1: 比較検証の背骨 (1-3か月)

- FACETS比較fixtureを定義し、Python出力、FACETS制御ファイル、FACETS出力、
  parameterization map、tolerance policy、version metadata、解釈メモを保存する。
- check種別を、strict numeric checkではなく、tolerance、rank、directional、
  non-comparableに分ける。
- 差分をbugと呼ぶ前に、期待される差分とparameterization上の理由を文書化する。

完了条件:

- 各scenarioについて、何を比較でき、何が比較対象外で、なぜそう扱うかを説明できる。

### Phase 2: FACETS風出力の読み替え (2-4か月)

- Crosswalkを実装し、FACETS table/graph/file名をアプリのtab、table、figure、
  download、validation statusへ対応づける。
- Measure/fit tableで、measure、SE、infit、outfit、ZSTD convention、
  displacement、count、observed/expected agreement、separation、strata、
  reliabilityのうち、実装・検証できる列を明示する。
- "FACETS-style packet" downloadにcrosswalk、measure/fit tables、Wright Map、
  CCC、expected score curve、residual exports、validation caveatsを束ねる。
- 外部FACETS runのmeasure/fit tableをimportし、Python runと横並びで確認できる
  補助機能を検討する。

完了条件:

- FACETSユーザーが、慣れた出力をcrosswalkで一つずつ確認できる。

### Phase 3: 実データ、linking、performance (3-6か月)

- 大規模upload、rerun、figure generation、ZIP exportのmemory/timeを測る。
- connectedness、subset、anchor drift、group-anchor、linking warningを強化する。
- 反復実施向けのtwo-run/multi-run calibration comparison exportを追加する。
- person、rater、institution、subgroup labelを含み得るexportにprivacy確認を追加する。

完了条件:

- データサイズの目安、fallback、linking evidence不足時の警告が文書化される。

### Phase 4: 報告packageとrelease gate (6-12か月)

- Reviewer-facing evidence binderとFACETS comparison appendix templateを作る。
- MML prior-SD sensitivity、latent regression diagnostics、plausible-value summaries、
  design simulation、posterior viewerを報告workflowに接続する。
- FACETS経験者とMFRM初学者のusability reviewを行い、supported/unsupported
  features、validated scenarios、tolerance results、high-stakes cautionsを含む
  release-readiness matrixを公開する。

完了条件:

- Release labelが、ユーザーテスト、比較fixture、supported/unsupported matrixで
  裏付けられる。

## 優先バックログ

P0:

- FACETS output crosswalk: FACETS output name、user question、app location、
  download file、status、validation evidence、caveat、next actionを持つ表を作る。
- Result、Figure、CCC、Wright Map、download routesをTutorial/Guided view/Helpで
  同じ表現にする。
- FACETS comparison fixtureとparameterization noteを追加する。
- Claim boundaryをREADME、画面内help、export caveatで同期する。

P1:

- Crosswalkを"FACETS-style packet" downloadに含める。
- FACETS table importとside-by-side comparisonを追加する。
- 大規模data/export performance、anchor/linking evidence gate、日英alias testを
  改善する。

P2:

- Design simulation workflowを拡張する。
- Posterior viewerとreport integrationを強化する。
- Reviewer binder templateとrelease-grade render checkを追加する。

## Claim policy

現時点で許容できる表現:

- "Standalone Python beta for exploratory MFRM workflows."
- "Targets functional parity for key MFRM workflows."
- "Can generate validation fixtures and optional external handoff artifacts."
- "Should be cross-checked before high-stakes decisions."

比較reportなしでは避ける表現:

- "FACETSと同じ結果を出す。"
- "FACETSと数値的に一致する。"
- "Operational high-stakes scoringで検証済み。"
- "Bias/DIFを因果的なfairness findingとして確定する。"
- "すべてのFACETS機能と出力慣行に対応している。"

## 参照情報

- FACETS product page: https://www.winsteps.com/facets.htm
- FACETS introduction/help: https://www.winsteps.com/facetman/introduction.htm
- FACETS Rasch model examples: https://www.winsteps.com/facetman64/raschmodels.htm
- FACETS output tables/files/graphs: https://www.winsteps.com/facetman64/outputtableindex.htm
- FACETS category probability curves: https://www.winsteps.com/facetman64/table8curves.htm
- Winsteps/FACETS support and learning-curve note: https://www.winsteps.com/comments.htm
- App validation stance: `../validation/README.md`
- App public positioning: `../README.md`
