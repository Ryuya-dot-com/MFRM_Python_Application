# APA manuscript template

Use **Report & Export → Start a paper** after fitting. The English Word
template opens with a title page and abstract, then guides the author through
the introduction, Method, Results, Discussion, and References. No research
context form is required before downloading. A reusable blank copy is available
as [Word](templates/mfrm_apa_manuscript.docx) or
[Markdown](templates/mfrm_apa_manuscript.md).

For a first worked example, see [作文評価から最初の報告を作る](#writing-assessment-first-report) (Japanese).

<a id="report-tables-and-notes"></a>

## 表と図から論文を書き始める

分析後の「レポートと出力」では、最初に基本レポートを作れます。入力データが残る結果では、分析設計、
ファセットの推定値・適合度、カテゴリー度数、信頼性の表と、ファセット図・カテゴリー確率曲線を収録します。
計算できない出力は含まれません。設計表の評定数と人数は別の単位であり、評定行を独立した参加者数とは扱いません。
重み合計も評定数と分けます。共有Personを持つ評価者ペアの集計は、観測された配置の記述であり、
それだけで連結性・同一課題への一致・測定精度を保証しません。カテゴリー表には記録された得点範囲を含め、未使用のカテゴリーも0件で示します。
個人別の行、追加の診断図、交互作用表、方法・執筆の補足資料は、詳細版を選ぶと確認できます。
保留の注意は基本版にも残り、詳細版を選んでも推論に使える範囲は変わりません。

これは論文の必須項目を自動的に満たすセットではありません。まず研究の問いを決め、方法には
対象者・評価者の割り当て・尺度・モデル・推定法・同定制約・欠測等を記します。結果は推定値だけでなく、
適合度、不確実性、カテゴリーの使用状況、残差の点検と併せて解釈します。
下の作文例の「方法・推定値・カテゴリーの3表」は、目的に応じて組み立てる候補です。
基本レポートにも設計とカテゴリー度数を収録するため、同じ分析の値から表を選び、研究目的に合わせて編集できます。
入力が残らない古い結果では設計・カテゴリー表を捏造せず、利用できる表だけを収録します。

「報告書を作る」の**「論文にする前に確認すること」**を開くと、分析設計・カテゴリー・残差の
資料を探すための表名と、次の画面へ進むボタンがあります。「確認事項」で保留の理由を読み、
「論文を書き始める」で英語の原稿テンプレートを保存できます。テンプレートは記録済みの人数や
モデル、記録済みの推定設定、観測された配置とカテゴリー度数を埋めますが、研究の問い・設計の理由・倫理的配慮・解釈は著者が記入します。
最適化器が終了した事実と、停留性・積分精度・区間の妥当性は分けて記載します。
現在の校正に解釈の保留がある場合、結果の下書きにも反映します。
結果の文章案は、移動先で必要なときに表示してください。

判断に閾値を使う場合は、その根拠と比較対象を方法に記し、丸め前の値で判断します。
表示桁を揃えた表の一致だけでは、境界付近の判断が同じとは限りません。2人の評価者の差のSEには
両者の共分散も必要です。各SEだけから差のCIを作ったり、別々のCIの重なりで差を判定したりしません。
推定の終了が失敗・不明の場合、CSVの技術的SE・CIが残っていても、現在の保留理由を優先します。
古い保存結果のReadyやConvergedという表示だけでは利用可能とは判断できません。

作文評価では、各評価者が担当した人数だけでなく、どの能力層・課題を担当したか、共通の受験者や
答案でどう接続したかも報告します。同じ答案を2人が評定しても、独立した答案が2つ増えるわけではありません。
集計した担当人数や重なりの表だけでは、能力に関連する割当や、Person×Taskに共有される変動は判別できません。
TAMとの数値的一致や小さいSEも、それらの仮定の適切さを保証しません。
今回の限定的な[仮定・境界の点検](mml_quadrature_cross_engine_review.md#robustness-boundary-review-20261023)は、
こうした点を調べるための例であり、実データの仮定違反を自動検出する機能ではありません。

詳細版の方法・執筆の補足では、強調や箇条書きの階層を保って表示します。箇条書きの
「文章案」「避ける表現」「次の確認」は執筆を助ける資料です。すべてを論文へ転記する必要はありません。
校正全体の解釈が保留されている場合は、信頼性・マップなどの関連する文章案も保留になります。
確認用の数値が残っていても、確定的な結論には使わず、示された理由と入力・設定を点検してください。

Word・PDF・HTMLの結果表は、次の構造に揃えています。

1. **Table 1**のように、番号を独立した行に太字で置きます。
2. 次の行には、*Facet Estimates and Fit*のような短い説明的な表題を斜体で置きます。
3. 表本体は縦罫線やセルごとの囲みを使わず、上端・見出しの下・下端の水平線で区切ります。
4. 表の下の*Note.*に、単位、略語、区間の水準、推定法の制約、解釈保留の理由を記します。

例えば表題を「推定法の注意をすべて含む長文」にせず、Noteに
「推定値の単位はロジット。SEは他の推定値を固定した条件付き近似である」のように説明します。
実際の方法に応じて書き分け、個人得点の事後SDを標本SEと呼ばないでください。
全行に共通する区間水準や不確実性の種類はNoteに移し、行によって異なる種類は列に残します。
区間の上下限は角括弧で1列にまとめます。50行を超える表は表示範囲をNoteに明記し、全行はCSVで確認できます。

この表の構成は[APA第7版の公式Table Setup](https://apastyle.apa.org/style-grammar-guidelines/tables-figures/tables)に基づきます
（確認日: 2026-10-16）。区間の列の構成例は[公式Sample Tables](https://apastyle.apa.org/style-grammar-guidelines/tables-figures/sample-tables)も参照してください。
表の体裁を整えることは、論文全体のAPA準拠、研究の妥当性、投稿先の要件を保証するものではありません。
表題は研究の対象が伝わる表現へ調整し、本文での言及順に番号を付け直してください。

<a id="report-faq-and-references"></a>

## 出力の末尾にあるQ&Aと参考文献

Word・PDF・HTMLの基本版・詳細版は、結果と補足資料の後に**Questions and Answers → References**を収録します。
文書は英語です。画面の「報告書を作る」では「結果の読み方：4つのQ&A」を開くと、選択中の言語で読めます。
4項目は、スピーチ評価のoverfit、PCMとRSM、Wright mapの向き、バイアス分析です。
それぞれに定義、具体例、次に確認することを示しています。

Q&Aは結果の読み方を補う資料です。全文を論文へ転載する必要はありません。
診断が実行されたか、推定値の解釈が保留されていないかは、実際の結果とNoteで確認してください。
評価者の厳しさが似ていることとoverfit、得点の一致と分離信頼性は区別します。
RSMでは段階構造を共通にし、PCMでは指定した項目やファセットの水準ごとに異なる構造を認めます。
同じ得点ラベルだけでモデルを選ばず、ルーブリック・カテゴリーの使用状況・適合度などを点検します。
評価者のロジットの向きは設定に依存し、バイアスの有無についての結論は引き続き保留します。

Q&Aの引用に対応する文献は、重複を除いてReferencesへ収録します。
定義は[FACETS 64-bitの不適合の点検](https://www.winsteps.com/facetman64/diagnosingmisfit.htm)、
[分離信頼性](https://www.winsteps.com/facetman64/reliability.htm)、
[ファセットの向き](https://www.winsteps.com/facetman64/positive.htm)、
[バイアス出力](https://www.winsteps.com/facetman64/table13.htm)、
[WinstepsのRSM・PCMの構造](https://www.winsteps.com/winman/isgroups.htm)を照合しました（確認日：2026-10-16）。
日付のないオンライン資料は`n.d.`と参照日を記し、参照年を発行年に置き換えません。
古典文献はモデルの根拠として残し、現在のアプリの計算法や検証状況とは区別します。
他のソフトウェアのマニュアルは、このアプリの検定・区間の妥当性や、出力の数値的一致の証明にはなりません。

<a id="strata-and-fit-decisions"></a>

## Strata・underfitを意思決定につなげる前に

表のNoteで定義を確認し、画面では「Strata・分離度の定義と、判断前の確認」などを開くと、
式と具体例を読めます。「報告書を作る」にも「表の注釈を読む：Strataとunderfit」があります。

| 指標 | 意味 | その値だけでは決められないこと |
|---|---|---|
| Strata（H） | 分離度GからH = (4G+1)/3で求める、区別できる水準数の目安 | 実際の集団数、群を分ける境界、順位・受賞区分 |
| Underfit | MnSq > 1。残差の変動がモデルの想定より大きい | 対象者の能力不足、評価者の不正、除外の必要性 |
| Overfit | MnSq < 1。残差の変動がモデルの想定より小さい | 優れた成績、採点の正確さ、模倣の有無 |
| 分離信頼性 | 推定値の広がりと不確実性から求める、違いを区別する程度 | 評価者間の得点の一致、相互交換可能性 |

Strataは仮想的な分布と精度についての指標で、小数や実際の対象数を超える値にもなります。
例えばG = 2、H = 3でも、**学生を上・中・下の3群に分類できると確定したわけではありません**。
推定値の分布・不確実性・評定数を確認し、群分けが必要なら目的と境界の根拠を別に示します。
定義は[FACETS 4の信頼性・分離度・Strata](https://www.winsteps.com/facetman64/reliability.htm)を参照してください。

**1との大小と、注意表示の閾値を区別します。** MnSq = 1はモデルのもとでの期待値であり、
全評定が予測と一致したことを意味しません。1.2は期待値を上回りますが、アプリの目安の範囲内です。
注意表示は丸める前の値で、0.50未満、0.50以上1.50以下、1.50を超え2.00以下、2.00超に分かれます。
この範囲は点検の優先順位のための目安で、合否・除外の基準ではありません。
[FACETS 4の適合度](https://www.winsteps.com/facetman64/table7.htm)と
[不適合の点検](https://www.winsteps.com/facetman64/diagnosingmisfit.htm)の説明を、目的・設計と併せて使います。

例えば、ピアノ審査員のInfit = 1.1、Outfit = 2.1なら、予想外の評定が一部に集中していないかを確認します。
元の得点、課題・演奏者の担当、評価基準の適用を点検し、値だけで審査員を除外しません。
スケート評価でも、Strataや適合度から公式順位の変更や不正を結論づけることはできません。
これらは説明用の仮定例であり、実際の審査結果や競技規則を検証した例ではありません。

**計算法も併記します。** SEを使う集計では、G = √max(SD²−RMSE², 0) / RMSEです。
SDは標本標準偏差、RMSEはSEの二乗平均平方根です。MMLの個人EAPでは、
G = SD(EAP) / √mean(事後SD²)という別の集計を使うため、FACETSの値と直接同一視しません。
両方ともH = (4G+1)/3ですが、同じ式でHを作ることは、元のGの意味が同じである証明ではありません。
保留の値は確認用とし、空欄を0に置き換えないでください。

## From output to an argument

Start by explaining the assessment problem to a reader outside the narrow
specialty. For each question, connect the design and outcome to the answer it
can support. In Results, give the answer, magnitude, uncertainty, and conditions,
then point to a table or figure. In Discussion, return to the problem and
distinguish evidence from possible explanations. State the research question
and the evidence supporting the interpretation before refining the prose.

The downloaded template uses retained rating and person counts, category use,
facet names, model, recorded estimation settings, and app version. Unknown fields remain bracketed. It does
not copy person records, invent ethics approval, assert convergence, or write
substantive conclusions. The optional **Results draft and checks** reuses the
existing sentence audit and carries its critical-review, bias, population-scale,
and simulation conditions. The article template and that editing worksheet
have different purposes. Existing detailed manuscript notes remain available
in the broader export tools.

## Word and citation workflow

For the current development review of report coverage, see the
[reporting requirements and source review](mfrm_reporting_requirements.md)
(Japanese, October 16, 2026). It maps eight reporting areas to existing outputs
and distinguishes conditional SEs, posterior SDs, contrasts, and interval
methods. The Japanese worked example below is an initial guide; the piano and
skating examples and the broader acceptance reviews remain planned. The guide
does not qualify statistical inference.

The Word baseline uses US Letter, one-inch margins, 12-point Times New Roman,
double spacing, a running-head placeholder, page numbers, APA heading levels,
and a hanging-indent reference placeholder. The introduction starts below the
repeated title. Replace prompts with prose and follow the target journal's
title-page, abstract, anonymization, and disclosure requirements.

Keep only tables and figures needed to answer the questions. Use a bold table
or figure number, an italic descriptive title, and explanatory notes for units,
uncertainty, sample sizes, or adjustments. Define the contrast and its direction;
do not manufacture a contrast SE from separate SEs without covariance. Use the
existing table/figure exports as evidence and check the final journal layout.

In Word, insert citations and the bibliography with the Zotero plugin; choose
APA 7th edition in Document Preferences. Citation fields are not preinstalled
in this template, and the app does not connect to a user's library at runtime.
See [Zotero's Word documentation](https://www.zotero.org/support/word_processor_plugin_usage).
Method-literature candidates are shown separately from the manuscript's
References section so that uncited works are not inserted automatically.

<a id="writing-assessment-first-report"></a>

## 作文評価のサンプルから、最初の報告を作る

**問いは「評価者の厳しさはどの程度異なり、その評定はモデルとどう適合するか」です。**
厳しさと適合度を別々に読み、研究に必要な表と図を選んでみましょう。
ここでは、アプリに含まれる合成データ「作文評価（30×4×2×4、960件）」を使います。
実際の受験者や採点者についての研究結果ではありません。

### 入力から結果まで

1. 「サンプルで学ぶ」または「サンプルデータで実行」から始めます。ほかのサンプルを選んでいる場合は、作文評価の960件の例に戻します。
2. 1行は、ある受験者のある課題について、1人の評価者が1つの観点に付けた得点です。`Person`は30人、`Rater`は4人、`Task`は2課題、`Criterion`は4観点、`Score`は0–5の6カテゴリーです。960件は受験者数ではなく評定数です。
3. 分析条件を確認します。下の数値例は、乱数シード20240101の組込みデータを、RSM・JMLE、最大100反復、アンカーなし、単位重みで分析したものです。評価者・課題・観点は各ファセット内で総和0とし、正の値ほど得点を下げる向きです。自分の分析の条件が異なる場合は、その出力で読み進めます。
4. 「結果」の「最初に確認」を読みます。計算が終わっていても、解釈が保留されている場合は、表示された理由と必要な確認を先に読みます。
5. 「レポートと出力 → 表・図を保存」を開きます。「データ表」では、内容名で表を選ぶと読み方が表示されます。「図」では保存する図を選べます。診断や図を省略して分析した場合は、必要な項目を含めて分析し直します。

ここでいうRSMは、カテゴリー間の段階を共通とするモデルです。JMLEでは、対象者の能力と
評価者・課題などのパラメータを同時に推定します。モデルと推定方法は別の選択です。

### 最初に選ぶ3表と2図

次の5点は、この問いに取り組むための基本セットです。論文での表番号・図番号や分割方法は、投稿先の規定と研究目的に合わせて変えられます。
CSVをそのまま全列掲載する必要はありません。

| 報告物 | 画面で選ぶもの／元データ | 読む箇所と、説明する内容 |
|---|---|---|
| 表1：分析設計 | 「分析の概要」「分析に含めた評定・除外した評定」「評価者の割付と評定の重なり」 | 30人と960評定を区別し、人数・課題・観点・カテゴリー、採用／除外、誰が誰を評価したかをまとめます。この例は全員を4人が評価します。接続していることだけで十分な精度を保証しません |
| 表2：評価者の厳しさと適合度 | 「推定値・適合度・不確実性」から`Facet = Rater`の行を読む | `Estimate`の向き、`Infit`・`Outfit`、`SE_Method`・`SE_Status`を確認します。評定数は「要素別の適合度」の`N_Fit`でも確認できます。厳しいことと、モデルに適合しないことは別です |
| 表3：カテゴリーの使われ方 | 「カテゴリーの使用度数と適合度」、必要に応じて「段階の順序と間隔」 | 0–5の全カテゴリーの`Count`・`Percent`を確認し、平均推定値・適合度・段階を合わせて読みます。未使用カテゴリーがあれば残して説明し、警告一つで統合を決めません |
| 図1：Wright Map | 「図」のWright Map。座標は「Wright Mapの座標」 | 対象者・評価者等を共通のロジット尺度で見ます。高い位置の意味はファセットによって異なります。この例では、評価者の正の推定値は厳しい方向です |
| 図2：カテゴリー確率曲線 | 「図」のカテゴリー確率曲線。数値は「カテゴリー確率曲線のデータ」 | 横軸上の各位置で、どのカテゴリーが選ばれやすいかを読みます。このRSM例では、他のファセット効果を0とした共通尺度です。確率の合計は1です。各評価者の実測得点分布とは区別してください |

個別保存のCSV名には`mfrm_`が付きます。例えば表2は`mfrm_measures.csv`、表3は
`mfrm_category_diagnostics.csv`です。ZIP内ではそれぞれ`measures.csv`、
`category_diagnostics.csv`になります。「詳細な表・技術的な記録も一覧に含める」を選ぶと、
主要な表以外も探せます。一括出力には、現在の共有設定で保存できるすべての表が含まれます。

### 数値を短い説明にする

上記の条件で確認した出力では、R1の厳しさは−0.726、R4は0.467ロジットでした。
Infit／OutfitはR1が1.006／0.995、R4が1.052／1.060です。例えば次のように、
符号・範囲・適合度を分けて説明できます。

> この合成データのRSM・JMLEによる記述では、R4はR1より厳しい側に位置した。
> 両者のInfit・Outfitはモデルの期待値である1に近かった。
> この結果だけでは、評価者の正確さ・公平さや、厳しさの差の統計的な確実性までは判断しない。

この例のJMLEのSE・区間は、他の推定値を固定した条件付き近似です。
個別のSEから評価者間の差のSEを求めるには、共分散の検討も必要です。
MMLへ切り替えた場合は、個人得点の事後SDを同じ意味のSEとして説明しないでください。
評価者の分離信頼性を加える場合も、厳しさの差を区別できる程度であり、一致度ではないことを明記します。

### 解釈が保留されている例では

ある評価者がすべての評定に最高点を付けた例では、極端な得点パターンにより、モデル全体の推定（校正）の解釈が保留されます。
技術的な確認用CSVに値が残っていても、確定的な厳しさとして報告しません。
報告書では、その理由と再点検の必要性を記し、保留されたSE・区間の空欄を埋めないでください。
まず元データの入力・符号化、割付、実際の採点状況を確認します。警告を消す目的で評定を変更しません。

「レポートと出力 → 確認事項」で、現在の結果から説明できる範囲を確認します。
「報告書を作る」は表・図を含む出力の点検、「論文を書き始める」は研究の問いや方法を書く際に使えます。
通常例でも、報告書を生成できたことだけで投稿可能とは判断しません。保存した設定・表と文章を照合します。
操作中に用語が分からなくなったらHelpで確認し、同じ分析の「結果」または「レポートと出力」に戻ります。

### 保存前に、自分の言葉で説明できるか

- R1とR4のどちらが厳しい側にあり、Infit・Outfitからは何を読み取れますか？
- 評価者の分離信頼性が高い場合、評価者同士の得点が一致していると言えますか？
- JMLEの条件付きSEと、MMLの個人得点の事後SDは、同じ意味で説明できますか？
- SE・区間が空欄になった保留例では、何を確認してから報告しますか？

<details>
<summary>回答の着眼点</summary>

この例ではR4が厳しい側です。厳しさの推定値は得点を上げ下げする方向、Infit・Outfitは
モデルへの適合を表すため、厳しい評価者を不適合と決めつけません。分離信頼性は厳しさの違いを
区別できる程度であり、得点の一致度ではありません。

JMLEの条件付きSEは、他の推定値を固定した近似です。MMLの個人得点の事後SDは、
当てはめたモデルと観測された評定の下での個人の潜在特性の広がりを表し、EAP推定値の標本SEとは
区別します。どちらもNoteや元表の方法・状態を読み、検証されていない精度を断定しません。

空欄は誤差がゼロ、差がない、結果が問題なく使える、という意味ではありません。
保留の理由を読み、入力・符号化・割付・実際の採点状況を点検します。保存や詳細版への切り替えで
保留は解除されません。残差や相互作用は、研究の問いに必要な根拠を追加し、その解釈の条件も確認します。

</details>

## Literature checked on 2026-09-13

The first three records below were checked in the user's Zotero library.
Aryadoust et al. and Wind et al. were reviewed at metadata/abstract level;
Eckes was checked at bibliographic level. Their full papers were not retrieved
for this task. The resulting prompts are our application of these sources,
not a claim that every prompt is a requirement imposed by those papers.

| Source | Role in the template | Zotero item key |
|---|---|---|
| Aryadoust, Ng, and Sayama (2021) | Make dimensionality, local dependence, fit, and reliability reporting explicit | ZWCUBLMS |
| Wind, Jones, and Grajeda (2023) | Ask authors to explain sparse allocation and the limits of design-specific evidence | NBA428KJ |
| Eckes (2005) | Provide a relevant application to consult when framing rater-effect questions | WL25CIFE |

Aryadoust, V., Ng, L. Y., & Sayama, H. (2021). A comprehensive review of Rasch
measurement in language assessment: Recommendations and guidelines for research.
*Language Testing, 38*(1), 6–40. https://doi.org/10.1177/0265532220927487

Eckes, T. (2005). Examining rater effects in TestDaF writing and speaking
performance assessments: A many-facet Rasch analysis. *Language Assessment
Quarterly, 2*(3), 197–221. https://doi.org/10.1207/s15434311laq0203_2

Wind, S. A., Jones, E., & Grajeda, S. (2023). Does sparseness matter? Examining
the use of generalizability theory and many-facet Rasch measurement in sparse
rating designs. *Applied Psychological Measurement, 47*(5–6), 351–364.
https://doi.org/10.1177/01466216231182148

The general reporting scaffold also refers to APA's external primary source:
Appelbaum, M., Cooper, H., Kline, R. B., Mayo-Wilson, E., Nezu, A. M., & Rao,
S. M. (2018). Journal article reporting standards for quantitative research in
psychology: The APA Publications and Communications Board task force report.
*American Psychologist, 73*(1), 3–25. https://doi.org/10.1037/amp0000191
([APA source](https://www.apa.org/pubs/journals/releases/amp-amp0000191.pdf)).
The exact Appelbaum record was not found by the library search; semantic
fallback hits were not treated as matches. No Zotero records were changed.

This is an APA 7 writing scaffold with MFRM-specific prompts, not certification
that a study satisfies every applicable JARS module or journal requirement.
