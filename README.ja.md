# MFRM — 多相ラッシュモデル分析アプリ

[English](README.md) · [紹介サイト・ガイド](https://ryuya-dot-com.github.io/MFRM_Python_Application/?lang=ja) · [macOS試用版](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-macos-20261030) · [Windows試用版](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-windows-20261030) · [Web版を開く](https://mfrmpythonapplication-xodghmcghngx8cz5zmhfeg.streamlit.app/) · [デスクトップ版の使い方](docs/desktop_user_guide.ja.md)

作文・スピーキング・実技試験などの評定データから、対象者の測定値、評価者の厳しさ、課題の難しさ、得点カテゴリーの使われ方を分析するアプリです。結果を図表で確認し、レポートを保存できます。

最初は作文の人工データを使ったサンプルを試してください。FACETS・ConQuest・Rのインストールは不要です。

## 利用を始める

| 利用方法 | 入口 |
|---|---|
| ブラウザーで試す | [Web版](https://mfrmpythonapplication-xodghmcghngx8cz5zmhfeg.streamlit.app/)。公開サービスでは人工データや適切に匿名化したデータを使います。 |
| macOS（Apple Silicon）のデスクトップ版を試す | [macOS試用版をダウンロード](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-macos-20261030)し、SHA-256を確認して `MFRM-Desktop.app` を展開します。Pythonは同梱されています。 |
| Windowsのデスクトップ版を試す | [Windows試用版をダウンロード](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-windows-20261030)し、[デスクトップ版の使い方](docs/desktop_user_guide.ja.md)に沿って全体を展開してEXEを開きます。Pythonは同梱されています。 |
| Rスクリプトで分析する | CRANパッケージ [`mfrmr`](https://cran.r-project.org/package=mfrmr)（[ドキュメント](https://ryuya-dot-com.github.io/mfrmr/)）。Python版とは独立した実装です。 |
| 自分のパソコンでWeb版を動かす | [インストール](README.md#install)・[起動](README.md#run)の手順（英語） |
| 開発・ビルドに参加する | [開発者向け手順](CONTRIBUTING.md)。デスクトップ版をソースから起動する手順もあります。 |

**配布状況：macOS（Apple Silicon / M1以降）および Windows 11 x64 向けの試用版を公開しています。** 上のRelease案内からアプリZIPを取得できます（Windows版ZIPはDropbox経由、macOS版ZIPはRelease添付）。このページで説明するデスクトップ機能を一つにまとめた、アプリ表示版0.2.19-betaの未公証・未署名プレビューです。ファイル名・SHA-256・確認した環境・既知の制約は各Releaseに記載しています。Intel MacやWindows on ARM、他のOSバージョンでの網羅確認は未完了です。別の **v0.2.19-beta** ReleaseはWeb版とPythonソースであり、デスクトップ版のダウンロードとは異なります。

GitHubの **Source code (zip)** はソースコードであり、デスクトップ版のインストーラーではありません。Web版の実行中の版・ソースは **アプリ情報** で確認できます。ソースのブランチと公開版では、利用できる機能が異なる場合があります。

## できること

| 目的 | 機能 |
|---|---|
| 対象者・評価者・課題などの違いを調べる | RSM／PCM、JMLE／固定母集団SDのMML |
| 採点尺度と適合度を見る | Wright map、Pathway map、カテゴリー曲線、推定値・適合度表 |
| 要因の組合せを調べる | 条件付きバイアス診断。デスクトップ試用版には別モデルの共同交互作用分析もあります |
| 指定した効果や差の不確実性を見る | 条件を満たす固定SD MMLの水準・事前指定差の点ごとの区間 |
| 評価者数・課題数の計画を比較する | 対応する完全均衡設計のG/D Study |
| 結果を保存・報告する | Web版のレポート・ダウンロード、デスクトップ版のプロジェクト・HTML・CSV・PNG・SVG |

デスクトップ試用版にはこれらを一つのアプリにまとめ、英語／日本語、説明付き／論文掲載用の図の切替を用意しています。Web版には、bounded GPCM、潜在回帰、アンカー・リンキング、残差PCA、シミュレーションなどの追加機能があります。機能があることと、すべての出力から統計的な結論を出せることは別です。

## サンプルで試す

**Web版：** サンプルの学習ガイド、または **ガイドなしで続ける** から入り、組み込みの作文例を分析します。結果の最初の説明を読み、目的に合う図表と「レポート・出力」へ進みます。[作文サンプルの読み方](docs/apa_manuscript_template.md#writing-assessment-first-report)も参照してください。

**デスクトップ試用版（macOS・Windows）：** **データ → サンプル → 次へ → 分析** と進みます。英語表示のときは **Data → Sample → Next → Analyze** です。[デスクトップガイド](docs/desktop_user_guide.ja.md)に図表の確認と保存を説明しています。

## データの準備

1行を、1人に対する1人の評価者による1つの課題・観点の評定にします。[作文の見本CSV](examples/desktop_writing.csv)と[CSVの準備方法](docs/desktop_user_guide.ja.md#自分のcsvを使う)を参照してください。

MFRMでは整数得点を使い、今回出現しなかった得点も含めて採点基準の最低点・最高点を設定します。未採点を0点として扱わないでください。デスクトップ版はCSVを読み込み、Web版では表の貼付け、TSV/TXT、Excel、Parquet、JSONも利用できます。読込後は行数と列の役割を確認します。

## 結果を読むときに

- 一般のMML表、個人得点の不確実性、推定した母集団SDの不確実性、自動的な同等性判断には利用上の制限があります。専用の指定区間は、それぞれの条件を満たす対象に限って提供されます。
- 条件付きバイアス診断だけで、バイアスの有無や有意性は確定できません。共同交互作用は別モデルであり、その仮定と適用範囲を確認します。
- G/Dは完全均衡の1～2ランダム要因設計に対応します。多変量G理論とG/Phiの区間推定は未実装です。
- SEや区間の空欄は0ではなく、利用できないことを示します。実装された方法の検証範囲に由来する制限は、同じ分析を繰り返しても解除されません。

出力ごとの説明を読み、[Windowsガイドの対応方法](docs/desktop_user_guide.ja.md#困ったとき)を確認してください。[分析の適用範囲・推論の説明（英語）](docs/analysis_scope.md)では、近似区間の仮定、他ソフトとの違い、入力不備と方法上の制約を区別しています。

## 保存と共有

デスクトップ版では、後で開き直すためのプロジェクト、共有用のHTML、表のCSV、図のPNG／SVGを保存できます。[保存・再読込](docs/desktop_user_guide.ja.md#保存する)を参照してください。論文掲載用の図では、画像と一緒に保存される図注ファイルも使います。

Web版では、必要なレポートや図表を「レポート・出力」から保存します。公開WebのPNG生成が利用できない場合は、案内されるHTML形式を利用してください。ブラウザーのセッションはバックアップではありません。

プロジェクトや報告には入力データ・対象者や評価者のIDが含まれることがあります。共有前に内容を確認してください。機密データは管理されたローカル環境で扱い、所属機関の規定に従ってください。公開Webへのアップロードは遠隔のサーバーで処理されます。分析IDは暗号化・匿名化の仕組みではありません。

## 利用条件・引用・問い合わせ

[MITライセンス](LICENSE)で、研究・教育・商用に利用できます。再配布時は必要な著作権・ライセンス表示を保持してください。ソフトウェアは無保証であり、データの扱い、モデルの仮定、結果の解釈・利用は利用者の責任です。[ライセンスの説明](LICENSE_NOTICE.md)も参照してください。

研究で利用した場合は、実際に使用したバージョンとソースを記録してください。[引用情報](CITATION.cff)と[引用方法](README.md#citing-this-software)を用意しています。ソフトウェアの引用は統計手法の引用を代替せず、MITライセンスに追加する利用条件でもありません。

不具合や質問は[GitHub Issues](https://github.com/Ryuya-dot-com/MFRM_Python_Application/issues)へ、OS、版、操作手順、メッセージを添えてお知らせください。ファイルや画面には個人を特定できる情報を含めないでください。
