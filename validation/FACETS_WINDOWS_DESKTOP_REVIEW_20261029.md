# WindowsデスクトップとFACETSの追加点検 — 2026-10-29

Windows 11、Python 3.12.16、FACETS 4.5.0で実行した限定的な確認である。
推定の比較対象は `9b761c1605425153c264792c4825bf35983a9d44`。
本変更はScorefileの行対応と外部出力parserを修正するもので、推定式・既定の収束条件は変更していない。

## 修正した不具合

1. `compute_scorefile` が確率行列を新しいRangeIndexで作成していた。
   欠損得点の除外や絞り込みで入力の行ラベルに隙間があると、横結合が行を増やし、
   `Length of values ... does not match length of index ...` で停止した。
   行ラベルの集合が同じでも順序が異なる場合は、確率を別の観測へ結び付けるおそれがあった。
   確率表へ観測表のindexを付与し、行の順序と対応を保つよう修正した。
2. FACETS Scorefile parserが `PtBis` を `PtMeasure` として読み、実出力の `PtMea` を無視していた。
   保存した89水準の `PtMeasure` が全件欠損になることを再現した。
   `PtMea`、`PtBis`、`PtBisI` をそれぞれ `PtMeasure`、`PtBiserial`、`PtBiserialIncluded` に分離した。
   実出力の見出しは[公式Scorefile説明](https://www.winsteps.com/facetman64/scorefileoutput.htm)でも確認した。

修正前に行対応の5ケースと相関parserの4ケースが失敗し、修正後に通過した。
Windows CIにもこれらの検査を加えた。外部FACETSのライセンスをCIの前提にはしない。
関連156件を `-W error` で一括実行し、追加した欠損CSVの実spawn検査を含む12件も通過した。
重複を除く対象は157件である。全リポジトリのテスト成功を意味しない。

## 人物推定・SE・適合度の比較

既存のPCM既知真値生成器（seed 731101）の共有閾値データをRSM/JMLEで分析した。
80人・1,920評定、Person/Rater/Task/Criterionの4ファセット、四カテゴリ、単位重み、
アンカーなしの完全平衡データである。Pythonはデスクトップと同じ `run_basic_analysis` を使った。
FACETSはBATCH=NOの表示あり経路で出力安定を確認して通常終了した。

| 項目 | 通常設定での最大絶対差 | 厳しい停止条件での最大絶対差 |
| --- | ---: | ---: |
| 人物を含む89水準の推定値 | 0.00770325 logit | 0.0000286492 logit |
| 同じ89水準のSE | 0.000693982 | 0.00000303279 |
| Infit | 0.00868107 | 0.00494165 |
| Outfit | 0.0128127 | 0.00494853 |
| Infit ZSTD | 0.0708168 | 0.00492159 |
| Outfit ZSTD | 0.0684886 | 0.00496872 |

通常設定はPython `reltol=1e-6`、FACETS `Convergence=0.5,0.01,0,0`。
追加条件はPython `reltol=1e-12`、FACETS `Convergence=0.000001,0.000001,0,0`。
Pythonの上限はどちらも400反復である。
厳しい条件では適合度4列×89水準がすべてFACETSの表示値±0.005に入った。
これは今回の表示精度での整合であり、FACETSの丸め前の値との完全一致ではない。
既定条件での差が減ったことは数値的な停止条件への感度を支持するが、任意データへの保証ではない。
SEの近さから区間被覆率を認定することもしない。

## 相関は定義を分ける

現在のアプリの `PTMEA` は、非人物ファセット内の観測得点とPerson測定値の重み付き相関であり、
人物行は未計算である。FACETSの `PtMea` は、その観測を構成する符号付き測定値の和との相関である。
したがって、同じ意味の値として引き算して合否判定しない。
[FACETS公式の計算法](https://winsteps.com/facetman64/gstats.htm)に従い、今回の単位重みRSMデータで
後者を独立に計算すると、厳しい停止条件の89水準すべてが表示値±0.005に入った。
アプリ既存の `PTMEA` とFACETS `PtMea` の差は非人物9水準で最大0.0858129だった。
この違いを消す目的でアプリの統計量を無断で置換してはいない。

**未解決：** 公式説明の期待相関式を独立に評価した値と、人物の `PtMeExp` 表示には最大0.0197167の差が残った。
この追加計算は診断用であり、アプリの推定・既存PTMEA計算には使っていない。
期待相関の同等性は未確認のままとし、出力へ合わせるための補正式は導入しない。

## 欠測・非連結の境界

既存の `facets_pcm_boundary_pilot.py` を新しい出力先で再実行した。
完全／リング型の連結した計画欠測／非連結の負の対照、共有／異質閾値、2反復を用い、
同じ12データセットをPCMとRSMで計24組比較した。

- 24/24組の実行が完了し、比較可能な16/16組は既存の直接一致ゲートを通過した。
- 主効果144個の最大差0.0116210 logit、閾値72個の最大差0.00926145 logit。
- カテゴリ度数144/144が一致した。
- 非連結の8/8組で独立した設計検査のnullityが1以上となり、Pythonは推論を保留した。
  同じ8組でFACETSは `Subset connection O.K.` を表示した。
- 非連結の1例をデスクトップ経路でも実行し、保存・復元後も
  `stat.jmle.identification_unqualified` が残り、89水準すべてでSE/CI表示が保留された。

これらは2反復の境界確認であり、欠測一般、推定法の優劣、誤検出率や検出力の認定ではない。
FACETSの短い接続表示とアプリの全制約を含む識別検査は別の保証として扱う。

## 保持した記録と再実行

今回のローカル証拠は `C:\MFRM_Run` 配下に保持した。

| 場所 | 内容 |
| --- | --- |
| `fx26` | 前段のPCM/JMLE 4条件、同一入力、FACETS生出力、比較CSV |
| `fr26` / `fs26` | 通常／厳しい停止条件のRSM入力、FACETS生出力、デスクトッププロジェクト |
| `fb26` | 境界パイロット24組の入力・仕様・生出力・失敗台帳・集計 |
| `fd26final` | 修正後parserの比較CSV、相関定義検査、非連結プロジェクト、review.json |
| `facets_desktop_recheck.py` | RSM比較のローカル実行コード（--output、--reltol、--facets-convergence） |
| `extended_facets_review.py` | 保存済みRSMの再計算とデスクトップ保留確認 |

集計とソースhashは [facets_windows_review_20261029.json](facets_windows_review_20261029.json) に記録した。
生出力は本リポジトリへ同梱していない。PCM／境界比較は既存スクリプトから再生成できる。

```powershell
python validation/facets_pcm_known_truth_smoke.py --output-dir C:/MFRM_Run/fx-new
python validation/facets_pcm_boundary_pilot.py --output-dir C:/MFRM_Run/fb-new
python -m pytest -q -W error tests/test_scorefile_row_alignment.py tests/test_operating_characteristics_facets_adapter.py
```

出力先は未作成の短いパスを使う。FACETSは既存インストールとライセンスが必要である。
Mac/ConQuestの再実行は今回の対象外であり、過去のMML比較の資格をこの確認で広げない。

更新EXEは `C:\MFRM_Run\MFRM-Desktop-20261029-2\MFRM-Desktop.exe` に配置した。
SHA256: `e7bfba26414297df7c3dc21dfd704de1a614139dcd97bc6af6fe50b65396bc5f`。
同梱ソースhash、起動、実Worker、3図入りHTMLとPNG、非連結プロジェクトの再推定0の復元を確認した。
EXE検証のstderrは空だった。CSVダイアログ操作は自動操作ツールのウィンドウ活性化エラーで中止し、
欠損CSVの解析完了は実spawnテストで確認した。詳細は[デスクトップ受入記録](../docs/desktop_prototype.md)に記した。
