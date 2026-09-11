# MML求積の比較：Python・mfrmr・TAM・ConQuest

確認日: 2026-09-11。ソース・公式資料・保存済み検証のレビューと、求積関数だけの実行照合。
新しいモデル推定、ConQuestの実行、独立した統計資格検証は実施していない。

**同じMMLモデルでも、有限個の点で近似した尤度は求積法によって異なる。**
最適化の収束、積分精度、標準誤差・区間の妥当性は別々に確認する必要がある。
今回の重要な収穫は、mfrmrに既存の比較証拠と感度確認機能があること、そして
点数を増やすだけでは解決しない範囲・重みの問題も確認できたことである。

## 確認した実体

| 対象 | 今回確認した版・範囲 |
|---|---|
| Pythonアプリ | `integration/unified-app`、計算ソースは `49f63c411a7d23e535dfc44c1de0ee34ba174a7e` |
| 開発中mfrmr | 指定されたルートの `development/`。`0.2.4.9000`、HEAD `6dfba8258403cb0ab2f86979e78b31519753e82f`、開発状態 |
| インストール済みmfrmr | `0.2.3.1`。今回の開発版ソースと別物なので、通常の `library(mfrmr)` を開発版の実行確認には使っていない |
| R・TAM | 実機 R `4.6.1`、TAM `4.3.25`。インストール済み名前空間の関数本体も確認 |
| ConQuest | 保存記録の `5.47.5`。ローカル実行ファイルのSHA-256は拡張計画記載のものと一致。今回の起動は未確認 |

mfrmrのソース・記録の基準ディレクトリは
`/Users/tohokusla/Dropbox/MFRM_Application/mfrmr/development`。
以下のmfrmrリンクは同じ親ディレクトリにあるこの読み取り専用の調査対象を指す。
ソースと古い実行記録の版を混同せず、今回の照合結果も過去の凍結研究に上書きしない。

## 共通の積分と、異なる数値近似

単位重み・一つの潜在能力を持つ比較可能なモデルでは、Personごとの反応ベクトルをまとめて

\[
\ell(\alpha,\beta,\sigma)
=\sum_p\log\int_{-\infty}^{\infty}
\underbrace{\prod_{j\in O_p}P(Y_{pj}=y_{pj}\mid\theta,\alpha)}_{L_p(\theta;\alpha)}
\frac1\sigma\phi\!\left(\frac{\theta-\mu_p}{\sigma}\right)d\theta,
\qquad \mu_p=X_p\beta
\]

を最大化する。`alpha` はファセット効果・ステップなど。欠測は観測集合 `O_p` から外す。
評定行ごとに積分してから積を取る式とは異なり、同じPersonの能力を共有する。
行内の反応重みと積分外のPerson重みも同じ操作ではないため、重み付き比較は別途対応を決める。

標準正規Gauss–Hermite（GH）の点 `z_q` と重み `w_q` なら、

\[
\ell_Q=\sum_p\log\sum_{q=1}^Q w_q L_p(\mu_p+\sigma z_q;\alpha).
\]

変数変換により正規密度とJacobianは重みに含まれる。ここで `dnorm` や `1/sigma` を
もう一度掛けると二重計上になる。積の過小化を防ぐ計算では、カテゴリーの対数確率を
Person内で足し、`log(w_q)` を加えてから `logsumexp` を取る。

| 対象 | 今回確認した求積 | 最適化・母集団の扱い |
|---|---|---|
| Python | `hermgauss` を標準正規へ変換し、点を `sigma` 倍。低水準API既定15点、前回の開発監査は31・61点 | 現行の自由SD経路はEM。共同停留性を確認するv2は開発用。母集団SDが動くと点も動く |
| mfrmr 0.2.4開発版 | Golub–Welschの固有値・固有ベクトルから標準正規GHを生成。公開 `fit_mfrm()` 既定31点 | 既定 `direct` は勾配法と必要時のpolishing。通常RSM/PCMは固定N(0,1)、母集団モデル有効時は `mu_p + sigma*z_q`。限定したEM・hybridもある |
| TAM 4.3.25の決定論的経路 | 既定 `nodes=seq(-6,6,len=21), snodes=0`。保存比較では31・61点等を明示した等間隔グリッド | 固定点で正規密度を評価するEM系。母集団平均・分散を固定するか推定するかを明示する |
| ConQuest | `gauss` はGH、`quadrature` は別方式。回帰変数なしでは前者、回帰変数ありでは後者が既定。保存比較は **`method=quadrature` を明示** | MMLはEM系。今回の記録から既定 `gauss` の結果まで検証済みとは扱えない |

TAMは `snodes>0` ならQMC、さらに `QMC=FALSE` なら通常のMonte Carloを選べる。
ConQuestの `quadrature` は点数と `minnode/maxnode` を持ち、既定範囲は−6～6、
`gauss` ではこの範囲指定は無視される。公式の呼び名が同じ「Gaussian quadrature」でも、
TAMの既定等間隔点をHermite多項式の根と同一視しない。
[TAM公式マニュアル](https://alexanderrobitzsch.r-universe.dev/TAM/doc/manual.html)、
[ConQuest公式コマンド参照](https://conquestmanual.acer.org/s4-00.html#estimate)

実機のTAM `tam.mml` → `tam_stud_prior(normalize=FALSE)` → `tam_calc_posterior` →
`tam_mml_compute_deviance` を追うと、1次元の決定論的経路の周辺確率は
`Delta * sum_q L_p(t_q) * dnorm(t_q, mu_p, sigma)` として評価される。
事後確率の行正規化と、事前密度をグリッド上で総和1に再正規化する処理は区別する。
後者を無断で追加すると、有限グリッドの母集団パラメータ依存の目的関数を変え得る。

mfrmrのGHは母集団の位置・尺度へ変換するが、Personの**事後モード・曲率に合わせる
適応的GHではない**。ファセットが複数あっても現在のモデルの積分は1次元である。
多次元化した場合の積分点数の増大を、ファセット数の増加と混同しない。

## 今回実行した求積関数の照合

Pythonの `gauss_hermite_normal` はASTから当該関数だけを、mfrmrはRの `parse()` から
同名関数の代入式だけを取り出して実行した。アプリの推定やインストール済み旧Rパッケージは使っていない。
点を昇順にそろえ、R出力のテキスト転送も含めて比較した。

| Q | 点の最大絶対差 | 重みの最大絶対差 | mfrmrでゼロになった重みの数 |
|---|---:|---:|---:|
| 31 | 1.42e−14 | 2.11e−15 | 0 |
| 61 | 5.33e−14 | 5.27e−15 | 4 |
| 121 | 5.68e−14 | 5.79e−15 | 36 |
| 181 | 5.33e−14 | 4.45e−15 | 74 |

Q31・Q61について、SDが0.5・1・2の場合の点の尺度変換と、正規分布の
重み総和・1次・2次・4次モーメントも検査した。点の絶対差1e−12、重み1e−13、
モーメント2e−11という今回の実装照合用許容差を満たした。統計的な許容差ではない。
標準正規GHの最大正の点はQ31で9.8934、Q61で14.4985。これらはGHの最外点であり、
等間隔法のような明示的な積分打切り境界とは異なる。

**新しく確認した注意点：微小重みのゼロと保存校正の要件が整合していない。**
Rの現在の固有ベクトル方式では、Q61の4点で重みが正の微小値からゼロになる。
対応するNumPy側重みの合計は約7.42e−37で、Q121・Q181でも約2.18e−36・2.47e−37。
通常の単位重みの反応尤度積は0～1なので、これらの項を落とす効果は絶対確率では小さい。
しかし非常に小さい周辺確率に対する相対誤差・対数尤度誤差は、この値だけでは保証できない。
この観測は過去のTAM差やPythonの停留性問題の原因を特定したものではない。

加えて、mfrmrの保存校正は生成した重みをそのまま格納する一方、検証時は
`any(weights <= 0)` を拒否する。Q61・121・181の今回生成した重みはこの述語に抵触し、
Q31は通った。**fit時の点数を増やすことと、保存校正のscoring点数を増やすことは別**である。
これは生成関数と検証条件の直接照合であり、保存校正API一式の新規実行試験ではない。
根拠: [生成・格納・検証](../../mfrmr/development/R/core-fixed-calibration.R)、
[GH生成関数](../../mfrmr/development/R/mfrm_core.R)。
対応時は重みを任意の小数に置き換えず、共通生成関数の数値安定性と高次数の回帰検査を検討する。
今回mfrmrには変更を加えていない。

微小重みと正重み要件の不整合は、旧インストール版をロードせず次で再確認できる。
`Rscript --vanilla` に渡し、最後の表を確認する。

```r
src <- "/Users/tohokusla/Dropbox/MFRM_Application/mfrmr/development/R/mfrm_core.R"
env <- new.env(parent = baseenv())
for (e in parse(src)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      identical(e[[2]], as.name("gauss_hermite_normal"))) eval(e, env)
}
q <- c(31L, 61L, 121L, 181L)
zeros <- vapply(q, function(n) sum(env$gauss_hermite_normal(n)$weights == 0), integer(1))
stopifnot(identical(zeros, c(0L, 4L, 36L, 74L))) # 今回の実機結果の再現検査
data.frame(q = q, zero_weights = zeros, strictly_positive = zeros == 0L)
```

今回の実行に使ったファイルのSHA-256:

- Python `streamlit_app.py`: `747150dfcf1e2429bc6a95add53ed0d28d535691a4f03e24064e4b53e94d3c17`
- mfrmr `R/mfrm_core.R`: `b2f78488f251f80f3a3f705b10c41631cbe7e6ba5b3e1cca284cf5660a1000e1`

## mfrmrの既存検証から使える根拠

以下は**保存済みの過去の観測結果**。新規実行や現HEAD全体の合格宣言ではない。

| 記録 | 結果と意味 |
|---|---|
| [TAM基本条件の再照合](../../mfrmr/development/inst/validation/tam-mml-core-current-head-record-0.2.4.md) | RSM/PCMの単純な対応条件ではQ61のdeviance差が約1e−12。モデル・座標を対応できる限定的な根拠 |
| [RSMストレス比較](../../mfrmr/development/inst/validation/tam-mml-release-stress-record-0.2.4.md) | 21データセット・42比較。mfrmrは42/42で収束を報告したが、ソフト間基準は12/42、Q31→Q61の安定性基準は6/21だけが通過。TAMの42/42は反復上限未到達の記録であり、共同スコア検査とは異なる |
| [密度・裾範囲の追加診断](../../mfrmr/development/inst/validation/tam-mml-density-diagnostic-record-0.2.4.md) | Q121→181の変化は21/21で既存review基準内。それでもソフト間一致は17/21。残る例でTAMの範囲を広げると差が縮小 |
| [PCM条件付きストレス比較](../../mfrmr/development/inst/validation/tam-pcm-mml-conditional-stress-record-0.2.4.md) | 15データセット。ソフト間基準の通過はQ31で3/15、Q61で6/15、Q181で15/15。Q31→61のmfrmr最大EAP変化は0.0583、事後SD変化は0.0731。RSMだけの現象ではない |
| [ConQuest・TAMのアルゴリズム監査](../../mfrmr/development/inst/validation/external-mml-algorithm-correlation-audit-record-0.2.3.md) | 明示したConQuest `quadrature` とTAMは固定グリッドMML/EMの共通部分を持つ。mfrmr既定の直接最適化と同じアルゴリズムではない。高相関は同等性の基準ではない |
| [対数中心化した連続積分](../../mfrmr/development/inst/validation/conquest-p2-log-centered-continuous-oracle-observation-record-0.2.3.md) | 凍結された13条件でモード中心化・裾誤差を確認した数値基準にGHが一致。以前の連続積分との差が残った例もあるため、「連続積分」と呼ぶだけで正しい基準とはしない |

裾範囲の具体例では、Q301を保ったままTAMの範囲を−6～6から−8～8に広げると、
強制極端反応の1条件のdeviance差は約2.90e−4から1.92e−10になった。
点を密にすることと裾を広げることを分ける必要を示す。結果を見た後の機序診断なので、
元の12/42や6/21を合格に書き換えない。またQ181を万能の既定値にする根拠でもない。

連続積分の基準はモード探索、対数中心化、裾の評価を伴うが、確率計算をmfrmr本体と共有している。
独立ソフトによるモデル検証とは異なり、報告された積分誤差も区間演算による厳密証明ではない。
Pythonへの再利用時はR側で同じ数値を見るだけでなく、反応確率・同定・尤度の独立照合を残す。

## 0.2.4で実装済みの対応と、Python側への含意

mfrmrの [修正記録](../../mfrmr/development/inst/validation/mml-quadrature-remedy-record-0.2.4.md)
と現在の [公開API](../../mfrmr/development/R/api-quadrature-sensitivity.R) を照合した。
`mml_quadrature_sensitivity()` はRSM・PCM・限定GPCMについて、同じデータ・モデルの再推定から
尤度、構造パラメータ、確率、EAP、事後SD等の変化を示す。既定 `c(31,41)` は比較の出発点であり、
十分な精度の宣言ではない。自動の安定性閾値を持たず、推論資格も昇格させない。

[保存校正の抽出](../../mfrmr/development/R/api-calibration.R) はこのreviewと、その中の
厳密に同じ最高次数fitを要求する。これは確認手順の強制であり、数値安定性の自動認定ではない。
推定と得点計算の積分は別に評価する設計になっており、現在のscoring既定は31点である。

母集団分散の座標も対応させる必要がある。mfrmrは `zeta=log(sigma^2)`、Python v2は
`eta=log(sigma)` を使うため、`zeta=2*eta`。同じ他パラメータ座標ならスコアは
`d ell / d eta = 2 * d ell / d zeta` となり、Hessianにも変換が必要である。
mfrmrの感度出力のSD標準誤差は、共同共分散の該当成分から
`0.5*sigma*sqrt(Var(zeta))` を計算する診断値である。正則化状態も記録され、公開SE資格を
上書きしない。Pythonで撤去した「他パラメータ固定の曲率をProfile SEと呼ぶ」計算とは区別する。

Python側の次の作業は、既存v2・R照合・感度アダプターを再利用し、次の順に限定する。

1. **同じ有限求積目的関数での停留性**：構造効果とSDを同時に検査する。
   SD依存の移動点を含めて微分し、EMの終了だけでは合格にしない。
2. **積分精度**：既存Q31/Q61開発確認を出発点に、同じデータで尤度・確率・EAP・事後SD・
   母集団SDの変化を調べる。TAMは範囲と密度、GHは次数と微小重みを記録する。
   高次数の相互一致だけで十分とせず、必要な条件では監査した連続積分も使う。
3. **独立実装と推論**：同定・切片・平均・分散・ステップ・カテゴリー支持・重み・
   尤度定数を対応させる。EAP同士を比較し、WLEとは分ける。SE・区間・被覆率は別ゲートに残す。
4. **高速化とUI**：上記を保存した上でJAXを評価する。画面では「推定の終了」「積分点を
   変えた結果の動き」「推論に使える範囲」を順に示し、推定と得点計算の点数を区別する。

前回のPython RSM/PCMのQ31/Q61検査は、ゆるい開発用条件での確認であり、今回のmfrmrの
証拠を引用して本番資格に昇格させない。JAXの高速化は高精度の積分を実行しやすくする候補であり、
求積誤差を自動的に消すものではない。NumPyroの導入もこのMML検査の前提ではない。
