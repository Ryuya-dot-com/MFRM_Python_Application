# FACETS 4.5.0 対照 PCM 境界パイロット

実施日: 2026-08-11  
親証拠: `facets_pcm_pipeline_pilot20_20260811/`

## 結論

完全均衡データで確認した FACETS–Python PCM/JMLE の一致は、今回のリング型計画欠測でも維持された。一方、意図的に Person–Rater を2成分へ分断した条件では、FACETS 4.5.0 は収束して `Subset connection O.K.` と報告したが、独立に展開した隣接カテゴリ設計とアプリ内監査はいずれも nullity 1 を検出した。

したがって、FACETS の収束や当該文言を、アプリの全予測子設計に対する識別保証として代用してはならない。アプリの `InferenceReady` を、収束とは独立した必須ゲートとして維持する。

## 事前登録した分離

同じ人数とほぼ同じ観測行数でも、観測構造が異なれば推論可能性は異なる。

- `complete`: 80 Persons × 4 Raters × 3 Tasks × 2 Criteria = 1,920行、Person–Rater 1成分
- `planned_connected`: 各Personを隣接2 Ratersが担当するリング、960行、1成分
- `disconnected_negative_control`: Persons前半をR01/R02、後半をR03/R04だけが担当、960行、2成分

各設計で、共有閾値と基準別異質閾値を同一乱数で生成し、FACETS/Pythonの両方でPCMとRSMを当てた。合計12生成run、24 fitである。

異質PCM生成データにRSMを当てた場合、単一RSM閾値の「真値」は定義しなかった。その値は観測設計とPerson分布に依存するpseudo-true量であり、基準別真値の単純平均を真値と呼ぶことは避けた。

## 実機結果

### エンジン一致

- 24/24 fitが計算を返した
- rank-full条件16/16が直接比較適格
- rank-full条件16/16が事前許容差に合格
- FACETS Table 8の尺度×カテゴリ度数144/144が入力と完全一致
- 主効果144個: FACETS対Python weighted MAE **0.002036**、最大差 **0.011621** logits
- 閾値72個: weighted MAE **0.003472**、最大差 **0.009261** logits

計画欠測では完全データより差がやや増えたが、事前許容差内だった。

| 設計 | Model | 最大run内主効果MAE | 最大主効果差 | 最大run内閾値MAE | 最大閾値差 |
|---|---:|---:|---:|---:|---:|
| complete | PCM | 0.001283 | 0.002212 | 0.003250 | 0.005481 |
| complete | RSM | 0.001487 | 0.002716 | 0.003748 | 0.006022 |
| planned connected | PCM | 0.004661 | 0.011567 | 0.005587 | 0.009261 |
| planned connected | RSM | 0.004568 | 0.011621 | 0.005874 | 0.008911 |

### モデル誤指定の方向

rank-fullの4つのdesign×replicateペア全てで、PCMのRSMに対する1観測当たりlog-likelihood改善は、共有閾値条件より異質閾値条件で大きかった。

| 設計 | 閾値条件 | PCM−RSM log-likelihood/obs 平均 | RSM−PCM AIC 平均 | RSM−PCM BIC 平均 |
|---|---:|---:|---:|---:|
| complete | shared | 0.000351 | -2.653 | -13.773 |
| complete | heterogeneous | 0.030638 | 113.650 | 102.530 |
| planned connected | shared | 0.000259 | -3.503 | -13.237 |
| planned connected | heterogeneous | 0.026657 | 47.181 | 37.447 |

共有条件では追加PCM閾値の改善はほぼ0で、AIC/BICは簡潔なRSMを選んだ。異質条件ではPCM改善が大きく、AIC/BICもPCMを選んだ。ただし2反復なので、選択確率や誤選択率の証拠ではない。

### 非連結負の対照

- 独立展開設計 nullity 1: 8/8 fit
- アプリ内 `EtaStructuralNullity=1`: 8/8 fit
- `PythonInferenceReady=false`: 8/8 fit
- 直接比較から除外: 8/8 fit
- FACETS収束: 8/8 fit
- FACETS `Subset connection O.K.`: 8/8 fit

この非連結構造では、第一群のPersonとRaterを同じ量だけ移動し、第二群を逆方向へ移動しても、全観測の `Person − Rater − Task − Criterion` は変化しない。全Raterの和を0とする制約も保てるため、尤度が識別できない自由方向が1本残る。

FACETS公式のconnectivity graph資料は、非連結subsetを接続データまたはgroup anchorで解消することを示している。今回の結果はその原則と整合するが、Table 3の短い文言だけではアプリが必要とする全設計rankを判定できなかった。

## メタ認知上の含意

1. 「FACETSが数値を返した」は、再現可能性の必要条件だが推論可能性の十分条件ではない。
2. 欠測率50%そのものが問題なのではなく、欠測パターンが連結性を保つかが重要である。
3. エンジン一致とモデル妥当性は別である。FACETSとPythonは誤指定RSMでもよく一致しうる。
4. AIC/BICやlog-likelihood差はモデル誤指定の警報になりうるが、小規模pilotを選択ルールへ直結させない。
5. 今後のJMLE・MML・CMLE比較でも、構造的不識別runを失敗分母に残しつつ、bias/RMSEの分母からは除外する必要がある。

## 次の段階

次は、この境界設計を20反復へ拡張する前に、同一rank-fullデータをJMLE・MML・exact CMLEへ渡す小規模estimand bridgeを作る。比較対象を次のように分ける。

- FACETS対Python JMLE: 同一推定量の実装一致
- JMLE対MML: fixed-person対母集団分布仮定を含む異なるestimandの感度
- JMLE対exact CMLE: Personを条件消去したfacet推定との差
- RSM対PCM: response-model誤指定の感度

各推定法を勝者として順位付けするのではなく、どの設計・仮定でどの方向の誤差が生じるかを記録する。

## 成果物

- `facets_pcm_boundary_pilot_plan_20260811.json`
- `facets_pcm_boundary_pilot.py`
- `facets_pcm_boundary_pilot_20260811/boundary_metrics.json`
- `facets_pcm_boundary_pilot_20260811/boundary_fit_ledger.csv`
- `facets_pcm_boundary_pilot_20260811/boundary_model_selection.csv`
- `facets_pcm_boundary_pilot_20260811/boundary_direction_checks.csv`
- `facets_pcm_boundary_pilot_assessment_20260811.json`
- `tests/test_facets_pcm_boundary_pilot.py`

FACETS PCM model: <https://www.winsteps.com/facetman64/raschmodels.htm>  
FACETS model control characters: <https://www.winsteps.com/facetman64/models.htm>  
FACETS Table 8 per-scale output: <https://www.winsteps.com/facetman/table8_1ratingscale.htm>  
FACETS connectivity graph: <https://www.ww.winsteps.com/facetman64/connectivitygraph.htm>
