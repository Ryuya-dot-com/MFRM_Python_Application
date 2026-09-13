# MFRM estimand bridge pilot

実施日: 2026-08-11  
親データ: `facets_pcm_boundary_pilot_20260811/retained_input/`

## 結論

同一のrank-fullデータを、FACETS JMLE、Python JMLE、Python MML（Person SD固定／自由）、native exact CMLEへ接続するbridgeが成立した。

- 親FACETS/Python JMLE契約32件をhash固定で再利用
- 新規MML 32/32 fitが収束・`InferenceReady=true`
- 新規exact CMLE 16/16 fitが条件付きrank full・有限解存在・収束・`InferenceReady=true`
- main-facet recovery 720行、threshold recovery 360行を同一RunIdで結合
- 制約残差最大 `1.11e-16`
- joint／marginal／conditional likelihoodごとのRSM対PCM方向判定16/16合格

これは推定法の順位表ではない。三つの尤度はPersonの扱いが異なるため、JMLE対MML対CMLEの差は、実装誤差ではなくestimandと仮定への感度として読む。

## Estimand契約

| Mode | Personの扱い | likelihood basis | 直接比較できる範囲 |
|---|---|---|---|
| FACETS / Python JMLE | 各Personを固定効果として同時推定 | joint | FACETS対Pythonのみ同一estimandのengine比較 |
| MML fixed SD=1, Q31 | 正規母集団、SD=1固定 | marginal | 同じMML mode内のRSM対PCM |
| MML free SD, Q31 | 正規母集団、SDを推定 | marginal | 同じMML mode内のRSM対PCM |
| exact CMLE | Person総得点で条件化してPersonを消去 | conditional | CMLE内のRSM対PCM |

CMLEのconditional log-likelihood/AICを、JMLEやMMLの値と数値比較していない。

## モデル誤指定方向

全likelihood basisで、共有閾値ではPCMの追加改善がほぼ0、異質閾値ではPCMの改善が大きかった。

| Mode | Design | shared: PCM−RSM LL/obs | heterogeneous: PCM−RSM LL/obs |
|---|---:|---:|---:|
| Python JMLE | complete | 0.000351 | 0.030638 |
| Python JMLE | planned connected | 0.000259 | 0.026657 |
| MML fixed SD=1 | complete | 0.000347 | 0.030594 |
| MML fixed SD=1 | planned connected | 0.000247 | 0.025935 |
| MML free SD | complete | 0.000342 | 0.030365 |
| MML free SD | planned connected | 0.000243 | 0.025611 |
| exact CMLE | complete | 0.000346 | 0.030258 |
| exact CMLE | planned connected | 0.000255 | 0.025778 |

16個のdesign×replicate×mode方向チェックは全て合格した。この一貫性は、異質なcriterion別閾値がRSMの共通閾値制約に適合しないことを、Personの扱いが異なる三尤度でも検出できたという実装証拠である。2反復なので、model-selection accuracyではない。

## 記述的truth recovery

正しいPCMを当てたmain facetの全条件合算RMSEは次の通りだった。

| Mode | complete | planned connected |
|---|---:|---:|
| FACETS JMLE | 0.0334 | 0.1111 |
| Python JMLE | 0.0336 | 0.1136 |
| MML fixed SD=1 | 0.0331 | 0.0956 |
| MML free SD | 0.0334 | 0.0905 |
| exact CMLE | 0.0334 | 0.0987 |

PCM閾値RMSEは次の通りだった。

| Mode | complete | planned connected |
|---|---:|---:|
| FACETS JMLE | 0.0834 | 0.1566 |
| Python JMLE | 0.0846 | 0.1615 |
| MML fixed SD=1 | 0.0790 | 0.1038 |
| MML free SD | 0.0779 | 0.0762 |
| exact CMLE | 0.0773 | 0.0743 |

計画欠測でMML/CMLEのRMSEがJMLEより小さいという観測は、次段階の仮説にはなるが、推定法の優位性を示さない。2つのseedだけであり、固定SD=1は実現Person SD約0.760–0.773とも一致していない。それでも一部RMSEが小さいため、単一pilotのRMSE最小値を「正しい仮定」と読み替える危険性がよく表れている。

## Person SD感度

自由SD MMLの推定範囲は0.740–0.872だった。

- complete: 条件・model平均で約0.751–0.776
- planned connected: 約0.805–0.842
- 実現Person SD: replicate 1 = 0.773、replicate 2 = 0.760

同じPersonを使っても、観測デザインを半分にすると自由SD推定が上方へ動いた。この変化は、欠測下の情報量、有限標本、response modelの影響が母集団分散とfacet推定へ伝播する可能性を示す。20反復以上で分布を確認すべき仮説であり、現時点でbiasとは断定しない。

## JMLEからの推定移動

主効果のPython JMLEからの平均絶対移動は、完全データより計画欠測で拡大した。

| Mode | complete | planned connected | 最大移動 |
|---|---:|---:|---:|
| FACETS JMLE | 約0.0011–0.0013 | 約0.0046 | 0.0116 |
| MML fixed SD=1 | 約0.0028–0.0030 | 約0.022–0.023 | 0.0823 |
| MML free SD | 約0.0113–0.0114 | 約0.033–0.034 | 0.1151 |
| exact CMLE | 約0.0110–0.0112 | 約0.025 | 0.0613 |

FACETS対Pythonの移動はengine差、MML/CMLEの移動はdifferent-estimand sensitivityである。同じ表に置いても意味ラベルを混ぜない。

## 計算資源

- exact CMLE 16 fit: 合計31.1秒、中央値1.91秒
- MML fixed SD=1 Q31 16 fit: 合計210.8秒、中央値12.28秒
- MML free SD Q31 16 fit: 合計172.6秒、中央値10.92秒
- bridge全体: 約418秒

大規模反復ではMMLが計算予算を支配する。次段階はmanifest分割、再開可能なfit ledger、mode別並列化を前提とする。

## 大局的な判断

このアプリの強みは「FACETSと同じ数値を出す」ことだけではない。FACETSでJMLE実装を校正したうえで、同一データを異なるPerson処理へ渡し、設計と仮定が結果をどの方向へ動かすかを可視化できる点にある。

現時点の最も価値ある次段階は、条件を増やす前に次を事前登録することである。

1. 20反復pilotでPerson SD、main facet、thresholdの分布を推定する
2. Monte Carlo SEを付け、RMSE差ではなく各estimandの条件別挙動として報告する
3. MMLについて正規分布誤指定（混合・歪み・重尾）を追加する
4. CMLEについて極端Person率と条件付き情報の低下を追加する
5. その後にのみ500反復研究とUI上の選択ガイダンスを検討する

## 成果物

- `estimand_bridge_pilot_plan_20260811.json`
- `estimand_bridge_pilot.py`
- `estimand_bridge_pilot_20260811/estimand_bridge_metrics.json`
- `estimand_bridge_pilot_20260811/estimand_bridge_run_ledger.csv`
- `estimand_bridge_pilot_20260811/estimand_bridge_recovery.csv`
- `estimand_bridge_pilot_20260811/estimand_bridge_thresholds.csv`
- `estimand_bridge_pilot_20260811/estimand_bridge_model_pairs.csv`
- `estimand_bridge_pilot_20260811/estimand_bridge_direction_checks.csv`
- `estimand_bridge_pilot_assessment_20260811.json`
- `tests/test_estimand_bridge_pilot.py`
