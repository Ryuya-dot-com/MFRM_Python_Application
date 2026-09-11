# 既知の割付機構：厳密oracle資格試験

実施日: 2026-08-11  
判定: **PASS（割付生成層のみ）**

## 目的

情報的割付の影響を推定法間で点検する前に、割付生成器自体が意図した確率分布を生成するかを分離して検証した。対象は、Person次数、Rater次数、direct Rater overlapの連結性を固定した二部グラフ上の条件付き指数型分布である。

\[
P(G\mid d_P,d_R,\text{connected}) \propto
\exp\{\gamma T(G)\},\qquad
T(G)=\sum_{(p,r)\in G}z(\theta_p)z(\delta_r).
\]

`gamma` は既知の生成条件であり、観測された割付やScoreから推定しない。`gamma=0` は条件付き状態空間上の一様分布、正値は能力とRater severityの整列、負値は反整列を優先する。

## 凍結設計

- 6 Persons、3 Raters、各Person次数2、各Rater次数4
- direct Rater overlapを連結に限定
- `gamma = -0.8, 0, +0.8`
- 可能な90状態を全列挙して正規化
- 各gammaについて20,000標本の独立Metropolis 2-switch連鎖を2本
- burn-in 5,000、thin 3
- Scoreを含まない行materializationを各連鎖終点で監査

計画、seed、閾値は実行前に [known_assignment_mechanism_pilot_plan_20260811.json](known_assignment_mechanism_pilot_plan_20260811.json) に固定した。最初の起動は統計処理前にimport pathで停止し、計画や出力を変更せずrunnerの起動経路だけを修正した。

## 結果

| gamma | 厳密期待統計量 | 厳密期待割付相関 | 2連鎖pool TV距離 | pool平均統計量誤差 | oracle外確率 |
|---:|---:|---:|---:|---:|---:|
| -0.8 | -3.95466 | -0.32955 | 0.03280 | -0.00930 | 0 |
| 0.0 | 約0 | 約0 | 0.04006 | -0.02897 | 0 |
| +0.8 | +3.95466 | +0.32955 | 0.03483 | -0.01840 | 0 |

- 全gammaで厳密状態数は90、確率和は浮動小数誤差の範囲で1だった。
- `gamma=0` の90状態は数値精度内で一様だった。
- 単独連鎖の最大TV距離は0.05443（登録上限0.10）。
- 2連鎖poolの最大TV距離は0.04006（登録上限0.07）。
- 最小statistic ESSは739.87（登録下限500）。
- 全120,000標本がoracle状態空間内にあり、次数、連結性、統計量更新残差の全ゲートを通過した。
- 実データ行へのmaterializationは、Person/Rater次数、Rater別行数、Rater×Task×Criterionセル度数、連結性、source-row identityを保持した。
- Score列を任意に変更するnegative controlでも、生成designとassignment mapは完全に不変だった。

機械可読判定は [assessment.json](known_assignment_mechanism_pilot_20260811/assessment.json)、環境とSHA-256は [manifest.json](known_assignment_mechanism_pilot_20260811/manifest.json) に保持した。

## FACETSとの境界

この段階では応答がまだ存在しないため、FACETSを呼ばない。したがってPASSはFACETS parity、JMLE/MML/CMLEのbias、または実データのmissingness mechanismを証明しない。

次段階では、この合格済み割付機構からScoreを参照せずにグラフを生成し、その後で同一の明示的PCM DGMから応答を生成する。FACETS 4.5 JMLEとPython JMLEは同一estimandの外部engine校正として比較し、MMLとexact CMLEは異なるPerson処理のsensitivityとして分離する。joint／marginal／conditional likelihood、AIC、BICの横断ランキングは行わない。

FACETSの通常表示で小数第2位に丸められたfitをraw計算へ戻さない。measure/thresholdの高精度出力と、丸められたfit表示の証拠層は引き続き分離する。

## 大局的判断

従来の決定論的quartile割付は「強いストレス条件」の再現性には有用だが、機構強度を連続的に変えたdose-responseや生成確率を明示できない。今回の指数型機構はその欠点を補う。一方、実データから`gamma`を同定する方法ではないため、当面はrepository-onlyのsimulation infrastructureとし、UIへ昇格させない。
