# FACETS 4.5.0 Windows起動障害の訂正診断

実施日: 2026-08-11  
現在の判定: **FACETS全体の障害ではない。hidden-window起動経路とlegacy path長の二つの運用障害。**

## 訂正

以前の「FACETS 4.5.0のインストール／ランタイム全体がこのWindows環境で使用不能」という結論は撤回する。ユーザーが示した通り、FACETSの対話GUIは起動できた。その後の対照実験により、クラッシュはPython固有ではなく、FACETS/Xojoのwindowを不可視にする経路へ局在した。

保持済みProcDumpは誤りではない。`XojoGUIFramework64.dll + 0x234DD8`、`RAX=0`、read address `0x58` という事実は、hidden-window経路で発生したnative null dereferenceを示す。ただし、それをvisible/interactive FACETSにも一般化した解釈が誤りだった。

## 起動経路の実証表

| 経路 | Pythonの関与 | window状態 | 結果 |
|---|---|---|---|
| 引数なしGUI起動 | なし | visible | 起動成功 |
| direct CMD `BATCH=YES`、公式Kct | なし | FACETS hidden batch | `0xC0000005`、reportなし |
| direct CMD、spec/report後に `BATCH=YES` | なし | FACETS hidden batch | `0xC0000005`、reportなし |
| direct CMD `BATCH=NO` | なし | visible | report生成、exit 0 |
| Python `BATCH=NO` + hidden `STARTUPINFO` | あり | hidden | `0xC0000005`、reportなし |
| Python `BATCH=NO` + visible output gate | あり | visible | report/Scorefiles生成、graceful close、exit 0 |

この対照により、「PythonスクリプトがFACETSへ不適切なhidden起動を要求していた」という説明が成立する。同時に、direct CMD `BATCH=YES`も失敗したため、Pythonインタプリタ自体がnative crashを起こしたわけではない。

## qualified visible adapter

現在の共通呼出しはWindowsで次を既定とする。

1. `BATCH=NO`でwindowをvisibleのまま起動する。
2. reportと登録した全facet Scorefileが存在し、非空で、size/mtimeが2秒以上安定するまで待つ。
3. FACETSプロセス所有のtop-level windowへ `WM_CLOSE` を送る。
4. exit 0、全output継続存在、強制終了なしを成功条件にする。
5. `BATCH=YES`は、別hostで独立qualificationした場合だけ明示的に選べる。

公式 `Kct.txt` と過去の4-facet project controlの両方がこのshared helperで合格した。過去project controlではreportと4 Scorefileがすべて生成された。

## 第二の障害: 260文字path境界

visible adapterによる補完実行の最初のdatasetでは、primary `Umean=6` passが成功したが、auxiliary `Umean=2` passがTable 5後に待機した。

- primary `scores.4.txt`: 絶対パス257文字、生成成功
- auxiliary `scores_u2.4.txt`: 絶対パス260文字、未生成
- 同一spec/dataの短い診断path: 最長98文字、5.3秒でreportと4 Scorefileを生成

これはモデルや丸めではなく、legacy Windows path境界である。共通呼出し層は現在、spec、report、推定される全Scorefileについて220文字を超える場合、FACETS起動前に明示的に失敗する。補完studyは最悪219文字の短い保存layoutへ登録し直し、12/12を完了した。

## 科学的境界

- hidden crashはJMLE計算やFACETSの2桁丸めに到達する前の運用障害である。
- path stallも推定結果の差ではない。
- FACETSの通常fit表示は小数第2位までであり、未報告raw fitを再構成していない。
- measure/thresholdの校正は登録済み `Umean=6` primary passを用いる。
- visible補完はFACETS/Pythonの同一estimand JMLE校正であり、MMLやCMLEの正解判定ではない。

## 保持する診断証拠

クラッシュdump、repair、MINIFAC 4.5.1対照、vendor packetは削除しない。これらはhidden-window経路が共有Xojo GUI runtime内で同じnull dereferenceを起こす証拠として有用である。ただし、今後は「FACETS 4.5.0全体が起動不能」という主張の証拠には使わない。

## 主要artifact

- [direct BATCH=YES plan](facets_direct_cmd_batchyes_plan_20260811.json)
- [direct BATCH=NO plan](facets_direct_cmd_batchno_plan_20260811.json)
- [visible launcher qualification](facets_visible_launcher_qualification_v2_20260811/result.json)
- [shared helper qualification](facets_shared_helper_qualification_20260811/result.json)
- [short-path auxiliary probe](fx2_20260811/result.json)
- [path amendment](known_assignment_facets_visible_supplement_path_amendment_20260811.json)
- [12/12 supplemental assessment](known_assignment_multivector_preflight4_20260811/v/assessment.json)
- [retained crash analysis](kafp_crash_capture_20260811/crash_analysis.json)
- [retained MINIFAC crash analysis](minifac_451_crash_capture_20260811/crash_analysis.json)
