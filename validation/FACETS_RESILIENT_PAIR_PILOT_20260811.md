# FACETS/Python resilient pair remediation pilot（2026-08-11）

## 結論

fresh-100 confirmatoryで観測されたFACETS exit-0/missing-U6-report障害に対し、凍結済みv1を変更しないversioned componentを作成した。登録した4つの故障注入scenario・13チェックは全て合格した。

最も重要な変更は、Python JMLEをFACETSより先にfitし、artifact hashを含むcompletion markerを先に書くことである。その後のFACETSが全retryを失敗しても、Python行は `StatisticalEvidenceReady=True`、`IncludedInStudy=True` のまま残り、FACETS側だけが `CalibrationReady=False` になる。外部校正engineの運用失敗が、既に成功したPython統計証拠を消すv1の結合を解いた。

## 登録契約

- retryは初回＋最大2回、合計最大3 try
- retry対象は `FACETS primary/auxiliary failed ... exit=0` かつ該当report欠落だけ
- nonzero exit、timeout、parse/count/convergence/direct-agreement/Python failureはretryしない
- 全tryのspec・stdout・stderr・report・outcome JSONを保持
- cross-process file lockでlegacy pair callを直列化
- `StatisticalEvidenceReady`、`CalibrationReady`、`PairFullyQualified`を分離
- 成功時は独立Pythonとlegacy pair内Pythonが `1e-12` 以内で一致することを要求
- `mfrm_app/*.py`、`streamlit_app.py`、直接validation依存、requirements、component、plan、FACETS.exeをhash化

## 結果

### FI01: 初回だけexit-0/missing-report

1回目の障害を保持し、2回目は実FACETS 4.5.0で成功した。Python statistical evidence、FACETS calibration、pair full qualificationは全てtrueだった。独立Pythonと再計算Pythonの最大差は、主効果・閾値とも0だった。

### FI02: retry対象障害を3回連続注入

3 try全てを保持して終了した。FACETS calibrationとpair full qualificationはfalseだが、先にcompletion-markしたPython statistical evidenceはtrueのまま残った。

### FI03: nonzero exitを注入

retryせず1 tryで終了した。ここでもPython evidenceはtrue、FACETS calibrationはfalseだった。成功するまで統計・数値障害を再試行する挙動はない。

### FI04: lock競合

4 workerの同時要求に対してlock同時保持最大は1だった。

## 同一性

- dependency files: 57
- dependency manifest SHA-256: `d1e7e8a754b8288872873faf1a621a87f85881950f555fe3b1cfd67971f7dd67`
- FACETS.exe SHA-256: `dfb0afb0faa18f026d1b3b4175f22e42cc3764430eb83cbd368c7a572b3593a1`
- FACETS reported version: 4.5.0
- 元confirmatory attempt 01024 marker: 不変
- confirmatory evidence replacement: false

## メタ認知的評価

v1の問題はFACETS数値校正ではなく、Python fitと外部reportの運命を1つのattempt exceptionへ結合したことだった。v2 componentはこれを分離した。また、v1のtop-level script hashだけではimport先変更を検出できないため、依存コードとFACETS binaryまで同一性を拡張した。

ただし、現在のlockは既存monolithic pair関数全体を囲む。このためPythonを独立fitした後、lock内でもう一度Pythonを計算する。証拠分離を優先した保守的実装であり、性能上の最終形ではない。将来はFACETS spec生成・呼出し・parseをPython fitから低レベル分離し、lock範囲をFACETS processだけへ狭めるべきである。

## 判定

別plan・fresh replicateを用いるH6-only runnerへの利用gateは合格した。一方、今回のcomponentやdiagnostic outputを過去のconfirmatoryへ戻してH6を99から100へ増やすこと、strict batch資格を合格へ変更すること、推定法ランキングへ利用することは禁止する。

## 成果物

- `facets_resilient_pair_plan_20260811.json`
- `facets_resilient_pair.py`
- `facets_resilient_pair_pilot.py`
- `facets_resilient_pair_pilot_20260811/`
- `facets_resilient_pair_pilot_assessment_20260811.json`
- `tests/test_facets_resilient_pair.py`
