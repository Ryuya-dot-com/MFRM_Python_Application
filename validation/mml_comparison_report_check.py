#!/usr/bin/env python3
"""Build reviewable reports from eight hash-verified retained native comparisons."""
import argparse
import hashlib
from html import escape
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

from mfrm_app.evidence import canonical_json
from mfrm_app.exports import build_mixed_asset_zip, to_csv_bytes, to_html_report
from mfrm_app.mml_comparison_report import build_comparison_report, comparison_frames, comparison_html, STATUS
from mml_pcm_continuous_refit import dump, sha
from mml_pcm_conquest_mc import rows
from mml_regression_history_timing import history_point, app_point
from mml_unit_weight_conquest_mc import PRIOR, verify_prior, verify_hashes

TIMING = ROOT/'validation/generated/mml_regression_history_timing_20260914'
TIMING_SHA = '5065621df66180145ef78823f72d6e294d009fea77aad8c7406dfdec43ade717'
TARGETS = dict(gradient=1e-4, nll_per_person=1e-6, eap=1e-4, posterior_sd=1e-4)


def inputs():
    native, native_spec = verify_prior()
    assert sha(TIMING/'summary.json') == TIMING_SHA
    timing = json.loads((TIMING/'summary.json').read_text())
    spec = json.loads((TIMING/'protocol.json').read_text())
    assert sha(TIMING/'protocol.json') == timing['protocol_sha256']
    verify_hashes(TIMING, timing['artifact_sha256']); verify_hashes(ROOT, spec['source_sha256'])
    assert timing['all_before_compatible'] and timing['arithmetic_controls_pass']
    return native, native_spec, timing, spec


def make_report(cell, bound, native, native_spec, timing):
    key = f'{cell}_b{bound}'; run = timing['results'][key]
    pair = native['results'][cell][f'b{bound}_refit']
    end = pair['endpoints']['cq']; value = run['returned']; continuous = run['retained_continuous']
    assert np.array_equal(run['returned_app'], end['coordinates'])
    assert run['returned_equals_last'] and run['selected_iteration'] == run['rows']
    title = ('通常' if cell.startswith('ordinary') else '広い裾')+'・'+('完全' if 'complete' in cell else '欠測')+f'：TAM／ConQuest比較（±{bound}）'
    def metric(name, value, target=None, unit='logit'):
        return dict(name=name, value=float(value), upper_target=target, unit=unit)
    def check(scope, basis, metrics, missing=()):
        return dict(scope=scope, target_basis=basis, metrics=metrics, missing_evidence=list(missing))
    target_basis = '既存の開発目安を用いた事後の数値画面検査。事前登録された同等性・推論資格の基準ではない。'
    checks = {
        'calibration':check('制約と尺度を揃えたTAM／ConQuestの返却校正値。',
            'この比較には事前に定めた同等性の許容幅がない。差の数値だけを記録する。',
            [metric('自由座標の最大差', pair['native_coordinate_difference'], unit='app coordinates incl. log(SD)'),
             metric('共通求積によるEAP最大差', pair['deterministic_score_difference']['eap']),
             metric('共通求積による事後SD最大差', pair['deterministic_score_difference']['sd'])],
            ['同等性を判断する許容幅']),
        'stationarity':check('丸められたConQuest返却座標での合計NLL勾配。内部座標の停留性は未認定。',
            target_basis,
            [metric('未正規化グリッドの勾配最大値', max(abs(np.array(value['raw_app_gradient']))), TARGETS['gradient'], 'NLL / app coordinate'),
             metric('正規化済みグリッドの勾配最大値', max(abs(np.array(value['app_gradient']))), TARGETS['gradient'], 'NLL / app coordinate')],
            ['出力丸めが勾配へ及ぼす誤差の評価']),
        'integration':check('同じConQuest返却座標で有限グリッドと連続参照を比較する。', target_basis,
            [metric('未正規化NLL差／32人', abs(value['raw_nll']-continuous['nll'])/32, TARGETS['nll_per_person'], 'NLL / person'),
             metric('正規化済みNLL差／32人', abs(value['nll']-continuous['nll'])/32, TARGETS['nll_per_person'], 'NLL / person'),
             metric('EAPの積分差', max(abs(np.array(value['eap'])-continuous['eap'])), TARGETS['eap']),
             metric('事後SDの積分差', max(abs(np.array(value['sd'])-continuous['sd'])), TARGETS['posterior_sd'])]),
        'native_scores':check('返却校正値に対するConQuestの標準得点1回と、決定的求積による再計算。',
            '1 seedの差から乱数変動の安定性は判定しない。別の固定校正値での多seed検査を流用しない。',
            [metric('標準EAPと再計算の最大差', end['native_score_difference']['eap']),
             metric('標準事後SDと再計算の最大差', end['native_score_difference']['sd'])],
            ['今回の返却校正値を固定した複数seedでの得点検査', '用途に応じた標準得点の精度目安']),
    }
    folder = PRIOR/cell/f'cq_b{bound}_refit'; history = rows(folder/'history.csv')
    previous = app_point(history_point(history[-2]))
    points = dict(initial=app_point(np.array(run['initial_native'])).tolist(),
        before_returned=previous.tolist(), returned=run['returned_app'], last=end['last_coordinates'],
        tam_returned=pair['endpoints']['tam']['coordinates'])
    sources = [TIMING/'summary.json', PRIOR/'summary.json', PRIOR/cell/'input.json',
               folder/'history.csv', folder/'model.cqc', folder/'review.txt']
    details = dict(case_id=key, response_data='Two paired response datasets with complete/missing designs; 32 persons, unit response and person weights.',
        source_artifacts=[dict(path=str(p.relative_to(ROOT)), sha256=sha(p)) for p in sources],
        model='PCM, regression mean beta*x, free variance, intercept fixed at zero by CASES',
        coordinate_mapping=native_spec['mapping'], coordinate_order=['beta','rater2','task1','criterion1','criterion2','c1step1','c1step2','c2step1','c2step2','log(SD)'],
        versions=dict(TAM=native_spec['tam_version'], ConQuest=native_spec['conquest_version']),
        calibration_quadrature=dict(method='fixed physical theta grid', bound=bound, spacing=.1, points=bound*20+1),
        scoring=dict(method='native posterior Monte Carlo', p_nodes=2000, seed=2, native_score_replicates=1),
        points=points, point_digest={name:hashlib.sha256(canonical_json(p).encode()).hexdigest() for name,p in points.items()},
        stopping=dict(termination=run['termination'], executed_rows=run['rows'], selected_iteration=run['selected_iteration'], returned_equals_last=True),
        displayed_objectives=[dict(source='selected history row', nll=run['selected']['reported_nll'], definition='normalized finite prior',
                compatible_point='before_returned', timing='Compatible within export rounding; same-row point may also be compatible near convergence.'),
            dict(source='final native report', nll=run['final_display_nll'], definition='normalized finite prior',
                compatible_point='before_returned', timing='Matches selected history display within their distinct display precisions.')],
        recomputed_objectives=[dict(definition='raw density times spacing', nll=value['raw_nll'], point='returned', gradient=value['raw_app_gradient']),
            dict(definition='normalized finite prior', nll=value['nll'], point='returned', gradient=value['app_gradient']),
            dict(definition='continuous normal integral', nll=continuous['nll'], point='returned', gradient=continuous['gradient'])],
        comparison_targets=TARGETS, target_status='Post-hoc engineering screen; upper targets use strict <.',
        rounding=dict(native_xsi_beta_variance_half_width=5e-7, internal_gradient_precision_qualified=False),
        continuous_reference=dict(settings=native_spec['continuous'], numeric_relative_mass_error_sum=continuous['numeric_relative_mass_error_sum'],
            tail_relative_mass_bound_sum=continuous['tail_relative_mass_bound_sum']),
        separate_continuous_refit_comparison=dict(point='returned', reference='retained continuous local refit',
            total_eap_error=end['total_score_error']['eap'], total_posterior_sd_error=end['total_score_error']['sd'],
            structural_error=end['structural'], log_sigma_error=end['log_sigma'], continuous_nll_loss_per_person=end['continuous_nll_loss_per_person']),
        inference='No global optimum, response-weight equivalence, SE/CI/coverage or scientific inference qualification.')
    return build_comparison_report(title=title, checks=checks, details=details)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--audit', action='store_true')
    args = parser.parse_args(); out = args.output_dir.resolve()
    native, native_spec, timing, timing_spec = inputs()
    if args.audit:
        saved = json.loads((out/'summary.json').read_text()); spec = json.loads((out/'protocol.json').read_text())
        assert sha(out/'protocol.json') == saved['protocol_sha256']
        verify_hashes(ROOT, spec['source_sha256']); verify_hashes(out, saved['artifact_sha256'])
        for cell in native['results']:
            for bound in (12,20):
                name = f'{cell}_b{bound}'
                assert make_report(cell,bound,native,native_spec,timing) == json.loads((out/f'{name}.json').read_text())
        print('Eight reports replayed exactly; sources and',len(saved['artifact_sha256']),'artifacts verified')
        return
    sources = dict(timing_spec['source_sha256'])
    for p in [Path(__file__).resolve(), ROOT/'mfrm_app/mml_comparison_report.py', ROOT/'mfrm_app/evidence.py', ROOT/'mfrm_app/exports.py', ROOT/'tests/test_mml_comparison_report.py']:
        sources[str(p.relative_to(ROOT))] = sha(p)
    out.mkdir(parents=True, exist_ok=False)
    dump(out/'protocol.json', dict(classification='POST_HOC_REPORTING_CHECK', scientific_inference_ready=False,
        qualification_eligible=False, source_sha256=sources, timing_summary_sha256=TIMING_SHA,
        native_summary_sha256=sha(PRIOR/'summary.json'), targets=TARGETS,
        question='Can four independent comparison answers preserve the wide-range failures and missing evidence without emitting an overall equivalence result?',
        target_basis='Explicit post-hoc engineering targets, not prospective equivalence margins. No cross-engine margin or single-seed stability threshold supplied.',
        scope='Eight retained native comparisons, no new fitting or app UI changes. Initial report language is Japanese.'))
    overview, assets = [], {}
    for cell in native['results']:
        for bound in (12,20):
            name = f'{cell}_b{bound}'; report = make_report(cell,bound,native,native_spec,timing)
            expected_integration = 'exceeds_targets' if cell.startswith('wide') and bound==12 else 'within_targets'
            assert report['checks']['integration']['status'] == expected_integration
            assert report['checks']['stationarity']['status'] == ('exceeds_targets' if expected_integration=='exceeds_targets' else 'not_assessed')
            assert report['checks']['calibration']['status'] == report['checks']['native_scores']['status'] == 'not_assessed'
            dump(out/f'{name}.json', report)
            html = comparison_html(report); (out/f'{name}.html').write_bytes(html); assets[f'{name}.html'] = html
            assets[f'{name}.json'] = (out/f'{name}.json').read_bytes()
            for key, frame in comparison_frames(report).items():
                payload = to_csv_bytes(frame); filename = f'{name}_{key}.csv'
                (out/filename).write_bytes(payload); assets[filename] = payload
            overview.append(dict(条件=report['title'], **{key:STATUS[c['status']][0] for key,c in report['checks'].items()}))
    table = pd.DataFrame(overview).rename(columns=dict(calibration='校正値',stationarity='停留性',integration='積分精度',native_scores='標準得点'))
    (out/'overview.csv').write_bytes(to_csv_bytes(table))
    html = to_html_report({'4つの答え':table}, 'TAM／ConQuest 数値比較').decode()
    links = '<ul>'+''.join(f'<li><a href="{name}.html">{escape(row["条件"])}</a></li>' for name,row in zip([f'{c}_b{b}' for c in native['results'] for b in (12,20)],overview))+'</ul>'
    html = html.replace('<html>', '<html lang="ja">').replace('<body>', '<body><p>既存データの開発検証です。目安内は記載した数値検査の範囲に限ります。</p>').replace('</body>', '<h2>条件ごとの詳細</h2>'+links+'</body>')
    (out/'index.html').write_text(html); assets['index.html']=html; assets['overview.csv']=to_csv_bytes(table)
    (out/'comparison_reports.zip').write_bytes(build_mixed_asset_zip(assets))
    inputs(); verify_hashes(ROOT, sources)
    dump(out/'summary.json', dict(classification='POST_HOC_REPORTING_CHECK', scientific_inference_ready=False,
        qualification_eligible=False, protocol_sha256=sha(out/'protocol.json'), reports=8,
        checks=dict(wide_b12_failures_preserved=True, missing_evidence_not_promoted=True, no_overall_equivalence=True),
        artifact_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()}))
    print('Saved eight reports:', out/'index.html')


if __name__ == '__main__':
    main()
