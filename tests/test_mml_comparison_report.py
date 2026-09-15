"""Comparison failures and missing evidence must survive reporting and export."""
from copy import deepcopy
import math

import pytest

from mfrm_app.mml_comparison_report import (
    QUESTIONS, build_comparison_report, comparison_frames, comparison_html,
)


def report_inputs():
    checks = {key:dict(scope="Returned calibration only", target_basis="Explicit development target",
        metrics=[dict(name="Error", value=0., upper_target=1e-4, unit="logit")], missing_evidence=[])
        for key in QUESTIONS}
    # The observed wide +/-12 pattern: small calibration discrepancy coexists
    # with appreciable gradient/integration error and no matched-seed experiment.
    checks['calibration']['metrics'][0].update(value=1.7608567518090368e-5, upper_target=None)
    checks['calibration']['missing_evidence'] = ['Equivalence margin']
    checks['stationarity']['metrics'][0].update(value=.28924760278567907)
    checks['stationarity']['missing_evidence'] = ['Internal-coordinate rounding assessment']
    checks['integration']['metrics'][0].update(value=.08725662163698744)
    checks['native_scores']['metrics'][0].update(value=.09747691784269286, upper_target=None)
    checks['native_scores']['missing_evidence'] = ['Repeated seeds at the returned calibration']
    return dict(title='Comparison', checks=checks, details={
        'native_converged':True,
        'displayed_nll':dict(value=151.7227955, evaluation_point='before_update'),
        'recomputed_nll':dict(value=151.722795628768, evaluation_point='returned'),
    })


def test_four_answers_preserve_negative_control_and_evaluation_points():
    source = report_inputs(); original = deepcopy(source)
    report = build_comparison_report(**source)
    assert source == original
    assert {k:v['status'] for k,v in report['checks'].items()} == dict(
        calibration='not_assessed', stationarity='exceeds_targets',
        integration='exceeds_targets', native_scores='not_assessed')
    assert report['details'] == source['details']
    assert report['scientific_inference_ready'] is report['qualification_eligible'] is False
    assert 'overall' not in report and 'equivalent' not in report
    source['details']['native_converged'] = False
    assert report['details']['native_converged'] is True


def test_missing_evidence_or_measurement_does_not_pass():
    source = report_inputs()
    source['checks']['stationarity']['metrics'][0]['value'] = 0.
    source['checks']['integration']['metrics'][0]['value'] = None
    report = build_comparison_report(**source)
    assert report['checks']['stationarity']['status'] == 'not_assessed'
    assert report['checks']['integration']['status'] == 'not_assessed'
    source['checks']['integration']['metrics'] = []
    assert build_comparison_report(**source)['checks']['integration']['status'] == 'not_assessed'


def test_target_uses_unrounded_numbers_and_strict_boundary():
    source = report_inputs()
    for value, expected in [(math.nextafter(1e-4, 0.), 'within_targets'),
                            (1e-4, 'exceeds_targets'),
                            (math.nextafter(1e-4, math.inf), 'exceeds_targets')]:
        source['checks']['integration']['metrics'][0]['value'] = value
        assert build_comparison_report(**source)['checks']['integration']['status'] == expected


@pytest.mark.parametrize('value', [True, '0', -1., math.inf, math.nan])
def test_invalid_measurements_are_rejected(value):
    source = report_inputs(); source['checks']['integration']['metrics'][0]['value'] = value
    with pytest.raises(ValueError): build_comparison_report(**source)


def test_invalid_target_and_nonfinite_provenance_are_rejected():
    source = report_inputs(); source['checks']['integration']['metrics'][0]['upper_target'] = 0.
    with pytest.raises(ValueError): build_comparison_report(**source)
    source = report_inputs(); source['details']['reported_nll'] = math.nan
    with pytest.raises(ValueError): build_comparison_report(**source)
    source = report_inputs(); del source['checks']['native_scores']
    with pytest.raises(ValueError): build_comparison_report(**source)


def test_export_recomputes_decisions_and_html_escapes_provided_text():
    source = report_inputs(); source['title'] = '<script>alert(1)</script>'
    source['checks']['integration']['scope'] = '<img src=x onerror=alert(1)>'
    source['details']['note'] = '</pre><script>alert(2)</script>'
    report = build_comparison_report(**source)
    report['checks']['integration']['status'] = 'within_targets'
    report['checks']['integration']['metrics'][0]['target_met'] = True
    frames = comparison_frames(report)
    assert frames['comparison_summary'].iloc[2]['Status'] == '目安未達'
    assert frames['comparison_metrics'].iloc[2]['TargetMet'] is False or frames['comparison_metrics'].iloc[2]['TargetMet'] == False
    for language in ('ja', 'en'):
        html = comparison_html(report, language=language).decode()
        assert '<script>' not in html and '<img ' not in html
        assert '&lt;script&gt;' in html and '<details>' in html and '<summary>' in html
        assert '<details open' not in html
        assert f'lang="{language}"' in html
    report['scientific_inference_ready'] = True
    with pytest.raises(ValueError): comparison_html(report)
