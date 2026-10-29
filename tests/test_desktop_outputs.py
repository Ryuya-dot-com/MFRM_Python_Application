"""Native/export parity and qualification from saved extended analyses."""
from copy import deepcopy
import re

import numpy as np
import pandas as pd
import pytest

from mfrm_app import desktop_outputs as output, desktop_analysis as desktop, project, engine, structural_intervals as intervals
from tests.test_desktop_workflows import analysis, _forbid_calculation
from tests.test_g_d_study_boundaries import ratings


def test_all_outputs_preserve_saved_scope_and_never_refit(analysis, tmp_path, monkeypatch):
    _forbid_calculation(monkeypatch)
    choices = output.figure_choices(analysis)
    assert any(o['kind'] == 'bias' for _, o in choices)
    assert any(o['kind'] == 'gd' for _, o in choices)
    figures = [(o, output.analysis_figure(analysis, **o)) for _, o in choices]
    config = analysis['result']['config']
    for options, fig in figures:
        if options['kind'] == 'category':
            steps, _ = output._steps(analysis['result'], options['level'])
            curves = fig.axes[0].lines
            eta = curves[0].get_xdata()
            weights = np.exp(np.outer(eta, np.arange(len(steps)+1)) - np.r_[0., np.cumsum(steps)])
            expected = weights / weights.sum(axis=1, keepdims=True)
            np.testing.assert_allclose(np.array([line.get_ydata() for line in curves]).T, expected, rtol=1e-12)
        if options['kind'] == 'wright':
            steps, _ = output._steps(analysis['result'], options['level'])
            dotted = [line.get_ydata()[0] for line in fig.axes[1].lines if line.get_linestyle() == ':']
            np.testing.assert_array_equal(dotted, steps)
            person_values = analysis['result']['facets']['person']['Estimate'].to_numpy()
            if np.ptp(person_values) < 1e-9:
                assert sum(p.get_width() for p in fig.axes[0].patches) == len(person_values)
                assert all(p.get_height() >= .1 for p in fig.axes[0].patches)
        if options['kind'] == 'interval':
            table = intervals.read_structural_intervals(analysis['result'], analysis['interval_analyses'][options['index']]['record'])
            np.testing.assert_array_equal(fig.axes[0].lines[0].get_xdata(), table.loc[table.CIEligible, 'Estimate'])
        if options['kind'] == 'bias':
            raw = analysis['bias_analyses'][0]['result']['table']
            expected = raw.assign(value=raw['Bias Size'].where(raw['Cell Status'].eq('interior'))).pivot(index='FacetA_Level', columns='FacetB_Level', values='value')
            np.testing.assert_allclose(fig.axes[0].images[0].get_array().filled(np.nan), expected.to_numpy(), equal_nan=True)
            assert fig.axes[0].images[0].get_clim()[1] >= 1
    report = output.report_html(analysis, pathway_statistic='InfitZSTD', pathway_facet='Rater')
    assert config['model'] + ' / ' + config['method'] in report
    assert 'LR Prob.' not in report and 'S.E.' not in report
    assert report.count('<svg') == len(choices)
    ids = re.findall(r'\bid="([^"]+)"', report)
    assert len(ids) == len(set(ids))
    assert 'cdn.' not in report and '<script' not in report
    path = tmp_path / 'result.mfrm'; project.save_project(path, analysis)
    assert output.report_html(project.load_project(path), pathway_statistic='InfitZSTD', pathway_facet='Rater') == report
    # Routine tests must not overwrite the retained evidence from earlier runs.
    evidence = tmp_path / 'native_replay'
    evidence.mkdir(parents=True, exist_ok=True)
    stem = config['model'] + '_' + config['method']
    project.save_project(evidence / (stem + '.mfrm'), analysis)
    (evidence / (stem + '.html')).write_text(report, encoding='utf-8')
    for i, (options, fig) in enumerate(figures):
        (evidence / f'{stem}_{i}.png').write_bytes(output.basic.figure_bytes(analysis, str(i), figure=fig, format='png'))


def test_changed_interval_source_is_held_in_figures_and_reports(analysis, monkeypatch):
    if not analysis['interval_analyses']: return
    _forbid_calculation(monkeypatch)
    monkeypatch.setattr(intervals, '_source_hash', lambda: 'changed-source')
    fig = output.analysis_figure(analysis, 'interval')
    assert not fig.axes[0].lines and not fig.axes[0].collections
    assert 'Withheld' in fig.axes[0].texts[0].get_text()
    assert 'Withheld' in output.report_html(analysis)


@pytest.mark.parametrize('two_facets', [False, True])
def test_gd_output_available_unavailable_and_html_escaping(two_facets, monkeypatch):
    data = ratings(two_facets)
    settings = dict(object_facet='Person', random_facets=['Rater', 'Item'] if two_facets else ['Rater'], target_population='<script>population</script>')
    good = desktop.run_analysis(data, workflow='gd', gd_settings=settings)
    missing = desktop.run_analysis(data.iloc[1:], workflow='gd', gd_settings=settings)
    _forbid_calculation(monkeypatch)
    for a in (good, missing):
        report = output.report_html(a)
        assert '<script>' not in report and '&lt;script&gt;population' in report
        fig = output.analysis_figure(a, 'gd')
        if a['gd_analysis']['result']['available']:
            table = a['gd_analysis']['result']['d_study']
            first = table if not two_facets else table.loc[table.Item.eq(table.Item.min())]
            np.testing.assert_array_equal(fig.axes[0].lines[0].get_ydata(), first.sort_values('Rater').G)
        else:
            assert not fig.axes[0].lines and 'missing' in report.lower()


def test_unknown_pcm_target_is_rejected(analysis):
    if analysis['fit_settings']['model'] == 'PCM':
        with pytest.raises(ValueError, match='step-facet'):
            output.analysis_figure(analysis, 'category', level='foreign')


def test_malformed_view_is_rejected_before_replacing_native_state(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from tkinter import filedialog
    from desktop_app import DesktopApp
    a = desktop.run_analysis(ratings(False), workflow='gd', gd_settings=dict(object_facet='Person', random_facets=['Rater']))
    a['view_settings'] = dict(figure='D Study · G / Phi', facet='すべて', statistic='Infit MnSq')
    a['content_fingerprint'] = project.content_fingerprint(a)
    path = tmp_path / 'view.mfrm'; project.save_project(path, a)
    assert project.load_project(path)['view_settings'] == a['view_settings']
    errors = []
    ui = SimpleNamespace(tr=lambda value: value, root=None, busy=False, analysis={'existing': True}, discard=lambda: True, error=errors.append,
                         set_data=lambda *args: pytest.fail('Reject malformed view before replacing input'))
    monkeypatch.setattr(filedialog, 'askopenfilename', lambda **kwargs: str(path))
    for malformed in ('bad', {}, {'figure': []}):
        a['view_settings'] = malformed
        a['content_fingerprint'] = project.content_fingerprint(a)
        path.write_bytes(project._json_bytes(dict(format=project.FORMAT, schema_version=2, analysis=project._encode(a))))
        DesktopApp.open_project(ui)
        assert 'view settings' in str(errors[-1])
        assert ui.analysis == {'existing': True}
