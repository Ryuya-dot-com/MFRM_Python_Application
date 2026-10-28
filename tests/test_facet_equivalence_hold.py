"""The advanced report must not bypass uncertainty holds with equivalence claims."""
import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from mfrm_app import output_qualification as qualification


@pytest.mark.parametrize('lang', ['en', 'ja'])
@pytest.mark.parametrize('method,free_sd', [('JMLE', False), ('MML', False), ('MML', True)])
@pytest.mark.parametrize('uncertainty', ['legacy_ready', 'zero', 'missing'])
def test_equivalence_panel_withholds_all_claims_including_legacy_values(lang, method, free_sd, uncertainty):
    def view(lang, method, free_sd, uncertainty):
        from unittest.mock import patch
        import numpy as np
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        st.session_state['lang'] = lang
        st.session_state['eq_bound'] = .5
        se = {'legacy_ready': .084, 'zero': 0., 'missing': np.nan}[uncertainty]
        result = {'config': {'method': method, 'model': 'PCM', 'estimate_population_sd': free_sd, 'facet_names': ['Rater']},
                  'prep': {'n_obs': 127}}
        diagnostics = {'measures': pd.DataFrame({'Facet': ['Rater']*2, 'Element': ['R0', 'R1'],
            'Estimate': [.012, -.012], 'SE': [se, se], 'SE_Status': ['Ready']*2, 'CI_Status': ['Ready']*2})}
        with patch.object(app, '_offer_fig_download'):
            app._render_facet_equivalence(result, diagnostics)
    at = AppTest.from_function(view, args=(lang, method, free_sd, uncertainty)).run(timeout=45)
    assert not at.exception
    assert not at.success and not at.error and not at.metric and not at.dataframe
    locale = json.loads((Path(app.__file__).parent / 'locales' / f'{lang}.json').read_text())['report_top']
    assert [w.value for w in at.warning] == [locale['equivalence_hold']]
    assert not at.get('plotly_chart') and not at.get('download_button')
    assert not at.number_input and not at.selectbox
    assert all(not item.value.startswith('report_top.') for item in [*at.caption, *at.markdown])


@pytest.mark.parametrize('lang', ['en', 'ja'])
def test_full_report_routes_to_the_held_panel(lang):
    def view(lang):
        from contextlib import ExitStack
        from unittest.mock import patch
        import streamlit as st
        import streamlit_app as app
        st.session_state['lang'] = lang
        others = ['_render_apa_report', '_render_manuscript_template_section', '_render_method_appendix_section',
                  '_render_manuscript_claim_guide_section', '_render_report_tables', '_render_reporting_checklist',
                  '_render_final_readiness_section', '_render_publication_document_section']
        with ExitStack() as stack:
            for name in others:
                stack.enter_context(patch.object(app, name))
            app.show_report_section({'config': {'method': 'MML'}}, {}, force_full=True)
    at = AppTest.from_function(view, args=(lang,)).run(timeout=45)
    assert not at.exception
    at.selectbox(key='report_detail_panel').set_value('Facet Equivalence').run()
    assert not at.exception
    locale = json.loads((Path(app.__file__).parent / 'locales' / f'{lang}.json').read_text())['report_top']
    assert [w.value for w in at.warning] == [locale['equivalence_hold']]


@pytest.mark.parametrize('method', ['JMLE', 'MML'])
def test_equivalence_hold_is_exported_with_stable_reason(method):
    item = qualification.facet_equivalence_qualification()
    assert item.status is qualification.OutputUse.WITHHELD
    assert not item.public_conclusion_allowed and not item.raw_technical_export_allowed
    row = app.build_output_qualification_table({'config': {'method': method}}).set_index('OutputID').loc['facet.equivalence']
    assert row.ReasonCode == 'stat.stat_003.facet_equivalence_unvalidated'
    assert row.QualificationStatus == 'WITHHELD' and not row.PublicConclusionAllowed
