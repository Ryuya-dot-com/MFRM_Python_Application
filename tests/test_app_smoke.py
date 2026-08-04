from streamlit.testing.v1 import AppTest


def test_app_initial_render():
    at = AppTest.from_file("streamlit_app.py").run(timeout=30)
    assert not at.exception
    assert any("MFRM FACETS-mode" in title.value for title in at.title)
    assert any("standalone Python runtime" in caption.value for caption in at.caption)
    assert any("selected built-in or generated dataset is synthetic" in caption.value for caption in at.caption)
    assert not any("Data privacy" in warning.value for warning in at.warning)
    assert not any("Keyboard shortcuts" in str(expander.label) for expander in at.expander)


def test_user_data_source_promotes_privacy_warning():
    at = AppTest.from_file("streamlit_app.py").run(timeout=30)
    at.button(key="onboarding_dismiss").click()
    at.run(timeout=30)

    assert not at.exception
    assert at.session_state["data_source_flat"] == "paste"
    assert any("pasted and uploaded rating files" in warning.value for warning in at.warning)
