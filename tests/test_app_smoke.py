from streamlit.testing.v1 import AppTest


def test_app_initial_render():
    at = AppTest.from_file("streamlit_app.py").run(timeout=30)
    assert not at.exception
    assert any("MFRM estimation" in title.value for title in at.title)
    assert any("standalone Python runtime" in caption.value for caption in at.caption)
    assert any("Data privacy" in warning.value for warning in at.warning)
