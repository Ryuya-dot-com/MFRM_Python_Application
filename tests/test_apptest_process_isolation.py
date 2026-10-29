"""An AppTest must not replace the entrypoint used by subsequent spawn tests."""
from pathlib import Path
import subprocess
import sys


def test_apptest_then_spawn_in_one_pytest_process(tmp_path):
    # Use a separate pytest process so this regression controls the ordering
    # and verifies the real conftest teardown, rather than calling it by hand.
    (tmp_path / "conftest.py").write_bytes(Path(__file__).with_name("conftest.py").read_bytes())
    (tmp_path / "test_sequence.py").write_text(
        '''import multiprocessing as mp
import os
from streamlit.testing.v1 import AppTest


def view(label):
    import streamlit as st
    st.write(label)


def exit_seven():
    os._exit(7)


def test_01_streamlit_view():
    app = AppTest.from_function(view, args=("ready",)).run()
    assert not app.exception
    assert app.markdown[0].value == "ready"


def test_02_spawn_uses_original_entrypoint():
    process = mp.get_context("spawn").Process(target=exit_seven)
    process.start()
    try:
        process.join(timeout=30)
        assert not process.is_alive()
        assert process.exitcode == 7
    finally:
        if process.is_alive():
            process.terminate()
            process.join()
        process.close()
''',
        encoding="utf-8",
    )
    completed = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "test_sequence.py"],
        cwd=tmp_path, capture_output=True, text=True, timeout=90,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
