PYTHON ?= python3
PORT ?= 8501

.PHONY: compile doctor self-test apptest benchmark parity verify run clean

compile:
	$(PYTHON) -m py_compile streamlit_app.py

doctor:
	$(PYTHON) streamlit_app.py --doctor

self-test:
	$(PYTHON) streamlit_app.py --self-test

apptest:
	$(PYTHON) -c 'from streamlit.testing.v1 import AppTest; at = AppTest.from_file("streamlit_app.py").run(timeout=30); assert not at.exception; assert any("MFRM FACETS-mode" in title.value for title in at.title); assert any("standalone Python runtime" in caption.value for caption in at.caption); assert any("Data privacy" in warning.value for warning in at.warning); print("AppTest smoke passed")'

benchmark:
	$(PYTHON) streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv

parity:
	$(PYTHON) streamlit_app.py --export-parity-fixture validation/generated/parity_fixture

verify: compile doctor self-test apptest benchmark parity

run:
	$(PYTHON) -m streamlit run streamlit_app.py --server.port $(PORT)

clean:
	rm -rf __pycache__ .pytest_cache validation/generated
