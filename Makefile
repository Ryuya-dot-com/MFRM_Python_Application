PYTHON ?= python3
PYTEST ?= $(PYTHON) -m pytest
PORT ?= 8501

.PHONY: compile doctor self-test apptest benchmark parity verify run clean

compile:
	$(PYTHON) -m py_compile streamlit_app.py

doctor:
	$(PYTHON) streamlit_app.py --doctor

self-test:
	$(PYTHON) streamlit_app.py --self-test

apptest:
	$(PYTEST) tests/test_app_smoke.py

benchmark:
	$(PYTHON) streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv

parity:
	$(PYTHON) streamlit_app.py --export-parity-fixture validation/generated/parity_fixture

verify: compile doctor self-test apptest benchmark parity

run:
	$(PYTHON) -m streamlit run streamlit_app.py --server.port $(PORT)

clean:
	rm -rf .pytest_cache validation/generated
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
