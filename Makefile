PYTHON ?= python3
PYTEST ?= $(PYTHON) -m pytest
PORT ?= 8501

.PHONY: compile doctor release-check self-test apptest benchmark demo parity verify run clean

compile:
	$(PYTHON) -m py_compile streamlit_app.py

doctor:
	$(PYTHON) streamlit_app.py --doctor

release-check:
	$(PYTHON) streamlit_app.py --release-check

self-test:
	$(PYTHON) streamlit_app.py --self-test

apptest:
	$(PYTEST) tests

benchmark:
	$(PYTHON) streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv

demo:
	$(PYTHON) streamlit_app.py --export-demo-report validation/generated/demo_report

parity:
	$(PYTHON) streamlit_app.py --export-parity-fixture validation/generated/parity_fixture

verify: compile doctor release-check self-test apptest benchmark demo parity

run:
	$(PYTHON) -m streamlit run streamlit_app.py --server.port $(PORT)

clean:
	rm -rf .pytest_cache validation/generated
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
