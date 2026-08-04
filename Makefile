PYTHON ?= python3
PYTEST ?= $(PYTHON) -m pytest
NATIVE_PYTEST_ARGS ?= -m "not legacy_compat" --ignore=tests/test_cross_engine_bundle.py --deselect=tests/test_classical_dif.py::test_validation_bundle_contents
PORT ?= 8501

.PHONY: compile doctor release-check self-test apptest benchmark demo verify run clean

compile:
	$(PYTHON) -m py_compile streamlit_app.py

doctor:
	$(PYTHON) streamlit_app.py --doctor

release-check:
	$(PYTHON) streamlit_app.py --release-check

self-test:
	$(PYTHON) streamlit_app.py --self-test

apptest:
	$(PYTEST) $(NATIVE_PYTEST_ARGS) tests

benchmark:
	$(PYTHON) streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv

demo:
	$(PYTHON) streamlit_app.py --export-demo-report validation/generated/demo_report

verify: compile doctor release-check self-test apptest benchmark demo

run:
	$(PYTHON) -m streamlit run streamlit_app.py --server.port $(PORT)

clean:
	rm -rf .pytest_cache validation/generated
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
