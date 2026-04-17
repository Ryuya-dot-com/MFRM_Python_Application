"""Tests for the advanced-model Stan code generators.

The sidebar's "🧪 Advanced models (Stan, download only)" expander can
emit Stan programs for DINA, HRM, Testlet RI / Bifactor, Mixture
Rasch, 2PL Binary, and Pairwise BTL. These tests verify each
generator produces a syntactically-plausible Stan program (required
blocks present, brace balance, no generator-crash regressions).
"""

from __future__ import annotations

import pytest

import streamlit_app as app


def _braces_balanced(code: str) -> bool:
    depth = 0
    for ch in code:
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth < 0:
                return False
    return depth == 0


def _has_all_blocks(code: str, blocks: list[str]) -> bool:
    for b in blocks:
        if f"{b} " not in code and f"{b}{{" not in code:
            return False
    return True


# ---------------------------------------------------------------------------
# DINA
# ---------------------------------------------------------------------------

def test_dina_generator_produces_non_empty_code():
    code = app.generate_dina_stan_code(n_items=5, n_attributes=3)
    assert isinstance(code, str) and len(code) > 200


def test_dina_generator_has_standard_blocks():
    code = app.generate_dina_stan_code(n_items=5, n_attributes=3)
    assert _has_all_blocks(code, ["data", "parameters", "model"])


def test_dina_generator_has_priors_for_slip_and_guess():
    code = app.generate_dina_stan_code(n_items=5, n_attributes=3)
    assert "slip" in code
    assert "guess" in code


# ---------------------------------------------------------------------------
# HRM
# ---------------------------------------------------------------------------

def test_hrm_generator_produces_non_empty_code():
    code = app.generate_hrm_stan_code(n_categories=5)
    assert isinstance(code, str) and len(code) > 200


def test_hrm_generator_has_standard_blocks():
    code = app.generate_hrm_stan_code(n_categories=5)
    assert _has_all_blocks(code, ["data", "parameters", "model"])


# ---------------------------------------------------------------------------
# Testlet (RI and Bifactor)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bifactor", [False, True])
def test_testlet_generator_produces_non_empty_code(bifactor: bool):
    code = app.generate_testlet_stan_code(n_categories=5, bifactor=bifactor)
    assert isinstance(code, str) and len(code) > 200


@pytest.mark.parametrize("bifactor", [False, True])
def test_testlet_generator_has_standard_blocks(bifactor: bool):
    code = app.generate_testlet_stan_code(n_categories=5, bifactor=bifactor)
    assert _has_all_blocks(code, ["data", "parameters", "model"])


# ---------------------------------------------------------------------------
# Mixture Rasch
# ---------------------------------------------------------------------------

def test_mixture_rasch_generator_produces_non_empty_code():
    code = app.generate_mixture_rasch_stan_code(n_classes=2)
    assert isinstance(code, str) and len(code) > 200


def test_mixture_rasch_generator_has_class_probs():
    code = app.generate_mixture_rasch_stan_code(n_classes=2)
    # Class probability parameter must appear somewhere for the mixture
    # to be interpretable.
    lowered = code.lower()
    assert ("class_probs" in lowered or "pi " in lowered
            or "mixing" in lowered or "dirichlet" in lowered)


# ---------------------------------------------------------------------------
# 2PL binary
# ---------------------------------------------------------------------------

def test_irt_2pl_generator_produces_non_empty_code():
    code = app.generate_2pl_binary_stan_code()
    assert isinstance(code, str) and len(code) > 200


def test_irt_2pl_generator_has_discrimination():
    code = app.generate_2pl_binary_stan_code()
    lowered = code.lower()
    assert "alpha" in lowered or "discrim" in lowered


# ---------------------------------------------------------------------------
# Pairwise BTL
# ---------------------------------------------------------------------------

def test_pairwise_btl_generator_produces_non_empty_code():
    code = app.generate_pairwise_btl_stan_code()
    assert isinstance(code, str) and len(code) > 200


def test_pairwise_btl_generator_has_ability_parameter():
    code = app.generate_pairwise_btl_stan_code()
    assert "ability" in code.lower()


# ---------------------------------------------------------------------------
# Global: every generator output parses with balanced braces
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("generator,kwargs", [
    (app.generate_dina_stan_code, {"n_items": 5, "n_attributes": 3}),
    (app.generate_hrm_stan_code, {"n_categories": 5}),
    (app.generate_testlet_stan_code, {"n_categories": 5, "bifactor": False}),
    (app.generate_testlet_stan_code, {"n_categories": 5, "bifactor": True}),
    (app.generate_mixture_rasch_stan_code, {"n_classes": 2}),
    (app.generate_2pl_binary_stan_code, {}),
    (app.generate_pairwise_btl_stan_code, {}),
])
def test_generator_output_has_balanced_braces(generator, kwargs):
    code = generator(**kwargs)
    assert _braces_balanced(code), (
        f"{generator.__name__}: unbalanced braces in generated Stan code"
    )


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("model_name", [
    "DINA", "HRM", "TESTLET_RI", "TESTLET_BIFACTOR",
    "MIXTURE_RASCH", "IRT_2PL_BINARY", "PAIRWISE_BTL",
])
def test_dispatch_returns_non_empty_code(model_name: str):
    code = app.generate_advanced_model_stan_code(model_name)
    assert isinstance(code, str) and len(code) > 100
