"""Bilingual and exposure gates for the shadow Help content contracts."""

from __future__ import annotations

import json
from pathlib import Path
import re
import string

import pytest

from mfrm_app.help_topics import (
    HELP_ALL_LOCALE_KEYS,
    HELP_NAV_LOCALE_KEYS,
    HELP_REQUIRED_LOCALE_KEYS,
)
from mfrm_app.terminology import TERMINOLOGY_REGISTRY
from mfrm_app.user_problems import (
    USER_PROBLEM_SPECS,
    USER_PROBLEM_SUPPORT_REFERENCE_LABEL_KEY,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
LOCALE_PATHS = {
    language: REPO_ROOT / "locales" / f"{language}.json"
    for language in ("en", "ja")
}


def _flatten(node: dict[str, object], prefix: str = "") -> dict[str, str]:
    leaves: dict[str, str] = {}
    for key, value in node.items():
        dotted = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            leaves.update(_flatten(value, dotted))
        elif isinstance(value, str):
            leaves[dotted] = value
    return leaves


def _problem_locale_keys() -> frozenset[str]:
    return frozenset(
        key
        for spec in USER_PROBLEM_SPECS.values()
        for key in (spec.title_key, spec.body_key, *spec.action_keys)
    ) | {USER_PROBLEM_SUPPORT_REFERENCE_LABEL_KEY}


def _placeholders(value: str) -> set[str]:
    return {
        field_name
        for _, field_name, _, _ in string.Formatter().parse(value)
        if field_name
    }


@pytest.fixture(scope="module")
def locale_leaves() -> dict[str, dict[str, str]]:
    return {
        language: _flatten(json.loads(path.read_text(encoding="utf-8")))
        for language, path in LOCALE_PATHS.items()
    }


def test_every_registered_contract_key_has_nonempty_bilingual_copy(locale_leaves):
    required = (
        HELP_ALL_LOCALE_KEYS
        | TERMINOLOGY_REGISTRY.locale_keys
        | _problem_locale_keys()
    )

    for language, leaves in locale_leaves.items():
        missing = sorted(required.difference(leaves))
        assert not missing, f"Missing registered {language} locale keys: {missing}"
        empty = sorted(key for key in required if not leaves[key].strip())
        assert not empty, f"Empty registered {language} locale values: {empty}"

        topic_keys = {key for key in leaves if key.startswith("help_topics.")}
        nav_keys = {key for key in leaves if key.startswith("help_nav.")}
        assert topic_keys == {
            key for key in HELP_REQUIRED_LOCALE_KEYS if key.startswith("help_topics.")
        }
        assert nav_keys == set(HELP_NAV_LOCALE_KEYS)

        section_bodies = [
            leaves[key]
            for key in topic_keys
            if ".sections." in key and key.endswith(".body")
        ]
        assert len(section_bodies) == len(set(section_bodies))


def test_registered_copy_has_placeholder_parity(locale_leaves):
    en = locale_leaves["en"]
    ja = locale_leaves["ja"]
    required = (
        HELP_ALL_LOCALE_KEYS
        | TERMINOLOGY_REGISTRY.locale_keys
        | _problem_locale_keys()
    )

    mismatches = {
        key: (_placeholders(en[key]), _placeholders(ja[key]))
        for key in required
        if _placeholders(en[key]) != _placeholders(ja[key])
    }
    assert not mismatches


def test_help_copy_does_not_expose_implementation_or_route_vocabulary(locale_leaves):
    internal_patterns = {
        "raw contract label": re.compile(
            r"\b(?:AnalysisID|EvidenceRecord|ReasonCode|ComputationState|StabilityState)\b",
            re.IGNORECASE,
        ),
        "implementation noun": re.compile(
            r"\b(?:backend|payload|schema|fixture|pipeline|manifest|runner)\b",
            re.IGNORECASE,
        ),
        "implementation state": re.compile(
            r"\b(?:cache|session state|stack trace|environment variable)\b",
            re.IGNORECASE,
        ),
        "function or environment name": re.compile(
            r"\b(?:compute_pca(?:_[a-z0-9_]+)?|MFRM_[A-Z0-9_]+)\b"
        ),
        "route identity": re.compile(r"\b(?:help|target|link)\.[a-z0-9_.:-]+"),
        "raw file extension": re.compile(r"\.(?:json|csv|zip)\b", re.IGNORECASE),
    }

    violations: list[str] = []
    for language, leaves in locale_leaves.items():
        for key in sorted(HELP_REQUIRED_LOCALE_KEYS | _problem_locale_keys()):
            value = leaves[key]
            for label, pattern in internal_patterns.items():
                if pattern.search(value):
                    violations.append(f"{language}:{key}: {label}")
    assert not violations, "Internal vocabulary reached Help copy:\n" + "\n".join(
        violations
    )


def test_help_copy_has_no_external_runtime_call_to_action(locale_leaves):
    external_patterns = {
        "TAM": re.compile(r"\bTAM\b"),
        "FACETS product": re.compile(r"\bFACETS\b"),
        "ConQuest": re.compile(r"\bConQuest\b", re.IGNORECASE),
        "mirt": re.compile(r"\bmirt\b", re.IGNORECASE),
        "R runtime": re.compile(r"(?:^|[\s(])R(?:script)?(?:$|[\s),.])"),
        "Julia": re.compile(r"\bJulia\b", re.IGNORECASE),
        "Stan": re.compile(r"\b(?:Stan|CmdStan)\b", re.IGNORECASE),
        "Posterior Viewer": re.compile(r"\bPosterior Viewer\b", re.IGNORECASE),
        "external engine": re.compile(
            r"\b(?:external engine|cross[- ]engine)\b", re.IGNORECASE
        ),
    }

    violations: list[str] = []
    for language, leaves in locale_leaves.items():
        for key in sorted(HELP_REQUIRED_LOCALE_KEYS | _problem_locale_keys()):
            value = leaves[key]
            for label, pattern in external_patterns.items():
                if pattern.search(value):
                    violations.append(f"{language}:{key}: {label}")
    assert not violations, "External execution route reached Help copy:\n" + "\n".join(
        violations
    )


def test_claim_critical_topics_keep_their_interpretation_boundaries(locale_leaves):
    en = locale_leaves["en"]
    ja = locale_leaves["ja"]

    residual_en = en["help_topics.results.residual_structure.cannot_show"].casefold()
    residual_ja = ja["help_topics.results.residual_structure.cannot_show"]
    assert "unidimensional" in residual_en
    assert any(boundary in residual_en for boundary in ("cannot", "does not", "not "))
    assert "一次元" in residual_ja
    assert any(boundary in residual_ja for boundary in ("できません", "証明しません", "示すものではありません"))

    interaction_en = en[
        "help_topics.results.differential_interaction.cannot_show"
    ].casefold()
    interaction_ja = ja[
        "help_topics.results.differential_interaction.cannot_show"
    ]
    assert "fair" in interaction_en
    assert any(boundary in interaction_en for boundary in ("cannot", "does not", "not "))
    assert "公平" in interaction_ja
    assert any(boundary in interaction_ja for boundary in ("できません", "ではありません", "示しません"))

    fit_en = en["help_topics.results.fit.cannot_show"].casefold()
    fit_ja = ja["help_topics.results.fit.cannot_show"]
    assert any(
        word in fit_en
        for word in ("remove", "removal", "exclude", "exclusion", "quality")
    )
    assert any(boundary in fit_en for boundary in ("cannot", "does not", "not "))
    assert any(word in fit_ja for word in ("除外", "質", "良否"))
    assert any(boundary in fit_ja for boundary in ("できません", "ではありません", "決めません"))

    privacy_en = en["help_topics.downloads.privacy.cannot_show"].casefold()
    privacy_ja = ja["help_topics.downloads.privacy.cannot_show"]
    assert any(
        word in privacy_en for word in ("privacy", "anonymous", "anonymity", "confidential")
    )
    assert any(boundary in privacy_en for boundary in ("cannot", "does not", "not "))
    assert any(
        word in privacy_ja for word in ("プライバシー", "個人", "非識別", "匿名")
    )
    assert any(boundary in privacy_ja for boundary in ("できません", "保証しません", "保証ではありません"))
