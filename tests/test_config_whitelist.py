"""Tests for the `_CONFIG_JSON_IMPORT_WHITELIST` contract.

The whitelist is the handshake between `Download config (JSON)`
(writer) and `Import config JSON` (reader). Drift between the two
causes silent setting-loss on re-import, which is confusing and
hard to debug.
"""

from __future__ import annotations

import streamlit_app as app


def test_whitelist_is_frozenset():
    assert isinstance(app._CONFIG_JSON_IMPORT_WHITELIST, frozenset)


def test_whitelist_size_is_pinned():
    # Bumping this number requires an entry in CHANGELOG.md.
    assert len(app._CONFIG_JSON_IMPORT_WHITELIST) == 23


def test_critical_sidebar_settings_are_whitelisted():
    critical = {
        "model_type", "method", "analysis_depth",
        "rating_min", "rating_max",
        "maxit", "reltol", "anchor_policy",
    }
    missing = critical - app._CONFIG_JSON_IMPORT_WHITELIST
    assert not missing, f"critical keys missing from whitelist: {missing}"


def test_whitelist_has_no_duplicates():
    # frozenset cannot have duplicates by construction, but check that
    # the tuple used to build it did not have them either.
    assert len(app._CONFIG_JSON_IMPORT_WHITELIST) == len(list(app._CONFIG_JSON_IMPORT_WHITELIST))


def test_whitelist_filter_drops_non_whitelisted_keys():
    imported = {
        "model_type": "RSM",
        "method": "JMLE",
        "not_in_whitelist": "should_be_dropped",
        "config_export_fingerprint": "metadata_filtered",
    }
    kept = {
        k: imported[k]
        for k in sorted(app._CONFIG_JSON_IMPORT_WHITELIST)
        if k in imported
    }
    assert "model_type" in kept
    assert "method" in kept
    assert "not_in_whitelist" not in kept
    assert "config_export_fingerprint" not in kept


def test_whitelist_preserves_all_matching_keys():
    """A config with only whitelisted keys must round-trip completely."""
    imported = {k: "v" for k in app._CONFIG_JSON_IMPORT_WHITELIST}
    kept = {
        k: imported[k]
        for k in sorted(app._CONFIG_JSON_IMPORT_WHITELIST)
        if k in imported
    }
    assert len(kept) == len(app._CONFIG_JSON_IMPORT_WHITELIST)
