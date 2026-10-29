"""Test process isolation and selection of locally retained study artifacts."""

import sys
import pytest


@pytest.fixture(autouse=True)
def restore_main_module_after_test():
    """AppTest installs a temporary __main__; spawn must see pytest's entrypoint."""
    original_main = sys.modules["__main__"]
    try:
        yield
    finally:
        sys.modules["__main__"] = original_main


_RETAINED_STUDY_TESTS = {
    # These seals use the retained timestamp trust root/receipt and local verifier.
    "tests/test_bias_population_controller.py::test_absent_or_probe_receipt_never_reaches_generation",
    "tests/test_bias_population_controller.py::test_seed_inventory_cannot_be_changed_by_updating_unsigned_local_manifest",
    "tests/test_bias_population_controller.py::test_seal_cannot_omit_required_source_or_runtime_bindings",
    "tests/test_bias_population_controller.py::test_restore_requires_exact_inventory_and_no_extra_files",
    "tests/test_mml_population_controller.py::test_absent_or_wrong_receipt_cannot_authorize_a_frozen_plan",
    "tests/test_mml_population_controller.py::test_source_bundle_change_fails_before_receipt_verification",
    "tests/test_mml_population_controller.py::test_ordered_seed_inventory_is_checked_semantically_even_with_updated_local_hash",
    "tests/test_known_assignment_confirmatory_strict_audit.py::test_facets_pair_metrics_reject_fractional_tries_and_bad_replay",
    "tests/test_known_assignment_confirmatory_strict_audit.py::test_registered_kac200_bundle_passes_strict_audit_without_publication",
    "tests/test_known_assignment_confirmatory_strict_audit_publish_v2.py::test_file_atomic_transport_publishes_exact_v1_bundle",
    "tests/test_known_assignment_confirmatory_strict_audit_publish_v2.py::test_file_atomic_transport_failure_removes_completion_and_owned_directory",
}


def pytest_collection_modifyitems(items):
    # Preserve the registered test-file hashes; selection is not study evidence.
    for item in items:
        if item.nodeid in _RETAINED_STUDY_TESTS:
            item.add_marker(pytest.mark.retained_evidence)
