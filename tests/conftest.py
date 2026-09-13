"""Release selection for test files whose bytes are frozen in study records."""

import pytest


_RETAINED_STUDY_TESTS = {
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
