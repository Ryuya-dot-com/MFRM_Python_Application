import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import threading
import time

import pandas as pd

from validation.estimand_distribution_study import generate_study_bundle
from validation.facets_resilient_pair import (
    cross_process_file_lock,
    dependency_manifest,
    dependency_manifest_digest,
    execute_with_bounded_retry,
    fit_resilient_pair,
    retryable_missing_report_failure,
)


def test_dependency_manifest_binds_application_and_direct_validation_code():
    manifest = dependency_manifest()
    assert "streamlit_app.py" in manifest
    assert "validation/estimand_distribution_study.py" in manifest
    assert "validation/facets_resilient_pair.py" in manifest
    assert "requirements.txt" in manifest
    assert any(path.startswith("mfrm_app/") for path in manifest)
    assert all(len(value) == 64 for value in manifest.values())
    assert len(dependency_manifest_digest(manifest)) == 64


def test_retry_policy_retries_only_exit_zero_missing_report(tmp_path):
    calls = []

    def operation(try_dir, try_number):
        calls.append(try_number)
        if try_number == 1:
            raise RuntimeError("FACETS primary failed for injected::PCM: exit=0")
        return "qualified"

    result = execute_with_bounded_retry(
        operation,
        attempts_root=tmp_path / "retryable",
        maximum_retries=2,
    )
    assert result.result == "qualified"
    assert calls == [1, 2]
    assert result.records[0]["RetryEligible"] is True
    assert result.records[1]["Succeeded"] is True
    first = json.loads(
        (tmp_path / "retryable" / "pair_try_01" / "try_outcome.json").read_text()
    )
    assert first["RetryEligible"] is True

    calls.clear()

    def nonretryable(_try_dir, try_number):
        calls.append(try_number)
        raise RuntimeError("FACETS primary failed for injected::PCM: exit=5")

    result = execute_with_bounded_retry(
        nonretryable,
        attempts_root=tmp_path / "nonretryable",
        maximum_retries=2,
    )
    assert result.result is None
    assert calls == [1]
    assert result.records[0]["RetryEligible"] is False


def test_retry_policy_exhausts_exactly_three_retained_tries(tmp_path):
    def operation(_try_dir, _try_number):
        raise RuntimeError("FACETS auxiliary failed for injected::PCM: exit=0")

    result = execute_with_bounded_retry(
        operation,
        attempts_root=tmp_path / "exhaustion",
        maximum_retries=2,
    )
    assert result.result is None
    assert len(result.records) == 3
    assert all(record["RetryEligible"] for record in result.records)
    assert len(list((tmp_path / "exhaustion").glob("pair_try_*/try_outcome.json"))) == 3


def test_retryable_contract_requires_missing_named_report(tmp_path):
    try_dir = tmp_path / "try"
    try_dir.mkdir()
    error = RuntimeError("FACETS primary failed for injected::PCM: exit=0")
    assert retryable_missing_report_failure(error, try_dir)
    nested = try_dir / "facets_runs" / "run"
    nested.mkdir(parents=True)
    (nested / "report_u6.txt").write_text("present", encoding="utf-8")
    assert not retryable_missing_report_failure(error, try_dir)


def test_cross_process_lock_also_serializes_threads(tmp_path):
    active = 0
    maximum_active = 0
    guard = threading.Lock()

    def worker():
        nonlocal active, maximum_active
        with cross_process_file_lock(tmp_path / "facets.lock", timeout_seconds=5):
            with guard:
                active += 1
                maximum_active = max(maximum_active, active)
            time.sleep(0.05)
            with guard:
                active -= 1

    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(lambda _index: worker(), range(3)))
    assert maximum_active == 1


def test_python_evidence_survives_exhausted_facets_retries(tmp_path):
    bundle = generate_study_bundle(replicates=1)
    manifest_row = bundle["manifest.csv"].iloc[0]
    run_id = str(manifest_row["RunId"])
    ratings = bundle["generated_ratings.csv"]
    ratings = ratings[ratings["RunId"].astype(str).eq(run_id)].drop(columns="RunId")
    truth = bundle["generated_facet_truth.csv"]
    truth = truth[truth["RunId"].astype(str).eq(run_id)]
    thresholds = bundle["generated_pcm_threshold_truth.csv"]
    thresholds = thresholds[thresholds["RunId"].astype(str).eq(run_id)]
    fake_facets = tmp_path / "Facets.exe"
    fake_facets.write_bytes(b"diagnostic placeholder")

    def always_missing_report(*_args, **_kwargs):
        raise RuntimeError(f"FACETS primary failed for {run_id}::PCM: exit=0")

    metrics = fit_resilient_pair(
        manifest_row,
        ratings,
        truth,
        thresholds,
        facets_exe=fake_facets,
        output_dir=tmp_path / "resilient",
        lock_path=tmp_path / "facets.lock",
        pair_fit=always_missing_report,
    )
    assert metrics["statistical_evidence_ready"] is True
    assert metrics["calibration_ready"] is False
    assert metrics["pair_fully_qualified"] is False
    assert metrics["facets_tries"] == 3
    assert (tmp_path / "resilient" / "python_independent" / "python_completion.json").is_file()
    runs = pd.read_csv(tmp_path / "resilient" / "combined_run_ledger.csv")
    python = runs[runs["EstimatorMode"].eq("PYTHON_JMLE")].iloc[0]
    facets = runs[runs["EstimatorMode"].eq("FACETS_4_5_JMLE")].iloc[0]
    assert bool(python["IncludedInStudy"])
    assert not bool(facets["IncludedInStudy"])
