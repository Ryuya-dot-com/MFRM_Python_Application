"""Exercise real spawn/pipe lifecycle, including Windows pipe EOF."""
import json
import multiprocessing as mp
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest

from desktop_app import AnalysisWorker, smoke_check


def _exit_without_result(sender):
    sender.close()
    os._exit(7)


def test_spawned_worker_cancel_failure_restart_and_restore(tmp_path):
    smoke_check(tmp_path)
    record = json.loads((tmp_path / "smoke.json").read_text(encoding="utf-8"))
    assert record["worker_result"] == "complete"
    assert record["events"].count("cancelled") == 1
    assert record["events"].count("failed") == 1
    assert {(item["model"], item["method"]) for item in record["extended"]} == {
        ("RSM", "JMLE"), ("PCM", "JMLE"), ("RSM", "MML"), ("PCM", "MML"), ("G/D", "ANOVA")}
    assert len(record["extended"]) == 5
    assert record["events"].count("complete") == 1 + len(record["extended"])
    assert record["report_svg_count"] == 3
    assert (tmp_path / "wright.png").read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    restored = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parents[1] / "desktop_app.py"),
         "--check-project", str(tmp_path / "作文.mfrm"), "--output-dir", str(tmp_path / "reopened")],
        capture_output=True, text=True, timeout=60,
    )
    assert restored.returncode == 0, restored.stderr
    assert (tmp_path / "report.html").read_bytes() == (tmp_path / "reopened/report.html").read_bytes()


def test_worker_reports_exit_without_a_terminal_message():
    context = mp.get_context("spawn")
    worker = AnalysisWorker()
    worker.receiver, sender = context.Pipe(duplex=False)
    worker.process = context.Process(target=_exit_without_result, args=(sender,))
    worker.process.start()
    sender.close()
    events = []
    try:
        deadline = time.monotonic() + 30
        while worker.process is not None:
            events.extend(worker.poll())
            assert time.monotonic() < deadline, "Worker exit was never delivered"
            time.sleep(0.01)
        assert events == [("failed", "Worker exited without a result (exit=7)")]
        assert worker.receiver is None
        assert worker.poll() == []
    finally:
        worker.close()


def test_spawned_worker_completes_csv_with_missing_score(tmp_path):
    import pandas as pd
    from desktop_app import demonstration_data

    csv = tmp_path / "missing-score.csv"
    data = demonstration_data()
    data.loc[1, "Score"] = float("nan")
    data.to_csv(csv, index=False)
    worker = AnalysisWorker()
    events = []
    try:
        worker.start(pd.read_csv(csv), dict(person_col="Person", facet_cols=["Rater", "Task"],
                                          score_col="Score", rating_min=0, rating_max=3))
        deadline = time.monotonic() + 60
        while worker.process is not None:
            events.extend(worker.poll())
            assert time.monotonic() < deadline, "Missing-score analysis did not finish"
            time.sleep(0.01)
        terminal = [(status, detail) for status, detail in events if status != "progress"]
        assert len(terminal) == 1 and terminal[0][0] == "complete", terminal
        scorefile = terminal[0][1]["scorefile"]
        assert len(scorefile) == len(data) - 1
        assert 1 not in scorefile.index
        assert scorefile[["P_0", "P_1", "P_2", "P_3"]].notna().all().all()
    finally:
        worker.close()


@pytest.mark.parametrize("where", ["poll", "recv"])
@pytest.mark.parametrize("notice", ["complete", "failed"])
def test_pipe_eof_preserves_terminal_notice_until_child_exits(where, notice):
    """A sender can close before is_alive() observes the child's exit."""
    class ClosedSender:
        delivered = False

        def poll(self):
            if self.delivered and where == "poll":
                raise BrokenPipeError(109, "The pipe has been ended")
            return True

        def recv(self):
            if self.delivered:
                raise EOFError
            self.delivered = True
            return notice, "payload"

    class ExitingProcess:
        def is_alive(self):
            return True

    worker = AnalysisWorker()
    worker.receiver = ClosedSender()
    worker.process = ExitingProcess()
    assert worker.poll() == []
    assert worker.terminal == (notice, "payload")
    assert worker.poll() == []
    assert worker.terminal == (notice, "payload")
