"""Windows-only FACETS BATCH=NO launcher with output-gated graceful close.

FACETS 4.5.0/4.5.1 can crash in their internal hidden ``BATCH=YES`` path on
some Windows builds while visible ``BATCH=NO`` analysis succeeds.  This module
does not treat file creation alone as completion: every registered output must
be nonempty and size/mtime-stable before WM_CLOSE is posted to the FACETS
window.  Forced termination is deliberately not a success path.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import subprocess
import time
from typing import Iterable, Sequence


@dataclass(frozen=True)
class VisibleFacetsInvocation:
    """Result and operational evidence from a visible-mode FACETS call."""

    completed: subprocess.CompletedProcess[str]
    readiness_paths: tuple[Path, ...]
    stable_seconds: float
    windows_closed: int
    forced_termination: bool


def _output_state(paths: Sequence[Path]) -> tuple[tuple[str, int, int], ...] | None:
    rows: list[tuple[str, int, int]] = []
    for path in paths:
        try:
            stat = path.stat()
        except FileNotFoundError:
            return None
        if not path.is_file() or stat.st_size <= 0:
            return None
        rows.append((str(path.resolve()), int(stat.st_size), int(stat.st_mtime_ns)))
    return tuple(rows)


def _post_close_to_process_windows(process_id: int) -> int:
    """Post WM_CLOSE to every top-level window owned by ``process_id``."""

    if os.name != "nt":
        return 0
    import ctypes
    from ctypes import wintypes

    user32 = ctypes.WinDLL("user32", use_last_error=True)
    callback_type = ctypes.WINFUNCTYPE(wintypes.BOOL, wintypes.HWND, wintypes.LPARAM)
    window_process_id = wintypes.DWORD()
    closed = 0

    @callback_type
    def callback(window: int, _parameter: int) -> bool:
        nonlocal closed
        user32.GetWindowThreadProcessId(window, ctypes.byref(window_process_id))
        if int(window_process_id.value) == int(process_id):
            if user32.PostMessageW(window, 0x0010, 0, 0):  # WM_CLOSE
                closed += 1
        return True

    if not user32.EnumWindows(callback, 0):
        error = ctypes.get_last_error()
        if error:
            raise OSError(error, "EnumWindows failed")
    return closed


def invoke_facets_batch_no(
    command: Sequence[str],
    *,
    cwd: Path,
    readiness_paths: Iterable[Path],
    timeout_seconds: float,
    stable_seconds: float = 2.0,
    poll_seconds: float = 0.10,
    close_timeout_seconds: float = 10.0,
    hide_window: bool = True,
) -> VisibleFacetsInvocation:
    """Run FACETS in ``BATCH=NO`` mode and close it after durable outputs exist.

    A successful return requires all readiness files to be nonempty and stable,
    at least one FACETS-owned window to accept WM_CLOSE, and process exit code
    zero.  Timeouts and missing windows terminate the child and raise; a hard
    kill is never normalized into success.
    """

    if os.name != "nt":
        raise OSError("FACETS visible-mode automation is supported only on Windows")
    paths = tuple(Path(path).resolve() for path in readiness_paths)
    if not paths:
        raise ValueError("At least one readiness path is required")
    if timeout_seconds <= 0 or stable_seconds <= 0 or poll_seconds <= 0:
        raise ValueError("Timeout and polling values must be positive")

    startupinfo = subprocess.STARTUPINFO()
    if hide_window:
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        startupinfo.wShowWindow = getattr(subprocess, "SW_HIDE", 0)
    process = subprocess.Popen(
        [str(part) for part in command],
        cwd=Path(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        creationflags=(
            getattr(subprocess, "CREATE_NO_WINDOW", 0) if hide_window else 0
        ),
        startupinfo=startupinfo,
    )
    started = time.monotonic()
    last_state: tuple[tuple[str, int, int], ...] | None = None
    stable_since: float | None = None
    windows_closed = 0
    forced_termination = False
    try:
        while True:
            return_code = process.poll()
            if return_code is not None:
                stdout, stderr = process.communicate()
                return VisibleFacetsInvocation(
                    completed=subprocess.CompletedProcess(
                        [str(part) for part in command], return_code, stdout, stderr
                    ),
                    readiness_paths=paths,
                    stable_seconds=0.0,
                    windows_closed=0,
                    forced_termination=False,
                )
            elapsed = time.monotonic() - started
            if elapsed >= timeout_seconds:
                raise subprocess.TimeoutExpired(command, timeout_seconds)
            state = _output_state(paths)
            now = time.monotonic()
            if state is None:
                last_state = None
                stable_since = None
            elif state != last_state:
                last_state = state
                stable_since = now
            elif stable_since is not None and now - stable_since >= stable_seconds:
                windows_closed = _post_close_to_process_windows(process.pid)
                if windows_closed <= 0:
                    raise RuntimeError(
                        "FACETS outputs became stable but no process window accepted WM_CLOSE"
                    )
                try:
                    return_code = process.wait(timeout=close_timeout_seconds)
                except subprocess.TimeoutExpired as exc:
                    raise RuntimeError(
                        "FACETS did not exit after its window accepted WM_CLOSE"
                    ) from exc
                stdout, stderr = process.communicate()
                final_state = _output_state(paths)
                if final_state is None:
                    raise RuntimeError("FACETS readiness outputs disappeared during close")
                return VisibleFacetsInvocation(
                    completed=subprocess.CompletedProcess(
                        [str(part) for part in command], return_code, stdout, stderr
                    ),
                    readiness_paths=paths,
                    stable_seconds=now - stable_since,
                    windows_closed=windows_closed,
                    forced_termination=False,
                )
            time.sleep(poll_seconds)
    except BaseException:
        if process.poll() is None:
            forced_termination = True
            process.terminate()
            try:
                process.wait(timeout=5.0)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=5.0)
        process.communicate()
        raise
