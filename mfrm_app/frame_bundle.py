from __future__ import annotations

from collections.abc import Iterable, MutableMapping
import re

import pandas as pd


def safe_frame_key(prefix: str, label: object, *, max_label_length: int | None = None) -> str:
    """Build a filesystem-friendly frame key from a prefix and arbitrary label."""
    safe_label = re.sub(r"[^A-Za-z0-9]+", "_", str(label)).strip("_")
    if max_label_length is not None:
        safe_label = safe_label[:max_label_length]
    return f"{prefix}_{safe_label}" if safe_label else prefix


def is_nonempty_frame(value: object) -> bool:
    """Return True for non-empty pandas DataFrames only."""
    return isinstance(value, pd.DataFrame) and not value.empty


def add_frame(
    frames: MutableMapping[str, pd.DataFrame],
    name: str,
    value: object,
    *,
    overwrite: bool = True,
    allow_empty: bool = False,
) -> bool:
    """Add a DataFrame to a frame bundle and report whether it was added."""
    if not isinstance(value, pd.DataFrame):
        return False
    if value.empty and not allow_empty:
        return False
    if not overwrite and name in frames:
        return False
    frames[name] = value
    return True


def add_frames(
    frames: MutableMapping[str, pd.DataFrame],
    items: Iterable[tuple[str, object]],
    *,
    overwrite: bool = True,
    allow_empty: bool = False,
) -> int:
    """Add all DataFrames from ``items`` and return the count added."""
    added = 0
    for name, value in items:
        if add_frame(frames, name, value, overwrite=overwrite, allow_empty=allow_empty):
            added += 1
    return added


def add_iterable_frame(
    frames: MutableMapping[str, pd.DataFrame],
    name: str,
    value: object,
    *,
    column_name: str,
    overwrite: bool = True,
    allow_empty: bool = False,
) -> bool:
    """Add a DataFrame from an iterable value, while leaving strings/mappings alone."""
    if isinstance(value, pd.DataFrame):
        return add_frame(frames, name, value, overwrite=overwrite, allow_empty=allow_empty)
    if isinstance(value, (str, bytes, dict)) or not isinstance(value, Iterable):
        return False
    try:
        frame = pd.DataFrame({column_name: list(value)})
    except Exception:
        return False
    return add_frame(frames, name, frame, overwrite=overwrite, allow_empty=allow_empty)
