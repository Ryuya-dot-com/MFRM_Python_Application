"""Hygiene gate that blocks internal phase / sprint section labels.

The repository must not ship internal phase / sprint codenames in any
artefact a reader of the source can see: code comments, docstrings,
UI strings, i18n locale files, documentation, and tests. Internal
labels are useful while a feature is being planned but they look
process-leaky once they land in the public codebase, and they make
diffs harder to follow when section organisers shift over time.

This gate scans the working tree for the forbidden patterns and
fails the test suite on any hit. Add new descriptive prose instead.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest


# NATO phonetic alphabet + number words + a small set of common milestone
# words. Using a finite vocabulary rather than ``\b[A-Z][a-z]+\b`` keeps prose
# like "phase transition" or "phase diagram" out of scope; capitalized common
# nouns that happen to follow "Phase " do not trigger the gate unless they
# match the vocabulary below.
_MULTI_LETTER_CODENAMES = (
    # NATO phonetic alphabet
    "Alpha", "Bravo", "Charlie", "Delta", "Echo", "Foxtrot", "Golf", "Hotel",
    "India", "Juliet", "Kilo", "Lima", "Mike", "November", "Oscar", "Papa",
    "Quebec", "Romeo", "Sierra", "Tango", "Uniform", "Victor", "Whiskey",
    "Yankee", "Zulu",
    # Spelled-out numbers
    "Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight",
    "Nine", "Ten",
    # Greek letters frequently used as milestone codenames
    "Beta", "Gamma", "Epsilon", "Theta", "Lambda", "Sigma", "Omega",
    # Common release / sprint vocabulary
    "Initial", "Final", "Pilot", "Cleanup", "Polish", "Rollout",
)


FORBIDDEN = re.compile(
    # Numeric milestone codenames: the label keyword followed by a number,
    # optionally with a sub-component (1-5, 2.0) and/or a trailing lowercase
    # suffix letter (7g).
    r"\b[Pp]hase[ -][0-9]+(?:[-.][0-9]+)?[a-z]?\b"
    # Alphabetic single-letter codenames: the label keyword followed by an
    # isolated capital letter at a word boundary. The trailing word boundary
    # keeps prose like "the alpha phase" or all-caps acronyms ("Phase UX")
    # out of scope -- only an isolated capital letter qualifies.
    r"|\b[Pp]hase[ -][A-Z]\b"
    # Multi-letter codenames drawn from a finite vocabulary so prose like
    # "phase transition" cannot self-trip the gate.
    r"|\b[Pp]hase[ -](?:" + "|".join(_MULTI_LETTER_CODENAMES) + r")\b"
    # Sprint codenames, in case sprint-style numbering ever reappears.
    r"|\b[Ss]print[ -][0-9]+\b"
)


# The hygiene test cannot scan its own source: literal examples of the
# forbidden patterns appear in this file by design (the regex above must
# pin them) and would otherwise self-trip the gate. The test for the gate
# itself is exercised by every other file in the working tree.
SELF_FILE = Path(__file__).resolve()

REPO_ROOT = Path(__file__).resolve().parent.parent

SCAN_PATTERNS: tuple[str, ...] = (
    "streamlit_app.py",
    "tests/*.py",
    "locales/*.json",
    "validation/*.md",
    "validation/**/*.md",
    "docs/*.md",
    "docs/**/*.md",
    "*.md",
    # CI workflow files: Phase labels in workflow names or comments would be
    # equally visible to anyone browsing the repository.
    ".github/workflows/*.yml",
    ".github/workflows/*.yaml",
    ".github/*.md",
    ".github/**/*.md",
    # Project metadata files.
    "*.toml",
    "*.cfg",
    "Makefile",
)


def _candidate_files() -> list[Path]:
    seen: set[Path] = set()
    for pattern in SCAN_PATTERNS:
        for path in REPO_ROOT.glob(pattern):
            if path.is_file() and path.resolve() != SELF_FILE:
                seen.add(path)
    return sorted(seen)


@pytest.mark.parametrize("path", _candidate_files(), ids=lambda p: str(p.relative_to(REPO_ROOT)))
def test_file_contains_no_internal_phase_or_sprint_labels(path: Path) -> None:
    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            match = FORBIDDEN.search(line)
            if match is not None:
                rel = path.relative_to(REPO_ROOT)
                raise AssertionError(
                    f"Internal label '{match.group()}' found in "
                    f"{rel}:{line_no} -- rewrite using a descriptive phrase "
                    f"(e.g. the actual section content, not the milestone "
                    f"codename). Offending line:\n  {line.rstrip()}"
                )
