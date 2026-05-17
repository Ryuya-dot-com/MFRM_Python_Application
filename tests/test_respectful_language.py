from __future__ import annotations

import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

TEXT_SUFFIXES = {
    ".md",
    ".py",
    ".json",
    ".toml",
    ".txt",
    ".csv",
    ".yml",
    ".yaml",
}

SKIP_DIRS = {
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    "__pycache__",
    ".ruff_cache",
    ".venv",
    "venv",
    "node_modules",
    "validation/generated",
}

RESPECTFUL_LANGUAGE_PATTERNS = [
    re.compile(r"\bbeginners?(?:[- ](?:friendly|oriented|facing))?\b", re.IGNORECASE),
    re.compile(r"\bnovices?\b", re.IGNORECASE),
    re.compile(r"\bnon[- ]?experts?\b", re.IGNORECASE),
    re.compile(r"\bnon[- ]?technical(?:\s+(?:users?|audience|readers?))?\b", re.IGNORECASE),
    re.compile(r"\blay(?:persons?|people|users?|audience|readers?)\b", re.IGNORECASE),
    re.compile(r"\binexperienced\b", re.IGNORECASE),
    re.compile(r"\bless experienced\b", re.IGNORECASE),
    re.compile(r"初心者|初学者|ビギナー|素人|非専門家"),
    re.compile(r"BeginnerAction|BeginnerInterpretation|beginner_case_guidance"),
]

ALLOWED_QUOTED_REFERENCES = (
    "layperson ratings of divergent",
)


def _is_scanned_path(path: Path) -> bool:
    rel = path.relative_to(REPO_ROOT).as_posix()
    if path.name == "test_respectful_language.py":
        return False
    if path.suffix.lower() not in TEXT_SUFFIXES:
        return False
    return not any(rel == skipped or rel.startswith(f"{skipped}/") for skipped in SKIP_DIRS)


def _is_allowed_reference(path: Path, line: str) -> bool:
    rel = path.relative_to(REPO_ROOT).as_posix()
    lowered = line.lower()
    return rel in {"locales/en.json", "locales/ja.json"} and any(
        reference in lowered for reference in ALLOWED_QUOTED_REFERENCES
    )


def test_respectful_language_avoids_user_ability_labels():
    violations: list[str] = []
    for path in sorted(REPO_ROOT.rglob("*")):
        if not path.is_file() or not _is_scanned_path(path):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for line_number, line in enumerate(text.splitlines(), start=1):
            if _is_allowed_reference(path, line):
                continue
            for pattern in RESPECTFUL_LANGUAGE_PATTERNS:
                match = pattern.search(line)
                if match:
                    rel = path.relative_to(REPO_ROOT).as_posix()
                    violations.append(f"{rel}:{line_number}: {match.group(0)!r}")

    assert not violations, (
        "Use task-centered wording such as 'guided', 'first-read', "
        "'interpretation support', or 'review checklist' instead of labeling "
        "users by presumed ability or experience:\n" + "\n".join(violations)
    )
