"""Distribution-readiness helpers for reproducible public releases."""

from __future__ import annotations

import re
from collections.abc import Mapping
from pathlib import Path

import pandas as pd


_REQUIREMENT_FLOOR_RE = re.compile(
    r"^\s*([A-Za-z0-9_.-]+)(?:\[[^\]]+\])?\s*>=\s*([^,\s;#]+)"
)


def package_key(package_name: str) -> str:
    """Return a PEP 503-style key for comparing package names."""
    return re.sub(r"[-_.]+", "-", str(package_name).strip()).lower()


def version_tuple(version_text: str) -> tuple[int, ...]:
    parts: list[int] = []
    for token in str(version_text).replace("-", ".").split("."):
        if token.isdigit():
            parts.append(int(token))
            continue
        digits = "".join(ch for ch in token if ch.isdigit())
        if digits:
            parts.append(int(digits))
        break
    return tuple(parts or [0])


def compare_versions(left: str, right: str) -> int:
    """Compare dotted numeric version floors."""
    left_tuple = version_tuple(left)
    right_tuple = version_tuple(right)
    max_len = max(len(left_tuple), len(right_tuple))
    left_tuple = left_tuple + (0,) * (max_len - len(left_tuple))
    right_tuple = right_tuple + (0,) * (max_len - len(right_tuple))
    if left_tuple == right_tuple:
        return 0
    return 1 if left_tuple > right_tuple else -1


def parse_requirement_floors(requirements_text: str) -> dict[str, dict[str, str]]:
    """Parse package lower bounds from requirements.txt-style content."""
    floors: dict[str, dict[str, str]] = {}
    for line_number, raw_line in enumerate(str(requirements_text).splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#") or line.startswith("-"):
            continue
        match = _REQUIREMENT_FLOOR_RE.match(line)
        if not match:
            continue
        package_name, floor = match.groups()
        floors[package_key(package_name)] = {
            "RequirementPackage": package_name,
            "RequirementFloor": floor,
            "RequirementLine": str(line_number),
        }
    return floors


def requirements_floor_contract(
    requirements_text: str,
    runtime_package_floors: Mapping[str, str],
) -> pd.DataFrame:
    """Compare requirements.txt floors with the runtime doctor contract."""
    requirement_floors = parse_requirement_floors(requirements_text)
    doctor_floors = {
        package_key(package_name): {
            "DoctorPackage": str(package_name),
            "DoctorFloor": str(floor),
        }
        for package_name, floor in runtime_package_floors.items()
    }
    keys = sorted(set(requirement_floors) | set(doctor_floors))
    rows: list[dict[str, str]] = []
    for key in keys:
        requirement = requirement_floors.get(key, {})
        doctor = doctor_floors.get(key, {})
        requirement_floor = requirement.get("RequirementFloor", "")
        doctor_floor = doctor.get("DoctorFloor", "")
        display_name = doctor.get("DoctorPackage") or requirement.get("RequirementPackage") or key

        if not requirement:
            status = "Blocker"
            detail = "Doctor checks this runtime package, but requirements.txt does not publish a floor."
        elif not doctor:
            status = "Review"
            detail = "requirements.txt publishes this floor, but the doctor does not check it."
        else:
            comparison = compare_versions(requirement_floor, doctor_floor)
            if comparison == 0:
                status = "Ready"
                detail = "requirements.txt and doctor use the same runtime floor."
            elif comparison > 0:
                status = "Review"
                detail = "requirements.txt is stricter than the doctor; update the doctor floor."
            else:
                status = "Blocker"
                detail = "requirements.txt is looser than the doctor; update the published floor."

        rows.append({
            "Package": display_name,
            "RequirementPackage": requirement.get("RequirementPackage", ""),
            "RequirementFloor": requirement_floor,
            "RequirementLine": requirement.get("RequirementLine", ""),
            "DoctorPackage": doctor.get("DoctorPackage", ""),
            "DoctorFloor": doctor_floor,
            "Status": status,
            "Detail": detail,
        })
    return pd.DataFrame(rows)


def requirements_floor_contract_from_file(
    requirements_path: str | Path,
    runtime_package_floors: Mapping[str, str],
) -> pd.DataFrame:
    """Load requirements.txt and compare it to the runtime doctor contract."""
    path = Path(requirements_path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        return pd.DataFrame([
            {
                "Package": package_name,
                "RequirementPackage": "",
                "RequirementFloor": "",
                "RequirementLine": "",
                "DoctorPackage": package_name,
                "DoctorFloor": str(floor),
                "Status": "Blocker",
                "Detail": f"Could not read {path.name}: {exc}",
            }
            for package_name, floor in runtime_package_floors.items()
        ])
    return requirements_floor_contract(text, runtime_package_floors)

