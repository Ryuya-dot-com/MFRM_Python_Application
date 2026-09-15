"""Software citation from CITATION.cff (JSON syntax is valid YAML 1.2)."""

import json
from pathlib import Path


CFF_PATH = Path(__file__).resolve().parents[1] / "CITATION.cff"
CFF_TEXT = CFF_PATH.read_text(encoding="utf-8")
METADATA = json.loads(CFF_TEXT)
AUTHOR = METADATA["authors"][0]["family-names"]
GIVEN_NAMES = METADATA["authors"][0]["given-names"]
INITIALS = " ".join(f"{part[0]}." for part in GIVEN_NAMES.split())
YEAR = METADATA["date-released"][:4]
IN_TEXT = f"({AUTHOR}, {YEAR})"
APA = (
    f"{AUTHOR}, {INITIALS} ({YEAR}). {METADATA['title']} "
    f"(Version {METADATA['version']}) [Computer software]. {METADATA['url']}"
)
BIBTEX = (
    "@software{mfrm_streamlit,\n"
    f"  author = {{{AUTHOR}, {GIVEN_NAMES}}},\n"
    f"  title = {{{{{METADATA['title']}}}}},\n"
    f"  year = {{{YEAR}}},\n"
    f"  version = {{{METADATA['version']}}},\n"
    f"  url = {{{METADATA['url']}}}\n"
    "}\n"
)
APA_MARKDOWN = APA.replace(METADATA["title"], f"*{METADATA['title']}*", 1)
UNAVAILABLE = (
    "Software citation metadata is not verified for this result's recorded version. "
    "Check the saved analysis configuration and the corresponding software release "
    "before completing the reference. Do not substitute the currently running app version."
)


def matches_result(result: object) -> bool:
    """Only use release metadata when the fit records that exact release."""
    config = result.get("config") if isinstance(result, dict) else None
    return isinstance(config, dict) and config.get("app_version") == METADATA["version"]


def analysis_statement(result: object) -> str:
    if not matches_result(result):
        return UNAVAILABLE
    return (
        f"The recorded analysis used {METADATA['title']} {IN_TEXT}, "
        f"Version {METADATA['version']}."
    )


def result_assets(result: object) -> dict[str, str]:
    """Citation sidecars contain no response rows or user-supplied text."""
    text = "# Software citation\n\n" + analysis_statement(result) + "\n"
    if not matches_result(result):
        return {"software_citation.md": text}
    text += (
        f"\n{APA_MARKDOWN}\n\nIn-text citation: {IN_TEXT}.\n\n"
        "Also cite the statistical methods actually used and retain the recorded "
        "model, estimator, settings, and source commit where available. "
        "Citation does not establish validity or remove any inference hold.\n"
    )
    return {"software_citation.md": text, "mfrm_software.bib": BIBTEX, "CITATION.cff": CFF_TEXT}
