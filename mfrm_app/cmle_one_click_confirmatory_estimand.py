"""Private, non-recommending SCI-01 estimand decision support."""

from __future__ import annotations

import hashlib
from html import escape
import io
from typing import Mapping
import zipfile

import pandas as pd


SCI01_SCHEMA_VERSION = "cmle_confirmatory_sci01_estimand_memo_v1"
SCI01_OPTION_IDS = (
    "SCI-01-A_observed_valid_estimand",
    "SCI-01-B_scheduled_slot_composite",
    "SCI-01-C_dual_estimand_with_bounds",
)
SCI01_SELECTION_COLUMNS = (
    "SchemaVersion", "DecisionId", "DecisionStatus", "SelectedOptionId",
    "ScientificOwnerRole", "AssumptionSetReference", "RationaleReference",
    "OwnerAttestationReference", "InvalidityCodeMapReference",
    "DecisionRecordedBeforeOutcomes", "HumanOutcomesInspectedBeforeDecision",
    "RepositoryGeneratedSelectionAllowed",
)
SCI01_BUNDLE_MEMBERS = (
    "BUNDLE_MANIFEST.csv", "START_HERE.md",
    "sci01_downstream_decision_impacts.csv", "sci01_estimand_comparison.html",
    "sci01_estimand_definitions.csv", "sci01_same_data_projection.csv",
    "sci01_selection_template.csv",
)


def _truth(value: object) -> bool:
    return value if isinstance(value, bool) else str(value).strip().lower() in {"true", "1", "yes"}


def _csv_bytes(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(index=False, float_format="%.17g").encode("utf-8")


def _sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def build_sci01_estimand_definitions() -> pd.DataFrame:
    """Return the three registered alternatives in display order, never rank order."""
    rows = (
        {
            "OptionId": SCI01_OPTION_IDS[0], "DisplayOrder": 1,
            "ShortLabel": "Observed-valid estimand", "TargetPopulation": "Valid eligible records",
            "PointEstimand": "D / V", "UncertaintyRule": "Wilson upper 95% bound for D of V",
            "InvalidOutcomeTreatment": "Excluded from point denominator; retained in denominator audit",
            "IdentificationRequirement": "Validity is noninformative for the scoped claim, or claim is explicitly conditional on validity",
            "AllowedClaim": "Dangerous-error risk among valid eligible records only",
            "ForbiddenClaim": "Scheduled-user, invited-user, or all-eligible risk without further assumptions",
            "PrimaryRisk": "Informative invalidity can make the conditional risk falsely reassuring",
            "RequiredProspectiveMaterial": "Target-population definition and missingness justification",
            "IsDirectlyCalculableFromFrozenSurface": True,
        },
        {
            "OptionId": SCI01_OPTION_IDS[1], "DisplayOrder": 2,
            "ShortLabel": "Scheduled-slot composite", "TargetPopulation": "Prospectively defined scheduled slots",
            "PointEstimand": "(D + prespecified composite failures) / scheduled slots",
            "UncertaintyRule": "Must be registered after the composite event and denominator are fixed",
            "InvalidOutcomeTreatment": "Only prospectively adjudicated codes enter the composite",
            "IdentificationRequirement": "Complete scheduled-slot ledger and fixed invalidity-code-to-event map",
            "AllowedClaim": "Composite event risk among the registered scheduled-slot population",
            "ForbiddenClaim": "Interpreting every invalidity as a dangerous error or using V+I as scheduled slots without proof",
            "PrimaryRisk": "Conflates technology, access, withdrawal, comprehension, and safety mechanisms",
            "RequiredProspectiveMaterial": "Scheduled-slot denominator and stakeholder-justified composite event map",
            "IsDirectlyCalculableFromFrozenSurface": False,
        },
        {
            "OptionId": SCI01_OPTION_IDS[2], "DisplayOrder": 3,
            "ShortLabel": "Dual estimand with bounds", "TargetPopulation": "Valid eligible records plus eligible records with unresolved outcomes",
            "PointEstimand": "D / V, accompanied by [D/(V+I), (D+I)/(V+I)]",
            "UncertaintyRule": "Observed and all-invalid-dangerous Wilson upper 95% bounds shown separately",
            "InvalidOutcomeTreatment": "No untestable allocation; both all-safe and all-dangerous endpoints reported",
            "IdentificationRequirement": "Exact D, V, and I counts; a separate prospective governing decision rule",
            "AllowedClaim": "Observed-valid risk and a transparent nonparametric uncertainty region",
            "ForbiddenClaim": "A single identified scheduled-population risk or automatic pass from the observed-valid result",
            "PrimaryRisk": "More complex communication; bounds alone do not choose a governing conclusion",
            "RequiredProspectiveMaterial": "Dual-estimand protocol, invalidity scope, and adjudication rule",
            "IsDirectlyCalculableFromFrozenSurface": True,
        },
    )
    frame = pd.DataFrame(rows)
    frame.insert(0, "SchemaVersion", SCI01_SCHEMA_VERSION)
    frame["DisplayOrderIsRank"] = False
    frame["RepositoryRecommendation"] = "none"
    frame["AutoSelectable"] = False
    return frame


def project_sci01_estimands(surface: pd.DataFrame) -> pd.DataFrame:
    """Project identical integer-count scenarios under all three alternatives."""
    required = {
        "DangerousErrors", "ValidEligibleN", "InvalidEligibleN", "TotalEligibleIfObserved",
        "ObservedValidDangerousProportionRaw", "AllInvalidSafeLowerRiskRaw",
        "AllInvalidDangerousUpperRiskRaw", "ObservedValidWilsonUpper95Raw",
        "AllInvalidDangerousWilsonUpper95Raw", "DecisionThresholdRaw",
    }
    missing = required - set(surface.columns)
    if missing:
        raise ValueError(f"Missing sensitivity columns: {sorted(missing)}")
    rows: list[dict[str, object]] = []
    for source_index, source in surface.reset_index(drop=True).iterrows():
        d, v, i = (int(source[name]) for name in ("DangerousErrors", "ValidEligibleN", "InvalidEligibleN"))
        threshold = float(source["DecisionThresholdRaw"])
        observed_ucb = float(source["ObservedValidWilsonUpper95Raw"])
        worst_ucb = float(source["AllInvalidDangerousWilsonUpper95Raw"])
        base = {
            "SchemaVersion": SCI01_SCHEMA_VERSION, "SourceScenarioRow": source_index + 1,
            "DangerousErrors": d, "ValidEligibleN": v, "InvalidEligibleN": i,
            "EligibleObservedOrInvalidN": v + i, "DecisionThresholdRaw": threshold,
            "DecisionOperator": "<", "DecisionUsesDisplayedRounding": False,
            "SampleSizeSelected": False, "HumanParticipants": 0,
        }
        rows.append({**base, "OptionId": SCI01_OPTION_IDS[0],
            "ReportedPointRaw": d / v, "ReportedLowerRaw": d / v, "ReportedUpperRaw": d / v,
            "GoverningWilsonUpper95Raw": observed_ucb,
            "ProjectedDecision": "pass" if observed_ucb < threshold else "fail",
            "ProjectionRole": "direct_conditional_estimand",
            "FinalScheduledSlotEstimandCalculated": False,
            "InterpretationLimit": "Conditional on valid eligible records"})
        rows.append({**base, "OptionId": SCI01_OPTION_IDS[1],
            "ReportedPointRaw": (d + i) / (v + i),
            "ReportedLowerRaw": d / (v + i), "ReportedUpperRaw": (d + i) / (v + i),
            "GoverningWilsonUpper95Raw": worst_ucb,
            "ProjectedDecision": "proxy_pass" if worst_ucb < threshold else "proxy_fail",
            "ProjectionRole": "all_invalid_dangerous_diagnostic_proxy",
            "FinalScheduledSlotEstimandCalculated": False,
            "InterpretationLimit": "Not the scheduled-slot estimand: scheduled denominator and code map are absent"})
        if worst_ucb < threshold:
            dual_decision = "robust_pass"
        elif observed_ucb < threshold:
            dual_decision = "observed_pass_not_robust"
        else:
            dual_decision = "observed_gate_fail"
        rows.append({**base, "OptionId": SCI01_OPTION_IDS[2],
            "ReportedPointRaw": d / v,
            "ReportedLowerRaw": d / (v + i), "ReportedUpperRaw": (d + i) / (v + i),
            "GoverningWilsonUpper95Raw": worst_ucb,
            "ProjectedDecision": dual_decision,
            "ProjectionRole": "partial_identification_bounds",
            "FinalScheduledSlotEstimandCalculated": False,
            "InterpretationLimit": "Bounds do not supply the still-unresolved governing rule"})
    return pd.DataFrame(rows)


def build_sci01_downstream_impacts() -> pd.DataFrame:
    """Expose, without scoring, which unresolved decisions SCI-01 changes."""
    impacts = {
        "SCI-04": "Multiplicity family depends on whether one, composite, or dual claims govern",
        "SCI-05": "Missingness rule is part of, or must remain coherent with, the estimand",
        "SCI-06": "Cluster target and uncertainty must match the selected population and event",
        "SCI-07": "Mechanism pooling changes event interpretation and cell family",
        "SCI-08": "Language-specific claim and denominator must use the same estimand definition",
        "SCI-09": "Accessibility exclusions can change the target population and invalidate comparisons",
        "SCI-02": "Sample size cannot be selected before event, denominator, and governing rule are fixed",
        "SCI-03": "Recruitment inflation needs the valid/scheduled denominator target",
        "OPS-01": "Blinded monitoring metrics depend on allowed denominator and invalidity information",
        "GOV-01": "Data collection must retain the fields needed for the chosen denominator and audit",
        "ETH-01": "Reviewed scope must match target population and treatment of withdrawal/invalidity",
        "ETH-02": "Consent and withdrawal handling constrain which records may enter any composite",
    }
    rows = []
    for option in SCI01_OPTION_IDS:
        for order, (decision, reason) in enumerate(impacts.items(), start=1):
            rows.append({"SchemaVersion": SCI01_SCHEMA_VERSION, "OptionId": option,
                "DownstreamOrder": order, "DownstreamDecisionId": decision,
                "ImpactMechanism": reason, "ImpactMagnitudeEstimated": False,
                "DirectionKnown": False, "DecisionResolved": False})
    return pd.DataFrame(rows)


def build_sci01_selection_template() -> pd.DataFrame:
    return pd.DataFrame([{
        "SchemaVersion": SCI01_SCHEMA_VERSION, "DecisionId": "SCI-01",
        "DecisionStatus": "unresolved", "SelectedOptionId": "",
        "ScientificOwnerRole": "", "AssumptionSetReference": "",
        "RationaleReference": "", "OwnerAttestationReference": "",
        "InvalidityCodeMapReference": "", "DecisionRecordedBeforeOutcomes": False,
        "HumanOutcomesInspectedBeforeDecision": False,
        "RepositoryGeneratedSelectionAllowed": False,
    }], columns=SCI01_SELECTION_COLUMNS)


def validate_sci01_selection(selection: pd.DataFrame) -> dict[str, object]:
    failures: list[str] = []
    if tuple(selection.columns) != SCI01_SELECTION_COLUMNS or len(selection) != 1:
        return {"valid": False, "complete": False, "failure_codes": ("selection_schema_mismatch",), "RecruitmentReady": False}
    row = selection.iloc[0]
    if str(row["SchemaVersion"]) != SCI01_SCHEMA_VERSION or str(row["DecisionId"]) != "SCI-01":
        failures.append("selection_identity_mismatch")
    status, selected = str(row["DecisionStatus"]).strip(), str(row["SelectedOptionId"]).strip()
    if _truth(row["RepositoryGeneratedSelectionAllowed"]): failures.append("repository_selection_forbidden")
    if status == "unresolved":
        if selected or any(str(row[col]).strip() for col in ("ScientificOwnerRole", "AssumptionSetReference", "RationaleReference", "OwnerAttestationReference", "InvalidityCodeMapReference")):
            failures.append("unresolved_row_contains_selection_material")
        complete = False
    elif status == "resolved_prospectively":
        if selected not in SCI01_OPTION_IDS: failures.append("unknown_option")
        required = ("ScientificOwnerRole", "AssumptionSetReference", "RationaleReference", "OwnerAttestationReference")
        if any(not str(row[col]).strip() for col in required): failures.append("missing_owner_or_evidence_reference")
        if selected in SCI01_OPTION_IDS[1:] and not str(row["InvalidityCodeMapReference"]).strip():
            failures.append("missing_invalidity_code_map")
        if "repository" in str(row["ScientificOwnerRole"]).lower(): failures.append("repository_cannot_be_scientific_owner")
        if not _truth(row["DecisionRecordedBeforeOutcomes"]) or _truth(row["HumanOutcomesInspectedBeforeDecision"]):
            failures.append("decision_not_prospectively_locked")
        complete = not failures
    else:
        failures.append("unknown_decision_status"); complete = False
    return {"valid": not failures, "complete": complete,
        "failure_codes": tuple(dict.fromkeys(failures)), "SelectedOptionId": selected,
        "RecruitmentReady": False}


def compute_sci01_content_identity(definitions: pd.DataFrame, projection: pd.DataFrame,
                                   impacts: pd.DataFrame, selection: pd.DataFrame) -> str:
    payload = b"\0".join(_csv_bytes(frame) for frame in (definitions, projection, impacts, selection))
    return _sha(payload)


def render_sci01_estimand_html(definitions: pd.DataFrame, projection: pd.DataFrame,
                               impacts: pd.DataFrame, content_sha256: str) -> str:
    cards = []
    for _, row in definitions.sort_values("DisplayOrder").iterrows():
        option = str(row["OptionId"])
        subset = projection.loc[projection["OptionId"].eq(option)]
        counts = subset["ProjectedDecision"].value_counts().sort_index()
        decision_counts = ", ".join(f"{escape(str(k))}={int(v)}" for k, v in counts.items())
        cards.append(f"<article><h2>{escape(str(row['ShortLabel']))}</h2><code>{escape(option)}</code>"
            f"<dl><dt>対象 / Target</dt><dd>{escape(str(row['TargetPopulation']))}</dd>"
            f"<dt>式 / Formula</dt><dd>{escape(str(row['PointEstimand']))}</dd>"
            f"<dt>許容主張 / Allowed</dt><dd>{escape(str(row['AllowedClaim']))}</dd>"
            f"<dt>禁止主張 / Forbidden</dt><dd>{escape(str(row['ForbiddenClaim']))}</dd>"
            f"<dt>主要リスク / Risk</dt><dd>{escape(str(row['PrimaryRisk']))}</dd></dl>"
            f"<p class='count'>72 registered scenarios: {decision_counts}</p></article>")
    downstream = "".join(f"<tr><td>{escape(str(r.OptionId))}</td><td>{escape(str(r.DownstreamDecisionId))}</td><td>{escape(str(r.ImpactMechanism))}</td></tr>" for r in impacts.itertuples())
    return f"""<!doctype html><html lang='ja'><head><meta charset='utf-8'><title>SCI-01 estimand memo</title>
<style>body{{font-family:system-ui,sans-serif;max-width:1120px;margin:2rem auto;padding:0 1rem;color:#172033}}header{{background:#7f1d1d;color:white;padding:1.2rem;border-radius:12px}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:1rem;margin:1rem 0}}article{{border:1px solid #cbd5e1;border-radius:12px;padding:1rem}}dt{{font-weight:700;margin-top:.7rem}}dd{{margin-left:0}}code{{word-break:break-all}}table{{border-collapse:collapse;width:100%;font-size:.86rem}}th,td{{border:1px solid #cbd5e1;padding:.45rem;text-align:left}}.warning{{background:#fff7ed;border-left:5px solid #ea580c;padding:1rem}}.count{{font-family:ui-monospace,monospace;background:#f1f5f9;padding:.5rem}}</style></head>
<body><header><h1>BLOCKED / 募集不可 — SCI-01 推定対象比較</h1><p>OptionsSelected=0 · RecruitmentReady=false · RepositoryRecommendation=none</p></header>
<p class='warning'>同じデータを3通りに投影した意思決定資料です。表示順は順位ではありません。Option B の数値は全無効例を危険事象とした診断用代理値であり、予定枠分母と無効コード表がないため scheduled-slot estimand そのものではありません。</p>
<p>This private read-only memo compares consequences; it does not select an estimand, verify assumptions, approve ethics, or authorize recruitment.</p>
<div class='grid'>{''.join(cards)}</div>
<h2>下流への影響 / Downstream impacts</h2><table><thead><tr><th>Option</th><th>Decision</th><th>Why it changes</th></tr></thead><tbody>{downstream}</tbody></table>
<p>Integer counts and unrounded binary64 values drive every projected label. The registered operator is strict &lt; 0.10; displayed rounding never drives a decision.</p>
<footer><p>SCI01ContentSHA256: <code>{escape(content_sha256)}</code></p><p>HumanParticipants=0 · ConfirmatoryOutcomesAvailable=false · PublicSurfaceEnabled=false</p></footer></body></html>"""


def validate_sci01_html(html: str, content_sha256: str) -> dict[str, object]:
    failures = []
    for required in ("BLOCKED", "募集不可", "OptionsSelected=0", "RecruitmentReady=false", content_sha256, *SCI01_OPTION_IDS):
        if required not in html: failures.append("missing_required_content")
    lower = html.lower()
    if any(token in lower for token in ("<form", "<input", "<script", "http://", "https://")):
        failures.append("writable_or_external_html_content")
    if any(token in lower for token in ("recommended option", "推奨案", "best option")):
        failures.append("recommendation_language_present")
    unique = tuple(dict.fromkeys(failures))
    return {"valid": not unique, "failure_codes": unique}


def _zip(payloads: Mapping[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name in sorted(payloads):
            info = zipfile.ZipInfo(name, (1980, 1, 1, 0, 0, 0)); info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16; info.create_system = 3
            archive.writestr(info, payloads[name], compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    return buffer.getvalue()


def build_private_sci01_bundle(surface: pd.DataFrame) -> dict[str, object]:
    definitions = build_sci01_estimand_definitions(); projection = project_sci01_estimands(surface)
    impacts = build_sci01_downstream_impacts(); selection = build_sci01_selection_template()
    selection_status = validate_sci01_selection(selection)
    content_sha = compute_sci01_content_identity(definitions, projection, impacts, selection)
    html = render_sci01_estimand_html(definitions, projection, impacts, content_sha)
    start = f"""# BLOCKED — private SCI-01 estimand comparison

OptionsSelected=0. SCI01Complete=false. RecruitmentReady=false. PublicSurfaceEnabled=false.

The display order is not a rank. Option B rows are diagnostic all-invalid-dangerous proxies, not a calculated scheduled-slot estimand. External scientific ownership is required before outcomes are inspected.

募集不可。推定対象は未選択で、このリポジトリは推奨・倫理承認・募集許可を行いません。

SCI01ContentSHA256: `{content_sha}`
"""
    payloads = {"START_HERE.md": start.encode(),
        "sci01_estimand_definitions.csv": _csv_bytes(definitions),
        "sci01_same_data_projection.csv": _csv_bytes(projection),
        "sci01_downstream_decision_impacts.csv": _csv_bytes(impacts),
        "sci01_selection_template.csv": _csv_bytes(selection),
        "sci01_estimand_comparison.html": html.encode()}
    manifest = pd.DataFrame([{"Member": name, "ByteCount": len(payload), "SHA256": _sha(payload),
        "HumanRows": 0, "ContainsSelection": False, "PublicArtifact": False} for name, payload in sorted(payloads.items())])
    payloads["BUNDLE_MANIFEST.csv"] = _csv_bytes(manifest)
    return {"definitions": definitions, "projection": projection, "impacts": impacts,
        "selection": selection, "selection_status": selection_status, "content_sha256": content_sha,
        "html": html, "manifest": manifest, "bundle_bytes": _zip(payloads)}


def validate_private_sci01_bundle(bundle_bytes: bytes) -> dict[str, object]:
    failures = []
    try:
        with zipfile.ZipFile(io.BytesIO(bundle_bytes), "r") as archive:
            infos = archive.infolist(); names = tuple(info.filename for info in infos)
            if names != tuple(sorted(SCI01_BUNDLE_MEMBERS)): failures.append("member_order_or_membership_mismatch")
            for info in infos:
                if info.date_time != (1980, 1, 1, 0, 0, 0): failures.append("non_deterministic_timestamp")
                if (info.external_attr >> 16) & 0o777 != 0o644: failures.append("non_deterministic_permissions")
            manifest = pd.read_csv(io.BytesIO(archive.read("BUNDLE_MANIFEST.csv")), dtype=str)
            if set(manifest["Member"]) != set(names) - {"BUNDLE_MANIFEST.csv"}: failures.append("manifest_member_mismatch")
            for _, row in manifest.iterrows():
                payload = archive.read(str(row["Member"]))
                if len(payload) != int(row["ByteCount"]): failures.append("manifest_byte_count_mismatch")
                if _sha(payload) != str(row["SHA256"]): failures.append("manifest_hash_mismatch")
                if int(row["HumanRows"]) or _truth(row["ContainsSelection"]): failures.append("human_rows_or_selection_present")
                if _truth(row["PublicArtifact"]): failures.append("public_artifact_present")
            start = archive.read("START_HERE.md").decode(); html = archive.read("sci01_estimand_comparison.html").decode()
            marker = "SCI01ContentSHA256: `"; identity = start.split(marker, 1)[1].split("`", 1)[0] if marker in start else ""
            if any(value not in start for value in ("BLOCKED", "募集不可", "OptionsSelected=0", "RecruitmentReady=false")): failures.append("first_read_not_blocked")
            if not validate_sci01_html(html, identity)["valid"]: failures.append("invalid_sci01_html")
    except (zipfile.BadZipFile, KeyError, ValueError, UnicodeDecodeError, IndexError) as error:
        failures.append(f"invalid_zip:{type(error).__name__}")
    unique = tuple(dict.fromkeys(failures))
    return {"valid": not unique, "failure_codes": unique, "bundle_sha256": _sha(bundle_bytes)}


def sci01_human_gate_status(selection_status: Mapping[str, object]) -> pd.DataFrame:
    return pd.DataFrame([{"SchemaVersion": SCI01_SCHEMA_VERSION, "HumanParticipants": 0,
        "HumanStudyStatus": "not_started_no_human_data", "OptionsSelected": 0,
        "SCI01Complete": bool(selection_status["complete"]), "SubstantiveEvidenceVerified": False,
        "RecruitmentReady": False, "SampleSizeSelected": False,
        "ConfirmatoryOutcomesAvailable": False, "PublicSurfaceEnabled": False}])
