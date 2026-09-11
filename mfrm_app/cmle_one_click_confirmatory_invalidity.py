"""Private, fail-closed invalidity adjudication support for SCI-01/SCI-05."""

from __future__ import annotations

import hashlib
from html import escape
import io
import math
from typing import Mapping
import zipfile

import pandas as pd

from mfrm_app.cmle_one_click_confirmatory_gate import wilson_upper_bound
from mfrm_app.cmle_one_click_confirmatory_protocol import REGISTERED_INVALIDITY_CODES


SCHEMA_VERSION = "cmle_confirmatory_invalidity_adjudication_v1"
COMPOSITE_FRACTIONS = (0.0, 0.25, 0.5, 0.75, 1.0)
THRESHOLD = 0.10
SELECTION_COLUMNS = (
    "SchemaVersion", "InvalidityReasonCode", "AdjudicationStatus",
    "TargetPopulationDisposition", "CompositeEventDisposition",
    "PrimaryAnalysisDisposition", "SensitivitySetDisposition",
    "DataUseAuthorityDisposition", "ScientificOwnerRole",
    "StatisticalReviewReference", "DomainReviewReference",
    "EthicsAuthorityReference", "RationaleReference",
    "RecordedBeforeOutcomes", "HumanOutcomesInspectedBeforeDecision",
    "RepositoryAutoClassificationAllowed",
)
BUNDLE_MEMBERS = (
    "BUNDLE_MANIFEST.csv", "START_HERE.md", "invalidity_adjudication_template.csv",
    "invalidity_adjudication_workbench.html", "invalidity_allowed_state_catalog.csv",
    "invalidity_code_taxonomy.csv", "invalidity_composite_fraction_sensitivity.csv",
    "invalidity_owner_review_matrix.csv",
)

_FAMILIES = {
    "none": ("valid_observed", "scientific_owner", "Fixed parent state; no invalidity"),
    "declined_before_start": ("consent", "ethics_authority", "Data use and denominator interpretation require consent authority"),
    "withdrawn_before_response_lock": ("withdrawal", "ethics_authority", "Withdrawal cannot be converted automatically into a dangerous event"),
    "withdrawn_after_response_lock": ("withdrawal", "ethics_authority", "Post-lock data use depends on approved withdrawal policy"),
    "ineligible_preregistered_rule": ("eligibility", "scientific_owner", "Verify the rule predates outcomes and defines the target population"),
    "technical_failure_before_response_lock": ("technology", "operations_owner", "Separate system reliability from response danger"),
    "technical_failure_after_response_lock": ("technology", "operations_owner", "Preserve locked response if authorized and technically recoverable"),
    "assignment_or_version_mismatch": ("protocol_integrity", "operations_owner", "Do not pool incompatible frozen content"),
    "duplicate_study_record": ("record_integrity", "data_governance_owner", "Resolve record identity without double-counting; never an automatic event"),
    "unknown_response_code": ("response_integrity", "statistical_owner", "Frozen code mapping and blinded review required"),
    "incomplete_primary_items": ("item_missingness", "statistical_owner", "Potentially outcome-related missingness requires sensitivity"),
    "moderator_interruption": ("administration", "operations_owner", "Intervention mechanism must remain auditable"),
    "accessibility_barrier": ("accessibility", "accessibility_owner", "Report separately; silent pooling can conceal exclusion"),
    "protocol_deviation_blinded_review": ("protocol_integrity", "operations_owner", "Document the frozen blinded adjudication rule"),
    "other_requires_protocol_amendment": ("unregistered", "protocol_authority", "Cannot be resolved inside the frozen vocabulary"),
}


def _truth(value: object) -> bool:
    return value if isinstance(value, bool) else str(value).strip().lower() in {"true", "1", "yes"}


def _csv(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(index=False, float_format="%.17g").encode()


def _sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def build_invalidity_code_taxonomy(codebook: pd.DataFrame) -> pd.DataFrame:
    required = {"InvalidityReasonCode", "Description", "OutcomeTiming", "MayBeOutcomeRelated", "RequiresProtocolAmendment"}
    if required - set(codebook.columns):
        raise ValueError("Invalid parent codebook schema")
    if tuple(codebook["InvalidityReasonCode"]) != REGISTERED_INVALIDITY_CODES:
        raise ValueError("Invalidity code order/content mismatch")
    rows = []
    for _, source in codebook.iterrows():
        code = str(source["InvalidityReasonCode"]); family, owner, guardrail = _FAMILIES[code]
        rows.append({"SchemaVersion": SCHEMA_VERSION, "DisplayOrder": len(rows) + 1,
            "InvalidityReasonCode": code, "Description": source["Description"],
            "MechanismFamily": family, "PrimaryOwnerRole": owner,
            "OutcomeTiming": source["OutcomeTiming"], "MayBeOutcomeRelated": bool(source["MayBeOutcomeRelated"]),
            "RequiresProtocolAmendment": bool(source["RequiresProtocolAmendment"]),
            "HardGuardrail": guardrail, "AutomaticallyDangerous": False,
            "AutomaticallySafe": code == "none", "RepositoryRecommendation": "none"})
    return pd.DataFrame(rows)


def build_allowed_state_catalog() -> pd.DataFrame:
    states = {
        "TargetPopulationDisposition": ("in_target", "outside_target_by_preregistered_rule", "authority_withheld", "unresolved"),
        "CompositeEventDisposition": ("dangerous_composite_event", "noncomposite_event", "sensitivity_only", "not_applicable", "unresolved"),
        "PrimaryAnalysisDisposition": ("include_observed", "exclude_preregistered", "authority_withheld", "pending"),
        "SensitivitySetDisposition": ("separate_stratum_required", "code_specific_tipping", "all_safe_all_dangerous_bounds", "not_applicable", "unresolved"),
        "DataUseAuthorityDisposition": ("authorized_use", "authority_excluded", "not_applicable", "unresolved"),
    }
    rows = []
    for dimension, values in states.items():
        for order, value in enumerate(values, 1):
            rows.append({"SchemaVersion": SCHEMA_VERSION, "Dimension": dimension,
                "DisplayOrder": order, "State": value, "DisplayOrderIsRank": False,
                "RepositoryRecommended": False})
    return pd.DataFrame(rows)


def build_owner_review_matrix(taxonomy: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in taxonomy.itertuples():
        roles = {row.PrimaryOwnerRole, "statistical_owner", "scientific_owner"}
        if row.MechanismFamily in {"consent", "withdrawal"}: roles.add("ethics_authority")
        if row.MechanismFamily == "accessibility": roles.add("accessibility_owner")
        for role in sorted(roles):
            rows.append({"SchemaVersion": SCHEMA_VERSION, "InvalidityReasonCode": row.InvalidityReasonCode,
                "OwnerRole": role, "ReviewRequired": row.InvalidityReasonCode != "none",
                "ReviewCompleted": False, "EvidenceReference": "", "RepositoryMaySelfApprove": False})
    return pd.DataFrame(rows)


def build_adjudication_template() -> pd.DataFrame:
    rows = []
    for code in REGISTERED_INVALIDITY_CODES:
        fixed = code == "none"
        rows.append({"SchemaVersion": SCHEMA_VERSION, "InvalidityReasonCode": code,
            "AdjudicationStatus": "fixed_parent_contract" if fixed else "unresolved",
            "TargetPopulationDisposition": "in_target" if fixed else "",
            "CompositeEventDisposition": "not_applicable" if fixed else "",
            "PrimaryAnalysisDisposition": "include_observed" if fixed else "",
            "SensitivitySetDisposition": "not_applicable" if fixed else "",
            "DataUseAuthorityDisposition": "not_applicable" if fixed else "",
            "ScientificOwnerRole": "repository_contract" if fixed else "",
            "StatisticalReviewReference": "fixed_parent_contract" if fixed else "",
            "DomainReviewReference": "fixed_parent_contract" if fixed else "",
            "EthicsAuthorityReference": "not_applicable" if fixed else "",
            "RationaleReference": "fixed_parent_contract" if fixed else "",
            "RecordedBeforeOutcomes": fixed, "HumanOutcomesInspectedBeforeDecision": False,
            "RepositoryAutoClassificationAllowed": False})
    return pd.DataFrame(rows, columns=SELECTION_COLUMNS)


def validate_adjudication(frame: pd.DataFrame) -> dict[str, object]:
    failures = []
    if tuple(frame.columns) != SELECTION_COLUMNS or len(frame) != 15:
        return {"valid": False, "complete": False, "failure_codes": ("schema_or_row_count_mismatch",), "RecruitmentReady": False}
    if tuple(frame["InvalidityReasonCode"].astype(str)) != REGISTERED_INVALIDITY_CODES:
        failures.append("code_identity_or_order_mismatch")
    allowed = build_allowed_state_catalog().groupby("Dimension")["State"].apply(set).to_dict()
    resolved_count = 0
    for _, row in frame.iterrows():
        code, status = str(row["InvalidityReasonCode"]), str(row["AdjudicationStatus"])
        if _truth(row["RepositoryAutoClassificationAllowed"]): failures.append("repository_auto_classification_forbidden")
        if code == "none":
            expected = ("fixed_parent_contract", "in_target", "not_applicable", "include_observed", "not_applicable", "not_applicable")
            actual = tuple(str(row[c]) for c in ("AdjudicationStatus", *SELECTION_COLUMNS[3:8]))
            if actual != expected: failures.append("none_fixed_state_mutated")
            continue
        dimensions = SELECTION_COLUMNS[3:8]
        if status == "unresolved":
            if any(str(row[c]).strip() for c in dimensions): failures.append("unresolved_row_contains_disposition")
            continue
        if code == "other_requires_protocol_amendment":
            if status != "amendment_required": failures.append("other_reason_requires_amendment_status")
            continue
        if status != "resolved_prospectively": failures.append("unknown_adjudication_status"); continue
        resolved_count += 1
        for dimension in dimensions:
            if str(row[dimension]) not in allowed[dimension] or str(row[dimension]) in {"unresolved", "pending"}:
                failures.append("missing_or_invalid_disposition")
        for reference in ("ScientificOwnerRole", "StatisticalReviewReference", "DomainReviewReference", "RationaleReference"):
            if not str(row[reference]).strip(): failures.append("missing_owner_or_review_reference")
        if not _truth(row["RecordedBeforeOutcomes"]) or _truth(row["HumanOutcomesInspectedBeforeDecision"]):
            failures.append("adjudication_not_prospectively_blinded")
        family = _FAMILIES[code][0]
        if family in {"consent", "withdrawal"}:
            if str(row["DataUseAuthorityDisposition"]) not in {"authorized_use", "authority_excluded"} or not str(row["EthicsAuthorityReference"]).strip():
                failures.append("missing_ethics_data_use_authority")
        if code == "duplicate_study_record" and str(row["CompositeEventDisposition"]) == "dangerous_composite_event":
            failures.append("duplicate_cannot_be_dangerous_event")
        if code == "accessibility_barrier":
            if "accessibility" not in str(row["ScientificOwnerRole"]).lower(): failures.append("missing_accessibility_owner")
            if str(row["SensitivitySetDisposition"]) != "separate_stratum_required": failures.append("accessibility_silent_pooling_forbidden")
    complete = resolved_count == 13 and str(frame.iloc[-1]["AdjudicationStatus"]) == "amendment_required" and not failures
    return {"valid": not failures, "complete": complete, "resolved_code_count": resolved_count,
        "failure_codes": tuple(dict.fromkeys(failures)), "RecruitmentReady": False}


def build_composite_fraction_sensitivity(surface: pd.DataFrame) -> pd.DataFrame:
    required = {"DangerousErrors", "ValidEligibleN", "InvalidEligibleN"}
    if required - set(surface.columns): raise ValueError("Missing source count columns")
    rows = []
    for source_index, source in surface.reset_index(drop=True).iterrows():
        d, v, i = (int(source[c]) for c in ("DangerousErrors", "ValidEligibleN", "InvalidEligibleN"))
        total = v + i
        for fraction in COMPOSITE_FRACTIONS:
            classified = int(math.floor(i * fraction + 0.5))
            errors = d + classified
            ucb = wilson_upper_bound(errors, total, confidence=.95)
            rows.append({"SchemaVersion": SCHEMA_VERSION, "SourceScenarioRow": source_index + 1,
                "DangerousErrorsObservedValid": d, "ValidEligibleN": v, "InvalidEligibleN": i,
                "CompositeFractionOfInvalidRegistered": fraction,
                "InvalidClassifiedCompositeN": classified, "CompositeEventNDiagnostic": errors,
                "DiagnosticDenominatorN": total, "CompositeRiskDiagnosticRaw": errors / total,
                "CompositeWilsonUpper95Raw": ucb, "DecisionThresholdRaw": THRESHOLD,
                "ProjectedDecision": "proxy_pass" if ucb < THRESHOLD else "proxy_fail",
                "AllocationRule": "floor(I*fraction+0.5)", "DecisionUsesDisplayedRounding": False,
                "InfersCodeFrequency": False, "SelectsCompositeDefinition": False,
                "FinalScheduledSlotEstimandCalculated": False, "HumanParticipants": 0})
    return pd.DataFrame(rows)


def content_identity(*frames: pd.DataFrame) -> str:
    return _sha(b"\0".join(_csv(frame) for frame in frames))


def render_html(taxonomy: pd.DataFrame, sensitivity: pd.DataFrame, identity: str) -> str:
    rows = "".join(f"<tr><td>{escape(str(r.InvalidityReasonCode))}</td><td>{escape(str(r.MechanismFamily))}</td><td>{escape(str(r.PrimaryOwnerRole))}</td><td>{escape(str(r.HardGuardrail))}</td></tr>" for r in taxonomy.itertuples())
    counts = sensitivity["ProjectedDecision"].value_counts()
    return f"""<!doctype html><html lang='ja'><head><meta charset='utf-8'><title>Invalidity adjudication</title><style>body{{font-family:system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#172033}}header{{background:#7f1d1d;color:white;padding:1rem;border-radius:12px}}.warn{{background:#fff7ed;border-left:5px solid #ea580c;padding:1rem}}table{{border-collapse:collapse;width:100%;font-size:.88rem}}th,td{{border:1px solid #cbd5e1;padding:.45rem;text-align:left}}code{{word-break:break-all}}</style></head><body><header><h1>BLOCKED / 募集不可 — 無効理由の裁定</h1><p>ResolvedCodes=0 · RecruitmentReady=false · RepositoryRecommendation=none</p></header><p class='warn'>無効は危険事象と同義ではありません。撤回、アクセシビリティ障壁、技術障害、重複を自動分類しません。表示された感度は診断用で、実際のコード頻度や予定枠推定対象を表しません。</p><p>Invalid is not synonymous with dangerous. External scientific, statistical, domain, accessibility, and ethics ownership must resolve the blank rows before outcome inspection.</p><p>Diagnostic rows={len(sensitivity)}; proxy_pass={int(counts.get('proxy_pass',0))}; proxy_fail={int(counts.get('proxy_fail',0))}. Raw strict &lt; 0.10 decisions only; displayed rounding is never used.</p><table><thead><tr><th>Code</th><th>Family</th><th>Primary owner</th><th>Guardrail</th></tr></thead><tbody>{rows}</tbody></table><footer><p>InvalidityContentSHA256: <code>{identity}</code></p><p>HumanParticipants=0 · ConfirmatoryOutcomesAvailable=false · PublicSurfaceEnabled=false</p></footer></body></html>"""


def validate_html(html: str, identity: str) -> dict[str, object]:
    failures = []
    for value in ("BLOCKED", "募集不可", "ResolvedCodes=0", "RecruitmentReady=false", identity, *REGISTERED_INVALIDITY_CODES):
        if value not in html: failures.append("missing_required_content")
    lower = html.lower()
    if any(x in lower for x in ("<form", "<input", "<script", "http://", "https://")): failures.append("writable_or_external_content")
    if any(x in lower for x in ("recommended option", "推奨案", "best option")): failures.append("recommendation_language")
    unique = tuple(dict.fromkeys(failures)); return {"valid": not unique, "failure_codes": unique}


def _zip(payloads: Mapping[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name in sorted(payloads):
            info = zipfile.ZipInfo(name, (1980,1,1,0,0,0)); info.compress_type = zipfile.ZIP_DEFLATED; info.external_attr = 0o100644 << 16; info.create_system = 3
            archive.writestr(info, payloads[name], compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    return buffer.getvalue()


def build_private_bundle(codebook: pd.DataFrame, surface: pd.DataFrame) -> dict[str, object]:
    taxonomy = build_invalidity_code_taxonomy(codebook); states = build_allowed_state_catalog()
    template = build_adjudication_template(); owners = build_owner_review_matrix(taxonomy)
    sensitivity = build_composite_fraction_sensitivity(surface)
    identity = content_identity(taxonomy, states, template, owners, sensitivity); html = render_html(taxonomy, sensitivity, identity)
    start = f"# BLOCKED — private invalidity adjudication\n\nResolvedCodes=0. AdjudicationComplete=false. RecruitmentReady=false. PublicSurfaceEnabled=false.\n\nInvalidity is not automatically dangerous. 募集不可。撤回・アクセス障壁・技術障害・重複を自動分類しません。\n\nInvalidityContentSHA256: `{identity}`\n"
    payloads = {"START_HERE.md": start.encode(), "invalidity_code_taxonomy.csv": _csv(taxonomy),
        "invalidity_allowed_state_catalog.csv": _csv(states), "invalidity_adjudication_template.csv": _csv(template),
        "invalidity_owner_review_matrix.csv": _csv(owners), "invalidity_composite_fraction_sensitivity.csv": _csv(sensitivity),
        "invalidity_adjudication_workbench.html": html.encode()}
    manifest = pd.DataFrame([{"Member": n, "ByteCount": len(p), "SHA256": _sha(p), "HumanRows": 0, "ContainsAdjudication": False, "PublicArtifact": False} for n,p in sorted(payloads.items())])
    payloads["BUNDLE_MANIFEST.csv"] = _csv(manifest)
    return {"taxonomy":taxonomy,"states":states,"template":template,"owners":owners,"sensitivity":sensitivity,"identity":identity,"html":html,"manifest":manifest,"bundle_bytes":_zip(payloads)}


def validate_private_bundle(bundle: bytes) -> dict[str, object]:
    failures=[]
    try:
        with zipfile.ZipFile(io.BytesIO(bundle),"r") as archive:
            infos=archive.infolist(); names=tuple(x.filename for x in infos)
            if names != tuple(sorted(BUNDLE_MEMBERS)): failures.append("membership_or_order")
            for info in infos:
                if info.date_time != (1980,1,1,0,0,0): failures.append("timestamp")
                if (info.external_attr>>16)&0o777 != 0o644: failures.append("permissions")
            manifest=pd.read_csv(io.BytesIO(archive.read("BUNDLE_MANIFEST.csv")),dtype=str)
            if set(manifest["Member"]) != set(names)-{"BUNDLE_MANIFEST.csv"}: failures.append("manifest_members")
            for _,row in manifest.iterrows():
                payload=archive.read(str(row["Member"]))
                if len(payload)!=int(row["ByteCount"]): failures.append("manifest_size")
                if _sha(payload)!=str(row["SHA256"]): failures.append("manifest_hash")
                if int(row["HumanRows"]) or _truth(row["ContainsAdjudication"]): failures.append("human_or_adjudication")
                if _truth(row["PublicArtifact"]): failures.append("public_artifact")
            start=archive.read("START_HERE.md").decode(); html=archive.read("invalidity_adjudication_workbench.html").decode(); marker="InvalidityContentSHA256: `"; identity=start.split(marker,1)[1].split("`",1)[0]
            if not validate_html(html,identity)["valid"]: failures.append("invalid_html")
            if any(x not in start for x in ("BLOCKED","募集不可","ResolvedCodes=0","RecruitmentReady=false")): failures.append("first_read")
    except (zipfile.BadZipFile,KeyError,ValueError,UnicodeDecodeError,IndexError) as error: failures.append(f"invalid_zip:{type(error).__name__}")
    unique=tuple(dict.fromkeys(failures)); return {"valid":not unique,"failure_codes":unique,"bundle_sha256":_sha(bundle)}
