#!/usr/bin/env python3
"""Build the registered private, non-recommending decision workbench."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

_MPL_CONFIG = Path(tempfile.gettempdir()) / "mfrm_app_matplotlib"
_MPL_CONFIG.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CONFIG))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.cmle_one_click_confirmatory_workbench import (  # noqa: E402
    DECISION_IDS,
    PREREQUISITES,
    WORKBENCH_BUNDLE_MEMBERS,
    build_decision_dependency_table,
    build_decision_option_catalog,
    build_decision_selection_template,
    build_private_decision_workbench_bundle,
    compute_workbench_content_identity,
    decision_workbench_human_gate_status,
    validate_decision_selections,
    validate_dependency_table,
    validate_private_decision_workbench_bundle,
    validate_workbench_html,
)

PLAN = ROOT / "validation/cmle_one_click_confirmatory_decision_workbench_plan_20260810.json"
AMENDMENT = ROOT / "validation/cmle_one_click_confirmatory_decision_workbench_amendment_20260810.json"
OUTPUT = ROOT / "validation/cmle_one_click_confirmatory_decision_workbench_20260810"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, float_format="%.17g")


def json_safe(value):
    if isinstance(value, dict): return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [json_safe(v) for v in value]
    if isinstance(value, np.generic): return json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def validate_registration() -> tuple[dict[str, object], dict[str, object]]:
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    amendment = json.loads(AMENDMENT.read_text(encoding="utf-8"))
    if plan.get("study_id") != "cmle_one_click_confirmatory_decision_workbench_v1":
        raise ValueError("Unexpected workbench study identity.")
    if amendment.get("study_id") != plan["study_id"] or amendment.get("parent_plan_sha256") != sha256_file(PLAN):
        raise ValueError("Workbench amendment identity failed.")
    mismatches = [relative for relative, expected in plan["parent_identity"].items() if not (ROOT / relative).is_file() or sha256_file(ROOT / relative) != expected]
    if mismatches: raise ValueError(f"Workbench parent identity failed: {mismatches}")
    return plan, amendment


def registration_audit(plan: dict[str, object]) -> pd.DataFrame:
    return pd.DataFrame([{"Artifact": relative, "ExpectedSHA256": expected, "ActualSHA256": sha256_file(ROOT / relative) if (ROOT / relative).is_file() else "", "Passed": (ROOT / relative).is_file() and sha256_file(ROOT / relative) == expected} for relative, expected in plan["parent_identity"].items()])


def catalog_audit(catalog: pd.DataFrame, plan: dict[str, object]) -> pd.DataFrame:
    registered = plan["option_catalog_contract"]["option_ids"]
    rows = []
    for decision_id in DECISION_IDS:
        frame = catalog.loc[catalog["DecisionId"].eq(decision_id)].sort_values("OptionOrderWithinDecision")
        actual = tuple(frame["OptionId"])
        expected = tuple(registered[decision_id])
        tradeoffs = frame[["Description", "Strength", "PrimaryRisk", "RequiredEvidence", "SampleSizeImpact", "EstimandImpact", "ParticipantBurdenImpact", "AccessibilityImpact", "PublicClaimImpact"]].astype(str).apply(lambda col: col.str.strip().ne("").all()).all()
        rows.append({"DecisionId": decision_id, "OptionCount": len(frame), "ExpectedOptionIds": "|".join(expected), "ActualOptionIds": "|".join(actual), "TradeoffFieldsComplete": bool(tradeoffs), "RecommendedCount": int(~frame["RecommendationStatus"].eq("not_recommended_by_repository").sum()) if False else int(frame["RecommendationStatus"].ne("not_recommended_by_repository").sum()), "AutoSelectableCount": int(frame["AutoSelectable"].sum()), "RankedCount": int(frame["OptionOrderIsRank"].sum()), "Passed": bool(len(frame) == 3 and actual == expected and tradeoffs and not frame["AutoSelectable"].any() and not frame["OptionOrderIsRank"].any() and frame["RecommendationStatus"].eq("not_recommended_by_repository").all())})
    return pd.DataFrame(rows)


def _complete_synthetic(catalog: pd.DataFrame, dependencies: pd.DataFrame) -> pd.DataFrame:
    selections = build_decision_selection_template()
    for decision_id in dependencies.sort_values("TopologicalOrder")["DecisionId"]:
        option = catalog.loc[catalog["DecisionId"].eq(decision_id)].iloc[0]
        mask = selections["DecisionId"].eq(decision_id)
        selections.loc[mask, "DecisionStatus"] = "resolved_prospectively"
        selections.loc[mask, "SelectedOptionId"] = option["OptionId"]
        selections.loc[mask, ["AssumptionSetReference", "RationaleReference", "OwnerAttestationReference"]] = "SYNTHETIC"
        selections.loc[mask, "DecisionRecordedBeforeOutcomes"] = True
        if bool(option["RequiresNumericInput"]):
            selections.loc[mask, "NumericValue"] = "100"; selections.loc[mask, "NumericUnit"] = "synthetic_units"
        if bool(option["ExternalAuthorityRequired"]):
            selections.loc[mask, "ExternalEvidenceReference"] = "SYNTHETIC-EXTERNAL"; selections.loc[mask, "EvidenceSource"] = "synthetic_external_fixture"
    return selections


def selection_attack_audit(catalog: pd.DataFrame, dependencies: pd.DataFrame) -> pd.DataFrame:
    default = build_decision_selection_template()
    hidden = default.copy(); hidden.loc[0, "SelectedOptionId"] = "SCI-01-A_observed_valid_estimand"
    child = default.copy(); mask = child["DecisionId"].eq("SCI-04"); child.loc[mask, "DecisionStatus"] = "resolved_prospectively"; child.loc[mask, "SelectedOptionId"] = "SCI-04-A_primary_all_cells_sensitivity_only"; child.loc[mask, ["AssumptionSetReference", "RationaleReference", "OwnerAttestationReference"]] = "SYNTHETIC"; child.loc[mask, "DecisionRecordedBeforeOutcomes"] = True
    complete = _complete_synthetic(catalog, dependencies)
    numeric = complete.copy(); numeric.loc[numeric["DecisionId"].eq("SCI-02"), ["NumericValue", "NumericUnit"]] = ""
    external = complete.copy(); external.loc[external["DecisionId"].eq("ETH-01"), "EvidenceSource"] = "repository_generated_fake"
    prospective = complete.copy(); prospective.loc[prospective["DecisionId"].eq("SCI-01"), "HumanOutcomesInspectedBeforeDecision"] = True
    cases = (
        ("default_blank_blocked", default, True, False, None),
        ("hidden_selection", hidden, False, False, "unresolved_row_contains_selection"),
        ("child_before_prerequisite", child, False, False, "unresolved_prerequisite"),
        ("numeric_missing", numeric, False, False, "missing_or_invalid_numeric_input"),
        ("repository_external_evidence", external, False, False, "invalid_external_evidence_source"),
        ("outcomes_inspected", prospective, False, False, "decision_not_prospectively_locked"),
        ("synthetic_complete_still_not_ready", complete, True, True, None),
    )
    rows = []
    for name, frame, expected_valid, expected_complete, required in cases:
        result = validate_decision_selections(frame, catalog)
        rows.append({"Case": name, "ExpectedValid": expected_valid, "ActualValid": result["valid"], "ExpectedSelectionsComplete": expected_complete, "ActualSelectionsComplete": result["SelectionsComplete"], "OptionsSelected": result["OptionsSelected"], "RecruitmentReady": result["RecruitmentReady"], "RequiredFailureCode": required or "", "ActualFailureCodes": "|".join(result["failure_codes"]), "Passed": bool(result["valid"] == expected_valid and result["SelectionsComplete"] == expected_complete and not result["RecruitmentReady"] and (required is None or required in result["failure_codes"]))})
    return pd.DataFrame(rows)


def _tamper(bundle: bytes) -> bytes:
    with zipfile.ZipFile(io.BytesIO(bundle), "r") as source: payloads = {name: source.read(name) for name in source.namelist()}
    payloads["protocol_decision_workbench.html"] = payloads["protocol_decision_workbench.html"].replace(b"BLOCKED", b"READY  ", 1)
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(payloads):
            info = zipfile.ZipInfo(name, (1980,1,1,0,0,0)); info.compress_type = zipfile.ZIP_DEFLATED; info.external_attr = 0o100644 << 16; info.create_system = 3; archive.writestr(info, payloads[name])
    return buffer.getvalue()


def bundle_audit(first: dict[str, object], second: dict[str, object]) -> pd.DataFrame:
    original = validate_private_decision_workbench_bundle(first["bundle_bytes"]); tampered = validate_private_decision_workbench_bundle(_tamper(first["bundle_bytes"]))
    return pd.DataFrame([{"BundleSHA256": sha256_bytes(first["bundle_bytes"]), "RepeatedSHA256": sha256_bytes(second["bundle_bytes"]), "ByteIdentical": first["bundle_bytes"] == second["bundle_bytes"], "RegisteredMembers": len(WORKBENCH_BUNDLE_MEMBERS), "ManifestRows": len(first["manifest"]), "OriginalValid": original["valid"], "TamperedValid": tampered["valid"], "TamperedFailureCodes": "|".join(tampered["failure_codes"]), "Passed": bool(first["bundle_bytes"] == second["bundle_bytes"] and original["valid"] and not tampered["valid"] and "manifest_hash_mismatch" in tampered["failure_codes"] and "invalid_workbench_html" in tampered["failure_codes"] and first["manifest"]["HumanRows"].eq(0).all() and not first["manifest"]["ContainsSelection"].any())}])


def render_figures(dependencies: pd.DataFrame, summary: pd.DataFrame, output: Path) -> None:
    ordered = list(dependencies.sort_values("TopologicalOrder")["DecisionId"]); index = {value: i for i, value in enumerate(ordered)}
    matrix = np.zeros((13,13), dtype=int)
    for child, parents in PREREQUISITES.items():
        for parent in parents: matrix[index[parent], index[child]] = 1
    fig, ax = plt.subplots(figsize=(9,7)); ax.imshow(matrix, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(13), ordered, rotation=55, ha="right"); ax.set_yticks(range(13), ordered)
    ax.set(xlabel="Dependent decision", ylabel="Prerequisite decision", title="Registered acyclic decision dependencies (1 = required first)")
    for i in range(13):
        for j in range(13):
            if matrix[i,j]: ax.text(j,i,"1",ha="center",va="center",color="white",fontweight="bold")
    fig.tight_layout(); fig.savefig(output / "decision_dependency_matrix.png", dpi=180); plt.close(fig)

    colors = {"blocked":"#b91c1c", "catalog_ready":"#0f766e", "guarded":"#2563eb", "unverified":"#a16207", "withheld":"#7c3aed"}
    frame = summary.sort_values("CardOrder", ascending=False); y=np.arange(len(frame))
    fig, ax=plt.subplots(figsize=(10,5.5)); ax.barh(y,np.ones(len(frame)),color=[colors.get(v,"#6b7280") for v in frame["Status"]]); ax.set_yticks(y,frame["CardId"])
    for position,(_,row) in enumerate(frame.iterrows()): ax.text(.02,position,f"{row['Status']}: {row['FirstRead']}",va="center",color="white",fontsize=8.5)
    ax.set(xlim=(0,1),xticks=[],title="Private decision workbench: options visible, all selections blocked"); [ax.spines[s].set_visible(False) for s in ("top","right","bottom")]
    fig.tight_layout(); fig.savefig(output / "decision_workbench_status.png",dpi=180); plt.close(fig)


def run_tests(output: Path) -> dict[str, object]:
    command=[sys.executable,"-m","pytest","-q","tests/test_cmle_one_click_confirmatory_workbench.py","tests/test_cmle_one_click_confirmatory_protocol.py","tests/test_cmle_one_click_confirmatory_dependence.py","tests/test_cmle_one_click_confirmatory_planning.py","tests/test_cmle_one_click_confirmatory_gate.py","tests/test_cmle_one_click_cognitive_interview.py","tests/test_cmle_one_click_comprehension.py","tests/test_cmle_one_click.py","tests/test_cmle_one_click_archive.py","tests/test_decision_stability.py","tests/test_threshold_decision_integration.py"]
    completed=subprocess.run(command,cwd=ROOT,text=True,capture_output=True,check=False); (output/"selected_tests_stdout.txt").write_text(completed.stdout,encoding="utf-8"); (output/"selected_tests_stderr.txt").write_text(completed.stderr,encoding="utf-8")
    return {"passed":completed.returncode==0,"returncode":completed.returncode,"command":command}


def write_review(output: Path, results: dict[str, object]) -> None:
    text=f"""# Private confirmatory decision workbench critical review

## Decision

The workbench software contract **{'passed' if results['contract_passed'] else 'failed'}**. All 39 options are visible, but zero are selected and recruitment remains blocked.

## What the workbench adds

Each of 13 unresolved decisions has three non-ranked options with strength, primary risk, evidence requirement, and qualitative impacts. The 32 registered prerequisite edges are acyclic. Sample-size decision SCI-02 has seven scientific prerequisites; recruitment, ethics, and consent occur later. Broad stage labels do not permit dependency skipping.

## What it deliberately does not do

The repository assigns no recommendation or default. A blank worksheet is valid but incomplete. Unknown/cross-decision choices, premature child decisions, missing numeric assumptions, repository-generated external evidence, and decisions after outcome inspection fail closed. Even a fully populated synthetic worksheet can be `SelectionsComplete=true` while `SubstantiveEvidenceVerified=false` and `RecruitmentReady=false`.

## Private offline artifact

The self-contained bilingual HTML has no form, scripts, external network resources, or writable controls. Its seven-member deterministic ZIP reproduced byte-for-byte and rejected a synthetic `BLOCKED` to `READY` mutation. Bundle SHA-256: `{results['bundle_sha256']}`.

The catalog is not exhaustive and qualitative impact labels are not effect-size or cost estimates. External scientific and institutional owners must resolve and substantively verify decisions prospectively. Human participants, selected options, numeric n, recruitment authority, confirmatory results, and public UI remain zero/false.
"""
    (output/"CONFIRMATORY_DECISION_WORKBENCH_CRITICAL_REVIEW.md").write_text(text,encoding="utf-8")


def main() -> None:
    parser=argparse.ArgumentParser(); parser.add_argument("--output",type=Path,default=OUTPUT); args=parser.parse_args()
    plan,amendment=validate_registration()
    if args.output.exists(): raise FileExistsError(f"Evidence output exists: {args.output}")
    args.output.mkdir(parents=True)
    first=build_private_decision_workbench_bundle(); second=build_private_decision_workbench_bundle()
    catalog=first["catalog"]; dependencies=first["dependencies"]; selections=first["selections"]; summary=first["summary"]
    human=decision_workbench_human_gate_status(first["selection_status"])
    reg_audit=registration_audit(plan); cat_audit=catalog_audit(catalog,plan); dep_result=validate_dependency_table(dependencies)
    dep_audit=pd.DataFrame([{"NodeCount":len(dependencies),"EdgeCount":sum(len(v) for v in PREREQUISITES.values()),"Valid":dep_result["valid"],"FailureCodes":"|".join(dep_result["failure_codes"]),"Passed":dep_result["valid"]}])
    select_audit=selection_attack_audit(catalog,dependencies); html_result=validate_workbench_html(first["html"],catalog,first["content_sha256"])
    html_audit=pd.DataFrame([{"HTMLByteCount":len(first["html"].encode()),"OptionIdsPresent":sum(option in first["html"] for option in catalog["OptionId"]),"Valid":html_result["valid"],"FailureCodes":"|".join(html_result["failure_codes"]),"Passed":html_result["valid"] and sum(option in first["html"] for option in catalog["OptionId"])==39}])
    bun_audit=bundle_audit(first,second)
    changed=catalog.copy(); changed.loc[0,"PrimaryRisk"] += " synthetic mutation"; changed_sha=compute_workbench_content_identity(changed,dependencies,selections)
    identity_audit=pd.DataFrame([{"OriginalSHA256":first["content_sha256"],"MutatedSHA256":changed_sha,"Passed":changed_sha!=first["content_sha256"]}])
    tables={"decision_option_catalog.csv":catalog,"decision_dependency_table.csv":dependencies,"decision_selection_template.csv":selections,"decision_workbench_summary.csv":summary,"bundle_manifest.csv":first["manifest"],"human_gate_status.csv":human,"registration_identity_audit.csv":reg_audit,"option_catalog_contract_audit.csv":cat_audit,"dependency_contract_audit.csv":dep_audit,"selection_validator_synthetic_audit.csv":select_audit,"html_contract_audit.csv":html_audit,"private_bundle_contract_audit.csv":bun_audit,"workbench_content_identity_mutation_audit.csv":identity_audit}
    for name,frame in tables.items(): write_csv(frame,args.output/name)
    (args.output/"protocol_decision_workbench.html").write_text(first["html"],encoding="utf-8"); (args.output/"private_confirmatory_decision_workbench.zip").write_bytes(first["bundle_bytes"])
    render_figures(dependencies,summary,args.output); tests=run_tests(args.output)
    row=human.iloc[0]
    gates={"identity_passed":bool(reg_audit["Passed"].all() and amendment["parent_plan_sha256"]==sha256_file(PLAN)),"catalog_passed":bool(len(catalog)==39 and cat_audit["Passed"].all()),"dependencies_passed":bool(dep_audit["Passed"].all()),"selection_fail_closed_passed":bool(select_audit["Passed"].all()),"html_passed":bool(html_audit["Passed"].all()),"bundle_passed":bool(bun_audit["Passed"].all()),"content_identity_passed":bool(identity_audit["Passed"].all()),"human_gate_passed":bool(int(row["HumanParticipants"])==0 and int(row["OptionsSelected"])==0 and not bool(row["SelectionsComplete"]) and not bool(row["SubstantiveEvidenceVerified"]) and not bool(row["RecruitmentReady"]) and not bool(row["MinimumValidPerCellRegistered"]) and not bool(row["PlannedRecruitmentNSelected"]) and not bool(row["ConfirmatoryResultAvailable"]) and not bool(row["PublicSurfaceEnabled"])),"tests_passed":bool(tests["passed"])}
    contract=all(gates.values()); results={**gates,"contract_passed":contract,"decision_count":13,"option_count":39,"dependency_edge_count":sum(len(v) for v in PREREQUISITES.values()),"sample_size_prerequisite_count":len(PREREQUISITES["SCI-02"]),"selection_attack_cases":len(select_audit),"options_selected":0,"selections_complete":False,"substantive_evidence_verified":False,"recruitment_ready":False,"bundle_members":len(WORKBENCH_BUNDLE_MEMBERS),"bundle_sha256":sha256_bytes(first["bundle_bytes"]),"workbench_content_sha256":first["content_sha256"],"human_participants":0,"minimum_valid_per_cell_registered":False,"planned_recruitment_n_selected":False,"confirmatory_result_available":False,"public_surface_enabled":False,"selected_tests":tests}
    write_review(args.output,results)
    files=sorted(path.relative_to(args.output).as_posix() for path in args.output.rglob("*") if path.is_file() and path.name!="decision.json")
    decision=json_safe({"study_id":plan["study_id"],"plan_sha256":sha256_file(PLAN),"amendment_sha256":sha256_file(AMENDMENT),"contract_passed":contract,"contract_interpretation":"private_nonrecommending_workbench_all_selections_blocked","implementation_sha256":{"mfrm_app/cmle_one_click_confirmatory_workbench.py":sha256_file(ROOT/"mfrm_app/cmle_one_click_confirmatory_workbench.py"),"tests/test_cmle_one_click_confirmatory_workbench.py":sha256_file(ROOT/"tests/test_cmle_one_click_confirmatory_workbench.py"),"validation/cmle_one_click_confirmatory_decision_workbench.py":sha256_file(Path(__file__))},"results":results,"output_sha256":{name:sha256_file(args.output/name) for name in files},"interpretation":{"options":"not recommendations or defaults","selection":"none","evidence":"not substantively verified","human_study":"not ready","sample_size":"not selected","human_data":"none","public_ui":"withheld"}})
    (args.output/"decision.json").write_text(json.dumps(decision,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8"); print(json.dumps(decision,indent=2,ensure_ascii=False,allow_nan=False))
    if not contract: raise SystemExit(1)


if __name__=="__main__": main()
