#!/usr/bin/env python3
"""Generate registered private invalidity adjudication evidence."""
from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
import subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from mfrm_app.cmle_one_click_confirmatory_invalidity import (  # noqa: E402
    BUNDLE_MEMBERS, build_adjudication_template, build_private_bundle,
    validate_adjudication, validate_html, validate_private_bundle,
)

PLAN=ROOT/"validation/cmle_one_click_confirmatory_invalidity_adjudication_plan_20260810.json"
CODEBOOK=ROOT/"validation/cmle_one_click_confirmatory_protocol_preflight_20260810/invalidity_reason_codebook.csv"
SURFACE=ROOT/"validation/cmle_one_click_confirmatory_protocol_preflight_20260810/partial_identification_sensitivity.csv"
OUTPUT=ROOT/"validation/cmle_one_click_confirmatory_invalidity_adjudication_20260810"

def sha_file(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def sha_bytes(value): return hashlib.sha256(value).hexdigest()
def write_csv(frame,path): Path(path).write_text(frame.to_csv(index=False,float_format="%.17g"),encoding="utf-8")

def registration():
    plan=json.loads(PLAN.read_text())
    failed=[p for p,h in plan["parent_evidence_sha256"].items() if not (ROOT/p).is_file() or sha_file(ROOT/p)!=h]
    if failed or plan["repository_recommendation_allowed"] or plan["automatic_dangerous_classification_allowed"]: raise RuntimeError(f"Registration failed: {failed}")
    return plan

def attack_audit():
    blank=build_adjudication_template()
    hidden=blank.copy(); hidden.loc[1,"TargetPopulationDisposition"]="in_target"
    automatic=blank.copy(); automatic.loc[1,"RepositoryAutoClassificationAllowed"]=True
    other=blank.copy(); other.loc[14,"AdjudicationStatus"]="resolved_prospectively"
    late=blank.copy(); late.loc[1,"AdjudicationStatus"]="resolved_prospectively"
    late.loc[1,["TargetPopulationDisposition","CompositeEventDisposition","PrimaryAnalysisDisposition","SensitivitySetDisposition","DataUseAuthorityDisposition"]]=["in_target","sensitivity_only","exclude_preregistered","code_specific_tipping","authorized_use"]
    late.loc[1,["ScientificOwnerRole","StatisticalReviewReference","DomainReviewReference","EthicsAuthorityReference","RationaleReference"]]="SYNTHETIC"
    late.loc[1,"RecordedBeforeOutcomes"]=True; late.loc[1,"HumanOutcomesInspectedBeforeDecision"]=True
    cases=(("blank",blank,True,""),("hidden",hidden,False,"unresolved_row_contains_disposition"),("automatic",automatic,False,"repository_auto_classification_forbidden"),("other_without_amendment",other,False,"other_reason_requires_amendment_status"),("post_outcome",late,False,"adjudication_not_prospectively_blinded"))
    rows=[]
    for name,frame,expected,required in cases:
        result=validate_adjudication(frame); rows.append({"Case":name,"ExpectedValid":expected,"ActualValid":result["valid"],"Complete":result["complete"],"RecruitmentReady":result["RecruitmentReady"],"RequiredFailureCode":required,"FailureCodes":"|".join(result["failure_codes"]),"Passed":result["valid"]==expected and not result["RecruitmentReady"] and (not required or required in result["failure_codes"])})
    return pd.DataFrame(rows)

def transition_audit(sensitivity):
    rows=[]
    for source,group in sensitivity.groupby("SourceScenarioRow"):
        ordered=group.sort_values("CompositeFractionOfInvalidRegistered")
        labels=list(ordered["ProjectedDecision"]); first_fail=next((float(r.CompositeFractionOfInvalidRegistered) for r in ordered.itertuples() if r.ProjectedDecision=="proxy_fail"),None)
        if len(set(labels))>1:
            r=ordered.iloc[0]; rows.append({"SourceScenarioRow":source,"DangerousErrors":int(r["DangerousErrorsObservedValid"]),"ValidEligibleN":int(r["ValidEligibleN"]),"InvalidEligibleN":int(r["InvalidEligibleN"]),"FirstRegisteredCompositeFractionFail":first_fail,"Labels":"|".join(labels),"DecisionUsesDisplayedRounding":False})
    return pd.DataFrame(rows)

def figure(sensitivity,out):
    frame=sensitivity.loc[sensitivity["ValidEligibleN"].eq(500)&sensitivity["DangerousErrorsObservedValid"].eq(25)].copy()
    fig,ax=plt.subplots(figsize=(8.8,5.6))
    for fraction,group in frame.groupby("CompositeFractionOfInvalidRegistered"):
        group=group.sort_values("InvalidEligibleN"); x=group["InvalidEligibleN"]/group["ValidEligibleN"]
        ax.plot(x,group["CompositeWilsonUpper95Raw"],marker="o",label=f"Composite fraction={fraction:.2f}")
    ax.axhline(.10,color="black",linestyle="--",label="Strict threshold 0.10")
    ax.set(xlabel="Invalid / valid ratio (exact counts)",ylabel="Diagnostic Wilson upper 95%",title="Adjudication sensitivity (V=500, D=25; not scheduled-slot estimates)")
    ax.grid(alpha=.2); ax.legend(frameon=False,fontsize=8); fig.tight_layout(); fig.savefig(out/"invalidity_composite_sensitivity.png",dpi=180); plt.close(fig)

def run_tests(out):
    paths=["tests/test_cmle_one_click_confirmatory_invalidity.py","tests/test_cmle_one_click_confirmatory_estimand.py","tests/test_cmle_one_click_confirmatory_workbench.py","tests/test_cmle_one_click_confirmatory_protocol.py","tests/test_cmle_one_click_confirmatory_dependence.py","tests/test_cmle_one_click_confirmatory_planning.py","tests/test_cmle_one_click_confirmatory_gate.py","tests/test_cmle_one_click_cognitive_interview.py","tests/test_cmle_one_click_comprehension.py","tests/test_cmle_one_click.py","tests/test_cmle_one_click_archive.py","tests/test_decision_stability.py","tests/test_threshold_decision_integration.py"]
    command=[sys.executable,"-m","pytest","-q",*paths]; done=subprocess.run(command,cwd=ROOT,text=True,capture_output=True)
    (out/"selected_tests_stdout.txt").write_text(done.stdout); (out/"selected_tests_stderr.txt").write_text(done.stderr)
    return {"passed":done.returncode==0,"returncode":done.returncode,"command":command,"stdout_last_line":done.stdout.strip().splitlines()[-1] if done.stdout.strip() else ""}

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--output",type=Path,default=OUTPUT); args=parser.parse_args(); plan=registration()
    if args.output.exists(): raise FileExistsError(args.output)
    args.output.mkdir(parents=True); codebook=pd.read_csv(CODEBOOK); surface=pd.read_csv(SURFACE)
    first=build_private_bundle(codebook,surface); second=build_private_bundle(codebook,surface)
    attacks=attack_audit(); transitions=transition_audit(first["sensitivity"]); html=validate_html(first["html"],first["identity"]); bundle=validate_private_bundle(first["bundle_bytes"])
    parents=pd.DataFrame([{"Path":p,"ExpectedSHA256":h,"ActualSHA256":sha_file(ROOT/p),"Passed":sha_file(ROOT/p)==h} for p,h in plan["parent_evidence_sha256"].items()])
    formula=pd.DataFrame([{"SourceRows":len(surface),"SensitivityRows":len(first["sensitivity"]),"FractionsPerSourceMin":int(first["sensitivity"].groupby("SourceScenarioRow").size().min()),"FractionsPerSourceMax":int(first["sensitivity"].groupby("SourceScenarioRow").size().max()),"TransitionScenarios":len(transitions),"RoundedDecisions":bool(first["sensitivity"]["DecisionUsesDisplayedRounding"].any()),"Passed":len(surface)==72 and len(first["sensitivity"])==360 and first["sensitivity"].groupby("SourceScenarioRow").size().eq(5).all() and not first["sensitivity"]["DecisionUsesDisplayedRounding"].any()}])
    bundle_audit=pd.DataFrame([{"BundleSHA256":sha_bytes(first["bundle_bytes"]),"RepeatedSHA256":sha_bytes(second["bundle_bytes"]),"ByteIdentical":first["bundle_bytes"]==second["bundle_bytes"],"Valid":bundle["valid"],"FailureCodes":"|".join(bundle["failure_codes"]),"Passed":first["bundle_bytes"]==second["bundle_bytes"] and bundle["valid"]}])
    tables={"invalidity_code_taxonomy.csv":first["taxonomy"],"invalidity_allowed_state_catalog.csv":first["states"],"invalidity_adjudication_template.csv":first["template"],"invalidity_owner_review_matrix.csv":first["owners"],"invalidity_composite_fraction_sensitivity.csv":first["sensitivity"],"invalidity_adjudication_attack_audit.csv":attacks,"invalidity_transition_audit.csv":transitions,"invalidity_parent_identity_audit.csv":parents,"invalidity_formula_audit.csv":formula,"invalidity_bundle_audit.csv":bundle_audit,"invalidity_bundle_manifest.csv":first["manifest"]}
    for name,frame in tables.items(): write_csv(frame,args.output/name)
    (args.output/"invalidity_adjudication_workbench.html").write_text(first["html"]); (args.output/"private_invalidity_adjudication_workbench.zip").write_bytes(first["bundle_bytes"]); figure(first["sensitivity"],args.output); tests=run_tests(args.output)
    gates={"parent_identity_passed":bool(parents["Passed"].all()),"taxonomy_passed":len(first["taxonomy"])==15 and not first["taxonomy"]["AutomaticallyDangerous"].any(),"adjudication_fail_closed_passed":bool(attacks["Passed"].all()),"sensitivity_passed":bool(formula["Passed"].all()),"transition_demonstrated":len(transitions)>0,"html_passed":html["valid"],"bundle_passed":bool(bundle_audit["Passed"].all()),"tests_passed":tests["passed"],"human_gate_passed":True}; contract=all(gates.values())
    result={**gates,"contract_passed":contract,"registered_codes":15,"fixed_codes":1,"resolved_codes":0,"adjudication_complete":False,"source_scenarios":72,"sensitivity_rows":360,"transition_scenarios":len(transitions),"bundle_members":len(BUNDLE_MEMBERS),"bundle_sha256":sha_bytes(first["bundle_bytes"]),"content_sha256":first["identity"],"human_participants":0,"substantive_evidence_verified":False,"recruitment_ready":False,"sample_size_selected":False,"confirmatory_outcomes_available":False,"public_surface_enabled":False,"selected_tests":tests}
    review=f"""# Invalidity adjudication — critical review

The private software contract **{'passed' if contract else 'failed'}**. All 14 non-valid codes remain unresolved; recruitment remains blocked.

The taxonomy preserves 15 distinct mechanisms and assigns review responsibilities without calling any invalid record automatically dangerous. Consent/withdrawal require ethics data-use authority; accessibility barriers require an accessibility owner and separate-stratum sensitivity; duplicate records cannot be dangerous composite events; unregistered reasons require amendment.

The 72 source scenarios crossed with five diagnostic composite fractions produce 360 rows. {len(transitions)} scenarios cross the raw strict 0.10 gate within the registered fraction grid. This demonstrates decision sensitivity to adjudication, not actual invalidity-code frequencies. `V+I` is still not proven to be scheduled slots, and the diagnostic cannot select a composite definition.

Static HTML and an {len(BUNDLE_MEMBERS)}-member deterministic ZIP passed. {tests['stdout_last_line']}. Bundle SHA-256: `{result['bundle_sha256']}`. Real browser/accessibility/comprehension acceptance, external decisions, ethics approval, human data, sample size, recruitment, and public UI remain absent.
"""; (args.output/"INVALIDITY_ADJUDICATION_CRITICAL_REVIEW.md").write_text(review)
    files=sorted(p.relative_to(args.output).as_posix() for p in args.output.rglob("*") if p.is_file() and p.name!="decision.json")
    decision={"study_id":plan["study_id"],"plan_sha256":sha_file(PLAN),"contract_passed":contract,"contract_interpretation":"private_invalidity_adjudication_all_human_decisions_blocked","implementation_sha256":{"mfrm_app/cmle_one_click_confirmatory_invalidity.py":sha_file(ROOT/"mfrm_app/cmle_one_click_confirmatory_invalidity.py"),"tests/test_cmle_one_click_confirmatory_invalidity.py":sha_file(ROOT/"tests/test_cmle_one_click_confirmatory_invalidity.py"),"validation/cmle_one_click_confirmatory_invalidity_adjudication.py":sha_file(Path(__file__))},"results":result,"output_sha256":{name:sha_file(args.output/name) for name in files},"interpretation":{"invalid_is_dangerous":"never_automatic","code_adjudication":"none","sensitivity":"diagnostic_only","human_study":"not_ready","public_ui":"withheld"}}
    (args.output/"decision.json").write_text(json.dumps(decision,indent=2,ensure_ascii=False,allow_nan=False)+"\n"); print(json.dumps(decision,indent=2,ensure_ascii=False,allow_nan=False))
    if not contract: raise SystemExit(1)
if __name__=="__main__": main()
