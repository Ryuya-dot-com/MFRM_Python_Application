"""Execute the declared additive reference amendment; original driver is unchanged."""
import argparse
import json
from pathlib import Path
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    sys.path.insert(0, str(out / "frozen"))
    from validation import medium_writing_check as study
    from mml_rasch_block_reference import BlockRaschIntegral

    amendment_path = out / "reference_implementation_amendment.json"
    amendment = json.loads(amendment_path.read_text())
    amendment_hash = study.digest(amendment_path)
    plan = json.loads((out / "protocol.json").read_text())
    numerical_manifest = json.loads((out / "numerical_source_manifest.json").read_text())
    driver_manifest = json.loads((out / "driver_source_manifest.json").read_text())

    def verify():
        if study.digest(out / "protocol.json") != study.PROTOCOL_SHA256:
            raise ValueError("Prespecified study protocol changed")
        if study.digest(out / "numerical_source_manifest.json") != plan["source_manifest_sha256"]:
            raise ValueError("Original numerical source manifest changed")
        if Path(study.app.__file__).resolve() != out / "frozen" / "streamlit_app.py":
            raise ValueError("Reference runner must import the frozen application")
        if Path(__file__).resolve() != out / "reference_implementation" / "medium_writing_references.py":
            raise ValueError("Execute the retained additive reference runner")
        if study.digest(out / "frozen/validation/medium_writing_check.py") != driver_manifest["sha256"]:
            raise ValueError("Original study driver changed")
        if study.digest(amendment_path) != amendment_hash:
            raise ValueError("Reference implementation amendment changed during execution")
        for relative, expected in {**{str(Path('frozen') / k): v for k, v in numerical_manifest.items()},
                                   **amendment["source_sha256"]}.items():
            if study.digest(out / relative) != expected:
                raise ValueError("Retained source changed: " + relative)

    verify()
    provenance_path = out / "reference_execution_implementation.json"
    provenance = dict(amendment_sha256=amendment_hash, implementation="independent_rasch_person_blocks10",
                      scientific_inference_ready=False)
    if provenance_path.exists():
        if json.loads(provenance_path.read_text()) != provenance:
            raise ValueError("Existing references use a different implementation")
    else:
        for n, design, delta in study.cases(plan):
            path = out / "references" / (study.case_id(n, design, delta) + ".json")
            if path.exists() or path.with_suffix(".started.json").exists():
                raise ValueError("A planned reference was already attempted before this amendment")
        study.save(provenance_path, provenance)

    # Change only the reference constructor for this process. Point selection,
    # settings, structural-block covariance, checks, caps and failure retention
    # all remain in the independently retained original driver.
    study.RaschIntegral = BlockRaschIntegral
    started = time.monotonic()
    for n, design, delta in study.cases(plan):
        remaining = plan["resource_limits"]["reference_stage_wall_seconds"] - (time.monotonic() - started)
        if remaining <= 0:
            break
        name = study.case_id(n, design, delta)
        folder = out / "cases" / name / "rep_000"
        def evaluate():
            result = study.continuous_case(folder, n, plan)
            result.update(provenance)
            return result
        result = study.attempt(out / "references" / (name + ".json"),
            min(remaining, plan["continuous_reference"]["per_case_timeout_seconds"]), evaluate)
        print("reference", name, result["status"], result.get("all_checks_pass"), flush=True)
    verify()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
