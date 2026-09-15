#!/usr/bin/env python3
"""Run the prepared PCM comparison from a host where ConQuest starts normally.

Uses only the retained synthetic ratings. Does not change ConQuest, signatures,
OS settings, the licence, or old output. Every launch keeps its console output.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
PREPARED = ROOT / "validation/generated/mml_pcm_external_refit_20260913"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    with path.open("x") as f:
        json.dump(value, f, indent=2)
        f.write("\n")


def launch(directory, command, timeout):
    started = time.monotonic()
    code, timed_out = None, False
    with (directory / "console.log").open("x") as log:
        try:
            code = subprocess.run(["/usr/bin/arch", "-x86_64", "/Applications/ConQuest/ConQuest"],
                input=command, text=True, stdout=log, stderr=subprocess.STDOUT,
                cwd=directory, timeout=timeout, check=False).returncode
        except subprocess.TimeoutExpired:
            timed_out = True
    console = (directory / "console.log").read_text(errors="replace")
    result = dict(exit_code=code, timed_out=timed_out, elapsed=time.monotonic()-started,
                  end_of_program="End of Program" in console)
    write(directory / "launch.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    spec = json.loads((PREPARED / "input.json").read_text())
    assert sha(PREPARED / "summary.json") == "ff1984b687504942e1ca01dc544cbadac9eb9beecbe08e96a429ef0257e00f23"
    assert not spec["scientific_inference_ready"] and not spec["qualification_eligible"]
    for name, expected in spec["input_sha256"].items():
        assert sha(PREPARED / name) == expected, name
    binary = Path("/Applications/ConQuest/ConQuest")
    assert sha(binary) == spec["source_sha256"][str(binary)], "Executable changed: review the new version first"
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(PREPARED / "wide.csv", out / "wide.csv")
    for config in spec["conquest"]:
        directory = out / config["id"]
        directory.mkdir()
        shutil.copyfile(PREPARED / config["id"] / "model.cqc", directory / "model.cqc")
    write(out / "input.json", dict(classification="OBSERVED_DEVELOPMENT_ONLY",
        scientific_inference_ready=False, qualification_eligible=False, plan=spec["conquest"],
        source_script_sha256=sha(Path(__file__).resolve()), binary_sha256=sha(binary),
        prepared_input_sha256=sha(PREPARED / "input.json"),
        input_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file()}))
    sentinel = out / "sentinel"
    sentinel.mkdir()
    print("Checking ConQuest startup without data...", flush=True)
    result = launch(sentinel, "quit;\n", 60)
    if result["exit_code"] != 0 or not result["end_of_program"]:
        print("Startup did not complete. See", sentinel / "console.log", flush=True)
        return 1
    results = {}
    required = ("parameters.csv", "amatrix.csv", "reg_coefficients.csv", "covariance.csv",
                "cases.csv", "history.csv", "review.txt")
    for config in spec["conquest"]:
        name = config["id"]
        directory = out / name
        print("Running", name, flush=True)
        result = launch(directory, (directory / "model.cqc").read_text(), 600)
        result["missing_outputs"] = [n for n in required if not (directory/n).is_file() or not (directory/n).stat().st_size]
        result["execution_complete"] = result["exit_code"] == 0 and result["end_of_program"] and not result["missing_outputs"]
        results[name] = result
        print(name, "execution complete:", result["execution_complete"], flush=True)
        if not result["execution_complete"]:
            break  # Preserve the first failure before attempting dependent comparisons.
    write(out / "execution.json", dict(scientific_inference_ready=False, qualification_eligible=False,
        runs=results, remaining_unattempted=[c["id"] for c in spec["conquest"] if c["id"] not in results],
        artifact_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file()}))
    complete = len(results) == len(spec["conquest"]) and all(r["execution_complete"] for r in results.values())
    print("Saved:", out, "Execution complete:", complete, "Numerical review pending.", flush=True)
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
