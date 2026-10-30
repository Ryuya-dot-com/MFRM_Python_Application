"""Check the built macOS .app bundle, verify icon/signature/restore, and package with START_HERE.txt."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import plistlib
import shutil
import subprocess
import tempfile


def check_build(build: Path, output: Path) -> dict:
    build, output = build.resolve(), output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    app_bundle = build / "dist" / "MFRM-Desktop.app"
    exe = app_bundle / "Contents" / "MacOS" / "MFRM-Desktop"
    info_plist_path = app_bundle / "Contents" / "Info.plist"
    if not exe.is_file():
        raise FileNotFoundError(exe)
    if not info_plist_path.is_file():
        raise FileNotFoundError(info_plist_path)

    plist = plistlib.loads(info_plist_path.read_bytes())
    icon_file = str(plist.get("CFBundleIconFile", ""))
    if not icon_file:
        raise ValueError("CFBundleIconFile is missing from MFRM-Desktop.app Info.plist")
    icon_candidates = [
        app_bundle / "Contents" / "Resources" / icon_file,
        app_bundle / "Contents" / "Resources" / (icon_file if icon_file.endswith(".icns") else f"{icon_file}.icns"),
    ]
    if not any(p.is_file() for p in icon_candidates):
        raise FileNotFoundError(f"Bundled .icns resource not found for CFBundleIconFile={icon_file}")

    manifest = json.loads((build / "desktop_build.json").read_text(encoding="utf-8"))
    root = Path(__file__).resolve().parents[1]
    for relative, expected in manifest["source_sha256"].items():
        if hashlib.sha256((root / relative).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Source changed since build: {relative}")

    smoke = output / "analysis"
    with tempfile.TemporaryDirectory(prefix="mfrm-macos-frozen-check-") as working:
        environment = dict(os.environ)
        environment.pop("PYTHONPATH", None)
        environment.pop("PYTHONHOME", None)

        def run(*args):
            with (output / "process.log").open("ab") as log:
                subprocess.run(
                    [str(exe), *map(str, args)],
                    cwd=working,
                    env=environment,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )

        run("--smoke-test", smoke)
        record = json.loads((smoke / "smoke.json").read_text(encoding="utf-8"))
        expected = {("RSM", "JMLE"), ("PCM", "JMLE"), ("RSM", "MML"), ("PCM", "MML"), ("G/D", "ANOVA")}
        if {(r["model"], r["method"]) for r in record["extended"]} != expected:
            raise ValueError("Frozen checks did not cover all shipped workflows")
        reopened = []
        for proj in [smoke / "作文.mfrm", *sorted(smoke.glob("*/analysis.mfrm"))]:
            before = hashlib.sha256(proj.read_bytes()).hexdigest()
            destination = output / "reopened" / proj.parent.name
            run("--check-project", proj, "--output-dir", destination)
            res = json.loads((destination / "reopen.json").read_text(encoding="utf-8"))
            if res["estimation_calls"] != 0:
                raise ValueError("Restore re-estimated a model")
            if (destination / "report.html").read_bytes() != (proj.parent / "report.html").read_bytes():
                raise ValueError(f"Restored report differs: {proj}")
            if hashlib.sha256(proj.read_bytes()).hexdigest() != before:
                raise ValueError(f"Original project was changed: {proj}")
            reopened.append(dict(workflow=proj.parent.name, **res))
        if len(reopened) != 6:
            raise ValueError("Expected the basic project and five extended projects to restore independently")

    archive = build / "MFRM-Desktop-macOS-arm64.zip"
    if archive.exists():
        archive.unlink()
    with tempfile.TemporaryDirectory(prefix="mfrm-macos-stage-") as stage_raw:
        stage_dir = Path(stage_raw) / "MFRM-Desktop-macOS-arm64"
        stage_dir.mkdir(parents=True)
        subprocess.run(["ditto", str(app_bundle), str(stage_dir / "MFRM-Desktop.app")], check=True)
        shutil.copyfile(root / "docs" / "desktop_macos_installation.md", stage_dir / "START_HERE.txt")
        subprocess.run(
            ["ditto", "-c", "-k", "--sequesterRsrc", "--keepParent", str(stage_dir), str(archive)],
            check=True,
        )

    with tempfile.TemporaryDirectory(prefix="mfrm-macos-verify-") as verify_raw:
        verify_dir = Path(verify_raw)
        subprocess.run(["ditto", "-x", "-k", str(archive), str(verify_dir)], check=True)
        extracted_app = verify_dir / "MFRM-Desktop-macOS-arm64" / "MFRM-Desktop.app"
        extracted_guide = verify_dir / "MFRM-Desktop-macOS-arm64" / "START_HERE.txt"
        if not extracted_guide.is_file():
            raise FileNotFoundError("START_HERE.txt missing from macOS archive")
        subprocess.run(["codesign", "--verify", "--deep", "--strict", str(extracted_app)], check=True)

    sha256 = hashlib.sha256(archive.read_bytes()).hexdigest()
    (build / "SHA256SUMS.txt").write_text(f"{sha256}  {archive.name}\n", encoding="utf-8")
    result = dict(
        workflows=record["extended"],
        reopened=reopened,
        archive=archive.name,
        bytes=archive.stat().st_size,
        sha256=sha256,
        icon_file=icon_file,
        codesign_adhoc_verified=True,
        apple_notarized=False,
    )
    (output / "checks.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(check_build(args.build_dir, args.output), indent=2))
