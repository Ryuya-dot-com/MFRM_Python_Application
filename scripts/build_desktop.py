"""Build a local prototype; does not sign, notarize, publish, or install it."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (sys.platform, platform.machine().lower()) not in {("darwin", "arm64"), ("win32", "amd64"), ("win32", "x86_64")}:
        parser.error("Initial targets are macOS arm64 and Windows x64; build on the target OS/CPU")
    root = Path(__file__).resolve().parents[1]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source_paths = [root / "desktop_app.py", root / "examples" / "desktop_writing.csv",
                    root / "LICENSE", root / "LICENSE_NOTICE.md", root / "locales" / "en.json",
                    Path(__file__).resolve(), root / "requirements-desktop-build.txt",
                    *sorted((root / "mfrm_app").glob("*.py"))]
    record = {
        "scope": "local prototype; developer signing/notarization and platform acceptance not completed",
        "platform": platform.platform(), "machine": platform.machine(), "python": sys.version,
        "source_sha256": {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths},
        "packages": {name: importlib.metadata.version(name) for name in
                     ("pyinstaller", "numpy", "pandas", "scipy", "matplotlib", "pyarrow")},
    }
    manifest = output / "desktop_build.json"
    manifest.write_text(json.dumps(record, indent=2), encoding="utf-8")
    command = [sys.executable, "-m", "PyInstaller", "--noconfirm", "--windowed", "--onedir",
               "--name", "MFRM-Desktop", "--distpath", str(output / "dist"),
               "--workpath", str(output / "work"), "--specpath", str(output),
               "--add-data", f"{manifest}:mfrm_app", "--exclude-module", "streamlit",
               "--add-data", f"{root / 'locales' / 'en.json'}:locales",
               "--add-data", f"{root / 'examples' / 'desktop_writing.csv'}:examples",
               "--add-data", f"{root / 'LICENSE'}:.",
               "--add-data", f"{root / 'LICENSE_NOTICE.md'}:.",
               "--hidden-import", "matplotlib.backends.backend_svg",
               "--hidden-import", "matplotlib.backends.backend_tkagg",
               "--exclude-module", "plotly", "--exclude-module", "pytest"]
    for name in ("numpy", "pandas", "scipy"):
        command.extend(["--copy-metadata", name])
    if sys.platform == "darwin":
        command.extend(["--target-arch", "arm64"])
    command.append(str(root / "desktop_app.py"))
    environment = dict(os.environ, PYINSTALLER_CONFIG_DIR=str(output / "cache"))
    environment.setdefault("MPLCONFIGDIR", str(output / "matplotlib-cache"))
    subprocess.run(command, cwd=root, env=environment, check=True)


if __name__ == "__main__":
    main()
