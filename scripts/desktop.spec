# -*- mode: python ; coding: utf-8 -*-
import json
import os
from pathlib import Path
import sys

from PyInstaller.utils.hooks import copy_metadata

root = Path(SPECPATH).parent
manifest = Path(os.environ["MFRM_DESKTOP_BUILD_MANIFEST"])
record = json.loads(manifest.read_text(encoding="utf-8"))
datas = [(str(manifest), "mfrm_app"),
         (str(root / "locales" / "en.json"), "locales"),
         (str(root / "examples" / "desktop_writing.csv"), "examples")]
datas.append((str(manifest.parent / 'third_party_notices'), 'third_party_notices'))
for name in ("LICENSE", "LICENSE_NOTICE.md", "CITATION.cff"):
    datas.append((str(root / name), "."))
for name in ("numpy", "pandas", "scipy"):
    datas += copy_metadata(name)

a = Analysis(
    [str(root / "desktop_app.py")], pathex=[str(root)], binaries=[], datas=datas,
    hiddenimports=["matplotlib.backends.backend_svg", "matplotlib.backends.backend_tkagg"],
    hookspath=[], hooksconfig={}, runtime_hooks=[],
    excludes=["streamlit", "plotly", "pytest"], noarchive=False, optimize=0,
    # Runtime qualification hashes the actual sources; load those same files.
    module_collection_mode={"mfrm_app": "py"},
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [], exclude_binaries=True, name="MFRM-Desktop",
    debug=False, bootloader_ignore_signals=False, strip=False, upx=True,
    console=False, disable_windowed_traceback=False, argv_emulation=False,
    target_arch="arm64" if sys.platform == "darwin" else None,
    codesign_identity=None, entitlements_file=None,
)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=True, upx_exclude=[], name="MFRM-Desktop")
if sys.platform == "darwin":
    # Set metadata before PyInstaller signs the bundle.
    app = BUNDLE(
        coll, name="MFRM-Desktop.app", icon=None,
        bundle_identifier="io.github.ryuya-dot-com.mfrm.desktop",
        version=record["bundle_version"],
        info_plist={"CFBundleVersion": record["bundle_version"],
                    "MFRMApplicationVersion": record["app_version"],
                    "CFBundleDocumentTypes": [{
                        "CFBundleTypeName": "MFRM project",
                        "CFBundleTypeRole": "Editor",
                        "LSHandlerRank": "Owner",
                        "LSItemContentTypes": ["io.github.ryuya-dot-com.mfrm.project"],
                    }],
                    "UTExportedTypeDeclarations": [{
                        "UTTypeIdentifier": "io.github.ryuya-dot-com.mfrm.project",
                        "UTTypeDescription": "MFRM project",
                        "UTTypeConformsTo": ["public.data"],
                        "UTTypeTagSpecification": {"public.filename-extension": ["mfrm"]},
                    }]},
    )
