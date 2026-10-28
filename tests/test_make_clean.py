"""Housekeeping must preserve retained research evidence and environments."""

from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import unittest


class MakeCleanTests(unittest.TestCase):
    @unittest.skipUnless(os.name == 'posix' and shutil.which('make'), 'POSIX make required')
    def test_clean_preserves_evidence_and_does_not_follow_symlinks(self):
        makefile = Path(__file__).resolve().parents[1] / "Makefile"
        caches = (
            ".pytest_cache/v/cache/nodeids", "__pycache__/app.pyc",
            "tests/__pycache__/test_app.pyc", "validation/__pycache__/runner.pyc",
            "mfrm_app/__pycache__/module.pyc",
            "mfrm_app/simulation/__pycache__/module.pyc",
        )
        protected = (
            "validation/generated/study/ledger.sqlite",
            "validation/generated/study/__pycache__/retained.pyc",
            "validation/source_snapshots/source.zip",
            "validation/older_study/__pycache__/retained.pyc",
            "backups/retained/__pycache__/module.pyc",
            ".venv/lib/__pycache__/module.pyc",
            ".venv-win/Lib/__pycache__/module.pyc",
            ".git/objects/retained", "tests/fixtures/__pycache__/fixture.pyc",
            "mfrm_app/module.py", "outside/__pycache__/keep.pyc",
        )
        with tempfile.TemporaryDirectory(prefix="mfrm-clean-test-") as tmp:
            root = Path(tmp)
            for name in caches + protected:
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"preserve this content\n")
            (root / "mfrm_app/linked").symlink_to(root / "outside", target_is_directory=True)
            subprocess.run(
                ["make", "-f", str(makefile), "clean"], cwd=root,
                check=True, capture_output=True, text=True,
            )
            for name in caches:
                self.assertFalse((root / name).exists(), name)
            for name in protected:
                self.assertEqual((root / name).read_bytes(), b"preserve this content\n", name)


if __name__ == "__main__":
    unittest.main()
