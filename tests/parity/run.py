"""The local parity run: R against the Python ssdiff release pinned in parity.toml.

Needs Python 3.11+ (standard library only), R with testthat, and network for
the first run. From anywhere:
    python tests/parity/run.py

1. Reads parity.toml.
2. Creates or reuses a virtual environment for the pinned ssdiff under
   ~/.cache/ssdiff-parity/ and installs that release into it; plskit and glmm
   come from its dependencies.
3. Runs generate.py in it: the committed golden files on the test pack, and
   the reference cases on the real pack into the git-ignored out/.
4. Fetches the pinned test pack by name and compares it with the committed
   test.wordpack.
5. Runs test-parity.R on the golden cases and on the reference cases; when
   both pass, writes the R plskit and fastglmm versions they ran with into the
   goldens' versions.json.
6. Compares inst/page_template.html with the installed ssdiff's.

A difference in step 4 or 6 prints the cp command that takes Python's file.
Exits with 1 when step 4, 5 or 6 fails.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TEST_PACK = Path("tests/testthat/fixtures/wordpack/test.wordpack")
TEMPLATE = Path("inst/page_template.html")
VERSIONS = ROOT / "tests" / "testthat" / "fixtures" / "golden" / "versions.json"


def same_bytes(python_file: str, copy: Path) -> bool:
    """Whether the R copy (a path from the ssdiff-r folder) equals Python's file; prints the cp command if not."""
    if Path(python_file).read_bytes() == (ROOT / copy).read_bytes():
        print(f"{copy}: same as {python_file}")
        return True
    print(f"{copy} differs from {python_file}. To take Python's file, from the ssdiff-r folder run:\n"
          f"  cp '{python_file}' '{copy}'")
    return False


def main() -> int:
    pins = tomllib.loads((HERE / "parity.toml").read_text(encoding="utf-8"))
    venv = Path.home() / ".cache" / "ssdiff-parity" / f"r-ssdiff-{pins['ssdiff']}"
    py = str(venv / "bin" / "python")
    if not os.path.exists(py):
        subprocess.run([sys.executable, "-m", "venv", str(venv)], check=True)
    subprocess.run([py, "-m", "pip", "install", "--quiet", f"ssdiff=={pins['ssdiff']}"], check=True)
    # -I: no PYTHONPATH, user site-packages or current folder, so only the venv's ssdiff is importable.
    subprocess.run([py, "-I", str(HERE / "generate.py")], check=True)
    ok = True

    resolve = f"from ssdiff.wordpack.fetch import resolve; print(resolve({pins['testpack']!r}))"
    fetch = subprocess.run([py, "-I", "-c", resolve], capture_output=True, text=True)
    if fetch.returncode:
        print(f"Cannot fetch the test pack {pins['testpack']}:\n{fetch.stderr}")
        ok = False
    else:
        ok &= same_bytes(fetch.stdout.strip(), TEST_PACK)

    golden_env = {k: v for k, v in os.environ.items() if k != "SSDIFF_PARITY_DIR"}
    passed = True
    test_parity = 'testthat::test_local(".", filter = "^parity$", stop_on_failure = TRUE)'
    for env in (golden_env, golden_env | {"SSDIFF_PARITY_DIR": str(HERE / "out")}):
        tests = subprocess.run(["Rscript", "-e", test_parity], cwd=ROOT, env=env)
        passed &= tests.returncode == 0
    ok &= passed
    if passed:
        r_versions = subprocess.run(
            ["Rscript", "-e", 'cat(format(packageVersion("plskit")), format(packageVersion("fastglmm")))'],
            capture_output=True, text=True, check=True).stdout.split()
        versions = json.loads(VERSIONS.read_text(encoding="utf-8"))
        versions |= dict(zip(("r_plskit", "r_fastglmm"), r_versions))
        VERSIONS.write_text(json.dumps(versions, indent=1) + "\n", encoding="utf-8")

    where = "import pathlib, ssdiff; print(pathlib.Path(ssdiff.__file__).parent / 'results' / 'page_template.html')"
    template = subprocess.run([py, "-I", "-c", where], capture_output=True, text=True, check=True).stdout.strip()
    ok &= same_bytes(template, TEMPLATE)
    print("Parity run: " + ("all steps passed" if ok else "a step failed (see above)"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
