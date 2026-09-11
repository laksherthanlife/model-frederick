"""Test-only subprocess probes for materialized frozen source trees.

Use the production frozen_runtime executors for historical evidence receipts.
This launcher instead permits white-box patches and deliberately tampered temporary
copies so the original release gates can be tested without altering parent imports.
It never substitutes module identities or treats a copied tree as execution.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import textwrap


ROOT = Path(__file__).resolve().parents[1]
_BOOTSTRAP = """
import json
import os
from pathlib import Path
import sys

root = Path(sys.argv[1]).resolve(strict=True)
original = Path(sys.argv.pop(2)).resolve(strict=True)
source = root / "src"
assert root != original
assert sys.flags.isolated and sys.dont_write_bytecode
assert os.environ["PYTHONPATH"] == str(source)
sys.path.insert(0, str(source))

def verify_test_modules():
    for name, module in tuple(sys.modules.items()):
        if name == "ystwin" or name.startswith("ystwin."):
            path = Path(module.__file__)
            assert path.resolve() == path and path.is_relative_to(source), name

verify_test_modules()

def forbid_original_checkout(event, arguments):
    if event == "open" and isinstance(arguments[0], (str, bytes)):
        candidate = Path(os.fsdecode(arguments[0])).resolve()
        if candidate.is_relative_to(original) and not candidate.is_relative_to(root):
            raise PermissionError("The current checkout is unavailable to this frozen-source probe")

sys.addaudithook(forbid_original_checkout)
"""


def run_frozen_python(root, program, *arguments, payload=None, timeout=120):
    """Run an explicit probe against real imports from a temporary source tree."""
    root = Path(root).resolve(strict=True)
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith(("PYTHON", "YSTWIN_", "GIT_"))}
    environment.update(PYTHONPATH=str(root / "src"), PYTHONDONTWRITEBYTECODE="1",
                       CI="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    # Probe locals (for example an "original" function saved before patching) must
    # not overwrite the bootstrap's import verifier or filesystem guard state.
    script = (_BOOTSTRAP
              + f"\nexec(compile({textwrap.dedent(program)!r}, '<frozen-source-probe>', 'exec'), "
                "{'__name__': '__main__', 'root': root, 'sys': sys, 'json': json, 'Path': Path})\n"
                "verify_test_modules()\n")
    return subprocess.run(
        [sys.executable, "-I", "-B", "-c", script,
         str(root), str(ROOT), *map(os.fspath, arguments)],
        input=json.dumps(payload, allow_nan=False), cwd=root, env=environment,
        capture_output=True, text=True, timeout=timeout,
    )
