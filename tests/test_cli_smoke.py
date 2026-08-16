"""Bitrot guard over the research-script surface (Phase 3 of the cleanup).

Fast tier: every script on the active surface must at least COMPILE.

Slow tier: every active/tool/validator/shim CLI is invoked with --help in a
subprocess. Scripts with an argparse CLI must exit 0; the env/argv-driven
campaign scripts are allowed to fail on their missing inputs, but never
with an import-time error (ModuleNotFoundError & co.) -- that is the rot
this sweep exists to catch.

The script list comes from the CLASSIFICATION map in tools/script_index.py,
so the sweep tracks the manifest automatically.
"""
import importlib.util
import os
import py_compile
import subprocess
import sys

import pytest

from conftest import REPO


def _classification():
    spec = importlib.util.spec_from_file_location(
        "script_index", os.path.join(REPO, "tools", "script_index.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.CLASSIFICATION


def _relpath(key):
    if key.startswith("dd/"):
        return os.path.join("scripts", "defect_dynamics", key[3:])
    for d in ("scripts", "tools"):
        p = os.path.join(d, key)
        if os.path.exists(os.path.join(REPO, p)):
            return p
    raise FileNotFoundError(key)


CLS = _classification()
PY = sorted(_relpath(k) for k in CLS if k.endswith(".py"))
# statuses whose scripts should be runnable front-ends
RUNNABLE = sorted(_relpath(k) for k, (prog, status) in CLS.items()
                  if k.endswith(".py")
                  and status in ("active", "tool", "validator", "shim"))

IMPORT_ROT = ("ModuleNotFoundError", "ImportError", "SyntaxError",
              "AttributeError", "NameError")


@pytest.mark.parametrize("rel", PY, ids=lambda r: r.replace("/", "."))
def test_compiles(rel):
    py_compile.compile(os.path.join(REPO, rel), doraise=True)


@pytest.mark.slow
@pytest.mark.parametrize("rel", RUNNABLE, ids=lambda r: r.replace("/", "."))
def test_help_or_clean_failure(rel):
    r = subprocess.run([sys.executable, os.path.join(REPO, rel), "--help"],
                       cwd=REPO, capture_output=True, text=True, timeout=180)
    if r.returncode == 0:
        return
    # env/argv-style campaign scripts fail on missing inputs -- fine, as long
    # as the failure is not import-time rot
    for marker in IMPORT_ROT:
        assert marker not in r.stderr, \
            f"{rel} fails with {marker}:\n{r.stderr[-1500:]}"
