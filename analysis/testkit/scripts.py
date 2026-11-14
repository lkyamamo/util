"""Run or import analysis scripts by path.

Many analysis folders are not importable package names (20260608_GrNrBaSqw,
lammps_GrBa_ovito), and most scripts are meant to be run, so tests either run
them as a subprocess or load them from their file.
"""

import importlib.util
import subprocess
import sys
from pathlib import Path


def run_script(path, *args, cwd=None):
    """Run a script with this interpreter; raise with its output if it fails."""
    proc = subprocess.run(
        [sys.executable, str(path), *map(str, args)],
        cwd=cwd, capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise AssertionError(
            f"{Path(path).name} exited {proc.returncode}\n"
            f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
        )
    return proc


def load_script(path, name=None):
    """Import a script from its file. Its folder goes on sys.path so sibling imports work."""
    path = Path(path).resolve()
    folder = str(path.parent)
    if folder not in sys.path:
        sys.path.insert(0, folder)
    spec = importlib.util.spec_from_file_location(name or path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
