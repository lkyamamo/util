"""SLURM job scripts and scheduler stdout (STREAM_OUTPUT)."""

from __future__ import annotations

import re
from datetime import datetime
from typing import Dict, List, Optional

SBATCH_RE = re.compile(r"^#SBATCH\s+(\S.*)$")
PBS_RE = re.compile(r"^#PBS\s+(\S.*)$")
MODULE_RE = re.compile(r"^\s*module\s+(?:load|add)\s+(.+)$")
EXPORT_EXE_RE = re.compile(r"^\s*(?:export\s+)?(\w+)=[\"']?(\S*(?:lmp|vasp|lammps)\S*?)[\"']?\s*$", re.I)
LAUNCH_RE = re.compile(r"^\s*(?:srun|mpirun|mpiexec|ibrun)\b.*$")
# `date` output, e.g. "Wed Aug 26 17:54:44 PDT 2026"
DATE_RE = re.compile(r"^[A-Z][a-z]{2} [A-Z][a-z]{2} +\d{1,2} \d\d:\d\d:\d\d [A-Z]{2,5} \d{4}$")

FAILURE_PATTERNS = [
    ("cancelled", re.compile(r"CANCELLED", re.I)),
    ("time_limit", re.compile(r"DUE TO TIME LIMIT")),
    ("oom", re.compile(r"oom[-_ ]kill|Out Of Memory", re.I)),
    ("segfault", re.compile(r"Segmentation fault|SIGSEGV")),
    ("mpi_abort", re.compile(r"MPI_ABORT|application called MPI_Abort")),
    ("error", re.compile(r"^ERROR\b|^\s*ERROR:", re.M)),
]


def _option(rest: str):
    rest = rest.split("#", 1)[0].strip()  # drop trailing comments
    if rest.startswith("--"):
        token = rest[2:]
        if "=" in token:
            key, _, val = token.partition("=")
        else:
            key, _, val = token.partition(" ")
        return key.strip(), val.strip()
    if rest.startswith("-"):
        key, _, val = rest[1:].partition(" ")
        return key.strip(), val.strip()
    return None


SHORT_TO_LONG = {"N": "nodes", "n": "ntasks", "p": "partition", "t": "time", "J": "job-name",
                 "A": "account", "o": "output", "e": "error", "c": "cpus-per-task", "C": "constraint",
                 "G": "gpus"}


def parse_script(text: str) -> Dict:
    directives: Dict[str, str] = {}
    modules: List[str] = []
    launch: List[str] = []
    executables: Dict[str, str] = {}
    for line in text.splitlines():
        s = line.strip()
        m = SBATCH_RE.match(s) or PBS_RE.match(s)
        if m:
            opt = _option(m.group(1))
            if opt:
                key, val = opt
                directives[SHORT_TO_LONG.get(key, key)] = val
            continue
        if s.startswith("#"):
            continue
        m = MODULE_RE.match(s)
        if m:
            modules.extend(m.group(1).split())
            continue
        m = EXPORT_EXE_RE.match(s)
        if m:
            executables[m.group(1)] = m.group(2)
        if LAUNCH_RE.match(s):
            launch.append(s)
    return {
        "directives": directives,
        "modules": modules,
        "launch": launch,
        "executables": executables,
        "scheduler": "slurm" if any(SBATCH_RE.match(l.strip()) for l in text.splitlines()) else
                     ("pbs" if "#PBS" in text else None),
    }


def _parse_date(line: str) -> Optional[str]:
    s = line.strip()
    if not DATE_RE.match(s):
        return None
    parts = s.split()
    # Drop the timezone token; strptime's %Z is unreliable across platforms.
    try:
        dt = datetime.strptime(" ".join(parts[:4] + parts[5:]), "%a %b %d %H:%M:%S %Y")
    except ValueError:
        return None
    return dt.isoformat()


def parse_stdout(head: str, tail: str) -> Dict:
    """Start/end timestamps (from `date` lines) and failure signatures."""
    start = next((d for d in map(_parse_date, head.splitlines()) if d), None)
    end = next((d for d in map(_parse_date, reversed(tail.splitlines())) if d), None)
    failures = [name for name, rx in FAILURE_PATTERNS if rx.search(tail)]
    errors = [l.strip() for l in tail.splitlines() if re.match(r"^\s*ERROR", l)][:10]
    return {"start": start, "end": end if end != start else None, "failures": failures, "errors": errors}
