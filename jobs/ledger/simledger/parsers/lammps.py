"""LAMMPS input scripts and log.lammps (head + tail only)."""

from __future__ import annotations

import re
from typing import Dict, List

from . import to_num

# Settings stored as single params (last occurrence wins).
SINGLE = ("units", "atom_style", "boundary", "dimension", "pair_style", "kspace_style",
          "bond_style", "angle_style", "dihedral_style", "timestep", "read_data",
          "read_restart", "processors", "neighbor")

INTEGRATORS = {"nve", "nvt", "npt", "nph", "langevin", "berendsen", "temp/berendsen",
               "temp/rescale", "temp/csvr", "nvt/sllod", "rigid/nvt", "rigid/nve"}

LOG_VERSION_RE = re.compile(r"^LAMMPS \((.+)\)")
ATOMS_RE = re.compile(r"^\s*(\d+) atoms\s*$")
LOOP_RE = re.compile(r"^Loop time of ([\d.eE+-]+) on (\d+) procs for (\d+) steps with (\d+) atoms")
PERF_RE = re.compile(r"^Performance:\s*(.*)$")
WALL_RE = re.compile(r"^Total wall time: (\d+):(\d\d):(\d\d)")


def _strip(line: str) -> str:
    return line.split("#", 1)[0].strip()


def parse_input(text: str) -> Dict:
    """Commands of interest, without evaluating variables.

    Values containing ${...} or v_ references are kept verbatim and flagged."""
    settings: Dict[str, str] = {}
    fixes: List[Dict] = []
    runs: List[int] = []
    minimize = 0
    pair_coeff: List[str] = []
    dumps: List[Dict] = []
    variables: Dict[str, str] = {}
    unresolved = set()
    joined = re.sub(r"&\s*\n", " ", text)  # continuation lines
    for raw in joined.splitlines():
        line = _strip(raw)
        if not line:
            continue
        cmd, _, rest = line.partition(" ")
        rest = rest.strip()
        for ref in re.findall(r"\$\{(\w+)\}|\$(\w)\b", rest):
            unresolved.add(ref[0] or ref[1])
        if cmd in SINGLE:
            settings[cmd] = rest
        elif cmd == "pair_coeff":
            pair_coeff.append(rest)
        elif cmd == "variable":
            parts = rest.split(None, 2)
            if len(parts) >= 3:
                variables[parts[0]] = f"{parts[1]} {parts[2]}"
        elif cmd == "fix":
            parts = rest.split()
            if len(parts) >= 3:
                fix = {"id": parts[0], "group": parts[1], "style": parts[2], "args": " ".join(parts[3:])}
                if parts[2] in INTEGRATORS:
                    args = parts[3:]
                    for kw in ("temp", "iso", "aniso", "x", "y", "z"):
                        if kw in args:
                            i = args.index(kw)
                            fix[kw] = " ".join(args[i + 1:i + 4])
                    if parts[2] in ("langevin", "temp/berendsen", "temp/rescale", "temp/csvr") and len(args) >= 3:
                        fix["temp"] = " ".join(args[:3])
                fixes.append(fix)
        elif cmd == "run":
            n = to_num(rest.split()[0]) if rest else None
            if n is not None:
                runs.append(int(n))
            else:
                unresolved.add(f"run {rest}")
        elif cmd == "minimize":
            minimize += 1
        elif cmd == "dump":
            parts = rest.split()
            if len(parts) >= 5:
                dumps.append({"id": parts[0], "group": parts[1], "style": parts[2],
                              "every": parts[3], "file": parts[4], "columns": " ".join(parts[5:])})
    # Variables defined in the script itself are not "unresolved".
    unresolved -= set(variables)
    total_steps = sum(runs) if runs else None
    dt = to_num(settings.get("timestep", ""))
    return {
        "settings": settings,
        "pair_coeff": pair_coeff,
        "fixes": fixes,
        "runs": runs,
        "total_steps": total_steps,
        "minimize": minimize,
        "dumps": dumps,
        "variables": variables,
        "unresolved": sorted(unresolved),
        "data_files": [settings[k].split()[0] for k in ("read_data", "read_restart") if k in settings],
        "simulated_time": (total_steps * dt) if (total_steps and dt) else None,
        "time_unit": {"metal": "ps", "real": "fs", "lj": "tau", "si": "s"}.get(settings.get("units", "")),
    }


def parse_log(head: str, tail: str) -> Dict:
    out: Dict = {"version": None, "n_atoms": None, "segments": [], "performance": None,
                 "wall_time_s": None, "errors": [], "warnings": 0, "finished": False}
    for line in head.splitlines():
        if out["version"] is None:
            m = LOG_VERSION_RE.match(line)
            if m:
                out["version"] = m.group(1)
        if out["n_atoms"] is None:
            m = ATOMS_RE.match(line)
            if m:
                out["n_atoms"] = int(m.group(1))
    for line in (head + "\n" + tail).splitlines():
        if line.startswith("WARNING"):
            out["warnings"] += 1
    lines = tail.splitlines()
    for i, line in enumerate(lines):
        m = LOOP_RE.match(line)
        if m:
            out["segments"].append({"loop_time_s": float(m.group(1)), "procs": int(m.group(2)),
                                    "steps": int(m.group(3)), "atoms": int(m.group(4))})
            out["n_atoms"] = int(m.group(4))
            continue
        m = PERF_RE.match(line)
        if m:
            out["performance"] = m.group(1).strip()
            continue
        m = WALL_RE.match(line)
        if m:
            h, mi, s = map(int, m.groups())
            out["wall_time_s"] = h * 3600 + mi * 60 + s
            out["finished"] = True
            continue
        if line.startswith("ERROR"):
            out["errors"].append(line.strip())
    if out["segments"]:
        out["procs"] = out["segments"][-1]["procs"]
    return out
