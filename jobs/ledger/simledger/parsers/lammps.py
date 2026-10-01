"""LAMMPS input scripts, log.lammps (head/tail for status, full stream for the protocol) and data files."""

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


# ---------------------------------------------------------------------------
# Full log stream: the protocol as it actually ran
# ---------------------------------------------------------------------------
#
# With `echo log/both` LAMMPS writes every command to the log, and when a line
# contains ${var} / $(expr) it writes the substituted line right after it. So
# lines still holding a $-reference are skipped and the substituted copy is
# used: these are the values that ran, including ones passed with -var.

COMMANDS = {"units", "atom_style", "boundary", "dimension", "pair_style", "pair_coeff", "mass",
            "timestep", "fix", "unfix", "run", "minimize", "velocity", "read_data", "read_restart",
            "write_data", "reset_timestep", "thermo_style", "kspace_style", "bond_style",
            "angle_style", "dihedral_style", "dump", "include", "replicate", "create_atoms", "delete_atoms"}
DOLLAR_RE = re.compile(r"\$(\{|\(|[A-Za-z])")
BOX_RE = re.compile(r"^\s*(orthogonal|triclinic) box = \(([^)]*)\) to \(([^)]*)\)(?: with tilt \(([^)]*)\))?")
THERMO_KEEP = ("Temp", "Press", "Volume", "Density", "PotEng", "KinEng", "TotEng", "Lx", "Ly", "Lz",
               "Pxx", "Pyy", "Pzz")

THERMOSTATS = {"langevin", "temp/berendsen", "temp/rescale", "temp/csvr", "temp/csld"}
BAROSTATS = {"press/berendsen"}


def _floats(text: str) -> List[float]:
    out = []
    for tok in text.replace(",", " ").split():
        v = to_num(tok)
        if v is not None:
            out.append(v)
    return out


def _fix_conditions(style: str, args: List[str]) -> Dict:
    """Target T/P from an integrator/thermostat/barostat fix."""
    out: Dict = {}

    def after(kw, n):
        if kw in args:
            i = args.index(kw)
            vals = [to_num(a) for a in args[i + 1:i + 1 + n]]
            return vals if all(v is not None for v in vals) else None
        return None

    if style in ("nvt", "npt", "nvt/sllod") or style.startswith("rigid/n"):
        t = after("temp", 3)
        if t:
            out["T"] = (t[0], t[1])
    if style in ("npt", "nph", "press/berendsen") or style.startswith("rigid/np"):
        for kw in ("iso", "aniso", "tri", "x", "y", "z"):
            p = after(kw, 3 if style != "press/berendsen" else 2)
            if p:
                out["P"] = (p[0], p[1])
                out["P_coupling"] = kw
                break
    if style in THERMOSTATS and len(args) >= 2:
        t0, t1 = to_num(args[0]), to_num(args[1])
        if t0 is not None and t1 is not None:
            out["T"] = (t0, t1)
    return out


def _ensemble(fixes: Dict[str, Dict]) -> str:
    styles = {f["style"] for f in fixes.values()}
    base = None
    if styles & {"npt"} or any(s.startswith("rigid/np") for s in styles):
        base = "NPT"
    elif "nph" in styles:
        base = "NPH" if not styles & THERMOSTATS else "NPT"
    elif styles & {"nvt", "nvt/sllod"} or any(s.startswith("rigid/nvt") for s in styles):
        base = "NVT"
    elif styles & {"nve", "nve/limit", "rigid/nve", "rigid"}:
        thermo = sorted(styles & THERMOSTATS)
        base = f"NVT ({thermo[0]})" if thermo else "NVE"
        if styles & BAROSTATS:
            base = base.replace("NVT", "NPT") + " + press/berendsen"
    extra = [s for s in ("deform", "wall/reflect", "spring", "efield", "momentum") if s in styles]
    if base is None:
        base = "no integrator" if not styles else "+".join(sorted(styles))
    return base + ("".join(f" + {e}" for e in extra) if extra else "")


def parse_log_stream(lines, substituted: bool = True) -> Dict:
    """Walk log (or input) lines: settings, box, masses and the run protocol.

    substituted=True for logs (skip lines that still hold $-references, their
    substituted copy follows). For an input script pass False: references are
    kept and reported as unresolved."""
    out: Dict = {"version": None, "units": None, "atom_style": None, "n_atoms": None, "box": None,
                 "masses": {}, "pair_style": None, "pair_coeff": [], "data_files": [], "write_data": [],
                 "timestep": None, "segments": [], "unresolved": set(), "has_commands": False,
                 "velocity_T": None}
    fixes: Dict[str, Dict] = {}
    seg = None          # segment being run (after its run/minimize command)
    header: List[str] = []
    rows: List[List[float]] = []

    def close_segment(loop_line=None):
        nonlocal seg, header, rows
        if seg is None:
            return
        if rows and header:
            half = rows[len(rows) // 2:]
            stats = {}
            for i, col in enumerate(header):
                if col in THERMO_KEEP or col.startswith(("v_", "c_")):
                    vals = [r[i] for r in half if i < len(r)]
                    if vals:
                        mean = sum(vals) / len(vals)
                        std = (sum((v - mean) ** 2 for v in vals) / len(vals)) ** 0.5
                        stats[col] = {"mean": mean, "std": std, "first": rows[0][i], "last": rows[-1][i]}
            seg["thermo"] = stats
            seg["thermo_rows"] = len(rows)
        if loop_line:
            m = LOOP_RE.match(loop_line)
            if m:
                seg["loop_time_s"] = float(m.group(1))
                seg["steps_done"] = int(m.group(3))
                out["n_atoms"] = int(m.group(4))
        out["segments"].append(seg)
        seg, header, rows = None, [], []

    for raw in lines:
        line = raw.rstrip("\n")
        s = line.strip()
        if not s:
            continue
        # thermo block
        if seg is not None:
            if s.startswith("Step ") or s == "Step":
                header, rows = s.split(), []
                continue
            if header and (s[0].isdigit() or s[0] == "-"):
                vals = _floats(s)
                if len(vals) == len(header):
                    rows.append(vals)
                    continue
            if s.startswith("Loop time"):
                close_segment(s)
                continue
        if out["version"] is None:
            m = LOG_VERSION_RE.match(s)
            if m:
                out["version"] = m.group(1)
                continue
        m = BOX_RE.match(line)
        if m:
            lo, hi = _floats(m.group(2)), _floats(m.group(3))
            tilt = _floats(m.group(4)) if m.group(4) else None
            if len(lo) == 3 and len(hi) == 3:
                box = {"lengths": [h - l for l, h in zip(lo, hi)], "tilt": tilt, "kind": m.group(1)}
                if out["box"] is None:
                    out["box"] = box
                out["box_last"] = box
            continue
        m = ATOMS_RE.match(line)
        if m and out["n_atoms"] is None:
            out["n_atoms"] = int(m.group(1))
            continue
        code, _, comment = s.partition("#")
        tokens = code.split()
        if not tokens or tokens[0] not in COMMANDS:
            continue
        if DOLLAR_RE.search(code):
            if substituted:
                continue
            out["unresolved"].update(a or b for a, b in re.findall(r"\$\{(\w+)\}|\$(\w)", code))
        cmd, args = tokens[0], tokens[1:]
        if args and (args[0] == "CPU" or (len(args) > 1 and args[1] == "=")):
            continue   # LAMMPS's own report, e.g. "  read_data CPU = 0.13 secs", not a command
        out["has_commands"] = True
        if cmd in ("units", "atom_style", "pair_style") and args:
            out[cmd] = " ".join(args) if cmd == "pair_style" else args[0]
        elif cmd == "pair_coeff":
            out["pair_coeff"].append(" ".join(args))
        elif cmd == "mass" and len(args) >= 2:
            m_val = to_num(args[1])
            if m_val is not None:
                out["masses"][args[0]] = {"mass": m_val, "label": comment.strip() or None}
        elif cmd == "timestep" and args:
            out["timestep"] = to_num(args[0])
        elif cmd == "read_data" and args:
            out["data_files"].append(args[0])
        elif cmd == "read_restart" and args:   # continues the same system; not another structure
            out.setdefault("restart_files", []).append(args[0])
        elif cmd == "write_data" and args:
            out["write_data"].append(args[0])
        elif cmd == "replicate" and len(args) >= 3:
            f = [to_num(a) for a in args[:3]]
            if all(f):
                out["replicate"] = (out.get("replicate") or 1) * int(f[0] * f[1] * f[2])
        elif cmd in ("create_atoms", "delete_atoms"):
            out["atoms_changed"] = True
        elif cmd == "velocity" and "create" in args:
            i = args.index("create")
            if i + 1 < len(args):
                out["velocity_T"] = to_num(args[i + 1])
        elif cmd == "fix" and len(args) >= 3:
            fixes[args[0]] = {"style": args[2], "group": args[1], "args": args[3:],
                              **_fix_conditions(args[2], args[3:])}
        elif cmd == "unfix" and args:
            fixes.pop(args[0], None)
        elif cmd in ("run", "minimize"):
            close_segment()
            if "box_at_run" not in out:     # the box the simulation starts from (after replicate etc.)
                out["box_at_run"] = out.get("box_last") or out.get("box")
            steps = to_num(args[0]) if (cmd == "run" and args) else (to_num(args[3]) if len(args) >= 4 else None)
            T = next((f["T"] for f in fixes.values() if "T" in f), None)
            P = next((f["P"] for f in fixes.values() if "P" in f), None)
            seg = {"kind": cmd, "steps": int(steps) if steps is not None else None,
                   "timestep": out["timestep"], "ensemble": "minimize" if cmd == "minimize" else _ensemble(fixes),
                   "T": T, "P": P, "fixes": sorted(f"{k}:{v['style']}" for k, v in fixes.items())}
            if not substituted:   # input script: no thermo/loop lines follow
                close_segment()
    close_segment()
    out["unresolved"] = sorted(out["unresolved"])
    return out


def collapse_segments(segments: List[Dict]) -> List[Dict]:
    """Merge consecutive identical segments (e.g. 128 chunks of the same NVT run)."""
    out: List[Dict] = []
    for s in segments:
        key = (s["kind"], s["ensemble"], s.get("T"), s.get("P"), s.get("timestep"), s.get("steps"))
        if out and out[-1]["_key"] == key:
            prev = out[-1]
            prev["count"] += 1
            for col, st in (s.get("thermo") or {}).items():
                acc = prev["_acc"].setdefault(col, [])
                acc.append(st["mean"])
            continue
        row = dict(s, count=1, _key=key, _acc={c: [st["mean"]] for c, st in (s.get("thermo") or {}).items()})
        out.append(row)
    for row in out:
        thermo = dict(row.get("thermo") or {})
        for col, means in row.pop("_acc").items():
            if col in thermo:
                thermo[col] = dict(thermo[col], mean=sum(means) / len(means))
        row["thermo"] = thermo
        row.pop("_key")
    return out


# ---------------------------------------------------------------------------
# Data files: header, masses, and atoms per type
# ---------------------------------------------------------------------------

TYPE_COLUMN = {"atomic": 1, "charge": 1, "dipole": 1, "sphere": 1, "ellipsoid": 1, "line": 1, "tri": 1,
               "full": 2, "molecular": 2, "bond": 2, "angle": 2, "template": 2}
HEADER_RE = re.compile(r"^\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([xyz])lo\s+[xyz]hi")


def parse_data_file(lines, atom_style: str = None, count_types: bool = True) -> Dict:
    """Header (atoms, types, box, masses) and, if count_types, atoms of each type."""
    out: Dict = {"n_atoms": None, "n_types": None, "box": {}, "tilt": None, "masses": {},
                 "atoms_style": None, "type_counts": None}
    section = None
    counts: Dict[str, int] = {}
    col = None
    remaining = None
    for raw in lines:
        s = raw.strip()
        if section is None or section == "header":
            if not s:
                continue
            m = re.match(r"^(\d+)\s+atoms\s*$", s)
            if m:
                out["n_atoms"] = int(m.group(1))
                continue
            m = re.match(r"^(\d+)\s+atom types\s*$", s)
            if m:
                out["n_types"] = int(m.group(1))
                continue
            m = HEADER_RE.match(s)
            if m:
                out["box"][m.group(3)] = float(m.group(2)) - float(m.group(1))
                continue
            if s.endswith("xy xz yz"):
                out["tilt"] = _floats(s)[:3]
                continue
        word = s.split("#")[0].strip()
        if word in ("Masses", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers",
                    "Pair Coeffs", "Bond Coeffs", "Angle Coeffs"):
            if section == "Atoms" and not count_types:
                break
            section = word
            if word == "Atoms":
                out["atoms_style"] = s.partition("#")[2].strip() or atom_style
                col = TYPE_COLUMN.get((out["atoms_style"] or "atomic").split()[0], 1)
                remaining = out["n_atoms"]
                if not count_types:
                    break
            elif section != "Masses" and out["type_counts"] is None and counts:
                break
            continue
        if section == "Masses" and s:
            parts = s.partition("#")
            toks = parts[0].split()
            if len(toks) >= 2 and to_num(toks[1]) is not None:
                out["masses"][toks[0]] = {"mass": to_num(toks[1]), "label": parts[2].strip() or None}
        elif section == "Atoms" and s:
            toks = s.split()
            if len(toks) > col:
                counts[toks[col]] = counts.get(toks[col], 0) + 1
                if remaining is not None:
                    remaining -= 1
                    if remaining == 0:
                        break
    if counts and (out["n_atoms"] is None or sum(counts.values()) == out["n_atoms"]):
        out["type_counts"] = counts
    return out


# Standard atomic masses, for naming types that carry no element label.
ELEMENT_MASSES = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007,
    "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085,
    "P": 30.974, "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078, "Ti": 47.867,
    "Cr": 51.996, "Fe": 55.845, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Ge": 72.630, "Zr": 91.224,
    "Ag": 107.87, "Sn": 118.71, "Pt": 195.08, "Au": 196.97, "Pb": 207.2,
}


def guess_element(mass: float, tol: float = 0.1):
    best = min(ELEMENT_MASSES.items(), key=lambda kv: abs(kv[1] - mass))
    return best[0] if abs(best[1] - mass) <= tol else None
