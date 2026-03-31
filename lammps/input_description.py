#!/usr/bin/env python3
"""
lammps_description.py
==================
Reads a LAMMPS input script and writes a human-readable description
covering the simulation setup and each logical step/phase.

Usage:
    python lammps_description.py <input_file> [output_file]

If no output_file is given, the description is written to
<input_file>.description.txt
"""

import sys
import re
from pathlib import Path
from datetime import datetime


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

SECTION_KEYWORDS = {
    # Initialisation
    "units", "dimension", "boundary", "atom_style", "atom_modify",
    "newton", "processors",
    # Geometry / topology
    "read_data", "read_restart", "read_dump", "lattice", "region",
    "create_box", "create_atoms", "delete_atoms", "replicate",
    # Force-field
    "pair_style", "pair_coeff", "pair_modify",
    "bond_style", "bond_coeff",
    "angle_style", "angle_coeff",
    "dihedral_style", "dihedral_coeff",
    "improper_style", "improper_coeff",
    "kspace_style", "kspace_modify",
    "special_bonds",
    # Groups / computes / fixes
    "group",
    "compute",
    "fix",
    "unfix",
    "fix_modify",
    # Variables & thermodynamics
    "variable",
    "thermo", "thermo_style", "thermo_modify",
    # Output
    "dump", "dump_modify", "undump",
    "restart", "write_restart", "write_data", "write_dump",
    "log",
    # Run control
    "minimize", "run", "rerun",
    "timestep", "reset_timestep",
    "velocity",
    "displace_atoms", "change_box", "set",
    # Misc
    "include", "jump", "label", "next", "print", "quit", "echo",
    "shell", "plugin", "package", "suffix", "if",
    "mass", "dielectric",
}

# Keywords that mark a new simulation "step" or phase
STEP_MARKERS = {"minimize", "run", "rerun"}

# Friendly descriptions for common commands
FRIENDLY = {
    "units":          lambda v: f"Unit system: {v}",
    "dimension":      lambda v: f"Spatial dimension: {v}D",
    "boundary":       lambda v: f"Boundary conditions: {v}",
    "atom_style":     lambda v: f"Atom style: {v}",
    "newton":         lambda v: f"Newton's 3rd-law flag: {v}",
    "processors":     lambda v: f"MPI processor grid: {v}",
    "read_data":      lambda v: f"Reads topology/coordinates from: {v.split()[0]}",
    "read_restart":   lambda v: f"Restarts from checkpoint: {v.split()[0]}",
    "lattice":        lambda v: f"Lattice: {v}",
    "region":         lambda v: f"Defines region '{v.split()[0]}' ({' '.join(v.split()[1:])})",
    "create_box":     lambda v: f"Creates simulation box with {v.split()[0]} atom type(s)",
    "create_atoms":   lambda v: f"Populates atoms of type {v.split()[0]}",
    "replicate":      lambda v: f"Replicates the box {v} times",
    "pair_style":     lambda v: f"Pairwise interaction: {v}",
    "pair_coeff":     lambda v: f"Pair coefficients: {v}",
    "bond_style":     lambda v: f"Bond style: {v}",
    "angle_style":    lambda v: f"Angle style: {v}",
    "dihedral_style": lambda v: f"Dihedral style: {v}",
    "improper_style": lambda v: f"Improper style: {v}",
    "kspace_style":   lambda v: f"Long-range electrostatics: {v}",
    "special_bonds":  lambda v: f"Special bond scaling: {v}",
    "mass":           lambda v: f"Mass assignment: {v}",
    "group":          lambda v: f"Group '{v.split()[0]}': {' '.join(v.split()[1:])}",
    "compute":        lambda v: f"Compute '{v.split()[0]}': {' '.join(v.split()[1:])}",
    "fix":            lambda v: _describe_fix(v),
    "unfix":          lambda v: f"Removes fix '{v}'",
    "variable":       lambda v: f"Variable '{v.split()[0]}' = {' '.join(v.split()[1:])}",
    "thermo":         lambda v: f"Thermo output every {v} steps",
    "thermo_style":   lambda v: f"Thermo columns: {v}",
    "dump":           lambda v: _describe_dump(v),
    "undump":         lambda v: f"Removes dump '{v}'",
    "restart":        lambda v: f"Writes restart every {v.split()[0]} step(s)",
    "write_restart":  lambda v: f"Writes restart file: {v}",
    "write_data":     lambda v: f"Writes data file: {v}",
    "timestep":       lambda v: f"Integration timestep: {v}",
    "reset_timestep": lambda v: f"Resets step counter to {v}",
    "velocity":       lambda v: f"Initialises velocities: {v}",
    "minimize":       lambda v: f"Energy minimisation ({v})",
    "run":            lambda v: f"MD run for {v} steps",
    "rerun":          lambda v: f"Re-runs trajectory: {v}",
    "log":            lambda v: f"Log file: {v}",
    "include":        lambda v: f"Includes sub-script: {v}",
    "print":          lambda v: f'Prints message: {v.strip(chr(34)).strip(chr(39))}',
    "echo":           lambda v: f"Echo mode: {v}",
    "dielectric":     lambda v: f"Dielectric constant: {v}",
    "change_box":     lambda v: f"Modifies box: {v}",
    "displace_atoms": lambda v: f"Displaces atoms: {v}",
}


def _describe_fix(args: str) -> str:
    parts = args.split()
    if len(parts) < 3:
        return f"Fix: {args}"
    fid, grp, style = parts[0], parts[1], parts[2]
    rest = " ".join(parts[3:])
    style_notes = {
        "nvt":   "NVT thermostat (constant N, V, T)",
        "npt":   "NPT barostat/thermostat (constant N, P, T)",
        "nve":   "NVE integrator (microcanonical)",
        "nvt/sllod": "NVT-SLLOD for shear flow",
        "langevin": "Langevin thermostat",
        "nose-hoover": "Nosé–Hoover thermostat",
        "shake": "SHAKE bond/angle constraint",
        "rigid": "Rigid-body integrator",
        "deform": "Box deformation",
        "ave/time": "Time-averaged output",
        "ave/chunk": "Chunk-averaged output",
        "ave/atom": "Per-atom averaged output",
        "spring": "Harmonic spring restraint",
        "wall/lj93": "LJ 9-3 wall",
        "momentum": "Zero linear/angular momentum",
        "recenter": "Re-centers the simulation",
    }
    note = style_notes.get(style.lower(), style)
    return f"Fix '{fid}' on group '{grp}': {note}" + (f" ({rest})" if rest else "")


def _describe_dump(args: str) -> str:
    parts = args.split()
    if len(parts) < 4:
        return f"Dump: {args}"
    did, grp, style, freq = parts[0], parts[1], parts[2], parts[3]
    fname = parts[4] if len(parts) > 4 else "?"
    return (f"Dump '{did}' — group '{grp}', format {style}, "
            f"every {freq} steps → {fname}")


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

class LammpsLine:
    """Represents one logical (possibly continued) LAMMPS command line."""
    def __init__(self, lineno: int, raw: str, keyword: str, args: str):
        self.lineno  = lineno
        self.raw     = raw
        self.keyword = keyword.lower()
        self.args    = args.strip()

    def describe(self) -> str:
        fn = FRIENDLY.get(self.keyword)
        if fn:
            try:
                return fn(self.args)
            except Exception:
                pass
        return f"{self.keyword.capitalize()}: {self.args}"


def strip_comment(line: str) -> str:
    """Remove inline # comments (respecting quoted strings is hard; we keep it simple)."""
    idx = line.find("#")
    if idx >= 0:
        return line[:idx]
    return line


def parse_lammps(path: Path) -> list[LammpsLine]:
    """Read a LAMMPS script and return a list of LammpsLine objects."""
    raw_lines = path.read_text(errors="replace").splitlines()

    logical_lines: list[LammpsLine] = []
    buffer = ""
    start_lineno = 1

    for i, raw in enumerate(raw_lines, start=1):
        stripped = strip_comment(raw).rstrip()
        if stripped.endswith("&"):           # line continuation
            buffer += " " + stripped[:-1].strip()
            if not buffer.strip():
                start_lineno = i + 1
        else:
            buffer += " " + stripped
            combined = buffer.strip()
            buffer = ""
            if not combined:
                start_lineno = i + 1
                continue
            tokens = combined.split()
            if not tokens:
                start_lineno = i + 1
                continue
            kw = tokens[0].lower()
            args = " ".join(tokens[1:])
            logical_lines.append(LammpsLine(start_lineno, combined, kw, args))
            start_lineno = i + 1

    return logical_lines


# ---------------------------------------------------------------------------
# Categoriser
# ---------------------------------------------------------------------------

class Section:
    def __init__(self, name: str):
        self.name  = name
        self.lines: list[LammpsLine] = []

    def add(self, line: LammpsLine):
        self.lines.append(line)


def categorise(parsed: list[LammpsLine]) -> list[Section]:
    """Split parsed commands into logical sections."""
    sections: list[Section] = []
    current = Section("Initialisation & Global Settings")
    step_count = 0

    for line in parsed:
        kw = line.keyword

        if kw in STEP_MARKERS:
            # Close prior section if non-empty
            if current.lines:
                sections.append(current)
            step_count += 1
            label = {
                "minimize": f"Step {step_count}: Energy Minimisation",
                "run":      f"Step {step_count}: MD Run",
                "rerun":    f"Step {step_count}: Trajectory Re-run",
            }[kw]
            current = Section(label)
            current.add(line)
            sections.append(current)
            current = Section(f"Post-Step {step_count} Setup")
        else:
            current.add(line)

    if current.lines:
        sections.append(current)

    # Drop empty trailing setup sections
    return [s for s in sections if s.lines]


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

def write_description(parsed: list[LammpsLine],
                      sections: list[Section],
                      src_path: Path,
                      out_path: Path):

    lines_out: list[str] = []
    W = lines_out.append

    def rule(char="=", width=72):
        W(char * width)

    def header(text: str):
        rule()
        W(text)
        rule()
        W("")

    def sub_header(text: str):
        W(text)
        rule("-", 40)

    # ---- Title block ----
    W(f"LAMMPS Input File Description")
    W(f"Generated : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    W(f"Source    : {src_path.resolve()}")
    W(f"Lines     : {sum(1 for _ in src_path.open(errors='replace'))}")
    W(f"Commands  : {len(parsed)}")
    W("")

    # ---- Quick summary ----
    header("QUICK SUMMARY")
    run_cmds    = [l for l in parsed if l.keyword == "run"]
    min_cmds    = [l for l in parsed if l.keyword == "minimize"]
    fix_cmds    = [l for l in parsed if l.keyword == "fix"]
    dump_cmds   = [l for l in parsed if l.keyword == "dump"]
    read_cmds   = [l for l in parsed if l.keyword in ("read_data", "read_restart", "read_dump")]
    ts_cmds     = [l for l in parsed if l.keyword == "timestep"]

    if read_cmds:
        W(f"  Input geometry : {read_cmds[0].args.split()[0]}")
    else:
        W("  Input geometry : (created internally)")

    units_cmds = [l for l in parsed if l.keyword == "units"]
    if units_cmds:
        W(f"  Unit system    : {units_cmds[0].args}")

    if ts_cmds:
        ts = ts_cmds[-1].args
        W(f"  Timestep       : {ts}")
        if run_cmds:
            try:
                total_steps = sum(int(r.args.split()[0]) for r in run_cmds)
                W(f"  Total MD steps : {total_steps:,}  "
                  f"({total_steps * float(ts):.4g} time units)")
            except ValueError:
                pass

    W(f"  MD run blocks  : {len(run_cmds)}")
    W(f"  Minimisations  : {len(min_cmds)}")
    W(f"  Active fixes   : {len(fix_cmds)}")
    W(f"  Dump streams   : {len(dump_cmds)}")
    W("")

    # ---- Section-by-section ----
    header("DETAILED WALKTHROUGH")

    for sec in sections:
        sub_header(sec.name)
        for line in sec.lines:
            desc = line.describe()
            W(f"  [line {line.lineno:>4}]  {desc}")
        W("")

    # ---- Force-field summary ----
    ff_kws = {"pair_style", "pair_coeff", "bond_style", "bond_coeff",
              "angle_style", "angle_coeff", "dihedral_style", "dihedral_coeff",
              "kspace_style", "special_bonds"}
    ff_lines = [l for l in parsed if l.keyword in ff_kws]
    if ff_lines:
        header("FORCE-FIELD OVERVIEW")
        for line in ff_lines:
            W(f"  {line.describe()}")
        W("")

    # ---- Fix summary ----
    if fix_cmds:
        header("FIX SUMMARY")
        for line in fix_cmds:
            W(f"  {line.describe()}")
        W("")

    # ---- Dump / output summary ----
    output_kws = {"dump", "restart", "write_restart", "write_data", "log"}
    output_lines = [l for l in parsed if l.keyword in output_kws]
    if output_lines:
        header("OUTPUT FILES")
        for line in output_lines:
            W(f"  {line.describe()}")
        W("")

    # ---- Raw script echo ----
    header("RAW LAMMPS SCRIPT")
    raw_text = src_path.read_text(errors="replace")
    for i, raw_line in enumerate(raw_text.splitlines(), start=1):
        W(f"  {i:>4}  {raw_line}")
    W("")

    out_path.write_text("\n".join(lines_out))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    src = Path(sys.argv[1])
    if not src.is_file():
        print(f"Error: '{src}' not found.")
        sys.exit(1)

    out = Path(sys.argv[2]) if len(sys.argv) > 2 else src.with_suffix(".description.txt")

    print(f"Parsing  : {src}")
    parsed   = parse_lammps(src)
    sections = categorise(parsed)
    write_description(parsed, sections, src, out)
    print(f"Written  : {out}")
    print(f"Commands : {len(parsed)}")
    print(f"Sections : {len(sections)}")


if __name__ == "__main__":
    main()