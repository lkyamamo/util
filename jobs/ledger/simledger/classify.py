"""File categorisation by name. Content is never read here."""

from __future__ import annotations

import fnmatch
import re

VASP_INPUTS = {"INCAR", "POSCAR", "KPOINTS", "POTCAR", "ICONST", "vdw_kernel.bindat"}
VASP_OUTPUTS = {"OUTCAR", "OSZICAR", "vasprun.xml", "CONTCAR", "REPORT", "EIGENVAL",
                "IBZKPT", "PCDAT", "vaspout.h5", "vasp.out"}
VASP_BINARIES = {"WAVECAR", "CHGCAR", "CHG", "LOCPOT", "ELFCAR", "PROCAR", "DOSCAR", "WAVEDER"}

# (category, glob patterns) checked in order; first match wins.
RULES = [
    ("readme", ["README*", "readme*", "notes*", "NOTES*", "*.md"]),
    ("scheduler_out", ["STREAM_OUTPUT", "slurm-*.out", "*.o[0-9]*", "*.e[0-9]*"]),
    ("job_script", ["*.slurm", "*.pbs", "*.sbatch", "submit*.sh", "job*.sh"]),
    ("vasp_input", sorted(VASP_INPUTS) + ["*POTCAR*"]),
    ("vasp_output", sorted(VASP_OUTPUTS)),
    ("trajectory", ["XDATCAR", "*.dump", "dump.*", "*.lammpstrj", "*.custom", "*.dcd", "*.xtc", "*.nc"]),
    ("large_binary", sorted(VASP_BINARIES) + ["restart.*", "*.restart", "*.rst"]),
    ("lammps_log", ["log.lammps", "log.*.lammps", "log.lammps.*"]),
    ("lammps_input", ["in.*", "*.in", "*.lmp", "*.lammps", "*.input"]),
    ("structure", ["*.data", "data.*", "*.xyz", "*.extxyz", "*.cif", "*.vasp", "POSCAR*", "CONTCAR*"]),
    ("potential", ["*.usc", "*.eam", "*.eam.alloy", "*.eam.fs", "*.meam", "*.tersoff", "*.sw",
                   "ffield*", "*.reax", "*.pb", "*.pt", "*.yace", "*.mtp", "*.snapcoeff", "*.snapparam"]),
    ("analysis_script", ["*.py", "*.ipynb", "*.sh", "*.m", "*.jl", "*.R", "*.gnu", "*.plt"]),
    ("figure", ["*.png", "*.pdf", "*.svg", "*.jpg", "*.jpeg", "*.gif", "*.mp4"]),
    ("data", ["*.csv", "*.dat", "*.txt", "*.npy", "*.npz", "*.h5", "*.hdf5", "*.json", "*.conf"]),
    ("log", ["*.log", "*.out", "*.err"]),
]

_COMPILED = [(cat, [re.compile(fnmatch.translate(p)) for p in pats]) for cat, pats in RULES]


def classify(name: str) -> str:
    for cat, regexes in _COMPILED:
        if any(r.match(name) for r in regexes):
            return cat
    return "other"


def code_from_categories(categories) -> str:
    cats = set(categories)
    lammps = bool(cats & {"lammps_input", "lammps_log"})
    vasp = bool(cats & {"vasp_input", "vasp_output"})
    if lammps and vasp:
        return "both"
    if lammps:
        return "lammps"
    if vasp:
        return "vasp"
    return "unknown"
