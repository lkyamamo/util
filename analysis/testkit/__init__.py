"""Shared helpers for analysis tests: synthetic systems with known answers,
LAMMPS file writers, and runners for scripts that are not importable modules.

Verified by analysis/testkit/test_testkit.py; see analysis/testkit/README.md.
"""

from .checks import Checker
from .lammps_io import read_dump, write_data, write_dump
from .scripts import load_script, run_script
from .systems import (
    dipole_series,
    harmonic_oscillators,
    ideal_gas,
    lattice_from_data,
    random_walk,
)

__all__ = [
    "Checker",
    "dipole_series",
    "harmonic_oscillators",
    "ideal_gas",
    "lattice_from_data",
    "load_script",
    "random_walk",
    "read_dump",
    "run_script",
    "write_data",
    "write_dump",
]
