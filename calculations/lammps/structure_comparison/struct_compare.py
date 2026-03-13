"""
compare_structures_pymatgen.py

Compares an original structure and a CG-relaxed structure, both in LAMMPS
data file format, using pymatgen. Computes:
  1. ΔV/V  — volumetric change
  2. Lattice parameter comparison (a, b, c, α, β, γ)
  3. MSD in fractional coordinates (volume-agnostic internal relaxation)
  4. StructureMatcher — symmetry-aware structural similarity
  5. Per-species MSD breakdown
  6. Top 5 most displaced atoms

Usage:
  python compare_structures_pymatgen.py                        # default output file
  python compare_structures_pymatgen.py --output my_results.txt
  python compare_structures_pymatgen.py --no-export            # print only
"""

import argparse
import sys
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.lammps.data import LammpsData
from pymatgen.analysis.structure_matcher import StructureMatcher


# =============================================================================
# ARGUMENT PARSING
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare original and CG-relaxed LAMMPS structures."
    )
    parser.add_argument(
        "orig_file",
        nargs="?",
        default="original.lammps",
        help="Path to original LAMMPS data file (default: original.lammps)"
    )
    parser.add_argument(
        "cg_file",
        nargs="?",
        default="cg_relaxed.lammps",
        help="Path to CG-relaxed LAMMPS data file (default: cg_relaxed.lammps)"
    )
    parser.add_argument(
        "--output", "-o",
        default="structure_comparison.txt",
        help="Output file path (default: structure_comparison.txt)"
    )
    parser.add_argument(
        "--no-export",
        action="store_true",
        default=False,
        help="Disable file export, print to terminal only"
    )
    return parser.parse_args()


# =============================================================================
# OUTPUT HANDLER
# =============================================================================
# Wraps both terminal and file output so the analysis code only
# calls one print function regardless of export settings.

class OutputHandler:
    def __init__(self, export: bool, filepath: str):
        self.export   = export
        self.filepath = filepath
        self.file     = None

    def __enter__(self):
        if self.export:
            self.file = open(self.filepath, "w")
            print(f"Exporting results to: {self.filepath}")
        return self

    def __exit__(self, *args):
        if self.file:
            self.file.close()

    def print(self, *args, **kwargs):
        # Always print to terminal
        print(*args, **kwargs)
        # Also write to file if export is on
        if self.file:
            print(*args, **kwargs, file=self.file)


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_analysis(orig_file: str, cg_file: str, out: OutputHandler):

    # --- Load structures ---
    out.print("=" * 55)
    out.print("  Structure Comparison: Original vs CG Relaxed")
    out.print("=" * 55)

    orig_lmp = LammpsData.from_file(orig_file, atom_style="atomic")
    cg_lmp   = LammpsData.from_file(cg_file,   atom_style="atomic")

    orig = orig_lmp.structure
    cg   = cg_lmp.structure

    N = len(orig)
    assert len(orig) == len(cg), \
        f"Atom count mismatch: {len(orig)} vs {len(cg)}"

    out.print(f"\nOriginal file   : {orig_file}")
    out.print(f"CG file         : {cg_file}")
    out.print(f"Number of atoms : {N}")
    out.print(f"Species present : {set(str(s.specie) for s in orig)}")

    # --- Volume ---
    V_orig    = orig.volume
    V_cg      = cg.volume
    dV_over_V = (V_cg - V_orig) / V_orig

    out.print(f"\n--- Volume ---")
    out.print(f"  Original : {V_orig:.4f} Å³")
    out.print(f"  CG       : {V_cg:.4f} Å³")
    out.print(f"  ΔV/V     : {dV_over_V:+.4f}  ({dV_over_V * 100:+.2f}%)")

    # --- Lattice parameters ---
    lo = orig.lattice
    lc = cg.lattice

    out.print(f"\n--- Lattice Parameters ---")
    out.print(f"  {'Param':<8} {'Original':>12} {'CG':>12} {'Δ':>12}")
    out.print(f"  {'-'*46}")
    for name, vo, vc in [
        ("a (Å)",   lo.a,     lc.a),
        ("b (Å)",   lo.b,     lc.b),
        ("c (Å)",   lo.c,     lc.c),
        ("α (°)",   lo.alpha, lc.alpha),
        ("β (°)",   lo.beta,  lc.beta),
        ("γ (°)",   lo.gamma, lc.gamma),
    ]:
        out.print(f"  {name:<8} {vo:>12.5f} {vc:>12.5f} {vc - vo:>+12.5f}")

    # --- Fractional coordinate MSD with COM correction ---
    frac_orig = np.array([site.frac_coords for site in orig])
    frac_cg   = np.array([site.frac_coords for site in cg])

    delta = frac_cg - frac_orig
    delta -= np.round(delta)

    masses     = np.array([site.specie.atomic_mass for site in orig])
    com_shift  = np.average(delta, axis=0, weights=masses)

    out.print(f"\n--- Center of Mass Shift (fractional) ---")
    out.print(f"  shift     : {com_shift}")
    out.print(f"  magnitude : {np.linalg.norm(com_shift):.6f}")

    delta -= com_shift

    per_atom_disp = np.sqrt(np.sum(delta**2, axis=1))
    msd           = np.mean(per_atom_disp**2)
    mean_disp     = np.mean(per_atom_disp)

    out.print(f"\n--- Fractional Coordinate Displacements (COM-corrected) ---")
    out.print(f"  MSD (fractional)         : {msd:.6f}")
    out.print(f"  sqrt(MSD) (fractional)   : {np.sqrt(msd):.6f}")
    out.print(f"  Mean displacement (frac) : {mean_disp:.6f}")

    # --- Per-species MSD ---
    species        = [str(site.specie) for site in orig]
    unique_species = sorted(set(species))

    if len(unique_species) > 1:
        out.print(f"\n--- Per-Species MSD (fractional, COM-corrected) ---")
        for sp in unique_species:
            idx    = [i for i, s in enumerate(species) if s == sp]
            sp_msd = np.mean(per_atom_disp[idx]**2)
            out.print(f"  {sp:<6}  MSD = {sp_msd:.6f}  (n={len(idx)})")

    # --- StructureMatcher ---
    out.print(f"\n--- StructureMatcher (symmetry-aware) ---")
    matcher = StructureMatcher()

    fit = matcher.fit(orig, cg)
    out.print(f"  Structures match : {fit}")

    rms_result = matcher.get_rms_dist(orig, cg)
    if rms_result is not None:
        rms, max_dist = rms_result
        out.print(f"  RMS distance (normalized)  : {rms:.6f}")
        out.print(f"  Max distance (normalized)  : {max_dist:.6f}")
    else:
        out.print("  RMS distance : could not be computed (structures too dissimilar)")

    # --- Top 5 most displaced atoms ---
    top5_idx = np.argsort(per_atom_disp)[::-1][:5]

    out.print(f"\n--- Top 5 Most Displaced Atoms (fractional, COM-corrected) ---")
    out.print(f"  {'Index':>6}  {'Species':>8}  {'|Δs|':>10}")
    for idx in top5_idx:
        out.print(f"  {idx:>6}  {species[idx]:>8}  {per_atom_disp[idx]:>10.6f}")

    out.print("\nDone.")


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    args = parse_args()

    with OutputHandler(
        export   = not args.no_export,
        filepath = args.output
    ) as out:
        run_analysis(args.orig_file, args.cg_file, out)