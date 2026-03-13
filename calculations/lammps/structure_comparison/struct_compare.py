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
"""

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.lammps.data import LammpsData
from pymatgen.analysis.structure_matcher import StructureMatcher


# =============================================================================
# FILE PATHS — edit these
# =============================================================================

ORIG_FILE = "original.lammps"
CG_FILE   = "cg_relaxed.lammps"


# =============================================================================
# STEP 1: LOAD STRUCTURES
# =============================================================================
# LammpsData.from_file() parses the LAMMPS data file, then .structure converts
# it into a pymatgen Structure object which gives us the lattice, species, and
# fractional coordinates cleanly.

print("=" * 55)
print("  Structure Comparison: Original vs CG Relaxed")
print("=" * 55)

orig_lmp = LammpsData.from_file(ORIG_FILE, atom_style="atomic")
cg_lmp   = LammpsData.from_file(CG_FILE,   atom_style="atomic")

orig = orig_lmp.structure
cg   = cg_lmp.structure

N = len(orig)
assert len(orig) == len(cg), \
    f"Atom count mismatch: {len(orig)} vs {len(cg)}"

print(f"\nNumber of atoms : {N}")
print(f"Species present : {set(str(s.specie) for s in orig)}")


# =============================================================================
# STEP 2: VOLUME COMPARISON
# =============================================================================
# pymatgen Structure objects expose .volume directly (in Å³).
# ΔV/V tells you whether the potential over- or under-binds volumetrically.

V_orig    = orig.volume
V_cg      = cg.volume
dV_over_V = (V_cg - V_orig) / V_orig

print(f"\n--- Volume ---")
print(f"  Original : {V_orig:.4f} Å³")
print(f"  CG       : {V_cg:.4f} Å³")
print(f"  ΔV/V     : {dV_over_V:+.4f}  ({dV_over_V * 100:+.2f}%)")


# =============================================================================
# STEP 3: LATTICE PARAMETER COMPARISON
# =============================================================================
# The Lattice object exposes a, b, c (lengths in Å) and alpha, beta, gamma
# (angles in degrees). Comparing these individually tells you whether the
# volume change is isotropic or if the cell is also shearing/distorting.

lo = orig.lattice
lc = cg.lattice

print(f"\n--- Lattice Parameters ---")
print(f"  {'Param':<8} {'Original':>12} {'CG':>12} {'Δ':>12}")
print(f"  {'-'*46}")
for name, vo, vc in [
    ("a (Å)",   lo.a,     lc.a),
    ("b (Å)",   lo.b,     lc.b),
    ("c (Å)",   lo.c,     lc.c),
    ("α (°)",   lo.alpha, lc.alpha),
    ("β (°)",   lo.beta,  lc.beta),
    ("γ (°)",   lo.gamma, lc.gamma),
]:
    print(f"  {name:<8} {vo:>12.5f} {vc:>12.5f} {vc - vo:>+12.5f}")


# =============================================================================
# STEP 4: FRACTIONAL COORDINATE MSD
# =============================================================================
# site.frac_coords gives the fractional (scaled) coordinates for each atom,
# which are dimensionless and independent of cell size — this is what makes
# the MSD fairly comparable even though the volume has changed.
#
# Minimum image convention: subtracting np.round(delta) ensures we always
# measure the shortest displacement across periodic boundaries.

frac_orig = np.array([site.frac_coords for site in orig])
frac_cg   = np.array([site.frac_coords for site in cg])

delta = frac_cg - frac_orig
delta -= np.round(delta)  # minimum image convention

# --- add this ---
masses = np.array([site.specie.atomic_mass for site in orig])
com_shift = np.average(delta, axis=0, weights=masses)

print(f"\n--- Center of Mass Shift (fractional) ---")
print(f"  {com_shift}")
print(f"  magnitude: {np.linalg.norm(com_shift):.6f}")

delta -= com_shift
# ----------------

per_atom_disp = np.sqrt(np.sum(delta**2, axis=1))
msd           = np.mean(per_atom_disp**2)
mean_disp     = np.mean(per_atom_disp)

print(f"\n--- Fractional Coordinate Displacements ---")
print(f"  MSD (fractional)         : {msd:.6f}")
print(f"  sqrt(MSD) (fractional)   : {np.sqrt(msd):.6f}")
print(f"  Mean displacement (frac) : {mean_disp:.6f}")


# =============================================================================
# STEP 5: PER-SPECIES MSD BREAKDOWN
# =============================================================================
# If your structure has multiple species, it's useful to see whether one
# sublattice is more distorted than the other — e.g. the anion cage may
# relax more than the cation sites.

species = [str(site.specie) for site in orig]
unique_species = sorted(set(species))

if len(unique_species) > 1:
    print(f"\n--- Per-Species MSD (fractional) ---")
    for sp in unique_species:
        idx = [i for i, s in enumerate(species) if s == sp]
        sp_msd = np.mean(per_atom_disp[idx]**2)
        print(f"  {sp:<6}  MSD = {sp_msd:.6f}  (n={len(idx)})")


# =============================================================================
# STEP 6: STRUCTURE MATCHER
# =============================================================================
# StructureMatcher is pymatgen's symmetry-aware comparison tool. It finds the
# best mapping between the two structures (accounting for translations,
# rotations, and periodic images) and returns:
#   - fit():         True/False — are they the same structure type?
#   - get_rms_dist(): (rms, max_dist) normalized by the average nn distance,
#                     so it's comparable across different materials/cell sizes.

print(f"\n--- StructureMatcher (symmetry-aware) ---")
matcher = StructureMatcher()

fit = matcher.fit(orig, cg)
print(f"  Structures match : {fit}")

rms_result = matcher.get_rms_dist(orig, cg)
if rms_result is not None:
    rms, max_dist = rms_result
    print(f"  RMS distance (normalized)  : {rms:.6f}")
    print(f"  Max distance (normalized)  : {max_dist:.6f}")
else:
    print("  RMS distance : could not be computed (structures too dissimilar)")


# =============================================================================
# STEP 7: TOP 5 MOST DISPLACED ATOMS
# =============================================================================
# Useful for spotting localised failure — e.g. if the potential is getting
# one particular Wyckoff site badly wrong.

top5_idx = np.argsort(per_atom_disp)[::-1][:5]

print(f"\n--- Top 5 Most Displaced Atoms (fractional) ---")
print(f"  {'Index':>6}  {'Species':>8}  {'|Δs|':>10}")
for idx in top5_idx:
    print(f"  {idx:>6}  {species[idx]:>8}  {per_atom_disp[idx]:>10.6f}")

print("\nDone.")