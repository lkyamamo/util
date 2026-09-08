#!/usr/bin/env python3
"""
structure_io.py

Reading, editing, and writing the two structure formats md_setup's path
generators work in - LAMMPS data files and VASP POSCARs - through a single ASE
Atoms object, so a caller can stay format-blind and pick the format at run time.

This module lives at the md_setup root so scripts in the per-tool subdirectories
can import it:

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import structure_io

Why ASE and not pymatgen
------------------------
lammps_data.py (still used by hand_place/) is built on pymatgen's LammpsData and
measures distances per-axis, which is wrong on a skewed cell. ASE reads both
formats into one object and its minimum-image helpers go through the cell, so a
single implementation here is correct for any lattice, orthogonal or triclinic.

Conventions shared by both formats
----------------------------------
- An atom "id" is its 1-based position in the file's coordinate/Atoms block.
  Adding atoms can renumber a POSCAR (see `add_water`), so ids returned for a
  written frame are not necessarily the ids the input used; `add_water` hands
  back the map between them.
- Positions are Cartesian everywhere in this module's interface. Fractional
  coordinates are an implementation detail of the cell operations.
- Nothing about the motion of the input is carried into the output. Velocities
  are dropped and LAMMPS image flags are neither read nor written: every frame is
  an independent starting point for its own relaxation, so a velocity or an image
  count inherited from the run that produced the input has no meaning in it.
  Positions are used exactly as the file states them, and new atoms are wrapped
  into the cell.

What each format keeps that the other has no notion of
------------------------------------------------------
  LAMMPS   numeric atom types (with the Masses table), per-atom charges,
           molecule ids, and the "Atoms # <style>" comment.
  POSCAR   selective dynamics, the Direct/Cartesian coordinate mode, and a
           species line whose order the POTCAR has to match.
Both are carried through; the format-specific bits are confined to `load`,
`add_water`, and `write_frame`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from ase import Atom, Atoms
from ase.constraints import FixAtoms
from ase.geometry import find_mic, get_distances
from ase.io import read as ase_read
from ase.io import write as ase_write

LAMMPS = "lammps"
POSCAR = "poscar"
FORMATS = (LAMMPS, POSCAR)

DESCRIPTION = {LAMMPS: "LAMMPS data file", POSCAR: "VASP POSCAR"}
# The knob for overriding the atom style; None in a format with no equivalent,
# which is how a caller knows to reject that flag.
STYLE_OPTION = {LAMMPS: "--atom-style", POSCAR: None}
# Whether the format can carry per-atom charges at all.
SUPPORTS_CHARGES = {LAMMPS: True, POSCAR: False}
# Whether the format can express per-atom freezing in the structure file itself.
SUPPORTS_FREEZING = {LAMMPS: False, POSCAR: True}

# Atom styles this module knows how to hand to ASE. Others are refused up front
# rather than written wrongly.
SUPPORTED_STYLES = ("atomic", "charge", "full")

# LAMMPS sections that describe bonded topology. Adding atoms to a file that has
# any of them would leave them inconsistent, and ASE's reader drops them without
# a word, so they have to be caught by reading the file directly.
TOPOLOGY_SECTIONS = ("Bonds", "Angles", "Dihedrals", "Impropers")

ELEMENT_WORDS = {"H": "hydrogen", "O": "oxygen", "Si": "silicon"}


@dataclass
class Frame:
    """
    A structure held open for editing, plus the per-format settings ASE does not
    keep on the Atoms object itself.
    """
    atoms: Atoms
    fmt: str
    n_input: int                       # atoms that came from the input file
    atom_style: str | None = None      # lammps
    specorder: list[str] = field(default_factory=list)   # lammps: type i+1 is specorder[i]
    direct: bool = True                # poscar
    source: str = ""
    dropped_velocities: bool = False   # the input carried velocities; see _drop_velocities


# -- format detection ---------------------------------------------------------

def detect_atom_style(input_file: str) -> str | None:
    """The style named in the "Atoms # <style>" comment, or None if absent."""
    with open(input_file) as handle:
        for line in handle:
            keyword, _, comment = line.partition("#")
            if keyword.strip() != "Atoms":
                continue
            words = comment.split()
            return words[0] if words else None
    raise ValueError(f"{input_file}: no Atoms section found")


def looks_like(input_file: str, fmt: str) -> bool:
    """
    Cheap check that a file is the format it is claimed to be, so the wrong
    --format is caught with a clear message instead of a parser traceback.
    """
    try:
        if fmt == LAMMPS:
            with open(input_file) as handle:
                return any(line.partition("#")[0].strip() == "Atoms" for line in handle)
        # A POSCAR's second line is a lone scale factor and the next three are a
        # 3x3 lattice.
        with open(input_file) as handle:
            lines = [handle.readline() for _ in range(5)]
        if len(lines[1].split()) != 1:
            return False
        float(lines[1].split()[0])
        for row in lines[2:5]:
            if len(row.split()) != 3:
                return False
            [float(value) for value in row.split()]
        return True
    except (OSError, ValueError, IndexError):
        return False


def find_topology(input_file: str, fmt: str) -> list[str]:
    """
    Topology sections present in the file. ASE's LAMMPS reader ignores them
    silently, so a caller that must refuse them has to be told here. A POSCAR has
    no such thing.
    """
    if fmt != LAMMPS:
        return []
    found = []
    with open(input_file) as handle:
        for line in handle:
            keyword = line.partition("#")[0].strip()
            if keyword in TOPOLOGY_SECTIONS and keyword not in found:
                found.append(keyword)
    return sorted(found)


# -- reading ------------------------------------------------------------------

def _derive_specorder(atoms: Atoms, input_file: str) -> list[str]:
    """
    The element of each LAMMPS atom type, indexed so specorder[i] is type i+1.

    Taken from the file's own type ids rather than re-derived, so writing the
    frame back preserves the numbering the potential's pair_coeff line depends
    on. ASE infers the elements themselves from the Masses table.
    """
    types = atoms.arrays["type"]
    symbols = np.array(atoms.get_chemical_symbols())
    mapping = {int(t): str(symbols[types == t][0]) for t in sorted(set(types.tolist()))}
    highest = max(mapping)
    missing = [t for t in range(1, highest + 1) if t not in mapping]
    if missing:
        raise ValueError(
            f"{input_file}: atom type(s) {missing} are declared but have no atoms, "
            f"so their element cannot be determined; this script needs every type "
            f"in the file to be in use"
        )
    return [mapping[t] for t in range(1, highest + 1)]


def load(input_file: str, fmt: str, atom_style: str | None = None,
         quiet: bool = False) -> Frame:
    """Read a structure. `atom_style` applies to LAMMPS only."""
    if fmt == LAMMPS:
        declared = detect_atom_style(input_file)
        style = atom_style or declared
        if style is None:
            raise ValueError(
                f"{input_file}: the Atoms section has no '# <style>' comment, so "
                f"the atom style cannot be detected; set --atom-style"
            )
        if atom_style is not None and declared is not None and declared != atom_style:
            raise ValueError(
                f"{input_file}: --atom-style is '{atom_style}' but the Atoms section "
                f"says '{declared}'; fix whichever is wrong (leave it unset to trust "
                f"the file)"
            )
        if style not in SUPPORTED_STYLES:
            raise ValueError(
                f"{input_file}: atom style '{style}' is not supported "
                f"(supported: {', '.join(SUPPORTED_STYLES)})"
            )
        # Image flags are deliberately not read; see the module docstring.
        atoms = ase_read(input_file, format="lammps-data", atom_style=style,
                         read_image_flags=False, sort_by_id=True)
        dropped = _drop_velocities(atoms)
        specorder = _derive_specorder(atoms, input_file)
        if not quiet:
            print(f"atom style: {style}"
                  f"{' (from the Atoms section comment)' if atom_style is None else ' (from --atom-style)'}"
                  f"; atom types " + ", ".join(f"{i + 1}={s}" for i, s in enumerate(specorder)))
        return Frame(atoms=atoms, fmt=fmt, n_input=len(atoms), atom_style=style,
                     specorder=specorder, source=input_file,
                     dropped_velocities=dropped)

    atoms = ase_read(input_file, format="vasp")
    frame = Frame(atoms=atoms, fmt=fmt, n_input=len(atoms),
                  direct=_coordinate_mode_is_direct(input_file), source=input_file,
                  dropped_velocities=_drop_velocities(atoms))
    if not quiet:
        frozen = _frozen_indices(atoms)
        extras = []
        if frozen:
            extras.append(f"{len(frozen)} atoms frozen by selective dynamics")
        if frame.dropped_velocities:
            extras.append("velocity block dropped")
        print(f"format: POSCAR ({len(atoms)} sites, {' '.join(species_order(frame))}"
              + (", " + ", ".join(extras) if extras else "") + ")")
    return frame


def _coordinate_mode(input_file: str) -> tuple[int, bool]:
    """
    The line index of the Direct/Cartesian keyword and which one it is. ASE does
    not record the mode it read, and flipping a file's mode on rewrite is a
    gratuitous diff.
    """
    with open(input_file) as handle:
        lines = [line.strip() for line in handle.read().splitlines()]
    for index, line in enumerate(lines[5:10], start=5):
        initial = line[:1].upper()
        if initial == "S":          # Selective dynamics
            continue
        if initial == "D":
            return index, True
        if initial in ("C", "K"):   # Cartesian
            return index, False
    raise ValueError(f"{input_file}: no Direct/Cartesian line found")


def _coordinate_mode_is_direct(input_file: str) -> bool:
    return _coordinate_mode(input_file)[1]


def _drop_velocities(atoms: Atoms) -> bool:
    """
    Discard any velocities the input carried. Returns whether there were any.

    Every frame is a fresh starting point for its own minimization, so a velocity
    inherited from whatever run produced the input is meaningless in it - and
    carrying velocities through ASE is treacherous anyway, because each rebuild
    of an Atoms re-applies its constraints and FixAtoms zeroes the momenta of the
    atoms it fixes. Dropping the block outright is both simpler and more honest
    than writing a half-preserved one.
    """
    if not atoms.has("momenta"):
        return False
    del atoms.arrays["momenta"]
    return True


# -- inspection ---------------------------------------------------------------

def atom_count(frame: Frame) -> int:
    return len(frame.atoms)


def species_order(frame: Frame) -> list[str]:
    """The species line's order: each element in the order it first appears."""
    order: list[str] = []
    for symbol in frame.atoms.get_chemical_symbols():
        if symbol not in order:
            order.append(symbol)
    return order


def is_triclinic(frame: Frame) -> bool:
    return not all(abs(angle - 90.0) < 1e-6 for angle in frame.atoms.cell.angles())


def species_label(frame: Frame, atom_id: int):
    """
    How an atom is named in reports and the manifest's *_type columns: the
    integer atom type for LAMMPS, the element symbol for a POSCAR.
    """
    symbol = frame.atoms[atom_id - 1].symbol
    if frame.fmt == LAMMPS:
        return frame.specorder.index(symbol) + 1
    return symbol


def element_of(frame: Frame, atom_id: int) -> str:
    return frame.atoms[atom_id - 1].symbol


def check_target(frame: Frame, atom_id: int, element: str) -> str | None:
    """A note if `atom_id` is not the element it is supposed to be, else None."""
    if not 1 <= atom_id <= len(frame.atoms):
        raise ValueError(
            f"atom id {atom_id} is not in {frame.source} (it has "
            f"{len(frame.atoms)} atoms, numbered from 1)"
        )
    actual = element_of(frame, atom_id)
    if actual == element:
        return None
    word = ELEMENT_WORDS.get(element, element)
    if frame.fmt == LAMMPS:
        return (f"atom {atom_id} has type {species_label(frame, atom_id)} ({actual}), "
                f"which is not {word}")
    return f"atom {atom_id} is {actual}, which is not {word}"


def identify_water(frame: Frame) -> tuple[int, list[int]]:
    """
    Pick the oxygen and the two hydrogens out of a three-atom structure, for
    reading O-H offsets off a water template.
    """
    symbols = [site.symbol for site in frame.atoms]
    oxygens = [i + 1 for i, s in enumerate(symbols) if s == "O"]
    hydrogens = [i + 1 for i, s in enumerate(symbols) if s == "H"]
    if len(oxygens) != 1 or len(hydrogens) != 2:
        raise ValueError(
            f"expected one O and two H in the water template, found "
            f"{' '.join(symbols)}"
        )
    return oxygens[0], hydrogens


def next_molecule_id(frame: Frame) -> int:
    if "mol-id" not in frame.atoms.arrays:
        return 0
    return int(frame.atoms.arrays["mol-id"].max()) + 1


def frozen_count(frame: Frame) -> int:
    """How many atoms selective dynamics holds fixed. Always 0 for LAMMPS."""
    return len(_frozen_indices(frame.atoms))


def _frozen_indices(atoms: Atoms) -> set[int]:
    frozen: set[int] = set()
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            frozen.update(int(i) for i in constraint.index)
    return frozen


# -- geometry -----------------------------------------------------------------
#
# All of this goes through the cell, so it is correct for a skewed lattice as
# well as an orthogonal box - in both formats.

def position_of(frame: Frame, atom_id: int) -> np.ndarray:
    if not 1 <= atom_id <= len(frame.atoms):
        raise ValueError(
            f"atom id {atom_id} is not in {frame.source} (it has "
            f"{len(frame.atoms)} atoms, numbered from 1)"
        )
    return np.asarray(frame.atoms.positions[atom_id - 1], dtype=float)


def displacement(frame: Frame, start, end) -> np.ndarray:
    """The shortest vector from `start` to `end` under PBC, for any cell."""
    delta = np.asarray(end, dtype=float) - np.asarray(start, dtype=float)
    vectors, _ = find_mic(np.atleast_2d(delta), frame.atoms.cell, pbc=True)
    return np.asarray(vectors[0], dtype=float)


def distance(frame: Frame, a, b) -> float:
    return float(np.linalg.norm(displacement(frame, a, b)))


def wrap_position(frame: Frame, position: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Wrap a position into the cell. Returns the wrapped position and the image it
    came from, so a caller can count how many atoms crossed a boundary.
    """
    cell = frame.atoms.cell
    frac = cell.scaled_positions(np.asarray(position, dtype=float).reshape(1, 3))[0]
    image = np.floor(frac).astype(int)
    return cell.cartesian_positions((frac - image).reshape(1, 3))[0], image


def nearest_existing(frame: Frame, position: np.ndarray, count: int = 1,
                     exclude: set[int] | None = None) -> list[tuple[int, object, float]]:
    """
    The `count` nearest atoms to a position under PBC, as a list of
    (atom_id, species label, distance), closest first. `exclude` drops atom ids
    from the search - an atom being deliberately walked toward a target is not a
    clash against that target.
    """
    atoms = frame.atoms
    if len(atoms) == 0:
        return []
    _, lengths = get_distances(np.asarray(position, dtype=float).reshape(1, 3),
                               atoms.positions, cell=atoms.cell, pbc=True)
    distances = np.asarray(lengths).ravel()
    ids = np.arange(1, len(atoms) + 1)
    if exclude:
        keep = ~np.isin(ids, list(exclude))
        ids, distances = ids[keep], distances[keep]
    if not len(ids):
        return []
    order = np.argsort(distances)[:count]
    return [(int(ids[i]), species_label(frame, int(ids[i])), float(distances[i]))
            for i in order]


# -- adding the water ---------------------------------------------------------

def _row_plan(frame: Frame, added: list[tuple[int, str]]) -> list[int]:
    """
    The order the atoms go in the output, as a permutation of indices into the
    already-extended Atoms.

    A POSCAR groups its sites by species and states the counts in a header, so a
    new atom has to join the end of its own species' run to keep that grouping -
    which is what renumbers everything after it. The input's species order is
    preserved rather than re-sorted, because the POTCAR has to match it. A
    species the structure does not have yet becomes a new trailing group, and the
    POTCAR has to grow to match.

    A LAMMPS file needs none of this: its types are per-atom, so new atoms simply
    go on the end.
    """
    if frame.fmt != POSCAR:
        return list(range(frame.n_input + len(added)))

    symbols = frame.atoms.get_chemical_symbols()
    last_index: dict[str, int] = {}
    for index in range(frame.n_input):
        last_index[symbols[index]] = index

    groups: dict[str, list[int]] = {}
    for index, species in added:
        groups.setdefault(species, []).append(index)

    plan: list[int] = []
    for index in range(frame.n_input):
        plan.append(index)
        for species, members in groups.items():
            if last_index.get(species) == index:
                plan.extend(members)
    for species, members in groups.items():
        if species not in last_index:
            plan.extend(members)
    return plan


def new_species(frame: Frame, spec: dict, labels) -> list[str]:
    """
    Species the water introduces that the input does not already contain. These
    become new groups on a POSCAR's species line (so the POTCAR must be extended)
    or a new LAMMPS atom type (so the pair_coeff line must be).
    """
    present = set(frame.atoms.get_chemical_symbols()[:frame.n_input])
    wanted: list[str] = []
    for label in labels:
        species = spec["species"][label]
        if species not in present and species not in wanted:
            wanted.append(species)
    return wanted


def _added_atoms(frame: Frame, spec: dict, labels) -> list[tuple[int, str]]:
    """Where the water atoms sit before the reindex, and what species each is."""
    return [(frame.n_input + offset, spec["species"][label])
            for offset, label in enumerate(labels)]


def index_map(frame: Frame, spec: dict, labels) -> dict[int, int]:
    """
    Input atom id -> id in the frames that will be written. Only the layout of
    the input decides this, so it is the same for every frame and can be worked
    out before any of them are built. The identity for LAMMPS, which appends.
    """
    plan = _row_plan(frame, _added_atoms(frame, spec, labels))
    where = {old: new for new, old in enumerate(plan)}
    return {old + 1: where[old] + 1 for old in range(frame.n_input)}


def add_water(frame: Frame, positions: dict, labels, spec: dict,
              freeze_substrate: bool = False) -> tuple[Frame, dict, dict, int]:
    """
    Copy the frame and add the water atoms, wrapping each into the cell.
    Returns (new frame, {label: atom id}, {input id: output id}, how many wrapped).

    The ids are indices into the frame as it will be written, which is what a
    caller needs to name the new atoms in a constraint.
    """
    atoms = frame.atoms.copy()
    wrapped_count = 0
    added: list[tuple[int, str]] = []

    for label in labels:
        wrapped, image = wrap_position(frame, positions[label])
        if image.any():
            wrapped_count += 1
        species = spec["species"][label]
        atom = Atom(species, wrapped)
        if frame.fmt == LAMMPS and "initial_charges" in atoms.arrays:
            atom.charge = float(spec["charges"][label])
        added.append((len(atoms), species))
        atoms.append(atom)

    if frame.fmt == LAMMPS and "mol-id" in atoms.arrays:
        atoms.arrays["mol-id"][frame.n_input:] = int(spec["molecule_id"])

    # Applied before the reindex so ASE remaps it along with everything else.
    if freeze_substrate:
        atoms.set_constraint(FixAtoms(indices=list(range(frame.n_input))))

    plan = _row_plan(frame, added)
    atoms = atoms[plan]

    # ASE keeps the type ids it read, and atoms appended above have none, so the
    # writer would emit them as type 0 - an invalid LAMMPS type in a file whose
    # Masses table still looks right. Dropping the array makes the writer assign
    # types from `specorder` instead, which is derived from the input.
    if "type" in atoms.arrays:
        del atoms.arrays["type"]

    where = {old: new for new, old in enumerate(plan)}
    ids = {label: where[index] + 1 for (index, _), label in zip(added, labels)}
    index_map = {old + 1: where[old] + 1 for old in range(frame.n_input)}

    result = Frame(**{**frame.__dict__, "atoms": atoms})
    # The caller settles the type->element list, because which type a new element
    # gets is a question about the potential's pair_coeff line, not about the
    # file. It must cover every element present.
    result.specorder = list(spec.get("specorder") or frame.specorder)
    if frame.fmt == LAMMPS:
        missing = sorted(set(atoms.get_chemical_symbols()) - set(result.specorder))
        if missing:
            raise ValueError(
                f"no atom type assigned to {', '.join(missing)}; the specorder "
                f"passed in does not cover every element in the frame"
            )
    return result, ids, index_map, wrapped_count


# -- writing ------------------------------------------------------------------

def frame_path(outdir, index: int, fmt: str) -> Path:
    """
    Where frame `index` is written, relative to `outdir`. VASP reads a file named
    literally POSCAR, so each frame gets its own directory and is ready to run
    once an INCAR, KPOINTS, and POTCAR are dropped in beside it.
    """
    if fmt == POSCAR:
        return Path(outdir) / str(index) / "POSCAR"
    return Path(outdir) / f"{index}.data"


def write_frame(frame: Frame, path, header_comment: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if frame.fmt == LAMMPS:
        ase_write(str(path), frame.atoms, format="lammps-data",
                  atom_style=frame.atom_style, specorder=frame.specorder,
                  masses=True, velocities=False, write_image_flags=False)
        # ASE's first line is a fixed banner; replace it with something that says
        # where the frame came from, as LAMMPS treats line 1 as a comment.
        lines = path.read_text().splitlines(keepends=True)
        lines[0] = f"{header_comment}\n"
        path.write_text("".join(lines))
        return

    ase_write(str(path), frame.atoms, format="vasp", direct=frame.direct,
              sort=False, vasp5=True)
    # ASE writes the species line as the comment; make it say where this came
    # from instead. VASP only reads line 1 as a label.
    lines = path.read_text().splitlines(keepends=True)
    lines[0] = f"{header_comment}\n"
    path.write_text("".join(lines))


def restraint_advice(fmt: str) -> str:
    """How to tell the user to relax a frame in this format."""
    if fmt == POSCAR:
        return ("relax each one with the O-Si distance held by an ICONST constraint "
                "and the hydrogens free before reading any energy off it")
    return ("relax each one with the O-Si distance restrained and the hydrogens "
            "free before reading any energy off it")
