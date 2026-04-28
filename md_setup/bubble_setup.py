#!/usr/bin/env python3
"""
bubble_setup.py

Create spherical bubble(s) in a LAMMPS data file (water medium) and
optionally fill with atoms carved from a second LAMMPS data file.

Both data files are assumed to share the same coordinate system.

Usage: edit the CONFIG section below and run as a script,
       or import create_bubble() directly. By default the output path is beside
       the input file, named ``<stem>_bubble<suffix>``.
"""

import numpy as np
from collections import defaultdict
from pathlib import Path

# ─── CONFIG ──────────────────────────────────────────────────────────────────
INPUT_FILE = "water.data"
# None = same directory as INPUT_FILE, filename stem + "_bubble" + extension (e.g. water.data → water_bubble.data)
OUTPUT_FILE: str | None = None


def default_bubble_output_path(input_file: str) -> str:
    """Output path beside input: ``<stem>_bubble<suffix>``."""
    p = Path(input_file)
    return str(p.with_name(f"{p.stem}_bubble{p.suffix}"))

# (x, y, z) bubble center positions in Angstroms
BUBBLE_CENTERS: list[tuple[float, float, float]] = [
    (0.0, 0.0, 0.0),
]

BUBBLE_RADIUS = 10.0  # Angstroms

# Atoms per water molecule — used for broken-molecule detection
WATER_MOL_SIZE = 3

# Path to filler LAMMPS data file, or None to leave bubble empty
FILLER_FILE = None  # e.g. "ethanol.data"

# Atoms per filler molecule; None = auto-detect from most common mol size
FILLER_MOL_SIZE = None

# Silica atom types (list of ints) for slab detection.
# None = auto-detect: any atom not belonging to a complete water molecule is treated as silica.
SILICA_ATOM_TYPES: list[int] | None = None
# ─────────────────────────────────────────────────────────────────────────────


# ── parsing ───────────────────────────────────────────────────────────────────

def parse_lammps_data(filename: str) -> dict:
    """
    Parse a LAMMPS data file.

    Returns a dict with:
        header, box, masses, atom_style,
        atoms (list of dicts), bonds, angles, dihedrals, impropers,
        n_atom_types, n_bond_types, n_angle_types, n_dihedral_types, n_improper_types
    """
    data = {
        'header': '',
        'box': {},
        'masses': {},
        'atom_style': 'full',
        'atoms': [],
        'bonds': [],
        'angles': [],
        'dihedrals': [],
        'impropers': [],
        'n_atoms': 0,
        'n_bonds': 0,
        'n_angles': 0,
        'n_dihedrals': 0,
        'n_impropers': 0,
        'n_atom_types': 0,
        'n_bond_types': 0,
        'n_angle_types': 0,
        'n_dihedral_types': 0,
        'n_improper_types': 0,
    }

    with open(filename) as f:
        lines = f.readlines()

    data['header'] = lines[0].strip() if lines else 'LAMMPS data file'

    section = None
    for line in lines[1:]:
        stripped = line.split('#')[0].strip()
        if not stripped:
            continue

        # ── section headers ──
        if stripped.startswith('Masses'):
            section = 'Masses'; continue
        if stripped.startswith('Atoms'):
            comment = line.split('#')
            if len(comment) > 1:
                data['atom_style'] = comment[1].strip()
            section = 'Atoms'; continue
        if stripped.startswith('Bonds'):
            section = 'Bonds'; continue
        if stripped.startswith('Angles'):
            section = 'Angles'; continue
        if stripped.startswith('Dihedrals'):
            section = 'Dihedrals'; continue
        if stripped.startswith('Impropers'):
            section = 'Impropers'; continue
        if stripped[0].isalpha():
            section = 'unknown'; continue

        parts = stripped.split()
        if not parts:
            continue

        # ── header block (section is None) ──
        if section is None:
            _parse_header_line(stripped, data)
            continue

        if section == 'Masses':
            if len(parts) >= 2:
                data['masses'][int(parts[0])] = float(parts[1])

        elif section == 'Atoms':
            atom = _parse_atom_line(parts, data['atom_style'])
            if atom:
                data['atoms'].append(atom)

        elif section == 'Bonds':
            if len(parts) >= 4:
                data['bonds'].append({
                    'id': int(parts[0]), 'type': int(parts[1]),
                    'atom1': int(parts[2]), 'atom2': int(parts[3]),
                })

        elif section == 'Angles':
            if len(parts) >= 5:
                data['angles'].append({
                    'id': int(parts[0]), 'type': int(parts[1]),
                    'atom1': int(parts[2]), 'atom2': int(parts[3]), 'atom3': int(parts[4]),
                })

        elif section == 'Dihedrals':
            if len(parts) >= 6:
                data['dihedrals'].append({
                    'id': int(parts[0]), 'type': int(parts[1]),
                    'atom1': int(parts[2]), 'atom2': int(parts[3]),
                    'atom3': int(parts[4]), 'atom4': int(parts[5]),
                })

        elif section == 'Impropers':
            if len(parts) >= 6:
                data['impropers'].append({
                    'id': int(parts[0]), 'type': int(parts[1]),
                    'atom1': int(parts[2]), 'atom2': int(parts[3]),
                    'atom3': int(parts[4]), 'atom4': int(parts[5]),
                })

    return data


def _parse_header_line(line: str, data: dict) -> None:
    if   'atom types'     in line: data['n_atom_types']     = int(line.split()[0])
    elif 'bond types'     in line: data['n_bond_types']     = int(line.split()[0])
    elif 'angle types'    in line: data['n_angle_types']    = int(line.split()[0])
    elif 'dihedral types' in line: data['n_dihedral_types'] = int(line.split()[0])
    elif 'improper types' in line: data['n_improper_types'] = int(line.split()[0])
    elif 'atoms'          in line: data['n_atoms']          = int(line.split()[0])
    elif 'bonds'          in line: data['n_bonds']          = int(line.split()[0])
    elif 'angles'         in line: data['n_angles']         = int(line.split()[0])
    elif 'dihedrals'      in line: data['n_dihedrals']      = int(line.split()[0])
    elif 'impropers'      in line: data['n_impropers']      = int(line.split()[0])
    elif 'xlo xhi'        in line:
        p = line.split(); data['box']['xlo'] = float(p[0]); data['box']['xhi'] = float(p[1])
    elif 'ylo yhi'        in line:
        p = line.split(); data['box']['ylo'] = float(p[0]); data['box']['yhi'] = float(p[1])
    elif 'zlo zhi'        in line:
        p = line.split(); data['box']['zlo'] = float(p[0]); data['box']['zhi'] = float(p[1])
    elif 'xy xz yz'       in line:
        p = line.split()
        data['box']['xy'] = float(p[0]); data['box']['xz'] = float(p[1]); data['box']['yz'] = float(p[2])


def _parse_atom_line(parts: list[str], style: str) -> dict | None:
    try:
        if style in ('full',):
            # atom_id mol_id type charge x y z [ix iy iz]
            a = dict(id=int(parts[0]), mol=int(parts[1]), type=int(parts[2]),
                     charge=float(parts[3]), x=float(parts[4]), y=float(parts[5]), z=float(parts[6]))
            if len(parts) > 9:
                a['ix'] = int(parts[7]); a['iy'] = int(parts[8]); a['iz'] = int(parts[9])
        elif style == 'molecular':
            # atom_id mol_id type x y z [ix iy iz]
            a = dict(id=int(parts[0]), mol=int(parts[1]), type=int(parts[2]),
                     charge=0.0, x=float(parts[3]), y=float(parts[4]), z=float(parts[5]))
            if len(parts) > 8:
                a['ix'] = int(parts[6]); a['iy'] = int(parts[7]); a['iz'] = int(parts[8])
        elif style == 'charge':
            # atom_id type charge x y z
            a = dict(id=int(parts[0]), mol=0, type=int(parts[1]),
                     charge=float(parts[2]), x=float(parts[3]), y=float(parts[4]), z=float(parts[5]))
        elif style == 'atomic':
            # atom_id type x y z
            a = dict(id=int(parts[0]), mol=0, type=int(parts[1]),
                     charge=0.0, x=float(parts[2]), y=float(parts[3]), z=float(parts[4]))
        else:
            # fall back to full
            a = dict(id=int(parts[0]), mol=int(parts[1]), type=int(parts[2]),
                     charge=float(parts[3]), x=float(parts[4]), y=float(parts[5]), z=float(parts[6]))
        return a
    except (ValueError, IndexError):
        return None


# ── geometry ──────────────────────────────────────────────────────────────────

def atoms_in_bubble(atoms: list[dict],
                    centers: list[tuple],
                    radius: float) -> set[int]:
    """Return atom IDs within radius of any center (Euclidean distance)."""
    if not atoms or not centers:
        return set()
    pos = np.array([[a['x'], a['y'], a['z']] for a in atoms])
    ids = np.array([a['id'] for a in atoms])
    inside = np.zeros(len(atoms), dtype=bool)
    for c in centers:
        inside |= np.linalg.norm(pos - np.array(c), axis=1) <= radius
    return set(ids[inside].tolist())


# ── silica detection & validation ────────────────────────────────────────────

def identify_silica_atoms(atoms: list[dict],
                          water_mol_size: int,
                          silica_atom_types: list[int] | None = None) -> list[dict]:
    """
    Return atoms belonging to the silica slab.

    If silica_atom_types is given, select by type.
    Otherwise auto-detect: atoms whose molecule has a size != water_mol_size
    (mol_id 0 atoms are always included as silica).
    """
    if silica_atom_types is not None:
        return [a for a in atoms if a['type'] in silica_atom_types]

    mol_sizes: dict[int, int] = defaultdict(int)
    for a in atoms:
        mol_sizes[a['mol']] += 1

    return [a for a in atoms if mol_sizes[a['mol']] != water_mol_size or a['mol'] == 0]


def validate_bubble_placement(centers: list[tuple],
                               radius: float,
                               silica_atoms: list[dict]) -> None:
    """
    Raise ValueError if any bubble center is inside the silica slab or if the
    bubble sphere would overlap it.

    The check is per-atom: the minimum distance from each center to any silica
    atom must be strictly greater than radius.
    """
    if not silica_atoms:
        return

    sil_pos = np.array([[a['x'], a['y'], a['z']] for a in silica_atoms])

    for i, c in enumerate(centers):
        dists = np.linalg.norm(sil_pos - np.array(c), axis=1)
        min_dist = float(dists.min())
        nearest = silica_atoms[int(dists.argmin())]

        if min_dist <= 0:
            raise ValueError(
                f"Bubble center {i} {c} is inside the silica slab "
                f"(nearest silica atom id={nearest['id']} at distance {min_dist:.3f} Å)."
            )
        if min_dist <= radius:
            raise ValueError(
                f"Bubble {i} (center={c}, r={radius} Å) overlaps the silica slab: "
                f"nearest silica atom id={nearest['id']} is only {min_dist:.3f} Å away. "
                f"Reduce the radius or move the center."
            )


# ── molecule filtering ────────────────────────────────────────────────────────

def _most_common_mol_size(atoms: list[dict]) -> int:
    sizes = defaultdict(int)
    for a in atoms:
        sizes[a['mol']] += 1
    counts = list(sizes.values())
    return max(set(counts), key=counts.count)


def remove_broken_molecules(atoms: list[dict], mol_size: int | None = None) -> list[dict]:
    """
    Keep only complete molecules. A molecule is complete if its atom count
    equals mol_size (or the most common size if mol_size is None).
    """
    if not atoms:
        return atoms
    expected = mol_size if mol_size is not None else _most_common_mol_size(atoms)
    mol_atoms: dict[int, list] = defaultdict(list)
    for a in atoms:
        mol_atoms[a['mol']].append(a)
    return [a for mol in mol_atoms.values() if len(mol) == expected for a in mol]


def _filter_bonds(bonds, valid_ids):
    return [b for b in bonds if b['atom1'] in valid_ids and b['atom2'] in valid_ids]

def _filter_angles(angles, valid_ids):
    return [a for a in angles
            if a['atom1'] in valid_ids and a['atom2'] in valid_ids and a['atom3'] in valid_ids]

def _filter_dihedrals(dihedrals, valid_ids):
    return [d for d in dihedrals
            if d['atom1'] in valid_ids and d['atom2'] in valid_ids
            and d['atom3'] in valid_ids and d['atom4'] in valid_ids]

def _filter_impropers(impropers, valid_ids):
    return _filter_dihedrals(impropers, valid_ids)


# ── renumbering ───────────────────────────────────────────────────────────────

def renumber(data: dict) -> dict:
    """Renumber all atom IDs, mol IDs, and topology IDs sequentially from 1."""
    atoms = data['atoms']

    # atom ID map
    id_map = {a['id']: new for new, a in enumerate(atoms, 1)}
    for new, a in enumerate(atoms, 1):
        a['id'] = new

    # mol ID map (preserve grouping, assign in first-seen order)
    seen_mols: dict[int, int] = {}
    mol_counter = 1
    mol_map: dict[int, int] = {}
    for a in atoms:
        old_mol = a['mol']
        if old_mol not in mol_map:
            mol_map[old_mol] = mol_counter
            mol_counter += 1

    for a in atoms:
        a['mol'] = mol_map[a['mol']]

    def remap_topo(items, keys):
        out = []
        for item in items:
            if all(item[k] in id_map for k in keys):
                for k in keys:
                    item[k] = id_map[item[k]]
                out.append(item)
        for new_id, item in enumerate(out, 1):
            item['id'] = new_id
        return out

    data['bonds']     = remap_topo(data['bonds'],     ['atom1', 'atom2'])
    data['angles']    = remap_topo(data['angles'],     ['atom1', 'atom2', 'atom3'])
    data['dihedrals'] = remap_topo(data['dihedrals'],  ['atom1', 'atom2', 'atom3', 'atom4'])
    data['impropers'] = remap_topo(data['impropers'],  ['atom1', 'atom2', 'atom3', 'atom4'])

    return data


# ── filler merge ──────────────────────────────────────────────────────────────

def merge_filler(base: dict, filler: dict,
                 centers: list[tuple], radius: float,
                 filler_mol_size: int | None = None) -> dict:
    """
    Carve a sphere of filler atoms from filler data, remove partial molecules
    (distance-based), and merge into base system with offset IDs and types.
    """
    # 1. atoms within bubble in filler coordinate system
    in_bubble = atoms_in_bubble(filler['atoms'], centers, radius)
    carved = [a for a in filler['atoms'] if a['id'] in in_bubble]

    # 2. remove partial filler molecules
    expected_size = filler_mol_size if filler_mol_size is not None \
                    else _most_common_mol_size(filler['atoms'])
    carved = remove_broken_molecules(carved, expected_size)
    if not carved:
        print("  Warning: no complete filler molecules found within bubble radius.")
        return base

    valid_filler = {a['id'] for a in carved}
    f_bonds      = _filter_bonds(filler['bonds'], valid_filler)
    f_angles     = _filter_angles(filler['angles'], valid_filler)
    f_dihedrals  = _filter_dihedrals(filler['dihedrals'], valid_filler)

    # 3. offsets
    atom_id_off  = max((a['id']  for a in base['atoms']), default=0)
    mol_id_off   = max((a['mol'] for a in base['atoms']), default=0)
    bond_id_off  = max((b['id']  for b in base['bonds']), default=0)
    angle_id_off = max((a['id']  for a in base['angles']), default=0)
    dih_id_off   = max((d['id']  for d in base['dihedrals']), default=0)
    type_off     = base['n_atom_types']
    btype_off    = base['n_bond_types']
    atype_off    = base['n_angle_types']
    dtype_off    = base['n_dihedral_types']

    # 4. build filler atom-ID → new-ID and mol → new-mol maps
    unique_mols = list(dict.fromkeys(a['mol'] for a in carved))
    fmol_map = {old: mol_id_off + new for new, old in enumerate(unique_mols, 1)}
    fid_map  = {}
    for new_idx, a in enumerate(carved, 1):
        fid_map[a['id']] = atom_id_off + new_idx

    # 5. append atoms
    for a in carved:
        base['atoms'].append({**a,
            'id':   fid_map[a['id']],
            'mol':  fmol_map[a['mol']],
            'type': a['type'] + type_off,
        })

    # 6. append bonds
    for i, b in enumerate(f_bonds, bond_id_off + 1):
        if b['atom1'] in fid_map and b['atom2'] in fid_map:
            base['bonds'].append({'id': i, 'type': b['type'] + btype_off,
                                  'atom1': fid_map[b['atom1']], 'atom2': fid_map[b['atom2']]})

    # 7. append angles
    for i, a in enumerate(f_angles, angle_id_off + 1):
        if a['atom1'] in fid_map and a['atom2'] in fid_map and a['atom3'] in fid_map:
            base['angles'].append({'id': i, 'type': a['type'] + atype_off,
                                   'atom1': fid_map[a['atom1']], 'atom2': fid_map[a['atom2']],
                                   'atom3': fid_map[a['atom3']]})

    # 8. append dihedrals
    for i, d in enumerate(f_dihedrals, dih_id_off + 1):
        if all(d[k] in fid_map for k in ['atom1', 'atom2', 'atom3', 'atom4']):
            base['dihedrals'].append({'id': i, 'type': d['type'] + dtype_off,
                                      'atom1': fid_map[d['atom1']], 'atom2': fid_map[d['atom2']],
                                      'atom3': fid_map[d['atom3']], 'atom4': fid_map[d['atom4']]})

    # 9. update type counts and masses
    base['n_atom_types']     += filler['n_atom_types']
    base['n_bond_types']     += filler['n_bond_types']
    base['n_angle_types']    += filler['n_angle_types']
    base['n_dihedral_types'] += filler['n_dihedral_types']
    for ftype, mass in filler['masses'].items():
        base['masses'][ftype + type_off] = mass

    print(f"  Merged {len(carved)} filler atoms ({len(unique_mols)} molecules)")
    return base


# ── writing ───────────────────────────────────────────────────────────────────

def write_lammps_data(data: dict, filename: str) -> None:
    style = data.get('atom_style', 'full')
    lines = [data['header'], '']

    lines.append(f"{len(data['atoms'])} atoms")
    if data['bonds']:     lines.append(f"{len(data['bonds'])} bonds")
    if data['angles']:    lines.append(f"{len(data['angles'])} angles")
    if data['dihedrals']: lines.append(f"{len(data['dihedrals'])} dihedrals")
    if data['impropers']: lines.append(f"{len(data['impropers'])} impropers")
    lines.append('')

    lines.append(f"{data['n_atom_types']} atom types")
    if data['n_bond_types']:     lines.append(f"{data['n_bond_types']} bond types")
    if data['n_angle_types']:    lines.append(f"{data['n_angle_types']} angle types")
    if data['n_dihedral_types']: lines.append(f"{data['n_dihedral_types']} dihedral types")
    if data['n_improper_types']: lines.append(f"{data['n_improper_types']} improper types")
    lines.append('')

    box = data['box']
    lines.append(f"{box['xlo']:.6f} {box['xhi']:.6f} xlo xhi")
    lines.append(f"{box['ylo']:.6f} {box['yhi']:.6f} ylo yhi")
    lines.append(f"{box['zlo']:.6f} {box['zhi']:.6f} zlo zhi")
    if 'xy' in box:
        lines.append(f"{box['xy']:.6f} {box['xz']:.6f} {box['yz']:.6f} xy xz yz")
    lines.append('')

    lines.append('Masses'); lines.append('')
    for t in sorted(data['masses']):
        lines.append(f"{t} {data['masses'][t]}")
    lines.append('')

    lines.append(f'Atoms # {style}'); lines.append('')
    for a in sorted(data['atoms'], key=lambda x: x['id']):
        ix, iy, iz = a.get('ix', 0), a.get('iy', 0), a.get('iz', 0)
        if style == 'full':
            lines.append(f"{a['id']} {a['mol']} {a['type']} {a['charge']:.6f} "
                         f"{a['x']:.6f} {a['y']:.6f} {a['z']:.6f} {ix} {iy} {iz}")
        elif style == 'molecular':
            lines.append(f"{a['id']} {a['mol']} {a['type']} "
                         f"{a['x']:.6f} {a['y']:.6f} {a['z']:.6f} {ix} {iy} {iz}")
        elif style in ('charge',):
            lines.append(f"{a['id']} {a['type']} {a['charge']:.6f} "
                         f"{a['x']:.6f} {a['y']:.6f} {a['z']:.6f}")
        else:  # atomic or unknown
            lines.append(f"{a['id']} {a['type']} {a['x']:.6f} {a['y']:.6f} {a['z']:.6f}")

    if data['bonds']:
        lines.append(''); lines.append('Bonds'); lines.append('')
        for b in sorted(data['bonds'], key=lambda x: x['id']):
            lines.append(f"{b['id']} {b['type']} {b['atom1']} {b['atom2']}")

    if data['angles']:
        lines.append(''); lines.append('Angles'); lines.append('')
        for a in sorted(data['angles'], key=lambda x: x['id']):
            lines.append(f"{a['id']} {a['type']} {a['atom1']} {a['atom2']} {a['atom3']}")

    if data['dihedrals']:
        lines.append(''); lines.append('Dihedrals'); lines.append('')
        for d in sorted(data['dihedrals'], key=lambda x: x['id']):
            lines.append(f"{d['id']} {d['type']} {d['atom1']} {d['atom2']} {d['atom3']} {d['atom4']}")

    if data['impropers']:
        lines.append(''); lines.append('Impropers'); lines.append('')
        for imp in sorted(data['impropers'], key=lambda x: x['id']):
            lines.append(f"{imp['id']} {imp['type']} {imp['atom1']} {imp['atom2']} {imp['atom3']} {imp['atom4']}")

    with open(filename, 'w') as f:
        f.write('\n'.join(lines) + '\n')


# ── main API ──────────────────────────────────────────────────────────────────

def create_bubble(
    input_file:          str,
    centers:             list[tuple[float, float, float]],
    radius:              float,
    output_file:         str | None = None,
    water_mol_size:      int            = 3,
    filler_file:         str | None     = None,
    filler_mol_size:     int | None     = None,
    silica_atom_types:   list[int] | None = None,
) -> None:
    """
    Create spherical bubble(s) in a LAMMPS silica/water data file.

    Args:
        input_file:        LAMMPS data file (silica slab + water).
        centers:           List of (x, y, z) bubble center positions in Angstroms.
                           Each center must lie in the water medium.
        radius:            Bubble radius in Angstroms.
                           Must not reach the silica slab from any center.
        output_file:       Path for the modified output file. If omitted, writes beside
                           ``input_file`` as ``<stem>_bubble<suffix>``.
        water_mol_size:    Atoms per water molecule (default 3 for SPC/TIP models).
        filler_file:       Optional LAMMPS data file to carve filler atoms from.
                           Must share the same coordinate system as input_file.
        filler_mol_size:   Atoms per filler molecule (None = auto-detect).
        silica_atom_types: Explicit list of silica atom types for slab detection.
                           None = auto-detect (any atom not in a complete water
                           molecule is treated as silica).
    """
    if output_file is None:
        output_file = default_bubble_output_path(input_file)

    print(f"Reading {input_file} ...")
    data = parse_lammps_data(input_file)
    print(f"  {len(data['atoms'])} atoms  |  {len(data['bonds'])} bonds  |  "
          f"{len(data['angles'])} angles")

    # ── validate bubble placement against silica slab ──
    silica_atoms = identify_silica_atoms(data['atoms'], water_mol_size, silica_atom_types)
    print(f"  {len(silica_atoms)} silica/substrate atoms detected")
    validate_bubble_placement(centers, radius, silica_atoms)
    print(f"  Bubble placement OK — all centers in water, radius clear of silica")

    # ── hollow out bubble(s) ──
    print(f"Removing atoms within bubble(s) (r={radius} Å, {len(centers)} center(s)) ...")
    remove_ids = atoms_in_bubble(data['atoms'], centers, radius)
    print(f"  {len(remove_ids)} atoms directly inside bubble(s)")
    data['atoms'] = [a for a in data['atoms'] if a['id'] not in remove_ids]

    # ── clean up broken water molecules ──
    n_before = len(data['atoms'])
    data['atoms'] = remove_broken_molecules(data['atoms'], water_mol_size)
    print(f"  {n_before - len(data['atoms'])} additional atoms removed (broken water molecules)")
    print(f"  {len(data['atoms'])} water atoms remaining")

    valid_ids = {a['id'] for a in data['atoms']}
    data['bonds']     = _filter_bonds(data['bonds'], valid_ids)
    data['angles']    = _filter_angles(data['angles'], valid_ids)
    data['dihedrals'] = _filter_dihedrals(data['dihedrals'], valid_ids)
    data['impropers'] = _filter_impropers(data['impropers'], valid_ids)

    # ── optional filler ──
    if filler_file:
        print(f"Reading filler file {filler_file} ...")
        filler = parse_lammps_data(filler_file)
        print(f"  {len(filler['atoms'])} atoms in filler")
        print("Carving and merging filler ...")
        data = merge_filler(data, filler, centers, radius, filler_mol_size)

    # ── renumber and write ──
    print("Renumbering IDs ...")
    data = renumber(data)

    print(f"Writing {output_file} ...")
    write_lammps_data(data, output_file)
    print(f"Done. Final system: {len(data['atoms'])} atoms, "
          f"{len(data['bonds'])} bonds, {len(data['angles'])} angles.")


# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    create_bubble(
        input_file=INPUT_FILE,
        centers=BUBBLE_CENTERS,
        radius=BUBBLE_RADIUS,
        output_file=OUTPUT_FILE,
        water_mol_size=WATER_MOL_SIZE,
        filler_file=FILLER_FILE,
        filler_mol_size=FILLER_MOL_SIZE,
        silica_atom_types=SILICA_ATOM_TYPES,
    )
