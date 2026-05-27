import argparse
from ase.io import read
from collections import defaultdict
from itertools import combinations

FORMAT_MAP = {
    '.dump':    'lammps-dump-text',
    '.lammpstrj': 'lammps-dump-text',
    '.xyz':     'xyz',
    '.extxyz':  'extxyz',
}

DEFAULT_TYPE_MAP = {1: 'Si', 2: 'O', 3: 'H'}

def parse_type_map(s):
    """Parse '1=Si,2=O,3=H' into {1: 'Si', 2: 'O', 3: 'H'}."""
    mapping = {}
    for pair in s.split(','):
        k, v = pair.split('=')
        mapping[int(k.strip())] = v.strip()
    return mapping

def remap_types(atoms, type_map):
    from ase import Atoms
    numbers = atoms.arrays.get('type', atoms.get_atomic_numbers())
    symbols = [type_map.get(int(n), f'X{n}') for n in numbers]
    remapped = Atoms(symbols=symbols, positions=atoms.positions,
                     cell=atoms.cell, pbc=atoms.pbc)
    return remapped

def load_frame(filepath, frame=-1, fmt=None, type_map=None):
    from pathlib import Path
    if fmt is None:
        suffix = Path(filepath).suffix.lower()
        name = Path(filepath).name.upper()
        if 'XDATCAR' in name:
            fmt = 'vasp-xdatcar'
        else:
            fmt = FORMAT_MAP.get(suffix)

    is_lammps = fmt and 'lammps' in fmt

    kwargs = {'filename': filepath, 'index': frame}
    if fmt:
        kwargs['format'] = fmt
    atoms = read(**kwargs)
    if isinstance(atoms, list):
        atoms = atoms[0]

    if is_lammps and type_map:
        atoms = remap_types(atoms, type_map)

    return atoms

def analyze(atoms):
    pair_dists = defaultdict(list)
    for i, j in combinations(range(len(atoms)), 2):
        d = atoms.get_distance(i, j, mic=True)
        pair = tuple(sorted([atoms[i].symbol, atoms[j].symbol]))
        pair_dists[pair].append((d, i, j))

    for pair in sorted(pair_dists):
        pair_dists[pair].sort()
        print(f"\n{pair[0]}-{pair[1]}:")
        for d, i, j in pair_dists[pair][:5]:
            print(f"  {atoms[i].symbol}({i})-{atoms[j].symbol}({j}): {d:.3f} Å")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Shortest bond lengths per atom-type pair.')
    parser.add_argument('file', help='Path to structure/trajectory file')
    parser.add_argument('-f', '--frame', type=int, default=-1,
                        help='Frame index to analyze (default: -1, last frame)')
    parser.add_argument('--format', dest='fmt', default=None,
                        help='ASE format string (auto-detected if omitted)')
    parser.add_argument('--types', default=None,
                        help='Type mapping for LAMMPS dumps, e.g. "1=Si,2=O,3=H" '
                             '(default: 1=Si,2=O,3=H)')
    args = parser.parse_args()

    type_map = parse_type_map(args.types) if args.types else DEFAULT_TYPE_MAP
    atoms = load_frame(args.file, frame=args.frame, fmt=args.fmt, type_map=type_map)
    print(f"Frame {args.frame} | {len(atoms)} atoms | {atoms.get_chemical_formula()}")
    analyze(atoms)