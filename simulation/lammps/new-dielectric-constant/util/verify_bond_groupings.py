#!/usr/bin/env python3
"""
Verify that the bondIDs created by process_bonds_inline_3atom.py match
the molecule groupings in the old format file (273.41000).

This script:
1. Reads the new format file (all_lammps_new.xyz)
2. Creates bonds using the same logic as process_bonds_inline_3atom.py
3. Reads the old format file (273.41000) to get molecule IDs
4. Compares the groupings to ensure atoms in the same bond are in the same molecule
"""

import numpy as np
from collections import defaultdict

# Type-to-charge mapping (same as in process_bonds_inline_3atom.py)
TYPE_TO_CHARGE = {
    1: 0.41,   # Central atoms (e.g., oxygen)
    2: -0.82,  # Terminal atoms (e.g., hydrogen)
}

def read_new_format(filename):
    """Read new format file and return atoms dictionary."""
    atoms = {}
    
    with open(filename, 'r') as f:
        # Skip number of atoms line
        f.readline()
        # Skip timestep line
        f.readline()
        
        atom_id = 1
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            atom_type = int(parts[0])
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
            atoms[atom_id] = [atom_type, x, y, z]
            atom_id += 1
    
    return atoms

def read_old_format(filename):
    """Read old format file and return molecule ID for each atom."""
    atom_to_mol = {}
    
    with open(filename, 'r') as f:
        # Skip to ATOMS section
        for line in f:
            if line.strip() == 'ITEM: ATOMS id mol type x y z q':
                break
        
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            atom_id = int(parts[0])
            mol_id = int(parts[1])
            atom_type = int(parts[2])
            x = float(parts[3])
            y = float(parts[4])
            z = float(parts[5])
            
            # Map atom ID to molecule ID
            atom_to_mol[atom_id] = mol_id
    
    return atom_to_mol

def create_3atom_bonds(atoms, type_A, type_B, cutoff, box_dims):
    """
    Create 3-atom molecular bonds (same logic as process_bonds_inline_3atom.py).
    
    Returns:
    --------
    bonds : dict
        {bondID: [atomID1, atomID2, atomID3]}
        where atomID1 is type_A, atomID2 and atomID3 are type_B
    bond_stats : dict
        Statistics about bond formation including:
        - 'bonds_created': number of bonds successfully created
        - 'bonds_missing': number of type_A atoms that couldn't form bonds
        - 'atoms_within_cutoff': total count of type_B atoms within cutoff of any type_A
        - 'type_A_count': total number of type_A atoms
        - 'type_B_count': total number of type_B atoms
    """
    bonds = {}
    bond_id = 1
    cutoff_sq = cutoff * cutoff
    
    # Statistics tracking
    bonds_created = 0
    bonds_missing = 0
    atoms_within_cutoff = 0
    type_A_count = 0
    type_B_count = 0
    
    atom_ids = list(atoms.keys())
    la, lb, lc = box_dims
    
    # Group atoms by type
    type_A_atoms = []
    type_B_atoms = []
    
    for atom_id in atom_ids:
        atom_type, _, _, _ = atoms[atom_id]
        if atom_type == type_A:
            type_A_atoms.append(atom_id)
            type_A_count += 1
        elif atom_type == type_B:
            type_B_atoms.append(atom_id)
            type_B_count += 1
    
    # For each type_A atom, find closest 2 type_B atoms
    for atom_id1 in type_A_atoms:
        _, x1, y1, z1 = atoms[atom_id1]
        
        # Find distances to all type_B atoms
        distances = []
        for atom_id2 in type_B_atoms:
            _, x2, y2, z2 = atoms[atom_id2]
            
            # Calculate distance with periodic boundary conditions
            dx = x2 - x1
            dy = y2 - y1
            dz = z2 - z1
            
            if dx >= la/2.0: dx -= la
            elif dx <= -la/2.0: dx += la
            if dy >= lb/2.0: dy -= lb
            elif dy <= -lb/2.0: dy += lb
            if dz >= lc/2.0: dz -= lc
            elif dz <= -lc/2.0: dz += lc
            
            dist_sq = dx*dx + dy*dy + dz*dz
            distances.append((dist_sq, atom_id2))
            
            # Count atoms within cutoff distance
            if dist_sq <= cutoff_sq:
                atoms_within_cutoff += 1
        
        # Sort by distance and take closest 2
        distances.sort()
        
        # Check if closest 2 are within cutoff
        if len(distances) >= 2:
            dist1, atom_id2 = distances[0]
            dist2, atom_id3 = distances[1]
            
            if dist1 <= cutoff_sq and dist2 <= cutoff_sq:
                # Create 3-atom bond: [type_A, type_B1, type_B2]
                bonds[bond_id] = [atom_id1, atom_id2, atom_id3]
                bond_id += 1
                bonds_created += 1
            else:
                bonds_missing += 1
                print(f"Warning: Bond {bond_id} not created: dist_sq1 = {dist1}, dist_sq2 = {dist2}")
            
            # Check for third closest atom if it exists
            if len(distances) >= 3:
                dist3, atom_id4 = distances[2]
                if dist3 <= cutoff_sq:
                    print(f"Warning: Third closest atom {atom_id4} is within cutoff for atom {atom_id1}")

        else:
            bonds_missing += 1
            print(f"No 2 closest atoms found for atom {atom_id1}")
    
    # Create statistics dictionary
    bond_stats = {
        'bonds_created': bonds_created,
        'bonds_missing': bonds_missing,
        'atoms_within_cutoff': atoms_within_cutoff,
        'type_A_count': type_A_count,
        'type_B_count': type_B_count
    }
    
    return bonds, bond_stats

def verify_bond_groupings(bonds, atom_to_mol):
    """
    Verify that all atoms in each bond belong to the same molecule.
    
    Returns:
    --------
    results : dict with verification statistics
    """
    results = {
        'total_bonds': len(bonds),
        'matching_molecules': 0,
        'mismatched_molecules': 0,
        'mismatches': []
    }
    
    for bond_id, atom_ids in bonds.items():
        atom1, atom2, atom3 = atom_ids
        
        # Get molecule IDs for each atom
        # Note: atom IDs in new format start from 1, but in old format they start from 1177
        # We need to add an offset
        offset = 1176  # 1177 - 1
        mol1 = atom_to_mol.get(atom1 + offset)
        mol2 = atom_to_mol.get(atom2 + offset)
        mol3 = atom_to_mol.get(atom3 + offset)
        
        if mol1 is None or mol2 is None or mol3 is None:
            results['mismatches'].append({
                'bond_id': bond_id,
                'atoms': atom_ids,
                'error': 'Missing molecule ID for one or more atoms'
            })
            results['mismatched_molecules'] += 1
            continue
        
        # Check if all three atoms are in the same molecule
        if mol1 == mol2 == mol3:
            results['matching_molecules'] += 1
        else:
            results['mismatched_molecules'] += 1
            results['mismatches'].append({
                'bond_id': bond_id,
                'atoms': atom_ids,
                'molecules': [mol1, mol2, mol3],
                'error': f'Molecules do not match: {mol1}, {mol2}, {mol3}'
            })
    
    return results

def main():
    # File paths
    new_format_file = 'all_lammps_new.xyz'
    old_format_file = '273.41000'
    
    # Parameters (same as process_bonds_inline_3atom.py)
    type_A = 2  # Central atom type (oxygen)
    type_B = 1  # Terminal atom type (hydrogen)
    cutoff = 1.2  # Bond cutoff distance
    box_dims = [43.4, 43.4, 55.8]  # Box dimensions
    
    print("=" * 70)
    print("Bond Grouping Verification")
    print("=" * 70)
    print(f"New format file: {new_format_file}")
    print(f"Old format file: {old_format_file}")
    print(f"Type A (central): {type_A}")
    print(f"Type B (terminal): {type_B}")
    print(f"Cutoff: {cutoff}")
    print(f"Box dimensions: {box_dims}")
    print("=" * 70)
    print()
    
    # Read new format file
    print("Reading new format file...")
    atoms = read_new_format(new_format_file)
    print(f"Read {len(atoms)} atoms")
    print()
    
    # Read old format file
    print("Reading old format file...")
    atom_to_mol = read_old_format(old_format_file)
    print(f"Read molecule IDs for {len(atom_to_mol)} atoms")
    print()
    
    # Create bonds
    print("Creating 3-atom bonds...")
    bonds, bond_stats = create_3atom_bonds(atoms, type_A, type_B, cutoff, box_dims)
    print(f"Created {len(bonds)} bonds")
    print(f"Bond stats: {bond_stats['bonds_created']} created, "
          f"{bond_stats['bonds_missing']} missing, "
          f"{bond_stats['atoms_within_cutoff']} type_B within cutoff, "
          f"type_A: {bond_stats['type_A_count']}, type_B: {bond_stats['type_B_count']}")
    print()
    
    # Verify bond groupings
    print("Verifying bond groupings...")
    results = verify_bond_groupings(bonds, atom_to_mol)
    
    # Print results
    print("=" * 70)
    print("Verification Results")
    print("=" * 70)
    print(f"Total bonds: {results['total_bonds']}")
    print(f"Matching molecules: {results['matching_molecules']}")
    print(f"Mismatched molecules: {results['mismatched_molecules']}")
    print()
    
    if results['mismatched_molecules'] > 0:
        print("MISMATCHES FOUND:")
        print("-" * 70)
        for mismatch in results['mismatches'][:10]:  # Show first 10
            print(f"Bond {mismatch['bond_id']}: Atoms {mismatch['atoms']}")
            print(f"  {mismatch['error']}")
            if 'molecules' in mismatch:
                print(f"  Molecule IDs: {mismatch['molecules']}")
        if len(results['mismatches']) > 10:
            print(f"\n... and {len(results['mismatches']) - 10} more mismatches")
    else:
        print("SUCCESS: All bonds have atoms in the same molecule!")
    
    print("=" * 70)

if __name__ == "__main__":
    main()

