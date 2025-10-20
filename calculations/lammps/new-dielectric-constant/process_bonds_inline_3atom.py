#!/usr/bin/env python3
"""
Process bonds inline from concatenated XYZ format - 3-atom molecular bonds.

Modified to align with old script's method:
- Creates 3-atom bonds (molecular bonds like water)
- Calculates molecular dipole moments
- Uses formula: dipole = (dx12)*q - (dx31)*q
- Uses type-specific charges from dictionary
"""

import sys
import numpy as np
import os

# Type-to-charge mapping (hardcoded)
# Based on old format: type 1 = -2.0, type 2 = 1.0
TYPE_TO_CHARGE = {
    1: -0.65966,  # Central atoms (e.g., oxygen)
    2: 0.32983,   # Terminal atoms (e.g., hydrogen)
}

def read_frame(f):
    """Read a single frame from file handle. Returns (atoms, num_atoms, timestep) or (None, None, None) at EOF."""
    atoms = {}
    
    line = f.readline()
    if not line or not line.strip():
        return None, None, None
    
    num_atoms = int(line.strip())
    timestep_line = f.readline().strip()
    timestep = int(timestep_line.split(':')[1].strip())
    
    for i in range(1, num_atoms + 1):
        line = f.readline()
        if not line:
            raise ValueError(f"Unexpected end of file while reading atom {i}")
        
        parts = line.strip().split()
        atoms[i] = [int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])]
    
    return atoms, num_atoms, timestep


def create_3atom_bonds(atoms, type_A, type_B, cutoff, box_dims):
    """
    Create 3-atom molecular bonds (like water molecules).
    
    Algorithm:
    1. Find all type_A atoms (e.g., oxygen)
    2. For each type_A, find closest 2 type_B atoms (e.g., hydrogens)
    3. Create a bond with these 3 atoms
    
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


def calculate_molecular_dipole(atoms, bonds, box_dims, type_to_charge):
    """
    Calculate total molecular dipole moment using old script's formula.
    
    Formula: dipole = (dx12)*q2 - (dx31)*q3
    Where:
    - atom1 = type_A (e.g., oxygen) with charge q1
    - atom2, atom3 = type_B (e.g., hydrogens) with charges q2, q3
    - dx12 = distance from atom1 to atom2
    - dx31 = distance from atom3 to atom1
    - Uses type-specific charges from dictionary
    """
    total = [0.0, 0.0, 0.0]
    dipole_mags = []
    
    la, lb, lc = box_dims
    
    for bond_id, atom_ids in bonds.items():
        atom_id1, atom_id2, atom_id3 = atom_ids
        
        # Get atom positions and types
        type1, x1, y1, z1 = atoms[atom_id1]
        type2, x2, y2, z2 = atoms[atom_id2]
        type3, x3, y3, z3 = atoms[atom_id3]
        
        # Get charges for each atom type
        q1 = type_to_charge.get(type1, 0.0)
        q2 = type_to_charge.get(type2, 0.0)
        q3 = type_to_charge.get(type3, 0.0)
        
        # Calculate dx12 (atom1 to atom2) with periodic boundary conditions
        # Note: old script uses atom1 - atom2 (not atom2 - atom1)
        dx12 = x1 - x2
        dy12 = y1 - y2
        dz12 = z1 - z2
        
        if dx12 >= la/2.0: dx12 -= la
        elif dx12 <= -la/2.0: dx12 += la
        if dy12 >= lb/2.0: dy12 -= lb
        elif dy12 <= -lb/2.0: dy12 += lb
        if dz12 >= lc/2.0: dz12 -= lc
        elif dz12 <= -lc/2.0: dz12 += lc
        
        # Calculate dx31 (atom3 to atom1) with periodic boundary conditions
        dx31 = x3 - x1
        dy31 = y3 - y1
        dz31 = z3 - z1
        
        if dx31 >= la/2.0: dx31 -= la
        elif dx31 <= -la/2.0: dx31 += la
        if dy31 >= lb/2.0: dy31 -= lb
        elif dy31 <= -lb/2.0: dy31 += lb
        if dz31 >= lc/2.0: dz31 -= lc
        elif dz31 <= -lc/2.0: dz31 += lc
        
        # Calculate dipole moment using old script's formula with type-specific charges
        dipole_moment = [
            (dx12) * q2 - (dx31) * q3,
            (dy12) * q2 - (dy31) * q3,
            (dz12) * q2 - (dz31) * q3
        ]
        
        dipole_mag = np.linalg.norm(dipole_moment)
        dipole_mags.append(dipole_mag)
        
        total[0] += dipole_moment[0]
        total[1] += dipole_moment[1]
        total[2] += dipole_moment[2]
    
    return total, dipole_mags


def load_processed_timesteps(output_file):
    """
    Load already processed timesteps from output file.
    Returns set of timesteps that have been completed.
    """
    processed = set()
    if not os.path.exists(output_file):
        return processed
    
    try:
        with open(output_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    # Try to parse timestep from line
                    # Format: "timestep  dipole_x  dipole_y  dipole_z" or "dipole_x  dipole_y  dipole_z"
                    parts = line.split()
                    if len(parts) >= 3:
                        # Check if first part is a timestep (integer)
                        try:
                            timestep = int(parts[0])
                            processed.add(timestep)
                        except ValueError:
                            # Not a timestep, just dipole values
                            pass
    except Exception as e:
        print(f"Warning: Could not read existing output file: {e}")
    
    return processed


def process_concatenated_file(filename, start, end, increment, type_A, type_B, 
                              cutoff, box_dims, type_to_charge, output_file):
    """
    Process concatenated frames and write dipole moments to output file.
    
    Features:
    - Incremental writing: writes each result immediately
    - Restart capability: skips already processed timesteps
    """
    # Load already processed timesteps for restart capability
    processed_timesteps = load_processed_timesteps(output_file)
    
    if processed_timesteps:
        print(f"Restart mode: Found {len(processed_timesteps)} already processed timesteps")
        print(f"Will skip: {sorted(list(processed_timesteps))[:10]}..." if len(processed_timesteps) > 10 else f"Will skip: {sorted(list(processed_timesteps))}")
    else:
        print("Starting fresh - no previous results found")
    
    # Open output file in append mode for incremental writing
    f_out = open(output_file, 'a')
    
    try:
        with open(filename, 'r') as f_in:
            counter = 0
            skipped = 0
            
            while True:
                atoms, num_atoms, timestep = read_frame(f_in)
                if atoms is None:
                    break
                
                # Check if this timestep should be processed
                if timestep < start or timestep >= end:
                    continue
                if (timestep - start) % increment != 0:
                    continue
                
                # Skip if already processed (restart capability)
                if timestep in processed_timesteps:
                    skipped += 1
                    print(f"Skipping timestep {timestep} (already processed)")
                    continue
                
                try:
                    bonds, bond_stats = create_3atom_bonds(atoms, type_A, type_B, cutoff, box_dims)
                    dipole, dipole_mags = calculate_molecular_dipole(atoms, bonds, box_dims, type_to_charge)
                    
                    mean_dipole = np.mean(dipole_mags) if dipole_mags else 0.0
                    
                    print(f"Timestep {timestep}: {len(bonds)} bonds, "
                          f"dipole = [{dipole[0]:14.6f}, {dipole[1]:14.6f}, {dipole[2]:14.6f}], "
                          f"mean_mag = {mean_dipole:.6f}")
                    print(f"Bond stats: {bond_stats['bonds_created']} created, "
                          f"{bond_stats['bonds_missing']} missing, "
                          f"{bond_stats['atoms_within_cutoff']} type_B within cutoff, "
                          f"type_A: {bond_stats['type_A_count']}, type_B: {bond_stats['type_B_count']}")
                    
                    # Write immediately with timestep for restart capability
                    f_out.write(f"{timestep}  {dipole[0]:14.6f}  {dipole[1]:14.6f}  {dipole[2]:14.6f}\n")
                    f_out.flush()  # Ensure data is written to disk
                    counter += 1
                    
                except Exception as e:
                    print(f"Error processing timestep {timestep}: {e}")
                    import traceback
                    traceback.print_exc()
                    continue
            
            print(f"\nProcessed {counter} new frames")
            if skipped > 0:
                print(f"Skipped {skipped} already processed frames")
            print(f"Output written to: {output_file}")
    
    finally:
        f_out.close()


if __name__ == "__main__":
    if len(sys.argv) < 12:
        print("Usage:")
        print("  python process_bonds_inline_3atom.py <input_file> <start> <end> <increment> "
              "<type_A> <type_B> <cutoff> <la> <lb> <lc> <output_file>")
        print("\nExample:")
        print("  python process_bonds_inline_3atom.py trajectory.dat 0 1000 10 1 2 1.2 "
              "35.0 35.0 35.0 dipole_output.txt")
        print("\nNote:")
        print("  - Creates 3-atom molecular bonds (type_A + 2*type_B)")
        print("  - Uses molecular dipole formula: (dx12)*q2 - (dx31)*q3")
        print("  - Uses type-specific charges from hardcoded dictionary:")
        for atom_type, charge in sorted(TYPE_TO_CHARGE.items()):
            print(f"    Type {atom_type}: {charge}")
        sys.exit(1)
    
    filename = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    increment = int(sys.argv[4])
    type_A = int(sys.argv[5])
    type_B = int(sys.argv[6])
    cutoff = float(sys.argv[7])
    la = float(sys.argv[8])
    lb = float(sys.argv[9])
    lc = float(sys.argv[10])
    output_file = sys.argv[11]
    
    print("=" * 70)
    print("Inline Bond Processing - 3-Atom Molecular Bonds")
    print("=" * 70)
    print(f"Input file: {filename}")
    print(f"Timesteps: {start} to {end} (step {increment})")
    print(f"Bond types: {type_A} (central) + 2*{type_B} (terminal)")
    print(f"Cutoff: {cutoff}")
    print(f"Box dimensions: {la} x {lb} x {lc}")
    print(f"Type-specific charges:")
    for atom_type, charge in sorted(TYPE_TO_CHARGE.items()):
        print(f"  Type {atom_type}: {charge}")
    print(f"Output file: {output_file}")
    print("=" * 70)
    print()
    
    process_concatenated_file(filename, start, end, increment, type_A, type_B,
                              cutoff, [la, lb, lc], TYPE_TO_CHARGE, output_file)

