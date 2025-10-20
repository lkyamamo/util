#!/usr/bin/env python3
"""
Test script to compare outputs of old and new dipole calculation methods.
This script implements the exact same logic as the old script but uses the new format input.
"""

import numpy as np

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

def create_3atom_bonds(atoms, type_A, type_B, cutoff, box_dims):
    """
    Create 3-atom molecular bonds (same logic as process_bonds_inline_3atom.py).
    """
    bonds = {}
    bond_id = 1
    cutoff_sq = cutoff * cutoff
    
    atom_ids = list(atoms.keys())
    la, lb, lc = box_dims
    
    # Group atoms by type
    type_A_atoms = []
    type_B_atoms = []
    
    for atom_id in atom_ids:
        atom_type, _, _, _ = atoms[atom_id]
        if atom_type == type_A:
            type_A_atoms.append(atom_id)
        elif atom_type == type_B:
            type_B_atoms.append(atom_id)
    
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
    
    return bonds

def calculate_dipole_old_method(atoms, bonds, box_dims):
    """
    Calculate dipole moment using the EXACT same method as the old script.
    Uses hardcoded charge of 0.41 for all atoms.
    """
    totalMoment = [0.0, 0.0, 0.0]
    dipole = []
    
    la, lb, lc = box_dims
    
    for bond_id, atom_ids in bonds.items():
        atomID1, atomID2, atomID3 = sorted(atom_ids)
        
        # Get positions
        _, x1, y1, z1 = atoms[atomID1]
        _, x2, y2, z2 = atoms[atomID2] 
        _, x3, y3, z3 = atoms[atomID3]
        
        # Calculate distances with periodic boundary conditions (same as old script)
        dx12 = x1 - x2
        if dx12 >= la/2.0:
            dx12 -= la
        elif dx12 <= -la/2.0:
            dx12 += la
            
        dx31 = x3 - x1
        if dx31 >= la/2.0:
            dx31 -= la
        elif dx31 <= -la/2.0:
            dx31 += la
            
        dy12 = y1 - y2
        if dy12 >= lb/2.0:
            dy12 -= lb
        elif dy12 <= -lb/2.0:
            dy12 += lb
            
        dy31 = y3 - y1
        if dy31 >= lb/2.0:
            dy31 -= lb
        elif dy31 <= -lb/2.0:
            dy31 += lb
            
        dz12 = z1 - z2
        if dz12 >= lc/2.0:
            dz12 -= lc
        elif dz12 <= -lc/2.0:
            dz12 += lc
            
        dz31 = z3 - z1
        if dz31 >= lc/2.0:
            dz31 -= lc
        elif dz31 <= -lc/2.0:
            dz31 += lc
        
        # Calculate dipole moment using EXACT same formula as old script
        dipoleMoment = [0, 0, 0]
        q = 0.41  # Hardcoded charge from old script
        
        dipoleMoment[0] += (dx12) * q
        dipoleMoment[0] -= (dx31) * q
        
        dipoleMoment[1] += (dy12) * q
        dipoleMoment[1] -= (dy31) * q
        
        dipoleMoment[2] += (dz12) * q
        dipoleMoment[2] -= (dz31) * q
        
        dipole_mag = np.linalg.norm(np.array(dipoleMoment))
        dipole.append(dipole_mag)
        
        totalMoment[0] += dipoleMoment[0]
        totalMoment[1] += dipoleMoment[1]
        totalMoment[2] += dipoleMoment[2]
    
    return totalMoment, dipole

def main():
    # Parameters matching the old script
    filename = "all_lammps_new.xyz"
    type_A = 2  # Central atom type (oxygen)
    type_B = 1  # Terminal atom type (hydrogen)  
    cutoff = 1.2  # Bond cutoff distance
    box_dims = [43.4, 43.4, 55.8]  # Box dimensions
    
    print("=" * 70)
    print("Testing New Format with Old Script Logic")
    print("=" * 70)
    print(f"Input file: {filename}")
    print(f"Type A (central): {type_A}")
    print(f"Type B (terminal): {type_B}")
    print(f"Cutoff: {cutoff}")
    print(f"Box dimensions: {box_dims}")
    print("=" * 70)
    print()
    
    # Read new format file
    print("Reading new format file...")
    atoms = read_new_format(filename)
    print(f"Read {len(atoms)} atoms")
    print()
    
    # Create bonds
    print("Creating 3-atom bonds...")
    bonds = create_3atom_bonds(atoms, type_A, type_B, cutoff, box_dims)
    print(f"Created {len(bonds)} bonds")
    print()
    
    # Calculate dipole using old method
    print("Calculating dipole moment using old script logic...")
    totalMoment, dipole_mags = calculate_dipole_old_method(atoms, bonds, box_dims)
    
    print("Results:")
    print("-" * 70)
    print(f"Mean dipole magnitude: {np.mean(np.array(dipole_mags))}")
    print(f"Box dimensions: {box_dims[0]} {box_dims[1]} {box_dims[2]}")
    print(f"Total dipole moment: [{totalMoment[0]:14.6f}, {totalMoment[1]:14.6f}, {totalMoment[2]:14.6f}]")
    print("=" * 70)
    
    # Write output in same format as old script
    with open("test_comparison_output.txt", "w") as f:
        f.write(f"{totalMoment[0]:14.6f}  {totalMoment[1]:14.6f}  {totalMoment[2]:14.6f}\n")

if __name__ == "__main__":
    main()

