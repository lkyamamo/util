#!/usr/bin/env python3
"""
Create example files in both formats for testing.

Creates:
1. example_concatenated.dat - Concatenated file (new format)
2. example_old_format/ - Directory with separate files (old format)
"""

import os

# Create directory for old format files
os.makedirs('example_old_format', exist_ok=True)

# Box dimensions
la = lb = lc = 10.0

# Charges
charge_A = -2.0  # Type 1 (central atom, e.g., oxygen)
charge_B = 1.0   # Type 2 (terminal atom, e.g., hydrogen)

# Create 3 asymmetric water-like molecules
# Molecule 1: O at (0, 0, 0), H1 at (0.8, 0, 0), H2 at (-0.3, 0, 0)
# Molecule 2: O at (3, 0, 0), H1 at (3.9, 0, 0), H2 at (2.6, 0, 0)
# Molecule 3: O at (6, 0, 0), H1 at (6.7, 0, 0), H2 at (5.4, 0, 0)

atoms_data = [
    # Molecule 1 (asymmetric)
    (1, 1, 0.0, 0.0, 0.0),      # O1
    (2, 1, 0.8, 0.0, 0.0),      # H1 (far)
    (3, 1, -0.3, 0.0, 0.0),     # H2 (near)
    # Molecule 2 (asymmetric)
    (4, 2, 3.0, 0.0, 0.0),      # O2
    (5, 2, 3.9, 0.0, 0.0),      # H3 (far)
    (6, 2, 2.6, 0.0, 0.0),      # H4 (near)
    # Molecule 3 (asymmetric)
    (7, 3, 6.0, 0.0, 0.0),      # O3
    (8, 3, 6.7, 0.0, 0.0),      # H5 (far)
    (9, 3, 5.4, 0.0, 0.0),      # H6 (near)
]

# Bond mapping: each molecule gets one bond ID
bond_mapping = {
    1: 1, 2: 1, 3: 1,  # Molecule 1
    4: 2, 5: 2, 6: 2,  # Molecule 2
    7: 3, 8: 3, 9: 3,  # Molecule 3
}

print("Creating example files...")
print(f"Total atoms: {len(atoms_data)}")
print(f"Total molecules: 3")
print("\nMolecule structure (asymmetric):")
print("  O at center, H1 at +0.8, H2 at -0.3")
print("  (H atoms at different distances from O)")
print()

# Create concatenated file (new format)
with open('example_concatenated.dat', 'w') as f:
    for frame in range(6):
        timestep = frame * 10
        shift = frame * 0.1
        
        f.write(f"{len(atoms_data)}\n")
        f.write(f"Atoms. Timestep: {timestep}\n")
        
        for atom_id, bond_id, x, y, z in atoms_data:
            if atom_id in [1, 4, 7]:
                atom_type = 1
            else:
                atom_type = 2
            
            f.write(f"{atom_type} {x+shift:.6f} {y:.6f} {z:.6f}\n")

print("✓ Created: example_concatenated.dat")

# Create old format files (separate files)
for frame in range(6):
    timestep = frame * 10
    shift = frame * 0.1
    
    with open(f'example_old_format/{timestep}', 'w') as f:
        # Header
        f.write("# Comment line 1\n")
        f.write("# Comment line 2\n")
        f.write(f"{len(atoms_data)}\n")
        f.write("# Comment line 4\n")
        f.write(f"0.0 {la} 0.0 {lb} 0.0 {lc}\n")
        f.write(f"0.0 {la} 0.0 {lb} 0.0 {lc}\n")
        f.write(f"0.0 {la} 0.0 {lb} 0.0 {lc}\n")
        f.write("# Comment line 8\n")
        
        # Atoms
        for atom_id, bond_id, x, y, z in atoms_data:
            if atom_id in [1, 4, 7]:
                atom_type = 1
                charge = charge_A
            else:
                atom_type = 2
                charge = charge_B
            
            f.write(f"{atom_id} {bond_id} {atom_type} {x+shift:.6f} {y:.6f} {z:.6f} {charge:.6f}\n")
    
    print(f"✓ Created: example_old_format/{timestep}")

print("\nExample files created successfully!")
print("\nCharges:")
print(f"  Type A (central): {charge_A}")
print(f"  Type B (terminal): {charge_B}")
print("\nExpected dipole per molecule:")
print("  dx12 = 0.0 - 0.8 = -0.8")
print("  dx31 = -0.3 - 0.0 = -0.3")
print(f"  dipole = (-0.8)*{charge_B} - (-0.3)*{charge_B} = {-0.8*charge_B} - {-0.3*charge_B} = {-0.5*charge_B}")
print(f"  Total (3 molecules): {-1.5*charge_B}")

