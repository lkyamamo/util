#!/usr/bin/env python3
"""
Script to extract timestep and thermodynamic data from VASP OSZICAR file.

The OSZICAR file contains timestep information in lines like:
   391 T=    63. E= -.27130923E+04 F= -.27141250E+04 E0= -.27141250E+04  EK= 0.11649E+01 SP= -.14E+00 SK= 0.46E-02

This script extracts:
- Timestep number
- T (Temperature)
- E (Energy)
- F (Free Energy)
- E0 (Energy at 0 K)
- EK (Kinetic Energy)
- SP (Entropy Pressure)
- SK (Entropy Kinetic)
"""

import re
import sys
from pathlib import Path


def parse_oszicar(filepath):
    """
    Parse OSZICAR file and extract timestep data.
    
    Parameters:
    -----------
    filepath : str or Path
        Path to the OSZICAR file
        
    Returns:
    --------
    list of dict
        List of dictionaries containing timestep data
    """
    # Pattern to match timestep lines
    # Format: "   391 T=    63. E= -.27130923E+04 F= -.27141250E+04 E0= -.27141250E+04  EK= 0.11649E+01 SP= -.14E+00 SK= 0.46E-02"
    pattern = r'^\s*(\d+)\s+T=\s*([\d.]+)\s+E=\s*([\d.E+-]+)\s+F=\s*([\d.E+-]+)\s+E0=\s*([\d.E+-]+)\s+EK=\s*([\d.E+-]+)\s+SP=\s*([\d.E+-]+)\s+SK=\s*([\d.E+-]+)'
    
    timesteps = []
    
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            match = re.match(pattern, line)
            if match:
                timestep_data = {
                    'timestep': int(match.group(1)),
                    'T': float(match.group(2)),
                    'E': float(match.group(3)),
                    'F': float(match.group(4)),
                    'E0': float(match.group(5)),
                    'EK': float(match.group(6)),
                    'SP': float(match.group(7)),
                    'SK': float(match.group(8))
                }
                timesteps.append(timestep_data)
    
    return timesteps


def main():
    """Main function to run the parser."""
    # Default file path
    default_file = Path(__file__).parent / 'output' / 'OSZICAR'
    
    # Allow command line argument for file path
    if len(sys.argv) > 1:
        filepath = Path(sys.argv[1])
    else:
        filepath = default_file
    
    if not filepath.exists():
        print(f"Error: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)
    
    # Parse the file
    timesteps = parse_oszicar(filepath)
    
    if not timesteps:
        print("Warning: No timestep data found in OSZICAR file.", file=sys.stderr)
        sys.exit(1)
    
    # Print results
    print(f"Found {len(timesteps)} timesteps")
    print("\nTimestep Data:")
    print("-" * 100)
    print(f"{'Step':<6} {'T':<8} {'E':<15} {'F':<15} {'E0':<15} {'EK':<15} {'SP':<12} {'SK':<12}")
    print("-" * 100)
    
    for data in timesteps:
        print(f"{data['timestep']:<6} "
              f"{data['T']:<8.2f} "
              f"{data['E']:<15.6E} "
              f"{data['F']:<15.6E} "
              f"{data['E0']:<15.6E} "
              f"{data['EK']:<15.6E} "
              f"{data['SP']:<12.6E} "
              f"{data['SK']:<12.6E}")
    
    # Optionally save to CSV
    if len(sys.argv) > 2 and sys.argv[2] == '--csv':
        csv_file = filepath.parent / 'oszicar_data.csv'
        import csv
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=timesteps[0].keys())
            writer.writeheader()
            writer.writerows(timesteps)
        print(f"\nData saved to: {csv_file}")


if __name__ == '__main__':
    main()

