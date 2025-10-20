#!/usr/bin/env python3
"""
Generate multiple unit test cases with different configurations.
"""

import os
import shutil
from pathlib import Path

def create_test_case(base_dir, test_name, config):
    """Create a test case directory with input files."""
    test_dir = base_dir / test_name
    test_dir.mkdir(exist_ok=True)
    
    # Create old format directory
    old_format_dir = test_dir / "old_format"
    old_format_dir.mkdir(exist_ok=True)
    
    # Create new format file
    new_format_file = test_dir / "concatenated.dat"
    
    # Generate data for each timestep
    timesteps = config['timesteps']
    num_molecules = config['num_molecules']
    box_size = config['box_size']
    
    # Generate old format files
    for timestep in timesteps:
        old_file = old_format_dir / str(timestep)
        with open(old_file, 'w') as f:
            f.write("# Comment line 1\n")
            f.write("# Comment line 2\n")
            f.write(f"{num_molecules * 3}\n")
            f.write("# Comment line 4\n")
            f.write(f"0.0 {box_size} 0.0 {box_size} 0.0 {box_size}\n")
            f.write(f"0.0 {box_size} 0.0 {box_size} 0.0 {box_size}\n")
            f.write(f"0.0 {box_size} 0.0 {box_size} 0.0 {box_size}\n")
            f.write("# Comment line 8\n")
            
            atom_id = 1
            for mol_id in range(1, num_molecules + 1):
                # Central atom (type 1)
                x = mol_id * 3.0 + timestep * 0.1
                y = 0.0
                z = 0.0
                f.write(f"{atom_id} {mol_id} 1 {x:.6f} {y:.6f} {z:.6f} -0.8476\n")
                atom_id += 1
                
                # Terminal atom 1 (type 2)
                x2 = x + 0.8
                f.write(f"{atom_id} {mol_id} 2 {x2:.6f} {y:.6f} {z:.6f} 0.4238\n")
                atom_id += 1
                
                # Terminal atom 2 (type 2)
                x3 = x - 0.3
                f.write(f"{atom_id} {mol_id} 2 {x3:.6f} {y:.6f} {z:.6f} 0.4238\n")
                atom_id += 1
    
    # Generate new format file
    with open(new_format_file, 'w') as f:
        for timestep in timesteps:
            f.write(f"{num_molecules * 3}\n")
            f.write(f"Atoms. Timestep: {timestep}\n")
            
            for mol_id in range(1, num_molecules + 1):
                # Central atom (type 1)
                x = mol_id * 3.0 + timestep * 0.1
                y = 0.0
                z = 0.0
                f.write(f"1 {x:.6f} {y:.6f} {z:.6f}\n")
                
                # Terminal atom 1 (type 2)
                x2 = x + 0.8
                f.write(f"2 {x2:.6f} {y:.6f} {z:.6f}\n")
                
                # Terminal atom 2 (type 2)
                x3 = x - 0.3
                f.write(f"2 {x3:.6f} {y:.6f} {z:.6f}\n")
    
        # Create config file
    config_file = test_dir / "config.txt"
    increment = timesteps[1] - timesteps[0] if len(timesteps) > 1 else 10
    with open(config_file, 'w') as f:
        f.write(f"Test Name: {test_name}\n")
        f.write(f"Number of Molecules: {num_molecules}\n")
        f.write(f"Box Size: {box_size}\n")
        f.write(f"Timesteps: {timesteps}\n")
        f.write(f"Start: {timesteps[0]}\n")
        f.write(f"End: {timesteps[-1] + increment}\n")
        f.write(f"Increment: {increment}\n")
    
    print(f"Created test case: {test_name}")
    return test_dir


def main():
    base_dir = Path("unit_tests")
    
    # Define test configurations
    test_configs = [
        {
            'name': 'test_01_small_system',
            'config': {
                'num_molecules': 3,
                'box_size': 10.0,
                'timesteps': [0, 10, 20, 30, 40, 50]
            }
        },
        {
            'name': 'test_02_medium_system',
            'config': {
                'num_molecules': 5,
                'box_size': 15.0,
                'timesteps': [0, 5, 10, 15, 20]
            }
        },
        {
            'name': 'test_03_large_system',
            'config': {
                'num_molecules': 10,
                'box_size': 30.0,
                'timesteps': [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
            }
        },
        {
            'name': 'test_04_small_box',
            'config': {
                'num_molecules': 2,
                'box_size': 5.0,
                'timesteps': [0, 10, 20, 30, 40]
            }
        },
        {
            'name': 'test_05_large_box',
            'config': {
                'num_molecules': 7,
                'box_size': 50.0,
                'timesteps': [0, 5, 10, 15, 20, 25, 30]
            }
        },
        {
            'name': 'test_06_few_timesteps',
            'config': {
                'num_molecules': 4,
                'box_size': 12.0,
                'timesteps': [0, 20, 40]
            }
        },
        {
            'name': 'test_07_many_timesteps',
            'config': {
                'num_molecules': 3,
                'box_size': 10.0,
                'timesteps': [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
            }
        },
        {
            'name': 'test_08_single_molecule',
            'config': {
                'num_molecules': 1,
                'box_size': 5.0,
                'timesteps': [0, 10, 20, 30, 40, 50]
            }
        },
        {
            'name': 'test_09_irregular_timesteps',
            'config': {
                'num_molecules': 4,
                'box_size': 15.0,
                'timesteps': [0, 3, 7, 12, 18, 25]
            }
        },
        {
            'name': 'test_10_dense_packing',
            'config': {
                'num_molecules': 8,
                'box_size': 10.0,
                'timesteps': [0, 5, 10, 15, 20]
            }
        }
    ]
    
    print(f"Generating {len(test_configs)} test cases...")
    print("=" * 70)
    
    for test_config in test_configs:
        create_test_case(
            base_dir,
            test_config['name'],
            test_config['config']
        )
    
    print("=" * 70)
    print(f"Successfully generated {len(test_configs)} test cases in {base_dir}/")
    
    # Create README
    readme_file = base_dir / "README.md"
    with open(readme_file, 'w') as f:
        f.write("# Unit Test Cases\n\n")
        f.write("This directory contains multiple unit test cases for validating the dipole moment calculation scripts.\n\n")
        f.write("## Test Cases\n\n")
        f.write("| Test | Molecules | Box Size | Timesteps | Description |\n")
        f.write("|------|-----------|----------|-----------|-------------|\n")
        
        for test_config in test_configs:
            config = test_config['config']
            f.write(f"| {test_config['name']} | {config['num_molecules']} | {config['box_size']} | {len(config['timesteps'])} | {test_config['name'].replace('_', ' ').title()} |\n")
        
        f.write("\n## Directory Structure\n\n")
        f.write("Each test case contains:\n")
        f.write("- `old_format/` - Old format input files (one per timestep)\n")
        f.write("- `concatenated.dat` - New format concatenated file\n")
        f.write("- `config.txt` - Test configuration\n\n")
        f.write("## Running Tests\n\n")
        f.write("Use the `run_all_unit_tests.py` script to run all tests:\n")
        f.write("```bash\n")
        f.write("python run_all_unit_tests.py\n")
        f.write("```\n")


if __name__ == "__main__":
    main()

