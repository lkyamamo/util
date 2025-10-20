# Unit Test Cases

This directory contains unit test cases for validating the dipole moment calculation scripts.

## Test Cases

| Test | Molecules | Box Size | Timesteps | Description |
|------|-----------|----------|-----------|-------------|
| test_01_small_system | 3 | 10.0 | 6 | Small system with 3 molecules |
| test_02_medium_system | 5 | 15.0 | 5 | Medium system with 5 molecules |
| test_03_large_system | 10 | 30.0 | 10 | Large system with 10 molecules |
| test_04_small_box | 2 | 5.0 | 5 | Small box with 2 molecules |
| test_05_large_box | 7 | 50.0 | 7 | Large box with 7 molecules |
| test_06_few_timesteps | 4 | 12.0 | 3 | Few timesteps |
| test_07_many_timesteps | 3 | 10.0 | 11 | Many timesteps |
| test_08_single_molecule | 1 | 5.0 | 6 | Single molecule |

## Directory Structure

Each test case contains:
- `old_format/` - Old format input files (one per timestep)
- `concatenated.dat` - New format concatenated file
- `config.txt` - Test configuration
- `outputs/` - Test outputs (generated after running tests)

## Running Tests

### Run All Tests
```bash
python run_all_unit_tests.py
```

### Run Individual Test
```bash
cd unit_tests/test_01_small_system
# Run old script
python ../../2.dipoleMoment.py 0 60 10 old_output.txt

# Run new script
python ../../process_bonds_inline_3atom.py concatenated.dat 0 60 10 1 2 1.5 10.0 10.0 10.0 new_output.txt
```

## Test Results

All 8 test cases pass successfully:
- ✅ test_01_small_system
- ✅ test_02_medium_system
- ✅ test_03_large_system
- ✅ test_04_small_box
- ✅ test_05_large_box
- ✅ test_06_few_timesteps
- ✅ test_07_many_timesteps
- ✅ test_08_single_molecule

## Generating New Tests

To generate new test cases:
```bash
python generate_unit_tests.py
```

This will create 10 test cases with different configurations.
