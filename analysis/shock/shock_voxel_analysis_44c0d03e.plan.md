---
name: Shock Voxel Analysis
overview: "Implement a streaming LAMMPS dump parser that bins atoms into 1 nm³ voxels, computes per-voxel physics quantities, and writes results to HDF5. Single file: `analysis/shock/shock_analysis.py`."
todos:
  - id: parse-header
    content: "Implement parse_header: reads LAMMPS dump header, returns metadata dict with box bounds, voxel grid dims, type_to_mass array, Si_type and O_type integers"
    status: pending
  - id: preallocate-hdf5
    content: "Implement preallocate_hdf5: creates all 8 datasets with correct shapes, dtypes, chunk shapes, fillvalues, and writes metadata as file attributes"
    status: pending
  - id: process-voxel
    content: "Implement process_voxel: takes float64 numpy array, computes density, full pressure, virial pressure, temperature, avg_speed, avg_O_speed, voxel_type, v_COM using vectorized numpy ops"
    status: pending
  - id: streaming-loop
    content: "Implement main streaming loop: reads atom lines into layer buffer dict keyed by (iy,iz), triggers flush_layer on ix change and at EOF"
    status: pending
  - id: flush-layer
    content: "Implement flush_layer: iterates all (iy,iz) in layer, calls process_voxel or writes NaN for empty voxels, writes float32-cast results into HDF5 slices"
    status: pending
  - id: cli-and-units
    content: Wire up argparse CLI (dump_file, output.h5, --voxel-size, --type-map), add unit conversion constants for pressure (to GPa), temperature (k_B in kcal/mol/K), and density (amu/A^3 to g/cm^3)
    status: pending
isProject: false
---

# Shock Voxel Analysis Plan

## Output file
`analysis/shock/shock_analysis.py`

## Architecture

```mermaid
flowchart TD
    A[LAMMPS dump file] --> B[parse_header]
    B --> C["metadata: timestep, N, Lx/Ly/Lz, nx/ny/nz, type_to_mass"]
    C --> D[preallocate HDF5 datasets]
    A --> E[stream atom lines]
    E --> F["layer buffer: dict[(iy,iz)] -> list of rows"]
    F --> G{new ix layer?}
    G -->|yes| H[flush previous layer]
    H --> I[for each voxel in layer]
    I --> J{len == 0?}
    J -->|yes| K[write NaN slice to HDF5]
    J -->|no| L[process_voxel]
    L --> M[write float32 slice to HDF5]
    G -->|no| E
```

## Entry point
```python
python shock_analysis.py <dump_file> <output.h5> --voxel-size 1.0 --type-map 1:H 2:O 3:Si
```

## Functions

### `parse_header(f) -> dict`
Read lines until `ITEM: ATOMS` line; return:
- `timestep`, `N`, `xlo/xhi/ylo/yhi/zlo/zhi`, `Lx/Ly/Lz`
- `nx, ny, nz = ceil(Lx/vs), ceil(Ly/vs), ceil(Lz/vs)`
- `type_to_mass`: numpy array indexed by int type, values in amu
- `Si_type`, `O_type`: int type IDs for Si and O from type mapping

### `preallocate_hdf5(h5file, nx, ny, nz, attrs) -> h5py.File`
Create all 8 datasets with `fillvalue=np.nan` (or 0 for `voxel_type`), chunk shape `(1, ny, nz)` or `(1, ny, nz, 3)` for v_COM. Write metadata as file attributes.

Datasets:
- `density, pressure, virial_pressure, temperature, avg_speed, avg_O_speed` — `(nx,ny,nz) float32`
- `voxel_type` — `(nx,ny,nz) uint8`, 0=empty, 1=water, 2=silica
- `v_COM` — `(nx,ny,nz,3) float32`

### `process_voxel(arr, masses, V, Si_type, O_type) -> dict`
Inputs: `arr` shape `(N,11)` float64, precomputed `masses` array.

Key computations (all float64 until final cast):
- `r_COM = np.average(arr[:,2:5], axis=0, weights=masses)`
- `v_COM = np.average(arr[:,5:8], axis=0, weights=masses)`
- `r_prime = arr[:,2:5] - r_COM`
- `virial_sum = np.sum(r_prime * arr[:,8:11])`
- `kinetic_term = np.dot(masses, np.sum(arr[:,5:8]**2, axis=1))`
- `v_th = arr[:,5:8] - v_COM`
- `thermal_KE = np.dot(masses, np.sum(v_th**2, axis=1))`
- `T = thermal_KE / (3 * N * k_B)`
- `density = total_mass / V` (V in Å³, convert to g/cm³)
- `voxel_type = 2 if Si_type in types else 1`
- `avg_O_speed`: only if `voxel_type == 1`, else NaN

Returns dict of float32-cast values ready to write.

### `flush_layer(layer_buf, ix, h5file, nx, ny, nz, ...)` 
Iterate over all `(iy, iz)` in current layer, call `process_voxel` or write NaN, write each result into `h5file['density'][ix, iy, iz]` etc.

### `main()`
- Parse args
- Open dump, call `parse_header`
- Open HDF5, call `preallocate_hdf5`
- Stream atom lines, accumulate into `layer_buf = defaultdict(list)`
- On ix change: call `flush_layer`, clear buffer
- After EOF: flush final layer

## Units note
LAMMPS real units: positions in Å, velocities in Å/ps, forces in kcal/mol/Å, masses in amu.
- Pressure: result of `(kinetic + virial) / (3V)` will be in kcal/mol/Å³ — convert to GPa (1 kcal/mol/Å³ ≈ 6.9479 GPa)
- Temperature: `k_B = 8.617e-5 eV/K = 1.9872e-3 kcal/mol/K`
- Density: `total_mass_amu * 1.6605e-27 kg / (V_Å³ * 1e-30 m³)` → kg/m³, or convert amu/Å³ to g/cm³ (factor: 1.6605)
