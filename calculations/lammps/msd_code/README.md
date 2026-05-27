# MSD input converters

Scripts here convert MD trajectories into `all.xyz`, the input format read by `calculations/codes/msd.cpp`.

## Output format (`all.xyz`)

Each frame:

```
<num_atoms>
<a> <b> <c> <alpha> <beta> <gamma>
<element> <x> <y> <z> <vx> <vy> <vz>
...
```

Frames are concatenated in one file. Both scripts skip the first `10000` frames by default.

---

## `lammps_to_msd_input_xyz.py`

**Input:** `all_lammps.xyz` — lab-specific concatenated trajectory (not standard LAMMPS dump).

Per frame: 9-line header, then atom lines with LAMMPS **type ID** as the first character.

```
(line 4: num_atoms)
(box bounds, blank lines, etc.)
<type> ...remaining columns...
```

**Configure** at top of script: `a/b/c/alpha/beta/gamma`, `start_timestep` (frame index), `type_dict`, `filename`.

```bash
python lammps_to_msd_input_xyz.py
```

---

## `lammps_custom_to_msd_input_xyz.py`

**Input:** standard LAMMPS **custom** dump with velocities and element names.

```
ITEM: TIMESTEP
<step>
ITEM: NUMBER OF ATOMS
<n>
ITEM: BOX BOUNDS ...
<3 box lines>
ITEM: ATOMS element x y z vx vy vz
<element> <x> <y> <z> <vx> <vy> <vz>
...
```

**Configure** at top of script: `input_filename`, `output_filename`, `start_frame`, `end_frame`, `type_map` (only if dump uses `type` instead of `element`).

```bash
python lammps_custom_to_msd_input_xyz.py
```

**LAMMPS dump command:**

```lammps
dump 1 all custom 100 traj.dump element x y z vx vy vz
dump_modify 1 element O H
```

Lattice parameters are taken from `ITEM: BOX BOUNDS` each frame.
