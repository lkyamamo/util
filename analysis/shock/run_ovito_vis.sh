#!/bin/bash

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Directory containing LAMMPS dump files
DUMP_DIR="/path/to/dump/dir"

# Glob pattern matching dump files — must contain a numeric timestep
# e.g. "full.*.lammpstrj" for files like full.0.lammpstrj, full.4000.lammpstrj
DUMP_GLOB="full.*.lammpstrj"

# Merged trajectory h5 produced by merge_h5.py
TRAJECTORY_H5="$SCRIPT_DIR/trajectory.h5"
# ---------------------------------------------------------------------------

module load ovito/3.11.2

ovito --script "$SCRIPT_DIR/ovito_vis.py" -- "$DUMP_DIR" "$DUMP_GLOB" "$TRAJECTORY_H5"
