import ovito
from ovito.io import import_file, export_file
from ovito.modifiers import AffineTransformationModifier
import numpy as np

pipeline = import_file("input.data", atom_style="full")

# Rotation matrix for +90° about +Y axis
R = np.array([
    [ 0, 0, 1, 0],
    [ 0, 1, 0, 0],
    [-1, 0, 0, 0]
])

pipeline.modifiers.append(
    AffineTransformationModifier(
        transformation=R,
        transform_particles=True,
        transform_velocities=True,   # handles vx vy vz automatically
        transform_box=True,          # rotates and recenters the simulation cell
    )
)

export_file(pipeline, "output.data", format="lammps/data", atom_style="full")