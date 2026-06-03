import ovito
from ovito.io import import_file, export_file
from ovito.modifiers import AffineTransformationModifier
import numpy as np

pipeline = import_file("test.data", atom_style="atomic")

# Rotation matrix for +90° about +Y axis
R = np.array([
    [ 0, 0, 1, 0],
    [ 0, 1, 0, 0],
    [-1, 0, 0, 0]
])

pipeline.modifiers.append(
    AffineTransformationModifier(
        transformation=R,
        operate_on={'particles', 'vector_properties', 'cell'},
    )
)

def shift_to_positive(frame, data):
    """Shift positions so the bounding box minimum corner is at the origin.
    Cell vectors whose primary axis component is negative are flipped so the
    bounding box has only positive extents."""
    cell = data.cell.matrix          # 3x4: columns 0-2 are cell vectors, column 3 is origin
    origin = cell[:, 3]
    a, b, c = cell[:, 0].copy(), cell[:, 1].copy(), cell[:, 2].copy()

    corners = np.array([
        origin, origin+a, origin+b, origin+c,
        origin+a+b, origin+a+c, origin+b+c, origin+a+b+c,
    ])
    lo = corners.min(axis=0)

    # Shift positions so bounding box minimum is at the origin.
    # Access the property as a numpy view (no copy) and subtract in-place.
    pos = data.particles_['Position_']
    pos[:] -= lo

    # Flip any cell vector whose primary-axis component is negative
    for v in (a, b, c):
        primary = np.argmax(np.abs(v))
        if v[primary] < 0:
            v *= -1

    data.cell_.matrix = np.column_stack([a, b, c, np.zeros(3)])

pipeline.modifiers.append(shift_to_positive)

export_file(pipeline, "test_rotated.data", format="lammps/data", atom_style="atomic")