"""
Quick test: does VTK render NaN cells as transparent in a volume?
A 4x4x4 grid where the center 2x2x2 cells have values and the rest are NaN.
If NaN = transparent, you should see a small cube floating in empty space.
If NaN = opaque/colored, you'll see the full 4x4x4 box.
"""
import numpy as np
import pyvista as pv

nx = ny = nz = 4
N = nx * ny * nz

grid = pv.ImageData()
grid.dimensions = (nx+1, ny+1, nz+1)
grid.spacing    = (1.0, 1.0, 1.0)
grid.origin     = (0.0, 0.0, 0.0)

# center 2x2x2 = value 1.0, rest = NaN
data = np.full((nx, ny, nz), np.nan, dtype=np.float32)
data[1:3, 1:3, 1:3] = 1.0

grid.cell_data['values'] = data.flatten(order='F')

p = pv.Plotter()
p.set_background('black')
p.add_volume(grid, scalars='values', opacity='linear', cmap='viridis')
p.add_text('NaN = transparent → small cube visible\nNaN = opaque → full box visible',
           position='upper_left', font_size=10, color='white')
p.show()
