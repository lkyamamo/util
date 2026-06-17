import pyvista as pv
import h5py
import numpy as np

h5_file = '/Users/loganyamamoto/Desktop/Research/grants/geo_sciences/finalized-bubble-collapse/analysis/water/voxel_r0128/trajectory.h5'

with h5py.File(h5_file, 'r') as f:
    print(list(f.keys()))
    print(f['density'].shape)
    density = np.array(f['density'])

#frame = density[0]

#test frame
frame = np.array([[[5,5,1],[5,5,5],[5,5,5]],[[5,5,5],[5,5,5],[5,5,5]],[[5,5,5],[5,5,5],[5,5,5]]])

# 2. Initialize a PyVista ImageData object
grid = pv.ImageData()

# 3. Match the grid dimensions to your matrix dimensions
grid.dimensions = frame.shape

flat = frame.flatten(order="F")
print(flat.shape)

# 4. Flatten the matrix with Fortran 'F' ordering to map it to the grid correctly
grid.point_data["My Values"] = frame.flatten(order="F")

# 5. Visualize the data as a solid volume block
grid.plot(volume=True, cmap="viridis")
grid.show(jupyter_backend='client')