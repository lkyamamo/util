import pyvista as pv
import h5py
import numpy as np

h5_file = 'trajectory.h5'

with h5py.File(h5_file, 'r') as f:
    data = f['data'][:]

print(data.shape)

pv.plot(data)