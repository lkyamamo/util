"""
Generate an XDMF descriptor for a shock voxel HDF5 trajectory file.

Usage:
    python generate_xdmf.py trajectory.h5 [output.xdmf]

The resulting .xdmf file can be opened directly in ParaView. Keep it in the
same directory as the .h5 file, or use absolute paths.
"""

import os
import sys
import textwrap
import warnings
import h5py

QUANTITIES = [
    'density', 'pressure', 'virial_pressure',
    'temperature', 'avg_speed', 'avg_O_speed', 'voxel_type',
]


def read_meta(h5file):
    attrs = dict(h5file.attrs) if h5file.attrs else {}
    T, nx, ny, nz = h5file[QUANTITIES[0]].shape

    if 'voxel_size' not in attrs:
        warnings.warn(
            "HDF5 file has no 'voxel_size' attribute; defaulting to 10.0 Å. "
            "Box bounds may be incorrect.",
        )
    voxel_size = float(attrs.get('voxel_size', 10.0))
    xlo = float(attrs.get('xlo', 0.0))
    ylo = float(attrs.get('ylo', 0.0))
    zlo = float(attrs.get('zlo', 0.0))
    xhi = float(attrs.get('xhi', xlo + nx * voxel_size))
    yhi = float(attrs.get('yhi', ylo + ny * voxel_size))
    zhi = float(attrs.get('zhi', zlo + nz * voxel_size))
    dx  = (xhi - xlo) / nx
    dy  = (yhi - ylo) / ny
    dz  = (zhi - zlo) / nz

    # Per-timestep step labels (optional)
    if 'timesteps' in h5file:
        timestep_labels = h5file['timesteps'][:]
    elif 'timestep_labels' in h5file:
        timestep_labels = h5file['timestep_labels'][:]
    else:
        timestep_labels = None

    return {
        'T': T, 'nx': nx, 'ny': ny, 'nz': nz,
        'xlo': xlo, 'ylo': ylo, 'zlo': zlo,
        'dx': dx, 'dy': dy, 'dz': dz,
        'timestep_labels': timestep_labels,
    }


def _h5_dtype_to_xdmf(dtype):
    """Map numpy dtype to XDMF NumberType + Precision strings."""
    import numpy as np
    dtype = np.dtype(dtype)
    if np.issubdtype(dtype, np.floating):
        return 'Float', str(dtype.itemsize)
    if np.issubdtype(dtype, np.signedinteger):
        return 'Int', str(dtype.itemsize)
    if np.issubdtype(dtype, np.unsignedinteger):
        return 'UInt', str(dtype.itemsize)
    return 'Float', '8'


def build_xdmf(h5_path, meta):
    T  = meta['T']
    nx, ny, nz = meta['nx'], meta['ny'], meta['nz']
    xlo, ylo, zlo = meta['xlo'], meta['ylo'], meta['zlo']
    dx, dy, dz = meta['dx'], meta['dy'], meta['dz']
    labels = meta['timestep_labels']

    # Detect which quantities actually exist and their dtypes
    with h5py.File(h5_path, 'r') as f:
        present = [(q, f[q].dtype) for q in QUANTITIES if q in f]

    # XDMF3 with 3DCoRectMesh uses ZYX dimension ordering throughout.
    # The H5 datasets are stored as (T, nx, ny, nz) — C-order, so when
    # XDMF reads the flat HyperSlab slice it gets nx*ny*nz values which
    # it then maps into the ZYX cell grid. The outer Dimensions on the
    # Attribute and HyperSlab must be in ZYX order to match the topology.
    timestep_blocks = []
    for t in range(T):
        time_val = int(labels[t]) if labels is not None else t

        attr_blocks = []
        for qty, dtype in present:
            num_type, precision = _h5_dtype_to_xdmf(dtype)
            # HyperSlab rows in XDMF3: start | stride | count (one row each)
            attr_blocks.append(f"""\
      <Attribute Name="{qty}" AttributeType="Scalar" Center="Cell">
        <DataItem ItemType="HyperSlab" Dimensions="{nz} {ny} {nx}">
          <DataItem Dimensions="3 4" Format="XML">
            {t} 0 0 0
            1  1 1 1
            1  {nx} {ny} {nz}
          </DataItem>
          <DataItem Format="HDF" NumberType="{num_type}" Precision="{precision}"
                    Dimensions="{T} {nx} {ny} {nz}">
            {h5_path}:/{qty}
          </DataItem>
        </DataItem>
      </Attribute>""")

        timestep_blocks.append(f"""\
    <Grid Name="t{t}" GridType="Uniform">
      <Time Value="{time_val}"/>
      <Topology TopologyType="3DCoRectMesh" Dimensions="{nz+1} {ny+1} {nx+1}"/>
      <Geometry Type="ORIGIN_DXDYDZ">
        <DataItem Dimensions="3" Format="XML">{xlo} {ylo} {zlo}</DataItem>
        <DataItem Dimensions="3" Format="XML">{dx} {dy} {dz}</DataItem>
      </Geometry>
{"".join(attr_blocks)}
    </Grid>""")

    return textwrap.dedent(f"""\
<?xml version="1.0" ?>
<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="3.0">
  <Domain>
    <Grid Name="ShockTrajectory" GridType="Collection" CollectionType="Temporal">
{"".join(timestep_blocks)}
    </Grid>
  </Domain>
</Xdmf>
""")


def write_clean_h5(src_path, dst_path):
    """Copy H5, replacing NaN with 0.0 in all float datasets so ParaView
    can compute a valid color range without manual rescaling."""
    import numpy as np
    print(f"Writing NaN-cleaned H5 → {dst_path}")
    with h5py.File(src_path, 'r') as src, h5py.File(dst_path, 'w') as dst:
        for k, v in src.attrs.items():
            dst.attrs[k] = v
        for key in src.keys():
            data = src[key][:]
            if np.issubdtype(data.dtype, np.floating):
                data = np.nan_to_num(data, nan=0.0)
            dst.create_dataset(key, data=data, compression='gzip', compression_opts=4)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('h5_file', help='Input trajectory HDF5 file')
    parser.add_argument('xdmf_file', nargs='?', help='Output XDMF path (default: alongside H5)')
    parser.add_argument('--clean', action='store_true',
                        help='Write a NaN-replaced H5 copy and point the XDMF at it')
    args = parser.parse_args()

    h5_path  = os.path.abspath(args.h5_file)
    h5_dir   = os.path.dirname(h5_path)
    h5_stem  = os.path.splitext(os.path.basename(h5_path))[0]
    xdmf_path = os.path.abspath(args.xdmf_file) if args.xdmf_file \
                else os.path.join(h5_dir, h5_stem + '.xdmf')

    if args.clean:
        clean_h5_path = os.path.join(h5_dir, h5_stem + '_clean.h5')
        write_clean_h5(h5_path, clean_h5_path)
        h5_path = clean_h5_path  # XDMF will point at the cleaned copy

    with h5py.File(h5_path, 'r') as f:
        meta = read_meta(f)

    print(f"Trajectory: {meta['T']} timesteps, "
          f"grid {meta['nx']}×{meta['ny']}×{meta['nz']}")
    print(f"Writing XDMF → {xdmf_path}")

    xdmf = build_xdmf(h5_path, meta)
    with open(xdmf_path, 'w') as f:
        f.write(xdmf)

    print("Done. Open the .xdmf file in ParaView.")
    print()
    print("Suggested ParaView workflow:")
    print("  1. File → Open → trajectory.xdmf")
    print("  2. Select 'Xdmf3ReaderT' when prompted")
    print("  3. Apply")
    print("  4. Representation → Volume")
    print("  5. Coloring → choose quantity")
    print("  6. Filters → Clip  (for cut plane)")
    print("  7. File → Save Animation  (for video export)")


if __name__ == "__main__":
    main()
