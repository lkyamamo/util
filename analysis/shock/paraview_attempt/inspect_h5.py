"""
Quick diagnostic: print dataset shapes, dtypes, and attributes for a shock H5 file.

Usage:
    python inspect_h5.py trajectory.h5
"""
import sys
import h5py

def main():
    if len(sys.argv) < 2:
        print("Usage: python inspect_h5.py trajectory.h5")
        sys.exit(1)

    path = sys.argv[1]
    with h5py.File(path, 'r') as f:
        print("=== Attributes ===")
        for k, v in f.attrs.items():
            print(f"  {k}: {v}")

        print("\n=== Datasets ===")
        for key in f.keys():
            ds = f[key]
            print(f"  {key:25s}  shape={ds.shape}  dtype={ds.dtype}")

if __name__ == "__main__":
    main()
