import sys
import glob
import re
import numpy as np
import h5py

def merge(output_dir, final_h5):
    files = sorted(glob.glob(f"{output_dir}/output_*.h5"),
                   key=lambda f: int(re.search(r'output_(\d+)\.h5$', f).group(1)))
    if not files:
        raise FileNotFoundError(f"No output_*.h5 files found in {output_dir}")

    T = len(files)
    print(f"Merging {T} files -> {final_h5}")

    with h5py.File(final_h5, 'w') as out:
        for i, path in enumerate(files):
            with h5py.File(path, 'r') as src:
                for name in src:
                    data = src[name][:]
                    if i == 0:
                        shape = (T,) + data.shape
                        out.create_dataset(name, shape=shape, dtype=data.dtype)
                        for key, val in src.attrs.items():
                            out.attrs[key] = val
                    out[name][i] = data
            if (i + 1) % 50 == 0 or (i + 1) == T:
                print(f"  {i + 1}/{T}")

    print("Done.")

if __name__ == "__main__":
    output_dir = sys.argv[1]
    final_h5   = sys.argv[2]
    merge(output_dir, final_h5)
