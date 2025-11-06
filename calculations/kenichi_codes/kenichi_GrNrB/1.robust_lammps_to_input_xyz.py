from pathlib import Path

a = b = c = 37.2514
alpha = beta = gamma = 90.0
start_frame = 0
end_frame = -1
skip_header = 2
type_map = {"1": "O", "2": "H"}

# sanity checks
if any(x == -1 for x in (a, b, c, alpha, beta, gamma)):
    raise ValueError("Must set lattice parameters (magnitudes and angles)")
if "CHANGE1" in type_map.values():
    raise ValueError("Must change type dictionary to element names")

filename = "../all_lammps_last500.xyz"
with open(filename, "r") as f:
    lines = f.readlines()

try:
    num_atoms = int(lines[0].strip())
except Exception as e:
    raise ValueError(f"First line must be integer atom count, got {lines[0]!r}") from e

start_line = start_frame*(num_atoms + skip_header) + skip_header

if end_frame == -1:
    end_frame = len(lines)
end_line = end_frame*(num_atoms + skip_header) + skip_header + 1    # to include last frame

box_dim_string = f"{a:.6f} {b:.6f} {c:.6f} {alpha:.6f} {beta:.6f} {gamma:.6f}\n"

with open("calculation.xyz", "w") as out:
    for j in range(start_line, end_line, num_atoms + skip_header):
        frame = lines[j:j+num_atoms]
        if len(frame) < num_atoms:
            break  # ignore a truncated last frame
        out.write(f"{num_atoms}\n")
        out.write(box_dim_string)
        for ln in frame:
            if not ln.strip():
                continue
            parts = ln.split()  # whitespace-delimited
            typ = parts[0]
            sym = type_map.get(typ, typ)  # leave as-is if not in map
            out.write(sym + " " + " ".join(parts[1:]) + "\n")
