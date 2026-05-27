import math

start_frame = 0
end_frame = -1

input_filename = "traj.dump"
output_filename = "all.xyz"

# Used only when the dump lists "type" instead of "element" in ITEM: ATOMS.
type_map = {"1": "O", "2": "H"}


def box_bounds_to_lattice(line1, line2, line3):
    xlo, xhi, *x_tilt = line1.split()
    ylo, yhi, *y_tilt = line2.split()
    zlo, zhi, *z_tilt = line3.split()

    xy = float(x_tilt[0]) if x_tilt else 0.0
    xz = float(y_tilt[0]) if y_tilt else 0.0
    yz = float(z_tilt[0]) if z_tilt else 0.0

    lx = float(xhi) - float(xlo)
    ly = float(yhi) - float(ylo)
    lz = float(zhi) - float(zlo)

    a = lx
    b = math.sqrt(ly * ly + xy * xy)
    c = math.sqrt(lz * lz + xz * xz + yz * yz)

    if b == 0.0 or c == 0.0:
        return a, b, c, 90.0, 90.0, 90.0

    gamma = math.degrees(math.acos(max(-1.0, min(1.0, xy / b))))
    beta = math.degrees(math.acos(max(-1.0, min(1.0, xz / c))))
    alpha = math.degrees(math.acos(max(-1.0, min(1.0, (xy * xz + ly * yz) / (b * c)))))

    return a, b, c, alpha, beta, gamma


def format_lattice(a, b, c, alpha, beta, gamma):
    return f"{a:.6f} {b:.6f} {c:.6f} {alpha:.6f} {beta:.6f} {gamma:.6f}\n"


def parse_atoms_header(header):
    if not header.startswith("ITEM: ATOMS"):
        raise ValueError(f"Expected 'ITEM: ATOMS', got {header!r}")

    columns = header.split()[2:]
    required = ("x", "y", "z", "vx", "vy", "vz")
    missing = [name for name in required if name not in columns]
    if missing:
        raise ValueError(f"Dump must include {missing} in ITEM: ATOMS line, got {columns}")

    if "element" in columns:
        symbol_col = "element"
    elif "type" in columns:
        symbol_col = "type"
    else:
        raise ValueError("Dump must include 'element' or 'type' in ITEM: ATOMS line")

    col_index = {name: columns.index(name) for name in columns}
    return symbol_col, col_index


def format_atom_line(parts, symbol_col, col_index):
    symbol = parts[col_index[symbol_col]]
    if symbol_col == "type":
        symbol = type_map.get(symbol, symbol)

    values = [parts[col_index[name]] for name in ("x", "y", "z", "vx", "vy", "vz")]
    return f"{symbol} {' '.join(values)}\n"


def read_frame(f):
    header = f.readline()
    if header == "":
        return None
    if not header.startswith("ITEM: TIMESTEP"):
        raise ValueError(f"Expected 'ITEM: TIMESTEP', got {header!r}")

    timestep = int(f.readline().strip())

    if f.readline().strip() != "ITEM: NUMBER OF ATOMS":
        raise ValueError("Expected 'ITEM: NUMBER OF ATOMS'")

    natoms = int(f.readline().strip())

    bounds_header = f.readline().strip()
    if not bounds_header.startswith("ITEM: BOX BOUNDS"):
        raise ValueError(f"Expected 'ITEM: BOX BOUNDS', got {bounds_header!r}")

    box_lines = [f.readline(), f.readline(), f.readline()]
    if any(line == "" for line in box_lines):
        raise ValueError("Unexpected EOF while reading box bounds")

    atoms_header = f.readline().strip()
    symbol_col, col_index = parse_atoms_header(atoms_header)

    atom_lines = []
    for _ in range(natoms):
        line = f.readline()
        if line == "":
            raise ValueError(f"Unexpected EOF while reading atoms (expected {natoms})")
        atom_lines.append(line)

    lattice = box_bounds_to_lattice(box_lines[0], box_lines[1], box_lines[2])
    return timestep, natoms, lattice, symbol_col, col_index, atom_lines


def main():
    frames_written = 0

    with open(input_filename, "r") as infile, open(output_filename, "w") as outfile:
        frame_idx = 0
        while True:
            frame = read_frame(infile)
            if frame is None:
                break

            timestep, natoms, lattice, symbol_col, col_index, atom_lines = frame

            if frame_idx >= start_frame and (end_frame == -1 or frame_idx < end_frame):
                outfile.write(f"{natoms}\n")
                outfile.write(format_lattice(*lattice))
                for line in atom_lines:
                    parts = line.split()
                    outfile.write(format_atom_line(parts, symbol_col, col_index))
                frames_written += 1

            frame_idx += 1

    print(f"Wrote {frames_written} frames to {output_filename}")


if __name__ == "__main__":
    main()
