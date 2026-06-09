import math

# =============================================================================
# CONFIGURATION
# =============================================================================

DUMP_FILE         = "../dump.lammpstrj"
OUTPUT_FILE       = "calculation.xyz"
WRITE_BUFFER_SIZE = 1_000_000   # ~1 MB mid-frame flush threshold; 0 = end-of-frame only

start_timestep    = 0           # first timestep to include (inclusive)
end_timestep      = -1          # last timestep to include (inclusive); -1 means all
type_map          = {"1": "O", "2": "H"}

# Fallback box parameters — used only when the input format does not include box bounds
# (e.g. plain XYZ dump).  Set to None to require the format to supply them.
FALLBACK_A     = None
FALLBACK_B     = None
FALLBACK_C     = None
FALLBACK_ALPHA = 90.0
FALLBACK_BETA  = 90.0
FALLBACK_GAMMA = 90.0

# =============================================================================
# END CONFIGURATION
# =============================================================================


# --- per-atom line mapping ---------------------------------------------------

def _map_atom_line_element(parts, col_sym, col_x, col_y, col_z):
    """Used when 'element' column is present — type_map is ignored."""
    return f"{parts[col_sym]} {parts[col_x]} {parts[col_y]} {parts[col_z]}\n"


def _map_atom_line_type(parts, col_sym, col_x, col_y, col_z, tm):
    """Used when only 'type' column is present — applies type_map."""
    element = tm.get(parts[col_sym], parts[col_sym])
    return f"{element} {parts[col_x]} {parts[col_y]} {parts[col_z]}\n"


# --- column detection --------------------------------------------------------

def _detect_columns(atoms_header):
    """Parse an 'ITEM: ATOMS ...' header line.

    Slices away the 'ITEM:' and 'ATOMS' tokens so that the returned indices
    align with actual atom data line splits.  Prefers 'element' over 'type'.
    """
    cols = {name: idx for idx, name in enumerate(atoms_header.split()[2:])}
    col_x = cols['x']
    col_y = cols['y']
    col_z = cols['z']
    if 'element' in cols:
        return ('element', cols['element'], col_x, col_y, col_z)
    return ('type', cols['type'], col_x, col_y, col_z)


# --- box conversion ----------------------------------------------------------

def _box_to_lattice(box, triclinic):
    """Convert raw LAMMPS box bounds to (a, b, c, alpha, beta, gamma).

    box = (xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz)
    xy/xz/yz are 0.0 for orthorhombic boxes.
    """
    xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz = box
    Lx = xhi - xlo
    Ly = yhi - ylo
    Lz = zhi - zlo
    if not triclinic:
        return Lx, Ly, Lz, 90.0, 90.0, 90.0
    a = Lx
    b = math.sqrt(Ly**2 + xy**2)
    c = math.sqrt(Lz**2 + xz**2 + yz**2)
    alpha = math.degrees(math.acos((xy * xz + Ly * yz) / (b * c)))
    beta  = math.degrees(math.acos(xz / c))
    gamma = math.degrees(math.acos(xy / b))
    return a, b, c, alpha, beta, gamma


def _lattice_line(frame):
    """Build the lattice parameter string for the XYZ comment line."""
    box = frame['box']
    if box is None:
        a, b, c = FALLBACK_A, FALLBACK_B, FALLBACK_C
        alpha, beta, gamma = FALLBACK_ALPHA, FALLBACK_BETA, FALLBACK_GAMMA
        if None in (a, b, c):
            raise ValueError(
                "Input format has no box bounds — set FALLBACK_A/B/C in config"
            )
    else:
        a, b, c, alpha, beta, gamma = _box_to_lattice(box, frame['triclinic'])
    return f"{a:.6f} {b:.6f} {c:.6f} {alpha:.6f} {beta:.6f} {gamma:.6f}"


# =============================================================================
# FORMAT REGISTRY
# =============================================================================

class CustomDumpHandler:
    """Standard LAMMPS custom dump.

    First line: 'ITEM: TIMESTEP'
    Sections:   TIMESTEP / NUMBER OF ATOMS / BOX BOUNDS / ATOMS
    """

    def detect(self, first_line):
        return first_line.strip() == 'ITEM: TIMESTEP'

    def probe(self, f):
        """Read first frame header to determine triclinic flag and column layout.
        Resets file to position 0 before returning.
        """
        for _ in range(4):           # ITEM: TIMESTEP / value / ITEM: NUMBER OF ATOMS / value
            f.readline()
        bounds_header = f.readline()
        triclinic = 'xy' in bounds_header
        n_box_lines = 3
        for _ in range(n_box_lines):
            f.readline()
        atoms_header = f.readline()
        mode, col_sym, col_x, col_y, col_z = _detect_columns(atoms_header)
        f.seek(0)
        return {
            'triclinic':   triclinic,
            'n_box_lines': n_box_lines,
            'has_box':     True,
            'mode':        mode,
            'col_sym':     col_sym,
            'col_x':       col_x,
            'col_y':       col_y,
            'col_z':       col_z,
        }

    def read_frame_id(self, f):
        """Stage 1: read ITEM:TIMESTEP block + ITEM:NUMBER OF ATOMS block.
        Returns None at EOF.
        """
        line = f.readline()
        if not line:
            return None                   # EOF
        # line == "ITEM: TIMESTEP\n"
        timestep = int(f.readline())
        f.readline()                      # "ITEM: NUMBER OF ATOMS"
        n_atoms = int(f.readline())
        return {'timestep': timestep, 'n_atoms': n_atoms}

    def parse_header_rest(self, f, probe_data):
        """Stage 2: parse box bounds; consume ATOMS header.
        f advances to the first atom data line.
        """
        f.readline()   # "ITEM: BOX BOUNDS ..."
        if probe_data['triclinic']:
            r0 = f.readline().split()
            r1 = f.readline().split()
            r2 = f.readline().split()
            xlo, xhi, xy = float(r0[0]), float(r0[1]), float(r0[2])
            ylo, yhi, xz = float(r1[0]), float(r1[1]), float(r1[2])
            zlo, zhi, yz = float(r2[0]), float(r2[1]), float(r2[2])
            box = (xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz)
        else:
            r0 = f.readline().split()
            r1 = f.readline().split()
            r2 = f.readline().split()
            xlo, xhi = float(r0[0]), float(r0[1])
            ylo, yhi = float(r1[0]), float(r1[1])
            zlo, zhi = float(r2[0]), float(r2[1])
            box = (xlo, xhi, ylo, yhi, zlo, zhi, 0.0, 0.0, 0.0)
        f.readline()   # "ITEM: ATOMS ..."  (columns already detected in probe)
        return {'box': box}

    def skip_rest(self, f, n_atoms, probe_data):
        """Skip remaining header lines (box + atoms header) + all atom lines."""
        n_skip = 1 + probe_data['n_box_lines'] + 1 + n_atoms
        for _ in range(n_skip):
            f.readline()


class XyzDumpHandler:
    """Plain XYZ dump (LAMMPS xyz style).

    First line: bare integer (n_atoms)
    Second line: comment (may contain 'Timestep: N')
    Remaining lines: element x y z ...
    No box bounds in file.
    """

    def detect(self, first_line):
        try:
            int(first_line.strip())
            return True
        except ValueError:
            return False

    def probe(self, f):
        f.seek(0)
        return {
            'has_box': False,
            'mode':    'element',
            'col_sym': 0,
            'col_x':   1,
            'col_y':   2,
            'col_z':   3,
        }

    def read_frame_id(self, f):
        """Read n_atoms line + comment line.
        Tries to extract timestep from 'Atoms. Timestep: N' pattern.
        """
        line = f.readline()
        if not line:
            return None          # EOF
        n_atoms  = int(line.strip())
        comment  = f.readline()
        timestep = None
        if 'Timestep:' in comment:
            try:
                timestep = int(comment.split('Timestep:')[1].strip().split()[0])
            except (ValueError, IndexError):
                pass
        return {'timestep': timestep, 'n_atoms': n_atoms}

    def parse_header_rest(self, f, probe_data):
        """Comment already consumed in read_frame_id; f is at first atom line."""
        return {'box': None}

    def skip_rest(self, f, n_atoms, probe_data):
        for _ in range(n_atoms):
            f.readline()


HANDLERS = [CustomDumpHandler(), XyzDumpHandler()]


# =============================================================================
# FRAME GENERATOR
# =============================================================================

def _parse_frames(filename, start_ts, end_ts):
    """Yield one dict per in-range frame with file handle positioned at first atom line.

    Two-stage approach:
      Stage 1 (always): read_frame_id — timestep + n_atoms only
      Stage 2 (in-range only): parse_header_rest — box bounds + advance to atom lines
    Out-of-range frames are skipped via skip_rest with raw readline() calls only.
    """
    with open(filename) as f:
        first_line = f.readline().strip()
        f.seek(0)

        handler = next((h for h in HANDLERS if h.detect(first_line)), None)
        if handler is None:
            raise ValueError(f"Unrecognised dump format (first line: {first_line!r})")

        probe_data = handler.probe(f)   # reads ahead then resets to 0

        frame_idx = 0
        while True:
            # --- stage 1: always — minimum read to decide skip vs parse ---
            frame_id = handler.read_frame_id(f)
            if frame_id is None:
                break                   # EOF

            ts      = frame_id['timestep'] if frame_id['timestep'] is not None else frame_idx
            n_atoms = frame_id['n_atoms']

            if end_ts != -1 and ts > end_ts:
                break                   # past the requested range; stop immediately

            if ts < start_ts:
                handler.skip_rest(f, n_atoms, probe_data)
                frame_idx += 1
                continue

            # --- stage 2: in-range — parse box bounds, advance to atom lines ---
            rest = handler.parse_header_rest(f, probe_data)

            yield {
                'timestep': ts,
                'n_atoms':  n_atoms,
                'box':      rest['box'],
                'triclinic': probe_data.get('triclinic', False),
                'mode':     probe_data['mode'],
                'col_sym':  probe_data['col_sym'],
                'col_x':    probe_data['col_x'],
                'col_y':    probe_data['col_y'],
                'col_z':    probe_data['col_z'],
                'f':        f,
            }
            frame_idx += 1


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    # Peek at the first in-range frame to build map_line once before any I/O.
    first_frame = next(_parse_frames(DUMP_FILE, start_timestep, end_timestep), None)
    if first_frame is None:
        raise ValueError("No frames found in the specified timestep range")

    mode    = first_frame['mode']
    col_sym = first_frame['col_sym']
    col_x   = first_frame['col_x']
    col_y   = first_frame['col_y']
    col_z   = first_frame['col_z']

    if mode == 'type' and not type_map:
        raise ValueError(
            "type_map is empty but dump uses numeric atom types — set type_map in config"
        )

    if mode == 'element':
        map_line = lambda parts: _map_atom_line_element(parts, col_sym, col_x, col_y, col_z)
    else:
        map_line = lambda parts: _map_atom_line_type(parts, col_sym, col_x, col_y, col_z, type_map)

    buf      = []
    buf_size = 0

    with open(OUTPUT_FILE, 'w') as out:
        for frame in _parse_frames(DUMP_FILE, start_timestep, end_timestep):
            print(f"Writing timestep {frame['timestep']} ({frame['n_atoms']} atoms)")

            # frame header is never buffered
            out.write(f"{frame['n_atoms']}\n{_lattice_line(frame)}\n")

            # pipeline: read one atom line → map → buffer
            f = frame['f']
            for _ in range(frame['n_atoms']):
                atom_line  = map_line(f.readline().split())
                buf.append(atom_line)
                buf_size  += len(atom_line)
                if WRITE_BUFFER_SIZE > 0 and buf_size >= WRITE_BUFFER_SIZE:
                    out.write(''.join(buf))
                    buf.clear()
                    buf_size = 0

            # unconditional end-of-frame flush
            if buf:
                out.write(''.join(buf))
                buf.clear()
                buf_size = 0

    print(f"Done. Output written to {OUTPUT_FILE}")
