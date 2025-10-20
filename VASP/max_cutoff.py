import sys
import numpy as np

def parse_vectors_from_text(text):
    """
    Parse 3 lines of 'x y z' into numpy vectors a, b, c.
    Accepts whitespace or comma separators; ignores blank lines.
    """
    lines = [ln.strip() for ln in text.strip().splitlines() if ln.strip()]
    if len(lines) < 3:
        raise ValueError("Need at least 3 non-empty lines with 'x y z' each.")
    def parse_line(ln):
        parts = ln.replace(",", " ").split()
        if len(parts) < 3:
            raise ValueError(f"Line has fewer than 3 numbers: {ln!r}")
        return np.array([float(parts[0]), float(parts[1]), float(parts[2])], dtype=float)
    a = parse_line(lines[0])
    b = parse_line(lines[1])
    c = parse_line(lines[2])
    return a, b, c

def inradius_parallelepiped(a, b, c):
    """
    Return (inradius, details) for the parallelepiped spanned by a,b,c.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    c = np.asarray(c, dtype=float)

    vol = abs(np.dot(a, np.cross(b, c)))
    if vol <= 0:
        raise ValueError("Degenerate cell (zero volume).")

    area_bc = np.linalg.norm(np.cross(b, c))
    area_ca = np.linalg.norm(np.cross(c, a))
    area_ab = np.linalg.norm(np.cross(a, b))
    if min(area_bc, area_ca, area_ab) <= 0:
        raise ValueError("Degenerate face area detected.")

    h_a = vol / area_bc  # separation of faces spanned by b,c
    h_b = vol / area_ca  # separation of faces spanned by c,a
    h_c = vol / area_ab  # separation of faces spanned by a,b
    r = 0.5 * min(h_a, h_b, h_c)

    return r, {
        "volume": vol,
        "h_a": h_a, "h_b": h_b, "h_c": h_c,
        "center": 0.5 * (a + b + c),  # center of insphere if origin is at 0
    }

if __name__ == "__main__":
    # Usage:
    #   python script.py path/to/file_with_3_lines.txt
    #   or: echo -e "x1 y1 z1\nx2 y2 z2\nx3 y3 z3" | python script.py
    if len(sys.argv) == 2 and sys.argv[1] != "-":
        with open(sys.argv[1], "r") as f:
            text = f.read()
    else:
        text = sys.stdin.read()

    a, b, c = parse_vectors_from_text(text)
    r, info = inradius_parallelepiped(a, b, c)
    print("a:", a)
    print("b:", b)
    print("c:", c)
    print(f"Inradius: {r:.10f}")
    print(f"Face separations: h_a={info['h_a']:.10f}, h_b={info['h_b']:.10f}, h_c={info['h_c']:.10f}")
    print("Center (relative to origin):", info["center"])