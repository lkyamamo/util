# Performance Notes

## Current implementation

The Python cleanup step (`cleanup_orphan_atoms.py`) has three memory-bound operations:

| Operation | Implementation | Memory | Time (12M atoms) |
|---|---|---|---|
| `count_atoms_by_type` | Streaming, line-by-line | O(1) | ~10–20 s per file |
| `_parse_atoms_section` | `list[_Atom]` with `__slots__` | ~1 GB | ~30–60 s |
| KDTree on all O or H | `scipy.spatial.KDTree` on full type set | ~200–300 MB | ~1–3 s |
| `_rewrite_data_file` | `path.read_text()` loads entire file | ~1.5 GB (briefly doubled) | ~30–60 s |

**Estimated total for 12M atoms: 2–5 minutes, ~2–4 GB RAM peak.** Acceptable for a one-time setup job.

---

## Hard limits

At 1 billion atoms the current implementation fails at three points:

- `path.read_text()` in `_rewrite_data_file` — ~125 GB file loaded into a Python string
- `list[_Atom]` in `_parse_atoms_section` — ~100 GB of Python objects
- KDTree on all O atoms (~300M points) — ~24 GB numpy array, prohibitively slow to build

---

## High-leverage improvements (ordered by impact)

### 1. Stream `_rewrite_data_file` instead of `read_text()`

Replace the full-file load with a line-by-line streaming rewrite to a temp file.
Eliminates the largest single memory spike and makes the rewrite scale to any file size.

```python
# current (bad for large files)
original_text = path.read_text()
lines = original_text.splitlines(keepends=True)

# replacement: stream input → write output line by line
with open(data_file) as fin, tempfile.NamedTemporaryFile(...) as fout:
    for line in fin:
        # process and write line-by-line
```

### 2. Replace `list[_Atom]` with numpy structured arrays in `_parse_atoms_section`

Storing all atom data in contiguous numpy arrays instead of Python objects reduces
memory by ~10× and speeds up subsequent numpy/KDTree operations.

```python
# pre-allocate or use np.fromiter / np.loadtxt approach
ids    = np.empty(n_atoms, dtype=np.int64)
types  = np.empty(n_atoms, dtype=np.int32)
coords = np.empty((n_atoms, 3), dtype=np.float64)
```

Memory: ~200 MB for 12M atoms vs ~1 GB with Python objects.

### 3. Spatially pre-filter the KDTree partner set

For the nearest-O (or nearest-H) query, only partner atoms within
`radius + shell_thickness + bond_cutoff` of the bubble center can ever be the
nearest partner to a shell atom. Filtering the KDTree input to this spatial
region shrinks it from millions of points to thousands.

```python
partner_cutoff = radius + shell_thickness + 5.0  # Å, conservative bond buffer
partner_dists  = np.linalg.norm(partner_coords - center, axis=1)
partner_coords = partner_coords[partner_dists <= partner_cutoff]
```

---

## For 1 billion atoms — fundamentally different approach required

File I/O alone becomes a 50–125 GB problem that cannot be held in memory.
Options at that scale:

- **OVITO Python API** — C++ backend, can stream and spatially filter without loading
  everything into Python memory.
- **MDAnalysis** — designed for large trajectory/topology files, has spatial selection.
- **LAMMPS binary restart files** — 3–5× smaller than text data files and faster to parse;
  avoids the text I/O bottleneck entirely.
- **Chunked streaming parser** — read and write the data file in fixed-size line batches,
  never holding more than one chunk in memory at once.
