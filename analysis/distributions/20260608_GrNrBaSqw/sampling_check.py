"""
sampling_check.py — size a run BEFORE submitting it.

Answers the two questions the parameter tables do not: will this configuration
be converged, and will it finish. Reads the SAME environment variables the
pipeline uses, so a planned configuration is checked without restating it:

    TRAJ / DYNAMICS_TRAJ, R_MAX, RDF_BINS, DSF_Q_MAX, DSF_N_Q_BINS,
    DSF_N_FRAMES, DSF_STRIDE, DYNAMICS_DT, VDOS_CORR_LENGTH,
    VDOS_CORR_INTERVAL, MSD_CORR_LENGTH, MSD_CORR_INTERVAL

    python sampling_check.py [--threads N]

WHAT IT PREDICTS
----------------
Sampling number M — the count of independent contributions behind the noisiest
output of each analysis. Relative statistical error is ~1/sqrt(M):

    M = 1e2   10%    exploratory only
    M = 1e3    3%    PASS MARK
    M = 1e4    1%    publication / comparison against measured data
    M = 1e6  0.1%    diminishing returns

M = (contributors per configuration) x (independent configurations), so it is
size-invariant: double the atoms and you can halve the frames.

Predictions here assume g(r) = 1 at a reference radius, i.e. a structureless
fluid. They are a rough LOWER BOUND, good to about an order of magnitude, not a
precise forecast — a real peak sits at some r_p with height g_p and holds
(r_p/r_ref)^2 * g_p times more, which is often 10-100x.

Measured BELOW predicted is not a bug either: it means g < 1 at that pair's
apparent peak, i.e. the two species avoid each other. Verified on a real
hydroxylated-silica frame, where H-Si comes back under prediction because H and
Si are never bonded — they are separated by the bridging O, so the argmax is
just picking the least-empty bin of a depleted curve. Read a LOW verdict on such
a pair as "this pair has no first shell", not as "add frames".

COST CONSTANTS ARE MEASURED, AND ONLY ROUGHLY TRANSFERABLE
----------------------------------------------------------
Calibrated on a 12-core Apple-silicon mac (analysis/.venv), single-threaded,
against these fixtures:
    rdf   134784-atom frame, 6 partials, R_MAX 4/8/16
    bad   5184 atoms x 10 frames, R_CUTOFF 2/4/6
    dsf   5184 atoms x 10 frames, Q_MAX 5/10
They reproduce the swept points to about +/-20%, but degrade away from them: at
R_MAX=6 the rdf constant (calibrated at 8) under-predicts compute by ~2.5x, and
small cutoffs are worse still because fixed overhead dominates. Treat a
prediction as the right order of magnitude and the SCALING as the reliable part
— the ratio between two candidate settings is much better than either absolute. Treat them as order-of-
magnitude on other hardware; the wall/core-seconds line each analysis now logs
is what lets them be recalibrated. vdos/msd have no absolute constant here —
their scaling law is printed instead of a fabricated time.
"""

import os
import sys
import numpy as np

# --- measured cost constants (see the docstring for provenance) --------------
K_RDF = 2.9e-8    # s per (atom * A^3) at 1 thread:  t ~ K * frames * N * R_MAX^3
K_BAD = 4.8e-7    # s per (atom * (rho R^3)^2):      t ~ K * frames * N * (rho R^3)^2
K_DSF = 1.3e-9    # s per (atom * q-vector):         t ~ K * frames * N * N_q

SAMPLING_TARGET = 1e3
# Measured parallel efficiency for freud (rdf/bad): 66% at 4 threads, 28% at 12.
FREUD_SPEEDUP = {1: 1.00, 2: 1.68, 4: 2.65, 6: 3.05, 8: 2.70, 12: 3.30}


def _env(name, default=None, cast=float):
    raw = os.environ.get(name, "").strip()
    if not raw:
        return default
    try:
        return cast(raw)
    except ValueError:
        return default


def human(seconds):
    """Readable across the six orders of magnitude these predictions span."""
    if seconds < 1:
        return f"{seconds*1000:,.0f} ms"
    if seconds < 120:
        return f"{seconds:,.1f} s"
    if seconds < 7200:
        return f"{seconds/60:,.1f} min"
    return f"{seconds/3600:,.1f} h"


def verdict(m):
    if m <= 0:
        return float('inf'), 'EMPTY'
    return 1.0 / np.sqrt(m), ('ok' if m >= SAMPLING_TARGET else 'LOW')


def read_first_frame(path):
    """Atom count, per-element counts and box volume from frame 0 only."""
    with open(path) as f:
        f.readline(); f.readline(); f.readline()
        n_atoms = int(f.readline().strip())
        f.readline()
        lo_hi = [tuple(map(float, f.readline().split()[:2])) for _ in range(3)]
        header = f.readline().split()[2:]
        col_el = header.index('element')
        counts = {}
        for _ in range(n_atoms):
            el = f.readline().split()[col_el]
            counts[el] = counts.get(el, 0) + 1
        bytes_per_frame = f.tell()
    volume = float(np.prod([hi - lo for lo, hi in lo_hi]))
    lengths = [hi - lo for lo, hi in lo_hi]
    return n_atoms, counts, volume, lengths, bytes_per_frame


def count_frames(path, bytes_per_frame):
    """Exact for small files, estimated from size for large ones; says which."""
    size = os.path.getsize(path)
    if size < 200e6:
        with open(path) as f:
            n = sum(1 for line in f if line.startswith('ITEM: TIMESTEP'))
        return n, 'exact'
    return max(1, int(round(size / bytes_per_frame))), 'estimated from file size'


def predict_rdf(n_atoms, counts, volume, frames, threads):
    r_max = _env("R_MAX", 20.0)
    bins = int(_env("RDF_BINS", 2000, float))
    dr = r_max / bins
    r_ref = 2.0
    print(f"\nrdf_freud.py   R_MAX={r_max} Å, {bins} bins (Δr={dr:.4f} Å), {frames} frames")
    print(f"  M below is for one bin at r = {r_ref:g} Å with g = 1. A real first peak at "
          f"r_p with height g_p\n  holds (r_p/{r_ref:g})^2 x g_p times more — e.g. a peak at "
          f"3.1 Å with g = 4 is ~10x these numbers.")
    els = sorted(counts)
    rows = []
    for i, a in enumerate(els):
        for b in els[i:]:
            # Pairs of A around B in a shell at r, for a structureless fluid.
            # Evaluated at a plausible first-peak radius; the real peak has
            # g ~ 2-4, so this is the conservative end.
            per_frame = counts[a] * (counts[b] / volume) * 4 * np.pi * r_ref ** 2 * dr
            if a == b:
                per_frame *= 0.5
            rows.append((f"{a}-{b}", per_frame * frames))
    width = max(len(k) for k, _ in rows)
    for label, m in sorted(rows, key=lambda t: t[1]):
        err, v = verdict(m)
        print(f"  {label.ljust(width)}  M ~ {m:9.3g}  {100*err:6.2f}%  {v}"
              f"{'   <-- limiting' if (label, m) == min(rows, key=lambda t: t[1]) else ''}")
    t1 = K_RDF * frames * n_atoms * r_max ** 3
    speedup = FREUD_SPEEDUP.get(threads, 3.3)
    print(f"  cost ~ frames x N x R_MAX^3  ->  ~{human(t1/speedup)} at {threads} threads "
          f"({human(t1)} at 1)")
    print(f"        R_MAX^3 dominates: halving R_MAX to {r_max/2:g} would be ~8x cheaper")


def predict_bad(n_atoms, volume, frames, threads):
    cutoff = 2.2   # no single env var; the per-pair map is parsed by bad_freud itself
    rho = n_atoms / volume
    coord = rho * (4 / 3) * np.pi * cutoff ** 3
    triplets = n_atoms * coord * (coord - 1) / 2 * frames
    print(f"\nbad_freud.py   assuming a ~{cutoff} Å cutoff (set per pair in R_CUTOFF)")
    print(f"  ~{coord:.1f} neighbours per atom  ->  {triplets:9.3g} triplets total")
    print(f"  M is the count in the PEAK bin, a few % of that; watch the reported value")
    t1 = K_BAD * frames * n_atoms * (rho * cutoff ** 3) ** 2
    speedup = FREUD_SPEEDUP.get(threads, 3.3)
    print(f"  cost ~ frames x N x (rho R^3)^2  ->  ~{human(t1/speedup)} at {threads} threads")
    print(f"        measured ~R^6: raising the cutoff 4->6 Å cost 16x. This is the most "
          f"expensive knob\n        in the pipeline per Å — set it from the first RDF "
          f"minimum, never generously")


def predict_dsf(n_atoms, volume, lengths, threads):
    q_max = _env("DSF_Q_MAX", 20.0)
    n_bins = int(_env("DSF_N_Q_BINS", 130, float))
    stride = int(_env("DSF_STRIDE", 600, float))
    n_frames_env = _env("DSF_N_FRAMES", None, float)
    frames = int((n_frames_env or 0) / stride) or 1
    n_q = q_max ** 3 * volume / (6 * np.pi ** 2)
    dq_min = 2 * np.pi / max(lengths)
    ceiling = q_max / dq_min
    print(f"\ndsf.py         Q_MAX={q_max} Å⁻¹, {n_bins} q-bins, ~{frames} frames used")
    print(f"  ~{n_q:,.0f} q-vectors total; 2π/L = {dq_min:.4f} Å⁻¹")
    if n_bins > ceiling:
        print(f"  N_Q_BINS={n_bins} EXCEEDS Q_MAX/(2π/L) = {ceiling:.0f} — low-q bins will be empty")
    # Lowest populated shell holds only the handful of lattice vectors at 2*pi/L.
    m_low = 6 * frames
    err, v = verdict(m_low)
    print(f"  lowest-q bin  M ~ {m_low:9.3g}  {100*err:6.2f}%  {v}   <-- always the limiting end")
    t1 = K_DSF * frames * n_atoms * n_q
    print(f"  cost ~ frames x N x Q_MAX^3  ->  ~{human(t1)} at 1 thread")
    print(f"        Q_MAX={q_max:g} -> {q_max*0.6:g} would be ~4.6x cheaper; this is the "
          f"single most\n        expensive default in the pipeline")


def predict_origins(tag, counts, frames, dt):
    corr_len = _env(f"{tag}_CORR_LENGTH")
    corr_int = _env(f"{tag}_CORR_INTERVAL")
    if corr_len is None or corr_int is None or dt is None:
        print(f"\n{tag.lower()}.py        skipped: set DYNAMICS_DT, {tag}_CORR_LENGTH "
              f"and {tag}_CORR_INTERVAL")
        return
    span = frames * dt
    n_origins = max(1, int((span - corr_len) / corr_int) + 1)
    n_indep = max(1, int(span // corr_len))
    print(f"\n{tag.lower()}.py        span {span/1000:.1f} ps, CORR_LENGTH {corr_len/1000:.2f} ps")
    print(f"  {n_origins} origins, of which {n_indep} are independent "
          f"(span / CORR_LENGTH)")
    for el in sorted(counts):
        m_raw = counts[el] * n_origins
        m_ind = counts[el] * n_indep
        err, v = verdict(m_ind)
        print(f"  {el:<4} {counts[el]:6d} atoms   M_raw = {m_raw:9.3g}   "
              f"M_indep = {m_ind:9.3g}  {100*err:6.2f}%  {v}")
    if n_origins > n_indep:
        waste = n_origins / n_indep
        print(f"  -> CORR_INTERVAL gives {waste:.0f}x more origins than independent windows: "
              f"runtime scales\n     with the former, precision only with the latter. "
              f"Raising it to {span/n_indep:.0f} fs costs nothing.")
    print(f"  cost ~ n_origins x CORR_LENGTH x N (no absolute constant measured here)")


if __name__ == '__main__':
    threads = 4
    if '--threads' in sys.argv:
        threads = int(sys.argv[sys.argv.index('--threads') + 1])

    struct = os.environ.get("TRAJ", "").strip()
    dyn = os.environ.get("DYNAMICS_TRAJ", "").strip()
    if not struct and not dyn:
        raise SystemExit(
            "sampling_check.py: set TRAJ (structural) and/or DYNAMICS_TRAJ (dynamics),\n"
            "  the same variables the pipeline uses. Nothing to size otherwise."
        )

    print("=" * 78)
    print(f"Sampling and cost pre-flight — target M >= {SAMPLING_TARGET:.0e} (~3% error), "
          f"{threads} threads")
    print("Predictions assume g(r)=1; a real peak has g~2-4, so the run should report MORE.")
    print("=" * 78)

    if struct and os.path.exists(struct):
        n, counts, vol, lengths, bpf = read_first_frame(struct)
        frames, how = count_frames(struct, bpf)
        print(f"\n{struct}\n  {n} atoms, {frames} frames ({how}), "
              f"box {np.round(lengths,2)}, V = {vol:,.0f} Å³")
        print(f"  composition: " + ", ".join(f"{e}={c}" for e, c in sorted(counts.items())))
        predict_rdf(n, counts, vol, frames, threads)
        predict_bad(n, vol, frames, threads)

    if dyn and os.path.exists(dyn):
        n, counts, vol, lengths, bpf = read_first_frame(dyn)
        frames, how = count_frames(dyn, bpf)
        dt = _env("DYNAMICS_DT")
        print(f"\n{dyn}\n  {n} atoms, {frames} frames ({how}), V = {vol:,.0f} Å³")
        predict_dsf(n, vol, lengths, threads)
        predict_origins("VDOS", counts, frames, dt)
        predict_origins("MSD", counts, frames, dt)

    print("\n" + "=" * 78)
    print("Raise M by adding frames or atoms — both enter linearly. Lower cost by "
          "lowering the\ncubic knobs (R_MAX, Q_MAX) or the sixth-power one (BAD cutoff) "
          "before adding threads:\nfreud gives 2.65x at 4 threads and only 3.30x at 12.")
