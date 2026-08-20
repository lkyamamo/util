#!/usr/bin/env python3
"""
generate_cascade_input.py — emit the stage-1 LAMMPS input for the
temperature-sweep pipeline.

One LAMMPS run does all of the thermalization: start from ice, deform to the
box size for a target density, ramp to the highest requested temperature, then
descend through the sorted temperature list. At each stop it holds, writes
thermalized_T<C>.data, and — only at the temperatures where diffusion was
asked for — runs an NVE production block writing dynamics_T<C>.lammpstrj.

Why generate the input instead of driving a static template with -var, which
is what every other pipeline here does: the block STRUCTURE varies. There are
N temperature stops, each optionally carrying an NVE block, and LAMMPS index
variables + "jump SELF" express that only awkwardly (tracking the previous
stop's temperature across iterations is the ugly part). The generated file
still lands in <sweep_dir>/input_files/in.input, so it is readable, diffable,
and rerunnable by hand exactly like a copied template would be. Precedent:
md_setup/create_bubble/generate_lammps_script.py.

Usage (normally called by submit_temperature_sweep.sh, but standalone-friendly):

    python3 generate_cascade_input.py \
        --start-data starting-structures/ICE_CUBIC.data \
        --replicate "6 6 6" --density 1.0 \
        --diffusion-temps-c "45;5" --dielectric-temps-c "45;65" \
        --out in.input

DEPENDENCIES
------------
Imports analysis/general/{box_size,box_size_from_data}.py for the density
solve, which pulls in ase. Same venv the analysis scripts already run in.
"""

import argparse
import sys
from pathlib import Path

# Where this file sits in the checkout, used only as the default --repo-root
# for standalone runs. submit_temperature_sweep.sh passes --repo-root
# explicitly, so the pipeline never depends on this script's depth in the tree.
DEFAULT_REPO_ROOT = Path(__file__).resolve().parents[3]

DEFAULT_PREAMBLE = Path(__file__).resolve().parent / "OH-cascade-preamble.input"


def load_box_size_helpers(repo_root):
    """Import the density-solve helpers from <repo_root>/analysis/general/.

    Imported here rather than at module scope for two reasons: the repo root
    is a runtime argument (the driver knows it and passes it, instead of this
    script guessing from __file__), and a checkout without analysis/general/
    should say which directory it looked in rather than raise a bare
    ModuleNotFoundError naming a module the reader has never heard of.

    The density solve is not reimplemented here — analysis/general already has
    it, and box_size.py's constants are the ones the rest of the repo uses.
    """
    general = Path(repo_root) / "analysis" / "general"
    missing = [
        name for name in ("box_size.py", "box_size_from_data.py")
        if not (general / name).is_file()
    ]
    if missing:
        paths = " ".join(f"analysis/general/{name}" for name in missing)
        raise SystemExit(
            f"Error: cannot find the density-solve helpers under {general}\n"
            f"  missing: {', '.join(missing)}\n"
            f"  repo root: {repo_root}\n"
            f"\n"
            f"Both files are tracked, so a checkout missing one has lost it\n"
            f"locally. Restore it in that checkout with:\n"
            f"    cd {repo_root} && git checkout -- {paths}\n"
            f"\n"
            f"Note that box_size_from_data.py imports box_size.py, so a checkout\n"
            f"without box_size.py cannot run either of them — this is not specific\n"
            f"to the sweep pipeline. Or pass --repo-root pointing at a complete\n"
            f"checkout."
        )

    if str(general) not in sys.path:
        sys.path.insert(0, str(general))

    from box_size import AMU_TO_G, cubic_dimension_A
    from box_size_from_data import read_lammps_data

    return AMU_TO_G, cubic_dimension_A, read_lammps_data

# Fixed by OH-cascade-preamble.input's pair_coeff line. Not configurable: see
# the header of that file.
ELEMENTS = "O H"

KELVIN_OFFSET = 273.15

# Temperature the preamble's "velocity all create" seeds the system at, and so
# the temperature the first ramp starts from.
INITIAL_TEMP_K = 10.0


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate the temperature-cascade LAMMPS input.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--repo-root", default=str(DEFAULT_REPO_ROOT),
                   help="util checkout holding analysis/general/ (the density solve)")
    p.add_argument("--start-data", required=True,
                   help="Starting structure .data file (read for atom counts)")
    p.add_argument("--out", required=True, help="Where to write the input")
    p.add_argument("--preamble", default=str(DEFAULT_PREAMBLE),
                   help=argparse.SUPPRESS)  # fixed; overridable only for tests
    p.add_argument("--replicate", default="6 6 6",
                   help='Supercell replication, e.g. "6 6 6"')
    p.add_argument("--density", type=float, default=1.0,
                   help="Target mass density in g/cc for the deform")

    p.add_argument("--diffusion-temps-c", default="",
                   help="SEMICOLON-separated Celsius list for the NVE/MSD route")
    p.add_argument("--dielectric-temps-c", default="",
                   help="SEMICOLON-separated Celsius list for the dielectric route")
    p.add_argument("--melt-temperature", type=float, default=None,
                   help="Optional Celsius melt above the union max; held, not saved")

    p.add_argument("--timestep", type=float, default=0.00025,
                   help="LAMMPS timestep in PICOSECONDS (metal units)")
    p.add_argument("--seed", type=int, default=156467,
                   help="velocity create seed")

    p.add_argument("--deform-length", type=int, default=100000,
                   help="Steps for the 10 K deform-to-density block")
    p.add_argument("--ramp-length", type=int, default=20000,
                   help="Steps for each temperature ramp")
    p.add_argument("--hold-length", type=int, default=100000,
                   help="Steps held at each temperature before writing it out")
    p.add_argument("--nve-length", type=int, default=50000,
                   help="Steps of NVE production at each diffusion temperature")
    p.add_argument("--dynamics-dump-every", type=int, default=10,
                   help="Steps between frames in dynamics_T<C>.lammpstrj")

    p.add_argument("--print-box-length", action="store_true",
                   help="Print only the solved cube edge (Angstrom) and exit")
    return p.parse_args()


def parse_temp_list(raw):
    """Semicolon-separated Celsius list -> list of normalized float values.

    Semicolons, not commas: these lists reach the SLURM stages through
    sbatch --export, which is comma-delimited and truncates silently.
    """
    out = []
    for chunk in raw.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        try:
            out.append(float(chunk))
        except ValueError:
            raise SystemExit(f"Error: not a number in temperature list: {chunk!r}")
    return out


def celsius_label(t_c):
    """Path-safe Celsius label. '%g' collapses 25 and 25.0 to one spelling, so
    the same temperature written two ways cannot produce two cascade blocks."""
    return f"{t_c:g}"


def to_kelvin(t_c):
    return t_c + KELVIN_OFFSET


def parse_replicate(raw):
    parts = raw.split()
    if len(parts) != 3:
        raise SystemExit(f'Error: --replicate needs three integers, got: {raw!r}')
    try:
        nx, ny, nz = (int(p) for p in parts)
    except ValueError:
        raise SystemExit(f'Error: --replicate values must be integers: {raw!r}')
    if nx < 1 or ny < 1 or nz < 1:
        raise SystemExit(f'Error: --replicate values must be >= 1: {raw!r}')
    return nx, ny, nz


def parse_preamble_masses(preamble_text):
    """{type_id: mass_amu} from the preamble's own "mass N <amu>" lines.

    These, not start.data's Masses section, are what LAMMPS integrates and so
    what the density solve must use. ICE_CUBIC.data says H is 1.0; the preamble
    says 1.00784, which is the difference between a 37.2514 A and a 37.2406 A
    cube at 1 g/cc.
    """
    masses = {}
    for line in preamble_text.splitlines():
        stripped = line.split("#")[0].strip()
        if not stripped.startswith("mass "):
            continue
        parts = stripped.split()
        if len(parts) >= 3:
            masses[int(parts[1])] = float(parts[2])
    if not masses:
        raise SystemExit("Error: no 'mass' lines found in the preamble")
    return masses


def solve_box_length(start_data, replicate, density, preamble_masses, repo_root):
    """Cube edge (Angstrom) giving `density` g/cc for the replicated cell."""
    AMU_TO_G, cubic_dimension_A, read_lammps_data = load_box_size_helpers(repo_root)
    info = read_lammps_data(start_data)
    nx, ny, nz = replicate
    n_cells = nx * ny * nz

    total_amu = 0.0
    for type_id, _symbol, _mass_from_file, count in info["composition"]:
        if type_id not in preamble_masses:
            raise SystemExit(
                f"Error: {start_data} has atom type {type_id}, which the preamble "
                f"has no mass line for (preamble types: {sorted(preamble_masses)})"
            )
        total_amu += preamble_masses[type_id] * count * n_cells

    return cubic_dimension_A(total_amu * AMU_TO_G, density)


def render_preamble(preamble_text, replicate_str, seed):
    rendered = preamble_text.replace("@REPLICATE@", replicate_str)
    rendered = rendered.replace("@SEED@", str(seed))
    for placeholder in ("@REPLICATE@", "@SEED@"):
        if placeholder in rendered:
            raise SystemExit(f"Error: {placeholder} survived substitution")
    return rendered


def render_deform(box_length, timestep_ps, deform_length):
    """The 10 K deform-to-density prep, from OH-therm.input's own block."""
    L = f"{box_length:.4f}"
    return f"""
timestep        {timestep_ps}

################ deform to the target density at 10 K ##################
# Box edge solved from the replicated atom count, the preamble's masses and
# --density (see generate_cascade_input.py). The literal 37.2514 that
# OH-therm.input hardcodes is what this reproduces for ICE_CUBIC 6 6 6 at
# 1 g/cc.
fix             def all deform 10 x final 0.0 {L} y final 0.0 {L} z final 0.0 {L}
fix             1 all nvt temp 10 10 $(100*dt)
run             {deform_length}
unfix           def
unfix           1
############### end deform prep ###############
"""


def render_stop(t_c, prev_k, do_diffusion, args, is_melt=False):
    """One temperature stop: ramp from prev_k, hold, save, optionally NVE."""
    t_k = to_kelvin(t_c)
    label = celsius_label(t_c)
    header = "melt (not saved)" if is_melt else f"T = {label} C"

    block = f"""
############### {header} — {t_k:.2f} K ###############
fix             1 all nvt temp {prev_k:.2f} {t_k:.2f} $(100*dt)
run             {args.ramp_length}
unfix           1

fix             1 all nvt temp {t_k:.2f} {t_k:.2f} $(100*dt)
run             {args.hold_length}
unfix           1
"""

    if is_melt:
        # The melt exists only to erase the ice's memory before descending. It
        # is deliberately not written out and never gets a diffusion run.
        return block + f"############### end melt ###############\n"

    block += f"""
# Starting structure for this temperature's dielectric production stage.
write_data      thermalized_T{label}C.data
"""

    if do_diffusion:
        block += f"""
# --- NVE diffusion production ({label} C) ---
# NVE, matching OH-therm.input's production block: the thermostat is released
# so the dynamics msd.py measures are not coupled to it. The system carries on
# from here into the next ramp down, which is intended.
reset_timestep  0
dump            mydyn all custom {args.dynamics_dump_every} dynamics_T{label}C.lammpstrj id element x y z vx vy vz
dump_modify     mydyn element {ELEMENTS}
dump_modify     mydyn sort id
fix             1 all nve
run             {args.nve_length}
unfix           1
undump          mydyn
"""

    block += f"############### end {label} C ###############\n"
    return block


def build_schedule(args):
    """Descending list of (t_c, do_diffusion) over the union of both routes."""
    diffusion = parse_temp_list(args.diffusion_temps_c)
    dielectric = parse_temp_list(args.dielectric_temps_c)

    if not diffusion and not dielectric:
        raise SystemExit(
            "Error: at least one of --diffusion-temps-c / --dielectric-temps-c "
            "must be non-empty."
        )

    # Deduplicate on the normalized label, so 25 and 25.0 are one stop.
    diffusion_labels = {celsius_label(t) for t in diffusion}
    union = {}
    for t in diffusion + dielectric:
        union[celsius_label(t)] = t

    # Descending: the cascade ramps up once, then only ever cools.
    ordered = sorted(union.items(), key=lambda kv: kv[1], reverse=True)
    return [(t_c, label in diffusion_labels) for label, t_c in ordered]


def main():
    args = parse_args()

    preamble_text = Path(args.preamble).read_text()
    preamble_masses = parse_preamble_masses(preamble_text)
    replicate = parse_replicate(args.replicate)
    box_length = solve_box_length(
        args.start_data, replicate, args.density, preamble_masses, args.repo_root
    )

    if args.print_box_length:
        print(f"{box_length:.4f}")
        return

    schedule = build_schedule(args)

    if args.melt_temperature is not None:
        highest = schedule[0][0]
        if args.melt_temperature <= highest:
            raise SystemExit(
                f"Error: --melt-temperature ({args.melt_temperature} C) must be above "
                f"the highest requested temperature ({highest} C); otherwise it is "
                f"just a stop the cascade already makes."
            )

    parts = [render_preamble(preamble_text, args.replicate, args.seed)]
    parts.append(render_deform(box_length, args.timestep, args.deform_length))

    prev_k = INITIAL_TEMP_K
    if args.melt_temperature is not None:
        parts.append(render_stop(args.melt_temperature, prev_k, False, args,
                                 is_melt=True))
        prev_k = to_kelvin(args.melt_temperature)

    for t_c, do_diffusion in schedule:
        parts.append(render_stop(t_c, prev_k, do_diffusion, args))
        prev_k = to_kelvin(t_c)

    parts.append("\nwrite_restart   final.restart\nwrite_data      final.data\n")

    Path(args.out).write_text("".join(parts))

    n_diff = sum(1 for _, d in schedule if d)
    print(f"Wrote {args.out}")
    print(f"  box edge      : {box_length:.4f} A ({args.density} g/cc, "
          f"replicate {args.replicate})")
    print(f"  temperatures  : {', '.join(celsius_label(t) + ' C' for t, _ in schedule)}"
          f"{'' if args.melt_temperature is None else f' (melt at {args.melt_temperature} C first)'}")
    print(f"  diffusion NVE : {n_diff} of {len(schedule)} stops")


if __name__ == "__main__":
    main()
