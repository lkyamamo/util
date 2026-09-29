# Chunked dielectric constant — manual, single temperature

One temperature, two `sbatch` steps, everything configured by hardcoded blocks
at the top of the scripts. This is the same chunked dielectric route the sweep
pipeline (`jobs/pipeline/lammps_to_temperature_sweep/`) runs for each of its
dielectric temperatures, structured the way
`20260617_dielectric_multi_traj/` is: numbered scripts you run yourself, no
driver, no temperature list, no aggregation CSV.

The thermalized structure is assumed to exist already — there is no cascade
here. Point `THERMALIZED_STRUCTURE` at a cascade stop
(`cascade/thermalized_T<C>.data`), or at the `final.data` of any earlier run
equilibrated at this temperature and density.

## Files

| File | |
|---|---|
| `0.production_submit.slurm` | Stage 1. Chunked NVT production: `N_TIMES` chunks of `NVT_LENGTH` steps, one `dumps/dielectric.<i>.custom` per chunk, restartable. |
| `dielectric-production.input` | The LAMMPS input stage 1 runs: chunk loop, restart support, and the by-hand dipole cross-check. |
| `1.calc_submit.slurm` | Stage 2. Drives the two analysis scripts below. |
| `2.calc_mpi.py` | MPI over the chunk files → one dipole line per frame, checkpointed per rank. |
| `3.dipole_std.py` | Dipole series → `<M>`, `<M²>`, variance, dielectric constant, per-timestep CSV. |

## Running it

```bash
# copy this directory somewhere on /scratch1 and work there
cp -r <util>/simulation/lammps/20260929_dielectric_chunked  /scratch1/lkyamamo/.../diel-T30C
cd /scratch1/lkyamamo/.../diel-T30C
mkdir -p logs                      # SLURM opens logs/*.out before the job's own mkdir runs

$EDITOR 0.production_submit.slurm  # THERMALIZED_STRUCTURE, TARGET_TEMP, the chunking
sbatch 0.production_submit.slurm

$EDITOR 1.calc_submit.slurm        # TEMPERATURE, LA/LB/LC, --ntasks
sbatch 1.calc_submit.slurm         # once stage 1 has finished
```

Stage 2 can also be submitted with `--dependency=afterok:<stage-1 jobid>` to
queue both at once, which is what the pipeline does.

Results land in `dipole_output/`:

- `dipole_output.txt` — the concatenated per-frame dipole series (e·Å), built
  atomically from the per-rank files once every rank has reported done.
- `summary.txt` — what `3.dipole_std.py` printed, including
  `eps_x/eps_y/eps_z/eps_total`.
- `dipole_output_timestep_data_<method>.csv` — the running/binned averages and
  the dielectric constant against timestep.
- `dipole_by_hand.txt` (in this directory, not `dipole_output/`) — LAMMPS's own
  dipole at the same cadence, for comparison against `dipole_output.txt`.

## Values the two stages share

The pipeline's driver guarantees these agree across the two stages. Here nothing
does, so check them by hand whenever you change one:

| | `0.production_submit.slurm` | `1.calc_submit.slurm` |
|---|---|---|
| Temperature (K) | `TARGET_TEMP` | `TEMPERATURE` |
| Frame cadence | `DUMP_EVERY` | `DUMP_EVERY` |
| H partial charge | `CHARGE_H` | `CHARGE_H` |

Temperatures are in **Kelvin** in both, unlike the sweep pipeline's Celsius
lists. `LA/LB/LC` in stage 2 are the box edges the production run actually had —
read them off the dump's box bounds if the density solve's value (37.2514 Å for
ICE_CUBIC 6×6×6 at 1.0 g/cc) does not apply.

## Resuming

Both stages resume by being submitted again.

- **Stage 1** reads `chunks_done.txt`, which the input script appends to only
  after that chunk's `restart.<i>` is safely written, and restarts from the last
  chunk named there. At most the chunk that was in flight is redone, and its
  dump is rewritten from scratch.
- **Stage 2** counts the lines already in each `dipole_output/ranks/dipole_rank_<r>.txt`
  and picks up from there. The rank→file assignment depends on the rank count, so
  **resume with the same `--ntasks`** — `2.calc_mpi.py` records the count the run
  started with and refuses to continue under a different one rather than
  corrupt the output.

To start over instead of resuming:

- Stage 1: delete `chunks_done.txt`, `restart.*`, `dumps/` and
  `dipole_by_hand.txt` (the input script *appends* to that one).
- Stage 2: delete `dipole_output/` entirely — removing only
  `dipole_output.txt` leaves the per-rank checkpoints in place.

## Deleting dumps as they are reduced

`DELETE_PROCESSED_DUMPS=1` in stage 2 unlinks each chunk file once every frame in
it has become a dipole line and those lines are fsynced. Every line carries all
the dielectric constant needs from its frame, so the dump has no further use —
but this is irreversible, and redoing the calc from scratch afterwards means
re-running stage 1. It defaults to 0 here; the sweep pipeline defaults it to 1,
because a sweep writes far more trajectory at once than disk can hold.

## Keeping in sync with the pipeline

`2.calc_mpi.py`, `3.dipole_std.py` and `dielectric-production.input` are copies:

| Here | Original |
|---|---|
| `2.calc_mpi.py` | `20260617_dielectric_multi_traj/1.calc_mpi.py` |
| `3.dipole_std.py` | `20260617_dielectric_multi_traj/2.dipole_std.py` |
| `dielectric-production.input` | `jobs/pipeline/lammps_to_temperature_sweep/dielectric-production.input` |

Only comments differ (they point at this directory's script names). `diff` them
when either side changes — the pipeline copies its analysis scripts from
`20260617_dielectric_multi_traj/`, so a fix made there does not reach this
directory on its own.
