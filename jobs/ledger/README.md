# simledger: ledger of LAMMPS/VASP runs

`simledger` scans project directories on the cluster and records every run in a searchable
SQLite database, plus one Markdown card per run. For each run it records:

- inputs and parameters
- status, with the evidence behind it
- results
- sub-runs and frames
- file inventory
- linked analysis
- notes

It never writes inside the scanned tree. Python 3.8+, standard library only.

## Where the ledger lives

`$SIMLEDGER_HOME`, default `~/ledger` (on endeavour: `/home1/lkyamamo/ledger`):

```
ledger.db             SQLite database, the source of truth
inbox/                job events dropped by SLURM hooks; read in by the next scan
cards/INDEX.md        every run, one line each
cards/runs.csv        flat table for Excel/pandas
cards/last_scan.md    what the latest scan found, plus open issues
cards/<project>/README.md     project summary, potentials/structures, analysis
cards/<project>/<run_id>.md   one card per run
```

From a Mac, browse it over VS Code Remote-SSH, or copy it:
`rsync -a lkyamamo@endeavour.usc.edu:ledger/ ~/ledger/`.

## Use

```bash
SL=/home1/lkyamamo/util/jobs/ledger/bin/simledger

python3 $SL survey /scratch2/lkyamamo          # layout check only: no parsing, no writes
python3 $SL scan /scratch2/lkyamamo            # incremental: only changed runs are reparsed
python3 $SL scan /scratch2/lkyamamo --full     # reparse everything

python3 $SL show 0239                          # run card (id, project/id, or any path inside a run)
python3 $SL show . --files                     # ...plus every file
python3 $SL list --project 20260825_locked --status failed
python3 $SL search code=vasp "ENCUT>=520" status=completed
python3 $SL search group~ring "energy_span>1"
python3 $SL search --text "water dissociation"     # READMEs, inputs, notes, analysis names
python3 $SL note 0256 "reference energy for the silanol series" --tag reference
python3 $SL sql "SELECT run_key, value_num FROM results WHERE key='energy_per_atom'"
```

Search terms are ANDed together:
- `key=value`, `key!=value`, `key>=n`, `key<n` and `key~substring`.
- `key` is a run field (`code`, `status`, `type`, `group`, `atoms`, `cores`, ...) or any parameter/result key.
- A parameter key can be limited to one source: `incar.ENCUT`, `slurm.partition`, `lammps.pair_style`, `user.temperature_C`.
- Bare words are full-text terms.

## What a card shows

Each card opens with a one-line summary, for example *"NVT MD of 5184 atoms (O 1728, H 3456) in a
37.25 Å cubic box, 1 g/cm³ at 288.15 K; 15 ns (60,000,000 steps × 0.00025 ps); pair usc (OH.usc)."*
It then has these sections:

- **System**: atoms, composition, box, density, units, potential, starting structure (and,
  for a broken link, the run it came from), timestep, length, ensemble, target and
  measured T/P/density.
- **Protocol**: one row per `run`/`minimize` stage, with ensemble, target T/P, steps, dt and
  time, plus the mean T/P/density over the second half of that stage. Repeated
  identical stages are merged ("468,750 × 128").
- **Sub-runs**: a conditions column for each sub-run (e.g. "NVT, 288.15 K, 15 ns"), and a
  frame set's energy span.

**Sources by code:**

| Code | Source |
|---|---|
| LAMMPS | Mainly `log.lammps` (plus `log.lammps.resumeN`, `log_*.lammps`). With `echo both` the log holds every command with `${var}`/`$(...)` already substituted, including values passed with `-var`, so it records what actually ran. |
| LAMMPS | The data file from `read_data` (or the `write_data` end state when that is missing) adds atoms per element and masses. `replicate` is accounted for. |
| LAMMPS | The input script is used only when no log has echoed `run` commands. Its `${...}` values then show as unresolved. |
| VASP | POSCAR/CONTCAR, INCAR (functional, ENCUT, MD settings), KPOINTS, POTCAR titles and masses, OSZICAR (measured MD temperature). |

## Analysis results

Each scan reads small result files and attaches their values to the matching run or sub-run.
Every value records its source file, so it can be traced back:

| File | Values | Attached to |
|---|---|---|
| `*_vs_temperature.csv` (sweep top) | every column, e.g. `eps_total`, `D_total` | the sub-run at that temperature |
| `dipole_output/summary.txt` | `eps_x/y/z/total`, dipole deviation, frames | the run/sub-run the analysis is linked to |
| `*diffusion.csv` (`label,D_...`) | `D_H`, `D_O`, `D_total` | " |
| `eos_summary.csv` | `bulk_modulus_P` (−V₀·dP/dV), `bulk_modulus_E0K` (V₀·d²E/dV²), V₀ | " |
| `SUMMARY.txt` (ring barriers) | per ring size: count and p10/median/p90/min barrier; or, per case, `barrier`, `barrier_frame`, `tail_dE` | the run, or the case's frame-set sub-run |
| `output.txt` (one-row polars table) | `avg_<column>` | the linked run |
| `*rdf*.csv` | `rdf_peak_r_<pair>`, `rdf_peak_g_<pair>` (position and height of the g(r) maximum) | " |
| any other CSV | column names (searchable); its values if it has one row | " |

Values the ledger computes rather than reads (bulk moduli, RDF peaks) carry a "derived: ..."
note. A result file containing a Python traceback is reported as an open issue, naming the
error. Each analysis folder also gets a "what it is" description, taken from its scripts'
docstrings.

Analysis folders link to runs by their name, e.g. `0239_T60C_msd` → sub-run `T60C` of 0239.
A sweep's card has a **Results by sub-run** table (ε, D, ... per temperature). The project
summary shows each run's headline result, as a range across its sub-runs where there is one.
Search works on all of these values, e.g. `simledger search "eps_total>80"` or `"D_total>4"`.

## Expected layout

```
<project>/runs/[group/...]/<name ending in a 4-digit id>/{input_files,run,setup}
<project>/analysis/<id>_<description>/      also 0242_0243_x, 0242-0245_x, 025_6-8_x
<project>/potentials/   <project>/structures/   (README.md "## <file>" sections become descriptions)
```

**Sub-runs and frames.** Inside a run, any directory holding `log.lammps`, `OUTCAR`, `OSZICAR`,
`vasprun.xml`, `STREAM_OUTPUT` or `INCAR` is a calculation. Calculations other than `run/`
become sub-runs, for example the temperatures of a sweep (`T45C`). Numbered sibling
directories (`.../1`, `.../2`, ...) are grouped as frames of one sub-run, and that sub-run
gets energy min/max/span and the highest-energy frame. An analysis named `<id>_T60C_...`
links to sub-run `T60C`.

**Status values.** Each one is backed by stored evidence:

| Status | Evidence |
|---|---|
| `completed` | `Total wall time` in the log, or the OUTCAR timing block |
| `failed` | a LAMMPS `ERROR`, a cancelled/time-limit/OOM job in `STREAM_OUTPUT`, or a non-zero exit from a hook |
| `unconverged` | a relaxation that finished without reaching the required accuracy |
| `running` | the log was modified within the last hour and has no end marker |
| `incomplete` | no end marker and not running |
| `partial` | sub-runs are mixed |
| `not_started` | no calculations yet |

**Inputs.** The copy inside the calculation directory is what actually ran, so it is read
first; `input_files/` is the fallback. Symlinks are inventoried and read only when they
resolve. Broken ones, such as links still pointing at `/scratch1`, are counted.

POTCAR contents are never stored, only the TITEL, VRHFIN, ZVAL, POMASS and ENMAX values.

## SLURM hooks

At the end of a job, a hook drops a small JSON file into `inbox/`. It never touches the
database and always exits 0, so it cannot fail a job. The event records:
- the job id, node list and exit code
- the submitted script (from `scontrol write batch_script`)
- any `--field` values

The next scan reads the events in and matches them to runs through `<project>/runs/.../<id>`.
This still works after the data moves, for example from `/scratch1` to `/scratch2`.

```bash
python3 /home1/lkyamamo/util/jobs/ledger/bin/simledger hook run --dir "$RUN_DIR" --exit-code "$LAMMPS_STATUS" \
  || echo "WARNING: ledger hook failed" >&2

python3 /home1/lkyamamo/util/jobs/ledger/bin/simledger hook analysis --dir "$INPUT_DIR" --type msd \
    --output "${OUT_DIR}/msd.dat" --field temperature_C="${T}" \
  || echo "WARNING: ledger hook failed" >&2
```

The scripts in `jobs/slurm/`, `analysis/` and `simulation/lammps/` already call these.

## Nightly scan

`crontab -e` on a login node:

```
30 2 * * * python3 /home1/lkyamamo/util/jobs/ledger/bin/simledger scan /scratch2/lkyamamo >> /home1/lkyamamo/ledger/scan.log 2>&1
```

A full first scan of `/scratch2/lkyamamo` (103 runs) takes about 6 minutes. Incremental rescans
take about 10 seconds.

## Config

`<ledger home>/simledger.toml` (Python 3.11+) can override `run_id_regex`, `exclude`,
`max_walk_depth`, `tail_bytes`, `running_window_s` and the other fields in
[config.py](simledger/config.py).

## Tests

```bash
cd jobs/ledger && python3 -m unittest discover -s tests -v
```
