# Plotting pipeline

Turns an analysis directory's CSVs into figures. This is the only place in the
repo that draws distribution plots.

The analysis scripts in `analysis/distributions/20260608_GrNrBaSqw/` write CSVs
and nothing else. Keeping the stages apart means:

- a figure can be restyled, re-unit'd, reformatted or redrawn **without
  recomputing anything** — the expensive half of a run is never repeated to
  change how a curve looks;
- there is **one** implementation of each plot, not one per analysis script;
- a plot can never disagree with the table it came from, because if a curve is
  not in a CSV it cannot be drawn.

## Running it

```bash
# in an analysis directory
~/util/jobs/pipeline/plotting/plot_pipeline.sh

# or point at one
~/util/jobs/pipeline/plotting/plot_pipeline.sh /path/to/analysis_dir

# one-off overrides, without touching a config
plot_pipeline.sh . --style publication --date 20260813
```

It is also called automatically as the last stage of `distribution_run.sh` and
`distribution_submit.slurm` (`RUN_PLOTS=1`, the default), and through
`submit_pipeline.sh --run-plots 1`.

## Config resolution

The first of these that exists wins, and **the one used is printed on every
run**, so it is never a guess:

| order | file | when |
|---|---|---|
| 1 | `--config <file>` | forced explicitly |
| 2 | `<analysis dir>/plot_pipeline.conf` | this run carries its own look |
| 3 | `jobs/pipeline/plotting/plot_pipeline.conf` | the tracked default — what most runs use |

To make one analysis directory's figures differ, copy the default into it and
edit the copy. Nothing tracked has to change. `submit_pipeline.sh --plot-config
FILE` does that copy for you when setting a run up.

Unlike `submit_pipeline.conf`, the default here is **tracked, not gitignored** —
it holds no personal paths, so it is a working default rather than an example.

## Output

One PNG per quantity, in `<date>_<calc>/` beside the CSV it came from:

| directory | from | contents |
|---|---|---|
| `<date>_rdf/` | `rdfs.csv`, `nrs.csv`, `wright.csv` | each partial, each convention column and its `_broadened` twin, each n(r), the wright curves |
| `<date>_bad/` | `bads.csv` | one per A-B-C triplet |
| `<date>_dsf/` | `sq.csv`, `dsf.csv` | S(q) partials/total/neutron, the S(q,ω) heatmaps |
| `<date>_vdos/` | `vdos.csv` | one per species and per weighted total |
| `<date>_msd/` | `msd.csv` | one per species |
| `<date>_vdos_dynmat/` | `vdos_dynmat.csv`, `..._modes.csv` | DOS curves, stretch/bend/rock, participation ratio, reduced DOS |

**The date comes from the CSV, not from today.** Plotting `20260813_bads.csv`
writes `20260813_bad/`, so figures stay married to the data that produced them
and re-plotting an old run does not stamp it with the current date. When a
directory holds several runs, `PLOT_DATE` chooses: `latest` (default), `all`, or
an explicit `YYYYMMDD`.

## What is and is not drawn

**No reference lines.** An individual plot carries its curve and nothing else —
no `g(r) = 1` line, no baseline under `T(r)`. The asymptotic limit each
convention column should approach is printed by `rdf_freud.py`'s convention
table instead, which is the check that actually catches a wrong normalization.

**A composite only where the combination is itself the result**, named so it
reads as one:

- `wright_composite.png` — T(r) broadened, raw, and the T⁰(r) baseline together;
  the figure a published curve is overlaid on.
- `character_composite.png` — stretch/bend/rock against the total. A band
  assignment is made by reading them against each other.
- `all_curves.png` / `all_species.png` — a partial read against the total it
  sums into, or one species outrunning another.

`PLOT_COMPOSITES=0` writes only the individual curves.

## Styles

- `analysis` (default) — titled, fully labelled, y ticks. For reading a run.
- `publication` — heavy lines and spines, large bold labels, no y ticks. For
  dropping a single curve into a document.

`publication` is what the old `rdf_freud_plot.py` / `bad_freud_plot.py` second
pass produced. Those scripts are gone; the style is a config value now, and it
applies to all six calculations rather than the two they covered.

## Adding a calculation

1. Write a function `plot_<name>(target_dir, date, files)` in
   `plot_distributions.py` returning `(out_dir, n_plots)`.
2. Add its CSV stem to `discover()`.
3. Add a row to `CALCULATIONS` with its switch and the CSV kinds it needs.
4. Add the switch to `plot_pipeline.conf` and to the pre-declare list in
   `plot_pipeline.sh`.
