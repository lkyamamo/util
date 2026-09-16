# Bond breaking analysis

Detects and classifies O–H bond breaking in bulk water from a LAMMPS trajectory,
and writes an annotated trajectory OVITO reads natively.

```bash
python bond_events.py --traj dynamics.lammpstrj --dt 2.5
python test_bond_events.py          # synthetic ground-truth checks
```

| File | What it is |
|---|---|
| `bond_events.py` | The analyzer. Standalone argparse tool, not wired into the pipeline. |
| `test_bond_events.py` | Synthetic ground-truth tests. No fixtures — trajectories are built in the test. |
| `calibrate_cutoffs.py` | **Experimental sketch, not used by anything.** See its STATUS section. |

## Why not just a distance cutoff

The simulations use `pair_style usc` with `atom_style atomic` and `0 bonds`, so
connectivity is dynamic but LAMMPS never reports it. Every other bond
determination in this repo is a per-frame distance cutoff, which is fine for
equilibrium averages and wrong for counting events.

The O–H stretch is ~3700 cm⁻¹, period ≈ 9 fs. `OH-therm.input` sets
`timestep 0.00025` ps with `dynamics_dump_frequency 10`, so `dynamics.lammpstrj`
has **dt = 2.5 fs** — about 3.6 samples per stretch period. Any pair whose
vibration straddles a hard cutoff "breaks" and "re-forms" every couple of frames.
That flicker averages out of a g(r) or a bond-angle distribution; it completely
dominates an event count, which would then scale with your dump frequency rather
than with chemistry.

So bonds run through a Schmitt trigger with a persistence requirement:

```
BONDED   -> UNBONDED   when r > r_break for >= tau_persist consecutive frames
UNBONDED -> BONDED     when r < r_form  for >= tau_persist consecutive frames
```

Excursions shorter than `tau_persist` are counted as **recrossings** and
reported separately. That number is how many breaks a plain single-cutoff
analysis would have reported that this one rejects — it is a headline result,
not a footnote. Check it: if recrossings hugely outnumber breaks, the cutoffs are
sitting inside the vibrational amplitude.

## Parameters

| Flag | Default | Why |
|---|---|---|
| `--r-form` | 1.15 Å | Below the O–H g(r) first minimum |
| `--r-break` | 1.35 Å | Above it — the two bracket the minimum, and the gap is the hysteresis |
| `--tau-persist` | 4 frames | 4 × 2.5 fs = 10 fs ≈ one full O–H stretch period. A state change that does not survive one stretch period is a vibration, not chemistry. |
| `--transfer-window` | 20 frames | 50 fs to watch a freed H before classifying the break |
| `--hbond-cutoff` | 3.5 Å | O···O first solvation shell, for the coordination descriptors |
| `--dt` | 2.5 fs | Time between dumped frames |
| `--temperature` | 300 K | For the `ke_h / (3/2)kT` descriptor only |

These defaults are for water with the OH USC potential. **They are not derived
from your specific run.** Sanity-check them against `rdf_freud.py`'s O–H g(r):
`r_form` and `r_break` should straddle the first minimum. Automating that is
parked in `calibrate_cutoffs.py` and is not trustworthy yet.

`--tau-persist` is counted in *kept* frames, so `--stride` changes the physics.
The tool warns when you use it.

## Reading the results, in OVITO

Open `output/bond_events.lammpstrj`. The extra columns become particle
properties, so **one expression works across the whole trajectory** — no
per-frame editing.

```
EventType != 0                            all atoms involved in an event this frame
EventAge >= 0 && EventAge < 20            atoms that reacted within the last 50 fs
Mechanism == 1                            atoms whose last event was a proton transfer
Mechanism == 1 && EventAge < 20           transfers, held on screen for 50 fs
Mechanism == 3                            dissociations — usually means bad cutoffs
Species == 1                              every hydronium, oxygen and hydrogens
Species == 2                              every hydroxide
Coordination != 2 && ParticleType == 1    any oxygen that is not part of an H2O
PartnerID == 12345                        the other atom in a specific event
```

### Selecting transfers specifically

`EventType` cannot do it on its own: it says an H broke away, but a Grotthuss
hop and a dissociation look identical at that instant — the difference is what
the H does over the *next* 50 fs, and OVITO expressions have no lookahead. That
is what `Mechanism` is for. It is assigned once the event is classified and then
written back to the frames the event happened on:

```
Mechanism == 1 && EventType == 2          the proton, at the instant it leaves
Mechanism == 1 && EventType == 4          the proton, at the instant it arrives
Mechanism == 1 && EventAge < 20           donor, acceptor and proton, for 50 fs
```

Unlike `EventType`, `Mechanism` **persists** after the event rather than
resetting each frame, so it pairs with `EventAge` to control how long a hop
stays highlighted.

Add a **Color Coding** modifier on `EventAge` (range 0–20, reversed) to get a
fade-out trail behind each reaction, or on `Species` for a static species map.

### Column meanings

| Column | Values |
|---|---|
| `Coordination` | Bonded partners this frame, from the hysteretic state — **not** an instantaneous cutoff, so it does not flicker |
| `EventType` | `0` none · `1` O lost a bond · `2` H lost a bond · `3` O gained a bond · `4` H gained a bond |
| `EventAge` | Frames since this atom's last event; `0` at the event, `-1` if it has never had one |
| `PartnerID` | Atom id of the other atom in the most recent event, `-1` if none |
| `Species` | `0` H₂O · `1` H₃O⁺ · `2` OH⁻ · `3` free H · `4` free O · `5` other |
| `Mechanism` | Classification of this atom's most recent event: `0` none · `1` transfer · `2` recross · `3` dissociation · `4` truncated. **Persists**, unlike `EventType`. |

Hydrogens inherit their oxygen's `Species`, so `Species == 1` selects a whole
hydronium rather than just its oxygen.

### Caveat: OVITO's Create Bonds will disagree with these columns

If you switch on OVITO's **Create Bonds** modifier, it draws sticks from a single
distance cutoff — exactly the criterion this tool rejects. The sticks will blink
on and off with the O–H stretch, and at the moments you care about they will
contradict the `Coordination` column: the column says the bond is intact because
the excursion has not persisted, the stick says it is gone.

Mitigation: set Create Bonds' O–H cutoff to **`r_break` (1.35 Å), not 1.2 Å**, so
the drawn bonds are the outer envelope of the tracked ones and never vanish while
`Coordination` still counts them. They will still appear slightly early on
re-formation.

Getting the sticks to match exactly needs an OVITO Python modifier injecting the
tracked topology per frame — deliberately not built, since it adds an `ovito`
dependency and only works from the scripting interface, not the GUI. Everything
needed to write it later is in `events.jsonl`: every event carries `o_id`,
`h_id`, and `frame`, which is enough to replay the exact bond set.

## Understanding *why* bonds break

Every break in `events.jsonl` is classified by what the H does next:

| Class | Meaning |
|---|---|
| `transfer` | H bonds to a **different** O — a Grotthuss proton hop. Should dominate overwhelmingly in equilibrium bulk water. |
| `recross` | H re-bonds to the **same** O after a confirmed unbonded period |
| `dissociation` | H bonds to no O for the whole window. Should be rare. |
| `truncated` | The classification window ran past the end of the trajectory |

and carries the descriptors that explain it:

- **`v_radial`** — relative O/H velocity along the bond, Å/ps. Large and positive
  means the bond was kinetically driven apart. Uses the `vx vy vz` already in the
  dump, which nothing else in this repo reads for this purpose.
- **`ke_h_over_kT`** — the H's kinetic energy over (3/2)kT. Consistently ≫ 1 means
  breaks come from the hot tail rather than from the reaction coordinate.
- **`n_o_neighbors_donor`** — O···O neighbours of the donor. Was it over-solvated?
- **`oo_distance_at_break`**, **`ohO_angle_at_break`** — the Grotthuss reaction
  coordinate. A hop at O···O ≈ 2.4 Å with a near-linear angle is textbook.

`bond_events.png` plots the last of these as a 2D histogram; it is the panel that
actually shows the mechanism.

## Sanity checks on a production run

- `transfer` should dominate. A large `dissociation` fraction means the cutoffs
  are wrong, not that the water is reacting.
- H₃O⁺ and OH⁻ counts should track each other; the `charge imbalance` line in
  `summary.txt` should hover near zero. A steady drift means bonds are being
  lost or gained without a matching partner — again a cutoff problem.
- `oo_distance_at_break` well above ~3 Å means the tool is pairing an H with an
  acceptor it never actually reached.

## Known limitation

At 2.5 fs/frame there are only ~3.6 samples per O–H stretch period. That is
enough to *identify* events with the persistence criterion but not enough to
resolve the transition path of a single hop, so do not over-interpret
`r_at_break` or `v_radial` for individual events — they are distributions, not
trajectories. For barrier-crossing detail, `dynamics_dump_frequency` in
`OH-therm.input:73` would need to drop from 10 to ~2.

## Outputs

| File | Contents |
|---|---|
| `bond_events.lammpstrj` | Input trajectory plus the five annotation columns |
| `events.jsonl` | One JSON object per event with every descriptor |
| `events_per_frame.csv` | Per-frame event counts and species populations |
| `summary.json` | Machine-readable counters (what the tests assert on) |
| `summary.txt` | The aggregate report |
| `bond_events.png` | Four diagnostic panels |

Skip the annotated trajectory with `--no-ovito-dump`; it is the slow half of the
run and the only part that rewrites the whole trajectory.

## Tests

`test_bond_events.py` builds trajectories whose answer is known by construction,
runs the real script on them, and compares — the same approach as
`analysis/shock/etc/test_pipeline.py`. There is no multi-frame water fixture in
this repo to use instead.

| Case | Ground truth |
|---|---|
| `flicker` | A pair oscillating across `r_break` every few frames: **zero** breaks, one recrossing per cycle. This test is the argument for the hysteresis criterion, made executable — a plain cutoff scores ~10 breaks. |
| `transfer` | One break classified `transfer`, donor 2→1 (OH⁻), acceptor 2→3 (H₃O⁺), O···O = 2.6 Å, plus every annotation column checked |
| `dissociation` | H leaves and never returns: one break, classified `dissociation` |
| `pbc` | The transfer case across a periodic face: results **identical**. Catches the bug class in `voxel_analysis.py:294`, which builds a `cKDTree` with no `boxsize=` and is silently non-periodic. |
| `eof` | A break whose window runs off the end: reported as `truncated` and counted in `n_pending_at_eof`, not silently dropped |
