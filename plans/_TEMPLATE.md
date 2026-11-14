# <What this plan delivers, in a few words>

**State:** draft (only the owner sets `approved`; code may follow once it is committed approved)

## Goal

What will be true when this ships, and why it matters.

## Scope

Every file or folder the work creates or changes, in backticks: `analysis/<topic>/`,
`analysis/<topic>/new_script.py`, `analysis/<topic>/test_new_script.py`,
`analysis/<topic>/README.md`. Code outside this list is refused at commit.

## Approach

How it will be built, in the order it will be built, with the decisions and the
alternatives rejected.

## Verification

The verification plan. Every subsection needs content; "n/a" needs a reason.

### Known answer

Synthetic inputs whose result is known by construction (see `analysis/testkit/`),
and the expected values.

### Invariants

Sum rules, limits and symmetries the output must satisfy (e.g. g(r) → 1 at large
r, VDOS integrates to 3N, translation and atom-order invariance).

### Independent implementation

A second code computing the same quantity on the same input (e.g. Kenichi's
Fortran vs freud, LAMMPS `compute` output, simledger's parsers, ASE/OVITO).

### Regression

Outputs pinned on a small real input (`analysis/testkit/data/`), and what a change
in them means.

### Tolerances and sources

Each tolerance with its justification (sampling error, written precision,
truncation error), and where each expected value comes from (analytic,
literature with citation, the other code).

### Slow checks

Convergence and sampling studies, marked `@pytest.mark.slow`.

### Done when

The concrete conditions for shipping, including the README `## Verification`
sections for every script touched.

## Out of scope

What this plan deliberately does not do.
