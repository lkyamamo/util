"""simledger: a read-only ledger of LAMMPS/VASP runs on the HPC filesystem.

The scanner walks project trees laid out as

    <project>/runs/[group/...]/<run_id>/{input_files,run,setup}
    <project>/analysis/<run_id>_<description>/
    <project>/potentials/   <project>/structures/

and records what it finds in a SQLite database plus Markdown run cards under
the ledger home (``$SIMLEDGER_HOME``, default ``~/ledger``). It never writes
inside the scanned tree. SLURM job scripts add job-time facts (job id, exit
code, script text) by dropping small JSON events into ``<ledger home>/inbox``
with ``simledger hook``; the next scan ingests them.
"""

__version__ = "0.1.0"

# Bump whenever parsing changes in a way that should refresh stored runs.
PARSER_VERSION = 2
