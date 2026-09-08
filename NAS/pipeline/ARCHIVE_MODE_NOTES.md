# Archive mode — removed, and what it would take to bring back

`nas.sh` once had a second transfer mode, `archive`: compress a remote tree into
`.tar.zst` / `.zst` pieces with SLURM jobs on the HPC, then rsync the artifacts
down. It was removed because it had never produced a correct result — the tree
walk was silently broken and the script could not even load on the transfer Mac.

## Getting the code back

The removed implementation is in git, unmodified:

```bash
git show d44d55b:NAS/pipeline/nas.sh
```

`d44d55b` is pinned deliberately — it is the last commit before archive mode was
removed, so this command keeps working no matter what happens to branch tips. The
relevant functions are `to_bytes`, `remote_dir_size_bytes`,
`remote_file_size_bytes`, `walk_tree`, `generate_slurm_script`,
`submit_slurm_job`, `poll_job`, `check_job_log`, `rsync_target`,
`cleanup_remote`, and `process_archive`.

`nas.sh` still has `resolve_entry` and the `"dirname:mode"` entry syntax, so
re-enabling is a dispatch change in `cmd_sync` rather than a rewrite. Any mode
other than `direct` currently errors and points here.

## Defects to fix first

1. **bash 4+ syntax; only `/bin/bash` 3.2 exists on the transfer Mac.**
   `poll_job` uses `"${state^^}"`; `process_archive` uses `declare -A active_jobs`
   and `"${!active_jobs[@]}"`. Replace `^^` with `tr '[:lower:]' '[:upper:]'`, and
   the associative array with two parallel indexed arrays (`JOB_IDS[]` /
   `JOB_TARGET_IDX[]`). Note this was not the only such bug: `${unit^^}` in
   `to_bytes` ran at *top level*, so every invocation of the whole script died
   with `bad substitution` before reaching the subcommand dispatch.

2. **`walk_tree` loses its separators.** `find ... -print0` is captured through
   `$( )`, and bash discards NUL bytes, so the entire listing collapses into one
   string and the `read -r -d ''` loop sees a single bogus item. Use
   `find ... -printf '%y\t%p\n'`: one round trip returning type *and* path
   together. That also removes the per-item
   `robust_ssh "[ -d '$item' ] && echo dir || echo file"` call, which is one SSH
   round trip **per file** and is what made the walk unusably slow on a big tree.

3. **The bundle target in `generate_slurm_script` is wrong twice.**
   `tar ... -C "$remote_source" $(find . -maxdepth 1 -type f)` — the `$(find .)`
   is evaluated in the SLURM script's cwd (`$HOME`), not `remote_source`, so it
   lists the wrong directory; and even pointed correctly it re-tars the large
   files `walk_tree` already emitted as individual `.zst` targets. Pass the
   explicit small-file list `walk_tree` computed, via `tar -T`.

4. **`cleanup_remote` deletes before it pulls.** The comment says "pull log to
   local before deleting", but the `rm -f` of the tarball and job log runs first,
   so the following rsync of that log always fails (swallowed by `|| true`).
   Swap the order.

5. **Manifest rows are keyed on the bare target name.** Names are only unique
   within a directory — two subdirectories both named `data` both yield
   `data.tar.zst`, and `manifest_update`'s `$2==name` then rewrites both rows.
   Key on `path/name`.

6. **`check_job_log` emits `0\n0`.** `grep -c 'NAS_SUCCESS' "$log" || echo 0`
   makes grep print its own `0` *and* exit 1, triggering the `|| echo 0`; the
   multiline result then breaks the `-gt` test. Use `grep -q`.

## Design note — it must rejoin the two-phase ordering

`nas.sh` now transfers in two global phases: every directory completes its main
transfer before any directory begins its deferred (trajectory) transfer. Archive
mode predates that and would ignore it.

When it returns, `walk_tree` should tag each target as main or deferred by
matching the source path against `RSYNC_PRIORITY_LAST` and `SKIP_PATTERNS`, so
phase 1 compresses and pulls every directory's main targets before phase 2 touches
any trajectory target.

Two things also changed underneath the old code, so it cannot be pasted back
as-is:

- `REMOTE_LOG_BASE` was removed from `nas.sh` (it was archive-only).
- The manifest schema gained a sixth column (`skipped`), and the old
  `COMPRESSING` / `COMPRESSED` / `UPLOADING` / `UPLOADED` statuses were replaced
  by the phase vocabulary (`UPLOADING_MAIN` → `MAIN_DONE` → `UPLOADING_DEFERRED`
  → `UPLOADED`, plus `FAILED_MAIN` / `FAILED_DEFERRED`). Archive mode needs its
  own status names, or a `type` column the status logic branches on.

## Config keys it needs

These were removed from `nas.config`; restore them when the mode returns.

```bash
# === Compression ===
COMPRESSION_LEVEL=3        # 1 (fastest) → 19 (smallest)
COMPRESSION_THREADS=16     # sets both zstd -T and SLURM --cpus-per-task

# === Thresholds ===
# Dirs larger than MAX_ARCHIVE_SIZE are not bundled — files compressed individually
MAX_ARCHIVE_SIZE=100G
# Files smaller than MIN_FILE_SIZE are bundled into a single tar.zst instead of individual .zst
MIN_FILE_SIZE=10M

# === SLURM ===
SLURM_PARTITION=priya
SLURM_TIME=4:00:00
MAX_CONCURRENT_JOBS=4      # max simultaneous SLURM compression jobs per directory
```
