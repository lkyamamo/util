#!/usr/bin/env bash
# =============================================================================
# nas.sh — central NAS transfer script
#
# Pulls directories from the HPC down to local NAS storage with rsync, in two
# global phases so that every directory's small/important files land before any
# bulk trajectory data starts moving:
#
#   PHASE 1 (main)     for every directory: priority files, then everything
#                      except the deferred patterns
#   PHASE 2 (deferred) for every directory: the deferred patterns
#                      (trajectories) — skippable entirely
#
# Usage:
#   ./nas.sh sync [--force] [--skip-trajectory|--no-skip-trajectory] [--dry-run] [dirname]
#   ./nas.sh status
#   ./nas.sh manifest [--force] <dirname>
#
# Requires bash 3.2 (the macOS system bash) — do not introduce bash 4+ syntax
# such as ${var^^}, ${var,,}, declare -A, or mapfile.
#
# SYMLINKS
#   `rsync -a` implies -l, so symlinks are copied *as symlinks*. The links under
#   runs/ are absolute /scratch1/lkyamamo paths that do not exist locally, so they
#   arrive dangling — the listing looks complete while the input it points at is
#   missing.
#
#   The fix is a relink pass, not --copy-unsafe-links. These targets are not
#   really lost: they are inside the transfer, just written in the *remote*
#   namespace. relink translates REMOTE_BASE/... to LOCAL_BASE/... and rewrites
#   the link as a RELATIVE one, which costs nothing and keeps the tree portable
#   if the drive is renamed.
#
#   Why not --copy-unsafe-links: rsync counts every absolute symlink as unsafe,
#   so it would replace each one with a full copy of its target. Most are small
#   (potentials ~2 KB, start.data ~800 KB) but some are DIRECTORY links into
#   another run's output — analysis/small_interface/voxel_0147/dumps ->
#   runs/small_interface/0147/run/full is ~2.1 TB, copied a second time. It also
#   destroys the sharing: runs/0177/run/final.data is the start.data for six
#   other runs, and as six independent copies that provenance is gone.
#
#   relink also works before the target exists, so a link pointing at trajectory
#   data can be fixed in phase 1 and simply starts resolving when phase 2 lands
#   the file. That is why it runs after *each* phase, not only at the end.
#
#   Links whose target is not under REMOTE_BASE cannot be translated. Those are
#   left untouched and recorded in LOG_DIR/unresolved_symlinks_<dirname>.tsv.
#
# DELETION HAZARD
#   Symlink targets are frequently *other runs'* outputs — runs/0177/run/final.data
#   is the start.data for six other runs. Pruning old runs by age or size will
#   silently break newer runs' provenance. Check inbound links before deleting.
#
# Archive mode (SLURM-side compression) was removed; see ARCHIVE_MODE_NOTES.md.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The config lives on the NAS itself, so every checkout and git worktree of this
# repo shares one config rather than each carrying a divergent copy. Fall back to
# the copy beside the script when the drive is not mounted (or on Linux, where
# /Volumes does not exist), and let NAS_CONFIG override both.
NAS_CONFIG_DEFAULT=/Volumes/Elements/nas.config
CONFIG="${NAS_CONFIG:-}"
if [[ -z "$CONFIG" ]]; then
    if [[ -f "$NAS_CONFIG_DEFAULT" ]]; then
        CONFIG="$NAS_CONFIG_DEFAULT"
    else
        CONFIG="${SCRIPT_DIR}/nas.config"
    fi
fi

SSH_CONTROL="${TMPDIR:-/tmp}/nas_ssh_ctl_$$"

# =============================================================================
# Load config
# =============================================================================

[[ -f "$CONFIG" ]] || {
    echo "ERROR: nas.config not found at $CONFIG"
    echo "       Expected ${NAS_CONFIG_DEFAULT} (is the drive mounted?)"
    echo "       or ${SCRIPT_DIR}/nas.config; override with NAS_CONFIG=<path>."
    exit 1
}
echo "  [config] $CONFIG"
source "$CONFIG"

# Transfer state (per-directory manifests, rsync logs) lives NEXT TO THE CONFIG,
# not next to the script. The repo has many git worktrees, and keying state to
# SCRIPT_DIR meant each of them tracked its own progress: running `sync` from one
# and `status` from another reported "(no manifest)" for everything while the real
# state sat elsewhere. Anchoring to the config means one drive, one set of state,
# whichever copy of the script you invoke. nas.config may override either.
CONFIG_DIR="$(cd "$(dirname "$CONFIG")" && pwd)"
MANIFEST_DIR="${MANIFEST_DIR:-${CONFIG_DIR}/manifests}"
LOG_DIR="${LOG_DIR:-${CONFIG_DIR}/logs}"

# rsync binary. Do not hardcode a path: this script runs on macOS (where the
# system rsync is 2.6.9 and too old for --info=progress2, so a Homebrew build is
# wanted) and on Linux hosts where rsync is simply /bin/rsync or /usr/bin/rsync.
# An explicit RSYNC in nas.config or the environment always wins.
if [[ -z "${RSYNC:-}" ]]; then
    for _candidate in /opt/homebrew/bin/rsync /usr/local/bin/rsync; do
        [[ -x "$_candidate" ]] && { RSYNC="$_candidate"; break; }
    done
    unset _candidate
    : "${RSYNC:=$(command -v rsync 2>/dev/null || true)}"
fi

HPC="${HPC_USER}@${HPC_HOST}"
FORCE="false"
DRY_RUN="false"
# SKIP_TRAJECTORY comes from nas.config; --skip-trajectory/--no-skip-trajectory override it
SKIP_TRAJECTORY="${SKIP_TRAJECTORY:-false}"
# Repair remote-namespace symlinks after each phase; --no-relink turns it off
RELINK="${RELINK:-true}"

mkdir -p "$MANIFEST_DIR" "$LOG_DIR"

# =============================================================================
# Preflight
# =============================================================================

preflight_local() {
    if [[ -z "${RSYNC:-}" || ! -x "$RSYNC" ]]; then
        echo "ERROR: no usable rsync found${RSYNC:+ at $RSYNC}."
        echo "       Set RSYNC in ${CONFIG} (or the environment) to its full path."
        echo "       On macOS the system rsync (2.6.9) is too old — brew install rsync."
        exit 1
    fi
}

preflight_remote() {
    if ! robust_ssh true; then
        echo "ERROR: cannot reach ${HPC} over ssh — aborting before any transfer."
        exit 1
    fi
}

# =============================================================================
# SSH helpers
# =============================================================================

setup_ssh_control() {
    ssh -fNM \
        -o ControlMaster=yes \
        -o ControlPath="$SSH_CONTROL" \
        -o ConnectTimeout=60 \
        -o ServerAliveInterval=30 \
        -o ServerAliveCountMax=5 \
        "$HPC" 2>/dev/null || true
}

cleanup_ssh_control() {
    ssh -O exit -o ControlPath="$SSH_CONTROL" "$HPC" 2>/dev/null || true
}

robust_ssh() {
    local attempt=0 delay=10
    while (( attempt < 5 )); do
        ssh -o ControlPath="$SSH_CONTROL" \
            -o ControlMaster=no \
            -o ConnectTimeout=60 \
            "$HPC" "$@" && return 0
        attempt=$(( attempt + 1 ))
        echo "  [ssh] attempt $attempt failed, retrying in ${delay}s..."
        sleep "$delay"
        delay=$(( delay * 2 ))
    done
    return 1
}

robust_rsync() {
    # Usage: robust_rsync [--logfile <path>] [extra rsync flags...] <src> <dest>
    local logfile=""
    if [[ "${1:-}" == "--logfile" ]]; then
        logfile="$2"
        shift 2
    fi

    local dry_flag=()
    [[ "$DRY_RUN" == "true" ]] && dry_flag=(--dry-run)

    local attempt=0 delay=10 rc=0
    while (( attempt < 5 )); do
        # NOTE: rsync must not run as a bare command here. Under `set -e` a bare
        # failing command aborts the whole script before the retry logic runs;
        # that only happens to work today because every caller sits inside an
        # `if` condition. Guard it explicitly so the retries are unconditional.
        #
        # NOTE: bash 3.2 + `set -u` treats an empty array expansion as unbound,
        # so every optional array below must use the "${arr[@]:+...}" form.
        if [[ -n "$logfile" ]]; then
            set +e
            $RSYNC \
                -az \
                --partial \
                --partial-dir=.rsync-partial \
                --info=progress2 \
                "${dry_flag[@]:+${dry_flag[@]}}" \
                -e "ssh -o ControlPath=$SSH_CONTROL -o ControlMaster=no" \
                "$@" 2>&1 | tee -a "$logfile"
            rc="${PIPESTATUS[0]}"
            set -e
        else
            set +e
            $RSYNC \
                -az \
                --partial \
                --partial-dir=.rsync-partial \
                --info=progress2 \
                "${dry_flag[@]:+${dry_flag[@]}}" \
                -e "ssh -o ControlPath=$SSH_CONTROL -o ControlMaster=no" \
                "$@"
            rc=$?
            set -e
        fi
        [[ $rc -eq 0 ]] && return 0

        # rc 1 (syntax/usage) and 2 (protocol incompatibility) are our fault, not
        # the network's — retrying them just burns 150s of backoff before failing.
        if (( rc == 1 || rc == 2 )); then
            echo "  [rsync] rc=$rc is not transient (usage/protocol error) — not retrying"
            return 1
        fi

        attempt=$(( attempt + 1 ))
        echo "  [rsync] attempt $attempt failed (rc=$rc), retrying in ${delay}s..."
        sleep "$delay"
        delay=$(( delay * 2 ))
    done
    return 1
}

# =============================================================================
# Pattern helpers
# =============================================================================

# Two lists govern what transfers, and they answer different questions:
#
#   RSYNC_PRIORITY_LAST — WHEN. These are held back to phase 2, so every
#       directory's main files land before any bulk data starts. Without
#       --skip-trajectory they all still transfer, just last. With
#       --skip-trajectory they are not transferred at all, and the directory
#       rests at MAIN_DONE so a later run without the flag completes it.
#
#   ALWAYS_EXCLUDE — WHETHER, unconditionally. Never transferred, in any pass,
#       whatever --skip-trajectory is set to. This is the list for data that is
#       derived and disposable: the dielectric dumps, for instance, are reduced
#       to dipole lines on the HPC and the sweep pipeline deletes them there.
#
# The distinction matters because ALWAYS_EXCLUDE entries need not be trajectories
# at all — dielectric.*.custom matches neither *.lammpstrj nor *.dump, so without
# this list it would transfer in phase 1, ahead of everything it dwarfs.
ALWAYS_EXCLUDES=()
build_always_excludes() {
    ALWAYS_EXCLUDES=()
    local pat
    for pat in "${ALWAYS_EXCLUDE[@]:-}"; do
        [[ -n "$pat" ]] && ALWAYS_EXCLUDES+=(--exclude="$pat")
    done
}

# Comma-joined list of what this directory is still missing, for the manifest's
# `skipped` column: the permanent exclusions, plus the deferred globs while
# --skip-trajectory is suppressing them. "-" when nothing is outstanding.
skipped_patterns_field() {
    local out="" pat
    for pat in "${ALWAYS_EXCLUDE[@]:-}"; do
        [[ -z "$pat" ]] && continue
        out="${out:+${out},}${pat}"
    done
    if [[ "$SKIP_TRAJECTORY" == "true" ]]; then
        for pat in "${RSYNC_PRIORITY_LAST[@]:-}"; do
            [[ -z "$pat" ]] && continue
            out="${out:+${out},}${pat}"
        done
    fi
    echo "${out:--}"
}

has_deferred_patterns() {
    local pat
    for pat in "${RSYNC_PRIORITY_LAST[@]:-}"; do
        [[ -n "$pat" ]] && return 0
    done
    return 1
}

# =============================================================================
# Relink — repair symlinks written in the remote namespace
#
# rsync -a copies symlinks as symlinks, so absolute REMOTE_BASE/... targets land
# pointing at paths that do not exist here. Translate them to LOCAL_BASE and
# rewrite as relative links. See the SYMLINKS note in the header for why this is
# preferred over --copy-unsafe-links.
# =============================================================================

# Relative path from a directory to a target, both absolute. Pure bash so it
# works with the macOS 3.2 shell and without GNU realpath.
relpath() {
    local from="$1" to="$2" common up=""
    common="$from"
    while [[ "$to" != "$common/"* && "$common" != "/" ]]; do
        common="$(dirname "$common")"
        up="../$up"
    done
    [[ "$common" == "/" ]] && { echo "$to"; return; }
    echo "${up}${to#$common/}"
}

relink_dir() {
    local dirname="$1"
    local root="${LOCAL_BASE}/${dirname}"
    [[ -d "$root" ]] || return 0

    local report="${LOG_DIR}/unresolved_symlinks_${dirname}.tsv"
    mkdir -p "$LOG_DIR"
    printf "link\ttarget\treason\n" > "${report}.tmp"

    local rewritten=0 already=0 unresolved=0 pending=0
    local link target mapped linkdir rel

    while IFS= read -r link; do
        target=$(readlink "$link") || continue

        # already relative — portable, leave it alone (this is what makes the
        # pass idempotent)
        if [[ "$target" != /* ]]; then
            already=$(( already + 1 ))
            continue
        fi

        if [[ "$target" == "${REMOTE_BASE}/"* ]]; then
            mapped="${LOCAL_BASE}/${target#${REMOTE_BASE}/}"
            linkdir=$(cd "$(dirname "$link")" && pwd)
            rel=$(relpath "$linkdir" "$mapped")

            if [[ "$DRY_RUN" == "true" ]]; then
                echo "    [would relink] ${link#$LOCAL_BASE/} -> $rel"
            else
                ln -sfn "$rel" "$link"
            fi
            rewritten=$(( rewritten + 1 ))

            # a rewritten link may still not resolve yet: its target can be
            # deferred to phase 2, or live in a directory not listed in
            # DIRECTORIES. Record it so it is visible either way.
            if [[ ! -e "$mapped" ]]; then
                pending=$(( pending + 1 ))
                printf "%s\t%s\ttarget-not-present-yet\n" \
                    "${link#$LOCAL_BASE/}" "$target" >> "${report}.tmp"
            fi
        else
            # outside REMOTE_BASE — nothing here can map it, so leave the link
            # exactly as it is and report it rather than papering over it
            unresolved=$(( unresolved + 1 ))
            printf "%s\t%s\toutside-REMOTE_BASE\n" \
                "${link#$LOCAL_BASE/}" "$target" >> "${report}.tmp"
        fi
    done < <(find "$root" -type l 2>/dev/null)

    mv "${report}.tmp" "$report"

    local verb="rewritten"
    [[ "$DRY_RUN" == "true" ]] && verb="would rewrite"
    echo "  [relink] ${dirname}: ${rewritten} ${verb}, ${already} already relative, ${pending} awaiting data, ${unresolved} unresolvable"
    if (( unresolved > 0 || pending > 0 )); then
        echo "  [relink] see $report"
    fi
}

relink_all() {
    local phase="$1"; shift
    echo ""
    echo "--- relink after ${phase} ---"
    local d
    for d in "$@"; do
        relink_dir "$d"
    done
}

# =============================================================================
# Manifest helpers
#
# Schema (TSV): idx  name  path  type  status  skipped
# One row per directory in direct mode.
#
# Status flow:
#   PENDING → UPLOADING_MAIN → MAIN_DONE → UPLOADING_DEFERRED → UPLOADED
#                  ↓                             ↓
#             FAILED_MAIN                  FAILED_DEFERRED
# =============================================================================

manifest_path() {
    echo "${MANIFEST_DIR}/manifest_${1}.tsv"
}

manifest_init() {
    local manifest="$1"
    if [[ ! -f "$manifest" ]]; then
        printf "idx\tname\tpath\ttype\tstatus\tskipped\n" > "$manifest"
    fi
}

manifest_has() {
    local manifest="$1" name="$2"
    [[ -f "$manifest" ]] || return 1
    awk -F'\t' -v n="$name" 'NR>1 && $2==n {f=1} END{exit !f}' "$manifest"
}

manifest_add() {
    local manifest="$1" name="$2" path="$3" type="$4"
    # exact match on the name column — a substring match would false-positive
    # against longer names and against the path column
    manifest_has "$manifest" "$name" && return 0
    local idx
    idx=$(( $(wc -l < "$manifest") ))
    printf "%d\t%s\t%s\t%s\tPENDING\t-\n" "$idx" "$name" "$path" "$type" >> "$manifest"
}

manifest_update() {
    local manifest="$1" name="$2" status="$3"
    awk -v name="$name" -v status="$status" \
        'BEGIN{FS=OFS="\t"} $2==name{$5=status} {print}' \
        "$manifest" > "${manifest}.tmp" && mv "${manifest}.tmp" "$manifest"
}

manifest_set_skipped() {
    local manifest="$1" name="$2" skipped="$3"
    awk -v name="$name" -v s="$skipped" \
        'BEGIN{FS=OFS="\t"} $2==name{$6=s} {print}' \
        "$manifest" > "${manifest}.tmp" && mv "${manifest}.tmp" "$manifest"
}

manifest_status() {
    local manifest="$1" name="$2"
    awk -v name="$name" 'BEGIN{FS="\t"} $2==name{print $5}' "$manifest"
}

manifest_skipped() {
    local manifest="$1" name="$2"
    awk -v name="$name" 'BEGIN{FS="\t"} $2==name{print $6}' "$manifest"
}

# =============================================================================
# Rsync passes
#
#   first    — RSYNC_PRIORITY_FIRST only
#   main     — everything except RSYNC_PRIORITY_LAST
#   deferred — RSYNC_PRIORITY_LAST only
#
# Pass `main` deliberately re-offers the priority-first files: -a makes that a
# cheap no-op, and it means a failure during `first` is still recoverable.
# =============================================================================

rsync_dir_pass() {
    local dirname="$1" pass="$2"

    local local_log_dir="${LOG_DIR}/${dirname}"
    mkdir -p "$local_log_dir" "${LOCAL_BASE}/${dirname}"
    local logfile="${local_log_dir}/rsync_${dirname}_${pass}_$(date +%Y%m%d_%H%M%S).log"

    build_always_excludes

    local filter_args=()
    local pat
    case "$pass" in
        first)
            for pat in "${RSYNC_PRIORITY_FIRST[@]:-}"; do
                [[ -n "$pat" ]] && filter_args+=(--include="$pat")
            done
            (( ${#filter_args[@]} == 0 )) && return 0
            # --include="*/" so rsync descends; --exclude="*" drops everything else
            filter_args=(--include="*/" "${filter_args[@]}" --exclude="*")
            ;;
        main)
            for pat in "${RSYNC_PRIORITY_LAST[@]:-}"; do
                [[ -n "$pat" ]] && filter_args+=(--exclude="$pat")
            done
            ;;
        deferred)
            for pat in "${RSYNC_PRIORITY_LAST[@]:-}"; do
                [[ -n "$pat" ]] && filter_args+=(--include="$pat")
            done
            (( ${#filter_args[@]} == 0 )) && return 0
            filter_args=(--include="*/" "${filter_args[@]}" --exclude="*")
            ;;
        *)
            echo "ERROR: unknown rsync pass '$pass'"; return 1 ;;
    esac

    echo "  [rsync:${pass}] log → $logfile"
    # rsync is first-match-wins, so the unconditional excludes go first: they
    # then beat any --include a pass adds later.
    robust_rsync --logfile "$logfile" \
        "${ALWAYS_EXCLUDES[@]:+${ALWAYS_EXCLUDES[@]}}" \
        "${filter_args[@]:+${filter_args[@]}}" \
        "${HPC}:${REMOTE_BASE}/${dirname}/" \
        "${LOCAL_BASE}/${dirname}/"
}

# =============================================================================
# Phase 1 — main transfer for one directory
# =============================================================================

sync_main() {
    local dirname="$1"
    local manifest
    manifest=$(manifest_path "$dirname")
    manifest_init "$manifest"
    manifest_add "$manifest" "$dirname" "$dirname" "dir"

    local cur_status
    cur_status=$(manifest_status "$manifest" "$dirname")
    if [[ "$FORCE" != "true" ]]; then
        case "$cur_status" in
            MAIN_DONE|UPLOADING_DEFERRED|UPLOADED)
                echo "  [skip] $dirname main phase already done ($cur_status)"
                return 0
                ;;
        esac
    fi

    manifest_update "$manifest" "$dirname" "UPLOADING_MAIN"
    manifest_set_skipped "$manifest" "$dirname" "$(skipped_patterns_field)"

    if rsync_dir_pass "$dirname" first && rsync_dir_pass "$dirname" main; then
        manifest_update "$manifest" "$dirname" "MAIN_DONE"
        echo "  [main done] $dirname"
        return 0
    else
        manifest_update "$manifest" "$dirname" "FAILED_MAIN"
        echo "  [failed main] $dirname"
        return 1
    fi
}

# =============================================================================
# Phase 2 — deferred (trajectory) transfer for one directory
# =============================================================================

sync_deferred() {
    local dirname="$1"
    local manifest
    manifest=$(manifest_path "$dirname")
    manifest_init "$manifest"

    local cur_status
    cur_status=$(manifest_status "$manifest" "$dirname")

    if [[ "$FORCE" != "true" ]]; then
        case "$cur_status" in
            MAIN_DONE|FAILED_DEFERRED) ;;
            UPLOADED)
                echo "  [skip] $dirname already UPLOADED"
                return 0
                ;;
            *)
                # never completed its main phase — not eligible for deferred
                echo "  [skip] $dirname main phase incomplete ($cur_status)"
                return 0
                ;;
        esac
    fi

    if ! has_deferred_patterns; then
        manifest_update "$manifest" "$dirname" "UPLOADED"
        echo "  [uploaded] $dirname (no deferred patterns configured)"
        return 0
    fi

    manifest_update "$manifest" "$dirname" "UPLOADING_DEFERRED"

    if rsync_dir_pass "$dirname" deferred; then
        manifest_update "$manifest" "$dirname" "UPLOADED"
        # the deferred files are down now, so whatever the main phase skipped is
        # no longer outstanding — don't leave a stale pattern list on an
        # UPLOADED row
        manifest_set_skipped "$manifest" "$dirname" "$(skipped_patterns_field)"
        echo "  [uploaded] $dirname"
        return 0
    else
        manifest_update "$manifest" "$dirname" "FAILED_DEFERRED"
        echo "  [failed deferred] $dirname"
        return 1
    fi
}

# =============================================================================
# Resolve mode for a directory entry
# =============================================================================

# Entry syntax "dirname" or "dirname:mode" is kept so archive mode can be
# re-enabled as a dispatch change; see ARCHIVE_MODE_NOTES.md.
resolve_entry() {
    local entry="$1"
    ENTRY_DIRNAME="${entry%%:*}"
    ENTRY_MODE="${entry##*:}"
    if [[ "$ENTRY_DIRNAME" == "$ENTRY_MODE" ]]; then
        ENTRY_MODE="$DEFAULT_MODE"
    fi
}

check_mode() {
    local dirname="$1" mode="$2"
    if [[ "$mode" != "direct" ]]; then
        echo "ERROR: unsupported mode '$mode' for ${dirname}."
        echo "       Only 'direct' is supported. Archive mode was removed; see"
        echo "       ${SCRIPT_DIR}/ARCHIVE_MODE_NOTES.md to bring it back."
        return 1
    fi
    return 0
}

# =============================================================================
# Subcommands
# =============================================================================

cmd_sync() {
    local target=""
    for arg in "$@"; do
        case "$arg" in
            --force)                FORCE="true" ;;
            --skip-trajectory)      SKIP_TRAJECTORY="true" ;;
            --no-skip-trajectory)   SKIP_TRAJECTORY="false" ;;
            --dry-run)              DRY_RUN="true" ;;
            --no-relink)            RELINK="false" ;;
            -*) echo "ERROR: unknown flag '$arg'"; usage ;;
            *)  target="$arg" ;;
        esac
    done

    preflight_local
    setup_ssh_control
    trap cleanup_ssh_control EXIT
    preflight_remote

    [[ "$DRY_RUN" == "true" ]] && echo "*** DRY RUN — rsync runs with --dry-run, nothing is written ***"
    [[ "$SKIP_TRAJECTORY" == "true" ]] && echo "*** --skip-trajectory active: not syncing ${RSYNC_PRIORITY_LAST[*]:-} ***"
    if (( ${#ALWAYS_EXCLUDE[@]} > 0 )); then
        echo "*** never transferred (ALWAYS_EXCLUDE): ${ALWAYS_EXCLUDE[*]:-} ***"
    fi

    # collect the directories this run applies to
    local -a run_dirs=()
    local entry
    for entry in "${DIRECTORIES[@]}"; do
        resolve_entry "$entry"
        if [[ -n "$target" && "$ENTRY_DIRNAME" != "$target" ]]; then continue; fi
        check_mode "$ENTRY_DIRNAME" "$ENTRY_MODE" || continue
        run_dirs+=("$ENTRY_DIRNAME")
    done

    if (( ${#run_dirs[@]} == 0 )); then
        echo "ERROR: no directories to sync${target:+ matching '$target'}"
        exit 1
    fi

    local d
    local main_failed=0 deferred_failed=0

    echo ""
    echo "############################################################"
    echo "# PHASE 1/2 — main transfer (priority + everything deferred-excluded)"
    echo "############################################################"
    for d in "${run_dirs[@]}"; do
        echo ""
        echo "=== ${d} [main] ==="
        sync_main "$d" || main_failed=$(( main_failed + 1 ))
    done

    # Relink after the phase, not after each directory: a link in one directory
    # often points into another, so waiting until the whole phase is down means
    # far fewer targets are still missing.
    if [[ "$RELINK" == "true" ]]; then
        relink_all "phase 1 (main)" "${run_dirs[@]}"
    fi

    echo ""
    echo "############################################################"
    if [[ "$SKIP_TRAJECTORY" == "true" ]]; then
        echo "# PHASE 2/2 — deferred transfer SKIPPED (--skip-trajectory)"
        echo "############################################################"
        echo ""
        echo "Directories remain at MAIN_DONE. Re-run without --skip-trajectory"
        echo "to pull the deferred files and complete them."
    else
        echo "# PHASE 2/2 — deferred transfer (${RSYNC_PRIORITY_LAST[*]:-none})"
        echo "############################################################"
        for d in "${run_dirs[@]}"; do
            echo ""
            echo "=== ${d} [deferred] ==="
            sync_deferred "$d" || deferred_failed=$(( deferred_failed + 1 ))
        done
        if [[ "$RELINK" == "true" ]]; then
            relink_all "phase 2 (deferred)" "${run_dirs[@]}"
        fi
    fi

    echo ""
    echo "=== sync complete ==="
    if (( main_failed > 0 || deferred_failed > 0 )); then
        echo "    ${main_failed} directory(ies) failed the main phase, ${deferred_failed} failed the deferred phase."
        echo "    Run './nas.sh status' for detail; re-running resumes where it left off."
        return 1
    fi
}

cmd_status() {
    printf "%-30s %-8s %-20s %s\n" "DIRECTORY" "MODE" "STATUS" "SKIPPED"
    printf '%.0s-' {1..90}; echo

    local entry
    for entry in "${DIRECTORIES[@]}"; do
        resolve_entry "$entry"
        local manifest
        manifest=$(manifest_path "$ENTRY_DIRNAME")
        if [[ ! -f "$manifest" ]]; then
            printf "%-30s %-8s %-20s %s\n" "$ENTRY_DIRNAME" "$ENTRY_MODE" "(no manifest)" "-"
            continue
        fi
        local st sk
        st=$(manifest_status "$manifest" "$ENTRY_DIRNAME")
        sk=$(manifest_skipped "$manifest" "$ENTRY_DIRNAME")
        printf "%-30s %-8s %-20s %s\n" \
            "$ENTRY_DIRNAME" "$ENTRY_MODE" "${st:-PENDING}" "${sk:--}"
    done
}

cmd_relink() {
    local target=""
    for arg in "$@"; do
        case "$arg" in
            --dry-run) DRY_RUN="true" ;;
            -*) echo "ERROR: unknown flag '$arg'"; usage ;;
            *)  target="$arg" ;;
        esac
    done

    [[ "$DRY_RUN" == "true" ]] && echo "*** DRY RUN — showing rewrites, changing nothing ***"

    local -a run_dirs=()
    local entry
    for entry in "${DIRECTORIES[@]}"; do
        resolve_entry "$entry"
        if [[ -n "$target" && "$ENTRY_DIRNAME" != "$target" ]]; then continue; fi
        run_dirs+=("$ENTRY_DIRNAME")
    done

    if (( ${#run_dirs[@]} == 0 )); then
        echo "ERROR: no directories to relink${target:+ matching '$target'}"
        exit 1
    fi

    relink_all "manual run" "${run_dirs[@]}"
}

cmd_manifest() {
    local dirname=""
    for arg in "$@"; do
        if [[ "$arg" == "--force" ]]; then
            FORCE="true"
        else
            dirname="$arg"
        fi
    done
    [[ -z "$dirname" ]] && { echo "Usage: ./nas.sh manifest [--force] <dirname>"; exit 1; }

    local mode="$DEFAULT_MODE"
    local entry
    for entry in "${DIRECTORIES[@]}"; do
        resolve_entry "$entry"
        if [[ "$ENTRY_DIRNAME" == "$dirname" ]]; then
            mode="$ENTRY_MODE"
            break
        fi
    done
    check_mode "$dirname" "$mode" || exit 1

    local manifest
    manifest=$(manifest_path "$dirname")
    if [[ "$FORCE" == "true" && -f "$manifest" ]]; then
        echo "  [force] removing existing manifest: $manifest"
        rm -f "$manifest"
    fi
    manifest_init "$manifest"
    manifest_add "$manifest" "$dirname" "$dirname" "dir"
    echo "Manifest written: $manifest"
}

# =============================================================================
# Entry point
# =============================================================================

usage() {
    echo "Usage:"
    echo "  ./nas.sh sync [flags] [dirname]        sync all or one directory"
    echo "  ./nas.sh status                        show per-directory status"
    echo "  ./nas.sh manifest [--force] <dirname>  create/reset a manifest
  ./nas.sh relink [--dry-run] [dirname]  repair remote-namespace symlinks"
    echo ""
    echo "Sync flags:"
    echo "  --force                 ignore recorded status and re-run transfers"
    echo "  --skip-trajectory       do not sync the RSYNC_PRIORITY_LAST globs at all"
    echo "  --no-skip-trajectory    force the deferred phase on (overrides nas.config)"
    echo "  --dry-run               pass --dry-run to rsync; transfer nothing
  --no-relink             skip the post-phase symlink repair"
    echo ""
    echo "Transfers run in two global phases: every directory completes its main"
    echo "transfer before any directory starts its deferred (trajectory) transfer."
    exit 1
}

case "${1:-}" in
    sync)     shift; cmd_sync "$@" ;;
    relink)   shift; cmd_relink "$@" ;;
    status)   cmd_status ;;
    manifest) shift; cmd_manifest "$@" ;;
    *)        usage ;;
esac
