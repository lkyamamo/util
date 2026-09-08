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
# Archive mode (SLURM-side compression) was removed; see ARCHIVE_MODE_NOTES.md.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSYNC=/opt/homebrew/bin/rsync
CONFIG="${SCRIPT_DIR}/nas.config"
MANIFEST_DIR="${SCRIPT_DIR}/manifests"
LOG_DIR="${SCRIPT_DIR}/logs"
SSH_CONTROL="${TMPDIR:-/tmp}/nas_ssh_ctl_$$"

# =============================================================================
# Load config
# =============================================================================

[[ -f "$CONFIG" ]] || { echo "ERROR: nas.config not found at $CONFIG"; exit 1; }
source "$CONFIG"

HPC="${HPC_USER}@${HPC_HOST}"
FORCE="false"
DRY_RUN="false"
# SKIP_TRAJECTORY comes from nas.config; --skip-trajectory/--no-skip-trajectory override it
SKIP_TRAJECTORY="${SKIP_TRAJECTORY:-false}"

mkdir -p "$MANIFEST_DIR" "$LOG_DIR"

# =============================================================================
# Preflight
# =============================================================================

preflight_local() {
    if [[ ! -x "$RSYNC" ]]; then
        echo "ERROR: rsync not found at $RSYNC"
        echo "       The macOS system rsync is too old for --info=progress2."
        echo "       Install a current one:  brew install rsync"
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

# Populates SKIP_EXCLUDES with --exclude args for every SKIP_PATTERNS entry,
# but only when trajectory skipping is active. Applied to *every* pass, so a
# skip pattern that is not also a deferred pattern is still honored.
SKIP_EXCLUDES=()
build_skip_excludes() {
    SKIP_EXCLUDES=()
    [[ "$SKIP_TRAJECTORY" != "true" ]] && return 0
    local pat
    for pat in "${SKIP_PATTERNS[@]:-}"; do
        [[ -n "$pat" ]] && SKIP_EXCLUDES+=(--exclude="$pat")
    done
}

# Comma-joined list of patterns actually suppressed on this run, for the
# manifest's `skipped` column. "-" when nothing was suppressed.
skipped_patterns_field() {
    if [[ "$SKIP_TRAJECTORY" != "true" ]]; then echo "-"; return; fi
    local out="" pat
    for pat in "${SKIP_PATTERNS[@]:-}"; do
        [[ -z "$pat" ]] && continue
        out="${out:+${out},}${pat}"
    done
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

    build_skip_excludes

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
    robust_rsync --logfile "$logfile" \
        "${SKIP_EXCLUDES[@]:+${SKIP_EXCLUDES[@]}}" \
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
            -*) echo "ERROR: unknown flag '$arg'"; usage ;;
            *)  target="$arg" ;;
        esac
    done

    preflight_local
    setup_ssh_control
    trap cleanup_ssh_control EXIT
    preflight_remote

    [[ "$DRY_RUN" == "true" ]] && echo "*** DRY RUN — rsync runs with --dry-run, nothing is written ***"
    [[ "$SKIP_TRAJECTORY" == "true" ]] && echo "*** --skip-trajectory active: excluding ${SKIP_PATTERNS[*]:-} ***"

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
    echo "  ./nas.sh manifest [--force] <dirname>  create/reset a manifest"
    echo ""
    echo "Sync flags:"
    echo "  --force                 ignore recorded status and re-run transfers"
    echo "  --skip-trajectory       exclude SKIP_PATTERNS and skip the deferred phase"
    echo "  --no-skip-trajectory    force the deferred phase on (overrides nas.config)"
    echo "  --dry-run               pass --dry-run to rsync; transfer nothing"
    echo ""
    echo "Transfers run in two global phases: every directory completes its main"
    echo "transfer before any directory starts its deferred (trajectory) transfer."
    exit 1
}

case "${1:-}" in
    sync)     shift; cmd_sync "$@" ;;
    status)   cmd_status ;;
    manifest) shift; cmd_manifest "$@" ;;
    *)        usage ;;
esac
