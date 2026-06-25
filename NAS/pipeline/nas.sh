#!/usr/bin/env bash
# =============================================================================
# nas.sh — central NAS transfer script
# Usage:
#   ./nas.sh sync [dirname]     process all or one directory
#   ./nas.sh status             show manifest summary
#   ./nas.sh manifest <dirname> rebuild manifest for a directory
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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
REMOTE_LOG_BASE="${REMOTE_BASE}/logs"

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
    # Usage: robust_rsync [--logfile <path>] <src> <dest>
    local logfile=""
    if [[ "${1:-}" == "--logfile" ]]; then
        logfile="$2"
        shift 2
    fi

    local attempt=0 delay=10
    while (( attempt < 5 )); do
        if [[ -n "$logfile" ]]; then
            rsync \
                -az \
                --partial \
                --partial-dir=.rsync-partial \
                --info=progress2 \
                -e "ssh -o ControlPath=$SSH_CONTROL -o ControlMaster=no" \
                "$@" 2>&1 | tee -a "$logfile"
            local rc="${PIPESTATUS[0]}"
        else
            rsync \
                -az \
                --partial \
                --partial-dir=.rsync-partial \
                --info=progress2 \
                -e "ssh -o ControlPath=$SSH_CONTROL -o ControlMaster=no" \
                "$@"
            local rc=$?
        fi
        [[ $rc -eq 0 ]] && return 0
        attempt=$(( attempt + 1 ))
        echo "  [rsync] attempt $attempt failed, retrying in ${delay}s..."
        sleep "$delay"
        delay=$(( delay * 2 ))
    done
    return 1
}

# =============================================================================
# Manifest helpers
# =============================================================================

manifest_path() {
    echo "${MANIFEST_DIR}/manifest_${1}.tsv"
}

manifest_init() {
    local manifest="$1"
    if [[ ! -f "$manifest" ]]; then
        printf "idx\tname\tpath\ttype\tstatus\n" > "$manifest"
    fi
}

manifest_add() {
    local manifest="$1" name="$2" path="$3" type="$4"
    # skip if already tracked
    if grep -qF "$name" "$manifest" 2>/dev/null; then return; fi
    local idx
    idx=$(( $(wc -l < "$manifest") ))
    printf "%d\t%s\t%s\t%s\tPENDING\n" "$idx" "$name" "$path" "$type" >> "$manifest"
}

manifest_update() {
    local manifest="$1" name="$2" status="$3"
    # use awk to rewrite the status column for the matching name
    awk -v name="$name" -v status="$status" \
        'BEGIN{FS=OFS="\t"} $2==name{$5=status} {print}' \
        "$manifest" > "${manifest}.tmp" && mv "${manifest}.tmp" "$manifest"
}

manifest_status() {
    local manifest="$1" name="$2"
    awk -v name="$name" 'BEGIN{FS="\t"} $2==name{print $5}' "$manifest"
}

manifest_get_field() {
    local manifest="$1" name="$2" field="$3"
    awk -v name="$name" -v f="$field" 'BEGIN{FS="\t"} $2==name{print $f}' "$manifest"
}

manifest_rows_by_status() {
    local manifest="$1" status="$2"
    awk -v s="$status" 'BEGIN{FS="\t"} NR>1 && $5==s {print $2}' "$manifest"
}

# =============================================================================
# Size helpers
# =============================================================================

# Convert human size (100G, 10M, etc.) to bytes
to_bytes() {
    local val="${1%[KMGTPE]*}" unit="${1##*[0-9]}"
    case "${unit^^}" in
        K) echo $(( val * 1024 )) ;;
        M) echo $(( val * 1024 * 1024 )) ;;
        G) echo $(( val * 1024 * 1024 * 1024 )) ;;
        T) echo $(( val * 1024 * 1024 * 1024 * 1024 )) ;;
        *) echo "$val" ;;
    esac
}

MAX_ARCHIVE_BYTES=$(to_bytes "$MAX_ARCHIVE_SIZE")
MIN_FILE_BYTES=$(to_bytes "$MIN_FILE_SIZE")

remote_dir_size_bytes() {
    # returns size in bytes of a remote directory
    robust_ssh "du -sb '$1' 2>/dev/null | cut -f1" || echo 0
}

remote_file_size_bytes() {
    robust_ssh "stat -c%s '$1' 2>/dev/null || stat -f%z '$1' 2>/dev/null || echo 0"
}

# =============================================================================
# Tree walk — find compress targets
# =============================================================================

# Populates global arrays:
#   TARGETS_NAME   — tarball/file name (e.g. child1.tar.zst)
#   TARGETS_PATH   — relative path of the source within REMOTE_BASE/dirname
#   TARGETS_TYPE   — dir | file | bundle
declare -a TARGETS_NAME TARGETS_PATH TARGETS_TYPE

walk_tree() {
    local remote_dir="$1"   # absolute remote path of dir to walk
    local rel_parent="$2"   # relative path from dirname root to parent of remote_dir
    local dirname="$3"      # top-level directory name (for manifest)

    # list direct children, separated by null
    local listing
    listing=$(robust_ssh "find '$remote_dir' -mindepth 1 -maxdepth 1 -print0 2>/dev/null")

    local has_files=false
    local -a direct_files direct_dirs
    direct_files=()
    direct_dirs=()

    while IFS= read -r -d '' item; do
        local item_type
        item_type=$(robust_ssh "[ -d '$item' ] && echo dir || echo file")
        if [[ "$item_type" == "dir" ]]; then
            direct_dirs+=("$item")
        else
            has_files=true
            direct_files+=("$item")
        fi
    done <<< "$listing"

    local base
    base=$(basename "$remote_dir")
    local rel_self="${rel_parent:+${rel_parent}/}${base}"

    if $has_files; then
        local dir_size
        dir_size=$(remote_dir_size_bytes "$remote_dir")

        if (( dir_size <= MAX_ARCHIVE_BYTES )); then
            # whole dir as one tar.zst
            TARGETS_NAME+=("${base}.tar.zst")
            TARGETS_PATH+=("$rel_self")
            TARGETS_TYPE+=("dir")
        else
            # dir too large — handle files individually and recurse into subdirs
            local -a small_files=()
            for f in "${direct_files[@]}"; do
                local fname fsize
                fname=$(basename "$f")
                fsize=$(remote_file_size_bytes "$f")
                if (( fsize >= MIN_FILE_BYTES )); then
                    TARGETS_NAME+=("${fname}.zst")
                    TARGETS_PATH+=("$rel_self")
                    TARGETS_TYPE+=("file")
                else
                    small_files+=("$fname")
                fi
            done
            if (( ${#small_files[@]} > 0 )); then
                TARGETS_NAME+=("${base}_bundle.tar.zst")
                TARGETS_PATH+=("$rel_self")
                TARGETS_TYPE+=("bundle")
            fi
            for subdir in "${direct_dirs[@]}"; do
                walk_tree "$subdir" "$rel_self" "$dirname"
            done
        fi
    else
        # no files here — recurse into subdirs
        for subdir in "${direct_dirs[@]}"; do
            walk_tree "$subdir" "$rel_parent" "$dirname"
        done
    fi
}

# =============================================================================
# SLURM job generation and submission
# =============================================================================

generate_slurm_script() {
    local dirname="$1"
    local target_name="$2"   # e.g. child1.tar.zst or file.zst
    local target_path="$3"   # relative path within dirname (parent of target)
    local target_type="$4"   # dir | file | bundle
    local job_name="nas_${dirname}_${target_name%%.*}"
    local remote_source="${REMOTE_BASE}/${dirname}/${target_path}"
    local remote_log_dir="${REMOTE_LOG_BASE}/${dirname}"
    local remote_tarball="${remote_log_dir}/${target_name}"

    local compress_cmd
    case "$target_type" in
        dir)
            local src_base
            src_base=$(basename "$target_path")
            local src_parent
            src_parent=$(dirname "${REMOTE_BASE}/${dirname}/${target_path}")
            compress_cmd="tar --use-compress-program=\"zstd -T${COMPRESSION_THREADS} -${COMPRESSION_LEVEL}\" -cf \"${remote_tarball}\" -C \"${src_parent}\" \"${src_base}\""
            ;;
        file)
            local fname="${target_name%.zst}"
            compress_cmd="zstd -T${COMPRESSION_THREADS} -${COMPRESSION_LEVEL} \"${remote_source}/${fname}\" -o \"${remote_tarball}\""
            ;;
        bundle)
            # compress all small files in the directory into one tar.zst
            compress_cmd="tar --use-compress-program=\"zstd -T${COMPRESSION_THREADS} -${COMPRESSION_LEVEL}\" -cf \"${remote_tarball}\" -C \"${remote_source}\" \$(find . -maxdepth 1 -type f)"
            ;;
    esac

    cat <<SLURM
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${COMPRESSION_THREADS}
#SBATCH --time=${SLURM_TIME}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --output=${remote_log_dir}/${target_name}_%j.log

module load zstd/1.5.6

mkdir -p "${remote_log_dir}"

${compress_cmd}

if [[ \$? -eq 0 ]]; then
    echo "NAS_SUCCESS"
else
    echo "NAS_FAILED"
    exit 1
fi
SLURM
}

submit_slurm_job() {
    local dirname="$1" target_name="$2" target_path="$3" target_type="$4"

    local script
    script=$(generate_slurm_script "$dirname" "$target_name" "$target_path" "$target_type")

    local remote_script="/tmp/nas_job_${dirname}_${target_name%%.*}_$$.sh"
    echo "$script" | robust_ssh "cat > '$remote_script'"

    local jobid
    jobid=$(robust_ssh "sbatch '$remote_script' | awk '{print \$NF}'")
    robust_ssh "rm -f '$remote_script'"

    echo "$jobid"
}

poll_job() {
    # returns: RUNNING | COMPLETED | FAILED
    local jobid="$1"
    local state
    state=$(robust_ssh "squeue -j '$jobid' -h -o '%T' 2>/dev/null || echo NOTFOUND")
    if [[ -z "$state" || "$state" == "NOTFOUND" ]]; then
        # job no longer in queue — check sacct
        state=$(robust_ssh "sacct -j '$jobid' --format=State --noheader 2>/dev/null | head -1 | tr -d ' '")
    fi
    case "${state^^}" in
        COMPLETED) echo "COMPLETED" ;;
        FAILED|CANCELLED|TIMEOUT|NODE_FAIL|OUT_OF_MEMORY) echo "FAILED" ;;
        *) echo "RUNNING" ;;
    esac
}

check_job_log() {
    # returns: SUCCESS | FAILED
    local dirname="$1" target_name="$2" jobid="$3"
    local log="${REMOTE_LOG_BASE}/${dirname}/${target_name}_${jobid}.log"
    local result
    result=$(robust_ssh "grep -c 'NAS_SUCCESS' '$log' 2>/dev/null || echo 0")
    if [[ "$result" -gt 0 ]]; then echo "SUCCESS"; else echo "FAILED"; fi
}

# =============================================================================
# Rsync helpers
# =============================================================================

rsync_target() {
    local dirname="$1" target_name="$2" target_path="$3" jobid="$4"
    local remote_log_dir="${REMOTE_LOG_BASE}/${dirname}"
    local remote_file="${remote_log_dir}/${target_name}"
    local local_dest="${LOCAL_BASE}/${dirname}/${target_path}"
    local local_log_dir="${LOG_DIR}/${dirname}"
    local logfile="${local_log_dir}/rsync_${target_name}_$(date +%Y%m%d_%H%M%S).log"

    mkdir -p "$local_dest" "$local_log_dir"
    echo "  [rsync log] $logfile"

    robust_rsync --logfile "$logfile" \
        "${HPC}:${remote_file}" \
        "${local_dest}/" && return 0 || return 1
}

cleanup_remote() {
    local dirname="$1" target_name="$2" jobid="$3"
    local remote_log_dir="${REMOTE_LOG_BASE}/${dirname}"
    robust_ssh "rm -f '${remote_log_dir}/${target_name}' '${remote_log_dir}/${target_name}_${jobid}.log'"
    # pull log to local before deleting
    local local_log="${LOG_DIR}/${dirname}"
    mkdir -p "$local_log"
    robust_rsync \
        "${HPC}:${remote_log_dir}/${target_name}_${jobid}.log" \
        "${local_log}/" 2>/dev/null || true
}

# =============================================================================
# Archive mode — process one directory
# =============================================================================

process_archive() {
    local dirname="$1"
    local manifest
    manifest=$(manifest_path "$dirname")
    manifest_init "$manifest"

    echo "  [archive] walking remote tree: ${REMOTE_BASE}/${dirname}"
    TARGETS_NAME=()
    TARGETS_PATH=()
    TARGETS_TYPE=()
    walk_tree "${REMOTE_BASE}/${dirname}" "" "$dirname"

    local total=${#TARGETS_NAME[@]}
    echo "  [archive] found $total compress target(s)"

    # add all targets to manifest
    for (( i=0; i<total; i++ )); do
        manifest_add "$manifest" "${TARGETS_NAME[$i]}" "${TARGETS_PATH[$i]}" "${TARGETS_TYPE[$i]}"
    done

    # sliding window: track active jobs as associative array jobid→target_index
    declare -A active_jobs   # jobid → index into TARGETS arrays
    local -a rsync_queue=()  # indices ready to rsync, in completion order

    local next_target=0

    while true; do
        # submit jobs to fill window
        while (( ${#active_jobs[@]} < MAX_CONCURRENT_JOBS && next_target < total )); do
            local idx=$next_target
            next_target=$(( next_target + 1 ))
            local tname="${TARGETS_NAME[$idx]}"
            local tpath="${TARGETS_PATH[$idx]}"
            local ttype="${TARGETS_TYPE[$idx]}"
            local cur_status
            cur_status=$(manifest_status "$manifest" "$tname")

            if [[ "$cur_status" == "UPLOADED" ]]; then
                echo "  [skip] $tname already UPLOADED"
                continue
            fi

            echo "  [submit] $tname ($ttype)"
            local jobid
            jobid=$(submit_slurm_job "$dirname" "$tname" "$tpath" "$ttype")
            manifest_update "$manifest" "$tname" "COMPRESSING"
            active_jobs["$jobid"]=$idx
            echo "  [submitted] $tname → job $jobid"
        done

        # break if nothing active and nothing pending
        if (( ${#active_jobs[@]} == 0 && next_target >= total && ${#rsync_queue[@]} == 0 )); then
            break
        fi

        # poll active jobs
        for jobid in "${!active_jobs[@]}"; do
            local state
            state=$(poll_job "$jobid")
            local idx="${active_jobs[$jobid]}"
            local tname="${TARGETS_NAME[$idx]}"

            if [[ "$state" == "COMPLETED" ]]; then
                local result
                result=$(check_job_log "$dirname" "$tname" "$jobid")
                unset "active_jobs[$jobid]"
                if [[ "$result" == "SUCCESS" ]]; then
                    manifest_update "$manifest" "$tname" "COMPRESSED"
                    rsync_queue+=("${idx}:${jobid}")
                    echo "  [compressed] $tname → queued for rsync"
                else
                    manifest_update "$manifest" "$tname" "FAILED_COMPRESSION"
                    echo "  [failed compression] $tname (job $jobid)"
                fi
            elif [[ "$state" == "FAILED" ]]; then
                unset "active_jobs[$jobid]"
                manifest_update "$manifest" "$tname" "FAILED_COMPRESSION"
                echo "  [failed compression] $tname (job $jobid)"
            fi
        done

        # process front of rsync queue
        if (( ${#rsync_queue[@]} > 0 )); then
            local entry="${rsync_queue[0]}"
            rsync_queue=("${rsync_queue[@]:1}")
            local idx="${entry%%:*}"
            local jobid="${entry##*:}"
            local tname="${TARGETS_NAME[$idx]}"
            local tpath="${TARGETS_PATH[$idx]}"

            echo "  [rsync] $tname"
            manifest_update "$manifest" "$tname" "UPLOADING"
            if rsync_target "$dirname" "$tname" "$tpath" "$jobid"; then
                cleanup_remote "$dirname" "$tname" "$jobid"
                manifest_update "$manifest" "$tname" "UPLOADED"
                echo "  [uploaded] $tname"
            else
                manifest_update "$manifest" "$tname" "FAILED_UPLOAD"
                echo "  [failed upload] $tname"
            fi
        fi

        sleep 30
    done
}

# =============================================================================
# Direct mode — process one directory
# =============================================================================

process_direct() {
    local dirname="$1"
    local manifest
    manifest=$(manifest_path "$dirname")
    manifest_init "$manifest"

    manifest_add "$manifest" "$dirname" "$dirname" "dir"

    local cur_status
    cur_status=$(manifest_status "$manifest" "$dirname")
    if [[ "$cur_status" == "UPLOADED" ]]; then
        echo "  [skip] $dirname already UPLOADED"
        return 0
    fi

    local local_log_dir="${LOG_DIR}/${dirname}"
    local logfile="${local_log_dir}/rsync_${dirname}_$(date +%Y%m%d_%H%M%S).log"
    mkdir -p "$local_log_dir"

    echo "  [direct] rsyncing ${dirname}..."
    echo "  [rsync log] $logfile"
    manifest_update "$manifest" "$dirname" "UPLOADING"
    if robust_rsync --logfile "$logfile" \
        "${HPC}:${REMOTE_BASE}/${dirname}/" \
        "${LOCAL_BASE}/${dirname}/"; then
        manifest_update "$manifest" "$dirname" "UPLOADED"
        echo "  [uploaded] $dirname"
    else
        manifest_update "$manifest" "$dirname" "FAILED_UPLOAD"
        echo "  [failed upload] $dirname"
        return 1
    fi
}

# =============================================================================
# Resolve mode for a directory entry
# =============================================================================

resolve_entry() {
    local entry="$1"
    ENTRY_DIRNAME="${entry%%:*}"
    ENTRY_MODE="${entry##*:}"
    if [[ "$ENTRY_DIRNAME" == "$ENTRY_MODE" ]]; then
        ENTRY_MODE="$DEFAULT_MODE"
    fi
}

# =============================================================================
# Subcommands
# =============================================================================

cmd_sync() {
    local target="${1:-}"
    setup_ssh_control
    trap cleanup_ssh_control EXIT

    # ensure remote log base exists
    robust_ssh "mkdir -p '${REMOTE_LOG_BASE}'"

    for entry in "${DIRECTORIES[@]}"; do
        resolve_entry "$entry"
        if [[ -n "$target" && "$ENTRY_DIRNAME" != "$target" ]]; then continue; fi

        echo ""
        echo "=== ${ENTRY_DIRNAME} [${ENTRY_MODE}] ==="
        mkdir -p "${LOCAL_BASE}/${ENTRY_DIRNAME}"

        case "$ENTRY_MODE" in
            archive) process_archive "$ENTRY_DIRNAME" ;;
            direct)  process_direct  "$ENTRY_DIRNAME" ;;
            *) echo "ERROR: unknown mode '$ENTRY_MODE' for $ENTRY_DIRNAME"; continue ;;
        esac
    done

    echo ""
    echo "=== sync complete ==="
}

cmd_status() {
    printf "%-30s %-12s %-12s %-12s %-12s %-12s\n" \
        "DIRECTORY" "PENDING" "COMPRESSING" "COMPRESSED" "UPLOADING" "UPLOADED" "FAILED"
    echo "$(printf '%.0s-' {1..100})"

    for entry in "${DIRECTORIES[@]}"; do
        resolve_entry "$entry"
        local manifest
        manifest=$(manifest_path "$ENTRY_DIRNAME")
        if [[ ! -f "$manifest" ]]; then
            printf "%-30s %s\n" "$ENTRY_DIRNAME" "(no manifest)"
            continue
        fi
        local pending compressing compressed uploading uploaded failed
        pending=$(awk 'BEGIN{FS="\t"} NR>1 && $5=="PENDING"{c++} END{print c+0}' "$manifest")
        compressing=$(awk 'BEGIN{FS="\t"} NR>1 && $5=="COMPRESSING"{c++} END{print c+0}' "$manifest")
        compressed=$(awk 'BEGIN{FS="\t"} NR>1 && $5=="COMPRESSED"{c++} END{print c+0}' "$manifest")
        uploading=$(awk 'BEGIN{FS="\t"} NR>1 && $5=="UPLOADING"{c++} END{print c+0}' "$manifest")
        uploaded=$(awk 'BEGIN{FS="\t"} NR>1 && $5=="UPLOADED"{c++} END{print c+0}' "$manifest")
        failed=$(awk 'BEGIN{FS="\t"} NR>1 && ($5=="FAILED_COMPRESSION" || $5=="FAILED_UPLOAD"){c++} END{print c+0}' "$manifest")
        printf "%-30s %-12s %-12s %-12s %-12s %-12s %-12s\n" \
            "$ENTRY_DIRNAME" "$pending" "$compressing" "$compressed" "$uploading" "$uploaded" "$failed"
    done
}

cmd_manifest() {
    local dirname="${1:-}"
    [[ -z "$dirname" ]] && { echo "Usage: ./nas.sh manifest <dirname>"; exit 1; }
    setup_ssh_control
    trap cleanup_ssh_control EXIT

    echo "Rebuilding manifest for $dirname..."
    TARGETS_NAME=()
    TARGETS_PATH=()
    TARGETS_TYPE=()

    # resolve mode to determine if we walk tree or treat as single dir
    local mode="$DEFAULT_MODE"
    for entry in "${DIRECTORIES[@]}"; do
        resolve_entry "$entry"
        if [[ "$ENTRY_DIRNAME" == "$dirname" ]]; then
            mode="$ENTRY_MODE"
            break
        fi
    done

    local manifest
    manifest=$(manifest_path "$dirname")
    manifest_init "$manifest"

    if [[ "$mode" == "archive" ]]; then
        walk_tree "${REMOTE_BASE}/${dirname}" "" "$dirname"
        for (( i=0; i<${#TARGETS_NAME[@]}; i++ )); do
            manifest_add "$manifest" "${TARGETS_NAME[$i]}" "${TARGETS_PATH[$i]}" "${TARGETS_TYPE[$i]}"
        done
        echo "Manifest written: $manifest (${#TARGETS_NAME[@]} targets)"
    else
        manifest_add "$manifest" "$dirname" "$dirname" "dir"
        echo "Manifest written: $manifest (direct mode, 1 entry)"
    fi
}

# =============================================================================
# Entry point
# =============================================================================

usage() {
    echo "Usage:"
    echo "  ./nas.sh sync [dirname]     sync all or one directory"
    echo "  ./nas.sh status             show manifest summary"
    echo "  ./nas.sh manifest <dirname> rebuild manifest for a directory"
    exit 1
}

case "${1:-}" in
    sync)     cmd_sync "${2:-}" ;;
    status)   cmd_status ;;
    manifest) cmd_manifest "${2:-}" ;;
    *)        usage ;;
esac
