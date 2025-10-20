#!/usr/bin/env bash
# sync-from-remote-list.sh — rsync directories *and* tarballs into the current directory
# Remote base: lkyamamo@endeavour.usc.edu:/project/priyav_216/Tians/<name>
# Local dest:  ./<name>   (dir)  or  ./<name> (file)
# Logs:        update_log.txt (truncated each run)
#
# Toggle console logging:
#   LOG_TO_CONSOLE=true  # log to console + file
#   LOG_TO_CONSOLE=false # log to file only

set -uo pipefail

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
  source "${SCRIPT_DIR}/config.sh"
fi

LIST_FILE="${LIST_FILE:-directory_list.txt}"
REMOTE="$HPC"
REMOTE_BASE="${REMOTE_BASE:-/project/priyav_216/Tians}"

# Use robust SSH configuration from config.sh
RSYNC_PATH="${RSYNC_PATH:-/usr/bin/rsync}"
RSYNC_SSH=(-e "$(get_rsync_ssh)")
RSYNC_EXTRA=(--timeout="${RSYNC_TIMEOUT:-600}")  # common extras

# IMPORTANT: Do NOT mix --inplace with --partial-dir.
RSYNC_DIR_EXTRA=(--partial-dir=.rsync-partial)   # for directories only
RSYNC_FILE_EXTRA=(--inplace --no-whole-file)     # for single large files (tarballs)

# Optional: enable compression over slow links by adding: RSYNC_EXTRA+=( -z )

# --- Logging ---
LOG="${LOG_FILE:-update_log.txt}"
: > "$LOG"

log() {
  if $LOG_TO_CONSOLE; then
    printf '%s\n' "$*" | tee -a "$LOG" >/dev/null
  else
    printf '%s\n' "$*" >>"$LOG"
  fi
}

run() {
  if $LOG_TO_CONSOLE; then
    "$@" 2>&1 | tee -a "$LOG"
  else
    "$@" >>"$LOG" 2>&1
  fi
}

[[ -f "$LIST_FILE" ]] || { log "Error: '$LIST_FILE' not found in $(pwd)"; exit 1; }

# Test SSH connection first
if ! test_ssh_connection; then
  log "Error: Cannot establish SSH connection to HPC"
  log "Please check your SSH configuration and network connectivity"
  exit 1
fi

# Preflight: confirm remote rsync exists & reachable
log "Checking remote rsync availability..."
if ! robust_ssh "command -v '$RSYNC_PATH' >/dev/null"; then
  log "Error: cannot find '$RSYNC_PATH' on $REMOTE"
  log "Please ensure rsync is installed on the remote system"
  exit 1
fi

processed=0 ok=0 warned=0 skipped=0 failed=0

rsync_with_retry() {
  local src="$1" dst="$2" type="$3"; shift 3
  local attempt=1 max=3 rc
  local extra=("$@")
  while (( attempt <= max )); do
    if run rsync -avP "${RSYNC_SSH[@]}" "${RSYNC_EXTRA[@]}" "${extra[@]}" \
           --rsync-path="$RSYNC_PATH" -- "$src" "$dst" 2>/dev/null; then
      return 0
    fi
    rc=$?
    # Treat partial/vanished as warning and continue
    if [[ $rc -eq 23 || $rc -eq 24 ]]; then
      log "   [warn] $type partial transfer (rsync exit $rc). Continuing..."
      return 0
    fi
    if (( attempt < max )); then
      log "   [retry] $type rsync failed (exit $rc). Retrying in $((attempt*5))s... [$attempt/$max]"
      sleep $((attempt*5))
    fi
    ((attempt++))
  done
  # Log the final failure but don't exit
  log "   [fail] $type rsync failed after $max attempts (exit $rc). Continuing with next item..."
  return "$rc"
}

# Read the list on FD 3 so ssh/rsync never consume it
exec 3<"$LIST_FILE"
while IFS= read -r raw <&3 || [[ -n "${raw:-}" ]]; do
  line="${raw%%#*}"
  line="${line%$'\r'}"
  line="${line#"${line%%[![:space:]]*}"}"
  line="${line%"${line##*[![:space:]]}"}"
  [[ -z "$line" ]] && continue

  name="$line"
  ((processed++)) || true

  remote_path="$REMOTE_BASE/$name"
  dest="$PWD/$name"

  log "==> [$processed] $REMOTE:$remote_path -> $dest"

  if robust_ssh "[ -d \"$remote_path\" ]" 2>/dev/null; then
    # Directory: use partial-dir for efficient resumes across many files
    mkdir -p "$dest"
    mkdir -p "$dest/.rsync-partial"
    if rsync_with_retry "$REMOTE:$remote_path/" "$dest/" "dir" "${RSYNC_DIR_EXTRA[@]}"; then
      ((ok++)) || true
    else
      ((failed++)) || true
    fi

  elif robust_ssh "[ -f \"$remote_path\" ]" 2>/dev/null; then
    # File (tarball): use in-place delta updates
    mkdir -p "$(dirname -- "$dest")"
    if rsync_with_retry "$REMOTE:$remote_path" "$dest" "file" "${RSYNC_FILE_EXTRA[@]}"; then
      ((ok++)) || true
    else
      ((failed++)) || true
    fi

  else
    log "   [skip] remote path not found (neither dir nor file): $REMOTE:$remote_path"
    ((skipped++)) || true
    continue
  fi
done
exec 3<&-

log ""
log "--- Summary ---"
log "Processed: $processed"
log "OK:        $ok"
log "Warnings:  $warned   (partial transfers)"
log "Skipped:   $skipped  (missing/inaccessible)"
log "Failed:    $failed   (non-transient errors)"
