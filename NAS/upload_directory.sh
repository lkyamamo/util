#!/usr/bin/env bash
# stripped-down: make a manifest on HPC, then upload items and update status in the manifest

set -uo pipefail
# Note: We don't use 'set -e' because we want to handle individual upload failures
# gracefully and continue processing other items

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
  source "${SCRIPT_DIR}/config.sh"
fi

# Override with environment variables if needed
MANIFEST="${MANIFEST:-${MANIFEST_FILE:-manifest.tsv}}"

# Function to update status in manifest
update_status() {
  local idx="$1" status="$2"
  local temp_manifest="${MANIFEST}.tmp"
  
  # Use || true to prevent script exit if update fails
  if awk -v idx="$idx" -v status="$status" -F'\t' '
    $1 == idx { $4 = status }
    { print }
  ' "$MANIFEST" > "$temp_manifest" 2>/dev/null; then
    mv "$temp_manifest" "$MANIFEST" 2>/dev/null || true
  else
    echo "Warning: Failed to update status for item $idx in manifest" >&2
    rm -f "$temp_manifest" 2>/dev/null || true
  fi
}

# Function to process a single directory (for use by other scripts)
process_single_directory() {
  local remote_dir="$1"
  local local_dest="$2"
  local manifest_file="${3:-manifest.tsv}"
  
  # Create manifest
  if ! "$(dirname "${BASH_SOURCE[0]}")/create_manifest.sh" "$remote_dir" "$manifest_file"; then
    echo "ERROR: Failed to create manifest for $remote_dir" >&2
    return 1
  fi
  
  # Upload directory
  if ! "$0" "$remote_dir" "$local_dest" "$manifest_file"; then
    echo "ERROR: Failed to upload directory $remote_dir" >&2
    return 1
  fi
  
  return 0
}

if (( $# < 2 )); then
  echo "Usage: $0 REMOTE_DIR NAS_DEST [MANIFEST_FILE]" >&2
  echo "Environment variables:" >&2
  echo "  DRY_RUN=true    - Show what would be done without executing" >&2
  echo "  MANIFEST=file   - Override manifest filename" >&2
  exit 2
fi

REMOTE_DIR="$1"        # directory on the HPC to scan (top-level entries)
NAS_DEST="$2"          # local/mounted NAS destination
[[ $# -ge 3 ]] && MANIFEST="$3"

# Check if manifest exists
[[ -f "$MANIFEST" ]] || { echo "Error: Manifest file '$MANIFEST' not found" >&2; exit 1; }

if [[ "$DRY_RUN" == "true" ]]; then
  echo "DRY RUN MODE - No actual transfers will be performed"
  echo "Manifest: $MANIFEST"
  echo "Remote dir: $REMOTE_DIR"
  echo "NAS dest: $NAS_DEST"
  echo
fi

mkdir -p -- "$NAS_DEST"
# Create partial directory for rsync to store incomplete transfers
mkdir -p -- "$NAS_DEST/.rsync-partial"

# Check destination space before starting
echo "Checking destination space..."
if ! df -h "$NAS_DEST" >/dev/null 2>&1; then
  echo "Warning: Could not check destination space" >&2
else
  df -h "$NAS_DEST"
fi

# Setup persistent SSH connection
echo "[1/3] Establishing persistent SSH connection..."
if ! setup_ssh_control; then
  echo "ERROR: Failed to establish SSH connection" >&2
  exit 1
fi

# Cleanup SSH connection on exit
trap 'cleanup_ssh_control' EXIT

echo "[2/3] Uploading items from manifest..."
# Read the manifest via a dedicated FD so we can rewrite it during the loop
exec 3<"$MANIFEST"
# skip header
IFS=$'\t' read -r _i _n _t _s <&3

processed=0
ok=0
failed=0
skipped=0

while IFS=$'\t' read -r idx name type status <&3; do
  [[ -z "${idx:-}" ]] && continue
  
  # Normalize manifest name: strip leading ./ and trailing /
  norm_name="${name#./}"     # drop leading ./
  norm_name="${norm_name%/}" # drop trailing /
  
  # Skip items that are already uploaded
  if [[ "$status" == "UPLOADED" ]]; then
    echo "  [SKIP] $norm_name  -> already uploaded"
    ((skipped++)) || true
    continue
  fi
  
  ((processed++)) || true

  if [[ "$type" == "d" ]]; then
    tarname="${norm_name}.tar.gz"
    echo "  [DIR]  $norm_name  -> tar on HPC, then rsync"
    
    if [[ "$DRY_RUN" == "true" ]]; then
      echo "    DRY RUN: Would check for existing tarball $tarname and rsync to $NAS_DEST/"
      ((ok++)) || true
      continue
    fi
    
    # Check if tarball already exists on remote and has non-zero size
    # Use stat to get file size (works on both BSD and GNU stat)
    # Note: No 'local' keyword here - we're in main script, not a function
    file_size=$(quick_ssh "stat -f%z '$REMOTE_DIR/$tarname' 2>/dev/null || stat -c%s '$REMOTE_DIR/$tarname' 2>/dev/null || echo 0" 2>/dev/null)
    if [[ -n "$file_size" && "$file_size" != "0" && "$file_size" =~ ^[0-9]+$ ]]; then
      echo "    Using existing tarball: $tarname (size: ${file_size} bytes)"
    else
      # File doesn't exist or is empty - clean up and create new
      quick_ssh "rm -f '$REMOTE_DIR/$tarname' '$REMOTE_DIR/${tarname}.partial.'*" 2>/dev/null || true
      echo "    Creating tarball: $tarname (this may take 1-4 hours for large directories)"
      # Use a much longer timeout for tar operations - 4 hours should be enough for even very large directories
      # Tar to temp file first, then move when complete (avoids partial files)
      # The persistent SSH connection will be reused for this operation
      # Use $$ in remote command to get remote PID for unique temp filename
      if ! ssh $(get_ssh_opts) -o ServerAliveInterval=30 -o ServerAliveCountMax=480 "$HPC" \
        "cd '$REMOTE_DIR' && tmp='${tarname}.partial.\$\$' && tar -${TAR_COMPRESSION:-czf} \"\$tmp\" '$norm_name' && mv \"\$tmp\" '$tarname'"; then
        echo "    ✗ Failed: tar creation error"
        # Clean up any partial tar files on failure
        quick_ssh "rm -f '$REMOTE_DIR/${tarname}.partial.'*" 2>/dev/null || true
        update_status "$idx" "FAILED"
        ((failed++)) || true
        continue
      fi
      echo "    Tar creation completed"
    fi
    
    # Now rsync the tarball (whether existing or newly created)
    if robust_rsync -- "$HPC:$REMOTE_DIR/$tarname" "$NAS_DEST/"; then
      # Cleanup remote tarball
      quick_ssh "rm -f '$REMOTE_DIR/$tarname'" || true
      update_status "$idx" "UPLOADED"
      ((ok++)) || true
      echo "    ✓ Success"
    else
      echo "    ✗ Failed: rsync error after retries"
      update_status "$idx" "FAILED"
      ((failed++)) || true
    fi
  else
    echo "  [FILE] $norm_name  -> rsync"
    
    if [[ "$DRY_RUN" == "true" ]]; then
      echo "    DRY RUN: Would rsync file to $NAS_DEST/"
      ((ok++)) || true
      continue
    fi
    
    # Use normalized name for rsync
    if robust_rsync -- "$HPC:$REMOTE_DIR/$norm_name" "$NAS_DEST/"; then
      update_status "$idx" "UPLOADED"
      ((ok++)) || true
      echo "    ✓ Success"
    else
      echo "    ✗ Failed: rsync error after retries"
      update_status "$idx" "FAILED"
      ((failed++)) || true
    fi
  fi
done
exec 3>&-

echo ""
echo "--- Summary ---"
echo "Processed: $processed"
echo "OK:        $ok"
echo "Failed:    $failed"
echo "Skipped:   $skipped"
echo ""
echo "[3/3] Done. See $MANIFEST for statuses."

# Cleanup is handled by trap on EXIT

# Exit with success (0) even if some uploads failed - this allows
# the parent script to continue processing other directories.
# Individual failures are tracked in the manifest status field.
exit 0
