#!/usr/bin/env bash
# stripped-down: make a manifest on HPC, then upload items and update status in the manifest

set -euo pipefail

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
  
  awk -v idx="$idx" -v status="$status" -F'\t' '
    $1 == idx { $4 = status }
    { print }
  ' "$MANIFEST" > "$temp_manifest" && mv "$temp_manifest" "$MANIFEST"
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
  
  # Skip items that are already uploaded
  if [[ "$status" == "UPLOADED" ]]; then
    echo "  [SKIP] $name  -> already uploaded"
    ((skipped++)) || true
    continue
  fi
  
  ((processed++)) || true

  if [[ "$type" == "d" ]]; then
    tarname="${name%/}.tar.gz"
    echo "  [DIR]  $name  -> tar on HPC, then rsync"
    
    if [[ "$DRY_RUN" == "true" ]]; then
      echo "    DRY RUN: Would check for existing tarball $tarname and rsync to $NAS_DEST/"
      ((ok++)) || true
      continue
    fi
    
    # Check if tarball already exists on remote
    if robust_ssh "test -f '$REMOTE_DIR/$tarname'" 2>/dev/null; then
      echo "    Using existing tarball: $tarname"
    else
      echo "    Creating tarball: $tarname"
      if ! robust_ssh "cd '$REMOTE_DIR' && tar -czf '$tarname' '$name'" 2>/dev/null; then
        echo "    ✗ Failed: tar creation error"
        update_status "$idx" "FAILED"
        ((failed++)) || true
        continue
      fi
    fi
    
    # Now rsync the tarball (whether existing or newly created)
    if rsync -a -e "$(get_rsync_ssh)" --timeout="${RSYNC_TIMEOUT:-600}" -- "$HPC:$REMOTE_DIR/$tarname" "$NAS_DEST/" 2>/dev/null; then
      # Cleanup remote tarball
      robust_ssh "rm -f '$REMOTE_DIR/$tarname'" >/dev/null 2>&1 || true
      update_status "$idx" "UPLOADED"
      ((ok++)) || true
      echo "    ✓ Success"
    else
      echo "    ✗ Failed: rsync error"
      update_status "$idx" "FAILED"
      ((failed++)) || true
    fi
  else
    echo "  [FILE] $name  -> rsync"
    
    if [[ "$DRY_RUN" == "true" ]]; then
      echo "    DRY RUN: Would rsync file to $NAS_DEST/"
      ((ok++)) || true
      continue
    fi
    
    if rsync -a -e "$(get_rsync_ssh)" --timeout="${RSYNC_TIMEOUT:-600}" -- "$HPC:$REMOTE_DIR/$name" "$NAS_DEST/" 2>/dev/null; then
      update_status "$idx" "UPLOADED"
      ((ok++)) || true
      echo "    ✓ Success"
    else
      echo "    ✗ Failed: rsync error"
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
