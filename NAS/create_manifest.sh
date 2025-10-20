#!/usr/bin/env bash
# Usage:
#   ./create_manifest.sh /remote/path/on/hpc [manifest_file]
#
# Output:
#   manifest_<dirname>.tsv  (columns: idx<TAB>name<TAB>type<TAB>status)
#
# Notes:
#   - Top-level only (no recursion). Requires GNU find on the HPC for -printf.
#   - Can override manifest filename via second argument or MANIFEST_FILE env var

set -euo pipefail

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
  source "${SCRIPT_DIR}/config.sh"
fi

# Override with command line argument if provided
[[ $# -ge 2 ]] && MANIFEST_FILE="$2"

if (( $# < 1 )); then
  echo "Usage: $0 REMOTE_DIR [MANIFEST_FILE]" >&2
  echo "Environment variables:" >&2
  echo "  MANIFEST_FILE=file  - Override manifest filename" >&2
  echo "  HPC_USER=user       - Override HPC username" >&2
  echo "  HPC_HOST=host       - Override HPC hostname" >&2
  echo "  REMOTE_BASE=path    - Override remote base path" >&2
  exit 2
fi

REMOTE_DIR="$1"

# Use provided manifest file or derive from directory name
if [[ -n "${MANIFEST_FILE:-}" ]]; then
  MANIFEST="$MANIFEST_FILE"
else
  # Derive a safe manifest filename from the dir being scanned
  DIR_CLEAN="${REMOTE_DIR%/}"
  DIR_LABEL="$(basename "$DIR_CLEAN")"
  if [[ -z "$DIR_LABEL" || "$DIR_LABEL" == "/" || "$DIR_LABEL" == "." ]]; then
    DIR_LABEL="${DIR_CLEAN//\//_}"
    [[ -z "$DIR_LABEL" ]] && DIR_LABEL="root"
  fi
  DIR_SAFE="$(printf '%s' "$DIR_LABEL" | tr -cs 'A-Za-z0-9._-' '_')"
  MANIFEST="manifest_${DIR_SAFE}.tsv"
fi

TMP_LIST="$(mktemp)"; trap 'rm -f "$TMP_LIST"' EXIT

echo "[1/2] Scanning directory: ${HPC}:${REMOTE_DIR}"

# Test SSH connection first
if ! test_ssh_connection; then
  echo "Error: Cannot establish SSH connection to HPC" >&2
  echo "Please check your SSH configuration and network connectivity" >&2
  exit 1
fi

# Build list of immediate children on the HPC: <type>\t<name>
echo "Scanning remote directory contents..."
if ! robust_ssh "cd '$REMOTE_DIR' && find . -mindepth 1 -maxdepth 1 -printf '%y\t%P\n'" > "$TMP_LIST" 2>/dev/null; then
  echo "Error: Failed to access directory '$REMOTE_DIR' on HPC" >&2
  echo "Please check that the directory exists and you have read permissions" >&2
  exit 1
fi

echo "[2/2] Writing manifest to $MANIFEST"

# Preserve previous statuses (join by name); new items => PENDING
if [[ -f "$MANIFEST" ]]; then
  awk -F'\t' '
    NR==FNR { if (FNR>1) st[$2]=$4; next }                 # old manifest -> map name->status
    BEGIN { print "idx\tname\ttype\tstatus"; i=0 }
    { i++; name=$2; type=$1; status = (name in st ? st[name] : "PENDING");
      printf "%d\t%s\t%s\t%s\n", i, name, type, status
    }
  ' "$MANIFEST" "$TMP_LIST" > "${MANIFEST}.tmp"
else
  awk -F'\t' '
    BEGIN { print "idx\tname\ttype\tstatus"; i=0 }
    { i++; printf "%d\t%s\t%s\tPENDING\n", i, $2, $1 }
  ' "$TMP_LIST" > "${MANIFEST}.tmp"
fi

mv "${MANIFEST}.tmp" "$MANIFEST"
echo "Done. $(wc -l < "$MANIFEST") lines (incl. header)."

# Function to create manifests for multiple directories (for use by other scripts)
create_manifests_for_directories() {
  local dir_list_file="$1"
  local remote_base="${2:-$REMOTE_BASE}"
  local local_base="${3:-$(pwd)}"
  
  if [[ ! -f "$dir_list_file" ]]; then
    echo "ERROR: Directory list file '$dir_list_file' not found" >&2
    return 1
  fi
  
  local processed=0
  local ok=0
  local failed=0
  
  while IFS= read -r dir_name || [[ -n "${dir_name:-}" ]]; do
    # Skip empty lines and comments
    dir_name="${dir_name%%#*}"
    dir_name="${dir_name%"${dir_name##*[![:space:]]}"}"
    [[ -z "$dir_name" ]] && continue
    
    ((processed++)) || true
    
    local remote_path="$remote_base/$dir_name"
    local manifest_file="manifest_${dir_name//[^A-Za-z0-9._-]/_}.tsv"
    
    echo "Creating manifest for: $dir_name"
    if "$0" "$remote_path" "$manifest_file"; then
      ((ok++)) || true
    else
      ((failed++)) || true
    fi
    
  done < "$dir_list_file"
  
  echo "Manifest creation summary: $processed processed, $ok successful, $failed failed"
  return $((failed > 0 ? 1 : 0))
}
