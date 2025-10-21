#!/usr/bin/env bash
# Process multiple directories from large_directory_list.txt
# Creates tarballs of each directory on HPC and syncs them to local directories

set -euo pipefail

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/config.sh" ]]; then
  source "${SCRIPT_DIR}/config.sh"
fi

# Configuration
LARGE_DIR_LIST="${LARGE_DIR_LIST:-large_directory_list.txt}"
REMOTE_BASE="${REMOTE_BASE:-/project/priyav_216/Tians}"
LOCAL_BASE="${LOCAL_BASE:-$(pwd)}"

# Function to log messages
log() {
  local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[$timestamp] $*"
}

# Function to check if all items in manifest are uploaded
check_manifest_status() {
  local manifest_file="$1"
  
  if [[ ! -f "$manifest_file" ]]; then
    return 1  # No manifest file exists
  fi
  
  # Check if all non-header lines have status "UPLOADED"
  local total_items=0
  local uploaded_items=0
  
  while IFS=$'\t' read -r idx name type status; do
    [[ "$idx" == "idx" ]] && continue  # Skip header
    [[ -z "${idx:-}" ]] && continue    # Skip empty lines
    
    ((total_items++)) || true
    if [[ "$status" == "UPLOADED" ]]; then
      ((uploaded_items++)) || true
    fi
  done < "$manifest_file"
  
  # Return 0 if all items are uploaded, 1 otherwise
  [[ $total_items -gt 0 && $uploaded_items -eq $total_items ]]
}

# Function to process a single directory
process_directory() {
  local remote_dir="$1"
  local local_dir="$2"
  local manifest_file="$3"
  
  log "Processing directory: $remote_dir"
  
  # Create local directory
  mkdir -p "$local_dir"
  
  # Check if manifest already exists and all items are uploaded
  if [[ -f "$manifest_file" ]]; then
    log "Found existing manifest: $manifest_file"
    if check_manifest_status "$manifest_file"; then
      log "All items in manifest are already UPLOADED - skipping directory"
      return 0
    else
      log "Manifest exists but has incomplete uploads - will retry"
    fi
  else
    log "No existing manifest found - creating new one"
  fi
  
  # Create or update manifest for this directory
  if ! "$SCRIPT_DIR/create_manifest.sh" "$remote_dir" "$manifest_file"; then
    log "ERROR: Failed to create manifest for $remote_dir"
    return 1
  fi
  
  # Upload directory contents
  if ! "$SCRIPT_DIR/upload_directory.sh" "$remote_dir" "$local_dir" "$manifest_file"; then
    log "ERROR: Failed to upload directory $remote_dir"
    return 1
  fi
  
  log "Completed processing: $remote_dir -> $local_dir"
  return 0
}

# Main processing function
main() {
  if (( $# > 0 )); then
    echo "Usage: $0" >&2
    echo "Environment variables:" >&2
    echo "  LARGE_DIR_LIST=file     - Override directory list file (default: large_directory_list.txt)" >&2
    echo "  REMOTE_BASE=path        - Override remote base path" >&2
    echo "  LOCAL_BASE=path         - Override local base path" >&2
    echo "  DRY_RUN=true            - Show what would be done without executing" >&2
    exit 2
  fi
  
  # Check if directory list exists
  if [[ ! -f "$LARGE_DIR_LIST" ]]; then
    log "ERROR: Directory list file '$LARGE_DIR_LIST' not found"
    log "Please create a file with one directory name per line"
    exit 1
  fi
  
  log "Starting large directory processing"
  log "Directory list: $LARGE_DIR_LIST"
  log "Remote base: $REMOTE_BASE"
  log "Local base: $LOCAL_BASE"
  
  # Test SSH connection first
  if ! test_ssh_connection; then
    log "ERROR: Cannot establish SSH connection to HPC"
    log "Please check your SSH configuration and network connectivity"
    exit 1
  fi
  
  if [[ "$DRY_RUN" == "true" ]]; then
    log "DRY RUN MODE - No actual transfers will be performed"
  fi
  
  # Process each directory
  local processed=0
  local ok=0
  local failed=0
  
  while IFS= read -r dir_name || [[ -n "${dir_name:-}" ]]; do
    # Skip empty lines and comments
    dir_name="${dir_name%%#*}"  # Remove comments
    dir_name="${dir_name%"${dir_name##*[![:space:]]}"}"  # Trim trailing whitespace
    [[ -z "$dir_name" ]] && continue
    
    ((processed++)) || true
    
    local remote_path="$REMOTE_BASE/$dir_name"
    local local_path="$LOCAL_BASE/$dir_name"
    local manifest_file="manifest_${dir_name//[^A-Za-z0-9._-]/_}.tsv"
    
    log ""
    log "=== [$processed] Processing: $dir_name ==="
    
    if process_directory "$remote_path" "$local_path" "$manifest_file"; then
      ((ok++)) || true
    else
      ((failed++)) || true
    fi
    
  done < "$LARGE_DIR_LIST"
  
  # Summary
  log ""
  log "=== FINAL SUMMARY ==="
  log "Total directories processed: $processed"
  log "Successfully completed: $ok"
  log "Failed: $failed"
  
  if [[ $failed -gt 0 ]]; then
    log "Some directories failed to process. Check the logs above for details."
    exit 1
  else
    log "All directories processed successfully!"
  fi
}

# Run main function
main "$@"
