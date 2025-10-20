#!/usr/bin/env bash
# Configuration file for NAS transfer scripts
# Source this file in your scripts or set these as environment variables

# HPC Connection Settings
export HPC_USER="${HPC_USER:-lkyamamo}"
export HPC_HOST="${HPC_HOST:-endeavour.usc.edu}"
export HPC="${HPC_USER}@${HPC_HOST}"

# Remote paths
export REMOTE_BASE="${REMOTE_BASE:-/project/priyav_216/Tians}"

# Local settings
export MANIFEST_FILE="${MANIFEST_FILE:-manifest.tsv}"
export LOG_FILE="${LOG_FILE:-update_log.txt}"
export LOCAL_BASE="${LOCAL_BASE:-$(pwd)}"

# Large directory processing
export LARGE_DIR_LIST="${LARGE_DIR_LIST:-large_directory_list.txt}"

# Transfer settings
export DRY_RUN="${DRY_RUN:-false}"
export LOG_TO_CONSOLE="${LOG_TO_CONSOLE:-false}"
export CONTINUE_ON_ERROR="${CONTINUE_ON_ERROR:-true}"  # Continue processing even if individual items fail

# SSH Options (array format for easy use)
export SSH_OPTS_STR="-n -T -o BatchMode=yes -o ControlMaster=no -o ConnectTimeout=60 -o ServerAliveInterval=30 -o ServerAliveCountMax=10 -o StrictHostKeyChecking=accept-new -o LogLevel=ERROR -o TCPKeepAlive=yes -o Compression=no -o ForwardAgent=no -o IdentitiesOnly=yes"

# Connection retry settings
export SSH_MAX_RETRIES="${SSH_MAX_RETRIES:-5}"
export SSH_RETRY_DELAY="${SSH_RETRY_DELAY:-10}"
export SSH_CONNECTION_TIMEOUT="${SSH_CONNECTION_TIMEOUT:-60}"

# Rsync settings
export RSYNC_PATH="${RSYNC_PATH:-/usr/bin/rsync}"
export RSYNC_TIMEOUT="${RSYNC_TIMEOUT:-600}"

# Helper function to get SSH options as array
get_ssh_opts() {
  echo $SSH_OPTS_STR
}

# Helper function to get rsync SSH command
get_rsync_ssh() {
  echo "ssh $(get_ssh_opts)"
}

# Robust SSH connection function with retries
robust_ssh() {
  local max_retries="${SSH_MAX_RETRIES:-5}"
  local retry_delay="${SSH_RETRY_DELAY:-10}"
  local attempt=1
  local last_exit_code=0
  
  while (( attempt <= max_retries )); do
    if ssh $(get_ssh_opts) "$HPC" "$@"; then
      return 0
    fi
    last_exit_code=$?
    
    if (( attempt < max_retries )); then
      echo "SSH attempt $attempt failed (exit code $last_exit_code). Retrying in ${retry_delay}s..." >&2
      sleep "$retry_delay"
      # Increase delay for subsequent attempts
      retry_delay=$((retry_delay * 2))
    fi
    ((attempt++))
  done
  
  echo "SSH failed after $max_retries attempts (final exit code $last_exit_code)" >&2
  return $last_exit_code
}

# Test SSH connection
test_ssh_connection() {
  echo "Testing SSH connection to $HPC..."
  if robust_ssh "echo 'SSH connection test successful'"; then
    echo "✓ SSH connection is working"
    return 0
  else
    echo "✗ SSH connection failed"
    return 1
  fi
}
