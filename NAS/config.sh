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

# Directory list files
export LARGE_DIR_LIST="${LARGE_DIR_LIST:-large_directory_list.txt}"
export LIST_FILE="${LIST_FILE:-directory_list.txt}"

# Transfer settings
export DRY_RUN="${DRY_RUN:-false}"
export LOG_TO_CONSOLE="${LOG_TO_CONSOLE:-false}"
export CONTINUE_ON_ERROR="${CONTINUE_ON_ERROR:-true}"  # Continue processing even if individual items fail
export DEBUG="${DEBUG:-false}"  # Enable debug mode for more verbose SSH output

# SSH ControlMaster socket path (for persistent connections)
export SSH_CONTROL_DIR="${SSH_CONTROL_DIR:-${HOME}/.ssh/nas_control}"
export SSH_CONTROL_SOCKET="${SSH_CONTROL_SOCKET:-${SSH_CONTROL_DIR}/${HPC_USER}@${HPC_HOST}:22}"

# SSH Options (array format for easy use)
# ServerAliveInterval=30 sends keepalives every 30s, ServerAliveCountMax=200 allows 100 minutes of silence
# Use ERROR loglevel by default (quiet), can enable DEBUG=true for verbose SSH output
# ControlMaster=auto enables connection multiplexing for persistent SSH connections
if [[ "${DEBUG:-false}" == "true" ]]; then
  export SSH_OPTS_STR="-n -T -o BatchMode=yes -o ControlMaster=auto -o ControlPath=${SSH_CONTROL_SOCKET} -o ControlPersist=10m -o ConnectTimeout=60 -o ServerAliveInterval=30 -o ServerAliveCountMax=200 -o StrictHostKeyChecking=accept-new -o LogLevel=INFO -o TCPKeepAlive=yes -o Compression=no -o ForwardAgent=no -o IdentitiesOnly=yes"
else
  export SSH_OPTS_STR="-n -T -o BatchMode=yes -o ControlMaster=auto -o ControlPath=${SSH_CONTROL_SOCKET} -o ControlPersist=10m -o ConnectTimeout=60 -o ServerAliveInterval=30 -o ServerAliveCountMax=200 -o StrictHostKeyChecking=accept-new -o LogLevel=ERROR -o TCPKeepAlive=yes -o Compression=no -o ForwardAgent=no -o IdentitiesOnly=yes"
fi

# Connection retry settings
export SSH_MAX_RETRIES="${SSH_MAX_RETRIES:-5}"
export SSH_RETRY_DELAY="${SSH_RETRY_DELAY:-10}"
export SSH_CONNECTION_TIMEOUT="${SSH_CONNECTION_TIMEOUT:-60}"
export SSH_QUICK_RETRIES="${SSH_QUICK_RETRIES:-2}"  # Fewer retries for quick checks
export SSH_QUICK_DELAY="${SSH_QUICK_DELAY:-2}"  # Shorter delay for quick checks

# Rsync settings
export RSYNC_PATH="${RSYNC_PATH:-/usr/bin/rsync}"
export RSYNC_TIMEOUT="${RSYNC_TIMEOUT:-600}"
export RSYNC_MAX_RETRIES="${RSYNC_MAX_RETRIES:-3}"
export RSYNC_RETRY_DELAY="${RSYNC_RETRY_DELAY:-5}"

# Tar settings
export TAR_COMPRESSION="${TAR_COMPRESSION:-czf}"  # czf = gzip, cf = no compression
export TAR_SSH_TIMEOUT="${TAR_SSH_TIMEOUT:-14400}"  # 4 hours for large directories

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
      echo "SSH retry ($attempt/$max_retries) - waiting ${retry_delay}s..." >&2
      sleep "$retry_delay"
      # Increase delay for subsequent attempts
      retry_delay=$((retry_delay * 2))
    fi
    ((attempt++))
  done
  
  # Only print final failure message if we've retried
  if (( max_retries > 1 )); then
    echo "SSH command failed after $max_retries attempts" >&2
  fi
  return $last_exit_code
}

# Quick SSH for fast checks (fewer retries, shorter delays)
quick_ssh() {
  local max_retries="${SSH_QUICK_RETRIES:-2}"
  local retry_delay="${SSH_QUICK_DELAY:-2}"
  local attempt=1
  local last_exit_code=0
  
  while (( attempt <= max_retries )); do
    if ssh $(get_ssh_opts) "$HPC" "$@"; then
      return 0
    fi
    last_exit_code=$?
    
    if (( attempt < max_retries )); then
      sleep "$retry_delay"
    fi
    ((attempt++))
  done
  
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

# Setup SSH master connection
setup_ssh_control() {
  # Create control directory if it doesn't exist
  mkdir -p "$SSH_CONTROL_DIR"
  
  # Test the connection to establish master
  if ssh $(get_ssh_opts) "$HPC" "echo 'SSH ControlMaster established'" >/dev/null 2>&1; then
    if [[ "${DEBUG:-false}" == "true" ]]; then
      echo "SSH ControlMaster connection established"
    fi
    return 0
  else
    echo "Failed to establish SSH ControlMaster connection" >&2
    return 1
  fi
}

# Cleanup SSH master connection
cleanup_ssh_control() {
  # Close the master connection
  ssh -O exit -o ControlPath="$SSH_CONTROL_SOCKET" "$HPC" >/dev/null 2>&1 || true
  
  # Clean up any stale socket files
  if [[ -S "$SSH_CONTROL_SOCKET" ]]; then
    rm -f "$SSH_CONTROL_SOCKET" || true
  fi
  
  if [[ "${DEBUG:-false}" == "true" ]]; then
    echo "SSH ControlMaster connection closed"
  fi
}

# Robust rsync function with retries
robust_rsync() {
  local max_retries="${RSYNC_MAX_RETRIES:-3}"
  local retry_delay="${RSYNC_RETRY_DELAY:-5}"
  local attempt=1
  local last_exit_code=0
  
  while (( attempt <= max_retries )); do
    # Use -avz for verbose output, compression, progress; --partial keeps partial transfers
    # --partial-dir stores partials in hidden directory to allow resume
    if rsync -avz --progress --partial --partial-dir=.rsync-partial --timeout="${RSYNC_TIMEOUT:-600}" \
         -e "$(get_rsync_ssh)" "$@"; then
      return 0
    fi
    last_exit_code=$?
    
    if (( attempt < max_retries )); then
      echo "Rsync attempt $attempt failed (exit code $last_exit_code). Retrying in ${retry_delay}s..." >&2
      sleep "$retry_delay"
    fi
    ((attempt++))
  done
  
  echo "Rsync failed after $max_retries attempts (final exit code $last_exit_code)" >&2
  return $last_exit_code
}
