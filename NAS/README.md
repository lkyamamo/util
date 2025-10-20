# NAS Transfer Scripts

A collection of bash scripts for transferring data between USC HPC cluster and local NAS storage.

## Scripts Overview

### 1. `create_manifest.sh`
Creates a manifest file listing all items in a remote HPC directory.

**Usage:**
```bash
./create_manifest.sh /remote/path/on/hpc [manifest_file]
```

**Features:**
- Scans top-level directory contents on HPC
- Creates TSV manifest with columns: idx, name, type, status
- Preserves existing status information when updating
- Safe filename generation for manifest files

### 2. `rsync_list.sh`
Syncs directories and files from HPC to local storage based on a list file.

**Usage:**
```bash
./rsync_list.sh
```

**Features:**
- Reads from `directory_list.txt` by default
- Handles both directories and files differently
- Retry logic with exponential backoff
- Comprehensive logging
- Partial transfer recovery

### 3. `upload_directory.sh`
Uploads items from a manifest, creating tarballs for directories.

**Usage:**
```bash
./upload_directory.sh REMOTE_DIR NAS_DEST [MANIFEST_FILE]
```

**Features:**
- Creates tarballs for directories on HPC
- Syncs files and tarballs to local NAS
- Updates manifest with transfer status
- Dry-run mode for testing
- Progress tracking and summary

### 4. `process_large_directories.sh`
Processes multiple directories from a list file, creating tarballs and syncing each to local directories.

**Usage:**
```bash
./process_large_directories.sh
```

**Features:**
- Reads directory list from `large_directory_list.txt`
- Creates individual manifests for each directory
- Processes each directory independently
- Creates local directories for each remote directory
- Comprehensive logging and progress tracking
- Dry-run mode support

## Configuration

All scripts use `config.sh` for centralized configuration:

```bash
# HPC Connection Settings
export HPC_USER="lkyamamo"
export HPC_HOST="endeavour.usc.edu"

# Remote paths
export REMOTE_BASE="/project/priyav_216/Tians"

# Local settings
export MANIFEST_FILE="manifest.tsv"
export LOG_FILE="update_log.txt"

# Transfer settings
export DRY_RUN="false"
export LOG_TO_CONSOLE="false"
export CONTINUE_ON_ERROR="true"  # Continue processing even if individual items fail
```

## Environment Variables

You can override any configuration by setting environment variables:

```bash
# Override HPC credentials
export HPC_USER="your_username"
export HPC_HOST="your_hpc_host"

# Enable dry-run mode
export DRY_RUN="true"

# Enable console logging
export LOG_TO_CONSOLE="true"

# Override manifest file
export MANIFEST_FILE="custom_manifest.tsv"

# Stop on first error (default is to continue)
export CONTINUE_ON_ERROR="false"
```

## Data Files

- `directory_list.txt`: List of directories to sync from HPC
- `file_list.txt`: List of tarball files
- `manifest_*.tsv`: Generated manifest files with transfer status

## Example Workflow

1. **Create a manifest:**
   ```bash
   ./create_manifest.sh /project/priyav_216/Tians/data
   ```

2. **Test with dry-run:**
   ```bash
   DRY_RUN=true ./upload_directory.sh /project/priyav_216/Tians/data /local/nas/dest
   ```

3. **Perform actual transfer:**
   ```bash
   ./upload_directory.sh /project/priyav_216/Tians/data /local/nas/dest
   ```

4. **Sync specific items:**
   ```bash
   ./rsync_list.sh
   ```

5. **Test SSH connection robustness:**
   ```bash
   ./test_ssh_robustness.sh
   ```

## Error Handling

All scripts include:
- **Resilient processing**: By default, scripts continue processing even if individual items fail
- Retry mechanisms for network operations (3 attempts with exponential backoff)
- Comprehensive logging with clear success/failure indicators
- Status tracking and reporting in manifests
- Graceful handling of partial transfers
- Clear error messages with ✓/✗ indicators

**Error Behavior:**
- Individual item failures are logged but don't stop the entire process
- Failed items are marked as "FAILED" in the manifest
- A summary shows total processed, successful, and failed items
- **Robust SSH connections** with automatic retry on failures
- **Connection validation** before major operations
- **Graceful handling** of network timeouts and disconnections

**SSH Robustness Features:**
- Automatic retry with exponential backoff (5 attempts by default)
- Connection keepalive settings to prevent timeouts
- Pre-flight connection testing before operations
- Enhanced timeout handling (60s connection, 600s rsync)
- TCP keepalive and compression optimization

## Requirements

- GNU `find` on the HPC (for `-printf` option)
- SSH key-based authentication to HPC
- `rsync` on both local and remote systems
- `awk` for manifest processing
