# NAS Transfer Scripts - Quick Tutorial

A simple guide to transferring data from HPC to local NAS storage.

## Overview

These scripts help you transfer data from USC HPC (`endeavour.usc.edu`) to your local NAS. They handle:
- Creating manifests (lists of files/directories to transfer)
- Uploading files and directories (creating tarballs for large directories)
- Syncing from remote using rsync
- Processing multiple directories automatically

---

## Quick Start

### 1. Configuration (`config.sh`)

**What it does:** Central configuration file that all scripts use. You usually don't need to edit this unless you want to change defaults.

**Key settings:**
- HPC connection: `lkyamamo@endeavour.usc.edu`
- Remote base path: `/project/priyav_216/Tians`
- SSH retry settings, rsync timeouts, etc.

**To customize:**
```bash
# Edit config.sh or set environment variables:
export HPC_USER="your_username"
export REMOTE_BASE="/your/remote/path"
export DRY_RUN="true"  # Test mode - no actual transfers
```

---

### 2. Simple File Sync (`rsync_list.sh`)

**What it does:** Syncs directories and files from HPC to local based on a list.

**Usage:**
```bash
# 1. Create a list file: directory_list.txt
#    Put one directory or file name per line:
echo "my_directory" >> directory_list.txt
echo "another_dir" >> directory_list.txt
echo "large_file.tar.gz" >> directory_list.txt

# 2. Run the script
./rsync_list.sh
```

**What happens:**
- Reads `directory_list.txt` (one item per line)
- For each item, checks if it's a directory or file on HPC
- Syncs directories with `rsync` (uses partial transfers for resume)
- Syncs files with `rsync` (uses in-place updates)
- Logs to `update_log.txt`

**Example:**
```bash
# Create your list
cat > directory_list.txt << EOF
data_2024
results
output.tar.gz
EOF

# Sync everything
./rsync_list.sh
```

**Note**
This script only reads from the remote HPC and writes to local storage. It does not create any files on the remote, so it's safe to use on read-only directories.

---

### 3. Upload with Manifest (`upload_directory.sh`)

**What it does:** Uploads a single directory from HPC to local NAS, creating tarballs for directories.

**Usage:**
```bash
./upload_directory.sh REMOTE_DIR LOCAL_DEST [MANIFEST_FILE]
```

**Example:**
```bash
# Upload a directory
./upload_directory.sh \
  /path/to/my_data \
  /local/nas/destination \
  manifest_my_data.tsv
```

**What happens:**
1. **Requires a manifest file** (created by `create_manifest.sh`)
2. Reads the manifest to see what needs uploading
3. For **directories**: Creates a tarball on HPC, then rsyncs it to local
4. For **files**: Directly rsyncs from HPC to local
5. Updates manifest status: `PENDING` → `UPLOADED` or `FAILED`
6. Skips items already marked `UPLOADED`

**Workflow:**
```bash
# Step 1: Create manifest
./create_manifest.sh /path/to/my_data manifest_my_data.tsv

# Step 2: Upload (test first with dry-run)
DRY_RUN=true ./upload_directory.sh \
  /path/to/my_data \
  /local/nas/dest \
  manifest_my_data.tsv

# Step 3: Actually upload
./upload_directory.sh \
  /path/to/my_data \
  /local/nas/dest \
  manifest_my_data.tsv
```

**Manifest format:**
```
idx	name	type	status
1	subdir1	d	PENDING
2	file1.txt	f	PENDING
3	subdir2	d	UPLOADED
```

---

**Note**
Since this script creates tarballs on the remote HPC, it cannot be used on read-only directories like project1.

### 4. Process Multiple Directories (`process_large_directories.sh`)

**What it does:** Automatically processes multiple directories from a list - creates manifests and uploads each one.

**Usage:**
```bash
# 1. Create a list of directories to process
cat > large_directory_list.txt << EOF
big_data_2024
simulation_results
experiment_1
EOF

# 2. Run the script
./process_large_directories.sh
```

**What happens:**
- Reads `large_directory_list.txt`
- For each directory:
  1. Creates a manifest (e.g., `manifest_big_data_2024.tsv`)
  2. Uploads all items from that manifest
  3. Skips directories where all items are already `UPLOADED`
- Creates local directories matching remote structure
- Shows progress and summary

**Example:**
```bash
# Create directory list
echo "1020-rxmd-algo" > large_directory_list.txt
echo "calculations" >> large_directory_list.txt

# Process everything
./process_large_directories.sh

# Output:
# manifest_1020-rxmd-algo.tsv
# manifest_calculations.tsv
# Local directories: ./1020-rxmd-algo/, ./calculations/
```

**Note**
Since this script uses `upload_directory.sh` which creates tarballs on the remote HPC, it cannot be used on read-only directories like project1.

---

## Common Workflows

### Workflow 1: Sync a few specific directories/files

```bash
# Create list
cat > directory_list.txt << EOF
dir1
dir2
file.tar.gz
EOF

# Sync
./rsync_list.sh
```

### Workflow 2: Upload one large directory with tracking

```bash
# Create manifest
./create_manifest.sh \
  /path/to/large_dir \
  manifest_large_dir.tsv

# Upload (can resume if interrupted)
./upload_directory.sh \
  /path/to/large_dir \
  /local/nas/large_dir \
  manifest_large_dir.tsv
```

### Workflow 3: Process many directories automatically

```bash
# Create list
cat > large_directory_list.txt << EOF
project1
project2
project3
EOF

# Process all
./process_large_directories.sh
```

---

## Useful Options

### Dry Run (Test Mode)
```bash
DRY_RUN=true ./upload_directory.sh REMOTE_DIR LOCAL_DEST MANIFEST
```

### Enable Console Logging
```bash
LOG_TO_CONSOLE=true ./rsync_list.sh
```

### Debug Mode (Verbose SSH)
```bash
DEBUG=true ./upload_directory.sh REMOTE_DIR LOCAL_DEST MANIFEST
```

### Uncompressed Tarballs (Faster)
```bash
TAR_COMPRESSION=cf ./upload_directory.sh REMOTE_DIR LOCAL_DEST MANIFEST
```

---

## File Structure

```
NAS/
├── config.sh                          # Configuration
├── create_manifest.sh                  # Create manifest from remote dir
├── upload_directory.sh                 # Upload one directory
├── process_large_directories.sh        # Process multiple directories
├── rsync_list.sh                      # Simple sync from list
├── directory_list.txt                  # List for rsync_list.sh
├── large_directory_list.txt            # List for process_large_directories.sh
└── manifest_*.tsv                      # Generated manifests
```

---

## Tips

1. **Always test first:** Use `DRY_RUN=true` to see what would happen
2. **Manifests are resumable:** If a transfer fails, rerun the script - it will skip `UPLOADED` items
3. **Large directories:** Creating tarballs can take 1-4 hours for very large dirs (TB scale)
4. **Partial transfers:** rsync automatically resumes interrupted transfers
5. **SSH connections:** Scripts use persistent SSH connections for efficiency

---

## Troubleshooting

**SSH connection fails:**
- Check your SSH key is set up: `ssh lkyamamo@endeavour.usc.edu`
- Verify network connectivity

**Manifest not found:**
- Run `create_manifest.sh` first to create the manifest

**Large directory takes too long:**
- Use `TAR_COMPRESSION=cf` for uncompressed (faster) tarballs
- Check remote disk space and I/O

**Transfer fails:**
- Check disk space: `df -h LOCAL_DEST`
- Check permissions
- Rerun the script - it will resume from where it left off

---

## Summary

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `rsync_list.sh` | Simple sync | `directory_list.txt` | Synced files/dirs |
| `create_manifest.sh` | Create tracking file | Remote dir | `manifest_*.tsv` |
| `upload_directory.sh` | Upload one dir | Remote dir + manifest | Uploaded files + updated manifest |
| `process_large_directories.sh` | Process many dirs | `large_directory_list.txt` | Multiple manifests + uploaded dirs |

**Quick decision tree:**
- Need to sync a few items? → `rsync_list.sh`
- Need to track upload progress? → `create_manifest.sh` + `upload_directory.sh`
- Need to process many directories? → `process_large_directories.sh`

