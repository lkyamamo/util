# NAS Transfer Scripts - Clean File Structure

## Essential Files

### Core Scripts
- **`config.sh`** - Central configuration file with SSH settings and robust connection functions
- **`create_manifest.sh`** - Creates manifest files listing remote directory contents
- **`upload_directory.sh`** - Uploads items from manifest, creating tarballs for directories
- **`rsync_list.sh`** - Syncs directories and files from HPC based on list files
- **`process_large_directories.sh`** - Main script for processing multiple directories from a list

### Configuration Files
- **`large_directory_list.txt`** - List of directories to process in batch
- **`directory_list.txt`** - List of specific directories/files to sync
- **`file_list.txt`** - List of tarball files to sync

### Documentation
- **`README.md`** - Complete usage documentation and examples

## File Purposes

### Scripts
1. **config.sh** - Contains all configuration variables, SSH settings, and robust connection functions
2. **create_manifest.sh** - Scans remote directories and creates TSV manifests with file/directory metadata
3. **upload_directory.sh** - Processes manifests to create tarballs and sync them locally
4. **rsync_list.sh** - Syncs specific items from list files using optimized rsync settings
5. **process_large_directories.sh** - Orchestrates batch processing of multiple directories

### Data Files
1. **large_directory_list.txt** - Contains directory names (one per line) for batch processing
2. **directory_list.txt** - Contains specific paths for targeted syncing
3. **file_list.txt** - Contains tarball filenames for file-based syncing

## Usage Workflow

### Single Directory Processing
```bash
# 1. Create manifest
./create_manifest.sh /remote/path

# 2. Upload directory
./upload_directory.sh /remote/path /local/dest
```

### Batch Directory Processing
```bash
# 1. Edit large_directory_list.txt with directory names
# 2. Process all directories
./process_large_directories.sh
```

### Targeted Syncing
```bash
# 1. Edit directory_list.txt with specific paths
# 2. Sync selected items
./rsync_list.sh
```

## Generated Files (Runtime)
- **`manifest_*.tsv`** - Generated manifest files with transfer status
- **`update_log.txt`** - Log file for rsync operations
- **`*.tar.gz`** - Tarballs created during directory processing

## Cleanup Status
✅ Removed all test files and temporary artifacts
✅ Kept only essential production files
✅ All scripts validated and working
✅ Clean, organized structure ready for production use
