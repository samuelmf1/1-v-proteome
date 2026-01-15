#!/bin/bash

# Load environment
module load singularity

# Define paths
TEMP_BUILD_DIR="/scratch/$USER/singularity_tmp"
CACHE_DIR="/scratch/$USER/singularity_cache"
FINAL_DEST="./singularity_cache"

mkdir -p "$TEMP_BUILD_DIR" "$CACHE_DIR" "$FINAL_DEST"

# Export variables
export SINGULARITY_TMPDIR="$TEMP_BUILD_DIR"
export SINGULARITY_CACHEDIR="$CACHE_DIR"

echo "Starting pull of Boltz2..."

# 1. Removed -v from pull (or you can use 'singularity -v pull')
# 2. Added -F to overwrite any failed/corrupted partial downloads
if singularity -v pull -F "$TEMP_BUILD_DIR/boltz2.sif" docker://nvcr.io/nim/mit/boltz2:1.0.0; then
    mv "$TEMP_BUILD_DIR/boltz2.sif" "$FINAL_DEST/nvcr.io-nim-mit-boltz2-1.0.0.img"
    echo "Pull and move complete."
else
    echo "Error: Singularity pull failed."
    exit 1
fi

echo "Cleaning up build artifacts..."
rm -rf "$TEMP_BUILD_DIR"