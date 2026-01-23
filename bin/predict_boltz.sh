#!/bin/bash
# Run Boltz predictions in batch mode
set -e

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <cache_dir> <filetype> <final_output_dir> <yaml_file1> [yaml_file2 ...]"
    exit 1
fi

CACHE_DIR=$1
FILETYPE=$2
FINAL_OUTPUT_DIR=$3
shift 3
YAML_FILES="$@"

# Strict GPU Check
nvidia-smi > /dev/null 2>&1 || { echo "No GPU found! Exiting to prevent CPU fallback."; exit 1; }

# Log GPU info to a file in the results directory for verification
mkdir -p "${FINAL_OUTPUT_DIR}/logs"
nvidia-smi --query-gpu=timestamp,name,utilization.gpu,utilization.memory --format=csv > "${FINAL_OUTPUT_DIR}/logs/gpu_usage_${SLURM_JOB_ID:-unknown}.log"

# Fix for "OSError: AF_UNIX path too long"
# Python multiprocessing uses TMPDIR for sockets. GPFS paths are too deep.
# /tmp inside the container is short and sufficient for sockets.
# SAFETY: Use a subdirectory to avoid wiping /tmp if bound from host.
export TMPDIR=$(mktemp -d -p /tmp 2>/dev/null || mktemp -d)

# Performance Optimization for L40S (Tensor Cores)
export TORCH_FLOAT32_MATMUL_PRECISION=medium # if lots of fails, switch to high

if [ "$FILETYPE" != "yaml" ]; then
    echo "Error: filetype must be 'yaml' for Boltz-2"
    exit 1
fi

# Set up local cache directories for the container environment
mkdir -p home_dir numba_cache mpl_config xdg_cache torch_home hf_home results tmp_dir

# Get absolute path of the staged cache directory
CACHE_PATH=$(readlink -f "$CACHE_DIR")

# Link the FULL cache (containing mols and weights)
ln -s "$CACHE_PATH" home_dir/.boltz

export HOME="$PWD/home_dir"
export NUMBA_CACHE_DIR="$PWD/numba_cache"
export MPLCONFIGDIR="$PWD/mpl_config"
export XDG_CACHE_HOME="$PWD/xdg_cache"
export TORCH_HOME="$PWD/torch_home"
export HF_HOME="$PWD/hf_home"

export OMPI_MCA_tmpdir_base="$PWD/tmp_dir"

for yaml_file in $YAML_FILES; do
    [ -e "$yaml_file" ] || continue
    base_name=$(basename "$yaml_file" .yaml)
    
    # Deduplication check
    if [ -d "${FINAL_OUTPUT_DIR}/${base_name}" ]; then
        echo "SKIPPING: ${base_name} (Already exists in ${FINAL_OUTPUT_DIR})"
        continue
    fi
    
    echo "Running Boltz-2 on $base_name"
    
    # Boltz-2 prediction command
    boltz predict "$yaml_file" \
        --cache "home_dir/.boltz" \
        --use_msa_server \
        --num_workers 23 \
        --out_dir "results/${base_name}"

    # Immediate Save: Copy results to final destination to save progress
    if [ -d "results/${base_name}" ]; then
        echo "Copying results for ${base_name} to ${FINAL_OUTPUT_DIR}..."
        cp -r "results/${base_name}" "${FINAL_OUTPUT_DIR}/"
    fi
        
    # Cleanup temporary directory to prevent 'successive task' accumulation
    rm -rf "$TMPDIR"/*
done
