#!/bin/bash
# Run Boltz predictions in batch mode
set -e

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <cache_dir> <filetype> <yaml_file1> [yaml_file2 ...]"
    exit 1
fi

CACHE_DIR=$1
FILETYPE=$2
shift 2
YAML_FILES="$@"

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

# Fix for Read-only /scratch error
export TMPDIR="$PWD/tmp_dir"
export TEMP="$PWD/tmp_dir"
export TMP="$PWD/tmp_dir"
export OMPI_MCA_tmpdir_base="$PWD/tmp_dir"

for yaml_file in $YAML_FILES; do
    [ -e "$yaml_file" ] || continue
    base_name=$(basename "$yaml_file" .yaml)
    
    echo "Running Boltz-2 on $base_name"
    
    # Boltz-2 prediction command
    boltz predict "$yaml_file" \
        --cache "home_dir/.boltz" \
        --use_msa_server \
        --num_workers 0 \
        --out_dir "results/${base_name}"
done
