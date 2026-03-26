#!/bin/bash
#SBATCH --job-name=af_ppi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --output=logs/04_%x_%j.out
#SBATCH --account=torch_pr_138_general

TARGET=$1
TARGET_DIR="$SCRATCH/alphafast/$TARGET"
PPI_SCRIPT="$HOME/alphafast_ppi.py"
VENV_PYTHON="$HOME/.af3/bin/python"

if [ -z "$TARGET" ]; then
    echo "Usage: sbatch alphafast_ppi.slurm <TARGET_DIR>"
    exit 1
fi

if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: $TARGET_DIR does not exist."
    exit 1
fi

if [ ! -f "$VENV_PYTHON" ]; then
    echo "Error: Python venv not found at $VENV_PYTHON"
    exit 1
fi

if [ ! -f "$PPI_SCRIPT" ]; then
    echo "Error: ppi.py not found at $PPI_SCRIPT"
    exit 1
fi

mkdir -p logs

echo "Starting PPI analysis for $TARGET at: $(date)"
echo "Input dir: $TARGET_DIR"
echo "Workers:   $SLURM_CPUS_PER_TASK"
echo "Output:    written alongside each *_model.cif"

# Layout: $TARGET_DIR/outputs_task_<N>/<structure_name>/{completion.txt,*_model.cif,...}
# Find every structure folder that has a completion.txt (quality-passed, CIF intact)
# and locate its *_model.cif sibling. Skip already-consolidated output_hits/output_nada.
# -maxdepth 1 ensures we stay in each structure folder.
mapfile -t CIF_FILES < <(
    find "$TARGET_DIR" \
        \( -path "$TARGET_DIR/output_hits" -o -path "$TARGET_DIR/output_nada" \) -prune -o \
        -name "completion.txt" -type f -print0 \
        | xargs -0 -I {} dirname {} \
        | xargs -I {} find {} -maxdepth 1 -name "*_model.cif" -type f 2>/dev/null \
        | sort
)

COUNT=${#CIF_FILES[@]}
echo "Found $COUNT completed, quality-passing CIF files to analyze."

if (( COUNT == 0 )); then
    echo "No CIF files found. Exiting."
    exit 0
fi

run_ppi() {
    local cif="$1"
    local script="$2"
    local python="$3"
    "$python" "$script" "$cif" --outdir "$(dirname "$cif")"
}
export -f run_ppi
export PPI_SCRIPT VENV_PYTHON

printf '%s\n' "${CIF_FILES[@]}" \
    | xargs -P "$SLURM_CPUS_PER_TASK" -I {} \
        bash -c 'run_ppi "$@"' _ {} "$PPI_SCRIPT" "$VENV_PYTHON"

echo "PPI analysis complete at: $(date)"
RESULT_COUNT=$(find "$TARGET_DIR" -name "*_interfaces.csv" | wc -l)
echo "$RESULT_COUNT interface CSV files written alongside their CIF files."


# --- Consolidation & Cleanup ---
echo "Consolidating results into hits and nada..."

HITS_DIR="$TARGET_DIR/output_hits"
NADA_DIR="$TARGET_DIR/output_nada"

mkdir -p "$HITS_DIR" "$NADA_DIR"

# Remove outputs_task_* wrapper from relative paths during consolidation.
strip_task_wrapper() {
    local rel_path="$1"
    echo "${rel_path#outputs_task_*/}"
}

# 1. Move folders with hits.
# Preserve relative path under TARGET_DIR to avoid basename collisions.
find "$TARGET_DIR" -name "*_interfaces.csv" -type f -print0 | while IFS= read -r -d '' csv_path; do
    PARENT_DIR=$(dirname "$csv_path")

    # Only process if not already in a consolidation folder
    if [[ "$PARENT_DIR" != "$HITS_DIR"* && "$PARENT_DIR" != "$NADA_DIR"* ]]; then
        REL_PARENT="$(strip_task_wrapper "${PARENT_DIR#$TARGET_DIR/}")"
        DEST_DIR="$HITS_DIR/$REL_PARENT"

        mkdir -p "$DEST_DIR"
        if rsync -a "$PARENT_DIR/" "$DEST_DIR/"; then
            rm -rf "$PARENT_DIR"
        else
            echo "[WARN] rsync failed; source preserved: $PARENT_DIR"
        fi
    fi
done

# 2. Move remaining completed folders.
# *_interfaces.csv is the only hit marker; otherwise route to nada.
find "$TARGET_DIR" -name "completion.txt" -type f -print0 | while IFS= read -r -d '' comp_path; do
    PARENT_DIR=$(dirname "$comp_path")

    if [[ "$PARENT_DIR" != "$HITS_DIR"* && "$PARENT_DIR" != "$NADA_DIR"* ]]; then
        REL_PARENT="$(strip_task_wrapper "${PARENT_DIR#$TARGET_DIR/}")"

        if find "$PARENT_DIR" -maxdepth 1 -type f -name "*_interfaces.csv" | grep -q .; then
            DEST_DIR="$HITS_DIR/$REL_PARENT"
        else
            DEST_DIR="$NADA_DIR/$REL_PARENT"
        fi

        mkdir -p "$DEST_DIR"
        if rsync -a "$PARENT_DIR/" "$DEST_DIR/"; then
            rm -rf "$PARENT_DIR"
        else
            echo "[WARN] rsync failed; source preserved: $PARENT_DIR"
        fi
    fi
done

echo "Cleaning up parent task directories..."

# 3. Safe Wrapper Cleanup: Only delete outputs_task_* if they are empty of subdirectories
find "$TARGET_DIR" -maxdepth 1 -type d -name "outputs_task_*" | while read -r task_dir; do
    # Count subdirectories inside (excluding the task_dir itself)
    sub_count=$(find "$task_dir" -mindepth 1 -maxdepth 1 -type d | wc -l)
    
    if [ "$sub_count" -eq 0 ]; then
        echo "Removing empty task wrapper: $task_dir"
        rm -rf "$task_dir"
    else
        echo "CAUTION: $task_dir still contains $sub_count folders. Skipping deletion."
    fi
done

# 4. Final Cleanup: Truly empty folders and input task dirs
find "$TARGET_DIR" -type d -empty -delete
find "$TARGET_DIR" -maxdepth 1 -type d -name "inputs_task_*" -exec rm -rf {} +

echo "Consolidation and cleanup complete at: $(date)"

