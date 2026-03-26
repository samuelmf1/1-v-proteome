#!/bin/bash
#SBATCH --job-name=af_cleanup
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16GB
#SBATCH --output=logs/03_%x_%j.out
#SBATCH --time=04:00:00 
#SBATCH --account=torch_pr_138_general

TARGET=$1
CURRENT_BATCH=${2:-48}
TARGET_DIR="$SCRATCH/alphafast/$TARGET"
THRESHOLD=0.58

# Basic safety guards before any destructive operation.
if [[ -z "$TARGET" ]]; then
    echo "ERROR: Missing target. Usage: $0 <TARGET> [CURRENT_BATCH]"
    exit 1
fi

if [[ -z "${SCRATCH:-}" ]]; then
    echo "ERROR: SCRATCH is not set. Refusing to continue."
    exit 1
fi

if [[ ! -d "$TARGET_DIR" ]]; then
    echo "ERROR: Target directory does not exist: $TARGET_DIR"
    exit 1
fi

case "$TARGET_DIR" in
    "$SCRATCH"/alphafast/*) ;;
    *)
        echo "ERROR: Unsafe TARGET_DIR outside expected root: $TARGET_DIR"
        exit 1
        ;;
esac

for cmd in jq bc rsync sbatch; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: $cmd"
        exit 1
    fi
done

# Dynamically count the target number of JSON inputs
GOAL=$(find "$TARGET_DIR/inputs" -maxdepth 1 -type f 2>/dev/null | wc -l | tr -d ' ')
if [ -z "$GOAL" ] || [ "$GOAL" -eq 0 ]; then GOAL=19699; fi # Fallback

echo "Starting cleanup for $TARGET (Current Batch: $CURRENT_BATCH) at: $(date)"

export THRESHOLD

process_folder() {
    local marker="$1"
    local dir=$(dirname "$marker")
    
    # 0. Skip if already verified
    if [[ -f "$dir/completion.txt" ]]; then 
        return 0
    fi
    
    # 1. Delete temp seed folders
    find "$dir" -maxdepth 1 -type d -name "seed*" -exec rm -rf {} +

    # 2. Quality Filter
    local json_file=$(find "$dir" -maxdepth 1 -name "*summary_confidences.json" | head -n 1)
    if [[ -f "$json_file" ]]; then
        local iptm=$(jq -r '.iptm' "$json_file")
        if (( $(echo "$iptm >= $THRESHOLD" | bc -l) )); then
            echo "[PASS] $dir - iPTM: $iptm"
            touch "$dir/completion.txt"
        else
	    echo "[REJECT] $dir - iPTM $iptm is below threshold. Wiping results to save space, but marking as complete."
            find "$dir" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
            touch "$dir/completion.txt"
        fi
    else
        echo "[FAIL] $dir - No summary JSON found. Wiping for retry."
        rm -rf "$dir"
    fi
}
export -f process_folder

# Find ALL completed run folders. The process_folder function already early-exits for verified ones.
find "$TARGET_DIR" -type f -name "TERMS_OF_USE.md" -print0 | \
    xargs -0 -P $SLURM_CPUS_PER_TASK -I {} bash -c 'process_folder "$@"' _ {}

# Catch-all: Remove any stray directories from failed jobs (e.g. OOM/timeout before TERMS_OF_USE.md)
# This finds any 'outputs_task_' directory and checks its subdirectories.
find "$TARGET_DIR" -maxdepth 1 -name "outputs_task_*" -type d 2>/dev/null | while read -r task_dir; do
    for struct_dir in "$task_dir"/*; do
        if [[ -d "$struct_dir" ]] && [[ ! -f "$struct_dir/completion.txt" ]]; then
            echo "[WIPE] Cleaning up incomplete output directory: $struct_dir"
            rm -rf "$struct_dir"
        fi
    done
done

# --- Run PPI analysis ---
echo "Running PPI analysis on completed structures..."
JOB_ID=$(sbatch --parsable --job-name="ppi_${TARGET}" $HOME/alphafast_ppi.slurm $TARGET)

# --- Count verified completions ---
count=$(find "$TARGET_DIR" -name "completion.txt" | wc -l)

if (( count >= GOAL )); then
    echo "SUCCESS: $count/$GOAL verified folders complete."
    exit 0
fi

# --- Re-submission Logic ---

remain=$(( GOAL - count ))

# Batch size step-down logic
if [ "$remain" -ge 10000 ]; then
    NEXT_BATCH=24
elif [ "$remain" -ge 5000 ]; then
    NEXT_BATCH=16
elif [ "$remain" -ge 1000 ]; then
    NEXT_BATCH=8
elif [ "$remain" -ge 200 ]; then
    NEXT_BATCH=2
else
    NEXT_BATCH=1
fi

echo "Progress: $count/$GOAL, $remain remain. Re-submitting with Batch Size: $NEXT_BATCH"
sbatch --dependency=afterany:$JOB_ID --job-name="sub_${TARGET}" $HOME/alphafast_submit.slurm "$TARGET" "$NEXT_BATCH"

rm -f "$HOME"/logs/02* 
rm -f "$HOME"/logs/pipeline* 
rm -f "$HOME"/logs/inference*
