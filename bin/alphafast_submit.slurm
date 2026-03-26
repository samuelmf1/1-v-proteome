#!/bin/bash
#SBATCH --job-name=af_submit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --output=logs/01_%x_%j.out
#SBATCH --account=torch_pr_138_general

DIR=$1
BATCH_SIZE=${2:-32}

if [ -z "$DIR" ]; then
    echo "Usage: sbatch alphafast_submit.slurm <DIR_NAME> [BATCH_SIZE]"
    exit 1
fi

BASE_DIR="$SCRATCH/alphafast/$DIR"
BASE_INPUT_DIR="$BASE_DIR/inputs"
MANIFEST="$BASE_DIR/todo_manifest.txt"
TAR_SOURCE="/home/sf4479/alphafast/af3_inputs_json_tar/${DIR}_json.tar.gz"

mkdir -p "$BASE_DIR" "$BASE_INPUT_DIR" "logs"

# 1. Extraction (Only if input dir is empty)
if [ ! -d "$BASE_INPUT_DIR" ] || [ -z "$(ls -A "$BASE_INPUT_DIR" 2>/dev/null)" ]; then
    echo "Extracting inputs from $TAR_SOURCE..."
    tar -zxvf "$TAR_SOURCE" -C "$BASE_INPUT_DIR"
fi

# 2. Manifest Generation 
echo "Scanning for completed structures..."
ls "$BASE_INPUT_DIR" | sort > all_files.tmp

# Find folders with completion.txt, get basename, and append .json
find "$BASE_DIR" -name "completion.txt" 2>/dev/null | xargs -I{} dirname {} | xargs -I{} basename {} | sort | sed 's/$/.json/' > finished_files.tmp

# Filter: keep files in 'all' that are NOT in 'finished'
grep -v -F -x -f finished_files.tmp all_files.tmp | grep -v af3_stats.txt | sort -R > "$MANIFEST"

COUNT=$(wc -l < "$MANIFEST")
rm all_files.tmp finished_files.tmp

if [ "$COUNT" -eq 0 ]; then
    echo "All structures completed and verified! No jobs to submit."
    exit 0
fi

echo "Prep complete. $COUNT files remaining. Current Batch: $BATCH_SIZE."

# Guard against bad input that can break arithmetic (e.g., 0 or non-numeric).
if ! [[ "$BATCH_SIZE" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: BATCH_SIZE must be a positive integer (got '$BATCH_SIZE')."
    exit 1
fi

# Force a single array task when batch size is 1.
if [ "$BATCH_SIZE" -eq 1 ]; then
    NUM_JOBS=0
else
    # Zero-based max array index for --array=0-N when COUNT > 0.
    NUM_JOBS=$(( (COUNT - 1) / BATCH_SIZE ))
fi

# Launch the Worker Array
JOB_ID=$(sbatch --parsable --job-name="${DIR}_task" --array=0-$NUM_JOBS alphafast_task.slurm "$DIR" "$BATCH_SIZE")

echo "Worker array submitted with Job ID: $JOB_ID"

# Schedule Cleanup
sbatch --dependency=afterany:$JOB_ID --job-name="cln_${DIR}" alphafast_cleanup.slurm "$DIR"
