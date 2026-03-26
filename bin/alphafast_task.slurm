#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48GB
#SBATCH --gres=gpu:1
#SBATCH --constraint=a100|h200|l40s|rtx6000|h100
#SBATCH --comment="preemption=yes;requeue=false;preemption_partitions_only=yes"
#SBATCH --time=3:00:00
#SBATCH --output=logs/02_%x_%A_%a.out
#SBATCH --account=torch_pr_138_general

DIR=$1
BATCH_SIZE=$2
BASE_DIR="$SCRATCH/alphafast/$DIR"
BASE_INPUT_DIR="$BASE_DIR/inputs"
MANIFEST="$BASE_DIR/todo_manifest.txt"

# Task-specific folders
TASK_INPUT_DIR="$BASE_DIR/inputs_task_${SLURM_ARRAY_TASK_ID}"
TASK_OUTPUT_DIR="$BASE_DIR/outputs_task_${SLURM_ARRAY_TASK_ID}"

# 1. Slice the manifest
START=$((SLURM_ARRAY_TASK_ID * BATCH_SIZE + 1))
END=$(((SLURM_ARRAY_TASK_ID + 1) * BATCH_SIZE))
FILES_TO_PROCESS=$(sed -n "${START},${END}p" "$MANIFEST")

if [ -z "$FILES_TO_PROCESS" ]; then
    echo "No work assigned for Task ID $SLURM_ARRAY_TASK_ID. Exiting."
    exit 0
fi

rm -rf "$TASK_INPUT_DIR"
mkdir -p "$TASK_INPUT_DIR" "$TASK_OUTPUT_DIR"

for f in $FILES_TO_PROCESS; do
    cp "$BASE_INPUT_DIR/$f" "$TASK_INPUT_DIR/"
done

# 2. Run AlphaFast
AF_EXE="/projects/work/public/apps/alphafast/20260221/scripts/run_alphafast.sh"
$AF_EXE \
    --input_dir "$TASK_INPUT_DIR" \
    --output_dir "$TASK_OUTPUT_DIR" \
    --num_gpus 1 \
    --gpu_devices 0

# 3. Clean up scratch inputs (Leave outputs for cleanup script)
rm -rf "$TASK_INPUT_DIR"
