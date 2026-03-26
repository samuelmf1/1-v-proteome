#!/bin/bash
#SBATCH --job-name=nf-VLP                      # Job name
#SBATCH --mem=8G                               # Job memory request
#SBATCH --cpus-per-task=1                      # number of cpu per task
#SBATCH --time=298:00:00

# Load required modules
module load nextflow
module load singularity
export NXF_SINGULARITY_CACHEDIR="${PWD}/singularity_cache"

# Run the pipeline
rm -rf boltz_cache
# nextflow run main.nf -profile slurm --targets example.tsv --fold boltz --filetype yaml
nextflow run main.nf -profile slurm --targets targets.tsv --filetype json --squashfs