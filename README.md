# Co-Folding Against Human Proteome; Protein-Protein Interactions (PPI)

This pipeline will download all human protein sequences from [STRING](https://string-db.org/cgi/download?sessionId=bVGGWIlTNLo8&species_text=Homo+sapiens&settings_expanded=0&min_download_score=0&filter_redundant_pairs=0&delimiter_type=txt) and create json files for each target sequence vs each protein sequence in the proteome as input for tools such as [AlphaFold](https://github.com/google-deepmind/alphafold3/tree/main) or [Boltz](https://github.com/jwohlwend/boltz).

Create a targets.tsv file with the following columns:

- id (name of your target sequence)
- type (protein, DNA, RNA, ligand)
- sequence

Example:

```{tsv}
id	type	sequence
p1	protein	ILAMKALGA
d1	DNA	GATTACA
```

```{bash}
module load nextflow # or use conda

nextflow run main.nf \
    -profile slurm \            # use this only if on a SLURM HPC
    --targets all_targets.tsv \ # custom input, default: targets.tsv
    --filetype yaml             # json or yaml, default: json
```

The result: numerous json/yaml files in a `inputs.tar.gz` file (or direct results if `--fold boltz` is used).

## Boltz Genetic/Structural Prediction Setup

If you are using `--fold boltz`, follows these steps to ensure the environment is ready:

### 1. Container Registry Login
Boltz-2 images are hosted on NVIDIA's registry. You must log in via Singularity/Apptainer:

```bash
# Get an API key from NVIDIA NGC and use it as the token
singularity registry login --username '$oauthtoken' docker://nvcr.io
pull_boltz_image.sh
```

### 2. Environment Variables
Ensure you have a designated cache directory for Singularity images:

```bash
export NXF_SINGULARITY_CACHEDIR="${PWD}/singularity_cache"
```

### 3. Running the Pipeline
The `run_nf.sh` script is provided for easy execution on SLURM clusters. It handles cleaning the cache and submitting the job:

```bash
sbatch run_nf.sh
```

### 4. Automated Cache Management
The pipeline automatically manages the Boltz cache (`boltz_cache/`). On the first run, it will:
- Download the correct Boltz-2 weights.
- Download and manually unpack the Chemical Component Dictionary (CCD) into `mols/`.
- Ensure RDKit atom properties are preserved to avoid structural featurization errors.

If you encounter issues with "Missing ALA" or "KeyError: 'name'", verify that `boltz_cache/mols` is populated with `.pkl` files. You can force a rebuild by deleting the `boltz_cache` directory.

## Requirements
- Nextflow (DSL2)
- Singularity / Apptainer
- NVIDIA GPUs with compatible drivers