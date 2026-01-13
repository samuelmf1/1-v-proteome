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

The result: numerous json files in a `inputs.tar.gz` file, ready to be run with downstream tools.