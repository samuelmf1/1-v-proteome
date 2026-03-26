#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process DL_STRING_SEQS {
    output:
    path "string_seqs.fa.gz"

    script:
    """
    wget -O string_seqs.fa.gz ${params.string_seq_url}
    """
}

process GENERATE_INPUTS {
    tag "${target_id}"
    cache false
    publishDir "${params.outdir}/inputs_${params.filetype}", mode: 'copy'

    input:
    tuple val(target_id), val(target_type), val(target_seq), path(seqs_file)

    output:
    tuple val(target_id), path("${target_id}_${params.filetype}.tar.gz"), emit: inputs_archive
    path "seq_lengths.tsv", emit: lengths

    script:
    """
    mkdir -p ${target_id}_inputs
    generate_seq_inputs.py \\
        --target_id "${target_id}" \\
        --target_type "${target_type}" \\
        --target_seq "${target_seq}" \\
        --seqs ${seqs_file} \\
        --filetype ${params.filetype} \\
        --outdir ${target_id}_inputs
    
    # Move the lengths file out and rename if necessary to match expectation
    # Python script generates 'af3_lengths.tsv'
    mv ${target_id}_inputs/af3_lengths.tsv seq_lengths.tsv
    
    # Tar the input directory
    # -C changes to directory so files are at root of tarball
    tar -czf ${target_id}_${params.filetype}.tar.gz -C ${target_id}_inputs .
    """
}

process MERGE_LENGTHS {
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path "lengths_*"

    output:
    path "seq_lengths.tsv"
    path "seq_stats.txt"

    script:
    """
    merge_seq_stats.py --inputs lengths_* --output_tsv seq_lengths.tsv --output_stats seq_stats.txt
    """
}


process GENERATE_MONOMER_INPUTS {
    tag "monomers"
    publishDir "${params.outdir}/monomers_${params.filetype}", mode: 'copy'

    input:
    path targets_file

    output:
    path "monomers_${params.filetype}.tar.gz", emit: monomer_archive
    path "monomers/*.${params.filetype}", emit: monomer_yamls, optional: true

    script:
    """
    mkdir -p monomers
    generate_seq_inputs.py \\
        --targets ${targets_file} \\
        --filetype ${params.filetype} \\
        --outdir monomers

    tar -czf monomers_${params.filetype}.tar.gz -C monomers .
    """
}

process SPARROW_FEATURES {
    tag "${yaml_file.baseName}"
    publishDir "${params.outdir}/sparrow", mode: 'copy'

    input:
    path yaml_file

    output:
    path "${yaml_file.baseName}_sparrow_features.tsv"
    path "${yaml_file.baseName}_sparrow_features_residues.tsv"

    script:
    """
    predict_sparrow.py --yaml ${yaml_file} --out ${yaml_file.baseName}_sparrow_features.tsv
    """
}


process PREPARE_INPUT {
    input:
    path input_file
    
    output:
    path "input.fasta"
    
    script:
    if (input_file.name.endsWith('.tsv') || input_file.name.endsWith('.csv')) {
        def sep = input_file.name.endsWith('.tsv') ? '\t' : ','
        """
        python3 -c "
        import csv
        with open('${input_file}', 'r') as f:
            reader = csv.DictReader(f, delimiter='${sep}')
            with open('input.fasta', 'w') as out:
                for row in reader:
                    out.write('>{}\\n{}\\n'.format(row['id'], row['sequence']))
        "
        """
    } else if (input_file.name.endsWith('.gz')) {
        """
        gunzip -c ${input_file} > input.fasta
        """
    } else {
        """
        cp ${input_file} input.fasta
        """
    }
}

process SPLIT_INPUT {
    input:
    path input_fasta
    
    output:
    path "chunk_*.fasta"
    
    script:
    """
    total=\$(grep -c '^>' ${input_fasta})
    # Calculate chunk size to get at most 3 chunks
    chunk_size=\$(( (total + 2) / 3 ))
    if [ \$chunk_size -le 0 ]; then chunk_size=1; fi
    
    # Split using awk
    awk -v sz="\$chunk_size" '
        /^>/ {
            if (n % sz == 0) {
                i++;
                if (out) close(out);
                out = "chunk_" i ".fasta";
            }
            n++;
        }
        { if (out) print > out; }
    ' ${input_fasta}
    """
}

process GENERATE_MSA {
    tag "${meta.id}"
    label 'process_high_memory'
    
    publishDir "${params.msa_out}", mode: 'copy'

    container 'docker://ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2'

    input:
    tuple val(meta), path(fasta)
    path db_dir

    output:
    path "*.a3m", emit: a3m

    script:
    """
    colabfold_search \
        ${fasta} \
        ${db_dir} \
        . \
        --threads ${task.cpus} \
        --db-load-mode 2 \
        --db1 uniref30_2302_db \
        --db3 colabfold_envdb_202108_db
    """
}

process CREATE_SQUASHFS {
    tag "${target_id}"
    publishDir "${params.outdir}/squashfs_${params.filetype}", mode: 'copy'

    input:
    tuple val(target_id), path(tarball)

    output:
    path "${target_id}-af3-pipeline-input-${params.filetype}.sqf"

    script:
    """
    mkdir -p input_${params.filetype}
    tar -xzf ${tarball} -C input_${params.filetype}
    mksquashfs input_${params.filetype} ${target_id}-af3-pipeline-input-${params.filetype}.sqf -all-root -keep-as-directory -no-xattrs
    """
}

process SUBMIT_ALPHAFAST {
    tag "serialized"

    input:
    path sqf_files

    output:
    path "submit_alphafast.done"

    script:
    """
    first_sqf=\$(printf "%s\n" ${sqf_files} | sort | head -n 1)
    if [ -z "\$first_sqf" ]; then
        echo "No squashfs files found for AlphaFast submission" >&2
        exit 1
    fi

    structure=\$(basename "\$first_sqf" ".sqf")
    structure=\${structure%-af3-pipeline-input-*}
    echo "Submitting initial AlphaFast structure: \$structure"
    sbatch ${projectDir}/bin/alphafast_submit.slurm "\$structure"

    touch submit_alphafast.done
    """
}

workflow {
    // 1. Prepare targets channel from TSV
    targets_ch = Channel.fromPath(params.targets)
        .splitCsv(sep: '\t', header: true)
        .map { row -> tuple(row.id, row.type, row.sequence) }

    // Generate Inputs for Monomers (always runs)
    GENERATE_MONOMER_INPUTS(file(params.targets))

    // Run SPARROW feature prediction for monomers (always runs)
    SPARROW_FEATURES(GENERATE_MONOMER_INPUTS.out.monomer_yamls.flatten())

    if (!params.skip_tarbell) {
        // 2. Prepare sequences channel (full file)
        seqs_ch = DL_STRING_SEQS()
        
        // 3. Parallel execution: Each target x Full sequence file
        GENERATE_INPUTS(targets_ch.combine(seqs_ch))
    
        GENERATE_INPUTS.out.lengths.collect().set { all_lengths }
        MERGE_LENGTHS(all_lengths)
        
        GENERATE_INPUTS.out.inputs_archive.set { inputs_archive_ch }
    } else {
        // Create inputs_archive_ch from existing files
        targets_ch
            .map { id, type, seq -> 
                def tarball_path = file("${params.outdir}/${id}_${params.filetype}.tar.gz")
                if (!tarball_path.exists()) {
                    error "Input archive not found: ${tarball_path}"
                }
                return tuple(id, tarball_path)
            }
            .set { inputs_archive_ch }
    }
    
    if (params.squashfs) {
        def monomer_archive_ch = GENERATE_MONOMER_INPUTS.out.monomer_archive
            .map { archive -> tuple("monomers", archive) }
        CREATE_SQUASHFS(inputs_archive_ch.mix(monomer_archive_ch)).set { squashfs_ch }
        if (params.alphafast) {
            SUBMIT_ALPHAFAST(squashfs_ch.collect())
        }
    }
}