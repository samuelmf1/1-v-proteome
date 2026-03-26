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

process UNPACK_INPUTS {
    tag "${target_id}"

    input:
    tuple val(target_id), path(inputs_archive)

    output:
    tuple val(target_id), path("inputs/*.yaml"), emit: input_files, optional: true

    script:
    """
    mkdir -p inputs
    tar -xzf ${inputs_archive} -C inputs
    
    # Filter existing results
    RESULTS_DIR="${params.outdir}/${target_id}_results"
    if [ -d "\$RESULTS_DIR" ]; then
        find inputs -name "*.yaml" | while read f; do
            BASE=\$(basename "\$f" .yaml)
            if [ -d "\$RESULTS_DIR/\$BASE" ]; then
                rm "\$f"
            fi
        done
    fi
    """
}

process SETUP_BOLTZ_CACHE {
    // Run once to populate cache (weights, CCD)
    
    input:
    path cache_dir
    
    output:
    path "boltz_ready", emit: ready
    
    script:
    """
    setup_boltz.sh ${cache_dir}
    """
}

process BOLTZ_PREDICT {
    tag "${target_id}_part${task.index}"

    
    input:
    tuple val(target_id), path(yaml_files)
    path cache_dir
    val ready // Dependency signal
    
    output:
    path "results/**/*"
    
    script:
    """
    predict_boltz.sh ${cache_dir} ${params.filetype} "${params.outdir}/${target_id}_results" ${yaml_files}
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

process BOLTZ_PREDICT_MONOMERS {
    tag "monomers"
    
    input:
    path yaml_files
    path cache_dir
    val ready // Dependency signal
    
    output:
    path "results/**/*"
    
    script:
    """
    predict_boltz.sh ${cache_dir} yaml "${params.outdir}/monomers/results" ${yaml_files}
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

workflow {
    // 1. Prepare targets channel from TSV
    targets_ch = Channel.fromPath(params.targets)
        .splitCsv(sep: '\t', header: true)
        .map { row -> tuple(row.id, row.type, row.sequence) }

    // Generate Inputs for Monomers (always runs)
    GENERATE_MONOMER_INPUTS(file(params.targets))

    // Run SPARROW feature prediction for monomers (always runs)
    SPARROW_FEATURES(GENERATE_MONOMER_INPUTS.out.monomer_yamls.flatten())

    // Run Monomer Fold Predictions (Boltz only)
    if (params.fold == "boltz") {
         // Define central cache directory
        def boltz_cache = file("${workflow.launchDir}/boltz_cache")
        
        // Run setup (download weights) ONCE
        SETUP_BOLTZ_CACHE(boltz_cache)

        // Run Prediction for Monomers
        BOLTZ_PREDICT_MONOMERS(GENERATE_MONOMER_INPUTS.out.monomer_yamls, boltz_cache, SETUP_BOLTZ_CACHE.out.ready)
    }

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
        CREATE_SQUASHFS(inputs_archive_ch.mix(monomer_archive_ch))
    }
    
    if (params.fold == "boltz") {
        UNPACK_INPUTS(inputs_archive_ch)
        
        // Batch files into chunks of 2000 (~20 jobs for 20k files)
        UNPACK_INPUTS.out.input_files
            .flatMap { id, files -> 
                if (files) {
                    // Ensure files is a list
                    def fileList = new ArrayList(files instanceof List ? files : [files])
                    java.util.Collections.shuffle(fileList)

                    if (fileList.isEmpty()) return []
                    
                    def batchSize = (int) Math.ceil(fileList.size() / 10.0)
                    def batches = fileList.collate(batchSize)
                    batches.collect { batch -> tuple(id, batch) }
                } else {
                    return []
                }
            }
            .set { boltz_batches_ch }
        
        // Run predictions using the populated cache
        // Note: SETUP_BOLTZ_CACHE and cache dir definition moved up
        def boltz_cache = file("${workflow.launchDir}/boltz_cache")
        BOLTZ_PREDICT(boltz_batches_ch, boltz_cache, SETUP_BOLTZ_CACHE.out.ready)
    }
}

workflow msa {
    if (params.msa_input) {
        msa_input_ch = Channel.fromPath(params.msa_input)
    } else if (params.input) {
        msa_input_ch = Channel.fromPath(params.input)
    } else {
        msa_input_ch = DL_STRING_SEQS()
    }

    db_path = file(params.msa_db)

    fasta_file = PREPARE_INPUT(msa_input_ch)
    chunks_ch = SPLIT_INPUT(fasta_file).flatten()
    
    input_ch_ready = chunks_ch.map { f -> tuple([id: f.baseName], f) }

    GENERATE_MSA(input_ch_ready, db_path)
}