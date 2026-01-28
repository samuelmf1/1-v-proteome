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
    publishDir "${params.outdir}", mode: 'copy'

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
    tuple val(target_id), path("inputs/*.yaml"), emit: input_files

    script:
    """
    mkdir -p inputs
    tar -xzf ${inputs_archive} -C inputs
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
    publishDir "${params.outdir}/monomers/inputs", mode: 'copy'

    input:
    path targets_file

    output:
    path "monomers/*.yaml", emit: monomer_yamls

    script:
    """
    mkdir -p monomers
    generate_seq_inputs.py \\
        --targets ${targets_file} \\
        --filetype yaml \\
        --outdir monomers
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

workflow {
    // 1. Prepare targets channel from TSV
    targets_ch = Channel.fromPath(params.targets)
        .splitCsv(sep: '\t', header: true)
        .map { row -> tuple(row.id, row.type, row.sequence) }

    // Run Monomer Predictions (Targets only)
    if (params.fold == "boltz") {
         // Define central cache directory
        def boltz_cache = file("${workflow.launchDir}/boltz_cache")
        
        // Run setup (download weights) ONCE
        SETUP_BOLTZ_CACHE(boltz_cache)

        // Generate Inputs for Monomers
        GENERATE_MONOMER_INPUTS(file(params.targets))

        // Run Prediction for Monomers
        // Split list of files into chunks if needed, but for now passing all
        // predict_boltz.sh can handle multiple files
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
    
    if (params.fold == "boltz") {
        UNPACK_INPUTS(inputs_archive_ch)
        
        // Batch files into chunks of 2000 (~20 jobs for 20k files)
        UNPACK_INPUTS.out.input_files
            .flatMap { id, files -> 
                def batches = files.collate(1600)
                batches.collect { batch -> tuple(id, batch) }
            }
            .set { boltz_batches_ch }
        
        // Run predictions using the populated cache
        // Note: SETUP_BOLTZ_CACHE and cache dir definition moved up
        def boltz_cache = file("${workflow.launchDir}/boltz_cache")
        BOLTZ_PREDICT(boltz_batches_ch, boltz_cache, SETUP_BOLTZ_CACHE.out.ready)
    }
}
