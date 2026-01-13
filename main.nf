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
    path "${target_id}_${params.filetype}.tar.gz", emit: inputs_archive
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

workflow {
    // 1. Prepare targets channel from TSV
    targets_ch = Channel.fromPath(params.targets)
        .splitCsv(sep: '\t', header: true)
        .map { row -> tuple(row.id, row.type, row.sequence) }

    // 2. Prepare sequences channel (full file)
    seqs_ch = DL_STRING_SEQS()
    
    // 3. Parallel execution: Each target x Full sequence file
    GENERATE_INPUTS(targets_ch.combine(seqs_ch))

    GENERATE_INPUTS.out.lengths.collect().set { all_lengths }
    MERGE_LENGTHS(all_lengths)
}
