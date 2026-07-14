// ============================================================
// Trimmomatic — 单端模式
// 输入: tuple(sample_id, read) + adapter_fasta
// 输出: tuple(sample_id, trimmed_read)
// ============================================================

process TRIMMOMATIC_SE {
    label 'process_low'
    tag "${sample_id}"
    publishDir "${params.outdir}/02_trimmomatic", mode: 'copy', pattern: "*.trimmed.fastq.gz"

    input:
    tuple val(sample_id), path(read)
    path adapter_fasta

    output:
    tuple val(sample_id),
          path("${sample_id}.trimmed.fastq.gz"), emit: trimmed_fq

    script:
    def phred_flag = params.trimmomatic.phred ? "-phred${params.trimmomatic.phred}" : ""
    """
    set -euo pipefail
    trimmomatic SE -threads ${task.cpus} ${phred_flag} \
        ${read} ${sample_id}.trimmed.fastq.gz \
        ILLUMINACLIP:${adapter_fasta}:${params.trimmomatic.illuminaclip} \
        LEADING:${params.trimmomatic.leading} TRAILING:${params.trimmomatic.trailing} \
        SLIDINGWINDOW:${params.trimmomatic.slidingwindow} MINLEN:${params.trimmomatic.minlen}
    """

    stub:
    """
    touch ${sample_id}.trimmed.fastq.gz
    """
}
