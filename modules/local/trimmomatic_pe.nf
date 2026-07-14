// ============================================================
// Trimmomatic — 双端模式
// 输入: tuple(sample_id, r1, r2) + adapter_fasta
// 输出: tuple(sample_id, trimmed_r1, trimmed_r2)
//
// 注意: unpaired 文件不 emit（下游不需要），Nextflow 自动清理
// ============================================================

process TRIMMOMATIC_PE {
    label 'process_low'
    tag "${sample_id}"
    publishDir "${params.outdir}/02_trimmomatic", mode: 'copy', pattern: "*.trimmed.fastq.gz"

    input:
    tuple val(sample_id), path(r1), path(r2)
    path adapter_fasta

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.trimmed.fastq.gz"),
          path("${sample_id}_R2.trimmed.fastq.gz"), emit: trimmed_fq

    script:
    def phred_flag = params.trimmomatic.phred ? "-phred${params.trimmomatic.phred}" : ""
    """
    set -euo pipefail
    trimmomatic PE -threads ${task.cpus} ${phred_flag} \
        ${r1} ${r2} \
        ${sample_id}_R1.trimmed.fastq.gz ${sample_id}_R1.unpaired.fastq.gz \
        ${sample_id}_R2.trimmed.fastq.gz ${sample_id}_R2.unpaired.fastq.gz \
        ILLUMINACLIP:${adapter_fasta}:${params.trimmomatic.illuminaclip} \
        LEADING:${params.trimmomatic.leading} TRAILING:${params.trimmomatic.trailing} \
        SLIDINGWINDOW:${params.trimmomatic.slidingwindow} MINLEN:${params.trimmomatic.minlen}
    """

    stub:
    """
    touch ${sample_id}_R1.trimmed.fastq.gz ${sample_id}_R2.trimmed.fastq.gz
    """
}
