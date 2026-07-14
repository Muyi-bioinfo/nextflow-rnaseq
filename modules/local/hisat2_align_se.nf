// ============================================================
// HISAT2 比对 — 单端模式
// 输入: tuple(sample_id, trimmed_read) + tuple(prefix, idx_files)
// 输出: tuple(sample_id, sorted.bam, sorted.bam.bai) + hisat2.log
//
// 注: samtools 使用一半线程，避免与 hisat2 管道两端合计超配
// ============================================================

process HISAT2_ALIGN_SE {
    label 'process_high'
    tag "${sample_id}"
    publishDir "${params.outdir}/03_hisat2", mode: 'copy'

    input:
    tuple val(sample_id), path(read), val(prefix), path(idx_files)

    output:
    tuple val(sample_id),
          path("${sample_id}.sorted.bam"),
          path("${sample_id}.sorted.bam.bai"), emit: bam_bai
    path("${sample_id}.hisat2.log"),           emit: hisat2_log

    script:
    // samtools 固定小线程数，hisat2 用满配（管道两端不翻倍）
    def samtools_cpus = (task.cpus as int) > 2 ? 2 : 1
    """
    set -euo pipefail
    hisat2 -x ${prefix} \
        -U ${read} \
        -p ${task.cpus} ${params.hisat2.extra} \
        2> ${sample_id}.hisat2.log \
        | samtools sort -@ ${samtools_cpus} -o ${sample_id}.sorted.bam -
    samtools index ${sample_id}.sorted.bam
    """

    stub:
    """
    touch ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai ${sample_id}.hisat2.log
    """
}
