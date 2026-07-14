// ============================================================
// featureCounts — 基因表达定量（汇总所有 BAM）
// 输入: 所有样本的 BAM + GTF 注释文件
// 输出: featurecounts.txt + featurecounts.summary.txt
// ============================================================

process FEATURECOUNTS {
    label 'process_medium'
    publishDir "${params.outdir}/04_featurecounts", mode: 'copy'

    input:
    path(bams)          // 所有 BAM (collect 后传入)
    path(gtf)

    output:
    path("featurecounts.txt"),          emit: counts
    path("featurecounts.summary.txt"),  emit: summary

    script:
    def bam_list = (bams instanceof List ? bams : [bams]).join(' ')
    def paired_flag = params.read_mode == 'paired' ? '-p' : ''
    def extra = params.featurecounts.extra ?: ''
    """
    set -euo pipefail
    featureCounts -T ${task.cpus} ${paired_flag} \
        -t ${params.featurecounts.feature_type} \
        -g ${params.featurecounts.attr_type} \
        -s ${params.featurecounts.strandedness} \
        -a ${gtf} \
        -o featurecounts.txt \
        ${extra} \
        ${bam_list}
    mv featurecounts.txt.summary featurecounts.summary.txt
    """

    stub:
    """
    touch featurecounts.txt featurecounts.summary.txt
    """
}
