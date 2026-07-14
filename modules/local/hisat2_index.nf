// ============================================================
// HISAT2 索引构建
// 输入: 参考基因组 FASTA
// 输出: .ht2 索引文件 (8 个) + 索引 basename（供下游 -x 参数）
// ============================================================

process HISAT2_INDEX {
    label 'process_medium'
    tag "${genome_fasta.simpleName}"
    cache 'deep'

    input:
    path genome_fasta

    output:
    tuple val("${genome_fasta.simpleName}"), path("*.ht2"), emit: index

    script:
    idx = genome_fasta.simpleName
    """
    set -euo pipefail
    hisat2-build -p ${task.cpus} ${genome_fasta} ${idx}
    """

    stub:
    idx = genome_fasta.simpleName
    """
    touch ${idx}.1.ht2 ${idx}.2.ht2 ${idx}.3.ht2 ${idx}.4.ht2 \
          ${idx}.5.ht2 ${idx}.6.ht2 ${idx}.7.ht2 ${idx}.8.ht2
    """
}
