// ============================================================
// DESeq2 — 差异表达分析（下游）
// 输入: featureCounts 计数矩阵 + GTF + 分组/对比 JSON + R 脚本
// 输出: *_all_results.csv, *_significant.csv, *.pdf
// ============================================================

process DESEQ2 {
    label 'process_R'
    publishDir "${params.outdir}/06_deseq2", mode: 'copy'

    input:
    path(counts)
    path(gtf)
    val(groups_json)
    val(contrasts_json)
    path(common_r)
    path(deseq2_r)
    path(gene2symbol_r)

    output:
    path("*.csv"), emit: tables
    path("*.pdf"), emit: plots
    path("*.png"), optional: true

    script:
    """
    set -euo pipefail
    Rscript ${deseq2_r} \\
        --counts ${counts} \\
        --groups '${groups_json}' \\
        --contrasts '${contrasts_json}' \\
        --gtf ${gtf} \\
        --outdir . \\
        --padj ${params.deseq2.padj_threshold} \\
        --log2fc ${params.deseq2.log2fc_threshold} \\
        --top_n ${params.deseq2.top_n_genes}
    """

    stub:
    """
    touch PCA_plot.pdf sample_distance_heatmap.pdf DEG_heatmap.pdf
${params.deseq2.contrasts.collect { c -> "    touch ${c.name}_all_results.csv ${c.name}_significant.csv ${c.name}_MA_plot.pdf ${c.name}_volcano_plot.pdf" }.join('\n')}
    """
}
