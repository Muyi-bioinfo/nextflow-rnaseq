// ============================================================
// clusterProfiler — GO + KEGG 功能富集分析（下游）
// 输入: DESeq2 显著差异基因 CSV + 对比名称列表 + R 脚本
// 输出: *_GO_enrichment.csv, *_GO_dotplot.pdf, *_KEGG_enrichment.csv, *_KEGG_dotplot.pdf
// ============================================================

process CLUSTERPROFILER {
    label 'process_R'
    publishDir "${params.outdir}/07_clusterprofiler", mode: 'copy'

    input:
    path(sig_files)            // 所有 _significant.csv 文件
    val(contrast_names_str)    // 逗号分隔的对比名称
    path(common_r)
    path(clusterprofiler_r)

    output:
    path("*.csv"), emit: tables, optional: true
    path("*.pdf"), emit: plots, optional: true

    script:
    def sig_str = (sig_files instanceof List ? sig_files : [sig_files]).join(',')
    """
    set -euo pipefail
    Rscript ${clusterprofiler_r} \\
        --sig_files '${sig_str}' \\
        --contrast_names '${contrast_names_str}' \\
        --org_db ${params.clusterprofiler.org_db} \\
        --kegg_org ${params.clusterprofiler.kegg_organism} \\
        --from_type ${params.clusterprofiler.from_type} \\
        --pval ${params.clusterprofiler.pvalue_cutoff} \\
        --qval ${params.clusterprofiler.qvalue_cutoff} \\
        --gene_id_col ${params.clusterprofiler.gene_id_col} \\
        --show_cat ${params.clusterprofiler.show_category} \\
        --outdir .
    """

    stub:
    """
${contrast_names_str.split(',').collect { name -> "    touch ${name.trim()}_GO_enrichment.csv ${name.trim()}_GO_dotplot.pdf ${name.trim()}_KEGG_enrichment.csv ${name.trim()}_KEGG_dotplot.pdf" }.join('\n')}
    """
}
