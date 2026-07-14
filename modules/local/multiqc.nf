// ============================================================
// MultiQC — 汇总质控报告
// 输入: 所有 FastQC 报告的合集 (raw + trimmed)
// 输出: multiqc_report.html + multiqc_data/
// ============================================================

process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/05_multiqc", mode: 'copy'

    input:
    path(input_files)          // 所有 QC 报告（raw + trimmed 混合）

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data/"),       emit: data, optional: true

    script:
    """
    set -euo pipefail
    multiqc . --force --outdir .
    """

    stub:
    """
    mkdir -p multiqc_data && touch multiqc_report.html
    """
}
