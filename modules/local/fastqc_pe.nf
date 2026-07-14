// ============================================================
// FastQC — 双端模式
// 输入: tuple(sample_id, r1, r2)
// 输出: tuple(sample_id, [html_files], [zip_files])
//
// 注: 输出用 glob pattern，自动适配 raw (*_R1_fastqc.html) 和
//     trimmed (*_R1.trimmed_fastqc.html) 两种命名
// ============================================================

process FASTQC_PE {
    label 'process_low'
    tag "${sample_id}"
    publishDir "${params.outdir}/01_fastqc", mode: 'copy', pattern: "*_fastqc.{html,zip}"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("*_fastqc.html"), path("*_fastqc.zip"), emit: reports

    script:
    """
    set -euo pipefail
    fastqc ${r1} ${r2} --outdir . --threads ${task.cpus} ${params.fastqc.extra}
    """

    stub:
    def r1_base = r1.name.replaceFirst(/\.fastq\.gz\$/, '')
    def r2_base = r2.name.replaceFirst(/\.fastq\.gz\$/, '')
    """
    touch ${r1_base}_fastqc.html ${r1_base}_fastqc.zip
    touch ${r2_base}_fastqc.html ${r2_base}_fastqc.zip
    """
}
