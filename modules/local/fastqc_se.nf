// ============================================================
// FastQC — 单端模式
// 输入: tuple(sample_id, read)
// 输出: tuple(sample_id, [html_files], [zip_files])
// ============================================================

process FASTQC_SE {
    label 'process_low'
    tag "${sample_id}"
    publishDir "${params.outdir}/01_fastqc", mode: 'copy', pattern: "*_fastqc.{html,zip}"

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id), path("*_fastqc.html"), path("*_fastqc.zip"), emit: reports

    script:
    """
    set -euo pipefail
    fastqc ${read} --outdir . --threads ${task.cpus} ${params.fastqc.extra}
    """

    stub:
    def read_base = read.name.replaceFirst(/\.fastq\.gz\$/, '')
    """
    touch ${read_base}_fastqc.html ${read_base}_fastqc.zip
    """
}
