// ============================================================
// 子流程 — SE 模式：QC → 剪切 → QC → 比对
//
// take: sample_ch      (sample_id, read)
//       adapter_ch     (adapter_fasta)
//       index_ch       (prefix, idx_files) — 单元素 channel
//
// emit: bam_bai        (sample_id, bam, bai)
//       qc_files       (all QC reports, flat list)
// ============================================================

include { FASTQC_SE as FASTQC_RAW_SE      } from '../../modules/local/fastqc_se'
include { FASTQC_SE as FASTQC_TRIMMED_SE  } from '../../modules/local/fastqc_se'
include { TRIMMOMATIC_SE     } from '../../modules/local/trimmomatic_se'
include { HISAT2_ALIGN_SE    } from '../../modules/local/hisat2_align_se'

workflow QC_TRIM_ALIGN_SE {
    take:
    sample_ch
    adapter_ch
    index_ch

    main:
    // 1. 原始 FastQC
    FASTQC_RAW_SE(sample_ch)

    // 2. Trimmomatic
    TRIMMOMATIC_SE(sample_ch, adapter_ch)

    // 3. 剪切后 FastQC
    FASTQC_TRIMMED_SE(TRIMMOMATIC_SE.out.trimmed_fq)

    // 4. HISAT2 比对 — use combine() so every sample pairs with the single index
    TRIMMOMATIC_SE.out.trimmed_fq
        .combine(index_ch)
        .map { sid, read, pfx, idx -> tuple(sid, read, pfx, idx) }
        .set { align_in }

    HISAT2_ALIGN_SE(align_in)

    // 汇总所有 QC 文件
    qc_files = FASTQC_RAW_SE.out.reports
        .mix(FASTQC_TRIMMED_SE.out.reports)
        .map { _sample, htmls, zips -> htmls + zips }
        .flatten()
        .collect()

    emit:
    bam_bai  = HISAT2_ALIGN_SE.out.bam_bai
    qc_files = qc_files
}
