// ============================================================
// 子流程 — PE 模式：QC → 剪切 → QC → 比对
//
// take: sample_ch      (sample_id, r1, r2)
//       adapter_ch     (adapter_fasta)
//       index_ch       (prefix, idx_files) — 单元素 channel
//
// emit: bam_bai        (sample_id, bam, bai)
//       qc_files       (all QC reports, flat list)
// ============================================================

include { FASTQC_PE as FASTQC_RAW_PE      } from '../../modules/local/fastqc_pe'
include { FASTQC_PE as FASTQC_TRIMMED_PE  } from '../../modules/local/fastqc_pe'
include { TRIMMOMATIC_PE     } from '../../modules/local/trimmomatic_pe'
include { HISAT2_ALIGN_PE    } from '../../modules/local/hisat2_align_pe'

workflow QC_TRIM_ALIGN_PE {
    take:
    sample_ch
    adapter_ch
    index_ch

    main:
    // 1. 原始 FastQC
    FASTQC_RAW_PE(sample_ch)

    // 2. Trimmomatic
    TRIMMOMATIC_PE(sample_ch, adapter_ch)

    // 3. 剪切后 FastQC
    FASTQC_TRIMMED_PE(TRIMMOMATIC_PE.out.trimmed_fq)

    // 4. HISAT2 比对 — use combine() so every sample pairs with the single index
    TRIMMOMATIC_PE.out.trimmed_fq
        .combine(index_ch)
        .map { sid, r1, r2, pfx, idx -> tuple(sid, r1, r2, pfx, idx) }
        .set { align_in }

    HISAT2_ALIGN_PE(align_in)

    // 汇总所有 QC 文件
    qc_files = FASTQC_RAW_PE.out.reports
        .mix(FASTQC_TRIMMED_PE.out.reports)
        .map { _sample, htmls, zips -> htmls + zips }
        .flatten()
        .collect()

    emit:
    bam_bai  = HISAT2_ALIGN_PE.out.bam_bai
    qc_files = qc_files
}
