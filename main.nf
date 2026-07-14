#!/usr/bin/env nextflow
// ============================================================
// RNA-seq 上游分析流程 — Nextflow DSL2 主文件
//
// DAG:
//   Layer 0: hisat2_index
//   Layer 1-3: qc_trim_align (subworkflow: FASTQC→TRIMMOMATIC→FASTQC→HISAT2)
//   Layer 4: featureCounts
//   Layer 5: MultiQC
//   Layer 6: DESeq2
//   Layer 7: clusterProfiler
//
// 用法:
//   nextflow run main.nf -profile docker \
//       --input samplesheet.csv \
//       --genome_fasta genome.fa \
//       --annotation_gtf genes.gtf
// ============================================================

nextflow.enable.dsl=2

// ============================================================
// Module includes
// ============================================================

// Layer 0: 索引
include { HISAT2_INDEX          } from './modules/local/hisat2_index'

// Layer 1-3: per-sample QC + 剪切 + 比对 (封装为 subworkflow)
include { QC_TRIM_ALIGN_PE      } from './subworkflows/local/qc_trim_align_pe'
include { QC_TRIM_ALIGN_SE      } from './subworkflows/local/qc_trim_align_se'

// Layer 4-5: 汇总
include { FEATURECOUNTS         } from './modules/local/featurecounts'
include { MULTIQC               } from './modules/local/multiqc'

// Layer 6-7: 下游分析
include { DESEQ2                } from './modules/local/deseq2'
include { CLUSTERPROFILER       } from './modules/local/clusterprofiler'

// ============================================================
// Main workflow
// ============================================================

workflow {

    // ---- 启动前校验 ----
    if (!params.input)           { exit 1, "Missing required parameter: --input (samplesheet CSV)" }
    if (!params.genome_fasta)    { exit 1, "Missing required parameter: --genome_fasta" }
    if (!params.annotation_gtf)  { exit 1, "Missing required parameter: --annotation_gtf" }

    log.info "[main] Read mode: ${params.read_mode}"
    log.info "[main] Output dir: ${params.outdir}"

    // ---- 参考文件 channels ----
    genome_fasta_ch      = channel.value(file(params.genome_fasta))
    gtf_ch               = channel.value(file(params.annotation_gtf))
    adapter_ch           = params.adapter_fasta
                           ? channel.value(file(params.adapter_fasta))
                           : channel.value(file("/dev/null"))

    // ---- R 脚本 channels (作为 process input 传入，避免 process 内引用 projectDir) ----
    common_r_ch           = channel.value(file("${projectDir}/bin/common.R"))
    deseq2_r_ch           = channel.value(file("${projectDir}/bin/deseq2.R"))
    gene2symbol_r_ch      = channel.value(file("${projectDir}/bin/gene2symbol.R"))
    clusterprofiler_r_ch  = channel.value(file("${projectDir}/bin/clusterprofiler.R"))

    // ---- 样本输入 channel (单次解析，.multiMap() 分流) ----
    // SE 模式用 null 补齐为 4 元素 tuple，统一 multiMap 索引
    sample_input_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            if (params.read_mode == 'paired')
                tuple(row.sample, row.group, file(row.fastq_1), file(row.fastq_2))
            else
                tuple(row.sample, row.group, file(row.fastq_1), null)
        }
        .ifEmpty { exit 1, "No samples found in ${params.input}" }

    log.info "[main] Samplesheet loaded: ${params.input}"

    // 单次解析 → multiMap 分流 (每个 item 同时进入三个下游通道)
    sample_input_ch
        .multiMap { it ->
            fq:         tuple(it[0], it[2], it[3])   // sample, fastq_1, fastq_2
            group:      tuple(it[0], it[1])           // sample, group
            validation: it[1]                          // group (for contrast validation)
        }
        .set { forked }

    // 对比校验: 确保 case/control 在样本表中存在（通过 multiMap 分支取 group，无需重复读文件）
    forked.validation
        .unique()
        .collect()
        .map { valid_groups ->
            def group_set = valid_groups as Set
            params.deseq2.contrasts.each { c ->
                if (!group_set.contains(c['case']))
                    exit 1, "Contrast '${c.name}': case group '${c['case']}' not found in samplesheet. Available: ${group_set}"
                if (!group_set.contains(c.control))
                    exit 1, "Contrast '${c.name}': control group '${c.control}' not found in samplesheet. Available: ${group_set}"
            }
            log.info "[main] Contrast validation passed: groups=${group_set}, contrasts=${params.deseq2.contrasts*.name}"
        }

    // ---- Layer 0: HISAT2 基因组索引 ----
    // 用户提供了已建好的索引路径 → 直接复用；未提供 → 从 FASTA 构建
    if (params.hisat2_index_prefix && file("${params.hisat2_index_prefix}.1.ht2").exists()) {
        def idx_name = params.hisat2_index_prefix.tokenize('/')[-1]
        index_ch = channel.value(tuple(idx_name, files("${params.hisat2_index_prefix}.*.ht2")))
        log.info "[main] Using pre-built HISAT2 index: ${params.hisat2_index_prefix}"
    } else {
        HISAT2_INDEX(genome_fasta_ch)
        index_ch = HISAT2_INDEX.out.index
    }

    // ---- Layer 1-3: QC → 剪切 → QC → 比对 (PE/SE 分支) ----
    if (params.read_mode == 'paired') {
        result = QC_TRIM_ALIGN_PE(
            forked.fq.map { sid, fq1, fq2 -> tuple(sid, fq1, fq2) },
            adapter_ch, index_ch)
    } else {
        result = QC_TRIM_ALIGN_SE(
            forked.fq.map { sid, fq, _null -> tuple(sid, fq) },
            adapter_ch, index_ch)
    }
    bam_ch    = result.bam_bai
    qc_all_ch = result.qc_files

    // ---- Layer 4: featureCounts ----
    all_bams_ch = bam_ch
        .map { _sid, bam, _bai -> bam }
        .collect()
    FEATURECOUNTS(all_bams_ch, gtf_ch)

    // ---- Layer 5: MultiQC ----
    MULTIQC(qc_all_ch)

    // ---- Layer 6: DESeq2 ----
    // 分组 JSON — Groovy 读取样本表构建 (可靠，已验证)
    def _delim = params.input.endsWith('.csv') ? ',' : '\t'
    def _sample_groups = [:] as LinkedHashMap
    new File(params.input).withReader { r ->
        r.eachLine { line, idx ->
            if (idx == 1) return              // skip header
            def cols = line.split(_delim)
            def sid  = cols[0].trim()
            def grp  = cols[1].trim()
            if (!_sample_groups.containsKey(grp)) _sample_groups[grp] = []
            _sample_groups[grp] << sid
        }
    }
    groups_ch = channel.value(new groovy.json.JsonBuilder(_sample_groups).toString())

    contrasts_ch = channel.value(
        new groovy.json.JsonBuilder(params.deseq2.contrasts).toString()
    )

    DESEQ2(FEATURECOUNTS.out.counts, gtf_ch, groups_ch, contrasts_ch,
           common_r_ch, deseq2_r_ch, gene2symbol_r_ch)

    // ---- Layer 7: clusterProfiler ----
    contrast_names_ch = channel.value(
        params.deseq2.contrasts*.name.join(',')
    )

    sig_csvs_ch = DESEQ2.out.tables
        .flatten()
        .filter { csv -> csv.name =~ /_significant\.csv$/ }
        .collect()
        .ifEmpty { log.warn "[main] No significant DEGs found in any contrast — skipping clusterProfiler"; [] }

    CLUSTERPROFILER(sig_csvs_ch, contrast_names_ch, common_r_ch, clusterprofiler_r_ch)
}
