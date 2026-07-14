# CLAUDE.md

> Nextflow RNA-seq 上游分析流程 — AI 辅助开发上下文

## 项目概述

基于 Nextflow DSL2 的 RNA-seq 双端/单端上游分析流程，覆盖原始 FASTQ → 差异基因列表 + GO/KEGG 富集的完整分析链。

**工具链**: `FastQC → Trimmomatic → HISAT2 → featureCounts → MultiQC → DESeq2 → clusterProfiler`  
**输入**: FASTQ (PE `{sample}_R1/R2.fastq.gz` 或 SE `{sample}.fastq.gz`) + samplesheet CSV  
**输出**: 差异基因 CSV + PCA/MA/Volcano/Heatmap + GO/KEGG 富集分析 + 基因 Symbol 注释

## 目录结构

```
nextflow-rnaseq/
├── main.nf                        # 主工作流入口 (DSL2)
├── nextflow.config                # 全局配置 — 参数/容器/Profile/report
├── nextflow_schema.json           # JSON Schema 参数校验 (nf-validation 插件)
├── modules/
│   └── local/
│       ├── hisat2_index.nf        #   HISAT2 基因组索引构建
│       ├── hisat2_align_pe.nf     #   HISAT2 比对 (双端)
│       ├── hisat2_align_se.nf     #   HISAT2 比对 (单端)
│       ├── fastqc_pe.nf           #   FastQC 质控 (双端)
│       ├── fastqc_se.nf           #   FastQC 质控 (单端)
│       ├── trimmomatic_pe.nf      #   Trimmomatic 剪切 (双端)
│       ├── trimmomatic_se.nf      #   Trimmomatic 剪切 (单端)
│       ├── featurecounts.nf       #   featureCounts 基因定量
│       ├── multiqc.nf             #   MultiQC 汇总报告
│       ├── deseq2.nf              #   DESeq2 差异表达
│       └── clusterprofiler.nf     #   clusterProfiler GO/KEGG 富集
├── subworkflows/
│   └── local/
│       ├── qc_trim_align_pe.nf    #   子流程: QC→剪切→QC→比对 (PE)
│       └── qc_trim_align_se.nf    #   子流程: QC→剪切→QC→比对 (SE)
├── bin/
│   ├── common.R                   #   R 公共函数 (parse_arg + R_LIBS_ONLY)
│   ├── deseq2.R                   #   DESeq2 差异分析主脚本
│   ├── gene2symbol.R              #   基因 ID→Symbol 注释库
│   └── clusterprofiler.R          #   GO/KEGG 富集分析脚本
├── conf/
│   ├── base.config                #   基础资源配置 (CPU/内存/重试策略)
│   ├── docker.config              #   Docker 执行器
│   ├── apptainer.config           #   Apptainer 执行器
│   ├── test.config                #   E. coli 端到端完整测试
│   └── test_stub.config           #   stub-run 免容器 DAG 快速校验
├── containers/
│   ├── Dockerfile                 #   统合镜像 (全部工具)
│   ├── apptainer.def              #   Apptainer 定义
│   └── build.sh                   #   双引擎构建脚本
├── test/                          # E2E 测试套件 (E. coli)
│   ├── samplesheet.csv            #   4 样本分组表
│   ├── adapters.fa                #   空 adapter 文件
│   ├── ecoli.fa                   #   NCBI E. coli K-12 MG1655 参考基因组
│   ├── ecoli.gtf                  #   NCBI E. coli K-12 MG1655 注释
│   ├── ecoli_ko.fa                #   400 基因敲除版基因组 (生成模拟 reads)
│   └── hisat2_index/              #   预建 HISAT2 索引 (gitignore)
├── .gitignore
├── CLAUDE.md
├── README.md
└── README_CN.md
```

## DAG 拓扑

```
Layer 0: hisat2_index      一次性建索引 (用户提供预建索引则跳过)
Layer 1: fastqc_raw        每个 sample 一次
Layer 2: trimmomatic       每个 sample 一次，依赖 Layer 1
Layer 3: fastqc_trimmed    每个 sample 一次，依赖 Layer 2
Layer 3: hisat2_align      每个 sample 一次，依赖 Layer 2 + Layer 0
Layer 4: featurecounts     汇总所有 BAM，单次运行
Layer 5: multiqc           汇总所有 QC 报告，单次运行
Layer 6: deseq2            读取计数矩阵，单次运行，依赖 Layer 4
Layer 7: clusterprofiler   GO + KEGG 富集，单次运行，依赖 Layer 6
```

Layer 1-3 封装在 `qc_trim_align_{pe,se}` 子流程中，持有 `{sample}` 分组键，Nextflow 自动并行执行所有样本。

## 架构设计要点

### Process 统一模板

每个 process 遵循相同结构：`label` + `tag` + `publishDir` + `input:` + `output:` + `script:` + `stub:`。所有 `script:` 块首行加 `set -euo pipefail`。

### PE/SE 模式切换

`nextflow.config` 中 `params.read_mode: "paired"|"single"` 控制全局。`main.nf` 通过 `if/else` 分支选择子流程版本。SE 模式自动用 `null` 补齐 tuple 元素以统一 `multiMap` 索引。

### main.nf 架构

`main.nf` 的 workflow body 负责：
- 参数校验 + 对比 group 存在性验证
- 样本表单次解析 → `.multiMap()` 分流为 FASTQ / group / validation 三个通道
- HISAT2 索引：用户提供预建索引路径则直接用 `channel.value()` 传递，否则调用 `HISAT2_INDEX` 从 FASTA 构建
- DESeq2 分组 JSON：从 Groovy 文件读取构建
- CLUSTERPROFILER 输入：`.flatten().filter().collect().ifEmpty{...}` 保证空输入时优雅跳过

### R 脚本协作

```
deseq2.nf  --gtf {input.gtf} →
  deseq2.R:
    source("common.R")              # parse_arg + .libPaths()
    if (has_gtf) {
      source("gene2symbol.R")
      anno <- load_gtf_annotation(gtf)     # rtracklayer::import()
      res_df <- add_gene_symbol(res_df, anno)
      sig    <- add_gene_symbol(sig, anno)
    }

clusterprofiler.nf  --sig_files {params.sig_files_str} →
  clusterprofiler.R:
    source("common.R")
    bitr(gene_ids, fromType = from_type, toType = "ENTREZID", OrgDb = org_db)
    enrichGO(gene = genes_entrez, OrgDb = org_db, ont = ...)
    enrichKEGG(gene = genes_entrez, organism = kegg_org)
```

`gene2symbol.R` 兼容两种 GTF 格式：
- **Ensembl GTF**: 使用 `gene_name` 列
- **NCBI GTF**: 自动将 `gene` 列重命名为 `gene_name`
- `rtracklayer::import()` 参数 `features` 已废弃，改用 `gtf[gtf$type == "gene"]` 过滤

### R 脚本隔离

`common.R` 中 `R_LIBS_ONLY` 环境变量支持 conda/容器环境隔离。当 `R_LIBS_ONLY` 设置时，`.libPaths()` 只使用该路径，避免系统 R 库与 conda R 库冲突。

### R 脚本日志规范

三份 R 脚本统一使用 `[模块名]` 前缀消息格式（`message("[DESeq2] ...")` / `[clusterProfiler]` / `[gene2symbol]`），便于 `grep` 从 Nextflow 日志中过滤特定模块输出。

### `.multiMap()` 频道分流

Nextflow 26 中 `Channel.create()` 和 `.into{}` 闭包语法已禁用。使用 `.multiMap()` 将单次 CSV 解析结果广播到多个下游通道：

```groovy
sample_input_ch
    .multiMap { it ->
        fq:         tuple(it[0], it[2], it[3])   // → 比对流程
        group:      tuple(it[0], it[1])           // → DESeq2 分组
        validation: it[1]                          // → 对比校验
    }
    .set { forked }
```

### hisat2_index 预建索引复用

不支持 `cache 'deep'`（跨 `work/` 清理不持久）。改为判断 `--hisat2_index_prefix` 提供的 `.1.ht2` 文件是否存在：

- **存在** → `channel.value()` 直接传递给下游，跳过 `HISAT2_INDEX`
- **不存在** → 从 genome FASTA 自动构建

## 配置约定

[nextflow.config](nextflow.config) 是唯一配置入口。关键段：

| 段 | 内容 |
|---|---|
| `params.input` | 样本表 CSV 路径（必需） |
| `params.read_mode` | `"paired"` / `"single"` |
| `params.outdir` | 输出根目录 (默认 `./results`) |
| `params.genome_fasta` | 参考基因组 FASTA（必需） |
| `params.annotation_gtf` | 基因注释 GTF（必需） |
| `params.hisat2_index_prefix` | HISAT2 索引 basename（可不提供，自动构建） |
| `params.adapter_fasta` | Trimmomatic adapter FASTA（可不提供） |
| `params.fastqc` | extra |
| `params.trimmomatic` | phred, illuminaclip, leading, trailing, slidingwindow, minlen |
| `params.hisat2` | extra |
| `params.featurecounts` | feature_type, attr_type, strandedness, extra |
| `params.deseq2` | padj_threshold, log2fc_threshold, contrasts, top_n_genes |
| `params.clusterprofiler` | org_db, kegg_organism, from_type, pvalue_cutoff, qvalue_cutoff, show_category, gene_id_col |

### 参数校验

[nextflow_schema.json](nextflow_schema.json) 通过 `nf-validation` 插件自动校验参数类型和格式：
- `genome_fasta` 只接受 `.fa/.fasta/.fna(.gz)` 后缀
- `annotation_gtf` 只接受 `.gtf(.gz)` 后缀
- `strandedness` 只允许 0/1/2
- `contrasts` 数组中每项必须有 `name`/`case`/`control`

插件在运行 `-profile docker/apptainer/test` 时自动加载。

### 新增样本

仅修改 samplesheet CSV：

```csv
sample,group,fastq_1,fastq_2
WT_1,WT,/path/to/WT_1_R1.fastq.gz,/path/to/WT_1_R2.fastq.gz
New_1,New,/path/to/New_1_R1.fastq.gz,/path/to/New_1_R2.fastq.gz   # ← 添加
```

SE 模式省略 `fastq_2` 列。对比校验会自动检查 `case`/`control` 是否在样本表中存在，拼写错误会在启动时立即报错。

### 新增对比

```groovy
// nextflow.config 或命令行 -c
params.deseq2.contrasts = [
    [name: "New_vs_Ctrl", 'case': "New", control: "Ctrl"]   // ← 添加
]
```

### 切换 PE/SE

```bash
nextflow run main.nf --read_mode "single" ...
```

无需修改任何 module 或 subworkflow。

## 运行命令

```bash
# 从 nextflow-rnaseq/ 目录执行

# 语法校验 (免容器)
nextflow run main.nf -profile test_stub -stub-run \
    --input test/samplesheet.csv \
    --genome_fasta test/ecoli.fa \
    --annotation_gtf test/ecoli.gtf

# Docker 模式 (生产)
bash containers/build.sh docker
nextflow run main.nf -profile docker \
    --input samplesheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf genes.gtf

# Apptainer 模式 (HPC)
bash containers/build.sh all
nextflow run main.nf -profile apptainer \
    --input samplesheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf genes.gtf

# 使用预建 HISAT2 索引 (跳过索引构建)
hisat2-build genome.fa /data/ref/hisat2/genome
nextflow run main.nf -profile docker \
    --genome_fasta genome.fa \
    --annotation_gtf genes.gtf \
    --hisat2_index_prefix /data/ref/hisat2/genome ...

# 中断恢复
nextflow run main.nf -profile docker ... -resume

# 查看参数帮助 (需要 nf-validation 插件)
nextflow run main.nf --help
```

## Profile

| Profile | 用途 | 关键配置 |
|---------|------|---------|
| `docker` | 生产 — Docker 容器 | `process.container = 'rnaseq-pipeline:1.0.0'` |
| `apptainer` | HPC — Apptainer 容器 | `apptainer.enabled = true` |
| `test` | 测试 — E. coli 端到端 | mamba 本地执行, NCBI GTF, gene-level 计数 |
| `test_stub` | 开发 — stub-run 免容器 | 极低资源, `process.container = ''` |

## 错误处理

| 退出码 | 策略 | 说明 |
|--------|------|------|
| 137 | retry (2次) | OOM kill — 重试可能分配到有更多内存的节点 |
| 143 | retry (2次) | SIGTERM — HPC 调度器抢占，换节点重试 |
| 其他 | terminate | 工具 bug / segfault — 重试无用，直接终止 |

配置在 [conf/base.config](conf/base.config) 的 `errorStrategy` 闭包中。

## 工具版本

| 工具 | 版本 | 安装方式 |
|------|------|---------|
| FastQC | 0.12.1 | Docker / conda |
| Trimmomatic | 0.39 | Docker / conda |
| HISAT2 | 2.2.1 | Docker / conda |
| Samtools | 1.18 | Docker / conda |
| featureCounts | 2.0.6 | Docker / conda |
| MultiQC | 1.21 | Docker / conda |
| R | 4.3.2 | Docker / conda |
| DESeq2 | 1.42.0 | Docker / conda |
| rtracklayer | 1.62.0 | Docker / conda |
| clusterProfiler | 4.10.0 | Docker / conda |
| enrichplot | 1.22.0 | Docker / conda |

## 容器化

`containers/build.sh` 是容器操作的统一入口：

```bash
bash containers/build.sh docker         # 构建 Docker 镜像
bash containers/build.sh all            # Docker + 转 SIF (推荐)
bash containers/build.sh test           # 验证镜像内所有工具
bash containers/build.sh clean          # 清理
```

镜像基于 `continuumio/miniconda3:24.1.2-0`，所有工具通过 conda 安装。

## Nextflow 注意事项

1. **Shell 健壮性**: 所有 process `script:` 块首行必须加 `set -euo pipefail`，确保管道错误不会静默吞噬失败
2. **Shell 中 JSON 参数**: 传给 R 脚本的 JSON 用单引号包裹 `'${groups_json}'`，防止 shell 解析 JSON 双引号
3. **`.multiMap()` 替代 `.into{}`**: Nextflow 26 禁用了 `Channel.create()` 和 `.into{}` 闭包语法。使用 `.multiMap()` 将单通道广播到多个分支
4. **`.flatten()` 修复 List emit**: `path("*.csv")` 在 output 中匹配多个文件时以单个 List 形式 emit。下游 `.filter{}` 前需 `.flatten()` 拆分为独立 item
5. **Stub 用 Groovy 迭代**: process `stub:` 块中避免 bash 解析 JSON，改用 Groovy `${params.xxx.collect{...}.join()}` 直接生成文件名
6. **`params` 不在 process 内使用**: process 定义中不引用 `$projectDir` 或 `params.hisat2_index_prefix`。脚本文件通过 `path` input 传入
7. **`optional: true` 处理空输出**: `CLUSTERPROFILER` 的 `path("*.csv")` 和 `path("*.pdf")` 设为 `optional: true`，因为无显著基因时不产生文件
8. **hisat2 | samtools 线程分配**: hisat2 使用 `${task.cpus}`，samtools 固定 2 线程。避免管道两端都占满导致 `2× cpus` 超配
9. **对比校验在启动时完成**: `main.nf` 从 `multiMap` 分支提取 unique groups 并与 `params.deseq2.contrasts` 比对。拼写错误立即 `exit 1` 而非等到 DESeq2 R 运行时
10. **`R_LIBS_ONLY` 隔离 conda R 库**: `common.R` 检查该环境变量，若设置则用 `.libPaths()` 替换整个搜索路径，避免系统 R 库与 conda R 库的二进制不兼容导致 segfault
11. **`test/` 下的模拟 reads 不提交**: `.gitignore` 中 `test/*.fastq.gz` 和 `test/hisat2_index/` 不纳入版本控制
12. **featureCounts `extra` 参数空字符串**: 当 `extra` 为空时，`def extra = params.featurecounts.extra ?: ''` 配合 shell 续行符处理空行，避免空参数导致命令行解析错误
13. **GTF 兼容 NCBI/Ensembl**: `gene2symbol.R` 自动检测 `gene_name` vs `gene` 列名，NCBI GTF 用 `gene` 列时自动重命名
14. **clusterProfiler `bitr()` 容错**: `bitr()` 用 `tryCatch` 包裹，keytype 不匹配时输出 warning 并优雅跳过，不中断整个流程

## 测试

### stub-run (语法/DAG 校验)

```bash
nextflow run main.nf -profile test_stub -stub-run \
    --input test/samplesheet.csv \
    --genome_fasta test/ecoli.fa \
    --annotation_gtf test/ecoli.gtf
# 预期: SUCCESS 21/21 (全部 process 执行 stub)
```

### E2E 端到端 (E. coli 真实分析)

```bash
# 1. 激活 conda/mamba 环境（保证命令行工具在 PATH 中）
conda activate nextflow-rnaseq

# 2. 设置 R 包隔离路径（触发 common.R 中的 .libPaths() 替换逻辑）
export R_LIBS_ONLY="$(conda info --base)/envs/nextflow-rnaseq/lib/R/library"

# 3. 运行
nextflow run main.nf -profile test \
    --input test/samplesheet.csv \
    --genome_fasta test/ecoli.fa \
    --annotation_gtf test/ecoli.gtf \
    --adapter_fasta test/adapters.fa \
    --hisat2_index_prefix test/hisat2_index/ecoli \
    --outdir test_results
# 预期: SUCCESS 20/20 (跳过索引构建)
#       344 DEGs, GO/KEGG 富集结果
```

E2E 测试用 E. coli K-12 MG1655 (NCBI GCF_000005845.2)：
- WT 样本 reads 从全基因组模拟，KO 样本从 400 基因敲除版基因组模拟
- featureCounts 用 `-t gene -g gene` 在基因层计数
- clusterProfiler 用 `org.EcK12.eg.db` + `from_type = "SYMBOL"`
- 由于 KO 基因是随机选择，GO 富集结果可能稀疏 — 属正常，真实数据会有实质通路聚集

### 数据生成

模拟 reads 用 `wgsim`（samtools 附带）生成：

```bash
# WT: 全基因组 100K PE reads
wgsim -N 100000 -1 100 -2 100 -d 300 -s 20 -e 0.01 \
    -S 42 ecoli.fa WT_1_R1.fastq WT_1_R2.fastq

# KO: 敲除版基因组 100K PE reads
wgsim -N 100000 -1 100 -2 100 -d 300 -s 20 -e 0.01 \
    -S 177 ecoli_ko.fa KO_1_R1.fastq KO_1_R2.fastq
```

## 不属于本项目的文件

- `../01.RNA_seq/` — Snakemake 版 RNA-seq 流程
- `test_results/` — 测试运行时产出 (gitignore)
- `work/` — Nextflow 工作目录 (gitignore)
