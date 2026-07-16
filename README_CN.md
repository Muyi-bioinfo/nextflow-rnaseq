# RNA-seq 分析流程 (Nextflow)

[English](README.md) | [中文](README_CN.md)

[![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-0dc09d)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Docker-支持-blue)](https://www.docker.com/)
[![Apptainer](https://img.shields.io/badge/Apptainer-支持-orange)](https://apptainer.org/)

基于 Nextflow DSL2 构建的 RNA-seq 上游分析流程，覆盖原始 FASTQ 质控至差异表达基因列表及功能富集的完整分析链，支持容器化部署。

**分析工具链**：FastQC → Trimmomatic → HISAT2 → featureCounts → MultiQC → DESeq2 → clusterProfiler（含基因 Symbol 注释 + GO/KEGG 富集分析）

## 特性

- **PE/SE 自动适配** — 一个 `--read_mode` 开关，全部 module 和 subworkflow 自动跟随
- **DSL2 模块化设计** — 11 个 process module + 2 个 subworkflow，统一的 `label`/`tag`/`publishDir`/`stub` 模式
- **编号输出目录** — `01_fastqc/` → … → `07_clusterprofiler/`，`ls` 即见执行顺序
- **多 Profile 运行** — `docker`、`apptainer`、`test`（E2E 测试）、`test_stub`（DAG 空跑校验）
- **预建索引复用** — 提供 `--hisat2_index_prefix` 指向已有索引路径，自动跳过索引构建
- **基因 Symbol 注释** — DESeq2 输出自动从同一 GTF 追加 gene symbol 列
- **GO/KEGG 功能富集** — clusterProfiler 对差异基因进行 GO + KEGG 通路富集分析
- **对比校验** — 启动时检查 `case`/`control` 是否存在于样本表中，拼写错误立即报错
- **Schema 校验** — `nextflow_schema.json` + `nf-validation` 插件自动校验参数类型和格式
- **样本表驱动** — 单个 CSV 定义样本、分组和 FASTQ 路径
- **可复现** — 容器镜像锁定全部工具版本
- **易于扩展** — 新增样本或对比仅修改外部配置，零 workflow 代码改动

📖 **[Nextflow 学习笔记](LEARNING_NEXTFLOW.md)** — 基于本项目代码系统学习 Nextflow 核心概念

## 快速开始

```bash
# 1. 进入项目
cd nextflow-rnaseq

# 2. 准备样本表 CSV
cat > samplesheet.csv << EOF
sample,group,fastq_1,fastq_2
WT_1,WT,/data/WT_1_R1.fastq.gz,/data/WT_1_R2.fastq.gz
WT_2,WT,/data/WT_2_R1.fastq.gz,/data/WT_2_R2.fastq.gz
KO_1,KO,/data/KO_1_R1.fastq.gz,/data/KO_1_R2.fastq.gz
KO_2,KO,/data/KO_2_R1.fastq.gz,/data/KO_2_R2.fastq.gz
EOF

# 3. 空跑校验（stub — 不调用真实工具）
nextflow run main.nf -profile test_stub -stub-run \
    --input samplesheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf genes.gtf

# 4. 正式运行
nextflow run main.nf -profile docker \
    --input samplesheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf genes.gtf
```

## 环境依赖

| 依赖 | 说明 |
|---|---|
| [Nextflow](https://www.nextflow.io/) ≥24.0 | 工作流引擎 |
| [Docker](https://www.docker.com/) ≥20.04 | Docker 容器模式 |
| [Apptainer](https://apptainer.org/) / [Singularity](https://sylabs.io/) ≥1.0 | HPC 容器模式 |
| [Conda](https://docs.conda.io/) / [Mamba](https://mamba.readthedocs.io/) | 本地测试模式（可选） |

> 先构建容器镜像：`bash containers/build.sh docker`

## 项目结构

```
nextflow-rnaseq/
├── main.nf                        # 主工作流 (DSL2)
├── nextflow.config                # 全局配置 — 参数 / 容器 / Profile
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
│   ├── deseq2.R                   #   DESeq2 差异分析 + 诊断图
│   ├── gene2symbol.R              #   基因 ID → Symbol 注释库 (rtracklayer)
│   └── clusterprofiler.R          #   GO/KEGG 富集分析脚本
├── conf/
│   ├── base.config                #   基础资源配置 + 错误重试策略
│   ├── docker.config              #   Docker 执行器
│   ├── apptainer.config           #   Apptainer 执行器
│   ├── test.config                #   E. coli 端到端测试
│   └── test_stub.config           #   stub-run DAG 验证
├── containers/
│   ├── Dockerfile                 #   统合容器镜像
│   ├── apptainer.def              #   Apptainer 定义
│   └── build.sh                   #   双引擎构建脚本
├── test/                          # E2E 测试套件 (E. coli)
│   ├── samplesheet.csv            #   4 样本测试表
│   ├── adapters.fa                #   空 adapter 文件
│   ├── ecoli.fa                   #   NCBI E. coli K-12 MG1655 基因组
│   ├── ecoli.gtf                  #   NCBI E. coli K-12 MG1655 注释
│   ├── ecoli_ko.fa                #   400 基因敲除版基因组 (生成模拟 reads)
│   └── hisat2_index/              #   预建 HISAT2 索引 (gitignore)
├── .gitignore
├── README.md
└── README_CN.md
```

## 流程 DAG

```
Layer 0: hisat2_index      一次性建基因组索引（提供预建索引则跳过）
Layer 1: fastqc_raw        每个样本独立执行，自动并行
Layer 2: trimmomatic        每个样本独立执行，PE/SE 自动适配
Layer 3: fastqc_trimmed    每个样本独立执行
Layer 3: hisat2_align      每个样本独立执行，PE/SE 自动适配
Layer 4: featurecounts     汇总所有 BAM，单次运行
Layer 5: multiqc           汇总 QC 报告，单次运行
Layer 6: deseq2            差异表达分析 + 基因 Symbol 注释
Layer 7: clusterprofiler   GO + KEGG 功能富集，依赖 Layer 6
```

Layer 1–3 封装在 `qc_trim_align_{pe,se}` 子流程中。持有样本分组键的 process 由 Nextflow 自动并行执行。

## 配置说明

[`nextflow.config`](nextflow.config) 是所有可调参数的统一入口。所有参数可通过命令行 `--param value` 覆盖或在配置文件中修改。

### 样本表

CSV 格式（逗号分隔；`.tsv` 和 `.txt` 也支持）：

```csv
sample,group,fastq_1,fastq_2
WT_1,WT,/data/WT_1_R1.fastq.gz,/data/WT_1_R2.fastq.gz
WT_2,WT,/data/WT_2_R1.fastq.gz,/data/WT_2_R2.fastq.gz
KO_1,KO,/data/KO_1_R1.fastq.gz,/data/KO_1_R2.fastq.gz
KO_2,KO,/data/KO_2_R1.fastq.gz,/data/KO_2_R2.fastq.gz
```

单端数据：省略 `fastq_2` 列。

### 双端 → 单端切换

```bash
nextflow run main.nf --read_mode "single" ...
```

无需修改任何 module 或 subworkflow。

### 差异比较设计

在 `nextflow.config` 中设置或通过 `-c` 传入：

```groovy
params.deseq2.contrasts = [
    [name: "KO_vs_WT", 'case': "KO", control: "WT"]
]
```

校验机制：启动时自动检查 `case` 和 `control` 分组是否存在于样本表中，拼写错误会立即给出明确错误信息。

### clusterProfiler 富集分析

```groovy
params.clusterprofiler {
    org_db        = "org.Hs.eg.db"    // 物种注释库（小鼠用 org.Mm.eg.db）
    kegg_organism = "hsa"             // KEGG 物种代码（mmu, rno, ...）
    from_type     = "ENSEMBL"         // DEG 表中的基因 ID 类型 → ENTREZID
    pvalue_cutoff = 0.05
    qvalue_cutoff = 0.2
    show_category = 15                // 气泡图展示的 top N 条目
    gene_id_col   = "gene_id"
}
```

### 预建 HISAT2 索引（跳过索引构建）

```bash
hisat2-build genome.fa /ref/hisat2/genome
nextflow run main.nf \
    --hisat2_index_prefix /ref/hisat2/genome ...
# → HISAT2_INDEX 步骤直接跳过
```

## 运行方式

### Stub 模式（DAG 校验）

```bash
nextflow run main.nf -profile test_stub -stub-run \
    --input test/samplesheet.csv \
    --genome_fasta test/ecoli.fa \
    --annotation_gtf test/ecoli.gtf
```

### E2E 端到端测试（E. coli, 本地 conda）

测试数据：E. coli K-12 MG1655 的 4 个模拟双端样本（WT vs KO，各 2 个
生物学重复，每样本 10 万对 reads）。基因组索引已预建，自动跳过构建步骤。

```bash
# 前置条件：激活安装了全部工具的 conda 环境，并设置 R_LIBS_ONLY
# 以隔离 R 包搜索路径，避免系统 R 库与 conda R 库冲突。
conda activate nextflow-rnaseq
export R_LIBS_ONLY="$(conda info --base)/envs/nextflow-rnaseq/lib/R/library"

nextflow run main.nf -profile test \
    --input test/samplesheet.csv \
    --genome_fasta test/ecoli.fa \
    --annotation_gtf test/ecoli.gtf \
    --adapter_fasta test/adapters.fa \
    --hisat2_index_prefix test/hisat2_index/ecoli \
    --outdir test_results
```

> `R_LIBS_ONLY` 将 R 的库搜索路径限定在 conda 环境内，防止系统 R 包与
> conda 安装的 Bioconductor 包之间二进制不兼容导致段错误。

### Docker 模式（生产环境）

```bash
bash containers/build.sh docker
nextflow run main.nf -profile docker \
    --input samplesheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf genes.gtf
```

### Apptainer 模式（HPC 集群）

```bash
bash containers/build.sh all
nextflow run main.nf -profile apptainer \
    --input samplesheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf genes.gtf
```

### 中断恢复

```bash
nextflow run main.nf -profile docker ... -resume
```

## Profile

| Profile | 用途 | 关键配置 |
|---------|------|---------|
| `docker` | 生产 — Docker 容器 | `process.container = 'rnaseq-pipeline:1.0.0'` |
| `apptainer` | HPC — Apptainer 容器 | `apptainer.enabled = true` |
| `test` | 测试 — E. coli 端到端 | 本地 conda, NCBI GTF, gene-level 计数 |
| `test_stub` | 开发 — stub-run DAG 校验 | 极低资源, 无容器 |

## 输出文件说明

### 各步骤输出

| 目录 | 关键文件 |
|---|---|
| `01_fastqc/` | `{read}_fastqc.{html,zip}`（原始 + 剪切后） |
| `02_trimmomatic/` | `{sample}_R{1,2}.trimmed.fastq.gz` |
| `03_hisat2/` | `{sample}.sorted.bam` + `.bam.bai` + `.hisat2.log` |
| `04_featurecounts/` | `featurecounts.txt`、`featurecounts.summary.txt` |
| `05_multiqc/` | `multiqc_report.html` |

### DESeq2 输出（`06_deseq2/`）

| 文件 | 内容 |
|---|---|
| `{contrast}_all_results.csv` | 全部基因统计表，含 gene_name 列 |
| `{contrast}_significant.csv` | 显著差异基因列表，含 gene_name 列 |
| `{contrast}_MA_plot.pdf` | MA 图（平均表达量 vs. log2FC） |
| `{contrast}_volcano_plot.pdf` | 火山图（有 GTF 时标注 gene symbol） |
| `PCA_plot.pdf` / `.png` | 样本 PCA 降维图 |
| `sample_distance_heatmap.pdf` | 样本间表达距离矩阵 |
| `DEG_heatmap.pdf` | Top N 差异基因表达热图 |

### clusterProfiler 输出（`07_clusterprofiler/`）

| 文件 | 内容 |
|---|---|
| `{contrast}_GO_enrichment.csv` | GO 富集分析结果（BP/CC/MF） |
| `{contrast}_GO_dotplot.pdf` | GO 富集气泡图 |
| `{contrast}_KEGG_enrichment.csv` | KEGG 通路富集结果 |
| `{contrast}_KEGG_dotplot.pdf` | KEGG 通路气泡图 |

## 基因 Symbol 注释

流程自动使用与比对和计数同一份 GTF 文件，为 DEG 输出追加基因 Symbol 列，确保上游到下游注释版本一致。

```bash
# 也可独立使用
Rscript bin/gene2symbol.R --input DEG.csv --gtf genes.gtf --output DEG_anno.csv
```

实现方式：`rtracklayer::import()` 标准解析 GTF（非正则）。同时兼容 NCBI GTF（`gene` 列）和 Ensembl GTF（`gene_name` 列）。

## 容器化支持

```bash
bash containers/build.sh docker      # 本地 Docker 镜像
bash containers/build.sh all         # Docker → SIF（HPC 推荐）
bash containers/build.sh test        # 验证镜像内全部工具
bash containers/build.sh clean       # 清理
```

统合镜像包含：FastQC 0.12.1、Trimmomatic 0.39、HISAT2 2.2.1、Samtools 1.18、featureCounts 2.0.6、MultiQC 1.21、R 4.3.2 + DESeq2 1.42.0 + rtracklayer 1.62.0 + clusterProfiler 4.10.0 + enrichplot 1.22.0。

## Nextflow 特性展示

| 特性 | 应用 | 说明 |
|---|---|---|
| **DSL2 模块** | 11 process + 2 subworkflow | `include { } from '...'`，结构清晰 |
| **`.multiMap()`** | 通道分流 | 单次 CSV 解析广播到 3 个下游消费者 |
| **`stub:` 块** | 空跑验证 | 所有 process 含 stub，`-stub-run` 校验 DAG |
| **Profile 体系** | 4 个 profile | docker / apptainer / test / test_stub |
| **Schema 校验** | `nextflow_schema.json` | nf-validation 插件自动参数检查 |
| **错误重试策略** | 按退出码区分 | OOM(137) 和抢占(143) → 重试；其他 → 终止 |
| **Groovy 配置逻辑** | 条件判断 | 预建索引检测、对比校验 |
| **`publishDir`** | 逐 process 输出 | 统一的编号目录布局 |
| **`cache 'deep'`** | 内容寻址缓存 | 索引跨批次复用 |
| **PE/SE 批处理** | `--read_mode` 单开关 | main.nf 中 `if/else` 分支自动切换 |
| **R 脚本隔离** | `R_LIBS_ONLY` 环境变量 | conda/容器环境中隔离 R 库路径 |
| **可选输出** | `optional: true` | clusterProfiler 无显著基因时优雅跳过 |

## 常见问题

**如何新增样本？** — 在 samplesheet CSV 中追加一行即可。

**如何用单端数据？** — `nextflow run main.nf --read_mode "single" ...`，无需改代码。

**如何跳过索引构建？** — 一次性建索引：`hisat2-build genome.fa /path/genome`，运行时指定 `--hisat2_index_prefix /path/genome`。

**如何彻底重跑？** — `rm -rf work/ test_results/ && nextflow run main.nf ...`

**如何中断恢复？** — `nextflow run main.nf -profile docker ... -resume`

**集群用 SGE 不是 Slurm？** — Nextflow 原生支持 SGE、LSF、PBS/Torque 和 Slurm，修改 executor 配置即可。

**如何修改集群资源申请？** — 编辑 `conf/base.config` 中的 process label 配置。

**如何为其他物种做富集分析？** — 设置 `--clusterprofiler.org_db` 和 `--clusterprofiler.kegg_organism`。Bioconductor 提供 100+ 物种 OrgDb 包。

## 工具版本清单

| 工具 | 版本 |
|---|---|
| FastQC | 0.12.1 |
| Trimmomatic | 0.39 |
| HISAT2 | 2.2.1 |
| Samtools | 1.18 |
| featureCounts | 2.0.6 |
| MultiQC | 1.21 |
| R | 4.3.2 |
| DESeq2 | 1.42.0 |
| rtracklayer | 1.62.0 |
| clusterProfiler | 4.10.0 |
| enrichplot | 1.22.0 |

## 许可证

本项目仅供教育和研究目的使用。
