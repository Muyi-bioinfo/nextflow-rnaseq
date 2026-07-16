# Nextflow 学习笔记 — 基于 nextflow-rnaseq 项目

> 本文档以 `nextflow-rnaseq` 项目的实际代码为例，按"运行 → 概念 → 实战"的路径，系统学习 Nextflow 核心知识。
>
> 建议阅读顺序：Part 一（跑起来）→ Part 二（理解概念）→ Part 三（实战串联）。

## 前置约定

- 终端命令以 `$` 开头，Groovy/Nextflow 代码使用 ` ```groovy ` 标记
- 示例代码来自项目实际文件，代码块上方标注了可点击的文件链接
- 本文档假设 samplesheet 内容如下（4 个样本，双端测序）：

```csv
sample,group,fastq_1,fastq_2
WT_1,WT,data/WT_1_R1.fastq.gz,data/WT_1_R2.fastq.gz
WT_2,WT,data/WT_2_R1.fastq.gz,data/WT_2_R2.fastq.gz
KO_1,KO,data/KO_1_R1.fastq.gz,data/KO_1_R2.fastq.gz
KO_2,KO,data/KO_2_R1.fastq.gz,data/KO_2_R2.fastq.gz
```

---

# 第一部分：快速开始

## 1. 安装与环境

### 安装 Nextflow

Nextflow 是一个独立的命令行工具，无需 root 权限，只需要 Java 运行环境。

```bash
# 方式一：官方推荐
curl -s https://get.nextflow.io | bash

# 方式二：conda / mamba
conda install -c bioconda nextflow

# 验证安装
./nextflow info
```

安装后把 `nextflow` 可执行文件移到 `PATH` 中的目录（如 `~/bin/`）即可全局使用。

### Java 环境

Nextflow 基于 Java 运行，需要 **Java 11 或更高版本**（推荐 Java 17/21）。

```bash
# 检查 Java 版本
java -version

# conda 安装
conda install -c conda-forge openjdk=17
```

> Nextflow 本身不需要 `JAVA_HOME` 环境变量，但如果系统中多个 Java 版本共存，可以通过 `JAVA_HOME` 指定使用的版本。

---

## 2. 基础命令

### `nextflow run` — 执行工作流

```bash
# 从本地文件运行
nextflow run main.nf

# 传递参数（双横线 --）
nextflow run main.nf --input samplesheet.csv --outdir ./results

# 指定 profile（单横线 -）
nextflow run main.nf -profile docker

# 并行进程数
nextflow run main.nf -profile docker --max_cpus 8
```

**关键区别**：
- `--param value`（双横线）：工作流**参数**，传给流程内部（对应 `params.xxx`）
- `-flag`（单横线）：Nextflow **运行选项**，控制执行器行为（如 `-resume`、`-profile`）

### `-resume` — 中断恢复 / 缓存复用

这是 Nextflow 最强大的特性之一。

```bash
# 首次运行
nextflow run main.nf -profile docker --input samplesheet.csv

# 如果运行中断（如 OOM、断电），加 -resume 恢复
nextflow run main.nf -profile docker --input samplesheet.csv -resume
```

**工作原理**：每个 process 的 input + script + 容器/conda → 生成唯一哈希。哈希匹配的 process 直接跳过（从 `work/` 取缓存），只重跑失败或内容变化的 process 及其下游。

### `-profile` — 执行环境切换

`-profile` 用于选择一组预定义的配置组合，定义在 [nextflow.config](nextflow.config) 中。

```bash
# 单一 profile
nextflow run main.nf -profile docker

# 同时使用多个 profile
nextflow run main.nf -profile docker,test
```

本项目中的 profile（详见 [§6.3](#63-profile--执行环境切换)）：

| Profile | 用途 | 行为 |
|---------|------|------|
| `docker` | 生产 — Docker 容器 | 启用 Docker，镜像内包含全部工具 |
| `apptainer` | HPC — Apptainer | 启用 Apptainer/Singularity |
| `test` | E. coli 端到端测试 | mamba 本地执行，真实分析 |
| `test_stub` | 开发 — 语法/DAG 校验 | 免容器，极低资源，配合 `-stub-run` |

---

## 3. 缓存机制深入

### 缓存哈希如何生成

→ [hisat2_index.nf L10](modules/local/hisat2_index.nf#L10)

```
缓存 key = hash(input 内容 + script 代码 + 容器/conda 环境 + cpus/memory 等 directives)


```
每个 task 在 `work/` 下的目录名就是这个哈希值。Nextflow 启动 task 前先查哈希是否命中，命中则跳过。
### 什么操作会破坏缓存
| 操作 | 缓存？ |
|------|:---:|
| 修改 process 的 `script:` 块 | ❌ 全部失效 |
| 修改 input 文件**内容**（文件名相同但内容变了） | ❌ 全部失效 |
| 更换容器镜像 | ❌ 全部失效 |
| 修改 `cpus` / `memory` 等 directives | ❌ 全部失效 |
| 改 `tag` 或 `publishDir` | ✅ 不影响 |
| 修改 params 中**未参与 process 计算**的值 | ✅ 不影响 |
| 增加新样本（旧样本的 task 仍命中） | ✅ 旧样本缓存保留 |
### `cache 'deep'`
本项目 使用了 `cache 'deep'`：
```groovy
process HISAT2_INDEX {
    cache 'deep'
    // ...
}
```

普通缓存仅在当前运行中保留。`cache 'deep'` 让结果**跨运行持久化**——即使 `work/` 被清理，索引也不会丢失。适合计算代价极高且结果很少变化的 task（如基因组索引）。

### 开发习惯

```
改完代码 → stub-run 验 DAG（几秒）→ 小数据测试 → 真实数据
```

不要在真实数据上边改边跑——每次改 `script:` 都会让缓存全部失效。

---

## 4. 运行产物

### 4.1 `work/` — Nextflow 内部执行目录

`work/` 是 Nextflow 的核心执行目录，**每个 process 的每次调用都会在其中创建一个独立子目录**。

```
work/
├── 3f/                          # 哈希前2位
│   └── a1b2c3d4e5f6.../         # 完整哈希 = task 唯一标识
│       ├── .command.sh           # 实际执行的 shell 脚本
│       ├── .command.run          # Nextflow task wrapper
│       ├── .command.log          # stdout（.command.sh 的输出）
│       ├── .command.err          # stderr
│       ├── .command.trace        # 资源使用记录
│       ├── .exitcode             # 退出码（0 = 成功）
│       ├── WT_1.sorted.bam       # task 产出文件
│       ├── WT_1.sorted.bam.bai
│       └── WT_1_R1.fastq.gz -> /original/path/  # 输入文件的符号链接
```

### 4.2 `outdir` — 最终结果目录

通过 `publishDir` 指令，从 `work/` 中把需要的文件复制到用户可见的结果目录。

### 4.3 `work/` 与 `outdir` 的关系

| 维度 | `work/` | `outdir` (`results/`) |
|------|---------|----------------------|
| **定位** | Nextflow 内部执行缓存 | 面向用户的结果目录 |
| **内容** | 所有中间产物 + 执行元数据 | 只包含 `publishDir` 发布的文件 |
| **结构** | 哈希目录，不可读 | 按分析步骤分目录，可读 |
| **可否删除** | 成功后可删（`nextflow clean -f`） | 永久保留 |
| **`-resume` 依赖** | 依赖此目录 | 不依赖 |

### 4.4 查看日志

#### 实时运行信息

运行时终端输出就是实时日志：

```
[1c/3e8a1b] process > HISAT2_ALIGN_PE (WT_1)   ← [已完成/总数] process名 (tag)
```

#### 单独 task 的执行日志

```bash
# 查找失败 task 的日志
find work -name ".command.log" -exec grep -l "ERROR" {} \;

# 查看某个 task 的完整输出
cat work/3f/a1b2c3d4/.command.log
cat work/3f/a1b2c3d4/.command.err
```

#### `.command.run` 中的 tag 信息

`.command.run` 第 3 行的 `### name:` 字段包含了 tag：

```
### name: 'QC_TRIM_ALIGN_PE:FASTQC_TRIMMED_PE (KO_1)'
         └─ 子流程:进程名 ──────┘          └─ tag 值 ─┘
```

```bash
# 查找所有 task 的 process 名和 tag
find work/ -name ".command.run" -exec grep "### name:" {} \; | sort
```

### 4.5 Report / Trace / Timeline — 运行报告

配置在 [nextflow.config L109-126](nextflow.config#L109)：
```groovy
dag {
    enabled   = true
    file      = "${params.outdir}/pipeline_dag.svg"
}
timeline {
    enabled   = true
    file      = "${params.outdir}/pipeline_timeline.html"
}
report {
    enabled   = true
    file      = "${params.outdir}/pipeline_report.html"
}
```

每次运行自动产出 3 个文件到 outdir：

| 文件 | 内容 | 回答的问题 |
|------|------|-----------|
| `pipeline_dag.svg` | 工作流 DAG 图 | 流程结构对不对？ |
| `pipeline_timeline.html` | 每个 task 的时间线甘特图 | 哪个步骤最慢？并行效率如何？ |
| `pipeline_report.html` | 综合运行摘要 | 这次运行整体怎么样？ |
| `pipeline_trace.txt` | 每个 task 的资源记录（TSV 格式） | 某个 task 用了多少 CPU？内存峰值？ |

**Trace 关键字段**：

| 字段 | 含义 | 实际用途 |
|------|------|---------|
| `name` | process 名 + tag | 定位是哪个样本的哪个步骤 |
| `status` | COMPLETED / FAILED / CACHED | 哪些重跑了，哪些走了缓存 |
| `%cpu` | CPU 使用率 | `780%` ≈ 8 核用满，`100%` ≈ 单线程瓶颈 |
| `peak_rss` | 内存峰值 | 接近 `memory` 限制就危险 |
| `exit` | 退出码 | 0=成功，137=OOM，143=SIGTERM |

```bash
# 找出耗时最长的 5 个 task
sort -t$'\t' -k9 -nr results/pipeline_trace.txt | head -5

# 找出所有失败的 task
grep "FAILED" results/pipeline_trace.txt
```

---

## 5. 开发工作流

推荐的操作闭环：

```bash
# 1. stub-run 校验（语法 + DAG，几秒完成）
nextflow run main.nf -profile test_stub -stub-run \
    --input test/samplesheet.csv \
    --genome_fasta test/ecoli.fa \
    --annotation_gtf test/ecoli.gtf

# 2. 小数据测试（E. coli 端到端，几分钟完成）
conda activate nextflow-rnaseq
nextflow run main.nf -profile test \
    --input test/samplesheet.csv \
    --genome_fasta test/ecoli.fa \
    --annotation_gtf test/ecoli.gtf \
    --adapter_fasta test/adapters.fa \
    --hisat2_index_prefix test/hisat2_index/ecoli \
    --outdir test_results

# 3. 真实数据运行（生产）
nextflow run main.nf -profile docker \
    --input samplesheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf genes.gtf

# 4. 查看报告（定位瓶颈和问题）
# 打开 results/pipeline_report.html
# 打开 results/pipeline_timeline.html
```

---

# 第二部分：核心概念

## 6. `nextflow.config` — 配置中心

[nextflow.config](nextflow.config) 是流程的**唯一配置入口**——所有参数默认值、资源分配、容器设置、Profile、报告都在这里定义。

### 6.1 段落总览

```groovy
// [nextflow.config](nextflow.config)
manifest { name = 'nextflow-rnaseq'; version = '1.0.0'; ... }   // ① 管线元数据

plugins { id 'nf-validation' }          // ② 插件加载：自动校验参数类型和格式

params {                                // ③ 用户参数默认值
    input      = null                   //    无默认值 → 用户必须提供
    read_mode  = "paired"               //    有默认值 → 可省略
    outdir     = "./results"
    trimmomatic { phred = "33"; ... }   //    嵌套参数
    deseq2 { contrasts = [...] }
}

process.container = 'rnaseq-pipeline:1.0.0'  // ④ 容器镜像

includeConfig 'conf/base.config'        // ⑤ 拆分配置（资源 + label + 错误处理）

profiles {                              // ⑥ Profile 定义
    docker    { includeConfig 'conf/docker.config' }
    apptainer { includeConfig 'conf/apptainer.config' }
    test      { includeConfig 'conf/test.config' }
    test_stub { includeConfig 'conf/test_stub.config' }
}

dag { enabled = true; ... }             // ⑦ 报告
timeline { enabled = true; ... }
report { enabled = true; ... }
```

### 6.2 params — 参数系统

#### 三个来源（优先级从高到低）

```
命令行 --xxx        最高优先级 → 覆盖配置文件
nextflow.config     默认值
params.xxx 未定义   访问时 → null
```

#### 嵌套参数覆盖

```groovy
// [nextflow.config](nextflow.config) 中定义
params {
    trimmomatic {
        phred = "33"
    }
}
```

```bash
# 命令行覆盖（点号变下划线）
nextflow run main.nf --trimmomatic_phred 64
```

#### 在项目中的使用场景

| 使用场景 | 在哪 | 示例 |
|---------|------|------|
| 定义默认值 | [nextflow.config](nextflow.config) | `params.outdir = "./results"` |
| 命令行覆盖 | 终端 | `--outdir /data/out` |
| 校验必填项 | [main.nf L48-50](main.nf#L48) | `if (!params.input) exit 1` |
| 控制流程分支 | [main.nf L108-115](main.nf#L108) | `if (params.read_mode == 'paired')` |
| 转成 channel | [main.nf L56-57](main.nf#L56) | `channel.value(file(params.genome_fasta))` |
| shell 脚本插值 | [trimmomatic_pe.nf L34](modules/local/trimmomatic_pe.nf#L34) | `${params.trimmomatic.minlen}` |
| Groovy 前导代码 | [featurecounts.nf L22](modules/local/featurecounts.nf#L22) | `def flag = params.read_mode == 'paired'` |

> **重要约束**：process 定义中不能直接引用 `$projectDir` 或 `params.xxx` 做文件路径 input。脚本文件必须通过 channel 以 `path` 类型传入。

### 6.3 Profile — 执行环境切换

**Profile 就是一组预设配置，通过 `-profile` 一键切换。**

定义在 [nextflow.config L94-107](nextflow.config#L94)，每个 profile 加载一个配置文件：
```groovy
profiles {
    docker    { includeConfig 'conf/docker.config' }
    apptainer { includeConfig 'conf/apptainer.config' }
    test      { includeConfig 'conf/test.config' }
    test_stub { includeConfig 'conf/test_stub.config' }
}
```

#### 4 个 Profile 对比

| Profile | 容器引擎 | process.container | 适用场景 |
|---------|---------|-------------------|------|
| `docker` | Docker | `rnaseq-pipeline:1.0.0` | 本地/服务器生产 |
| `apptainer` | Apptainer | `rnaseq-pipeline:1.0.0` | HPC 集群生产 |
| `test` | 无（本地 conda） | `""` | E. coli 端到端测试 |
| `test_stub` | 无 | `""` | 语法/DAG 快速校验 |

**Profile 的覆盖原理**：后加载的覆盖先加载的。例如 `-profile test` 会加载 [test.config](conf/test.config)，其中 `process.container = ''` 覆盖了主配置的 `'rnaseq-pipeline:1.0.0'`。

```bash
# 同一条命令，改 profile 就能在不同环境跑
nextflow run main.nf -profile docker ...
nextflow run main.nf -profile apptainer ...
nextflow run main.nf -profile test ...
nextflow run main.nf -profile test_stub -stub-run ...
```

### 6.4 Container — 容器镜像

本项目采用**统合镜像**策略：一个镜像包含全部 10+ 工具。

```dockerfile
# [containers/Dockerfile](containers/Dockerfile)
FROM continuumio/miniconda3:24.1.2-0

# 分两层安装（缓存友好：改命令行工具不影响 R 层重建）
RUN conda install ...  # Layer 1: R + DESeq2 + clusterProfiler + 绘图库
RUN conda install ...  # Layer 2: FastQC + Trimmomatic + HISAT2 + Samtools + ...
```

→ [nextflow.config L88](nextflow.config#L88)

```bash
# 构建
bash containers/build.sh docker      # Docker 镜像
bash containers/build.sh all         # Docker + Apptainer SIF
bash containers/build.sh test        # 验证所有工具可用


```
**容器与 process 的连接**：
```groovy
process.container = 'rnaseq-pipeline:1.0.0'
//                   ↑ 所有 process 默认用这个镜像
```

Nextflow 启动每个 task 时，自动用这个镜像创建容器，task 在容器内执行。

### 6.5 配置文件拆分

`includeConfig` 将配置按职责拆到多个文件（文本级插入，类似 C 的 `#include`）：

| 文件 | 职责 |
|------|------|
| [nextflow.config](nextflow.config) | 参数定义 + Profile 声明 + 报告 |
| [conf/base.config](conf/base.config) | 资源分配 + label + 错误处理 |
| [conf/docker.config](conf/docker.config) | Docker 引擎开关 |
| [conf/apptainer.config](conf/apptainer.config) | Apptainer 引擎开关 |
| [conf/test.config](conf/test.config) | 测试：压低资源 + 覆盖物种参数 |
| [conf/test_stub.config](conf/test_stub.config) | stub-run：极低资源 + 免容器 |

---

## 7. Process — 最小执行单元

Process 是 Nextflow 中最核心的概念——每个 process 定义了一个独立的计算步骤，每个 step 的一个输入 item 对应一个并行 task。

### 7.1 完整模板

以 [hisat2_align_pe.nf](modules/local/hisat2_align_pe.nf) 为例：

```groovy
// [hisat2_align_pe.nf](modules/local/hisat2_align_pe.nf)
process HISAT2_ALIGN_PE {                    // ① process 名称（全局唯一）
    label 'process_high'                      // ② 资源标签
    tag "${sample_id}"                        // ③ 任务标签
    publishDir "${params.outdir}/03_hisat2", mode: 'copy'  // ④ 输出发布

    input:                                    // ⑤ 输入定义
    tuple val(sample_id), path(r1), path(r2), val(prefix), path(idx_files)

    output:                                   // ⑥ 输出定义
    tuple val(sample_id),
          path("${sample_id}.sorted.bam"),
          path("${sample_id}.sorted.bam.bai"), emit: bam_bai
    path("${sample_id}.hisat2.log"),           emit: hisat2_log

    script:                                   // ⑦ 执行脚本
    def samtools_cpus = (task.cpus as int) > 2 ? 2 : 1
    """
    set -euo pipefail
    hisat2 -x ${prefix} \
        -1 ${r1} -2 ${r2} \
        -p ${task.cpus} ...
        | samtools sort -@ ${samtools_cpus} -o ${sample_id}.sorted.bam -
    samtools index ${sample_id}.sorted.bam
    """

    stub:                                     // ⑧ 测试桩（可选）
    """
    touch ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai ${sample_id}.hisat2.log
    """
}
```

### 7.2 `input:` — 输入定义

定义 process 需要哪些数据，以及如何从 channel 中接收。

#### 输入限定符

| 限定符 | 含义 | 示例 |
|--------|------|------|
| `val(x)` | 接收一个普通值 | `val(sample_id)` ← `"WT_1"` |
| `path(x)` | 接收一个文件路径（Nextflow 自动 staging 到 task 目录） | `path(r1)` ← `WT_1_R1.fastq.gz` |
| `tuple(...)` | 把多个值/文件打包成一个元组 | `tuple val(id), path(fq1), path(fq2)` |
| `each(x)` | 对集合中每个元素分别启动一个 task | `each samples` |

#### channel 数据如何匹配 input

```
channel 发送                              input 接收
─────────────────────────────────────────────────────────────
["WT_1", R1.fq, R2.fq, "genome", [*.ht2]]
   │      │      │       │         │
   ▼      ▼      ▼       ▼         ▼
tuple val(sample_id), path(r1), path(r2), val(prefix), path(idx_files)
```

**按位置匹配，不看名字。**

### 7.3 `output:` — 输出定义

定义 process 产生哪些文件/值，声明后传给下游 channel。

```groovy
output:
tuple val(sample_id),
      path("${sample_id}.sorted.bam"),
      path("${sample_id}.sorted.bam.bai"), emit: bam_bai   // ← 命名通道
path("${sample_id}.hisat2.log"),           emit: hisat2_log
```

| 机制 | 说明 |
|------|------|
| `emit: name` | 给输出通道命名，下游通过 `PROCESS.out.name` 引用 |
| `optional: true` | 文件可能不存在时不报错（如无显著基因时 clusterProfiler 不跑） |

**下游引用**：

```groovy
HISAT2_ALIGN_PE.out.bam_bai      // → channel: [sample_id, bam, bai]
HISAT2_ALIGN_PE.out.hisat2_log   // → channel: [hisat2.log]
```

### 7.4 `script:` — 执行脚本

每个 task 在隔离的工作目录中运行的 shell 脚本。

→ [hisat2_align_pe.nf L24](modules/local/hisat2_align_pe.nf#L24)

```groovy
script:
def samtools_cpus = (task.cpus as int) > 2 ? 2 : 1   // Groovy 前导代码
"""
set -euo pipefail                                      // 必须：管道错误不静默吞噬
hisat2 -x ${prefix} -1 ${r1} -2 ${r2} -p ${task.cpus} | samtools sort ...
"""


```
| 特点 | 说明 |
|------|------|
| 多行字符串 | 使用 `"""..."""` 三引号 |
| 变量插值 | `${r1}`, `${task.cpus}` 等 Groovy 变量直接插值 |
| `set -euo pipefail` | **强制要求** |
| 当前目录 | task 的独立工作目录 |
**HISAT2 管道中的 CPU 分配**（[hisat2_align_pe.nf L24](modules/local/hisat2_align_pe.nf#L24)）：
```groovy
def samtools_cpus = (task.cpus as int) > 2 ? 2 : 1
// 如果 task.cpus > 2 → samtools 用 2 核，否则用 1 核
```

为什么？管道 `hisat2 | samtools sort` 两端同时运行。hisat2 计算密集需要满配 CPU，samtools 是 IO 密集 2 核足够。避免 `2 × cpus` 超配被调度器 kill。

### 7.5 `stub:` — 测试桩

```groovy
stub:
"""
touch ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai
"""
```

配合 `-stub-run` 使用：只创建空文件，不执行真实工具。用于快速校验语法和 DAG 连线。

→ [hisat2_align_pe.nf L13](modules/local/hisat2_align_pe.nf#L13)
→ [trimmomatic_pe.nf L12](modules/local/trimmomatic_pe.nf#L12)
→ [fastqc_pe.nf L13](modules/local/fastqc_pe.nf#L13)

```bash
nextflow run main.nf -profile test_stub -stub-run ...
# 预期: SUCCESS 21/21（全部 process 执行 stub，几秒完成）


```
### 7.6 `publishDir` — 输出发布
task 成功后，将 output 中的文件从 `work/` 拷贝到用户可见的 outdir。
```groovy
// — 全部发布
publishDir "${params.outdir}/03_hisat2", mode: 'copy'

// — 过滤发布（只发布 trimmed reads）
publishDir "${params.outdir}/02_trimmomatic", mode: 'copy', pattern: "*.trimmed.fastq.gz"

// — 过滤发布（只发布报告）
publishDir "${params.outdir}/01_fastqc", mode: 'copy', pattern: "*_fastqc.{html,zip}"
```

| mode | 说明 |
|------|------|
| `'copy'` | 复制，work/ 保留原件（本项目默认，最安全） |
| `'link'` | 硬链接（同文件系统），节省磁盘 |
| `'symlink'` | 符号链接 outdir → work/ |
| `'move'` | 移动，work/ 中原件移走（影响 resume） |

### 7.7 `tag` — 任务标签

```groovy
tag "${sample_id}"
```

运行时终端显示：

```
[1c/3e8a1b] process > HISAT2_ALIGN_PE (WT_1)   ← tag 在括号里
[3f/7b2c4d] process > HISAT2_ALIGN_PE (KO_1)
```

不加 tag 时括号里是空白的。并行 task 多时，tag 让你一眼分辨谁是谁。

**tag 信息存在哪？** `.command.run` 第 3 行：

```
### name: 'HISAT2_ALIGN_PE (WT_1)'
```

```bash
# 查看所有 task 的 tag
find work/ -name ".command.run" -exec grep "### name:" {} \; | sort
```

### 7.8 cpus / memory / time — 资源标签

本项目不在 process 中直接写资源值，而是通过 **label 标签** + **[conf/base.config](conf/base.config)** 统一管理：

→ [hisat2_align_pe.nf L10](modules/local/hisat2_align_pe.nf#L10)
→ [conf/base.config](conf/base.config)

```groovy
// process 中声明 label
process HISAT2_ALIGN_PE {
    label 'process_high'
}
// conf/base.config 中定义 label 对应的资源
process {
 cpus = 2 // 默认值（未设 label 的 process 用这个）
 memory = 2.GB
 time = 4.h
 withLabel: process_low { cpus = 2; memory = 4.GB; time = 2.h }
 withLabel: process_medium { cpus = 4; memory = 16.GB; time = 6.h }
 withLabel: process_high { cpus = 8; memory = 32.GB; time = 12.h }
 withLabel: process_R { cpus = 2; memory = 8.GB; time = 4.h }
 // 错误处理：OOM(137) / HPC抢占(143) 重试，其他直接终止
 errorStrategy = { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
 maxRetries = 2
}
```

#### Label 体系总览

| Label | cpus | memory | time | 适用场景 |
|-------|:----:|:------:|:----:|------|
| `process_low` | 2 | 4 GB | 2 h | FastQC |
| `process_medium` | 4 | 16 GB | 6 h | Trimmomatic, featureCounts |
| `process_high` | 8 | 32 GB | 12 h | HISAT2 align, HISAT2 index |
| `process_R` | 2 | 8 GB | 4 h | DESeq2, clusterProfiler |
| 默认（无 label） | 2 | 2 GB | 4 h | — |

#### 错误处理

| 退出码 | 含义 | 策略 |
|--------|------|------|
| `137` | OOM kill（128+9） | retry（可能分配到更多内存的节点） |
| `143` | SIGTERM（128+15） | retry（HPC 调度器抢占） |
| 其他 | 工具 bug / segfault | terminate（重试无意义） |

#### 时间单位

| 写法 | 含义 |
|------|------|
| `2.h` | 2 小时 |
| `30.m` | 30 分钟 |
| `45.s` | 45 秒 |
| `3.d` | 3 天 |
| `4.GB` / `512.MB` | 内存单位 |

---

## 8. Channel — 数据管道

Channel 是 Nextflow 的"管道"——数据在 process 之间通过 channel 流动。**Channel 是异步队列，先进先出，只能消费一次**（`channel.value()` 除外）。

**核心认知：channel 有几个 item，process 就并行几个 task。**

```
[输入文件] → channel → process A → channel → process B → channel → [输出]
```

### 8.1 工厂方法：创建 channel

#### `channel.fromPath()` — 从文件路径创建

→ [main.nf L70](main.nf#L70)

```groovy
channel.fromPath(params.input)
// → 发送 1 个 Path 对象: samplesheet.csv
```

#### `channel.value()` — 创建可重复消费的单值 channel

→ [main.nf L56-57](main.nf#L56)

```groovy
genome_fasta_ch = channel.value(file(params.genome_fasta))
// → 发送 1 个 item: genome.fa (Path 对象)
// → 可被多个 process 重复消费，不会耗尽
```

本项目大量使用：参考基因组、GTF、R 脚本、JSON 字符串都通过它传递。

#### `channel.of()` — 手动创建 channel

```groovy
channel.of("A", "B", "C")
// → 依次发送 "A", "B", "C"
```

本项目未使用，常用在测试场景。

#### `channel.fromFilePairs()` — 自动配对 FASTQ

```groovy
channel.fromFilePairs("raw/*_{R1,R2}.fastq.gz")
// → 按 {sample}_R{1,2}.fastq.gz 模式匹配
// → tuple("WT_1", [WT_1_R1.fq, WT_1_R2.fq])
```

本项目用 CSV samplesheet 管理配对，所以没用这个。

### 8.2 转换操作符

#### `.map()` — 转换每个 item

→ [main.nf L72](main.nf#L72)
→ [main.nf L132](main.nf#L132)

```groovy
// main.nf L72 — 把 CSV 行转成 tuple
.map { row -> tuple(row.sample, row.group, file(row.fastq_1), file(row.fastq_2)) }
// 输入: {sample: "WT_1", group: "WT", ...}
// 输出: ["WT_1", "WT", WT_1_R1.fq, WT_1_R2.fq]

// main.nf L132 — 丢掉不需要的元素
.map { _sid, bam, _bai -> bam }
// 输入: ["WT_1", WT_1.bam, WT_1.bam.bai]
// 输出: WT_1.bam
```

#### `.filter()` — 保留满足条件的 item

→ [main.nf L169](main.nf#L169)

```groovy
.filter { csv -> csv.name =~ /_significant\.csv$/ }
// 输入: KO_vs_WT_all_results.csv, KO_vs_WT_significant.csv, KO_vs_WT_MA_plot.pdf
// 输出: KO_vs_WT_significant.csv
```

#### `.flatten()` — 拆开嵌套列表

→ [main.nf L168](main.nf#L168)

```groovy
.flatten()
// 输入: [a.csv, b.csv, c.csv]  ← 1 个 item（是个 List）
// 输出: a.csv, b.csv, c.csv     ← 3 个独立 item
```

当 process 的 `output` 使用 `path("*.csv")` 匹配多个文件时，这些文件会以 **单个 List 形式 emit**。下游 `.filter{}` 前必须先 `.flatten()` 拆开。

#### `.collect()` — 收集所有 item 成一个列表

→ [main.nf L133](main.nf#L133)

```groovy
.collect()
// 输入: WT_1.bam, WT_2.bam, KO_1.bam, KO_2.bam  (4 个 item)
// 输出: [WT_1.bam, WT_2.bam, KO_1.bam, KO_2.bam]  (1 个 item)
```

**`.collect()` 是 Scatter → Gather 转换的关键操作**——之前每个样本并行，之后只有一个 task。

### 8.3 合并/组合操作符

#### `.mix()` — 混合两个 channel 的 item

→ [qc_trim_align_pe.nf L42-43](subworkflows/local/qc_trim_align_pe.nf#L42)

```groovy
FASTQC_RAW_PE.out.reports
    .mix(FASTQC_TRIMMED_PE.out.reports)
// channel A: [WT_1, [raw.html, raw.zip]], [WT_2, ...], ...
// channel B: [WT_1, [trimmed.html, trimmed.zip]], [WT_2, ...], ...
// .mix() 后 — 8 个独立 item 交错输出:
//   [WT_1, [raw.html, raw.zip]]      ← FASTQC_RAW 第 1 个
//   [WT_2, [raw.html, raw.zip]]      ← FASTQC_RAW 第 2 个
//   [KO_1, [raw.html, raw.zip]]      ← FASTQC_RAW 第 3 个
//   [KO_2, [raw.html, raw.zip]]      ← FASTQC_RAW 第 4 个
//   [WT_1, [trimmed.html, trimmed.zip]]  ← FASTQC_TRIMMED 第 1 个
//   ...                              （共 8 个 item）
```

不要求 key 匹配，就是简单地把两份数据汇入同一条 channel。

#### `.combine()` — 笛卡尔积

→ [qc_trim_align_pe.nf L34-35](subworkflows/local/qc_trim_align_pe.nf#L34)

```groovy
TRIMMOMATIC_PE.out.trimmed_fq    // 4 个: [WT_1, R1, R2], [WT_2, ...], ...
    .combine(index_ch)            // 1 个: [ecoli, [*.ht2]]
// 输出 (4×1=4):
//   [WT_1, R1, R2, ecoli, [*.ht2]]
//   [WT_2, R1, R2, ecoli, [*.ht2]]
//   [KO_1, R1, R2, ecoli, [*.ht2]]
//   [KO_2, R1, R2, ecoli, [*.ht2]]
```

每个样本 × 每个索引。本项目索引只有 1 个，实际效果就是"给每个样本附上索引信息"。

#### `.join()` — 按 key 关联（SQL JOIN）

```groovy
// channel A: [sample_id, fastq_1, fastq_2]
// channel B: [sample_id, group]
// .join() → [sample_id, fastq_1, fastq_2, group]
// 只保留两边都存在的 key
```

本项目未直接使用（用 `.multiMap()` 替代了 join 的需求）。

### 8.4 分流操作符

#### `.multiMap()` — 一个 channel 广播到多个下游

本项目最核心的分流操作（[main.nf L84-89](main.nf#L84)）：

```groovy
sample_input_ch
    .multiMap { it ->
        fq:         tuple(it[0], it[2], it[3])   // 给比对流程
        group:      tuple(it[0], it[1])           // 预留（未使用）
        validation: it[1]                          // 给对比校验
    }
    .set { forked }

// 输入: ["WT_1", "WT", R1.fq, R2.fq]
// 输出三个分支:
//   forked.fq         → ["WT_1", R1.fq, R2.fq]
//   forked.group      → ["WT_1", "WT"]
//   forked.validation → "WT"
```

一次解析 samplesheet，三路同时使用，无需重复读文件。

### 8.5 其他操作符

#### `.ifEmpty()` — channel 为空时的处理

→ [main.nf L78](main.nf#L78)
→ [main.nf L171](main.nf#L171)

```groovy
// main.nf L78
.ifEmpty { exit 1, "No samples found in ${params.input}" }

// main.nf L171
.ifEmpty { log.warn "[main] No significant DEGs found — skipping clusterProfiler"; [] }
```

#### `.unique()` — 去重

→ [main.nf L93](main.nf#L93)

```groovy
forked.validation       // "WT", "WT", "KO", "KO"
    .unique()           // "WT", "KO"
```

#### `.view()` / `.dump()` — 调试打印

```groovy
// 看看 channel 里有什么，不改变数据
my_ch.view { "DEBUG: $it" }
my_ch.dump("tag_name")
```

### 8.6 本项目 Channel 操作符使用情况

| 操作符 | [main.nf](main.nf) | [qc_trim_align_pe.nf](subworkflows/local/qc_trim_align_pe.nf) |
|--------|:---:|:---:|
| `channel.fromPath()` | ✅ L70 | |
| `channel.value()` | ✅ (大量) | |
| `.splitCsv()` | ✅ L71 | |
| `.map()` | ✅ L72,95,120,124,132 | ✅ L36,44 |
| `.filter()` | ✅ L169 | |
| `.flatten()` | ✅ L168 | ✅ L45 |
| `.collect()` | ✅ L94,133,170 | ✅ L46 |
| `.combine()` | | ✅ L35 |
| `.mix()` | | ✅ L43 |
| `.multiMap()` | ✅ L84 | |
| `.set()` | ✅ L89 | ✅ L37 |
| `.ifEmpty()` | ✅ L78,171 | |
| `.unique()` | ✅ L93 | |
| `.set()` | ✅ L89 | ✅ L37 |
| `.ifEmpty()` | ✅ L78,171 | |
| `.unique()` | ✅ L93 | |
---
## 9. 多样本处理 — Scatter-Gather 模式
Nextflow 不需要显式写 `for` 循环——**一个 channel 有几个 item，process 就执行几次 task。**
### 9.1 Scatter：每个样本一个 task

→ [qc_trim_align_pe.nf L25](subworkflows/local/qc_trim_align_pe.nf#L25)

```groovy
FASTQC_RAW_PE(sample_ch)     // sample_ch 有 4 个 item
```

Nextflow 自动启动 4 个并行 task：

```
item 1: ["WT_1", R1, R2] → task 1: fastqc WT_1_R1.fq WT_1_R2.fq  ┐
item 2: ["WT_2", R1, R2] → task 2: fastqc WT_2_R1.fq WT_2_R2.fq  │ 4 task
item 3: ["KO_1", R1, R2] → task 3: fastqc KO_1_R1.fq KO_1_R2.fq  │ 同时
item 4: ["KO_2", R1, R2] → task 4: fastqc KO_2_R1.fq KO_2_R2.fq  ┘
```

### 9.2 Gather：`.collect()` 汇总

→ [main.nf L131-134](main.nf#L131)

```groovy
all_bams_ch = bam_ch
    .map { _sid, bam, _bai -> bam }   // 丢 sample_id 和 bai
    .collect()                          // → [WT_1.bam, WT_2.bam, KO_1.bam, KO_2.bam]

FEATURECOUNTS(all_bams_ch, gtf_ch)
// → 只跑 1 个 task: featureCounts -a gtf -o counts.txt WT_1.bam WT_2.bam KO_1.bam KO_2.bam
```

### 9.3 本项目的分层

```
Layer 1-3: SCATTER（并行）          Layer 4-7: GATHER（单任务）
─────────────────────────          ─────────────────────────
FASTQC      4 task 并行
TRIMMOMATIC 4 task 并行
FASTQC      4 task 并行
HISAT2      4 task 并行
         ─── .collect() ───→       FEATURECOUNTS    1 task
                                   MULTIQC          1 task
                                   DESEQ2           1 task
                                   CLUSTERPROFILER  1 task
```

### 9.4 时间轴可视化

```
时间 ───────────────────────────────────────────────→

FASTQC_RAW     ██ WT_1 ██ WT_2 ██ KO_1 ██ KO_2
TRIMMOMATIC    ████ WT_1 ████ WT_2 ████ KO_1 ████ KO_2
FASTQC_TRIMMED ██ WT_1 ██ WT_2 ██ KO_1 ██ KO_2
HISAT2_ALIGN   ████████ WT_1 ████████ WT_2 ████████ KO_1 ████████ KO_2
               ──────────── 全部完成，.collect() 汇集 ────────────
FEATURECOUNTS  ████████████████████████████████  1 task
MULTIQC        ████                              1 task
DESEQ2         ██████                            1 task
CLUSTERPROFILER ████                             1 task


```
### 9.5 特殊处理：`.combine()` 给每个样本补上共用数据

→ [qc_trim_align_pe.nf L34-37](subworkflows/local/qc_trim_align_pe.nf#L34)

```groovy
TRIMMOMATIC_PE.out.trimmed_fq       // 4 个样本
    .combine(index_ch)               // 1 个索引
    .map { sid, r1, r2, pfx, idx -> tuple(sid, r1, r2, pfx, idx) }
```

每个人拿到一份索引，4 个 task 可以并行比对。

---

## 10. Workflow — 编排器

Workflow 是 process 的编排器——定义哪些 process 按什么顺序执行，它们之间的 channel 怎么连接。它自身不执行计算，只负责"连线"。

### 10.1 结构：三部分

```groovy
workflow 名字 {
    take:    // 从外部接收什么 channel（函数参数）
    参数名

    main:   // 内部逻辑：调用 process，串联 channel（函数体）
    PROCESS_A(input_ch)
    PROCESS_B(PROCESS_A.out.xxx)

    emit:   // 向外部暴露什么 channel（返回值）
    输出名 = PROCESS_B.out.yyy
}


```
### 10.2 主流程（无名 workflow）

[main.nf L45](main.nf#L45) — 没有名字的顶层 workflow 就是入口：

```groovy
workflow {              // ← 无名 = 主入口
    // 无 take — 直接从 params 和 channel factory 读取
    gtf_ch = channel.value(file(params.annotation_gtf))

    // 校验 + channel 准备 + 流程编排
    HISAT2_INDEX(genome_fasta_ch)
    result = QC_TRIM_ALIGN_PE(fq_ch, adapter_ch, index_ch)
    FEATURECOUNTS(all_bams_ch, gtf_ch)
    DESEQ2(...)
    CLUSTERPROFILER(...)
}


```
主流程不支持 `take:` 和 `emit:`——它是入口，输入来自 `params` 和 channel factory，输出通过 `publishDir` 写盘。
### 10.3 子流程（命名 workflow）

[qc_trim_align_pe.nf L17-51](subworkflows/local/qc_trim_align_pe.nf#L17) — 把 4 个上游步骤封装成一个可复用单元：

```groovy
workflow QC_TRIM_ALIGN_PE {
    take:                              // ← 入参：3 个 channel
    sample_ch                          //    (sample_id, fastq_1, fastq_2)
    adapter_ch                         //    adapter Fasta 文件
    index_ch                           //    HISAT2 索引

    main:                              // ← 内部连线
    FASTQC_RAW_PE(sample_ch)
    TRIMMOMATIC_PE(sample_ch, adapter_ch)
    FASTQC_TRIMMED_PE(TRIMMOMATIC_PE.out.trimmed_fq)
    TRIMMOMATIC_PE.out.trimmed_fq
        .combine(index_ch)
        .map { sid, r1, r2, pfx, idx -> tuple(sid, r1, r2, pfx, idx) }
        .set { align_in }
    HISAT2_ALIGN_PE(align_in)
    qc_files = FASTQC_RAW_PE.out.reports
        .mix(FASTQC_TRIMMED_PE.out.reports)
        .map { _sample, htmls, zips -> htmls + zips }
        .flatten().collect()

    emit:                              // ← 出参：2 个 channel
    bam_bai  = HISAT2_ALIGN_PE.out.bam_bai
    qc_files = qc_files
}


```
**为什么要封装成子流程？** PE 和 SE 模式的 4 个步骤相同但 process 不同。封装后 main.nf 只需一行调用（[main.nf L118-126](main.nf#L118)）：

```groovy
if (params.read_mode == 'paired') {
    result = QC_TRIM_ALIGN_PE(fq_ch, adapter_ch, index_ch)
} else {
    result = QC_TRIM_ALIGN_SE(fq_ch, adapter_ch, index_ch)
}
bam_ch    = result.bam_bai       // ← 通过 .属性名 取子流程 emit
qc_all_ch = result.qc_files


```
调用方只需要知道：传 3 个 channel 进去，拿 2 个 channel 出来——内部 4 个 process 怎么连线的，外部完全不需要关心。
### 10.4 Workflow vs Process
| | process | workflow |
|------|---------|----------|
| 做什么 | 执行一个计算步骤 | 编排多个 step |
| 内容 | `input:` + `output:` + `script:` | `take:` + `main:` + `emit:` |
| 内部 | shell 脚本 | 调用 process / 子 workflow |
| 产出 | 文件 | channel |
| 并行 | 每个 input item 一个 task | 不直接产生 task |
---
## 11. DSL2 模块系统
DSL2 是 Nextflow 的模块化语法（对比 DSL1 的单文件脚本模式）。本项目第一行就声明了（[main.nf L20](main.nf#L20)）：

```groovy
nextflow.enable.dsl=2
```

### 11.1 Module — 一个 .nf 文件 = 一个 process

```
modules/local/
├── hisat2_index.nf          ← process HISAT2_INDEX
├── hisat2_align_pe.nf       ← process HISAT2_ALIGN_PE
├── hisat2_align_se.nf       ← process HISAT2_ALIGN_SE
├── fastqc_pe.nf             ← process FASTQC_PE
├── fastqc_se.nf             ← process FASTQC_SE
├── trimmomatic_pe.nf        ← process TRIMMOMATIC_PE
├── trimmomatic_se.nf        ← process TRIMMOMATIC_SE
├── featurecounts.nf         ← process FEATURECOUNTS
├── multiqc.nf               ← process MULTIQC
├── deseq2.nf                ← process DESEQ2
└── clusterprofiler.nf       ← process CLUSTERPROFILER
```

一个文件 = 一个 process = 一个 Module。每个 Module 结构完全统一：`label` + `tag` + `publishDir` + `input:` + `output:` + `script:` + `stub:`。

### 11.2 `include {}` — 导入语句

#### 基本导入

→ [main.nf L27](main.nf#L27)

```groovy
include { HISAT2_INDEX } from './modules/local/hisat2_index'
```


#### 别名导入：同一个 Module 实例化两次

→ [qc_trim_align_pe.nf L12-13](subworkflows/local/qc_trim_align_pe.nf#L12)

```groovy
include { FASTQC_PE as FASTQC_RAW_PE      } from '../../modules/local/fastqc_pe'
include { FASTQC_PE as FASTQC_TRIMMED_PE  } from '../../modules/local/fastqc_pe'
```

`fastqc_pe.nf` 只有一个 process，但流程中需要用两次——原始 reads 和剪切后 reads。`as` 别名让同一个 Module 产生两个独立实例：

```
fastqc_pe.nf
  ├──→ FASTQC_RAW_PE      输入: 原始 FASTQ,  tag: WT_1
  └──→ FASTQC_TRIMMED_PE  输入: 剪切后 FASTQ, tag: WT_1


```
两个实例互不影响，有自己的 input channel、tag、task 哈希。
#### 导入 Subworkflow
```groovy
include { QC_TRIM_ALIGN_PE } from './subworkflows/local/qc_trim_align_pe'
include { QC_TRIM_ALIGN_SE } from './subworkflows/local/qc_trim_align_se'
```

#### PE/SE 独立 Module

fastqc、trimmomatic、hisat2_align 各有两个文件（PE/SE）。因为 `input:`/`output:`/`script:` 完全不同（文件数量、参数都不同），拆成独立 Module 比在 process 内部 `if/else` 更清晰。

### 11.3 DSL2 的关键规则

1. **所有 process 和 workflow 必须通过 `include {}` 导入后才能调用**
2. **channel 只能消费一次**（除 `channel.value()` 可重复消费）
3. **调用 process 时传 channel，参数数量和位置与 `input:` 对应**
4. **`as` 只支持 process，workflow 不支持**——所以 PE/SE 需要两个独立 subworkflow 文件

### 11.4 本项目 Include 完整树

```
[main.nf](main.nf)
├── include { HISAT2_INDEX      } from modules/local/hisat2_index.nf
├── include { QC_TRIM_ALIGN_PE  } from subworkflows/.../qc_trim_align_pe.nf
│   └── 内部:
│       ├── include { FASTQC_PE as FASTQC_RAW_PE      } from fastqc_pe.nf
│       ├── include { FASTQC_PE as FASTQC_TRIMMED_PE  } from fastqc_pe.nf
│       ├── include { TRIMMOMATIC_PE } from trimmomatic_pe.nf
│       └── include { HISAT2_ALIGN_PE } from hisat2_align_pe.nf
├── include { QC_TRIM_ALIGN_SE  } from subworkflows/.../qc_trim_align_se.nf
├── include { FEATURECOUNTS     } from modules/local/featurecounts.nf
├── include { MULTIQC           } from modules/local/multiqc.nf
├── include { DESEQ2            } from modules/local/deseq2.nf
└── include { CLUSTERPROFILER   } from modules/local/clusterprofiler.nf


```
---
# 第三部分：实战串联
## 12. main.nf 逐段解读
以 4 个样本为例，按执行顺序逐步拆解完整数据流。
### 12.1 启动校验 + 参考文件 channel

→ [main.nf L48-66](main.nf#L48)

```groovy
if (!params.input)           { exit 1, "Missing required parameter: --input" }
if (!params.genome_fasta)    { exit 1, "Missing required parameter: --genome_fasta" }
if (!params.annotation_gtf)  { exit 1, "Missing required parameter: --annotation_gtf" }

genome_fasta_ch  = channel.value(file(params.genome_fasta))
gtf_ch           = channel.value(file(params.annotation_gtf))
deseq2_r_ch      = channel.value(file("${projectDir}/bin/deseq2.R"))


```
`channel.value()` 创建可重复消费的 channel。基因组文件要同时给 HISAT2 建索引、给 DESeq2 做注释，每个 process 都能取到。
### 12.2 样本输入 channel

→ [main.nf L70-78](main.nf#L70)

```groovy
sample_input_ch = channel.fromPath(params.input)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample, row.group, file(row.fastq_1), file(row.fastq_2)) }
    .ifEmpty { exit 1, "No samples found in ${params.input}" }
```

> 假设 `samplesheet.csv` 内容：
>
> sample,group,fastq_1,fastq_2
> WT_1,WT,data/WT_1_R1.fastq.gz,data/WT_1_R2.fastq.gz
> WT_2,WT,data/WT_2_R1.fastq.gz,data/WT_2_R2.fastq.gz
> KO_1,KO,data/KO_1_R1.fastq.gz,data/KO_1_R2.fastq.gz
> KO_2,KO,data/KO_2_R1.fastq.gz,data/KO_2_R2.fastq.gz

**每一步的数据形态**：

```
① channel.fromPath("samplesheet.csv")
   → [samplesheet.csv]  1 个 Path 对象

② .splitCsv(header: true)
   → {sample: "WT_1", group: "WT", fastq_1: "data/WT_1_R1.fastq.gz", ...}
   → {sample: "WT_2", group: "WT", ...}
   → {sample: "KO_1", group: "KO", ...}
   → {sample: "KO_2", group: "KO", ...}    4 个 Map

③ .map { row -> tuple(...) }
   → ["WT_1", "WT", WT_1_R1.fastq.gz, WT_1_R2.fastq.gz]
   → ["WT_2", "WT", WT_2_R1.fastq.gz, WT_2_R2.fastq.gz]
   → ["KO_1", "KO", KO_1_R1.fastq.gz, KO_1_R2.fastq.gz]
   → ["KO_2", "KO", KO_2_R1.fastq.gz, KO_2_R2.fastq.gz]  4 个 tuple


```
### 12.3 multiMap 三路分流

→ [main.nf L83-89](main.nf#L83)

```groovy
sample_input_ch
    .multiMap { it ->
        fq:         tuple(it[0], it[2], it[3])   // [sample, R1, R2]
        group:      tuple(it[0], it[1])           // [sample, group]
        validation: it[1]                          // group 名
    }
    .set { forked }
```

**以第一个 item `["WT_1", "WT", WT_1_R1.fastq.gz, WT_1_R2.fastq.gz]` 为例**：

```
    输入: ["WT_1", "WT", WT_1_R1.fastq.gz, WT_1_R2.fastq.gz]
              │
       multiMap 分三路
              │
     ├── forked.fq:         ["WT_1", WT_1_R1.fastq.gz, WT_1_R2.fastq.gz]
     ├── forked.group:      ["WT_1", "WT"]
     └── forked.validation: "WT"


```
### 12.4 对比校验 — unique + collect + map 联合使用

→ [main.nf L92-104](main.nf#L92)

```groovy
forked.validation
    .unique()        // 去重: "WT", "WT", "KO", "KO" → "WT", "KO"
    .collect()       // 汇集: "WT", "KO" → ["WT", "KO"]
    .map { valid_groups ->
        def group_set = valid_groups as Set   // → {"WT", "KO"}
        params.deseq2.contrasts.each { c ->
            if (!group_set.contains(c['case']))
                exit 1, "Contrast '${c.name}': case group '${c['case']}' not found in samplesheet"
            if (!group_set.contains(c.control))
                exit 1, "Contrast '${c.name}': control group '${c.control}' not found in samplesheet"
        }
        log.info "[main] Contrast validation passed: groups=${group_set}"
    }


```
**为什么 `.unique()` + `.collect()` + `.map()` 要配合使用？**
- `.unique()` 去重 → 得到不重复的 group 名
- `.collect()` 汇集 → 把多个 item 变成一个 List
- `.map()` 校验 → 拿到完整 List 后，一次性检查所有 contrasts
如果不 `.collect()`，`.map()` 会**对每个 item 单独执行一次**——每来一个 group 都跑一遍校验循环（不合理）。汇集后只执行一次。
**为什么不在每样本并行步骤后做？** 在流程启动时就验证，group 拼写错误 2 秒内报错，不用等到 40 分钟上游跑完到 DESeq2 才崩。
### 12.5 HISAT2 索引

→ [main.nf L108-115](main.nf#L108)

```groovy
if (params.hisat2_index_prefix && file("${params.hisat2_index_prefix}.1.ht2").exists()) {
    def idx_name = params.hisat2_index_prefix.tokenize('/')[-1]
    index_ch = channel.value(tuple(idx_name, files("${params.hisat2_index_prefix}.*.ht2")))
    log.info "[main] Using pre-built HISAT2 index: ${params.hisat2_index_prefix}"
} else {
    HISAT2_INDEX(genome_fasta_ch)
    index_ch = HISAT2_INDEX.out.index
}


```
用户提供了预建索引 → 直接复用；未提供 → 从 FASTA 自动构建。
### 12.6 Layer 1-3：Scatter 阶段

→ [main.nf L118-126](main.nf#L118)

```groovy
if (params.read_mode == 'paired') {
    result = QC_TRIM_ALIGN_PE(
        forked.fq.map { sid, fq1, fq2 -> tuple(sid, fq1, fq2) },
        adapter_ch, index_ch)
} else {
    result = QC_TRIM_ALIGN_SE(...)
}
bam_ch    = result.bam_bai
qc_all_ch = result.qc_files


```
子流程内部的 `.combine()`（[qc_trim_align_pe.nf L34-35](subworkflows/local/qc_trim_align_pe.nf#L34)）给每个样本附上索引信息：

```
TRIMMOMATIC_PE.out.trimmed_fq (4 个):
  [WT_1, R1.trimmed, R2.trimmed], [WT_2, ...], [KO_1, ...], [KO_2, ...]

index_ch (1 个):
  [ecoli, [*.ht2]]

.combine() → 4×1=4:
  [WT_1, R1.trimmed, R2.trimmed, ecoli, [*.ht2]]
  [WT_2, R1.trimmed, R2.trimmed, ecoli, [*.ht2]]
  [KO_1, R1.trimmed, R2.trimmed, ecoli, [*.ht2]]
  [KO_2, R1.trimmed, R2.trimmed, ecoli, [*.ht2]]


```
子流程尾部 mix + flatten + collect 汇总所有 QC 文件（[qc_trim_align_pe.nf L42-46](subworkflows/local/qc_trim_align_pe.nf#L42)）：

```
FASTQC_RAW.reports (4) .mix() FASTQC_TRIMMED.reports (4) → 8 个
  .map { 丢 sample_id, htmls+zips 合并 } → 8 个文件列表（每个 4 文件）
  .flatten()  → 8×4=32 个独立文件
  .collect()  → 1 个包含所有 32 个文件的列表 → 传给 MULTIQC


```
### 12.7 Layer 4-6：Gather 阶段

→ [main.nf L131-160](main.nf#L131)

```groovy
// featureCounts
all_bams_ch = bam_ch
    .map { _sid, bam, _bai -> bam }   // 丢 sample_id 和 bai
    .collect()                          // → [WT_1.bam, WT_2.bam, KO_1.bam, KO_2.bam]
FEATURECOUNTS(all_bams_ch, gtf_ch)     // 1 个 task，所有 BAM 一起计数

// MultiQC
MULTIQC(qc_all_ch)                     // 1 个 task

// DESeq2
// 分组 JSON — Groovy 读取 CSV 构建
def _sample_groups = [:] as LinkedHashMap
new File(params.input).withReader { r ->
    r.eachLine { line, idx ->
        if (idx == 1) return              // skip header
        def cols = line.split(',')
        def sid = cols[0].trim()
        def grp = cols[1].trim()
        if (!_sample_groups.containsKey(grp))
            _sample_groups[grp] = []
        _sample_groups[grp] << sid
    }
}
// → {"WT":["WT_1","WT_2"],"KO":["KO_1","KO_2"]}
groups_ch = channel.value(new groovy.json.JsonBuilder(_sample_groups).toString())
contrasts_ch = channel.value(
    new groovy.json.JsonBuilder(params.deseq2.contrasts).toString()
)
DESEQ2(FEATURECOUNTS.out.counts, gtf_ch, groups_ch, contrasts_ch,
       common_r_ch, deseq2_r_ch, gene2symbol_r_ch)


```
> **为什么 DESeq2 的分组 JSON 要重新读文件，而不是用 `forked.group` channel 构建？** 因为 `forked.group` 转成 `{group: [samples]}` 需要 `groupTuple + collect + map` 链式操作，代码繁琐。samplesheet 只有几十行，Groovy 直接读文件几行代码意图更清晰。
### 12.8 Layer 7：clusterProfiler 的 flatten+filter+collect+ifEmpty 链路

→ [main.nf L167-173](main.nf#L167)

```groovy
sig_csvs_ch = DESEQ2.out.tables
    .flatten()                                           // 拆开 List
    .filter { csv -> csv.name =~ /_significant\.csv$/ }  // 只保留显著基因文件
    .collect()                                           // 汇集
    .ifEmpty { log.warn "[main] No significant DEGs found in any contrast — skipping clusterProfiler"; [] }
// 如果所有对比都没有显著基因 → 打 warning → 给空列表 → CLUSTERPROFILER 不会启动

CLUSTERPROFILER(sig_csvs_ch, contrast_names_ch, common_r_ch, clusterprofiler_r_ch)
```

**每一步的数据变化**（假设 KO_vs_WT 有显著基因）：

```
① DESEQ2.out.tables
   → [KO_vs_WT_all_results.csv, KO_vs_WT_significant.csv, KO_vs_WT_MA_plot.pdf, ...]
   （1 个 item，是 List）

② .flatten()
   → KO_vs_WT_all_results.csv       (独立 item)
   → KO_vs_WT_significant.csv       (独立 item)
   → KO_vs_WT_MA_plot.pdf           (独立 item)
   → ...

③ .filter { csv -> csv.name =~ /_significant\.csv$/ }
   → KO_vs_WT_significant.csv  ✅ 匹配正则
   → 其他文件被过滤 ❌

④ .collect()
   → [KO_vs_WT_significant.csv]     (1 个 item，List)

⑤ .ifEmpty { ...; [] }
   → 不空，原样传递
```

### 12.9 完整数据流图


```
samplesheet.csv
 │ fromPath → splitCsv → map → ifEmpty
 │ │
 │ sample_input_ch (4 个 tuple)
 │ │
 │ .multiMap { fq / group / validation }
 │ │ │
 │ │ └─→ unique → collect → map (校验)
 │ │
 │ ▼
 │ forked.fq (给比对)
 │ │
 │ ▼
 │ 子流程
 │ │ │
 │ │ .combine(index_ch) │ .mix(raw + trimmed)
 │ │ .map() 重组 │ .map() → .flatten() → .collect()
 │ │ │
 │ ▼ ▼
 │ HISAT2_ALIGN_PE MULTIQC ← 所有 QC 文件
 │ │
 │ ▼
 │ bam_ch (4 个 BAM)
 │ │
 │ ├─ .map { 只留 bam } → .collect { 所有 BAM }
 │ ▼
 │ FEATURECOUNTS → counts.txt
 │ │
 │ ▼
 │ DESEQ2 → all_results.csv + significant.csv + 图
 │ │
 │ └─ .flatten() → .filter{significant} → .collect() → .ifEmpty{}
 │ │
 │ ▼
 │ CLUSTERPROFILER
 │
 └─ channel.value() × N (参考文件、R 脚本、JSON — 直接送入各 process)
```

---

# 第四部分：附录

## 13. 错误处理与调试

### 13.1 如何定位失败的 task

```bash
# 方法一：从 trace.txt 筛选
grep "FAILED" results/pipeline_trace.txt
awk -F'\t' '$5=="FAILED"' results/pipeline_trace.txt

# 方法二：查找非 0 退出码
find work/ -name ".exitcode" -exec grep -lv "^0$" {} \;

# 方法三：从 .nextflow.log 查找
grep -i "error\|failed\|killed" .nextflow.log | tail -20
```

### 13.2 查看失败 task 的错误信息

```bash
# 假设 task 目录是 work/3f/a1b2c3d4/

# stderr（工具报错通常在这里）
cat work/3f/a1b2c3d4/.command.err

# stdout + Nextflow stage/unstage 日志
cat work/3f/a1b2c3d4/.command.out

# 完整输出（.command.err + .command.out 合并）
cat work/3f/a1b2c3d4/.command.log

# 退出码
cat work/3f/a1b2c3d4/.exitcode
```

### 13.3 本项目错误策略

配置在 [conf/base.config](conf/base.config)：

```groovy
process {
    errorStrategy = { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries    = 2
}
```

| 退出码 | 含义 | 策略 | 原因 |
|--------|------|------|------|
| `137` | OOM kill（128+9） | retry (2次) | 换节点可能有更多内存 |
| `143` | SIGTERM（128+15） | retry (2次) | HPC 调度器抢占，换节点 |
| 其他 | 工具 bug / segfault | terminate | 重试无意义 |

### 13.4 常见失败排查

| 症状 | 排查路径 |
|------|---------|
| OOM | 看 trace.txt 的 `peak_rss` vs 配置的 `memory`；调大 memory 或降 cpu |
| 工具报错 | 看 `.command.err`；检查输入文件是否完整 |
| 输入文件不存在 | 看 `.command.run` 的 stage 日志；检查路径 |
| 任务超时 | 看 trace.txt 的 `duration` vs 配置的 `time`；调大 time |
| 缓存不命中 | 检查是否修改了 script / input 内容 / 容器镜像 |

---

## 14. Groovy 速查

### 14.1 三元表达式

```groovy
condition ? true_value : false_value

// ⚠️ 坑：数字 0 是 falsy！
def x = 0 ? "yes" : "no"    // → "no"

// ✅ 安全做法：显式 as int
def cpus = (task.cpus as int) > 2 ? 2 : 1
```

### 14.2 Elvis 运算符 `?:`

→ [featurecounts.nf L22](modules/local/featurecounts.nf#L22)

```groovy
def result = value ?: defaultValue
// 等价于
def result = value ? value : defaultValue
// 本项目使用：
def extra = params.featurecounts.extra ?: ''
```

### 14.3 `.each{}` 遍历

```groovy
// Python: for item in list
// Groovy:
["A", "B", "C"].each { item -> println item }

// 不写参数名，默认叫 it
["A", "B"].each { println it }
```

### 14.4 `.collect{}` 遍历+收集返回值

```groovy
// .each{} — 不收集返回值
["a", "b"].each { it.toUpperCase() }    // → ["a", "b"]  原样！

// .collect{} — 收集返回值
["a", "b"].collect { it.toUpperCase() } // → ["A", "B"]  新列表！


```

### 14.5 与 Java 互操作

→ [main.nf L96](main.nf#L96)
→ [main.nf L143-152](main.nf#L143)
→ [main.nf L153-156](main.nf#L153)
→ [main.nf L163-165](main.nf#L163)

```groovy
// 转 Set
def s = ["WT", "KO"] as Set
s.contains("WT")  // → true

// Map 字面量
def m = [:] as LinkedHashMap     // LinkedHashMap 保持插入顺序
m["key"] = []
m["key"] << "value"              // << = 追加到列表

// 读文件
new File("samplesheet.csv").withReader { r ->   // withReader 自动关闭
    r.eachLine { line, idx ->                   // idx 从 1 开始
        if (idx == 1) return                    // 跳过 header
        def cols = line.split(',')
    }
}

// 构建 JSON
new groovy.json.JsonBuilder([WT: ["WT_1", "WT_2"]]).toString()
// → '{"WT":["WT_1","WT_2"]}'

// spread 运算符
params.deseq2.contrasts*.name
// 等价于 params.deseq2.contrasts.collect { it.name }
```

### 14.6 常见陷阱

| 陷阱 | 说明 | 对策 |
|------|------|------|
| `0` 是 falsy | `def x = 0 ? "a" : "b"` → `"b"` | 用 `!= null` 显式判空 |
| `""` 是 falsy | 空串被当作 false | 同上 |
| `task.cpus` 可能是 String | `"8"` 和 `2` 比较不可靠 | 用 `as int` 显式转换 |
| `==` 调的是 `equals()` | Groovy 中 `==` 比较值而非引用 | 引用比较用 `is()` |

---

## 15. 速查表

### 15.1 Channel 操作符速查

| 操作符 | 作用 | 位置（示例） |
|--------|------|------|
| `channel.fromPath()` | 从文件路径创建 channel | [main.nf L70](main.nf#L70) |
| `channel.value()` | 创建可重复消费的单值 channel | [main.nf L56-57](main.nf#L56) |
| `.splitCsv()` | 解析 CSV 文件 | [main.nf L71](main.nf#L71) |
| `.map()` | 转换每个 item | [main.nf L72](main.nf#L72) |
| `.filter()` | 保留满足条件的 item | [main.nf L169](main.nf#L169) |
| `.flatten()` | 拆开嵌套列表 | [main.nf L168](main.nf#L168) |
| `.collect()` | 收集所有 item 成列表 | [main.nf L133](main.nf#L133) |
| `.mix()` | 混合两个 channel | [qc_trim_align_pe.nf L43](subworkflows/local/qc_trim_align_pe.nf#L43) |
| `.combine()` | 笛卡尔积 | [qc_trim_align_pe.nf L35](subworkflows/local/qc_trim_align_pe.nf#L35) |
| `.multiMap()` | 一个 channel 广播到多个下游 | [main.nf L84](main.nf#L84) |
| `.unique()` | 去重 | [main.nf L93](main.nf#L93) |
| `.ifEmpty()` | channel 为空时的处理 | [main.nf L78](main.nf#L78) |
| `.set()` | 给 channel 命名 | [main.nf L89](main.nf#L89) |
| `.view()` / `.dump()` | 调试打印 | — |

### 15.2 Process Directive 速查

| directive | 作用 | 定义位置 |
|-----------|------|------|
| `label` | 资源标签 | process 定义中（如 [hisat2_align_pe.nf L10](modules/local/hisat2_align_pe.nf#L10)） |
| `tag` | 任务标签 | process 定义中（如 [hisat2_align_pe.nf L11](modules/local/hisat2_align_pe.nf#L11)） |
| `publishDir` | 输出发布路径 + 模式 | process 定义中（如 [hisat2_align_pe.nf L12](modules/local/hisat2_align_pe.nf#L12)） |
| `cpus` | CPU 核数 | [conf/base.config](conf/base.config) (`withLabel`) |
| `memory` | 内存限制 | [conf/base.config](conf/base.config) |
| `time` | 时间限制 | [conf/base.config](conf/base.config) |
| `errorStrategy` | 失败策略 | [conf/base.config L11](conf/base.config#L11) |
| `maxRetries` | 最大重试次数 | [conf/base.config L12](conf/base.config#L12) |
| `cache` | 缓存策略 | [hisat2_index.nf L10](modules/local/hisat2_index.nf#L10) (`'deep'`) |
| `container` | 容器镜像 | [nextflow.config L88](nextflow.config#L88) (全局) |

### 15.3 常用命令速查

```bash
# 语法校验
nextflow run main.nf -profile test_stub -stub-run ...

# 小数据测试
nextflow run main.nf -profile test ...

# 生产运行
nextflow run main.nf -profile docker --input samplesheet.csv ...

# 中断恢复
nextflow run main.nf -profile docker ... -resume

# 清理 work 目录
nextflow clean -f

# 查看参数帮助
nextflow run main.nf --help

# 查看所有 task 的 tag
find work/ -name ".command.run" -exec grep "### name:" {} \; | sort

# 查找失败 task
grep "FAILED" results/pipeline_trace.txt

# 查看某个 task 的错误输出
cat work/xx/hash/.command.err
```

### 15.4 项目文件导航

| 文件 | 内容 |
|------|------|
| [main.nf](main.nf) | 主工作流入口 |
| [nextflow.config](nextflow.config) | 全局配置中心 |
| [conf/base.config](conf/base.config) | 资源分配 + label + 错误处理 |
| [conf/docker.config](conf/docker.config) | Docker 引擎开关 |
| [conf/apptainer.config](conf/apptainer.config) | Apptainer 引擎开关 |
| [conf/test.config](conf/test.config) | E. coli 测试配置 |
| [conf/test_stub.config](conf/test_stub.config) | stub-run 配置 |
| [modules/local/](modules/local/) | 11 个 process Module |
| [subworkflows/local/](subworkflows/local/) | 2 个子流程（PE/SE） |
| [bin/](bin/) | 4 个 R 脚本 |
| [containers/](containers/) | Dockerfile + Apptainer.def + build.sh |
| [test/](test/) | E. coli 测试数据 |
