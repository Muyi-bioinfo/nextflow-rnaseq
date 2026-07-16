# RNA-seq Pipeline (Nextflow)

[English](README.md) | [中文](README_CN.md)

[![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-0dc09d)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Docker-Supported-blue)](https://www.docker.com/)
[![Apptainer](https://img.shields.io/badge/Apptainer-Supported-orange)](https://apptainer.org/)

A Nextflow DSL2 RNA-seq analysis pipeline covering raw FASTQ QC through differential expression and functional enrichment, with full containerization support.

**Toolchain**: FastQC → Trimmomatic → HISAT2 → featureCounts → MultiQC → DESeq2 → clusterProfiler (+ gene symbol annotation, GO/KEGG enrichment)

## Features

- **PE/SE auto-adapt** — one `--read_mode` switch, all modules and subworkflows adapt
- **Modular DSL2 design** — 11 process modules + 2 subworkflows, each with consistent `label`/`tag`/`publishDir`/`stub` pattern
- **Numbered outputs** — `01_fastqc/` → … → `07_clusterprofiler/`, execution order reflected in `ls`
- **Multi-profile runtime** — `docker`, `apptainer`, `test` (E2E), `test_stub` (dry-run DAG validation)
- **Pre-built index reuse** — provide `--hisat2_index_prefix` pointing to existing index to skip index building
- **Gene symbol annotation** — DESeq2 output auto-annotated with gene names from the same GTF
- **GO/KEGG enrichment** — clusterProfiler for functional enrichment analysis on DEGs
- **Contrast validation** — pre-flight check catches mismatched group names before execution
- **Schema validation** — `nextflow_schema.json` with `nf-validation` plugin for parameter type/format checks
- **Samplesheet-driven** — single CSV defines samples, groups, and FASTQ paths
- **Reproducible** — pinned tool versions in unified container image
- **Extensible** — add samples or contrasts without touching workflow code

📖 **[Nextflow Learning Guide](LEARNING_NEXTFLOW.md)** — Learn Nextflow core concepts from this project's actual code (in Chinese)

## Quick Start

```bash
# 1. Enter project
cd nextflow-rnaseq

# 2. Prepare samplesheet CSV
cat > samplesheet.csv << EOF
sample,group,fastq_1,fastq_2
WT_1,WT,/data/WT_1_R1.fastq.gz,/data/WT_1_R2.fastq.gz
WT_2,WT,/data/WT_2_R1.fastq.gz,/data/WT_2_R2.fastq.gz
KO_1,KO,/data/KO_1_R1.fastq.gz,/data/KO_1_R2.fastq.gz
KO_2,KO,/data/KO_2_R1.fastq.gz,/data/KO_2_R2.fastq.gz
EOF

# 3. Dry run (stub — no real tools)
nextflow run main.nf -profile test_stub -stub-run \
    --input samplesheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf genes.gtf

# 4. Real run
nextflow run main.nf -profile docker \
    --input samplesheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf genes.gtf
```

## Requirements

| Dependency | Notes |
|---|---|
| [Nextflow](https://www.nextflow.io/) ≥24.0 | Workflow engine |
| [Docker](https://www.docker.com/) ≥20.04 | For Docker mode |
| [Apptainer](https://apptainer.org/) / [Singularity](https://sylabs.io/) ≥1.0 | For HPC mode |
| [Conda](https://docs.conda.io/) / [Mamba](https://mamba.readthedocs.io/) | For local test mode (optional) |

> Build the container image first: `bash containers/build.sh docker`

## Project Structure

```
nextflow-rnaseq/
├── main.nf                        # Main workflow (DSL2)
├── nextflow.config                # Global config — params / container / profiles
├── nextflow_schema.json           # JSON Schema validation (nf-validation plugin)
├── modules/
│   └── local/
│       ├── hisat2_index.nf        #   HISAT2 genome indexing
│       ├── hisat2_align_pe.nf     #   HISAT2 alignment (PE)
│       ├── hisat2_align_se.nf     #   HISAT2 alignment (SE)
│       ├── fastqc_pe.nf           #   FastQC quality control (PE)
│       ├── fastqc_se.nf           #   FastQC quality control (SE)
│       ├── trimmomatic_pe.nf      #   Trimmomatic trimming (PE)
│       ├── trimmomatic_se.nf      #   Trimmomatic trimming (SE)
│       ├── featurecounts.nf       #   featureCounts gene quantification
│       ├── multiqc.nf             #   MultiQC report aggregation
│       ├── deseq2.nf              #   DESeq2 differential expression
│       └── clusterprofiler.nf     #   clusterProfiler GO/KEGG enrichment
├── subworkflows/
│   └── local/
│       ├── qc_trim_align_pe.nf    #   Subworkflow: QC→Trim→QC→Align (PE)
│       └── qc_trim_align_se.nf    #   Subworkflow: QC→Trim→QC→Align (SE)
├── bin/
│   ├── common.R                   #   R shared utility (parse_arg + R_LIBS_ONLY)
│   ├── deseq2.R                   #   DESeq2 + diagnostics
│   ├── gene2symbol.R              #   Gene ID → Symbol annotation
│   └── clusterprofiler.R          #   GO + KEGG enrichment
├── conf/
│   ├── base.config                #   Base resource labels + error strategy
│   ├── docker.config              #   Docker executor
│   ├── apptainer.config           #   Apptainer executor
│   ├── test.config                #   E. coli end-to-end test
│   └── test_stub.config           #   stub-run DAG validation
├── containers/
│   ├── Dockerfile                 #   Unified container image
│   ├── apptainer.def              #   Apptainer definition
│   └── build.sh                   #   Dual-engine build script
├── test/                          # E2E test suite (E. coli)
│   ├── samplesheet.csv            #   4-sample test sheet
│   ├── adapters.fa                #   Empty adapter file
│   ├── ecoli.fa                   #   NCBI E. coli K-12 MG1655 genome
│   ├── ecoli.gtf                  #   NCBI E. coli K-12 MG1655 annotation
│   ├── ecoli_ko.fa                #   400-gene knockout genome (read simulation)
│   └── hisat2_index/              #   Pre-built HISAT2 index (gitignored)
├── .gitignore
├── README.md
└── README_CN.md
```

## Pipeline DAG

```
Layer 0: hisat2_index       One-time genome index (skipped if pre-built index provided)
Layer 1: fastqc_raw         Per-sample raw QC
Layer 2: trimmomatic         Per-sample trim  (PE/SE auto)
Layer 3: fastqc_trimmed     Per-sample post-trim QC
Layer 3: hisat2_align       Per-sample alignment (PE/SE auto)
Layer 4: featurecounts      Aggregate gene counts
Layer 5: multiqc            Aggregate QC report
Layer 6: deseq2              Differential expression + gene symbols
Layer 7: clusterprofiler    GO + KEGG enrichment
```

Layers 1–3 are encapsulated in the `qc_trim_align_{pe,se}` subworkflows. Per-sample processes run in parallel automatically.

## Configuration

[`nextflow.config`](nextflow.config) is the single point of configuration. All parameters settable via `--param value` on the command line or by editing the config file.

### Samplesheet

CSV format (comma-separated; `.tsv` and `.txt` also supported):

```csv
sample,group,fastq_1,fastq_2
WT_1,WT,/data/WT_1_R1.fastq.gz,/data/WT_1_R2.fastq.gz
WT_2,WT,/data/WT_2_R1.fastq.gz,/data/WT_2_R2.fastq.gz
KO_1,KO,/data/KO_1_R1.fastq.gz,/data/KO_1_R2.fastq.gz
KO_2,KO,/data/KO_2_R1.fastq.gz,/data/KO_2_R2.fastq.gz
```

Single-end: omit the `fastq_2` column.

### PE → SE switch

```bash
nextflow run main.nf --read_mode "single" ...
```

No module or subworkflow changes needed.

### Contrasts

Set in `nextflow.config` or via `-c`:

```groovy
params.deseq2.contrasts = [
    [name: "KO_vs_WT", 'case': "KO", control: "WT"]
]
```

Validation: the pipeline checks that `case` and `control` groups exist in the samplesheet before execution. Mismatches produce an immediate, clear error.

### ClusterProfiler (enrichment)

```groovy
params.clusterprofiler {
    org_db        = "org.Hs.eg.db"    // species OrgDb (org.Mm.eg.db for mouse)
    kegg_organism = "hsa"             // KEGG code (mmu, rno, ...)
    from_type     = "ENSEMBL"         // gene ID type → ENTREZID
    pvalue_cutoff = 0.05
    qvalue_cutoff = 0.2
    show_category = 15                // top N terms in dotplots
    gene_id_col   = "gene_id"
}
```

### Pre-built HISAT2 index (skip index building)

```bash
hisat2-build genome.fa /ref/hisat2/genome
nextflow run main.nf \
    --hisat2_index_prefix /ref/hisat2/genome ...
# → HISAT2_INDEX step is skipped entirely
```

## Usage

### Stub mode (DAG validation)

```bash
nextflow run main.nf -profile test_stub -stub-run \
    --input test/samplesheet.csv \
    --genome_fasta test/ecoli.fa \
    --annotation_gtf test/ecoli.gtf
```

### E2E test (E. coli, local conda)

Test data consists of 4 simulated paired-end samples (WT vs KO, 2 replicates
each, 100K read pairs per sample) for E. coli K-12 MG1655. The pre-built
HISAT2 index is detected and index building is skipped automatically.

```bash
# Prerequisites: activate your conda environment with all tools installed,
# and set R_LIBS_ONLY to isolate R packages from system libraries.
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

> `R_LIBS_ONLY` isolates R's library search path to the conda environment,
> preventing binary incompatibility between system R packages and
> conda-installed Bioconductor packages.

### Docker mode (production)

```bash
bash containers/build.sh docker
nextflow run main.nf -profile docker \
    --input samplesheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf genes.gtf
```

### Apptainer mode (HPC)

```bash
bash containers/build.sh all
nextflow run main.nf -profile apptainer \
    --input samplesheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf genes.gtf
```

### Resume interrupted runs

```bash
nextflow run main.nf -profile docker ... -resume
```

## Profiles

| Profile | Use Case | Key Config |
|---------|----------|------------|
| `docker` | Production — Docker containers | `process.container = 'rnaseq-pipeline:1.0.0'` |
| `apptainer` | HPC — Apptainer containers | `apptainer.enabled = true` |
| `test` | Test — E. coli end-to-end | Local conda, NCBI GTF, gene-level counting |
| `test_stub` | Dev — stub-run DAG check | Minimal resources, no container |

## Output Files

### Per-step outputs

| Directory | Key Files |
|---|---|
| `01_fastqc/` | `{read}_fastqc.{html,zip}` (raw + trimmed) |
| `02_trimmomatic/` | `{sample}_R{1,2}.trimmed.fastq.gz` |
| `03_hisat2/` | `{sample}.sorted.bam` + `.bam.bai` + `.hisat2.log` |
| `04_featurecounts/` | `featurecounts.txt`, `featurecounts.summary.txt` |
| `05_multiqc/` | `multiqc_report.html` |

### DESeq2 output (`06_deseq2/`)

| File | Description |
|---|---|
| `{contrast}_all_results.csv` | All genes with statistics + gene_name column |
| `{contrast}_significant.csv` | Significant DEGs + gene_name column |
| `{contrast}_MA_plot.pdf` | MA plot |
| `{contrast}_volcano_plot.pdf` | Volcano plot (labelled with gene symbols when GTF provided) |
| `PCA_plot.pdf` / `.png` | Sample PCA |
| `sample_distance_heatmap.pdf` | Sample-to-sample distance |
| `DEG_heatmap.pdf` | Top N DEG expression heatmap |

### ClusterProfiler output (`07_clusterprofiler/`)

| File | Description |
|---|---|
| `{contrast}_GO_enrichment.csv` | GO enrichment (BP/CC/MF) results |
| `{contrast}_GO_dotplot.pdf` | GO dotplot — top enriched terms |
| `{contrast}_KEGG_enrichment.csv` | KEGG pathway enrichment results |
| `{contrast}_KEGG_dotplot.pdf` | KEGG pathway dotplot |

## Gene Symbol Annotation

The pipeline automatically annotates DEG output with gene symbols using the same GTF file used for counting — ensuring complete consistency across the analysis chain.

```bash
# Also usable standalone:
Rscript bin/gene2symbol.R --input DEG.csv --gtf genes.gtf --output DEG_anno.csv
```

Implementation: `rtracklayer::import()` (not regex). Compatible with both NCBI (`gene` column) and Ensembl (`gene_name` column) GTF formats.

## Container Support

```bash
bash containers/build.sh docker      # Local Docker image
bash containers/build.sh all         # Docker → SIF (recommended for HPC)
bash containers/build.sh test        # Verify all tools
bash containers/build.sh clean       # Remove images
```

Unified image contains: FastQC 0.12.1, Trimmomatic 0.39, HISAT2 2.2.1, Samtools 1.18, featureCounts 2.0.6, MultiQC 1.21, R 4.3.2 + DESeq2 1.42.0 + rtracklayer 1.62.0 + clusterProfiler 4.10.0 + enrichplot 1.22.0.

## Adding Samples or Contrasts

**New sample** — only the samplesheet CSV:

```csv
sample,group,fastq_1,fastq_2
...
New_1,New,/data/New_1_R1.fastq.gz,/data/New_1_R2.fastq.gz   # ← append row
```

**New contrast** — only `nextflow.config` or `-c` override:

```groovy
params.deseq2.contrasts = [
    [name: "KO_vs_WT", 'case': "KO", control: "WT"],
    [name: "New_vs_Ctrl", 'case': "New", control: "Ctrl"]   // ← add
]
```

## Nextflow Features

This pipeline demonstrates:

| Feature | Usage | Notes |
|---|---|---|
| **DSL2 modules** | 11 processes + 2 subworkflows | `include { } from '...'`, clean separation |
| **`.multiMap()`** | Channel fan-out | Broadcasting single CSV parse to 3 consumers |
| **`stub:` blocks** | Dry-run validation | All processes have stubs, `-stub-run` for DAG check |
| **Profile system** | 4 profiles | docker / apptainer / test / test_stub |
| **Schema validation** | `nextflow_schema.json` | nf-validation plugin for param checks |
| **Error strategy** | Exit-code aware retry | OOM (137) and preemption (143) → retry; others → terminate |
| **Groovy in config** | Conditional logic | Pre-built index detection, contrast validation |
| **`publishDir`** | Per-process output | Consistent numbered directory layout |
| **`cache 'deep'`** | Content-addressable | HISAT2 index reuse across runs |
| **PE/SE batching** | Single `--read_mode` switch | All affected modules adapt via `if/else` in main.nf |
| **R script isolation** | `R_LIBS_ONLY` env var | Conda/container library path isolation |
| **Optional outputs** | `optional: true` | clusterProfiler handles empty DEGs gracefully |

## FAQ

**How to add new samples?** — Append rows to your samplesheet CSV.

**How to use single-end reads?** — `nextflow run main.nf --read_mode "single" ...`. No code changes.

**How to skip index building?** — Build index once: `hisat2-build genome.fa /path/genome`, then pass `--hisat2_index_prefix /path/genome`.

**How to clean restart?** — `rm -rf work/ test_results/ && nextflow run main.nf ...`

**How to resume after interruption?** — `nextflow run main.nf -profile docker ... -resume`

**My cluster uses SGE not Slurm?** — Nextflow natively supports SGE, LSF, PBS/Torque, and Slurm via executor configuration.

**How to change cluster resource requests?** — Edit `conf/base.config` process labels.

**How to use a different species for enrichment?** — Set `--clusterprofiler.org_db` and `--clusterprofiler.kegg_organism`. 100+ OrgDb packages available on Bioconductor.

## Tool Versions

| Tool | Version |
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

## License

For educational and research purposes.
