# Advanced RNA-seq Analysis Pipeline
A comprehensive Snakemake workflow for RNA-seq data analysis that combines alignment-based quantification (HISAT2/featureCounts), alignment-free quantification (Salmon), and extensive quality control metrics in a single, easy-to-use pipeline.

## 📚 Table of Contents

- [Overview](#overview)
- [Pipeline Components](#pipeline-components)
- [Workflow Diagram](#workflow-diagram)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Input Requirements](#input-requirements)
- [Output Description](#output-description)
- [Parameters](#parameters)
- [Conda Environments](#conda-environments)
- [License](#license)
- [Contact](#contact)

## 🔍 Overview

This pipeline integrates three complementary workflows for comprehensive RNA-seq analysis: (run them sequentially!)

1. **Core RNA-seq Processing**: Raw read QC, trimming, alignment, and gene-level quantification
2. **Advanced RNA-seq QC**: Multiple QC tools to evaluate alignment quality and RNA integrity
3. **Transcript Quantification**: Alignment-free transcript-level quantification with standard and decoy-aware indices

By combining these approaches, the pipeline delivers robust gene expression estimates while providing detailed quality metrics to ensure reliable results.

## 🧩 Pipeline Components

The workflow consists of three main component pipelines:

### 1. Core RNA-seq Processing

- **Quality Control**: FastQC for raw reads
- **Read Preprocessing**: Fastp for adapter trimming and quality filtering
- **Alignment**: HISAT2 for spliced read alignment to the reference genome (you can also use other aligners)
- **BAM Processing**: Samtools for filtering, sorting, and indexing
- **Expression Quantification**: featureCounts for gene-level counting
- **Results Summary**: MultiQC for aggregated quality reporting

### 2. Advanced RNA-seq QC

- **Fragment Size Analysis**: Picard CollectInsertSizeMetrics
- **Alignment QC**: Qualimap BAM QC for alignment quality metrics
- **RNA-seq Specific QC**: Qualimap RNA-seq for gene model coverage statistics
- **Transcript Analysis**: RSeQC for advanced RNA-seq QC metrics:
  - Read distribution across genomic features
  - GC content analysis
  - Transcript integrity number (TIN) calculation

### 3. Transcript Quantification

- **Alignment-free Quantification**: Salmon for direct transcript quantification
- **Multiple Index Support**:
  - Standard transcriptome index
  - Decoy-aware transcriptome index (reduces mapping bias)

## 📊 Workflow Diagram

```
                         ┌─────────────┐
                         │   Raw FASTQ │
                         └──────┬──────┘
                                │
                         ┌──────▼──────┐
                     ┌───┤   FastQC    ├───┐
                     │   └──────┬──────┘   │
                     │          │          │
                     │   ┌──────▼──────┐   │
                     │   │    Fastp    │   │
                     │   └──────┬──────┘   │
                     │          │          │
          ┌──────────▼──┐       │       ┌──▼───────────┐
          │   HISAT2    │◄──────┴───────►    Salmon    │
          └──────┬──────┘                └──────┬──────┘
                 │                              │
        ┌────────▼─────────┐              ┌─────▼──────┐
        │ Samtools process │              │  Standard  │
        └────────┬─────────┘              │   Index    │
                 │                        └────────────┘
     ┌───────────┼───────────┐                  │
     │           │           │            ┌─────▼──────┐
┌────▼────┐ ┌────▼────┐ ┌────▼─────┐      │Decoy-aware │
│ Picard  │ │Qualimap │ │  RSeQC   │      │   Index    │
└────┬────┘ └────┬────┘ └────┬─────┘      └────────────┘
     │           │           │
     └───────────┼───────────┘
                 │
           ┌─────▼─────┐
           │ featureC. │
           └─────┬─────┘
                 │
           ┌─────▼─────┐
           │  MultiQC  │
           └───────────┘
```

## 🛠️ Requirements

- Snakemake ≥ 7.0
- Conda (for environment management)
- 64GB+ RAM recommended
- 20+ CPU cores recommended for optimal performance

## 📥 Installation

Clone the repository:

```bash
git clone https://github.com/gynecoloji/SnakeMake_RNAseq.git
cd SnakeMake_RNAseq
```

Create conda environments:

```bash
# Main snakemake environment
conda env create -f envs/snakemake.yaml

# QC tools environment
conda env create -f envs/qualimap.yaml

# RSeQC environment
conda env create -f envs/RSeQC.yaml

# Salmon environment
conda env create -f envs/salmon.yaml
```

## 🚀 Usage

1. Place paired-end FASTQ files in the `data/` directory following the naming convention:
   - `{sample}_R1_001.fastq.gz`
   - `{sample}_R2_001.fastq.gz`

2. Configure reference paths in the Snakefile:
   - HISAT2 index: `ref/ENSEMBL/genome`
   - GTF annotation: `ref/Homo_sapiens.GRCh38.102.gtf`
   - Salmon indices: `ref/Salmon_index_Grch38` and `ref/Salmon_index_decoy_Grch38`
   - RSeQC bed file: `ref/ENSEMBL_hg38.bed`

3. Run the workflow:

```bash
# Dry run to verify
snakemake -n -s snakefile_name

# Run with 20 cores
snakemake --use-conda --cores 20 -s snakefile_name -p

# Run on cluster with SLURM (Not tested yet)
snakemake --use-conda --cluster "sbatch --mem=64g --cpus-per-task={threads}" --jobs 20 -s snakefile_name
```

## 📁 Input Requirements

- **FASTQ files**: Paired-end reads with naming pattern `{sample}_R1_001.fastq.gz` and `{sample}_R2_001.fastq.gz`
- **Reference files**:
  - HISAT2 genome index
  - GTF gene annotation file (ENSEMBL format)
  - BED file for RSeQC tools
  - Salmon indices (standard and decoy-aware)

## 📊 Output Description

The pipeline generates organized outputs in the `results/` directory:

```
results/
├── fastqc/                # Raw read quality reports
├── trimmed/               # Trimmed reads and QC reports
├── hisat2/                # Alignment files and summaries
├── samtools/              # Processed BAM files and flagstat reports
├── featurecounts/         # Gene-level counts
├── picard/                # Insert size metrics and plots
├── qualimap_bamqc/        # General alignment QC
├── qualimap_rnaseq/       # RNA-specific QC
├── rseqc/                 # Advanced RNA-seq QC metrics
├── quants/                # Salmon standard index quantification
├── quants_decoy/          # Salmon decoy-aware quantification
└── multiqc_report.html    # Combined QC report
```

## ⚙️ Parameters

Key configurable parameters in the workflow:

- **HISAT2**:
  - Max fragment length: 3000bp
  - Mixed/discordant alignments: disabled

- **Samtools**:
  - Filtering: Properly paired reads (0x2) and primary alignments (-F 0x100)
  - Uniquely mapped: Only reads with NH:i:1 tag

- **featureCounts**:
  - Count paired reads (-p)
  - Feature type: exon
  - ID attribute: gene_id
  - Multi-mapping: count all (-M)
  - Strand-specificity: unstranded (-s 0)

- **Salmon**:
  - Library type: automatic (-l A)
  - Mapping validation: enabled (--validateMappings)

## 🧪 Conda Environments

The workflow uses four Conda environments:

1. **snakemake.yaml**: Core tools (FastQC, fastp, HISAT2, Samtools, featureCounts, MultiQC)
2. **qualimap.yaml**: Qualimap and Picard for alignment QC
3. **RSeQC.yaml**: RSeQC tools for RNA-specific QC
4. **salmon.yaml**: Salmon for transcript quantification

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 📧 Contact

For questions or feedback, please open an issue on the GitHub repository or contact the author.

---

Last updated: April 27, 2025  
Created by: gynecoloji
