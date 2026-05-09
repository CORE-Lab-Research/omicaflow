# OmicaFlow: Modular Multi-Omics Cancer Analysis Pipeline

[![Build Status](https://img.shields.io/badge/Snakemake-7.32.4-brightgreen)](https://snakemake.readthedocs.io)
[![R](https://img.shields.io/badge/R-4.3.1-blue)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.11-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**OmicaFlow** is a modular, reproducible multi-omics pipeline for cancer research, designed to integrate DNA (SNV/CNV), RNA (differential expression), and methylation data to identify genes with converging alterations and assess their prognostic value.

## Overview

OmicaFlow addresses a key challenge in cancer genomics: **identifying genes that are simultaneously altered at multiple molecular levels** (mutation, expression, methylation) and determining whether these converging alterations predict patient outcomes. The pipeline implements a complete workflow from TCGA data acquisition to survival analysis, with all steps configurable via a single YAML file.

### Key Features

- **Multi-Omics Integration**: Identifies genes with somatic mutations + transcriptional over-expression + promoter hypomethylation
- **Clinical Validation**: Survival analysis (Kaplan-Meier, Cox proportional hazards) to test prognostic value
- **Single Configuration**: All parameters adjustable via `config/base.yaml`
- **Modular Design**: Enable/disable modules independently via config toggles
- **HPC-Ready**: Pre-configured profiles for SLURM, PBS/Torque, and interactive execution
- **Reproducible**: Conda environments + renv for complete version pinning
- **Extensible**: Add new omics modules without modifying existing code
- **Comprehensive Logging**: Terminal + file logs for all modules with input validation

## Scientific Rationale

Cancer is driven by alterations at multiple molecular levels. Genes that show **converging alterations** across DNA, RNA, and epigenome are likely to be key drivers with strong functional impact. OmicaFlow systematically identifies such genes and validates their clinical relevance through survival analysis.

**Research Question**: Are there genes that simultaneously experience somatic mutations, transcriptional over-expression, and promoter hypomethylation — and do these genes identify patient subgroups with worse prognosis?

## Quick Start

### Prerequisites
- Conda/Miniconda installed
- R 4.3+ (optional if using conda R environment)
- Git

### Installation
```bash
# Clone the repository
git clone https://github.com/CORE-Lab-Research/omicaflow.git
cd omicaflow

# Create and activate Snakemake environment
conda env create -f envs/snakemake.yml
conda activate omicaflow-snakemake

# Create and activate R environment
conda env create -f envs/r_base.yml
conda activate omicaflow-r-base

# (Optional) Initialize R packages via renv
Rscript -e "renv::restore()"
```

### Configuration
Edit `config/base.yaml` to adjust parameters:
```yaml
project:
  cancer_type: "LUAD"  # TCGA project code (e.g., LUAD, BRCA, COAD)

modules:
  acquisition: { enabled: true }
  qc: { enabled: true }
  dna: { enabled: true }
  rna: { enabled: true }
  methylation: { enabled: true }
  integration: { enabled: true }  # Multi-omics integration
  survival: { enabled: false }    # Requires clinical data
  reporting: { enabled: true }
```

### Run Pipeline
```bash
# Dry-run to validate workflow
snakemake -n

# Run locally with all available cores
snakemake --cores all

# Run on HPC with SLURM
snakemake --profile workflow/profiles/slurm

# Run on HPC with PBS/Torque
snakemake --profile workflow/profiles/pbs

# Run specific module only
snakemake results/integration/LUAD/converging_genes.tsv
```

## Pipeline Modules

| Module | Description | Key Tools | Output |
|--------|-------------|-----------|--------|
| **M00 Acquisition** | TCGA data download via GDC API | TCGAbiolinks | Raw MAF, CNV, RNA, Methylation |
| **M01 QC** | Quality control filtering | R (dplyr, readr) | Filtered data, sample list |
| **M02 DNA Analysis** | SNV/CNV analysis, driver prediction | maftools | Driver genes, mutational burden |
| **M03 RNA Analysis** | Differential expression, pathway enrichment | DESeq2, clusterProfiler | DEGs, enriched pathways |
| **M04 Methylation** | Differential methylation, functional annotation | limma, minfi | DMPs, promoter annotations |
| **M06 Integration** | Multi-omics convergence analysis | R (dplyr) | Genes with 3-omics alterations |
| **M07 Survival** | Kaplan-Meier, Cox PH survival analysis | survival, survminer | Prognostic genes, KM plots |
| **M05 Reporting** | HTML report generation | Jinja2, Python | Summary report |

## Expected Results (TCGA-LUAD Example)

Based on published TCGA studies, typical results for Lung Adenocarcinoma (LUAD):

- **Driver genes**: 50-200 genes (TP53, KRAS, EGFR expected in top 10)
- **DEGs**: 2,000-5,000 genes (padj ≤ 0.05, |LFC| ≥ 1.5)
- **DMPs**: 10,000-50,000 probes (padj ≤ 0.05)
- **Converging genes** (mutation + overexpression + hypomethylation): **5-20 genes**
- **Prognostic genes**: 1-5 genes with significant survival association (HR > 1.5, p < 0.05)

See `docs/EXPECTED_RESULTS.md` for detailed interpretation guide.

## Documentation

| Document | Purpose |
|----------|---------|
| [USAGE.md](docs/USAGE.md) | How to run on local/HPC (SLURM/PBS/interactive) |
| [DEVELOPER_GUIDE.md](docs/DEVELOPER_GUIDE.md) | Architecture, algorithms, tools, extensibility |
| [EXPECTED_RESULTS.md](docs/EXPECTED_RESULTS.md) | What to expect from each module, interpretation |
| [IMPLEMENTATION_SUMMARY.md](docs/IMPLEMENTATION_SUMMARY.md) | Research question, next steps, publication strategy |
| [CONTRIBUTING.md](CONTRIBUTING.md) | How to contribute, coding standards |
| [OMICAFLOW_MVP_PLAN.md](OMICAFLOW_MVP_PLAN.md) | Progress tracking (T001-T027) |

## Project Structure
```
OmicaFlow/
├── config/              # Configuration files
│   └── base.yaml        # Single source of truth for all parameters
├── workflow/            # Snakemake orchestration
│   ├── Snakefile        # Main workflow with conditional module inclusion
│   ├── rules/           # Per-module Snakemake rules
│   └── profiles/        # HPC scheduler profiles (SLURM, PBS)
├── modules/             # Analysis scripts (R/Python)
│   ├── acquisition/     # M00: TCGA data download
│   ├── qc/              # M01: Quality control
│   ├── dna/             # M02: SNV/CNV analysis
│   ├── rna/             # M03: DEG analysis
│   ├── methylation/     # M04: DMP analysis
│   ├── integration/     # M06: Multi-omics integration
│   ├── survival/        # M07: Survival analysis
│   └── reporting/       # M05: Jinja2 report generation
├── templates/           # Jinja2 HTML templates
├── envs/                # Conda environment specifications
├── data/                # Raw/intermediate data (gitignored except test/)
├── results/             # Analysis outputs (gitignored)
├── logs/                # Per-module log files (gitignored)
└── docs/                # Documentation
```

## Resource Requirements

### Recommended (HPC)
- **Cores**: 32-64
- **Memory**: 64 GB RAM
- **Runtime**: 8-15 hours for full TCGA cohort (~500 samples)
- **Storage**: ~50 GB for TCGA-LUAD data + results

### Minimum (Local)
- **Cores**: 4-8
- **Memory**: 16 GB RAM
- **Runtime**: 24-48 hours for subset (50-100 samples)
- **Storage**: ~10 GB for subset

## Publication Strategy

### Target Journals
- **Mid-Impact (IF 5-7)**: *Molecular Cancer*, *Cancers*, *Frontiers in Oncology*
  - Requirements: Clear 3-omics integration + survival validation
- **High-Impact (IF 10+)**: *Nature Communications*, *Genome Medicine*
  - Additional requirements: External validation, novel genes with strong survival signal

### Key Contributions
1. Modular, reproducible multi-omics pipeline (methods paper)
2. Identification of genes with converging 3-omics alterations (application paper)
3. Clinical validation through survival analysis
4. Open-source tool for cancer genomics community

## Citations & Acknowledgments

OmicaFlow builds upon the following key tools and methods:

### Workflow Management
- **Snakemake**: Mölder, F., Jablonski, K.P., Letcher, B., et al. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10:33. DOI: [10.12688/f1000research.29032.3](https://doi.org/10.12688/f1000research.29032.3)

### Data Acquisition
- **TCGAbiolinks**: Colaprico, A., Silva, T.C., Olsen, C., et al. (2015). TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. *Nucleic Acids Research*, 44(8):e71. DOI: [10.1093/nar/gkv1507](https://doi.org/10.1093/nar/gkv1507)

### DNA Analysis
- **maftools**: Mayakonda, A., Lin, D.C., Assenov, Y., Plass, C., & Koeffler, H.P. (2018). Maftools: efficient and comprehensive analysis of somatic variants in cancer. *Genome Research*, 28(11):1747-1756. DOI: [10.1101/gr.239244.118](https://doi.org/10.1101/gr.239244.118)

### RNA Analysis
- **DESeq2**: Love, M.I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12):550. DOI: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

### Methylation & Statistical Analysis
- **limma**: Ritchie, M.E., Phipson, B., Wu, D., et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43(7):e47. DOI: [10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)

### Data Source
- **TCGA**: The Cancer Genome Atlas Research Network. Available at: https://portal.gdc.cancer.gov

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on:
- Code style and standards
- How to add new modules
- Testing requirements
- Pull request process

## Support

- **Documentation**: See `docs/` directory for comprehensive guides
- **Issues**: Report bugs or request features via [GitHub Issues](https://github.com/CORE-Lab-Research/omicaflow/issues)
- **Questions**: Open a discussion on the repository

## Acknowledgments

This pipeline was developed to address the need for integrated multi-omics analysis in cancer research. We thank the TCGA Research Network and the Bioconductor community for providing the data and tools that make this work possible.

---

**Citation**: If you use OmicaFlow in your research, please cite this repository and the key tools listed above.

**Version**: 1.0.0  
**Last Updated**: 2026-05-09