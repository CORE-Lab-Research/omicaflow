# OmicaFlow: Modular Multi-Omics Cancer Analysis Pipeline

[![Build Status](https://img.shields.io/badge/Snakemake-7.32.4-brightgreen)](https://snakemake.readthedocs.io)
[![R](https://img.shields.io/badge/R-4.3.1-blue)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.11-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**OmicaFlow** is a modular, reproducible multi-omics pipeline framework for cancer research. It is designed to integrate DNA (SNV/CNV), RNA (differential expression), and methylation data to identify genes with converging alterations across multiple molecular levels.

## Overview

OmicaFlow provides a standardized workflow to address the challenge of **multi-omics data integration**. By automating the process from data acquisition to clinical validation, it enables researchers to systematically identify molecular drivers and assess their prognostic value across different cancer types.

### Key Features

- **Multi-Omics Integration**: Identifies genes with somatic mutations + transcriptional over-expression + promoter hypomethylation.
- **Clinical Validation**: Survival analysis (Kaplan-Meier, Cox proportional hazards) to test prognostic value.
- **Single Configuration**: All parameters adjustable via `config/base.yaml`.
- **Modular Design**: Enable/disable modules independently via config toggles.
- **HPC-Ready**: Pre-configured profiles for SLURM, PBS/Torque, and interactive execution.
- **Reproducible**: Conda environments + renv for complete version pinning.
- **Extensible**: Add new omics modules without modifying existing code.
- **Comprehensive Logging**: Terminal + file logs for all modules with input validation.

## Quick Start

### Prerequisites
- **Conda/Miniconda** installed
- **R 4.3+** (optional if using conda R environment)
- **Git**

### Installation
```bash
# Clone the repository
git clone https://github.com/your-org/omicaflow.git
cd omicaflow

# Create and activate the Snakemake environment
conda env create -f envs/snakemake.yml
conda activate omicaflow-snakemake

# Create and activate the R environment
conda env create -f envs/r_base.yml
```

### Configuration
Edit `config/base.yaml` to adjust parameters:
```yaml
project:
  cancer_type: "BRCA"  # TCGA project code (e.g., BRCA, LUAD, COAD)

modules:
  acquisition: { enabled: true }
  qc: { enabled: true }
  dna: { enabled: true }
  rna: { enabled: true }
  methylation: { enabled: true }
  integration: { enabled: true }
  survival: { enabled: false }
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
```

## Pipeline Modules

| Module | Description | Key Tools | Output |
|--------|-------------|-----------|--------|
| **Acquisition** | TCGA data download via GDC API | TCGAbiolinks | Raw MAF, CNV, RNA, Methylation |
| **QC** | Quality control filtering | R (dplyr, readr) | Filtered data, sample list |
| **DNA Analysis** | SNV/CNV analysis, driver prediction | maftools | Driver genes, mutational burden |
| **RNA Analysis** | Differential expression, pathway enrichment | DESeq2, clusterProfiler | DEGs, enriched pathways |
| **Methylation** | Differential methylation, functional annotation | limma, minfi | DMPs, promoter annotations |
| **Integration** | Multi-omics convergence analysis | R (dplyr) | Genes with 3-omics alterations |
| **Survival** | Kaplan-Meier, Cox PH survival analysis | survival, survminer | Prognostic genes, KM plots |
| **Reporting** | HTML report generation | Jinja2, Python | Summary report |

## Resource Requirements

### Recommended (HPC)
- **Cores**: 32-64
- **Memory**: 64 GB RAM
- **Runtime**: 8-15 hours for full TCGA cohort (~500 samples)
- **Storage**: ~50 GB for cohort data + results

### Minimum (Local)
- **Cores**: 4-8
- **Memory**: 16 GB RAM
- **Runtime**: 24-48 hours for subset (50-100 samples)
- **Storage**: ~10 GB for subset

## Documentation

| Document | Purpose |
|----------|---------|
| [USAGE.md](docs/USAGE.md) | How to run on local/HPC (SLURM/PBS/interactive) |
| [DEVELOPER_GUIDE.md](docs/DEVELOPER_GUIDE.md) | Architecture, algorithms, tools, extensibility |
| [CONTRIBUTING.md](CONTRIBUTING.md) | How to contribute, coding standards |

## Project Structure
```
OmicaFlow/
├── config/              # Configuration files
├── workflow/            # Snakemake orchestration
│   ├── Snakefile        # Main workflow
│   ├── rules/           # Per-module Snakemake rules
│   └── profiles/        # HPC scheduler profiles (SLURM, PBS)
├── modules/             # Analysis scripts (R, Python)
├── templates/           # Jinja2 HTML templates
├── envs/                # Conda environment specifications
├── data/                # Data directory (gitignored except test/)
├── results/             # Analysis outputs (gitignored)
└── docs/                # Documentation
```

## Citations & References

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
- **TCGA**: The Cancer Genome Atlas Research Network. Available at: [https://portal.gdc.cancer.gov](https://portal.gdc.cancer.gov)

## Support

- **Documentation**: See `docs/` directory for comprehensive guides.
- **Issues**: Report bugs or request features via [GitHub Issues](https://github.com/your-org/omicaflow/issues).
- **Questions**: Open a discussion on the repository.

## Acknowledgments

This pipeline was developed to address the need for integrated multi-omics analysis in cancer research. We thank the TCGA Research Network and the Bioconductor community for providing the data and tools that make this work possible.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**Citation**: If you use OmicaFlow in your research, please cite this repository and the key tools listed above.

**Version**: 1.0.0  
**Last Updated**: 2026-05-09
