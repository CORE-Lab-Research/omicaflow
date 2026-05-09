# OmicaFlow: Modular Multi-Omics Cancer Analysis Pipeline

[![Build Status](https://img.shields.io/badge/Snakemake-7.32.4-brightgreen)](https://snakemake.readthedocs.io)
[![R](https://img.shields.io/badge/R-4.3.1-blue)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.11-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**OmicaFlow** is a modular, reproducible multi-omics pipeline framework for cancer research. It is designed to integrate DNA (SNV/CNV), RNA (differential expression), and methylation data to identify genes with converging alterations across multiple molecular levels.

## Overview

OmicaFlow provides a standardized workflow to address the challenge of **multi-omics data integration**. By automating the process from data acquisition to clinical validation, it enables researchers to systematically identify molecular drivers and assess their prognostic value across different cancer types.

### Key Features

- **Multi-Omics Integration**: Identifies genes with somatic mutations + transcriptional alterations + promoter hypomethylation.
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
Adjust parameters in `config/base.yaml` to define your project:
```yaml
project:
  cancer_type: "BRCA"  # TCGA project code (e.g., BRCA, LUAD, COAD)

modules:
  acquisition: { enabled: true }
  dna: { enabled: true }
  rna: { enabled: true }
  methylation: { enabled: true }
  integration: { enabled: true }
  survival: { enabled: false }
  reporting: { enabled: true }
```

### Execution
```bash
# Dry-run
snakemake -n

# Run with local cores
snakemake --cores all

# Run on HPC (SLURM example)
snakemake --profile workflow/profiles/slurm
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
- **Storage**: ~50 GB for full cohort data + results

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

## References

Snakemake: Mölder F, Jablonski KP, Letcher B, Hall MB, Tomkins-Tinch CH, Sochat V, Forster J, Lee S, Twardziok SO, Kanitz A, Wilm A, Holtgrewe M, Rahmann S, Nahnsen S, Köster J. Sustainable data analysis with Snakemake. F1000Research. 2021;10:33. doi: 10.12688/f1000research.29032.3

TCGAbiolinks: Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D, Sabedot TS, Malta TM, Pagnotta SM, Castiglioni I, Ceccarelli M, Bontempi G, Noushmehr H. TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. Nucleic Acids Research. 2015;44(8):e71. doi: 10.1093/nar/gkv1507

maftools: Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. Maftools: efficient and comprehensive analysis of somatic variants in cancer. Genome Research. 2018;28(11):1747-1756. doi: 10.1101/gr.239244.118

DESeq2: Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15(12):550. doi: 10.1186/s13059-014-0550-8

limma: Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research. 2015;43(7):e47. doi: 10.1093/nar/gkv007

minfi: Fortin JP, Triche TJ, Hansen KD. Preprocessing, normalization and integration of the Illumina HumanMethylationEPIC array with minfi. Bioinformatics. 2016;33(4):558-560. doi: 10.1093/bioinformatics/btw691

clusterProfiler: Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012;16(5):284-287. doi: 10.1089/omi.2011.0118

survival: Therneau TM, Grambsch PM. Modeling Survival Data: Extending the Cox Model. Springer; 2000.

survminer: Kassambara A, Kosinski M, Biecek P. survminer: Drawing Survival Curves using 'ggplot2'. R package version 0.4.9. 2021.

TCGA Data: The Cancer Genome Atlas Research Network. Available at: https://portal.gdc.cancer.gov

## Support

- **Documentation**: See `docs/` directory for comprehensive guides.
- **Issues**: Report bugs or request features via [GitHub Issues](https://github.com/your-org/omicaflow/issues).
- **Questions**: Open a discussion on the repository.

## Acknowledgments

This pipeline was developed to address the need for integrated multi-omics analysis in cancer research. We thank the TCGA Research Network and the Bioconductor community for providing the data and tools that make this work possible.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
