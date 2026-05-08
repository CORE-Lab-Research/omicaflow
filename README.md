# OmicaFlow: Modular Multi-Omics Cancer Analysis Pipeline

[![Build Status](https://img.shields.io/badge/Snakemake-7.32.4-brightgreen)](https://snakemake.readthedocs.io)
[![R](https://img.shields.io/badge/R-4.3.1-blue)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.11-blue)](https://www.python.org/)

Modular, reproducible multi-omics pipeline for cancer research with MVP implementation focusing on 3-omics: **DNA (SNV/CNV) → RNA (DEG) → Methylation**.

## Features
- **Single Configuration**: All parameters adjustable via `config/base.yaml`
- **Modular Design**: Enable/disable modules via config toggles
- **Snakemake Orchestration**: Reproducible workflow with conda environment management
- **Jinja2 Reporting**: Parameterized HTML reports
- **TCGA Integration**: Direct download via TCGAbiolinks
- **Extensible**: Easily add new omics modules without modifying existing code

## Quick Start

### Prerequisites
- Conda/Miniconda installed
- R 4.3+ with renv (for R package management)
- Snakemake 7.32+ (installed via conda)

### Installation
```bash
# Clone the repository
git clone <repo-url>
cd omicaflow

# Create and activate Snakemake environment
conda env create -f envs/snakemake.yml
conda activate omicaflow-snakemake

# Create and activate R environment
conda env create -f envs/r_base.yml
conda activate omicaflow-r-base

# Initialize R packages via renv
Rscript -e "renv::restore()"
```

### Configuration
Edit `config/base.yaml` to adjust parameters:
```yaml
project:
  cancer_type: "LUAD"  # TCGA project code

modules:
  acquisition: { enabled: true }
  qc: { enabled: true }
  dna: { enabled: true }
  rna: { enabled: true }
  methylation: { enabled: true }
  reporting: { enabled: true }
```

### Run Pipeline
```bash
# Dry-run to check workflow
snakemake -n

# Run the full pipeline
snakemake --cores all

# Run specific module only
snakemake results/dna/LUAD/Driver_genes.tsv
```

## Module Structure
| Module | Description | Dependencies |
|--------|-------------|--------------|
| M00 Acquisition | TCGA data download | None |
| M01 QC | Quality control filtering | M00 |
| M02 DNA Analysis | SNV/CNV analysis with maftools | M01 |
| M03 RNA Analysis | DEG analysis with DESeq2 | M01 |
| M04 Methylation | DMP detection with limma | M01 |
| M05 Reporting | Jinja2 HTML report generation | M02, M03, M04 |

## Project Structure
```
OmicaFlow/
├── config/          # Configuration files
├── workflow/        # Snakemake workflow and rules
├── modules/         # Analysis scripts (R/Python)
├── templates/       # Jinja2 report templates
├── envs/            # Conda environment specs
├── results/         # Generated outputs (gitignored)
└── data/            # Raw/intermediate data (gitignored)
```

## Publication Strategy
- **Methods Paper**: Modular multi-omics pipeline with Snakemake orchestration
- **Application Paper**: Integrated 3-omics signature for cancer prognosis

## Progress Tracking
See [OMICAFLOW_MVP_PLAN.md](OMICAFLOW_MVP_PLAN.md) for detailed progress tracking and implementation plan.

## License
MIT License (to be determined)

## Contact
For questions or contributions, please open an issue on the repository.