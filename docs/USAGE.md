# OmicaFlow Usage Documentation

## Prerequisites
- Conda/Miniconda installed
- (Optional) R 4.3+ if not using conda R environment
- Git

## Local Setup

### Windows
1. Install Miniconda: https://docs.conda.io/en/latest/miniconda.html
2. Open Anaconda Prompt or PowerShell
3. Clone repo:
   ```bash
   git clone https://github.com/CORE-Lab-Research/omicaflow.git
   cd omicaflow
   ```
4. Create and activate Snakemake environment:
   ```bash
   conda env create -f envs/snakemake.yml
   conda activate omicaflow-snakemake
   ```
5. Create and activate R environment:
   ```bash
   conda env create -f envs/r_base.yml
   conda activate omicaflow-r-base
   ```
6. (Optional) Initialize R packages via renv:
   ```bash
   Rscript -e "renv::restore()"
   ```

### macOS/Linux
Same as Windows but use terminal.

## HPC/SSH Setup
1. SSH into HPC cluster:
   ```bash
   ssh username@hpc-cluster.edu
   ```
2. Load conda module (if required by cluster):
   ```bash
   module load conda
   ```
3. Clone repo and set up environments as in Local Setup steps 3-5
4. Create a job submission script (e.g., `run_omicaflow.sh` for SLURM):
   ```bash
   #!/bin/bash
   #SBATCH -n 16
   #SBATCH --mem 32G
   #SBATCH -t 24:00:00
   conda activate omicaflow-snakemake
   snakemake --cores 16
   ```
5. Submit job:
   ```bash
   sbatch run_omicaflow.sh
   ```

## Common Commands
- Dry-run (validate workflow without execution):
  ```bash
  snakemake -n
  ```
- Run all enabled modules:
  ```bash
  snakemake --cores all
  ```
- Run specific module output (e.g., DNA analysis results):
  ```bash
  snakemake results/dna/LUAD/Driver_genes.tsv
  ```
- Disable a module: Edit `config/base.yaml` set `enabled: false` for that module

## Configuration
Adjust all parameters in `config/base.yaml`:
- Change cancer type: `project.cancer_type: "LUAD"`
- Toggle modules: `modules.dna.enabled: false`
- Adjust analysis thresholds: `rna.padj_threshold: 0.05`

## Troubleshooting
- **Rscript not found**: Install R and add to PATH, or use conda R environment
- **Conda environment creation fails**: Try `conda config --set channel_priority strict`
- **Module not running**: Check `config/base.yaml` toggle is `enabled: true`
- **Snakemake rule missing**: Ensure module toggle is enabled in config (rules are conditionally included)