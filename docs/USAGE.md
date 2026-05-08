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

### Option A: SLURM Profile (Default)
Run pipeline using pre-configured SLURM profile:
```bash
snakemake --profile workflow/profiles/slurm
```
This uses `workflow/profiles/slurm/config.yaml` which sets default threads, memory, and time limits.

For custom job submission, create a script (e.g., `run_omicaflow.sh`):
```bash
#!/bin/bash
#SBATCH -n 16
#SBATCH --mem 32G
#SBATCH -t 24:00:00
conda activate omicaflow-snakemake
snakemake --profile workflow/profiles/slurm
```
Submit with: `sbatch run_omicaflow.sh`

### Option B: PBS/Torque Profile
If your HPC uses PBS/Torque, use the pre-configured PBS profile:
```bash
snakemake --profile workflow/profiles/pbs
```
This uses `qsub` for job submission. Edit `workflow/profiles/pbs/config.yaml` to adjust resources.

### Option C: Interactive Mode (No Scheduler)
Run Snakemake directly in an interactive shell (most HPCs support `qsub -I` for PBS or `srun -I` for SLURM):
```bash
# Start interactive session (example for PBS)
qsub -I -l nodes=1:ppn=8,mem=32gb,walltime=24:00:00

# Inside interactive session:
conda activate omicaflow-snakemake
snakemake --cores 8  # Use available cores
```

### Option D: Other Schedulers (LSF, SGE)
Create custom profiles in `workflow/profiles/` following the SLURM/PBS examples:
- LSF: Use `bsub` in the `cluster` directive of your profile's config.yaml
- SGE: Use `qsub` with SGE-specific flags

### Option E: Simple Batch Script (Any Scheduler)
Create a generic submission script:
```bash
#!/bin/bash
# For PBS:
#PBS -l nodes=1:ppn=16,mem=32gb,walltime=24:00:00
# For LSF: #BSUB -n 16 -M 32768 -W 24:00

conda activate omicaflow-snakemake
snakemake --cores 16
```
Submit with your scheduler's native command (`qsub`, `bsub`, etc.).

## Logging System
Each module generates logs in `logs/` directory:
- Terminal output: Progress messages print to terminal in real-time
- Log files: Per-module logs saved as `logs/{module}_{cancer_type}.log`
- Example: `logs/dna_LUAD.log` contains all DNA module output and errors

View logs in real-time during run:
```bash
tail -f logs/dna_LUAD.log
```

## Test Dataset (Quick Testing)
A minimal fake dataset is available in `data/test/` for quick pipeline validation without downloading TCGA data:
- `data/test/SNV.maf`
- `data/test/CNV.tsv`
- `data/test/RNA_FPKM.tsv`
- `data/test/Methylation_Beta.tsv`

To use test data, edit `config/base.yaml`:
```yaml
project:
  cancer_type: "test"
modules:
  acquisition: { enabled: false }  # Skip download, use test data
  qc: { enabled: true }
  dna: { enabled: true }
  # ... other modules
```
Then create symlinks or copy test data to expected paths.

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
- View log of specific module:
  ```bash
  cat logs/dna_LUAD.log
  ```

## Configuration
Adjust all parameters in `config/base.yaml`:
- Change cancer type: `project.cancer_type: "LUAD"`
- Toggle modules: `modules.dna.enabled: false`
- Adjust analysis thresholds: `rna.padj_threshold: 0.05`
- Set log directory: `project.log_dir: "logs/"` (default)

## Troubleshooting
- **Rscript not found**: Install R and add to PATH, or use conda R environment
- **Conda environment creation fails**: Try `conda config --set channel_priority strict`
- **Module not running**: Check `config/base.yaml` toggle is `enabled: true`
- **Snakemake rule missing**: Ensure module toggle is enabled in config (rules are conditionally included)
- **Log file not created**: Check `logs/` directory exists and is writable
- **Test data not found**: Ensure `data/test/` files are present (they are committed to repo)