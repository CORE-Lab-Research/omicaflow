# OmicaFlow Developer Guide

## Table of Contents
1. [Architecture Overview](#architecture-overview)
2. [Tools & Technologies](#tools--technologies)
3. [Module Design & Logic](#module-design--logic)
4. [Algorithms & Statistical Methods](#algorithms--statistical-methods)
5. [Data Flow & File Formats](#data-flow--file-formats)
6. [Extending the Pipeline](#extending-the-pipeline)
7. [Testing & Validation](#testing--validation)

---

## Architecture Overview

### High-Level Architecture
OmicaFlow follows a **modular, workflow-based architecture** using Snakemake for orchestration:

```
┌─────────────────────────────────────────────────────────────┐
│                    config/base.yaml                         │
│              (Single Source of Truth)                       │
└─────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────┐
│                  workflow/Snakefile                         │
│         (Conditional Module Inclusion)                      │
└─────────────────────────────────────────────────────────────┘
                            │
        ┌───────────────────┼───────────────────┐
        ▼                   ▼                   ▼
    ┌───────┐          ┌───────┐          ┌───────┐
    │  M00  │──────────│  M01  │──────────│  M02  │
    │ Acq   │          │  QC   │          │  DNA  │
    └───────┘          └───────┘          └───────┘
                            │
        ┌───────────────────┼───────────────────┐
        ▼                   ▼                   ▼
    ┌───────┐          ┌───────┐          ┌───────┐
    │  M03  │          │  M04  │          │  M06  │
    │  RNA  │          │ Meth  │          │ Integ │
    └───────┘          └───────┘          └───────┘
                            │
        ┌───────────────────┼───────────────────┐
        ▼                   ▼                   ▼
    ┌───────┐          ┌───────┐          ┌───────┐
    │  M07  │          │  M05  │          │ Logs  │
    │ Surv  │          │Report │          │       │
    └───────┘          └───────┘          └───────┘
```

### Design Principles
1. **Modularity**: Each module is independent, with clear input/output contracts
2. **Configurability**: All parameters adjustable via single YAML config
3. **Reproducibility**: Conda environments + renv for version pinning
4. **Extensibility**: Add new modules without modifying existing code
5. **Observability**: Comprehensive logging (terminal + file) for all modules

### Directory Structure
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

---

## Tools & Technologies

### Orchestration Layer
**Snakemake 7.32.4**
- **Purpose**: Workflow management, dependency resolution, parallel execution
- **Why Snakemake**: Python-native, conda integration, HPC scheduler support
- **Key Features Used**:
  - Conditional rule inclusion (module toggles)
  - Wildcard-based file patterns (`{cancer_type}`)
  - Per-rule conda environments
  - Log file directives
  - Cluster profiles (SLURM, PBS)

### R Ecosystem (Statistical Analysis)
**R 4.3.1** with Bioconductor packages:

| Package | Version | Purpose | Module |
|---------|---------|---------|--------|
| TCGAbiolinks | 2.30.0 | TCGA data download via GDC API | M00 |
| maftools | 2.16.0 | SNV annotation, driver prediction | M02 |
| DESeq2 | 1.40.2 | Differential expression analysis | M03 |
| edgeR | 3.42.4 | Alternative DEG method | M03 |
| limma | 3.56.2 | DMP detection (methylation) | M04 |
| minfi | 1.46.0 | Methylation array processing | M04 |
| clusterProfiler | 4.8.0 | Pathway enrichment (GO/KEGG) | M03 |
| survival | 3.5.5 | Kaplan-Meier, Cox PH models | M07 |
| survminer | 0.4.9 | Survival plot generation | M07 |

**Why R for Omics**: Bioconductor ecosystem is the gold standard for genomics analysis.

### Python Ecosystem (Utilities & Reporting)
**Python 3.11** with:

| Package | Version | Purpose | Module |
|---------|---------|---------|--------|
| Jinja2 | 3.1.2 | HTML report templating | M05 |
| pandas | 2.0 | Data manipulation | M05 |
| PyYAML | 6.0 | Config parsing | All |

### Environment Management
- **Conda**: Tool-level isolation (Snakemake, R, Python)
- **renv**: R package version pinning (optional, for reproducibility)

### Version Control & Collaboration
- **Git**: Version control
- **GitHub**: Remote repository hosting

---

## Module Design & Logic

### Module Interface Contract
Each module follows a standard interface:

```yaml
# Snakemake Rule Structure
rule module_name:
    input:
        file1="path/to/input1.tsv"
        file2="path/to/input2.tsv"
    output:
        result="path/to/output.tsv"
    log:
        "logs/module_{cancer_type}.log"
    params:
        threshold=config["module"]["threshold"]
    conda:
        "envs/r_base.yml"
    script:
        "modules/module/script.R"
```

**R Script Structure**:
```r
# 1. Load libraries
library(package)

# 2. Get Snakemake parameters
input_file <- snakemake@input$file1
param <- snakemake@params$threshold
log_file <- snakemake@log[[1]]

# 3. Setup logging
sink(log_file, type="output", append=TRUE)
sink(log_file, type="message", append=TRUE)
message("=== START MODULE: ", Sys.time(), " ===")

# 4. Input validation
if (!file.exists(input_file)) stop("Missing input")

# 5. Analysis logic
result <- analyze(input_file, param)

# 6. Save output
write_tsv(result, snakemake@output$result)

# 7. Close logging
message("=== MODULE COMPLETED: ", Sys.time(), " ===")
sink(type="output")
sink(type="message")
```

### Module-Specific Logic

#### M00: Acquisition (TCGA Data Download)
**Logic**:
1. Query TCGA GDC API for specified cancer type
2. Download 4 data types: SNV (MAF), CNV (segments), RNA (FPKM), Methylation (beta values)
3. Use TCGAbiolinks for API interaction and data preparation
4. Create sample barcode mapping table

**Key Functions**:
- `GDCquery()`: Query TCGA data catalog
- `GDCdownload()`: Download files via API
- `GDCprepare()`: Parse and format data

**Output Format**:
- SNV: MAF format (Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification, etc.)
- CNV: TSV (Sample, Chromosome, Start, End, segment_mean, gene_symbol)
- RNA: TSV (gene_id as rows, samples as columns, FPKM values)
- Methylation: TSV (probe_id, gene_symbol, samples as columns, beta values 0-1)

#### M01: QC (Quality Control)
**Logic**:
1. **Sample-level filtering**:
   - SNV: Call rate ≥ 95% (remove low-quality samples)
   - RNA: Missing values ≤ 20% per sample
   - CNV: Minimum 5 segments per sample
2. **Feature-level filtering**:
   - Methylation: Probes with NA ≤ 10%
3. **Cross-omics matching**: Keep only samples present in all omics layers

**Thresholds** (adjustable in config):
- `min_sample_callrate`: 0.95
- `min_rna_mapping_rate`: 0.8
- `max_methylation_na`: 0.1

**Output**: Filtered data files + final sample list (intersection of all omics)

#### M02: DNA Analysis (SNV/CNV)
**Logic**:
1. **SNV Analysis**:
   - Load MAF file with maftools
   - Annotate variants (VEP/ANNOVAR if available)
   - Identify driver genes: MAF ≥ 5% (mutated in ≥5% of samples)
   - Calculate mutational burden per sample
2. **CNV Analysis**:
   - Identify amplifications (segment_mean > 0.5) and deletions (segment_mean < -0.5)
   - Count events per gene

**Key Algorithms**:
- **Driver prediction**: Frequency-based (MAF threshold)
- **Mutational burden**: Total mutations per sample

**Output**:
- Driver genes: Hugo_Symbol, n_samples, total_mutations, MAF
- Mutational burden: Sample, total_mutations, missense, nonsense, frameshift

#### M03: RNA Analysis (DEG)
**Logic**:
1. **Normalization**: DESeq2 Variance Stabilizing Transformation (VST)
2. **Differential Expression**:
   - Model: `~ condition` (Tumor vs Normal)
   - Test: Wald test (DESeq2 default)
   - Shrinkage: apeglm for log2FoldChange
3. **Filtering**: padj ≤ 0.05, |log2FoldChange| ≥ 1.5
4. **Pathway Enrichment**: GSEA with GO/KEGG databases

**Key Algorithms**:
- **DESeq2**: Negative binomial GLM for count data
- **GSEA**: Gene Set Enrichment Analysis (hypergeometric test)

**Output**:
- DEG results: gene_id, log2FoldChange, padj
- Normalized expression matrix
- Enriched pathways: pathway_id, p_value, genes

#### M04: Methylation Analysis (DMP)
**Logic**:
1. **Normalization**: SWAN (Subset-quantile Within Array Normalization)
2. **Differential Methylation**:
   - Model: `~ condition` (Tumor vs Normal)
   - Test: limma (linear model + empirical Bayes)
3. **Filtering**: adj.P.Val ≤ 0.05
4. **Functional Annotation**: Classify probes by genomic location (TSS, promoter, CGI)

**Key Algorithms**:
- **limma**: Linear model with empirical Bayes moderation
- **SWAN**: Normalization for Illumina arrays

**Output**:
- DMP results: probe_id, gene_symbol, logFC, adj.P.Val
- Functional annotation: location (TSS/Exon/Intergenic), cpg_island (TRUE/FALSE)

#### M06: Multi-Omics Integration (Converging Genes)
**Logic**:
1. **Define criteria**:
   - DNA: Driver genes (MAF ≥ 5%)
   - RNA: Over-expressed (log2FoldChange ≥ 1.5, padj ≤ 0.05)
   - Methylation: Hypomethylated at promoter (logFC ≤ -0.5, padj ≤ 0.05)
2. **Find intersection**: Genes meeting ALL THREE criteria
3. **Rank by effect size**: Sort by MAF, then log2FoldChange

**Key Algorithm**:
- **Set intersection**: Inner join on gene symbol across 3 omics layers

**Output**:
- Converging genes: gene, maf, log2FoldChange, logFC (methylation)
- Integration summary: Counts at each level (DNA only, DNA+RNA, DNA+RNA+Methylation)
- Venn diagram data: Binary matrix (in_driver, in_deg, in_dmp)

#### M07: Survival Analysis (Prognosis)
**Logic**:
1. **Stratification**: For each converging gene, split patients by median expression (High vs Low)
2. **Kaplan-Meier**: Estimate survival curves for High vs Low groups
3. **Log-rank test**: Test for significant survival difference
4. **Cox Proportional Hazards**: Calculate Hazard Ratio (HR) and 95% CI
5. **Multiple testing correction**: Benjamini-Hochberg FDR

**Key Algorithms**:
- **Kaplan-Meier estimator**: Non-parametric survival curve
- **Log-rank test**: Compare survival distributions
- **Cox PH model**: `Surv(time, status) ~ risk_group`

**Output**:
- Survival results: gene, n_samples, n_events, HR, HR_lower, HR_upper, p_value, p_adj
- Kaplan-Meier plots (PDF) for top 3 significant genes

---

## Algorithms & Statistical Methods

### Statistical Tests Used

| Module | Test | Purpose | Null Hypothesis | Alternative |
|--------|------|---------|-----------------|-------------|
| M02 DNA | Frequency-based | Driver prediction | Gene mutated by chance | Gene is driver |
| M03 RNA | Wald test (DESeq2) | Differential expression | No expression change | Expression differs |
| M04 Methylation | Moderated t-test (limma) | Differential methylation | No methylation change | Methylation differs |
| M06 Integration | Set intersection | Find converging genes | N/A | N/A |
| M07 Survival | Log-rank test | Survival difference | No survival difference | Survival differs |
| M07 Survival | Cox PH | Hazard ratio | HR = 1 | HR ≠ 1 |

### Multiple Testing Correction
All modules use **Benjamini-Hochberg FDR** (False Discovery Rate) to control for multiple comparisons:
- M03 RNA: `p.adjust(method="BH")` on DEG p-values
- M04 Methylation: `p.adjust(method="BH")` on DMP p-values
- M07 Survival: `p.adjust(method="BH")` on log-rank p-values

**Why FDR**: More powerful than Bonferroni for exploratory omics studies, controls expected proportion of false positives.

### Normalization Methods

#### RNA-seq (M03)
**DESeq2 VST (Variance Stabilizing Transformation)**:
- Transforms count data to log2-like scale
- Stabilizes variance across expression range
- Accounts for library size differences

**Formula**: `vst(counts) = log2(counts + pseudocount) - size_factor_adjustment`

#### Methylation (M04)
**SWAN (Subset-quantile Within Array Normalization)**:
- Corrects for probe type bias (Type I vs Type II probes on Illumina arrays)
- Quantile normalization within probe type subsets

**Why SWAN**: Illumina 450K/EPIC arrays have two probe types with different intensity distributions.

### Thresholds & Rationale

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| DEG padj | 0.05 | Standard FDR cutoff (5% false positives) |
| DEG LFC | 1.5 | Biological significance (≥2.8-fold change) |
| DMP padj | 0.05 | Standard FDR cutoff |
| DMP logFC | -0.5 | Hypomethylation (≥30% decrease in beta value) |
| Driver MAF | 0.05 | Mutated in ≥5% of cohort (recurrent mutations) |
| Survival HR | 1.5 | Clinically meaningful (50% increased risk) |

---

## Data Flow & File Formats

### Data Flow Diagram
```
TCGA GDC API
     │
     ▼
[M00] Raw Data (MAF, CNV, FPKM, Beta)
     │
     ▼
[M01] Filtered Data + Sample List
     │
     ├──────────┬──────────┐
     ▼          ▼          ▼
  [M02]      [M03]      [M04]
   DNA        RNA        Meth
     │          │          │
     └──────────┼──────────┘
                ▼
            [M06] Converging Genes
                │
                ├──────────┐
                ▼          ▼
            [M07]      [M05]
            Survival   Report
```

### File Format Specifications

#### TSV (Tab-Separated Values)
**Standard format for all tabular data**:
- Delimiter: `\t` (tab)
- Header: First row contains column names
- Encoding: UTF-8
- Missing values: `NA`

**Example (DEG results)**:
```
gene_id	log2FoldChange	padj
ENSG00000141510	2.5	0.001
ENSG00000171862	-1.8	0.023
```

#### MAF (Mutation Annotation Format)
**Standard format for SNV data**:
- Required columns: Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification
- Format: TSV with specific column names
- Used by maftools

#### RDS (R Data Serialization)
**For complex R objects** (e.g., DESeqDataSet, SummarizedExperiment):
- Binary format
- Preserves R object structure
- Not human-readable (use for intermediate storage only)

### Sample Barcode Format (TCGA)
**Format**: `TCGA-XX-XXXX-XXA-XXX-XXXX-XX`
- First 12 characters: Patient ID (`TCGA-XX-XXXX`)
- Characters 14-15: Sample type (`01`=Primary Tumor, `11`=Normal)

**Example**: `TCGA-05-4250-01A-01D-2088-08`

---

## Extending the Pipeline

### Adding a New Module (Example: Proteomics)

#### Step 1: Create Module Directory
```bash
mkdir -p modules/proteomics
```

#### Step 2: Write Analysis Script
Create `modules/proteomics/analyze_proteomics.R`:
```r
# M08 Proteomics Module
library(dplyr)
library(readr)

# Get Snakemake parameters
input_file <- snakemake@input$proteomics_data
log_file <- snakemake@log[[1]]

# Setup logging
sink(log_file, type="output", append=TRUE)
message("=== START PROTEOMICS MODULE ===")

# Input validation
if (!file.exists(input_file)) stop("Missing input")

# Analysis logic
result <- analyze_proteomics(input_file)

# Save output
write_tsv(result, snakemake@output$result)

message("=== PROTEOMICS MODULE COMPLETED ===")
sink(type="output")
```

#### Step 3: Create Snakemake Rule
Create `workflow/rules/proteomics.smk`:
```python
rule analyze_proteomics:
    input:
        proteomics_data="data/proteomics/{cancer_type}/protein_abundance.tsv"
    output:
        result="results/proteomics/{cancer_type}/DEP_results.tsv"
    log:
        "logs/proteomics_{cancer_type}.log"
    params:
        padj_threshold=config["proteomics"]["padj_threshold"]
    conda:
        "envs/r_base.yml"
    script:
        "modules/proteomics/analyze_proteomics.R"
```

#### Step 4: Update Config
Add to `config/base.yaml`:
```yaml
modules:
  proteomics: { enabled: false }

proteomics:
  padj_threshold: 0.05
  lfc_threshold: 1.5
```

#### Step 5: Update Snakefile
Add to `workflow/Snakefile`:
```python
if config["modules"]["proteomics"]["enabled"]:
    include: "rules/proteomics.smk"

# Update rule all
rule all:
    input:
        proteomics=(
            f"results/proteomics/{config['project']['cancer_type']}/DEP_results.tsv"
            if config["modules"]["proteomics"]["enabled"]
            else []
        )
```

#### Step 6: Test
```bash
# Enable module
# Edit config/base.yaml: proteomics: { enabled: true }

# Dry-run
snakemake -n

# Run
snakemake --cores all
```

---

## Testing & Validation

### Unit Testing (Per Module)
**Test with minimal fake data**:
```bash
# Use test dataset
snakemake results/dna/test/Driver_genes.tsv --cores 1
```

**Expected behavior**:
- Module completes without errors
- Output files created
- Log file shows progress messages

### Integration Testing (Full Pipeline)
**Test with TCGA subset** (50-100 samples):
```bash
# Edit config to use subset
snakemake --cores all
```

**Validation checks**:
1. All output files created
2. Known biomarkers present (e.g., TP53)
3. Result counts in expected range

### Reproducibility Testing
**Run pipeline twice, compare results**:
```bash
# Run 1
snakemake --cores all

# Backup results
cp -r results results_run1

# Clean
snakemake --delete-all-output

# Run 2
snakemake --cores all

# Compare
diff -r results results_run1
```

**Expected**: Identical results (deterministic algorithms)

### Performance Benchmarking
**Add benchmark directive to rules**:
```python
rule analyze_rna:
    benchmark:
        "benchmarks/rna_{cancer_type}.txt"
```

**Output**: Runtime, memory usage per rule

---

## Troubleshooting

### Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| `Rscript not found` | R not in PATH | Use conda R env: `conda activate omicaflow-r-base` |
| `No DEGs found` | Thresholds too strict | Relax padj/LFC in config |
| `Memory error (DESeq2)` | Large RNA matrix | Increase memory allocation or use subset |
| `Missing input file` | Module disabled or failed | Check logs/, enable prerequisite modules |
| `Snakemake rule not found` | Module not enabled in config | Set `enabled: true` in config/base.yaml |

### Debug Mode
**Enable verbose logging**:
```bash
snakemake --cores all --verbose --printshellcmds
```

**Check module logs**:
```bash
tail -f logs/dna_BRCA.log
```

---

## Performance Optimization

### Parallelization
- **Snakemake**: Automatically parallelizes independent rules
- **Optimal cores**: 16-32 for TCGA cohorts (500 samples)
- **Bottleneck**: M03 RNA (DESeq2 is single-threaded)

### Memory Management
- **Peak usage**: M03 RNA (16-32 GB for 500 samples)
- **Optimization**: Use `--resources mem_mb=32000` to limit memory

### HPC Best Practices
- Use cluster profiles (`--profile workflow/profiles/slurm`)
- Set appropriate walltime (12-24 hours for full cohort)
- Monitor with `squeue` (SLURM) or `qstat` (PBS)

---

## References

### Key Publications
- **DESeq2**: Love et al. (2014) Genome Biology
- **limma**: Ritchie et al. (2015) Nucleic Acids Research
- **maftools**: Mayakonda et al. (2018) Genome Research
- **TCGAbiolinks**: Colaprico et al. (2016) Nucleic Acids Research

### External Resources
- TCGA Data Portal: https://portal.gdc.cancer.gov
- Snakemake Documentation: https://snakemake.readthedocs.io
- Bioconductor: https://bioconductor.org

---

*Last updated: 2026-05-09*