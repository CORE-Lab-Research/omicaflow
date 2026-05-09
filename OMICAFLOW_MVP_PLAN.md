# OmicaFlow MVP Comprehensive Plan: 3-Omics Multi-Cancer Pipeline
## Version: 1.0 | Status: Implementation Phase | Last Updated: 2026-05-09
*Self-contained plan for implementation and progress tracking. Reference this file if conversation context is lost.*

---

## 1. Executive Summary
Modular, reproducible multi-omics pipeline for cancer research targeting MVP 3-omics: **DNA (SNV/CNV) → RNA (DEG) → Methylation**. Uses Snakemake orchestration, R+Python hybrid stack, single `config/base.yaml` for all adjustable parameters, and Jinja2 for parameterized report generation. Designed for immediate MVP delivery and frictionless future module expansion.

---

## 2. Technology Stack Selection & Rationale
### Orchestration: Snakemake (over Nextflow)
| Feature | Snakemake | Nextflow | Rationale for Choice |
|---------|-----------|----------|----------------------|
| Syntax | Python-native | Groovy-based | Aligns with R+Python hybrid stack, lower learning curve |
| Single config support | Native `config.yaml` | Requires separate config | Meets requirement for all params in one file |
| Jinja2 integration | Native (Python import) | Requires custom plugins | Supports user-requested Jinja2 report templating |
| Conda support | Per-rule native support | Native support | Matches reproducibility requirements |
| Scale | Lab/small cluster | Cloud/HPC enterprise | Sufficient for MVP and mid-term growth |

### Full Tool Stack
| Category | Tools | Purpose |
|----------|-------|---------|
| Orchestration | Snakemake 7.30+ | Pipeline workflow management |
| Reporting | Jinja2 3.1+, Python 3.10+ | Parameterized HTML/PDF report generation |
| R Analysis | R 4.3+, TCGAbiolinks, maftools, DESeq2, minfi, limma | Omics-specific statistical analysis |
| R Environment | renv 1.0+ | R package version pinning |
| Python Utilities | pandas 2.0+, requests 2.31+ | Data parsing, API helpers |
| Conda | conda 23.0+ | Isolated tool environment management |
| QC Aggregation | MultiQC 1.19+ | Cross-module QC report generation |

---

## 3. Core Design Principles
1. **Single Source of Truth**: All adjustable parameters (module toggles, thresholds, tool paths, output dirs) live in one `config/base.yaml` file. No hidden parameters in scripts.
2. **File-Based Module Communication**: Modules exchange data via TSV (tables), RDS (R objects), and HDF5 (large matrices) in standardized directory paths.
3. **Toggle-First Module Activation**: Every module is enabled/disabled via `config.base.yaml` boolean flags. Snakemake skips disabled modules automatically.
4. **Reproducibility**: Per-module conda environments + renv snapshot for R packages. All tool versions pinned.
5. **Extensibility**: New modules only require adding a new directory, Snakemake rule, and config flag. No changes to existing module code.

---

## 4. Single Configuration File Specification
All parameters adjustable from `config/base.yaml` (full example already created in implementation):
- Path: `E:\Project\omicaflow\config\base.yaml`
- Contains all module toggles, tool parameters, thresholds, and environment settings

---

## 5. Pipeline Architecture
### 5.1 High-Level Data Flow
```
config/base.yaml → Acquisition → QC → DNA/RNA/Methylation (parallel) → Reporting → Final Report
```

### 5.2 Module Dependency Graph
| Module ID | Module Name | Dependencies | Output Directory |
|-----------|-------------|--------------|------------------|
| M00 | Acquisition | None | `data/raw/{cancer_type}/` |
| M01 | QC | M00 | `data/qc/{cancer_type}/` |
| M02 | DNA Analysis | M01 | `results/dna/{cancer_type}/` |
| M03 | RNA Analysis | M01 | `results/rna/{cancer_type}/` |
| M04 | Methylation Analysis | M01 | `results/methylation/{cancer_type}/` |
| M05 | Reporting | M02, M03, M04 | `results/reports/{cancer_type}/` |

### 5.3 Directory Structure (Implemented)
```
OmicaFlow/
├── config/
│   └── base.yaml          # Single adjustable config file (CREATED)
├── workflow/
│   ├── Snakefile          # Main Snakemake workflow (CREATED)
│   └── rules/             # Per-module Snakemake rules
│       ├── acquisition.smk (CREATED)
│       ├── qc.smk
│       ├── dna.smk
│       ├── rna.smk
│       ├── methylation.smk
│       └── reporting.smk
├── modules/               # Analysis scripts (R/Python)
│   ├── acquisition/       # (download_tcga.R CREATED)
│   ├── qc/
│   ├── dna/
│   ├── rna/
│   ├── methylation/
│   └── reporting/
├── templates/
│   └── jinja2/            # Jinja2 report templates
├── envs/                  # Conda environment specs (CREATED)
│   ├── snakemake.yml
│   └── r_base.yml
├── results/               # Generated outputs (gitignored)
├── data/                  # Raw/intermediate data (gitignored)
└── docs/                  # Plan, references, user guides
```

---

## 6. MVP Module Details
### 6.1 M00: Acquisition Module (IMPLEMENTED)
- **Toggle**: `modules.acquisition.enabled`
- **Input**: `acquisition.*` params from config
- **Process**: R script using TCGAbiolinks to download matched TCGA samples for DNA/SNV/CNV/RNA/Methylation
- **Output**: Raw MAF, CNV, FPKM, Beta value files + sample barcode mapping table
- **Tools**: TCGAbiolinks, R 4.3+
- **Files Created**:
  - `workflow/rules/acquisition.smk`
  - `modules/acquisition/download_tcga.R`

### 6.2 M01: QC Module (TODO)
- **Toggle**: `modules.qc.enabled`
- **Input**: Raw files from M00 + `qc.*` params
- **Process**: Sample-level filtering, feature-level filtering, MultiQC aggregation
- **Output**: Filtered omics files + QC summary report
- **Tools**: MultiQC, R (qc helper scripts)

### 6.3 M02: DNA Analysis Module (TODO)
- **Toggle**: `modules.dna.enabled`
- **Input**: Filtered MAF/CNV from M01 + `dna.*` params
- **Process**: SNV annotation/driver prediction (maftools), CNV peak calling (GISTIC2)
- **Output**: Annotated SNV table, GISTIC CNV results, driver gene list
- **Tools**: maftools, GISTIC2, R 4.3+

### 6.4 M03: RNA Analysis Module (TODO)
- **Toggle**: `modules.rna.enabled`
- **Input**: Filtered FPKM from M01 + `rna.*` params
- **Process**: Normalization (DESeq2 VST), DEG calling, pathway enrichment (GSEA)
- **Output**: DEG results table, normalized expression matrix, pathway enrichment results
- **Tools**: DESeq2, edgeR, clusterProfiler, R 4.3+

### 6.5 M04: Methylation Analysis Module (TODO)
- **Toggle**: `modules.methylation.enabled`
- **Input**: Filtered Beta values from M01 + `methylation.*` params
- **Process**: Normalization (SWAN), DMP calling (limma), promoter/CGI annotation
- **Output**: DMP results table, normalized methylation matrix, functional annotation
- **Tools**: minfi, limma, R 4.3+

### 6.6 M05: Reporting Module (TODO)
- **Toggle**: `modules.reporting.enabled`
- **Input**: All module outputs + `reporting.*` params
- **Process**: Python script loads Jinja2 templates, passes context from config and module outputs
- **Output**: Parameterized HTML report with embedded figures, tables, and analysis summaries
- **Tools**: Jinja2, Python 3.10+

---

## 7. Orchestration Workflow
### Snakemake Rule Logic for Module Toggles
```python
# Excerpt from workflow/Snakefile (CREATED)
configfile: "config/base.yaml"

# Only include rules if module is enabled
if config["modules"]["acquisition"]["enabled"]:
    include: "workflow/rules/acquisition.smk"
# ... repeat for all modules
```

### Jinja2 Integration Example
```python
# Excerpt from modules/reporting/render_report.py (TODO)
from jinja2 import Environment, FileSystemLoader
import yaml

# Load config
with open("config/base.yaml") as f:
    config = yaml.safe_load(f)

# Load Jinja2 template
env = FileSystemLoader(config["reporting"]["template_dir"])
template = env.get_template("summary_report.html.j2")

# Render with context from module outputs
html_out = template.render(
    cancer_type=config["project"]["cancer_type"],
    deg_count=len(pd.read_csv("results/rna/LUAD/DEG_results.tsv")),
    driver_genes=pd.read_csv("results/dna/LUAD/Driver_genes.tsv")["gene"].tolist(),
    # ... add all context variables
)
```

---

## 8. Extensibility Framework
To add a new module (e.g., Proteomics):
1. Create `modules/proteomics/` with analysis scripts
2. Create `workflow/rules/proteomics.smk` with Snakemake rule
3. Add toggle to `config/base.yaml`: `proteomics: { enabled: false }`
4. Add `proteomics` to `envs/` conda spec
5. Update `workflow/Snakefile` to include rule if enabled
6. (Optional) Add Jinja2 template for proteomics report section

No changes to existing modules required.

---

## 9. Progress Tracking Table
*Update `Status` column as tasks are completed*
| Task ID | Task Description | Module | Priority | Dependencies | Status | Notes |
|---------|------------------|--------|----------|--------------|--------|-------|
| T001 | Initialize repo structure | Infra | High | None | **COMPLETED** | Directories created, git init |
| T002 | Write `config/base.yaml` | Infra | High | T001 | **COMPLETED** | Full single config file |
| T003 | Pin conda envs per module | Infra | High | T001 | **COMPLETED** | snakemake.yml, r_base.yml |
| T004 | Pin R packages via renv | Infra | High | T001 | **SKIPPED** | Rscript not in PATH |
| T005 | Implement M00 Acquisition | M00 | High | T002, T003 | **COMPLETED** | Rule + R script created |
| T006 | Implement M01 QC | M01 | High | T005 | **COMPLETED** | Rule + R script created |
| T007 | Implement M02 DNA Analysis | M02 | High | T006 | **COMPLETED** | Rule + R script created |
| T008 | Implement M03 RNA Analysis | M03 | High | T006 | **COMPLETED** | Rule + R script created |
| T009 | Implement M04 Methylation | M04 | High | T006 | **COMPLETED** | Rule + R script created |
| T010 | Implement M05 Reporting | M05 | Medium | T007, T008, T009 | **COMPLETED** | Jinja2 template + Python script |
| T011 | End-to-end test with LUAD | All | High | T010 | **TODO** | Need R installed to run |
| T012 | Document module interfaces | Docs | Medium | T011 | **TODO** | For future extension |
| T013 | Add MIT License | Docs | Medium | T012 | **COMPLETED** | Committed locally |
| T014 | Add CONTRIBUTING.md | Docs | Medium | T013 | **COMPLETED** | Committed locally |
| T015 | Add usage docs (local + HPC) | Docs | Medium | T014 | **COMPLETED** | Committed locally |
| T016 | Push all local commits to remote | Infra | High | T015 | **TODO** | Push when ready |
| T017 | Add logging (terminal+file) to all modules | All | High | T010 | **COMPLETED** | Logs in logs/, sink() in R, logging in Python |
| T018 | Add input validation to R scripts | All | High | T017 | **COMPLETED** | Check file existence before processing |
| T019 | Create SLURM profile for HPC | Infra | Medium | T015 | **COMPLETED** | workflow/profiles/slurm/config.yaml |
| T020 | Create minimal test dataset | Test | Medium | T019 | **COMPLETED** | data/test/ with fake TCGA-like files |
| T021 | Update USAGE.md with new features | Docs | Medium | T020 | **COMPLETED** | Logging, SLURM, test data sections |
| T022 | Push all enhancements to remote | Infra | High | T021 | **TODO** | Push when ready |
| T023 | Add non-SLURM HPC support (PBS profile, docs) | Infra | Medium | T022 | **COMPLETED** | PBS profile + USAGE update |
| T024 | Implement M06 Multi-Omics Integration | M06 | High | T010, T017, T018 | **COMPLETED** | Rule + R script for converging genes |
| T025 | Implement M07 Survival Analysis | M07 | High | T024 | **COMPLETED** | Rule + R script for survival analysis |
| T026 | Update USAGE.md with new modules | Docs | Medium | T025 | **COMPLETED** | Added integration/survival sections |
| T027 | Push all modules to remote | Infra | High | T026 | **TODO** | Push when ready |
| T023 | Add non-SLURM HPC support (PBS profile, docs) | Infra | Medium | T022 | **COMPLETED** | PBS profile + USAGE update |

---

## 10. Development Roadmap (4-Week Sprint)
| Week | Deliverables |
|------|--------------|
| 1 | Infra setup (T001-T004), M00 Acquisition (T005) |
| 2 | M01 QC (T006), M02 DNA (T007), M03 RNA (T008) |
| 3 | M04 Methylation (T009), M05 Reporting (T010) |
| 4 | End-to-end testing (T011), documentation (T012), MVP validation |

---

## 11. Novelty & Publication Strategy
### MVP Methods Paper (Target: *Bioinformatics*)
- Focus: Modular, configurable 3-omics pipeline with single-config Snakemake orchestration, Jinja2 reporting
- Proof-of-concept: Apply to LUAD TCGA cohort

### Future Application Paper (Target: *Frontiers in Oncology*)
- Focus: Integrated 3-omics signature for LUAD prognosis prediction
- Novelty: Combined SNV driver + DEG + methylation biomarker panel

---

## 12. Appendix
### Key Tool Versions (Pinned)
- Snakemake: 7.32.4 (in snakemake.yml)
- R: 4.3.1 (in r_base.yml)
- Python: 3.11 (in snakemake.yml)
- Jinja2: 3.1.2 (in snakemake.yml)
- DESeq2: 1.40.2 (in r_base.yml)
- maftools: 2.16.0 (in r_base.yml)

### References
- TCGA Data Portal: https://portal.gdc.cancer.gov
- Snakemake Documentation: https://snakemake.readthedocs.io
- TCGAbiolinks: https://bioconductor.org/packages/TCGAbiolinks

---

## 13. Next Steps (Immediate)
1. Install R and add to PATH (to complete T004 renv init)
2. Implement M01 QC module (T006):
   - Create `workflow/rules/qc.smk`
   - Create `modules/qc/qc_processing.R`
   - Add MultiQC integration
3. Continue with M02-M05 modules sequentially

*Last implementation update: 2026-05-09 00:41 (GMT+7)*