# Survival Analysis Module Rule (M07)
rule analyze_survival:
    input:
        converging_genes="results/integration/{cancer_type}/converging_genes.tsv",
        norm_expr="results/rna/{cancer_type}/Normalized_expression.tsv",
        clinical_data="data/clinical/{cancer_type}/clinical_data.tsv"  # User needs to provide this
    output:
        survival_results="results/survival/{cancer_type}/survival_results.tsv"
    log:
        "logs/survival_{cancer_type}.log"
    params:
        survival_endpoint=config.get("survival", {}).get("endpoint", "OS")  # OS or PFS
    conda:
        "envs/r_base.yml"
    script:
        "modules/survival/analyze_survival.R"