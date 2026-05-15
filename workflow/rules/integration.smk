# Multi-Omics Integration Module Rule (M06)
rule integrate_omics:
    input:
        driver_genes="results/dna/{cancer_type}/Driver_genes.tsv",
        deg_results="results/rna/{cancer_type}/DEG_results.tsv",
        dmp_results="results/methylation/{cancer_type}/DMP_results.tsv"
    output:
        converging_genes="results/integration/{cancer_type}/converging_genes.tsv",
        integration_summary="results/integration/{cancer_type}/integration_summary.tsv",
        venn_data="results/integration/{cancer_type}/venn_data.tsv"
    params:
        maf_threshold=config["dna"]["driver_maf_threshold"],
        deg_lfc_threshold=config["rna"]["lfc_threshold"],
        deg_padj_threshold=config["rna"]["padj_threshold"],
        dmp_logfc_threshold=-0.5,  # Fixed threshold for hypomethylation
        dmp_padj_threshold=config["methylation"]["dmp_padj_threshold"]
    conda:
        "../../envs/r_base.yml"
    script:
        "modules/integration/find_converging_genes.R"