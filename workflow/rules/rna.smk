# RNA Analysis Module Rule (M03)
rule analyze_rna:
    input:
        filtered_rna="data/qc/{cancer_type}/filtered_RNA_FPKM.tsv",
        sample_list="data/qc/{cancer_type}/final_sample_list.tsv"
    output:
        deg_results="results/rna/{cancer_type}/DEG_results.tsv",
        norm_expr="results/rna/{cancer_type}/Normalized_expression.tsv",
        pathway_enrichment="results/rna/{cancer_type}/Pathway_enrichment.tsv",
        survival_assoc="results/rna/{cancer_type}/Survival_association.tsv",
        rna_qc="results/rna/{cancer_type}/Expression_QC.html"
    log:
        "logs/rna_{cancer_type}.log"
    params:
        deg_tool=config["rna"]["deg_tool"],
        padj_threshold=config["rna"]["padj_threshold"],
        lfc_threshold=config["rna"]["lfc_threshold"],
        enrichment_dbs=config["rna"]["enrichment_databases"]
    conda:
        "../envs/r_base.yml"
    script:
        "modules/rna/analyze_rna.R"