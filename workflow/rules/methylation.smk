# Methylation Analysis Module Rule (M04)
rule analyze_methylation:
    input:
        filtered_meth="data/qc/{cancer_type}/filtered_Methylation_Beta.tsv",
        sample_list="data/qc/{cancer_type}/final_sample_list.tsv"
    output:
        dmp_results="results/methylation/{cancer_type}/DMP_results.tsv",
        norm_meth="results/methylation/{cancer_type}/Normalized_methylation.tsv",
        func_annotation="results/methylation/{cancer_type}/Functional_annotation.tsv",
        meth_qc="results/methylation/{cancer_type}/Methylation_QC.html"
    log:
        "logs/methylation_{cancer_type}.log"
    params:
        array_type=config["methylation"]["array_type"],
        norm_method=config["methylation"]["normalization_method"],
        padj_threshold=config["methylation"]["dmp_padj_threshold"]
    conda:
        "../../envs/r_base.yml"
    script:
        "../../modules/methylation/analyze_methylation.R"