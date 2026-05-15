# QC Module Rule (M01)
rule run_qc:
    input:
        snv="data/raw/{cancer_type}/SNV.maf",
        cnv="data/raw/{cancer_type}/CNV.tsv",
        rna="data/raw/{cancer_type}/RNA_FPKM.tsv",
        methylation="data/raw/{cancer_type}/Methylation_Beta.tsv"
    output:
        filtered_snv="data/qc/{cancer_type}/filtered_SNV.maf",
        filtered_cnv="data/qc/{cancer_type}/filtered_CNV.tsv",
        filtered_rna="data/qc/{cancer_type}/filtered_RNA_FPKM.tsv",
        filtered_meth="data/qc/{cancer_type}/filtered_Methylation_Beta.tsv",
        sample_list="data/qc/{cancer_type}/final_sample_list.tsv",
        qc_report="data/qc/{cancer_type}/QC_report.html"
    log:
        "logs/qc_{cancer_type}.log"
    params:
        min_callrate=config["qc"]["min_sample_callrate"],
        min_map_rate=config["qc"]["min_rna_mapping_rate"],
        max_na=config["qc"]["max_methylation_na"]
    conda:
        "../../envs/r_base.yml"
    script:
        "../modules/qc/run_qc.R"