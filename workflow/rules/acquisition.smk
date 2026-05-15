# Acquisition Module Rule (M00)
rule download_tcga_data:
    output:
        snv="data/raw/{cancer_type}/SNV.maf",
        cnv="data/raw/{cancer_type}/CNV.tsv",
        rna="data/raw/{cancer_type}/RNA_FPKM.tsv",
        methylation="data/raw/{cancer_type}/Methylation_Beta.tsv",
        sample_map="data/raw/{cancer_type}/sample_map.tsv"
    log:
        "logs/acquisition_{cancer_type}.log"
    params:
        cancer_type="{cancer_type}",
        cache_dir=config["acquisition"]["cache_dir"],
        tcga_url=config["acquisition"]["tcga_api_url"]
    conda:
        "../../envs/r_base.yml"
    script:
        "modules/acquisition/download_tcga.R"