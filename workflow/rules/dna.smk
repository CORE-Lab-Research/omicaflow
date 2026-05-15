# DNA Analysis Module Rule (M02)
rule analyze_dna:
    input:
        filtered_snv="data/qc/{cancer_type}/filtered_SNV.maf",
        filtered_cnv="data/qc/{cancer_type}/filtered_CNV.tsv",
        sample_list="data/qc/{cancer_type}/final_sample_list.tsv"
    output:
        annotated_snv="results/dna/{cancer_type}/SNV_annotated.tsv",
        driver_genes="results/dna/{cancer_type}/Driver_genes.tsv",
        cnv_results="results/dna/{cancer_type}/CNV_GISTIC.tsv",
        mutational_burden="results/dna/{cancer_type}/Mutational_burden.tsv",
        dna_summary="results/dna/{cancer_type}/DNA_integration_summary.tsv"
    log:
        "logs/dna_{cancer_type}.log"
    params:
        snv_tool=config["dna"]["snv_tool"],
        cnv_tool=config["dna"]["cnv_tool"],
        driver_prediction=config["dna"]["driver_prediction"],
        maf_threshold=config["dna"]["driver_maf_threshold"]
    conda:
        "../../envs/r_base.yml"
    script:
        "../modules/dna/analyze_dna.R"