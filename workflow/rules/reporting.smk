# Reporting Module Rule (M05)
rule generate_report:
    input:
        dna_summary="results/dna/{cancer_type}/DNA_integration_summary.tsv",
        rna_deg="results/rna/{cancer_type}/DEG_results.tsv",
        meth_dmp="results/methylation/{cancer_type}/DMP_results.tsv",
        sample_list="data/qc/{cancer_type}/final_sample_list.tsv"
    output:
        report="results/reports/{cancer_type}/OmicaFlow_Report.html"
    log:
        "logs/reporting_{cancer_type}.log"
    params:
        template_dir=config["reporting"]["template_dir"],
        project_name=config["reporting"]["project_name"]
    conda:
        "envs/snakemake.yml"
    script:
        "modules/reporting/render_report.py"