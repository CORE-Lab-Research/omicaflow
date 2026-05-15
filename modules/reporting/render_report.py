#!/usr/bin/env python3
"""
OmicaFlow Reporting Module: Jinja2 Template Rendering
Generates parameterized HTML reports from analysis outputs
"""

import os
import yaml
from jinja2 import Environment, FileSystemLoader
from datetime import datetime
import pandas as pd
import sys

def main():
    # Setup logging using snakemake log if available
    if 'snakemake' in globals():
        log_file = snakemake.log[0]
        cancer_type = snakemake.wildcards.cancer_type
        project_name = snakemake.params.project_name
        template_dir = snakemake.params.template_dir
        
        # Inputs
        dna_summary_path = snakemake.input.dna_summary
        dna_driver_path = snakemake.input.dna_driver_genes
        rna_deg_path = snakemake.input.rna_deg
        meth_dmp_path = snakemake.input.meth_dmp
        conv_genes_path = snakemake.input.converging_genes
        sample_list_path = snakemake.input.sample_list
        
        output_path = snakemake.output.report
    else:
        # Fallback for manual testing
        log_file = "logs/reporting_manual.log"
        cancer_type = "BRCA"
        project_name = "OmicaFlow_Analysis"
        template_dir = "templates/jinja2/"
        dna_summary_path = f"results/dna/{cancer_type}/DNA_integration_summary.tsv"
        dna_driver_path = f"results/dna/{cancer_type}/Driver_genes.tsv"
        rna_deg_path = f"results/rna/{cancer_type}/DEG_results.tsv"
        meth_dmp_path = f"results/methylation/{cancer_type}/DMP_results.tsv"
        conv_genes_path = f"results/integration/{cancer_type}/converging_genes.tsv"
        sample_list_path = f"data/qc/{cancer_type}/final_sample_list.tsv"
        output_path = f"results/reports/{cancer_type}/OmicaFlow_Report.html"

    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    log_fh = open(log_file, "a")
    
    def log(msg):
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{ts}] {msg}"
        print(line)
        log_fh.write(line + "\n")
        log_fh.flush()
    
    log("=== START REPORTING MODULE ===")
    log(f"Project: {project_name} | Cancer: {cancer_type}")
    
    context = {
        "project_name": project_name,
        "cancer_type": cancer_type,
        "generation_date": datetime.now().strftime("%B %d, %Y %H:%M:%S"),
        "timestamp": datetime.now().isoformat(),
    }
    
    # Load DNA results
    log("Loading DNA results...")
    try:
        driver_genes = pd.read_csv(dna_driver_path, sep="\t")
        dna_summary = pd.read_csv(dna_summary_path, sep="\t")
        
        context["dna_total_snv"] = int(dna_summary[dna_summary["metric"] == "Total SNV mutations"]["value"].iloc[0])
        context["dna_driver_genes_count"] = len(driver_genes)
        context["dna_top_genes"] = driver_genes.head(10).to_dict("records")
    except Exception as e:
        log(f"ERROR DNA: {e}")
        context["dna_total_snv"] = 0
        context["dna_driver_genes_count"] = 0
        context["dna_top_genes"] = []
    
    # Load RNA results
    log("Loading RNA results...")
    try:
        deg_results = pd.read_csv(rna_deg_path, sep="\t")
        context["rna_deg_count"] = len(deg_results)
        context["rna_up_count"] = len(deg_results[deg_results["log2FoldChange"] > 0])
        context["rna_down_count"] = len(deg_results[deg_results["log2FoldChange"] < 0])
        context["rna_top_degs"] = deg_results.head(10).to_dict("records")
    except Exception as e:
        log(f"ERROR RNA: {e}")
        context["rna_deg_count"] = 0
        context["rna_up_count"] = 0
        context["rna_down_count"] = 0
        context["rna_top_degs"] = []
    
    # Load Methylation results
    log("Loading Methylation results...")
    try:
        dmp_results = pd.read_csv(meth_dmp_path, sep="\t")
        context["meth_dmp_count"] = len(dmp_results)
    except Exception as e:
        log(f"ERROR Meth: {e}")
        context["meth_dmp_count"] = 0
    
    # Load Converging Genes (Multi-Omics)
    log("Loading Integration results...")
    try:
        conv_genes = pd.read_csv(conv_genes_path, sep="\t")
        context["conv_genes_count"] = len(conv_genes)
        context["conv_genes_list"] = conv_genes.head(20).to_dict("records")
        log(f"Found {context['conv_genes_count']} converging genes")
    except Exception as e:
        log(f"ERROR Integration: {e}")
        context["conv_genes_count"] = 0
        context["conv_genes_list"] = []
    
    # Load sample list
    try:
        sample_list = pd.read_csv(sample_list_path, sep="\t")
        context["total_samples"] = len(sample_list)
    except Exception as e:
        log(f"ERROR Samples: {e}")
        context["total_samples"] = 0
    
    # Render template
    log(f"Rendering template from {template_dir}...")
    try:
        env = Environment(loader=FileSystemLoader(template_dir))
        template = env.get_template("summary_report.html.j2")
        html_output = template.render(**context)
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, "w") as f:
            f.write(html_output)
        log(f"Report generated: {output_path}")
    except Exception as e:
        log(f"ERROR Rendering: {e}")
    
    log("=== REPORTING MODULE COMPLETED ===")
    log_fh.close()
    
    print(f"Report generated: {output_path}")

if __name__ == "__main__":
    main()