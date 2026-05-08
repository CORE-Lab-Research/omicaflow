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
    # Setup logging to file and terminal
    log_file = "logs/reporting_main.log"
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    log_fh = open(log_file, "a")
    
    def log(msg):
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{ts}] {msg}"
        print(line)
        log_fh.write(line + "\n")
        log_fh.flush()
    
    log("=== START REPORTING MODULE ===")
    
    # Load config
    log("Loading config...")
    with open("config/base.yaml") as f:
        config = yaml.safe_load(f)
    
    cancer_type = config["project"]["cancer_type"]
    project_name = config["project"]["name"]
    template_dir = config["reporting"]["template_dir"]
    log(f"Cancer type: {cancer_type}")
    log(f"Project name: {project_name}")
    
    # Prepare context for template
    context = {
        "project_name": project_name,
        "cancer_type": cancer_type,
        "generation_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "timestamp": datetime.now().isoformat(),
    }
    
    # Load DNA results
    dna_dir = f"results/dna/{cancer_type}/"
    log(f"Loading DNA results from {dna_dir}...")
    try:
        driver_genes = pd.read_csv(os.path.join(dna_dir, "Driver_genes.tsv"), sep="\t")
        dna_summary = pd.read_csv(os.path.join(dna_dir, "DNA_integration_summary.tsv"), sep="\t")
        
        context["dna_total_snv"] = int(dna_summary[dna_summary["metric"] == "Total SNV mutations"]["value"].iloc[0])
        context["dna_driver_genes"] = len(driver_genes)
        context["dna_top_genes"] = driver_genes.head(5).to_dict("records")
        log(f"DNA results loaded: {context['dna_total_snv']} SNVs, {context['dna_driver_genes']} driver genes")
    except Exception as e:
        log(f"ERROR loading DNA results: {e}")
        context["dna_total_snv"] = 0
        context["dna_driver_genes"] = 0
        context["dna_top_genes"] = []
    
    # Load RNA results
    rna_dir = f"results/rna/{cancer_type}/"
    log(f"Loading RNA results from {rna_dir}...")
    try:
        deg_results = pd.read_csv(os.path.join(rna_dir, "DEG_results.tsv"), sep="\t")
        context["rna_deg_count"] = len(deg_results)
        context["rna_up_count"] = len(deg_results[deg_results["log2FoldChange"] > 0])
        context["rna_down_count"] = len(deg_results[deg_results["log2FoldChange"] < 0])
        log(f"RNA results loaded: {context['rna_deg_count']} DEGs")
    except Exception as e:
        log(f"ERROR loading RNA results: {e}")
        context["rna_deg_count"] = 0
        context["rna_up_count"] = 0
        context["rna_down_count"] = 0
    
    # Load Methylation results
    meth_dir = f"results/methylation/{cancer_type}/"
    log(f"Loading Methylation results from {meth_dir}...")
    try:
        dmp_results = pd.read_csv(os.path.join(meth_dir, "DMP_results.tsv"), sep="\t")
        context["meth_dmp_count"] = len(dmp_results)
        log(f"Methylation results loaded: {context['meth_dmp_count']} DMPs")
    except Exception as e:
        log(f"ERROR loading Methylation results: {e}")
        context["meth_dmp_count"] = 0
    
    # Load sample list
    try:
        sample_list = pd.read_csv(f"data/qc/{cancer_type}/final_sample_list.tsv", sep="\t")
        context["total_samples"] = len(sample_list)
        log(f"Sample list loaded: {context['total_samples']} samples")
    except Exception as e:
        log(f"ERROR loading sample list: {e}")
        context["total_samples"] = 0
    
    # Render template
    log(f"Rendering template from {template_dir}...")
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template("summary_report.html.j2")
    html_output = template.render(**context)
    
    # Write output
    output_dir = f"results/reports/{cancer_type}/"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "OmicaFlow_Report.html")
    
    with open(output_path, "w") as f:
        f.write(html_output)
    
    log(f"Report generated: {output_path}")
    log("=== REPORTING MODULE COMPLETED ===")
    log_fh.close()
    
    print(f"Report generated: {output_path}")

if __name__ == "__main__":
    main()