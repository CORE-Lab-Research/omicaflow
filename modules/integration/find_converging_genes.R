# M06 Multi-Omics Integration Module: Find Converging Genes
# Identifies genes with somatic mutation + over-expression + promoter hypomethylation

library(dplyr)
library(readr)
library(tidyr)

# Get parameters from Snakemake
driver_genes <- snakemake@input$driver_genes
deg_results <- snakemake@input$deg_results
dmp_results <- snakemake@input$dmp_results
log_file <- snakemake@log[[1]]

# Setup logging
log_dir <- dirname(log_file)
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
log_con <- file(log_file, open = "a")
sink(log_con, type = "output")
sink(log_con, type = "message")
message("=== START MULTI-OMICS INTEGRATION MODULE: ", Sys.time(), " ===")

# Input validation
input_files <- c(driver_genes, deg_results, dmp_results)
for (f in input_files) {
    if (!file.exists(f)) {
        message("ERROR: Input file not found: ", f)
        stop(paste("Missing input file:", f))
    }
}
message("All input files validated")

# Create output directory
output_dir <- dirname(snakemake@output$converging_genes)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
}

# --- Load Data ---
message("Loading driver genes...")
driver_df <- read_tsv(driver_genes, show_col_types = FALSE) %>%
    mutate(gene = Hugo_Symbol) %>%
    select(gene, n_samples, maf)
message("Driver genes loaded: ", nrow(driver_df), " genes")

message("Loading DEG results...")
deg_df <- read_tsv(deg_results, show_col_types = FALSE) %>%
    mutate(gene = gene_id) %>%
    select(gene, log2FoldChange, padj)
message("DEG results loaded: ", nrow(deg_df), " genes")

message("Loading DMP results...")
dmp_df <- read_tsv(dmp_results, show_col_types = FALSE) %>%
    mutate(gene = gene_symbol) %>%
    select(gene, logFC, adj.P.Val)
message("DMP results loaded: ", nrow(dmp_df), " genes")

# --- Define Thresholds (from Snakemake params) ---
maf_threshold <- snakemake@params$maf_threshold
deg_lfc_threshold <- snakemake@params$deg_lfc_threshold
deg_padj_threshold <- snakemake@params$deg_padj_threshold
dmp_logfc_threshold <- snakemake@params$dmp_logfc_threshold
dmp_padj_threshold <- snakemake@params$dmp_padj_threshold

message("Using thresholds from config:")
message("  - Driver genes: MAF >= ", maf_threshold)
message("  - DEG: LFC >= ", deg_lfc_threshold, " (over-expression), padj <= ", deg_padj_threshold)
message("  - DMP: logFC <= ", dmp_logfc_threshold, " (hypomethylation), padj <= ", dmp_padj_threshold)

# --- Find Converging Genes ---
message("Finding genes with mutation + over-expression...")
mut_expr_genes <- driver_df %>%
    filter(maf >= maf_threshold) %>%
    inner_join(deg_df %>%
                  filter(log2FoldChange >= deg_lfc_threshold, padj <= deg_padj_threshold),
              by = "gene")
message("Mutation + expression genes: ", nrow(mut_expr_genes))

# Add methylation (hypomethylation) filter
message("Adding methylation (hypomethylation) filter...")
converging_genes <- mut_expr_genes %>%
    inner_join(dmp_df %>%
                  filter(logFC <= dmp_logfc_threshold, adj.P.Val <= dmp_padj_threshold),
              by = "gene") %>%
    arrange(desc(maf), desc(log2FoldChange), logFC)  # Sort by mutation frequency, then expression, then methylation

message("Converging genes (3-omics): ", nrow(converging_genes))

# Save results
write_tsv(converging_genes, snakemake@output$converging_genes)
message("Converging genes saved: ", snakemake@output$converging_genes)

# --- Generate Summary Report ---
summary_df <- data.frame(
    metric = c("Driver genes (MAF >= threshold)", 
               "Over-expressed genes (LFC >= threshold)", 
               "Hypomethylated genes (logFC <= threshold)",
               "Genes with mutation + over-expression",
               "Genes with mutation + over-expression + hypomethylation"),
    count = c(
        nrow(driver_df %>% filter(maf >= maf_threshold)),
        nrow(deg_df %>% filter(log2FoldChange >= deg_lfc_threshold, padj <= deg_padj_threshold)),
        nrow(dmp_df %>% filter(logFC <= dmp_logfc_threshold, adj.P.Val <= dmp_padj_threshold)),
        nrow(mut_expr_genes),
        nrow(converging_genes)
    ),
    stringsAsFactors = FALSE
)

write_tsv(summary_df, snakemake@output$integration_summary)
message("Integration summary saved: ", snakemake@output$integration_summary)

# --- Generate Venn Diagram Data (for visualization) ---
venn_data <- list(
    driver_genes = driver_df %>% filter(maf >= maf_threshold) %>% pull(gene),
    over_expressed = deg_df %>% filter(log2FoldChange >= deg_lfc_threshold, padj <= deg_padj_threshold) %>% pull(gene),
    hypomethylated = dmp_df %>% filter(logFC <= dmp_logfc_threshold, adj.P.Val <= dmp_padj_threshold) %>% pull(gene)
)

venn_df <- data.frame(
    gene = unique(c(venn_data$driver_genes, venn_data$over_expressed, venn_data$hypomethylated)),
    in_driver = 0,
    in_deg = 0,
    in_dmp = 0,
    stringsAsFactors = FALSE
)

venn_df$in_driver <- as.integer(venn_df$gene %in% venn_data$driver_genes)
venn_df$in_deg <- as.integer(venn_df$gene %in% venn_data$over_expressed)
venn_df$in_dmp <- as.integer(venn_df$gene %in% venn_data$hypomethylated)

write_tsv(venn_df, snakemake@output$venn_data)
message("Venn diagram data saved: ", snakemake@output$venn_data)

message("=== MULTI-OMICS INTEGRATION MODULE COMPLETED: ", Sys.time(), " ===")
message("Converging genes found: ", nrow(converging_genes))
# Close sink connections
sink(type = "message")
sink(type = "output")
close(log_con)