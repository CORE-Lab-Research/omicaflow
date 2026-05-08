# M04 Methylation Analysis Module: DMP Detection and Annotation
# Uses minfi for methylation processing and limma for DMP calling

library(minfi)
library(limma)
library(dplyr)
library(readr)
library(ggplot2)

# Get parameters from Snakemake
filtered_meth <- snakemake@input$filtered_meth
sample_list <- snakemake@input$sample_list
padj_thresh <- snakemake@params$padj_threshold

# Create output directory
output_dir <- dirname(snakemake@output$dmp_results)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Load and prepare methylation data ---
meth_data <- read_tsv(filtered_meth, show_col_types = FALSE)

# For MVP, assume half are tumor, half are normal (simplified)
sample_df <- read_tsv(sample_list, show_col_types = FALSE)
n_samples <- ncol(meth_data) - 2  # minus probe_id and gene_symbol
tumor_samples <- colnames(meth_data)[3:(n_samples/2 + 2)]
normal_samples <- colnames(meth_data)[(n_samples/2 + 3):ncol(meth_data)]

# Create design matrix for limma
pheno <- data.frame(
    sample = colnames(meth_data)[!colnames(meth_data) %in% c("probe_id", "gene_symbol")],
    condition = factor(ifelse(colnames(meth_data)[!colnames(meth_data) %in% c("probe_id", "gene_symbol")] %in% tumor_samples, "Tumor", "Normal"))
)

meth_matrix <- meth_data %>%
    select(-c(probe_id, gene_symbol)) %>%
    as.matrix()

# --- DMP detection with limma ---
design <- model.matrix(~ condition, data = pheno)
fit <- lmFit(meth_matrix, design)
fit <- eBayes(fit)

dmp_results <- topTable(fit, coef = "conditionTumor", number = Inf, adjust.method = "BH") %>%
    rownames_to_column("probe_id") %>%
    filter(adj.P.Val <= padj_thresh) %>%
    left_join(meth_data %>% select(probe_id, gene_symbol), by = "probe_id") %>%
    arrange(adj.P.Val)

write_tsv(dmp_results, snakemake@output$dmp_results)

# --- Normalized methylation matrix (simplified: use beta values as-is) ---
norm_meth <- meth_data
write_tsv(norm_meth, snakemake@output$norm_meth)

# --- Functional annotation (simplified) ---
func_annotation <- dmp_results %>%
    mutate(
        location = case_when(
            grepl("TSS", gene_symbol) ~ "TSS",
            grepl("Exon", gene_symbol) ~ "Exon",
            TRUE ~ "Intergenic"
        ),
        cpg_island = grepl("CpG", probe_id)
    ) %>%
    select(probe_id, gene_symbol, location, cpg_island, logFC, adj.P.Val)

write_tsv(func_annotation, snakemake@output$func_annotation)

# --- Generate QC Report ---
# Simplified: create histogram of beta values
beta_long <- meth_matrix %>%
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "sample", values_to = "beta")

p <- ggplot(beta_long, aes(x = beta)) +
    geom_histogram(bins = 30, fill = "steelblue") +
    theme_minimal() +
    ggtitle("Distribution of Methylation Beta Values") +
    xlab("Beta Value (0=Unmethylated, 1=Methylated)")

ggsave(file.path(output_dir, "beta_hist.png"), p, width = 8, height = 6)

qc_html <- paste0(
    "<html><head><title>Methylation QC Report</title></head><body>",
    "<h1>Methylation Analysis QC Report</h1>",
    "<h2>DMP Summary</h2>",
    "<p>Total DMPs (padj <= ", padj_thresh, "): ", nrow(dmp_results), "</p>",
    "<h2>Beta Value Distribution</h2><img src='beta_hist.png' width='600px'>",
    "</body></html>"
)

writeLines(qc_html, snakemake@output$meth_qc)

message("Methylation analysis completed. DMPs found: ", nrow(dmp_results))