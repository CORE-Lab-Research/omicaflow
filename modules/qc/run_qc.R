# M01 QC Module: Quality Control Processing
# Filters raw omics data based on QC thresholds

library(dplyr)
library(readr)
library(tidyr)
library(tibble)

# Get parameters from Snakemake
min_callrate <- snakemake@params$min_callrate
min_map_rate <- snakemake@params$min_map_rate
max_na <- snakemake@params$max_na
log_file <- snakemake@log[[1]]

# Setup logging
log_dir <- dirname(log_file)
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
log_con <- file(log_file, open = "a")
sink(log_con, type = "output")
sink(log_con, type = "message")
message("=== START QC MODULE: ", Sys.time(), " ===")

# Input validation
input_files <- c(snakemake@input$snv, snakemake@input$cnv, snakemake@input$rna, snakemake@input$methylation)
for (f in input_files) {
    if (!file.exists(f)) {
        message("ERROR: Input file not found: ", f)
        stop(paste("Missing input file:", f))
    }
}
message("All input files validated")

# --- SNV QC ---
message("Processing SNV data...")
snv_data <- read_tsv(snakemake@input$snv, show_col_types = FALSE)
# Filter by callrate per sample
sample_callrates <- colMeans(!is.na(snv_data[,-1]))
valid_samples_snv <- names(sample_callrates)[sample_callrates >= min_callrate]
message("SNV QC: ", length(valid_samples_snv), "/", length(sample_callrates), " samples passed callrate threshold")

# --- CNV QC ---
message("Processing CNV data...")
cnv_data <- read_tsv(snakemake@input$cnv, show_col_types = FALSE)
# Filter samples with at least 5 segments (standard filter for TCGA-like data)
sample_segment_counts <- cnv_data %>% count(Sample)
valid_samples_cnv <- sample_segment_counts %>% filter(n >= 5) %>% pull(Sample)
message("CNV QC: ", length(valid_samples_cnv), "/", nrow(sample_segment_counts), " samples passed segment threshold")

# --- RNA QC ---
message("Processing RNA data...")
rna_data <- read_tsv(snakemake@input$rna, show_col_types = FALSE)
# Filter samples with high zero/NA counts
sample_na_counts <- colSums(is.na(rna_data[,-1]))
valid_samples_rna <- names(sample_na_counts)[sample_na_counts / nrow(rna_data) <= 0.2] # max 20% NA
message("RNA QC: ", length(valid_samples_rna), "/", length(sample_na_counts), " samples passed NA threshold")

# --- Methylation QC ---
message("Processing Methylation data...")
meth_data <- read_tsv(snakemake@input$methylation, show_col_types = FALSE)
# Filter probes with high NA rate
probe_na_rates <- rowMeans(is.na(meth_data %>% select(-c(probe_id, gene_symbol))))
valid_probes_idx <- which(probe_na_rates <= max_na)
meth_filtered <- meth_data[valid_probes_idx, ]
message("Methylation QC: ", nrow(meth_filtered), "/", nrow(meth_data), " probes passed NA threshold")

# --- Sample Matching ---
message("Matching samples across platforms...")
final_samples <- intersect(intersect(valid_samples_snv, valid_samples_cnv), valid_samples_rna)
# Subset methylation samples if they exist in others
meth_samples <- colnames(meth_filtered)[!colnames(meth_filtered) %in% c("probe_id", "gene_symbol")]
final_samples <- intersect(final_samples, meth_samples)

message("Total matched samples: ", length(final_samples))

if (length(final_samples) == 0) {
    message("ERROR: No samples match across all platforms!")
    stop("Sample matching failed: zero common samples")
}

# --- Save Filtered Data ---
message("Saving filtered data...")
# SNV
snv_filtered <- snv_data %>% select(Hugo_Symbol, all_of(final_samples))
write_tsv(snv_filtered, snakemake@output$filtered_snv)

# CNV
cnv_filtered <- cnv_data %>% filter(Sample %in% final_samples)
write_tsv(cnv_filtered, snakemake@output$filtered_cnv)

# RNA
rna_filtered <- rna_data %>% select(gene_id, all_of(final_samples))
write_tsv(rna_filtered, snakemake@output$filtered_rna)

# Methylation
meth_filtered_final <- meth_filtered %>% select(probe_id, gene_symbol, all_of(final_samples))
write_tsv(meth_filtered_final, snakemake@output$filtered_meth)

# Sample List
write_tsv(data.frame(sample = final_samples), snakemake@output$sample_list)

# --- QC Report (Simple HTML) ---
qc_html <- paste0(
    "<html><head><title>QC Report</title></head><body>",
    "<h1>OmicaFlow QC Report</h1>",
    "<h2>SNV QC</h2><p>Total mutations: ", nrow(snv_data), 
    " | Valid samples (callrate >= ", min_callrate, "): ", length(valid_samples_snv), "</p>",
    "<h2>CNV QC</h2><p>Valid samples (>=5 segments): ", length(valid_samples_cnv), "</p>",
    "<h2>RNA QC</h2><p>Valid samples (NA <=20%): ", length(valid_samples_rna), "</p>",
    "<h2>Methylation QC</h2><p>Valid probes (NA <=", max_na*100, "%): ", length(valid_probes_idx), "</p>",
    "<h2>Final Sample List</h2><p>Total matched samples: ", length(final_samples), "</p>",
    "</body></html>"
)
writeLines(qc_html, snakemake@output$qc_report)
message("QC report saved: ", snakemake@output$qc_report)

message("=== QC MODULE COMPLETED: ", Sys.time(), " ===")
# Close sink connections
sink(type = "message")
sink(type = "output")
close(log_con)