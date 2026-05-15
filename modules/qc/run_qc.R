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
sink(log_file, type = "output", append = TRUE)
sink(log_file, type = "message", append = TRUE)
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

# Create output directory
output_dir <- dirname(snakemake@output$qc_report)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
}

# --- DNA SNV QC ---
message("Processing SNV QC...")
snv_data <- read_tsv(snakemake@input$snv, show_col_types = FALSE)
message("SNV data loaded: ", nrow(snv_data), " rows, ", length(unique(snv_data$Tumor_Sample_Barcode)), " samples")

# Filter SNV samples by mutation count (proxy for call rate in MAF)
sample_counts <- snv_data %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(n_mut = n(), .groups = "drop")
avg_mut <- mean(sample_counts$n_mut)
message("Average mutations per sample: ", round(avg_mut, 2))

valid_samples_snv <- sample_counts %>%
    filter(n_mut >= avg_mut * 0.1) %>% # Filter extreme outliers with very few mutations
    pull(Tumor_Sample_Barcode)
message("Valid SNV samples (>= 10% of avg muts): ", length(valid_samples_snv))

filtered_snv <- snv_data %>%
    filter(Tumor_Sample_Barcode %in% valid_samples_snv)
write_tsv(filtered_snv, snakemake@output$filtered_snv)
message("Filtered SNV saved: ", snakemake@output$filtered_snv)

# --- DNA CNV QC ---
message("Processing CNV QC...")
cnv_data <- read_tsv(snakemake@input$cnv, show_col_types = FALSE)
message("CNV data loaded: ", nrow(cnv_data), " rows, ", length(unique(cnv_data$Sample)), " samples")

valid_samples_cnv <- cnv_data %>%
    group_by(Sample) %>%
    summarise(n_segments = n(), .groups = "drop") %>%
    filter(n_segments >= 5) %>%  # Minimum 5 CNV segments
    pull(Sample)
message("Valid CNV samples (>=5 segments): ", length(valid_samples_cnv))

filtered_cnv <- cnv_data %>%
    filter(Sample %in% valid_samples_cnv)
write_tsv(filtered_cnv, snakemake@output$filtered_cnv)
message("Filtered CNV saved: ", snakemake@output$filtered_cnv)

# --- RNA QC ---
message("Processing RNA QC...")
rna_data <- read_tsv(snakemake@input$rna, show_col_types = FALSE)
message("RNA data loaded: ", nrow(rna_data), " genes, ", ncol(rna_data)-1, " samples")

# Filter RNA by missingness (samples with >20% missing values)
sample_na_rates <- rna_data %>%
    select(-gene_id) %>%
    summarise(across(everything(), ~ mean(is.na(.))) * 100) %>%
    pivot_longer(everything(), names_to = "sample", values_to = "na_pct")
message("RNA NA rates calculated for ", nrow(sample_na_rates), " samples")

valid_samples_rna <- sample_na_rates %>%
    filter(na_pct <= 20) %>%
    pull(sample)
message("Valid RNA samples (NA <=20%): ", length(valid_samples_rna))

filtered_rna <- rna_data %>%
    select(gene_id, all_of(valid_samples_rna))
write_tsv(filtered_rna, snakemake@output$filtered_rna)
message("Filtered RNA saved: ", snakemake@output$filtered_rna)

# --- Methylation QC ---
message("Processing Methylation QC...")
meth_data <- read_tsv(snakemake@input$methylation, show_col_types = FALSE)
message("Methylation data loaded: ", nrow(meth_data), " probes, ", ncol(meth_data)-2, " samples")

# Filter probes with > threshold missing values (per row)
# Note: probes are rows, samples are columns
probe_na_rates <- rowMeans(is.na(meth_data %>% select(-c(probe_id, gene_symbol))))
message("Methylation probe NA rates calculated")

valid_probes_idx <- which(probe_na_rates <= max_na)
message("Valid methylation probes (NA <=", max_na*100, "%): ", length(valid_probes_idx))

filtered_meth <- meth_data[valid_probes_idx, ]
write_tsv(filtered_meth, snakemake@output$filtered_meth)
message("Filtered Methylation saved: ", snakemake@output$filtered_meth)

# --- Generate final sample list ---
message("Generating final sample list...")
final_samples <- intersect(
    intersect(valid_samples_snv, valid_samples_cnv),
    intersect(valid_samples_rna, colnames(filtered_meth)[!colnames(filtered_meth) %in% c("probe_id", "gene_symbol")])
)
message("Total matched samples: ", length(final_samples))

final_sample_df <- data.frame(
    sample = final_samples,
    patient_id = substr(final_samples, 1, 12),
    stringsAsFactors = FALSE
)
write_tsv(final_sample_df, snakemake@output$sample_list)
message("Sample list saved: ", snakemake@output$sample_list)

# --- Generate QC Report (simple HTML) ---
message("Generating QC report...")
qc_html <- paste0(
    "<html><head><title>QC Report - ", snakemake@params$cancer_type, "</title></head><body>",
    "<h1>OmicaFlow QC Report</h1>",
    "<h2>SNV QC</h2><p>Initial samples: ", length(unique(snv_data$Tumor_Sample_Barcode)),
    " | Valid samples (callrate >= ", min_callrate, "): ", length(valid_samples_snv), "</p>",
    "<h2>CNV QC</h2><p>Valid samples (>=5 segments): ", length(valid_samples_cnv), "</p>",
    "<h2>RNA QC</h2><p>Valid samples (NA <=20%): ", length(valid_samples_rna), "</p>",
    "<h2>Methylation QC</h2><p>Valid probes (NA <=", max_na*100, "%): ", length(valid_probes), "</p>",
    "<h2>Final Sample List</h2><p>Total matched samples: ", length(final_samples), "</p>",
    "</body></html>"
)
writeLines(qc_html, snakemake@output$qc_report)
message("QC report saved: ", snakemake@output$qc_report)

message("=== QC MODULE COMPLETED: ", Sys.time(), " ===")
# Close sink connections
sink(type = "output")
sink(type = "message")