# M01 QC Module: Quality Control Processing
# Filters raw omics data based on QC thresholds

library(dplyr)
library(readr)

# Get parameters from Snakemake
min_callrate <- snakemake@params$min_callrate
min_map_rate <- snakemake@params$min_map_rate
max_na <- snakemake@params$max_na

# Create output directory
output_dir <- dirname(snakemake@output$qc_report)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- DNA SNV QC ---
snv_data <- read_tsv(snakemake@input$snv, show_col_types = FALSE)

# Filter SNV by call rate (per sample)
sample_callrates <- snv_data %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(call_rate = n() / nrow(snv_data), .groups = "drop")

valid_samples_snv <- sample_callrates %>%
    filter(call_rate >= min_callrate) %>%
    pull(Tumor_Sample_Barcode)

filtered_snv <- snv_data %>%
    filter(Tumor_Sample_Barcode %in% valid_samples_snv)

write_tsv(filtered_snv, snakemake@output$filtered_snv)

# --- DNA CNV QC ---
cnv_data <- read_tsv(snakemake@input$cnv, show_col_types = FALSE)
valid_samples_cnv <- cnv_data %>%
    group_by(Sample) %>%
    summarise(n_segments = n(), .groups = "drop") %>%
    filter(n_segments >= 5) %>%  # Minimum 5 CNV segments
    pull(Sample)

filtered_cnv <- cnv_data %>%
    filter(Sample %in% valid_samples_cnv)

write_tsv(filtered_cnv, snakemake@output$filtered_cnv)

# --- RNA QC ---
rna_data <- read_tsv(snakemake@input$rna, show_col_types = FALSE)

# Filter RNA by missingness (samples with >20% missing values)
sample_na_rates <- rna_data %>%
    select(-gene_id) %>%
    summarise(across(everything(), ~ mean(is.na(.))) * 100) %>%
    pivot_longer(everything(), names_to = "sample", values_to = "na_pct")

valid_samples_rna <- sample_na_rates %>%
    filter(na_pct <= 20) %>%
    pull(sample)

filtered_rna <- rna_data %>%
    select(gene_id, all_of(valid_samples_rna))

write_tsv(filtered_rna, snakemake@output$filtered_rna)

# --- Methylation QC ---
meth_data <- read_tsv(snakemake@input$methylation, show_col_types = FALSE)

# Filter probes with >10% missing values
probe_na_rates <- meth_data %>%
    select(-c(probe_id, gene_symbol)) %>%
    summarise(across(everything(), ~ mean(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "probe", values_to = "na_rate")

valid_probes <- probe_na_rates %>%
    filter(na_rate <= max_na) %>%
    pull(probe)

filtered_meth <- meth_data %>%
    filter(probe_id %in% valid_probes)

write_tsv(filtered_meth, snakemake@output$filtered_meth)

# --- Generate final sample list ---
final_samples <- intersect(
    intersect(valid_samples_snv, valid_samples_cnv),
    intersect(valid_samples_rna, colnames(filtered_meth)[!colnames(filtered_meth) %in% c("probe_id", "gene_symbol")])
)

final_sample_df <- data.frame(
    sample = final_samples,
    patient_id = substr(final_samples, 1, 12),
    stringsAsFactors = FALSE
)

write_tsv(final_sample_df, snakemake@output$sample_list)

# --- Generate QC Report (simple HTML) ---
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

message("QC processing completed. Final matched samples: ", length(final_samples))