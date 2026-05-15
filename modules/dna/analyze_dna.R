# M02 DNA Analysis Module: SNV/CNV Analysis
# Uses maftools for SNV annotation and driver prediction

library(maftools)
library(readr)
library(dplyr)

# Get parameters from Snakemake
filtered_snv <- snakemake@input$filtered_snv
filtered_cnv <- snakemake@input$filtered_cnv
sample_list <- snakemake@input$sample_list
maf_threshold <- snakemake@params$maf_threshold
log_file <- snakemake@log[[1]]

# Setup logging
log_dir <- dirname(log_file)
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
sink(log_file, type = "output", append = TRUE)
sink(log_file, type = "message", append = TRUE)
message("=== START DNA ANALYSIS MODULE: ", Sys.time(), " ===")

# Input validation
input_files <- c(filtered_snv, filtered_cnv, sample_list)
for (f in input_files) {
    if (!file.exists(f)) {
        message("ERROR: Input file not found: ", f)
        stop(paste("Missing input file:", f))
    }
}
message("All input files validated")

# Create output directory
output_dir <- dirname(snakemake@output$annotated_snv)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
}

# --- SNV Analysis with maftools ---
message("Loading SNV data...")
maf_data <- read_tsv(filtered_snv, show_col_types = FALSE)
message("SNV data loaded: ", nrow(maf_data), " rows")

message("Loading sample metadata for maftools...")
sample_meta <- read_tsv(sample_list, show_col_types = FALSE) %>%
    rename(Tumor_Sample_Barcode = sample)

message("Creating MAF object with maftools...")
maf_obj <- read.maf(maf = filtered_snv, clinicalData = sample_meta)

# Get annotated SNV data
annotated_snv <- maf_obj@data
write_tsv(annotated_snv, snakemake@output$annotated_snv)
message("Annotated SNV saved: ", snakemake@output$annotated_snv, " (", nrow(annotated_snv), " rows)")

# Driver gene prediction (top mutated genes by MAF)
message("Identifying driver genes (MAF >= ", maf_threshold, ")...")
gene_summary <- annotated_snv %>%
    group_by(Hugo_Symbol) %>%
    summarise(
        n_samples = n_distinct(Tumor_Sample_Barcode),
        total_mutations = n(),
        .groups = "drop"
    ) %>%
    mutate(maf = n_samples / length(unique(sample_meta$Tumor_Sample_Barcode))) %>%
    filter(maf >= maf_threshold) %>%
    arrange(desc(maf))

write_tsv(gene_summary, snakemake@output$driver_genes)
message("Driver genes saved: ", snakemake@output$driver_genes, " (", nrow(gene_summary), " genes)")

# Mutational burden calculation
message("Calculating mutational burden...")
mut_burden <- annotated_snv %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(
        total_mutations = n(),
        missense = sum(Variant_Classification == "Missense_Mutation"),
        nonsense = sum(Variant_Classification == "Nonsense_Mutation"),
        frameshift = sum(grepl("frameshift", Variant_Classification, ignore.case = TRUE)),
        .groups = "drop"
    )

write_tsv(mut_burden, snakemake@output$mutational_burden)
message("Mutational burden saved: ", snakemake@output$mutational_burden, " (", nrow(mut_burden), " samples)")

# --- CNV Analysis (simplified for MVP) ---
message("Loading CNV data...")
cnv_data <- read_tsv(filtered_cnv, show_col_types = FALSE)
message("CNV data loaded: ", nrow(cnv_data), " rows, ", length(unique(cnv_data$Sample)), " samples")

# Thresholds for Amp/Del (could be moved to config)
amp_thresh <- 0.5
del_thresh <- -0.5

message("Summarizing CNV events (Amp > ", amp_thresh, ", Del < ", del_thresh, ")...")
cnv_summary <- cnv_data %>%
    filter(segment_mean > amp_thresh | segment_mean < del_thresh) %>%
    group_by(Sample, gene_symbol) %>%
    summarise(
        avg_segment_mean = mean(segment_mean),
        event_type = ifelse(avg_segment_mean > amp_thresh, "amplification", "deletion"),
        .groups = "drop"
    ) %>%
    group_by(gene_symbol, event_type) %>%
    summarise(
        n_samples = n_distinct(Sample),
        .groups = "drop"
    ) %>%
    arrange(desc(n_samples))

write_tsv(cnv_summary, snakemake@output$cnv_results)
message("CNV summary saved: ", snakemake@output$cnv_results, " (", nrow(cnv_summary), " entries)")

# --- DNA Integration Summary ---
message("Generating DNA integration summary...")
dna_summary <- data.frame(
    metric = c("Total SNV mutations", "Driver genes (MAF >= threshold)", "Samples with CNV events"),
    value = c(
        nrow(annotated_snv),
        nrow(gene_summary),
        length(unique(cnv_data$Sample))
    ),
    stringsAsFactors = FALSE
)

write_tsv(dna_summary, snakemake@output$dna_summary)
message("DNA summary saved: ", snakemake@output$dna_summary)

message("=== DNA ANALYSIS MODULE COMPLETED: ", Sys.time(), " ===")
message("Driver genes found: ", nrow(gene_summary))
# Close sink connections
sink(type = "output")
sink(type = "message")