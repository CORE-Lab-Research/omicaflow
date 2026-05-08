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

# Create output directory
output_dir <- dirname(snakemake@output$annotated_snv)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- SNV Analysis with maftools ---
maf_data <- read_tsv(filtered_snv, show_col_types = FALSE)

# Create MAF object for maftools
maf_obj <- read.maf(maf = filtered_snv)

# Get annotated SNV data
annotated_snv <- maf_obj@data
write_tsv(annotated_snv, snakemake@output$annotated_snv)

# Driver gene prediction (top mutated genes by MAF)
gene_summary <- annotated_snv %>%
    group_by(Hugo_Symbol) %>%
    summarise(
        n_samples = n_distinct(Tumor_Sample_Barcode),
        total_mutations = n(),
        .groups = "drop"
    ) %>%
    mutate(maf = n_samples / length(unique(annotated_snv$Tumor_Sample_Barcode))) %>%
    filter(maf >= maf_threshold) %>%
    arrange(desc(maf))

write_tsv(gene_summary, snakemake@output$driver_genes)

# Mutational burden calculation
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

# --- CNV Analysis (simplified for MVP) ---
cnv_data <- read_tsv(filtered_cnv, show_col_types = FALSE)

# Simple CNV summary: count amplifications/deletions per gene
cnv_summary <- cnv_data %>%
    filter(segment_mean > 0.5 | segment_mean < -0.5) %>%  # Threshold for amp/del
    group_by(Sample, gene_symbol) %>%
    summarise(
        avg_segment_mean = mean(segment_mean),
        event_type = ifelse(avg_segment_mean > 0.5, "amplification", "deletion"),
        .groups = "drop"
    ) %>%
    group_by(gene_symbol, event_type) %>%
    summarise(
        n_samples = n_distinct(Sample),
        .groups = "drop"
    ) %>%
    arrange(desc(n_samples))

write_tsv(cnv_summary, snakemake@output$cnv_results)

# --- DNA Integration Summary ---
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

message("DNA analysis completed. Driver genes found: ", nrow(gene_summary))