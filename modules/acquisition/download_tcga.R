# M00 Acquisition Module: TCGA Data Download
# Uses TCGAbiolinks to download multi-omics data for a given cancer type

library(TCGAbiolinks)
library(dplyr)

# Get parameters from Snakemake
cancer_type <- snakemake@params$cancer_type
cache_dir <- snakemake@params$cache_dir
tcga_url <- snakemake@params$tcga_url
log_file <- snakemake@log[[1]]

# Setup logging: create log directory, sink output to file + terminal
log_dir <- dirname(log_file)
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
sink(log_file, type = "output", append = TRUE)
sink(log_file, type = "message", append = TRUE)
message("=== START ACQUISITION MODULE: ", Sys.time(), " ===")
message("Cancer type: ", cancer_type)
message("Cache directory: ", cache_dir)

# Create cache directory if not exists
if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
    message("Created cache directory: ", cache_dir)
}

# Query TCGA data
message("Querying TCGA GDC API...")
query <- GDCquery(
    project = paste0("TCGA-", cancer_type),
    data.category = c("Simple Nucleotide Variation", "Copy Number Variation", "Transcriptome Profiling", "DNA Methylation"),
    data.type = c("Masked Somatic Mutation", "Copy Number Segment", "Gene Expression Quantification", "Methylation Beta Value"),
    workflow.type = c("MuTect2", "GISTIC2", "HTSeq - FPKM", "Illumina Human Methylation 450")
)
message("TCGA query completed")

# Download data
message("Downloading TCGA data...")
GDCdownload(query, method = "api", directory = cache_dir, files.per.chunk = 10)
message("Download completed")

# Prepare output directories
output_dir <- dirname(snakemake@output$snv)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
}

# Process and save SNV data
message("Processing SNV data...")
snv_data <- GDCprepare(query, directory = cache_dir) %>%
    filter(data.type == "Masked Somatic Mutation") %>%
    select(-starts_with("unused"))
write.table(snv_data, snakemake@output$snv, sep = "\t", quote = FALSE, row.names = FALSE)
message("SNV data saved: ", snakemake@output$snv, " (", nrow(snv_data), " rows)")

# Process and save CNV data
message("Processing CNV data...")
cnv_data <- GDCprepare(query, directory = cache_dir) %>%
    filter(data.type == "Copy Number Segment") %>%
    select(-starts_with("unused"))
write.table(cnv_data, snakemake@output$cnv, sep = "\t", quote = FALSE, row.names = FALSE)
message("CNV data saved: ", snakemake@output$cnv, " (", nrow(cnv_data), " rows)")

# Process and save RNA-seq data
message("Processing RNA-seq data...")
rna_data <- GDCprepare(query, directory = cache_dir) %>%
    filter(data.type == "Gene Expression Quantification") %>%
    select(-starts_with("unused"))
write.table(rna_data, snakemake@output$rna, sep = "\t", quote = FALSE, row.names = FALSE)
message("RNA-seq data saved: ", snakemake@output$rna, " (", nrow(rna_data), " rows)")

# Process and save Methylation data
message("Processing Methylation data...")
meth_data <- GDCprepare(query, directory = cache_dir) %>%
    filter(data.type == "DNA Methylation") %>%
    select(-starts_with("unused"))
write.table(meth_data, snakemake@output$methylation, sep = "\t", quote = FALSE, row.names = FALSE)
message("Methylation data saved: ", snakemake@output$methylation, " (", nrow(meth_data), " rows)")

# Create sample mapping table
message("Creating sample mapping table...")
sample_map <- data.frame(
    barcode = unique(c(snv_data$Tumor_Sample_Barcode, rna_data$barcode, meth_data$barcode)),
    patient_id = substr(barcode, 1, 12),
    stringsAsFactors = FALSE
)
write.table(sample_map, snakemake@output$sample_map, sep = "\t", quote = FALSE, row.names = FALSE)
message("Sample map saved: ", snakemake@output$sample_map, " (", nrow(sample_map), " samples)")

message("=== ACQUISITION MODULE COMPLETED: ", Sys.time(), " ===")
# Close sink connections
sink(type = "output")
sink(type = "message")