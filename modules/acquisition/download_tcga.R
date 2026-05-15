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

# Define data types to download
data_configs <- list(
    snv = list(
        category = "Simple Nucleotide Variation",
        type = "Masked Somatic Mutation",
        workflow = "MuTect2",
        output = snakemake@output$snv
    ),
    cnv = list(
        category = "Copy Number Variation",
        type = "Copy Number Segment",
        workflow = "GISTIC2",
        output = snakemake@output$cnv
    ),
    rna = list(
        category = "Transcriptome Profiling",
        type = "Gene Expression Quantification",
        workflow = "HTSeq - FPKM",
        output = snakemake@output$rna
    ),
    methylation = list(
        category = "DNA Methylation",
        type = "Methylation Beta Value",
        workflow = "Illumina Human Methylation 450",
        output = snakemake@output$methylation
    )
)

all_barcodes <- c()

# Loop through and process each data type
for (name in names(data_configs)) {
    conf <- data_configs[[name]]
    message("--- Processing ", toupper(name), " ---")
    
    tryCatch({
        # Query
        query <- GDCquery(
            project = paste0("TCGA-", cancer_type),
            data.category = conf$category,
            data.type = conf$type,
            workflow.type = conf$workflow
        )
        
        # Download
        GDCdownload(query, method = "api", directory = cache_dir, files.per.chunk = 10)
        
        # Prepare
        data <- GDCprepare(query, directory = cache_dir)
        
        # Save output
        if (name == "snv") {
            write.table(data, conf$output, sep = "\t", quote = FALSE, row.names = FALSE)
            all_barcodes <- c(all_barcodes, data$Tumor_Sample_Barcode)
        } else if (name == "rna") {
            # Handle SummarizedExperiment or data.frame
            if (is(data, "SummarizedExperiment")) {
                expr_data <- as.data.frame(assay(data))
                expr_data <- expr_data %>% rownames_to_column("gene_id")
                write.table(expr_data, conf$output, sep = "\t", quote = FALSE, row.names = FALSE)
                all_barcodes <- c(all_barcodes, colnames(expr_data)[-1])
            } else {
                write.table(data, conf$output, sep = "\t", quote = FALSE, row.names = FALSE)
                all_barcodes <- c(all_barcodes, data$barcode)
            }
        } else {
            write.table(data, conf$output, sep = "\t", quote = FALSE, row.names = FALSE)
            if ("barcode" %in% colnames(data)) all_barcodes <- c(all_barcodes, data$barcode)
            if ("Sample" %in% colnames(data)) all_barcodes <- c(all_barcodes, data$Sample)
        }
        
        message(toupper(name), " data saved to ", conf$output)
        
    }, error = function(e) {
        message("ERROR processing ", name, ": ", e$message)
    })
}

# Create sample mapping table
message("Creating sample mapping table...")
sample_map <- data.frame(
    barcode = unique(all_barcodes),
    patient_id = substr(unique(all_barcodes), 1, 12),
    stringsAsFactors = FALSE
)
write.table(sample_map, snakemake@output$sample_map, sep = "\t", quote = FALSE, row.names = FALSE)
message("Sample map saved: ", snakemake@output$sample_map, " (", nrow(sample_map), " samples)")

message("=== ACQUISITION MODULE COMPLETED: ", Sys.time(), " ===")
# Close sink connections
sink(type = "output")
sink(type = "message")