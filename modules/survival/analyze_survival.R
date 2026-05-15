# M07 Survival Analysis Module: Test Prognostic Value of Genes
# Tests if gene expression stratifies patients into significant survival groups

library(survival)
library(survminer)
library(dplyr)
library(readr)

# Get parameters from Snakemake
converging_genes <- snakemake@input$converging_genes
norm_expr <- snakemake@input$norm_expr
clinical_data <- snakemake@input$clinical_data  # Expected: barcode, OS.time, OS.status, PFS.time, PFS.status
log_file <- snakemake@log[[1]]

# Setup logging
log_dir <- dirname(log_file)
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
log_con <- file(log_file, open = "a")
sink(log_con, type = "output")
sink(log_con, type = "message")
message("=== START SURVIVAL ANALYSIS MODULE: ", Sys.time(), " ===")

# Input validation
input_files <- c(converging_genes, norm_expr, clinical_data)
for (f in input_files) {
    if (!file.exists(f)) {
        message("ERROR: Input file not found: ", f)
        stop(paste("Missing input file:", f))
    }
}
message("All input files validated")

# Create output directory
output_dir <- dirname(snakemake@output$survival_results)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
}

# --- Load Data ---
message("Loading converging genes...")
conv_genes_df <- read_tsv(converging_genes, show_col_types = FALSE)
message("Converging genes loaded: ", nrow(conv_genes_df), " genes")

message("Loading normalized expression matrix...")
expr_data <- read_tsv(norm_expr, show_col_types = FALSE)
message("Expression data loaded: ", nrow(expr_data), " genes, ", ncol(expr_data)-1, " samples")

message("Loading clinical data...")
clin_data <- read_tsv(clinical_data, show_col_types = FALSE)
message("Clinical data loaded: ", nrow(clin_data), " patients")
message("Available columns: ", paste(colnames(clin_data), collapse = ", "))

# --- Prepare Data for Survival Analysis ---
# For simplicity, use Overall Survival (OS) if available, else PFS
surv_time_col <- if ("OS.time" %in% colnames(clin_data)) "OS.time" else if ("PFS.time" %in% colnames(clin_data)) "PFS.time" else NULL
surv_status_col <- if ("OS.status" %in% colnames(clin_data)) "OS.status" else if ("PFS.status" %in% colnames(clin_data)) "PFS.status" else NULL

if (is.null(surv_time_col) || is.null(surv_status_col)) {
    message("WARNING: Standard survival columns not found. Using first two numeric columns as time/status.")
    # Fallback: use first two columns after barcode
    surv_cols <- clin_data %>%
        select(-barcode) %>%
        select(where(is.numeric)) %>%
        colnames()
    if (length(surv_cols) >= 2) {
        surv_time_col <- surv_cols[1]
        surv_status_col <- surv_cols[2]
        message("Using columns: time = ", surv_time_col, ", status = ", surv_status_col)
    } else {
        stop("Cannot find suitable survival time/status columns in clinical data")
    }
}

message("Using survival endpoint: time = ", surv_time_col, ", status = ", surv_status_col)

# Extract survival data
surv_data <- clin_data %>%
    select(barcode, all_of(c(surv_time_col, surv_status_col))) %>%
    rename(time = !!sym(surv_time_col), status = !!sym(surv_status_col)) %>%
    mutate(
        time = as.numeric(time),
        status = as.numeric(status)
    ) %>%
    filter(!is.na(time) & !is.na(status))

message("Valid survival samples: ", nrow(surv_data))

# --- Expression Data Preparation ---
# Ensure expression data has gene_id as first column
expr_matrix <- expr_data %>%
    column_to_rownames(var = "gene_id") %>%
    as.matrix()

# Get common samples between expression and clinical data
common_samples <- intersect(colnames(expr_matrix), surv_data$barcode)
message("Common samples (expression + clinical): ", length(common_samples))

if (length(common_samples) < 10) {
    stop("Too few common samples between expression and clinical data")
}

# Subset data to common samples
expr_common <- expr_matrix[, common_samples]
surv_common <- surv_data %>%
    filter(barcode %in% common_samples) %>%
    arrange(barcode)

# Ensure sample order matches
expr_common <- expr_common[, surv_common$barcode]

# --- Survival Analysis Per Gene ---
message("Starting survival analysis for ", nrow(conv_genes_df), " converging genes...")

survival_results <- data.frame(
    gene = character(),
    n_samples = integer(),
    n_events = integer(),
    hr = numeric(),
    hr_lower = numeric(),
    hr_upper = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
)

# For each converging gene, test if high vs low expression predicts survival
for (i in 1:nrow(conv_genes_df)) {
    gene_name <- conv_genes_df$gene[i]
    
    # Get expression values for this gene
    expr_values <- expr_common[gene_name, ]
    
    # Skip if all values are NA or constant
    if (all(is.na(expr_values)) || sd(expr_values, na.rm = TRUE) == 0) {
        next
    }
    
    # Stratify by median expression (high vs low)
    median_expr <- median(expr_values, na.rm = TRUE)
    risk_group <- ifelse(expr_values > median_expr, "High", "Low")
    
    # Create survival data frame for this gene
    gene_surv_data <- data.frame(
        time = surv_common$time,
        status = surv_common$status,
        risk = risk_group,
        stringsAsFactors = FALSE
    ) %>%
        filter(!is.na(time) & !is.na(status))
    
    # Skip if too few events
    if (sum(gene_surv_data$status) < 5) {
        next
    }
    
    # Fit Kaplan-Meier and log-rank test
    km_fit <- survfit(Surv(time, status) ~ risk, data = gene_surv_data)
    logrank_test <- survdiff(Surv(time, status) ~ risk, data = gene_surv_data)
    
    # Calculate hazard ratio from Cox PH model
    cox_model <- coxph(Surv(time, status) ~ risk, data = gene_surv_data)
    hr_summary <- summary(cox_model)
    
    # Extract results
    hr <- hr_summary$conf.int[1, "exp(coef)"]
    hr_lower <- hr_summary$conf.int[1, "lower .95"]
    hr_upper = hr_summary$conf.int[1, "upper .95"]
    p_value <- logrank_test$pvalue
    n_samples <- nrow(gene_surv_data)
    n_events <- sum(gene_surv_data$status)
    
    # Store results
    survival_results <- rbind(survival_results, data.frame(
        gene = gene_name,
        n_samples = n_samples,
        n_events = n_events,
        hr = hr,
        hr_lower = hr_lower,
        hr_upper = hr_upper,
        p_value = p_value,
        stringsAsFactors = FALSE
    ))
    
    # Progress indicator every 10 genes
    if (i %% 10 == 0) {
        message("Processed ", i, "/", nrow(conv_genes_df), " genes...")
    }
}

message("Survival analysis completed for ", nrow(survival_results), " genes")

# --- Adjust for Multiple Testing (FDR) ---
if (nrow(survival_results) > 0) {
    survival_results <- survival_results %>%
        mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
        arrange(p_adj)
    
    message("Significant genes (FDR < 0.05): ", sum(survival_results$p_adj < 0.05))
} else {
    survival_results$p_adj <- numeric()
}

# Save results
write_tsv(survival_results, snakemake@output$survival_results)
message("Survival results saved: ", snakemake@output$survival_results)

# --- Generate Kaplan-Meier Plots for Top Genes ---
if (nrow(survival_results) > 0) {
    # Get top 3 significant genes (or all if less than 3)
    top_n <- min(3, nrow(survival_results))
    top_genes <- survival_results[1:top_n, ]$gene
    
    for (gene_name in top_genes) {
        # Get expression values
        expr_values <- expr_common[gene_name, ]
        median_expr <- median(expr_values, na.rm = TRUE)
        risk_group <- ifelse(expr_values > median_expr, "High", "Low")
        
        # Create survival data
        gene_surv_data <- data.frame(
            time = surv_common$time,
            status = surv_common$status,
            risk = risk_group,
            stringsAsFactors = FALSE
        ) %>%
            filter(!is.na(time) & !is.na(status))
        
        # Skip if insufficient data
        if (nrow(gene_surv_data) < 10 || sum(gene_surv_data$status) < 3) {
            next
        }
        
        # Create Kaplan-Meier plot
        km_fit <- survfit(Surv(time, status) ~ risk, data = gene_surv_data)
        
        # Survival plot
        surv_plot <- ggsurvplot(
            km_fit,
            data = gene_surv_data,
            pval = TRUE,
            conf.int = TRUE,
            risk.table = TRUE,
            risk.table.col = "strata",
            ggtheme = theme_minimal(),
            palette = c("#E7B800", "#2E9FDF"),
            title = paste("Survival by", gene_name, "Expression"),
            subtitle = paste("HR =", round(survival_results[survival_results$gene == gene_name, "hr"], 2),
                           ", p =", signif(survival_results[survival_results$gene == gene_name, "p_value"], 3))
        )
        
        # Save plot
        plot_file <- file.path(output_dir, paste0("km_", gene_name, ".pdf"))
        ggsave(plot_file, surv_plot$plot, width = 8, height = 6)
        message("Kaplan-Meier plot saved: ", plot_file)
    }
}

message("=== SURVIVAL ANALYSIS MODULE COMPLETED: ", Sys.time(), " ===")
message("Significant prognostic genes (FDR < 0.05): ", 
        if (nrow(survival_results) > 0) sum(survival_results$p_adj < 0.05) else 0)
# Close sink connections
sink(type = "message")
sink(type = "output")
close(log_con)