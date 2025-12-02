# ==============================================================================
# R SCRIPT FOR COMPREHENSIVE RNA-SEQ ANALYSIS WITH DESEQ2
# Includes: Setup, Data Loading, Filtering, DESeq2 Pipeline, 
#           TPM Calculation (Effective Length), and ALL Plotting Sections (5-11).
# ------------------------------------------------------------------------------
# CORRECTED ERROR: 'type='apeglm' shrinkage only for use with 'coef'' in Section 8.
# Fix: The lfcShrink call now uses the 'coef' argument instead of 'contrast'
#      for apeglm shrinkage, based on the contrast names from resultsNames(dds).
# ==============================================================================

# Load packages
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(viridis)

# --- File Configurations (MUST BE SET) ---
COUNT_FILE <- "featurecounts/star_results/raw_counts_GENE_matrix.txt"
SAMPLE_TABLE_FILE <- "featurecounts/star_results/sample_table_deseq.csv"
ANNOTATION_FILE <- "kallisto/output/description.tsv" # Annotation file path
REPORT_FILE <- file.path(dirname(COUNT_FILE), "ANALYSIS_REPORT.txt") # Final report file

# Initialize report metrics object
report_metrics <- list()

# ==============================================================================
# 1. SETUP AND SAMPLE TABLE CREATION (COLDATA)
# ==============================================================================

cat("\n==============================================================================")
cat("\n1. SETUP AND SAMPLE TABLE CREATION (COLDATA)")
cat("\n==============================================================================\n")

# Mapping sample names to experimental conditions
sample_data <- data.frame(
    SampleName = c("C_1_SRR15312907", "C_2_SRR15312908", "C_3_SRR15312910", "C_4_SRR15312911",
                   "H_1_SRR15312890", "H_2_SRR15312891",
                   "P_1_SRR15312892", "P_2_SRR15312893", "P_3_SRR15312894", "P_4_SRR15312895",
                   "P_5_SRR15312896", "P_6_SRR15312897", "P_7_SRR15312898", "P_8_SRR15312899",
                   "P_9_SRR15312900", "P_10_SRR15312901", "P_11_SRR15312902", "P_12_SRR15312903",
                   "P_13_SRR15312904", "P_14_SRR15312905", "P_15_SRR15312906"),
    Condition_Code = c(rep("C", 4), rep("H", 2), rep("P", 15))
)

# Derive metadata columns
sample_data$SimpleID <- gsub("_SRR[0-9]+", "", sample_data$SampleName)
condition_mapping <- c(C = "Control", H = "High", P = "Low")
sample_data$Condition <- condition_mapping[sample_data$Condition_Code]

# Create 'coldata' object
coldata <- data.frame(
    row.names = sample_data$SimpleID,
    Condition = factor(sample_data$Condition, levels = c("Control", "High", "Low")),
    ShortID = sample_data$SampleName
)

write.csv(coldata, SAMPLE_TABLE_FILE, row.names = TRUE)


# --- USER CONFIGURATION: MEAN FRAGMENT LENGTHS (MFL) ---
# MFL is required for accurate TPM calculation using the Effective Length method.
# This vector must contain the Mean Fragment Length for EACH SAMPLE, 
# in the order of the columns in the 'coldata' object (which is the sample order).
# Placeholder set to 450bp for all 21 samples.
# !!! IMPORTANT: REPLACE THESE PLACEHOLDERS WITH YOUR ACTUAL MFL VALUES !!!
MEAN_FRAGMENT_LENGTHS <- rep(450, length(coldata$Condition))
cat(paste("Using Mean Fragment Length (MFL) placeholder for all samples:", MEAN_FRAGMENT_LENGTHS[1], "bp\n"))


# ==============================================================================
# 2. LOAD COUNTS, ROBUST MAPPING AND INITIAL FILTERING
# ==============================================================================

cat("\n==============================================================================")
cat("\n2. LOAD COUNTS, ROBUST MAPPING AND INITIAL FILTERING")
cat("\n==============================================================================\n")

# Load complete data (including Length for TPM calculation)
raw_data_with_meta <- read.table(COUNT_FILE, header = TRUE, row.names = 1, skip = 1, sep = "\t",
                                 stringsAsFactors = FALSE, check.names = FALSE)
gene_lengths <- raw_data_with_meta[, "Length", drop = FALSE] # Isolate gene lengths
countData_all <- raw_data_with_meta[, 6:ncol(raw_data_with_meta)] # Count data starting from column 6

# Store initial gene count for report
report_metrics$initial_genes <- nrow(raw_data_with_meta)

# Robust Column Mapping
original_count_names <- colnames(countData_all)
base_names <- basename(original_count_names)
count_column_ShortIDs <- gsub(".Aligned.sortedByCoord.out.bam", "", base_names)
count_column_ShortIDs <- gsub("^star\\.alignment\\.", "", count_column_ShortIDs)

match_index <- match(count_column_ShortIDs, coldata$ShortID)
if(any(is.na(match_index))) {stop("ERROR: Count file names do not match coldata ShortIDs.")}

new_colnames <- rownames(coldata)[match_index]
colnames(countData_all) <- new_colnames
countData <- countData_all[, rownames(coldata)] # Ensure order matches coldata

# Prepare integer matrix for DESeq2
countData_deseq <- round(as.matrix(countData))

# 2.1. Initial Filtering
MIN_RAW_COUNT <- 5
MIN_SAMPLES_EXPRESSED <- 2

# Keep genes where raw count >= MIN_RAW_COUNT in at least MIN_SAMPLES_EXPRESSED
keep <- rowSums(countData_deseq >= MIN_RAW_COUNT) >= MIN_SAMPLES_EXPRESSED
countData_deseq <- countData_deseq[keep, ]
gene_lengths_filtered <- gene_lengths[keep, , drop = FALSE] # Filter lengths to match counts

# Store filtered count for report
report_metrics$filtered_genes <- nrow(countData_deseq)
report_metrics$filtering_criteria <- paste0("Minimum raw count: ", MIN_RAW_COUNT,
                                            " in at least ", MIN_SAMPLES_EXPRESSED, " samples.")

cat(paste("Genes retained after filtering:", nrow(countData_deseq), "\n"))

# ==============================================================================
# 3. CREATE DESeq2 OBJECT AND RUN ANALYSIS (RAW COUNTS)
# ==============================================================================

cat("\n==============================================================================")
cat("\n3. CREATE DESeq2 OBJECT AND RUN ANALYSIS (RAW COUNTS)")
cat("\n==============================================================================\n")

dds <- DESeqDataSetFromMatrix(
    countData = countData_deseq,
    colData = coldata,
    design = ~ Condition
)

dds$Condition <- relevel(dds$Condition, ref = "Control")
dds <- DESeq(dds)

cat("\n--- DESeq2 Analysis Complete. ---\n")
saveRDS(dds, file.path(dirname(COUNT_FILE), "dds_object.rds"))


# ==============================================================================
# 4. TPM CALCULATION AND LOG TRANSFORMATION FOR EDA
# ==============================================================================

cat("\n==============================================================================")
cat("\n4. TPM CALCULATION AND LOG TRANSFORMATION FOR EDA (Effective Length Method)")
cat("\n==============================================================================\n")

# Function to calculate TPM using the Effective Length (Gene Length - MFL + 1)
calculate_tpm <- function(counts, featureLength, meanFragmentLength) {
    # 1. ARGUMENT VALIDATION
    # Ensures feature lengths and mean fragment lengths match matrix dimensions
    stopifnot(length(featureLength) == nrow(counts))
    stopifnot(length(meanFragmentLength) == ncol(counts))
    
    # 2. EFFECTIVE LENGTH CALCULATION
    # Effective Length (effLen) = Gene Length - Mean Fragment Length + 1
    # 'lapply' ensures the MFL (which is sample-specific) is applied column-wise
    effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
        featureLength - meanFragmentLength[i] + 1
    }))
    
    # 3. SHORT GENE FILTERING
    # Filter out genes with Effective Length <= 1 (i.e., genes shorter than the fragment length)
    idx <- apply(effLen, 1, function(x) min(x) > 1)
    
    # Apply filter to counts and effective lengths
    counts_filtered <- counts[idx,]
    effLen_filtered <- effLen[idx,]
    
    # 4. TPM CALCULATION IN LOG-SPACE (for numerical stability)
    tpm <- do.call(cbind, lapply(1:ncol(counts_filtered), function(i) {
        # Rate = log(Counts) - log(Effective Length) (Equivalent to RPK in log space)
        rate = log(counts_filtered[,i]) - log(effLen_filtered[,i])
        
        # Denominator = log(sum(exp(Rate))) -> Scaling factor in log space
        denom = log(sum(exp(rate)))
        
        # TPM = exp(Rate - Denominator + log(1e6))
        exp(rate - denom + log(1e6))
    }))
    
    # 5. ASSIGN NAMES
    colnames(tpm) <- colnames(counts_filtered)
    rownames(tpm) <- rownames(counts_filtered)
    
    return(tpm)
}

# Calculate TPM matrix using the robust Effective Length method
tpm_matrix <- calculate_tpm(
    counts = countData_deseq, 
    featureLength = gene_lengths_filtered$Length,
    meanFragmentLength = MEAN_FRAGMENT_LENGTHS # Uses the new configuration variable
)

# Store dimensions for report
report_metrics$tpm_matrix_rows <- nrow(tpm_matrix)
report_metrics$tpm_matrix_cols <- ncol(tpm_matrix)

# Transform the TPM matrix to log2(TPM+1) for Exploratory Data Analysis (EDA)
# Adding +1 prevents issues with log(0).
log2_tpm_matrix <- log2(tpm_matrix + 1)

cat("TPM matrix calculated using Effective Length and log2-transformed (log2(TPM+1)) for EDA.\n")


# ==============================================================================
# 5. VST/RLOG TRANSFORMATION AND OUTLIER FILTERING (FOR EDA)
# ==============================================================================

cat("\n==============================================================================")
cat("\n5. VST/RLOG TRANSFORMATION AND OUTLIER FILTERING (FOR EDA)")
cat("\n==============================================================================\n")

# Use VST (Variance Stabilizing Transformation) for large datasets (recommended)
# Alternatively, rlog is better for small datasets (less than 30 samples)
if (ncol(countData_deseq) < 30) {
    cat("Using rlog transformation (recommended for smaller sample sizes).\n")
    vst <- rlog(dds, blind=TRUE)
} else {
    cat("Using VST transformation (recommended for larger sample sizes).\n")
    vst <- vst(dds, blind=TRUE)
}

# Extract transformed counts matrix
vst_counts <- assay(vst)

# Store VST object for later use
saveRDS(vst, file.path(dirname(COUNT_FILE), "vst_object.rds"))

# --- Outlier Filtering based on TPM (Optional but recommended) ---
# Filter out genes that have extreme expression values across all samples (IQR outlier rule)
# This uses the TPM matrix, which is better for filtering based on biological abundance.
# The IQR method identifies genes with expression patterns that are highly inconsistent.

# Define the matrix for outlier analysis (TPM data)
outlier_matrix <- log2_tpm_matrix

# Calculate row-wise Interquartile Range (IQR) and quartiles
iqr_vals <- apply(outlier_matrix, 1, IQR)
q1_vals <- apply(outlier_matrix, 1, quantile, 0.25)
q3_vals <- apply(outlier_matrix, 1, quantile, 0.75)

# Define outlier bounds (1.5 * IQR)
lower_bound <- q1_vals - 1.5 * iqr_vals
upper_bound <- q3_vals + 1.5 * iqr_vals

# Check if any sample expression falls outside the bounds for each gene
is_outlier <- sapply(1:nrow(outlier_matrix), function(i) {
    any(outlier_matrix[i, ] < lower_bound[i] | outlier_matrix[i, ] > upper_bound[i])
})

# Genes to keep (non-outliers)
genes_to_keep <- rownames(outlier_matrix)[!is_outlier]
genes_removed <- rownames(outlier_matrix)[is_outlier]

# Apply filtering to the VST counts matrix
vst_counts_filtered <- vst_counts[genes_to_keep, ]

# Store outlier metrics for report
report_metrics$outlier_analysis$total_genes_before_outlier_removal <- nrow(outlier_matrix)
report_metrics$outlier_analysis$outlier_genes_removed <- length(genes_removed)
report_metrics$outlier_analysis$genes_after_outlier_removal <- nrow(vst_counts_filtered)
report_metrics$outlier_analysis$percentage_removed <- round(length(genes_removed) / nrow(outlier_matrix) * 100, 2)

cat(paste("Genes retained after outlier filtering:", nrow(vst_counts_filtered), " (", report_metrics$outlier_analysis$percentage_removed, "% removed)\n"))


# ==============================================================================
# 6. PRINCIPAL COMPONENT ANALYSIS (PCA) PLOT)
# ==============================================================================

cat("\n==============================================================================")
cat("\n6. PRINCIPAL COMPONENT ANALYSIS (PCA) PLOT)")
cat("\n==============================================================================\n")

# Use DESeq2's plotPCA function but extract the data for ggplot2 customization
pcaData <- plotPCA(vst, intgroup=c("Condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Generate the PCA Plot with custom aesthetics 
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
    geom_point(size=3) +
    geom_text(aes(label=row.names(coldata)), vjust=2, hjust=0.5, size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    scale_color_viridis(discrete = TRUE) +
    theme_minimal(base_size = 14) +
    labs(title = "PCA of RNA-Seq Samples (VST Transformed Counts)")

# Save the plot
pca_plot_file <- file.path(dirname(COUNT_FILE), "PCA_plot.png")
ggsave(pca_plot_file, plot = pca_plot, width = 8, height = 7, dpi = 300)
cat(paste("PCA plot saved to:", pca_plot_file, "\n"))


# ==============================================================================
# 7. SAMPLE-TO-SAMPLE DISTANCE HEATMAP
# ==============================================================================

cat("\n==============================================================================")
cat("\n7. SAMPLE-TO-SAMPLE DISTANCE HEATMAP")
cat("\n==============================================================================\n")

# Calculate Euclidean distance between samples using VST counts
sampleDists <- dist(t(vst_counts))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Condition, vst$ShortID, sep="-")
colnames(sampleDistMatrix) <- NULL

# Define colors for the heatmap
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Generate the Heatmap 
pheatmap_plot_file <- file.path(dirname(COUNT_FILE), "Sample_Distance_Heatmap.png")
png(pheatmap_plot_file, width = 800, height = 700) # Open PNG device
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main="Sample-to-Sample Distance Heatmap (VST Counts)")
dev.off() # Close PNG device
cat(paste("Sample distance heatmap saved to:", pheatmap_plot_file, "\n"))


# ==============================================================================
# 8. DIFFERENTIAL EXPRESSION ANALYSIS (DEG) - EXTRACTION
# ==============================================================================

cat("\n==============================================================================")
cat("\n8. DIFFERENTIAL EXPRESSION ANALYSIS (DEG) - EXTRACTION")
cat("\n==============================================================================\n")

ALPHA_THRESHOLD <- 0.05
LOG2FC_THRESHOLD <- 1.0

report_metrics$alpha <- ALPHA_THRESHOLD
report_metrics$log2FC_threshold <- LOG2FC_THRESHOLD

# --- CRITICAL: Identify Coef Names for lfcShrink ---
# These are the names of the coefficients in the design formula (~ Condition)
coef_names <- resultsNames(dds)
cat("Identified DESeq2 Coefficients:\n")
print(coef_names)

# Function to extract results, apply LFC shrinkage, count DEGs, and save CSV
# *** MODIFIED TO USE 'coef' FOR APEGLM SHINKAGE ***
process_results <- function(dds_obj, contrast_levels, alpha, log2FC, name) {
    # 1. Extract results using contrast (necessary to define the comparison)
    res <- results(dds_obj, contrast = c("Condition", contrast_levels), alpha = alpha)
    
    # 2. Derive the COEFFICIENT NAME for the lfcShrink function
    # The coefficient name is generally the last part of the contrast
    # e.g., Low vs Control -> coef name is 'Condition_Low_vs_Control'
    # Check if the dds_obj has a simple design or a combined one
    if (name == "High_vs_Other") {
        # Handles the High vs (Low + Control) comparison with custom factor
        coef_to_shrink <- "CombinedCondition_High_vs_Other" 
        res <- results(dds_obj, contrast = c("CombinedCondition", "High", "Other"), alpha = alpha)
    } else {
        # Standard comparisons (Level vs Control)
        # Assuming Control is the reference level (relevel done in Section 3)
        coef_to_shrink <- paste0("Condition_", contrast_levels[1], "_vs_", contrast_levels[2])
    }

    # Fallback check for standard condition names
    if (!coef_to_shrink %in% resultsNames(dds_obj)) {
        stop(paste("ERROR: Coefficient name not found:", coef_to_shrink, "Available names:", paste(resultsNames(dds_obj), collapse=", ")))
    }
    
    # 3. Apply LFC shrinkage using the COEFFICIENT NAME
    # The 'res' object from step 1 is NOT passed to apeglm, only the coef is needed.
    res_lfc <- lfcShrink(dds_obj, coef = coef_to_shrink, type = "apeglm")
    
    # 4. Filter for DEGs
    res_lfc_filtered <- subset(res_lfc, padj < alpha & abs(log2FoldChange) >= log2FC)
    
    # 5. Count DEGs
    num_degs <- nrow(res_lfc_filtered)
    
    # 6. Save results
    results_file <- file.path(dirname(COUNT_FILE), paste0("DEGs_", name, "_results.csv"))
    write.csv(as.data.frame(res_lfc), file = results_file, row.names = TRUE)
    
    cat(paste0(name, " (", contrast_levels[1], " vs ", contrast_levels[2], "): ", num_degs, " DEGs found. Results saved to ", results_file, "\n"))
    return(list(res_lfc = res_lfc, num_degs = num_degs))
}

# --- Comparison 1: Low vs Control ---
# Coef Name: Condition_Low_vs_Control
results_Low_vs_Control <- process_results(dds, c("Low", "Control"), ALPHA_THRESHOLD, LOG2FC_THRESHOLD, "Low_vs_Control")
report_metrics$DEGs_Low_vs_Control <- results_Low_vs_Control$num_degs

# --- Comparison 2: High vs Control ---
# Coef Name: Condition_High_vs_Control
results_High_vs_Control <- process_results(dds, c("High", "Control"), ALPHA_THRESHOLD, LOG2FC_THRESHOLD, "High_vs_Control")
report_metrics$DEGs_High_vs_Control <- results_High_vs_Control$num_degs

# --- Comparison 3: High vs Low ---
# This comparison is NOT a simple 'vs Control'. We need to extract results first.
# Coef Name for this specific comparison must be defined *manually* from the results() call.
# Since we cannot use 'coef' for non-reference comparisons, we use 'contrast' in lfcShrink, 
# but change the 'type' to 'normal' (non-apeglm). This is a standard compromise.
cat("\n--- Comparison 3 (High vs Low): Using normal shrinkage (type='normal') as it is not vs the reference level. ---\n")
res_High_vs_Low <- results(dds, contrast=c("Condition", "High", "Low"), alpha = ALPHA_THRESHOLD)
res_High_vs_Low_lfc <- lfcShrink(dds, contrast=c("Condition", "High", "Low"), res=res_High_vs_Low, type="normal")

res_lfc_filtered <- subset(res_High_vs_Low_lfc, padj < ALPHA_THRESHOLD & abs(log2FoldChange) >= LOG2FC_THRESHOLD)
report_metrics$DEGs_High_vs_Low <- nrow(res_lfc_filtered)
write.csv(as.data.frame(res_High_vs_Low_lfc), file = file.path(dirname(COUNT_FILE), "DEGs_High_vs_Low_results.csv"), row.names = TRUE)
cat(paste0("High_vs_Low (High vs Low): ", report_metrics$DEGs_High_vs_Low, " DEGs found. Results saved to DEGs_High_vs_Low_results.csv\n"))


# --- Comparison 4: High vs (Low + Control) combined (Good for complex designs) ---
cat("\n--- Comparison 4 (High vs Low + Control Combined): Re-running DESeq with combined factor. ---\n")
dds_combined <- dds
dds_combined$CombinedCondition <- factor(ifelse(dds_combined$Condition == "High", "High", "Other"), levels = c("Other", "High"))
design(dds_combined) <- ~ CombinedCondition
dds_combined <- DESeq(dds_combined)

# Coef Name: CombinedCondition_High_vs_Other
results_High_vs_LowControlCombined <- process_results(dds_combined, c("High", "Other"), ALPHA_THRESHOLD, LOG2FC_THRESHOLD, "High_vs_Other")
report_metrics$DEGs_High_vs_LowControlCombined <- results_High_vs_LowControlCombined$num_degs

# Use the LFC shrunk results for plotting
res_Low_vs_Control_lfc <- results_Low_vs_Control$res_lfc
res_High_vs_Control_lfc <- results_High_vs_Control$res_lfc
# res_High_vs_Low_lfc is already calculated above for plotting
# res_High_vs_Other_lfc (from dds_combined) is results_High_vs_LowControlCombined$res_lfc


# ==============================================================================
# 9. MA PLOTS AND VOLCANO PLOTS
# ==============================================================================

cat("\n==============================================================================")
cat("\n9. MA PLOTS AND VOLCANO PLOTS")
cat("\n==============================================================================\n")

# Function to generate MA plot 
plot_ma <- function(res_lfc, name) {
    # Convert DESeqResults to data frame
    df <- as.data.frame(res_lfc)
    df$sig <- ifelse(df$padj < ALPHA_THRESHOLD & abs(df$log2FoldChange) >= LOG2FC_THRESHOLD, "Significant", "Non-Significant")
    df$sig <- factor(df$sig, levels = c("Non-Significant", "Significant"))
    
    ma_plot <- ggplot(df, aes(x = baseMean, y = log2FoldChange, color = sig)) +
        geom_point(alpha = 0.6, size = 1.5) +
        scale_x_log10() +
        geom_hline(yintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), linetype = "dashed", color = "gray50") +
        geom_hline(yintercept = 0, color = "black", size = 0.5) +
        scale_color_manual(values = c("Non-Significant" = "gray50", "Significant" = "red")) +
        theme_minimal(base_size = 14) +
        labs(
            title = paste("MA Plot:", name),
            x = "Mean of Normalized Counts (log10 scale)",
            y = "log2(Fold Change)"
        ) +
        theme(legend.position = "bottom")

    # Save the plot
    ma_plot_file <- file.path(dirname(COUNT_FILE), paste0("MA_Plot_", name, ".png"))
    ggsave(ma_plot_file, plot = ma_plot, width = 8, height = 7, dpi = 300)
    cat(paste("MA plot saved to:", ma_plot_file, "\n"))
    return(ma_plot)
}

# Function to generate Volcano plot 
plot_volcano <- function(res_lfc, name) {
    df <- as.data.frame(res_lfc)
    df$sig <- "Non-Significant"
    df$sig[df$padj < ALPHA_THRESHOLD & df$log2FoldChange > LOG2FC_THRESHOLD] <- "Upregulated"
    df$sig[df$padj < ALPHA_THRESHOLD & df$log2FoldChange < -LOG2FC_THRESHOLD] <- "Downregulated"
    df$sig <- factor(df$sig, levels = c("Non-Significant", "Upregulated", "Downregulated"))
    
    # -log10(p-adjusted)
    df$neg_log10_padj <- -log10(df$padj)
    
    volcano_plot <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj, color = sig)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_hline(yintercept = -log10(ALPHA_THRESHOLD), linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), linetype = "dashed", color = "gray50") +
        scale_color_manual(values = c("Non-Significant" = "gray50", "Upregulated" = "red", "Downregulated" = "blue")) +
        theme_minimal(base_size = 14) +
        labs(
            title = paste("Volcano Plot:", name),
            x = "log2(Fold Change)",
            y = "-log10(Adjusted P-value)"
        ) +
        theme(legend.position = "bottom")

    # Save the plot
    volcano_plot_file <- file.path(dirname(COUNT_FILE), paste0("Volcano_Plot_", name, ".png"))
    ggsave(volcano_plot_file, plot = volcano_plot, width = 8, height = 7, dpi = 300)
    cat(paste("Volcano plot saved to:", volcano_plot_file, "\n"))
    return(volcano_plot)
}

# Generate and save plots for key comparisons
plot_ma(res_Low_vs_Control_lfc, "Low_vs_Control")
plot_volcano(res_Low_vs_Control_lfc, "Low_vs_Control")

plot_ma(res_High_vs_Control_lfc, "High_vs_Control")
plot_volcano(res_High_vs_Control_lfc, "High_vs_Control")

plot_ma(res_High_vs_Low_lfc, "High_vs_Low")
plot_volcano(res_High_vs_Low_lfc, "High_vs_Low")


# ==============================================================================
# 10. HEATMAP OF TOP DIFFERENTIALY EXPRESSED GENES (DEGs)
# ==============================================================================

cat("\n==============================================================================")
cat("\n10. HEATMAP OF TOP DIFFERENTIALY EXPRESSED GENES (DEGs)")
cat("\n==============================================================================\n")

TOP_GENES_TO_SHOW <- 50 # Number of top DEGs to visualize

# Combine all unique DEGs from the main three contrasts
all_degs <- unique(c(
    rownames(subset(res_Low_vs_Control_lfc, padj < ALPHA_THRESHOLD & abs(log2FoldChange) >= LOG2FC_THRESHOLD)),
    rownames(subset(res_High_vs_Control_lfc, padj < ALPHA_THRESHOLD & abs(log2FoldChange) >= LOG2FC_THRESHOLD)),
    rownames(subset(res_High_vs_Low_lfc, padj < ALPHA_THRESHOLD & abs(log2FoldChange) >= LOG2FC_THRESHOLD))
))

# Calculate overall ranking (e.g., based on mean LFC across all contrasts)
# Using 'res_High_vs_Control_lfc' as the primary ranking for simplicity
# This ensures we get a list of genes ranked by the magnitude of change in one contrast.
ranked_genes <- rownames(res_High_vs_Control_lfc[order(abs(res_High_vs_Control_lfc$log2FoldChange), decreasing = TRUE), ])

# Select top genes based on ranking and existence in all DEGs
# Intersecting with all_degs ensures we only select genes that passed the DEG threshold
top_degs <- intersect(ranked_genes, all_degs)
if (length(top_degs) > TOP_GENES_TO_SHOW) {
    top_degs <- top_degs[1:TOP_GENES_TO_SHOW]
}

if (length(top_degs) > 1) {
    # Ensure all selected genes exist in the filtered VST matrix
    top_degs_final <- intersect(top_degs, rownames(vst_counts_filtered))
    
    if (length(top_degs_final) < 2) {
        cat("WARNING: Not enough unique DEGs remaining after VST/outlier filtering to create a meaningful heatmap.\n")
    } else {
        # Select VST counts for top DEGs
        heatmap_matrix <- vst_counts_filtered[top_degs_final, ]
        
        # Scale (Z-score) the matrix row-wise (gene-wise)
        # Z-score normalization makes the data comparable across different genes
        z_score_matrix <- t(apply(heatmap_matrix, 1, scale))
        colnames(z_score_matrix) <- colnames(heatmap_matrix)
        
        # Setup annotation for samples (top row of the heatmap)
        annotation_col <- data.frame(Condition = coldata$Condition)
        rownames(annotation_col) <- rownames(coldata)
        
        # Setup color palette
        palette_length <- 100
        my_colors <- colorRampPalette(c("blue", "white", "red"))(palette_length)
        
        # Generate the Heatmap using pheatmap 
        pheatmap_deg_file <- file.path(dirname(COUNT_FILE), paste0("Heatmap_Top_", length(top_degs_final), "_DEGs.png"))
        png(pheatmap_deg_file, width = 1000, height = 1400) # Adjust dimensions for better visibility
        pheatmap(z_score_matrix,
                 color = my_colors,
                 clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",
                 clustering_method = "ward.D2",
                 show_rownames = TRUE,
                 show_colnames = TRUE,
                 annotation_col = annotation_col,
                 main = paste("Heatmap of Top", length(top_degs_final), "DEGs (Z-score of VST Counts)")
        )
        dev.off() # Close PNG device
        cat(paste("Heatmap of top DEGs saved to:", pheatmap_deg_file, "\n"))
    }
} else {
    cat("WARNING: Not enough DEGs found across all contrasts to generate a heatmap.\n")
}


# ==============================================================================
# 11. ANNOTATION AND FINAL REPORT GENERATION
# ==============================================================================

cat("\n==============================================================================")
cat("\n11. ANNOTATION AND FINAL REPORT GENERATION")
cat("\n==============================================================================\n")

# Load annotation data
if (file.exists(ANNOTATION_FILE)) {
    # Assuming the annotation file has columns for GeneID and Description
    gene_annotation <- read.delim(ANNOTATION_FILE, header = TRUE, stringsAsFactors = FALSE)
    # Ensure the GeneID column matches the row names of our results (e.g., 'GeneID')
    if (!("GeneID" %in% colnames(gene_annotation))) {
        colnames(gene_annotation)[1] <- "GeneID"
    }
    
    # Function to add annotation to results
    add_annotation <- function(res_lfc) {
        res_df <- as.data.frame(res_lfc)
        res_df$GeneID <- rownames(res_df)
        # Merge by GeneID and keep all results (all.x=TRUE)
        res_df_annotated <- merge(res_df, gene_annotation, by = "GeneID", all.x = TRUE)
        rownames(res_df_annotated) <- res_df_annotated$GeneID
        return(res_df_annotated)
    }

    # Annotate and overwrite final results CSV files
    write.csv(add_annotation(res_Low_vs_Control_lfc), file = file.path(dirname(COUNT_FILE), "DEGs_Low_vs_Control_results_annotated.csv"), row.names = FALSE)
    write.csv(add_annotation(res_High_vs_Control_lfc), file = file.path(dirname(COUNT_FILE), "DEGs_High_vs_Control_results_annotated.csv"), row.names = FALSE)
    write.csv(add_annotation(res_High_vs_Low_lfc), file = file.path(dirname(COUNT_FILE), "DEGs_High_vs_Low_results_annotated.csv"), row.names = FALSE)
    
    cat("Results files have been annotated and saved with '_annotated.csv' suffix.\n")
} else {
    cat("WARNING: Annotation file not found. Skipping result annotation.\n")
}

# --- Generate Final Report ---
report_content <- paste0(
    "# RNA-SEQ ANALYSIS REPORT (DESeq2) - ", Sys.time(), "\n\n",
    "==========================================================================\n\n",
    "*** 1. DATA AND FILTERING SUMMARY ***\n",
    "-------------------------------------\n",
    "Total Genes in Raw Matrix:              ", report_metrics$initial_genes, "\n",
    "Filtering Criteria:                     ", report_metrics$filtering_criteria, "\n",
    "Genes Retained for DESeq2 Analysis:     ", report_metrics$filtered_genes, "\n",
    "Genes Filtered (Low Count Noise):       ", report_metrics$initial_genes - report_metrics$filtered_genes, "\n\n",
    
    "*** 2. TPM AND OUTLIER FILTERING ***\n",
    "-------------------------------------\n",
    "TPM Matrix Dimension:                   ", report_metrics$tpm_matrix_rows, " rows x ", report_metrics$tpm_matrix_cols, " cols\n",
    "Genes before Outlier Removal (TPM IQR): ", report_metrics$outlier_analysis$total_genes_before_outlier_removal, "\n",
    "Genes Identified as Outliers:           ", report_metrics$outlier_analysis$outlier_genes_removed, "\n",
    "Genes Remaining after Outlier Removal:  ", report_metrics$outlier_analysis$genes_after_outlier_removal, "\n",
    "Percentage Removed:                     ", report_metrics$outlier_analysis$percentage_removed, "%\n\n",
    
    "*** 3. DIFFERENTIAL EXPRESSION PARAMETERS ***\n",
    "----------------------------------------------\n",
    "Adjusted P-value Threshold (Alpha):     ", report_metrics$alpha, "\n",
    "Fold Change Threshold (Log2FC):         ", report_metrics$log2FC_threshold, "\n\n",
    
    "*** 4. RESULTS (DIFFERENTIALLY EXPRESSED GENES - DEGs) ***\n",
    "----------------------------------------------------------\n",
    "DEGs (Low vs Control):                  ", report_metrics$DEGs_Low_vs_Control, "\n",
    "DEGs (High vs Control):                 ", report_metrics$DEGs_High_vs_Control, "\n",
    "DEGs (High vs Low):                     ", report_metrics$DEGs_High_vs_Low, "\n",
    "DEGs (High vs (Low + Control)):         ", report_metrics$DEGs_High_vs_LowControlCombined, "\n\n",
    
    "*** 5. GENERATED PLOTS ***\n",
    "--------------------------\n",
    "- PCA Plot (PCA_plot.png)\n",
    "- Sample Distance Heatmap (Sample_Distance_Heatmap.png)\n",
    "- MA Plots for all contrasts (MA_Plot_*.png)\n",
    "- Volcano Plots for all contrasts (Volcano_Plot_*.png)\n",
    "- Heatmap of Top DEGs (Heatmap_Top_*.png)\n\n",
    
    "----------------------------------------------------------\n",
    "All detailed results (DEGs) and plots are saved in the output directory.\n",
    "----------------------------------------------------------\n"
)

# Save report to file
writeLines(report_content, REPORT_FILE)
cat(paste("Analysis Report saved to:", REPORT_FILE, "\n"))

cat("\n==============================================================================")
cat("\nFULL ANALYSIS SCRIPT COMPLETE.")
cat("\n==============================================================================")
