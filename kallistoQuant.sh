#!/bin/bash

# ==============================================================================
# Kallisto Batch Quantification Script
# This script runs 'kallisto quant' for all paired samples
# found in the input reads directory.
# ==============================================================================

# --- Path Configurations ---
# The directory containing the paired and compressed FASTQ files (.fastq.gz)
INPUT_DIR="output/final_trimmed_reads"

# Path to the Kallisto index file generated in the previous step
INDEX_FILE="kallisto/index/kallisto_mRNA_index.idx"

# Directory where the quantification results will be saved
OUTPUT_BASE_DIR="kallisto/quant_results"

# Number of threads to use (adjust based on your system capacity)
THREADS=75

# Number of bootstrap samples (recommended for differential expression analysis)
BOOTSTRAP_SAMPLES=100

# Create the main output directory if it doesn't already exist
mkdir -p "$OUTPUT_BASE_DIR"

echo "Starting Kallisto quantification..."
echo "Index: $INDEX_FILE"
echo "Threads: $THREADS"
echo "Bootstrap Samples: $BOOTSTRAP_SAMPLES"
echo "--------------------------------------------------------"

# 1. Loop over all R1 read files
# The pattern *_1_paired.fastq.gz ensures we only process the initial pair
for R1_FILE in "$INPUT_DIR"/*_1_paired.fastq.gz; do
    
    # Check if the R1 file was found
    if [ ! -f "$R1_FILE" ]; then
        echo "WARNING: No R1 files found in $INPUT_DIR. Check the path."
        break
    fi

    # 2. Extract the base name of the sample (e.g., C_1_SRR15312907)
    # Remove the full path and the suffix '_1_paired.fastq.gz'
    BASE_NAME=$(basename "$R1_FILE" | sed 's/_1_paired.fastq.gz//')
    
    # 3. Build the R2 file name (by replacing '_1' with '_2')
    R2_FILE="${R1_FILE/_1_paired.fastq.gz/_2_paired.fastq.gz}"
    
    # 4. Define the output directory for this specific sample
    SAMPLE_OUTDIR="$OUTPUT_BASE_DIR/$BASE_NAME"

    echo "Processing sample: $BASE_NAME"
    echo "  R1: $(basename "$R1_FILE")"
    echo "  R2: $(basename "$R2_FILE")"
    echo "  Output: $SAMPLE_OUTDIR"
    
    # 5. Check if the R2 file exists
    if [ ! -f "$R2_FILE" ]; then
        echo "FATAL ERROR: Corresponding R2 file ($R2_FILE) not found. Skipping sample."
        continue
    fi
    
    # 6. Execute Kallisto quant
    # We use saveCommand to log the exact command
    saveCommand kallisto quant \
        -i "$INDEX_FILE" \
        -o "$SAMPLE_OUTDIR" \
        -b "$BOOTSTRAP_SAMPLES" \
        -t "$THREADS" \
        "$R1_FILE" \
        "$R2_FILE"
    
    echo "Quantification for $BASE_NAME complete."
    echo "---"

done

echo "--------------------------------------------------------"
echo "Kallisto quantification finished for all samples."
echo "Results available in: $OUTPUT_BASE_DIR"