#!/bin/bash

# ==============================================================================
# SCRIPT 04: FastQC & MultiQC for FILTERED Data (Post-Trimming/Cutadapt)
# Includes smart checking to avoid reprocessing.
# ==============================================================================

# --- CONDA ENVIRONMENT CONTROL ---
echo "Activating conda environment: multiqc_env"
# Ensures FastQC and MultiQC are accessible.
conda activate multiqc_env

# --- Path and Parameter Configuration ---
# INPUT DIRECTORY: processed and renamed files
INPUT_DIR="output/final_trimmed_reads"
FASTQC_OUTPUT_DIR="reports/fastqc_final_data"
LOG_DIR="log/fastqc_final" 
THREADS=75 # Use the maximum number of threads you configured

echo "--- Starting Smart FastQC and MultiQC Reporting for FINAL Reads ---"

# 1. Create necessary directories
mkdir -p "$FASTQC_OUTPUT_DIR"
mkdir -p "$LOG_DIR"

FASTQC_ERROR=0

# --- 2. FastQC Execution (Smart Check and Process) ---

echo "Checking existing FastQC reports and running only for missing files..."

# Loop through ALL .fastq.gz files in the FINAL input folder
for file_path in "$INPUT_DIR"/*.fastq.gz; do
    
    # Check if the .fastq.gz file exists
    if [ ! -f "$file_path" ]; then
        continue
    fi
    
    # Get the base filename (Ex: H_2_SRR15312891_1_paired.fastq.gz)
    SAMPLE_NAME=$(basename "$file_path")
    
    # **THE CHECK WORKS LIKE THIS:**
    # 3. Define the expected output report filename (FastQC generates a ZIP)
    # Remove the .fastq.gz extension and add the _fastqc.zip suffix
    EXPECTED_REPORT="${FASTQC_OUTPUT_DIR}/${SAMPLE_NAME%.fastq.gz}_fastqc.zip"
    
    # 4. Check if the report (ZIP) ALREADY EXISTS in the output directory
    if [ -f "$EXPECTED_REPORT" ]; then
        echo "Skipping $SAMPLE_NAME: FastQC report already found at $(basename "$EXPECTED_REPORT")."
        # If it exists, move to the next file in the loop. (THIS IS THE CHECK YOU NEED)
        continue
    fi
    
    # 5. If the report DOES NOT EXIST, processing starts:
    LOG_FILE_PATH="${LOG_DIR}/${SAMPLE_NAME%.fastq.gz}.fastqc.err"
    echo "Processing missing file: $SAMPLE_NAME | Log: $(basename "$LOG_FILE_PATH")"
    
    # Execute FastQC for the sample
    saveCommand fastqc -t $THREADS "$file_path" -o "$FASTQC_OUTPUT_DIR" 2> "$LOG_FILE_PATH"
    
    # 6. Check if the FastQC command was successful
    if [ $? -ne 0 ]; then
        echo "CRITICAL ERROR: FastQC failed on file $SAMPLE_NAME. Check the log file for details."
        FASTQC_ERROR=1
    fi

done

if [ "$FASTQC_ERROR" -ne 0 ]; then
    echo "--- FastQC finished, but with errors on some files. ---"
else
    echo "--- FastQC successfully completed for all required files. ---"
fi

echo "Reports are saved in: $FASTQC_OUTPUT_DIR"
echo "Individual error logs are saved in: $LOG_DIR/"


# --- 3. MultiQC Execution (Runs AFTER all samples are done) ---

echo ""
echo "Starting MultiQC to generate combined report..."

# The MultiQC output directory will be 'reports/'
REPORT_DIR=$(dirname "$FASTQC_OUTPUT_DIR") # reports/

# MultiQC scans the results from the new fastqc data folder
saveCommand multiqc "$FASTQC_OUTPUT_DIR" -o "$REPORT_DIR" -n "multiqc_final_report"

if [ $? -ne 0 ]; then
    echo "ERROR: MultiQC execution failed. Check if MultiQC is installed and accessible."
    exit 1
fi

echo "MultiQC report successfully generated: $REPORT_DIR/multiqc_final_report.html"
echo "--- Analysis Workflow Complete ---"

# --- CONDA DEACTIVATE ---
conda deactivate
echo "Conda environment deactivated."