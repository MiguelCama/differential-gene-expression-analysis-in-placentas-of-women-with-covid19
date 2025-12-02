#!/bin/bash

# Defines the input and output directories
INPUT_DIR="input"
FASTQC_OUTPUT_DIR="reports/fastqc_data"
LOG_DIR="log"

echo "--- Starting Smart FastQC and MultiQC Reporting ---"

# 1. Create necessary directories
mkdir -p "$FASTQC_OUTPUT_DIR"
mkdir -p "$LOG_DIR"

FASTQC_ERROR=0

# --- 2. FastQC Execution (Check and Process) ---

echo "Checking existing FastQC reports and running only for missing files..."

# Loop through ALL .fastq files in the input folder
for file_path in "$INPUT_DIR"/*.fastq; do
    
    # Check if the .fastq file exists (in case there are no files in INPUT_DIR)
    if [ ! -f "$file_path" ]; then
        continue
    fi
    
    # Get the base file name (Ex: SRR15312890_1.fastq)
    SAMPLE_NAME=$(basename "$file_path")
    
    # **EXISTENCE CHECK:**
    # 3. Define the expected output report file name (FastQC generates a ZIP)
    # Ex: reports/fastqc_data/SRR15312890_1_fastqc.zip
    EXPECTED_REPORT="${FASTQC_OUTPUT_DIR}/${SAMPLE_NAME%.fastq}_fastqc.zip"
    
    # 4. Check if the report (ZIP) ALREADY EXISTS in the output directory
    if [ -f "$EXPECTED_REPORT" ]; then
        echo "Skipping $SAMPLE_NAME: FastQC report already found at $(basename "$EXPECTED_REPORT")."
        # If it exists, proceed to the next file in the loop
        continue
    fi
    
    # 5. If the report DOES NOT EXIST, processing starts:
    LOG_FILE_PATH="${LOG_DIR}/${SAMPLE_NAME}.fastqc.err"
    echo "Processing missing file: $SAMPLE_NAME | Log: $LOG_FILE_PATH"
    
    # Execute FastQC for the sample
    saveCommand fastqc -t 50 "$file_path" -o "$FASTQC_OUTPUT_DIR" 2> "$LOG_FILE_PATH"
    
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


# --- 3. MultiQC Execution (Runs AFTER all samples are finished) ---

echo ""
echo "Starting MultiQC to generate combined report..."

# The MultiQC output directory will be 'reports/'
REPORT_DIR=$(dirname "$FASTQC_OUTPUT_DIR") 

# MultiQC scans the results of all files in the fastqc_data folder
saveCommand multiqc "$FASTQC_OUTPUT_DIR" -o "$REPORT_DIR"

if [ $? -ne 0 ]; then
    echo "ERROR: MultiQC execution failed. Check if MultiQC is installed and accessible."
    exit 1
fi

echo "MultiQC report successfully generated: $REPORT_DIR/multiqc_report.html"
echo "--- Analysis Workflow Complete ---"