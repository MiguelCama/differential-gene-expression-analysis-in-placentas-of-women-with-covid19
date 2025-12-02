#!/bin/bash

# ==============================================================================
# SCRIPT 03: FastQC & MultiQC for TRIMMED Data (Post-Trimmomatic)
# Includes smart checking to avoid reprocessing.
# ==============================================================================

# --- CONDA ENVIRONMENT CONTROL ---
# IMPORTANT: This line assumes you have a conda environment named 'multiqc_env' 
# with FastQC and MultiQC installed.
echo "Activating conda environment: multiqc_env"
conda activate multiqc_env

# --- Path and Parameter Configuration ---
# INPUT DIRECTORY: Onde o Trimmomatic salvou os arquivos emparelhados (.fastq.gz)
INPUT_DIR="output/trimmed_trimmomatic"

# OUTPUT DIRECTORIES
FASTQC_OUTPUT_DIR="reports/fastqc_trimmomatic"
LOG_DIR="log/fastqc_trimmomatic"
 
THREADS=55 # Use o número de threads adequado para sua máquina.

echo "--- Starting Smart FastQC and MultiQC Reporting for TRIMMED Reads ---"

# 1. Create necessary directories
mkdir -p "$FASTQC_OUTPUT_DIR"
mkdir -p "$LOG_DIR"

FASTQC_ERROR=0

# --- 2. FastQC Execution (Smart Check and Process) ---

echo "Checking existing FastQC reports and running only for missing files in $INPUT_DIR..."

# Loop through paired .fastq.gz files in the input folder
# Usamos o padrão *_paired.fastq.gz para garantir que processamos os reads que serão usados.
for file_path in "$INPUT_DIR"/*_paired.fastq.gz; do
    
    # Check if the .fastq.gz file exists (handle case where no files match the pattern)
    if [ ! -f "$file_path" ]; then
        continue
    fi
    
    # Get the base filename (Ex: SRR123456_1_paired.fastq.gz)
    SAMPLE_NAME=$(basename "$file_path")
    
    # 3. Define the expected output report filename (FastQC generates a ZIP)
    # Remove the .fastq.gz extension and add the _fastqc.zip suffix
    EXPECTED_REPORT="${FASTQC_OUTPUT_DIR}/${SAMPLE_NAME%.fastq.gz}_fastqc.zip"
    
    # 4. Check if the report (ZIP) ALREADY EXISTS
    if [ -f "$EXPECTED_REPORT" ]; then
        echo "Skipping $SAMPLE_NAME: FastQC report already found at $(basename "$EXPECTED_REPORT")."
        continue
    fi
    
    # 5. If the report DOES NOT EXIST, processing starts:
    LOG_FILE_PATH="${LOG_DIR}/${SAMPLE_NAME%.fastq.gz}.fastqc.err"
    echo "Processing missing file: $SAMPLE_NAME | Log: $(basename "$LOG_FILE_PATH")"
    
    # Execute FastQC for the sample
    # Note: If 'saveCommand' is a custom function in your workflow, replace it 
    # with the plain command: fastqc
    fastqc -t $THREADS "$file_path" -o "$FASTQC_OUTPUT_DIR" 2> "$LOG_FILE_PATH"
    
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

# The MultiQC output directory will be the parent folder of FASTQC_OUTPUT_DIR (reports/)
REPORT_DIR=$(dirname "$FASTQC_OUTPUT_DIR") # reports/

# MultiQC scans the results from the new fastqc data folder
# Note: If 'saveCommand' is a custom function in your workflow, replace it 
# with the plain command: multiqc
multiqc "$FASTQC_OUTPUT_DIR" -o "$REPORT_DIR" -n "multiqc_trimmomatic_report"

if [ $? -ne 0 ]; then
    echo "ERROR: MultiQC execution failed. Check if MultiQC is installed and accessible."
    exit 1
fi

echo "MultiQC report successfully generated: $REPORT_DIR/multiqc_trimmomatic_report.html"
echo "--- Analysis Workflow Complete ---"

# --- CONDA DEACTIVATE ---
conda deactivate
echo "Conda environment deactivated."
