#!/bin/bash

# ==============================================================================
# SCRIPT 00: Download Raw Data (fastq-dump)
# DESCRIPTION: Downloads paired-end sequencing data from SRA using fastq-dump,
#              generating compressed files (.fastq.gz).
# USAGE: ./run_fastqdump.sh
# ==============================================================================

# --- Path and Parameter Configuration (Must match the Trimmomatic script) ---
INPUT_DIR="input"
METADATA_FILE="config/samples.txt"
THREADS=50 # Number of threads for faster download and processing by SRA Toolkit.

echo "--- Starting raw data download with fastq-dump ---"

# --- 1. Tool Verification ---
if ! command -v fastq-dump &> /dev/null
then
    echo "CRITICAL ERROR: 'fastq-dump' (SRA Toolkit) was not found."
    echo "Please ensure the SRA Toolkit is installed and in your PATH."
    exit 1
fi

# 2. Create output directory
mkdir -p "$INPUT_DIR"
echo "Download destination directory: $INPUT_DIR"
echo "Reading SRR numbers from: $METADATA_FILE"
echo "---"

# --- Processing Loop ---
# awk reads column 1 of samples.txt, ignoring the first line (header)
awk 'NR > 1 {print $1}' "$METADATA_FILE" | while read SRR_NUMBER
do
    SRR_NUMBER=$(echo "$SRR_NUMBER" | tr -d '[:space:]') # Clean up spaces
    
    # Adds a check to ensure SRR_NUMBER is not empty.
    if [ -z "$SRR_NUMBER" ]; then
        continue
    fi

    echo "====================================================="
    echo "Checking and downloading SRR ${SRR_NUMBER}"
    echo "====================================================="

    # Defines expected output file names (now compressed in .fastq.gz)
    R1_OUT="${INPUT_DIR}/${SRR_NUMBER}_1.fastq.gz"
    R2_OUT="${INPUT_DIR}/${SRR_NUMBER}_2.fastq.gz"

    # --- 1. Smart Check for Existing Output Files ---
    # If both compressed files exist, skip the download.
    if [[ -f "$R1_OUT" && -f "$R2_OUT" ]]; then
        echo "Skipping ${SRR_NUMBER}: Both output files (.fastq.gz) already exist in ${INPUT_DIR}."
        continue 
    fi

    # --- 2. fastq-dump Execution ---
    echo "  -> Starting fastq-dump for ${SRR_NUMBER}. Generating .fastq.gz files..."
    
    # Used options:
    # --split-files: Splits paired-end reads into *_1.fastq.gz and *_2.fastq.gz
    # --gzip: Compresses the output
    # --skip-technical: Ignores technical reads
    # --outdir: Specifies the output directory
    
    saveCommand fastq-dump --split-files --gzip --skip-technical \
               --outdir "$INPUT_DIR" \
               --threads "$THREADS" \
               "$SRR_NUMBER"
    
    # Checks the return status
    if [ $? -eq 0 ]; then
        echo "  -> Successful download, splitting, and compression for ${SRR_NUMBER}."
    else
        echo "FATAL ERROR in fastq-dump for ${SRR_NUMBER}. Continuing to next sample."
    fi
    echo ""

done

echo "--- fastq-dump processing finished. Raw data in ${INPUT_DIR}/ ---"
echo "NOTE: The generated files are now compressed (.fastq.gz)."