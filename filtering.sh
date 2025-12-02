#!/bin/bash

# ==============================================================================
# BATCH TRIMMING PIPELINE: Trimmomatic (Quality/Illumina) -> Cutadapt (Poly-A/G)
# ==============================================================================

# --- Path and Parameter Configuration ---
INPUT_DIR="input"
# Intermediate output directory (Trimmomatic results, paired)
TRIMMO_OUT_DIR="output/trimmed_trimmomatic"
# Final output directory (Cutadapt results, cleaned and renamed)
FINAL_OUT_DIR="output/final_trimmed_reads"

METADATA_FILE="config/samples.txt"
ADAPTER_FASTA="config/adapters.fasta" # Should contain Illumina adapters only.

THREADS=25 # Number of threads to use (Adjust based on your server resources)

# Create output directories if they do not exist
mkdir -p $TRIMMO_OUT_DIR
mkdir -p $FINAL_OUT_DIR

# --- Processing Functions ---

# NOTE: Assuming 'saveCommand' is a user-defined function for logging commands.

# Function to run Trimmomatic (Quality and Illumina Adapters)
run_trimmomatic() {
    local R1_IN=$1
    local R2_IN=$2
    local R1_PAIRED_OUT=$3
    local R2_PAIRED_OUT=$4
    local R1_UNPAIRED_OUT=$5
    local R2_UNPAIRED_OUT=$6
    local BASE_NAME=$7

    echo "  -> 1/2. Running Trimmomatic (Quality/Illumina)..."

    saveCommand trimmomatic PE -threads $THREADS -phred33 \
        -trimlog "${TRIMMO_OUT_DIR}/${BASE_NAME}.trimmomatic.log" \
        $R1_IN \
        $R2_IN \
        $R1_PAIRED_OUT \
        $R1_UNPAIRED_OUT \
        $R2_PAIRED_OUT \
        $R2_UNPAIRED_OUT \
        ILLUMINACLIP:$ADAPTER_FASTA:2:30:10:8:true \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:20 \
        MINLEN:36
}

# Function to run Cutadapt (Poly-A/G Removal and Final Renaming)
run_cutadapt() {
    local R1_TRIMMO=$1
    local R2_TRIMMO=$2
    local R1_FINAL=$3
    local R2_FINAL=$4

    echo "  -> 2/2. Running Cutadapt (Poly-A/G Removal)..."

    # Cutadapt Parameters for Poly-A/G/C/T cleanup:
    # -a "A{100}" / -A "T{100}": Remove Poly-A/T from 3' ends (R1/R2)
    # -g "G{100}" / -G "C{100}": Remove Poly-G/C from 5' ends (R1/R2)
    saveCommand cutadapt -j $THREADS -pair-adapters \
        -a "A{100}" -A "T{100}" \
        -g "G{100}" -G "C{100}" \
        --minimum-length 36 \
        -o $R1_FINAL -p $R2_FINAL \
        $R1_TRIMMO $R2_TRIMMO
}


# --- Metadata Processing Loop ---
# Uses 'awk' to safely read the SRR number and Sample Name from config/samples.txt
awk 'NR > 1 {print $1, $2}' $METADATA_FILE | while read SRR_NUMBER SAMPLE_NAME
do
    echo "====================================================="
    echo "Starting: SRR ${SRR_NUMBER} | Sample: ${SAMPLE_NAME}"
    echo "====================================================="

    # --- 1. Define File Paths ---
    R1_IN="${INPUT_DIR}/${SRR_NUMBER}_1.fastq"
    R2_IN="${INPUT_DIR}/${SRR_NUMBER}_2.fastq"
    
    # Trimmomatic Paired Output (Intermediate)
    R1_TRIMMO_P="${TRIMMO_OUT_DIR}/${SRR_NUMBER}_1_paired.fastq.gz"
    R2_TRIMMO_P="${TRIMMO_OUT_DIR}/${SRR_NUMBER}_2_paired.fastq.gz"
    # Trimmomatic Unpaired Output (Intermediate, not processed further)
    R1_TRIMMO_U="${TRIMMO_OUT_DIR}/${SRR_NUMBER}_1_unpaired.fastq.gz"
    R2_TRIMMO_U="${TRIMMO_OUT_DIR}/${SRR_NUMBER}_2_unpaired.fastq.gz"
    
    # Final Cutadapt Output (Cleaned and Renamed)
    R1_FINAL="${FINAL_OUT_DIR}/${SAMPLE_NAME}_${SRR_NUMBER}_1_paired.fastq.gz"
    R2_FINAL="${FINAL_OUT_DIR}/${SAMPLE_NAME}_${SRR_NUMBER}_2_paired.fastq.gz"

    # Check for input files
    if [[ ! -f "$R1_IN" || ! -f "$R2_IN" ]]; then
        echo "ERROR: Input files ($R1_IN and $R2_IN) not found. Skipping sample."
        continue
    fi

    # 1. Execute Trimmomatic (Quality Trimming and Illumina Adapter Removal)
    run_trimmomatic "$R1_IN" "$R2_IN" "$R1_TRIMMO_P" "$R1_TRIMMO_U" "$R2_TRIMMO_P" "$R2_TRIMMO_U" "$SRR_NUMBER"
    
    if [ $? -ne 0 ]; then
        echo "FATAL ERROR in Trimmomatic. Stopping sample processing."
        continue
    fi

    # 2. Execute Cutadapt (Poly-A/G Removal and Final Renaming)
    # Reads the paired files from Trimmomatic and writes the final renamed files.
    run_cutadapt "$R1_TRIMMO_P" "$R2_TRIMMO_P" "$R1_FINAL" "$R2_FINAL"

    if [ $? -ne 0 ]; then
        echo "FATAL ERROR in Cutadapt. Stopping sample processing."
        continue
    fi
    
    echo "Processing COMPLETE for ${SAMPLE_NAME}. Final files saved to ${FINAL_OUT_DIR}/"
    echo ""

done

echo "--- All samples processing finished. ---"