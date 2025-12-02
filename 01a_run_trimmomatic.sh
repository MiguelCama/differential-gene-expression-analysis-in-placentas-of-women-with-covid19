#!/bin/bash

# ==============================================================================
# SCRIPT 01: Trimmomatic (Quality Control and Illumina Adapter Removal)
# ==============================================================================

# --- Path and Parameter Configuration ---
INPUT_DIR="input"
TRIMMO_OUT_DIR="output/trimmed_trimmomatic"
METADATA_FILE="config/samples.txt"
ADAPTER_FASTA="config/adapters.fasta" 

THREADS=50

mkdir -p $TRIMMO_OUT_DIR

# Function to run Trimmomatic 
run_trimmomatic() {
    local R1_IN=$1
    local R2_IN=$2
    local BASE_NAME=$3 
    
    local R1_PAIRED_OUT="${TRIMMO_OUT_DIR}/${BASE_NAME}_1_paired.fastq.gz"
    local R2_PAIRED_OUT="${TRIMMO_OUT_DIR}/${BASE_NAME}_2_paired.fastq.gz"
    local R1_UNPAIRED_OUT="${TRIMMO_OUT_DIR}/${BASE_NAME}_1_unpaired.fastq.gz"
    local R2_UNPAIRED_OUT="${TRIMMO_OUT_DIR}/${BASE_NAME}_2_unpaired.fastq.gz"

    echo "    -> Running Trimmomatic for ${BASE_NAME}..."

    # Comando Trimmomatic
    saveCommand trimmomatic PE -threads $THREADS -phred33 \
        -trimlog "${TRIMMO_OUT_DIR}/${BASE_NAME}.trimmomatic.log" \
        "$R1_IN" \
        "$R2_IN" \
        "$R1_PAIRED_OUT" \
        "$R1_UNPAIRED_OUT" \
        "$R2_PAIRED_OUT" \
        "$R2_UNPAIRED_OUT" \
        ILLUMINACLIP:"$ADAPTER_FASTA":2:30:10:8:true \
        MINLEN:36
}


# --- Processing Loop ---
# O awk lê a coluna 1 do samples.txt
awk 'NR > 1 {print $1}' $METADATA_FILE | while read SRR_NUMBER
do
    echo "====================================================="
    echo "Checking and running Trimmomatic for SRR ${SRR_NUMBER}"
    echo "====================================================="

    R1_IN="${INPUT_DIR}/${SRR_NUMBER}_1.fastq"
    R2_IN="${INPUT_DIR}/${SRR_NUMBER}_2.fastq"
    
    R1_PAIRED_OUT_EXPECTED="${TRIMMO_OUT_DIR}/${SRR_NUMBER}_1_paired.fastq.gz"
    R2_PAIRED_OUT_EXPECTED="${TRIMMO_OUT_DIR}/${SRR_NUMBER}_2_paired.fastq.gz"

    # --- 1. Checagem de Arquivos de Saída JÁ EXISTENTES (LINHA 59 CORRIGIDA) ---
    if [[ -f "$R1_PAIRED_OUT_EXPECTED" && -f "$R2_PAIRED_OUT_EXPECTED" ]]; then
        echo "Skipping ${SRR_NUMBER}: Paired output files already exist in ${TRIMMO_OUT_DIR}."
        continue 
    fi

    # --- 2. Checagem de Arquivos de Entrada ---
    if [[ ! -f "$R1_IN" || ! -f "$R2_IN" ]]; then
        echo "ERROR: Input files ($R1_IN and $R2_IN) not found. Skipping sample."
        continue 
    fi

    # --- 3. Execução (Se as saídas não existirem e as entradas existirem) ---
    run_trimmomatic "$R1_IN" "$R2_IN" "$SRR_NUMBER"
    
    if [ $? -ne 0 ]; then
        echo "FATAL ERROR in Trimmomatic for ${SRR_NUMBER}. Continuing to next sample."
    fi
    echo ""

done

echo "--- Trimmomatic processing finished. Output in ${TRIMMO_OUT_DIR}/ ---"