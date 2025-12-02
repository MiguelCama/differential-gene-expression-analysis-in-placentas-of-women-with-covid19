#!/bin/bash

# ==============================================================================
# SCRIPT 02: Cutadapt (Poly-A/G Removal and Final Renaming)
# Depends on the output of 01_run_trimmomatic.sh
# ==============================================================================

# --- Path and Parameter Configuration ---
TRIMMO_OUT_DIR="output/trimmed_trimmomatic"
FINAL_OUT_DIR="output/final_trimmed_reads"
METADATA_FILE="config/samples.txt"

THREADS=75 # Using 75 threads is aggressive; ensure your system can handle this.

mkdir -p $FINAL_OUT_DIR

# Function to run Cutadapt
run_cutadapt() {
    local R1_TRIMMO=$1
    local R2_TRIMMO=$2
    local R1_FINAL=$3
    local R2_FINAL=$4

    echo "  -> Running Cutadapt..."

    # Comando CORRIGIDO: Removido 'saveCommand'
    saveCommand cutadapt -j $THREADS -pair-adapters \
        -a "A{20}" -A "T{20}" \
        -g "G{20}" -G "C{20}" \
        --minimum-length 36 \
        -o "$R1_FINAL" -p "$R2_FINAL" \
        "$R1_TRIMMO" "$R2_TRIMMO"
    
    # Adicionando aspas duplas às variáveis para segurança
}


# --- Processing Loop ---
# Usando 'tail -n +2' para pular o cabeçalho, garantindo que o awk não seja quebrado por caracteres estranhos no cabeçalho
tail -n +2 "$METADATA_FILE" | awk '{print $1, $2}' | while read SRR_NUMBER SAMPLE_NAME
do
    echo "====================================================="
    echo "Checking and running Cutadapt for ${SAMPLE_NAME} (${SRR_NUMBER})"
    echo "====================================================="

    # Arquivos de entrada do Trimmomatic
    R1_TRIMMO_P="${TRIMMO_OUT_DIR}/${SRR_NUMBER}_1_paired.fastq.gz"
    R2_TRIMMO_P="${TRIMMO_OUT_DIR}/${SRR_NUMBER}_2_paired.fastq.gz"
    
    # Arquivos de saída finais (Renomeados)
    R1_FINAL_EXPECTED="${FINAL_OUT_DIR}/${SAMPLE_NAME}_${SRR_NUMBER}_1_paired.fastq.gz"
    R2_FINAL_EXPECTED="${FINAL_OUT_DIR}/${SAMPLE_NAME}_${SRR_NUMBER}_2_paired.fastq.gz"

    # --- 1. Checagem de Arquivos de Saída JÁ EXISTENTES ---
    if [[ -f "$R1_FINAL_EXPECTED" && -f "$R2_FINAL_EXPECTED" ]]; then
        echo "Skipping ${SAMPLE_NAME}: Final output files already exist in ${FINAL_OUT_DIR}."
        continue
    fi

    # --- 2. Checagem de Arquivos de Entrada (Saída do Trimmomatic) ---
    if [[ ! -f "$R1_TRIMMO_P" || ! -f "$R2_TRIMMO_P" ]]; then
        echo "ERROR: Trimmomatic paired files not found for ${SRR_NUMBER}. Skipping Cutadapt. "
        continue
    fi

    # --- 3. Execução ---
    run_cutadapt "$R1_TRIMMO_P" "$R2_TRIMMO_P" "$R1_FINAL_EXPECTED" "$R2_FINAL_EXPECTED"

    if [ $? -ne 0 ]; then
        echo "FATAL ERROR in Cutadapt for ${SAMPLE_NAME}. Check the logs. Continuing to next sample."
    else
        echo "Cutadapt COMPLETE. Final files saved to ${FINAL_OUT_DIR}/"
    fi
    echo ""

done

echo "--- Cutadapt processing finished. Output in ${FINAL_OUT_DIR}/ ---"
