#!/bin/bash

# ==============================================================================
# SCRIPT 03: Cleanup (Remove initial FASTQ files)
# Deve ser executado APÓS o Cutadapt (Script 02).
# ==============================================================================

# --- Path Configuration ---
INPUT_DIR="input"
FINAL_OUT_DIR="output/final_trimmed_reads"
METADATA_FILE="config/samples.txt"

echo "--- Starting Smart Cleanup: Deleting initial FASTQ files ---"

CLEANUP_COUNT=0

# Loop através de todas as amostras
# Lê as colunas SRR e Nome_da_Amostra
awk 'NR > 1 {print $1, $2}' $METADATA_FILE | while read SRR_NUMBER SAMPLE_NAME
do
    echo "Checking sample: ${SRR_NUMBER}"

    # Arquivos originais descompactados que devem ser removidos
    R1_IN="${INPUT_DIR}/${SRR_NUMBER}_1.fastq"
    R2_IN="${INPUT_DIR}/${SRR_NUMBER}_2.fastq"
    
    # Arquivos finais filtrados e compactados (Verificação de segurança)
    R1_FINAL_EXPECTED="${FINAL_OUT_DIR}/${SAMPLE_NAME}_${SRR_NUMBER}_1_paired.fastq.gz"
    R2_FINAL_EXPECTED="${FINAL_OUT_DIR}/${SAMPLE_NAME}_${SRR_NUMBER}_2_paired.fastq.gz"

    # --- Lógica de Verificação de Segurança ---
    
    # 1. Checa se as cópias finais JÁ EXISTEM
    if [[ -f "$R1_FINAL_EXPECTED" && -f "$R2_FINAL_EXPECTED" ]]; then
        
        # 2. Checa se os arquivos originais AINDA EXISTEM no INPUT_DIR para remoção
        if [[ -f "$R1_IN" && -f "$R2_IN" ]]; then
            
            echo "    -> SUCCESS: Final files exist. Removing originals for ${SRR_NUMBER}."
            
            # Comando de remoção (use 'rm -f' para forçar e evitar prompts)
            rm -f "$R1_IN"
            rm -f "$R2_IN"

            if [ $? -eq 0 ]; then
                echo "    -> Removed: $R1_IN and $R2_IN"
                CLEANUP_COUNT=$((CLEANUP_COUNT + 2))
            else
                echo "    -> WARNING: Failed to remove files for ${SRR_NUMBER}."
            fi

        elif [[ -f "$R1_IN" || -f "$R2_IN" ]]; then
            # Caso em que um arquivo existe e o outro não (Limpeza parcial)
            echo "    -> INFO: Removing any remaining file in $INPUT_DIR for ${SRR_NUMBER}."
            rm -f "$R1_IN" "$R2_IN"
            
        else
            echo "    -> SKIP: Original FASTQ files already removed for ${SRR_NUMBER}."
        fi

    else
        echo "    -> WARNING: Final processed files for ${SRR_NUMBER} NOT FOUND in $FINAL_OUT_DIR."
        echo "    -> ABORTING REMOVAL for this sample to prevent data loss."
    fi
    echo "-----------------------------------------------------"
done

echo "--- Cleanup Complete ---"
echo "Total files successfully removed from '$INPUT_DIR': $CLEANUP_COUNT"
echo "All quality-controlled reads remain in '$FINAL_OUT_DIR'."