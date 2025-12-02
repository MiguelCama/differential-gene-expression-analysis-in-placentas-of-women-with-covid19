#!/bin/bash

# ==============================================================================
# SCRIPT 03_FIX_RENAME: Corrige nomenclatura de arquivos antigos em FINAL_OUT_DIR
# ==============================================================================

# --- Path and Parameter Configuration ---
FINAL_OUT_DIR="output/final_trimmed_reads"
METADATA_FILE="config/samples.txt"

echo "--- Starting Naming Correction in $FINAL_OUT_DIR ---"

FIXED_COUNT=0

# Loop através de todas as amostras
# Lê as colunas SRR, Nome_da_Amostra e Grupo (o awk lê o grupo, mas não o usamos)
awk 'NR > 1 {print $1, $2, $3}' $METADATA_FILE | while read SRR_NUMBER NEW_SAMPLE_NAME GROUP
do
    # 1. Tenta inferir o NOME ANTIGO baseado no GRUPO e no número SRR (Lógica de Renomeação)

    # Exemplo: SRR15312891 (H_2) foi nomeado como 'Controle_2'
    # Exemplo: SRR15312894 (P_3) foi nomeado como 'Infectada_3'
    
    # Tentativa de Nome Antigo (usamos o sed para remover o prefixo de Controles ou Infectadas)
    # Procuramos o arquivo baseado no SRR, que é a parte constante
    
    # Padrão para encontrar o arquivo antigo (ex: Controle_2_SRR15312891_1_paired.fastq.gz)
    # A parte 'Controle_2_' ou 'Infectada_3_' é o que precisa ser substituído.
    
    # Encontra todos os arquivos que contêm o SRR na pasta final
    FILES_TO_CHECK=$(find "$FINAL_OUT_DIR" -maxdepth 1 -name "*${SRR_NUMBER}_[12]_paired.fastq.gz" 2>/dev/null)

    if [ -z "$FILES_TO_CHECK" ]; then
        # Se não houver arquivos antigos, continua para o próximo
        continue
    fi
    
    # 2. Define o NOME NOVO (Correto, vindo do samples.txt atual)
    NEW_PREFIX="${NEW_SAMPLE_NAME}_${SRR_NUMBER}"
    
    echo "Processing SRR ${SRR_NUMBER}: Target name is ${NEW_SAMPLE_NAME}"

    # 3. Itera sobre os arquivos encontrados e renomeia
    for OLD_FILE_PATH in $FILES_TO_CHECK; do
        
        # Extrai apenas o nome do arquivo
        OLD_FILE_BASENAME=$(basename "$OLD_FILE_PATH")
        
        # Cria o NOVO nome do arquivo: 
        # Mantém apenas o sufixo '_1_paired.fastq.gz' ou '_2_paired.fastq.gz'
        # e substitui o prefixo antigo pelo novo (${NEW_PREFIX}).
        # O 'sed' remove tudo ANTES do SRR, e o 'OLD_PREFIX' é o que sobra
        
        # Ex: OLD_FILE_BASENAME = Controle_2_SRR15312891_1_paired.fastq.gz
        # Novo Sufixo: _1_paired.fastq.gz
        
        # Tentativa de extração do sufixo R1/R2
        # Exemplo: _1_paired.fastq.gz
        SUFFIX=$(echo "$OLD_FILE_BASENAME" | sed 's/.*SRR[0-9]*//') 
        
        # Cria o nome final correto
        NEW_FILE_BASENAME="${NEW_PREFIX}${SUFFIX}"
        NEW_FILE_PATH="${FINAL_OUT_DIR}/${NEW_FILE_BASENAME}"

        # Verifica se o arquivo já está nomeado corretamente
        if [[ "$OLD_FILE_BASENAME" == "${NEW_FILE_BASENAME}" ]]; then
            echo "    -> SKIP: $OLD_FILE_BASENAME already has the correct name. "
            continue
        fi
        
        # Executa a renomeação (mv)
        echo "    -> RENAMING: $OLD_FILE_BASENAME -> $NEW_FILE_BASENAME"
        mv "$OLD_FILE_PATH" "$NEW_FILE_PATH"
        
        if [ $? -eq 0 ]; then
            FIXED_COUNT=$((FIXED_COUNT + 1))
        else
            echo "    -> ERROR: Failed to rename $OLD_FILE_BASENAME."
        fi

    done
    echo "-----------------------------------------------------"
done

echo "--- Renaming Cleanup Complete ---"
echo "Total files successfully renamed: $FIXED_COUNT"
echo "Agora você pode prosseguir com o processamento das amostras faltantes."