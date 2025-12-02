#!/bin/bash

# --- 01_rename_srr.sh: Remove o código SRR do nome dos arquivos SAM ---

# --- Configuration ---
ALIGNMENT_DIR="bowtie2/alignment"

echo "--- Iniciando a Renomeação dos Arquivos SAM ---"

# Verifica se o diretório existe
if [ ! -d "$ALIGNMENT_DIR" ]; then
    echo "!!! Diretório de alinhamento não encontrado: $ALIGNMENT_DIR. Abortando. !!!"
    exit 1
fi

# Itera sobre todos os arquivos .sam
find "$ALIGNMENT_DIR" -maxdepth 1 -name "*.sam" | while read FULL_PATH_SAM; do
    # Extrai o nome do arquivo (ex: C_1_SRR15312907.sam)
    FILENAME=$(basename "$FULL_PATH_SAM")
    
    # Define o padrão a ser removido: _SRR e 8 ou mais dígitos até .sam
    # O novo nome será a string antes de '_SRR' mais '.sam'
    # Exemplo: C_1_SRR15312907.sam -> C_1.sam
    
    # Usa o comando 'sed' para remover o padrão '_SRR[dígitos]'
    NEW_FILENAME=$(echo "$FILENAME" | sed -E 's/_[S]RR[0-9]{8,}\.sam/\.sam/')
    
    # Garante que o nome não é o mesmo (para evitar renomear se o padrão não for encontrado)
    if [ "$FILENAME" != "$NEW_FILENAME" ]; then
        NEW_FULL_PATH="${ALIGNMENT_DIR}/${NEW_FILENAME}"
        echo "Renomeando: $FILENAME -> $NEW_FILENAME"
        mv "$FULL_PATH_SAM" "$NEW_FULL_PATH"
    else
        echo "Aviso: Nome de arquivo não alterado (Padrão SRR não encontrado ou erro): $FILENAME"
    fi
done

echo "--- Renomeação Concluída ---"

# Após a execução, os arquivos serão:
# C_1_SRR15312907.sam -> C_1.sam
# P_15_SRR15312906.sam -> P_15.sam