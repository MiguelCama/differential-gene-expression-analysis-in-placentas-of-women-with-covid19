#!/bin/bash

# Define o diretório de entrada
INPUT_DIR="input"
REPORT_SCRIPT="script/00b_run_multiqc_pre.sh"

# --- 1. Verificação de Espaço e Cálculo ---

echo "--- Disk Space Verification for Decompression ---"

# 1.1 Calcula o tamanho total (compactado) dos arquivos .fastq.gz
TOTAL_COMPACTED_BYTES=$(du -cb "$INPUT_DIR"/*.fastq.gz 2>/dev/null | tail -n 1 | awk '{print $1}')

# Verifica se os arquivos foram encontrados
if [ -z "$TOTAL_COMPACTED_BYTES" ] || [ "$TOTAL_COMPACTED_BYTES" == "0" ]; then
    echo "ERROR: No '.fastq.gz' files found in the directory '$INPUT_DIR'."
    exit 1
fi

# 1.2 Estima o tamanho necessário (descompactado)
EXPANSION_RATE=4
NEEDED_BYTES=$(echo "$TOTAL_COMPACTED_BYTES * $EXPANSION_RATE" | bc)

# 1.3 Obtém o espaço livre em bytes no sistema de arquivos atual
FREE_BYTES=$(df -B 1 . | awk 'NR==2 {print $4}')

# --- 2. Conversão para Formato Legível (GB) ---

# Função para converter bytes para Gigabytes (com duas casas decimais)
bytes_to_gb() {
    # Divide por (1024^3) para obter GB
    echo "scale=2; $1 / 1073741824" | bc
}

TOTAL_COMPACTED_GB=$(bytes_to_gb "$TOTAL_COMPACTED_BYTES")
NEEDED_GB=$(bytes_to_gb "$NEEDED_BYTES")
FREE_GB=$(bytes_to_gb "$FREE_BYTES")

# --- 3. Imprime Resultados ---

echo "------------------------------------------------------"
# LINHAS CORRIGIDAS: Removido alinhamento por espaços invisíveis
echo "Total size of compressed files (.fastq.gz): ${TOTAL_COMPACTED_GB} GB"
echo "Required disk space (4x estimate): ${NEEDED_GB} GB"
echo "Available free space on disk: ${FREE_GB} GB"
echo "------------------------------------------------------"

# --- 4. Verificação de Suficiência de Espaço ---

if [ "$FREE_BYTES" -gt "$NEEDED_BYTES" ]; then
    
    # Garante uma margem de segurança de 10%
    MARGIN_CHECK=$(echo "scale=0; $NEEDED_BYTES * 1.1 / 1" | bc) 
    
    # LINHA CORRIGIDA: Removido o espaço invisível
    if [ "$FREE_BYTES" -gt "$MARGIN_CHECK" ]; then
        echo "SPACE SUFFICIENT. Safety margin OK."
    else
        echo "ATTENTION: Space is sufficient, but the margin is tight."
    fi
else
    echo "CRITICAL ERROR: Free space (${FREE_GB} GB) is insufficient for decompression (${NEEDED_GB} GB)."
    echo "Free up space and try again. Aborting."
    exit 1
fi

# --- 5. User Confirmation ---

echo ""
read -p "Press [ENTER] to start SMART DECOMPRESSION, or Ctrl+C to cancel."

# --- 6. Decompression Execution (MODIFICADO para ignorar arquivos existentes) ---

echo ""
echo "Starting smart decompression in the '$INPUT_DIR' folder..."

# Inicializa o sinalizador de erro
DECOMPRESSION_FAILED=0

# Loop através de todos os arquivos compactados
for GZ_FILE in "$INPUT_DIR"/*.fastq.gz; do
    # Verifica se encontrou arquivos .gz
    if [ ! -f "$GZ_FILE" ]; then
        continue
    fi

    # Determina o nome do arquivo descompactado (removendo .gz)
    FASTQ_FILE="${GZ_FILE%.gz}"
    
    # Verifica se o arquivo descompactado JÁ EXISTE
    if [ -f "$FASTQ_FILE" ]; then
        echo "Skipping $(basename "$GZ_FILE"): '$FASTQ_FILE' already exists."
        continue
    fi

    echo "Decompressing: $(basename "$GZ_FILE")"
    
    # Executa o gunzip -k para descompactar, MANTENDO o .gz original
    gunzip -k "$GZ_FILE"
    
    # Verifica se houve erro no gunzip
    if [ $? -ne 0 ]; then
        echo "ERROR: Decompression failed for $GZ_FILE."
        DECOMPRESSION_FAILED=1
    fi
done

if [ "$DECOMPRESSION_FAILED" -ne 0 ]; then
    echo "--- Decompression Completed with Errors ---"
else
    echo "--- Decompression Complete ---"
fi

# --- 7. Run MultiQC Report Generation ---

echo "Generating MultiQC report..."

# Verifica se houve falha de descompactação antes de continuar
if [ "$DECOMPRESSION_FAILED" -ne 0 ]; then
    echo "Skipping MultiQC report due to prior decompression errors."
    exit 1
fi

# Verifica se o script MultiQC existe e é executável
if [ -x "$REPORT_SCRIPT" ]; then
    # Executa o script de relatório
    "$REPORT_SCRIPT" "$INPUT_DIR"
    echo "MultiQC report generation requested. Check 'reports/' folder for results."
else
    echo "WARNING: MultiQC script ($REPORT_SCRIPT) not found or not executable. Skipping report generation."
    echo "Please ensure the next script is created and set to executable (chmod +x $REPORT_SCRIPT)."
fi

echo "--- Analysis Workflow Complete ---"