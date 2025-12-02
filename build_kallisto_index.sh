#!/bin/bash

# ==============================================================================
# Kallisto Index Preparation Script (CORRIGIDO)
# Focando exclusivamente na feature 'mRNA' para quantificação.
# ==============================================================================

# --- File and Path Configuration (Read-Only) ---
GENOME_FASTA="kallisto/input/genoma.fna"
ANNOTATION_GFF="kallisto/input/anotacao_genes.gff"

# Features de interesse (Hardcoded conforme solicitado)
FEATURE_TYPE="mRNA"
FILE_SUFFIX="mRNA"

# Intermediate Files
NORMALIZED_GFF="kallisto/input/anotacao_genes.normalized.gff"
GFF_FOR_BEDTOOLS="kallisto/input/transcripts_bedtools.gff"
OUTPUT_FASTA="kallisto/index/transcripts_raw.fasta"

# Dynamic Output Files
CLEANED_FASTA="kallisto/index/transcripts_${FILE_SUFFIX}_cleaned.fasta"
KALLISTO_INDEX="kallisto/index/kallisto_${FILE_SUFFIX}_index.idx"

# --- Setup Directories ---
mkdir -p kallisto/input kallisto/index

echo "Starting Kallisto index preparation (Target: $FEATURE_TYPE)..."
echo "Output index will be: $KALLISTO_INDEX"
echo "--------------------------------------------------------"

# --- 1. GFF Normalization (Ensuring Tab Delimiters) ---
echo "1/4: Normalizando GFF (convertendo espaços para TABs)..."

# CRITICAL FIX: Use sed in three steps for maximum robustness against mixed whitespace:
# 1. Remove leading/trailing whitespace.
# 2. Convert all sequences of whitespace (spaces, tabs, etc.) into a single space.
# 3. Convert all remaining single spaces into a single TAB delimiter.
cat "$ANNOTATION_GFF" | \
    sed -E 's/^[[:space:]]*//; s/[[:space:]]*$//; s/[[:space:]]+/ /g; s/ /\t/g' > "$NORMALIZED_GFF"

if [ ! -s "$NORMALIZED_GFF" ]; then
    echo "ERROR: GFF normalization failed. The file is empty."
    exit 1
fi
echo "Normalização concluída: $(wc -l < "$NORMALIZED_GFF") linhas."


# --- 2. Filter GFF and Prepare for bedtools (Injecting Name attribute) ---
echo "2/4: Filtrando por '$FEATURE_TYPE' e injetando 'Name=' com 'transcript_id'..."

# AWK: 
# 1. Filtra na coluna 3 (FEATURE_TYPE).
# 2. Extrai o transcript_id (o mais confiável) ou o ID.
# 3. Garante que o atributo Name= seja adicionado (necessário para bedtools -name).
awk -v FEATURE="$FEATURE_TYPE" '
BEGIN {OFS = "\t"} 
/^#/ {print; next} # Print headers

{
    # 1. Filtra a feature na coluna 3
    if ($3 == FEATURE) {
        transcript_id = ""
        # 2. Tenta extrair o transcript_id
        if (match($9, /transcript_id=([^;]+)/, arr)) {
            transcript_id = arr[1]
        }
        # Fallback: Tenta extrair o ID (pode ser usado se o transcript_id estiver faltando)
        else if (match($9, /ID=([^;]+)/, arr)) {
            transcript_id = arr[1]
        }

        # 3. Injete o atributo Name
        if (transcript_id != "") {
            # Apenas adiciona Name= se não existir, para bedtools getfasta
            if ($9 !~ /Name=/) {
                # Usa o transcript_id para o Name
                $9 = $9 ";Name=" transcript_id
            }
            print $0
        } else {
            print "WARNING: Linha de " FEATURE " pulada por falta de ID válido: " $0 > "/dev/stderr"
        }
    }
}
' "$NORMALIZED_GFF" > "$GFF_FOR_BEDTOOLS"


if [ ! -s "$GFF_FOR_BEDTOOLS" ]; then
    echo "ERROR: O arquivo GFF filtrado ($GFF_FOR_BEDTOOLS) está vazio. Cheque o arquivo de entrada ou a extração de ID."
    exit 1
fi
echo "Total de linhas $FEATURE_TYPE prontas: $(wc -l < "$GFF_FOR_BEDTOOLS")"


# --- 3. Generate Transcript FASTA (using bedtools getfasta) ---
echo "3/4: Extraindo sequências FASTA (getfasta)..."
# bedtools getfasta usa o atributo Name na coluna 9 como cabeçalho FASTA.
saveCommand bedtools getfasta -fi "$GENOME_FASTA" -bed "$GFF_FOR_BEDTOOLS" -fo "$OUTPUT_FASTA" -name

if [ ! -s "$OUTPUT_FASTA" ]; then
    echo "ERROR: A extração FASTA falhou. Verifique o caminho do genoma e as coordenadas GFF."
    exit 1
fi
echo "Arquivo FASTA bruto gerado: $OUTPUT_FASTA"


# --- 4. Clean FASTA Headers and Index ---
echo "4/4: Limpando cabeçalhos FASTA e criando o índice Kallisto..."

# Remove a sujeira PÓS-ID gerada por bedtools (::chr:start-end)
sed '/>/ s/::.*//' "$OUTPUT_FASTA" > "$CLEANED_FASTA"

if [ ! -s "$CLEANED_FASTA" ]; then
    echo "ERROR: A limpeza do FASTA falhou."
    exit 1
fi

# Criação do índice Kallisto (Usando o FASTA limpo)
saveCommand kallisto index -i "$KALLISTO_INDEX" "$CLEANED_FASTA" --make-unique

if [ $? -ne 0 ]; then
    echo "Kallisto indexing failed. Review the output above."
    exit 1
fi

echo "--------------------------------------------------------"
echo "Kallisto indexing process completed successfully!"
echo "Index File: $KALLISTO_INDEX"
echo "--------------------------------------------------------"