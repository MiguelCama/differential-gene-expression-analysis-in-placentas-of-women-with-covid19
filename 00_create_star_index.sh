#!/bin/bash

# --- 00_create_star_index_human.sh: STAR Index Generation (GRCh38 ONLY) ---

# --- Configuration ---
# Directory where the index will be saved.
STAR_INDEX_DIR="genome/index/star_human"

# Arquivos de entrada corrigidos para o caminho exato que você forneceu:
# ATENÇÃO: Os caminhos são relativos ao diretório onde o script será executado.
GENOME_FASTA="genome/GCA_000001405.15_GRCh38_full_analysis_set.fna" 
ANNOTATION_GFF="annotation/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sorted.gff"

# Number of threads for index generation (Can be high, e.g., 20-30)
THREADS=45 

mkdir -p "$STAR_INDEX_DIR"

echo "--- Starting STAR Index Generation (Human ONLY) ---"

# Verifica a existência dos arquivos de entrada antes de prosseguir
if [ ! -f "$GENOME_FASTA" ]; then
    echo "!!! ERRO: Arquivo FASTA não encontrado em: $GENOME_FASTA"
    exit 1
fi
if [ ! -f "$ANNOTATION_GFF" ]; then
    echo "!!! ERRO: Arquivo GFF não encontrado em: $ANNOTATION_GFF"
    exit 1
fi

# --- STAR Genome Generate Command ---
saveCommand STAR --runThreadN $THREADS \
     --runMode genomeGenerate \
     --genomeDir "$STAR_INDEX_DIR" \
     --genomeFastaFiles "$GENOME_FASTA" \
     --sjdbGTFfile "$ANNOTATION_GFF" \
     --sjdbOverhang 100 

if [ $? -ne 0 ]; then
    echo "!!! ERRO FATAL na geração do índice STAR. Verifique os logs e arquivos."
    exit 1
fi

echo "--- STAR Index Generation COMPLETED successfully ---"