#!/bin/bash

# --- 02_run_featurecounts_quantification.sh: Count Gene features using featureCounts ---

# --- Configuration ---
# Directory where STAR placed the final sorted BAM files
BAM_DIR="star/alignment"  # MUDEI: agora aponta para os resultados do STAR
# Annotation file: Using the full GTF path you provided
ANNOTATION_FILE="annotation/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf" 
# Output directory for the final count matrix
OUTPUT_DIR="featurecounts/star_results" 
# Number of threads for internal parallelization
THREADS=55 

# --- Files and Paths ---
OUTPUT_MATRIX="$OUTPUT_DIR/raw_counts_GENE_matrix.txt"
LOG_FILE="$OUTPUT_DIR/featurecounts_log.txt"

mkdir -p "$OUTPUT_DIR"

echo "--- Starting featureCounts Quantification (GTF, Exon, Stranded -s 2) ---"
echo "Processing all samples simultaneously using $THREADS threads."
echo "Using annotation file: $ANNOTATION_FILE"
echo "BAM files directory: $BAM_DIR"
echo "---"

# Find all STAR sorted BAM files (*Aligned.sortedByCoord.out.bam is the STAR output naming)
BAM_FILES=$(find "$BAM_DIR" -maxdepth 1 -name "*Aligned.sortedByCoord.out.bam" | sort | tr '\n' ' ')

if [ -z "$BAM_FILES" ]; then
    echo "!!! ERROR: No BAM files found in $BAM_DIR. Aborting quantification. !!!"
    exit 1
fi

echo "Found $(echo $BAM_FILES | wc -w) BAM files to process:"
for bam in $BAM_FILES; do
    echo "  - $(basename $bam)"
done

# --- featureCounts Command (Optimized for STAR BAMs) ---
# -F GTF: Specifies the format.
# -t 'exon': Specifies the feature type to count.
# -g 'gene_id': Specifies the attribute to group features by gene.
# -s 2: REVERSELY STRANDED (TruSeq Stranded kit) - IMPORTANTE: verifique se é o correto para seu protocolo
# -p: count fragments (para paired-end)
# --countReadPairs: count read pairs instead of reads
# -B: only count read pairs that have both ends aligned
# -C: do not count read pairs that have ends mapping to different chromosomes or same chr but different strands

echo "Starting featureCounts at $(date)"
saveCommand featureCounts -a "$ANNOTATION_FILE" \
             -o "$OUTPUT_MATRIX" \
             -F GTF -t 'exon' -g 'gene_id' \
             -p --countReadPairs \
             -B -C \
             -T "$THREADS" \
             --primary \
             --ignoreDup \
             -s 2 \
             $BAM_FILES 2>> "$LOG_FILE" 

FC_EXIT_CODE=$?

if [ $FC_EXIT_CODE -ne 0 ]; then
    echo "!!! FATAL ERROR in featureCounts. Check $LOG_FILE for details."
    # If the error is about the attribute, suggest alternatives
    if grep -q "failed to find the gene identifier attribute" "$LOG_FILE"; then
        echo "HINT: The attribute 'gene_id' might be incorrect. Try changing '-g 'gene_id'' to '-g 'gene'' or '-g 'Parent''."
    fi
    exit 1
else
    echo "Quantification COMPLETED successfully at $(date)."
    echo "Count matrix saved to: $OUTPUT_MATRIX"
    
    # Print summary of results
    if [ -f "$OUTPUT_MATRIX" ]; then
        echo "--- Summary of counts ---"
        echo "Total genes counted: $(tail -n +3 "$OUTPUT_MATRIX" | wc -l)"
        echo "Matrix header:"
        head -n 1 "$OUTPUT_MATRIX" | cut -f1-5
    fi
fi

echo "--- featureCounts process finished. ---"