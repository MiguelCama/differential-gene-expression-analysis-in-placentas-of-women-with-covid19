#!/bin/bash

# --- 01_run_star_parallel_human.sh: Parallel STAR Spliced Mapping (Human ONLY) ---

# --- Configuration ---
READS_DIR="output/final_trimmed_reads"
STAR_INDEX_DIR="genome/index/star_human" # Novo índice humano
OUTPUT_DIR="star/alignment"
LOGS_DIR="star/logs"
THREADS_PER_JOB=8 

# Create output directories
mkdir -p "$OUTPUT_DIR" "$LOGS_DIR"

echo "--- Starting Parallel STAR Mapping (Human ONLY) ---"
# Utiliza a variável SLURM_NTASKS ou um valor padrão para paralelização
SLURM_NTASKS=${SLURM_NTASKS:-75} 
MAX_JOBS=$(($SLURM_NTASKS / $THREADS_PER_JOB))
echo "Threads per Job: $THREADS_PER_JOB. Max Concurrent Jobs: $MAX_JOBS"

# --- Function for STAR ---
run_star() {
    local SAMPLE_NAME="$1"
    
    local R1="${READS_DIR}/${SAMPLE_NAME}_1_paired.fastq.gz"
    local R2="${READS_DIR}/${SAMPLE_NAME}_2_paired.fastq.gz"
    local STAR_OUTPUT_PREFIX="${OUTPUT_DIR}/${SAMPLE_NAME}."
    local LOG_FILE="${LOGS_DIR}/${SAMPLE_NAME}.star.log"
    local THREADS="${THREADS_PER_JOB}"

    echo "[$(date +%H:%M:%S)] Starting STAR: $SAMPLE_NAME"
    echo "STAR Command Started at $(date)" > "$LOG_FILE"
    
    # STAR command (Outputting BAM, SortedByCoordinate)
    saveCommand STAR --genomeDir "$STAR_INDEX_DIR" \
         --readFilesIn "$R1" "$R2" \
         --runThreadN "$THREADS" \
         --outFileNamePrefix "$STAR_OUTPUT_PREFIX" \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterMultimapNmax 20 \
         --outFilterScoreMinOverLread 0.5 \
         --outFilterMatchNminOverLread 0.5 \
         --outFilterType BySJout \
         --twopassMode Basic \
         2>> "$LOG_FILE" 

    if [ $? -ne 0 ]; then
        echo "!!! ERROR in STAR for $SAMPLE_NAME. Check $LOG_FILE" >> "$LOG_FILE"
        return 1
    fi
    
    echo "[$(date +%H:%M:%S)] STAR COMPLETED successfully for $SAMPLE_NAME."
    return 0
}

# --- Parallel Execution ---
export -f run_star
export READS_DIR STAR_INDEX_DIR OUTPUT_DIR LOGS_DIR THREADS_PER_JOB

# Generate the list of sample base names (Ex: C_1_SRR15312907, P_1_SRR15312892, etc.)
SAMPLE_LIST=$(find "$READS_DIR" -maxdepth 1 -name "*_1_paired.fastq.gz" -printf "%f\n" | sed 's/_1_paired\.fastq\.gz//')

# Run STAR in parallel
echo "$SAMPLE_LIST" | parallel \
    --jobs "$MAX_JOBS" \
    --joblog "$LOGS_DIR/star_job_status.log" \
    run_star {}

echo "--- End of Parallel STAR Mapping ---"