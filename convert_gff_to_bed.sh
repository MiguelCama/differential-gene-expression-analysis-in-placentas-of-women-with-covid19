#!/bin/bash

# ==============================================================================
# SCRIPT: convert_gff_to_bed.sh
# DESCRIPTION: Converts all GFF files in the 'gff_split/' directory 
#              into BED format using the 'gff2bed' tool.
# USAGE: ./convert_gff_to_bed.sh
# ==============================================================================

# --- Configuration ---
INPUT_DIR="gff_split"
OUTPUT_DIR="bedFiles"
GFF2BED_TOOL="gff2bed" # Assumes 'gff2bed' is in your PATH

echo "--- Starting GFF to BED Conversion ---"
echo "Input Directory: $INPUT_DIR"
echo "Output Directory: $OUTPUT_DIR"

# Check if the gff_split directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory '$INPUT_DIR' not found. Please run the GFF splitting script first."
    exit 1
fi

# Check if the gff2bed tool is available
if ! command -v "$GFF2BED_TOOL" &> /dev/null
then
    echo "CRITICAL ERROR: The required tool '$GFF2BED_TOOL' could not be found."
    echo "Please ensure 'gff2bed' (from the BEDOPS suite) is installed and available in your system's PATH."
    exit 1
fi

# 1. Create the output directory
mkdir -p "$OUTPUT_DIR"
echo "Created output directory: $OUTPUT_DIR"
echo "---"

TOTAL_CONVERTED=0
TOTAL_SKIPPED=0

# 2. Iterate through all GFF files in the input directory
for GFF_FILE in "$INPUT_DIR"/*.gff; do
    
    # Check if the pattern matched any files (prevents error if dir is empty)
    if [ ! -f "$GFF_FILE" ]; then
        continue
    fi
    
    # Extract the base file name (e.g., gene.gff)
    BASE_NAME=$(basename "$GFF_FILE")
    
    # Define the output BED file name (e.g., bedFiles/gene.bed)
    OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME%.gff}.bed"

    echo "Converting $BASE_NAME -> $(basename "$OUTPUT_FILE")"
    
    # Check if the BED file already exists
    if [ -f "$OUTPUT_FILE" ]; then
        echo "  -> Skipping: Output file already exists. Delete it to rerun."
        TOTAL_SKIPPED=$((TOTAL_SKIPPED + 1))
        continue
    fi
    
    # Execute the conversion
    # gff2bed reads from stdin (<) and writes to stdout (>), which we redirect to the output file
    "$GFF2BED_TOOL" < "$GFF_FILE" > "$OUTPUT_FILE"
    
    # Check the exit status of the conversion command
    if [ $? -eq 0 ]; then
        echo "  -> Success."
        TOTAL_CONVERTED=$((TOTAL_CONVERTED + 1))
    else
        echo "  -> ERROR: Conversion failed for $BASE_NAME. Check input file format."
    fi
    
done

echo "---"
echo "Conversion Summary:"
echo "  Total files converted successfully: $TOTAL_CONVERTED"
echo "  Total files skipped (already existed): $TOTAL_SKIPPED"
echo "All BED files are located in the '$OUTPUT_DIR' directory."
echo "--- Conversion Finished ---"