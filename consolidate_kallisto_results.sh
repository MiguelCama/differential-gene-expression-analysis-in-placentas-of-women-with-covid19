#!/bin/bash

# --- 03_consolidate_kallisto_results.sh: Combine all Kallisto results and add gene IDs ---

# --- Configuration ---
KALLISTO_DIR="kallisto/quant_results"
OUTPUT_DIR="kallisto/consolidated"
DESCRIPTION_TSV="kallisto/output/description.tsv" 
SAMPLE_FILE="config/samples.txt"

mkdir -p "$OUTPUT_DIR"

echo "--- Starting Kallisto Results Consolidation ---"

# --- Step 1: Create the corrected Python script ---
cat > "$OUTPUT_DIR/combine_kallisto.py" << 'EOF'
# -*- coding: utf-8 -*-
import os
import sys
import glob
from collections import defaultdict

# Configuration
kallisto_dir = "kallisto/quant_results"
output_dir = "kallisto/consolidated"
description_file = "kallisto/output/description.tsv" 
sample_file = "config/samples.txt"

def read_sample_metadata(sample_file):
    """Read sample metadata from samples.txt file"""
    sample_mapping = {}
    try:
        with open(sample_file, 'r') as f:
            header = f.readline().strip().split('\t')
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    srr = parts[0]
                    sample_name = parts[1]
                    sample_mapping[srr] = sample_name
        print(f"Successfully loaded metadata for {len(sample_mapping)} samples")
        return sample_mapping
    except Exception as e:
        print(f"Error reading sample file: {e}")
        return {}

def read_description_tsv(description_file):
    """Read the pre-processed description TSV for annotations"""
    print(f"Reading pre-processed descriptions from {description_file}...")
    annotations = {}
    try:
        with open(description_file, 'r') as f:
            header = f.readline().strip().split('\t')
            if header != ['target_id', 'gene_id', 'product_description']:
                print(f"Warning: Description file header mismatch: {header}")
                
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 3:
                    target_id, gene_id, product_description = parts
                    annotations[target_id] = (gene_id, product_description)
            
        print(f"Successfully loaded {len(annotations)} annotations.")
        return annotations
    except Exception as e:
        print(f"Error reading description file {description_file}: {e}")
        return {}


def read_kallisto_abundance(abundance_file):
    """Read Kallisto abundance.tsv file without pandas"""
    data = {}
    try:
        with open(abundance_file, 'r') as f:
            header = f.readline().strip().split('\t')
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    target_id = parts[0]
                    est_counts = float(parts[3])
                    tpm = float(parts[4])
                    data[target_id] = (est_counts, tpm)
        return data
    except Exception as e:
        print(f"Error reading abundance file: {e}")
        return {}

def combine_kallisto_results():
    """Combine all Kallisto abundance files without pandas"""
    
    # Read sample metadata
    print("Reading sample metadata...")
    sample_mapping = read_sample_metadata(sample_file)
    if not sample_mapping:
        print("Error: No sample metadata found!")
        return False
    
    # Read annotations from the pre-processed TSV file
    annotations = read_description_tsv(description_file)
    if not annotations:
        print("Error: No annotations loaded!")
        return False
        
    # Find all Kallisto result directories
    print("Finding Kallisto result directories...")
    sample_dirs = {}
    missing_samples = []
    
    for srr, sample_name in sample_mapping.items():
        dir_pattern1 = os.path.join(kallisto_dir, f"{sample_name}_{srr}")
        dir_pattern2 = os.path.join(kallisto_dir, f"{srr}")
        
        if os.path.exists(dir_pattern1) and os.path.exists(os.path.join(dir_pattern1, "abundance.tsv")):
            sample_dirs[sample_name] = dir_pattern1
        elif os.path.exists(dir_pattern2) and os.path.exists(os.path.join(dir_pattern2, "abundance.tsv")):
            sample_dirs[sample_name] = dir_pattern2
        else:
            missing_samples.append((sample_name, srr))
    
    if missing_samples:
        print(f"Warning: {len(missing_samples)} samples not found: (Showing first 5)")
        for sample_name, srr in missing_samples[:5]:
            print(f"  - {sample_name} (SRR: {srr})")
        
    if not sample_dirs:
        print("Error: No Kallisto results found!")
        return False
    
    # Initialize data structures
    all_target_ids = set()
    est_counts_data = defaultdict(dict)
    tpm_data = defaultdict(dict)
    
    # Process each sample
    print(f"Processing {len(sample_dirs)} samples...")
    for sample_name, dir_path in sample_dirs.items():
        abundance_file = os.path.join(dir_path, "abundance.tsv")
        
        print(f"  Processing {sample_name}...")
        data = read_kallisto_abundance(abundance_file)
        if not data:
            print(f"    Warning: No data found for {sample_name}")
            continue
        
        # Store data
        for target_id, (est_count, tpm) in data.items():
            all_target_ids.add(target_id)
            est_counts_data[target_id][sample_name] = est_count
            tpm_data[target_id][sample_name] = tpm
    
    # Create sorted lists
    sorted_target_ids = sorted(all_target_ids)
    sorted_samples = sorted(sample_dirs.keys())
    
    # DEBUG: Check annotation matching
    print("\n=== DEBUG: Checking annotation matching ===")
    matched_count = 0
    for target_id in sorted_target_ids:
        if target_id in annotations:
            matched_count += 1
    
    print(f"Overall matching (Kallisto IDs in description TSV): {matched_count}/{len(sorted_target_ids)} ({matched_count/len(sorted_target_ids)*100:.2f}%)")
    
    # Check first 5 matches
    print("\nFirst 5 annotated transcripts:")
    found_preview = 0
    for target_id in sorted_target_ids:
        if target_id in annotations and found_preview < 5:
            gene_id, product = annotations[target_id]
            print(f"  {target_id} -> Gene: {gene_id}, Product: {product[:50]}...")
            found_preview += 1
    print("=== END DEBUG ===\n")
    
    # Write est_counts file
    print("Writing est_counts file with annotations...")
    with open(os.path.join(output_dir, "kallisto_est_counts.tsv"), 'w') as f:
        # Header: target_id, gene_id, product_description, sample counts...
        f.write("target_id\tgene_id\tproduct_description\t" + "\t".join(sorted_samples) + "\n")
        
        # Write data
        for target_id in sorted_target_ids:
            gene_id, product = annotations.get(target_id, ("unknown", "unknown"))
            counts = [str(est_counts_data[target_id].get(sample, 0)) for sample in sorted_samples]
            f.write(f"{target_id}\t{gene_id}\t{product}\t" + "\t".join(counts) + "\n")
    
    # Write TPM file
    print("Writing TPM file with annotations...")
    with open(os.path.join(output_dir, "kallisto_tpm.tsv"), 'w') as f:
        # Header: target_id, gene_id, product_description, sample TPMs...
        f.write("target_id\tgene_id\tproduct_description\t" + "\t".join(sorted_samples) + "\n")
        
        # Write data
        for target_id in sorted_target_ids:
            gene_id, product = annotations.get(target_id, ("unknown", "unknown"))
            tpms = [str(tpm_data[target_id].get(sample, 0)) for sample in sorted_samples]
            f.write(f"{target_id}\t{gene_id}\t{product}\t" + "\t".join(tpms) + "\n")
    
    # Create detailed summary
    total_transcripts = len(sorted_target_ids)
    annotated_transcripts = sum(1 for target_id in sorted_target_ids if target_id in annotations)
    
    summary_stats = {
        'total_transcripts': total_transcripts,
        'samples_processed': len(sample_dirs),
        'annotated_transcripts': annotated_transcripts,
        'samples_found': sorted_samples
    }
    
    with open(os.path.join(output_dir, "consolidation_summary.txt"), 'w') as f:
        f.write("Kallisto Results Consolidation Summary\n")
        f.write("=====================================\n\n")
        f.write("Annotation Statistics:\n")
        f.write(f"  Total transcripts: {summary_stats['total_transcripts']}\n")
        f.write(f"  Transcripts Annotated (matched in description TSV): {summary_stats['annotated_transcripts']} ({summary_stats['annotated_transcripts']/summary_stats['total_transcripts']*100:.1f}%)\n")
        f.write(f"  Samples processed: {summary_stats['samples_processed']}\n\n")
        f.write("Samples found:\n")
        for sample in summary_stats['samples_found']:
            f.write(f"  - {sample}\n")
            
    print("\nConsolidation completed successfully!")
    print(f"Total transcripts: {summary_stats['total_transcripts']}")
    print(f"Samples processed: {summary_stats['samples_processed']}")
    print(f"Annotated transcripts: {summary_stats['annotated_transcripts']} ({summary_stats['annotated_transcripts']/summary_stats['total_transcripts']*100:.1f}%)")
    
    return True

if __name__ == "__main__":
    success = combine_kallisto_results()
    sys.exit(0 if success else 1)
EOF

# --- Step 2: Run the corrected Python script ---
echo "Running consolidation script..."
python3 "$OUTPUT_DIR/combine_kallisto.py"

if [ $? -eq 0 ]; then
    echo "--- Consolidation COMPLETED ---"
    echo "Output files:"
    echo "  - $OUTPUT_DIR/kallisto_est_counts.tsv (with enriched annotations)"
    echo "  - $OUTPUT_DIR/kallisto_tpm.tsv (with enriched annotations)"  
    echo "  - $OUTPUT_DIR/consolidation_summary.txt"
    
    # Show preview of results
    echo ""
    echo "--- Preview of est_counts file (with annotations) ---"
    head -3 "$OUTPUT_DIR/kallisto_est_counts.tsv" | cut -f1-8
    echo "..."
    
    echo ""
    echo "--- Preview of TPM file (with annotations) ---"
    head -3 "$OUTPUT_DIR/kallisto_tpm.tsv" | cut -f1-8
    echo "..."
    
else
    echo "!!! ERROR in consolidation script !!!"
    exit 1
fi

echo "--- Kallisto results consolidation finished ---"