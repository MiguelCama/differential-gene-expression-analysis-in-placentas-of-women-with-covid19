import sys
import os

def split_gff_by_feature(input_gff_path):
    """
    Reads a GFF file and splits it into multiple files based on the feature type (column 3).
    Header lines (starting with #) are preserved and included in every output file.
    """
    # --- Configuration ---
    OUTPUT_DIR = "gff_split"
    TAB = '\t'
    
    # Check if input file exists
    if not os.path.exists(input_gff_path):
        print(f"ERROR: Input file '{input_gff_path}' not found.")
        sys.exit(1)

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output files will be saved in: {OUTPUT_DIR}/")
    print("---")

    headers = []
    feature_types = set()

    print("--- 1. Identifying unique feature types and extracting headers ---")

    try:
        with open(input_gff_path, 'r') as infile:
            for line in infile:
                # Remove trailing whitespace
                line = line.rstrip('\n')
                
                if line.startswith('##') or line.startswith('#!'):
                    # Collect header lines (lines starting with #, including comments and directives)
                    headers.append(line)
                elif not line.startswith('#'):
                    # Process feature lines
                    try:
                        # Split line by tab. GFF is strictly tab-separated.
                        fields = line.split(TAB)
                        # Column 3 (index 2) is the feature type
                        if len(fields) >= 3:
                            feature_type = fields[2]
                            feature_types.add(feature_type)
                    except IndexError:
                        # Skip lines that might be malformed but don't start with '#'
                        continue

    except Exception as e:
        print(f"CRITICAL ERROR while reading input file: {e}")
        sys.exit(1)

    # Check if any feature types were found
    if not feature_types:
        print("WARNING: No features found in the GFF file after skipping headers. Exiting.")
        sys.exit(0)

    print(f"Found the following features: {', '.join(sorted(feature_types))}")
    print("---")
    
    # Pre-open output files dictionary for writing
    output_handles = {}

    try:
        # --- 2. Write headers and open file handles ---
        for feature in feature_types:
            output_file = os.path.join(OUTPUT_DIR, f"{feature}.gff")
            output_handles[feature] = open(output_file, 'w')
            
            # Write all headers to the new file
            for header_line in headers:
                output_handles[feature].write(header_line + '\n')
            
            print(f"Prepared file for feature: '{feature}' -> {output_file}")

        # --- 3. Process the file again and distribute feature lines ---
        print("\n--- Distributing feature lines to respective files ---")
        
        with open(input_gff_path, 'r') as infile:
            for line in infile:
                line = line.rstrip('\n')
                
                # We only care about non-header lines now
                if not line.startswith('#'):
                    try:
                        fields = line.split(TAB)
                        if len(fields) >= 3:
                            feature_type = fields[2]
                            
                            # Write line to the corresponding file handle
                            if feature_type in output_handles:
                                output_handles[feature_type].write(line + '\n')
                            
                    except IndexError:
                        # Ignore malformed lines
                        continue
    
    except Exception as e:
        print(f"ERROR during file distribution: {e}")
    
    finally:
        # --- 4. Clean up and close all file handles ---
        for handle in output_handles.values():
            handle.close()

    print("---")
    print(f"Processing complete. Check the '{OUTPUT_DIR}' folder for results.")

if __name__ == "__main__":
    # Check for command line argument
    if len(sys.argv) != 2:
        print("ERROR: Please provide the input GFF file path.")
        print(f"Usage: python {sys.argv[0]} <input_gff_file>")
        sys.exit(1)
        
    gff_file_path = sys.argv[1]
    split_gff_by_feature(gff_file_path)