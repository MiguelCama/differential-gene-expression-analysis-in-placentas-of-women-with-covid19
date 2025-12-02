# ==============================================================================
# R SCRIPT FOR PROCESSING BED FILE (mRNA)
# Objective: Extract Gene IDs (Symbol), Product Descriptions, and GenBank IDs
# to create an annotation file 'description.tsv'.
# CORRECTION: Uses the Gene Symbol (from the 'gene=' field in BED) as 'gene_id'
# for compatibility with DESeq2 results.
# ==============================================================================

library(tidyverse)
library(readr)
library(stringr)

# --- 2. Path and Variable Configuration ---
INPUT_FILE <- "bedFiles/mRNA.bed"
OUTPUT_FILE <- "kallisto/output/description.tsv"

# Column 10 (X10) contains all BED attributes
ATTRIB_COL_INDEX <- 10 

# --- 3. Data Reading ---
# Reads the BED file (no header, all columns as text)
data_raw <- read_tsv(
  INPUT_FILE, 
  col_names = FALSE, 
  col_types = cols(.default = "c"), 
  comment = "#" 
)

# --- 4. Processing, Attribute Extraction, and Selection ---
descriptions_df <- data_raw %>%
  # 4.1. Create base columns
  mutate(target_id = paste0(X1, ":", X2, "-", X3),
         ATTRIBUTES = !!sym(paste0("X", ATTRIB_COL_INDEX))) %>%
  
  # 4.2. Field Extraction and Cleaning
  mutate(
    # Extract the Transcript ID (for reference)
    TRANSCRIPT_ID = str_extract(ATTRIBUTES, "transcript_id=[^;]+") %>%
      str_replace("transcript_id=", "") %>% 
      str_remove_all("\"|;"),
      
    # Extract the Gene Symbol (GENE) - THIS IS THE FIELD WE NEED!
    GENE_SYMBOL = str_extract(ATTRIBUTES, "gene=[^;]+") %>%
      str_replace("gene=", "") %>% 
      str_remove_all("\"|;"),
      
    # Extract the Product Description (PRODUCT)
    PRODUCT = str_extract(ATTRIBUTES, "product=[^;]+") %>%
      str_replace("product=", "") %>% 
      str_remove_all("\"|;"),
      
    # NEW: Extract the GenBank ID, looking for the explicit 'GenBank:' pattern
    GENBANK_ID = str_extract(ATTRIBUTES, "GenBank:[^,;]+") %>%
      str_replace("GenBank:", "") %>% 
      str_remove_all("\"|;")
  ) %>%
  
  # 4.3. Final Definition and Selection
  select(
    target_id,
    # <--- CRUCIAL CORRECTION: Use Gene Symbol as 'gene_id' --->
    gene_id = GENE_SYMBOL, 
    product_description = PRODUCT,
    genbank_id = GENBANK_ID # NEW: Add the extracted GenBank ID
  ) %>%
  # Remove rows where the gene symbol was not found (if any)
  filter(!is.na(gene_id) & gene_id != "")

# --- 5. Folder Creation and File Writing ---
dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)
write_tsv(descriptions_df, OUTPUT_FILE)

cat("Processing successful. Description file SAVED.\n")
cat("The 'description.tsv' file now contains 'target_id', 'gene_id', 'product_description', and 'genbank_id'.\n")
cat("Description file path:", OUTPUT_FILE, "\n")