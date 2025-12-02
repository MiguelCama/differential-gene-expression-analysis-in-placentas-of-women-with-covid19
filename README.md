# Differential Gene Expression Analysis in Placentas from COVID-19–Positive Pregnancies

This repository contains all scripts used to perform a full RNA-seq pipeline for differential gene expression analysis in placental samples from pregnant women infected with SARS-CoV-2 compared with healthy controls.

COVID-19 is not only a respiratory condition but also a vascular disease capable of causing significant inflammation and circulatory impairment. During pregnancy, SARS-CoV-2 infection may alter placental physiology and influence fetal development. In this study, we analyzed placental transcriptomes and identified distinct gene-expression patterns between infected and control groups, suggesting potential placental dysfunction associated with COVID-19. These results highlight the importance of transcriptomics for understanding pregnancy-related alterations linked to viral infection.

---

## 📁 Repository Structure

This project includes scripts for **quality control**, **preprocessing**, **alignment**, **quantification**, and **downstream differential expression and functional analysis**.

### **1. Preprocessing & Quality Control**
- `00a_downloadAndDecompress.sh` – Download and decompress raw data  
- `decompress.sh` – Additional decompression utilities  
- `getSRR.sh` – Retrieve SRR data  
- `01a_run_trimmomatic.sh` – Read trimming  
- `02a_run_cutadapt.sh` – Adapter trimming  
- `00b_run_multiqc_pre.sh` – QC before processing  
- `01b_run_multiqc_post.sh` / `02b_run_multiqc_post.sh` – QC after processing  

### **2. Alignment (STAR)**
- `00_create_star_index.sh` – STAR genome index creation  
- `01_run_star_parallel.sh` – Parallelized alignment  
- `01_renameMappedSamples_srr.sh` – Sample renaming  
- `02_run_samtools_parallel.sh` – BAM processing  
- `03_run_featurecounts_quantificationStar.sh` – Gene quantification via FeatureCounts  

### **3. Kallisto-based Quantification**
- `build_kallisto_index.sh` / `kallisto_index.sh` – Kallisto index creation  
- `kallistoQuant.sh` – Quantification  
- `consolidate_kallisto_results.sh` – Merge Kallisto outputs  

### **4. Annotation Processing**
- `convert_gff_to_bed.sh` – Convert GFF to BED  
- `split_gff_by_feature.py` – Split GFF by feature type  
- `process_bed_attributes.R` – Parse BED attributes  

### **5. Cleanup & Utilities**
- `03_cleanUp.sh` – Remove intermediate files  
- `filtering.sh` – Filtering utilities  
- `fixSampleNames.sh` – Standardize sample naming  

### **6. Downstream Analysis**
- `DESeq2.R` – Differential gene expression analysis  
- `GO.R` – Gene Ontology enrichment  

---

## 📄 Dataset File: *dados-projeto.xlsx*

The repository includes an Excel file containing metadata for all samples used in the analysis:

👉 **Download:** `dados-projeto.xlsx`

| SRR ID        | Sample | Condition                     |
|---------------|--------|-------------------------------|
| SRR15312890   | H_1    | INFECTED (High Viral Load)    |
| SRR15312891   | H_2    | INFECTED (High Viral Load)    |
| SRR15312892   | P_1    | INFECTED                      |
| SRR15312893   | P_2    | INFECTED                      |
| SRR15312894   | P_3    | INFECTED                      |
| SRR15312895   | P_4    | INFECTED                      |
| SRR15312896   | P_5    | INFECTED                      |
| SRR15312897   | P_6    | INFECTED                      |
| SRR15312898   | P_7    | INFECTED                      |
| SRR15312899   | P_8    | INFECTED                      |
| SRR15312900   | P_9    | INFECTED                      |
| SRR15312901   | P_10   | INFECTED                      |
| SRR15312902   | P_11   | INFECTED                      |
| SRR15312903   | P_12   | INFECTED                      |
| SRR15312904   | P_13   | INFECTED                      |
| SRR15312905   | P_14   | INFECTED                      |
| SRR15312906   | P_15   | INFECTED                      |
| SRR15312907   | C_1    | CONTROL (Non-Infected)        |
| SRR15312908   | C_2    | CONTROL (Non-Infected)        |
| SRR15312910   | C_3    | CONTROL (Non-Infected)        |
| SRR15312911   | C_4    | CONTROL (Non-Infected)        |

---

## 🚀 Requirements

Depending on pipeline components:

- **Linux environment**
- **R** (DESeq2, tidyverse, clusterProfiler, etc.)
- **Python 3** (pandas, argparse)
- **STAR**
- **samtools**
- **Kallisto**
- **MultiQC**
- Standard bioinformatics tools available in `$PATH`

---

## 🎯 Purpose

This repository ensures full reproducibility of our transcriptomic analysis and supports research on how SARS-CoV-2 infection may impact placental biology and pregnancy outcomes.

---


