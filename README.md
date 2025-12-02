Differential Gene Expression Analysis in Placentas from COVID-19–Positive Pregnancies

This repository contains all scripts used to perform a full RNA-seq pipeline for differential gene expression analysis in placental samples from pregnant women infected with SARS-CoV-2 compared with healthy controls.
COVID-19 is not only a respiratory condition but also a vascular disease capable of causing significant inflammation and circulatory impairment. During pregnancy, SARS-CoV-2 infection may alter placental physiology and influence fetal development. In this study, we analyzed placental transcriptomes and identified distinct gene-expression patterns between infected and control groups, suggesting potential placental dysfunction associated with COVID-19. These results highlight the importance of transcriptomics for understanding pregnancy-related alterations linked to viral infection.

📁 Repository Structure
This project includes scripts for quality control, preprocessing, alignment, quantification, and downstream differential expression and functional analysis.
1. Preprocessing & Quality Control
  00a_downloadAndDecompress.sh – Download and decompress raw data
  decompress.sh – Additional decompression utilities
  getSRR.sh – Retrieve SRR data
  01a_run_trimmomatic.sh – Read trimming
  02a_run_cutadapt.sh – Adapter trimming
  00b_run_multiqc_pre.sh – QC before processing
  01b_run_multiqc_post.sh / 02b_run_multiqc_post.sh – QC after processing

2. Alignment (STAR)
  00_create_star_index.sh – STAR genome index creation
  01_run_star_parallel.sh – Parallelized alignment
  01_renameMappedSamples_srr.sh – Sample renaming
  02_run_samtools_parallel.sh – BAM processing
  03_run_featurecounts_quantificationStar.sh – Gene quantification via FeatureCounts

3. Kallisto-based Quantification
  build_kallisto_index.sh / kallisto_index.sh – Index creation
  kallistoQuant.sh – Quantification
  consolidate_kallisto_results.sh – Merge Kallisto outputs

4. Annotation Processing
  convert_gff_to_bed.sh – Convert GFF to BED
  split_gff_by_feature.py – Split GFF by feature type
  process_bed_attributes.R – Parse BED attributes

5. Cleanup & Utilities
  03_cleanUp.sh – Remove intermediate files
  filtering.sh – Filtering utilities
  fixSampleNames.sh – Standardize sample naming

6. Downstream Analysis
  DESeq2.R – Differential gene expression analysis
  GO.R – Gene Ontology enrichment

🚀 Requirements
Depending on pipeline components:
  Linux environment
  R (DESeq2, tidyverse, clusterProfiler, etc.)
  Python 3 (pandas, argparse)
  STAR, samtools, Kallisto, MultiQC

Standard bioinformatics tools in $PATH

🎯 Purpose

This repository ensures full reproducibility of our transcriptomic analysis and supports research on how SARS-CoV-2 infection may impact placental biology and pregnancy outcomes.
