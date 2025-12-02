#!/bin/bash

saveCommand kallisto index \
    -i kallisto/index/kallisto_mRNA_index.idx \
    kallisto/input/mRNA.fasta --make-unique

echo "Indexing complete."