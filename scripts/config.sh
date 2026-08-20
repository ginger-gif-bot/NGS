#!/bin/bash

# Reference genome
REFERENCE="data/reference/e_coli_ref_genome/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"

# Directories
TRIMMED_READS="results/trimmed"
SAM_FILES="results/alignment"
SORTED_BAM="results/sorted"
VCF_DIR="results/vcf"
CONSENSUS="results/consensus"
CHAIN="results/chain"

# Threads
THREADS=8
