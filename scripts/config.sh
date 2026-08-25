#!/bin/bash

# Reference genome
REFERENCE="data/reference/tb_ref_genome/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna"

# Directories
TRIMMED_READS="results/trimmed"
SAM_FILES="results/alignment"
SORTED_BAM="results/sorted"
VCF_DIR="results/vcf"
CONSENSUS="results/consensus"
CHAIN="results/chain"

# Threads
THREADS=8

# GETTING DATA 
METADATA_FILE="metadata/CRyPTIC_reuse_table_20240917.csv"
SENSITIVE_IDS="metadata/sensitive_ids.txt"
RESISTANT_IDS="metadata/resistant_ids.txt"
TOTAL_IDS="metadata/all_ids.txt"

# Batches
BATCH_1_S="metadata/batch_1_S.txt"
BATCH_2_S="metadata/batch_2_S.txt"
BATCH_3_R="metadata/batch_3_R.txt"
BATCH_4_R="metadata/batch_4_R.txt"

