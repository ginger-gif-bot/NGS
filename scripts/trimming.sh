#!/bin/bash
set -e

for forward in data/raw_reads/raw_fastq/*_1.fastq
do
reverse="$(forward/_1.fastq/_2.fastq)"
extract_id="$(forward%_1.fastq)"
extract_id="$(extract_id##*/)"  # removes the file path
extract_id="$(extract_id%.sra)" 
trimmomatic PE \
    "$forward" \
    "$reverse" \
    "results/fastqc_trimmed/${extract_id}_1_paired.fastq" \
    "results/fastqc_trimmed/${extract_id}_1_unpaired.fastq" \
    "results/fastqc_trimmed/${extract_id}_2_paired.fastq" \
    "results/fastqc_trimmed/${extract_id}_2_unpaired.fastq" \
