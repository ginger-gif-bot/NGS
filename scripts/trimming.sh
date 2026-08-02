#!/bin/bash
set -e

echo "Starting the trimming ..."
total=$(ls data/raw_reads/raw_fastq | wc -l )
# count=0
# for forward in data/raw_reads/raw_fastq/*_1.fastq
# do
# ((count+=1))
# reverse="${forward/_1.fastq/_2.fastq}"
# extract_id="${forward%_1.fastq}"
# extract_id="${extract_id##*/}"  # removes the file path
# extract_id="${extract_id%.sra}"
# echo -e "\nTrimming [$count/$total]\n" 
# trimmomatic PE \
#     "$forward" \
#     "$reverse" \
#     "results/fastqc_trimmed/${extract_id}_1_paired.fastq" \
#     "results/fastqc_trimmed/${extract_id}_1_unpaired.fastq" \
#     "results/fastqc_trimmed/${extract_id}_2_paired.fastq" \
#     "results/fastqc_trimmed/${extract_id}_2_unpaired.fastq" \
#     ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 \
#     LEADING:3 \
#     TRAILING:3 \
#     SLIDINGWINDOW:4:15 \
#     MINLEN:36
# done
# echo "Done" 

trimmomatic PE \
    data/raw_reads/raw_fastq/ERR13230455.sra_1.fastq \
    data/raw_reads/raw_fastq/ERR13230455.sra_2.fastq \
    out1.fastq out2.fastq out3.fastq out4.fastq \
    ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36




