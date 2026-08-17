#!/bin/bash
exec > >(tee -a "logs/trimming.log") 2>&1
echo "=== Run started $(date) ==="

echo "Starting the trimming ..."
total=$(ls data/raw_reads/raw_fastq | wc -l )
count=0

> results/failed_samples.txt

for forward in data/raw_reads/raw_fastq/*_1.fastq
do
((count+=1))
reverse="${forward/_1.fastq/_2.fastq}"
extract_id="${forward%_1.fastq}"
extract_id="${extract_id##*/}"  # removes the file path
extract_id="${extract_id%.sra}"
echo -e "\nTrimming [$count/$total]\n" 
start=$SECONDS
trimmomatic PE \
    -phred33 \
    "$forward" \
    "$reverse" \
    "results/trimmed/${extract_id}_1_paired.fastq" \
    "results/trimmed/${extract_id}_1_unpaired.fastq" \
    "results/trimmed/${extract_id}_2_paired.fastq" \
    "results/trimmed/${extract_id}_2_unpaired.fastq" \
    ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36 || {
        echo "WARNING: Trimmomatic failed on $extract_id - skipping"
        echo "$extract_id" >> results/failed_samples.txt
        continue 
    }
    echo "done in $((SECONDS - start))s"
done
echo "Done" 












