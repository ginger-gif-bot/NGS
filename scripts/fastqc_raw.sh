#!/bin/bash
exec > >(tee -a "logs/fastqc_raw.log") 2>&1
echo "=== Run started $(date) ==="

set -e

echo "Starting FastQC ..."
count=0
total=$(ls data/raw_reads/raw_fastq/ | wc -l)

for file in data/raw_reads/raw_fastq/*.fastq
do 
    ((count+=1))
    base=$(basename "$file" .fastq)
    if [ ! -f "results/fastqc_raw/${base}_fastqc.zip" ] || \
       [ ! -f "results/fastqc_raw/${base}_fastqc.html" ]
    then
       echo -e "\n[$count/$total] Running FastQC on "$file""
       start=$SECONDS
       fastqc "$file" -o results/fastqc_raw || \
       { rm -f results/fastqc_raw/"${base}"_fastqc.*; exit 1; }
       echo -e "\ndone in $((SECONDS - start))s\n"
    else
       echo "[$count/$total] already exists, skipping "$file""
    fi
done
echo -e "\n\nAnalysis complete" 






