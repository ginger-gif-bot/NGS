#!/bin/bash
set -e

echo "Starting FastQC ..."
count=0
total=$(ls results/trimmed| wc -l)
for file in results/trimmed/*.fastq
do 
    ((count+=1))
    base=$(basename "$file" .fastq)
    if [ ! -f "results/fastqc_trimmed/${base}_fastqc.zip" ] || \
       [ ! -f "results/fastqc_trimmed/${base}_fastqc.html" ]
    then
       echo -e "\n[$count/$total] Running FastQC on "$file""
       fastqc "$file" -o results/fastqc_trimmed
    else
       echo "[$count/$total] already exists, skipping "$file""
    fi
done
echo -e "\n\nAnalysis complete" 








