#!/bin/bash
exec > >(tee -a "logs/fastqc_trim.log") 2>&1
echo "=== Run started $(date) ==="
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
       start=$SECONDS
       echo -e "\n[$count/$total] Running FastQC on "$file""
       fastqc "$file" -o results/fastqc_trimmed || \
       { rm -f results/fastqc_trimmed/"${base}"_fastqc.*; exit 1; }
       echo "done in $((SECONDS - start))s"
    else
       echo "[$count/$total] already exists, skipping "$file""
    fi
done
echo -e "\n\nAnalysis complete" 








