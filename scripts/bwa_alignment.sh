#!/bin/bash
set -e

source scripts/config.sh

LOG_DIR="logs/bwa"
OUTPUT_DIR="results/alignment"

total=$(ls "$TRIMMED_READS"/*_1_paired.fastq | wc -l)
count=0

if [ ! -f "${REFERENCE}.bwt" ]
then
    echo "Index not found -> run scripts/reference.sh first"
    exit 1
fi
    
for reads in "$TRIMMED_READS"/*_1_paired.fastq
do 
((count+=1))
R1="$reads"
R2="${reads/_1_paired.fastq/_2_paired.fastq}"
id="${reads%_1_paired.fastq}"
id="${id##*/}"

if [[ -f "$R1" && -f "$R2" ]]
then
    if [ -s "$OUTPUT_DIR/${id}.sam" ]
    then
        echo "Output already exists -> skipping ${id}..."
    else
        echo "Processing [${count}/${total}]..."
        start=$SECONDS
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${id}\tSM:${id}\tPL:ILLUMINA" \
            "$REFERENCE" "$R1" "$R2" \
            > "$OUTPUT_DIR/${id}.sam.tmp" \
            2> "$LOG_DIR/${id}.log"

        mv "$OUTPUT_DIR/${id}.sam.tmp" "$OUTPUT_DIR/${id}.sam"
        echo "done in $((SECONDS - start))s"
    fi
else
    echo "R1 or R2 is missing -> skipping ${id}..."
fi
done 

