#!/bin/bash
exec > >(tee -a "logs/sorting_sam.log") 2>&1
echo "=== Run started $(date) ==="
set -e

source scripts/config.sh

total=$(ls "$SAM_FILES"/*.sam|wc -l )
count=0

for files in "$SAM_FILES"/*.sam
do
((count+=1))
id="${files%.sam}"
id="${id##*/}"
    if [ -s "$SORTED_BAM/${id}.sorted.bam" ]
    then
    echo "[${count}/${total}] Sorted file exist, skipping ${id} ..."
    else
        echo "Processing [${count}/${total}]..."
        start=$SECONDS
        samtools sort -@ "$THREADS" -o "$SORTED_BAM/${id}.sorted.bam.tmp" "$files" 
        mv "$SORTED_BAM/${id}.sorted.bam.tmp" "$SORTED_BAM/${id}.sorted.bam"
        samtools index "$SORTED_BAM/${id}.sorted.bam"
        echo "done in $((SECONDS - start))s"
    fi
done
