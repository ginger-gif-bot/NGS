#!/bin/bash
exec > >(tee -a "logs/consenus.log") 2>&1
echo "=== Run started $(date) ==="
set -e 

source scripts/config.sh

total=$(ls "$VCF_DIR"/*.gz | wc -l)
count=0

for files in "$VCF_DIR"/*.gz
do
id="${files%.gz}"
id="${id%.vcf}"
id="${id##*/}"
((count +=1))
    if [ -s "$CONSENSUS/${id}.fasta" ]
    then
        echo "Consensus exists, skipping ${id} ..."
    else
        start=$SECONDS
        echo "Processing [${count}/${total}] ... "
        bcftools consensus -f "$REFERENCE" --chain "$CHAIN/${id}.chain" "$files" > "$CONSENSUS/${id}.fasta.tmp"
        mv "$CONSENSUS/${id}.fasta.tmp" "$CONSENSUS/${id}.fasta"
        echo -e "\ndone in $((SECONDS - start))s\n"
    fi
done
