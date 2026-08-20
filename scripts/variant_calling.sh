#!/bin/bash
exec > >(tee -a "logs/variant_calling.log") 2>&1
echo "=== Run started $(date) ==="
set -e 
set -o pipefail

source scripts/config.sh

total=$(ls "$SORTED_BAM"/*.bam | wc -l)
count=0
for files in "$SORTED_BAM"/*.bam
do 
id="${files%.bam}"
id="${id%.sorted}"
id="${id##*/}"
((count +=1))
    if [ -s "$VCF_DIR/${id}.vcf.gz" ]
    then
        echo "VCF exists, skipping ${id} ..."
    else
        echo "Processing [${count}/${total}] ... "
        start=$SECONDS
        bcftools mpileup -f "$REFERENCE" -O u --threads "$THREADS" "$files" \
        | bcftools call -m -v -O z --threads "$THREADS" -o "$VCF_DIR/${id}.vcf.gz.tmp"
        mv "$VCF_DIR/${id}.vcf.gz.tmp" "$VCF_DIR/${id}.vcf.gz"
        bcftools index "$VCF_DIR/${id}.vcf.gz"
        echo -e "\ndone in $((SECONDS - start))s\n"
    fi
done
