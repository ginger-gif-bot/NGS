#!/bin/bash
exec > >(tee -a "logs/getting_data.log") 2>&1
echo "=== Run started $(date) ==="

set -e

source scripts/config.sh
# 1. Filtering the data

echo -e "Isolating the SRA ids...\n"
awk -F "," 'NR > 1 && $9 == "S" && $14 == "S" {print $1}' "$METADATA_FILE" | head -50 > "$SENSITIVE_IDS"
echo "Saved the Sensitive accession ids in the "$SENSITIVE_IDS" file."
awk -F "," 'NR > 1 && $9 == "R" && $14 == "R"  {print $1}' "$METADATA_FILE" | head -50 > "$RESISTANT_IDS"
echo "Saved the Resistant accession ids in the "$RESISTANT_IDS" file."

cat "$SENSITIVE_IDS" "$RESISTANT_IDS" > "$TOTAL_IDS"
echo "Saved all the accession ids in the "$TOTAL_IDS" file."

# 2. Prefetch
echo -e "\nStarting prefetch\n"

CURRENT_BATCH="$BATCH_1_S"
count=0
total=$(wc -l < "$CURRENT_BATCH")

while read id
do 
    ((count +=1))
    if [ ! -f "data/raw_reads/raw_sra/${id}/${id}.sra" ]
    then 
        start=$SECONDS
        echo -e "\n[$count/$total] Downloading $id ...\n"
        prefetch "$id" -O data/raw_reads/raw_sra
        echo -e "\ndone in $((SECONDS - start))s\n"
    else
    echo "[$count/$total] file already exists, skipping "$id"" 
    fi
done < "$CURRENT_BATCH"

echo "All the sra sequences have been downloaded"

# 3. Fasterq dump
echo -e "\nStarting fasterq-dump\n"

count=0

while read id
do
    ((count +=1))
    if [ ! -f "data/raw_reads/raw_fastq/${id}.sra_1.fastq" ] || \
       [ ! -f "data/raw_reads/raw_fastq/${id}.sra_2.fastq" ]
    then
        start=$SECONDS
        echo "[$count/$total] Downloading fastq $id ..."
        fasterq-dump data/raw_reads/raw_sra/"$id"/"$id".sra \
        --split-files \
        -O data/raw_reads/raw_fastq/ || \
        { rm -f data/raw_reads/raw_fastq/"${id}".sra_*.fastq; exit 1; }
        echo -e "\ndone in $((SECONDS - start))s\n"
        rm -rf "data/raw_reads/raw_sra/${id}/"
        echo -e "\nSRA files deleted for ${id}\n"
    else
        echo "[$count/$total] already exists, skipping "$id""
    fi
done < "$CURRENT_BATCH"
echo "All the fastq sequences have been downloaded"

# To check the column numbers
# head -1 metadata/CRyPTIC_reuse_table_20240917.csv | tr ',' '\n' | grep -n "BINARY"

# To check if some ids have same accesion ids
# comm -12 <(sort metadata/sensitive_ids.txt) <(sort metadata/resistant_ids.txt)
# echo "No same ids found"




