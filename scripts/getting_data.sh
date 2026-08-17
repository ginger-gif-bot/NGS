#!/bin/bash
exec > >(tee -a "logs/getting_data.log") 2>&1
echo "=== Run started $(date) ==="

set -e
# 1. Filtering the data
Metadata_file="metadata/SraRuntable.csv"
saved_acc_id="metadata/selected_sra_run_ids"
echo -e "Isolating the SRA ids...\n"
awk -F "," 'NR > 1 {print $1}' "$Metadata_file" > "$saved_acc_id"
echo "Saved the accession ids in the "$saved_acc_id" file."

# 2. Prefetch
echo -e "\nStarting prefetch\n"

count=0
total=$(wc -l < "$saved_acc_id")

while read id
do 
    ((count +=1))
    if [ ! -f "data/raw_reads/raw_sra/${id}/${id}.sra" ]
    then 
        start=$SECONDS
        echo -e "\n[$count/$total] Downloading $id ...\n"
        prefetch "$id" -O data/raw_reads/raw_sra
        echo "done in $((SECONDS - start))s"
    else
    echo "[$count/$total] file already exists, skipping "$id"" 
    fi
done < $saved_acc_id

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
        echo "done in $((SECONDS - start))s"
    else
        echo "[$count/$total] already exists, skipping "$id""
    fi
done < $saved_acc_id
echo "All the fastq sequences have been downloaded"










