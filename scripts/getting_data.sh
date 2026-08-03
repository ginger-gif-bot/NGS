#!/bin/bash
set -e
# 1. Filtering the data
Metadata_file=$1
saved_acc_id=$2 
echo "Filtering the data"
awk -F "," 'NR>1 && $6 < 80000 {print $1}' "$Metadata_file" > "$saved_acc_id"
echo "Saved the accession ids in the "$saved_acc_id""

# 2. Prefetch
echo "Starting prefetch"

count=0
total=$(wc -l < "$saved_acc_id")

while read id
do 
    ((count +=1))
    if [ ! -f "data/raw_reads/raw_sra/${id}/${id}.sra" ]
    then 
        echo -e "\n[$count/$total] Downloading $id ...\n"
        prefetch "$id" -O data/raw_reads/raw_sra
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
        echo "[$count/$total] Downloading fastq $id ..."
        fasterq-dump data/raw_reads/raw_sra/"$id"/"$id".sra \
        --split-files \
        -O data/raw_reads/raw_fastq
    else
        echo "[$count/$total] already exists, skipping "$id""
    fi
done < $saved_acc_id
echo "All the fastq sequences have been downloaded"










