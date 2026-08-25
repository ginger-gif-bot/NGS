#!/bin/bash
set -e
source scripts/config.sh

head -25 "$SENSITIVE_IDS" > "metadata/batch_1_S.txt"
tail -25 "$SENSITIVE_IDS" > "metadata/batch_2_S.txt"

head -25 "$RESISTANT_IDS" > "metadata/batch_3_R.txt"
tail -25 "$RESISTANT_IDS" > "metadata/batch_4_R.txt"

duplicates_S=$(comm -12 <(sort metadata/batch_1_S.txt) <(sort metadata/batch_2_S.txt))
if [ -z "$duplicates_S" ]
then
    echo "No same ids found for sensitive batches"
else
    echo "WARNING!! Duplicates ids found"
    echo "$duplicates_S"
fi

duplicates_R=$(comm -12 <(sort metadata/batch_3_R.txt) <(sort metadata/batch_4_R.txt))
if [ -z "$duplicates_R" ]
then
    echo "No same ids found for resistant batches"
else
    echo "WARNING!! Duplicates ids found"
    echo "$duplicates_R"
fi

while read id
do
    if [ -f "data/raw_reads/raw_fastq/${id}.sra_1.fastq" ]
    then 
        echo "The file with id ${id} exists, deleting it ..."
        rm -f "data/raw_reads/raw_fastq"/${id}*.fastq
        echo "The file with id ${id} has been deleted."
    else
        echo "File with id ${id} does not exist."
    fi
done < "$RESISTANT_IDS"