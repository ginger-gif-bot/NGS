#!/bin/bash
set -e 

SAMPLE=$1

if [ -z "$1" ]
then 
echo "Usage: "$0" <filename>"
exit 1
fi 
echo "Starting the NGS pipeline"

# 1.
echo "Step 0: Obtain the data"
# Raw sequencing genome(FASTQ) --> NCBI SRA
# Reference genome (FASTA) --> NCBI Assembly
# Both require different tools
# esearch --> efetch --> prefetch --> fasterq-dump
# 1. 
# esearch -db sra -query "Escherichia coli[Organism]
# AND ILLUMINA[Platform] AND \
# PAIRED[Layout] AND \
# public[Access]" -retmax 20 \
# | efetch -format runinfo \
# | awk -F "," 'NR > 1 && 490 < $8 && $8 < 500 {print $1, $8}'
# | wc -l
# search not doing, because its taking too long and its crashing

echo "Step 1: Reading the FASTQ file"
# Input --> fastq file
# FastQC produces an HTML report and a zip file
# save the report and the zip file in qc folder 
echo "Saving the Report and zip file in qc folder"

