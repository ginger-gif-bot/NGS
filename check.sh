#!/bin/bash

comm -12 <(sort metadata/sensitive_ids.txt) \
  <(ls data/raw_reads/raw_fastq/ | grep "_1" | sed 's/\.sra_1.*//' | sort) | wc -l

comm -12 <(sort metadata/resistant_ids.txt) \
  <(ls data/raw_reads/raw_fastq/ | grep "_1" | sed 's/\.sra_1.*//' | sort) | wc -l

mkdir -p data/raw_reads/holding
while read id
do
  mv data/raw_reads/raw_fastq/${id}.sra_*.fastq data/raw_reads/holding/
done < "metadata/batch_2_S.txt"