#!bin/bash
exec > >(tee -a "logs/msa_tree.log") 2>&1
echo "=== Run started $(date) ==="
set -e

for gene in rpoB katG gyrA
do
mkdir -p results/tree/${gene}
echo -e "\nProcessing ${gene}...\n"
mafft --auto results/cds_seq/combined_files/${gene}_combined.fasta > results/msa/${gene}_aligned.fasta
iqtree -s results/msa/${gene}_aligned.fasta -m MFP -bb 1000 -nt AUTO -pre results/tree/${gene}/${gene}_tree
done