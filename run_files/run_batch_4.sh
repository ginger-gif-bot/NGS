#!/bin/bash
set -e
# Change the current batch in getting_data.sh before running
exec > >(tee -a "logs/batch_2.log") 2>&1

echo "==== BATCH 4 STARTED $(date) ===="

echo -e "\n--- Starting: getting_data.sh ---\n"
conda run -n ngs_env bash scripts/getting_data.sh
echo "--- Done: getting_data.sh ---"

echo -e "\n--- Starting: fastqc_raw.sh ---\n"
conda run -n qc_env bash scripts/fastqc_raw.sh
echo "--- Done: fastqc_raw.sh ---"

echo -e "\n--- Starting: multiqc_raw.sh ---\n"
conda run -n qc_env bash scripts/multiqc_raw.sh
echo "--- Done: multiqc_raw.sh ---"

echo -e "\n--- Starting: trimming.sh ---\n"
conda run -n trim_env bash scripts/trimming.sh
echo "--- Done: trimming.sh ---"

echo -e "\n--- Starting: fastqc_trim.sh ---\n"
conda run -n qc_env bash scripts/fastqc_trim.sh
echo "--- Done: fastqc_trim.sh ---"

echo -e "\n--- Starting: multiqc_trim.sh ---\n"
conda run -n qc_env bash scripts/multiqc_trim.sh
echo "--- Done: multiqc_trim.sh ---"

echo -e "\n--- Starting: bwa_alignment.sh ---\n"
conda run -n ngs_env bash scripts/bwa_alignment.sh
echo "--- Done: bwa_alignment.sh ---"

echo -e "\n--- Starting: sorting_sam.sh ---\n"
conda run -n ngs_env bash scripts/sorting_sam.sh
echo "--- Done: sorting_sam.sh ---"

echo -e "\n--- Starting: variant_calling.sh ---\n"
conda run -n bcftools_env bash scripts/variant_calling.sh
echo "--- Done: variant_calling.sh ---"

echo -e "\n--- Starting: consensus.sh ---\n"
conda run -n bcftools_env bash scripts/consensus.sh
echo "--- Done: consensus.sh ---"

echo "==== BATCH 4 COMPLETE $(date) ===="