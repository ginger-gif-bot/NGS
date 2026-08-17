#!/bin/bash
exec > >(tee -a "logs/multiqc_raw.log") 2>&1
echo "=== Run started $(date) ==="
set -e

echo "Starting MultiQC..."
multiqc -o results/multiqc_raw results/fastqc_raw
echo "Done"



