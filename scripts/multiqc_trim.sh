#!/bin/bash
exec > >(tee -a "logs/multiqc_trim.log") 2>&1
echo "=== Run started $(date) ==="
set -e

echo "Starting MultiQC..."
multiqc -o results/multiqc_trimmed results/fastqc_trimmed
echo "Done"




