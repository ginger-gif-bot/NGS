#!/bin/bash

set -e

echo "Starting MultiQC..."
multiqc -o results/multiqc_trimmed results/fastqc_trimmed
echo "Done"




