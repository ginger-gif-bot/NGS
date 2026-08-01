#!/bin/bash

set -e

echo "Starting MultiQC..."
multiqc -o results/multiqc_raw results/fastqc_raw
echo "Done"