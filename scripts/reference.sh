#!/bin/bash
exec > >(tee -a "logs/reference.log") 2>&1
echo "=== Run started $(date) ==="
set -e

# Finding the reference genome ... -> fetch directly from NCBI
# Reference genome retrieval:
# 1. Go to the NCBI Genome database and search for "Mycobacterium tuberculosis".
# 2. Apply the filters:
#       - Reference genomes
#       - Annotated genomes
#       - Annotated by NCBI RefSeq
# 3. Select the resulting reference assembly:
#       - Mycobacterium tuberculosis H37Rv (strain)
#       - Assembly accession: GCF_000195955.2
#       - Chromosome accession: NC_000962.3
# 4. Select "Download Package".
# 5. In the download options, select:
#       - RefSeq only
#       - Genome sequences (FASTA)
#       - Sequence and annotation (GBFF)
#    and leave unnecessary reports/files unselected.
# 6. Download the package and extract/unzip it.
# 7. Use the genomic FASTA file:
#       GCF_000195955.2_ASM19595v2_genomic.fna
#    as the reference genome for the downstream BWA analysis.

source scripts/config.sh

echo -e "Preparing the reference for BWA...\n"
if [ -f "$REFERENCE" ]
then 
    echo -e "Reference file is present\n"
    if [ -f "$REFERENCE.bwt" ] && \
       [ -f "$REFERENCE.pac" ] && \
       [ -f "$REFERENCE.ann" ] && \
       [ -f "$REFERENCE.amb" ] && \
       [ -f "$REFERENCE.sa" ]
    then
        echo "BWA Index files are present. Skipping creation ..."
    else
        echo "BWA Index files are absent, creating them..."
        bwa index "$REFERENCE" || { rm -f "$REFERENCE".{amb,ann,bwt,pac,sa}; exit 1; }
    fi
else
    echo "Reference file is absent"
    exit 1
fi





