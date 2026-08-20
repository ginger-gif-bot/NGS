#!/bin/bash

if samtools quickcheck -v results/sorted/*.sorted.bam 
then
    echo "ALL BAMs OK -> removing SAM files"
    read -p "Delete SAM files? (y/n)" confirm
    if [ "$confirm" = "y" ]
    then
        rm results/alignment/*.sam
        echo "Deleted"
    else
        echo "Skipped"
    fi
else 
    echo "BAM check failed -> keeping SAM files"
    exit 1
fi

# Don't run it if you're unsure about the process
# Delete only when you are way futher in the pipeline
# It deletes everything 
