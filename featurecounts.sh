#!/bin/bash

# Go to your aligned reads folder
cd /home/shobita/Bulk_RNA_Seq_Analysis/alignedreads

# Loop over all BAM files
for bam in *.bam; do
    start=$(date +%s)  # start time

    echo "Processing $bam ..."
    featureCounts -s 0 -a /home/shobita/Bulk_RNA_Seq_Analysis/Homo_sapiens.GRCh38.114.gtf \
        -o /home/shobita/Bulk_RNA_Seq_Analysis/quants/${bam%.bam}_featurecounts.txt \
        "$bam"

    end=$(date +%s)  # end time
    runtime=$(( (end - start) / 60 ))  # in minutes

    echo "✅ Completed $bam in $runtime minutes."
    echo "------------------------------------"
done
