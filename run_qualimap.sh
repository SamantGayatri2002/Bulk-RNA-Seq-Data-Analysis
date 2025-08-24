#!/bin/bash
for bamfile in alignedreads/*.bam
do
    sample=$(basename "$bamfile" .bam)
      qualimap rnaseq \
        -bam "$bamfile" \
        -gtf Homo_sapiens.GRCh38.114.gtf \
        -outdir rnaseq_qc_"$sample" \
        --java-mem-size=8G
done
