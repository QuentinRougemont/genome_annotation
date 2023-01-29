#!/bin/bash

#Author: QR
#Date 23-01-23
#Purposes: script to merge all bam from the gsnap run

cd 04_mapped
samtools merge all.merged.bam *.sorted.bam
samtools sort -n all.merged.bam -o all.merged.sorted.bam

filterBam --uniq --paired --pairwiseAlignment --in all.merged.sorted.bam  --out all.merged.uniq.bam

samtools sort all.merged.uniq.bam -o Aligned.out.ss.bam  

#cleanup:
rm all.merged.bam
rm all.merged.sorted.bam 
rm all.merged.uniq.bam
