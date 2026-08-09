#!/bin/bash
###
# This script filter out multi-mapping reads and remove duplicate reads
####

for i in B73.*.bam; do
Out=${i//B73_Mo17.Bowtie2./} 

samtools view -@ 50 -h -b -q 30 $i | samtools sort -@ 50 - -o Q30.$Out;
#samtools view -@ 30 -h $i | LC_ALL=C grep -v 'XS:i' - | samtools view -@ 30 -h -bS - | samtools sort -@ 30 - -o Uniq.$Out;
echo "... Done unique filter  ...."

picard MarkDuplicates I=Q30.$Out O=DeDup.Q30.$Out M=Metrics.Q30.${Out//.bam/.txt} REMOVE_DUPLICATES=true;
#picard MarkDuplicates I=Uniq.$Out O=DeDup.Uniq.$Out M=Metrics.Uniq.${Out//.bam/.txt} REMOVE_DUPLICATES=true;

echo "... Done Remove Duplicates ...."

done;
