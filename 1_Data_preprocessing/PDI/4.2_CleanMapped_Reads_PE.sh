#!/bin/bash
###
# This script filter out low quality MAPQ<30 (multi-mapping) reads and remove duplicate reads
####

for i in *.bam; do
	Out=${i//B73.Clean./B73.} 

	samtools view -@ 30 -h -b -q 30 $i | samtools sort -@ 30 - -o Q30.$Out;
	echo "... Done unique filter  ...."

	picard MarkDuplicates I=Q30.$Out O=DeDup.Q30.$Out M=Metrics.Q30.${Out//.bam/.txt} REMOVE_DUPLICATES=true;
	echo "... Done Remove Duplicates ...."
done;
