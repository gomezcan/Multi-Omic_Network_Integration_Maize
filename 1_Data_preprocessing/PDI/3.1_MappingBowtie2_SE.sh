#!/bin/bash

module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1


Index=Index_B73v4.dna_bowtie2

for i in Clean*.fastq.gz; do
	out=${i//.fastq.gz/.sam};
	name=${i//.gz/}
	# Mapping 0 mismacht in seed
	echo "Mapping $i to B73 with default "
	bowtie2 -p 25 --un B73.un.$name --no-unal -x $Index -U $i -S B73.$out;
	samtools view -@ 25 -h -bS B73.$out | samtools sort -@ 25 - -o B73.${out//.sam/.bam};
	rm B73.$out;
	echo ".... Donde ... "
done;
