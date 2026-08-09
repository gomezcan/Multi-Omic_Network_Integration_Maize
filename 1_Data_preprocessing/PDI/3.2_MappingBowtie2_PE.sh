#!/bin/bash

module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1

for i in Clean.*_1.fastq.gz; do
	# set inputs
	Input_1=$i
	Input_2=${i//_1.fastq.gz/_2.fastq.gz}
	# set outputs
	out=${i//_1.fastq.gz/.sam};
	out=${out//Clean./};
	# Un aling reads
	name=${i//_1.fastq.gz/}

	# Mapping 0 mismacht in seed
	echo "Mapping $i to B73 with default "
	bowtie2 -p 50 --al-conc-gz B73.un.$name.fq.gz --no-mixed --no-discordant --no-unal \
	-x ../../Index_B73v4.dna_bowtie2 -1 $Input_1 -2 $Input_2 -S B73.$out;
	
	# compress sam to bam and sort bam
	samtools view -@ 30 -h -bS B73.$out | samtools sort -@ 30 - -o B73.${out//.sam/.bam};
	rm B73.$out;
	echo ".... Donde ... "
done;
