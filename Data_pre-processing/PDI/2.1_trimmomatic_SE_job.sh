#!/bin/bash
# module load Trimmomatic 

for i in *.gz; do
	trimmomatic SE -threads 30 $i Clean.$i ILLUMINACLIP:Adapter.fastq:2:40:15 SLIDINGWINDOW:4:20 MINLEN:27;
	echo " ... Done $i  ...."
done;
