#!/bin/bash

for i in *1.fastq.gz; do
	File1=$i;
	File2=${i//1.fastq.gz/2.fastq.gz};	
	trimmomatic PE -threads 30 $File1 $File2 Clean.$File1 UnPaired.$File1 Clean.$File2 UnPaired.$File2 ILLUMINACLIP:Adapter.fastq:2:40:15:1:TRUE SLIDINGWINDOW:4:20 MINLEN:27

done;
