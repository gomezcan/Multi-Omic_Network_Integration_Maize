#!/bin/bash

for i in *.gz; do
	trimmomatic SE -threads 50 $i Clean.$i ILLUMINACLIP:Adapter.fastq:2:40:15 SLIDINGWINDOW:4:20 MINLEN:30;
done;
