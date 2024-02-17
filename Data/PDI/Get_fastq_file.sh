#!/bin/bash

#################################
#####   Requare modules:     ####
# module load sra-toolkit       #
# module load parallel          #
#################################

## Input file requare at least two columns: SRR if in first column and Lybrary Layout in second (Single or Paired):
#  SRA	LibraryLayout
#  SRR8524980	Paired
#  SRR8524981	Paired
#  SRR504832	Single
#  SRR504833	Single

Input=$1

donwload_single() {
  # download SE files
  i=$(echo $1 | cut -f1)
  echo "$i"
  fastq-dump --gzip $Input;
  mv ${1}.fastq.gz RawReads/SE/;
}

donwload_paired() {
  # download PE files
  i=$(echo $1 | cut -f1)
  echo "$i"
  fastq-dump --gzip --split-files $Input;
  mv ${1}_1.fastq.gz RawReads/PE/;
  mv ${1}_2.fastq.gz RawReads/PE/;
}

export -f donwload_single
export -f donwload_paired

# Create out directory
makedir RawReads
makedir RawReads/PE
makedir RawReads/SE

# Subset SE runs to into tem file
grep -w 'Single' $1 | cut -f1 - > tem_SE;

# Subset PE runs to into tem file
grep -w 'Paired' $1 | cut -f1 - > tem_PE;

# Call downloading function
parallel -j 50 donwload_single :::: tem_SE;

# Call downloading function
parallel -j 50 donwload_paired :::: tem_PE;

# remove tem file
rm tem_SE;
rm tem_PE;




