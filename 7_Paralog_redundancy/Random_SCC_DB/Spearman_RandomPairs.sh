#!/bin/bash

# conda activate /maindisk/fabio/miniconda3

Call_sccPair() {
  
  Rscript Spearman_RandomPairs.R $1;
  #echo " ..  Donde $1 .."
}

export -f Call_sccPair

# Call downloading function
parallel -j 50 Call_sccPair :::: $1;

#conda deactivate
