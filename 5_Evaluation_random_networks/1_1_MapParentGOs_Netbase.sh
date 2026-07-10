#!/bin/bash

# conda activate /maindisk/fabio/miniconda3

Call_1_1() {
  
  Rscript 1_1_MapParentGOs_Netbase.R $1;
  echo " ..  Donde $1 .."
}

export -f Call_1_1

# Call downloading function
parallel -j 32 Call_1_1 :::: $1;

#conda deactivate
