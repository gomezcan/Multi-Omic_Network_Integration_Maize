#!/bin/bash

# conda activate /maindisk/fabio/miniconda3

Call_2_2() {
  
  Rscript 2_2_MapParentGOs_CommonFunct.R $1;
  echo " ..  Donde $1 .."
}

export -f Call_2_2

# Call downloading function
parallel -j 32 Call_2_2 :::: $1;

#conda deactivate
