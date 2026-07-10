#!/bin/bash

# conda activate /maindisk/fabio/miniconda3

Call_3_1() {
  
  Rscript 3_1_MapParentGOs_Common_targets.R $1;
  echo " ..  Donde $1 .."
}

export -f Call_3_1

# Call downloading function
parallel -j 40 Call_3_1 :::: $1;

#conda deactivate
