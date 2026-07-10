#!/bin/bash

# conda activate /maindisk/fabio/miniconda3

Call_pepDist() {
  
  Rscript 2_Pep_distance_Calculation.R $1;
  echo " ..  Donde $1 .."
}

export -f Call_pepDist

# Call downloading function
parallel -j 50 Call_pepDist :::: $1;

#conda deactivate
