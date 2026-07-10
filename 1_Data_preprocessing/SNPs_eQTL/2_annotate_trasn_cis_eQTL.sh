#!/bin/bash


# -a eQTL_FarmCPU.bed
# $2 eQTL_FarmCPU_Annotated_trans.bed
# $3 eQTL_FarmCPU_Annotated_cis.bed

intersectBed -a $1 -b Zea_mays.B73_RefGen_v4.46.bed -wao > $2;
closestBed -a   $1 -b All.Summit_10.2020.bed -D 'ref' > $3;

