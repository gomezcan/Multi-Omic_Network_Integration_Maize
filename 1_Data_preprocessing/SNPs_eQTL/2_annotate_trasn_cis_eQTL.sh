#!/bin/bash

# Usage: 2_annotate_trasn_cis_eQTL.sh <eQTL.bed> <out_gene_overlap.bed> <out_peak_distance.bed>
#   $1  eQTL BED (chr, start, end, snp, ...)
#   $2  output: eQTLs annotated by gene-body overlap (intersectBed vs. gene BED)
#   $3  output: eQTLs annotated by distance to the nearest PDI peak summit (closestBed)

intersectBed -a $1 -b Zea_mays.B73_RefGen_v4.46.bed -wao > $2;
closestBed -a   $1 -b All.Summit_10.2020.bed -D 'ref' > $3;
