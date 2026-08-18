#!/bin/bash

# 1: keep eQTLs supported by all 8 models (col 7 = n_models_support)
# 2: keep eQTLs within 50 kb of the target TSS (col 24 = snp_distance_to_TSS)
# NOTE 2026-08-18: restored the missing `print $0` in the distance filter — the
# previously committed awk had an empty action and emitted zero rows.
cat 00_filtered_cis_trans_all_eQTL_results.txt | awk -v OFS='\t' '{if($7==8) print $0}' | awk -v OFS='\t' '{if($24 <= 50000) print $0}' | cut -f1,2,3,4,8,20,22,24 \
	| awk -v OFS='\t' '{print $3,$4,$4,$2,$5,$6,$1}' > 00_clean.cis.eQTLs.bed

closestBed -a 00_clean.cis.eQTLs.bed -b All.Summit_10.2020.bed -D 'ref' | awk -v OFS='\t' '{if(sqrt($12*$12)<=20) print $0}'  | cut -f1,2,4,6,7,9,11,12 > cis_eQTL.pdi.network.txt
