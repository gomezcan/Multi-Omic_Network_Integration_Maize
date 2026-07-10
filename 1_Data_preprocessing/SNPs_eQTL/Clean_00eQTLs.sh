#!/bin/bash

# 1: remove eQTLs without support in all 8 models
# 2: keep eQTLs presented 50 k around TSS
cat 00_filtered_cis_trans_all_eQTL_results.txt | awk -v OFS='\t' '{if($7==8) print $0}' | awk -v OFS='\t' '{if($24 <=50000)}' | cut -f1,2,3,4,8,20,22,24 \
	| awk -v OFS='\t' '{print $3,$4,$4,$2,$5,$6,$1}' > 00_clean.cis.eQTLs.bed

closestBed -a 00_clean.cis.eQTLs.bed -b All.Summit_10.2020.bed -D 'ref' | awk -v OFS='\t' '{if(sqrt($12*$12)<=20) print $0}'  | cut -f1,2,4,6,7,9,11,12 > cis_eQTL.pdi.network.txt
