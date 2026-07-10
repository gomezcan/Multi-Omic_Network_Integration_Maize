#!/bin/bash

#==> 01_genome_wide_eQTL_FarmCPU_variable_snpPCs_2exprPCs_top100_per_gene.txt <==
#snp	snp_chr	snp_pos	ref	alt	effect	std_err	p.value	pass_bonf	gene	gene_chr	gene_start	gene_stop	snp_distance_to_gene	snp_relative_to_gene	snp_loc_relative_to_gene_loc

# ==> 01_phenolic_and_tf_genes_eQTL_FarmCPU_variable_snpPCs_2exprPCs_top100_per_gene.txt <==
#snp	tvalue	p.value	pass_bonf	stderr	estimate	snp_chr	snp_pos	gene	gene_start	gene_stop	gene_chr	snp_distance_to_gene	snp_relative_to_gene	snp_loc_relative_to_gene_loc

## Target columns: snp_chr snp_pos snp_pos snp  
#
# line for 01_phenolic_and_tf_genes 
# cat $1 | awk  -v OFS='\t' '{print $7, $8, $8, $1}' | grep -v 'snp' | sort -k1,1 -k2,2n > $2
# line for 01_genome_wide 
cat $1 | awk  -v OFS='\t' '{print $2, $3, $3, $1}' | grep -v 'snp' | sort -k1,1 -k2,2n > $2

