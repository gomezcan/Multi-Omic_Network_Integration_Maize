# =============================================================================
# EDITORIAL NOTE (added 2026-08-18 when archiving for publication; not by the
# original author). This file was run interactively, block by block, and is not
# runnable end-to-end as written: genotypeio.PlinkReader is called only for
# datasets 10 and 09, yet the later per-tissue blocks call
# trans.map_permutations(genotype_df, ...) without reloading genotype_df. Reload
# the matching PLINK prefix in each block before re-running. The permutation
# settings themselves (nperms=10000, maf_threshold=0.05, seed=123456) are the
# production ones, and the resulting per-tissue threshold is the 5th percentile
# of the permuted minimum-p distribution.
# =============================================================================
# !/bin/bash
# 
# ## Deactivate conda env
# conda deactivate
# 
# ## Activate environment
# source venv/bin/activate
# 
# ## Look for help page
# python3 -m tensorqtl --help
# 
# ## Run your commands within python environment
# cd tensorqtl/tensorqtl

import torch
from torch.utils import data
import numpy as np
import pandas as pd
import scipy.stats as stats
from collections import OrderedDict
import sys
import os
import time
import tensorqtl
import trans

# sys.path.insert(1, os.path.dirname(__file__))
import genotypeio
from core import *
  
## Run permutations for eQTL without including any covariates for trans associations
plink_prefix_path='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/10_Seedling_V1_602genos'
expression_bed='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/10_Seedling_V1_602genos.expression.bed.gz'
covariates_file='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/10_Seedling_V1_602genos.PEER_covariates_null.txt'

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# All genes Emperical p-values based on permutations Note: Some do not have SNPs in the window
trans_df = trans.map_permutations(genotype_df, covariates_df, permutations=None,chr_s=None, nperms=10000, maf_threshold=0.05, batch_size=20000, logger=None, seed=123456, verbose=True)
nperms = len(trans_df['minp_true_df'])
pval_perm = trans_df['minp_empirical']

# Write the dataframe to a text file
pd.DataFrame(pval_perm).to_csv('/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/10_Seedling_V1_602genos_trans_qtl.txt', sep='\t',index=True)


## Now the other datasets
plink_prefix_path='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/09_Seedling_V1_304genos'
expression_bed='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/09_Seedling_V1_304genos.expression.bed.gz'
covariates_file='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/09_Seedling_V1_304genos.PEER_covariates_null.txt'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
trans_df = trans.map_permutations(genotype_df, covariates_df, permutations=None,chr_s=None, nperms=10000, maf_threshold=0.05, batch_size=20000, logger=None, seed=123456, verbose=True)
nperms = len(trans_df['minp_true_df'])
pval_perm = trans_df['minp_empirical']
pd.DataFrame(pval_perm).to_csv('/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/09_Seedling_V1_304genos_trans_qtl.txt', sep='\t',index=True)

expression_bed='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/08_Seedling_V1_187genos.expression.bed.gz'
covariates_file='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/08_Seedling_V1_187genos.PEER_covariates_null.txt'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
trans_df = trans.map_permutations(genotype_df, covariates_df, permutations=None,chr_s=None, nperms=10000, maf_threshold=0.05, batch_size=20000, logger=None, seed=123456, verbose=True)
nperms = len(trans_df['minp_true_df'])
pval_perm = trans_df['minp_empirical']
pd.DataFrame(pval_perm).to_csv('/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/08_Seedling_V1_187genos_trans_qtl.txt', sep='\t',index=True)


expression_bed='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/07_LMAD_131genos.expression.bed.gz'
covariates_file='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/07_LMAD_131genos.PEER_covariates_null.txt'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
trans_df = trans.map_permutations(genotype_df, covariates_df, permutations=None,chr_s=None, nperms=10000, maf_threshold=0.05, batch_size=20000, logger=None, seed=123456, verbose=True)
nperms = len(trans_df['minp_true_df'])
pval_perm = trans_df['minp_empirical']
pd.DataFrame(pval_perm).to_csv('/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/07_LMAD_131genos_trans_qtl.txt', sep='\t',index=True)

expression_bed='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/06_Kern_169genos.expression.bed.gz'
covariates_file='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/06_Kern_169genos.PEER_covariates_null.txt'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
trans_df = trans.map_permutations(genotype_df, covariates_df, permutations=None,chr_s=None, nperms=10000, maf_threshold=0.05, batch_size=20000, logger=None, seed=123456, verbose=True)
nperms = len(trans_df['minp_true_df'])
pval_perm = trans_df['minp_empirical']
pd.DataFrame(pval_perm).to_csv('/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/06_Kern_169genos_trans_qtl.txt', sep='\t',index=True)

expression_bed='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/05_Gshoot_178genos.expression.bed.gz'
covariates_file='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/05_Gshoot_178genos.PEER_covariates_null.txt'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
trans_df = trans.map_permutations(genotype_df, covariates_df, permutations=None,chr_s=None, nperms=10000, maf_threshold=0.05, batch_size=20000, logger=None, seed=123456, verbose=True)
nperms = len(trans_df['minp_true_df'])
pval_perm = trans_df['minp_empirical']
pd.DataFrame(pval_perm).to_csv('/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/05_Gshoot_178genos_trans_qtl.txt', sep='\t',index=True)

expression_bed='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/04_L3Base_178genos.expression.bed.gz'
covariates_file='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/04_L3Base_178genos.PEER_covariates_null.txt'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
trans_df = trans.map_permutations(genotype_df, covariates_df, permutations=None,chr_s=None, nperms=10000, maf_threshold=0.05, batch_size=20000, logger=None, seed=123456, verbose=True)
nperms = len(trans_df['minp_true_df'])
pval_perm = trans_df['minp_empirical']
pd.DataFrame(pval_perm).to_csv('/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/04_L3Base_178genos_trans_qtl.txt', sep='\t',index=True)

expression_bed='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/03_L3Tip_180genos.expression.bed.gz'
covariates_file='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/03_L3Tip_180genos.PEER_covariates_null.txt'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
trans_df = trans.map_permutations(genotype_df, covariates_df, permutations=None,chr_s=None, nperms=10000, maf_threshold=0.05, batch_size=20000, logger=None, seed=123456, verbose=True)
nperms = len(trans_df['minp_true_df'])
pval_perm = trans_df['minp_empirical']
pd.DataFrame(pval_perm).to_csv('/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/03_L3Tip_180genos_trans_qtl.txt', sep='\t',index=True)

expression_bed='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/02_LMAN_183genos.expression.bed.gz'
covariates_file='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/02_LMAN_183genos.PEER_covariates_null.txt'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
trans_df = trans.map_permutations(genotype_df, covariates_df, permutations=None,chr_s=None, nperms=10000, maf_threshold=0.05, batch_size=20000, logger=None, seed=123456, verbose=True)
nperms = len(trans_df['minp_true_df'])
pval_perm = trans_df['minp_empirical']
pd.DataFrame(pval_perm).to_csv('/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/02_LMAN_183genos_trans_qtl.txt', sep='\t',index=True)

expression_bed='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/01_GRoot_175genos.expression.bed.gz'
covariates_file='/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/01_GRoot_175genos.PEER_covariates_null.txt'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
trans_df = trans.map_permutations(genotype_df, covariates_df, permutations=None,chr_s=None, nperms=10000, maf_threshold=0.05, batch_size=20000, logger=None, seed=123456, verbose=True)
nperms = len(trans_df['minp_true_df'])
pval_perm = trans_df['minp_empirical']
pd.DataFrame(pval_perm).to_csv('/home/<chtc-user>/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/01_GRoot_175genos_trans_qtl.txt', sep='\t',index=True)
