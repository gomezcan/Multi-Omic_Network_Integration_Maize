#!/bin/bash
# 
# ## Deactivate conda env
# conda deactivate
# 
# ## Activate environment
# source venv/bin/activate
# 
# # Watch the GPUs
# watch -n 1 nvidia-smi

## Run permutations for eQTL without including any covariates
tensorqtl \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/10_Seedling_V1_602genos \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/10_Seedling_V1_602genos.expression.bed.gz \
'10_Seedling_V1_602genos' \
--mode cis \
--window 100000 \
--seed 123456 \
--output_text \
--load_split \
--qvalue_lambda 0 \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/

tensorqtl \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/08_Seedling_V1_187genos.bed \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/08_Seedling_V1_187genos.expression.bed.gz \
'08_Seedling_V1_187genos' \
--mode cis \
--window 100000 \
--seed 123456 \
--output_text \
--load_split \
--qvalue_lambda 0 \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/

tensorqtl \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/09_Seedling_V1_304genos.bed \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/09_Seedling_V1_304genos.expression.bed.gz \
'09_Seedling_V1_304genos' \
--mode cis \
--window 100000 \
--seed 123456 \
--output_text \
--load_split \
--qvalue_lambda 0 \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/

tensorqtl \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/01_GRoot_175genos.bed \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/01_GRoot_175genos.expression.bed.gz \
'01_GRoot_175genos' \
--mode cis \
--window 100000 \
--seed 123456 \
--output_text \
--load_split \
--qvalue_lambda 0 \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/

tensorqtl \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/02_LMAN_183genos.bed \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/02_LMAN_183genos.expression.bed.gz \
'02_LMAN_183genos' \
--mode cis \
--window 100000 \
--seed 123456 \
--output_text \
--load_split \
--qvalue_lambda 0 \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/

tensorqtl \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/03_L3Tip_180genos.bed \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/03_L3Tip_180genos.expression.bed.gz \
'03_L3Tip_180genos' \
--mode cis \
--window 100000 \
--seed 123456 \
--output_text \
--load_split \
--qvalue_lambda 0 \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/

tensorqtl \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/04_L3Base_178genos.bed \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/04_L3Base_178genos.expression.bed.gz \
'04_L3Base_178genos' \
--mode cis \
--window 100000 \
--seed 123456 \
--output_text \
--load_split \
--qvalue_lambda 0 \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/

tensorqtl \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/05_Gshoot_178genos.bed \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/05_Gshoot_178genos.expression.bed.gz \
'05_Gshoot_178genos' \
--mode cis \
--window 100000 \
--seed 123456 \
--output_text \
--load_split \
--qvalue_lambda 0 \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/

tensorqtl \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/06_Kern_169genos.bed \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/06_Kern_169genos.expression.bed.gz \
'06_Kern_169genos' \
--mode cis \
--window 100000 \
--seed 123456 \
--output_text \
--load_split \
--qvalue_lambda 0 \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/

tensorqtl \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/07_LMAD_131genos.bed \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/07_LMAD_131genos.expression.bed.gz \
'07_LMAD_131genos' \
--mode cis \
--window 100000 \
--seed 123456 \
--output_text \
--load_split \
--qvalue_lambda 0 \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/results/02_eQTL_permutations/


