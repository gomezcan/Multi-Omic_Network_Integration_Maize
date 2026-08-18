#!/bin/bash

## Deactivate conda env
conda deactivate

## Activate environment
source venv/bin/activate

## Make a temporary directory
mkdir tmp

## Download and GFF and convert to GTF file using cufflinks
wget ftp://ftp.ensemblgenomes.org/pub/plants/current/gtf/zea_mays/Zea_mays.B73_RefGen_v4.50.gtf.gz -P tmp
gunzip tmp/Zea_mays.B73_RefGen_v4.50.gtf.gz
mv tmp/Zea_mays.B73_RefGen_v4.50.gtf $PWD/../data/06_annotation_data/
rm -rf tmp

## Activate the virtual envirionment (On maching with GPUs) set up using instructions here: https://github.com/broadinstitute/gtex-pipeline
source $HOME/venv/bin/activate

## Collapse gtf file to have gene level expression data
python3 $HOME/gtex-pipeline/gene_model/collapse_annotation.py \
$PWD/../data/06_annotation_data/Zea_mays.B73_RefGen_v4.50.gtf \
$PWD/../data/06_annotation_data/collapsed_Zea_mays.B73_RefGen_v4.50.gtf

## Convert expression raw TPM to a .gct file for various datasets
Rscript 02.1_make_gct_files_mazaheri.R
Rscript 02.1_make_gct_files_kremling.R

## Make Index chromosome file for vcf...Note this is going to run on machine with the GPUs
tabix -p vcf $HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_Missing_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz
tabix --list-chroms $HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_Missing_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz \
> PANZEA_WIDIV_chr_list.txt

## Apply filtering and normalization to expression datasets using gtex-pipeline
python3 $HOME/gtex-pipeline/qtl/src/eqtl_prepare_expression.py \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/001.rTPM_Kremling_GRoot_expression_175g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/001.ReadCount_Kremling_GRoot_expression_175g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/06_annotation_data/collapsed_Zea_mays.B73_RefGen_v4.50.gtf \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/001.taxa_lookup_GRoot_175g_agpv4.txt \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_chr_list.txt \
'01_GRoot_175genos' \
--tpm_threshold 0.1 \
--count_threshold 6 \
--sample_frac_threshold 0.2 \
--normalization_method tmm \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data

python3 $HOME/gtex-pipeline/qtl/src/eqtl_prepare_expression.py \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/002.rTPM_Kremling_LMAN_expression_183g_agpv4.gct  \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/002.ReadCount_Kremling_LMAN_expression_183g_agpv4.gct  \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/06_annotation_data/collapsed_Zea_mays.B73_RefGen_v4.50.gtf \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/002.taxa_lookup_LMAN_183g_agpv4.txt  \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_chr_list.txt \
'02_LMAN_183genos' \
--tpm_threshold 0.1 \
--count_threshold 6 \
--sample_frac_threshold 0.2 \
--normalization_method tmm \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data

python3 $HOME/gtex-pipeline/qtl/src/eqtl_prepare_expression.py \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/003.rTPM_Kremling_L3Tip_expression_180g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/003.ReadCount_Kremling_L3Tip_expression_180g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/06_annotation_data/collapsed_Zea_mays.B73_RefGen_v4.50.gtf \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/003.taxa_lookup_L3Tip_180g_agpv4.txt \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_chr_list.txt \
'03_L3Tip_180genos' \
--tpm_threshold 0.1 \
--count_threshold 6 \
--sample_frac_threshold 0.2 \
--normalization_method tmm \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data

python3 $HOME/gtex-pipeline/qtl/src/eqtl_prepare_expression.py \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/004.rTPM_Kremling_L3Base_expression_178g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/004.ReadCount_Kremling_L3Base_expression_178g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/06_annotation_data/collapsed_Zea_mays.B73_RefGen_v4.50.gtf \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/004.taxa_lookup_L3Base_178g_agpv4.txt \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_chr_list.txt \
'04_L3Base_178genos' \
--tpm_threshold 0.1 \
--count_threshold 6 \
--sample_frac_threshold 0.2 \
--normalization_method tmm \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data

python3 $HOME/gtex-pipeline/qtl/src/eqtl_prepare_expression.py \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/005.rTPM_Kremling_GShoot_expression_178g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/005.ReadCount_Kremling_GShoot_expression_178g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/06_annotation_data/collapsed_Zea_mays.B73_RefGen_v4.50.gtf \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/005.taxa_lookup_GShoot_178g_agpv4.txt \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_chr_list.txt \
'05_Gshoot_178genos' \
--tpm_threshold 0.1 \
--count_threshold 6 \
--sample_frac_threshold 0.2 \
--normalization_method tmm \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data

python3 $HOME/gtex-pipeline/qtl/src/eqtl_prepare_expression.py \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/006.rTPM_Kremling_Kern_expression_169g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/006.ReadCount_Kremling_Kern_expression_169g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/06_annotation_data/collapsed_Zea_mays.B73_RefGen_v4.50.gtf \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/006.taxa_lookup_Kern_169g_agpv4.txt \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_chr_list.txt \
'06_Kern_169genos' \
--tpm_threshold 0.1 \
--count_threshold 6 \
--sample_frac_threshold 0.2 \
--normalization_method tmm \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data

python3 $HOME/gtex-pipeline/qtl/src/eqtl_prepare_expression.py \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/007.rTPM_Kremling_LMAD_expression_131g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/007.ReadCount_Kremling_LMAD_expression_131g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/06_annotation_data/collapsed_Zea_mays.B73_RefGen_v4.50.gtf \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/007.taxa_lookup_LMAD_131g_agpv4.txt \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_chr_list.txt \
'07_LMAD_131genos' \
--tpm_threshold 0.1 \
--count_threshold 6 \
--sample_frac_threshold 0.2 \
--normalization_method tmm \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data

python3 $HOME/gtex-pipeline/qtl/src/eqtl_prepare_expression.py \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/02_mazaheri_2019/01_gct_files/001.rTPM_WIDIV_seedling_V1_expression_187g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/02_mazaheri_2019/01_gct_files/001.ReadCount_WIDIV_seedling_V1_expression_187g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/06_annotation_data/collapsed_Zea_mays.B73_RefGen_v4.50.gtf \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/02_mazaheri_2019/01_gct_files/001.taxa_lookup_187g_agpv4.txt \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_chr_list.txt \
'08_Seedling_V1_187genos' \
--tpm_threshold 0.1 \
--count_threshold 6 \
--sample_frac_threshold 0.2 \
--normalization_method tmm \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data

python3 $HOME/gtex-pipeline/qtl/src/eqtl_prepare_expression.py \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/02_mazaheri_2019/01_gct_files/002.rTPM_WIDIV_seedling_V1_expression_304g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/02_mazaheri_2019/01_gct_files/002.ReadCount_WIDIV_seedling_V1_expression_304g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/06_annotation_data/collapsed_Zea_mays.B73_RefGen_v4.50.gtf \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/02_mazaheri_2019/01_gct_files/002.taxa_lookup_304g_agpv4.txt \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_chr_list.txt \
'09_Seedling_V1_304genos' \
--tpm_threshold 0.1 \
--count_threshold 6 \
--sample_frac_threshold 0.2 \
--normalization_method tmm \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data

python3 $HOME/gtex-pipeline/qtl/src/eqtl_prepare_expression.py \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/02_mazaheri_2019/01_gct_files/003.rTPM_WIDIV_seedling_V1_expression_602g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/02_mazaheri_2019/01_gct_files/003.ReadCount_WIDIV_seedling_V1_expression_602g_agpv4.gct \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/06_annotation_data/collapsed_Zea_mays.B73_RefGen_v4.50.gtf \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/02_mazaheri_2019/01_gct_files/003.taxa_lookup_602g_agpv4.txt \
$HOME/Beast_mount/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_chr_list.txt \
'10_Seedling_V1_602genos' \
--tpm_threshold 0.1 \
--count_threshold 6 \
--sample_frac_threshold 0.2 \
--normalization_method tmm \
-o $HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data

## Calculate PEER factors for each dataset
for prefix in 01_GRoot_175genos 02_LMAN_183genos 03_L3Tip_180genos 04_L3Base_178genos 05_Gshoot_178genos 06_Kern_169genos 07_LMAD_131genos 08_Seedling_V1_187genos 09_Seedling_V1_304genos 10_Seedling_V1_602genos
do
num_peer=60
docker run --rm -v $HOME/Beast_mount/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data:/data -t broadinstitute/gtex_eqtl:V8 /bin/bash \
-c "Rscript /src/run_PEER.R /data/${prefix}.expression.bed.gz ${prefix} ${num_peer} --output_dir '/data/' "
done

## Subset individual datasets and filter for MAF and missing data for each  (On Beast server for this session)
Rscript 02.3_make_individual_beds.R

plink2 --make-bed \
--output-chr 26 \
--keep $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/003.taxa_lookup_602g_agpv4.txt \
--vcf $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/ref_fixed_widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.vcf.gz \
--maf 0.05 \
--geno 0.1 \
--max-alleles 2 \
--threads 18 \
--out $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/10_Seedling_V1_602genos

plink2 --make-bed \
--output-chr 26 \
--keep $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/001.taxa_lookup_187g_agpv4.txt \
--vcf $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_Missing_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz \
--maf 0.05 \
--geno 0.1 \
--threads 18 \
--out $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/08_Seedling_V1_187genos

plink2 --make-bed \
--output-chr 26 \
--keep $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/002.taxa_lookup_304g_agpv4.txt \
--vcf $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_Missing_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz \
--maf 0.05 \
--geno 0.1 \
--threads 18 \
--out $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/09_Seedling_V1_304genos

plink2 --make-bed \
--output-chr 26 \
--keep $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/001.taxa_lookup_GRoot_175g_agpv4.txt \
--vcf $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_Missing_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz \
--maf 0.05 \
--geno 0.1 \
--threads 18 \
--out $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/01_GRoot_175genos

plink2 --make-bed \
--output-chr 26 \
--keep $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/002.taxa_lookup_LMAN_183g_agpv4.txt \
--vcf $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_Missing_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz \
--maf 0.05 \
--geno 0.1 \
--threads 18 \
--out $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/02_LMAN_183genos

plink2 --make-bed \
--output-chr 26 \
--keep $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/003.taxa_lookup_L3Tip_180g_agpv4.txt \
--vcf $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_Missing_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz \
--maf 0.05 \
--geno 0.1 \
--threads 18 \
--out $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/03_L3Tip_180genos

plink2 --make-bed \
--output-chr 26 \
--keep $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/004.taxa_lookup_L3Base_178g_agpv4.txt \
--vcf $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_Missing_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz \
--maf 0.05 \
--geno 0.1 \
--threads 18 \
--out $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/04_L3Base_178genos

plink2 --make-bed \
--output-chr 26 \
--keep $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/005.taxa_lookup_GShoot_178g_agpv4.txt \
--vcf $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_Missing_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz \
--maf 0.05 \
--geno 0.1 \
--threads 18 \
--out $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/05_Gshoot_178genos

plink2 --make-bed \
--output-chr 26 \
--keep $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/006.taxa_lookup_Kern_169g_agpv4.txt \
--vcf $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_Missing_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz \
--maf 0.05 \
--geno 0.1 \
--threads 18 \
--out $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/06_Kern_169genos

plink2 --make-bed \
--output-chr 26 \
--keep $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/007.taxa_lookup_LMAD_131g_agpv4.txt \
--vcf $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_Missing_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz \
--maf 0.05 \
--geno 0.1 \
--threads 18 \
--out $HOME/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/07_LMAD_131genos

## Calculate MAF for each dataset
plink2 --bfile 01_GRoot_175genos --freq --out 01_GRoot_175genos
plink2 --bfile 02_LMAN_183genos --freq --out 02_LMAN_183genos
plink2 --bfile 03_L3Tip_180genos --freq --out 03_L3Tip_180genos
plink2 --bfile 04_L3Base_178genos --freq --out 04_L3Base_178genos
plink2 --bfile 05_Gshoot_178genos --freq --out 05_Gshoot_178genos
plink2 --bfile 06_Kern_169genos --freq --out 06_Kern_169genos
plink2 --bfile 07_LMAD_131genos --freq --out 07_LMAD_131genos
plink2 --bfile 09_Seedling_V1_304genos --freq --out 09_Seedling_V1_304genos

