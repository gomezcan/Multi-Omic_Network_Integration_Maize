#!/bin/bash

## Download Panzea SNPs 
iget -K -r -N 12 /iplant/home/shared/panzea/hapmap3/hmp321/imputed/uplifted_APGv4  $PWD/../data/

## Combine into one large vcf and keep only lines overlapping with WIDIV and Panzea and in metabolite profiling experiment.
## Keep only biallelic SNPs, Filter for MAF > 0.05 and set heterozygotes to NA.
bcftools concat $PWD/../data/uplifted_APGv4/hmp321_agpv4_chr1.vcf.gz \
$PWD/../data/uplifted_APGv4/hmp321_agpv4_chr2.vcf.gz \
$PWD/../data/uplifted_APGv4/hmp321_agpv4_chr3.vcf.gz \
$PWD/../data/uplifted_APGv4/hmp321_agpv4_chr4.vcf.gz \
$PWD/../data/uplifted_APGv4/hmp321_agpv4_chr5.vcf.gz \
$PWD/../data/uplifted_APGv4/hmp321_agpv4_chr6.vcf.gz \
$PWD/../data/uplifted_APGv4/hmp321_agpv4_chr7.vcf.gz \
$PWD/../data/uplifted_APGv4/hmp321_agpv4_chr8.vcf.gz \
$PWD/../data/uplifted_APGv4/hmp321_agpv4_chr9.vcf.gz \
$PWD/../data/uplifted_APGv4/hmp321_agpv4_chr10.vcf.gz --threads 15 \
| bcftools view -Oz -S $PWD/../data/01_metabolite_data/lines_with_phenolics_overlap_widiv_panzea_hmp_names.txt --threads 15 \
| bcftools view -Oz -m2 -M2 -v snps --threads 15 \
| bcftools view -Oz -q 0.05:minor --threads 15 \
| bcftools +setGT -Oz -- -t q -i 'GT="het"' -n "./." \
> $PWD/../data/Biallelic_MAF_filtered_Homozygous_hmp321_304g_agpv4_chr_all.vcf.gz

## Index the file
bcftools index ../data/Biallelic_MAF_filtered_Homozygous_hmp321_304g_agpv4_chr_all.vcf.gz --threads 20

## Download TASSEL
git clone https://bitbucket.org/tasseladmin/tassel-5-standalone.git

## Download data from Mazaheri et al. 2019 and unzip 
wget http://datadryad.org/api/v2/datasets/doi%253A10.5061%252Fdryad.n0m260p/download 
unzip download

## Make a temporary directory
mkdir tmp

## Convert .hmp.txt to VCF. Change the -Xms and -Xmx flags depending on your system hardware specs (RAM available min/max to be used by TASSEL)
perl tassel-5-standalone/run_pipeline.pl -Xms4G -Xmx64G \
-h widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.hmp.txt \
-export -exportType VCF

## Generate the site summary for SNPs in the file
perl tassel-5-standalone/run_pipeline.pl -Xms4G -Xmx64G \
-h widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.hmp.txt \
-GenotypeSummaryPlugin -endPlugin -export summary

## Move the .vcf file into the temporary directory
mv widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.vcf tmp

## Download B73 V4 reference sequence
wget https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz -P tmp
gunzip tmp/Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz
sed 's/Chr//g' tmp/Zm-B73-REFERENCE-GRAMENE-4.0.fa > tmp/no.chr.prefix.Zm-B73-REFERENCE-GRAMENE-4.0.fa

## Extract the reference allele in the vcf from the fasta file
(grep -v "#" tmp/widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.vcf \
| awk '{printf("%s:%s-%s\n",$1,$2,$2);}' | \
while read P; do samtools faidx tmp/no.chr.prefix.Zm-B73-REFERENCE-GRAMENE-4.0.fa ${P} ; done) \
> tmp/ref_allele_out.txt
sed '$!N;s/\n/ /' tmp/ref_allele_out.txt > tmp/ref_allele_out_wide.txt

## Fix the reference allele using R
Rscript 01.1_fix_REF_ALT_var.R

## Convert .hmp.txt with fixed ref to VCF
perl tassel-5-standalone/run_pipeline.pl -Xms4G -Xmx64G -h tmp/ref_fixed_widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.hmp.txt -export -exportType VCF

## Move the .vcf to the data directory
mv ref_fixed_widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.vcf $PWD/../data/

## bgzip and index the file
bgzip $PWD/../data/ref_fixed_widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.vcf
bcftools index $PWD/../data/ref_fixed_widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.vcf.gz --threads 20

## Subset the lines overlapping with WIDIV and Panzea and in metabolite profiling experiment
bcftools view -Oz -v snps $PWD/../data/ref_fixed_widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.vcf.gz \
-S $PWD/../data/01_metabolite_data/lines_with_phenolics_overlap_widiv_panzea_widiv_names.txt \
-o $PWD/../data/ref_fixed_widiv_304g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.vcf.gz

## Keep only biallelic SNPs, Filter for MAF > 0.05 and set heterozygotes to NA
bcftools view -Oz -m2 -M2 -v snps $PWD/../data/ref_fixed_widiv_304g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.vcf.gz \
| bcftools view -Oz -q 0.05:minor \
| bcftools +setGT -Oz -- -t q -i 'GT="het"' -n "./." \
| bcftools reheader --samples $PWD/../data/01_metabolite_data/widiv_panzea_name_converter.txt \
> $PWD/../data/Biallelic_MAF_filtered_Homozygous_widiv_304g_agpv4_chr_all.vcf.gz

## Index the file
bcftools index $PWD/../data/Biallelic_MAF_filtered_Homozygous_widiv_304g_agpv4_chr_all.vcf.gz --threads 20

## Merge the two vcf files 
bcftools concat --allow-overlaps --rm-dups all -Oz $PWD/../data/Biallelic_MAF_filtered_Homozygous_widiv_304g_agpv4_chr_all.vcf.gz \
$PWD/../data/Biallelic_MAF_filtered_Homozygous_hmp321_304g_agpv4_chr_all.vcf.gz --threads 20 \
> $PWD/../data/PANZEA_WIDIV_Biallelic_MAF_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz

## Do filtering with concated SNP data to make sure all SNPs are filtered equally
bcftools view -Oz -v snps $PWD/../data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz --threads 28 \
| bcftools view -Oz -m2 -M2 --threads 20 \
| bcftools +fill-tags -Oz --threads 20 \
| bcftools +setGT -Oz -- -t q -i 'GT="het"' -n "./." \
| bcftools view -i 'F_MISSING<0.1' --threads 20 \
| bcftools view -Oz -q 0.05:minor --threads 20 \
> $PWD/../data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_Missing_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz

## Make a SNPstats file 
bcftools stats $PWD/../data/02_genotypic_data/PANZEA_WIDIV_Biallelic_MAF_Missing_filtered_Homozygous_304g_agpv4_chr_all.vcf.gz > $PWD/../data/02_genotypic_data/combined_snp_file.stats

# ## Clean up the working directories
# mkdir $PWD/../data/02_genotypic_data
# rm -rf tmp/
# rm -rf tassel-5-standalone/
# rm -rf ../data/uplifted_APGv4/
# rm download \
# 62biomAP_v_B73_SNPMatrix.txt.gz \
# README_for_* \
# summary* \
# widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.hmp.txt
# mv $PWD/../data/combined_snp_file.stats \
# $PWD/../data/*.vcf* $PWD/../data/02_genotypic_data







