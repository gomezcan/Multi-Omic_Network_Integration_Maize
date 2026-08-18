## Fix REF/ALT variants in WIDIV (Mazaheri et al. 2019) Hapmap formatted file 

## Jonas Rodriguez ##
## 03-29-2021 ##

## Libraries
library(data.table)
library(tidyverse)

## Read REF allele extracted from the fasta file
r_vars = fread("tmp/ref_allele_out_wide.txt", header = FALSE)
colnames(r_vars) = c("query","REF")
r_vars = r_vars %>% 
  mutate(pos = str_replace(r_vars$query, pattern = ".*:([^-]*)-.*", replacement = "\\1")) %>% 
  mutate(chrom = str_replace(r_vars$query, pattern = ".*>([^-]*):.*", replacement = "\\1")) %>% 
  select(chrom,pos,REF)
r_vars$pos = as.integer(r_vars$pos)

## Read in Hapmap file from Mazaheri et al. 2019
pub_dat = fread("widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.hmp.txt")
pub_dat$chrom = as.character(pub_dat$chrom)

## Determine alleles for each SNP using Summary information generated from TASSEL
a_vars = fread("summary3.txt")
a_vars = a_vars %>% 
  select("Chromosome","Physical Position", "Major Allele","Minor Allele")
colnames(a_vars) = c("chrom","pos","major","minor")
a_vars$chrom = as.character(a_vars$chrom)
a_vars$pos = as.integer(a_vars$pos)

## Join the reference alleles and Major/Minor
j_vars = full_join(a_vars,r_vars)

## set REF/ALT 
j_vars$ALT = ifelse(j_vars$REF==j_vars$major,j_vars$minor, j_vars$major) 
j_vars[is.na(j_vars)] <- "N"

## Keep common column names from j_vars df to join with Hapmap file
j_vars = j_vars %>% 
  mutate(allele = paste0(REF,"/",ALT)) %>% 
  select(chrom,pos,allele)
  
## Join the dataframes
df_joined = left_join(j_vars, pub_dat, by = c("chrom","pos")) %>% 
  mutate(alleles = allele) %>% 
  select(-allele) %>% 
  select(rs, alleles,chrom,pos, everything())

## Keep the necessary columns and write out the fixed hapmap file
fwrite(df_joined, file  = "tmp/ref_fixed_widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.hmp.txt", col.names = TRUE, sep = "\t", quote = FALSE)
