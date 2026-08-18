# =============================================================================
# EDITORIAL NOTE (added 2026-08-18 when archiving for publication; not by the
# original author). RUN THIS AFTER 06.2_new_eQTL_analysis_complie_results.R,
# which produces the trans_results.txt this script's (commented) awk consumes.
# This is the LIBERAL branch: trans associations only, gated at
# n_models_support >= 2, feeding the GAN trans-eQTL layer.
# =============================================================================
## Deleting any files from trash
# find /home/<chtc-user>/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/09_Seedling_V1_304genos/my_results -maxdepth 1 -name '*' | xargs -P 10 rm -r

## After copying from condor make another copy on server and untar then combine to a single file
library(data.table)
library(tidyverse)

# DS9 ---------------------------------------------------------------------
#### Dataset number 9 (v1 seedling) do this for cis and trans files separately to prevent errors
# cd /home/<chtc-user>/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/09_Seedling_V1_304genos/
# awk '{ if ($6 >= 2) { print } }' trans_results.txt > trans_results_filtered_liberal.txt

setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/12_combined_results_liberal/")
header = fread("00_header.txt")

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/09_Seedling_V1_304genos/trans_results_filtered_liberal.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_seedling_v1_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS7 ---------------------------------------------------------------------
## Dataset number 7 (LMAD)
# cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/07_LMAD_131genos/
# awk '{ if ($6 >= 2) { print } }' trans_results.txt > trans_results_filtered_liberal.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/12_combined_results_liberal/")
header = fread("00_header.txt")

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/07_LMAD_131genos/trans_results_filtered_liberal.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_LMAD_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS6 ---------------------------------------------------------------------
## Dataset number 6 (Kernel)
# cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/06_Kern_169genos/
# awk '{ if ($6 >= 2) { print } }' trans_results.txt > trans_results_filtered_liberal.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/12_combined_results_liberal/")
header = fread("00_header.txt")

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/06_Kern_169genos/trans_results_filtered_liberal.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_Kern_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS5 ---------------------------------------------------------------------
## Dataset number 5 (Germinating shoot)
# cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/05_Gshoot_178genos/
# awk '{ if ($6 >= 2) { print } }' trans_results.txt > trans_results_filtered_liberal.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/12_combined_results_liberal/")
header = fread("00_header.txt")

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/05_Gshoot_178genos/trans_results_filtered_liberal.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_Gshoot_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS4 ---------------------------------------------------------------------
## Dataset number 4 (L3 base)
# cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/04_L3Base_178genos/
# awk '{ if ($6 >= 2) { print } }' trans_results.txt > trans_results_filtered_liberal.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/12_combined_results_liberal/")
header = fread("00_header.txt")

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/04_L3Base_178genos/trans_results_filtered_liberal.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_L3Base_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS3 ---------------------------------------------------------------------
## Dataset number 3 (L3 tip)
# cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/03_L3Tip_180genos/
# awk '{ if ($6 >= 2) { print } }' trans_results.txt > trans_results_filtered_liberal.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/12_combined_results_liberal/")
header = fread("00_header.txt")

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/03_L3Tip_180genos/trans_results_filtered_liberal.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_L3Tip_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS2 ---------------------------------------------------------------------
## Dataset number 2 (LMAN)
# cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/02_LMAN_183genos/
# awk '{ if ($6 >= 2) { print } }' trans_results.txt > trans_results_filtered_liberal.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/12_combined_results_liberal/")
header = fread("00_header.txt")

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/02_LMAN_183genos/trans_results_filtered_liberal.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_LMAN_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS1 ---------------------------------------------------------------------
## Dataset number 1 (GRoot)
# cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/01_GRoot_175genos/
# awk '{ if ($6 >= 2) { print } }' trans_results.txt > trans_results_filtered_liberal.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/12_combined_results_liberal/")
header = fread("00_header.txt")

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/01_GRoot_175genos/trans_results_filtered_liberal.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_GRoot_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# Combine all ---------------------------------------------------------------------
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/12_combined_results_liberal/")
my_files = list.files()
my_files = my_files[!my_files%in%"00_header.txt"]
my_tissues = gsub("^(?:[^_]+_){1}([^_]+).*", "\\1", my_files)

# my_files = my_files[1:2]
# my_tissues = my_tissues[1:2]

x = list()
for(i in seq_along(my_files)){
  x[[i]] = fread(my_files[[i]]) %>% 
    mutate(dataset = my_tissues[i], .before = 1)
}

x = bind_rows(x) %>% 
  distinct(.) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  filter(n_models_support>=2)
  
fwrite(x,"00_filtered_cis_trans_all_eQTL_results_liberal.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

