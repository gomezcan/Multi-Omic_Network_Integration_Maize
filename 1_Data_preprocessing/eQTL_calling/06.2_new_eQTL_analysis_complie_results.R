# =============================================================================
# EDITORIAL NOTE (added 2026-08-18 when archiving for publication; not by the
# original author). RUN THIS BEFORE the _liberal sibling: the shell block below
# is what creates trans_results.txt, which the liberal script's (commented) awk
# expects to already exist. This script is the STRICT branch and feeds the eGRN
# cis-eQTL set; the liberal branch feeds the GAN trans-eQTL set.
# Two hazards when re-running: the shell commands are inline and uncommented, so
# source() this file only if that is what you want; and the appends (>>) plus the
# directory glob make it non-idempotent -- clear prior outputs between runs.
# =============================================================================
## Deleting any files from trash
# find /home/<chtc-user>/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/09_Seedling_V1_304genos/my_results -maxdepth 1 -name '*' | xargs -P 10 rm -r

## After copying from condor make another copy on server and untar then combine to a single file
library(data.table)
library(tidyverse)

# DS9 ---------------------------------------------------------------------
#### Dataset number 9 (v1 seedling) do this for cis and trans files separately to prevent errors
cd /home/<chtc-user>/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/09_Seedling_V1_304genos/
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'cis*' --to-command="tail -n +2" \; >> cis_results.txt
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'trans*' --to-command="tail -n +2" \; >> trans_results.txt
awk '{ if ($6 >= 4) { print } }' trans_results.txt > trans_results_filtered.txt

setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/11_combined_results/")
header = fread("00_header.txt")
cis_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/09_Seedling_V1_304genos/cis_results.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos)
  # filter(n_models_support >=8)
length(unique(cis_tmp$gene_id))
fwrite(cis_tmp,"cis_seedling_v1_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/09_Seedling_V1_304genos/trans_results_filtered.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_seedling_v1_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS7 ---------------------------------------------------------------------
## Dataset number 7 (LMAD)
cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/07_LMAD_131genos/
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'cis*' --to-command="tail -n +2" \; >> cis_results.txt
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'trans*' --to-command="tail -n +2" \; >> trans_results.txt
awk '{ if ($6 >= 4) { print } }' trans_results.txt > trans_results_filtered.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/11_combined_results/")
header = fread("00_header.txt")
cis_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/07_LMAD_131genos/cis_results.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos)
# filter(n_models_support >=8)
length(unique(cis_tmp$gene_id))
fwrite(cis_tmp,"cis_LMAD_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/07_LMAD_131genos/trans_results_filtered.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_LMAD_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS6 ---------------------------------------------------------------------
## Dataset number 6 (Kernel)
cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/06_Kern_169genos/
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'cis*' --to-command="tail -n +2" \; >> cis_results.txt
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'trans*' --to-command="tail -n +2" \; >> trans_results.txt
awk '{ if ($6 >= 4) { print } }' trans_results.txt > trans_results_filtered.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/11_combined_results/")
header = fread("00_header.txt")
cis_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/06_Kern_169genos/cis_results.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos)
# filter(n_models_support >=8)
length(unique(cis_tmp$gene_id))
fwrite(cis_tmp,"cis_Kern_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/06_Kern_169genos/trans_results_filtered.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_Kern_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS5 ---------------------------------------------------------------------
## Dataset number 5 (Germinating shoot)
cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/05_Gshoot_178genos/
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'cis*' --to-command="tail -n +2" \; >> cis_results.txt
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'trans*' --to-command="tail -n +2" \; >> trans_results.txt
awk '{ if ($6 >= 4) { print } }' trans_results.txt > trans_results_filtered.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/11_combined_results/")
header = fread("00_header.txt")
cis_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/05_Gshoot_178genos/cis_results.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos)
# filter(n_models_support >=8)
length(unique(cis_tmp$gene_id))
fwrite(cis_tmp,"cis_Gshoot_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/05_Gshoot_178genos/trans_results_filtered.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_Gshoot_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS4 ---------------------------------------------------------------------
## Dataset number 4 (L3 base)
cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/04_L3Base_178genos/
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'cis*' --to-command="tail -n +2" \; >> cis_results.txt
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'trans*' --to-command="tail -n +2" \; >> trans_results.txt
awk '{ if ($6 >= 4) { print } }' trans_results.txt > trans_results_filtered.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/11_combined_results/")
header = fread("00_header.txt")
cis_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/04_L3Base_178genos/cis_results.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos)
# filter(n_models_support >=8)
length(unique(cis_tmp$gene_id))
fwrite(cis_tmp,"cis_L3Base_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/04_L3Base_178genos/trans_results_filtered.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_L3Base_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS3 ---------------------------------------------------------------------
## Dataset number 3 (L3 tip)
cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/03_L3Tip_180genos/
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'cis*' --to-command="tail -n +2" \; >> cis_results.txt
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'trans*' --to-command="tail -n +2" \; >> trans_results.txt
awk '{ if ($6 >= 4) { print } }' trans_results.txt > trans_results_filtered.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/11_combined_results/")
header = fread("00_header.txt")
cis_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/03_L3Tip_180genos/cis_results.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos)
# filter(n_models_support >=8)
length(unique(cis_tmp$gene_id))
fwrite(cis_tmp,"cis_L3Tip_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/03_L3Tip_180genos/trans_results_filtered.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_L3Tip_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS2 ---------------------------------------------------------------------
## Dataset number 2 (LMAN)
cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/02_LMAN_183genos/
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'cis*' --to-command="tail -n +2" \; >> cis_results.txt
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'trans*' --to-command="tail -n +2" \; >> trans_results.txt
awk '{ if ($6 >= 4) { print } }' trans_results.txt > trans_results_filtered.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/11_combined_results/")
header = fread("00_header.txt")
cis_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/02_LMAN_183genos/cis_results.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos)
# filter(n_models_support >=8)
length(unique(cis_tmp$gene_id))
fwrite(cis_tmp,"cis_LMAN_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/02_LMAN_183genos/trans_results_filtered.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_LMAN_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# DS1 ---------------------------------------------------------------------
## Dataset number 1 (GRoot)
cd ~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/01_GRoot_175genos/
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'cis*' --to-command="tail -n +2" \; >> cis_results.txt
find . -type f -name '*.tar.gz' -exec tar -xzf '{}' --wildcards --no-anchored 'trans*' --to-command="tail -n +2" \; >> trans_results.txt
awk '{ if ($6 >= 4) { print } }' trans_results.txt > trans_results_filtered.txt

rm(list=ls())
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/11_combined_results/")
header = fread("00_header.txt")
cis_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/01_GRoot_175genos/cis_results.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos)
# filter(n_models_support >=8)
length(unique(cis_tmp$gene_id))
fwrite(cis_tmp,"cis_GRoot_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

trans_tmp = fread("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/01_GRoot_175genos/trans_results_filtered.txt") %>% 
  setNames(., colnames(header)) %>% 
  arrange(snp_chrom, snp_pos) %>% 
  mutate(snp_up_down_stream = ifelse(snp_chrom==gene_chrom, snp_up_down_stream, NA))
# filter(n_models_support >=8)
length(unique(trans_tmp$gene_id))
fwrite(trans_tmp,"trans_GRoot_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

# Combine all ---------------------------------------------------------------------
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/11_combined_results/")
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
  
fwrite(x,"00_filtered_cis_trans_all_eQTL_results.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

