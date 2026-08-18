## Compute percent variance explained for eQTL results 

library(data.table)
library(tidyverse)

## Read in eQTL results
setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/12_combined_results_liberal/")
eqtls = fread("00_filtered_cis_trans_all_eQTL_results_liberal.txt") %>% 
  mutate(dataset = ifelse(dataset=="seedling","Seedling",dataset)) %>% 
  select(-c(gene_chrom,gene_TSS,gene_TTS,snp_distance_to_TSS,snp_up_down_stream))

## Read in gff annotation combined with syntenic information
gene_anno = fread("~/eqtl_pdi_coexpr/data/06_annotation_data/Zea_mays.B73_RefGen_v4.50.gff3.with.desc.synteny.txt")

## Join eQTL results with gene annotation and synteny then calculate distance of SNP to TSS of gene if on the same chromosome
eqtls = eqtls %>% 
  left_join(., gene_anno) %>% 
  mutate(snp_distance_to_TSS = ifelse(snp_chrom==gene_chrom, abs(snp_pos-gene_TSS), no = NA)) %>% 
  mutate(snp_up_down_stream = ifelse(is.na(snp_distance_to_TSS), NA,
                                     ifelse(snp_pos <= gene_TSS,paste0("upstream"),
                                                          ifelse(snp_pos >= gene_TSS & snp_pos <= gene_TTS,"within_gene_body",
                                                                 ifelse(snp_pos >= gene_TSS, paste0("downstream"),NA)))))
  
## Read in MAF tables for each dataset
setwd("~/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/")
maf_files = list.files(pattern = "*.afreq")
dsets = list()
for(i in maf_files){
  dsets[[i]] = fread(i) %>% 
    mutate(dataset = gsub(".afreq","", i)) %>% 
    mutate(dataset = gsub("^(?:[^_]+_){1}([^_]+).*", "\\1", dataset)) %>% 
    mutate(maf = pmin(ALT_FREQS, 1-(ALT_FREQS))) %>% 
    mutate(n_obs = OBS_CT/2) %>% 
    select(dataset,ID,maf, n_obs) %>% 
    setNames(., c("dataset","snp","maf","n_obs"))
}

## Loop through adding pve to a new dataframe
my_df = list()
for(i in seq_along(dsets)){
  
  ds_filter = unique(dsets[[i]]$dataset)
  
  my_df[[i]] = eqtls %>% 
    filter(dataset == all_of(ds_filter)) %>% 
    left_join(., dsets[[i]]) %>% 
    mutate(pve_naive_model = (2*(effect_naive_model^2)*maf*(1-maf))/(2*(effect_naive_model^2)*maf*(1-maf)+((std_err_naive_model)^2)*2*n_obs*maf*(1-maf)), .before = naive_model_pval)
  
}

my_df = bind_rows(my_df) %>% 
  relocate(snp_distance_to_TSS, .before = gene_type) %>% 
  relocate(snp_up_down_stream, .before = gene_type) %>% 
  relocate(maf, .before = n_models_support) %>% 
  relocate(n_obs, .before = snp)
colnames(my_df)[8] = "snp_maf"

setwd("~/eqtl_pdi_coexpr/results/05_eQTL_results_genome_wide/12_combined_results_liberal/")
fwrite(my_df,"00_filtered_cis_trans_all_eQTL_results_with_pve_and_synteny_liberal.txt",sep = "\t", row.names = FALSE, col.names = TRUE)


